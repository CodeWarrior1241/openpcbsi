"""G8 3D regions on openEMS (#9-#16, #18), after Backplane's em/m3d.py and its G8 recommendations,
generalised to geometry cut from any KiCad board.

    python3 -m openpcbsi region BOARD SPEC.toml [-c CONFIG]

    [region]
    name = "rx0_caps"
    nets = ["A_P", "A_N", "B_P", "B_N"]     # the nets that carry ports where they cross the box
    box = [x0, y0, x1, y1]                  # mm, board coordinates
    fmax_ghz = 20
    lossy = false                           # dispersive-free lossy dielectric (tan d) and copper conductivity
    max_timesteps = 400000
    end_criteria = 1e-4
    lead_mm = 2.0                           # synthetic straight lead-in at each port
    preflight_only = false

The model (built here in KiCad's python, run by em3d.py in the openEMS python):
- every copper polygon in the box on every layer (tracks with their thickness, pads, zone fills with
  holes and islands), vias as cylinders with their pads, solder mask and air above and below;
- two-pad parts between region nets (value from the KiCad value field): resistors as lumped resistors;
  a capacitor whose charging time through the ports exceeds the run (any AC-coupling capacitor) as its
  metal body bridging the pads, since in band it is a short; a small capacitor or an inductor as a
  series R-L-C lumped element;
- at each crossing of a port net with the box edge, the track is extended straight out of the box by
  `lead_mm` (zero-thickness, as openEMS's microstrip port expects) and fed by an MSLPort beyond the
  absorbing layer; P/N pairs on one side are driven +1/-1 in one run (mixed-mode parameters);
- PML on the four sides, MUR top and bottom with 1 mm of air; mesh lines on copper edges by the 1/3-2/3
  rule, >= 5 cells across the narrowest track and gap, graded outwards, lines closer than 5 um merged;
- pre-flight refuses a model with copper running along the absorbing layer, or ports shorter than the
  PML plus the feed-to-measurement distance;
- two grids; each run records whether it ended by the energy criterion or by the time-step limit; the
  result is accepted only if both ended by the criterion and agree on line impedance and TDR;
- the line impedance of the model is read from the TDR (renormalised to the port's measured impedance)
  and compared with the 2D solver for the same geometry (FS2); beyond 2 % the TDR is reported relative
  to the model's own line.
"""
import hashlib
import json
import math
import os
import re
import subprocess
import time

import numpy as np

from . import tdr as T

VENV_PY = os.environ.get("OPENPCBSI_OPENEMS_PYTHON", "/media/fpgadev/Dev_Tools/OpenEMS/install/venv/bin/python")
C0 = 299792458.0
# Grid-to-grid tolerances. The absolute impedance level of a run moves with the mesh by up to ~2.5 %
# (seen on three grids of one pair: 98.95, 101.09, 99.9 ohm against 99.6 ohm in 2D, the box-to-lead ratio
# constant to 0.1 %), so the line impedance may differ by 3 % between grids, the TDR's shape by 2 %, and
# the summed power of a column of S by 5 % above 1 before a result counts as not passive.
Z_GRID_TOL_PCT = 3.0
SHAPE_TOL_PCT = 2.0
PASSIVE_TOL = 1.05
T_RUN_S = 0.5e-9      # a series C whose 100 ohm time constant exceeds this is left out of the FDTD model
PML_CELLS = 8
AIR_UM = 1000.0


# ---- geometry helpers (um) --------------------------------------------------------------------
def clip_ring(ring, box):
    """Sutherland-Hodgman clip of a polygon against an axis-aligned box (x0, y0, x1, y1)."""
    x0, y0, x1, y1 = box
    pts = [tuple(p) for p in ring]
    for inside, inter in (
            (lambda p: p[0] >= x0, lambda a, b: (x0, a[1] + (b[1] - a[1]) * (x0 - a[0]) / (b[0] - a[0]))),
            (lambda p: p[0] <= x1, lambda a, b: (x1, a[1] + (b[1] - a[1]) * (x1 - a[0]) / (b[0] - a[0]))),
            (lambda p: p[1] >= y0, lambda a, b: (a[0] + (b[0] - a[0]) * (y0 - a[1]) / (b[1] - a[1]), y0)),
            (lambda p: p[1] <= y1, lambda a, b: (a[0] + (b[0] - a[0]) * (y1 - a[1]) / (b[1] - a[1]), y1))):
        out = []
        n = len(pts)
        for i in range(n):
            a, b = pts[i - 1], pts[i]
            ia, ib = inside(a), inside(b)
            if ib:
                if not ia:
                    out.append(inter(a, b))
                out.append(b)
            elif ia:
                out.append(inter(a, b))
        pts = out
        if not pts:
            return []
    return pts if len(pts) >= 3 else []



def grow_ring(pts, d):
    """Polygon grown outwards by d (vertex miter offset; pads and tracks are convex or nearly so)."""
    n = len(pts)
    if n < 3 or d <= 0:
        return pts
    area = sum(pts[i][0] * pts[(i + 1) % n][1] - pts[(i + 1) % n][0] * pts[i][1] for i in range(n))
    sgn = 1.0 if area > 0 else -1.0             # counter-clockwise: the outward normal of (dx, dy) is (dy, -dx)
    out = []
    for i in range(n):
        (xa, ya), (xb, yb), (xc, yc) = pts[i - 1], pts[i], pts[(i + 1) % n]
        n1 = (yb - ya, -(xb - xa))
        n2 = (yc - yb, -(xc - xb))
        l1, l2 = math.hypot(*n1) or 1.0, math.hypot(*n2) or 1.0
        n1 = (sgn * n1[0] / l1, sgn * n1[1] / l1)
        n2 = (sgn * n2[0] / l2, sgn * n2[1] / l2)
        bx, by = n1[0] + n2[0], n1[1] + n2[1]
        c = (bx * n1[0] + by * n1[1])
        k = d / c if c > 1e-3 else d
        out.append((xb + bx * k, yb + by * k))
    return out

def capsule(x1, y1, x2, y2, w, n=6):
    """Track segment as a polygon: rectangle with polygonal round caps."""
    dx, dy = x2 - x1, y2 - y1
    L = math.hypot(dx, dy)
    a = math.atan2(dy, dx) if L > 0 else 0.0
    r = w / 2
    pts = []
    for k in range(n + 1):
        t = a + math.pi / 2 + math.pi * k / n
        pts.append((x2 + r * math.cos(t - math.pi), y2 + r * math.sin(t - math.pi)))
    for k in range(n + 1):
        t = a - math.pi / 2 + math.pi * k / n
        pts.append((x1 + r * math.cos(t + math.pi), y1 + r * math.sin(t + math.pi)))
    return pts


def pip(x, y, ring):
    inside = False
    n = len(ring)
    for i in range(n):
        (xa, ya), (xb, yb) = ring[i - 1], ring[i]
        if (ya > y) != (yb > y) and x < xa + (y - ya) * (xb - xa) / (yb - ya):
            inside = not inside
    return inside


def seg_cross(x1, y1, x2, y2, side, c):
    """Crossing of a segment with the line x = c (side 'x') or y = c (side 'y'): parameter t or None."""
    a, b = (x1, x2) if side == "x" else (y1, y2)
    if (a - c) * (b - c) > 0 or a == b:
        return None
    return (c - a) / (b - a)


# ---- mesh --------------------------------------------------------------------------------------
def mesh_axis(lo, hi, fixed, zones, hmax, ratio=1.3, min_space=5.0):
    """Lines from lo to hi by a sweep: the cell at x is at most target(x), where target is h inside each
    of `zones` [(a, b, h)] and grows by `ratio` per cell away from it, up to hmax. `fixed` lines are
    kept (lines closer than min_space to them are dropped)."""
    zs = [(a, b, h) for a, b, h in zones if b >= lo and a <= hi]
    # every fixed line is a small zone of its own: the sweep grades towards it instead of jumping
    fl = sorted(set(round(f, 3) for f in fixed if lo < f < hi))
    for i, f in enumerate(fl):
        gaps = [abs(f - g) for g in (fl[i - 1] if i else None, fl[i + 1] if i + 1 < len(fl) else None) if g is not None]
        hf = max(min_space, min(gaps + [hmax]))
        zs.append((f, f, hf))

    def target(x):
        t = hmax
        for a, b, h in zs:
            d = max(a - x, 0.0, x - b)
            # geometric growth: n cells of h * ratio^k cover h (ratio^n - 1)/(ratio - 1) of distance
            t = min(t, h + (ratio - 1.0) * d)
        return max(t, min_space)
    fx = sorted(set(round(f, 3) for f in fixed if lo < f < hi))
    fset = set(fx) | {round(lo, 3), round(hi, 3)}
    lines, x, k = [lo], lo, 0
    while x < hi - 1e-6:
        s = target(x)
        for _ in range(4):                     # do not step over a finer region ahead
            s = min(s, target(min(x + s, hi)))
        nxt = x + s
        while k < len(fx) and fx[k] <= x + 1e-6:
            k += 1
        if k < len(fx) and fx[k] <= nxt + 1e-6:
            nxt = fx[k]
            k += 1
        if hi - nxt < 0.5 * target(hi):
            nxt = hi
        nxt = min(nxt, hi)
        if nxt - lines[-1] >= min_space:
            lines.append(nxt)
        elif (round(nxt, 3) in fset or nxt == hi) and len(lines) > 1 and round(lines[-1], 3) not in fset:
            lines[-1] = nxt                     # a fixed line (a copper sheet, a port) always wins
        elif round(nxt, 3) in fset and round(lines[-1], 3) in fset:
            lines.append(nxt)                   # two fixed lines closer than min_space: keep both
        x = nxt
    # smooth: no cell more than 1.5 x a neighbour (split the larger one)
    for _ in range(6):
        d = np.diff(lines)
        bad = [i for i in range(len(d)) if (i > 0 and d[i] > 2.0 * d[i - 1]) or (i < len(d) - 1 and d[i] > 2.0 * d[i + 1])]
        bad = [i for i in bad if d[i] / 2 >= min_space]
        if not bad:
            break
        for i in sorted(bad, reverse=True):
            lines.insert(i + 1, (lines[i] + lines[i + 1]) / 2)
    return [round(v, 3) for v in lines]


def mesh_axis_sym(lo, hi, c, fixed, zones, hmax, **kw):
    """mesh_axis mirrored about c (the centre of a differential pair): P and N see the same cells, so the
    mesh adds no mode conversion of its own. Lines beyond the mirrored span follow the domain edge."""
    R = max(hi - c, c - lo)
    fr = sorted({round(f if f >= c else 2 * c - f, 3) for f in fixed if lo <= f <= hi})
    zr = []
    for a, b_, h in zones:
        for aa, bb in ((a, b_), (2 * c - b_, 2 * c - a)):
            if bb >= c:
                zr.append((max(aa, c), bb, h))
    half = mesh_axis(c, c + R, fr, zr, hmax, **kw)
    full = sorted({round(2 * c - v, 3) for v in half} | set(half))
    full = [v for v in full if lo <= v <= hi]
    for v in (lo, hi):
        if not full or min(abs(v - x) for x in full) > 1e-3:
            full.append(round(v, 3))
    return sorted(full)


def thirds(edges, h):
    """1/3-2/3 rule: for an edge at e with metal on side s (+1 / -1), lines at e + s h/3 and e - s 2h/3."""
    out = []
    for e, s in edges:
        out += [e + s * h / 3.0, e - s * 2.0 * h / 3.0]
    return out


# ---- the model ---------------------------------------------------------------------------------
class Region:
    def __init__(self, ctx, spec, res_scale=1.0):
        self.ctx, self.spec, self.b = ctx, spec, ctx.b
        self.res_scale = res_scale
        self.um = 1000.0
        bx = [v * self.um for v in spec["box"]]
        self.box = (min(bx[0], bx[2]), min(bx[1], bx[3]), max(bx[0], bx[2]), max(bx[1], bx[3]))
        self.fmax = float(spec.get("fmax_ghz", 20.0)) * 1e9
        self.lead = float(spec.get("lead_mm", 2.0)) * self.um
        self.nets = list(spec["nets"])
        self.problems = []

    # z of each copper sheet (um, z up from the bottom of the board)
    def stack(self):
        st = self.b.stack
        th = st.thickness * self.um
        Z = {}
        # inner copper is a sheet: put it on the face towards the region's signal layers (a plane under
        # F.Cu on its top face), so the strip-to-plane spacing is the real dielectric thickness (the
        # middle of the copper made it half a copper thickness too far: ~4 % high impedance)
        sig = [L for L in getattr(self, "signal_layers", [])]
        sig_idx = [st.index[L] for L in sig if L in st.index]
        for i, c in enumerate(st.copper):
            if i == 0:
                Z[c["name"]] = th - c["z_bot"] * self.um                    # F.Cu: on the dielectric below
            elif i == len(st.copper) - 1:
                Z[c["name"]] = th - c["z_top"] * self.um                    # B.Cu: under the dielectric above
            elif sig_idx and all(k < i for k in sig_idx):
                Z[c["name"]] = th - c["z_top"] * self.um                    # signals above: the top face
            elif sig_idx and all(k > i for k in sig_idx):
                Z[c["name"]] = th - c["z_bot"] * self.um                    # signals below: the bottom face
            else:
                Z[c["name"]] = th - (c["z_top"] + c["z_bot"]) / 2 * self.um
        diel = [{"z0": th - d["z_bot"] * self.um, "z1": th - d["z_top"] * self.um, "er": d["er"], "tand": d["tand"]}
                for d in st.dielectrics]
        # inner copper is a sheet: the dielectrics on either side meet at it (no empty copper layer)
        inner = [Z[c["name"]] for c in st.copper[1:-1]]
        for d in diel:
            for zs in inner:
                if 0 <= zs - d["z1"] < 60:
                    d["z1"] = zs
                if 0 <= d["z0"] - zs < 60:
                    d["z0"] = zs
        mask = []
        for side, cu in (("F", st.copper[0]["name"]), ("B", st.copper[-1]["name"])):
            m = st.mask.get(side)
            if m and m.get("t"):
                t = m["t"] * self.um
                z = Z[cu]
                mask.append({"z0": z, "z1": z + t, "er": m["er"]} if side == "F" else {"z0": z - t, "z1": z, "er": m["er"]})
        thick = {c["name"]: c["t"] * self.um for c in st.copper}
        return Z, diel, mask, thick

    def build(self):
        b, ctx = self.b, self.ctx
        um = self.um
        x0, y0, x1, y1 = self.box
        # signal layers of the region's nets (for the inner sheet faces)
        self.signal_layers = sorted({b.layers[int(b.tracks["layer"][k])] for n in self.nets if n in b.net_id
                                     for k in b.tracks_of(n)})
        Z, diel, mask, thick = self.stack()
        top_cu, bot_cu = b.stack.copper[0]["name"], b.stack.copper[-1]["name"]
        # conformal mask over thick outer copper: it covers the sidewalls too (t_mask out from each edge),
        # as the 2D solver has it; the grown outline goes with the sheet and its edges onto the mesh
        mask_t = {} if self.spec.get("lossy") else {
            ln: next((m["z1"] - m["z0"] for m in mask if (m["z0"] >= Z[ln] - 1e-3) == (ln == top_cu)), 0.0)
            for ln in (top_cu, bot_cu) if thick.get(ln, 0) > 0}
        mask_edges_x, mask_edges_y = [], []
        L = b.layers
        T_ = b.tracks
        # ---- port crossings of the region nets --------------------------------------------------
        ports = []
        for net in self.nets:
            if net not in b.net_id:
                self.problems.append(f"net {net} not on the board")
                continue
            for k in b.tracks_of(net):
                xa, ya, xb, yb = (T_["x1"][k] * um, T_["y1"][k] * um, T_["x2"][k] * um, T_["y2"][k] * um)
                for side, c in (("xmin", x0), ("xmax", x1), ("ymin", y0), ("ymax", y1)):
                    ax = "x" if side[0] == "x" else "y"
                    t = seg_cross(xa, ya, xb, yb, ax, c)
                    if t is None:
                        continue
                    px, py = xa + t * (xb - xa), ya + t * (yb - ya)
                    lo, hi = (y0, y1) if ax == "x" else (x0, x1)
                    pv = py if ax == "x" else px
                    if not (lo < pv < hi):
                        continue
                    w_ = float(T_["w"][k]) * um
                    # the straight run of this track inside the box from the crossing: to the segment's
                    # inner end (less an allowance for the pad or via there) or across the whole box
                    inner = [(ex, ey) for ex, ey in ((xa, ya), (xb, yb)) if x0 < ex < x1 and y0 < ey < y1]
                    if inner:
                        run = math.hypot(inner[0][0] - px, inner[0][1] - py) - max(2 * w_, 300.0)
                    else:
                        run = (x1 - x0) if ax == "x" else (y1 - y0)
                    ports.append({"net": net, "side": side, "x": px, "y": py, "layer": L[int(T_["layer"][k])],
                                  "w": w_, "run_um": round(max(run, 0.0), 1)})
        # one port per (net, side): the crossing nearest the middle of the side
        uniq = {}
        for p in ports:
            k2 = (p["net"], p["side"], p["layer"])
            mid = ((y0 + y1) / 2) if p["side"][0] == "x" else ((x0 + x1) / 2)
            d = abs((p["y"] if p["side"][0] == "x" else p["x"]) - mid)
            if k2 not in uniq or d < uniq[k2][0]:
                uniq[k2] = (d, p)
        ports = [v[1] for v in uniq.values()]
        sides = sorted({p["side"] for p in ports})
        # domain: the box, extended by the lead on every side that carries a port
        dom = [x0 - (self.lead if "xmin" in sides else 0), y0 - (self.lead if "ymin" in sides else 0),
               x1 + (self.lead if "xmax" in sides else 0), y1 + (self.lead if "ymax" in sides else 0)]
        # reference plane of each port: the nearest layer whose copper (any zone fill) covers the point
        fills = {}
        for z in b.snap["zones"]:
            if z.get("net") and z.get("polys"):
                fills.setdefault(z["layer"], []).append(z)
        for p in ports:
            best = None
            for ln, zs in fills.items():
                if ln == p["layer"]:
                    continue
                for z in zs:
                    for poly in z["polys"]:
                        if pip(p["x"] / um, p["y"] / um, poly[0]) and not any(pip(p["x"] / um, p["y"] / um, h) for h in poly[1:]):
                            d = abs(Z[ln] - Z[p["layer"]])
                            if best is None or d < best[0]:
                                best = (d, ln)
            if best is None:
                self.problems.append(f"port {p['net']} {p['side']}: no reference plane under the port")
                p["ref"] = None
            else:
                p["ref"] = best[1]
        # ---- metal ------------------------------------------------------------------------------
        box = (x0, y0, x1, y1)
        sheets = []           # {"layer", "z", "pts", "metal", "prio", "thick"}
        sig_edges_x, sig_edges_y, fine_boxes = [], [], []
        w_min, gap_min = 1e9, 1e9
        for k in range(len(T_["x1"])):
            xa, ya, xb, yb = (T_["x1"][k] * um, T_["y1"][k] * um, T_["x2"][k] * um, T_["y2"][k] * um)
            w = float(T_["w"][k]) * um
            if max(xa, xb) + w < x0 or min(xa, xb) - w > x1 or max(ya, yb) + w < y0 or min(ya, yb) - w > y1:
                continue
            ln = L[int(T_["layer"][k])]
            # tracks along an axis: an exact rectangle (clean mesh edges); others: a capsule
            if abs(ya - yb) < 1e-6 or abs(xa - xb) < 1e-6:
                xa2, xb2, ya2, yb2 = min(xa, xb), max(xa, xb), min(ya, yb), max(ya, yb)
                if abs(ya - yb) < 1e-6:
                    ring = [(xa2, ya2 - w / 2), (xb2, ya2 - w / 2), (xb2, ya2 + w / 2), (xa2, ya2 + w / 2)]
                else:
                    ring = [(xa2 - w / 2, ya2), (xa2 + w / 2, ya2), (xa2 + w / 2, yb2), (xa2 - w / 2, yb2)]
            else:
                ring = capsule(xa, ya, xb, yb, w)
            cr = clip_ring(ring, box)
            if cr:
                sheets.append({"layer": ln, "z": Z[ln], "pts": cr, "metal": True, "prio": 50, "thick": thick[ln], "kind": "track"})
                tm = mask_t.get(ln, 0.0)
                if tm:
                    if abs(ya - yb) < 1e-6 or abs(xa - xb) < 1e-6:
                        gx = [p_[0] for p_ in ring]
                        gy = [p_[1] for p_ in ring]
                        gring = [(min(gx) - tm, min(gy) - tm), (max(gx) + tm, min(gy) - tm), (max(gx) + tm, max(gy) + tm), (min(gx) - tm, max(gy) + tm)]
                    else:
                        gring = capsule(xa, ya, xb, yb, w + 2 * tm)
                    gr = clip_ring(gring, box)
                    if gr:
                        sheets[-1]["mask_pts"] = gr
                        if abs(ya - yb) < 1e-6:
                            mask_edges_y += [min(p_[1] for p_ in gr), max(p_[1] for p_ in gr)]
                        elif abs(xa - xb) < 1e-6:
                            mask_edges_x += [min(p_[0] for p_ in gr), max(p_[0] for p_ in gr)]
                xs, ys = [p[0] for p in cr], [p[1] for p in cr]
                fine_boxes.append((min(xs), min(ys), max(xs), max(ys), w))
                w_min = min(w_min, w)
                if abs(ya - yb) < 1e-6:
                    sig_edges_y += [(min(ys), +1), (max(ys), -1)]
                elif abs(xa - xb) < 1e-6:
                    sig_edges_x += [(min(xs), +1), (max(xs), -1)]
        for p in b.snap["pads"]:
            for ln in p["layers"]:
                for poly in p.get("poly") or []:
                    ring = [(x * um, y * um) for x, y in poly[0]]
                    cr = clip_ring(ring, box)
                    if cr:
                        sheets.append({"layer": ln, "z": Z[ln], "pts": cr, "metal": True, "prio": 50, "thick": thick[ln], "kind": "pad"})
                        if mask_t.get(ln):
                            gr = clip_ring(grow_ring(ring, mask_t[ln]), box)
                            if gr:
                                sheets[-1]["mask_pts"] = gr
                        xs, ys = [q[0] for q in cr], [q[1] for q in cr]
                        fine_boxes.append((min(xs), min(ys), max(xs), max(ys), min(max(xs) - min(xs), max(ys) - min(ys))))
                        sig_edges_x += [(min(xs), +1), (max(xs), -1)]
                        sig_edges_y += [(min(ys), +1), (max(ys), -1)]
        for ln, zs in fills.items():
            outlines = []
            for z in zs:
                for poly in z["polys"]:
                    outlines.append(poly)
            for poly in outlines:
                ring0 = poly[0]
                depth = 0
                for other in outlines:
                    if other is poly:
                        continue
                    for h in other[1:]:
                        if pip(ring0[0][0], ring0[0][1], h):
                            depth += 1
                # planes run on under the port leads to the domain edge (into the PML): the ports'
                # reference; signal copper stays clipped to the box
                cr = clip_ring([(x * um, y * um) for x, y in ring0], tuple(dom))
                if not cr:
                    continue
                sheets.append({"layer": ln, "z": Z[ln], "pts": cr, "metal": True, "prio": 10 + 2 * depth, "thick": 0.0, "kind": "plane"})
                for h in poly[1:]:
                    ch = clip_ring([(x * um, y * um) for x, y in h], tuple(dom))
                    if ch:
                        sheets.append({"layer": ln, "z": Z[ln], "pts": ch, "metal": False, "prio": 11 + 2 * depth, "thick": 0.0, "kind": "hole"})
        # #19 variants: the standard fix applied to the model
        var = self.spec.get("variant") or {}
        names = b.layers
        for ref in var.get("void_under", []):
            for p in b.pads_by_ref.get(ref, []):
                if not p["layers"]:
                    continue
                li = b.L[p["layers"][0]]
                below = li + 1 if li < len(names) - 1 else li - 1          # the plane next to the pad's layer
                m = var.get("margin_mm", 0.15) * um
                hx, hy = p["sx"] * um / 2 + m, p["sy"] * um / 2 + m
                ring = [(p["x"] * um - hx, p["y"] * um - hy), (p["x"] * um + hx, p["y"] * um - hy),
                        (p["x"] * um + hx, p["y"] * um + hy), (p["x"] * um - hx, p["y"] * um + hy)]
                sheets.append({"layer": names[below], "z": Z[names[below]], "pts": ring, "metal": False, "prio": 40,
                               "thick": 0.0, "kind": "hole"})
        if var.get("fill_layer") in Z:
            fx0, fy0, fx1, fy1 = [v * um for v in var["fill_box"]]
            sheets.append({"layer": var["fill_layer"], "z": Z[var["fill_layer"]], "pts": [(fx0, fy0), (fx1, fy0), (fx1, fy1), (fx0, fy1)],
                           "metal": True, "prio": 45, "thick": 0.0, "kind": "plane"})
        self._extra_vias = [{"x": x * um, "y": y * um, "r": 75.0, "pad_r": 150.0, "z0": min(Z.values()), "z1": max(Z.values()),
                             "pads": sorted(set(Z.values())), "net": (self.ctx.cfg.reference_nets or ["GND"])[0]}
                            for x, y in var.get("add_vias", [])]
        # pair gaps: facing edges of different nets on one layer along the same axis
        for ax, edges in (("y", sig_edges_y), ("x", sig_edges_x)):
            ups = sorted(e for e, s in edges if s == -1)
            downs = sorted(e for e, s in edges if s == +1)
            for e in ups:
                nxt = [d for d in downs if d > e + 1.0]
                if nxt and nxt[0] - e < 500:
                    gap_min = min(gap_min, nxt[0] - e)
        # ---- vias ---------------------------------------------------------------------------------
        vias = []
        for v in b.snap["vias"]:
            x, y = v["x"] * um, v["y"] * um
            if not (x0 < x < x1 and y0 < y < y1):
                continue
            zs_ = [Z[l] for l in v["span"] if l in Z]
            vias.append({"x": x, "y": y, "r": v["drill"] * um / 2, "pad_r": v["d"] * um / 2, "z0": min(zs_), "z1": max(zs_),
                         "pads": [Z[l] for l in v["flashed"] if l in Z], "net": v["net"]})
        vias += getattr(self, "_extra_vias", [])
        # ---- series parts -------------------------------------------------------------------------
        from .pdn import parse_value
        lumped = []
        region_nets = set(self.nets)
        for f in b.snap["footprints"]:
            pads = [p for p in b.pads_by_ref.get(f["ref"], []) if p["net"]]
            if len(pads) != 2 or f.get("dnp"):
                continue
            if not ({pads[0]["net"], pads[1]["net"]} <= region_nets and pads[0]["net"] != pads[1]["net"]):
                continue
            (xa, ya), (xb, yb) = (pads[0]["x"] * um, pads[0]["y"] * um), (pads[1]["x"] * um, pads[1]["y"] * um)
            if not all(x0 < x < x1 and y0 < y < y1 for x, y in ((xa, ya), (xb, yb))):
                continue
            val = parse_value(f.get("value", ""))
            pre = re.match(r"[A-Za-z]+", f["ref"]).group(0).upper()
            unit = re.search(r"[0-9.]\s*[pnuµm]?\s*([FHR])", f.get("value", ""), re.I)
            if pre not in ("C", "L", "FB", "R"):  # an unusual prefix: the value's unit decides
                pre = {"F": "C", "H": "L", "R": "R"}.get(unit.group(1).upper(), pre) if unit else pre
            if val is None:
                self.problems.append(f"{f['ref']} value '{f.get('value')}' not understood: modelled as 0 ohm")
            ax = "x" if abs(xb - xa) >= abs(yb - ya) else "y"
            sx, sy = pads[0]["sx"] * um, pads[0]["sy"] * um
            if ax == "x":
                ia, ib = sorted([xa, xb])
                start = [ia + sx / 2, (ya + yb) / 2 - min(sx, sy) / 2]
                stop = [ib - sx / 2, (ya + yb) / 2 + min(sx, sy) / 2]
            else:
                ia, ib = sorted([ya, yb])
                start = [(xa + xb) / 2 - min(sx, sy) / 2, ia + sy / 2]
                stop = [(xa + xb) / 2 + min(sx, sy) / 2, ib - sy / 2]
            ln = pads[0]["layers"][0]
            el = {"ref": f["ref"], "value": f.get("value"), "dir": ax, "start": start, "stop": stop, "z": Z[ln]}
            if pre == "C":
                c = val or 1e-6
                el.update(C=c, R=0.03)
                # A series C blocks the excitation's DC: it charges through both ports' 50 ohm with
                # tau = 100 ohm * C, far longer than the run for any coupling capacitor, so the energy
                # never decays. When tau exceeds the run, C is left out of the FDTD model (series R-L
                # only); its reactance at the lowest reported frequency says whether that matters.
                if 100.0 * c > T_RUN_S:
                    f_lo = self.fmax / 40
                    xc = 1 / (2 * math.pi * f_lo * c)
                    el.pop("C")
                    # in band it is a short: model the part as what it is there, a metal body bridging its
                    # pads (pad centre to pad centre, pad wide, ~half the pad pitch tall: 0201 0.3 mm, 0402
                    # 0.5 mm, 0603 0.8 mm), so the mounting inductance comes from the geometry. (openEMS's
                    # lumped series R-L element forces a tiny time step and gave non-passive results.)
                    pitch = math.hypot(xb - xa, yb - ya)
                    wdt = (sy if ax == "x" else sx)
                    if ax == "x":
                        b0 = [min(xa, xb), (ya + yb) / 2 - wdt / 2]
                        b1 = [max(xa, xb), (ya + yb) / 2 + wdt / 2]
                    else:
                        b0 = [(xa + xb) / 2 - wdt / 2, min(ya, yb)]
                        b1 = [(xa + xb) / 2 + wdt / 2, max(ya, yb)]
                    h = min(1000.0, max(300.0, 0.5 * pitch))
                    up_ = ln == b.stack.copper[0]["name"]
                    zc0 = Z[ln] + thick.get(ln, 0.0) if up_ else Z[ln] - thick.get(ln, 0.0)
                    el.update(kind="body", body={"start": b0, "stop": b1, "z0": zc0 if up_ else zc0 - h, "z1": zc0 + h if up_ else zc0})
                    el.pop("R", None)
                    el.pop("L", None)
                    el["C_omitted"] = {"C": c, "xc_ohm_at_f_lo": round(xc, 4), "f_lo_ghz": round(f_lo / 1e9, 3)}
                    if xc > 0.05 * 50:
                        self.problems.append(f"{f['ref']} {f.get('value')}: modelled as a short (a series C this large "
                                             f"does not settle in the run): its reactance, {xc:.1f} ohm at {f_lo / 1e9:.2f} GHz, "
                                             "is missing from the result at the low end of the band")
            elif pre in ("L", "FB"):
                el.update(L=val if pre == "L" and val else 1e-9, R=0.05)
            else:
                el.update(R=max(val or 0.0, 1e-3))      # a plain lumped resistor (no series L: see above)
            lumped.append(el)
            fine_boxes.append((min(start[0], stop[0]), min(start[1], stop[1]), max(start[0], stop[0]), max(start[1], stop[1]), 50.0))
        # ---- leads and ports ------------------------------------------------------------------------
        leads = []
        for p in ports:
            w = p["w"]
            if p["side"] == "xmin":
                a, c = (dom[0], p["y"] - w / 2), (x0, p["y"] + w / 2)
            elif p["side"] == "xmax":
                a, c = (x1, p["y"] - w / 2), (dom[2], p["y"] + w / 2)
            elif p["side"] == "ymin":
                a, c = (p["x"] - w / 2, dom[1]), (p["x"] + w / 2, y0)
            else:
                a, c = (p["x"] - w / 2, y1), (p["x"] + w / 2, dom[3])
            p["lead"] = [a[0], a[1], c[0], c[1]]
            p["z"] = Z[p["layer"]]
            p["z_ref"] = Z[p["ref"]] if p.get("ref") else None
            leads.append(p)
        # pairs and excitation
        part = self.ctx.partner
        order = []
        for side in sides:
            here = [p for p in ports if p["side"] == side]
            done = set()
            for p in sorted(here, key=lambda q: q["net"]):
                if p["net"] in done:
                    continue
                q = next((r for r in here if r["net"] == part.get(p["net"])), None)
                if q:
                    pp, nn = sorted([p, q], key=lambda r: (not r["net"].endswith(("P", "+")), r["net"]))
                    order.append({"kind": "pair", "side": side, "ports": [pp, nn]})
                    done |= {p["net"], q["net"]}
                else:
                    order.append({"kind": "se", "side": side, "ports": [p]})
                    done.add(p["net"])
        exc_side = sides[0] if sides else None
        k = 1
        for grp in order:
            for i, p in enumerate(grp["ports"]):
                p["n"] = k
                k += 1
                if grp["side"] == exc_side and grp is next(g for g in order if g["side"] == exc_side):
                    common = self.spec.get("excite_mode") == "common"
                    p["excite"] = (1.0 if (i == 0 or common) else -1.0) if grp["kind"] == "pair" else 1.0
                else:
                    p["excite"] = 0.0
        # ---- mesh ---------------------------------------------------------------------------------
        er_max = max([d["er"] for d in diel] + [1.0])
        lam = C0 / self.fmax / math.sqrt(er_max) * 1e6
        hmax = min(lam / 20.0, 120.0)     # also keeps the PML (8 boundary cells) thin
        fine = float(self.spec["fine_um"]) if self.spec.get("fine_um") else \
            max(10.0, min(40.0, min(w_min, gap_min if gap_min < 1e8 else w_min) / 5.0))
        fine *= self.res_scale
        h_lead = min(50.0, hmax) * self.res_scale
        zx, zy = [], []
        for (a, b_, c, d, w) in fine_boxes:
            # fine across a straight track only; along it the cell follows the wavelength
            if (c - a) > 4 * (d - b_):
                zy.append((b_ - 3 * fine, d + 3 * fine, fine))
                zx.append((a, c, min(h_lead * 2, hmax)))
            elif (d - b_) > 4 * (c - a):
                zx.append((a - 3 * fine, c + 3 * fine, fine))
                zy.append((b_, d, min(h_lead * 2, hmax)))
            else:
                zx.append((a - 3 * fine, c + 3 * fine, fine))
                zy.append((b_ - 3 * fine, d + 3 * fine, fine))
        for p in leads:
            la = p["lead"]
            if p["side"][0] == "x":
                zx.append((min(la[0], la[2]), max(la[0], la[2]), h_lead))
            else:
                zy.append((min(la[1], la[3]), max(la[1], la[3]), h_lead))
        fx = thirds(sig_edges_x, fine) + mask_edges_x
        fy = thirds(sig_edges_y, fine) + mask_edges_y
        for v in vias:
            for r in (v["r"], v["pad_r"]):
                fx += [v["x"] - r, v["x"] + r]
                fy += [v["y"] - r, v["y"] + r]
            fx.append(v["x"])
            fy.append(v["y"])
            zx.append((v["x"] - v["pad_r"] - 3 * fine, v["x"] + v["pad_r"] + 3 * fine, fine * 1.5))
            zy.append((v["y"] - v["pad_r"] - 3 * fine, v["y"] + v["pad_r"] + 3 * fine, fine * 1.5))
        for el in lumped:
            bx_ = el.get("body") or el
            fx += [bx_["start"][0], bx_["stop"][0]]
            fy += [bx_["start"][1], bx_["stop"][1]]
        for p in leads:
            la = p["lead"]
            if p["side"][0] == "x":
                fy += thirds([(la[1], +1), (la[3], -1)], fine)
            else:
                fx += thirds([(la[0], +1), (la[2], -1)], fine)
        pair = next((g for g in order if g["kind"] == "pair"), None)
        if pair and pair["side"][0] == "x":           # pair along x: mirror the y mesh about its centre
            c = (pair["ports"][0]["y"] + pair["ports"][1]["y"]) / 2
            X = mesh_axis(dom[0], dom[2], fx, zx, hmax)
            Y = mesh_axis_sym(dom[1], dom[3], c, fy, zy, hmax)
        elif pair:                                    # pair along y: mirror the x mesh
            c = (pair["ports"][0]["x"] + pair["ports"][1]["x"]) / 2
            X = mesh_axis_sym(dom[0], dom[2], c, fx, zx, hmax)
            Y = mesh_axis(dom[1], dom[3], fy, zy, hmax)
        else:
            X = mesh_axis(dom[0], dom[2], fx, zx, hmax)
            Y = mesh_axis(dom[1], dom[3], fy, zy, hmax)
        zfix = sorted(set(Z.values()))
        zt = max(zfix)
        zb = min(zfix)
        zones_z = []
        z_sig = [Z[ln] for ln in self.signal_layers if ln in Z]
        for d in diel:
            zfix += [d["z0"], d["z1"]]
            # the dielectric between a region signal layer and its plane holds most of the line's field
            # (and sets its capacitance): finer there than elsewhere in the stack
            next_to_signal = any(min(abs(d["z0"] - z), abs(d["z1"] - z)) < 1.0 for z in z_sig)
            nd = float(self.spec.get("diel_cells_signal", 8) if next_to_signal else self.spec.get("diel_cells", 4))
            zones_z.append((d["z0"], d["z1"], max(5.0, (d["z1"] - d["z0"]) / nd) * self.res_scale))
        for m in mask:
            zfix += [m["z0"], m["z1"]]
        top_cu, bot_cu = b.stack.copper[0]["name"], b.stack.copper[-1]["name"]
        for ln, t in thick.items():           # faces of the thick outer copper (inner copper is a sheet: no ~9 um cells)
            if t > 0 and not self.spec.get("lossy") and ln in (top_cu, bot_cu):
                zfix += [Z[ln] + t] if ln == top_cu else [Z[ln] - t]
                mt = next((m["z1"] - m["z0"] for m in mask if (m["z0"] >= Z[ln] - 1e-3) == (ln == top_cu)), 0.0)
                if mt:                         # the mask's surface over the copper
                    zfix += [Z[ln] + t + mt] if ln == top_cu else [Z[ln] - t - mt]
        for el in lumped:
            if el.get("body"):
                zfix += [el["body"]["z0"], el["body"]["z1"]]
        Zl = mesh_axis(zb - AIR_UM, zt + AIR_UM, zfix, zones_z, hmax * 0.8, ratio=1.4, min_space=10.0)
        self.model = {"domain": dom, "box": list(box), "z_sheets": Z, "dielectrics": diel, "mask": mask,
                      "sheets": sheets, "vias": vias, "lumped": lumped, "ports": leads, "groups": order,
                      "mesh": {"x": X, "y": Y, "z": Zl}, "fmax": self.fmax, "lossy": bool(self.spec.get("lossy", False)),
                      "end_criteria": float(self.spec.get("end_criteria", 1e-4)),
                      "max_timesteps": int(self.spec.get("max_timesteps", 400000)),
                      "pml_cells": PML_CELLS, "fine_um": fine, "w_min_um": w_min, "gap_min_um": gap_min,
                      "excite_mode": self.spec.get("excite_mode", "odd"),
                      "problems": self.problems}
        return self.model

    # ---- checks on the model -----------------------------------------------------------------------
    def mesh_report(self):
        m = self.model
        X, Y = np.array(m["mesh"]["x"]), np.array(m["mesh"]["y"])
        def cells(lines, a, b):
            return int(((lines > a + 1e-6) & (lines < b - 1e-6)).sum()) + 1
        across_t, across_g = [], []
        for s in m["sheets"]:
            if s["kind"] != "track":
                continue
            xs, ys = [p[0] for p in s["pts"]], [p[1] for p in s["pts"]]
            wx, wy = max(xs) - min(xs), max(ys) - min(ys)
            if wy < wx:
                across_t.append(cells(Y, min(ys), max(ys)))
            else:
                across_t.append(cells(X, min(xs), max(xs)))
        tracks = [s for s in m["sheets"] if s["kind"] == "track"]
        for i, s in enumerate(tracks):
            ys = [p[1] for p in s["pts"]]
            xs = [p[0] for p in s["pts"]]
            for t in tracks[i + 1:]:
                if t["layer"] != s["layer"]:
                    continue
                ys2 = [p[1] for p in t["pts"]]
                xs2 = [p[0] for p in t["pts"]]
                if max(xs) - min(xs) > max(ys) - min(ys) and min(max(xs), max(xs2)) > max(min(xs), min(xs2)):
                    lo, hi = sorted([(max(ys), min(ys2)), (max(ys2), min(ys))], key=lambda q: abs(q[1] - q[0]))[0]
                    g = abs(hi - lo)
                    if 1 < g < 500:
                        across_g.append(cells(Y, min(lo, hi), max(lo, hi)))
        sp = np.concatenate([np.diff(X), np.diff(Y), np.diff(np.array(m["mesh"]["z"]))])
        return {"cells": int(len(X) * len(Y) * len(m["mesh"]["z"])), "nx": len(X), "ny": len(Y), "nz": len(m["mesh"]["z"]),
                "min_cells_across_track": min(across_t) if across_t else None,
                "min_cells_across_gap": min(across_g) if across_g else None,
                "min_spacing_um": round(float(sp.min()), 2), "fine_um": m["fine_um"]}

    def preflight(self):
        m = self.model
        X, Y = np.array(m["mesh"]["x"]), np.array(m["mesh"]["y"])
        dom = m["domain"]
        why = []
        pml = {"xmin": X[PML_CELLS] - X[0], "xmax": X[-1] - X[-1 - PML_CELLS],
               "ymin": Y[PML_CELLS] - Y[0], "ymax": Y[-1] - Y[-1 - PML_CELLS]}
        cell = {"xmin": X[1] - X[0], "xmax": X[-1] - X[-2], "ymin": Y[1] - Y[0], "ymax": Y[-1] - Y[-2]}
        keep = {s: pml[s] + 8 * cell[s] for s in pml}
        if not m["ports"]:
            why.append("no port: no region net crosses the box edge")
        for p in m["ports"]:
            if p.get("ref") is None:
                why.append(f"port {p['net']} ({p['side']}) has no reference plane")
            la = p["lead"]
            L = abs(la[2] - la[0]) if p["side"][0] == "x" else abs(la[3] - la[1])
            need = pml[p["side"]] + 4 * cell[p["side"]] + 300.0
            if L < need:
                why.append(f"port {p['net']} ({p['side']}): lead {L:.0f} um shorter than the PML plus the feed and measurement planes ({need:.0f} um)")
        # copper must cross a side (a line into the absorbing layer) or keep 8 cells clear of the PML
        for s in m["sheets"]:
            if not s["metal"] or s["kind"] == "plane":
                continue
            xs, ys = [p[0] for p in s["pts"]], [p[1] for p in s["pts"]]
            for side, inside, crosses in (
                    ("xmin", min(xs) - dom[0] < keep["xmin"], min(xs) <= m["box"][0] + 1e-3),
                    ("xmax", dom[2] - max(xs) < keep["xmax"], max(xs) >= m["box"][2] - 1e-3),
                    ("ymin", min(ys) - dom[1] < keep["ymin"], min(ys) <= m["box"][1] + 1e-3),
                    ("ymax", dom[3] - max(ys) < keep["ymax"], max(ys) >= m["box"][3] - 1e-3)):
                if not inside:
                    continue
                along = (max(ys) - min(ys)) if side[0] == "x" else (max(xs) - min(xs))
                if crosses and along < 2000.0:
                    continue
                why.append(f"{s['kind']} on {s['layer']} within {keep[side]:.0f} um of the {side} absorbing boundary")
                break
        for v in m["vias"]:
            for side, d in (("xmin", v["x"] - dom[0]), ("xmax", dom[2] - v["x"]), ("ymin", v["y"] - dom[1]), ("ymax", dom[3] - v["y"])):
                if d - v["pad_r"] < keep[side]:
                    why.append(f"via at ({v['x'] / 1000:.2f}, {v['y'] / 1000:.2f}) mm within {keep[side]:.0f} um of the {side} boundary")
        for g in m["groups"]:
            if g["kind"] == "pair" and self.spec.get("excite_mode") != "common":
                e = [p["excite"] for p in g["ports"]]
                if e[0] != -e[1]:
                    why.append(f"pair {g['ports'][0]['net']}/{g['ports'][1]['net']} not driven with opposite polarity")
        uniq = []
        for w in why:
            if w not in uniq:
                uniq.append(w)
        return {"ok": not uniq, "problems": uniq[:20], "pml_um": {k: round(v, 1) for k, v in pml.items()}}

    def path_connected(self):
        """Port-to-port connectivity of the model: the region nets' clipped tracks, joined through vias and
        series parts (an AC-coupled lane changes net at its capacitor). Each lane (a port and the ports it
        should reach: same polarity on the other sides) must end up in one piece."""
        b, um = self.b, self.um
        box = self.box
        T_ = b.tracks
        par = {}

        def find(a):
            while par.setdefault(a, a) != a:
                par[a] = par[par[a]]
                a = par[a]
            return a

        def union(a, c):
            par[find(a)] = find(c)
        pts = []                                   # (x, y, layer, node)
        nets = set(self.nets)
        for net in nets:
            if net not in b.net_id:
                continue
            for k in b.tracks_of(net):
                seg = _clip_seg(T_["x1"][k] * um, T_["y1"][k] * um, T_["x2"][k] * um, T_["y2"][k] * um, box)
                if seg is None:
                    continue
                l = b.layers[int(T_["layer"][k])]
                na, nb = ("t", k, 0), ("t", k, 1)
                union(na, nb)
                pts += [(seg[0], seg[1], l, na, float(T_["w"][k]) * um), (seg[2], seg[3], l, nb, float(T_["w"][k]) * um)]
        def near(x, y, l, r):
            return [nd for (px, py, pl, nd, w) in pts if (l is None or pl == l) and math.hypot(px - x, py - y) <= r + w / 2 + 5]
        for v in self.model["vias"]:
            ns = near(v["x"], v["y"], None, v["pad_r"])
            for nd in ns[1:]:
                union(ns[0], nd)
        for f in b.snap["footprints"]:
            pads = [p for p in b.pads_by_ref.get(f["ref"], []) if p["net"] in nets]
            for p in pads:                        # tracks ending on a pad join through the pad
                ns = near(p["x"] * um, p["y"] * um, None, max(p["sx"], p["sy"]) * um / 2)
                for nd in ns[1:]:
                    union(ns[0], nd)
        for e in self.model["lumped"]:            # series parts join their two pads' copper
            pads = [p for p in b.pads_by_ref.get(e["ref"], []) if p["net"]]
            groups = [near(p["x"] * um, p["y"] * um, None, max(p["sx"], p["sy"]) * um / 2) for p in pads]
            if len(groups) == 2 and groups[0] and groups[1]:
                union(groups[0][0], groups[1][0])
        node_of = {}
        for p in self.model["ports"]:
            ns = near(p["x"], p["y"], p["layer"], p["w"])
            node_of[p["n"]] = find(ns[0]) if ns else None
        res = {}
        groups = self.model["groups"]
        for g in groups:
            for i, p in enumerate(g["ports"]):
                others = [q for g2 in groups if g2 is not g and g2["kind"] == g["kind"] and g2["side"] != g["side"]
                          for j, q in enumerate(g2["ports"]) if j == i]
                key = f"{p['net']}@{p['side']}"
                res[key] = bool(others) and node_of.get(p["n"]) is not None and \
                    any(node_of.get(q["n"]) == node_of[p["n"]] for q in others)
        return (all(res.values()) if res else False), res


def _clip_seg(x1, y1, x2, y2, box):
    """Liang-Barsky clip of a segment to an axis-aligned box: (x1, y1, x2, y2) inside, or None."""
    x0, y0, xx, yy = box
    dx, dy = x2 - x1, y2 - y1
    t0, t1 = 0.0, 1.0
    for p, q in ((-dx, x1 - x0), (dx, xx - x1), (-dy, y1 - y0), (dy, yy - y1)):
        if p == 0:
            if q < 0:
                return None
            continue
        r = q / p
        if p < 0:
            t0 = max(t0, r)
        else:
            t1 = min(t1, r)
        if t0 > t1:
            return None
    return (x1 + t0 * dx, y1 + t0 * dy, x1 + t1 * dx, y1 + t1 * dy)


# ---- running ---------------------------------------------------------------------------------------
def _run_grid(model, od, gi, threads, log):
    gdir = os.path.join(od, f"grid{gi}")
    os.makedirs(gdir, exist_ok=True)
    mp = os.path.join(gdir, "model.json")
    json.dump(model, open(mp, "w"))
    runner = os.path.join(os.path.dirname(os.path.abspath(__file__)), "em3d.py")
    t0 = time.time()
    from . import jobs
    with jobs.slot("openems", log) as share:
        th = threads or share
        p = subprocess.run([VENV_PY, "-I", runner, mp, gdir, str(th)], capture_output=True, text=True,
                           env=dict(os.environ, OMP_NUM_THREADS=str(th)))
    logtxt = p.stdout + "\n" + p.stderr
    open(os.path.join(gdir, "runner.log"), "w").write(logtxt)
    if p.returncode != 0:
        raise RuntimeError(f"openEMS run failed (grid {gi}): {logtxt[-1500:]}")
    ports = json.load(open(os.path.join(gdir, "ports.json")))
    # openEMS prints "Max. number of timesteps: N" at the start of every run; only this warning means
    # the run stopped on the limit
    ended = "timesteps" if re.search(r"number of timesteps was reached", logtxt) else "end_criteria"
    log(f"  grid {gi}: {ports['cells']} cells, {time.time() - t0:.0f} s, ended by {ended}")
    return {"ports": ports, "ended_by": ended, "seconds": round(time.time() - t0, 1), "cells": ports["cells"]}


def _post(model, run):
    """S-parameters, line impedance and TDR from one grid's port waves."""
    P = {int(k): v for k, v in run["ports"]["ports"].items()}
    f = np.array(run["ports"]["freq_hz"])
    a = {n: np.array(v["a_re"]) + 1j * np.array(v["a_im"]) for n, v in P.items()}
    bb = {n: np.array(v["b_re"]) + 1j * np.array(v["b_im"]) for n, v in P.items()}
    zr = {n: np.array(v["zref_re"]) for n, v in P.items()}
    out = {"freq_ghz": (f / 1e9).round(4).tolist()}
    exc = [g for g in model["groups"] if any(p["excite"] for p in g["ports"])]
    if not exc:
        return out
    g = exc[0]
    other = [x for x in model["groups"] if x is not g]
    db = lambda s: (20 * np.log10(np.abs(s) + 1e-15)).round(3).tolist()
    band = (f > 0.05 * f.max()) & (f < 0.5 * f.max())
    if g["kind"] == "pair" and model.get("excite_mode") == "common":
        p1, n1 = g["ports"][0]["n"], g["ports"][1]["n"]
        ac = a[p1] + a[n1]
        scc11 = (bb[p1] + bb[n1]) / ac
        out.update(Scc11_db=db(scc11))
        tot = np.abs(scc11) ** 2
        far = next((x for x in other if x["kind"] == "pair"), None)
        if far:
            p2, n2 = far["ports"][0]["n"], far["ports"][1]["n"]
            scc21 = (bb[p2] + bb[n2]) / ac
            out.update(Scc21_db=db(scc21))
            tot = tot + np.abs(scc21) ** 2
        z_port = float(np.mean((1.0 / (1.0 / zr[p1] + 1.0 / zr[n1]))[band]))     # Zcomm = Zeven / 2
        s11, zref_s = scc11, 25.0
    elif g["kind"] == "pair":
        p1, n1 = g["ports"][0]["n"], g["ports"][1]["n"]
        ad = a[p1] - a[n1]
        sdd11, scd11 = (bb[p1] - bb[n1]) / ad, (bb[p1] + bb[n1]) / ad
        out.update(Sdd11_db=db(sdd11), Scd11_db=db(scd11))
        far = next((x for x in other if x["kind"] == "pair"), None)
        tot = np.abs(sdd11) ** 2 + np.abs(scd11) ** 2
        if far:
            p2, n2 = far["ports"][0]["n"], far["ports"][1]["n"]
            sdd21, scd21 = (bb[p2] - bb[n2]) / ad, (bb[p2] + bb[n2]) / ad
            out.update(Sdd21_db=db(sdd21), Scd21_db=db(scd21))
            tot = tot + np.abs(sdd21) ** 2 + np.abs(scd21) ** 2
        z_port = float(np.mean((zr[p1] + zr[n1])[band]))
        s11, zref_s = sdd11, 100.0
    else:
        p1 = g["ports"][0]["n"]
        s11 = bb[p1] / a[p1]
        out.update(S11_db=db(s11))
        tot = np.abs(s11) ** 2
        far = next((x for x in other if x["kind"] == "se"), None)
        if far:
            s21 = bb[far["ports"][0]["n"]] / a[p1]
            out.update(S21_db=db(s21))
            tot = tot + np.abs(s21) ** 2
        z_port = float(np.mean(zr[p1][band]))
        zref_s = 50.0
    out["passive"] = bool(np.max(tot[f > 0.02 * f.max()]) <= PASSIVE_TOL)
    out["z_port_ohm"] = round(z_port, 2)
    # TDR referenced to the port's own measured line impedance
    fu = np.linspace(0, f.max(), len(f))
    s_i = np.interp(fu, f, s11.real) + 1j * np.interp(fu, f, s11.imag)
    t, z = T.tdr(fu, T.renorm(s_i, z_port, zref_s), Z0=z_port, rise=tdr_rise_ps(f.max()) * 1e-12)
    out["tdr_t_ps"] = (t[t < 1.5e-9] * 1e12).round(2).tolist()
    out["tdr_z_ohm"] = z[t < 1.5e-9].round(3).tolist()
    return out


def tdr_rise_ps(fmax):
    return max(35e-12, 0.7 / fmax) * 1e12


def _line_z(model, post):
    """Impedance of the model's line: the TDR plateau on the uniform run of the excited port's track,
    from the box edge to its first discontinuity (pad, via, part, end), clear of both ends by 0.6 rise times.
    None when that run is too short to hold a plateau at this rise time (a part or via close to the box
    edge): the region then has no line reading to calibrate, rather than one taken on the discontinuity."""
    t, z = np.array(post.get("tdr_t_ps", [])), np.array(post.get("tdr_z_ohm", []))
    if not len(t):
        return None
    box = model["box"]
    p = next(q for q in model["ports"] if q["excite"])
    L = p.get("run_um") or ((box[2] - box[0]) if p["side"][0] == "x" else (box[3] - box[1]))
    ee = model.get("eps_eff_2d") or 3.0
    tl = 2 * L * 1e-6 * math.sqrt(ee) / C0 * 1e12            # round trip along the run, ps
    rise = tdr_rise_ps(model["fmax"])
    lo, hi = 0.6 * rise, tl - 0.6 * rise                      # the edge is centred on the discontinuity
    if hi - lo < max(20.0, 0.25 * tl):
        return None
    m = (t > lo) & (t < hi)
    return round(float(np.median(z[m])), 2) if m.any() else None


def run(ctx, spec, cache_dir, runs_dir, log=print):
    """Build, check and run a region; returns the result dict (also written to <dir>/result.json)."""
    t0 = time.time()
    reg = Region(ctx, spec, 1.0)
    model = reg.build()
    pf = reg.preflight()
    mesh = reg.mesh_report()
    conn, conn_by = reg.path_connected()
    key = hashlib.sha256(json.dumps({"m": model, "s": spec}, sort_keys=True, default=str).encode()).hexdigest()[:16]
    od = os.path.join(runs_dir, "em", f"{spec.get('name', 'region')}_{key}")
    os.makedirs(od, exist_ok=True)
    res = {"name": spec.get("name"), "board": ctx.b.path, "box_mm": spec["box"], "nets": reg.nets,
           "ports_kind": "line", "boundaries": ["PML_8"] * 4 + ["MUR", "MUR"], "preflight": pf, "mesh": mesh,
           "path_connected": conn, "path_by_net": conn_by, "series_parts": [{k: v for k, v in e.items() if k in ("ref", "value", "C", "R", "L", "C_omitted", "kind", "body")}
                                                                             for e in model["lumped"]],
           "problems": model["problems"], "lossy": model["lossy"], "dir": od}
    # 2D reference for the line at the excited port (FS2)
    try:
        res.update(_two_d(ctx, model))
        model["eps_eff_2d"] = res.get("eps_eff_2d")
    except Exception as e:
        res["two_d_error"] = str(e)
    if spec.get("preflight_only") or not pf["ok"]:
        res["ran"] = False
        json.dump(res, open(os.path.join(od, "result.json"), "w"), indent=1)
        return res
    threads = int(spec.get("threads", os.environ.get("OPENPCBSI_EM_THREADS", 0)))   # 0: the job slot's share
    runs = []
    for gi, scale in enumerate(spec.get("res_scales", [1.0, 0.7])):
        r2 = Region(ctx, spec, scale)
        m2 = r2.build()
        m2["eps_eff_2d"] = model.get("eps_eff_2d")
        rr = _run_grid(m2, od, gi, threads, log)
        rr["post"] = _post(m2, rr)
        rr["z_line_ohm"] = _line_z(m2, rr["post"])
        rr["res_scale"] = scale
        runs.append(rr)
    fine = runs[-1]
    res.update({k: v for k, v in fine["post"].items()})
    res["runs_per_grid"] = 1
    res["z_line_ohm"] = fine["z_line_ohm"]
    pe = next((q for q in model["ports"] if q["excite"]), None)
    if pe:                                     # the box length along the excited port's axis
        bx = model["box"]
        res["line_length_mm"] = round(((bx[2] - bx[0]) if pe["side"][0] == "x" else (bx[3] - bx[1])) / 1000.0, 3)
    res["eps_eff"] = model.get("eps_eff_2d")
    # convergence (#15): both runs ended by the criterion, line impedance and TDR agree
    conv = {"runs": [{"res_scale": r["res_scale"], "cells": r["cells"], "seconds": r["seconds"], "ended_by": r["ended_by"],
                      "z_line_ohm": r["z_line_ohm"]} for r in runs]}
    za, zb = runs[0]["z_line_ohm"], runs[-1]["z_line_ohm"]
    if za and zb:
        conv["z_line_diff_pct"] = round(abs(za - zb) / zb * 100, 2)
    ta, tb = np.array(runs[0]["post"].get("tdr_z_ohm", [])), np.array(runs[-1]["post"].get("tdr_z_ohm", []))
    if len(ta) and len(tb):
        n = min(len(ta), len(tb))
        conv["tdr_diff_ohm"] = round(float(np.max(np.abs(ta[:n] - tb[:n]))), 2)
        # The TDR's shape, each grid's curve over its own port reading: whether the structures in the box
        # (vias, transitions, parts) are resolved. The absolute level of a whole curve moves with the mesh
        # by up to ~2.5 % (the zero-thickness port lead's own discretisation, common to every point of the
        # curve); that level is what the 2D calibration checks and corrects, not the grid comparison.
        za_p, zb_p = runs[0]["post"].get("z_port_ohm"), runs[-1]["post"].get("z_port_ohm")
        if za_p and zb_p:
            conv["tdr_shape_diff_pct"] = round(float(np.max(np.abs(ta[:n] / za_p - tb[:n] / zb_p))) * 100, 2)
    if any(r["ended_by"] != "end_criteria" for r in runs):
        conv["converged"], conv["reason"] = False, "a run hit the time-step limit before the energy criterion"
    elif conv.get("z_line_diff_pct", 0) > Z_GRID_TOL_PCT or conv.get("tdr_shape_diff_pct", 0) > SHAPE_TOL_PCT:
        conv["converged"], conv["reason"] = False, "the two grids disagree on line impedance or TDR shape"
    else:
        conv["converged"], conv["reason"] = True, "both runs settled; line impedance and TDR shape agree"
    res["convergence"] = conv
    # calibration (#16): the model's line against 2D for the same geometry
    if res.get("z_2d_ohm") and res.get("z_line_ohm"):
        ag = (res["z_line_ohm"] - res["z_2d_ohm"]) / res["z_2d_ohm"] * 100
        res["calibration"] = {"z_model_ohm": res["z_line_ohm"], "z_2d_ohm": res["z_2d_ohm"], "agree_pct": round(ag, 2),
                              "renormalised": abs(ag) > 2.0,
                              "note": "TDR reported relative to the model's own line" if abs(ag) > 2.0 else "within 2 % of 2D"}
    elif res.get("z_2d_ohm"):
        pe = next((q for q in model["ports"] if q["excite"]), {})
        res["calibration"] = {"z_model_ohm": None, "z_2d_ohm": res["z_2d_ohm"], "agree_pct": None, "renormalised": False,
                              "note": f"not calibrated: the excited port's track runs {pe.get('run_um', 0) / 1000:.1f} mm "
                                      "uniform inside the box, too short for a TDR plateau at this rise time "
                                      "(move the box edge further from the first part or via)"}
    if res.get("passive") is False:
        res.setdefault("problems", []).append("result is not passive")
    res["ran"] = True
    res["seconds"] = round(time.time() - t0, 1)
    json.dump(res, open(os.path.join(od, "result.json"), "w"), indent=1)
    return res


def _two_d(ctx, model):
    from . import xsec
    p = next((q for q in model["ports"] if q["excite"]), None)
    if p is None:
        return {}
    b = ctx.b
    li = b.L[p["layer"]]
    refs = [b.layers[r] for r in ctx.ref_layers(li)]
    grp = next(g for g in model["groups"] if p in g["ports"])
    xc = ctx.cfg["xsec"]
    if grp["kind"] == "pair":
        q = grp["ports"][1]
        gap = abs((q["y"] - p["y"]) if p["side"][0] == "x" else (q["x"] - p["x"])) - p["w"]
        sec = xsec.build_section(b.stack, p["layer"], refs, [p["w"] / 1000, q["w"] / 1000], [gap / 1000], xc, xc["mask_er"])
        r = xsec.solve_lines(sec, 5.0, (1e9,))["freqs"][0]
        if model.get("excite_mode") == "common":
            return {"z_2d_ohm": round(r["Zcm"], 2), "eps_eff_2d": round(r.get("eeff_even", r["eeff_odd"]), 3), "line_kind": "pair, common mode",
                    "z_odd_2d": round(r["Zodd"], 2), "z_even_2d": round(r["Zeven"], 2)}
        return {"z_2d_ohm": round(r["Zdiff"], 2), "eps_eff_2d": round(r["eeff_odd"], 3), "line_kind": "pair",
                "z_odd_2d": round(r["Zodd"], 2), "z_even_2d": round(r["Zeven"], 2)}
    sec = xsec.build_section(b.stack, p["layer"], refs, [p["w"] / 1000], [], xc, xc["mask_er"])
    r = xsec.solve_lines(sec, 5.0, (1e9,))["freqs"][0]
    return {"z_2d_ohm": round(r["Z0"], 2), "eps_eff_2d": round(r["eeff"], 3), "line_kind": "se"}
