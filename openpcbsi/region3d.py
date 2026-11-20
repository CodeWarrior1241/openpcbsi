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
    rotate_deg = 0                          # turn the board about rotate_about (default: the box centre) first
    include_mm = 1.0                        # other nets' copper beyond this from the region's copper is left out
    term_ohm = 50                           # load at a pad where a region net ends in the box (terminate = false: none)
    max_hours = 6                           # refused when predicted to take longer; stopped there on the clock
    pml_cells = 8                           # absorbing layer thickness, in cells
    third_grid = true                       # a third, finer grid when two straddle the TDR limit (time allowing)
    rise_ps / nyquist_ghz                   # the link's edge and Nyquist (default: the nets' class)
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
- a region net that ends on a part's pad in the box: a lumped port of term_ohm from the pad to the
  nearest reference plane under it (the receiver's or cable's load, and a measured port), not an open end;
- the nets beyond the region's series parts (and their partners) join the region; other nets' copper
  farther than include_mm from the region's copper is left out, nearer copper meshed coarser;
- PML on the four sides, MUR top and bottom with 1 mm of air; mesh lines on the region copper's edges by
  the 1/3-2/3 rule, >= 5 cells across the narrowest track and gap, graded outwards; no two lines closer
  than 5 um (10 um in z), after mirroring about a pair too (merge_lines);
- pre-flight refuses a model with copper ending in or running along the absorbing layer, ports shorter
  than the PML plus the feed-to-measurement distance, a lane that reaches no other port or termination,
  a pair cut by a box corner, a port off square; run() then refuses a model predicted to stop on the
  step limit or to take longer than max_hours (run_estimate), all before anything is solved;
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

OPENEMS_ENV = "OPENPCBSI_OPENEMS_PYTHON"
C0 = 299792458.0


def openems_python():
    """The Python of an environment with openEMS and CSXCAD, from OPENPCBSI_OPENEMS_PYTHON. There is no
    default: openpcbsi never runs an interpreter it was not pointed at."""
    from .snapshot import InputError
    py = os.environ.get(OPENEMS_ENV, "").strip()
    if not py:
        raise InputError(f"3D regions need openEMS: set {OPENEMS_ENV} to the Python of an environment that has "
                         "openEMS and CSXCAD installed (e.g. <venv>/bin/python)")
    if not os.path.isfile(py) or not os.access(py, os.X_OK):
        raise InputError(f"{OPENEMS_ENV}={py} is not an executable Python")
    return py
# Grid-to-grid tolerances. The absolute impedance level of a run moves with the mesh by up to ~2.5 %
# (seen on three grids of one pair: 98.95, 101.09, 99.9 ohm against 99.6 ohm in 2D, the box-to-lead ratio
# constant to 0.1 %), so the line impedance may differ by 3 % between grids, the TDR's shape by 2 %, and
# the summed power of a column of S by 5 % above 1 before a result counts as not passive.
CLEAN_PCT = 5.0       # both grids within this in the box, and under the reflection floor: clean
TDR_LIMIT_PCT = 10.0  # the TDR excursion a region is judged against (emplan.judge, FS9)
PASSIVE_TOL = 1.05
T_RUN_S = 0.5e-9      # a series C whose 100 ohm time constant exceeds this is left out of the FDTD model
PML_CELLS = 8
KEEP_H = 5.0          # pre-flight clearance from the absorbing boundary: at least this many field heights
AIR_UM = 1000.0
ADD_CLEAR_UM = 100.0  # an added (variant) via keeps this much from other nets' copper
COARSE = 2.5          # other nets' copper: cells this many times the region lines' fine cell
TERM_OHM = 50.0       # a region net that ends on a pad in the box: its load to the reference plane
INCLUDE_MM = 1.0      # other nets' copper farther than this from the region's copper is left out (see _inclusion)
MAX_PORT_ANGLE = 10.0  # degrees between a port's track and the side's normal


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
    """Point in polygon (even-odd), vectorised over the ring's edges (zone outlines have thousands)."""
    r = np.asarray(ring, dtype=float)
    if len(r) < 3:
        return False
    xb, yb = r[:, 0], r[:, 1]
    xa, ya = np.roll(xb, 1), np.roll(yb, 1)
    span = (ya > y) != (yb > y)
    with np.errstate(divide="ignore", invalid="ignore"):
        cross = span & (x < xa + (y - ya) * (xb - xa) / (yb - ya))
    return bool(np.count_nonzero(cross) % 2)


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

    ZA = np.array([z[0] for z in zs], float)
    ZB = np.array([z[1] for z in zs], float)
    ZH = np.array([z[2] for z in zs], float)

    def target(x):
        # geometric growth: n cells of h * ratio^k cover h (ratio^n - 1)/(ratio - 1) of distance
        t = hmax
        if len(ZH):
            d = np.maximum(np.maximum(ZA - x, 0.0), x - ZB)
            t = min(t, float(np.min(ZH + (ratio - 1.0) * d)))
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
    return merge_lines([round(v, 3) for v in lines], fset, min_space)


def merge_lines(lines, fixed, min_space):
    """No two lines closer than min_space. The FDTD time step follows the smallest cell (Courant), so one
    pair of lines 1 nm apart makes every step ~1000 times shorter for the whole model. Of two lines too
    close, a fixed line (a copper sheet's face, a port, the 1/3-2/3 lines of an edge) wins over a graded
    one; two fixed lines become one between them (each moves by under min_space / 2); the domain's
    ends always stay."""
    lines = sorted(lines)
    if len(lines) < 2:
        return lines
    lo, hi = lines[0], lines[-1]
    isfix = lambda v: round(v, 3) in fixed
    out = [lo]
    for v in lines[1:]:
        if v - out[-1] >= min_space - 1e-9:
            out.append(v)
            continue
        if len(out) == 1:                          # out[-1] is the domain's first line: it stays
            if v == hi:
                out.append(v)                      # a domain thinner than min_space: keep both ends
            continue
        if v == hi:
            out[-1] = v                            # the domain's last line replaces the one before it
        elif isfix(v) and not isfix(out[-1]):
            out[-1] = v
        elif isfix(v) and isfix(out[-1]):
            out[-1] = round((out[-1] + v) / 2, 3)
            fixed = set(fixed) | {round(out[-1], 3)}
        # else: a graded line too close to the last one is dropped
    return out


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
    # mirrored as offsets from c: the centre line is one line (mirroring the rounded values gave it twice,
    # 1 nm apart: a 1 nm cell, and a time step ~1000 times too short)
    off = [v - half[0] for v in half]
    full = sorted({round(c - d, 3) for d in off} | {round(c + d, 3) for d in off})
    full = [v for v in full if lo + 1e-3 < v < hi - 1e-3] + [round(lo, 3), round(hi, 3)]
    ms = kw.get("min_space", 5.0)
    fixed = {round(f, 3) for f in fixed} | {round(2 * c - f, 3) for f in fixed} | {round(c, 3)}
    return merge_lines(sorted(set(full)), fixed, ms)


def thirds(edges, h):
    """1/3-2/3 rule: for an edge at e with metal on side s (+1 / -1), lines at e + s h/3 and e - s 2h/3."""
    out = []
    for e, s in edges:
        out += [e + s * h / 3.0, e - s * 2.0 * h / 3.0]
    return out


# ---- the model ---------------------------------------------------------------------------------
class Rotated:
    """The board seen in a frame turned by `deg` (counter-clockwise in board coordinates) about (cx, cy)
    mm: tracks, pads (centres and outlines), vias and zone fills turned; everything else the board's own.
    A route at 45 degrees becomes a route along an axis, so a box side can cut it square and its port
    lead continues it straight (an axis-aligned box cuts a diagonal pair at an angle, with a bend at each
    lead and the legs' spacing wrong along the edge). Planes and other copper at 45 degrees in the turned
    frame follow the rectilinear mesh as staircases at its resolution."""

    def __init__(self, b, deg, cx, cy):
        self._b = b
        c, s_ = math.cos(math.radians(deg)), math.sin(math.radians(deg))

        def R(x, y):
            x, y = np.asarray(x, float) - cx, np.asarray(y, float) - cy
            return cx + c * x - s_ * y, cy + s_ * x + c * y

        def ring(r):
            if not r:
                return r
            X, Y = R([q[0] for q in r], [q[1] for q in r])
            return [(float(a), float(b_)) for a, b_ in zip(X, Y)]
        T = dict(b.tracks)
        T["x1"], T["y1"] = R(b.tracks["x1"], b.tracks["y1"])
        T["x2"], T["y2"] = R(b.tracks["x2"], b.tracks["y2"])
        self.tracks = T
        pads = []
        for p in b.snap["pads"]:
            x, y = R(p["x"], p["y"])
            q = dict(p, x=float(x), y=float(y))
            if p.get("poly"):
                q["poly"] = [[ring(r) for r in poly] for poly in p["poly"]]
            if abs(s_) > 0.99:                        # a quarter turn: the pad's sides swap
                q["sx"], q["sy"] = p["sy"], p["sx"]
            pads.append(q)
        vias = []
        for v in b.snap["vias"]:
            x, y = R(v["x"], v["y"])
            vias.append(dict(v, x=float(x), y=float(y)))
        zones = [dict(z, polys=[[ring(r) for r in poly] for poly in z.get("polys") or []]) for z in b.snap["zones"]]
        self.snap = dict(b.snap, pads=pads, vias=vias, zones=zones)
        self.pads_by_ref = {}
        for p in pads:
            self.pads_by_ref.setdefault(p["ref"], []).append(p)
        self.deg, self.about = deg, (cx, cy)
        self._R = R

    def point(self, x, y):
        """A board point (mm) in the turned frame."""
        X, Y = self._R(x, y)
        return float(X), float(Y)

    def __getattr__(self, k):
        return getattr(self._b, k)


class Region:
    def __init__(self, ctx, spec, res_scale=1.0):
        self.ctx, self.spec, self.b = ctx, spec, ctx.b
        if spec.get("rotate_deg"):
            # the box is given in the turned frame (see Rotated), about rotate_about or the box's centre
            cx, cy = spec.get("rotate_about") or ((spec["box"][0] + spec["box"][2]) / 2, (spec["box"][1] + spec["box"][3]) / 2)
            self.b = Rotated(ctx.b, float(spec["rotate_deg"]), cx, cy)
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

    def _inclusion(self, Z, diel):
        """(test, distance in um): the test says whether a piece of other nets' copper comes within the
        inclusion distance of the region's own copper (its nets' tracks, pads and vias in the box, and the
        series parts between them). Copper farther out is left out of the model: its coupling to the
        region's lines falls as 1 / (1 + (d/h)^2) (under 1 % beyond 10 field heights h), and the absorbing
        sides already stand for the board beyond the box with a matched continuation, so a far stub or via
        adds nothing but an end near the absorber. Planes (zone fills) always go in. The distance is
        `include_mm` (the spec's, else OPENPCBSI_EM_INCLUDE_MM, else INCLUDE_MM), at least 2 KEEP_H field
        heights; `include_mm = 0` keeps all copper (the model as before the inclusion distance)."""
        b, um = self.b, self.um
        x0, y0, x1, y1 = self.box
        inc = self.spec.get("include_mm", os.environ.get("OPENPCBSI_EM_INCLUDE_MM", INCLUDE_MM))
        if not float(inc):
            return None, None
        from .parse import use
        use()                                          # shapely, from the pinned parsers' site
        from shapely.geometry import LineString, Point
        from shapely.ops import unary_union
        from shapely.prepared import prep
        sig = [Z[l] for l in self.signal_layers if l in Z]
        fh = max([min([d["z1"] - d["z0"] for d in diel if min(abs(d["z0"] - z), abs(d["z1"] - z)) < 1.0] or [0.0]) for z in sig] or [0.0])
        D = max(float(inc) * um, 2 * KEEP_H * fh)
        nets = set(self.nets)
        geo = []
        T_ = b.tracks
        for n in nets:
            if n not in b.net_id:
                continue
            for k in b.tracks_of(n):
                cl = _clip_seg(T_["x1"][k] * um, T_["y1"][k] * um, T_["x2"][k] * um, T_["y2"][k] * um, self.box)
                if cl is not None:
                    geo.append(LineString([cl[:2], cl[2:]]).buffer(float(T_["w"][k]) * um / 2, 8))
        for p in b.snap["pads"]:
            if p["net"] in nets and x0 < p["x"] * um < x1 and y0 < p["y"] * um < y1:
                geo.append(Point(p["x"] * um, p["y"] * um).buffer(max(p["sx"], p["sy"]) * um / 2, 8))
        for v in b.snap["vias"]:
            if v["net"] in nets and x0 < v["x"] * um < x1 and y0 < v["y"] * um < y1:
                geo.append(Point(v["x"] * um, v["y"] * um).buffer(v["d"] * um / 2, 8))
        for f in b.snap["footprints"]:                 # series parts between region nets: their pads
            ps = [p for p in b.pads_by_ref.get(f["ref"], []) if p["net"]]
            if len(ps) == 2 and {ps[0]["net"], ps[1]["net"]} <= nets:
                for p in ps:
                    geo.append(Point(p["x"] * um, p["y"] * um).buffer(max(p["sx"], p["sy"]) * um / 2, 8))
        if not geo:
            return None, None
        zone = prep(unary_union(geo).buffer(D, 8))

        def test(kind, *a):
            if kind == "seg":
                (xa, ya, xb, yb), w = a
                return zone.intersects(LineString([(xa, ya), (xb, yb)]).buffer(w / 2, 4))
            if kind == "poly":
                return zone.intersects(Polygon(a[0]).buffer(0))
            x, y, r = a
            return zone.intersects(Point(x, y).buffer(r, 8))
        from shapely.geometry import Polygon
        return test, round(D, 1)

    def _series_nets(self):
        """The region's nets with the nets beyond its series parts in the box (an AC-coupling capacitor,
        a series resistor or inductor) and their pair partners: a lane runs through its series parts to
        a port or a termination, never into a capacitor's pad as if that were its end."""
        b, um = self.b, self.um
        x0, y0, x1, y1 = self.box
        ref = set(self.ctx.cfg.reference_nets or []) | {z["net"] for z in b.snap["zones"] if z.get("net") and z.get("polys")}
        nets, added = list(self.nets), []
        for _ in range(4):
            grew = False
            for f in b.snap["footprints"]:
                ps = [p for p in b.pads_by_ref.get(f["ref"], []) if p["net"]]
                if len(ps) != 2 or f.get("dnp") or not re.match(r"(C|R|L|FB)\d", f["ref"]):
                    continue
                if not all(x0 < p["x"] * um < x1 and y0 < p["y"] * um < y1 for p in ps):
                    continue
                a, c = ps[0]["net"], ps[1]["net"]
                for u, v in ((a, c), (c, a)):
                    if u in nets and v not in nets and v not in ref:
                        for n in (v, self.ctx.partner.get(v)):
                            if n and n not in nets and n in b.net_id:
                                nets.append(n)
                                added.append(n)
                                grew = True
            if not grew:
                break
        self.nets, self.added_nets = nets, added

    def build(self):
        b, ctx = self.b, self.ctx
        um = self.um
        if not hasattr(self, "added_nets"):
            self._series_nets()
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
                    sl = math.hypot(xb - xa, yb - ya)
                    ang = math.degrees(math.acos(min(1.0, abs((xb - xa) if ax == "x" else (yb - ya)) / sl))) if sl > 0 else 90.0
                    ports.append({"net": net, "side": side, "x": px, "y": py, "layer": L[int(T_["layer"][k])],
                                  "w": w_, "run_um": round(max(run, 0.0), 1), "angle_deg": round(ang, 1)})
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
        # ports the model cannot measure correctly: refused in the pre-flight, (side, message)
        self.port_problems = []
        part = self.ctx.partner
        by_net = {}
        for p in ports:
            by_net.setdefault(p["net"], set()).add(p["side"])
        for n, sd in sorted(by_net.items()):
            q = part.get(n)
            if q in by_net and n < q and sd != by_net[q]:
                odd = sorted(sd ^ by_net[q])
                self.port_problems.append((odd[0], f"pair {n}/{q} leaves the box on different sides ({n}: {', '.join(sorted(sd))}; "
                                                   f"{q}: {', '.join(sorted(by_net[q]))}): a box corner cuts the pair, so its legs "
                                                   "would get single-ended ports"))
        for p in ports:
            if p["angle_deg"] > MAX_PORT_ANGLE:
                self.port_problems.append((p["side"], f"port {p['net']} ({p['side']}): the track crosses the box edge at "
                                                      f"{p['angle_deg']:.0f} degrees to its normal; the straight port lead meets it "
                                                      "with a bend at the edge (a pair's legs also get the wrong spacing): the side "
                                                      "must cut the route where it runs square to it"))
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
        # ---- which copper goes in -------------------------------------------------------------------
        box = (x0, y0, x1, y1)
        keep_in, inc_um = self._inclusion(Z, diel)
        # the region's own copper on inner layers with its real thickness (a zero-thickness sheet reads a
        # stripline's impedance ~9 % high against 2D); other nets' inner copper stays a sheet
        thick_inner = bool(self.spec.get("thick_inner", True)) and not self.spec.get("lossy")
        left_out = {"tracks": 0, "pads": 0, "vias": 0}
        # ---- metal ------------------------------------------------------------------------------
        sheets = []           # {"layer", "z", "pts", "metal", "prio", "thick"}
        sig_edges_x, sig_edges_y, fine_boxes = [], [], []
        w_min, gap_min = 1e9, 1e9
        for k in range(len(T_["x1"])):
            xa, ya, xb, yb = (T_["x1"][k] * um, T_["y1"][k] * um, T_["x2"][k] * um, T_["y2"][k] * um)
            w = float(T_["w"][k]) * um
            if max(xa, xb) + w < x0 or min(xa, xb) - w > x1 or max(ya, yb) + w < y0 or min(ya, yb) - w > y1:
                continue
            ln = L[int(T_["layer"][k])]
            own = b.net_names[int(T_["net"][k])] in self.nets
            if keep_in is not None and b.net_names[int(T_["net"][k])] not in self.nets and \
                    not keep_in("seg", (xa, ya, xb, yb), w):
                left_out["tracks"] += 1
                continue
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
                # the uncut centreline goes with the sheet: the pre-flight judges how the track meets the
                # absorbing boundary (crosses it, ends near it, runs along it), which the clipped outline hides
                sheets.append({"layer": ln, "z": Z[ln], "pts": cr, "metal": True, "prio": 50, "thick": thick[ln], "kind": "track",
                               "net": b.net_names[int(T_["net"][k])], "seg": [round(v, 3) for v in (xa, ya, xb, yb)], "w": w,
                               "thick_inner": own and ln not in (top_cu, bot_cu) and thick_inner})
                tm = mask_t.get(ln, 0.0)
                if tm and own:
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
                if not own:                    # other nets' copper: resolved coarser (it couples to the
                    fine_boxes.append((min(xs), min(ys), max(xs), max(ys), w, COARSE)) # region's lines only weakly)
                    continue
                fine_boxes.append((min(xs), min(ys), max(xs), max(ys), w))
                w_min = min(w_min, w)
                if abs(ya - yb) < 1e-6:
                    sig_edges_y += [(min(ys), +1), (max(ys), -1)]
                elif abs(xa - xb) < 1e-6:
                    sig_edges_x += [(min(xs), +1), (max(xs), -1)]
        for p in b.snap["pads"]:
            if keep_in is not None and p["net"] not in self.nets and p.get("poly") and \
                    not any(keep_in("poly", [(x * um, y * um) for x, y in poly[0]]) for poly in p["poly"] if len(poly[0]) >= 3):
                if x0 < p["x"] * um < x1 and y0 < p["y"] * um < y1:
                    left_out["pads"] += 1
                continue
            for ln in p["layers"]:
                for poly in p.get("poly") or []:
                    ring = [(x * um, y * um) for x, y in poly[0]]
                    cr = clip_ring(ring, box)
                    if cr:
                        sheets.append({"layer": ln, "z": Z[ln], "pts": cr, "metal": True, "prio": 50, "thick": thick[ln], "kind": "pad", "net": p["net"]})
                        if mask_t.get(ln):
                            gr = clip_ring(grow_ring(ring, mask_t[ln]), box)
                            if gr:
                                sheets[-1]["mask_pts"] = gr
                        xs, ys = [q[0] for q in cr], [q[1] for q in cr]
                        if p["net"] not in self.nets:
                            fine_boxes.append((min(xs), min(ys), max(xs), max(ys), min(max(xs) - min(xs), max(ys) - min(ys)), COARSE))
                            continue
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
        # added vias (the return-via fix): a through via of the reference net, with a pad only where its own
        # net's plane is, and a clearance hole in every other net's plane (as the board's own vias have): a
        # plain cylinder would short each power plane it passes through to ground
        ref_net = (self.ctx.cfg.reference_nets or ["GND"])[0]
        self._extra_vias = []
        for v in var.get("add_vias", []):
            x, y = float(v[0]), float(v[1])
            d, drill = (float(v[2]), float(v[3])) if len(v) >= 4 else (0.3, 0.15)
            if isinstance(b, Rotated):
                x, y = b.point(x, y)
            pads = []
            for ln, zs in fills.items():
                inside_own = inside_other = False
                for z in zs:
                    for poly in z["polys"]:
                        if pip(x, y, poly[0]) and not any(pip(x, y, h) for h in poly[1:]):
                            if z["net"] == ref_net:
                                inside_own = True
                            else:
                                inside_other = True
                if inside_own:
                    pads.append(Z[ln])
                if inside_other:
                    ring = [(x * um + (d / 2 + 0.1) * um * math.cos(a), y * um + (d / 2 + 0.1) * um * math.sin(a))
                            for a in np.linspace(0, 2 * math.pi, 24, endpoint=False)]
                    sheets.append({"layer": ln, "z": Z[ln], "pts": ring, "metal": False, "prio": 40, "thick": 0.0, "kind": "hole"})
            self._extra_vias.append({"x": x * um, "y": y * um, "r": drill * um / 2, "pad_r": d * um / 2, "z0": min(Z.values()),
                                     "z1": max(Z.values()), "pads": sorted(set(pads)), "net": ref_net, "added": True})
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
            if keep_in is not None and v["net"] not in self.nets and not keep_in("pt", x, y, v["d"] * um / 2):
                left_out["vias"] += 1
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
            if pre not in ("C", "L", "FB", "R"):  # a connector or IC with two pads on region nets: no series part
                continue
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
        # ---- terminations ----------------------------------------------------------------------------
        # A region net that ends on a part's pad inside the box (a connector pin, a BGA ball, a receiver)
        # would be an open end in a board-only model: |S11| = 1 by construction. Its pad gets a lumped port
        # of `term_ohm` (default TERM_OHM) to the nearest reference plane under it: the receiver's or the
        # mated cable's load, and a measured second port. Two legs of a pair: 100 ohm differential and
        # 25 ohm common mode, which is how a GTH, LVDS or PCIe receiver terminates on die.
        terms = []
        r_t = float(self.spec.get("term_ohm", TERM_OHM))
        series_refs = {e["ref"] for e in lumped}
        if self.spec.get("terminate", True):
            for p in b.snap["pads"]:
                if p["net"] not in region_nets or p.get("ref") in series_refs or not p["layers"]:
                    continue
                px, py = p["x"] * um, p["y"] * um
                if not (x0 < px < x1 and y0 < py < y1):
                    continue
                lns = [l for l in p["layers"] if l in Z]
                on = [l for l in lns if any(b.layers[int(T_["layer"][k])] == l for k in b.tracks_of(p["net"])
                                            if math.hypot(T_["x1"][k] * um - px, T_["y1"][k] * um - py) < max(p["sx"], p["sy"]) * um
                                            or math.hypot(T_["x2"][k] * um - px, T_["y2"][k] * um - py) < max(p["sx"], p["sy"]) * um)]
                ln = (on or lns or [None])[0]
                if ln is None:
                    continue
                best = None
                for lr, zs in fills.items():
                    if lr == ln:
                        continue
                    for z in zs:
                        for poly in z["polys"]:
                            if pip(p["x"], p["y"], poly[0]) and not any(pip(p["x"], p["y"], h) for h in poly[1:]):
                                d = abs(Z[lr] - Z[ln])
                                if best is None or d < best[0]:
                                    best = (d, lr)
                if best is None:
                    self.problems.append(f"{p.get('ref')} pad {p.get('num')} ({p['net']}): no reference plane under it for its termination")
                    continue
                h = min(p["sx"], p["sy"]) * um / 2
                terms.append({"net": p["net"], "ref": p.get("ref"), "pad": p.get("num"), "x": px, "y": py, "layer": ln,
                              "ref_layer": best[1], "R": r_t, "w": 2 * h,
                              "start": [px - h, py - h, Z[best[1]]], "stop": [px + h, py + h, Z[ln]]})
                fine_boxes.append((px - h, py - h, px + h, py + h, 2 * h))
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
        for side, lyr in sorted({(p["side"], p["layer"]) for p in ports}):
            # one group per side and layer: a route that leaves on the side it came in, on another layer
            # after a via, has its two ports on one side
            here = [p for p in ports if p["side"] == side and p["layer"] == lyr]
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
        # the excited group: a pair when the region has one (never one leg of a pair on its own), on the
        # side whose lanes run longest into the box (the line reading and calibration need that run)
        paired = [g for g in order if g["kind"] == "pair"]
        cand = paired or [g for g in order if part.get(g["ports"][0]["net"]) not in by_net]
        if self.spec.get("excite_side"):
            cand = [g for g in cand if g["side"] == self.spec["excite_side"]] or cand
        exc = max(cand, key=lambda g: (min(p["run_um"] for p in g["ports"]), g["side"] == (sides[0] if sides else None))) if cand else None
        k = 1
        for grp in order:
            for i, p in enumerate(grp["ports"]):
                p["n"] = k
                k += 1
                if grp is exc:
                    common = self.spec.get("excite_mode") == "common"
                    p["excite"] = (1.0 if (i == 0 or common) else -1.0) if grp["kind"] == "pair" else 1.0
                else:
                    p["excite"] = 0.0
        # the terminations as groups of their own (side "pad"), numbered after the line ports
        done = set()
        for t in sorted(terms, key=lambda q: (q["ref"] or "", q["net"])):
            if id(t) in done:
                continue
            q = min((r for r in terms if r is not t and id(r) not in done and r["net"] == part.get(t["net"])),
                    key=lambda r: (r["ref"] != t["ref"], math.hypot(r["x"] - t["x"], r["y"] - t["y"])), default=None)
            grp_t = [t, q] if q else [t]
            grp_t.sort(key=lambda r: (not r["net"].endswith(("P", "+")), r["net"]))
            for r in grp_t:
                r["n"], r["excite"], r["side"] = k, 0.0, "pad"
                k += 1
                done.add(id(r))
            order.append({"kind": "pair" if q else "se", "side": "pad", "ports": grp_t})
        # ---- mesh ---------------------------------------------------------------------------------
        er_max = max([d["er"] for d in diel] + [1.0])
        lam = C0 / self.fmax / math.sqrt(er_max) * 1e6
        hmax = min(lam / 20.0, 120.0)     # also keeps the PML (8 boundary cells) thin
        fine = float(self.spec["fine_um"]) if self.spec.get("fine_um") else \
            max(10.0, min(40.0, min(w_min, gap_min if gap_min < 1e8 else w_min) / 5.0))
        fine *= self.res_scale
        h_lead = min(50.0, hmax) * self.res_scale
        zx, zy = [], []
        for fb in fine_boxes:
            a, b_, c, d, w = fb[:5]
            hf = min(fine * fb[5], hmax) if len(fb) > 5 else fine
            # fine across a straight track only; along it the cell follows the wavelength
            if (c - a) > 4 * (d - b_):
                zy.append((b_ - 3 * hf, d + 3 * hf, hf))
                zx.append((a, c, min(h_lead * 2, hmax)))
            elif (d - b_) > 4 * (c - a):
                zx.append((a - 3 * hf, c + 3 * hf, hf))
                zy.append((b_, d, min(h_lead * 2, hmax)))
            else:
                zx.append((a - 3 * hf, c + 3 * hf, hf))
                zy.append((b_ - 3 * hf, d + 3 * hf, hf))
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
            hv = fine * 1.5 if v.get("net") in self.nets else min(fine * COARSE, hmax)
            zx.append((v["x"] - v["pad_r"] - 3 * fine, v["x"] + v["pad_r"] + 3 * fine, hv))
            zy.append((v["y"] - v["pad_r"] - 3 * fine, v["y"] + v["pad_r"] + 3 * fine, hv))
        for el in lumped:
            bx_ = el.get("body") or el
            fx += [bx_["start"][0], bx_["stop"][0]]
            fy += [bx_["start"][1], bx_["stop"][1]]
        for t in terms:
            fx += [t["start"][0], t["stop"][0]]
            fy += [t["start"][1], t["stop"][1]]
        for p in leads:
            la = p["lead"]
            if p["side"][0] == "x":
                fy += thirds([(la[1], +1), (la[3], -1)], fine)
            else:
                fx += thirds([(la[0], +1), (la[2], -1)], fine)
        pair = exc if (exc and exc["kind"] == "pair") else next((g for g in order if g["kind"] == "pair" and g["side"] != "pad"), None)
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
        z_ms = 10.0
        for ln in {s_["layer"] for s_ in sheets if s_.get("thick_inner")}:
            t = thick.get(ln, 0.0)
            if t > 0:                         # both faces of the inner copper, about its sheet's z
                zfix += [Z[ln] - t / 2, Z[ln] + t / 2]
                z_ms = min(z_ms, 0.45 * t)
        Zl = mesh_axis(zb - AIR_UM, zt + AIR_UM, zfix, zones_z, hmax * 0.8, ratio=1.4, min_space=z_ms)
        # field height: the largest spacing from a layer with signal copper in the box to its nearest
        # neighbouring copper (the dielectric that holds its field); the pre-flight's clearance floor
        field_h = 0.0
        for zs in {s["z"] for s in sheets if s["metal"] and s["kind"] in ("track", "pad")}:
            adj = [d["z1"] - d["z0"] for d in diel if min(abs(d["z0"] - zs), abs(d["z1"] - zs)) < 1.0]
            if adj:
                field_h = max(field_h, min(adj))
        self.model = {"domain": dom, "box": list(box), "z_sheets": Z, "dielectrics": diel, "mask": mask,
                      "sheets": sheets, "vias": vias, "lumped": lumped, "ports": leads, "terms": terms, "groups": order,
                      "mesh": {"x": X, "y": Y, "z": Zl}, "fmax": self.fmax, "lossy": bool(self.spec.get("lossy", False)),
                      "end_criteria": float(self.spec.get("end_criteria", 1e-4)),
                      "max_timesteps": int(self.spec.get("max_timesteps", 400000)),
                      "pml_cells": int(self.spec.get("pml_cells", PML_CELLS)), "fine_um": fine, "w_min_um": w_min, "gap_min_um": gap_min,
                      "field_h_um": round(field_h, 2), "include_um": inc_um, "left_out": left_out,
                      "port_problems": self.port_problems, "nets": list(self.nets), "added_nets": self.added_nets,
                      "rotate_deg": float(self.spec.get("rotate_deg") or 0.0),
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
        """Refuse a model whose result the absorbing boundaries would spoil.

        Ports need a reference plane and a lead longer than the PML plus the feed and measurement planes.
        Copper near the PML (closer than `keep`: the PML plus 8 cells, at least KEEP_H field heights): the
        PML absorbs a wave that runs INTO it, and a line that passes through it is terminated there like a
        matched load, so a track
        that crosses the absorbing side cleanly is benign, whatever its net and at any angle steeper than
        30 degrees to the side (a 45-degree track included; the PML's absorption falls off only at grazing
        incidence). What spoils a result is copper that ENDS near or inside the PML (an open end, a pad, a
        via: its fringing field reaches into the absorber, which then loads or reflects it unlike the real
        board) and copper that RUNS ALONG the PML (a line parallel to the absorber couples to it over its
        whole length: an unphysical loss). So, per side, every track's copper within `keep` must belong to
        a chain of segments (one net, one layer, joined end to end inside the band) that reaches an
        absorbing side of the box, has no free end inside the band, has no segment within 30 degrees of
        parallel to the side for more than 2 cells, and spans along the side no more than a 30-degree
        crossing would. Pads keep the rule they had (a pad is accepted only where it is cut by the box edge
        and spans under 2 mm along it); vias within `keep` are refused."""
        m = self.model
        X, Y = np.array(m["mesh"]["x"]), np.array(m["mesh"]["y"])
        dom = m["domain"]
        box = m["box"]
        why = []                                   # (side or None, message)
        nP = int(m.get("pml_cells", PML_CELLS))
        pml = {"xmin": X[nP] - X[0], "xmax": X[-1] - X[-1 - nP],
               "ymin": Y[nP] - Y[0], "ymax": Y[-1] - Y[-1 - nP]}
        cell = {"xmin": X[1] - X[0], "xmax": X[-1] - X[-2], "ymin": Y[1] - Y[0], "ymax": Y[-1] - Y[-2]}
        # the clearance: the PML plus 8 cells, and never less than KEEP_H field heights. The cell count
        # alone shrinks as copper comes closer to the edge (the mesh refines around it, so the PML and its
        # 8 cells get thinner): a box moved towards a track would pass for that reason only. A line's
        # fringing field reaches a few dielectric heights sideways (coupling ~ 1 / (1 + (d/h)^2): 4 % at 5 h).
        floor = KEEP_H * float(m.get("field_h_um") or 0.0)
        keep = {s: max(pml[s] + 8 * cell[s], floor) for s in pml}
        if not m["ports"]:
            why.append((None, "no port: no region net crosses the box edge"))
        for p in m["ports"]:
            if p.get("ref") is None:
                why.append((None, f"port {p['net']} ({p['side']}) has no reference plane"))
            la = p["lead"]
            L = abs(la[2] - la[0]) if p["side"][0] == "x" else abs(la[3] - la[1])
            need = pml[p["side"]] + 4 * cell[p["side"]] + 300.0
            if L < need:
                why.append((p["side"], f"port {p['net']} ({p['side']}): lead {L:.0f} um shorter than the PML plus the feed and measurement planes ({need:.0f} um)"))
        why += _track_band_problems(m, keep, cell)
        # pads (and any other signal copper that is not a track): cut by the box edge, or 8 cells clear of the PML
        for s in m["sheets"]:
            if not s["metal"] or s["kind"] in ("plane", "track"):
                continue
            xs, ys = [p[0] for p in s["pts"]], [p[1] for p in s["pts"]]
            for side, d, crosses in (
                    ("xmin", min(xs) - dom[0], min(xs) <= box[0] + 1e-3),
                    ("xmax", dom[2] - max(xs), max(xs) >= box[2] - 1e-3),
                    ("ymin", min(ys) - dom[1], min(ys) <= box[1] + 1e-3),
                    ("ymax", dom[3] - max(ys), max(ys) >= box[3] - 1e-3)):
                if d >= keep[side]:
                    continue
                along = (max(ys) - min(ys)) if side[0] == "x" else (max(xs) - min(xs))
                if crosses and along < 2000.0:
                    continue
                why.append((side, f"{s['kind']} on {s['layer']} at {(min(xs) + max(xs)) / 2000:.2f}, {(min(ys) + max(ys)) / 2000:.2f} mm: "
                                  f"{max(d, 0):.0f} um from the {side} absorbing boundary (keep {keep[side]:.0f} um)"))
                break
        for v in m["vias"]:
            for side, d in (("xmin", v["x"] - dom[0]), ("xmax", dom[2] - v["x"]), ("ymin", v["y"] - dom[1]), ("ymax", dom[3] - v["y"])):
                if d - v["pad_r"] < keep[side]:
                    why.append((side, f"via at ({v['x'] / 1000:.2f}, {v['y'] / 1000:.2f}) mm within {keep[side]:.0f} um of the {side} boundary"))
        why += [(side, w) for side, w in m.get("port_problems", [])]
        # an added via (a variant's return via) must clear every other net's copper: one that lands on a
        # signal via or a track shorts it to ground and the variant measures that short
        for v in m["vias"]:
            if not v.get("added"):
                continue
            for o in m["vias"]:
                if o is not v and o.get("net") != v["net"] and math.hypot(o["x"] - v["x"], o["y"] - v["y"]) < o["pad_r"] + v["pad_r"] + ADD_CLEAR_UM:
                    why.append((None, f"added via at ({v['x'] / 1000:.2f}, {v['y'] / 1000:.2f}) mm overlaps the {o.get('net')} via at "
                                      f"({o['x'] / 1000:.2f}, {o['y'] / 1000:.2f}) mm"))
            for sh in m["sheets"]:
                if not sh.get("metal") or sh["kind"] not in ("track", "pad") or sh.get("net") == v["net"]:
                    continue
                xs, ys = [p_[0] for p_ in sh["pts"]], [p_[1] for p_ in sh["pts"]]
                dx = max(min(xs) - v["x"], 0.0, v["x"] - max(xs))
                dy = max(min(ys) - v["y"], 0.0, v["y"] - max(ys))
                if math.hypot(dx, dy) < v["pad_r"] + ADD_CLEAR_UM and _ring_dist(sh["pts"], v["x"], v["y"]) < v["pad_r"] + ADD_CLEAR_UM:
                    why.append((None, f"added via at ({v['x'] / 1000:.2f}, {v['y'] / 1000:.2f}) mm touches {sh.get('net') or 'a pad'}'s "
                                      f"copper on {sh['layer']}"))
                    break
        # every lane must run through the box to a port on another side or to a termination: a lane that
        # does not is an open end (|S11| = 1 whatever the layout) or a break the solver would take as one
        if m["ports"]:
            conn, by = self.path_connected()
            if not conn:
                cut = [k for k, v in by.items() if not v]
                why.append((None, f"not connected through the box: {', '.join(cut)} reach{'es' if len(cut) == 1 else ''} no port on "
                                  "another side and no termination (an open end, or a break the model would solve as one)"))
        if not any(p.get("excite") for p in m["ports"]):
            why.append((None, "no port is driven"))
        for g in m["groups"]:
            if g["kind"] == "pair" and self.spec.get("excite_mode") != "common":
                e = [p["excite"] for p in g["ports"]]
                if e[0] != -e[1] or (g["side"] == "pad" and any(e)):
                    why.append((None, f"pair {g['ports'][0]['net']}/{g['ports'][1]['net']} not driven with opposite polarity"))
        uniq, sides = [], {}
        for side, w in why:
            if w not in uniq:
                uniq.append(w)
                sides.setdefault(side or "other", []).append(w)
        return {"ok": not uniq, "problems": uniq[:20], "n_problems": len(uniq), "by_side": sides,
                "keep_um": {k: round(v, 1) for k, v in keep.items()}, "pml_um": {k: round(v, 1) for k, v in pml.items()}}

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
                pts += [(seg[0], seg[1], l, na, float(T_["w"][k]) * um, net), (seg[2], seg[3], l, nb, float(T_["w"][k]) * um, net)]

        def near(x, y, l, r, net=None):
            return [nd for (px, py, pl, nd, w, n_) in pts if (l is None or pl == l) and (net is None or n_ == net)
                    and math.hypot(px - x, py - y) <= r + w / 2 + 5]
        # tracks that meet: an end of one on the copper of another (same net and layer), end to end or a T
        segs = []
        for net in nets:
            if net not in b.net_id:
                continue
            for k in b.tracks_of(net):
                seg = _clip_seg(T_["x1"][k] * um, T_["y1"][k] * um, T_["x2"][k] * um, T_["y2"][k] * um, box)
                if seg is not None:
                    segs.append((k, net, b.layers[int(T_["layer"][k])], seg, float(T_["w"][k]) * um))
        for i, (ka, na_, la, sa, wa) in enumerate(segs):
            for kb, nb_, lb, sb, wb in segs[i + 1:]:
                if la != lb or na_ != nb_:
                    continue
                d = min(_pt_seg(sa[0], sa[1], *sb), _pt_seg(sa[2], sa[3], *sb), _pt_seg(sb[0], sb[1], *sa), _pt_seg(sb[2], sb[3], *sa))
                if d <= max(wa, wb) / 2 + 1.0:
                    union(("t", ka, 0), ("t", kb, 0))
        for v in self.model["vias"]:
            ns = near(v["x"], v["y"], None, v["pad_r"], v.get("net"))
            for nd in ns[1:]:
                union(ns[0], nd)
        for f in b.snap["footprints"]:
            pads = [p for p in b.pads_by_ref.get(f["ref"], []) if p["net"] in nets]
            for p in pads:                        # tracks ending on a pad join through the pad
                ns = near(p["x"] * um, p["y"] * um, None, max(p["sx"], p["sy"]) * um / 2, p["net"])
                for nd in ns[1:]:
                    union(ns[0], nd)
        for e in self.model["lumped"]:            # series parts join their two pads' copper
            pads = [p for p in b.pads_by_ref.get(e["ref"], []) if p["net"]]
            groups = [near(p["x"] * um, p["y"] * um, None, max(p["sx"], p["sy"]) * um / 2, p["net"]) for p in pads]
            if len(groups) == 2 and groups[0] and groups[1]:
                union(groups[0][0], groups[1][0])
        node_of = {}
        for p in self.model["ports"] + self.model.get("terms", []):
            ns = near(p["x"], p["y"], p["layer"], p["w"], p["net"])
            node_of[p["n"]] = find(ns[0]) if ns else None
        res = {}
        groups = self.model["groups"]
        for g in groups:
            for i, p in enumerate(g["ports"]):
                others = [q for g2 in groups if g2 is not g and g2["kind"] == g["kind"]
                          and (g2["side"], g2["ports"][0]["layer"]) != (g["side"], g["ports"][0]["layer"])
                          for j, q in enumerate(g2["ports"]) if j == i]
                key = f"{p['net']}@{p['side']}" + (f"/{p['layer']}" if sum(q["net"] == p["net"] and q["side"] == p["side"]
                                                                             for g2 in groups for q in g2["ports"]) > 1 else "")
                res[key] = bool(others) and node_of.get(p["n"]) is not None and \
                    any(node_of.get(q["n"]) == node_of[p["n"]] for q in others)
        return (all(res.values()) if res else False), res


def _ring_dist(pts, x, y):
    """Distance from (x, y) to a polygon (0 inside)."""
    if pip(x, y, pts):
        return 0.0
    return min(_pt_seg(x, y, *pts[i], *pts[(i + 1) % len(pts)]) for i in range(len(pts)))


def _pt_seg(px, py, x1, y1, x2, y2):
    dx, dy = x2 - x1, y2 - y1
    L2 = dx * dx + dy * dy
    t = 0.0 if L2 <= 0 else max(0.0, min(1.0, ((px - x1) * dx + (py - y1) * dy) / L2))
    return math.hypot(px - x1 - t * dx, py - y1 - t * dy)


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


STEEP_SIN = 0.5           # a track meets a side at 30 degrees or more: it crosses it (below: it runs along it)


def _track_band_problems(m, keep, cell):
    """Tracks within `keep` of an absorbing side (see Region.preflight): [(side, message)].

    Per side, the tracks whose copper comes within `keep` of the domain edge are grouped into chains
    (one net, one layer, segments joined end to end inside the band). A chain is benign when it reaches a
    side of the box that is the domain edge (it passes into the PML and is terminated there), has no free
    end inside the band (an end joined to no other segment, or an end cut at a port-side box edge), has no
    segment within 30 degrees of parallel to the side for more than 2 cells, and spans along the side no
    more than a 30-degree crossing does."""
    box, dom = m["box"], m["domain"]
    sides = {"xmin": (0, dom[0], +1), "xmax": (0, dom[2], -1), "ymin": (1, dom[1], +1), "ymax": (1, dom[3], -1)}
    edge = {"xmin": box[0], "xmax": box[2], "ymin": box[1], "ymax": box[3]}
    absorbing = {s: abs(edge[s] - sides[s][1]) < 1e-6 for s in sides}   # no port lead: the box edge is the PML's
    tracks = [s for s in m["sheets"] if s.get("kind") == "track" and s.get("metal") and s.get("seg")]
    key = lambda s, x, y: (s["layer"], s["net"], round(x), round(y))
    joints = {}
    for s in tracks:
        xa, ya, xb, yb = s["seg"]
        if math.hypot(xb - xa, yb - ya) > 1e-6:
            for x, y in ((xa, ya), (xb, yb)):
                joints[key(s, x, y)] = joints.get(key(s, x, y), 0) + 1

    def on_edge(x, y):
        """The box side a clipped end lies on, or None (the end is the track's own)."""
        for s, (ax, _, _) in sides.items():
            if abs((x, y)[ax] - edge[s]) < 1e-3:
                return s
        return None
    out = []
    for side, (ax, d0, sg) in sides.items():
        k, c = keep[side], cell[side]
        dist = lambda p: (p[ax] - d0) * sg                # from the domain edge, inwards
        items = []
        for s in tracks:
            w = s["w"]
            cl = _clip_seg(*s["seg"], box)
            if cl is None:                                   # centreline outside, copper inside: a sliver
                dmin = min(dist(p) for p in s["pts"])        # cut lengthwise by the box edge
                if dmin < k:
                    out.append((side, f"track of {s['net'] or 'no net'} on {s['layer']} is cut lengthwise by the box edge, "
                                      f"{max(dmin, 0):.0f} um from the {side} absorbing boundary (keep {k:.0f} um)"))
                continue
            lim = d0 + sg * (k + w / 2)
            band = list(box)                                 # the box within keep (+ half the track) of this side
            if sg > 0:
                band[ax + 2] = min(band[ax + 2], lim)
            else:
                band[ax] = max(band[ax], lim)
            if band[0] > band[2] or band[1] > band[3]:
                continue
            part = _clip_seg(*cl, band)
            if part is None:
                continue
            xa, ya, xb, yb = s["seg"]
            L = math.hypot(xb - xa, yb - ya)
            steep = L > 1e-6 and abs((xb - xa, yb - ya)[ax]) / L >= STEEP_SIN
            items.append({"s": s, "cl": cl, "part": part, "steep": steep, "L": L})
        if not items:
            continue
        # chains: segments of one net and layer joined at an end inside the band
        par = list(range(len(items)))

        def find(i):
            while par[i] != i:
                par[i] = par[par[i]]
                i = par[i]
            return i
        at = {}
        for i, it in enumerate(items):
            xa, ya, xb, yb = it["s"]["seg"]
            for x, y in ((xa, ya), (xb, yb)):
                kk = key(it["s"], x, y)
                if kk in at:
                    par[find(i)] = find(at[kk])
                else:
                    at[kk] = i
        chains = {}
        for i in range(len(items)):
            chains.setdefault(find(i), []).append(items[i])
        for ch in chains.values():
            s0 = ch[0]["s"]
            wmax = max(it["s"]["w"] for it in ch)
            what, reaches = [], False
            near = min(min(dist(it["part"][:2]), dist(it["part"][2:])) - it["s"]["w"] / 2 for it in ch)
            for it in ch:
                s, w = it["s"], it["s"]["w"]
                if it["L"] <= 1e-6:
                    what.append("is a dot of copper")
                    continue
                xa, ya, xb, yb = s["seg"]
                cl = it["cl"]
                for (ox, oy), (cx, cy) in (((xa, ya), cl[:2]), ((xb, yb), cl[2:])):
                    if math.hypot(ox - cx, oy - cy) > 1e-3:            # cut at the box edge
                        e = on_edge(cx, cy)
                        if e and absorbing[e]:
                            reaches = True
                        elif dist((cx, cy)) < k + w / 2:
                            what.append("is cut at the box edge beside a port lead")
                    elif dist((ox, oy)) < k + w / 2 and joints.get(key(s, ox, oy), 0) < 2:
                        what.append("ends")
                pl = math.hypot(it["part"][2] - it["part"][0], it["part"][3] - it["part"][1])
                if not it["steep"] and pl > 2 * c:
                    what.append(f"runs along it for {pl:.0f} um")
            ta = 1 - ax
            span = [v for it in ch for v in (it["part"][ta], it["part"][ta + 2])]
            if max(span) - min(span) > 1.8 * k + 2 * wmax + c:
                what.append(f"spans {max(span) - min(span):.0f} um along it")
            if not reaches and not what:
                what.append("does not reach an absorbing side of the box")
            if what:
                uniq = list(dict.fromkeys(what))
                out.append((side, f"track of {s0['net'] or 'no net'} on {s0['layer']} {', '.join(uniq)}, "
                                  f"{max(near, 0):.0f} um from the {side} absorbing boundary (keep {k:.0f} um)"))
    return out


# ---- what a run will cost ---------------------------------------------------------------------------
SPEED_MCS = 350.0        # Mcells/s assumed for a machine that has finished no grid yet (measured runs: 300-900)
SETTLE = 2.0             # simulated time to the -40 dB end criterion, in excitation lengths (J1 launch: 1.6)
STEP_CAP = 4.0         # a run's step limit: this many times the steps its excitation needs (run_estimate)
TAIL_DB = -50.0        # a run cut on the step limit counts when every port is this far below the largest peak
MAX_HOURS = 6.0          # a region predicted to take longer (typically) is refused; on the clock, it is stopped there


def run_estimate(model, speed_mcs=None, slow_mcs=None):
    """The FDTD time step (Courant, from the smallest cell on each axis: openEMS's own step is at most
    ~1.3 times longer), the simulated time a run needs (openEMS's Gaussian excitation, 2 * 9 / (2 pi fc),
    then SETTLE excitation lengths to the end criterion), the time steps that takes, whether the step
    limit cuts the run short, and the hours at `speed_mcs` Mcells/s: typical (`hours`, a structure that
    settles) and at most (`hours_max`, a run that goes to its step limit, which is enforced)."""
    d = [float(np.min(np.diff(np.array(model["mesh"][a])))) * 1e-6 for a in "xyz"]
    dt = 1.0 / (C0 * math.sqrt(sum(1.0 / v ** 2 for v in d)))
    fc = model["fmax"] / 2
    exc = 2 * 9.0 / (2 * math.pi * fc)
    need = SETTLE * exc
    steps = need / dt
    cap = int(model.get("max_timesteps", 400000))
    cells = len(model["mesh"]["x"]) * len(model["mesh"]["y"]) * len(model["mesh"]["z"])
    sp = float(speed_mcs or SPEED_MCS)
    return {"dt_fs": round(dt * 1e15, 3), "min_cell_um": [round(v * 1e6, 3) for v in d], "sim_ps": round(need * 1e12, 1),
            "steps": int(steps), "step_limit": cap, "hits_step_limit": steps > cap, "cells": cells,
            "speed_mcs": round(sp, 1), "speed_measured": speed_mcs is not None,
            "hours": round(cells * min(steps, cap) / (sp * 1e6) / 3600, 2),
            "hours_max": round(cells * cap / (min(sp, float(slow_mcs or sp)) * 1e6) / 3600, 2)}


def measured_speed(cache_dir, slowest=False):
    """This machine's openEMS speed (Mcells/s) over the last grids it finished: the median, or with
    `slowest` the slowest of them (shared with other solves, a run can drop to half its speed), or None."""
    try:
        v = json.load(open(os.path.join(cache_dir, "em_speed.json"))).get("mcs", [])
        return (float(np.min(v)) if slowest else float(np.median(v))) if v else None
    except (OSError, ValueError, TypeError):
        return None


def _record_speed(cache_dir, logtxt):
    m = re.findall(r"Speed:\s*([0-9.]+)\s*MCells/s", logtxt)
    if not m or not cache_dir:
        return
    p = os.path.join(cache_dir, "em_speed.json")
    try:
        v = json.load(open(p)).get("mcs", [])
    except (OSError, ValueError):
        v = []
    v = (v + [float(m[-1])])[-10:]
    os.makedirs(cache_dir, exist_ok=True)
    tmp = p + f".{os.getpid()}.tmp"
    json.dump({"mcs": v}, open(tmp, "w"))
    os.replace(tmp, p)


# ---- running ---------------------------------------------------------------------------------------
RUNNER = os.path.join(os.path.dirname(os.path.abspath(__file__)), "em3d.py")


def grid_key(model):
    """Hash of everything one grid's solve depends on: the model (geometry, materials, mesh, ports,
    excitation, boundaries, end criterion and time-step limit) and the runner that turns it into openEMS."""
    h = hashlib.sha256(json.dumps(model, sort_keys=True, default=str).encode())
    h.update(open(RUNNER, "rb").read())
    return h.hexdigest()


def _reuse_grid(gdir, key):
    """A finished grid of an earlier (possibly interrupted) run of the same region: its record says it
    completed, its inputs hash the same, and its port waves are complete. None otherwise."""
    try:
        done = json.load(open(os.path.join(gdir, "done.json")))
        if done.get("key") != key:
            return None
        ports = json.load(open(os.path.join(gdir, "ports.json")))
        if hashlib.sha256(open(os.path.join(gdir, "ports.json"), "rb").read()).hexdigest() != done.get("ports_sha256"):
            return None
        if not ports.get("ports") or any(len(v.get("a_re", [])) != len(ports["freq_hz"]) for v in ports["ports"].values()):
            return None
        return {"ports": ports, "ended_by": done["ended_by"], "seconds": done["seconds"], "cells": ports["cells"], "reused": True,
                "tail_db": done.get("tail_db", tail_db(ports))}
    except (OSError, ValueError, KeyError, TypeError):
        return None


DIVERGE = 100.0        # energy this many times its lowest point after the excitation has passed: the run is unstable


def _solve(cmd, env, deadline):
    """Run openEMS, reading its progress as it goes: (returncode, log text, stopped) where `stopped` is
    None, "unstable" (the field energy, after falling at least 20 dB below its peak, grew DIVERGE times
    above its lowest point: an FDTD instability, which only grows) or "time" (past `deadline`)."""
    import threading
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, env=env, bufsize=1)
    lines, state = [], {"stopped": None, "peak": 0.0, "low": None}

    def watch():
        while p.poll() is None:
            if deadline and time.time() > deadline:
                state["stopped"] = "time"
                p.kill()
                return
            time.sleep(5)
    threading.Thread(target=watch, daemon=True).start()
    for line in p.stdout:
        lines.append(line)
        m = re.search(r"Energy: ~([0-9.eE+-]+)", line)
        if not m:
            continue
        try:
            e = float(m.group(1))
        except ValueError:
            continue
        state["peak"] = max(state["peak"], e)
        if e > 0 and e < 1e-2 * state["peak"]:
            state["low"] = e if state["low"] is None else min(state["low"], e)
        if state["low"] and e > DIVERGE * state["low"]:
            state["stopped"] = "unstable"
            p.kill()
            break
    p.wait()
    txt = "".join(lines)
    if state["stopped"] == "unstable":
        txt += (f"\nopenpcbsi: stopped, unstable: the field energy fell to {state['low']:.3g} and then grew past "
                f"{DIVERGE:g} times that\n")
    elif state["stopped"] == "time":
        txt += "\nopenpcbsi: stopped at the region's time limit (max_hours)\n"
    return p.returncode, txt, state["stopped"]


def _run_grid(model, od, gi, threads, log, cache_dir=None, deadline=None):
    from .si import inside
    gdir = inside(od, os.path.join(od, f"grid{gi}"))
    os.makedirs(gdir, exist_ok=True)
    key = grid_key(model)
    old = _reuse_grid(gdir, key)
    if old:
        log(f"  grid {gi}: reused the finished solve in {gdir} ({old['cells']} cells, ended by {old['ended_by']})")
        return old
    for f in ("done.json", "ports.json"):           # a stale or partial grid: solved again from scratch
        if os.path.exists(os.path.join(gdir, f)):
            os.remove(os.path.join(gdir, f))
    py = openems_python()
    mp = os.path.join(gdir, "model.json")
    json.dump(model, open(mp, "w"))
    t0 = time.time()
    from . import jobs
    t_wait = time.time()
    with jobs.slot("openems", log) as share:
        # time spent waiting for a job slot (other solves on the machine) does not count against max_hours
        if isinstance(deadline, dict):
            deadline["t"] += time.time() - t_wait
        t0 = time.time()
        th = threads or share
        rc, logtxt, stopped = _solve([py, "-I", RUNNER, mp, gdir, str(th)], dict(os.environ, OMP_NUM_THREADS=str(th)),
                                     deadline["t"] if isinstance(deadline, dict) else deadline)
    open(os.path.join(gdir, "runner.log"), "w").write(logtxt)
    if stopped:
        # no ports.json and no done.json: nothing to read, and a rerun solves it again
        secs = round(time.time() - t0, 1)
        log(f"  grid {gi}: stopped after {secs:.0f} s: " + ("unstable (the field energy grew after decaying)" if stopped == "unstable"
                                                          else "the region's time limit"))
        return {"ports": None, "ended_by": stopped, "seconds": secs, "cells": None, "tail_db": None}
    if rc == 0:
        _record_speed(cache_dir, logtxt)
    if rc != 0:
        raise RuntimeError(f"openEMS run failed (grid {gi}): {logtxt[-1500:]}")
    ports = json.load(open(os.path.join(gdir, "ports.json")))
    # openEMS prints "Max. number of timesteps: N" at the start of every run; only this warning means
    # the run stopped on the limit
    ended = "timesteps" if re.search(r"number of timesteps was reached", logtxt) else "end_criteria"
    secs = round(time.time() - t0, 1)
    tdb = tail_db(ports)
    log(f"  grid {gi}: {ports['cells']} cells, {secs:.0f} s, ended by {ended}" + (f" (port waves {tdb:.0f} dB below their peak)" if tdb is not None else ""))
    # written last, atomically: its presence means this grid finished (an interrupted run leaves none)
    rec = {"key": key, "ended_by": ended, "seconds": secs, "cells": ports["cells"], "tail_db": tdb,
           "ports_sha256": hashlib.sha256(open(os.path.join(gdir, "ports.json"), "rb").read()).hexdigest(),
           "date": time.strftime("%Y-%m-%d %H:%M:%S")}
    tmp = os.path.join(gdir, f"done.json.{os.getpid()}.tmp")
    with open(tmp, "w") as fh:
        json.dump(rec, fh, indent=1)
    os.replace(tmp, os.path.join(gdir, "done.json"))
    return {"ports": ports, "ended_by": ended, "seconds": secs, "cells": ports["cells"], "tail_db": tdb}


def tail_db(ports):
    """The ports' voltage over the last tenth of the run against the largest port peak, in dB; None
    for a grid solved before the runner recorded it."""
    P = list(ports.get("ports", {}).values())
    pk = max([v.get("u_peak") or 0.0 for v in P] or [0.0])
    tl = [v.get("u_tail") for v in P if v.get("u_tail") is not None]
    if not pk or not tl or any(t != t for t in tl):
        return None
    return round(20 * math.log10(max(max(tl), 1e-30) / pk), 1)


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
    # at the link's own rate: S-parameters at its Nyquist, and the in-box reflections (time-gated to the
    # box's window, so neither port lead counts) seen with the link's own edge
    fn = model.get("f_nyq_hz")
    if fn and fn <= f.max():
        i = int(np.argmin(np.abs(f - fn)))
        r_k = next((k for k in ("Sdd11_db", "S11_db", "Scc11_db") if k in out), None)
        t_k = next((k for k in ("Sdd21_db", "S21_db", "Scc21_db") if k in out), None)
        out["at_nyquist"] = {"f_ghz": round(float(f[i]) / 1e9, 3), "reflection_db": out[r_k][i] if r_k else None,
                             "transmission_db": out[t_k][i] if t_k else None}
    win = window_ps(model)
    if win:
        rise = tdr_rise_ps(f.max()) * 1e-12
        rn = max(float(model.get("rise_net_ps") or 0.0) * 1e-12, rise)
        tg, zg = T.gated_tdr(fu, T.renorm(s_i, z_port, zref_s), 1.0, win, rise, rn)
        m = tg > 0
        if m.any():
            out["tdr_net"] = {"min_pct": round(float(zg[m].min() - 1) * 100, 2), "max_pct": round(float(zg[m].max() - 1) * 100, 2),
                              "rise_ps": round(rn * 1e12, 1), "window_ps": [round(v, 1) for v in win]}
    return out


def net_rate(ctx, nets, spec=None):
    """(Nyquist Hz, 10-90 rise ps, where from) of the fastest of the region's nets: the spec's `rise_ps` /
    `nyquist_ghz`, else the net class's `rise_ps` and `nyquist_ghz`; a rise not given is 30 % of a bit,
    0.15 / Nyquist (an assumption, stated in the result). (None, None, None) without a class rate."""
    spec = spec or {}
    fn = float(spec["nyquist_ghz"]) * 1e9 if spec.get("nyquist_ghz") else None
    rise = float(spec["rise_ps"]) if spec.get("rise_ps") else None
    src = "region spec" if (fn or rise) else None
    for n in nets:
        c = ctx.cfg.class_of(n)
        if not c:
            continue
        q = ctx.cfg.class_param(c, "nyquist_ghz", 0) or 0
        r = ctx.cfg.class_param(c, "rise_ps", 0) or 0
        if q and not spec.get("nyquist_ghz") and (fn is None or q * 1e9 > fn):
            fn, src = q * 1e9, f"net class {c}"
        if r and not spec.get("rise_ps") and (rise is None or r < rise):
            rise, src = float(r), f"net class {c}"
    if fn and not rise:
        rise, src = 0.15 / fn * 1e12, (src or "") + " (rise taken as 30 % of a bit, 0.15 / Nyquist)"
    return fn, rise, src


GRID_PTS = 1.5         # grids agree on an excursion within this many percentage points ...
GRID_REL = 0.15        # ... or within this fraction of the excursion, whichever is larger
GRID_S21_DB = 0.1      # ... and on the transmission at Nyquist within this many dB (or GRID_REL of it)


def grid_agreement(a, b, limit=None):
    """Two grids' results, each {"exc": [min, max] % in the box against one reference, "net": [min, max] %
    at the link's own edge or None, "s21": dB at the link's Nyquist or None}: {converged, agreed_on,
    reason, undecided, spread}. They agree on the metric when every excursion matches within
    max(GRID_PTS, GRID_REL x its size) and the transmission within max(GRID_S21_DB, GRID_REL x it); on
    the decision when both put the judged excursion (at the link's edge when there is one) on the same
    side of `limit`. Either is converged; neither is undecided."""
    limit = TDR_LIMIT_PCT if limit is None else limit
    spread, ok = {}, True
    for k in ("exc", "net"):
        if a.get(k) is None or b.get(k) is None:
            continue
        size = max(abs(v) for v in list(a[k]) + list(b[k]))
        tol = max(GRID_PTS, GRID_REL * size)
        d = max(abs(a[k][0] - b[k][0]), abs(a[k][1] - b[k][1]))
        spread[k] = {"diff_pts": round(d, 2), "tol_pts": round(tol, 2)}
        ok &= d <= tol
    if a.get("s21") is not None and b.get("s21") is not None:
        d = abs(a["s21"] - b["s21"])
        tol = max(GRID_S21_DB, GRID_REL * max(abs(a["s21"]), abs(b["s21"])))
        spread["s21"] = {"diff_db": round(d, 3), "tol_db": round(tol, 3)}
        ok &= d <= tol
    key = "net" if (a.get("net") is not None and b.get("net") is not None) else "exc"
    if a.get(key) is None or b.get(key) is None:
        # a box too short for an in-box TDR window (the leads' steps cover the round trip at this edge):
        # the grids are compared on the transmission at the link's Nyquist alone, and the result says so
        if "s21" in spread:
            agree = spread["s21"]["diff_db"] <= spread["s21"]["tol_db"]
            return {"converged": agree, "agreed_on": "transmission" if agree else None, "undecided": False, "spread": spread,
                    "reason": ("the box is too short for an in-box TDR window, so the grids are compared on the transmission at "
                               f"the link's Nyquist alone: within {spread['s21']['diff_db']} of {spread['s21']['tol_db']} dB"
                               + ("" if agree else " (not within)"))}
        return {"converged": False, "undecided": False, "reason": "no in-box TDR to compare (box too short for a window) and no "
                "transmission at the link's Nyquist", "spread": spread}
    wa, wb = max(abs(v) for v in a[key]), max(abs(v) for v in b[key])
    same = (wa <= limit) == (wb <= limit)
    what = "at the link's own edge" if key == "net" else "in the box"
    if ok:
        return {"converged": True, "agreed_on": "result", "undecided": False, "spread": spread,
                "reason": f"the grids agree on the TDR {what} (" + ", ".join(
                    f"{k} within {v['diff_pts']} of {v['tol_pts']} points" if "diff_pts" in v else f"{k} within {v['diff_db']} of {v['tol_db']} dB"
                    for k, v in spread.items()) + ")"}
    if same:
        return {"converged": True, "agreed_on": "decision", "undecided": False, "spread": spread,
                "reason": f"the grids differ ({spread}), but both put the TDR {what} ({wa:.1f} and {wb:.1f} %) "
                          f"{'within' if wa <= limit else 'beyond'} the {limit:g} % limit: the same verdict"}
    return {"converged": False, "agreed_on": None, "undecided": True, "spread": spread,
            "reason": f"the grids straddle the {limit:g} % limit ({wa:.1f} and {wb:.1f} % {what})"}


def _convergence(grids, runs, zref):
    """The grid comparison of a solved region (see grid_agreement): the first and the last grid given."""
    conv = {"runs": [{"res_scale": r["res_scale"], "cells": r["cells"], "seconds": r["seconds"], "ended_by": r["ended_by"],
                      "z_line_ohm": r["z_line_ohm"], "tail_db": r.get("tail_db"), "z_port_ohm": r["post"].get("z_port_ohm")}
                     for r in runs], "reference_ohm": round(float(zref), 2) if zref else None}
    cut = [r for r in runs if r["ended_by"] != "end_criteria"]
    loud = [r for r in cut if r.get("tail_db") is None or r["tail_db"] > TAIL_DB]
    if loud:
        conv.update(converged=False, undecided=False, reason=("a run hit the time-step limit with its port waves still "
                    f"{max((r.get('tail_db') or 0) for r in loud):.0f} dB from their peak (limit {TAIL_DB:g} dB)"))
        return conv
    g = []
    for (_, m_), r in zip([grids[0], grids[-1]], [runs[0], runs[-1]]):
        e = tdr_excursion(m_, {"tdr_t_ps": r["post"].get("tdr_t_ps"), "tdr_z_ohm": r["post"].get("tdr_z_ohm"),
                               "z_line_ohm": None, "z_2d_ohm": zref})
        n = r["post"].get("tdr_net")
        ny = r["post"].get("at_nyquist") or {}
        g.append({"exc": e and [e["min_pct"], e["max_pct"]], "net": n and [n["min_pct"], n["max_pct"]], "s21": ny.get("transmission_db")})
    conv["per_grid"] = g
    # for the record: the in-box shape with each grid's own level removed (a level shift is the port lead's)
    w = window_ps(grids[-1][1])
    if w:
        ts = [np.array(r["post"].get("tdr_t_ps") or []) for r in (runs[0], runs[-1])]
        zs = [np.array(r["post"].get("tdr_z_ohm") or []) for r in (runs[0], runs[-1])]
        n = min(len(ts[0]), len(ts[1]))
        m = (ts[0][:n] > w[0]) & (ts[0][:n] < w[1])
        if n and m.any():
            q = [z_[:n][m] / np.mean(z_[:n][m]) for z_ in zs]
            conv["shape_diff_level_removed_pct"] = round(float(np.max(np.abs(q[0] - q[1]))) * 100, 2)
    conv.update(grid_agreement(g[0], g[-1]))
    if cut and conv["converged"]:
        conv["reason"] += (f"; {len(cut)} run(s) stopped on the step limit with the port waves {max(r['tail_db'] for r in cut):.0f} dB "
                           "below their peak (energy left ringing in the model's planes, which the board's decoupling would damp)")
    return conv


def window_ps(model):
    """The TDR window inside the box (ps): from 0.6 rise times after the port (past the near lead's step)
    to 0.6 short of the round trip to the far side (before the far lead's); to one rise time past it when a
    termination in the box ends the line (its pad is part of what is judged). None when the box is too
    short for one."""
    pe = excited_port(model)
    if not pe:
        return None
    box = model["box"]
    L = (box[2] - box[0]) if pe["side"][0] == "x" else (box[3] - box[1])
    tl = 2 * L * 1e-6 * math.sqrt(model.get("eps_eff_2d") or 3.0) / C0 * 1e12
    rise = tdr_rise_ps(model["fmax"])
    lo, hi = 0.6 * rise, tl - 0.6 * rise
    if any(g["side"] == "pad" for g in model["groups"]):
        hi = tl + rise
    return (lo, hi) if hi > lo else None


def excited_port(model):
    """The excited group's first port (P of a pair, driven +1), whatever order the nets were listed in; its
    run is the shorter of the group's (the line reading needs both legs uniform). None without one."""
    g = next((g for g in model["groups"] if any(p.get("excite") for p in g["ports"])), None)
    if g is None:
        return None
    p = dict(g["ports"][0])
    p["run_um"] = min(q.get("run_um") or 0.0 for q in g["ports"])
    return p


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
    p = excited_port(model)
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
    # the run folder: the region's name (from net and part names: through safe_name, so a hierarchical
    # "/sheet/NET" or a "../" in a name cannot leave runs/em) and a hash of what the result depends on
    # (not the thread count or whether this call only checks: a rerun finds the grids an earlier one finished)
    from .si import inside, safe_name
    hspec = {k: v for k, v in spec.items() if k not in ("preflight_only", "threads")}
    key = hashlib.sha256(json.dumps({"m": model, "s": hspec}, sort_keys=True, default=str).encode()).hexdigest()[:16]
    em_dir = os.path.join(runs_dir, "em")
    od = inside(em_dir, os.path.join(em_dir, f"{safe_name(spec.get('name') or 'region')}_{key}"))
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
    res["terminations"] = [{k: t[k] for k in ("ref", "pad", "net", "R", "ref_layer")} for t in model.get("terms", [])]
    res["left_out"] = dict(model.get("left_out") or {}, include_um=model.get("include_um"))
    # what the run will cost, every grid built first: refused before solving when a run would stop on the
    # step limit or take longer than max_hours (a 1 nm cell made one model need weeks)
    grids = []
    for scale in spec.get("res_scales", [1.0, 0.7]):
        m2 = Region(ctx, spec, scale).build()
        m2["eps_eff_2d"] = model.get("eps_eff_2d")
        # the step limit: STEP_CAP times the steps the excitation needs (the spec's max_timesteps at most).
        # A structure that keeps ringing (a plane fragment between ground planes, which the real board's
        # decoupling damps) stops there instead of running for hours; it is accepted when its port waves
        # have died away (TAIL_DB), whatever is left in the planes
        e0 = run_estimate(m2)
        m2["max_timesteps"] = int(min(m2["max_timesteps"], max(20000, STEP_CAP * e0["steps"])))
        grids.append((scale, m2))
    speed = measured_speed(cache_dir)
    slow = min(measured_speed(cache_dir, slowest=True) or SPEED_MCS, SPEED_MCS)
    est = [dict(run_estimate(m2, speed, slow), res_scale=sc) for sc, m2 in grids]
    hours = round(sum(e["hours"] for e in est), 2)
    hours_max = round(sum(e["hours_max"] for e in est), 2)
    max_h = float(spec.get("max_hours", os.environ.get("OPENPCBSI_EM_MAX_HOURS", MAX_HOURS)))
    res["estimate"] = {"grids": est, "hours": hours, "hours_max": hours_max, "max_hours": max_h}
    for e in est:
        if e["hits_step_limit"]:
            pf["problems"].append(f"grid {e['res_scale']}: a run needs {e['steps']} time steps of {e['dt_fs']} fs (smallest cells "
                                  f"{e['min_cell_um']} um) and would stop on the limit of {e['step_limit']} before the end criterion")
    if hours > max_h:
        pf["problems"].append(f"predicted {hours:.1f} h (up to {hours_max:.1f} h) for {len(est)} grids ({', '.join(str(e['cells']) for e in est)} cells at "
                              f"{est[0]['speed_mcs']:.0f} Mcells/s{'' if speed else ', not yet measured on this machine'}), over "
                              f"max_hours = {max_h:g}: make the box smaller or raise max_hours")
    # calibration: whether the excited lane's uniform run can hold a TDR plateau, said before solving
    need = plateau_um(model["fmax"], model.get("eps_eff_2d"))
    pe = excited_port(model)
    if pe:
        res["calibration_plan"] = {"run_um": pe.get("run_um"), "needs_um": round(need), "possible": (pe.get("run_um") or 0) >= need}
    pf["n_problems"] = len(pf["problems"])
    pf["ok"] = not pf["problems"]
    cp = res.get("calibration_plan") or {}
    if hours_max > max_h and hours <= max_h:
        log(f"  {spec.get('name')}: a run that goes to its step limit would take up to {hours_max:.1f} h; it is stopped at "
            f"max_hours = {max_h:g} (no answer then)")
    log(f"  {spec.get('name')}: predicted {hours:.1f} h, at most {hours_max:.1f} h ({', '.join(str(e['cells']) for e in est)} cells, step {est[0]['dt_fs']} fs, "
        f"{est[0]['speed_mcs']:.0f} Mcells/s); {len(res['terminations'])} terminations; line calibration "
        + ("possible" if cp.get("possible") else f"not possible in this box (uniform run {(cp.get('run_um') or 0) / 1000:.2f} mm, "
           f"needs {(cp.get('needs_um') or 0) / 1000:.1f} mm): TDR against the 2D impedance of the same geometry"))
    if spec.get("preflight_only") or not pf["ok"]:
        res["ran"] = False
        res["verdict"] = verdict(res)
        json.dump(res, open(os.path.join(od, "result.json"), "w"), indent=1)
        return res
    threads = int(spec.get("threads", os.environ.get("OPENPCBSI_EM_THREADS", 0)))   # 0: the job slot's share
    runs = []
    clock = {"t": time.time() + max_h * 3600.0}  # max_hours holds on the clock (waits for a job slot excepted)
    fnq, rnet, rsrc = net_rate(ctx, reg.nets, spec)
    res["link"] = {"nyquist_ghz": fnq and round(fnq / 1e9, 3), "rise_ps": rnet and round(rnet, 1), "from": rsrc}

    def solve_grid(gi, scale, m2):
        rr = _run_grid(m2, od, gi, threads, log, cache_dir, clock)
        if rr["ended_by"] == "unstable" and not spec.get("no_retry"):
            # an FDTD instability belongs to one mesh (rev D REFCLK0: unstable at res_scale 1.0, stable at 0.9
            # and in a larger box, with either absorber): the grid once more at 0.9 of the cell size
            s2 = round(scale * 0.9, 3)
            log(f"  grid {gi}: unstable at res_scale {scale}; solved again at {s2}")
            m2 = Region(ctx, spec, s2).build()
            m2["eps_eff_2d"] = model.get("eps_eff_2d")
            m2["max_timesteps"] = int(min(m2["max_timesteps"], max(20000, STEP_CAP * run_estimate(m2)["steps"])))
            grids[gi] = (s2, m2)
            rr = _run_grid(m2, od, f"{gi}r", threads, log, cache_dir, clock)
            rr["retried"] = {"unstable_at": scale, "res_scale": s2}
            scale = s2
        if rr["ports"] is None:
            return rr
        m2["f_nyq_hz"], m2["rise_net_ps"] = fnq, rnet       # post-processing only (after the grid's key)
        rr["post"] = _post(m2, rr)
        rr["z_line_ohm"] = _line_z(m2, rr["post"])
        rr["res_scale"] = scale
        return rr
    for gi, (scale, m2) in enumerate(grids):
        if runs and not runs[-1].get("reused"):
            # the next grid predicted from the one just solved (its measured time, scaled by cells and time
            # step): one that would pass max_hours is not started
            e_prev = est[gi - 1]
            e_next = est[gi]
            pred = runs[-1]["seconds"] * (e_next["cells"] / max(e_prev["cells"], 1)) * (e_prev["dt_fs"] / max(e_next["dt_fs"], 1e-9))
            if time.time() + pred > clock["t"]:
                res.update(ran=True, stopped={"grid": gi, "res_scale": scale, "why": "time", "seconds": 0.0,
                                              "predicted_s": round(pred)}, seconds=round(time.time() - t0, 1))
                res["convergence"] = {"converged": False, "reason": (f"grid {gi} not started: from grid {gi - 1}'s {runs[-1]['seconds']:.0f} s it "
                                                                     f"would take about {pred / 3600:.1f} h, past max_hours = {max_h:g}")}
                res["verdict"] = verdict(res)
                json.dump(res, open(os.path.join(od, "result.json"), "w"), indent=1)
                return res
        rr = solve_grid(gi, scale, m2)
        if rr["ports"] is None:                 # stopped: unstable, or out of time; no answer from this region
            res.update(ran=True, stopped={"grid": gi, "res_scale": scale, "why": rr["ended_by"], "seconds": rr["seconds"]},
                       seconds=round(time.time() - t0, 1))
            res["convergence"] = {"converged": False, "reason": (f"grid {gi} was unstable: the field energy grew after it had decayed "
                                                                 "(an FDTD instability), so it was stopped" if rr["ended_by"] == "unstable"
                                                                 else f"grid {gi} was stopped at max_hours = {max_h:g}")}
            res["verdict"] = verdict(res)
            json.dump(res, open(os.path.join(od, "result.json"), "w"), indent=1)
            return res
        runs.append(rr)
    fine = runs[-1]
    res.update({k: v for k, v in fine["post"].items()})
    res["runs_per_grid"] = 1
    res["z_line_ohm"] = fine["z_line_ohm"]
    pe = excited_port(model)
    if pe:                                     # the box length along the excited port's axis
        bx = model["box"]
        res["line_length_mm"] = round(((bx[2] - bx[0]) if pe["side"][0] == "x" else (bx[3] - bx[1])) / 1000.0, 3)
    res["eps_eff"] = model.get("eps_eff_2d")
    # convergence (#15): the grids compared on what is reported, the in-box TDR excursion against one
    # reference for both (the 2D line), the excursion at the link's own edge and the transmission at its
    # Nyquist; accepted when they agree within a tolerance scaled to the result, or when both land on the
    # same side of the limit the region is judged against. When they straddle it, one more grid decides.
    zr = res.get("z_2d_ohm") or float(np.mean([r["post"].get("z_port_ohm") or 0 for r in runs]))
    conv = _convergence(grids, runs, zr)
    if conv.get("undecided") and spec.get("third_grid", True) and len(grids) >= 2:
        s3 = round(grids[-1][0] * 0.75, 3)
        m3 = Region(ctx, spec, s3).build()
        m3["eps_eff_2d"] = model.get("eps_eff_2d")
        m3["max_timesteps"] = int(min(m3["max_timesteps"], max(20000, STEP_CAP * run_estimate(m3)["steps"])))
        e3 = run_estimate(m3, speed, slow)
        if time.time() + e3["hours"] * 3600 < clock["t"]:
            log(f"  {spec.get('name')}: the grids straddle the limit; a third grid (res_scale {s3}, {e3['cells']} cells, ~{e3['hours']:.1f} h) decides")
            rr3 = solve_grid(len(grids), s3, m3)
            if rr3["ports"] is not None:
                grids.append((s3, m3))
                runs.append(rr3)
                conv = _convergence(grids[-2:], runs[-2:], zr)
                conv["third_grid"] = f"res_scale {s3}: compared with the grid before it"
            else:
                conv["third_grid"] = f"res_scale {s3}: stopped ({rr3['ended_by']})"
        else:
            conv["third_grid"] = f"not run: ~{e3['hours']:.1f} h would pass max_hours = {max_h:g}"
        fine = runs[-1]
        res.update({k: v for k, v in fine["post"].items()})
        res["z_line_ohm"] = fine["z_line_ohm"]
    res["convergence"] = conv
    # calibration (#16): the model's line against 2D for the same geometry
    if res.get("z_2d_ohm") and res.get("z_line_ohm"):
        ag = (res["z_line_ohm"] - res["z_2d_ohm"]) / res["z_2d_ohm"] * 100
        res["calibration"] = {"z_model_ohm": res["z_line_ohm"], "z_2d_ohm": res["z_2d_ohm"], "agree_pct": round(ag, 2),
                              "renormalised": abs(ag) > 2.0,
                              "note": "TDR reported relative to the model's own line" if abs(ag) > 2.0 else "within 2 % of 2D"}
    elif res.get("z_2d_ohm"):
        pe = excited_port(model) or {}
        res["calibration"] = {"z_model_ohm": None, "z_2d_ohm": res["z_2d_ohm"], "agree_pct": None, "renormalised": False,
                              "note": f"not calibrated: the excited port's track runs {pe.get('run_um', 0) / 1000:.1f} mm "
                                      "uniform inside the box, too short for a TDR plateau at this rise time "
                                      "(move the box edge further from the first part or via)"}
    res["tdr_excursion"] = tdr_excursion(model, res)
    res["reflection_floor_db"] = reflection_floor(res)
    # clean, below the model's resolution: both grids within CLEAN_PCT in the box (at the link's edge when
    # there is one) and the return loss under the floor the model's own leads set. A positive answer about
    # the layout (FS9 met), not "not usable"
    g = (res.get("convergence") or {}).get("per_grid") or []
    fl = res.get("reflection_floor_db")
    rk = next((k for k in ("Sdd11_db", "S11_db", "Scc11_db") if res.get(k)), None)
    if res["convergence"].get("converged") and g and fl is not None and rk:
        key = "net" if all(x.get("net") for x in g) else "exc"
        f_ = np.array(res["freq_ghz"])
        band = f_ <= ((res.get("link") or {}).get("nyquist_ghz") or f_.max())
        worst_s11 = float(np.max(np.array(res[rk])[band & (f_ > 0.02 * f_.max())])) if (band & (f_ > 0.02 * f_.max())).any() else None
        if all(x.get(key) and max(abs(v) for v in x[key]) <= CLEAN_PCT for x in g) and worst_s11 is not None and worst_s11 <= fl:
            res["clean"] = {"tdr_pct": max(max(abs(v) for v in x[key]) for x in g), "on": key, "s11_db": round(worst_s11, 1), "floor_db": fl}
    if res.get("passive") is False:
        res.setdefault("problems", []).append("result is not passive")
    res["ran"] = True
    res["seconds"] = round(time.time() - t0, 1)
    res["verdict"] = verdict(res)
    json.dump(res, open(os.path.join(od, "result.json"), "w"), indent=1)
    return res


def tdr_excursion(model, res):
    """The TDR's lowest and highest point inside the box (from 0.6 rise times after the port to one rise
    time past the box's round trip along the excited port's axis), in % of the reference line: the
    model's own calibrated line, else the 2D impedance of the same geometry. (The port's own reading of
    its zero-thickness lead is not a reference: it reads ~7 % above the line on the calibration case.)"""
    t, z = np.array(res.get("tdr_t_ps") or []), np.array(res.get("tdr_z_ohm") or [])
    ref, what = (res.get("z_line_ohm"), "the model's calibrated line") if res.get("z_line_ohm") else (res.get("z_2d_ohm"), "2D")
    win = window_ps(model)
    if not len(t) or not ref or not win:
        return None
    lo, hi = win
    m = (t > lo) & (t < hi)
    if not m.any():
        return None
    return {"min_pct": round(float(z[m].min() / ref - 1) * 100, 2), "max_pct": round(float(z[m].max() / ref - 1) * 100, 2),
            "ref_ohm": round(float(ref), 2), "ref": what, "window_ps": [round(lo, 1), round(hi, 1)]}


def reflection_floor(res):
    """The reflection a uniform line shows in this model, in dB: the port lead's own reading against the
    S-parameters' reference (50 ohm a leg), plus the step from the zero-thickness lead to the board's
    track at each box edge (|z_port - z_line| / (z_port + z_line), twice), added in phase. On the
    validation lines this bounds the measured |S11| (-23.3 vs -23, -21.3 vs -17, -19.7 vs -18 dB): an
    S11 or Sdd11 near or below it is the model's, not the layout's."""
    zp, zt = res.get("z_port_ohm"), res.get("z_line_ohm") or res.get("z_2d_ohm")
    zr = {"pair": 100.0, "pair, common mode": 25.0}.get(res.get("line_kind"), 50.0)
    if not zp or not zt:
        return None
    g = abs(zp - zr) / (zp + zr) + 2 * abs(zp - zt) / (zp + zt)
    return round(20 * math.log10(max(g, 1e-6)), 1)


def plateau_um(fmax, eeff=None):
    """The uniform run a line needs inside the box for its TDR plateau (see _line_z) at fmax's rise time."""
    rise = tdr_rise_ps(fmax)
    tl = max(1.2 * rise + 20.0, (1.2 * rise) / 0.75)          # round trip, ps
    return tl * 1e-12 * C0 / (2 * math.sqrt(eeff or 3.0)) * 1e6


def verdict(res):
    """One line: what the region result is good for."""
    if not res.get("ran"):
        pf = res.get("preflight", {})
        if pf.get("ok"):
            return "checked only (not run)"
        return f"refused before solving: {'; '.join(pf.get('problems', [])[:2])}"
    conv, cal = res.get("convergence", {}), res.get("calibration", {})
    if res.get("stopped"):
        return f"no answer: {conv.get('reason')}"
    bad = []
    if not conv.get("converged"):
        bad.append(f"not converged ({conv.get('reason')})")
    if res.get("passive") is False:
        bad.append("not passive")
    if bad:
        return "not usable: " + "; ".join(bad)
    ex = res.get("tdr_excursion") or {}
    if res.get("clean"):
        c = res["clean"]
        return (f"pass: clean, below the model's resolution (TDR within {c['tdr_pct']:.1f} % "
                f"{'at the link edge' if c['on'] == 'net' else 'in the box'} on both grids, return loss {c['s11_db']} dB under "
                f"the model's {c['floor_db']} dB floor)")
    note = "" if cal.get("agree_pct") is not None else \
        f" (line not calibrated in this box: TDR against {ex.get('ref', '2D')}, {ex.get('ref_ohm', res.get('z_2d_ohm'))} ohm)"
    tdr = f"; TDR {ex['min_pct']:+.1f} / {ex['max_pct']:+.1f} %" if ex else ""
    fl = res.get("reflection_floor_db")
    fl_s = f"; S11/Sdd11 below {fl:.0f} dB is the model's own (port leads), not resolved" if fl is not None else ""
    return "usable: both grids agree" + tdr + fl_s + note


def _two_d(ctx, model):
    from . import xsec
    p = excited_port(model)
    if p is None:
        return {}
    b = ctx.b
    li = b.L[p["layer"]]
    refs = [b.layers[r] for r in ctx.ref_layers(li)]
    grp = next(g for g in model["groups"] if any(q["n"] == p["n"] for q in g["ports"]))
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
