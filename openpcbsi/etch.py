"""MF1 / MF2: places where etchant is held between copper of one net (acid traps, narrow same-net gaps).

Geometry. A ball of radius r (the reach of the etchant flow; by default half the fabricator's
same-net spacing) is rolled through the space between copper. What it cannot reach is a pocket.
Copper of a different net never closes a pocket: such gaps are ordinary clearance (KiCad's DRC), so a
pocket counts only where every wall is copper of one net. Two kinds are reported:

  wedge  two edges of one net meeting at an angle. With E the distance from the wedge tip to the
         nearest ball centre, the opening angle is 2 asin(r / E). Capillary criterion (Concus and
         Finn, 1969): liquid with contact angle theta is held in a wedge whose half-angle alpha
         satisfies alpha + theta < 90 deg, i.e. opening angle < 180 - 2 theta.
  gap    a slot or channel narrower than 2 r between copper of one net (constant width).

Flush time. For each pocket the mean time for a dissolved ion to diffuse out to the flowing etchant
is solved on the pocket: D lap(tau) = -1, tau = 0 at the mouth, no flux through the copper walls.
A dead-end slot of length L gives tau_max = L^2 / 2D; a sector of radius R gives R^2 / 4D. Exchange
through the top of the resist channel is ignored, so tau ranks pockets against each other; it is not a
prediction of the fabricator's process.

The candidates are found on the board raster (25 um), then measured on a raster `hires_div` times
finer, drawn from the vector data in 2 mm tiles.
"""
import math

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy import ndimage

from .drc import finding, result

EDT = ndimage.distance_transform_edt
S8 = np.ones((3, 3), bool)
PRE = 1.45          # a pocket is measured when its tip is further than PRE x r from the nearest ball centre
SLOT_RATIO = 0.7    # width at depth / widest: a wedge narrows towards its tip (0.4), a slot does not (1.0)


def pockets(lab, rp, want_masks=False):
    """Pockets of a net-label raster (0 = no copper) for a ball of radius `rp` pixels.

    Returns a list of dicts: i, j (tip, the pixel furthest from any ball centre), e (that distance,
    px), d_max and d_deep (largest distance to copper in the pocket, and in its deepest 40 %), area
    (px), net (label of the walls), sl (bounding slices) and, if asked, the label image."""
    cu = lab != 0
    if not cu.any() or cu.all():
        return ([], None, None) if want_masks else []
    D, ind = EDT(~cu, return_indices=True)
    n1 = lab[ind[0], ind[1]]
    reach = D >= rp + 0.5
    E = EDT(~reach) if reach.any() else np.full(lab.shape, 1e6)
    pocket = ~cu & (E > rp)
    # a pixel near the border between the nearest-copper regions of two nets has a wall of another net
    vb = np.zeros(lab.shape, bool)
    gap = ~cu
    h = (n1[:, 1:] != n1[:, :-1]) & gap[:, 1:] & gap[:, :-1]
    v = (n1[1:, :] != n1[:-1, :]) & gap[1:, :] & gap[:-1, :]
    vb[:, 1:] |= h
    vb[:, :-1] |= h
    vb[1:, :] |= v
    vb[:-1, :] |= v
    if vb.any():
        pocket &= EDT(~vb) > rp
    lbl, n = ndimage.label(pocket, structure=S8)
    if not n:
        return ([], lbl, None) if want_masks else []
    index = np.arange(1, n + 1)
    e_max = ndimage.maximum(E, lbl, index)
    keep = np.flatnonzero(e_max >= PRE * rp)
    sls = ndimage.find_objects(lbl)
    out = []
    for k in keep:
        sl = sls[k]
        m = lbl[sl] == k + 1
        e, d = E[sl][m], D[sl][m]
        ia = int(np.argmax(e))
        ii, jj = np.nonzero(m)
        deep = e >= rp + 0.6 * (e[ia] - rp)
        out.append({"i": sl[0].start + int(ii[ia]), "j": sl[1].start + int(jj[ia]), "e": float(e[ia]),
                    "d_max": float(d.max()), "d_deep": float(d[deep].max()), "area": int(m.sum()),
                    "net": int(n1[sl][m][ia]), "sl": sl, "label": k + 1,
                    "open": bool(_touches_open(lbl, sl, k + 1, gap)),
                    # cut by the border of the raster: the pocket goes on outside it
                    "cut": sl[0].start == 0 or sl[1].start == 0 or sl[0].stop == lab.shape[0] or sl[1].stop == lab.shape[1]})
    return (out, lbl, gap) if want_masks else out


def _touches_open(lbl, sl, k, gap):
    """Whether a pocket has a mouth: a neighbouring gap pixel that is not part of it."""
    i0, i1 = max(sl[0].start - 1, 0), min(sl[0].stop + 1, lbl.shape[0])
    j0, j1 = max(sl[1].start - 1, 0), min(sl[1].stop + 1, lbl.shape[1])
    m = lbl[i0:i1, j0:j1] == k
    ring = ndimage.binary_dilation(m, structure=S8) & ~m
    return (ring & gap[i0:i1, j0:j1]).any()


def classify(p, rp):
    """('gap', width_px) for a constant-width slot, ('wedge', opening angle in degrees) otherwise."""
    if p["d_max"] > 0 and p["d_deep"] / p["d_max"] >= SLOT_RATIO:
        return "gap", max(2.0 * p["d_deep"] - 0.5, 1.0)
    s = min(1.0, (rp + 0.5) / max(p["e"], rp + 0.5))
    return "wedge", 2.0 * math.degrees(math.asin(s))


def flush_time(lbl, gap, p, res_mm, diff):
    """Largest and mean diffusion time out of a pocket (s); None if the pocket has no mouth."""
    sl = p["sl"]
    i0, i1 = max(sl[0].start - 1, 0), min(sl[0].stop + 1, lbl.shape[0])
    j0, j1 = max(sl[1].start - 1, 0), min(sl[1].stop + 1, lbl.shape[1])
    m = lbl[i0:i1, j0:j1] == p["label"]
    g = gap[i0:i1, j0:j1]
    n = int(m.sum())
    idx = -np.ones(m.shape, np.int64)
    idx[m] = np.arange(n)
    diag = np.zeros(n)
    rows, cols = [], []
    mouth = False
    for di, dj in ((0, 1), (1, 0)):
        a_m, b_m = m[: m.shape[0] - di, : m.shape[1] - dj], m[di:, dj:]
        a_g, b_g = g[: m.shape[0] - di, : m.shape[1] - dj], g[di:, dj:]
        a_i, b_i = idx[: m.shape[0] - di, : m.shape[1] - dj], idx[di:, dj:]
        both = a_m & b_m
        rows += [a_i[both], b_i[both]]
        cols += [b_i[both], a_i[both]]
        np.add.at(diag, a_i[both], 1.0)
        np.add.at(diag, b_i[both], 1.0)
        for inside, other_gap, ii in ((a_m & ~b_m, b_g, a_i), (b_m & ~a_m, a_g, b_i)):
            d = inside & other_gap                # a gap pixel outside the pocket: the mouth (tau = 0)
            if d.any():
                mouth = True
                np.add.at(diag, ii[d], 1.0)
    if not mouth:
        return None, None
    rows, cols = np.concatenate(rows), np.concatenate(cols)
    A = sp.coo_matrix((-np.ones(len(rows)), (rows, cols)), shape=(n, n)).tocsr() + sp.diags(diag)
    lone = diag == 0                              # a pixel with no neighbour at all (cannot happen with a mouth)
    if lone.any():
        A = A + sp.diags(lone.astype(float))
    tau = spla.spsolve(A.tocsc(), np.full(n, (res_mm * 1e-3) ** 2 / diff))
    return float(tau.max()), float(tau.mean())


def _tiles(shape, size, halo):
    ny, nx = shape
    for i in range(0, ny, size):
        for j in range(0, nx, size):
            yield (i, min(i + size, ny), j, min(j + size, nx),
                   max(i - halo, 0), min(i + size + halo, ny), max(j - halo, 0), min(j + size + halo, nx))


def candidates(b, li, rp_px, tile=1024):
    """Stage 1 on the board raster: tip positions (mm) and extents of pockets worth measuring."""
    r = b.raster(li)
    g = b.grid()
    halo = int(math.ceil(6 * rp_px))
    out = []
    for i0, i1, j0, j1, a0, a1, c0, c1 in _tiles(r.shape, tile, halo):
        lab = r[a0:a1, c0:c1]
        if not lab.any():
            continue
        for p in pockets(lab, rp_px):
            i, j = p["i"] + a0, p["j"] + c0
            if not (i0 <= i < i1 and j0 <= j < j1):
                continue
            sl = p["sl"]
            x, y = g.xy(i, j)
            out.append({"x": float(x), "y": float(y),
                        "x0": g.x0 + (sl[1].start + c0) * g.res, "x1": g.x0 + (sl[1].stop + c0) * g.res,
                        "y0": g.y0 + (sl[0].start + a0) * g.res, "y1": g.y0 + (sl[0].stop + a0) * g.res})
    return out


def measure(b, li, cands, rp_mm, res, tile_mm, diff):
    """Stage 2: the tiles that hold a candidate, redrawn at `res` and measured."""
    halo = 6 * rp_mm
    rp = rp_mm / res
    tiles = set()
    for c in cands:
        for tx in range(int(math.floor(c["x0"] / tile_mm)), int(math.floor(c["x1"] / tile_mm)) + 1):
            for ty in range(int(math.floor(c["y0"] / tile_mm)), int(math.floor(c["y1"] / tile_mm)) + 1):
                tiles.add((tx, ty))
    out = []
    for tx, ty in sorted(tiles):
        x0, y0 = tx * tile_mm, ty * tile_mm
        g, lab = b.raster_window(li, x0 - halo, y0 - halo, x0 + tile_mm + halo, y0 + tile_mm + halo, res)
        ps, lbl, gap = pockets(lab, rp, want_masks=True)
        for p in ps:
            x, y = g.xy(p["i"], p["j"])
            if not (x0 <= x < x0 + tile_mm and y0 <= y < y0 + tile_mm):
                continue
            kind, val = classify(p, rp)
            tmax, tmean = flush_time(lbl, gap, p, res, diff) if p["open"] and not p["cut"] else (None, None)
            sl = p["sl"]
            out.append({"kind": kind, "x": float(x), "y": float(y), "net": p["net"],
                        "angle_deg": val if kind == "wedge" else None,
                        "width_mm": val * res if kind == "gap" else None,
                        "depth_mm": max(p["e"] - rp, 0.0) * res, "area_mm2": p["area"] * res * res,
                        "tau_max_s": tmax, "tau_mean_s": tmean, "enclosed": not p["open"] and not p["cut"], "cut": p["cut"],
                        "bbox": (g.x0 + sl[1].start * res, g.y0 + sl[0].start * res,
                                 g.x0 + sl[1].stop * res, g.y0 + sl[0].stop * res)})
    return out


def _merge_gaps(gaps, tol):
    """Pieces of one gap cut by the tile borders: join pieces of one net whose boxes touch."""
    n = len(gaps)
    parent = list(range(n))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a
    order = sorted(range(n), key=lambda k: gaps[k]["bbox"][0])
    for a_, ia in enumerate(order):
        A = gaps[ia]["bbox"]
        for ib in order[a_ + 1:]:
            B = gaps[ib]["bbox"]
            if B[0] > A[2] + tol:
                break
            if gaps[ia]["net"] == gaps[ib]["net"] and B[1] <= A[3] + tol and B[3] >= A[1] - tol:
                parent[find(ia)] = find(ib)
    groups = {}
    for k in range(n):
        groups.setdefault(find(k), []).append(gaps[k])
    out = []
    for gs in groups.values():
        w = min(g["width_mm"] for g in gs)
        area = sum(g["area_mm2"] for g in gs)
        worst = min(gs, key=lambda g: g["width_mm"])
        taus = [g["tau_max_s"] for g in gs if g["tau_max_s"] is not None]
        out.append(dict(worst, width_mm=w, area_mm2=area, length_mm=area / max(w, 1e-6),
                        tau_max_s=max(taus) if taus and len(gs) == 1 else None,
                        enclosed=all(g["enclosed"] for g in gs), cut=any(g["cut"] for g in gs), pieces=len(gs)))
    return out


def run(ctx, fabr, log=print):
    """Both rules over every copper layer. Returns [acid_trap result, same_net_gap result]."""
    b, cfg = ctx.b, ctx.cfg
    p = cfg["dfm"]["acid"]
    g = b.grid()
    snet = fabr.get("same_net_spacing_mm")
    rp_mm = p["ball_mm"] or (snet / 2 if snet else 0.125)
    res = g.res / max(1, int(p["hires_div"]))
    cf = 180.0 - 2.0 * p["contact_angle_deg"]
    limit = min(cf, p["max_angle_deg"])
    judge_at = limit - p["angle_tol_deg"]
    areas = ctx.areas()
    f_acid, f_gap = [], []
    m_acid, m_gap = {}, {}
    for li, lname in enumerate(b.layers):
        cands = candidates(b, li, rp_mm / g.res)
        found = measure(b, li, cands, rp_mm, res, p["tile_mm"], p["diffusivity_m2s"]) if cands else []
        wedges = [f for f in found if f["kind"] == "wedge" and f["angle_deg"] < judge_at]
        gaps = _merge_gaps([f for f in found if f["kind"] == "gap"], 2 * res)
        gaps = [f for f in gaps if f["length_mm"] >= p["min_gap_len_mm"] and
                (snet is None or f["width_mm"] < snet - res)]
        for f in wedges:
            net = b.net_names[f["net"]]
            area = ctx.area_name(g.sample(areas, [f["x"]], [f["y"]])[0])
            tau = f["tau_max_s"]
            msg = (f"acute copper corner of {f['angle_deg']:.0f} deg (limit {limit:.0f} deg: liquid with a "
                   f"{p['contact_angle_deg']:g} deg contact angle is held below it), {f['depth_mm'] * 1000:.0f} um deep"
                   + (f", flush time {tau:.2f} s" if tau is not None else ""))
            slow = tau is None or tau >= p["flush_warn_s"]
            if not slow:
                msg += f" (below the {p['flush_warn_s']:g} s from which a corner warns)"
            f_acid.append(finding("acid_trap", "warning" if slow else "info", msg, net, lname, f["x"], f["y"], value=round(f["angle_deg"], 1),
                                  limit=limit, unit="deg", area=area, depth_mm=round(f["depth_mm"], 4),
                                  area_mm2=round(f["area_mm2"], 5), tau_s=None if tau is None else round(tau, 3)))
        for f in gaps:
            net = b.net_names[f["net"]]
            area = ctx.area_name(g.sample(areas, [f["x"]], [f["y"]])[0])
            what = "enclosed opening" if f["enclosed"] else "gap"
            msg = (f"{what} of {f['width_mm'] * 1000:.0f} um in copper of one net, {'at least ' if f['cut'] else ''}{f['length_mm']:.2f} mm long"
                   + (f" (same-net spacing {snet:g} mm, {fabr.src('same_net_spacing_mm')})" if snet else ""))
            f_gap.append(finding("same_net_gap", "warning", msg, net, lname, f["x"], f["y"], value=round(f["width_mm"], 4),
                                 limit=snet, unit="mm", area=area, length_mm=round(f["length_mm"], 3),
                                 kind="enclosed" if f["enclosed"] else "gap"))
        taus = [f["tau_max_s"] for f in wedges if f["tau_max_s"] is not None]
        m_acid[lname] = {"candidates": len(cands), "wedges": len(wedges),
                         "slow": sum(1 for f in wedges if f["tau_max_s"] is None or f["tau_max_s"] >= p["flush_warn_s"]),
                         "min_angle_deg": round(min((f["angle_deg"] for f in wedges), default=0.0), 1) if wedges else None,
                         "max_flush_s": round(max(taus), 3) if taus else None}
        m_gap[lname] = {"gaps": len(gaps), "min_width_mm": round(min((f["width_mm"] for f in gaps), default=0.0), 4) if gaps else None}
        log(f"  {lname:8s} {len(cands):5d} candidates  {len(wedges):4d} acute corners  {len(gaps):4d} narrow same-net gaps")
    params = {"ball_mm": rp_mm, "res_mm": res, "limit_deg": limit, "judged_below_deg": judge_at, "flush_warn_s": p["flush_warn_s"],
              "contact_angle_deg": p["contact_angle_deg"], "diffusivity_m2s": p["diffusivity_m2s"],
              "same_net_spacing_mm": snet}
    return [result("acid_trap", "MF1", "Acute corners in copper of one net (acid traps)", params, m_acid, f_acid),
            result("same_net_gap", "MF2", "Gaps narrower than the same-net spacing", params, m_gap, f_gap)]


# ------------------------------------------------------------------------------------------------
def _wedge(angle_deg, rot_deg, n=420, apex=(210.3, 210.6)):
    """A label raster: copper of one net everywhere except a wedge of the given opening angle."""
    ii, jj = np.mgrid[0:n, 0:n]
    x, y = jj + 0.5 - apex[0], ii + 0.5 - apex[1]
    a = np.degrees(np.arctan2(y, x)) - rot_deg
    a = (a + 180.0) % 360.0 - 180.0
    return np.where(np.abs(a) < angle_deg / 2, 0, 1).astype(np.int32)


def validation_cases():
    """Opening angles of drawn wedges, a slot's width, and the two closed-form flush times."""
    from .validate import _case
    out = []
    rp = 20.0
    for ang, rot in ((30.0, 0.0), (45.0, 17.0), (60.0, 45.0), (75.0, 100.0), (90.0, 33.0)):
        lab = _wedge(ang, rot)
        lab[np.hypot(*(np.mgrid[0:420, 0:420] - 210.0)) > 200] = 0      # the wedge opens into free space
        ps = [p for p in pockets(lab, rp) if abs(p["i"] - 210) < 12 and abs(p["j"] - 210) < 12]
        val = classify(ps[0], rp)[1] if ps and classify(ps[0], rp)[0] == "wedge" else 0.0
        out.append(_case(f"wedge of {ang:g} deg drawn at {rot:g} deg: opening angle", "etch", val, ang, 6.0, "deg",
                         "2 asin(r / E), exact for straight edges"))
    # dead-end slot, 16 px wide and 300 px long, in copper of one net; diffusion time L^2 / 2D
    n, w, L = 520, 16, 300
    lab = np.ones((n, n), np.int32)
    lab[:, :100] = 0
    lab[252:252 + w, 100:100 + L] = 0
    ps, lbl, gap = pockets(lab, rp, want_masks=True)
    slot = max(ps, key=lambda p: p["area"]) if ps else None
    kind, val = classify(slot, rp) if slot else ("", 0.0)
    out.append(_case("dead-end slot 16 px wide: width", "etch", val if kind == "gap" else 0.0, float(w), 7.0, "px",
                     "drawn width"))
    res, diff = 0.00625, 1e-9
    tmax, _ = flush_time(lbl, gap, slot, res, diff) if slot else (0.0, 0.0)
    ref = (L * res * 1e-3) ** 2 / (2 * diff)
    out.append(_case("dead-end slot 300 px long: flush time", "etch", tmax or 0.0, ref, 5.0, "s",
                     "L^2 / 2D (one-dimensional diffusion to the mouth)"))
    # a quarter-disc pocket handed to the solver directly: tau_max = R^2 / 4D for any opening angle
    R = 120
    ii, jj = np.mgrid[0:R + 4, 0:R + 4]
    rr = np.hypot(ii + 0.5, jj + 0.5)
    lblq = (rr < R).astype(np.int32)
    gapq = np.ones(lblq.shape, bool)
    pq = {"sl": (slice(0, R + 4), slice(0, R + 4)), "label": 1}
    tq, _ = flush_time(lblq, gapq, pq, res, diff)
    out.append(_case("sector of radius 120 px: flush time", "etch", tq, (R * res * 1e-3) ** 2 / (4 * diff), 3.0, "s",
                     "R^2 / 4D (radial diffusion in a sector with insulating walls)"))
    return out
