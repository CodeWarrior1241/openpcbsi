"""Copper geometry for the manufacturing review: spacing between nets (MF12), gaps in copper of one
net (MF2) and acute corners (MF1, a ranking).

Everything here is measurement of drawn copper. Nothing models the etching process.

  spacing   Between the copper of two nets, the gap is the sum of the two distances from a point on
            the border between their nearest-copper regions. Judged against the fabricator's
            published spacing (track to track, pad to track, and the figures it accepts in a BGA).
  pocket    A ball of radius r is rolled through the space between copper; what it cannot reach,
            where every wall is copper of one net, is a pocket. r is a fixed measuring scale
            (`ball_mm`), not a property of any process. A pocket of constant width is a gap,
            judged against the fabricator's same-net spacing where a track bounds it. A pocket that
            narrows to a point is a wedge: with E the distance from its tip to the nearest ball
            centre, the opening angle is 2 asin(r / E), which does not depend on r.
  corners   Wedges with an opening angle below `acute_deg` are listed, sharpest and thinnest track
            first. That is the drafting rule of thumb about acute angles and nothing more: no
            fabricator in the profile publishes a limit for it, so it is ranked, not judged.

Candidates are found on the board raster (25 um) and measured on a raster `hires_div` times finer,
drawn from the vector data in tiles.
"""
import math

import numpy as np
from scipy import ndimage

from . import geom, par
from .drc import finding, result

EDT = ndimage.distance_transform_edt
S8 = np.ones((3, 3), bool)
PRE = 1.45          # a pocket is measured when its tip is further than PRE x r from the nearest ball centre
SLOT_RATIO = 0.7    # width at depth / widest: a wedge narrows towards its tip (0.4), a slot does not (1.0)


def fields(lab):
    """Distance of every pixel to the nearest copper (pixels), that copper's index and its label."""
    cu = lab != 0
    D, ind = EDT(~cu, return_indices=True)
    return D, ind, lab[ind[0], ind[1]]


def pockets(lab, rp, want_masks=False, f=None):
    """Pockets of a net-label raster (0 = no copper) for a ball of radius `rp` pixels.

    Returns a list of dicts: i, j (tip, the pixel furthest from any ball centre), e (that distance,
    px), d_max and d_deep (largest distance to copper in the pocket, and in its deepest 40 %), area
    (px), net (label of the walls), sl (bounding slices) and, if asked, the label image."""
    cu = lab != 0
    if not cu.any() or cu.all():
        return ([], None, None) if want_masks else []
    D, ind, n1 = f if f is not None else fields(lab)
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


def close_pairs(lab, f, limit_px, smooth=7):
    """Places where copper of two nets is closer than `limit_px`: arrays i, j (pixel on the border),
    gap (px), the two labels, and the index (2, n) of the nearest copper pixel on each side.

    At two adjacent pixels whose nearest copper belongs to different nets, the gap is the distance
    between those two copper pixels less one pixel's extent in that direction. The readings scatter
    by the pixel grid, so each is averaged with its neighbours along the border."""
    D, ind, n1 = f
    S, W = np.zeros(lab.shape), np.zeros(lab.shape)
    recs = []
    for a, b in (((slice(None), slice(None, -1)), (slice(None), slice(1, None))),
                 ((slice(None, -1), slice(None)), (slice(1, None), slice(None)))):
        la, lb = n1[a], n1[b]                    # on copper, n1 is the pixel's own label and D is 0
        m = (la != lb) & (D[a] + D[b] < limit_px + 3)
        if not m.any():
            continue
        ii, jj = np.nonzero(m)
        ia, ib = ind[:, a[0], a[1]][:, m], ind[:, b[0], b[1]][:, m]
        dy, dx = (ia[0] - ib[0]).astype(float), (ia[1] - ib[1]).astype(float)
        d = np.hypot(dx, dy)
        # a copper pixel's centre lies up to one pixel extent inside the true edge (a full extent for an
        # edge along the grid, less for an oblique one): 0.7 of it is taken off, which keeps the reading
        # within one pixel of the drawn gap at any angle (validation cases)
        s = np.where(d <= 1.0, 0.0, np.maximum(d - 0.7 * (np.abs(dx) + np.abs(dy)) / np.maximum(d, 1e-9), 0.0))
        np.add.at(S, (ii, jj), s)
        np.add.at(W, (ii, jj), 1.0)
        recs.append((ii, jj, la[m], lb[m], ia, ib))
    if not recs:
        z = np.zeros(0, int)
        return z, z, np.zeros(0), z, z, np.zeros((2, 0), int), np.zeros((2, 0), int)
    Ss, Ws = ndimage.uniform_filter(S, smooth), ndimage.uniform_filter(W, smooth)
    ii, jj, la, lb, ia, ib = (np.concatenate([r[k] for r in recs], axis=-1) for k in range(6))
    val = Ss[ii, jj] / np.maximum(Ws[ii, jj], 1e-12)
    keep = val < limit_px
    return ii[keep], jj[keep], val[keep], la[keep], lb[keep], ia[:, keep], ib[:, keep]


def _tiles(shape, size, halo):
    ny, nx = shape
    for i in range(0, ny, size):
        for j in range(0, nx, size):
            yield (i, min(i + size, ny), j, min(j + size, nx),
                   max(i - halo, 0), min(i + size + halo, ny), max(j - halo, 0), min(j + size + halo, nx))


def candidates(b, li, rp_px, space_px=None, tile=1024):
    """Stage 1 on the board raster: pockets worth measuring (tip and extent, mm) and, with `space_px`,
    points where two nets come closer than that (x, y, estimated gap in mm)."""
    r = b.raster(li)
    g = b.grid()
    halo = int(math.ceil(6 * rp_px))
    pk, sp = [], []
    for i0, i1, j0, j1, a0, a1, c0, c1 in _tiles(r.shape, tile, halo):
        lab = r[a0:a1, c0:c1]
        if not lab.any():
            continue
        f = fields(lab)
        for p in pockets(lab, rp_px, f=f):
            i, j = p["i"] + a0, p["j"] + c0
            if not (i0 <= i < i1 and j0 <= j < j1):
                continue
            sl = p["sl"]
            x, y = g.xy(i, j)
            pk.append({"x": float(x), "y": float(y),
                       "x0": g.x0 + (sl[1].start + c0) * g.res, "x1": g.x0 + (sl[1].stop + c0) * g.res,
                       "y0": g.y0 + (sl[0].start + a0) * g.res, "y1": g.y0 + (sl[0].stop + a0) * g.res})
        if space_px is not None:
            ii, jj, s = close_pairs(lab, f, space_px)[:3]
            ii, jj = ii + a0, jj + c0
            core = (ii >= i0) & (ii < i1) & (jj >= j0) & (jj < j1)
            if core.any():
                x, y = g.xy(ii[core], jj[core])
                sp.append(np.c_[x, y, s[core] * g.res])
    return pk, (np.vstack(sp) if sp else np.zeros((0, 3)))


def measure(b, li, cands, rp_mm, res, tile_mm):
    """Stage 2 for pockets: the tiles that hold a candidate, redrawn at `res` and measured."""
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
        for p in pockets(lab, rp):
            x, y = g.xy(p["i"], p["j"])
            if not (x0 <= x < x0 + tile_mm and y0 <= y < y0 + tile_mm):
                continue
            kind, val = classify(p, rp)
            sl = p["sl"]
            out.append({"kind": kind, "x": float(x), "y": float(y), "net": p["net"],
                        "angle_deg": val if kind == "wedge" else None,
                        "width_mm": val * res if kind == "gap" else None,
                        "area_mm2": p["area"] * res * res,
                        "enclosed": not p["open"] and not p["cut"], "cut": p["cut"],
                        "bbox": (g.x0 + sl[1].start * res, g.y0 + sl[0].start * res,
                                 g.x0 + sl[1].stop * res, g.y0 + sl[0].stop * res)})
    return out


def measure_spacing(b, li, pts, limit_of, res, tile_mm, max_tiles, coarse_res, reach_mm=0.2):
    """Stage 2 for spacing. `pts` are stage-1 points (x, y, estimated gap); `limit_of(x, y, pad_a,
    pad_b)` gives (limit mm, kind, limit key) or None when the pair is not judged here. The tiles
    with the smallest estimates are measured first; beyond `max_tiles` the estimate is kept.

    Returns ({(net a, net b, kind): {gap_mm, x, y, n, limit, key, res}}, tiles measured, tiles left)."""
    halo = reach_mm + 0.15
    by = {}
    for x, y, s in pts:
        k = (int(math.floor(x / tile_mm)), int(math.floor(y / tile_mm)))
        by[k] = min(by.get(k, 1e9), s)
    order = sorted(by, key=by.get)
    fine, rest = order[:max_tiles], set(order[max_tiles:])
    found = {}

    def add(na, nb, kind, gap, x, y, limit, key, r):
        k = (min(na, nb), max(na, nb), kind)
        e = found.setdefault(k, {"gap_mm": 1e9, "n": 0})
        e["n"] += 1
        if gap < e["gap_mm"]:
            e.update(gap_mm=gap, x=float(x), y=float(y), limit=limit, key=key, res=r)

    for tx, ty in fine:
        x0, y0 = tx * tile_mm, ty * tile_mm
        g, lab, pm = b.raster_window(li, x0 - halo, y0 - halo, x0 + tile_mm + halo, y0 + tile_mm + halo, res, pads=True)
        ii, jj, s, la, lb, ia, ib = close_pairs(lab, fields(lab), reach_mm / res)
        if not len(ii):
            continue
        x, y = g.xy(ii, jj)
        core = (x >= x0) & (x < x0 + tile_mm) & (y >= y0) & (y < y0 + tile_mm)
        pa, pb = pm[ia[0], ia[1]], pm[ib[0], ib[1]]
        for k in np.flatnonzero(core):
            lim = limit_of(x[k], y[k], bool(pa[k]), bool(pb[k]))
            if lim is None:
                continue
            gap = s[k] * res
            if gap < lim[0] + res:                  # one fine pixel either side of the limit is the reading's slack
                add(int(la[k]), int(lb[k]), lim[1], gap, x[k], y[k], lim[0], lim[2], res)
    if rest:                                            # not refined: only what is clear at the coarse resolution
        for x, y, s in pts:
            if (int(math.floor(x / tile_mm)), int(math.floor(y / tile_mm))) not in rest:
                continue
            lim = limit_of(x, y, False, False)
            if lim is not None and s < lim[0] - 1.5 * coarse_res:
                add(0, 0, lim[1] + " (coarse)", s, x, y, lim[0], lim[2], coarse_res)
    return found, len(fine), len(rest)


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
        out.append(dict(worst, width_mm=w, area_mm2=area, length_mm=area / max(w, 1e-6),
                        enclosed=all(g["enclosed"] for g in gs), cut=any(g["cut"] for g in gs), pieces=len(gs)))
    return out


def track_at(b, li, net, x, y, reach):
    """Width of the narrowest track of `net` on layer `li` whose edge passes within `reach` of a point
    (the track that forms a corner or bounds a gap there), or None if only pads and zones do."""
    T = b.tracks
    sel = np.flatnonzero((T["layer"] == li) & (T["net"] == net))
    if not len(sel):
        return None
    d = geom.point_seg_dist(x, y, T["x1"][sel], T["y1"][sel], T["x2"][sel], T["y2"][sel])
    near = np.abs(d - T["w"][sel] / 2) <= reach
    return float(T["w"][sel][near].min()) if near.any() else None


def run(ctx, fabr, in_bga, want=("gap", "acid", "space"), log=print):
    """Spacing rule, same-net gap rule and the list of acute corners, over every copper layer.
    `in_bga(x, y)` says whether a point is inside a BGA field (the fabricator's local figures apply).
    Returns {"fab_space": result, "same_net_gap": result, "acid_trap": ranking}."""
    b, cfg = ctx.b, ctx.cfg
    p = cfg["dfm"]["acid"]
    sp_cfg = cfg["dfm"]["space"]
    g = b.grid()
    snet = fabr.judge("same_net_spacing_mm", "same_net_gap")
    rp_mm = max(p["ball_mm"], (snet or 0.0) / 2)
    res = g.res / max(1, int(p["hires_div"]))
    judge_at = p["acute_deg"] - p["angle_tol_deg"]
    space = fabr.judge("track_space_mm", "fab_space")
    padtrk = fabr.judge("pad_to_track_mm", "fab_space")
    bga_sp = fabr.judge("bga_track_mm", "fab_space")
    bga_pt = fabr.judge("bga_pad_to_track_mm", "fab_space")
    do_space = "space" in want and (space is not None or padtrk is not None)

    cur = {"space": space, "key": "track_space_mm"}         # the track spacing of the layer being measured

    def limit_of(x, y, pad_a, pad_b):
        space = cur["space"]
        bga = in_bga(x, y)
        if pad_a and pad_b:
            return None                                    # pad to pad: MF4
        if pad_a or pad_b:
            if padtrk is None:
                return None
            return (bga_pt, "pad to track in a BGA field", "bga_pad_to_track_mm") if bga and bga_pt is not None else \
                   (padtrk, "pad to track", "pad_to_track_mm")
        if space is None:
            return None
        return (bga_sp, "track spacing in a BGA field", "bga_track_mm") if bga and bga_sp is not None else \
               (space, "track spacing", cur["key"])

    widest = max([v for v in (space, padtrk) if v is not None] +
                 [e["value"] for k, e in fabr.limits.items() if k.startswith("track_space_mm@")])  if do_space else None
    space_px = (widest + 1.5 * g.res) / g.res if do_space else None
    def layer(li):
        lname = b.layers[li]
        cur["key"] = fabr.layer_key("track_space_mm", lname)
        cur["space"] = fabr.judge(cur["key"], "fab_space") if space is not None else None
        f_gap, f_sp, corners, counts = [], [], [], {"on_track": 0, "wide_copper": 0}
        m_sp = None
        cands, pts = candidates(b, li, rp_mm / g.res, space_px)
        found = measure(b, li, cands, rp_mm, res, p["tile_mm"]) if cands and ("gap" in want or "acid" in want) else []
        # ---- acute corners (ranking)
        nw = 0
        if "acid" in want:
            for f in found:
                if f["kind"] != "wedge" or f["angle_deg"] >= judge_at:
                    continue
                nw += 1
                tw = track_at(b, li, f["net"], f["x"], f["y"], 0.04)
                counts["on_track" if tw is not None else "wide_copper"] += 1
                corners.append({"layer": lname, "x": round(f["x"], 3), "y": round(f["y"], 3), "net": b.net_names[f["net"]],
                                "angle_deg": round(f["angle_deg"], 1), "track_mm": tw, "bga": bool(in_bga(f["x"], f["y"]))})
        # ---- gaps in copper of one net
        gaps = []
        n_track = 0
        if "gap" in want and snet is not None:
            gaps = _merge_gaps([f for f in found if f["kind"] == "gap"], 2 * res)
            gaps = [f for f in gaps if f["length_mm"] >= p["min_gap_len_mm"] and f["width_mm"] < snet - res]
            for f in gaps:
                tw = track_at(b, li, f["net"], f["x"], f["y"], f["width_mm"] + 0.06)
                what = "enclosed opening" if f["enclosed"] else "gap"
                msg = (f"{what} of {f['width_mm'] * 1000:.0f} um in copper of one net, {'at least ' if f['cut'] else ''}"
                       f"{f['length_mm']:.2f} mm long")
                if tw is not None:
                    n_track += 1
                    msg += f", beside a {tw:g} mm track (same-net track spacing {snet:g} mm, {fabr.src('same_net_spacing_mm')})"
                else:
                    msg += (f", between pads or pour copper only: the fabricator's {snet:g} mm figure is for tracks, "
                            f"so this is listed, not judged")
                f_gap.append(finding("same_net_gap", "warning" if tw is not None else "info", msg, b.net_names[f["net"]], lname,
                                     f["x"], f["y"], value=round(f["width_mm"], 4), limit=snet, unit="mm",
                                     length_mm=round(f["length_mm"], 3), kind="track" if tw is not None else "copper",
                                     basis=fabr.src("same_net_spacing_mm"), field="BGA" if in_bga(f["x"], f["y"]) else ""))
        m_gap = {"gaps": len(gaps), "beside_a_track": n_track,
                        "min_width_mm": round(min(f["width_mm"] for f in gaps), 4) if gaps else None}
        # ---- spacing between nets
        n_fine = n_left = n_below = 0
        if do_space:
            fs, n_fine, n_left = measure_spacing(b, li, pts, limit_of, res, p["tile_mm"], sp_cfg["max_tiles"], g.res,
                                                 max(v for v in (cur["space"], padtrk) if v is not None) + 0.02)
            for (na, nb, kind), e in sorted(fs.items(), key=lambda kv: kv[1]["gap_mm"] - kv[1]["limit"]):
                names = f"{b.net_names[na]} / {b.net_names[nb]}" if na or nb else "nets not resolved"
                below = e["gap_mm"] < e["limit"] - e["res"]
                n_below += below
                how = (f"below the {e['limit'] * 1000:.0f} um minimum" if below else
                       f"within the measuring slack ({e['res'] * 1000:g} um) of the {e['limit'] * 1000:.0f} um minimum, so neither a pass nor a fail can be shown")
                f_sp.append(finding("fab_space", "error" if below else "warning", f"{names}: {e['gap_mm'] * 1000:.0f} um apart, {how} "
                                    f"for {kind.replace(' (coarse)', '')} ({fabr.src(e['key'])}); "
                                    f"{e['n'] * e['res']:.2f} mm of it, measured at {e['res'] * 1000:g} um",
                                    b.net_names[na] if na else "", lname, e["x"], e["y"], value=round(e["gap_mm"], 4), limit=e["limit"],
                                    unit="mm", aggressor=b.net_names[nb] if nb else "", kind=kind, basis=fabr.src(e["key"])))
            est = pts[:, 2].min() if len(pts) else None
            m_sp = {"min_gap_estimate_mm": None if est is None else round(float(est), 3), "tiles_measured": n_fine,
                           "tiles_not_refined": n_left, "below_limit": n_below}
        line = (f"  {lname:8s} {len(cands):5d} pockets  {nw:4d} acute corners  {len(gaps):4d} same-net gaps"
                + (f"  {n_below:4d} spacings below the limit ({n_fine} tiles measured"
                   + (f", {n_left} at 25 um only" if n_left else "") + ")" if do_space else ""))
        return {"f_gap": f_gap, "f_sp": f_sp, "corners": corners, "counts": counts, "m_gap": m_gap, "m_sp": m_sp, "line": line}

    for l in range(len(b.layers)):                 # loaded before the workers start, so they share them
        b.raster(l)
    b._window_index()
    f_gap, f_sp, m_gap, m_sp, corners = [], [], {}, {}, []
    counts = {"on_track": 0, "wide_copper": 0}
    for li, o in enumerate(par.pmap(layer, range(len(b.layers)), par.workers(cfg))):
        lname = b.layers[li]
        f_gap += o["f_gap"]
        f_sp += o["f_sp"]
        corners += o["corners"]
        for k in counts:
            counts[k] += o["counts"][k]
        m_gap[lname] = o["m_gap"]
        if o["m_sp"] is not None:
            m_sp[lname] = o["m_sp"]
        log(o["line"])
    corners.sort(key=lambda c: (c["track_mm"] is None, c["track_mm"] or 0.0, c["angle_deg"]))
    params = {"ball_mm": rp_mm, "res_mm": res, "same_net_spacing_mm": snet, "min_gap_len_mm": p["min_gap_len_mm"]}
    out = {}
    if "gap" in want:
        out["same_net_gap"] = result("same_net_gap", "MF2", "Gaps in copper of one net below the same-net track spacing", params, m_gap, f_gap)
    if "space" in want:
        out["fab_space"] = result("fab_space", "MF12", "Copper spacing between nets, measured", {
            "track_space_mm": space, "pad_to_track_mm": padtrk, "bga_track_mm": bga_sp, "bga_pad_to_track_mm": bga_pt,
            "res_mm": res, "tolerance_mm": round(res, 5)}, m_sp, f_sp)
    if "acid" in want:
        out["acid_trap"] = {"rule": "acid_trap", "spec": "MF1", "title": "Acute corners in copper of one net",
                            "params": {"acute_deg": p["acute_deg"], "listed_below_deg": judge_at, "ball_mm": rp_mm, "res_mm": res},
                            "summary": dict(counts, corners=len(corners)), "rows": corners, "findings": []}
    return out


# ------------------------------------------------------------------------------------------------
def _wedge(angle_deg, rot_deg, n=420, apex=(210.3, 210.6)):
    """A label raster: copper of one net everywhere except a wedge of the given opening angle."""
    ii, jj = np.mgrid[0:n, 0:n]
    x, y = jj + 0.5 - apex[0], ii + 0.5 - apex[1]
    a = np.degrees(np.arctan2(y, x)) - rot_deg
    a = (a + 180.0) % 360.0 - 180.0
    return np.where(np.abs(a) < angle_deg / 2, 0, 1).astype(np.int32)


def validation_cases():
    """Measurement of drawn shapes: wedge angles, a slot's width, and the gap between two nets drawn
    square-on and at an angle."""
    from .validate import _case
    out = []
    rp = 20.0
    for ang, rot in ((30.0, 0.0), (45.0, 17.0), (60.0, 45.0), (75.0, 100.0), (90.0, 33.0)):
        lab = _wedge(ang, rot)
        lab[np.hypot(*(np.mgrid[0:420, 0:420] - 210.0)) > 200] = 0      # the wedge opens into free space
        ps = [p for p in pockets(lab, rp) if abs(p["i"] - 210) < 12 and abs(p["j"] - 210) < 12]
        val = classify(ps[0], rp)[1] if ps and classify(ps[0], rp)[0] == "wedge" else 0.0
        out.append(_case(f"wedge of {ang:g} deg drawn at {rot:g} deg: opening angle", "etch", val, ang, 6.0, "deg",
                         "drawn angle (2 asin(r / E) is exact for straight edges)"))
    n, w, L = 520, 16, 300
    lab = np.ones((n, n), np.int32)
    lab[:, :100] = 0
    lab[252:252 + w, 100:100 + L] = 0
    ps = pockets(lab, rp)
    slot = max(ps, key=lambda p: p["area"]) if ps else None
    kind, val = classify(slot, rp) if slot else ("", 0.0)
    out.append(_case("dead-end slot 16 px wide: width", "etch", val if kind == "gap" else 0.0, float(w), 7.0, "px", "drawn width"))
    for rot, name in ((0.0, "square-on"), (17.0, "at 17 deg"), (30.0, "at 30 deg"), (45.0, "at 45 deg"), (68.0, "at 68 deg")):
        ii, jj = np.mgrid[0:300, 0:300]
        d = (jj - 150.0) * math.cos(math.radians(rot)) + (ii - 150.0) * math.sin(math.radians(rot))
        gap_w = 14.4
        lab = np.where(d < -gap_w / 2, 1, np.where(d > gap_w / 2, 2, 0)).astype(np.int32)
        s = close_pairs(lab, fields(lab), 40.0)[2]
        out.append(_case(f"gap of 14.4 px between two nets, {name}: narrowest reading", "etch", float(s.min()), gap_w, 7.0, "px",
                         "drawn gap, to one pixel"))
    return out
