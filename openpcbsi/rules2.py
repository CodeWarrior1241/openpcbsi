"""Board-wide stitching map (#23, RP3) and return proximity on boards without planes (#27).

stitching_map: for the ground planes of the board (layers with reference copper over at least 20 % of
the outline), the distance from every point where two planes overlap to the nearest via that joins
them; its 50 / 90 / 99 % percentiles, and the largest via-free area (the largest connected region
farther than `void_mm` from any stitching via). Warns when the 90th percentile exceeds `p90_mm` or a
void exceeds `void_warn_mm2`, naming the fast nets that cross it.

return_proximity: tracks of signal nets with no reference plane under them (two-layer boards, planeless
areas) are judged by their distance to the nearest reference copper on the same or an adjacent layer,
which carries their return current. Warns where a run longer than `min_len_mm` stays farther than
`max_mm` from it.
"""
import math

import numpy as np
from scipy import ndimage

from .drc import finding, result


def _down(m, k):
    """Downsample a boolean raster by k (any True in the block)."""
    ny, nx = m.shape
    ny2, nx2 = ny // k, nx // k
    return m[:ny2 * k, :nx2 * k].reshape(ny2, k, nx2, k).any(axis=(1, 3))


def stitching_map(ctx):
    b, cfg = ctx.b, ctx.cfg
    p = cfg["rules"]["stitching_map"]
    g = b.grid()
    k = max(1, int(round(p["cell_mm"] / g.res)))
    cell = g.res * k
    planes = [li for li in range(len(b.layers)) if ctx.refmask(li).mean() > 0.2 * _outline_frac(ctx)]
    metrics, finds = {}, []
    if len(planes) < 2:
        return result("stitching_map", "RP3", "Stitching map between ground planes", p, metrics, finds,
                      {"skipped": "fewer than two ground planes"})
    both = None
    for li in planes:
        m = _down(ctx.refmask(li), k)
        both = m if both is None else both + m.astype(np.int16)
    region = both >= 2                                   # at least two planes overlap here
    ref = set(ctx.ref_ids.tolist())
    V = b.vias
    via_mask = np.zeros_like(region)
    n_vias = 0
    for i in range(len(V["x"])):
        if int(V["net"][i]) not in ref:
            continue
        fl = [li for li in planes if li in b.via_flashed[i] or _spans(b, i, li)]
        if len(fl) < 2:
            continue
        ii, jj = g.ij(V["x"][i], V["y"][i])
        ii, jj = int(ii) // k, int(jj) // k
        if 0 <= ii < via_mask.shape[0] and 0 <= jj < via_mask.shape[1]:
            via_mask[ii, jj] = True
            n_vias += 1
    if not n_vias or not region.any():
        return result("stitching_map", "RP3", "Stitching map between ground planes", p, metrics, finds,
                      {"skipped": "no via joins two ground planes"})
    dist = ndimage.distance_transform_edt(~via_mask) * cell
    vals = dist[region]
    p50, p90, p99 = (float(np.percentile(vals, q)) for q in (50, 90, 99))
    lab, n = ndimage.label(region & (dist > p["void_mm"]))
    void = {"area_mm2": 0.0}
    if n:
        sizes = ndimage.sum(np.ones_like(lab), lab, index=range(1, n + 1))
        kk = int(np.argmax(sizes)) + 1
        ys, xs = np.nonzero(lab == kk)
        x0, y0 = g.x0 + (xs.min() * k) * g.res, g.y0 + (ys.min() * k) * g.res
        x1, y1 = g.x0 + ((xs.max() + 1) * k) * g.res, g.y0 + ((ys.max() + 1) * k) * g.res
        cx, cy = float(np.mean(xs)) * cell + g.x0, float(np.mean(ys)) * cell + g.y0
        # nets of fast classes with tracks crossing the void
        fast, T = set(), b.tracks
        for t in range(len(T["x1"])):
            net = b.net_names[T["net"][t]]
            if not ctx.net_class(net):
                continue
            for u in (0.0, 0.5, 1.0):
                x = T["x1"][t] + u * (T["x2"][t] - T["x1"][t])
                y = T["y1"][t] + u * (T["y2"][t] - T["y1"][t])
                i2, j2 = g.ij(x, y)
                i2, j2 = int(i2) // k, int(j2) // k
                if 0 <= i2 < lab.shape[0] and 0 <= j2 < lab.shape[1] and lab[i2, j2] == kk:
                    fast.add(net)
        void = {"area_mm2": round(float(sizes[kk - 1]) * cell * cell, 1), "bbox": [round(v, 2) for v in (x0, y0, x1, y1)],
                "x": round(cx, 2), "y": round(cy, 2), "fast_nets": sorted(fast)[:20]}
    metrics["_map"] = {"planes": [b.layers[li] for li in planes], "stitching_vias": n_vias,
                       "p50_mm": round(p50, 2), "p90_mm": round(p90, 2), "p99_mm": round(p99, 2), "void": void}
    if p90 > p["p90_mm"]:
        finds.append(finding("stitching_map", "warning", f"90 % of the plane area is within {p90:.1f} mm of a stitching via "
                             f"(limit {p['p90_mm']} mm); median {p50:.1f} mm", "", "", void.get("x"), void.get("y"),
                             value=round(p90, 2), limit=p["p90_mm"], unit="mm"))
    if void["area_mm2"] > p["void_warn_mm2"]:
        finds.append(finding("stitching_map", "warning", f"{void['area_mm2']:.0f} mm2 more than {p['void_mm']} mm from any stitching via"
                             + (f", crossed by {', '.join(void['fast_nets'][:6])}" if void.get("fast_nets") else ""),
                             (void.get("fast_nets") or [""])[0], "", void.get("x"), void.get("y"),
                             value=void["area_mm2"], limit=p["void_warn_mm2"], unit="mm2"))
    return result("stitching_map", "RP3", "Stitching map between ground planes", p, metrics, finds,
                  {"p50_mm": round(p50, 2), "p90_mm": round(p90, 2), "void_area_mm2": void["area_mm2"]})


def _outline_frac(ctx):
    """Share of the raster inside the board outline (planes are judged against the board, not the bbox)."""
    if not hasattr(ctx, "_outline_frac"):
        from . import geom
        g = ctx.b.grid()
        m = np.zeros(g.shape, bool)
        geom.raster_rings(g, ctx.b.outline_rings(), m)
        ctx._outline_frac = max(float(m.mean()), 0.05)
    return ctx._outline_frac


def _spans(b, i, li):
    span = b.snap["vias"][i].get("span") or []
    if len(span) == 2:
        lo, hi = sorted(b.L[x] for x in span)
        return lo <= li <= hi
    return False


def return_proximity(ctx):
    b, cfg = ctx.b, ctx.cfg
    p = cfg["rules"]["return_proximity"]
    g = b.grid()
    k = max(1, int(round(0.1 / g.res)))
    cell = g.res * k
    T = b.tracks
    nets = ctx.nets("return_proximity")
    if not nets:                                       # no classes: every signal net
        ref = set(cfg.reference_nets)
        nets = [n for n in b.net_names[1:] if n not in ref and b.kicad_class.get(n) != "PWR"
                and not n.startswith("unconnected-")]
    dmap = {}
    metrics, finds = {}, []
    for net in nets:
        worst, run_len, run_start = 0.0, 0.0, None
        for t in b.tracks_of(net):
            li = int(T["layer"][t])
            refs = ctx.ref_layers(li)
            planned = any(ctx.refmask(r).mean() > 0.2 * _outline_frac(ctx) for r in refs)
            if planned:
                continue                                # a plane is there: the plane rules judge it
            if li not in dmap:
                near = [x for x in (li - 1, li, li + 1) if 0 <= x < len(b.layers)]
                m = None
                for x in near:
                    mm = _down(ctx.refmask(x), k)
                    m = mm if m is None else (m | mm)
                dmap[li] = ndimage.distance_transform_edt(~m) * cell if m.any() else None
            L = float(T["len"][t])
            if dmap[li] is None:
                d = math.inf
            else:
                ds = []
                for u in (0.0, 0.25, 0.5, 0.75, 1.0):
                    x = T["x1"][t] + u * (T["x2"][t] - T["x1"][t])
                    y = T["y1"][t] + u * (T["y2"][t] - T["y1"][t])
                    i, j = g.ij(x, y)
                    i, j = min(max(int(i) // k, 0), dmap[li].shape[0] - 1), min(max(int(j) // k, 0), dmap[li].shape[1] - 1)
                    ds.append(dmap[li][i, j])
                d = float(np.median(ds))
            worst = max(worst, d)
            if d > p["max_mm"]:
                run_len += L
                run_start = run_start or (float(T["x1"][t]), float(T["y1"][t]), b.layers[li])
        if run_start is None and worst == 0.0:
            continue
        metrics[net] = {"max_distance_mm": None if math.isinf(worst) else round(worst, 2), "far_mm": round(run_len, 2)}
        if run_len > p["min_len_mm"]:
            finds.append(finding("return_proximity", "warning",
                                 f"{run_len:.1f} mm with no plane and the nearest reference copper "
                                 + ("absent" if math.isinf(worst) else f"{worst:.1f} mm away") + f" (limit {p['max_mm']} mm)",
                                 net, run_start[2], run_start[0], run_start[1], value=None if math.isinf(worst) else round(worst, 2),
                                 limit=p["max_mm"], unit="mm"))
    return result("return_proximity", "RP", "Return path on boards without planes", p, metrics, finds)


def breakout(ctx):
    """#24: breakout necks inside the fields (BGA and connector areas). For each pair (or single net) of a
    signal class: the tracks inside a field that are narrower than the net's width outside it, their
    length and their 2D impedance (with the local pair gap, or uncoupled), judged against the class
    target. These findings are about the fields, so they count (they are not set aside as area findings)."""
    from collections import defaultdict
    from . import xsec
    from .planner import pair_geometries
    b, cfg = ctx.b, ctx.cfg
    p = cfg["rules"]["breakout"]
    T = b.tracks
    metrics, finds = {}, []
    where = {}
    pair_geometries(ctx, where)
    gap_of = {}
    for (li, w, g), wh in where.items():
        for n in wh["nets"]:
            gap_of.setdefault((n, li, round(w, 4)), g)
    xc = cfg["xsec"]
    done = set()
    for net in ctx.nets("breakout"):
        if net in done:
            continue
        partner = ctx.partner.get(net)
        pos = min(net, partner) if partner else net
        done |= {net, partner} if partner else {net}
        cls = cfg.class_of(pos)
        ks = b.tracks_of(pos)
        if not len(ks):
            continue
        mx, my = (T["x1"][ks] + T["x2"][ks]) / 2, (T["y1"][ks] + T["y2"][ks]) / 2
        codes = ctx.area_at(mx, my)
        inside = codes > 0 if hasattr(codes, "__gt__") else np.array([c > 0 for c in codes])
        out_w = T["w"][ks][~inside]
        nominal = float(np.median(out_w)) if len(out_w) else float(np.max(T["w"][ks]))
        necks = defaultdict(lambda: {"len": 0.0, "x": None, "y": None, "area": ""})
        for k, ins, x, y, c in zip(ks, inside, mx, my, codes):
            w = float(T["w"][k])
            if not ins or w >= nominal - 1e-4:
                continue
            key = (int(T["layer"][k]), round(w, 4))
            nk = necks[key]
            nk["len"] += float(T["len"][k])
            if nk["x"] is None:
                nk["x"], nk["y"], nk["area"] = float(x), float(y), ctx.area_name(c)
        for (li, w), nk in necks.items():
            if nk["len"] < p["neck_min_mm"]:
                continue
            lname = b.layers[li]
            refs = [b.layers[r] for r in ctx.ref_layers(li)]
            g = gap_of.get((pos, li, w)) if partner else None
            try:
                if partner and g:
                    sec = xsec.build_section(b.stack, lname, refs, [w, w], [g], xc, xc["mask_er"])
                    z = xsec.solve_lines(sec, xc["res_um"], (1e9,))["freqs"][0]["Zdiff"]
                else:
                    sec = xsec.build_section(b.stack, lname, refs, [w], [], xc, xc["mask_er"])
                    z = xsec.solve_lines(sec, xc["res_um"], (1e9,))["freqs"][0]["Z0"] * (2 if partner else 1)
            except Exception:
                continue
            target = cfg.class_param(cls, "z_diff" if partner else "z_se") if cls else None
            tol = cfg.class_param(cls, "z_tol", 10.0) if cls else 10.0
            key = f"{pos}{'/' + partner if partner else ''} {lname} w={w}"
            metrics[key] = {"neck_mm": round(nk["len"], 2), "w_mm": w, "gap_mm": g, "Z_ohm": round(z, 1),
                            "nominal_w_mm": round(nominal, 4), "field": nk["area"], "target": target}
            if target and abs(z - target) > tol:
                finds.append(finding("breakout", "warning",
                                     f"{nk['len']:.1f} mm neck of w={w} mm on {lname} in field {nk['area']}: "
                                     f"{'Zdiff' if partner else 'Z0'} {z:.0f} ohm against {target:g} +/- {tol:g} ohm"
                                     + ("" if g or not partner else " (tracks not coupled there)"),
                                     pos, lname, nk["x"], nk["y"], value=round(z, 1), limit=tol, unit="ohm",
                                     neck_mm=round(nk["len"], 2)))
    return result("breakout", "breakout", "Breakout necks inside the fields", p, metrics, finds)
