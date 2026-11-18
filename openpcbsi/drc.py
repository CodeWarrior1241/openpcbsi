"""G1: electrical rule checks over every net and layer (HyperLynx DRC class).

Each rule returns {"rule", "spec", "title", "params", "metrics": {key: {...}}, "findings": [...]}.
Metrics are what G15 compares between two boards; findings are what G16 draws on the board.
Findings inside exclusion areas (BGA and connector fields) carry the area name and are counted
separately; they are listed, not failed.
"""
import math
from collections import defaultdict

import numpy as np
from scipy import ndimage
from scipy.spatial import cKDTree

from . import geom
from .context import runs

C0 = 299792458.0


def finding(rule, sev, msg, net="", layer="", x=None, y=None, x2=None, y2=None, value=None, limit=None,
            unit="", area="", **extra):
    f = {"rule": rule, "severity": sev, "msg": msg, "net": net, "layer": layer, "x": _r(x), "y": _r(y),
         "x2": _r(x2), "y2": _r(y2), "value": _r(value), "limit": limit, "unit": unit, "area": area}
    f.update(extra)
    return f


def _r(v, n=4):
    if v is None:
        return None
    return round(float(v), n)


def result(rule, spec, title, params, metrics, findings, summary=None):
    return {"rule": rule, "spec": spec, "title": title, "params": params, "metrics": metrics,
            "findings": findings, "summary": summary or {}}


# ==============================================================================================
# PC1 / PC2 / unreferenced: reference-plane continuity, split and gap crossings

PAD_VOID_MARGIN_MM = 0.15   # the void pad_plane asks for under a wide pad: the pad + 0.15 mm (its fix hint)


def _wide_pads(ctx, net):
    """The pads of `net` that pad_plane judges: SMD and connector pads wider than `ratio` x the track that
    enters them. [(pad, track width, pad width)]"""
    b, T = ctx.b, ctx.b.tracks
    p = ctx.cfg["rules"]["pad_plane"]
    ks = b.tracks_of(net)
    out = []
    for pad in b.pads_of(net):
        if pad["kind"] not in ("smd", "conn") or not pad["layer_ids"]:
            continue
        li = pad["layer_ids"][0]
        near = [k for k in ks if T["layer"][k] == li and min(
            math.hypot(T["x1"][k] - pad["x"], T["y1"][k] - pad["y"]),
            math.hypot(T["x2"][k] - pad["x"], T["y2"][k] - pad["y"])) <= max(pad["sx"], pad["sy"])]
        tw = float(np.median(T["w"][near])) if near else (float(np.median(T["w"][ks])) if len(ks) else 0)
        pw = min(pad["sx"], pad["sy"])
        if not tw or pw <= p["ratio"] * tw:
            continue
        out.append((pad, tw, pw))
    return out


def _own_pad_voids(ctx, net, li, pad_nets):
    """The voids ref_plane does not count against `net`'s track on layer `li`: under each of the net's own
    pads on that layer for which pad_plane asks for a void (a wide pad on a class it judges, Nyquist at or
    above its fast_ghz), the pad's outline grown by the margin pad_plane recommends plus one raster cell
    (the fill edge's quantisation). [(outer ring as an (n, 2) array, grow distance mm)]"""
    p = ctx.cfg["rules"]["pad_plane"]
    if not p.get("enabled", True) or net not in pad_nets:
        return []
    nyq = ctx.cparam(net, "nyquist_ghz", None)
    if nyq is None or nyq < p.get("fast_ghz", 3.0):
        return []
    d = p.get("void_margin_mm", PAD_VOID_MARGIN_MM) + ctx.g.res
    out = []
    for pad, tw, pw in _wide_pads(ctx, net):
        if li not in pad["layer_ids"]:
            continue
        for poly in pad["poly"] or []:
            if poly and len(poly[0]) >= 3:
                out.append((np.asarray(poly[0], float), d))
    return out


def _in_grown(px, py, rings):
    """Points (any shape) inside a ring grown by its distance, for [(ring, d)] (pads: small rings)."""
    hit = np.zeros(np.shape(px), bool)
    for ring, d in rings:
        x0, y0 = ring.min(axis=0) - d
        x1, y1 = ring.max(axis=0) + d
        cand = ~hit & (px >= x0) & (px <= x1) & (py >= y0) & (py <= y1)
        if not cand.any():
            continue
        x, y = px[cand], py[cand]
        inside = np.zeros(x.shape, bool)
        dmin = np.full(x.shape, np.inf)
        for (xa, ya), (xb, yb) in zip(ring, np.roll(ring, -1, axis=0)):
            if yb != ya:
                inside ^= ((ya > y) != (yb > y)) & (x < xa + (y - ya) * (xb - xa) / (yb - ya))
            ex, ey = xb - xa, yb - ya
            L2 = ex * ex + ey * ey
            t = np.clip(((x - xa) * ex + (y - ya) * ey) / L2, 0, 1) if L2 > 0 else np.zeros(x.shape)
            dmin = np.minimum(dmin, np.hypot(x - xa - t * ex, y - ya - t * ey))
        hit[cand] = inside | (dmin <= d)
    return hit


def rule_ref_plane(ctx):
    b, cfg, g = ctx.b, ctx.cfg, ctx.g
    p = cfg["rules"]["ref_plane"]
    ps = cfg["rules"]["split_crossing"]
    step = p["step"]
    metrics, finds = {}, []
    split_metrics, split_finds = {}, []
    pad_nets = set(ctx.nets("pad_plane"))
    for net in ctx.nets("ref_plane"):
        m = {"uncovered_mm": 0.0, "uncovered_in_areas_mm": 0.0, "length_mm": 0.0, "excused_own_pad_void_mm": 0.0}
        sm = {"gap_crossings": 0, "split_crossings": 0, "gap_crossings_in_areas": 0, "split_crossings_in_areas": 0}
        for li in sorted(set(b.tracks["layer"][b.tracks_of(net)].tolist())):
            refs = ctx.ref_layers(li)
            # pad_plane asks for the reference to be voided under the net's own wide pads: the track that
            # enters such a pad runs over that void by design, so it is not failed for it here (any other
            # void, and any part of this one beyond the pad + the recommended margin, still is)
            own_voids = _own_pad_voids(ctx, net, li, pad_nets)
            for path in ctx.path_samples(net, li, step):
                x, y, w, dl = path["x"], path["y"], path["w"], path["dl"]
                m["length_mm"] += path["length"]
                # direction per sample for the across-track offsets
                dx = np.gradient(x) if len(x) > 1 else np.array([1.0])
                dy = np.gradient(y) if len(y) > 1 else np.array([0.0])
                nn = np.hypot(dx, dy) + 1e-12
                nx_, ny_ = -dy / nn, dx / nn
                area = ctx.area_at(x, y)
                for rli in refs:
                    h = b.stack.gap(b.layers[li], b.layers[rli])
                    half = w / 2 + p["margin_h"] * h
                    nacross = 9
                    u = np.linspace(-1, 1, nacross)
                    px = x[:, None] + nx_[:, None] * half[:, None] * u[None, :]
                    py = y[:, None] + ny_[:, None] * half[:, None] * u[None, :]
                    rm = ctx.refmask(rli)
                    cov = g.sample(rm.view(np.uint8), px, py).astype(bool)
                    excused = np.zeros(cov.shape, bool)
                    if own_voids:
                        excused = ~cov & _in_grown(px, py, own_voids)
                        m["excused_own_pad_void_mm"] += float(dl[excused.any(axis=1)].sum())
                    centre_real = cov[:, nacross // 2]
                    cov = cov | excused
                    full = cov.all(axis=1)
                    centre = cov[:, nacross // 2]
                    rast = b.raster(rli)
                    under = g.sample(rast, x, y)
                    bad = ~full
                    out_len = float((dl[bad & (area == 0)]).sum())
                    in_len = float((dl[bad & (area != 0)]).sum())
                    m["uncovered_mm"] += out_len
                    m["uncovered_in_areas_mm"] += in_len
                    for a0, a1 in runs(bad):
                        L = float(dl[a0:a1].sum())
                        if L < 2 * step:
                            continue
                        kind = "margin" if centre[a0:a1].all() else (
                            "split" if (under[a0:a1][~centre[a0:a1]] != 0).any() else "void")
                        ar = ctx.area_name(np.bincount(area[a0:a1]).argmax())
                        sev = "info" if ar else ("warning" if kind == "margin" else "error")
                        finds.append(finding("ref_plane", sev,
                                             f"{L:.2f} mm without solid reference on {b.layers[rli]} ({kind})",
                                             net, b.layers[li], x[a0], y[a0], x[a1 - 1], y[a1 - 1], L,
                                             p["limit_mm"], "mm", ar, ref_layer=b.layers[rli], kind=kind))
                    # PC2: crossings of the centreline: ref -> (void|other) -> ref. Over the void pad_plane asks
                    # for under the net's own pad (state 3) the track is not on its plane: a run that ends
                    # there ends off the plane; a run of nothing else is no crossing
                    state = np.where(centre_real, 1, np.where(excused[:, nacross // 2], 3, np.where(under != 0, 2, 0)))
                    for a0, a1 in runs(state != 1):
                        if a0 == 0 or a1 >= len(state):
                            continue          # starts or ends off the plane: not a crossing
                        real = state[a0:a1] != 3
                        L = float(dl[a0:a1][real].sum())
                        if not real.any() or L < ps["min_gap_mm"]:
                            continue
                        kind = "split" if (state[a0:a1] == 2).any() else "gap"
                        ar = ctx.area_name(area[(a0 + a1) // 2])
                        key = f"{kind}_crossings" + ("_in_areas" if ar else "")
                        sm[key] += 1
                        split_finds.append(finding("split_crossing", "info" if ar else "error",
                                                   f"crosses a {L:.2f} mm {kind} in {b.layers[rli]}",
                                                   net, b.layers[li], x[(a0 + a1) // 2], y[(a0 + a1) // 2],
                                                   value=L, unit="mm", area=ar, ref_layer=b.layers[rli], kind=kind))
        metrics[net] = {k: round(v, 3) for k, v in m.items()}
        split_metrics[net] = sm
    return [result("ref_plane", "PC1", "Reference-plane continuity under the track (+ margin)",
                   p, metrics, finds),
            result("split_crossing", "PC2", "Crossings of gaps and splits in the reference plane",
                   ps, split_metrics, split_finds)]


def rule_unreferenced(ctx):
    b, cfg, g = ctx.b, ctx.cfg, ctx.g
    p = cfg["rules"]["unreferenced"]
    metrics, finds = {}, []
    step = 0.1
    T = b.tracks
    for net in ctx.nets("unreferenced"):
        ks = b.tracks_of(net)
        if not len(ks):
            continue
        tot = 0.0
        for k in ks:
            li = int(T["layer"][k])
            refs = [b.L[n] for n in b.stack.neighbours(b.layers[li])]
            x, y, dl = geom.sample_segment(T["x1"][k], T["y1"][k], T["x2"][k], T["y2"][k], step)
            ok = np.zeros(len(x), bool)
            for rli in refs:
                ok |= g.sample(ctx.anyplane(rli).view(np.uint8), x, y).astype(bool)
            bad = ~ok
            for a0, a1 in runs(bad):
                L = (a1 - a0) * dl
                tot += L
                if L >= p["min_len_mm"]:
                    ar = ctx.area_name(ctx.area_at(x[(a0 + a1) // 2], y[(a0 + a1) // 2]))
                    finds.append(finding("unreferenced", "info" if ar else "warning",
                                         f"{L:.2f} mm with no plane on either adjacent layer", net, b.layers[li],
                                         x[a0], y[a0], x[a1 - 1], y[a1 - 1], L, p["min_len_mm"], "mm", ar))
        if tot > 0:
            metrics[net] = {"unreferenced_mm": round(tot, 3)}
    return result("unreferenced", "G1", "Track segments with no plane on either adjacent layer", p, metrics, finds)


# ==============================================================================================
# PC3 return vias, PC4 via stubs, PC9 stitching

def _ref_via_layers(ctx):
    """Indices of reference-net vias and the set of layers on which each joins its plane."""
    b = ctx.b
    vzl = ctx.via_zone_layers()
    idx = np.flatnonzero(np.isin(b.vias["net"], ctx.ref_ids))
    return idx, [vzl[k] for k in idx]


def rule_return_via(ctx):
    b, cfg = ctx.b, ctx.cfg
    p = cfg["rules"]["return_via"]
    vtl = ctx.via_track_layers()
    ridx, rlay = _ref_via_layers(ctx)
    rx, ry = b.vias["x"][ridx], b.vias["y"][ridx]
    # reference-net pads (ball pads, pins), allowed inside exclusion areas
    rpads = [q for q in b.pads if q["net_id"] in set(ctx.ref_ids.tolist())]
    px = np.array([q["x"] for q in rpads]) if rpads else np.zeros(0)
    py = np.array([q["y"] for q in rpads]) if rpads else np.zeros(0)
    # caps between a power net and a reference net, for vias referenced to a power area
    caps = _decaps(ctx)
    metrics, finds = {}, []
    for net in ctx.nets("return_via"):
        worst = 0.0
        nvias = 0
        for k in b.vias_of(net):
            layers = vtl[k]
            if len(layers) < 2:
                continue
            nvias += 1
            need = set()
            for li in layers:
                need |= set(ctx.ref_layers(li))
            x, y = b.vias["x"][k], b.vias["y"][k]
            # which needed reference layers are actually power copper at this via?
            power_ref = {}
            for rli in need:
                r = b.raster(rli)
                # look 0.5 mm around (the via's own antipad is a void)
                ang = np.linspace(0, 2 * np.pi, 12, endpoint=False)
                s = ctx.g.sample(r, x + 0.5 * np.cos(ang), y + 0.5 * np.sin(ang))
                s = s[s != 0]
                if len(s) and not np.isin(s, ctx.ref_ids).any():
                    power_ref[rli] = b.net_names[np.bincount(s).argmax()]
            need_gnd = need - set(power_ref)
            ok = np.array([need_gnd <= L for L in rlay], bool) if len(rlay) else np.zeros(0, bool)
            d = np.hypot(rx[ok] - x, ry[ok] - y)
            dist = float(d.min()) if len(d) else math.inf
            area = ctx.area_name(ctx.area_at([x], [y])[0])
            if area and p["count_pads_in_areas"] and len(px):
                dp = np.hypot(px - x, py - y)
                dist = min(dist, float(dp.min()))
            msg = f"nearest return via joining {', '.join(b.layers[i] for i in sorted(need_gnd))} at {dist:.2f} mm"
            for rli, pn in power_ref.items():
                cd = [math.hypot(c["x"] - x, c["y"] - y) for c in caps if c["power"] == pn]
                cdist = min(cd) if cd else math.inf
                msg += f"; {b.layers[rli]} is {pn} here, nearest {pn}-GND capacitor {cdist:.1f} mm"
            worst = max(worst, dist if math.isfinite(dist) else 99.0)
            if dist > p["limit_mm"]:
                finds.append(finding("return_via", "info" if area else "error", msg, net,
                                     "/".join(b.layers[i] for i in sorted(layers)), x, y, value=dist,
                                     limit=p["limit_mm"], unit="mm", area=area))
        if nvias:
            metrics[net] = {"layer_change_vias": nvias, "worst_return_mm": round(worst, 3)}
    return result("return_via", "PC3", "Return vias at layer changes", p, metrics, finds)


def rule_via_stub(ctx):
    b, cfg = ctx.b, ctx.cfg
    p = cfg["rules"]["via_stub"]
    vtl, vzl = ctx.via_track_layers(), ctx.via_zone_layers()
    metrics, finds = {}, []
    st = b.stack
    for net in ctx.nets("via_stub"):
        nyq = ctx.cparam(net, "nyquist_ghz", None)
        rows = []
        for k in b.vias_of(net):
            conn = vtl[k] | vzl[k]
            if not conn:
                continue
            top, bot = int(b.vias["top"][k]), int(b.vias["bot"][k])
            bd = (b.snap["vias"][k].get("backdrill") or {}) if k < len(b.snap["vias"]) else {}
            if bd.get("bottom_to") in b.L:          # #26: back-drilled from the bottom up to that layer
                bot = min(bot, b.L[bd["bottom_to"]])
            if bd.get("top_to") in b.L:
                top = max(top, b.L[bd["top_to"]])
            lo, hi = min(conn), max(conn)
            stub_top = st.barrel(b.layers[top], b.layers[lo]) if lo > top else 0.0
            stub_bot = st.barrel(b.layers[hi], b.layers[bot]) if bot > hi else 0.0
            L = max(stub_top, stub_bot)
            er = st.mean_er()
            f = C0 / (4 * L * 1e-3 * math.sqrt(er)) / 1e9 if L > 0 else math.inf
            rows.append({"x": round(b.vias["x"][k], 3), "y": round(b.vias["y"][k], 3),
                         "layers": [b.layers[i] for i in sorted(conn)], "stub_mm": round(L, 4),
                         "backdrilled": bool(bd),
                         "f_res_ghz": round(f, 2) if math.isfinite(f) else None})
            if nyq and f < p["factor"] * nyq:
                area = ctx.area_name(ctx.area_at([b.vias["x"][k]], [b.vias["y"][k]])[0])
                finds.append(finding("via_stub", "info" if area else "error",
                                     f"stub {L:.3f} mm, quarter-wave resonance {f:.1f} GHz < {p['factor']:.0f} x Nyquist ({nyq} GHz)",
                                     net, "/".join(b.layers[i] for i in sorted(conn)), b.vias["x"][k], b.vias["y"][k],
                                     value=f, limit=p["factor"] * nyq, unit="GHz", area=area))
        if rows:
            fr = [r["f_res_ghz"] for r in rows if r["f_res_ghz"]]
            metrics[net] = {"vias": len(rows), "max_stub_mm": max(r["stub_mm"] for r in rows),
                            "min_f_res_ghz": min(fr) if fr else None, "list": rows}
    return result("via_stub", "PC4", "Via stub length and quarter-wave resonance", p, metrics, finds)


def rule_stitching(ctx):
    b, cfg = ctx.b, ctx.cfg
    p = cfg["rules"]["stitching"]
    ridx, rlay = _ref_via_layers(ctx)
    stitch = np.array([len(L) >= 2 for L in rlay], bool)
    sx, sy = b.vias["x"][ridx][stitch], b.vias["y"][ridx][stitch]
    metrics, finds = {}, []
    tree = cKDTree(np.c_[sx, sy]) if len(sx) else None
    for net in ctx.nets("stitching"):
        worst = 0.0
        for li in sorted(set(b.tracks["layer"][b.tracks_of(net)].tolist())):
            for path in ctx.path_samples(net, li, 0.1):
                if tree is None or path["length"] < 1.0:
                    continue
                pts = np.c_[path["x"], path["y"]]
                ptree = cKDTree(pts)
                cand = tree.query_ball_point(pts, r=p["band_mm"])
                vs = sorted(set(v for c in cand for v in c))
                if vs:
                    _, nearest = ptree.query(np.c_[sx[vs], sy[vs]])
                    s = np.sort(path["s"][nearest])
                else:
                    s = np.zeros(0)
                s = np.r_[0.0, s, path["length"]]
                gaps = np.diff(s)
                j = int(np.argmax(gaps))
                gmax = float(gaps[j])
                # skip parts of the corridor inside exclusion areas
                mid = (s[j] + s[j + 1]) / 2
                i_mid = int(np.argmin(abs(path["s"] - mid)))
                area = ctx.area_name(ctx.area_at([path["x"][i_mid]], [path["y"][i_mid]])[0])
                worst = max(worst, gmax if not area else 0)
                if gmax > p["max_gap_mm"]:
                    finds.append(finding("stitching", "info" if area else "warning",
                                         f"{gmax:.2f} mm along the corridor without a stitching via within {p['band_mm']} mm",
                                         net, b.layers[li], path["x"][i_mid], path["y"][i_mid], value=gmax,
                                         limit=p["max_gap_mm"], unit="mm", area=area))
        metrics[net] = {"max_stitch_gap_mm": round(worst, 3)}
    # board edge
    edge_gaps = []
    for ring in b.outline_rings():
        ex, ey, es, L = geom.resample_ring(ring, 0.1)
        if tree is None:
            break
        etree = cKDTree(np.c_[ex, ey])
        d, near = etree.query(np.c_[sx, sy])
        sel = d <= p["edge_band_mm"]
        s = np.sort(es[near[sel]])
        if len(s) < 2:
            edge_gaps.append(L)
            continue
        gaps = np.diff(np.r_[s, s[0] + L])
        for j in np.flatnonzero(gaps > p["board_edge_gap_mm"]):
            mid = (s[j] + gaps[j] / 2) % L
            i = int(np.argmin(abs(es - mid)))
            finds.append(finding("stitching", "warning",
                                 f"board edge: {gaps[j]:.1f} mm without a stitching via within {p['edge_band_mm']} mm",
                                 "", "edge", ex[i], ey[i], value=gaps[j], limit=p["board_edge_gap_mm"], unit="mm"))
        edge_gaps.append(float(gaps.max()))
    if edge_gaps:                       # only when an edge stretch was measured (no planes: nothing to check)
        metrics["_board_edge"] = {"max_edge_gap_mm": round(max(edge_gaps), 2), "stitching_vias": int(stitch.sum())}
    return result("stitching", "PC9", "Ground stitching along corridors and the board edge", p, metrics, finds)


# ==============================================================================================
# PC8 distance to plane and board edges

def rule_edge_distance(ctx):
    b, cfg, g = ctx.b, ctx.cfg, ctx.g
    p = cfg["rules"]["edge_distance"]
    rings = b.outline_rings()
    oseg = []
    for r in rings:
        a = np.asarray(r)
        oseg.append(np.c_[a, np.roll(a, -1, axis=0)])
    oseg = np.vstack(oseg) if oseg else np.zeros((0, 4))
    # plane extent per layer: reference copper with small holes (antipads) filled
    plane_dist = {}

    def pdist(rli):
        if rli not in plane_dist:
            m = ctx.refmask(rli)
            filled = ndimage.binary_fill_holes(m)
            holes = filled & ~m
            lab, n = ndimage.label(holes)
            if n:
                sizes = ndimage.sum(np.ones_like(lab), lab, index=np.arange(1, n + 1)) * g.res ** 2
                big = np.flatnonzero(sizes > 4.0) + 1
                ext = m | (holes & ~np.isin(lab, big))
            else:
                ext = m
            plane_dist[rli] = (ndimage.distance_transform_edt(ext) * g.res).astype(np.float32)
        return plane_dist[rli]

    metrics, finds = {}, []
    T = b.tracks
    for net in ctx.nets("edge_distance"):
        min_edge, min_plane = math.inf, math.inf
        loc_e = loc_p = None
        for k in b.tracks_of(net):
            li = int(T["layer"][k])
            x, y, dl = geom.sample_segment(T["x1"][k], T["y1"][k], T["x2"][k], T["y2"][k], 0.05)
            area = ctx.area_at(x, y)
            out = area == 0
            if not out.any():
                continue
            x, y = x[out], y[out]
            if len(oseg):
                d = np.min(np.stack([geom.point_seg_dist(x, y, *s) for s in oseg]), axis=0) - T["w"][k] / 2
                i = int(np.argmin(d))
                if d[i] < min_edge:
                    min_edge, loc_e = float(d[i]), (x[i], y[i], b.layers[li])
            for rli in ctx.ref_layers(li):
                pd = g.sample(pdist(rli), x, y) - T["w"][k] / 2
                on = pd > 0
                if on.any():
                    i = int(np.argmin(np.where(on, pd, np.inf)))
                    if pd[i] < min_plane:
                        min_plane, loc_p = float(pd[i]), (x[i], y[i], b.layers[li], b.layers[rli])
        if math.isfinite(min_edge) or math.isfinite(min_plane):
            metrics[net] = {"min_board_edge_mm": round(min_edge, 3) if math.isfinite(min_edge) else None,
                            "min_plane_edge_mm": round(min_plane, 3) if math.isfinite(min_plane) else None}
        if loc_e and min_edge < p["board_edge_mm"]:
            finds.append(finding("edge_distance", "warning", f"{min_edge:.2f} mm from the board edge", net, loc_e[2],
                                 loc_e[0], loc_e[1], value=min_edge, limit=p["board_edge_mm"], unit="mm"))
        if loc_p and min_plane < p["plane_edge_mm"]:
            finds.append(finding("edge_distance", "warning",
                                 f"{min_plane:.2f} mm from the edge of its reference copper on {loc_p[3]}", net,
                                 loc_p[2], loc_p[0], loc_p[1], value=min_plane, limit=p["plane_edge_mm"], unit="mm"))
    return result("edge_distance", "PC8", "Distance to plane and board edges (outside exclusion areas)",
                  p, metrics, finds)


# ==============================================================================================
# PC6 pair coupling, PC7 parallel length against spacing

def _window(g, x, y, margin):
    i0, j0 = g.ij(np.min(x) - margin, np.min(y) - margin)
    i1, j1 = g.ij(np.max(x) + margin, np.max(y) + margin)
    return max(int(i0), 0), min(int(i1) + 1, g.ny), max(int(j0), 0), min(int(j1) + 1, g.nx)


def rule_pair_coupling(ctx):
    b, cfg, g = ctx.b, ctx.cfg, ctx.g
    p = cfg["rules"]["pair_coupling"]
    metrics, finds = {}, []
    T = b.tracks
    done = set()
    for net in ctx.nets("pair_coupling"):
        other = ctx.partner.get(net)
        if not other or net in done:
            continue
        done |= {net, other}
        pos, neg = sorted([net, other])
        lim, lim_src = ctx.limit("pair_coupling", "limit_mm", pos)
        max_gap, _ = ctx.limit("pair_coupling", "max_gap_mm", pos)
        lp = float(T["len"][b.tracks_of(pos)].sum())
        ln = float(T["len"][b.tracks_of(neg)].sum())
        unc_out = unc_in = 0.0
        gaps = []
        per_layer = defaultdict(list)
        for li in sorted(set(T["layer"][b.tracks_of(pos)].tolist())):
            r = b.raster(li)
            nid = b.net_id[neg]
            for path in ctx.path_samples(pos, li, 0.025):
                x, y, w, dl = path["x"], path["y"], path["w"], path["dl"]
                i0, i1, j0, j1 = _window(g, x, y, 1.0)
                sub = r[i0:i1, j0:j1] != nid
                if sub.all():
                    d = np.full(len(x), np.inf)
                else:
                    edt = ndimage.distance_transform_edt(sub) * g.res
                    ii, jj = g.ij(x, y)
                    d = edt[np.clip(ii - i0, 0, i1 - i0 - 1), np.clip(jj - j0, 0, j1 - j0 - 1)]
                gap = d - w / 2 - g.res / 2     # cell centre to nearest partner cell centre
                area = ctx.area_at(x, y)
                bad = gap > max_gap
                unc_out += float(dl[bad & (area == 0)].sum())
                unc_in += float(dl[bad & (area != 0)].sum())
                good = ~bad & (area == 0)
                if good.any():
                    gaps.append(gap[good])
                    per_layer[b.layers[li]].append(gap[good])
                for a0, a1 in runs(bad & (area == 0)):
                    L = float(dl[a0:a1].sum())
                    if L >= 0.2:
                        finds.append(finding("pair_coupling", "info", f"{L:.2f} mm uncoupled (gap > {max_gap} mm)",
                                             pos, b.layers[li], x[a0], y[a0], x[a1 - 1], y[a1 - 1], L, unit="mm"))
        gl = np.concatenate(gaps) if gaps else np.zeros(0)
        key = f"{pos}/{neg}"
        metrics[key] = {"len_p_mm": round(lp, 3), "len_n_mm": round(ln, 3), "skew_mm": round(abs(lp - ln), 3),
                        "uncoupled_mm": round(unc_out, 3), "uncoupled_in_areas_mm": round(unc_in, 3),
                        "gap_median_mm": round(float(np.median(gl)), 4) if len(gl) else None,
                        "gap_by_layer": {k: round(float(np.median(np.concatenate(v))), 4) for k, v in per_layer.items()}}
        metrics[key]["limit_mm"], metrics[key]["limit_source"] = lim, lim_src
        if unc_out > lim:
            finds.append(finding("pair_coupling", "error", f"{unc_out:.2f} mm uncoupled outside the fields (limit {lim:g} mm, {lim_src})",
                                 pos, "", value=unc_out, limit=lim, unit="mm", limit_source=lim_src))
    return result("pair_coupling", "PC6", "Pair coupling along the length (P-to-N gap)", p, metrics, finds)


def coupling_profile(ctx, net, max_spacing, step=0.025, exclude=()):
    """For each layer and centreline sample of `net`: edge spacing to the nearest other-net copper and
    which net that is. Reference nets and `exclude` are not aggressors."""
    b, g = ctx.b, ctx.g
    T = b.tracks
    skip = set(ctx.ref_ids.tolist()) | {b.net_id[net]} | {b.net_id[e] for e in exclude if e in b.net_id} | {0}
    out = []
    for li in sorted(set(T["layer"][b.tracks_of(net)].tolist())):
        r = b.raster(li)
        for path in ctx.path_samples(net, li, step):
            x, y, w, dl = path["x"], path["y"], path["w"], path["dl"]
            i0, i1, j0, j1 = _window(g, x, y, max_spacing + 0.5)
            sub = r[i0:i1, j0:j1]
            agg = ~np.isin(sub, np.array(sorted(skip), dtype=np.int16))
            if not agg.any():
                continue
            edt, (ri, rj) = ndimage.distance_transform_edt(~agg, return_indices=True)
            ii, jj = g.ij(x, y)
            ii = np.clip(ii - i0, 0, i1 - i0 - 1)
            jj = np.clip(jj - j0, 0, j1 - j0 - 1)
            d = edt[ii, jj] * g.res - w / 2 - g.res / 2
            who = sub[ri[ii, jj], rj[ii, jj]]
            out.append({"layer": li, "x": x, "y": y, "dl": dl, "w": w, "spacing": d, "agg": who,
                        "area": ctx.area_at(x, y)})
    return out


def rule_parallel(ctx):
    b, cfg = ctx.b, ctx.cfg
    p = cfg["rules"]["parallel"]
    metrics, finds = {}, []
    pair_nets = set(ctx.partner)
    for net in ctx.nets("parallel"):
        partner = ctx.partner.get(net)
        prof = coupling_profile(ctx, net, p["max_spacing_mm"], exclude=[partner] if partner else [])
        lim_any = ctx.cparam(net, "min_spacing_mm", p["nonpair_min_mm"])
        lim_pairs = ctx.cparam(net, "min_spacing_pairs_mm", lim_any)
        coupled = defaultdict(lambda: {"length_mm": 0.0, "min_spacing_mm": math.inf})
        viol = defaultdict(float)
        worst = math.inf
        for pr in prof:
            sel = (pr["spacing"] <= p["max_spacing_mm"]) & (pr["area"] == 0)
            for a in np.unique(pr["agg"][sel]):
                ms = sel & (pr["agg"] == a)
                an = b.net_names[a]
                c = coupled[an]
                c["length_mm"] += float(pr["dl"][ms].sum())
                c["min_spacing_mm"] = min(c["min_spacing_mm"], float(pr["spacing"][ms].min()))
                lim = lim_pairs if an in pair_nets else lim_any
                if lim:
                    mv = ms & (pr["spacing"] < lim)
                    if mv.any():
                        viol[an] += float(pr["dl"][mv].sum())
                        i = int(np.flatnonzero(mv)[np.argmin(pr["spacing"][mv])])
                        if viol[an] == float(pr["dl"][mv].sum()):   # first time: one marker per aggressor
                            finds.append(finding("parallel", "warning",
                                                 f"runs at {pr['spacing'][i]:.3f} mm from {an} (limit {lim} mm)",
                                                 net, b.layers[pr["layer"]], pr["x"][i], pr["y"][i],
                                                 value=float(pr["spacing"][i]), limit=lim, unit="mm", aggressor=an))
            if sel.any():
                worst = min(worst, float(pr["spacing"][sel].min()))
        top = sorted(coupled.items(), key=lambda kv: -kv[1]["length_mm"])[:8]
        metrics[net] = {"min_spacing_mm": round(worst, 4) if math.isfinite(worst) else None,
                        "violation_mm": round(sum(viol.values()), 3),
                        "coupled": {k: {"length_mm": round(v["length_mm"], 3), "min_spacing_mm": round(v["min_spacing_mm"], 4)}
                                    for k, v in top}}
        for an, L in viol.items():
            for f in finds:
                if f["net"] == net and f.get("aggressor") == an and f["rule"] == "parallel":
                    f["msg"] += f"; {L:.2f} mm below the limit"
                    f["length_mm"] = round(L, 3)
    return result("parallel", "PC7", "Distance to other nets and parallel length (outside the fields)",
                  p, metrics, finds)


# ==============================================================================================
# PC5 pads over planes, PC11 decoupling, PC13 thermal vias, power overlap, stubs

def _pad_mask(g, pad):
    rings = [ring for poly in pad["poly"] for ring in poly]
    if not rings:
        return None, None
    allp = np.vstack([np.asarray(r) for r in rings])
    i0, j0 = g.ij(allp[:, 0].min(), allp[:, 1].min())
    i1, j1 = g.ij(allp[:, 0].max(), allp[:, 1].max())
    m = np.zeros(g.shape, bool)
    geom.raster_rings(g, rings, m)
    sl = (slice(max(int(i0) - 1, 0), int(i1) + 2), slice(max(int(j0) - 1, 0), int(j1) + 2))
    return m[sl], sl


def rule_pad_plane(ctx):
    b, cfg, g = ctx.b, ctx.cfg, ctx.g
    p = cfg["rules"]["pad_plane"]
    metrics, finds = {}, []
    for net in ctx.nets("pad_plane"):
        rows = []
        for pad, tw, pw in _wide_pads(ctx, net):
            li = pad["layer_ids"][0]
            m, sl = _pad_mask(g, pad)
            if m is None or not m.any():
                continue
            covs = {}
            for rli in ctx.ref_layers(li):
                covs[b.layers[rli]] = round(float(ctx.refmask(rli)[sl][m].mean()), 3)
            state = "solid" if all(v > 0.95 for v in covs.values()) else (
                "cut-out" if all(v < 0.05 for v in covs.values()) else "partial")
            rows.append({"pad": f"{pad['ref']}.{pad['num']}", "pad_w": round(pw, 3), "track_w": round(tw, 3),
                         "ref_cover": covs, "state": state})
            # RP1: judged on fast classes (Nyquist >= fast_ghz): solid reference -> error, partial -> warning
            nyq = ctx.cparam(net, "nyquist_ghz", None)
            fast = nyq is not None and nyq >= p.get("fast_ghz", 3.0)
            sev = ("error" if state == "solid" else "warning" if state == "partial" else "info") if fast else "info"
            finds.append(finding("pad_plane", sev,
                                 f"pad {pad['ref']}.{pad['num']} {pw:.2f} mm wide on a {tw:.2f} mm track; reference {state} {covs}"
                                 + (f" (class Nyquist {nyq:g} GHz)" if fast else ""),
                                 net, b.layers[li], pad["x"], pad["y"], value=pw / tw, unit="x", pad=f"{pad['ref']}.{pad['num']}",
                                 state=state, area=ctx.area_name(ctx.area_at([pad["x"]], [pad["y"]])[0])))
        if rows:
            metrics[net] = {"wide_pads": len(rows), "solid": sum(r["state"] == "solid" for r in rows),
                            "cut_out": sum(r["state"] == "cut-out" for r in rows), "list": rows}
    return result("pad_plane", "PC5", "Pads wider than the track and the reference under them", p, metrics, finds)


def _decaps(ctx):
    if hasattr(ctx, "_decap_list"):
        return ctx._decap_list
    b = ctx.b
    p = ctx.cfg["rules"]["decap"]
    refset = set(ctx.cfg.reference_nets)
    out = []
    for ref, pads in b.pads_by_ref.items():
        if not any(ref.startswith(pf) and ref[len(pf):len(pf) + 1].isdigit() for pf in p["prefixes"]):
            continue
        cp = [q for q in pads if q["net"]]
        if len(cp) != 2:
            continue
        nets = {q["net"] for q in cp}
        if len(nets & refset) != 1 or len(nets) != 2:
            continue
        power = (nets - refset).pop()
        if not (power in ctx.zone_nets or b.kicad_class.get(power) == "PWR" or power in ctx.cfg["rails"]):
            continue
        f = b.fps[ref]
        out.append({"ref": ref, "value": f["value"], "fpid": f["fpid"], "power": power, "x": f["x"], "y": f["y"],
                    "side": f["side"], "pads": cp, "dnp": f["dnp"]})
    ctx._decap_list = out
    return out


def _pad_via_dist(pad, vx, vy, vd):
    rings = [ring for poly in pad["poly"] for ring in poly]
    if not rings or not len(vx):
        return np.full(len(vx), np.inf)
    inside = geom.points_in_rings(vx, vy, rings)
    d = np.full(len(vx), np.inf)
    for r in rings:
        a = np.asarray(r)
        for s in np.c_[a, np.roll(a, -1, axis=0)]:
            d = np.minimum(d, geom.point_seg_dist(vx, vy, *s))
    d = np.maximum(d - vd / 2, 0)
    d[inside] = 0.0
    return d


def _poured(ctx, pad):
    """True if copper of the pad's own net surrounds most of the pad on its layer (a pour or plane)."""
    if not pad["layer_ids"]:
        return False
    r = ctx.b.raster(pad["layer_ids"][0])
    ang = np.linspace(0, 2 * np.pi, 24, endpoint=False)
    rr = max(pad["sx"], pad["sy"]) / 2 + 0.15
    s = ctx.g.sample(r, pad["x"] + rr * np.cos(ang), pad["y"] + rr * np.sin(ang))
    return bool((s == pad["net_id"]).mean() > 0.5)


def rule_decap(ctx):
    b, cfg = ctx.b, ctx.cfg
    p = cfg["rules"]["decap"]
    V = b.vias
    metrics, finds = {}, []
    loads = defaultdict(list)
    for ref in p["load_refs"]:
        for q in b.pads_by_ref.get(ref, []):
            loads[q["net"]].append((q["x"], q["y"]))
    for c in _decaps(ctx):
        if c["dnp"]:
            continue
        row = {"value": c["value"], "net": c["power"], "side": c["side"]}
        for q in c["pads"]:
            vk = b.vias_of(q["net"])
            d = _pad_via_dist(q, V["x"][vk], V["y"][vk], V["d"][vk]) if len(vk) else np.zeros(0)
            nearest = float(d.min()) if len(d) else math.inf
            n_close = int((d <= p["via_mm"]).sum()) if len(d) else 0
            tag = "gnd" if q["net"] in cfg.reference_nets else "pwr"
            row[f"{tag}_via_mm"] = round(nearest, 3) if math.isfinite(nearest) else None
            row[f"{tag}_vias"] = n_close
            poured = _poured(ctx, q)
            row[f"{tag}_poured"] = poured
            if nearest > p["via_mm"]:
                finds.append(finding("decap", "info" if poured else "warning",
                                     f"{c['ref']} ({c['value']}) {tag} pad: nearest {q['net']} via {nearest:.2f} mm"
                                     + (f" (pad sits in a {q['net']} pour on its layer)" if poured else ""),
                                     c["power"], q["layers"][0] if q["layers"] else "", q["x"], q["y"],
                                     value=nearest if math.isfinite(nearest) else 99, limit=p["via_mm"], unit="mm",
                                     part=c["ref"]))
        if loads.get(c["power"]):
            lx, ly = np.array(loads[c["power"]]).T
            dl_ = float(np.min(np.hypot(lx - c["x"], ly - c["y"])))
            row["load_mm"] = round(dl_, 2)
            if dl_ > p["load_mm"]:
                finds.append(finding("decap", "info", f"{c['ref']} ({c['value']}) is {dl_:.1f} mm from its nearest load pin",
                                     c["power"], "", c["x"], c["y"], value=dl_, limit=p["load_mm"], unit="mm", part=c["ref"]))
        metrics[c["ref"]] = row
    return result("decap", "PC11", "Decoupling capacitors: pad-to-via distance and distance to the load",
                  p, metrics, finds)


def rule_thermal_vias(ctx):
    b, cfg = ctx.b, ctx.cfg
    p = cfg["rules"]["thermal_vias"]
    therm = cfg["thermal"]["parts"]
    V = b.vias
    metrics, finds = {}, []
    parts = p["parts"] or [r for r, pads in b.pads_by_ref.items() if len(pads) >= 4]
    for ref in parts:
        best = None
        for q in b.pads_by_ref.get(ref, []):
            if q["net"] not in cfg.reference_nets or not q["poly"] or q["kind"] != "smd":
                continue
            rings = [ring for poly in q["poly"] for ring in poly]
            area = abs(sum(geom.ring_area(r) for r in rings))
            if area < 1.0 and ref not in p["parts"]:
                continue
            vk = b.vias_of(q["net"])
            n_in = int(geom.points_in_rings(V["x"][vk], V["y"][vk], rings).sum()) if len(vk) else 0
            if best is None or area > best[1]:
                best = (q, area, n_in)
        if best is None:
            continue
        q, area, n_in = best
        pw = therm.get(ref, {}).get("power_w")
        metrics[ref] = {"pad": q["num"], "pad_area_mm2": round(area, 2), "vias_in_pad": n_in, "power_w": pw}
        if n_in <= 1:
            finds.append(finding("thermal_vias", "info",
                                 f"{ref} thermal pad {q['num']} ({area:.1f} mm2): {n_in} via(s)" + (f", {pw} W" if pw else ""),
                                 q["net"], q["layers"][0] if q["layers"] else "", q["x"], q["y"], value=n_in, unit="vias",
                                 part=ref))
    return result("thermal_vias", "PC13", "Thermal vias in exposed pads", p, metrics, finds)


def rule_power_overlap(ctx):
    b, cfg, g = ctx.b, ctx.cfg, ctx.g
    p = cfg["rules"]["power_overlap"]
    refset = set(cfg.reference_nets)
    metrics, finds = {}, []
    for a, c in zip(b.layers[:-1], b.layers[1:]):
        h = b.stack.gap(a, c)
        if h > p["max_dielectric_mm"]:
            continue
        za = {z["net"] for z in b.zones if z["layer"] == a and z["net"] not in refset}
        zc = {z["net"] for z in b.zones if z["layer"] == c and z["net"] not in refset}
        if not za or not zc:
            continue
        ra, rc = b.raster(a), b.raster(c)
        for na in za:
            ma = ra == b.net_id[na]
            for nc in zc:
                if na == nc:
                    continue
                ov = ma & (rc == b.net_id[nc])
                A = float(ov.sum()) * g.res ** 2
                if A >= p["min_area_mm2"]:
                    ii, jj = np.nonzero(ov)
                    x, y = g.xy(ii.mean(), jj.mean())
                    key = f"{na}@{a}/{nc}@{c}"
                    metrics[key] = {"area_mm2": round(A, 2), "dielectric_mm": round(h, 4)}
                    finds.append(finding("power_overlap", "warning",
                                         f"{na} on {a} overlaps {nc} on {c} by {A:.1f} mm2 across {h*1000:.0f} um",
                                         na, a, x, y, value=A, unit="mm2"))
    return result("power_overlap", "G1", "Power areas of different nets overlapping across a thin dielectric",
                  p, metrics, finds)


def rule_stubs(ctx):
    """Dangling track ends (all nets) and T-branches (nets in scope)."""
    b, cfg = ctx.b, ctx.cfg
    p = cfg["rules"]["stubs"]
    T, V = b.tracks, b.vias
    tol = 0.002
    key = lambda li, x, y: (li, round(x / tol), round(y / tol))
    scope = set(ctx.nets("stubs"))
    metrics, finds = {}, []
    vtree = cKDTree(np.c_[V["x"], V["y"]]) if len(V["x"]) else None
    for nid in np.unique(T["net"]):
        if nid == 0:
            continue
        net = b.net_names[nid]
        ks = np.flatnonzero(T["net"] == nid)
        adj = defaultdict(list)
        n_ = len(ks)
        lab = geom.cluster(np.r_[T["x1"][ks], T["x2"][ks]] + 1000.0 * np.r_[T["layer"][ks], T["layer"][ks]],
                           np.r_[T["y1"][ks], T["y2"][ks]], 0.005)
        rep = {}
        for i, k in enumerate(ks):
            li = int(T["layer"][k])
            for e, (xx, yy) in enumerate(((T["x1"][k], T["y1"][k]), (T["x2"][k], T["y2"][k]))):
                c = lab[i + e * n_]
                if c not in rep:
                    rep[c] = key(li, xx, yy)
                adj[rep[c]].append(k)
        pads = [q for q in b.pads if q["net_id"] == nid]

        def terminal(node, hw):
            li, xi, yi = node
            x, y = xi * tol, yi * tol
            if vtree is not None:
                for v in vtree.query_ball_point([x, y], r=1.0):
                    if V["net"][v] == nid and math.hypot(V["x"][v] - x, V["y"][v] - y) <= V["d"][v] / 2 + hw:
                        return True
            for q in pads:
                if li in q["layer_ids"] and abs(q["x"] - x) <= q["sx"] + q["sy"] and abs(q["y"] - y) <= q["sx"] + q["sy"]:
                    rings = [ring for poly in q["poly"] for ring in poly]
                    if rings and (geom.points_in_rings([x], [y], rings)[0]
                                  or np.isfinite(_pad_via_dist(q, np.array([x]), np.array([y]), np.array([0.0]))[0])
                                  and _pad_via_dist(q, np.array([x]), np.array([y]), np.array([0.0]))[0] <= hw):
                        return True
            # a track end lying on the body of another track (T without a vertex)
            for k in ks:
                if T["layer"][k] == li and geom.point_seg_dist(x, y, T["x1"][k], T["y1"][k], T["x2"][k], T["y2"][k]) < T["w"][k] / 2 \
                        and min(math.hypot(T["x1"][k] - x, T["y1"][k] - y), math.hypot(T["x2"][k] - x, T["y2"][k] - y)) > 0.01:
                    return True
            # inside a zone of the same net on that layer
            r = b.raster(li)
            return ctx.g.sample(r, [x], [y])[0] == nid and any(z["net"] == net and z["layer"] == b.layers[li] for z in b.zones)

        dangling = []
        for node, kk in adj.items():
            if len(kk) == 1 and not terminal(node, float(T["w"][kk[0]]) / 2):
                k = kk[0]
                L = float(T["len"][k])
                if L >= p["dangling_min_mm"]:
                    dangling.append((node, L))
        branches = []
        if net in scope:
            for node, kk in adj.items():
                if len(kk) >= 3:
                    lens = sorted(float(T["len"][k]) for k in kk)
                    branches.append((node, lens))
        if dangling or branches:
            metrics[net] = {"dangling": len(dangling), "dangling_mm": round(sum(L for _, L in dangling), 3),
                            "branch_points": len(branches)}
        for (li, xi, yi), L in dangling:
            finds.append(finding("stubs", "warning", f"dangling track end ({L:.2f} mm segment)", net, b.layers[li],
                                 xi * tol, yi * tol, value=L, unit="mm"))
        for (li, xi, yi), lens in branches:
            if True:
                finds.append(finding("stubs", "info" if lens[0] <= p["max_stub_mm"] else "warning",
                                     f"T-branch: shortest branch segment {lens[0]:.2f} mm", net, b.layers[li],
                                     xi * tol, yi * tol, value=lens[0], limit=p["max_stub_mm"], unit="mm"))
    return result("stubs", "G1", "Dangling track ends and branches in nets in scope", p, metrics, finds)


# ==============================================================================================

RULES = [
    ("ref_plane", rule_ref_plane),          # also yields split_crossing
    ("unreferenced", rule_unreferenced),
    ("return_via", rule_return_via),
    ("via_stub", rule_via_stub),
    ("edge_distance", rule_edge_distance),
    ("stitching", rule_stitching),
    ("pair_coupling", rule_pair_coupling),
    ("parallel", rule_parallel),
    ("pad_plane", rule_pad_plane),
    ("decap", rule_decap),
    ("thermal_vias", rule_thermal_vias),
    ("power_overlap", rule_power_overlap),
    ("stubs", rule_stubs),
]


def _rules2():
    from . import rules2
    return [("stitching_map", rules2.stitching_map), ("return_proximity", rules2.return_proximity),
            ("breakout", rules2.breakout)]


def run(ctx, only=None, log=print):
    import time
    out = []
    for name, fn in RULES + _rules2():
        if only and name not in only:
            continue
        if not ctx.cfg["rules"][name].get("enabled", True):
            continue
        t = time.time()
        r = fn(ctx)
        rs = r if isinstance(r, list) else [r]
        for x in rs:
            x["seconds"] = round(time.time() - t, 2)
            sev = defaultdict(int)
            for f in x["findings"]:
                sev[f["severity"] + ("_in_areas" if f.get("area") else "")] += 1
            x["summary"].update(dict(sev))
            log(f"  {x['rule']:15s} {x['spec']:5s} {len(x['metrics']):4d} items  "
                + "  ".join(f"{k}={v}" for k, v in sorted(sev.items())) + f"  ({x['seconds']} s)")
        out += rs
    return out
