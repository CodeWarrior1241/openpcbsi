"""MF7 / MF8: the board in a convection reflow oven.

MF7, board level. The 2.5-D finite-volume model of thermal.py (one node per copper layer per cell)
with heat capacity added: the laminate, the copper and each part as a lump on the cells under it.
Both faces exchange heat with the oven air through a film coefficient h (forced convection with the
radiation linearised into it; about 50 W/m2K in the literature). The air temperature follows the
oven's zones in time. Backward Euler, one factorisation. Per part: peak temperature, time above
liquidus, steepest ramp, and the lag time tau (below). With `calibrate`, the spike zones are shifted
until the median peak over the board is the assembler's stated reflow temperature, which is what
an operator does with a profiling board.

MF8, pad level (tombstoning). In a ramp every node rises at the ramp rate r and lags the air by
r tau_i, where tau = (K + H)^-1 C 1 is the first moment of the node's response (its Elmore delay):
K conduction, H the film conductances, C the heat capacities. Two joints of one chip therefore reach
liquidus tau_1 - tau_2 seconds apart, whatever the ramp rate. tau is solved in a window around the
part: the part's own copper layer at 50 um, the layers below at 200 um, the chip as two nodes (its
halves) that take heat from the air over their own surface, pass it to the pads through the unmelted
paste, and are joined through the body.

Not modelled: latent heat of the paste, flux evaporation, shadowing by tall neighbours, the
non-uniform film coefficient of real ovens, infrared heating. The results compare places on one
board; they do not predict a measured profile.
"""
import math

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from . import geom, packages, thermal
from .drc import finding, result
from .thermal import K_CU, K_FR4_IN, K_FR4_THRU, RHOC_CU, RHOC_FR4

K_SOLDER_PASTE = 2.0


def _lu(A):
    return spla.splu(A.tocsc(), permc_spec="MMD_AT_PLUS_A", diag_pivot_thresh=0.0, options={"SymmetricMode": True})


def air_profile(r, shift=0.0):
    """Time (s) and air temperature (C) of the oven: zones held for zone_s each, then the cooling zones.
    `shift` is added to the last `spike_zones` heating zones."""
    dt = r["dt_s"]
    zones = list(r["zones_c"])
    for k in range(1, int(r["spike_zones"]) + 1):
        zones[-k] += shift
    seg = [(z, r["zone_s"]) for z in zones] + [(z, r["cool_s"]) for z in r["cool_c"]]
    t_end = sum(s for _, s in seg)
    t = np.arange(1, int(round(t_end / dt)) + 1) * dt
    ta = np.empty(len(t))
    edge = 0.0
    for z, s in seg:
        ta[(t > edge - 1e-9) & (t <= edge + s + 1e-9)] = z
        edge += s
    return t, ta


def _part_cover(mesh, fp):
    """Share of each cell (of the part's outer layer) under the part's body."""
    g = mesh["g"]
    k, ny, nx = mesh["k"], mesh["ny"], mesh["nx"]
    m = np.zeros(g.shape, bool)
    rings = [r for poly in fp["courtyard"] for r in poly]
    if not rings:
        x0, y0, x1, y1 = fp["bbox"]
        rings = [[(x0, y0), (x1, y0), (x1, y1), (x0, y1)]]
    geom.raster_rings(g, rings, m)
    return m[:ny * k, :nx * k].reshape(ny, k, nx, k).mean(axis=(1, 3)) * mesh["board"]


def transient(A, idx, gconv, cap, t, ta, t0, surf):
    """Backward-Euler march. Returns the temperature history of the nodes `surf` (len(t), len(surf))."""
    dt = t[1] - t[0] if len(t) > 1 else t[0]
    lu = _lu(A + sp.diags(cap / dt))
    T = np.full(A.shape[0], float(t0))
    hist = np.empty((len(t), len(surf)), dtype=np.float32)
    for n in range(len(t)):
        T = lu.solve(cap / dt * T + gconv * ta[n])
        hist[n] = T[surf]
    return hist


def board(ctx, fabr, out_dir, log=print):
    """MF7. Returns the rule result (findings, per-part rows, the calibrated profile)."""
    b, cfg = ctx.b, ctx.cfg
    r = cfg["reflow"]
    h_fab = fabr.get("h_w_m2k", 50.0)
    h_top, h_bot = r["h_top"] or h_fab, r["h_bot"] or h_fab
    mesh = thermal.build(ctx, r["cell_mm"])
    L, ny, nx = mesh["frac"].shape
    bm = mesh["board"]
    A, idx, gconv = thermal.assemble(b.stack, b.vias, mesh, h_top, h_bot)
    C = thermal.capacity(b.stack, mesh)
    parts = []
    extra = np.zeros(A.shape[0])                 # film conductance of the parts' side faces
    for ref, fp in b.fps.items():
        if fp.get("dnp") or not any(p["kind"] == "smd" for p in b.pads_by_ref.get(ref, [])):
            continue
        cov = _part_cover(mesh, fp)
        if cov.sum() <= 0:
            continue
        body = packages.body_of(fp, cfg)
        l = 0 if fp["side"] == "F" else L - 1
        w = cov / cov.sum()
        C[l] += body["C"] * w
        side = 4.0 * math.sqrt(max(body["area"], 0.0)) * body["height"] * 1e-6       # m2, a square body of that area
        extra[idx[l][bm]] += (h_top if l == 0 else h_bot) * side * w[bm]
        parts.append({"ref": ref, "layer": l, "cov": cov, "body": body, "x": fp["x"], "y": fp["y"]})
    A = (A + sp.diags(extra)).tocsc()
    gconv = gconv + extra
    cap = np.zeros(A.shape[0])
    for l in range(L):
        cap[idx[l][bm]] = C[l][bm]
    surf = np.concatenate([idx[0][bm], idx[L - 1][bm]])
    nb = int(bm.sum())
    t, ta0 = air_profile(r)
    _, ta1 = air_profile(r, 1.0)
    base = transient(A, idx, gconv, cap, t, ta0, r["start_c"], surf)
    target = r["peak_target_c"] or fabr.get("reflow_peak_c")
    shift = 0.0
    if r["calibrate"] and target:
        unit = transient(A, idx, gconv, cap, t, ta1 - ta0, 0.0, surf)       # response to +1 C on the spike zones
        lo, hi = -150.0, 150.0
        for _ in range(40):
            shift = (lo + hi) / 2
            if np.median((base + shift * unit).max(axis=0)) < target:
                lo = shift
            else:
                hi = shift
        hist = base + shift * unit
        del unit
    else:
        hist = base
    del base
    tau = _lu(A).solve(cap)                                                    # lag time of every node, s
    dt = r["dt_s"]
    liq = r["liquidus_c"]
    peak = hist.max(axis=0)
    rows = []
    for p in parts:
        w = p["cov"][bm]
        w = w / w.sum()
        off = 0 if p["layer"] == 0 else nb
        T = hist[:, off:off + nb] @ w
        above = T >= liq
        ramp = np.diff(T) / dt
        win = max(1, int(round(r["ramp_window_s"] / dt)))                    # steepest rise, averaged over the window
        ramp = np.convolve(ramp, np.ones(win) / win, mode="valid") if len(ramp) >= win else ramp
        rows.append({"ref": p["ref"], "x": p["x"], "y": p["y"], "side": "F" if p["layer"] == 0 else "B",
                     "peak_c": round(float(T.max()), 1), "tal_s": round(float(above.sum() * dt), 1),
                     "t_liquidus_s": round(float(t[np.argmax(above)]), 1) if above.any() else None,
                     "ramp_c_s": round(float(ramp.max()), 2),
                     "lag_s": round(float(tau[idx[p["layer"]][bm]] @ w), 1),
                     "body": p["body"]["source"], "assumed": p["body"]["assumed"]})
    lim = {k: r[k] for k in ("peak_min_c", "peak_max_c", "tal_min_s", "tal_max_s", "spread_warn_c", "ramp_max_c_s", "ramp_window_s")}
    fs = []
    lname = lambda row: b.layers[0] if row["side"] == "F" else b.layers[-1]
    for row in rows:
        why = []
        if row["peak_c"] < lim["peak_min_c"]:
            why.append((f"peaks at {row['peak_c']:.0f} C, below {lim['peak_min_c']:g} C: risk of an unwetted joint",
                        row["peak_c"], lim["peak_min_c"], "C", "peak_low"))
        elif row["peak_c"] > lim["peak_max_c"]:
            why.append((f"peaks at {row['peak_c']:.0f} C, above {lim['peak_max_c']:g} C", row["peak_c"], lim["peak_max_c"],
                        "C", "peak_high"))
        if row["peak_c"] >= liq and not (lim["tal_min_s"] <= row["tal_s"] <= lim["tal_max_s"]):
            why.append((f"{row['tal_s']:.0f} s above liquidus ({liq:g} C), outside {lim['tal_min_s']:g}-{lim['tal_max_s']:g} s",
                        row["tal_s"], lim["tal_max_s"] if row["tal_s"] > lim["tal_max_s"] else lim["tal_min_s"], "s", "tal"))
        if row["ramp_c_s"] > lim["ramp_max_c_s"]:
            why.append((f"heats at {row['ramp_c_s']:.1f} C/s, above {lim['ramp_max_c_s']:g} C/s", row["ramp_c_s"],
                        lim["ramp_max_c_s"], "C/s", "ramp"))
        for msg, val, limit, unit, kind in why:
            fs.append(finding("reflow", "warning", f"{row['ref']} {msg} (lag {row['lag_s']:.0f} s; {row['body']})",
                              layer=lname(row), x=row["x"], y=row["y"], value=val, limit=limit, unit=unit,
                              part=row["ref"], kind=kind))
    spread = None
    if rows:
        cold, hot = min(rows, key=lambda x: x["peak_c"]), max(rows, key=lambda x: x["peak_c"])
        spread = round(hot["peak_c"] - cold["peak_c"], 1)
        if spread > lim["spread_warn_c"]:
            fs.append(finding("reflow", "warning", f"peak temperature differs by {spread:.0f} C across the parts: {cold['ref']} "
                              f"{cold['peak_c']:.0f} C, {hot['ref']} {hot['peak_c']:.0f} C (limit {lim['spread_warn_c']:g} C)",
                              layer=lname(cold), x=cold["x"], y=cold["y"], value=spread, limit=lim["spread_warn_c"], unit="C",
                              part=cold["ref"], kind="spread"))
    img = None
    try:
        from . import report
        pk = np.full((ny, nx), np.nan)
        pk[bm] = peak[:nb]
        img = report.heatmap(pk, mesh, b, f"{out_dir}/reflow_peak_top.png", "Peak temperature in reflow, top layer (C)")
    except Exception as e:
        log(f"  (no image: {e})")
    n_assumed = sum(1 for x in rows if x["assumed"])
    zones = list(r["zones_c"])
    for k in range(1, int(r["spike_zones"]) + 1):
        zones[-k] = round(zones[-k] + shift, 1)
    log(f"  {len(rows)} parts, h = {h_top:g}/{h_bot:g} W/m2K, spike zones shifted {shift:+.1f} C for a median peak of "
        f"{float(np.median(peak)):.1f} C" + (f" (target {target:g} C)" if target else ""))
    if rows:
        log(f"  peak {min(x['peak_c'] for x in rows):.1f} to {max(x['peak_c'] for x in rows):.1f} C over the parts; "
            f"lag {min(x['lag_s'] for x in rows):.0f} to {max(x['lag_s'] for x in rows):.0f} s; {n_assumed} parts with an assumed body")
    params = {"cell_mm": r["cell_mm"], "dt_s": dt, "h_top": h_top, "h_bot": h_bot, "h_source": fabr.src("h_w_m2k") if not r["h_top"] else "board config",
              "zones_c": zones, "zone_s": r["zone_s"], "cool_c": r["cool_c"], "cool_s": r["cool_s"], "shift_c": round(shift, 1),
              "target_c": target, "target_source": "board config" if r["peak_target_c"] else fabr.src("reflow_peak_c"),
              "liquidus_c": liq, "limits": lim, "limits_source": "proposed defaults ([reflow] in the board config)"}
    return dict(result("reflow", "MF7", "Reflow: peak temperature, time above liquidus and lag of every part", params,
                       {x["ref"]: {"peak_c": x["peak_c"], "tal_s": x["tal_s"], "lag_s": x["lag_s"]} for x in rows}, fs,
                       {"spread_c": spread, "median_peak_c": round(float(np.median(peak)), 1), "parts": len(rows),
                        "assumed_bodies": n_assumed}),
                parts=rows, image=img)


# ------------------------------------------------------------------------------------------------
# MF8: the window model

def _blockmean(a, k):
    ny, nx = a.shape[0] // k, a.shape[1] // k
    return a[:ny * k, :nx * k].reshape(ny, k, nx, k).mean(axis=(1, 3))


def window_matrix(stack, cu, inside, vias, lf, kf, m, res, h_fine, h_far):
    """Conduction + film matrix and heat capacities of a window.

    cu: list (one per copper layer) of boolean rasters of the window at `res` mm, sized a multiple of
    kf * m; inside: the board area; vias: dict x, y (pixel coordinates in the window), top, bot, drill.
    Layer `lf` is meshed at kf pixels per cell, every other layer at kf * m. Returns A (CSC), the
    capacity vector, the fine-layer node index (ny, nx; -1 outside the board) and the cell size (m)."""
    L = len(cu)
    kc = kf * m
    bf, bc = _blockmean(inside, kf) > 0.5, _blockmean(inside, kc) > 0.5
    nyf, nxf = bf.shape
    nyc, nxc = bc.shape
    half = thermal._half_dielectric(stack)
    cs_f, cs_c = kf * res * 1e-3, kc * res * 1e-3
    idx, n = [], 0
    for l in range(L):
        bm = bf if l == lf else bc
        ix = -np.ones(bm.shape, np.int64)
        ix[bm] = n + np.arange(int(bm.sum()))
        n += int(bm.sum())
        idx.append(ix)
    rows, cols, vals = [], [], []
    diag = np.zeros(n)
    cap = np.zeros(n)

    def link(a, c, G):
        rows.extend([a, c]); cols.extend([c, a]); vals.extend([-G, -G])
        np.add.at(diag, a, G)
        np.add.at(diag, c, G)

    for l in range(L):
        k = kf if l == lf else kc
        cs = cs_f if l == lf else cs_c
        frac = _blockmean(cu[l], k)
        t = stack.copper[l]["t"] * 1e-3
        Gs = K_CU * t * frac + K_FR4_IN * (half[l] * 1e-3 + t * (1 - frac))
        ix = idx[l]
        ny, nx = ix.shape
        for di, dj in ((0, 1), (1, 0)):
            a, c = ix[: ny - di, : nx - dj], ix[di:, dj:]
            ok = (a >= 0) & (c >= 0)
            link(a[ok], c[ok], 2 / (1 / Gs[: ny - di, : nx - dj][ok] + 1 / Gs[di:, dj:][ok]))
        ok = ix >= 0
        cap[ix[ok]] = cs * cs * (RHOC_CU * t * frac[ok] + RHOC_FR4 * (half[l] * 1e-3 + t * (1 - frac[ok])))
    have_v = vias is not None and len(vias["x"]) > 0
    for l in range(L - 1):
        dz = (stack.copper[l + 1]["z_top"] - stack.copper[l]["z_bot"]) * 1e-3
        fine = l == lf or l + 1 == lf
        k = kf if fine else kc
        cs = cs_f if fine else cs_c
        shape = (nyf, nxf) if fine else (nyc, nxc)
        G = np.full(shape, K_FR4_THRU * cs * cs / dz)
        if have_v:
            span = (vias["top"] <= l) & (vias["bot"] >= l + 1)
            i, j = np.floor(vias["y"] / k).astype(int), np.floor(vias["x"] / k).astype(int)
            ok = span & (i >= 0) & (i < shape[0]) & (j >= 0) & (j < shape[1])
            ann = math.pi * ((vias["drill"] / 2 + 0.02) ** 2 - (vias["drill"] / 2) ** 2) * 1e-6
            np.add.at(G, (i[ok], j[ok]), K_CU * ann[ok] / dz)
        if fine:
            fi = idx[lf]
            ci = np.repeat(np.repeat(idx[l + 1 if l == lf else l], m, axis=0), m, axis=1)
            ok = (fi >= 0) & (ci >= 0)
            link(fi[ok], ci[ok], G[ok])
        else:
            a, c = idx[l], idx[l + 1]
            ok = (a >= 0) & (c >= 0)
            link(a[ok], c[ok], G[ok])
    far = L - 1 if lf == 0 else 0
    gconv = np.zeros(n)
    gconv[idx[lf][idx[lf] >= 0]] += h_fine * cs_f * cs_f
    if far != lf:
        gconv[idx[far][idx[far] >= 0]] += h_far * cs_c * cs_c
    rows = np.concatenate([np.atleast_1d(x) for x in rows])
    cols = np.concatenate([np.atleast_1d(x) for x in cols])
    vals = np.concatenate([np.atleast_1d(x) for x in vals])
    A = sp.coo_matrix((vals, (rows, cols)), shape=(n, n)).tocsr() + sp.diags(diag + gconv)
    return A, cap, gconv, idx[lf], cs_f


def chip_lag(ctx, fp, chip, pads, h_part, h_far, p):
    """Lag times (s) of the two joints of a chip: pad, termination and their mean, per pad."""
    b = ctx.b
    g = b.grid()
    L = len(b.layers)
    lf = 0 if fp["side"] == "F" else L - 1
    kf = max(1, int(round(p["cell_mm"] / g.res)))
    m = max(1, int(p["coarse"]))
    kc = kf * m
    pts = np.vstack([np.asarray(r, dtype=float) for pad in pads for poly in pad["poly"] for r in poly])
    mar = p["margin_mm"]
    j0 = int(math.floor((pts[:, 0].min() - mar - g.x0) / g.res / kc)) * kc
    j1 = int(math.ceil((pts[:, 0].max() + mar - g.x0) / g.res / kc)) * kc
    i0 = int(math.floor((pts[:, 1].min() - mar - g.y0) / g.res / kc)) * kc
    i1 = int(math.ceil((pts[:, 1].max() + mar - g.y0) / g.res / kc)) * kc
    j0, i0 = max(j0, 0), max(i0, 0)
    j1, i1 = min(j1, (g.nx // kc) * kc), min(i1, (g.ny // kc) * kc)
    sl = (slice(i0, i1), slice(j0, j1))
    cu = [b.raster(l)[sl] != 0 for l in range(L)]
    V = b.vias
    vx, vy = (V["x"] - g.x0) / g.res - j0, (V["y"] - g.y0) / g.res - i0
    near = (vx >= 0) & (vx < j1 - j0) & (vy >= 0) & (vy < i1 - i0)
    vias = {"x": vx[near], "y": vy[near], "top": V["top"][near], "bot": V["bot"][near], "drill": V["drill"][near]}
    A, cap, gconv, fidx, cs = window_matrix(b.stack, cu, thermal._inside(ctx)[sl], vias, lf, kf, m, g.res, h_part, h_far)
    n = A.shape[0]
    sub = geom.Grid.__new__(geom.Grid)
    sub.res, sub.x0, sub.y0, sub.nx, sub.ny = g.res, g.x0 + j0 * g.res, g.y0 + i0 * g.res, j1 - j0, i1 - i0
    covs = []
    for pad in pads:
        mk = np.zeros(sub.shape, bool)
        geom.raster_rings(sub, [ring for poly in pad["poly"] for ring in poly], mk)
        covs.append(_blockmean(mk, kf) * (fidx >= 0))
    # two nodes for the halves of the chip
    surf = (chip["L"] * chip["W"] + 2 * chip["H"] * (chip["L"] + chip["W"])) * 1e-6
    gap = max(math.hypot(pads[0]["x"] - pads[1]["x"], pads[0]["y"] - pads[1]["y"]), 0.1) * 1e-3
    g_body = chip["k"] * chip["W"] * chip["H"] * 1e-6 / gap
    g_paste = p["paste_k"] / (p["paste_um"] * 1e-6) * cs * cs            # per fully covered cell
    rows, cols, vals = [n, n + 1], [n + 1, n], [-g_body, -g_body]
    dg = np.zeros(n + 2)
    dg[n] = dg[n + 1] = g_body + h_part * surf / 2
    for k, cov in enumerate(covs):
        cells = fidx[cov > 0]
        gk = g_paste * cov[cov > 0]
        rows += [np.full(len(cells), n + k), cells]
        cols += [cells, np.full(len(cells), n + k)]
        vals += [-gk, -gk]
        np.add.at(dg, cells, gk)
        dg[n + k] += gk.sum()
    rows = np.concatenate([np.atleast_1d(x) for x in rows])
    cols = np.concatenate([np.atleast_1d(x) for x in cols])
    vals = np.concatenate([np.atleast_1d(x) for x in vals])
    A2 = sp.bmat([[A, None], [None, sp.csr_matrix((2, 2))]]).tocsr() + \
        sp.coo_matrix((vals, (rows, cols)), shape=(n + 2, n + 2)).tocsr() + sp.diags(dg)
    c2 = np.concatenate([cap, [chip["C"] / 2, chip["C"] / 2]])
    tau = _lu(A2).solve(c2)
    out = []
    for k, cov in enumerate(covs):
        w = cov[cov > 0]
        t_pad = float(tau[fidx[cov > 0]] @ w / w.sum()) if w.sum() > 0 else float("nan")
        out.append({"pad": pads[k]["num"], "net": pads[k]["net"], "tau_pad_s": t_pad, "tau_term_s": float(tau[n + k]),
                    "tau_joint_s": (t_pad + float(tau[n + k])) / 2})
    return out, n + 2


def chips_in_scope(ctx):
    """Two-terminal chips to judge: [(footprint, chip, [pad a, pad b])]."""
    b, cfg = ctx.b, ctx.cfg
    top = packages.size_rank(cfg["dfm"]["pad_balance"]["max_package"])
    out = []
    for ref, fp in b.fps.items():
        pads = b.pads_by_ref.get(ref, [])
        if fp.get("dnp") or len(pads) != 2 or any(p["kind"] != "smd" or not p["poly"] or len(p["layer_ids"]) != 1 for p in pads):
            continue
        if pads[0]["layer_ids"] != pads[1]["layer_ids"] or pads[0]["net"] == pads[1]["net"]:
            continue
        chip = packages.chip_of(fp)
        if chip is None or packages.size_rank(chip["code"]) > top:
            continue
        out.append((fp, chip, pads))
    return out


def tombstone(ctx, fabr, order=None, log=print):
    """MF8 over the chips in scope (`order`: references to solve first when max_parts limits the run)."""
    b, cfg = ctx.b, ctx.cfg
    p = cfg["dfm"]["tombstone"]
    r = cfg["reflow"]
    h_fab = fabr.get("h_w_m2k", 50.0)
    h_top, h_bot = r["h_top"] or h_fab, r["h_bot"] or h_fab
    small = packages.size_rank(cfg["dfm"]["pad_balance"]["small_package"])
    chips = chips_in_scope(ctx)
    if order:
        rank = {ref: k for k, ref in enumerate(order)}
        chips.sort(key=lambda c: rank.get(c[0]["ref"], len(rank)))
    skipped = 0
    if p["max_parts"] and len(chips) > p["max_parts"]:
        skipped = len(chips) - p["max_parts"]
        chips = chips[: p["max_parts"]]
    rows, fs = [], []
    for fp, chip, pads in chips:
        front = fp["side"] == "F"
        try:
            j, n = chip_lag(ctx, fp, chip, pads, h_top if front else h_bot, h_bot if front else h_top, p)
        except Exception as e:                       # a window that cannot be solved is reported, not hidden
            fs.append(finding("tombstone", "warning", f"{fp['ref']}: pad-level model failed: {e}", part=fp["ref"],
                              layer=b.layers[pads[0]["layer_ids"][0]], x=fp["x"], y=fp["y"], kind="problem"))
            continue
        d = j[0]["tau_joint_s"] - j[1]["tau_joint_s"]
        first, last = (j[1], j[0]) if d > 0 else (j[0], j[1])
        row = {"ref": fp["ref"], "package": chip["code"], "kind": chip["kind"], "x": fp["x"], "y": fp["y"],
               "dt_s": round(abs(d), 3), "first_pad": first["pad"], "first_net": first["net"], "last_pad": last["pad"],
               "last_net": last["net"], "tau_s": [round(j[0]["tau_joint_s"], 2), round(j[1]["tau_joint_s"], 2)], "nodes": n}
        rows.append(row)
        if abs(d) >= p["dt_warn_s"]:
            sev = "warning" if packages.size_rank(chip["code"]) <= small else "info"
            fs.append(finding("tombstone", sev, f"{fp['ref']} ({chip['code']}): pad {first['pad']} ({first['net']}) reaches liquidus "
                              f"{abs(d):.2f} s before pad {last['pad']} ({last['net']}) (limit {p['dt_warn_s']:g} s, proposed)",
                              net=last["net"], layer=b.layers[pads[0]["layer_ids"][0]], x=fp["x"], y=fp["y"],
                              value=round(abs(d), 3), limit=p["dt_warn_s"], unit="s", part=fp["ref"], kind="dt"))
    rows.sort(key=lambda x: -x["dt_s"])
    if rows:
        ds = np.array([x["dt_s"] for x in rows])
        log(f"  {len(rows)} chips solved" + (f" ({skipped} left out by max_parts)" if skipped else "")
            + f": difference between the two joints median {np.median(ds):.2f} s, 90 % below {np.percentile(ds, 90):.2f} s, "
              f"largest {ds.max():.2f} s ({rows[0]['ref']})")
    params = dict(p, h_top=h_top, h_bot=h_bot, limit_source="proposed default ([dfm.tombstone] dt_warn_s)")
    return dict(result("tombstone", "MF8", "Tombstoning: time between the two joints of a chip reaching liquidus", params,
                       {x["ref"]: {"dt_s": x["dt_s"]} for x in rows}, fs, {"chips": len(rows), "not_solved": skipped}), parts=rows)


# ------------------------------------------------------------------------------------------------
def validation_cases():
    """Closed forms on a uniform two-layer board: step response, ramp lag, lag time, window model."""
    from .model import Stackup
    from .validate import _case
    st = Stackup([], ["F.Cu", "B.Cu"], 1.6)
    ny = nx = 12
    mesh = {"k": 1, "cs": 1.0, "ny": ny, "nx": nx, "board": np.ones((ny, nx), bool),
            "frac": np.ones((2, ny, nx)), "g": geom.window_grid(0, 0, nx, ny, 1.0)}
    h = 50.0
    A, idx, gconv = thermal.assemble(st, None, mesh, h, h)
    C = thermal.capacity(st, mesh)
    cap = np.zeros(A.shape[0])
    for l in range(2):
        cap[idx[l].ravel()] = C[l].ravel()
    tau_ref = float(C[:, 0, 0].sum() / (2 * h * 1e-6))
    out = []
    dt = 0.5
    t = np.arange(1, int(6 * tau_ref / dt)) * dt
    hist = transient(A, idx, gconv, cap, t, np.full(len(t), 200.0), 25.0, np.array([0]))[:, 0]
    k = int(round(tau_ref / dt)) - 1
    out.append(_case("step of the air temperature: rise after one time constant", "reflow", float(hist[k] - 25.0),
                     175.0 * (1 - math.exp(-t[k] / tau_ref)), 1.0, "C", "lumped body, tau = C / (h_top + h_bot)"))
    rate = 1.5
    hist = transient(A, idx, gconv, cap, t, 25.0 + rate * t, 25.0, np.array([0]))[:, 0]
    out.append(_case("ramp of the air temperature: steady lag", "reflow", float(25.0 + rate * t[-1] - hist[-1]),
                     rate * tau_ref * (1 - math.exp(-t[-1] / tau_ref)), 1.0, "C", "lag = ramp rate x tau"))
    out.append(_case("lag time from the steady solve", "reflow", float(_lu(A).solve(cap)[0]), tau_ref, 0.5, "s",
                     "(K + H)^-1 C 1 = C / (h_top + h_bot) for a uniform board"))
    res, kf, m = 0.025, 2, 4
    npx = 160
    cu = [np.ones((npx, npx), bool), np.ones((npx, npx), bool)]
    A, cap, gconv, fidx, cs = window_matrix(st, cu, np.ones((npx, npx), bool), None, 0, kf, m, res, h, h)
    tw = _lu(A).solve(cap)
    out.append(_case("window model (50 um over 200 um cells): lag time", "reflow", float(tw[fidx[40, 40]]), tau_ref, 0.5, "s",
                     "same uniform board: the two mesh sizes must agree with the lumped value"))
    return out
