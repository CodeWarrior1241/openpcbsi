"""MF7 / MF8: thermal rankings of the board in a convection reflow oven.

These are models, not measurements, and nothing here has been compared with a profiled board or
with an assembler's defect data. They order the parts of one board; they do not judge them. A limit
applies only when the user supplies one with its source (`[dfm.judge]`).

MF7, board level. The 2.5-D finite-volume model of thermal.py (one node per copper layer per cell)
with heat capacity: the laminate, the copper and each part as a lump on the cells under it. Both
faces, and the side faces of each part, exchange heat with the oven air through a film coefficient.
The air temperature follows the oven's zones, and each point of the board enters a zone when the
belt carries it there. With `calibrate`, the spike zones are shifted until the median peak over the
parts is the assembler's stated reflow temperature; the absolute peak is therefore an input, and
only the differences between parts are a result. The run is repeated with the film coefficient and
the assumed part bodies varied, and each part's difference from the median is given as a range.

MF8, chip level. A window around each two-terminal chip: its own copper layer at 50 um, the other
layers at 200 um, the chip as two nodes heated over their own surface and joined to the pads
through the unmelted paste and to each other through the body. The window's edge takes its values
from the board model, so the planes do not end at the window. One number per chip, the time between
its two joints reaching liquidus:
  lag difference   tau = (K + H)^-1 C 1 is the first moment of a node's response; in a slow ramp
                   two joints reach liquidus tau_1 - tau_2 apart. Solved for every chip.
  belt term        the belt carries one pad into each zone before the other: their distance along
                   the belt over its speed. Added for the stated direction ([reflow] travel); with
                   none stated the value excludes it and its range covers the four directions.
The range adds the sensitivity to the assumed paste, film coefficient and window size. For the
layouts at the top of the list, the same window marched through the oven profile is a consistency
check of the value (two computations of one quantity); a disagreement beyond the sensitivity is
reported as a problem.
Timing is one input to tombstoning. The wetting forces (pad geometry, paste volume and offset,
placement, alloy) are not modelled; the land-pattern geometry is listed beside the timing.

Not modelled: latent heat, flux, shadowing by tall parts, the non-uniform film coefficient of real
ovens, infrared heating, the temperature inside a part (a part is as warm as the board under it),
panels, and the second pass of a double-sided board.
"""
import math
import os

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy import ndimage

from . import geom, packages, par, thermal
from .drc import finding
from .thermal import K_CU, K_FR4_IN, K_FR4_THRU, RHOC_CU, RHOC_FR4

# Variants of the model's parameters (not of its form): the spread of a result over them is its
# sensitivity to those parameters, not an error bar. A part without body data gets an estimated body
# (packages.body_of: its fabrication outline and a typical height for its class), marked "estimated",
# and the estimated bodies are varied x0.5 and x2.
VARIANTS = (("film coefficient x 0.6", 0.6, 1.0), ("film coefficient x 1.6", 1.6, 1.0),
            ("model and estimated bodies x 0.5", 1.0, 0.5), ("model and estimated bodies x 2", 1.0, 2.0))


def _lu(A):
    return spla.splu(A.tocsc(), permc_spec="MMD_AT_PLUS_A", diag_pivot_thresh=0.0, options={"SymmetricMode": True})


def profile_fn(r, shift=0.0, unit=False):
    """Air temperature against the time since a point entered the oven. `shift` is added to the last
    `spike_zones` heating zones; `unit` gives the change for a shift of +1 C instead."""
    zones = list(r["zones_c"])
    ks = range(len(zones) - int(r["spike_zones"]), len(zones))
    for k in ks:
        zones[k] += shift
    temps = zones + list(r["cool_c"])
    edges = np.cumsum([r["zone_s"]] * len(zones) + [r["cool_s"]] * len(r["cool_c"]))
    if unit:
        temps = [1.0 if k in ks else 0.0 for k in range(len(zones))] + [0.0] * len(r["cool_c"])
    vals = np.array([0.0 if unit else r["start_c"]] + temps)

    def fn(tau):
        tau = np.asarray(tau, dtype=float)
        k = np.searchsorted(edges, tau, side="left") + 1
        k = np.where(tau <= 0, 0, np.minimum(k, len(vals) - 1))
        return vals[k]
    return fn, float(edges[-1])


def times(r):
    _, end = profile_fn(r)
    return np.arange(1, int(round(end / r["dt_s"])) + 1) * r["dt_s"]


def march(A, cap, gconv, off, pfn, t, t0, project=None, keep=False):
    """Backward Euler through the oven. `off` is each node's delay behind the leading edge (s).
    Returns (projected history, full history or None): `project` (rows x nodes) picks or averages."""
    dt = t[0]
    lu = _lu(A + sp.diags(cap / dt))
    T = np.full(A.shape[0], float(t0))
    act = np.flatnonzero(gconv > 0)
    hp = np.empty((len(t), project.shape[0]), dtype=np.float32) if project is not None else None
    full = np.empty((len(t), A.shape[0]), dtype=np.float32) if keep else None
    for n in range(len(t)):
        rhs = cap / dt * T
        rhs[act] += gconv[act] * pfn(t[n] - off[act])
        T = lu.solve(rhs)
        if hp is not None:
            hp[n] = project @ T
        if keep:
            full[n] = T
    return hp, full


def _part_cover(ctx, mesh, fp):
    """Share of each cell (of the part's outer layer) under the part's body."""
    cache = ctx.__dict__.setdefault("_reflow_cov", {})
    if fp["ref"] not in cache:
        g = mesh["g"]
        k, ny, nx = mesh["k"], mesh["ny"], mesh["nx"]
        m = np.zeros(g.shape, bool)
        rings = [r for poly in fp["courtyard"] for r in poly]
        if not rings:
            x0, y0, x1, y1 = fp["bbox"]
            rings = [[(x0, y0), (x1, y0), (x1, y1), (x0, y1)]]
        geom.raster_rings(g, rings, m)
        cache[fp["ref"]] = m[:ny * k, :nx * k].reshape(ny, k, nx, k).mean(axis=(1, 3)) * mesh["board"]
    return cache[fp["ref"]]


def build_model(ctx, fabr, h_scale=1.0, body_scale=1.0):
    """The board with its parts: matrix, capacities, film conductances, belt delays, part weights."""
    b, cfg = ctx.b, ctx.cfg
    r = cfg["reflow"]
    h_fab = fabr.get("h_w_m2k", 50.0)
    h_top, h_bot = (r["h_top"] or h_fab) * h_scale, (r["h_bot"] or h_fab) * h_scale
    mesh = ctx.__dict__.get("_reflow_mesh")
    if mesh is None:
        mesh = ctx._reflow_mesh = thermal.build(ctx, r["cell_mm"])
    L, ny, nx = mesh["frac"].shape
    bm = mesh["board"]
    A, idx, gconv = thermal.assemble(b.stack, b.vias, mesh, h_top, h_bot)
    C = thermal.capacity(b.stack, mesh)
    N = A.shape[0]
    extra = np.zeros(N)
    parts, rows, cols, vals = [], [], [], []
    left_out = []
    models = ctx.__dict__.get("_reflow_models")
    if models is None:
        from . import bodies
        from .run import CACHE
        try:
            models = ctx._reflow_models = bodies.bodies(b.path, b.cache_dir or CACHE)
        except Exception:
            models = ctx._reflow_models = {}
    for ref, fp in b.fps.items():
        if fp.get("dnp") or not any(p["kind"] == "smd" for p in b.pads_by_ref.get(ref, [])):
            continue
        why = packages.not_a_part(fp, b.pads_by_ref.get(ref, []))
        if why:
            left_out.append({"ref": ref, "why": why})
            continue
        cov = _part_cover(ctx, mesh, fp)
        if cov.sum() <= 0:
            continue
        body = packages.body_of(fp, cfg, models.get(ref))
        if body is None:
            left_out.append({"ref": ref, "why": "no body data: no 3D model (or its file is missing), no Mass or Height field, and no "
                                                "package class to estimate one; give its mass or height in [reflow.parts." + ref + "]"})
            continue
        # the body variants vary every body not measured: estimates by class and 3D-model volumes (whose
        # material is not known)
        scale = body_scale if (body["assumed"] or body.get("model")) else 1.0
        l = 0 if fp["side"] == "F" else L - 1
        w = cov / cov.sum()
        C[l] += body["C"] * scale * w
        side = 4.0 * math.sqrt(max(body["area"], 0.0)) * body["height"] * scale * 1e-6      # m2, a square body of that area
        extra[idx[l][bm]] += (h_top if l == 0 else h_bot) * side * w[bm]
        nodes = idx[l][bm][w[bm] > 0]
        rows += [len(parts)] * len(nodes)
        cols += nodes.tolist()
        vals += w[bm][w[bm] > 0].tolist()
        parts.append({"ref": ref, "layer": l, "body": body, "x": fp["x"], "y": fp["y"], "nodes": nodes})
    A = (A + sp.diags(extra)).tocsc()
    gconv = gconv + extra
    cap = np.zeros(N)
    for l in range(L):
        cap[idx[l][bm]] = C[l][bm]
    # delay of each node behind the leading edge of the board
    v = r["zone_mm"] / r["zone_s"]
    jj, ii = np.meshgrid(np.arange(nx), np.arange(ny))
    cx, cy = mesh["g"].x0 + (jj + 0.5) * mesh["cs"], mesh["g"].y0 + (ii + 0.5) * mesh["cs"]
    trav = r["travel"] or "x"
    axis = cx if "x" in trav else cy
    lead = axis[bm].min() if trav.startswith("-") else axis[bm].max()
    delay = np.abs(lead - axis) / v
    off = np.zeros(N)
    for l in range(L):
        off[idx[l][bm]] = delay[bm]
    W = sp.csr_matrix((vals, (rows, cols)), shape=(len(parts), N)) if parts else sp.csr_matrix((0, N))
    return {"mesh": mesh, "A": A, "idx": idx, "gconv": gconv, "cap": cap, "off": off, "parts": parts, "W": W, "left_out": left_out,
            "h_top": h_top, "h_bot": h_bot, "belt_mm_s": v, "delay_s": float(delay[bm].max()), "cx": cx, "cy": cy}


def run_profile(M, r, target, keep=False):
    """March the model through the oven and, with a target, shift the spike zones until the median
    peak over the parts meets it. Returns (part history, full history or None, shift)."""
    t = times(r)
    pfn, _ = profile_fn(r)
    hp, full = march(M["A"], M["cap"], M["gconv"], M["off"], pfn, t, r["start_c"], M["W"], keep)
    shift = 0.0
    if r["calibrate"] and target and hp.shape[1]:
        ufn, _ = profile_fn(r, unit=True)
        up, ufull = march(M["A"], M["cap"], M["gconv"], M["off"], ufn, t, 0.0, M["W"], keep)
        lo, hi = -150.0, 150.0
        for _ in range(40):
            shift = (lo + hi) / 2
            if np.median((hp + shift * up).max(axis=0)) < target:
                lo = shift
            else:
                hi = shift
        hp = hp + shift * up
        if keep:
            full += shift * ufull
    return hp, full, shift


def board(ctx, fabr, out_dir, log=print):
    """MF7. A ranking of the parts by how far their peak lies from the board's median."""
    b, cfg = ctx.b, ctx.cfg
    r = cfg["reflow"]
    target = r["peak_target_c"] or fabr.get("reflow_peak_c")
    M = build_model(ctx, fabr)
    mesh, bm = M["mesh"], M["mesh"]["board"]
    L, ny, nx = mesh["frac"].shape
    t = times(r)
    keep = M["A"].shape[0] * len(t) * 4 < 400e6
    hp, full, shift = run_profile(M, r, target, keep)
    tau = _lu(M["A"]).solve(M["cap"])
    liq, dt = r["liquidus_c"], r["dt_s"]
    peaks = hp.max(axis=0) if hp.shape[1] else np.zeros(0)
    med = float(np.median(peaks)) if len(peaks) else None
    dev = {"base": peaks - med} if len(peaks) else {}
    lags = {"base": tau}
    if r["sensitivity"] and len(peaks):
        def variant(v):
            name, hs, bs = v
            Mv = build_model(ctx, fabr, hs, bs)
            hv, _, _ = run_profile(Mv, r, target)
            pv = hv.max(axis=0)
            return name, pv - np.median(pv), _lu(Mv["A"]).solve(Mv["cap"]) if bs == 1.0 else None
        for name, d, tv in par.pmap(variant, VARIANTS, par.workers(cfg)):
            dev[name] = d
            if tv is not None:
                lags[name] = tv
    node_peak = full.max(axis=0) if full is not None else None
    rows = []
    n_cold = max(1, int(round(0.1 * len(peaks)))) if len(peaks) else 0
    for k, p in enumerate(M["parts"]):
        T = hp[:, k]
        above = T >= liq
        ds = [float(d[k]) for d in dev.values()]
        row = {"ref": p["ref"], "x": p["x"], "y": p["y"], "side": "F" if p["layer"] == 0 else "B",
               "peak_c": round(float(T.max()), 1), "d_peak_c": round(ds[0], 1), "d_range_c": [round(min(ds), 1), round(max(ds), 1)],
               "tal_s": round(float(above.sum() * dt), 0), "lag_s": round(float(M["W"][k] @ tau), 1),
               "body": p["body"]["source"], "assumed": p["body"]["assumed"], "estimated": p["body"]["assumed"],
               "from_model": bool(p["body"].get("model"))}
        if node_peak is not None and len(p["nodes"]) >= 4:
            row["within_c"] = round(float(node_peak[p["nodes"]].max() - node_peak[p["nodes"]].min()), 1)
        rows.append(row)
    rows.sort(key=lambda x: x["d_peak_c"])
    fs = []
    j = cfg["dfm"]["judge"]
    if j.get("source") and j.get("reflow_cold_c") is not None:
        for row in rows:
            if -row["d_peak_c"] > j["reflow_cold_c"]:
                fs.append(finding("reflow", "warning", f"{row['ref']} peaks {-row['d_peak_c']:.1f} C below the board median "
                                  f"(sensitivity {-row['d_range_c'][1]:.1f} to {-row['d_range_c'][0]:.1f} C over the parameter variants), more than "
                                  f"{j['reflow_cold_c']:g} C (board config: {j['source']})", layer=b.layers[0 if row["side"] == "F" else -1],
                                  x=row["x"], y=row["y"], value=-row["d_peak_c"], limit=j["reflow_cold_c"], unit="C", part=row["ref"],
                                  kind="cold", basis=f"board config: {j['source']}"))
    img = None
    if node_peak is not None:
        try:
            from . import report
            pk = np.full((ny, nx), np.nan)
            pk[bm] = node_peak[M["idx"][0][bm]]
            img = report.heatmap(pk, mesh, b, f"{out_dir}/reflow_peak_top.png", "Modelled peak temperature in reflow, top layer (C)")
        except Exception as e:
            log(f"  (no image: {e})")
    n_assumed = sum(1 for x in rows if x["assumed"])
    zones = list(r["zones_c"])
    for k in range(len(zones) - int(r["spike_zones"]), len(zones)):
        zones[k] = round(zones[k] + shift, 1)
    if rows:
        log(f"  {len(rows)} parts; spike zones shifted {shift:+.1f} C so that the median peak is {med:.1f} C"
            + (f" (target {target:g} C)" if target else "") + f"; belt {M['belt_mm_s']:.1f} mm/s, {M['delay_s']:.0f} s from leading to trailing edge")
        log(f"  coldest {rows[0]['ref']} {rows[0]['d_peak_c']:+.1f} C from the median (sensitivity {rows[0]['d_range_c'][0]:+.1f} to "
            f"{rows[0]['d_range_c'][1]:+.1f}), hottest {rows[-1]['ref']} {rows[-1]['d_peak_c']:+.1f} C; {n_assumed} parts with estimated bodies "
            f"(give a mass or height in [reflow.parts.REF] or a Mass / Height field); {len(M.get('left_out', []))} footprints left out "
            f"(not a soldered part); {len(dev)} parameter variants")
    params = {"cell_mm": r["cell_mm"], "dt_s": dt, "h_top": M["h_top"], "h_bot": M["h_bot"],
              "h_source": fabr.src("h_w_m2k") if not r["h_top"] else "board config", "zones_c": zones, "zone_s": r["zone_s"],
              "cool_c": r["cool_c"], "cool_s": r["cool_s"], "shift_c": round(shift, 1), "target_c": target,
              "target_source": "board config" if r["peak_target_c"] else fabr.src("reflow_peak_c"), "liquidus_c": liq,
              "belt_mm_s": round(M["belt_mm_s"], 2), "travel": r["travel"] or "x (assumed: [reflow] travel not set)", "variants": list(dev)}
    return {"rule": "reflow", "spec": "MF7", "title": "Reflow: how far each part's peak lies from the board median (model)",
            "params": params, "findings": fs, "rows": rows, "image": img,
            "summary": {"parts": len(rows), "assumed_bodies": n_assumed, "median_peak_c": None if med is None else round(med, 1),
                        "spread_c": round(rows[-1]["d_peak_c"] - rows[0]["d_peak_c"], 1) if rows else None,
                        "spread_range_c": [round(float(min(d.max() - d.min() for d in dev.values())), 1),
                                           round(float(max(d.max() - d.min() for d in dev.values())), 1)] if dev else None,
                        "left_out": M.get("left_out", []),
                        "note": "the range is the spread over variants of the model's parameters: its sensitivity to them, not an error bar"},
            "_internal": {"M": M, "tau": lags, "hist": full, "t": t, "r": r, "shift": shift}}


# ------------------------------------------------------------------------------------------------
# MF8: the window model

def _blockmean(a, k):
    ny, nx = a.shape[0] // k, a.shape[1] // k
    return a[:ny * k, :nx * k].reshape(ny, k, nx, k).mean(axis=(1, 3))


def window_matrix(stack, cu, inside, vias, lf, kf, m, res, h_fine, h_far):
    """Conduction + film matrix and heat capacities of a window.

    cu: list (one per copper layer) of boolean rasters of the window at `res` mm, sized a multiple of
    kf * m; inside: the board area; vias: dict x, y (pixel coordinates in the window), top, bot, drill.
    Layer `lf` is meshed at kf pixels per cell, every other layer at kf * m. Returns A (CSR, film on the
    diagonal), the capacities, the film conductances, the node index of every layer (-1 outside the
    board) and the fine cell size (m)."""
    L = len(cu)
    kc = kf * m
    bf, bc = _blockmean(inside, kf) > 0.5, _blockmean(inside, kc) > 0.5
    nyf, nxf = bf.shape
    nyc, nxc = bc.shape
    half = thermal._half_dielectric(stack)
    cs_f, cs_c = kf * res * 1e-3, kc * res * 1e-3
    idx, n = [], 0
    for l in range(L):
        bmk = bf if l == lf else bc
        ix = -np.ones(bmk.shape, np.int64)
        ix[bmk] = n + np.arange(int(bmk.sum()))
        n += int(bmk.sum())
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
    return A, cap, gconv, idx, cs_f


def chip_window(ctx, fp, chip, pads, h_part, h_far, p, paste_scale=1.0, margin=None):
    """The window model of one chip: matrix with the two chip nodes appended, capacities, film
    conductances, the nodes on the window's edge with their layer and position, and the pad weights."""
    b = ctx.b
    g = b.grid()
    L = len(b.layers)
    lf = 0 if fp["side"] == "F" else L - 1
    kf = max(1, int(round(p["cell_mm"] / g.res)))
    m = max(1, int(p["coarse"]))
    kc = kf * m
    pts = np.vstack([np.asarray(r, dtype=float) for pad in pads for poly in pad["poly"] for r in poly])
    mar = p["margin_mm"] if margin is None else margin
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
    A, cap, gconv, idx, cs = window_matrix(b.stack, cu, thermal._inside(ctx)[sl], vias, lf, kf, m, g.res, h_part, h_far)
    n = A.shape[0]
    fidx = idx[lf]
    x0, y0 = g.x0 + j0 * g.res, g.y0 + i0 * g.res
    sub = geom.Grid.__new__(geom.Grid)
    sub.res, sub.x0, sub.y0, sub.nx, sub.ny = g.res, x0, y0, j1 - j0, i1 - i0
    covs = []
    for pad in pads:
        mk = np.zeros(sub.shape, bool)
        geom.raster_rings(sub, [ring for poly in pad["poly"] for ring in poly], mk)
        covs.append(_blockmean(mk, kf) * (fidx >= 0))
    surf = (chip["L"] * chip["W"] + 2 * chip["H"] * (chip["L"] + chip["W"])) * 1e-6
    gap = max(math.hypot(pads[0]["x"] - pads[1]["x"], pads[0]["y"] - pads[1]["y"]), 0.1) * 1e-3
    g_body = chip["k"] * chip["W"] * chip["H"] * 1e-6 / gap
    g_paste = p["paste_k"] * paste_scale / (p["paste_um"] * 1e-6) * cs * cs            # per fully covered cell
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
    A2 = (sp.bmat([[A, None], [None, sp.csr_matrix((2, 2))]]).tocsr() +
          sp.coo_matrix((vals, (rows, cols)), shape=(n + 2, n + 2)).tocsr() + sp.diags(dg)).tocsr()
    g2 = np.concatenate([gconv, [h_part * surf / 2, h_part * surf / 2]])
    c2 = np.concatenate([cap, [chip["C"] / 2, chip["C"] / 2]])
    # nodes on the edge of the window: one ring of coarse cells in every layer
    bnodes, blayer, bx, by = [], [], [], []
    for l, ix in enumerate(idx):
        k = kf if l == lf else kc
        w = m if l == lf else 1
        edge = np.zeros(ix.shape, bool)
        edge[:w, :] = edge[-w:, :] = edge[:, :w] = edge[:, -w:] = True
        ii, jj = np.nonzero(edge & (ix >= 0))
        bnodes.append(ix[ii, jj])
        blayer.append(np.full(len(ii), l))
        bx.append(x0 + (jj + 0.5) * k * g.res)
        by.append(y0 + (ii + 0.5) * k * g.res)
    wts = []
    for cov in covs:
        w = np.zeros(n + 2)
        w[fidx[cov > 0]] = cov[cov > 0] / max(cov.sum(), 1e-12)
        wts.append(w)
    return {"A": A2, "cap": c2, "gconv": g2, "n": n, "bnodes": np.concatenate(bnodes), "blayer": np.concatenate(blayer),
            "bx": np.concatenate(bx), "by": np.concatenate(by), "wts": wts, "x": float(pts[:, 0].mean()), "y": float(pts[:, 1].mean())}


def _board_weights(M, layer, x, y):
    """Sparse matrix (points x board nodes): bilinear weights of the board model's nodes at points of
    one copper layer each; cells outside the board are left out and the rest renormalised."""
    mesh = M["mesh"]
    cs = mesh["cs"]
    L, ny, nx = mesh["frac"].shape
    fx, fy = (x - mesh["g"].x0) / cs - 0.5, (y - mesh["g"].y0) / cs - 0.5
    j0, i0 = np.floor(fx).astype(int), np.floor(fy).astype(int)
    tx, ty = fx - j0, fy - i0
    rows, cols, vals = [], [], []
    tot = np.zeros(len(x))
    for di, dj, w in ((0, 0, (1 - tx) * (1 - ty)), (0, 1, tx * (1 - ty)), (1, 0, (1 - tx) * ty), (1, 1, tx * ty)):
        i, j = np.clip(i0 + di, 0, ny - 1), np.clip(j0 + dj, 0, nx - 1)
        node = M["idx"][layer, i, j]
        ok = (node >= 0) & (w > 0)
        rows.append(np.flatnonzero(ok))
        cols.append(node[ok])
        vals.append(w[ok])
        np.add.at(tot, np.flatnonzero(ok), w[ok])
    rows, cols, vals = np.concatenate(rows), np.concatenate(cols), np.concatenate(vals)
    lone = tot <= 0                                   # a point with no board cell around it: the nearest board cell
    if lone.any():
        bm = mesh["board"]
        bi, bj = np.nonzero(bm)
        for k in np.flatnonzero(lone):
            q = int(np.argmin((bi - fy[k]) ** 2 + (bj - fx[k]) ** 2))
            rows = np.r_[rows, k]
            cols = np.r_[cols, M["idx"][layer[k]][bi[q], bj[q]]]
            vals = np.r_[vals, 1.0]
            tot[k] = 1.0
    return sp.csr_matrix((vals / tot[rows], (rows, cols)), shape=(len(x), M["A"].shape[0]))


def chip_lag(Wm, M=None, tau_board=None):
    """Lag times of the two joints (pad, termination, mean per pad). With a board model the window's
    edge is held at the board's lag time there; without one the edge is insulated."""
    n2 = Wm["A"].shape[0]
    if M is not None and tau_board is not None and len(Wm["bnodes"]):
        tb = _board_weights(M, Wm["blayer"], Wm["bx"], Wm["by"]) @ tau_board
        free = np.ones(n2, bool)
        free[Wm["bnodes"]] = False
        A = Wm["A"]
        rhs = Wm["cap"][free] - A[free][:, Wm["bnodes"]] @ tb
        tau = np.empty(n2)
        tau[Wm["bnodes"]] = tb
        tau[free] = _lu(A[free][:, free]).solve(rhs)
    else:
        tau = _lu(Wm["A"]).solve(Wm["cap"])
    n = Wm["n"]
    return [{"tau_pad_s": float(w @ tau), "tau_term_s": float(tau[n + k]), "tau_joint_s": (float(w @ tau) + float(tau[n + k])) / 2}
            for k, w in enumerate(Wm["wts"])]


def chip_transient(Wm, I, liquidus):
    """The window marched through the oven profile of the board run: time (s) at which each joint
    crosses liquidus, or None. The window's edge follows the board model's temperatures."""
    M, t, r = I["M"], I["t"], I["r"]
    dt = t[0]
    n2 = Wm["A"].shape[0]
    free = np.ones(n2, bool)
    free[Wm["bnodes"]] = False
    A = Wm["A"]
    Aff, Afb = A[free][:, free], A[free][:, Wm["bnodes"]]
    Wb = _board_weights(M, Wm["blayer"], Wm["bx"], Wm["by"])
    cap, gc = Wm["cap"][free], Wm["gconv"][free]
    lu = _lu(Aff + sp.diags(cap / dt))
    pfn, _ = profile_fn(r, I["shift"])
    v = M["belt_mm_s"]
    cx, cy, bm = M["cx"], M["cy"], M["mesh"]["board"]
    trav = r["travel"] or "x"
    axis, here = (cx, Wm["x"]) if "x" in trav else (cy, Wm["y"])
    lead = axis[bm].min() if trav.startswith("-") else axis[bm].max()
    off = abs(lead - here) / v
    T = np.full(int(free.sum()), float(r["start_c"]))
    sel = np.cumsum(free) - 1
    n = Wm["n"]
    probes = []
    for k, w in enumerate(Wm["wts"]):
        wf = w[free]
        e = np.zeros(len(T))
        e[sel[n + k]] = 1.0
        probes.append((wf + e) / 2)
    P = np.array(probes)
    hist = np.empty((len(t), len(probes)))
    for s in range(len(t)):
        T = lu.solve(cap / dt * T + gc * pfn(t[s] - off) - Afb @ (Wb @ I["hist"][s]))
        hist[s] = P @ T
    out = []
    for k in range(len(probes)):
        y = hist[:, k]
        up = np.flatnonzero((y[1:] >= liquidus) & (y[:-1] < liquidus))
        if not len(up):
            out.append(None)
            continue
        s = int(up[0])
        out.append(float(t[s] + (liquidus - y[s]) / (y[s + 1] - y[s]) * dt))
    return out


def chips_in_scope(ctx):
    """Two-terminal chips to rank: [(footprint, chip, [pad a, pad b])]."""
    b, cfg = ctx.b, ctx.cfg
    top = packages.size_rank(cfg["dfm"]["chips"]["max_package"])
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


def _rings(pad):
    return [ring for poly in pad["poly"] for ring in poly]


def pad_copper(ctx, pad, plane_cache):
    """What leaves a pad: width of copper of its net crossing a ring just outside it (mm), the share of
    the ring in copper, vias in the pad, and vias of its net into a plane within reach."""
    b = ctx.b
    p = ctx.cfg["dfm"]["chips"]
    g = b.grid()
    ring_px = max(2, int(round(p["ring_mm"] / g.res)))
    rings = _rings(pad)
    a = np.vstack([np.asarray(r, dtype=float) for r in rings])
    x0, y0, x1, y1 = a[:, 0].min(), a[:, 1].min(), a[:, 0].max(), a[:, 1].max()
    mm = (ring_px + 3) * g.res
    i0, j0 = g.ij(x0 - mm, y0 - mm)
    i1, j1 = g.ij(x1 + mm, y1 + mm)
    i0, j0 = max(int(i0), 0), max(int(j0), 0)
    i1, j1 = min(int(i1) + 1, g.ny), min(int(j1) + 1, g.nx)
    sub = geom.Grid.__new__(geom.Grid)
    sub.res, sub.x0, sub.y0, sub.nx, sub.ny = g.res, g.x0 + j0 * g.res, g.y0 + i0 * g.res, j1 - j0, i1 - i0
    M = geom.raster_rings(sub, rings)
    d = ndimage.distance_transform_edt(~M)
    ring = (d > ring_px - 0.5) & (d <= ring_px + 0.5)
    lab = b.raster(pad["layer_ids"][0])[i0:i1, j0:j1]
    n = int(ring.sum())
    cu = int((ring & (lab == pad["net_id"])).sum()) if pad["net_id"] else 0
    V = b.vias
    reach = p["via_reach_mm"]
    same = V["net"] == pad["net_id"] if pad["net_id"] else np.zeros(len(V["x"]), bool)
    inbox = np.flatnonzero(same & (V["x"] >= x0) & (V["x"] <= x1) & (V["y"] >= y0) & (V["y"] <= y1))
    in_pad = int(geom.points_in_rings(V["x"][inbox], V["y"][inbox], rings).sum()) if len(inbox) else 0
    close = np.flatnonzero(same & (V["x"] >= x0 - reach) & (V["x"] <= x1 + reach) & (V["y"] >= y0 - reach) & (V["y"] <= y1 + reach))
    planes = 0
    for k in close:
        k = int(k)
        if k not in plane_cache:
            plane_cache[k] = bool(b.zone_connected_layers(k))
        planes += plane_cache[k]
    return {"pad": pad["num"], "net": pad["net"], "entry_mm": round(cu * g.res, 3), "share": round(cu / n if n else 0.0, 3),
            "vias_in_pad": in_pad, "plane_vias": planes}


def _sig(x, n=2):
    """Two significant digits; anything below 5 ms is 0 (the model is not good to that)."""
    return 0.0 if abs(x) < 0.005 else round(x, n - 1 - int(math.floor(math.log10(abs(x)))))


def chips(ctx, fabr, internal=None, exposure=None, log=print):
    """MF8. A ranking of the two-terminal chips by the time between their two joints reaching liquidus,
    with the copper and land-pattern geometry of each beside it."""
    b, cfg = ctx.b, ctx.cfg
    p = cfg["dfm"]["chips"]
    r = cfg["reflow"]
    h_fab = fabr.get("h_w_m2k", 50.0)
    h_top, h_bot = r["h_top"] or h_fab, r["h_bot"] or h_fab
    M = internal["M"] if internal else None
    tau_b = internal["tau"]["base"] if internal else None
    scope = chips_in_scope(ctx)
    skipped = 0
    if p["max_parts"] and len(scope) > p["max_parts"]:
        skipped = len(scope) - p["max_parts"]
        scope = scope[: p["max_parts"]]
    for l in range(len(b.layers)):                 # loaded before the workers start, so they share them
        b.raster(l)
    thermal._inside(ctx)
    nw = par.workers(cfg)

    def first(k):
        fp, chip, pads = scope[k]
        front = fp["side"] == "F"
        hp, hf = (h_top, h_bot) if front else (h_bot, h_top)
        try:
            j = chip_lag(chip_window(ctx, fp, chip, pads, hp, hf, p), M, tau_b)
        except Exception as e:                       # a window that cannot be solved is reported, not hidden
            return {"problem": f"{fp['ref']}: {e}"}
        d = j[0]["tau_joint_s"] - j[1]["tau_joint_s"]
        a, z = (1, 0) if d > 0 else (0, 1)
        cache = {}
        geo = [pad_copper(ctx, pad, cache) for pad in pads]
        for q in geo:
            q["exposure"] = round((exposure or {}).get((fp["ref"], q["pad"]), 1.0), 3)
        # the belt carries one pad into each zone before the other: by the pads' distance along the belt
        # over its speed. Added to the lag difference for the direction the board travels; with the
        # direction not stated, the ranking value excludes it and its range covers all four
        v = r["zone_mm"] / r["zone_s"]
        belts = {dname: sgn * (pads[a][ax] - pads[z][ax]) / v for dname, ax, sgn in (("x", "x", 1), ("-x", "x", -1), ("y", "y", 1), ("-y", "y", -1))}
        belt = belts[r["travel"]] if r["travel"] else 0.0
        return {"k": k, "ref": fp["ref"], "package": chip["code"], "kind": chip["kind"], "x": fp["x"], "y": fp["y"], "belt_s": round(belt, 3),
                "belt_max_s": round(max(abs(x_) for x_ in belts.values()), 3), "lag_s": abs(d), "_belts": belts,
                "layer": b.layers[pads[0]["layer_ids"][0]], "dt_s": abs(d) + belt, "first_pad": pads[a]["num"],
                "first_net": pads[a]["net"], "last_pad": pads[z]["num"], "last_net": pads[z]["net"],
                "first": a, "last": z, "pads": geo, "land": packages.land_pattern(chip, pads)}
    rows, problems = [], []
    for x in par.pmap(first, range(len(scope)), nw):
        (problems if "problem" in x else rows).append(x["problem"] if "problem" in x else x)
    # parts with the same layout and nearly the same result share a row
    rows.sort(key=lambda x: -x["dt_s"])
    out = []
    for x in rows:
        sig = (x["package"], x["kind"], x["first_net"] == "GND", x["last_net"] == "GND",
               tuple((q["entry_mm"], q["vias_in_pad"], q["plane_vias"]) for q in sorted(x["pads"], key=lambda q: q["pad"] != x["first_pad"])))
        for gq in out:
            if gq["_sig"] == sig and abs(gq["dt_s"] - x["dt_s"]) <= 0.15 * gq["dt_s"] + 0.02:
                gq["refs"].append(x["ref"])
                gq["n"] += 1
                break
        else:
            out.append(dict(x, refs=[x["ref"]], n=1, _sig=sig))
    liq = r["liquidus_c"]

    def detail(i):
        x = out[i]
        fp, chip, pads = scope[x["k"]]
        hp, hf = (h_top, h_bot) if fp["side"] == "F" else (h_bot, h_top)
        res, vals, probs = {"variants": {}}, [x["lag_s"]], []
        for name, kw, tkey, hs in (("paste x 0.5", {"paste_scale": 0.5}, "base", 1.0), ("paste x 2", {"paste_scale": 2.0}, "base", 1.0),
                                   ("window x 2", {"margin": 2 * p["margin_mm"]}, "base", 1.0),
                                   ("film coefficient x 0.6", {}, "film coefficient x 0.6", 0.6),
                                   ("film coefficient x 1.6", {}, "film coefficient x 1.6", 1.6)):
            tb = internal["tau"].get(tkey) if internal else None
            if internal and tb is None:
                continue
            try:
                jv = chip_lag(chip_window(ctx, fp, chip, pads, hp * hs, hf * hs, p, **kw), M, tb)
                dv = jv[x["last"]]["tau_joint_s"] - jv[x["first"]]["tau_joint_s"]
                res["variants"][name] = round(dv, 3)
                vals.append(dv)
            except Exception as e:
                probs.append(f"{fp['ref']} ({name}): {e}")
        # the range: the parameter variants, plus the belt term -- for the stated direction, or over all
        # four when the direction is not known
        if r["travel"]:
            res["dt_range_s"] = [_sig(min(vals) + x["belt_s"]), _sig(max(vals) + x["belt_s"])]
        else:
            res["dt_range_s"] = [_sig(max(min(vals) - x["belt_max_s"], 0.0)), _sig(max(vals) + x["belt_max_s"])]
        if internal and internal.get("hist") is not None:
            # a consistency check of the ranking value, not a second figure: the same window marched
            # through the oven (with the board model's direction) against the lag difference plus that
            # direction's belt term. The two are different computations of one quantity (a first moment
            # and a transient); a disagreement larger than the parameter sensitivity is reported
            trav = r["travel"] or "x"
            try:
                tc = chip_transient(chip_window(ctx, fp, chip, pads, hp, hf, p), internal, liq)
                if None not in tc:
                    oven = tc[x["last"]] - tc[x["first"]]
                    est = x["lag_s"] + x["_belts"][trav]
                    spread = max(vals) - min(vals)
                    res["oven_check"] = {"direction": trav, "oven_s": _sig(oven), "lag_plus_belt_s": _sig(est),
                                         "differ_s": _sig(oven - est), "sensitivity_s": _sig(spread)}
                    if abs(oven - est) > max(spread, 0.05):
                        probs.append(f"{fp['ref']}: the oven transient gives {oven:.3g} s, the lag difference plus the belt term "
                                     f"{est:.3g} s (travel {trav}): they differ by more than the parameter sensitivity ({spread:.2g} s)")
                        # the two computations disagree: the value is not resolved; its range covers both
                        res["oven_check"]["agrees"] = False
                        d_ = oven - est
                        res["dt_range_s"] = [_sig(min(res["dt_range_s"][0], res["dt_range_s"][0] + d_)),
                                             _sig(max(res["dt_range_s"][1], res["dt_range_s"][1] + d_))]
                    else:
                        res["oven_check"]["agrees"] = True
            except Exception as e:
                probs.append(f"{fp['ref']} (transient): {e}")
        return res, probs
    top = range(min(p["detail"], len(out)))
    for i, (res, probs) in zip(top, par.pmap(detail, top, nw)):
        out[i].update(res)
        problems += probs
    for x in out:
        for k in ("_sig", "k", "first", "last", "_belts"):
            x.pop(k, None)
    for x in out:
        x["dt_s"] = _sig(x["dt_s"])
        x["lag_s"] = _sig(x["lag_s"])
    # the order the model can resolve: each row's range of ranks given its sensitivity (a row without a
    # range uses its value); rows whose ranges overlap are not ordered by this ranking
    lo = np.array([(x.get("dt_range_s") or [x["dt_s"]] * 2)[0] for x in out])
    hi = np.array([(x.get("dt_range_s") or [x["dt_s"]] * 2)[1] for x in out])
    for i, x in enumerate(out):
        x["rank"] = i + 1
        x["rank_range"] = [int(1 + (lo > hi[i]).sum()), int((hi >= lo[i]).sum())]
    resolved = sum(1 for x in out[:p["detail"]] if x["rank_range"][0] == x["rank_range"][1])
    # what explains the ranking: the concordance of each listed geometry feature with dt (0.5 = none,
    # 1 = the feature alone orders the chips as the model does)
    feats = {"a plane via at one pad only": [float(bool(q["pads"][0]["plane_vias"]) != bool(q["pads"][1]["plane_vias"])) for q in out],
             "copper leaving the pads, ratio": [max(e["entry_mm"] for e in q["pads"]) / max(min(e["entry_mm"] for e in q["pads"]), 0.1) for q in out],
             "toe of the pads": [q["land"]["toe_mm"] for q in out]}
    dts = np.array([x["dt_s"] for x in out])

    def conc(f):
        f = np.asarray(f, float)
        num = den = 0.0
        for i in range(len(f)):
            dd, ff = dts[i + 1:] - dts[i], f[i + 1:] - f[i]
            m = dd != 0
            num += float(((dd[m] * ff[m]) > 0).sum() + 0.5 * (ff[m] == 0).sum())
            den += float(m.sum())
        return round(num / den, 2) if den else None
    explained = {k: conc(v) for k, v in feats.items()} if len(out) > 2 else {}
    fs = []
    j = cfg["dfm"]["judge"]
    if j.get("source") and j.get("joint_dt_s") is not None:
        for x in out:
            if x["dt_s"] >= j["joint_dt_s"]:
                rng = f" (sensitivity {x['dt_range_s'][0]:g} to {x['dt_range_s'][1]:g} s)" if x.get("dt_range_s") else ""
                fs.append(finding("chips", "warning", f"{', '.join(x['refs'][:6])}{' ...' if x['n'] > 6 else ''} ({x['package']}): pad "
                                  f"{x['first_pad']} ({x['first_net']}) reaches liquidus {x['dt_s']:g} s before pad {x['last_pad']} "
                                  f"({x['last_net']}){rng}, at least {j['joint_dt_s']:g} s (board config: {j['source']})",
                                  net=x["last_net"], layer=x["layer"], x=x["x"], y=x["y"], value=x["dt_s"], limit=j["joint_dt_s"],
                                  unit="s", part=x["ref"], kind="dt", basis=f"board config: {j['source']}"))
    if out:
        ds = np.array([x["dt_s"] for x in out for _ in range(x["n"])])
        chk = [(x["oven_check"]["lag_plus_belt_s"], x["oven_check"]["oven_s"]) for x in out if x.get("oven_check")]
        log(f"  {len(rows)} chips in {len(out)} layouts" + (f" ({skipped} left out by max_parts)" if skipped else "")
            + f": median {np.median(ds):.2f} s, largest {ds.max():.2g} s ({out[0]['refs'][0]}"
            + (f" and {out[0]['n'] - 1} like it" if out[0]["n"] > 1 else "") + ")")
        if chk:
            log("  check: lag difference plus belt term against the same window marched through the oven: "
                + ", ".join(f"{a:g} / {c:g} s" for a, c in chk[:6]))
        if not r["travel"]:
            log("  the belt direction is not stated ([reflow] travel): the time between the joints excludes the belt term, and its "
                "range covers the four directions")
        log(f"  {resolved} of the first {min(p['detail'], len(out))} rows have a rank of their own; what explains the order: "
            + ", ".join(f"{k} {v}" for k, v in explained.items() if v is not None))
    for msg in problems[:5]:
        log(f"  problem: {msg}")
    params = {"cell_mm": p["cell_mm"], "coarse": p["coarse"], "margin_mm": p["margin_mm"], "paste_k": p["paste_k"],
              "paste_um": p["paste_um"], "h_top": h_top, "h_bot": h_bot, "edge": "board model" if M is not None else "insulated",
              "detail": p["detail"], "liquidus_c": liq,
              "travel": r["travel"] or "not stated: the value excludes the belt term; the range covers the four directions"}
    return {"rule": "chips", "spec": "MF8", "title": "Chips: time between the two joints reaching liquidus (model), with the land geometry",
            "params": params, "findings": fs, "rows": out, "problems": problems,
            "summary": {"chips": len(rows), "layouts": len(out), "not_solved": skipped + len(problems),
                        "top_resolved": f"{resolved} of the first {min(p['detail'], len(out))} rows have a rank of their own; the others "
                                        f"overlap within their sensitivity (rank_range)",
                        "explained_by": explained,
                        "explained_note": "concordance of each feature with dt over all layouts (0.5 none, 1 the feature alone gives the "
                                          "model's order): where a simple feature explains it, the model adds little beyond it",
                        "unresolved": [x["refs"][0] for x in out if x.get("oven_check") and x["oven_check"].get("agrees") is False]}}


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
    off = np.zeros(A.shape[0])
    pick = sp.csr_matrix(([1.0], ([0], [0])), shape=(1, A.shape[0]))
    hist = march(A, cap, gconv, off, lambda x: np.full(np.shape(x), 200.0), t, 25.0, pick)[0][:, 0]
    k = int(round(tau_ref / dt)) - 1
    out.append(_case("step of the air temperature: rise after one time constant", "reflow", float(hist[k] - 25.0),
                     175.0 * (1 - math.exp(-t[k] / tau_ref)), 1.0, "C", "lumped body, tau = C / (h_top + h_bot)"))
    rate = 1.5
    hist = march(A, cap, gconv, off, lambda x: 25.0 + rate * np.asarray(x), t, 25.0, pick)[0][:, 0]
    out.append(_case("ramp of the air temperature: steady lag", "reflow", float(25.0 + rate * t[-1] - hist[-1]),
                     rate * tau_ref * (1 - math.exp(-t[-1] / tau_ref)), 1.0, "C", "lag = ramp rate x tau"))
    out.append(_case("lag time from the steady solve", "reflow", float(_lu(A).solve(cap)[0]), tau_ref, 0.5, "s",
                     "(K + H)^-1 C 1 = C / (h_top + h_bot) for a uniform board"))
    res, kf, m = 0.025, 2, 4
    npx = 160
    cu = [np.ones((npx, npx), bool), np.ones((npx, npx), bool)]
    A, cap, gconv, widx, cs = window_matrix(st, cu, np.ones((npx, npx), bool), None, 0, kf, m, res, h, h)
    tw = _lu(A).solve(cap)
    out.append(_case("window model (50 um over 200 um cells): lag time", "reflow", float(tw[widx[0][40, 40]]), tau_ref, 0.5, "s",
                     "same uniform board: the two mesh sizes must agree with the lumped value"))
    return out
