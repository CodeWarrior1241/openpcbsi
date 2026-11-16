"""G13: steady-state board temperature (2.5-D finite volumes).

One node per copper layer per square cell. In-plane: copper (fraction of the cell covered, from
the raster) in parallel with the adjacent half-dielectrics (FR4 in-plane conductivity). Through
the stack: dielectric conduction in parallel with the via barrels in the cell. Top and bottom
surfaces lose heat to ambient through h (convection + linearised radiation). Heat sources: the
parts' dissipation spread over their courtyards on their side, plus the I^2 R of each rail from
the G10 DC solution. Junction estimate: Tj = T(board under the part) + P theta_jb."""
import math

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from . import geom
from .validate import register

K_CU, K_FR4_IN, K_FR4_THRU = 390.0, 0.8, 0.3


def build(ctx, cell_mm):
    b = ctx.b
    g = b.grid()
    k = max(1, int(round(cell_mm / g.res)))
    cs = k * g.res
    ny, nx = g.shape[0] // k, g.shape[1] // k
    inside = np.zeros(g.shape, bool)
    geom.raster_rings(g, b.outline_rings(), inside)
    board = inside[:ny * k, :nx * k].reshape(ny, k, nx, k).mean(axis=(1, 3)) > 0.5
    frac = []
    for li in range(len(b.layers)):
        r = b.raster(li) != 0
        frac.append(r[:ny * k, :nx * k].reshape(ny, k, nx, k).mean(axis=(1, 3)))
    return {"k": k, "cs": cs, "ny": ny, "nx": nx, "board": board, "frac": np.array(frac), "g": g}


def solve(ctx, mesh, heat, h_top, h_bot, ambient):
    """heat: (L, ny, nx) watts per node. Returns T (L, ny, nx) in C."""
    b = ctx.b
    st = b.stack
    L, ny, nx = mesh["frac"].shape
    cs = mesh["cs"] * 1e-3
    board = mesh["board"]
    idx = -np.ones((L, ny, nx), np.int64)
    nb = int(board.sum())
    for l in range(L):
        idx[l][board] = l * nb + np.arange(nb)
    N = L * nb
    rows, cols, vals = [], [], []

    def add(a, c, G):
        rows.extend([a, c]); cols.extend([c, a]); vals.extend([-G, -G])

    # half dielectric thickness belonging to each copper layer
    half = np.zeros(L)
    for d in st.dielectrics:
        if d["above"] in st.index:
            half[st.index[d["above"]]] += d["t"] / 2
        if d.get("below") in st.index:
            half[st.index[d["below"]]] += d["t"] / 2
    for l in range(L):
        t = st.copper[l]["t"] * 1e-3
        Gs = K_CU * t * mesh["frac"][l] + K_FR4_IN * (half[l] * 1e-3 + t * (1 - mesh["frac"][l]))
        for di, dj in ((0, 1), (1, 0)):
            a = idx[l, : ny - di, : nx - dj]
            c = idx[l, di:, dj:]
            ok = (a >= 0) & (c >= 0)
            Gab = 2 / (1 / Gs[: ny - di, : nx - dj][ok] + 1 / Gs[di:, dj:][ok])
            add(a[ok], c[ok], Gab)
    # vias per cell
    V = b.vias
    j = ((V["x"] - mesh["g"].x0) / mesh["cs"]).astype(int)
    i = ((V["y"] - mesh["g"].y0) / mesh["cs"]).astype(int)
    for l in range(L - 1):
        dz = (st.copper[l + 1]["z_top"] - st.copper[l]["z_bot"]) * 1e-3
        G = np.full((ny, nx), K_FR4_THRU * cs * cs / dz)
        span = (V["top"] <= l) & (V["bot"] >= l + 1)
        ok = span & (i >= 0) & (i < ny) & (j >= 0) & (j < nx)
        ann = math.pi * ((V["drill"] / 2 + 0.02) ** 2 - (V["drill"] / 2) ** 2) * 1e-6
        np.add.at(G, (i[ok], j[ok]), K_CU * ann[ok] / dz)
        a, c = idx[l][board], idx[l + 1][board]
        add(a, c, G[board])
    rows = np.concatenate([np.atleast_1d(r) for r in rows])
    cols = np.concatenate([np.atleast_1d(c) for c in cols])
    vals = np.concatenate([np.atleast_1d(v) for v in vals])
    A = sp.coo_matrix((vals, (rows, cols)), shape=(N, N)).tocsr()
    diag = -np.asarray(A.sum(axis=1)).ravel()
    gconv = np.zeros(N)
    gconv[idx[0][board]] += h_top * cs * cs
    gconv[idx[L - 1][board]] += h_bot * cs * cs
    A = A + sp.diags(diag + gconv)
    q = np.zeros(N)
    for l in range(L):
        q[idx[l][board]] += heat[l][board]
    T = spla.spsolve(A.tocsc(), q + gconv * ambient)
    out = np.full((L, ny, nx), np.nan)
    for l in range(L):
        out[l][board] = T[idx[l][board]]
    return out


def run(ctx, cache_dir, out_dir, log=print, dc=None):
    b, cfg = ctx.b, ctx.cfg
    p = cfg["thermal"]
    mesh = build(ctx, p["cell_mm"])
    L, ny, nx = mesh["frac"].shape
    heat = np.zeros((L, ny, nx))
    g = mesh["g"]
    parts = []
    for ref, spec in p["parts"].items():
        f = b.fps.get(ref)
        if not f:
            continue
        m = np.zeros(g.shape, bool)
        rings = [r for poly in f["courtyard"] for r in poly]
        if rings:
            geom.raster_rings(g, rings, m)
        else:
            x0, y0, x1, y1 = f["bbox"]
            geom.raster_rings(g, [[(x0, y0), (x1, y0), (x1, y1), (x0, y1)]], m)
        k = mesh["k"]
        cov = m[:ny * k, :nx * k].reshape(ny, k, nx, k).mean(axis=(1, 3))
        cov = cov * mesh["board"]
        if cov.sum() == 0:
            continue
        l = 0 if f["side"] == "F" else L - 1
        heat[l] += spec["power_w"] * cov / cov.sum()
        parts.append((ref, spec, l, cov > 0))
    # Joule heating of the rails (G10)
    joule = 0.0
    if p.get("joule", True) and cfg["rails"]:
        from . import pdn
        for name, rail in cfg["rails"].items():
            r, sol = pdn.dc_drop(ctx, name, rail)
            if sol is None:
                continue
            net, Vf = sol
            for e in net.br:
                va, vb = Vf[e["a"]], Vf[e["b"]]
                if np.isnan(va) or np.isnan(vb):
                    continue
                P = (va - vb) ** 2 / e["R"]
                if P <= 0:
                    continue
                na, nbn = net.nodes[e["a"]], net.nodes[e["b"]]
                x, y = (na[2] + nbn[2]) / 2, (na[3] + nbn[3]) / 2
                jj, ii = int((x - g.x0) / mesh["cs"]), int((y - g.y0) / mesh["cs"])
                if 0 <= ii < ny and 0 <= jj < nx:
                    heat[na[1], ii, jj] += P
                    joule += P
    T = solve(ctx, mesh, heat, p["h_top"], p["h_bot"], p["ambient_c"])
    rows = []
    for ref, spec, l, m in parts:
        tb = float(np.nanmax(T[l][m]))
        tj = tb + spec["power_w"] * spec.get("theta_jb", 0.0)
        fp = b.fps.get(ref, {})
        rows.append({"ref": ref, "power_w": spec["power_w"], "T_board_c": round(tb, 1),
                     "theta_jb": spec.get("theta_jb"), "Tj_est_c": round(tj, 1), "x": fp.get("x"), "y": fp.get("y")})
        log(f"  {ref:5s} {spec['power_w']:5.2f} W  board under part {tb:6.1f} C  Tj ~ {tj:6.1f} C")
    tot = float(heat.sum())
    log(f"  total {tot:.2f} W (of which I^2R of the rails {joule:.3f} W); board max {np.nanmax(T):.1f} C, "
        f"mean top {np.nanmean(T[0]):.1f} C at h = {p['h_top']}/{p['h_bot']} W/m2K, {p['ambient_c']} C ambient")
    img = None
    try:
        from . import report
        img = report.heatmap(T[0], mesh, b, f"{out_dir}/thermal_top.png", "Temperature, top layer (C)")
    except Exception as e:
        log(f"  (no image: {e})")
    return {"rule": "thermal", "spec": "G13", "title": "Board temperature", "validation": _status(cache_dir),
            "params": p, "parts": rows, "total_w": round(tot, 3), "joule_w": round(joule, 4),
            "T_max_c": round(float(np.nanmax(T)), 2), "T_mean_top_c": round(float(np.nanmean(T[0])), 2),
            "image": img, "metrics": {r["ref"]: {"Tj_est_c": r["Tj_est_c"], "T_board_c": r["T_board_c"]} for r in rows},
            "findings": _judge(rows, ctx), "summary": {"T_max_c": round(float(np.nanmax(T)), 1)}}


def _judge(rows, ctx):
    from . import judge
    return judge.thermal(rows, ctx.cfg, lambda ref: ctx.b.fps.get(ref, {}).get("fields", {}))


def _status(cache_dir):
    from .validate import status
    return status(cache_dir, "thermal")


def validation_cases(ctx):
    """Uniform heat on the whole top layer of the real board mesh: T - Ta = q / (h_top + h_bot)
    (thin board, Biot << 1, adiabatic edges; error from through-thickness conduction < 0.5 %)."""
    from .validate import _case
    mesh = build(ctx, 2.0)
    L, ny, nx = mesh["frac"].shape
    heat = np.zeros((L, ny, nx))
    q = 100.0  # W/m^2
    A = (mesh["cs"] * 1e-3) ** 2
    heat[0][mesh["board"]] = q * A
    T = solve(ctx, mesh, heat, 10.0, 10.0, 25.0)
    ref = 25.0 + q / 20.0
    return [_case("uniform 100 W/m2 on the board, h = 10 + 10 W/m2K: mean temperature", "thermal",
                  float(np.nanmean(T)), ref, 0.5, "C", "energy balance q = (h_top + h_bot) dT (exact for Bi << 1)")]
