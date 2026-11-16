"""G10 DC drop, G11 PDN impedance with a decoupling optimiser, G12 plane resonances and transfer
impedance. One network per supply net, shared by all three:

  - zone copper on each layer as square cells (fraction of copper per cell from the fine raster),
    joined by sheet conductances (DC) or R + jwL plane-pair branches (AC, cavity per square
    L = mu0 h_eff, shunt C = eps A / h to the adjacent reference planes),
  - tracks as series branches R = L / (sigma w t) (+ jwL from the height over the nearest plane),
  - via barrels between the layers where the via is connected (R of the plated annulus, partial L),
  - pads as nodes; track ends, vias and pads joined where they touch.

Ground is ideal (the return plane is the AC reference; DC drop is for the supply net only).
Capacitors: series C-ESR-ESL by value and package, plus GND-side mounting inductance.
"""
import math
import re
from collections import defaultdict

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.spatial import cKDTree

from . import geom

MU0 = 4e-7 * math.pi
EPS0 = 8.8541878128e-12
C0 = 299792458.0


def _pad_dist(p, x, y):
    """Distance from a point to a pad's copper (0 inside)."""
    rings = [r for poly in p["poly"] for r in poly]
    if not rings:
        return math.hypot(p["x"] - x, p["y"] - y)
    if geom.points_in_rings([x], [y], rings)[0]:
        return 0.0
    d = math.inf
    for r in rings:
        a = np.asarray(r)
        for s in np.c_[a, np.roll(a, -1, axis=0)]:
            d = min(d, float(geom.point_seg_dist(x, y, *s)))
    return d


def parse_value(s):
    """'100n' -> 1e-7, '4.7u 6.3V' -> 4.7e-6, '10mR' -> 0.01, '1R' -> 1."""
    if not s:
        return None
    m = re.match(r"\s*([0-9.]+)\s*([pnumkMR]?)([0-9]*)", s.replace("µ", "u"))
    if not m:
        return None
    v = float(m.group(1))
    unit = m.group(2)
    if m.group(3):  # 4u7 style
        v = float(f"{m.group(1)}.{m.group(3)}")
    mult = {"p": 1e-12, "n": 1e-9, "u": 1e-6, "m": 1e-3, "k": 1e3, "M": 1e6, "R": 1.0, "": 1.0}[unit]
    if unit == "m" and "R" in s:
        mult = 1e-3
    return v * mult


def package_of(fpid):
    m = re.search(r"_(0201|0402|0603|0805|1206|1210|2220)_", fpid)
    if m:
        return m.group(1)
    if "Elec" in fpid or "CP_" in fpid:
        return "elec"
    return "0603"


def esr_default(c, pkg):
    if pkg == "elec":
        return 0.1
    if c >= 10e-6:
        return 0.003
    if c >= 1e-6:
        return 0.006
    if c >= 1e-7:
        return 0.015
    return 0.04


class Network:
    """Nodes and branches of one supply net. Branch kinds: 'sheet', 'track', 'via', 'link'."""

    def __init__(self, ctx, net, cell_mm, sigma=5.8e7, plating_um=20.0, ac=False):
        self.ctx, self.b, self.net = ctx, ctx.b, net
        self.nid = self.b.net_id[net]
        self.cell = cell_mm
        self.sigma = sigma
        self.plating = plating_um * 1e-3
        self.nodes = []          # (kind, layer, x, y, extra)
        self.br = []             # dict(a, b, kind, R, L, info)
        self.shunt = defaultdict(lambda: [0.0, 0.0])   # node -> [C, G/omega factor (tand)]
        self.ac = ac
        self._build()

    def add_node(self, kind, li, x, y, extra=None):
        self.nodes.append((kind, li, x, y, extra))
        return len(self.nodes) - 1

    def add_branch(self, a, b, kind, R, L=0.0, info=None):
        if a == b:
            return
        self.br.append({"a": a, "b": b, "kind": kind, "R": max(R, 1e-9), "L": L, "info": info})

    # ------------------------------------------------------------------------------------------
    def _build(self):
        b, ctx, g = self.b, self.ctx, self.b.grid()
        st = b.stack
        k = max(1, int(round(self.cell / g.res)))
        self.k = k
        cs = k * g.res
        self.cell = cs
        self.cell_index = {}
        zl = sorted({b.L[z["layer"]] for z in b.zones if z["net"] == self.net})
        # --- zone cells
        for li in zl:
            m = np.zeros(g.shape, bool)
            for z in b.zones:
                if z["net"] == self.net and z["layer"] == b.layers[li]:
                    geom.raster_rings(g, [r for poly in z["polys"] for r in poly], m)
            ny, nx = g.shape[0] // k, g.shape[1] // k
            frac = m[:ny * k, :nx * k].reshape(ny, k, nx, k).mean(axis=(1, 3))
            idx = -np.ones((ny, nx), np.int64)
            ii, jj = np.nonzero(frac > 0.2)
            t = st.cu(b.layers[li])["t"] * 1e-3
            heff, cpa, tand = self._plane_pair(li)
            for i, j in zip(ii, jj):
                x = g.x0 + (j + 0.5) * cs
                y = g.y0 + (i + 0.5) * cs
                idx[i, j] = self.add_node("cell", li, x, y, frac[i, j])
                if cpa:
                    sh = self.shunt[idx[i, j]]
                    sh[0] += cpa * (cs * 1e-3) ** 2 * frac[i, j]
                    sh[1] = tand
            Gs = self.sigma * t
            for di, dj in ((0, 1), (1, 0)):
                a = idx[: ny - di, : nx - dj]
                c = idx[di:, dj:]
                ok = (a >= 0) & (c >= 0)
                fa = frac[: ny - di, : nx - dj][ok]
                fc = frac[di:, dj:][ok]
                for na, nc, f1, f2 in zip(a[ok], c[ok], fa, fc):
                    ff = math.sqrt(f1 * f2)
                    self.add_branch(int(na), int(nc), "sheet", 1.0 / (Gs * ff), MU0 * heff / ff if heff else 0.0,
                                    {"layer": li})
            self.cell_index[li] = (idx, frac)
        # --- pads
        self.pad_node = {}
        for p in b.pads_of(self.net):
            if not p["layer_ids"]:
                continue
            n = self.add_node("pad", p["layer_ids"][0], p["x"], p["y"], f"{p['ref']}.{p['num']}")
            self.pad_node[(p["ref"], p["num"], p["x"], p["y"])] = n
            for li in p["layer_ids"]:
                cn = self.cell_at(li, p["x"], p["y"])
                if cn is not None:
                    self.add_branch(n, cn, "link", self._spread(li, max(p["sx"], p["sy"]) / 2))
        # --- vias: a node per connected layer
        vtl, vzl = ctx.via_track_layers(), ctx.via_zone_layers()
        V = b.vias
        self.via_nodes = {}
        for v in b.vias_of(self.net):
            layers = sorted(vtl[v] | vzl[v])
            if not layers:
                continue
            ns = {}
            for li in layers:
                ns[li] = self.add_node("via", li, V["x"][v], V["y"][v], int(v))
                cn = self.cell_at(li, V["x"][v], V["y"][v])
                if cn is not None and li in vzl[v]:
                    self.add_branch(ns[li], cn, "link", self._spread(li, V["d"][v] / 2))
            r_out = V["drill"][v] / 2 + self.plating
            area = math.pi * (r_out ** 2 - (V["drill"][v] / 2) ** 2) * 1e-6
            for la, lb in zip(layers[:-1], layers[1:]):
                dz = st.barrel(b.layers[la], b.layers[lb]) * 1e-3
                Lv = MU0 * dz / (2 * math.pi) * math.log(max(self._return_dist(V["x"][v], V["y"][v]), 2 * r_out) / (r_out * 1e-0)) \
                    if self.ac else 0.0
                self.add_branch(ns[la], ns[lb], "via", dz / (self.sigma * area), Lv, {"via": int(v), "layers": (la, lb)})
            self.via_nodes[int(v)] = ns
        # --- tracks (split at T junctions)
        T = b.tracks
        ks = b.tracks_of(self.net)
        ends = defaultdict(list)
        for t in ks:
            li = int(T["layer"][t])
            for (x, y) in ((T["x1"][t], T["y1"][t]), (T["x2"][t], T["y2"][t])):
                ends[(li, round(x, 4), round(y, 4))].append(t)
        node_at = {}
        for key in ends:
            li, x, y = key
            node_at[key] = self.add_node("end", li, x, y)
        # extra split points: an end lying on the body of another track
        splits = defaultdict(list)
        for key in ends:
            li, x, y = key
            for t in ks:
                if T["layer"][t] != li:
                    continue
                if (round(T["x1"][t], 4), round(T["y1"][t], 4)) == (x, y) or (round(T["x2"][t], 4), round(T["y2"][t], 4)) == (x, y):
                    continue
                if geom.point_seg_dist(x, y, T["x1"][t], T["y1"][t], T["x2"][t], T["y2"][t]) < T["w"][t] / 2:
                    L2 = T["len"][t] ** 2 or 1
                    s = ((x - T["x1"][t]) * (T["x2"][t] - T["x1"][t]) + (y - T["y1"][t]) * (T["y2"][t] - T["y1"][t])) / L2
                    splits[t].append((min(max(s, 0), 1), node_at[key]))
        for t in ks:
            li = int(T["layer"][t])
            a = node_at[(li, round(T["x1"][t], 4), round(T["y1"][t], 4))]
            c = node_at[(li, round(T["x2"][t], 4), round(T["y2"][t], 4))]
            pts = [(0.0, a)] + sorted(splits.get(t, [])) + [(1.0, c)]
            w = T["w"][t] * 1e-3
            th = st.cu(b.layers[li])["t"] * 1e-3
            Lp = self._track_L(li, T["w"][t]) if self.ac else 0.0
            for (s0, n0), (s1, n1) in zip(pts[:-1], pts[1:]):
                L = (s1 - s0) * T["len"][t] * 1e-3
                self.add_branch(n0, n1, "track", L / (self.sigma * w * th), Lp * L, {"track": int(t)})
        # join track ends to pads, vias and zone cells (a track end touches what lies within w/2)
        hw = defaultdict(float)
        for key, tl in ends.items():
            hw[key] = max(T["w"][t] for t in tl) / 2
        pads_on = defaultdict(list)
        for (ref, num, px, py), pn in self.pad_node.items():
            p = next(q for q in b.pads_by_ref[ref] if q["num"] == num and q["x"] == px and q["y"] == py)
            for li in p["layer_ids"]:
                pads_on[li].append((p, pn))
        for (li, x, y), n in node_at.items():
            w2 = hw[(li, x, y)]
            for p, pn in pads_on.get(li, []):
                if abs(p["x"] - x) <= p["sx"] + p["sy"] and abs(p["y"] - y) <= p["sx"] + p["sy"]:
                    if _pad_dist(p, x, y) <= w2 + 1e-4:
                        self.add_branch(n, pn, "link", 1e-6)
            for v, ns in self.via_nodes.items():
                if li in ns and math.hypot(V["x"][v] - x, V["y"][v] - y) <= V["d"][v] / 2 + w2:
                    self.add_branch(n, ns[li], "link", 1e-6)
            cn = self.cell_at(li, x, y)
            if cn is not None:
                self.add_branch(n, cn, "link", self._spread(li, 0.1))
        # vias inside pads (via in pad, bottom-side caps on vias)
        for v, ns in self.via_nodes.items():
            for (ref, num, px, py), pn in self.pad_node.items():
                p = next(q for q in b.pads_by_ref[ref] if q["num"] == num and q["x"] == px and q["y"] == py)
                for li in p["layer_ids"]:
                    if li in ns and abs(px - V["x"][v]) < p["sx"] + p["sy"] and abs(py - V["y"][v]) < p["sx"] + p["sy"]:
                        rings = [r for poly in p["poly"] for r in poly]
                        if rings and geom.points_in_rings([V["x"][v]], [V["y"][v]], rings)[0]:
                            self.add_branch(ns[li], pn, "link", 1e-6)
        self.N = len(self.nodes)

    # ------------------------------------------------------------------------------------------
    def cell_at(self, li, x, y):
        if li not in self.cell_index:
            return None
        idx, frac = self.cell_index[li]
        g = self.b.grid()
        j = int((x - g.x0) // self.cell)
        i = int((y - g.y0) // self.cell)
        if 0 <= i < idx.shape[0] and 0 <= j < idx.shape[1] and idx[i, j] >= 0:
            return int(idx[i, j])
        # nearest copper cell within one cell
        best = None
        for di in (-1, 0, 1):
            for dj in (-1, 0, 1):
                ii, jj = i + di, j + dj
                if 0 <= ii < idx.shape[0] and 0 <= jj < idx.shape[1] and idx[ii, jj] >= 0:
                    best = int(idx[ii, jj])
        return best

    def _spread(self, li, r_mm):
        t = self.b.stack.cu(self.b.layers[li])["t"] * 1e-3
        return max(math.log(max(self.cell / 2, r_mm * 1.01) / r_mm), 0.05) / (2 * math.pi * self.sigma * t)

    def _plane_pair(self, li):
        """Cavity of a zone on layer li with the adjacent reference planes: (h_eff for L, C per m^2, tand)."""
        b, ctx = self.b, self.ctx
        hs, cpa, td = [], 0.0, []
        for n in b.stack.neighbours(b.layers[li]):
            nli = b.L[n]
            if ctx.refmask(nli).mean() > 0.3:
                h = b.stack.gap(b.layers[li], n) * 1e-3
                er = b.stack.er_between(b.layers[li], n)
                hs.append(h)
                cpa += EPS0 * er / h
                td.append(b.stack.tand_between(b.layers[li], n))
        if not hs:
            return None, 0.0, 0.02
        heff = 1.0 / sum(1.0 / h for h in hs)
        return heff, cpa, float(np.mean(td))

    def _return_dist(self, x, y):
        if not hasattr(self, "_gtree"):
            V = self.b.vias
            m = np.isin(V["net"], self.ctx.ref_ids)
            self._gtree = cKDTree(np.c_[V["x"][m], V["y"][m]]) if m.any() else None
        if self._gtree is None:
            return 1.0
        d, _ = self._gtree.query([x, y])
        return max(d, 0.2)

    def _track_L(self, li, w):
        """Per-metre inductance of a track over its nearest plane (Wheeler/Hammerstad, air)."""
        b = self.b
        hs = [b.stack.gap(b.layers[li], n) for n in b.stack.neighbours(b.layers[li]) if self.ctx.refmask(b.L[n]).mean() > 0.3]
        h = min(hs) if hs else 0.5
        from .xsec import hammerstad_jensen
        z_air, _ = hammerstad_jensen(w, h, 1.0)
        return z_air / C0

    # ------------------------------------------------------------------------------------------
    def matrix(self, f=0.0, extra=None):
        """Nodal admittance matrix (complex if f > 0). extra: list of (node, Y) shunts to ground."""
        if not hasattr(self, "_arr"):
            self._arr = (np.array([e["a"] for e in self.br]), np.array([e["b"] for e in self.br]),
                         np.array([e["R"] for e in self.br]), np.array([e["L"] for e in self.br]))
            sh = list(self.shunt.items())
            self._sh = (np.array([n for n, _ in sh], dtype=np.int64), np.array([v[0] for _, v in sh]),
                        np.array([v[1] for _, v in sh]))
        a, c, R, L = self._arr
        N = self.N
        if f > 0:
            w = 2 * math.pi * f
            y = 1.0 / (R + 1j * w * L)
        else:
            y = 1.0 / R
        Y = sp.coo_matrix((np.r_[-y, -y], (np.r_[a, c], np.r_[c, a])), shape=(N, N)).tocsr()
        diag = -np.asarray(Y.sum(axis=1)).ravel()
        if f > 0:
            diag = diag.astype(complex)
            n, C, td = self._sh
            if len(n):
                np.add.at(diag, n, 1j * w * C * (1 - 1j * td))
        if extra:
            for n, yv in extra:
                diag[n] += yv
        return Y + sp.diags(diag)

    def components(self):
        a = np.array([e["a"] for e in self.br])
        c = np.array([e["b"] for e in self.br])
        A = sp.coo_matrix((np.ones(len(a)), (a, c)), shape=(self.N, self.N))
        return sp.csgraph.connected_components(A, directed=False)[1]


# ==============================================================================================

def rail_nodes(net, rail, ctx):
    """Source, sense and load nodes of a rail."""
    b = ctx.b
    src_pads = [p for p in b.pads_by_ref.get(rail.get("source", ""), []) if p["net"] == rail["net"]]
    loads = [p for r in rail.get("loads", []) for p in b.pads_by_ref.get(r, []) if p["net"] == rail["net"]]
    sense = []
    if rail.get("sense"):
        # the sense net itself may be the rail (pads of parts on both), or a net one part away
        nets_of = {ref: {p["net"] for p in pads} for ref, pads in b.pads_by_ref.items() if not b.fps[ref]["dnp"]}
        for ref, nets in nets_of.items():
            if rail["sense"] in nets and rail["net"] in nets:
                sense += [p for p in b.pads_by_ref[ref] if p["net"] == rail["net"]]
        if not sense:
            # sense net -> series part -> rail: e.g. FB -> divider -> CORE_SNS -> track/zero-ohm -> VCCINT
            for ref, nets in nets_of.items():
                if rail["sense"] in nets and len(nets) == 2:
                    other = (nets - {rail["sense"]}).pop()
                    for ref2, nets2 in nets_of.items():
                        if other in nets2 and rail["net"] in nets2 and ref2 != ref:
                            sense += [p for p in b.pads_by_ref[ref2] if p["net"] == rail["net"]]
    node = lambda p: net.pad_node.get((p["ref"], p["num"], p["x"], p["y"]))
    return ([node(p) for p in src_pads if node(p) is not None], [node(p) for p in sense if node(p) is not None],
            [(node(p), p) for p in loads if node(p) is not None])


def dc_drop(ctx, rail_name, rail, log=print):
    cfg = ctx.cfg["dc"]
    net = Network(ctx, rail["net"], cfg["cell_mm"], cfg["sigma"], cfg["plating_um"])
    src, sense, loads = rail_nodes(net, rail, ctx)
    out = {"rail": rail_name, "net": rail["net"], "nodes": net.N, "branches": len(net.br),
           "v_nom": rail["v"], "current_a": rail["current_a"], "problems": []}
    if not src:
        out["problems"].append(f"no source pad of {rail.get('source')} on {rail['net']}")
        return out, None
    if not loads:
        out["problems"].append("no load pads")
        return out, None
    comp = net.components()
    main = comp[src[0]]
    floating = [p["ref"] + "." + p["num"] for n, p in loads if comp[n] != main]
    if floating:
        out["problems"].append(f"{len(floating)} load pads not connected to the source in copper: {floating[:6]}")
    loads = [(n, p) for n, p in loads if comp[n] == main]
    if not loads:
        out["problems"].append("no load pad is connected to the source")
        return out, None
    fixed = sense[0] if sense and comp[sense[0]] == main else src[0]
    out["regulation"] = "remote sense" if sense and fixed == sense[0] else "local"
    keep = np.flatnonzero(comp == main)
    remap = -np.ones(net.N, np.int64)
    remap[keep] = np.arange(len(keep))
    Y = net.matrix(0.0)[keep][:, keep].tocsc()
    I = np.zeros(len(keep))
    for s in src:
        I[remap[s]] += rail["current_a"] / len(src)
    for n, p in loads:
        I[remap[n]] -= rail["current_a"] / len(loads)
    fi = remap[fixed]
    free = np.ones(len(keep), bool)
    free[fi] = False
    fidx = np.flatnonzero(free)
    V = np.zeros(len(keep))
    V[fi] = rail["v"]
    rhs = I[fidx] - Y[fidx][:, [fi]].toarray().ravel() * rail["v"]
    V[fidx] = spla.spsolve(Y[fidx][:, fidx].tocsc(), rhs)
    Vfull = np.full(net.N, np.nan)
    Vfull[keep] = V
    win = rail["v"] * rail.get("tol_pct", 3.0) / 100
    lv = np.array([Vfull[n] for n, p in loads])
    worst_i = int(np.argmin(lv))
    out.update({
        "v_source": round(float(np.mean([Vfull[s] for s in src])), 6),
        "v_load_min": round(float(lv.min()), 6), "v_load_max": round(float(lv.max()), 6),
        "drop_mv": round(1000 * (rail["v"] - float(lv.min())), 3),
        "drop_pct_of_window": round(100 * (rail["v"] - float(lv.min())) / win, 2),
        "worst_ball": loads[worst_i][1]["ref"] + "." + loads[worst_i][1]["num"],
        "loads": len(loads),
    })
    # branch currents
    Ib = []
    for e in net.br:
        va, vb = Vfull[e["a"]], Vfull[e["b"]]
        if np.isnan(va) or np.isnan(vb):
            continue
        Ib.append((e, (va - vb) / e["R"]))
    b = ctx.b
    vias = defaultdict(float)
    tracks = []
    for e, i in Ib:
        if e["kind"] == "via":
            vias[e["info"]["via"]] = max(vias[e["info"]["via"]], abs(i))
        elif e["kind"] == "track":
            t = e["info"]["track"]
            li = int(b.tracks["layer"][t])
            w = b.tracks["w"][t]
            th = b.stack.cu(b.layers[li])["t"]
            A = w * th  # mm^2
            J = abs(i) / A
            A_mil2 = A / (0.0254 ** 2)
            kk = 0.048 if b.stack.is_outer(b.layers[li]) else 0.024
            dT = (abs(i) / (kk * A_mil2 ** 0.725)) ** (1 / 0.44) if abs(i) > 0 else 0.0
            tracks.append({"track": t, "layer": b.layers[li], "w": round(w, 3), "I": round(abs(i), 4),
                           "J_A_mm2": round(J, 2), "dT_ipc2221_C": round(dT, 2),
                           "x": round((b.tracks["x1"][t] + b.tracks["x2"][t]) / 2, 3),
                           "y": round((b.tracks["y1"][t] + b.tracks["y2"][t]) / 2, 3)})
    tracks.sort(key=lambda r: -r["J_A_mm2"])
    # sheet current density in cells (A/mm of width / thickness)
    hot = []
    for li, (idx, frac) in net.cell_index.items():
        th = b.stack.cu(b.layers[li])["t"]
        Vc = np.where(idx >= 0, Vfull[np.clip(idx, 0, None)], np.nan)
        gy, gx = np.gradient(Vc, net.cell * 1e-3)
        Js = np.hypot(gx, gy) * net.sigma * th * 1e-3 * 1e-3 / th      # A/mm^2
        Js = np.nan_to_num(Js)
        top = np.argsort(Js.ravel())[::-1][:5]
        for q in top:
            i, j = np.unravel_index(q, Js.shape)
            if Js[i, j] > 0:
                g = b.grid()
                hot.append({"layer": b.layers[li], "x": round(g.x0 + (j + 0.5) * net.cell, 2),
                            "y": round(g.y0 + (i + 0.5) * net.cell, 2), "J_A_mm2": round(float(Js[i, j]), 2)})
    hot.sort(key=lambda r: -r["J_A_mm2"])
    vl = sorted(vias.items(), key=lambda kv: -kv[1])
    out["max_via_current_a"] = round(vl[0][1], 4) if vl else 0.0
    out["vias_carrying"] = len(vl)
    out["top_vias"] = [{"x": round(b.vias["x"][v], 3), "y": round(b.vias["y"][v], 3), "I": round(i, 4),
                        "drill_mm": round(float(b.vias["drill"][v]), 3)} for v, i in vl[:20]]
    out["top_tracks"] = tracks[:20]
    out["max_track_J_A_mm2"] = tracks[0]["J_A_mm2"] if tracks else 0.0
    out["max_track_dT_C"] = max((t["dT_ipc2221_C"] for t in tracks), default=0.0)
    out["hot_spots"] = hot[:5]
    out["ball_voltages"] = {p["ref"] + "." + p["num"]: round(float(Vfull[n]), 6) for n, p in loads}
    out["budget"], out["neck"] = drop_budget(net, Ib, fixed, loads[worst_i][0], b)
    if ctx.cfg["dc"].get("solve_ground", True):
        try:
            out.update(ground_rise(ctx, rail, rail["current_a"]))
            if "ground_rise_mv" in out:
                out["drop_with_ground_mv"] = round(out["drop_mv"] + out["ground_rise_mv"], 4)
        except Exception as e:
            out["ground"] = f"not solved: {e}"
    return out, (net, Vfull)


def ground_rise(ctx, rail, current):
    """#27: the ground return, solved like the rail (not assumed ideal). The rail's current enters the
    ground net at the loads' ground pads and leaves at the regulator's ground pads; the result is the
    rise of the worst load's ground above the regulator's. Needs `regulator` in the rail (config or
    detection); a coarse cell keeps whole-board ground planes affordable."""
    reg = rail.get("regulator")
    gnd = (ctx.cfg.reference_nets or ["GND"])[0]
    b = ctx.b
    if not reg or reg not in b.pads_by_ref or gnd not in b.net_id:
        return {"ground": "not solved: the rail names no regulator" if not reg else f"not solved: no {gnd} on {reg}"}
    cell = ctx.cfg["dc"].get("ground_cell_mm", 0.5)
    net = Network(ctx, gnd, cell, ctx.cfg["dc"]["sigma"], ctx.cfg["dc"]["plating_um"])
    node = lambda p: net.pad_node.get((p["ref"], p["num"], p["x"], p["y"]))
    src = [node(p) for p in b.pads_by_ref[reg] if p["net"] == gnd]
    snk = [node(p) for r in rail.get("loads", []) for p in b.pads_by_ref.get(r, []) if p["net"] == gnd]
    src, snk = [n for n in src if n is not None], [n for n in snk if n is not None]
    if not src or not snk:
        return {"ground": "not solved: no ground pads reached"}
    comp = net.components()
    main = comp[src[0]]
    snk = [n for n in snk if comp[n] == main]
    if not snk:
        return {"ground": "not solved: the loads' ground pads are not connected to the regulator's"}
    keep = np.flatnonzero(comp == main)
    remap = -np.ones(net.N, np.int64)
    remap[keep] = np.arange(len(keep))
    Y = net.matrix(0.0)[keep][:, keep].tocsc()
    I = np.zeros(len(keep))
    for n in snk:
        I[remap[n]] += current / len(snk)
    free = np.ones(len(keep), bool)
    for n in src:                              # every ground pad of the regulator is the 0 V reference
        if comp[n] == main:
            free[remap[n]] = False
    fidx = np.flatnonzero(free)
    V = np.zeros(len(keep))
    V[fidx] = spla.spsolve(Y[fidx][:, fidx].tocsc(), I[fidx])
    rise = max(V[remap[n]] for n in snk)
    return {"ground": "solved", "ground_rise_mv": round(float(rise) * 1000, 4), "ground_net": gnd, "ground_cell_mm": cell}


def drop_budget(net, Ib, src, dst, b):
    """#31: where the drop is spent. Along any path from the regulated node to the worst load the
    branch drops sum to the total drop exactly; the path taken is the one carrying the most current
    (widest-path search). Drops are grouped by element (track / via / plane per layer) in path order;
    `neck` is the track on that path with the highest current per mm of width."""
    import heapq
    adj = defaultdict(list)
    for e, i in Ib:
        if i > 0:
            adj[e["a"]].append((e["b"], e, i))
        elif i < 0:
            adj[e["b"]].append((e["a"], e, -i))
    best = {src: math.inf}
    prev = {}
    h = [(-math.inf, src)]
    while h:
        negw, u = heapq.heappop(h)
        if u == dst:
            break
        if -negw < best.get(u, 0):
            continue
        for v, e, i in adj[u]:
            w = min(-negw, i)
            if w > best.get(v, 0):
                best[v] = w
                prev[v] = (u, e, i)
                heapq.heappush(h, (-w, v))
    if dst not in prev:
        return [], None
    path = []
    v = dst
    while v != src:
        u, e, i = prev[v]
        path.append((e, i))
        v = u
    path.reverse()
    groups, neck = [], None
    for e, i in path:
        li = net.nodes[e["a"]][1]
        what = e["kind"] if e["kind"] == "via" else f"{e['kind']} {b.layers[li] if li is not None and li >= 0 else ''}".strip()
        mv = 1000 * i * e["R"]
        if groups and groups[-1]["what"] == what:
            groups[-1]["mv"] += mv
            groups[-1]["elements"] += 1
        else:
            groups.append({"what": what, "mv": mv, "elements": 1})
        if e["kind"] == "track":
            t = e["info"]["track"]
            w = float(b.tracks["w"][t])
            if neck is None or i / w > neck["A_per_mm"]:
                neck = {"layer": b.layers[int(b.tracks["layer"][t])], "w_mm": round(w, 3), "I_a": round(i, 4),
                        "A_per_mm": round(i / w, 3), "x": round(float((b.tracks["x1"][t] + b.tracks["x2"][t]) / 2), 3),
                        "y": round(float((b.tracks["y1"][t] + b.tracks["y2"][t]) / 2), 3)}
    for g in groups:
        g["mv"] = round(g["mv"], 4)
    return groups, neck


# ==============================================================================================
# G11 / G12

def cap_models(ctx, rail_net):
    b, cfg = ctx.b, ctx.cfg["pdn"]
    caps = []
    for ref, pads in b.pads_by_ref.items():
        if not ref.startswith("C") or b.fps[ref]["dnp"]:
            continue
        nets = {p["net"] for p in pads}
        if rail_net not in nets or not (nets & set(ctx.cfg.reference_nets)):
            continue
        f = b.fps[ref]
        c = parse_value(f["value"])
        if not c:
            continue
        pkg = package_of(f["fpid"])
        esl = cfg["esl_nh"].get(pkg, 1.5 if pkg == "elec" else 0.6) * 1e-9
        pp = [p for p in pads if p["net"] == rail_net][0]
        gp = [p for p in pads if p["net"] in ctx.cfg.reference_nets][0]
        # GND-side mounting: nearest GND via to the GND pad, barrel to the nearest GND plane
        V = b.vias
        gv = np.flatnonzero(np.isin(V["net"], ctx.ref_ids))
        if len(gv):
            d = np.hypot(V["x"][gv] - gp["x"], V["y"][gv] - gp["y"])
            k = gv[int(np.argmin(d))]
            li = gp["layer_ids"][0]
            planes = [n for n in range(len(b.layers)) if ctx.refmask(n).mean() > 0.3]
            hz = min(b.stack.barrel(b.layers[li], b.layers[n]) for n in planes) if planes else 0.2
            s = max(math.hypot(pp["x"] - gp["x"], pp["y"] - gp["y"]), 0.3)
            l_mount = MU0 * hz * 1e-3 / (2 * math.pi) * math.log(s / (V["drill"][k] / 2)) + MU0 * float(d.min()) * 1e-3 * 0.1
        else:
            l_mount = 0.5e-9
        caps.append({"ref": ref, "value": f["value"], "C": c, "pkg": pkg, "ESR": esr_default(c, pkg),
                     "ESL": esl + l_mount, "pad": pp})
    return caps


def zcap(cap, w):
    return cap["ESR"] + 1j * w * cap["ESL"] + 1.0 / (1j * w * cap["C"])


def load_step(freqs, Z, step_a, rise_ns=1.0, n=1 << 15):
    """Voltage droop for a current step: v(t) = IFFT(Z(f) I(f)) on a uniform grid; Z between the log
    samples by interpolation of magnitude and phase, below the first sample held at its value (the
    regulator's DC regulation is not modelled, so the droop is the transient part)."""
    fmax = freqs[-1]
    fu = np.linspace(0, fmax, n)
    lm = np.interp(np.log(np.maximum(fu, freqs[0])), np.log(freqs), np.log(np.abs(Z)))
    ph = np.interp(np.log(np.maximum(fu, freqs[0])), np.log(freqs), np.unwrap(np.angle(Z)))
    Zu = np.exp(lm) * np.exp(1j * ph)
    dt = 1 / (2 * fmax)
    t = np.arange(2 * (n - 1)) * dt
    sig = rise_ns * 1e-9 / 2.563                     # Gaussian edge with the given 10-90 % rise time
    th = t[-1] / 2                                   # step up, then back down half way: a periodic record
    i_t = step_a * ((t >= 0) & (t < th)).astype(float)
    fr = np.fft.rfftfreq(len(t), dt)
    I = np.fft.rfft(i_t) * np.exp(-0.5 * (2 * np.pi * fr * sig) ** 2) * np.exp(-2j * np.pi * fr * 4 * sig)
    v = np.fft.irfft(Zu * I, n=len(t))
    m = t < t[-1] / 2
    return {"load_step_mv": round(float(np.max(np.abs(v[m]))) * 1000, 3), "step_a": step_a, "step_rise_ns": rise_ns}


def pdn_impedance(ctx, rail_name, rail, log=print, optimise=True):
    """Z(f) at the load port of a rail, transfer impedance from the source, resonances, optimiser.

    The bare board (planes, tracks, vias, regulator) is reduced once per frequency to an impedance
    matrix over the port nodes (load balls, capacitor pads, source); capacitors are then shunt
    admittances on that small matrix, so the optimiser does not re-solve the mesh."""
    cfg = ctx.cfg["pdn"]
    net = Network(ctx, rail["net"], cfg["cell_mm"], ac=True)
    src, sense, loads = rail_nodes(net, rail, ctx)
    out = {"rail": rail_name, "net": rail["net"], "nodes": net.N}
    if not src or not loads:
        out["problem"] = "no source or no load pads"
        return out
    comp = net.components()
    main = comp[src[0]]
    keep = np.flatnonzero(comp == main)
    remap = -np.ones(net.N, np.int64)
    remap[keep] = np.arange(len(keep))
    loads = [(n, p) for n, p in loads if comp[n] == main]
    caps = cap_models(ctx, rail["net"])
    for c in caps:
        n = net.pad_node.get((c["pad"]["ref"], c["pad"]["num"], c["pad"]["x"], c["pad"]["y"]))
        c["node"] = n if n is not None and comp[n] == main else None
    caps_on = [c for c in caps if c["node"] is not None]
    shunt_r = parse_value(ctx.b.fps.get(rail.get("source", ""), {}).get("value", "")) or 0.0
    zr = rail.get("zout", {})
    r_reg, l_reg = zr.get("r", 0.002), zr.get("l_nh", 50.0) * 1e-9
    freqs = np.logspace(math.log10(cfg["fmin_hz"]), math.log10(cfg["fmax_hz"]), int(cfg["points"]))
    fmax_t = rail.get("target_fmax_hz", cfg.get("target_fmax_hz", 1e8))
    zt = rail["v"] * rail.get("ripple_pct", 3.0) / 100 / rail.get("step_a", rail["current_a"] / 2)
    ports = sorted(set([n for n, p in loads] + [c["node"] for c in caps_on] + [src[0]]))
    pidx = {n: k for k, n in enumerate(ports)}
    P = np.array([remap[n] for n in ports])
    u = np.zeros(len(ports))
    for n, p in loads:
        u[pidx[n]] += 1.0 / len(loads)
    Zp, modes = [], []
    for f in freqs:
        w = 2 * math.pi * f
        Y = net.matrix(f, [(src[0], 1.0 / (shunt_r + r_reg + 1j * w * l_reg))])[keep][:, keep].tocsc()
        lu = spla.splu(Y)
        rhs = np.zeros((len(keep), len(ports)), complex)
        rhs[P, np.arange(len(ports))] = 1.0
        X = lu.solve(rhs)
        Zp.append(X[P, :])
        modes.append(np.abs(X @ u))
    Zp = np.array(Zp)

    def port_z(capset):
        """|Z| at the load port and |Z21| (source -> load) with the given capacitors."""
        Zl, Z21 = np.zeros(len(freqs), complex), np.zeros(len(freqs), complex)
        for k, f in enumerate(freqs):
            w = 2 * math.pi * f
            Yc = np.zeros(len(ports), complex)
            for c in capset:
                Yc[pidx[c["node"]]] += 1.0 / zcap(c, w)
            Zk = Zp[k]
            # (Z^-1 + Yc)^-1 = Z (I + Yc Z)^-1
            Zn = Zk @ np.linalg.inv(np.eye(len(ports)) + Yc[:, None] * Zk)
            Zl[k] = u @ Zn @ u
            Z21[k] = Zn[pidx[src[0]], :] @ u
        return Zl, Z21

    def with_added(added):
        extra = []
        for k, (cname, count) in enumerate(added.items()):
            cv = parse_value(cname.split()[0])
            pkg = cname.split()[1] if len(cname.split()) > 1 else "0402"
            for m in range(count):
                n = loads[(k * 7 + m * 3) % len(loads)][0]
                extra.append({"ref": f"+{cname}#{m}", "C": cv, "ESR": esr_default(cv, pkg),
                              "ESL": cfg["esl_nh"].get(pkg, 0.5) * 1e-9 + 0.25e-9, "node": n})
        return extra

    # #22: the rail's own band (vendor noise spec), else from fmin to target_fmax_hz
    bh = rail.get("band_hz") or [freqs[0], fmax_t]
    band = (freqs >= bh[0] * 0.999) & (freqs <= bh[1] * 1.001)
    if not band.any():
        band = freqs <= fmax_t
    Z, Z21 = port_z(caps_on)
    mag = np.abs(Z)
    viol = (mag > zt) & band
    out["band_hz"] = [float(bh[0]), float(bh[1])]
    out["band_source"] = "rail band_hz" if rail.get("band_hz") else "default (fmin to target_fmax_hz)"
    # package and die (#22): Z at the die = (Z_board + R_pkg + jw L_pkg) || 1/(jw C_die)
    pk = rail.get("package")
    if pk:
        w = 2 * np.pi * freqs
        zpk = Z + pk.get("r_mohm", 0.0) * 1e-3 + 1j * w * pk.get("l_nh", 0.0) * 1e-9
        cd = pk.get("c_die_nf", 0.0) * 1e-9
        zdie = zpk / (1 + 1j * w * cd * zpk) if cd > 0 else zpk
        out["package"] = {**pk, "Z_die_ohm": np.abs(zdie).tolist(), "max_Z_die_ohm": round(float(np.abs(zdie)[band].max()), 5)}
    # load step (#22): droop of the rail voltage at the loads for a step of step_a with rise step_rise_ns
    try:
        out.update(load_step(freqs, Z, rail.get("step_a", rail["current_a"] / 2), rail.get("step_rise_ns", 1.0)))
    except Exception as e:
        out["load_step_error"] = str(e)
    out.update({"target_ohm": round(zt, 5), "target_fmax_hz": fmax_t, "caps": len(caps_on), "shunt_ohm": shunt_r,
                "freqs_hz": freqs.tolist(), "Z_ohm": mag.tolist(), "Z21_ohm": np.abs(Z21).tolist(),
                "max_Z_ohm": round(float(mag[band].max()), 5), "f_at_max_hz": float(freqs[band][int(np.argmax(mag[band]))]),
                "max_ratio": round(float((mag[band] / zt).max()), 3),
                "violation_band_hz": [float(freqs[viol].min()), float(freqs[viol].max())] if viol.any() else None,
                "cap_list": [{"ref": c["ref"], "value": c["value"], "ESL_nH": round(c["ESL"] * 1e9, 3)} for c in caps_on]})
    # G12: resonances = local maxima of |Z| at the load; location = node with the largest response
    res = []
    for i in range(1, len(freqs) - 1):
        if mag[i] > mag[i - 1] and mag[i] > mag[i + 1] and mag[i] > 0.3 * zt:
            vmax = int(np.argmax(modes[i]))
            node = net.nodes[keep[vmax]]
            res.append({"f_hz": float(freqs[i]), "Z_ohm": round(float(mag[i]), 5),
                        "x": round(node[2], 2), "y": round(node[3], 2), "layer": ctx.b.layers[node[1]]})
    out["resonances"] = res
    if optimise:
        out["optimiser"] = optimise_caps(caps_on, zt, lambda cs: np.abs(port_z(cs)[0])[band], with_added,
                                         cfg["candidates"])
    return out


def optimise_caps(caps_on, zt, zfun, with_added, candidates, max_add=12):
    """Greedy: add the candidate that lowers the worst |Z|/Zt most, until it meets the target;
    then try removing existing capacitors one at a time without breaking the target."""
    added = defaultdict(int)
    base = zfun(caps_on)
    cur = base
    steps = []
    for _ in range(max_add):
        if (cur / zt).max() <= 1.0:
            break
        best = None
        for cand in candidates:
            trial = dict(added)
            trial[cand] = trial.get(cand, 0) + 1
            z = zfun(caps_on + with_added(trial))
            score = (z / zt).max()
            if best is None or score < best[0]:
                best = (score, cand, z)
        if best is None or best[0] >= (cur / zt).max() - 1e-3:
            break
        added[best[1]] += 1
        cur = best[2]
        steps.append({"add": best[1], "max_ratio": round(float(best[0]), 3)})
    removable = []
    if (cur / zt).max() <= 1.0:
        kept = list(caps_on)
        for c in sorted(caps_on, key=lambda c: c["C"]):
            trial = [x for x in kept if x is not c]
            z = zfun(trial + with_added(added))
            if (z / zt).max() <= 1.0:
                kept = trial
                removable.append(c["ref"])
    return {"start_ratio": round(float((base / zt).max()), 3), "end_ratio": round(float((cur / zt).max()), 3),
            "add": dict(added), "steps": steps, "removable": removable,
            "met": bool((cur / zt).max() <= 1.0)}


# ==============================================================================================
# validation (G17) for the network solvers

def validation_cases():
    """Closed-form cases: rectangular sheet resistance and the series resonance of a capacitor."""
    from .validate import _case
    out = []
    # a bare rectangular sheet of copper L x W at thickness t: R = L / (sigma W t)
    sigma, t, L, W, cell = 5.8e7, 35e-6, 10.0, 2.0, 0.1
    nx, ny = int(L / cell), int(W / cell)
    N = nx * ny
    idx = np.arange(N).reshape(ny, nx)
    Gs = sigma * t
    rows, cols, vals = [], [], []
    for a, c in ((idx[:, :-1], idx[:, 1:]), (idx[:-1, :], idx[1:, :])):
        rows += [a.ravel(), c.ravel()]
        cols += [c.ravel(), a.ravel()]
        vals += [np.full(a.size, -Gs)] * 2
    Y = sp.coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(N, N)).tocsr()
    Y = Y + sp.diags(-np.asarray(Y.sum(axis=1)).ravel())
    # left column at 1 V, right column at 0 V
    fixed = np.r_[idx[:, 0], idx[:, -1]]
    Vf = np.r_[np.ones(ny), np.zeros(ny)]
    free = np.setdiff1d(np.arange(N), fixed)
    V = np.zeros(N)
    V[fixed] = Vf
    V[free] = spla.spsolve(Y[free][:, free].tocsc(), -Y[free][:, fixed] @ Vf)
    I = (Y @ V)[idx[:, 0]].sum()
    R_num = 1.0 / I
    R_ref = ((nx - 1) * cell * 1e-3) / (sigma * W * 1e-3 * t)
    out.append(_case("DC: 10 x 2 mm sheet, 35 um copper, resistance", "dc", R_num, R_ref, 0.5, "ohm",
                     "R = L / (sigma W t), exact"))
    # capacitor model: |Z| at self-resonance equals ESR, f_srf = 1 / (2 pi sqrt(L C))
    cap = {"C": 1e-7, "ESR": 0.015, "ESL": 0.6e-9}
    fs = 1 / (2 * math.pi * math.sqrt(cap["ESL"] * cap["C"]))
    fr = np.logspace(math.log10(fs) - 0.5, math.log10(fs) + 0.5, 20001)
    zz = np.abs([zcap(cap, 2 * math.pi * f) for f in fr])
    out.append(_case("capacitor model self-resonance frequency", "pdn", float(fr[int(np.argmin(zz))]), fs, 0.1, "Hz",
                     "1 / (2 pi sqrt(LC)), exact"))
    # parallel-plate cavity: first resonance of a rectangular plane pair a x b (open edges): f10 = c / (2 a sqrt(er))
    out += cavity_case()
    return out


def cavity_case():
    """Lumped-mesh cavity of a 40 x 20 mm plane pair (h = 0.1 mm, er 4.3): first mode f = c / (2 a sqrt(er))."""
    from .validate import _case
    a, bb, h, er, cell = 40e-3, 20e-3, 0.1e-3, 4.3, 0.5e-3
    nx, ny = int(a / cell), int(bb / cell)
    N = nx * ny
    idx = np.arange(N).reshape(ny, nx)
    Lsq = MU0 * h
    Cc = EPS0 * er * cell * cell / h
    rows, cols = [], []
    for p, q in ((idx[:, :-1], idx[:, 1:]), (idx[:-1, :], idx[1:, :])):
        rows += [p.ravel(), q.ravel()]
        cols += [q.ravel(), p.ravel()]
    rows, cols = np.concatenate(rows), np.concatenate(cols)
    Gm = sp.coo_matrix((-np.ones(len(rows)), (rows, cols)), shape=(N, N)).tocsr()
    Gm = Gm + sp.diags(-np.asarray(Gm.sum(axis=1)).ravel())
    # generalized eigenproblem (1/L) G v = w^2 C v
    vals = spla.eigsh(Gm / Lsq, k=4, M=sp.identity(N) * Cc, sigma=-1e15, which="LM")[0]
    fs = np.sort(np.sqrt(np.abs(vals)) / (2 * math.pi))
    f10 = fs[fs > 1e6][0]
    ref = C0 / (2 * a * math.sqrt(er))
    return [_case("cavity mesh 40 x 20 mm: first resonance", "pdn", float(f10), ref, 2.0, "Hz",
                  "rectangular cavity, open edges: f10 = c / (2 a sqrt(er))")]
