"""Shared state for checks: exclusion areas, reference masks, pairs, ordered track paths."""
import math
from collections import defaultdict

import numpy as np
from scipy import ndimage

from . import geom


class Context:
    def __init__(self, board, config):
        self.b = board
        self.cfg = config
        self.g = board.grid()
        self.ref_ids = np.array(sorted(board.net_id[n] for n in config.reference_nets if n in board.net_id),
                                dtype=np.int16)
        self._refmask = {}
        self._anyplane = {}
        self._areas = None
        self._area_names = []
        names = set(board.net_names[1:])
        self.partner = {}
        for n in board.net_names[1:]:
            p = config.partner(n, names)
            if p:
                self.partner[n] = p
                self.partner[p] = n
        self.zone_nets = set(z["net"] for z in board.zones)

    # ------------------------------------------------------------------------------------------
    def refmask(self, li):
        if li not in self._refmask:
            r = self.b.raster(li)
            m = np.isin(r, self.ref_ids)
            if self.cfg["board"]["power_as_reference"]:
                pw = [self.b.net_id[n] for n in self.zone_nets if n in self.b.net_id]
                m |= np.isin(r, pw)
            self._refmask[li] = m
        return self._refmask[li]

    def anyplane(self, li):
        """Reference nets or any net that has a zone (a plane) on this layer."""
        if li not in self._anyplane:
            r = self.b.raster(li)
            lname = self.b.layers[li]
            ids = set(self.ref_ids.tolist())
            ids |= {self.b.net_id[z["net"]] for z in self.b.zones if z["layer"] == lname and z["net"] in self.b.net_id}
            self._anyplane[li] = np.isin(r, np.array(sorted(ids), dtype=np.int16))
        return self._anyplane[li]

    def ref_layers(self, li):
        lname = self.b.layers[li]
        over = self.cfg["reference_layers"].get(lname)
        if over:
            return [self.b.L[x] for x in over]
        out = []
        for n in self.b.stack.neighbours(lname):
            nli = self.b.L[n]
            # a neighbour counts as a reference layer if it carries reference copper at all
            if self.refmask(nli).mean() > 0.02:
                out.append(nli)
        return out or [self.b.L[n] for n in self.b.stack.neighbours(lname)]

    # ------------------------------------------------------------------------------------------
    def areas(self):
        """int8 raster: 0 outside any exclusion area, k for area k (names in self.area_names)."""
        if self._areas is not None:
            return self._areas
        a = np.zeros(self.g.shape, dtype=np.int8)
        self._area_names = [""]
        for k, (name, spec) in enumerate(self.cfg["areas"].items(), start=1):
            m = np.zeros(self.g.shape, dtype=bool)
            for ref in spec.get("footprints", []):
                f = self.b.fps.get(ref)
                if not f:
                    continue
                rings = [ring for poly in f["courtyard"] for ring in poly]
                if rings:
                    geom.raster_rings(self.g, rings, m)
                else:
                    x0, y0, x1, y1 = f["bbox"]
                    geom.raster_rings(self.g, [[(x0, y0), (x1, y0), (x1, y1), (x0, y1)]], m)
            for rect in spec.get("rects", []):
                x0, y0, x1, y1 = rect
                geom.raster_rings(self.g, [[(x0, y0), (x1, y0), (x1, y1), (x0, y1)]], m)
            margin = spec.get("margin", 0.0)
            if margin > 0 and m.any():
                d = ndimage.distance_transform_edt(~m) * self.g.res
                m = d <= margin
            a[m & (a == 0)] = k
            self._area_names.append(name)
        self._areas = a
        return a

    @property
    def area_names(self):
        self.areas()
        return self._area_names

    def area_at(self, x, y):
        a = self.g.sample(self.areas(), x, y)
        return a

    def area_name(self, code):
        return self.area_names[int(code)] if code else ""

    # ------------------------------------------------------------------------------------------
    def nets(self, rule):
        return self.cfg.nets_in_rule(self.b, rule)

    def net_class(self, net):
        return self.cfg.class_of(net)

    def cparam(self, net, key, default=None):
        c = self.cfg.class_of(net)
        return self.cfg.class_param(c, key, default) if c else default

    # ------------------------------------------------------------------------------------------
    def layer_paths(self, net, li, tol=0.005):
        """Ordered chains of the net's tracks on one layer: [[(k, flipped), ...], ...]."""
        b = self.b
        ks = b.tracks_of(net, li)
        if len(ks) == 0:
            return []
        T = b.tracks
        lab = geom.cluster(np.r_[T["x1"][ks], T["x2"][ks]], np.r_[T["y1"][ks], T["y2"][ks]], 0.005)
        n = len(ks)
        cid = {}
        for i, k in enumerate(ks):
            cid[(k, 0)] = lab[i]
            cid[(k, 1)] = lab[n + i]
        adj = defaultdict(list)
        for k in ks:
            adj[cid[(k, 0)]].append((k, 0))
            adj[cid[(k, 1)]].append((k, 1))
        used = set()
        chains = []
        starts = [n for n, v in adj.items() if len(v) != 2] + list(adj.keys())
        for s in starts:
            for (k0, end0) in adj[s]:
                if k0 in used:
                    continue
                chain = []
                node, k, end = s, k0, end0
                while k is not None and k not in used:
                    used.add(k)
                    flipped = end == 1          # we enter at end -> walk from x2 to x1
                    chain.append((k, flipped))
                    node = cid[(k, 0)] if flipped else cid[(k, 1)]
                    nxt = [(kk, ee) for kk, ee in adj[node] if kk not in used]
                    if len(adj[node]) == 2 and nxt:
                        k, end = nxt[0]
                    else:
                        k = None
                chains.append(chain)
        return chains

    def path_samples(self, net, li, step):
        """Samples along each ordered chain: list of dicts with x, y, s (arc length), w, k, dl."""
        T = self.b.tracks
        out = []
        for chain in self.layer_paths(net, li):
            xs, ys, ss, ws, kk, dls = [], [], [], [], [], []
            s0 = 0.0
            for k, fl in chain:
                x1, y1, x2, y2 = T["x1"][k], T["y1"][k], T["x2"][k], T["y2"][k]
                if fl:
                    x1, y1, x2, y2 = x2, y2, x1, y1
                px, py, dl = geom.sample_segment(x1, y1, x2, y2, step)
                n = len(px)
                xs.append(px); ys.append(py)
                ss.append(s0 + (np.arange(n) + 0.5) * dl)
                ws.append(np.full(n, T["w"][k])); kk.append(np.full(n, k)); dls.append(np.full(n, dl))
                s0 += T["len"][k]
            out.append({"x": np.concatenate(xs), "y": np.concatenate(ys), "s": np.concatenate(ss),
                        "w": np.concatenate(ws), "k": np.concatenate(kk), "dl": np.concatenate(dls),
                        "length": s0})
        return out

    # ------------------------------------------------------------------------------------------
    def via_zone_layers(self):
        """For every via: set of layers on which its own net's copper surrounds it (zone connection)."""
        if hasattr(self, "_vzl"):
            return self._vzl
        b, g = self.b, self.g
        V = b.vias
        n = len(V["x"])
        out = [set() for _ in range(n)]
        ang = np.linspace(0, 2 * np.pi, 8, endpoint=False)
        for li in range(len(b.layers)):
            r = b.raster(li)
            rr = V["d"] / 2 + 2 * g.res
            px = V["x"][:, None] + rr[:, None] * np.cos(ang)[None, :]
            py = V["y"][:, None] + rr[:, None] * np.sin(ang)[None, :]
            s = g.sample(r, px, py)
            hit = (s == V["net"][:, None]).sum(axis=1) >= 2
            for k in np.flatnonzero(hit):
                out[k].add(li)
        self._vzl = out
        return out

    def via_track_layers(self):
        """For every via: set of layers with a track of its net ending on it, or a pad of its net under it."""
        if hasattr(self, "_vtl"):
            return self._vtl
        b = self.b
        V, T = b.vias, b.tracks
        n = len(V["x"])
        out = [set() for _ in range(n)]
        from scipy.spatial import cKDTree
        if n:
            tree = cKDTree(np.c_[V["x"], V["y"]])
            for (xa, ya) in (("x1", "y1"), ("x2", "y2")):
                pts = np.c_[T[xa], T[ya]]
                res = tree.query_ball_point(pts, r=float(V["d"].max() / 2 + 0.01))
                for t, cand in enumerate(res):
                    for k in cand:
                        if V["net"][k] == T["net"][t] and math.hypot(V["x"][k] - pts[t, 0], V["y"][k] - pts[t, 1]) <= V["d"][k] / 2 + 0.01:
                            out[k].add(int(T["layer"][t]))
            for p in b.pads:
                if p["kind"] not in ("smd", "conn") or not p["poly"]:
                    continue
                cand = tree.query_ball_point([p["x"], p["y"]], r=max(p["sx"], p["sy"]))
                cand = [k for k in cand if V["net"][k] == p["net_id"]]
                if not cand:
                    continue
                rings = [ring for poly in p["poly"] for ring in poly]
                inside = geom.points_in_rings(V["x"][cand], V["y"][cand], rings)
                for k, ins in zip(cand, inside):
                    if ins:
                        out[k].update(p["layer_ids"])
        self._vtl = out
        return out


def runs(mask):
    """Start/stop index pairs of True runs in a 1-D bool array."""
    m = np.r_[False, np.asarray(mask, bool), False]
    d = np.diff(m.astype(np.int8))
    return list(zip(np.flatnonzero(d == 1), np.flatnonzero(d == -1)))
