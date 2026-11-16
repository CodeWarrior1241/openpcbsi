"""Board model built from a snapshot: stackup geometry, numpy arrays, per-layer copper rasters."""
import fnmatch
import math
import os
import re

import numpy as np

from . import geom

DEFAULT_MASK_ER = 3.8
C0 = 299792458.0


class Stackup:
    """Copper layers top to bottom with z positions (mm from the board top surface, mask excluded)."""

    def __init__(self, snap_stackup, copper_names, board_thickness):
        self.copper = []      # [{name, t, z_top, z_bot}]
        self.dielectrics = []  # [{name, t, er, tand, z_top, z_bot, above, below, material}]
        self.mask = {}         # side -> {t, er}
        z = 0.0
        last_cu = None
        items = [s for s in snap_stackup if s.get("type") not in (None, "finish")]
        if not any(s.get("type") == "copper" for s in items):
            # no stackup in the file: assume equal dielectrics, 35 um outer, 18 um inner copper
            n = len(copper_names)
            items = []
            for k, nm in enumerate(copper_names):
                items.append({"name": nm, "type": "copper", "thickness": 0.035 if k in (0, n - 1) else 0.018})
                if k < n - 1:
                    items.append({"name": f"dielectric {k+1}", "type": "core", "thickness": board_thickness / (n - 1),
                                  "epsilon_r": 4.3, "loss_tangent": 0.02})
        for s in items:
            ty = (s.get("type") or "").lower()
            t = s.get("thickness") or 0.0
            if "mask" in ty:
                side = "F" if s["name"].startswith("F.") else "B"
                self.mask[side] = {"t": t, "er": s.get("epsilon_r") or DEFAULT_MASK_ER}
                continue
            if ty == "copper":
                self.copper.append({"name": s["name"], "t": t, "z_top": z, "z_bot": z + t})
                last_cu = s["name"]
                z += t
            elif ty in ("core", "prepreg", "dielectric") or s["name"].startswith("dielectric"):
                self.dielectrics.append({"name": s["name"], "t": t, "er": s.get("epsilon_r") or 4.3,
                                         "tand": s.get("loss_tangent") or 0.02, "z_top": z, "z_bot": z + t,
                                         "above": last_cu, "material": s.get("material")})
                z += t
        for d in self.dielectrics:
            below = [c for c in self.copper if c["z_top"] >= d["z_bot"] - 1e-9]
            d["below"] = below[0]["name"] if below else None
        self.thickness = z
        self.names = [c["name"] for c in self.copper]
        self.index = {n: i for i, n in enumerate(self.names)}

    def cu(self, name):
        return self.copper[self.index[name]]

    def gap(self, a, b):
        """Dielectric distance between the facing surfaces of copper layers a and b (mm)."""
        ca, cb = self.cu(a), self.cu(b)
        if ca["z_top"] > cb["z_top"]:
            ca, cb = cb, ca
        return cb["z_top"] - ca["z_bot"]

    def er_between(self, a, b):
        """Thickness-weighted permittivity of the dielectric between two copper layers."""
        ca, cb = self.cu(a), self.cu(b)
        z0, z1 = sorted((ca["z_bot"], cb["z_top"]))
        tot, acc = 0.0, 0.0
        for d in self.dielectrics:
            ov = min(z1, d["z_bot"]) - max(z0, d["z_top"])
            if ov > 0:
                tot += ov
                acc += ov * d["er"]
        return acc / tot if tot else 4.3

    def tand_between(self, a, b):
        ca, cb = self.cu(a), self.cu(b)
        z0, z1 = sorted((ca["z_bot"], cb["z_top"]))
        tot, acc = 0.0, 0.0
        for d in self.dielectrics:
            ov = min(z1, d["z_bot"]) - max(z0, d["z_top"])
            if ov > 0:
                tot += ov
                acc += ov * d["tand"]
        return acc / tot if tot else 0.02

    def neighbours(self, name):
        i = self.index[name]
        out = []
        if i > 0:
            out.append(self.names[i - 1])
        if i < len(self.names) - 1:
            out.append(self.names[i + 1])
        return out

    def is_outer(self, name):
        return self.index[name] in (0, len(self.names) - 1)

    def barrel(self, a, b):
        """Via barrel length between the centres of copper layers a and b."""
        ca, cb = self.cu(a), self.cu(b)
        return abs((ca["z_top"] + ca["z_bot"]) / 2 - (cb["z_top"] + cb["z_bot"]) / 2)

    def mean_er(self):
        tot = sum(d["t"] for d in self.dielectrics)
        return sum(d["t"] * d["er"] for d in self.dielectrics) / tot if tot else 4.3

    def describe(self):
        rows = []
        for c in self.copper:
            rows.append({"layer": c["name"], "t_um": round(c["t"] * 1000, 1), "z_mm": round(c["z_top"], 4)})
        return rows


class Board:
    def __init__(self, snap, config=None, cache_dir=None):
        self.snap = snap
        self.path = snap["board"]
        self.sha = snap["sha256"]
        self.name = os.path.splitext(os.path.basename(self.path))[0]
        self.cache_dir = cache_dir
        self.config = config
        self.stack = Stackup(snap["stackup"], snap["copper_layers"], snap["board_thickness"])
        self.layers = self.stack.names
        self.L = self.stack.index
        self.net_names = [""] + sorted(snap["nets"].keys())
        self.net_id = {n: i for i, n in enumerate(self.net_names)}
        self.kicad_class = {n: v["netclass"] for n, v in snap["nets"].items()}

        T = snap["tracks"]
        self.tracks = {
            "x1": np.array([t["x1"] for t in T]), "y1": np.array([t["y1"] for t in T]),
            "x2": np.array([t["x2"] for t in T]), "y2": np.array([t["y2"] for t in T]),
            "w": np.array([t["w"] for t in T]),
            "layer": np.array([self.L.get(t["layer"], -1) for t in T], dtype=np.int64),
            "net": np.array([self.net_id.get(t["net"], 0) for t in T], dtype=np.int64),
        }
        self.tracks["len"] = np.hypot(self.tracks["x2"] - self.tracks["x1"], self.tracks["y2"] - self.tracks["y1"])
        V = snap["vias"]
        self.vias = {
            "x": np.array([v["x"] for v in V]), "y": np.array([v["y"] for v in V]),
            "d": np.array([v["d"] for v in V]), "drill": np.array([v["drill"] for v in V]),
            "net": np.array([self.net_id.get(v["net"], 0) for v in V], dtype=np.int64),
            "top": np.array([self.L.get(v["span"][0], 0) for v in V], dtype=np.int64),
            "bot": np.array([self.L.get(v["span"][1], len(self.layers) - 1) for v in V], dtype=np.int64),
        }
        self.via_flashed = [set(self.L[l] for l in v["flashed"]) for v in V]
        self.pads = snap["pads"]
        for p in self.pads:
            p["net_id"] = self.net_id.get(p["net"], 0)
            p["layer_ids"] = [self.L[l] for l in p["layers"] if l in self.L]
        self.fps = {f["ref"]: f for f in snap["footprints"]}
        self.pads_by_ref = {}
        for p in self.pads:
            self.pads_by_ref.setdefault(p["ref"], []).append(p)
        self.zones = [z for z in snap["zones"] if not z.get("rule_area")]
        ol = snap["outline"]
        pts = np.vstack([np.asarray(r) for poly in ol for r in poly]) if ol else np.zeros((0, 2))
        if len(pts):
            self.bbox = (pts[:, 0].min(), pts[:, 1].min(), pts[:, 0].max(), pts[:, 1].max())
        else:
            allx = np.r_[self.tracks["x1"], self.vias["x"]]
            ally = np.r_[self.tracks["y1"], self.vias["y"]]
            self.bbox = (allx.min() - 1, ally.min() - 1, allx.max() + 1, ally.max() + 1)
        self._rasters = {}
        self._grid = None

    # -------------------------------------------------------------------------------------------
    def nets_matching(self, patterns):
        out = []
        for n in self.net_names[1:]:
            for p in patterns:
                if fnmatch.fnmatchcase(n, p) or (p.startswith("re:") and re.search(p[3:], n)):
                    out.append(n)
                    break
        return out

    def grid(self, res=None):
        res = res or (self.config.raster_res if self.config else 0.025)
        if self._grid is None or abs(self._grid.res - res) > 1e-12:
            x0, y0, x1, y1 = self.bbox
            self._grid = geom.Grid(x0 - 0.5, y0 - 0.5, x1 + 0.5, y1 + 0.5, res)
        return self._grid

    def raster(self, layer):
        """int16 raster of net ids of all copper on a layer (zones, pads, tracks, flashed vias)."""
        li = self.L[layer] if isinstance(layer, str) else layer
        if li in self._rasters:
            return self._rasters[li]
        g = self.grid()
        key = f"raster_{self.sha[:16]}_{g.res:.4f}_{li}.npy"
        cpath = os.path.join(self.cache_dir, key) if self.cache_dir else None
        if cpath and os.path.exists(cpath):
            r = np.load(cpath)
            if r.shape == g.shape:
                self._rasters[li] = r
                return r
        lname = self.layers[li]
        r = np.zeros(g.shape, dtype=np.int16)
        zs = sorted([z for z in self.zones if z["layer"] == lname], key=lambda z: z.get("priority", 0))
        for z in zs:
            nid = self.net_id.get(z["net"], 0)
            rings = [ring for poly in z["polys"] for ring in poly]
            m = geom.raster_rings(g, rings)
            r[m] = nid
        for p in self.pads:
            if li in p["layer_ids"] and p["poly"]:
                rings = [ring for poly in p["poly"] for ring in poly]
                m = geom.raster_rings(g, rings)
                r[m] = p["net_id"]
        T = self.tracks
        sel = np.flatnonzero(T["layer"] == li)
        for k in sel:
            geom.raster_capsule(g, r, T["x1"][k], T["y1"][k], T["x2"][k], T["y2"][k], T["w"][k] / 2, T["net"][k])
        V = self.vias
        for k in range(len(V["x"])):
            if li in self.via_flashed[k]:
                geom.raster_capsule(g, r, V["x"][k], V["y"][k], V["x"][k], V["y"][k], V["d"][k] / 2, V["net"][k])
        self._rasters[li] = r
        if cpath:                      # atomic: parallel runs never read a half-written cache (#33)
            tmp = f"{cpath}.{os.getpid()}.tmp.npy"
            np.save(tmp, r)
            os.replace(tmp, cpath)
        return r

    def outline_rings(self):
        return [ring for poly in self.snap["outline"] for ring in poly]

    # -------------------------------------------------------------------------------------------
    def tracks_of(self, net, layer=None):
        nid = self.net_id[net] if isinstance(net, str) else net
        m = self.tracks["net"] == nid
        if layer is not None:
            m &= self.tracks["layer"] == (self.L[layer] if isinstance(layer, str) else layer)
        return np.flatnonzero(m)

    def vias_of(self, net):
        nid = self.net_id[net] if isinstance(net, str) else net
        return np.flatnonzero(self.vias["net"] == nid)

    def pads_of(self, net):
        nid = self.net_id[net] if isinstance(net, str) else net
        return [p for p in self.pads if p["net_id"] == nid]

    def via_attached_layers(self, k, tol=0.01):
        """Copper layers on which via k is connected: a track end on it, a pad on it, or flashed in a zone."""
        x, y, nid = self.vias["x"][k], self.vias["y"][k], self.vias["net"][k]
        T = self.tracks
        m = T["net"] == nid
        r = self.vias["d"][k] / 2 + tol
        near = m & ((np.hypot(T["x1"] - x, T["y1"] - y) <= r) | (np.hypot(T["x2"] - x, T["y2"] - y) <= r))
        layers = set(T["layer"][near].tolist())
        for p in self.pads:
            if p["net_id"] == nid and abs(p["x"] - x) < max(p["sx"], p["sy"]) and abs(p["y"] - y) < max(p["sx"], p["sy"]):
                if p["poly"] and geom.points_in_rings([x], [y], [ring for poly in p["poly"] for ring in poly])[0]:
                    layers.update(p["layer_ids"])
        return layers

    def zone_connected_layers(self, k):
        """Layers on which a zone of the via's own net covers the via centre."""
        g = self.grid()
        out = set()
        nid = self.vias["net"][k]
        for li, lname in enumerate(self.layers):
            if not any(z["layer"] == lname and self.net_id.get(z["net"]) == nid for z in self.zones):
                continue
            r = self.raster(li)
            if g.sample(r, [self.vias["x"][k]], [self.vias["y"][k]])[0] == nid:
                # copper of the net around the via, beyond its own pad: test a ring just outside the pad
                ang = np.linspace(0, 2 * np.pi, 8, endpoint=False)
                rr = self.vias["d"][k] / 2 + 2 * g.res
                s = g.sample(r, self.vias["x"][k] + rr * np.cos(ang), self.vias["y"][k] + rr * np.sin(ang))
                if (s == nid).any():
                    out.add(li)
        return out
