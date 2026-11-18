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


MAX_RASTER_CELLS = 3e8       # cells per layer raster (int16: 0.6 GB); a 430 x 430 mm board at 25 um


def pad_poly(p, layer):
    """A pad's polygons on one copper layer (by name): its padstack's shape there, else its one shape."""
    return (p.get("layer_poly") or {}).get(layer, p["poly"])


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
        self.tracks["arc"] = np.array([bool(t.get("arc")) for t in T], dtype=bool)          # chords of an arc track
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
        for z in self.zones:
            z["net_id"] = self.net_id.get(z["net"], 0)
        # copper that is not a track, pad, via or zone: graphics and text on copper layers. Strokes (lines,
        # arc chords with their circle, text strokes, discs as zero-length strokes) as arrays like the
        # tracks'; filled shapes as polygons. Kept apart from the tracks: they are no part of a net's routing.
        G = [(g, sg, None) for g in snap.get("copper_gfx", []) for sg in g["segs"]] + \
            [(g, a[:5], a[5:8]) for g in snap.get("copper_gfx", []) for a in g["arcs"]]
        G = [x for x in G if x[0]["layer"] in self.L]
        self.gfx = {"x1": np.array([s_[0] for _, s_, _ in G], float), "y1": np.array([s_[1] for _, s_, _ in G], float),
                    "x2": np.array([s_[2] for _, s_, _ in G], float), "y2": np.array([s_[3] for _, s_, _ in G], float),
                    "w": np.array([s_[4] for _, s_, _ in G], float),
                    "layer": np.array([self.L[g["layer"]] for g, _, _ in G], dtype=np.int64),
                    "net": np.array([self.net_id.get(g["net"], 0) for g, _, _ in G], dtype=np.int64),
                    "arc": np.array([c if c else (np.nan, np.nan, np.nan) for _, _, c in G], float).reshape(-1, 3),
                    "exact": np.array([bool(g.get("exact", True)) for g, _, _ in G], bool),
                    "item": np.array([id(g) for g, _, _ in G], dtype=np.int64)}
        self.gfx_polys = [{"layer": self.L[g["layer"]], "net_id": self.net_id.get(g["net"], 0), "polys": g["polys"],
                           "exact": bool(g.get("exact", True)), "item": id(g), "src": g.get("src", "")}
                          for g in snap.get("copper_gfx", []) if g["polys"] and g["layer"] in self.L]
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
            cells = ((x1 - x0 + 1) / res) * ((y1 - y0 + 1) / res)
            cap = float(os.environ.get("OPENPCBSI_MAX_RASTER_CELLS", MAX_RASTER_CELLS))
            if cells > cap:
                from .snapshot import InputError
                raise InputError(f"the board's extent ({x1 - x0:.0f} x {y1 - y0:.0f} mm) at {res * 1000:g} um is {cells / 1e6:.0f} million "
                                 f"cells per layer, above the {cap / 1e6:.0f} million this run allows (about {2 * cap / 1e9:.1f} GB per "
                                 f"layer): check the board outline, or raise OPENPCBSI_MAX_RASTER_CELLS")
            self._grid = geom.Grid(x0 - 0.5, y0 - 0.5, x1 + 0.5, y1 + 0.5, res)
        return self._grid

    def raster(self, layer):
        """int16 raster of net ids of all copper on a layer (zones, pads, tracks, flashed vias)."""
        li = self.L[layer] if isinstance(layer, str) else layer
        if li in self._rasters:
            return self._rasters[li]
        g = self.grid()
        key = f"raster_{self.sha[:16]}_v{self.snap.get('version', 0)}_{g.res:.4f}_{li}{'_nn' if getattr(self, '_nonet', False) else ''}.npy"
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
            nid = z["net_id"]
            rings = [ring for poly in z["polys"] for ring in poly]
            m = geom.raster_rings(g, rings)
            r[m] = nid
        for gp in self.gfx_polys:
            if gp["layer"] == li:
                r[geom.raster_rings(g, [ring for poly in gp["polys"] for ring in poly])] = gp["net_id"]
        for p in self.pads:
            if li in p["layer_ids"] and p["poly"]:
                rings = [ring for poly in pad_poly(p, lname) for ring in poly]
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
        F = self.gfx
        for k in np.flatnonzero(F["layer"] == li):
            geom.raster_capsule(g, r, F["x1"][k], F["y1"][k], F["x2"][k], F["y2"][k], F["w"][k] / 2, F["net"][k])
        self._rasters[li] = r
        if cpath:                      # atomic: parallel runs never read a half-written cache (#33)
            tmp = f"{cpath}.{os.getpid()}.tmp.npy"
            np.save(tmp, r)
            os.replace(tmp, cpath)
        return r

    def outline_rings(self):
        return [ring for poly in self.snap["outline"] for ring in poly]

    def label_nonet_pads(self):
        """Give copper without a net a label of its own, so that it is copper for every check: pads
        (fiducials, mounting pads), tracks, vias, zones and graphics or text on copper layers. Items
        without a net that touch (on a layer, or through a pad or via) share one label, so they are not
        judged against each other. Labels are appended to net_names: '<no net> REF.PAD' for a group with
        a pad, else '<no net> <what> on <layer>'."""
        if getattr(self, "_nonet", False):
            return
        items = []                                         # (kind, index, layers, shapes per layer, name)
        for k, p in enumerate(self.pads):
            if p["net_id"] == 0 and p["poly"] and p["layer_ids"] and p["kind"] != "npth":
                items.append(("pad", k, p["layer_ids"], f"{p['ref']}.{p['num']}"))
        T = self.tracks
        for k in np.flatnonzero(T["net"] == 0):
            if T["layer"][k] >= 0:
                items.append(("track", int(k), [int(T["layer"][k])],
                              f"track on {self.layers[T['layer'][k]]} at ({T['x1'][k]:.2f}, {T['y1'][k]:.2f})"))
        V = self.vias
        for k in np.flatnonzero(V["net"] == 0):
            items.append(("via", int(k), sorted(self.via_flashed[k]), f"via at ({V['x'][k]:.2f}, {V['y'][k]:.2f})"))
        for k, z in enumerate(self.zones):
            if z["net_id"] == 0 and z["layer"] in self.L:
                items.append(("zone", k, [self.L[z["layer"]]], f"zone on {z['layer']}"))
        F = self.gfx
        by_item = {}
        for k in np.flatnonzero(F["net"] == 0):
            by_item.setdefault(int(F["item"][k]), []).append(int(k))
        src = {gp["item"]: gp["src"] for gp in self.gfx_polys}
        src.update({id(g): g.get("src", "") for g in self.snap.get("copper_gfx", [])})
        for it, ks in by_item.items():
            items.append(("gfx", it, [int(F["layer"][ks[0]])], f"{src.get(it, 'graphic')} on {self.layers[F['layer'][ks[0]]]} "
                                                              f"at ({F['x1'][ks[0]]:.2f}, {F['y1'][ks[0]]:.2f})"))
        for k, gp in enumerate(self.gfx_polys):
            if gp["net_id"] == 0 and gp["item"] not in by_item:
                items.append(("gpoly", k, [gp["layer"]], f"{gp['src'] or 'graphic'} on {self.layers[gp['layer']]}"))
        parent = list(range(len(items)))

        def find(a):
            while parent[a] != a:
                parent[a] = parent[parent[a]]
                a = parent[a]
            return a
        if len(items) > 1:
            try:
                from .parse import use
                use()
                from shapely import geometry, ops, strtree
                shapes = {}

                def geo(kind, k, li):
                    if kind == "pad":
                        polys = pad_poly(self.pads[k], self.layers[li])
                        return ops.unary_union([geometry.Polygon(pl[0], pl[1:]).buffer(0) for pl in polys if len(pl[0]) >= 3])
                    if kind == "track":
                        return geometry.LineString([(T["x1"][k], T["y1"][k]), (T["x2"][k], T["y2"][k])]).buffer(T["w"][k] / 2)
                    if kind == "via":
                        return geometry.Point(V["x"][k], V["y"][k]).buffer(V["d"][k] / 2)
                    if kind == "zone":
                        return ops.unary_union([geometry.Polygon(pl[0], pl[1:]).buffer(0) for pl in self.zones[k]["polys"] if len(pl[0]) >= 3])
                    if kind == "gfx":
                        parts = [geometry.LineString([(F["x1"][j], F["y1"][j]), (F["x2"][j], F["y2"][j])]).buffer(F["w"][j] / 2)
                                 if (F["x1"][j], F["y1"][j]) != (F["x2"][j], F["y2"][j]) else
                                 geometry.Point(F["x1"][j], F["y1"][j]).buffer(F["w"][j] / 2) for j in by_item[k]]
                        parts += [geometry.Polygon(pl[0], pl[1:]).buffer(0) for gp in self.gfx_polys if gp["item"] == k for pl in gp["polys"]]
                        return ops.unary_union(parts)
                    return ops.unary_union([geometry.Polygon(pl[0], pl[1:]).buffer(0) for pl in self.gfx_polys[k]["polys"] if len(pl[0]) >= 3])
                for li in range(len(self.layers)):
                    idx = [i for i, it in enumerate(items) if li in it[2]]
                    if len(idx) < 2:
                        continue
                    gs = [geo(items[i][0], items[i][1], li) for i in idx]
                    tree = strtree.STRtree(gs)
                    for a_, ga in enumerate(gs):
                        for b_ in tree.query(ga):
                            b_ = int(b_)
                            if b_ > a_ and ga.intersects(gs[b_]):
                                parent[find(idx[a_])] = find(idx[b_])
            except ImportError:
                pass                                       # no geometry library: every item is its own label
        groups = {}
        for i in range(len(items)):
            groups.setdefault(find(i), []).append(items[i])
        for members in groups.values():
            pads_ = sorted(m[3] for m in members if m[0] == "pad")
            name = f"<no net> {pads_[0] if pads_ else members[0][3]}"
            while name in self.net_id:
                name += "'"
            nid = len(self.net_names)
            self.net_id[name] = nid
            self.net_names.append(name)
            for kind, k, _, _ in members:
                if kind == "pad":
                    self.pads[k]["net_id"] = nid
                elif kind == "track":
                    T["net"][k] = nid
                elif kind == "via":
                    V["net"][k] = nid
                elif kind == "zone":
                    self.zones[k]["net_id"] = nid
                elif kind == "gfx":
                    F["net"][np.array(by_item[k])] = nid
                    for gp in self.gfx_polys:
                        if gp["item"] == k:
                            gp["net_id"] = nid
                else:
                    self.gfx_polys[k]["net_id"] = nid
        self._nonet = True
        self._rasters = {}
        self._win = None

    def _window_index(self):
        """Per layer: zone fills as edge arrays with a bounding box, and the bounding boxes of the pads."""
        if getattr(self, "_win", None) is None:
            zs = {li: [] for li in range(len(self.layers))}
            for z in sorted(self.zones, key=lambda z: z.get("priority", 0)):
                li = self.L.get(z["layer"])
                segs = geom.ring_segments([ring for poly in z["polys"] for ring in poly])
                if li is None or not len(segs):
                    continue
                bb = (segs[:, [0, 2]].min(), segs[:, [1, 3]].min(), segs[:, [0, 2]].max(), segs[:, [1, 3]].max())
                zs[li].append((z["net_id"], segs, bb))
            for gp in self.gfx_polys:
                segs = geom.ring_segments([ring for poly in gp["polys"] for ring in poly])
                if len(segs):
                    bb = (segs[:, [0, 2]].min(), segs[:, [1, 3]].min(), segs[:, [0, 2]].max(), segs[:, [1, 3]].max())
                    zs[gp["layer"]].append((gp["net_id"], segs, bb))
            pb = np.full((len(self.pads), 4), np.nan)
            for k, p in enumerate(self.pads):
                if p["poly"]:
                    a = np.vstack([np.asarray(r, dtype=float) for pl in [p["poly"]] + list((p.get("layer_poly") or {}).values())
                                   for poly in pl for r in poly])
                    pb[k] = (a[:, 0].min(), a[:, 1].min(), a[:, 0].max(), a[:, 1].max())
            on = np.zeros((len(self.layers), len(self.pads)), bool)
            for k, p in enumerate(self.pads):
                on[p["layer_ids"], k] = True
            self._win = {"zones": zs, "pad_bb": pb, "pad_on": on}
        return self._win

    def raster_window(self, layer, x0, y0, x1, y1, res, pads=False):
        """Net-id raster (int32) of the copper of one layer inside a window, at any resolution: the same
        drawing order as raster(), from the vector data. Returns (grid, raster), and with `pads` also a
        mask of the pixels that are pad copper."""
        li = self.L[layer] if isinstance(layer, str) else layer
        w = self._window_index()
        g = geom.window_grid(x0, y0, x1, y1, res)
        r = np.zeros(g.shape, dtype=np.int32)
        for nid, segs, bb in w["zones"][li]:
            if bb[0] > x1 or bb[2] < x0 or bb[1] > y1 or bb[3] < y0:
                continue
            r[geom.raster_segs(g, segs)] = nid
        pb = w["pad_bb"]
        sel = np.flatnonzero(w["pad_on"][li] & (pb[:, 0] <= x1) & (pb[:, 2] >= x0) & (pb[:, 1] <= y1) & (pb[:, 3] >= y0))
        pm = np.zeros(g.shape, dtype=bool) if pads else None
        for k in sel:
            p = self.pads[k]
            m = geom.raster_rings(g, [ring for poly in pad_poly(p, self.layers[li]) for ring in poly])
            r[m] = p["net_id"]
            if pads:
                pm |= m
        T = self.tracks
        hw = T["w"] / 2
        sel = np.flatnonzero((T["layer"] == li) & (np.minimum(T["x1"], T["x2"]) - hw <= x1) & (np.maximum(T["x1"], T["x2"]) + hw >= x0)
                             & (np.minimum(T["y1"], T["y2"]) - hw <= y1) & (np.maximum(T["y1"], T["y2"]) + hw >= y0))
        for k in sel:
            geom.raster_capsule(g, r, T["x1"][k], T["y1"][k], T["x2"][k], T["y2"][k], hw[k], T["net"][k])
        V = self.vias
        vr = V["d"] / 2
        sel = np.flatnonzero((V["x"] - vr <= x1) & (V["x"] + vr >= x0) & (V["y"] - vr <= y1) & (V["y"] + vr >= y0))
        for k in sel:
            if li in self.via_flashed[k]:
                geom.raster_capsule(g, r, V["x"][k], V["y"][k], V["x"][k], V["y"][k], vr[k], V["net"][k])
        F = self.gfx
        hw = F["w"] / 2
        for k in np.flatnonzero((F["layer"] == li) & (np.minimum(F["x1"], F["x2"]) - hw <= x1) & (np.maximum(F["x1"], F["x2"]) + hw >= x0)
                                & (np.minimum(F["y1"], F["y2"]) - hw <= y1) & (np.maximum(F["y1"], F["y2"]) + hw >= y0)):
            geom.raster_capsule(g, r, F["x1"][k], F["y1"][k], F["x2"][k], F["y2"][k], hw[k], F["net"][k])
        return (g, r, pm) if pads else (g, r)

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
