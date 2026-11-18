"""Exact clearances from the board's vector data: holes to copper, copper to the board edge.

Every copper item is kept as KiCad describes it, as a core and a radius (the copper is every point
within the radius of the core):

    straight track     its centre line, half its width
    arc track          its chords, each with the arc's circle: the distance is to the arc itself
    via                its centre, half its diameter (on the layers where it is flashed)
    standard pad       padgeom.core: a point, a segment or a convex polygon, and a radius
    zone fill          its polygons' edges, as written in the file (KiCad's own fill, exact)

A distance to these is exact to the coordinates' precision (KiCad's unit, 1 nm); a limit is met when
the distance is within TOL (10 nm, the rounding of coordinates drawn at a limit) of it. Two items are not
exact this way, and carry an error bound instead: a custom pad, or a chamfered pad with round corners
(its polygon lies up to POLY_ERR inside its true shape), and the board outline where it is an arc (each
chord lies its sagitta inside the arc: well under 1 um). A measurement then gives a pair (lo, hi): the true distance
lies between them, and lo == hi when every item involved is exact.
"""
import math

import numpy as np
from scipy.spatial import cKDTree

from . import geom

TOL = 1e-5           # mm: the tolerance of an exact comparison, 10 nm. KiCad stores coordinates as whole nanometres, so
                     # copper drawn at a limit (a via placed at a computed point beside a turned pad) can sit a few nm
                     # short of it from rounding alone; no fabricator works to that
POLY_ERR = 0.002     # mm: how far inside its true shape a pad polygon with arcs may lie (KiCad: 2 um, openpcbsi: 1 um)
EDGE_ERR = 0.002     # mm: the same for the board outline's arcs
PIECE = 0.25         # mm: long edges are indexed as pieces no longer than this


# ------------------------------------------------------------------------------------------------
# distances (vectorised over the second argument)

def _cross(ax, ay, bx, by):
    return ax * by - ay * bx


def seg_segs(ax1, ay1, ax2, ay2, bx1, by1, bx2, by2):
    """Distance from one segment A (a point when its ends coincide) to arrays of segments B; 0 where
    they cross."""
    d = np.minimum.reduce([geom.point_seg_dist(ax1, ay1, bx1, by1, bx2, by2), geom.point_seg_dist(ax2, ay2, bx1, by1, bx2, by2),
                           geom.point_seg_dist(bx1, by1, ax1, ay1, ax2, ay2), geom.point_seg_dist(bx2, by2, ax1, ay1, ax2, ay2)])
    if ax1 == ax2 and ay1 == ay2:
        return d
    o1 = _cross(ax2 - ax1, ay2 - ay1, bx1 - ax1, by1 - ay1)
    o2 = _cross(ax2 - ax1, ay2 - ay1, bx2 - ax1, by2 - ay1)
    o3 = _cross(bx2 - bx1, by2 - by1, ax1 - bx1, ay1 - by1)
    o4 = _cross(bx2 - bx1, by2 - by1, ax2 - bx1, ay2 - by1)
    return np.where((o1 * o2 < 0) & (o3 * o4 < 0), 0.0, d)


def _in_wedge(vx, vy, ax, ay, bx, by):
    """Whether direction v lies within the minor arc's angle from a to b (all from the centre)."""
    s = _cross(ax, ay, bx, by)
    return (_cross(ax, ay, vx, vy) * s >= 0) & (_cross(vx, vy, bx, by) * s >= 0)


def point_arcs(px, py, x1, y1, x2, y2, cx, cy, R):
    """Distance from a point to the minor arcs from (x1, y1) to (x2, y2) on circles (cx, cy, R)."""
    vx, vy = px - cx, py - cy
    on = _in_wedge(vx, vy, x1 - cx, y1 - cy, x2 - cx, y2 - cy)
    ends = np.minimum(np.hypot(px - x1, py - y1), np.hypot(px - x2, py - y2))
    return np.where(on, np.abs(np.hypot(vx, vy) - R), ends)


def seg_arcs(ax1, ay1, ax2, ay2, x1, y1, x2, y2, cx, cy, R):
    """Distance from a segment A to minor arcs: 0 where it crosses one; otherwise the least of the ends
    of each to the other and of the point of A nearest the centre, which is where a common normal (a
    radius) can meet A."""
    d = np.minimum.reduce([point_arcs(ax1, ay1, x1, y1, x2, y2, cx, cy, R), point_arcs(ax2, ay2, x1, y1, x2, y2, cx, cy, R),
                           geom.point_seg_dist(x1, y1, ax1, ay1, ax2, ay2), geom.point_seg_dist(x2, y2, ax1, ay1, ax2, ay2)])
    dx, dy = ax2 - ax1, ay2 - ay1
    L2 = dx * dx + dy * dy
    if L2 == 0:
        return d
    t = np.clip(((cx - ax1) * dx + (cy - ay1) * dy) / L2, 0, 1)
    fx, fy = ax1 + t * dx, ay1 + t * dy
    d = np.minimum(d, np.where(_in_wedge(fx - cx, fy - cy, x1 - cx, y1 - cy, x2 - cx, y2 - cy), np.abs(np.hypot(fx - cx, fy - cy) - R), np.inf))
    # crossings: |A(t) - C| = R
    b = 2 * ((ax1 - cx) * dx + (ay1 - cy) * dy)
    c = (ax1 - cx) ** 2 + (ay1 - cy) ** 2 - R * R
    disc = b * b - 4 * L2 * c
    ok = disc >= 0
    sq = np.sqrt(np.where(ok, disc, 0))
    for sgn in (-1, 1):
        tt = (-b + sgn * sq) / (2 * L2)
        qx, qy = ax1 + tt * dx, ay1 + tt * dy
        hit = ok & (tt >= 0) & (tt <= 1) & _in_wedge(qx - cx, qy - cy, x1 - cx, y1 - cy, x2 - cx, y2 - cy)
        d = np.where(hit, 0.0, d)
    return d


def seg_core(ax1, ay1, ax2, ay2, pts):
    """Distance from segment A to a core (a point, a segment, or a polygon: 0 inside it)."""
    P = np.asarray(pts, dtype=float)
    if len(P) == 1:
        return float(geom.point_seg_dist(P[0, 0], P[0, 1], ax1, ay1, ax2, ay2))
    if len(P) == 2:
        return float(seg_segs(ax1, ay1, ax2, ay2, P[:1, 0], P[:1, 1], P[1:, 0], P[1:, 1])[0])
    if geom.points_in_rings(np.array([ax1]), np.array([ay1]), [P])[0]:
        return 0.0
    Q = np.roll(P, -1, axis=0)
    return float(seg_segs(ax1, ay1, ax2, ay2, P[:, 0], P[:, 1], Q[:, 0], Q[:, 1]).min())


def _pieces(segs, piece=PIECE):
    """Edges cut into pieces no longer than `piece` (exact: the same lines), and the source index."""
    if not len(segs):
        return np.zeros((0, 4)), np.zeros(0, int)
    L = np.hypot(segs[:, 2] - segs[:, 0], segs[:, 3] - segs[:, 1])
    n = np.maximum(1, np.ceil(L / piece).astype(int))
    src = np.repeat(np.arange(len(segs)), n)
    k = np.concatenate([np.arange(m) for m in n])
    t0, t1 = k / n[src], (k + 1) / n[src]
    s = segs[src]
    out = np.c_[s[:, 0] + t0 * (s[:, 2] - s[:, 0]), s[:, 1] + t0 * (s[:, 3] - s[:, 1]),
                s[:, 0] + t1 * (s[:, 2] - s[:, 0]), s[:, 1] + t1 * (s[:, 3] - s[:, 1])]
    # the original end points exactly
    first = k == 0
    last = k == n[src] - 1
    out[first, 0], out[first, 1] = s[first, 0], s[first, 1]
    out[last, 2], out[last, 3] = s[last, 2], s[last, 3]
    return out, src


# ------------------------------------------------------------------------------------------------

class Copper:
    """The copper of one layer as exact items, indexed for nearness queries."""

    def __init__(self, b, li):
        self.b, self.li = b, li
        T = b.tracks
        sel = np.flatnonzero(T["layer"] == li)
        arcs = [b.snap["tracks"][k].get("arc") for k in sel]
        is_arc = np.array([isinstance(a, (list, tuple)) for a in arcs], bool)
        self.t = {k: T[k][sel] for k in ("x1", "y1", "x2", "y2", "net")}
        self.t["r"] = T["w"][sel] / 2
        self.t["arc"] = is_arc
        self.t["c"] = np.array([a if isinstance(a, (list, tuple)) else (0.0, 0.0, 0.0) for a in arcs], float).reshape(-1, 3)
        V = b.vias
        vs = np.array([k for k in range(len(V["x"])) if li in b.via_flashed[k]], int)
        self.v = {"x": V["x"][vs], "y": V["y"][vs], "r": V["d"][vs] / 2, "net": V["net"][vs]} if len(vs) else None
        self.pads = []
        for p in b.pads:
            if li not in p["layer_ids"]:
                continue
            poly = (p.get("layer_poly") or {}).get(b.layers[li], p["poly"])
            core = (p.get("layer_core") or {}).get(b.layers[li], p.get("core")) if p.get("layer_poly") else p.get("core")
            if not poly:
                continue
            a = np.vstack([np.asarray(r, float) for pl in poly for r in pl])
            bb = (a[:, 0].min(), a[:, 1].min(), a[:, 0].max(), a[:, 1].max())
            if core and core.get("pts") is not None:
                item = ("core", core["pts"], core["r"], 0.0)
            else:
                item = ("poly", [np.asarray(r, float) for pl in poly for r in pl], 0.0, 0.0 if core else POLY_ERR)
            self.pads.append((bb, p["net_id"], item, p))
        self.pad_bb = np.array([x[0] for x in self.pads]).reshape(-1, 4)
        segs, nets = [], []
        for z in b.zones:
            if z["layer"] != b.layers[li]:
                continue
            s_ = geom.ring_segments([ring for poly in z["polys"] for ring in poly])
            if len(s_):
                segs.append(s_)
                nets.append(np.full(len(s_), b.net_id.get(z["net"], 0)))
        segs = np.vstack(segs) if segs else np.zeros((0, 4))
        self.z, src = _pieces(segs)
        self.z_net = np.concatenate(nets)[src] if nets else np.zeros(0, int)
        self.z_tree = cKDTree(np.c_[(self.z[:, 0] + self.z[:, 2]) / 2, (self.z[:, 1] + self.z[:, 3]) / 2]) if len(self.z) else None

    def clearance(self, ax1, ay1, ax2, ay2, rh, nid, reach):
        """The least distance from a hole (its axis A and radius rh: a stadium) to copper of any other net
        than `nid` (any copper when nid is 0) within `reach` of its wall: (lo, hi, net id, x, y) of the
        nearest item, or None. lo < 0 means the copper cuts into the hole."""
        other = (lambda n: (n != nid) & (n != 0)) if nid else (lambda n: n != 0)
        R = rh + reach
        x0, y0, x1, y1 = min(ax1, ax2) - R, min(ay1, ay2) - R, max(ax1, ax2) + R, max(ay1, ay2) + R
        best = None

        def take(lo, hi, n, x, y):
            nonlocal best
            if best is None or lo < best[0] or (lo == best[0] and hi < best[1]):
                best = (float(lo), float(hi), int(n), float(x), float(y))
        # tracks
        t = self.t
        if len(t["x1"]):
            m = other(t["net"]) & (np.minimum(t["x1"], t["x2"]) - t["r"] <= x1) & (np.maximum(t["x1"], t["x2"]) + t["r"] >= x0) \
                & (np.minimum(t["y1"], t["y2"]) - t["r"] <= y1) & (np.maximum(t["y1"], t["y2"]) + t["r"] >= y0)
            k = np.flatnonzero(m)
            if len(k):
                d = seg_segs(ax1, ay1, ax2, ay2, t["x1"][k], t["y1"][k], t["x2"][k], t["y2"][k])
                ka = t["arc"][k]
                if ka.any():
                    c = t["c"][k[ka]]
                    d[ka] = seg_arcs(ax1, ay1, ax2, ay2, t["x1"][k[ka]], t["y1"][k[ka]], t["x2"][k[ka]], t["y2"][k[ka]], c[:, 0], c[:, 1], c[:, 2])
                d = d - rh - t["r"][k]
                j = int(np.argmin(d))
                take(d[j], d[j], t["net"][k[j]], (t["x1"][k[j]] + t["x2"][k[j]]) / 2, (t["y1"][k[j]] + t["y2"][k[j]]) / 2)
        # vias
        v = self.v
        if v is not None:
            m = other(v["net"]) & (v["x"] - v["r"] <= x1) & (v["x"] + v["r"] >= x0) & (v["y"] - v["r"] <= y1) & (v["y"] + v["r"] >= y0)
            k = np.flatnonzero(m)
            if len(k):
                d = geom.point_seg_dist(v["x"][k], v["y"][k], ax1, ay1, ax2, ay2) - rh - v["r"][k]
                j = int(np.argmin(d))
                take(d[j], d[j], v["net"][k[j]], v["x"][k[j]], v["y"][k[j]])
        # pads
        if len(self.pad_bb):
            pb = self.pad_bb
            for j in np.flatnonzero((pb[:, 0] <= x1) & (pb[:, 2] >= x0) & (pb[:, 1] <= y1) & (pb[:, 3] >= y0)):
                _, n, (kind, data, r, err), p = self.pads[j]
                if not other(np.array(n)):
                    continue
                if kind == "core":
                    d = seg_core(ax1, ay1, ax2, ay2, data) - r - rh
                else:
                    inside = geom.points_in_rings(np.array([ax1]), np.array([ay1]), data)[0]
                    d = -rh if inside else float(min(seg_segs(ax1, ay1, ax2, ay2, *np.roll(np.c_[q, np.roll(q, -1, axis=0)], 0, axis=0).T).min()
                                                     for q in data)) - rh
                take(d - err, d, n, p["x"], p["y"])
        # zone fills: the edges near the hole; a hole inside another net's fill is found on the raster
        # (a centre inside copper but further than a pixel from its edge is an overlap either way)
        if self.z_tree is not None:
            cx, cy = (ax1 + ax2) / 2, (ay1 + ay2) / 2
            half = math.hypot(ax2 - ax1, ay2 - ay1) / 2
            k = np.array(self.z_tree.query_ball_point((cx, cy), R + half + PIECE / 2 + 1e-9), int)
            if len(k):
                k = k[other(self.z_net[k])]
            if len(k):
                d = seg_segs(ax1, ay1, ax2, ay2, self.z[k, 0], self.z[k, 1], self.z[k, 2], self.z[k, 3]) - rh
                j = int(np.argmin(d))
                take(d[j], d[j], self.z_net[k[j]], (self.z[k[j], 0] + self.z[k[j], 2]) / 2, (self.z[k[j], 1] + self.z[k[j], 3]) / 2)
        g = self.b.grid()
        i, jj = g.ij((ax1 + ax2) / 2, (ay1 + ay2) / 2)
        r_ = self.b.raster(self.li)
        if 0 <= int(i) < r_.shape[0] and 0 <= int(jj) < r_.shape[1]:
            n = int(r_[int(i), int(jj)])
            if n and other(np.array(n)) and (best is None or best[0] > -rh):
                take(-rh, -rh, n, (ax1 + ax2) / 2, (ay1 + ay2) / 2)
        return best

    def edge(self, outline_segs, tree, limit, piece_err=None, skip=()):
        """Copper of this layer nearer the outline than `limit`: (lo, hi, net id, x, y) of the nearest,
        or None. Copper outside the outline is at distance 0. Items whose reference point lies in one of
        the `skip` circles (x, y, radius: castellated pads, which reach the edge by design) are not
        judged."""
        best = None
        sk = np.asarray(skip, float).reshape(-1, 3)
        if not len(outline_segs):
            return None
        reach = limit + PIECE / 2 + 1e-9
        rings = self.b.outline_rings()

        def near(px, py, extra):
            return tree.query_ball_point(np.c_[px, py], reach + extra)

        def take(lo, hi, n, x, y):
            nonlocal best
            if len(sk) and (np.hypot(sk[:, 0] - x, sk[:, 1] - y) < sk[:, 2]).any():
                return
            if best is None or lo < best[0]:
                best = (float(lo), float(hi), int(n), float(x), float(y))
        O = outline_segs

        E = piece_err if piece_err is not None else np.full(len(O), EDGE_ERR)
        last_err = [0.0]
        last_j = [0]

        def closest(ax1, ay1, ax2, ay2, ks):
            """The least distance to the outline pieces ks; its piece's error bound is kept in last_err."""
            ks = np.asarray(ks, int)
            if not len(ks):
                return np.inf
            d = seg_segs(ax1, ay1, ax2, ay2, O[ks, 0], O[ks, 1], O[ks, 2], O[ks, 3])
            j = int(np.argmin(d))
            last_err[0] = float(E[ks[j]])
            last_j[0] = int(ks[j])
            return float(d[j])

        def at(x1, y1, x2, y2):
            """The point of a segment nearest the outline piece found last (where a finding is placed)."""
            q = O[last_j[0]]
            dx, dy = x2 - x1, y2 - y1
            L2 = dx * dx + dy * dy
            best_p, best_d = (x1, y1), np.inf
            for px, py in ((q[0], q[1]), (q[2], q[3]), ((q[0] + q[2]) / 2, (q[1] + q[3]) / 2)):
                t = 0.0 if L2 == 0 else min(max(((px - x1) * dx + (py - y1) * dy) / L2, 0.0), 1.0)
                cx_, cy_ = x1 + t * dx, y1 + t * dy
                dd = math.hypot(cx_ - px, cy_ - py)
                if dd < best_d:
                    best_p, best_d = (cx_, cy_), dd
            return best_p
        t = self.t
        if len(t["x1"]):
            mx, my = (t["x1"] + t["x2"]) / 2, (t["y1"] + t["y2"]) / 2
            half = np.hypot(t["x2"] - t["x1"], t["y2"] - t["y1"]) / 2 + t["r"]
            for k, ks in enumerate(tree.query_ball_point(np.c_[mx, my], reach + half)):
                if ks and t["net"][k] != -1:
                    # an arc chord lies inside its arc by its sagitta: a bound, not exact
                    sag = 0.0
                    if t["arc"][k]:
                        R = t["c"][k, 2]
                        sag = R - math.sqrt(max(R * R - (half[k] - t["r"][k]) ** 2, 0.0))
                    d = closest(t["x1"][k], t["y1"][k], t["x2"][k], t["y2"][k], ks) - t["r"][k]
                    take(d - sag - last_err[0], d, t["net"][k], *at(t["x1"][k], t["y1"][k], t["x2"][k], t["y2"][k]))
        v = self.v
        if v is not None:
            for k, ks in enumerate(tree.query_ball_point(np.c_[v["x"], v["y"]], reach + v["r"])):
                if ks:
                    d = closest(v["x"][k], v["y"][k], v["x"][k], v["y"][k], ks) - v["r"][k]
                    take(d - last_err[0], d, v["net"][k], v["x"][k], v["y"][k])
        for bb, n, (kind, data, r, err), p in self.pads:
            cx, cy = (bb[0] + bb[2]) / 2, (bb[1] + bb[3]) / 2
            ks = tree.query_ball_point((cx, cy), reach + math.hypot(bb[2] - bb[0], bb[3] - bb[1]) / 2)
            if not ks:
                continue
            ks = np.asarray(ks, int)
            if kind == "core":
                P = np.asarray(data, float)
                segs_ = [(P[0], P[-1])] if len(P) <= 2 else [(P[i], P[(i + 1) % len(P)]) for i in range(len(P))]
            else:
                segs_ = [(q[i], q[(i + 1) % len(q)]) for q in data for i in range(len(q))]
            d, e_ = np.inf, 0.0
            for a_, c_ in segs_:
                dd = closest(a_[0], a_[1], c_[0], c_[1], ks)
                if dd < d:
                    d, e_ = dd, last_err[0]
            d -= r
            pts = np.asarray(data[0] if kind == "poly" else data, float)
            if rings and not geom.points_in_rings(pts[:1, 0], pts[:1, 1], rings)[0]:
                d = 0.0
            take(d - err - e_, d, n, p["x"], p["y"])
        if len(self.z):
            z = self.z
            mx, my = (z[:, 0] + z[:, 2]) / 2, (z[:, 1] + z[:, 3]) / 2
            cand = [k for k, ks in enumerate(tree.query_ball_point(np.c_[mx, my], reach + PIECE / 2)) if ks]
            if cand:
                cand = np.asarray(cand, int)
                outside = rings and ~geom.points_in_rings(mx[cand], my[cand], rings)
                for i, k in enumerate(cand):
                    ks = tree.query_ball_point((mx[k], my[k]), reach + PIECE / 2)
                    if rings is not None and len(rings) and outside[i]:
                        d, e_ = 0.0, 0.0
                    else:
                        d = closest(z[k, 0], z[k, 1], z[k, 2], z[k, 3], ks)
                        e_ = last_err[0]
                    take(d - e_, d, self.z_net[k], mx[k], my[k])
        return best


def outline_index(b):
    """The board outline's edges as pieces, with a nearness index and each piece's error bound: 0 for a
    straight edge, the sagitta of its chord where the edge is a chord of one of the outline's arcs
    (snapshot `outline_arcs`)."""
    s_ = geom.ring_segments(b.outline_rings())
    O, src = _pieces(s_)
    err = np.zeros(len(s_))
    for cx, cy, R in b.snap.get("outline_arcs") or []:
        on = (np.abs(np.hypot(s_[:, 0] - cx, s_[:, 1] - cy) - R) < 1e-5) & (np.abs(np.hypot(s_[:, 2] - cx, s_[:, 3] - cy) - R) < 1e-5)
        half = np.hypot(s_[:, 2] - s_[:, 0], s_[:, 3] - s_[:, 1]) / 2
        err = np.where(on, np.maximum(err, R - np.sqrt(np.maximum(R * R - half * half, 0.0))), err)
    if "outline_arcs" not in b.snap:
        err[:] = EDGE_ERR                          # an older snapshot: no arc information, the general bound
    return O, (cKDTree(np.c_[(O[:, 0] + O[:, 2]) / 2, (O[:, 1] + O[:, 3]) / 2]) if len(O) else None), err[src] if len(src) else np.zeros(0)
