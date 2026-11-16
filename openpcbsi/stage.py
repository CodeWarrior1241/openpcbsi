"""Layout stage (#29): placement, routing or final, from the board's own copper.

A net is incomplete when its pads fall into more than one copper island (pads, tracks, vias and zone
fills joined where they touch). The stage is `placement` when no signal track exists yet, `routing`
while any net with two or more pads is incomplete, `final` otherwise. During placement and routing,
findings that only say "this is not connected yet" become `pending` instead of `fail`.
"""
import math
from collections import defaultdict


class _UF:
    def __init__(self):
        self.p = {}

    def find(self, a):
        p = self.p
        p.setdefault(a, a)
        while p[a] != a:
            p[a] = p[p[a]]
            a = p[a]
        return a

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.p[ra] = rb


def _pip(x, y, ring):
    inside = False
    for i in range(len(ring)):
        (xa, ya), (xb, yb) = ring[i - 1], ring[i]
        if (ya > y) != (yb > y) and x < xa + (y - ya) * (xb - xa) / (yb - ya):
            inside = not inside
    return inside


def analyse(snap):
    by_net = defaultdict(lambda: {"pads": [], "tracks": [], "vias": [], "zones": []})
    for i, p in enumerate(snap["pads"]):
        if p["net"]:
            by_net[p["net"]]["pads"].append(i)
    for i, t in enumerate(snap["tracks"]):
        if t["net"]:
            by_net[t["net"]]["tracks"].append(i)
    for i, v in enumerate(snap["vias"]):
        if v["net"]:
            by_net[v["net"]]["vias"].append(i)
    for i, z in enumerate(snap["zones"]):
        if z.get("net") and z.get("polys"):
            by_net[z["net"]]["zones"].append(i)
    incomplete, open_pads = [], 0
    for net, it in by_net.items():
        if len(it["pads"]) < 2:
            continue
        uf = _UF()
        pts = []      # (x, y, layer or None for all layers, reach, node)
        for i in it["pads"]:
            p = snap["pads"][i]
            r = max(p.get("sx") or 0, p.get("sy") or 0) / 2
            for l in p["layers"] or [None]:
                pts.append((p["x"], p["y"], l, r, ("p", i)))
        for i in it["tracks"]:
            t = snap["tracks"][i]
            uf.union(("t", i, 0), ("t", i, 1))
            pts.append((t["x1"], t["y1"], t["layer"], t["w"] / 2, ("t", i, 0)))
            pts.append((t["x2"], t["y2"], t["layer"], t["w"] / 2, ("t", i, 1)))
        for i in it["vias"]:
            v = snap["vias"][i]
            for l in v["span"]:
                pts.append((v["x"], v["y"], l, v["d"] / 2, ("v", i)))
        # join points that touch (same layer, or a through pad on any layer)
        cell = 0.5
        grid = defaultdict(list)
        for k, (x, y, l, r, node) in enumerate(pts):
            grid[(int(x // cell), int(y // cell))].append(k)
        for k, (x, y, l, r, node) in enumerate(pts):
            gx, gy = int(x // cell), int(y // cell)
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    for k2 in grid.get((gx + dx, gy + dy), ()):
                        if k2 <= k:
                            continue
                        x2, y2, l2, r2, node2 = pts[k2]
                        if l is not None and l2 is not None and l != l2:
                            continue
                        if math.hypot(x - x2, y - y2) <= r + r2 + 1e-3:
                            uf.union(node, node2)
        # a point inside a zone fill of its layer joins that fill
        for i in it["zones"]:
            z = snap["zones"][i]
            for j, poly in enumerate(z["polys"]):
                for (x, y, l, r, node) in pts:
                    if l not in (None, z["layer"]):
                        continue
                    if _pip(x, y, poly[0]) and not any(_pip(x, y, h) for h in poly[1:]):
                        uf.union(node, ("z", i, j))
        # tracks crossing pads mid-segment (a pad on the track, not at its end)
        for i in it["tracks"]:
            t = snap["tracks"][i]
            for j in it["pads"]:
                p = snap["pads"][j]
                if p["layers"] and t["layer"] not in p["layers"]:
                    continue
                ax, ay, bx, by = t["x1"], t["y1"], t["x2"], t["y2"]
                L2 = (bx - ax) ** 2 + (by - ay) ** 2
                u = 0 if L2 == 0 else max(0, min(1, ((p["x"] - ax) * (bx - ax) + (p["y"] - ay) * (by - ay)) / L2))
                if math.hypot(ax + u * (bx - ax) - p["x"], ay + u * (by - ay) - p["y"]) <= max(p.get("sx") or 0, p.get("sy") or 0) / 2 + t["w"] / 2:
                    uf.union(("t", i, 0), ("p", j))
        roots = defaultdict(list)
        for i in it["pads"]:
            roots[uf.find(("p", i))].append(i)
        if len(roots) > 1:
            incomplete.append(net)
            open_pads += sum(len(v) for v in sorted(roots.values(), key=len)[:-1])
    signal_tracks = [t for t in snap["tracks"] if t["net"]]
    stage = "placement" if not signal_tracks else ("routing" if incomplete else "final")
    return {"stage": stage, "incomplete_nets": sorted(incomplete), "unrouted_pads": open_pads}
