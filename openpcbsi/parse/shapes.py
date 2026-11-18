"""Pad shapes, margins and board outlines as polygons, computed by openpcbsi from the file's own fields,
so that every reader (and so every check) gets the same geometry. Checked against KiCad's own polygons
(the pcbnew golden snapshots) and KiCad's IPC-2581 export on the reference corpus.

Coordinates are KiCad's: mm, y grows downwards, a positive angle turns counter-clockwise on screen.
Arcs are cut into chords no more than ARC_ERR from the true curve (KiCad's own polygons use 2 um).
"""
import math

ARC_ERR = 0.001          # mm: largest distance of a chord from its arc
JOIN = 0.005             # mm: outline ends this close are one point


def _qs(r):
    """Segments per quarter circle for a radius, so that the chord error stays below ARC_ERR."""
    if r <= ARC_ERR:
        return 4
    return max(4, int(math.ceil((math.pi / 2) / math.acos(max(-1.0, 1 - ARC_ERR / r)))))


def _sh():
    from shapely import affinity, geometry, ops  # noqa: the parsers' site is on the path (openpcbsi.parse.use)
    return geometry, affinity, ops


def to_board(g, cx, cy, deg):
    """A pad-frame geometry placed at (cx, cy) and turned by `deg` (KiCad's sense)."""
    geometry, affinity, _ = _sh()
    a = math.radians(deg or 0.0)
    c, s = math.cos(a), math.sin(a)
    return affinity.affine_transform(g, [c, s, -s, c, cx, cy])


def rect(w, h):
    geometry, _, _ = _sh()
    return geometry.box(-w / 2, -h / 2, w / 2, h / 2)


def roundrect(w, h, rratio=0.0, chamfer=0.0, corners=()):
    """KiCad's rounded rectangle (radius = ratio x the smaller side), with chamfered corners cut off."""
    geometry, _, _ = _sh()
    r = max(0.0, min(rratio, 0.5)) * min(w, h)
    if r > 1e-9:
        g = geometry.box(-w / 2 + r, -h / 2 + r, w / 2 - r, h / 2 - r).buffer(r, quad_segs=_qs(r))
    else:
        g = rect(w, h)
    c = max(0.0, min(chamfer, 0.5)) * min(w, h)
    if c > 1e-9:
        for name in corners or ():
            sx = 1 if "right" in name else -1
            sy = -1 if "top" in name else 1
            x0, y0 = sx * w / 2, sy * h / 2
            tri = geometry.Polygon([(x0, y0), (x0 - sx * c, y0), (x0, y0 - sy * c)])
            g = g.difference(tri.buffer(1e-9))
    return g


def oval(w, h):
    geometry, _, _ = _sh()
    if abs(w - h) < 1e-9:
        return geometry.Point(0, 0).buffer(w / 2, quad_segs=_qs(w / 2))
    r = min(w, h) / 2
    half = (max(w, h) - min(w, h)) / 2
    line = geometry.LineString([(-half, 0), (half, 0)]) if w > h else geometry.LineString([(0, -half), (0, half)])
    return line.buffer(r, quad_segs=_qs(r))


def trapezoid(w, h, dx, dy):
    """KiCad's trapezoid: `rect_delta` (dx, dy) widens one side and narrows the opposite one."""
    geometry, _, _ = _sh()
    ddx, ddy = dx / 2, dy / 2
    return geometry.Polygon([(-w / 2 - ddy, h / 2 + ddx), (w / 2 + ddy, h / 2 - ddx),
                             (w / 2 - ddy, -h / 2 + ddx), (-w / 2 + ddy, -h / 2 - ddx)])


def _arc_pts(s, m, e):
    (x1, y1), (x2, y2), (x3, y3) = s, m, e
    d = 2 * (x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2))
    if abs(d) < 1e-12:
        return [s, e]
    ux = ((x1 ** 2 + y1 ** 2) * (y2 - y3) + (x2 ** 2 + y2 ** 2) * (y3 - y1) + (x3 ** 2 + y3 ** 2) * (y1 - y2)) / d
    uy = ((x1 ** 2 + y1 ** 2) * (x3 - x2) + (x2 ** 2 + y2 ** 2) * (x1 - x3) + (x3 ** 2 + y3 ** 2) * (x2 - x1)) / d
    r = math.hypot(x1 - ux, y1 - uy)
    a1, a2, a3 = (math.atan2(p[1] - uy, p[0] - ux) for p in (s, m, e))
    norm = lambda a: (a + 2 * math.pi) % (2 * math.pi)
    sweep = norm(a3 - a1)
    if not norm(a2 - a1) < sweep:
        sweep -= 2 * math.pi
    step = 2 * math.acos(max(-1.0, 1 - ARC_ERR / r)) if r > ARC_ERR else math.pi / 8
    n = max(2, int(math.ceil(abs(sweep) / step)))
    pts = [(ux + r * math.cos(a1 + sweep * k / n), uy + r * math.sin(a1 + sweep * k / n)) for k in range(n + 1)]
    pts[0], pts[-1] = tuple(s), tuple(e)              # the given ends exactly: chains must join
    return pts


def _num(sx, key):
    """Numbers of a sub-expression like ['start', x, y] in a parsed s-expression list."""
    for it in sx[1:]:
        if isinstance(it, list) and it and str(it[0]) == key:
            return [float(v) for v in it[1:] if isinstance(v, (int, float))]
    return None


def _flag(sx, key):
    for it in sx[1:]:
        if isinstance(it, list) and it and str(it[0]) == key:
            return str(it[1]) if len(it) > 1 else ""
    return None


def primitive(sx):
    """A custom pad's primitive (gr_poly, gr_rect, gr_circle, gr_line, gr_arc), from its s-expression,
    in the pad's frame."""
    geometry, _, ops = _sh()
    kind = str(sx[0])
    width = (_num(sx, "width") or [None])[0]
    if width is None:                                  # KiCad 9/10: (stroke (width w) (type solid))
        st = next((it for it in sx[1:] if isinstance(it, list) and it and str(it[0]) == "stroke"), None)
        width = ((_num(st, "width") if st else None) or [0.0])[0]
    fill = _flag(sx, "fill")
    filled = fill not in (None, "no", "none")
    if kind == "gr_poly":
        pts_e = next((it for it in sx[1:] if isinstance(it, list) and it and str(it[0]) == "pts"), [])
        pts = []
        for p in pts_e[1:]:
            if isinstance(p, list) and str(p[0]) == "xy":
                pts.append((float(p[1]), float(p[2])))
            elif isinstance(p, list) and str(p[0]) == "arc":
                s_, m_, e_ = _num(p, "start"), _num(p, "mid"), _num(p, "end")
                pts += _arc_pts(tuple(s_), tuple(m_), tuple(e_))
        if len(pts) < 3:
            return None
        poly = geometry.Polygon(pts).buffer(0)
        if fill is None or filled:                 # a pad polygon is filled unless it says otherwise
            return poly.buffer(width / 2, quad_segs=_qs(width / 2)) if width > 0 else poly
        return geometry.LinearRing(pts).buffer(max(width, 1e-6) / 2, quad_segs=_qs(width / 2))
    if kind == "gr_rect":
        (x0, y0), (x1, y1) = _num(sx, "start"), _num(sx, "end")
        b = geometry.box(min(x0, x1), min(y0, y1), max(x0, x1), max(y0, y1))
        if filled:
            return b.buffer(width / 2, join_style=2) if width > 0 else b
        return b.exterior.buffer(max(width, 1e-6) / 2, join_style=2)
    if kind == "gr_circle":
        (cx, cy), (ex, ey) = _num(sx, "center"), _num(sx, "end")
        r = math.hypot(ex - cx, ey - cy)
        if filled:
            return geometry.Point(cx, cy).buffer(r + width / 2, quad_segs=_qs(r + width / 2))
        return geometry.Point(cx, cy).buffer(r, quad_segs=_qs(r)).exterior.buffer(max(width, 1e-6) / 2, quad_segs=_qs(width / 2))
    if kind == "gr_line":
        (x0, y0), (x1, y1) = _num(sx, "start"), _num(sx, "end")
        return geometry.LineString([(x0, y0), (x1, y1)]).buffer(max(width, 1e-6) / 2, quad_segs=_qs(width / 2))
    if kind == "gr_arc":
        pts = _arc_pts(tuple(_num(sx, "start")), tuple(_num(sx, "mid")), tuple(_num(sx, "end")))
        return geometry.LineString(pts).buffer(max(width, 1e-6) / 2, quad_segs=_qs(width / 2))
    return None


def pad_local(shape, sx, sy, rratio=0.0, chamfer=0.0, corners=(), delta=(0.0, 0.0), anchor="rect", prims=()):
    """A pad's shape in its own frame. `shape` is KiCad's name (rect, roundrect, circle, oval, trapezoid,
    custom); a custom pad is its anchor joined with its primitives (parsed s-expressions)."""
    shape = (shape or "rect").lower()
    if shape == "circle":
        return oval(sx, sx)
    if shape == "oval":
        return oval(sx, sy)
    if shape == "roundrect" or shape == "chamfered_rect":
        return roundrect(sx, sy, rratio, chamfer, corners)
    if shape == "trapezoid":
        return trapezoid(sx, sy, *delta) if any(abs(v) > 1e-9 for v in delta) else rect(sx, sy)
    if shape == "custom":
        _, _, ops = _sh()
        parts = [oval(sx, sx) if anchor == "circle" else rect(sx, sy)]
        parts += [g for g in (primitive(p) for p in prims) if g is not None and not g.is_empty]
        return ops.unary_union(parts)
    return rect(sx, sy)


def grow(g, amount):
    """A pad shape grown (or shrunk) by a margin, with round corners as KiCad grows pads for a mask opening."""
    if abs(amount) < 1e-12:
        return g
    if amount > 0:
        return g.buffer(amount, quad_segs=_qs(amount))
    return g.buffer(amount, join_style=2)


def rings(g):
    """A shapely (Multi)Polygon as the snapshot's polygons: [[outline, hole, ...], ...], points in mm."""
    if g is None or g.is_empty:
        return []
    geoms = getattr(g, "geoms", [g])
    out = []
    for p in geoms:
        if p.geom_type != "Polygon" or p.is_empty:
            continue
        rs = [[(round(x, 6), round(y, 6)) for x, y in list(p.exterior.coords)[:-1]]]
        rs += [[(round(x, 6), round(y, 6)) for x, y in list(h.coords)[:-1]] for h in p.interiors]
        out.append(rs)
    return out


def outline(items):
    """Board outline polygons from Edge.Cuts items given as point chains (closed or open) in board
    coordinates: [[outline, cutout, ...], ...] as KiCad's GetBoardPolygonOutlines gives them."""
    geometry, _, ops = _sh()
    # ends that almost meet are joined, as KiCad chains its outline (ends within JOIN of each other)
    chains = [[(float(x), float(y)) for x, y in pts] for pts in items if len(pts) >= 2]
    if chains:
        import numpy as np
        from scipy.spatial import cKDTree
        E = np.array([c[0] for c in chains] + [c[-1] for c in chains])
        rep = list(range(len(E)))
        for i, j in sorted(cKDTree(E).query_pairs(JOIN)):
            ri, rj = rep[i], rep[j]
            while rep[ri] != ri:
                ri = rep[ri]
            while rep[rj] != rj:
                rj = rep[rj]
            if ri != rj:
                rep[rj] = ri
        root = []
        for k in range(len(E)):
            r = k
            while rep[r] != r:
                r = rep[r]
            root.append(r)
        n = len(chains)
        for k, c in enumerate(chains):
            c[0] = tuple(E[root[k]])
            c[-1] = tuple(E[root[n + k]])
    lines = []
    for pts in chains:
        pts = [(round(x, 6), round(y, 6)) for x, y in pts]
        if len(pts) >= 2:
            lines.append(geometry.LineString(pts))
    if not lines:
        return []
    merged = ops.unary_union(lines)
    faces = [geometry.Polygon(p.exterior) for p in ops.polygonize(merged)]
    rep = [f.representative_point() for f in faces]
    depth = [sum(1 for j, g in enumerate(faces) if j != i and g.contains(rep[i])) for i in range(len(faces))]
    out = []
    for i, f in enumerate(faces):
        if depth[i] % 2:
            continue                                  # a cutout: subtracted from the face directly around it
        inner = [faces[j] for j in range(len(faces)) if depth[j] == depth[i] + 1 and f.contains(rep[j])]
        g = f
        for h in inner:
            g = g.difference(h)
        out += rings(g)
    return out


def padstack(raw):
    """(mode, {layer name: entry}) of a pad's (padstack ...) block, from its parsed s-expression; mode is
    None when the pad has none."""
    if not raw:
        return None, {}
    ps = next((it for it in raw[1:] if isinstance(it, list) and it and str(it[0]) == "padstack"), None)
    if ps is None:
        return None, {}
    mode = _flag(ps, "mode") or "normal"
    entries = {str(it[1]): it for it in ps[1:] if isinstance(it, list) and len(it) > 1 and str(it[0]) == "layer"}
    return mode, entries


def entry_params(e, base):
    """A padstack layer entry's shape parameters, the pad's own where the entry gives none."""
    q = dict(base)
    sh = _flag(e, "shape")
    if sh:
        q["shape"] = sh
    sz = _num(e, "size")
    if sz:
        q["sx"], q["sy"] = sz[0], sz[1] if len(sz) > 1 else sz[0]
    for key, name in (("rratio", "roundrect_rratio"), ("chamfer", "chamfer_ratio")):
        v = _num(e, name)
        if v:
            q[key] = v[0]
    ch = next((it for it in e[1:] if isinstance(it, list) and it and str(it[0]) == "chamfer"), None)
    if ch is not None:
        q["corners"] = [str(c).lower() for c in ch[1:]]
    d = _num(e, "rect_delta")
    if d:
        q["delta"] = (d[0], d[1] if len(d) > 1 else 0.0)
    opt = next((it for it in e[1:] if isinstance(it, list) and it and str(it[0]) == "options"), None)
    if opt is not None:
        q["anchor"] = _flag(opt, "anchor") or q["anchor"]
    pr = next((it for it in e[1:] if isinstance(it, list) and it and str(it[0]) == "primitives"), None)
    if pr is not None:
        q["prims"] = [it for it in pr[1:] if isinstance(it, list)]
    return q
