"""Pad sizes measured from the pad's own polygons, never from its anchor.

A custom pad in KiCad has an anchor (a small rectangle or circle, often 0.2 mm) and primitives that
make the real copper; its "size" is the anchor's. Every check that needs a pad's size reads these
measurements instead:

    w, h          extent of the polygon in the pad's own frame (rotated back by the pad's angle), mm
    area, perim   area and outer perimeter of the polygon (holes subtracted from the area), mm2, mm
    inscribed_d   diameter of the largest disc that fits inside (the aperture's narrowest "width" for a
                  rounded or L-shaped opening; min(w, h) for a rectangle), mm

Pure Python (numpy, scipy): the snapshot computes it from KiCad's polygons, and so does every file
parser, so the numbers do not depend on which reader produced the polygons.
"""
import math

import numpy as np


def _rings(polys):
    return [np.asarray(r, dtype=float) for poly in polys or [] for r in poly if len(r) >= 3]


def _ring_area(r):
    x, y = r[:, 0], r[:, 1]
    return 0.5 * float(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))


def _ring_perim(r):
    d = np.diff(np.vstack([r, r[:1]]), axis=0)
    return float(np.hypot(d[:, 0], d[:, 1]).sum())


def measure(polys, cx=0.0, cy=0.0, rot_deg=0.0, symmetric=False):
    """Measurements of one pad's polygons (the [outline, hole, ...] lists of a snapshot), or None
    when there is no polygon. `symmetric`: the polygon is a rectangle, rounded rectangle, oval or
    circle, whose largest inscribed disc is exactly min(w, h) (no raster needed)."""
    if not polys:
        return None
    a = math.radians(rot_deg or 0.0)
    ca, sa = math.cos(a), math.sin(a)
    area, perim, pts = 0.0, 0.0, []
    for poly in polys:
        rings = [np.asarray(r, dtype=float) for r in poly if len(r) >= 3]
        if not rings:
            continue
        area += abs(_ring_area(rings[0])) - sum(abs(_ring_area(h)) for h in rings[1:])
        perim += _ring_perim(rings[0])
        pts.append(rings[0])
    if not pts:
        return None
    P = np.vstack(pts)
    # KiCad: y grows downwards and a positive angle turns counter-clockwise on screen; turning the
    # points back by the angle about the centre gives the pad's own frame (its extent in x and y)
    dx, dy = P[:, 0] - cx, P[:, 1] - cy
    u = dx * ca - dy * sa
    v = dx * sa + dy * ca
    w, h = float(u.max() - u.min()), float(v.max() - v.min())
    return {"w": round(w, 5), "h": round(h, 5), "area": round(area, 6), "perim": round(perim, 5),
            "inscribed_d": round(min(w, h) if symmetric else inscribed_diameter(polys, min(w, h)), 4)}


def inscribed_diameter(polys, scale=None):
    """Diameter of the largest disc inside the polygons: twice the largest distance from an inside
    point to the boundary, on a raster of at most 2 um (1/100 of the smaller extent)."""
    from scipy import ndimage
    from matplotlib.path import Path
    rings = _rings(polys)
    if not rings:
        return 0.0
    P = np.vstack(rings)
    x0, y0 = P.min(axis=0)
    x1, y1 = P.max(axis=0)
    s = scale if scale and scale > 0 else min(x1 - x0, y1 - y0)
    res = min(0.002, max(s / 100.0, 1e-4))
    pad = 2 * res
    nx, ny = int(math.ceil((x1 - x0 + 2 * pad) / res)) + 1, int(math.ceil((y1 - y0 + 2 * pad) / res)) + 1
    if nx * ny > 4_000_000:                       # very large polygon: coarser raster
        f = math.sqrt(nx * ny / 4_000_000)
        res *= f
        nx, ny = int(nx / f) + 1, int(ny / f) + 1
    xs = x0 - pad + (np.arange(nx) + 0.5) * res
    ys = y0 - pad + (np.arange(ny) + 0.5) * res
    X, Y = np.meshgrid(xs, ys)
    pts = np.c_[X.ravel(), Y.ravel()]
    inside = np.zeros(len(pts), bool)
    for poly in polys:
        rs = [np.asarray(r, dtype=float) for r in poly if len(r) >= 3]
        if not rs:
            continue
        m = Path(rs[0]).contains_points(pts)
        for hole in rs[1:]:
            m &= ~Path(hole).contains_points(pts)
        inside |= m
    inside = inside.reshape(ny, nx)
    if not inside.any():
        return 0.0
    d = ndimage.distance_transform_edt(inside) * res
    # the transform measures to the centre of the first outside cell: half a cell beyond the edge
    return max(0.0, 2.0 * (float(d.max()) - 0.5 * res))


def core(shape, sx, sy, cx, cy, rot_deg=0.0, corner_r=0.0, delta=(0.0, 0.0), chamfered=False):
    """A standard pad's exact shape as a core and a radius: the copper is every point within `r` of the
    core, which is one point (circle), a segment (oval), or a convex polygon (rectangle, rounded
    rectangle's inner rectangle, trapezoid). Distances to it are then exact, with no arc approximation:
    distance(core) - r. In board coordinates (KiCad's: y down, a positive angle turns counter-clockwise
    on screen). None for shapes not exact this way (custom, chamfered with round corners); pts None for
    a shape whose polygon has straight edges only (the polygon itself is exact)."""
    shape = (shape or "rect").lower()
    a = math.radians(rot_deg or 0.0)
    c, s = math.cos(a), math.sin(a)
    place = lambda pts: [[cx + c * x + s * y, cy - s * x + c * y] for x, y in pts]
    if shape == "circle":
        return {"r": sx / 2, "pts": place([(0.0, 0.0)])}
    if shape == "oval":
        if abs(sx - sy) < 1e-9:
            return {"r": sx / 2, "pts": place([(0.0, 0.0)])}
        half = (max(sx, sy) - min(sx, sy)) / 2
        return {"r": min(sx, sy) / 2, "pts": place([(-half, 0.0), (half, 0.0)] if sx > sy else [(0.0, -half), (0.0, half)])}
    if shape in ("roundrect", "chamfered_rect"):
        r = max(0.0, min(corner_r, min(sx, sy) / 2))
        if chamfered:
            # with round corners: not exact this way; without: straight edges only, so the pad's own
            # polygon is exact (pts None tells the caller to use it)
            return None if r > 1e-9 else {"r": 0.0, "pts": None}
        w, h = sx / 2 - r, sy / 2 - r
        return {"r": r, "pts": place([(-w, -h), (w, -h), (w, h), (-w, h)])}
    if shape == "trapezoid":
        ddx, ddy = delta[0] / 2, delta[1] / 2
        return {"r": 0.0, "pts": place([(-sx / 2 - ddy, sy / 2 + ddx), (sx / 2 + ddy, sy / 2 - ddx),
                                         (sx / 2 - ddy, -sy / 2 + ddx), (-sx / 2 + ddy, -sy / 2 - ddx)])}
    if shape == "rect":
        return {"r": 0.0, "pts": place([(-sx / 2, -sy / 2), (sx / 2, -sy / 2), (sx / 2, sy / 2), (-sx / 2, sy / 2)])}
    return None
