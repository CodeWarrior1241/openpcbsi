"""Geometry helpers: even-odd polygon rasterisation, capsules, segment distances, polylines."""
import math

import numpy as np


class Grid:
    """A raster over the board: cell (i, j) has its centre at (x0 + (j+0.5) res, y0 + (i+0.5) res)."""

    def __init__(self, x0, y0, x1, y1, res):
        self.res = float(res)
        self.x0 = math.floor(x0 / res) * res
        self.y0 = math.floor(y0 / res) * res
        self.nx = int(math.ceil((x1 - self.x0) / res)) + 1
        self.ny = int(math.ceil((y1 - self.y0) / res)) + 1

    @property
    def shape(self):
        return (self.ny, self.nx)

    def ij(self, x, y):
        j = np.floor((np.asarray(x) - self.x0) / self.res).astype(np.int64)
        i = np.floor((np.asarray(y) - self.y0) / self.res).astype(np.int64)
        return i, j

    def xy(self, i, j):
        return self.x0 + (np.asarray(j) + 0.5) * self.res, self.y0 + (np.asarray(i) + 0.5) * self.res

    def inside(self, i, j):
        return (i >= 0) & (i < self.ny) & (j >= 0) & (j < self.nx)

    def sample(self, arr, x, y, fill=0):
        i, j = self.ij(x, y)
        ok = self.inside(i, j)
        out = np.full(np.shape(i), fill, dtype=arr.dtype)
        out[ok] = arr[i[ok], j[ok]]
        return out

    def extent(self):
        """matplotlib imshow extent (y down)."""
        return (self.x0, self.x0 + self.nx * self.res, self.y0 + self.ny * self.res, self.y0)


def raster_rings(grid, rings, out=None, value=True):
    """Even-odd fill of a set of rings (an outline with holes, or fractured polygons).

    Works in the bounding-box window of the rings, so small polygons cost little."""
    if out is None:
        out = np.zeros(grid.shape, dtype=bool)
    pts = [np.asarray(r, dtype=float) for r in rings if len(r) >= 3]
    if not pts:
        return out
    allp = np.vstack(pts)
    i0, j0 = grid.ij(allp[:, 0].min(), allp[:, 1].min())
    i1, j1 = grid.ij(allp[:, 0].max(), allp[:, 1].max())
    i0, j0 = max(int(i0) - 1, 0), max(int(j0) - 1, 0)
    i1, j1 = min(int(i1) + 2, grid.ny), min(int(j1) + 2, grid.nx)
    if i0 >= i1 or j0 >= j1:
        return out
    sub = Grid.__new__(Grid)
    sub.res = grid.res
    sub.x0 = grid.x0 + j0 * grid.res
    sub.y0 = grid.y0 + i0 * grid.res
    sub.nx, sub.ny = j1 - j0, i1 - i0
    m = _raster_rings_full(sub, pts)
    win = out[i0:i1, j0:j1]
    if win.dtype == bool:
        win |= m
    else:
        win[m] = value
    return out


def _raster_rings_full(grid, rings):
    ny, nx = grid.shape
    out = np.zeros((ny, nx), dtype=bool)
    segs = []
    for r in rings:
        a = np.asarray(r, dtype=float)
        if len(a) < 3:
            continue
        b = np.roll(a, -1, axis=0)
        segs.append(np.hstack([a, b]))
    if not segs:
        return out
    s = np.vstack(segs)
    x1, y1, x2, y2 = s.T
    # rows whose centre lies in [min(y), max(y))
    fy1 = (y1 - grid.y0) / grid.res - 0.5
    fy2 = (y2 - grid.y0) / grid.res - 0.5
    lo = np.ceil(np.minimum(fy1, fy2)).astype(np.int64)
    hi = np.ceil(np.maximum(fy1, fy2)).astype(np.int64)  # exclusive
    lo = np.clip(lo, 0, ny)
    hi = np.clip(hi, 0, ny)
    cnt = hi - lo
    keep = cnt > 0
    if not keep.any():
        return out
    x1, y1, x2, y2, lo, cnt, fy1, fy2 = x1[keep], y1[keep], x2[keep], y2[keep], lo[keep], cnt[keep], fy1[keep], fy2[keep]
    rep = np.repeat(np.arange(len(lo)), cnt)
    rows = lo[rep] + (np.arange(rep.size) - np.repeat(np.cumsum(cnt) - cnt, cnt))
    t = (rows - fy1[rep]) / (fy2[rep] - fy1[rep])
    xc = x1[rep] + t * (x2[rep] - x1[rep])
    col = np.ceil((xc - grid.x0) / grid.res - 0.5).astype(np.int64)  # first cell centre >= crossing
    col = np.clip(col, 0, nx)
    order = np.lexsort((col, rows))
    rows, col = rows[order], col[order]
    # pair consecutive crossings in each row (even-odd)
    n = rows.size
    if n % 2:
        # numerical oddity; drop the last crossing of an odd row
        pass
    # group by row and pair
    starts = np.flatnonzero(np.r_[True, rows[1:] != rows[:-1]])
    lens = np.diff(np.r_[starts, n])
    pos_in_row = np.arange(n) - np.repeat(starts, lens)
    even = pos_in_row % 2 == 0
    # each even crossing pairs with the next one in the same row
    nxt = np.arange(n) + 1
    valid = even & (nxt < n)
    valid[valid] &= rows[nxt[valid]] == rows[valid]
    r_ = rows[valid]
    c0 = col[valid]
    c1 = col[nxt[valid]]
    diff = np.zeros((ny, nx + 1), dtype=np.int32)
    np.add.at(diff, (r_, c0), 1)
    np.add.at(diff, (r_, c1), -1)
    return np.cumsum(diff, axis=1)[:, :nx] > 0


def raster_capsule(grid, out, x1, y1, x2, y2, r, value=True):
    """Fill a stadium (track) of half-width r into out."""
    res = grid.res
    i0, j0 = grid.ij(min(x1, x2) - r, min(y1, y2) - r)
    i1, j1 = grid.ij(max(x1, x2) + r, max(y1, y2) + r)
    i0, j0 = max(int(i0), 0), max(int(j0), 0)
    i1, j1 = min(int(i1) + 1, grid.ny), min(int(j1) + 1, grid.nx)
    if i0 >= i1 or j0 >= j1:
        return
    xs = grid.x0 + (np.arange(j0, j1) + 0.5) * res
    ys = grid.y0 + (np.arange(i0, i1) + 0.5) * res
    X, Y = np.meshgrid(xs, ys)
    d = point_seg_dist(X, Y, x1, y1, x2, y2)
    m = d <= r
    sub = out[i0:i1, j0:j1]
    if sub.dtype == bool:
        sub |= m
    else:
        sub[m] = value


def point_seg_dist(px, py, x1, y1, x2, y2):
    dx, dy = x2 - x1, y2 - y1
    L2 = dx * dx + dy * dy
    if np.isscalar(L2) and L2 == 0:
        return np.hypot(px - x1, py - y1)
    t = ((px - x1) * dx + (py - y1) * dy) / np.where(L2 == 0, 1, L2)
    t = np.clip(t, 0, 1)
    return np.hypot(px - (x1 + t * dx), py - (y1 + t * dy))


def seg_seg_dist(a, b):
    """Minimum distance between segments a=(x1,y1,x2,y2) and b."""
    ax1, ay1, ax2, ay2 = a
    bx1, by1, bx2, by2 = b
    if _seg_intersect(a, b):
        return 0.0
    return min(point_seg_dist(ax1, ay1, bx1, by1, bx2, by2), point_seg_dist(ax2, ay2, bx1, by1, bx2, by2),
               point_seg_dist(bx1, by1, ax1, ay1, ax2, ay2), point_seg_dist(bx2, by2, ax1, ay1, ax2, ay2))


def _seg_intersect(a, b):
    def orient(px, py, qx, qy, rx, ry):
        return (qx - px) * (ry - py) - (qy - py) * (rx - px)
    ax1, ay1, ax2, ay2 = a
    bx1, by1, bx2, by2 = b
    d1 = orient(bx1, by1, bx2, by2, ax1, ay1)
    d2 = orient(bx1, by1, bx2, by2, ax2, ay2)
    d3 = orient(ax1, ay1, ax2, ay2, bx1, by1)
    d4 = orient(ax1, ay1, ax2, ay2, bx2, by2)
    return (d1 * d2 < 0) and (d3 * d4 < 0)


def sample_segment(x1, y1, x2, y2, step):
    """Points along a segment at about `step` spacing (cell-centred), with the length per point."""
    L = math.hypot(x2 - x1, y2 - y1)
    n = max(1, int(math.ceil(L / step)))
    t = (np.arange(n) + 0.5) / n
    return x1 + t * (x2 - x1), y1 + t * (y2 - y1), L / n


def ring_area(r):
    a = np.asarray(r)
    x, y = a[:, 0], a[:, 1]
    return 0.5 * float(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))


def points_in_rings(px, py, rings):
    """Even-odd point-in-polygon for arrays of points."""
    px = np.asarray(px, dtype=float)
    py = np.asarray(py, dtype=float)
    inside = np.zeros(px.shape, dtype=bool)
    for r in rings:
        a = np.asarray(r, dtype=float)
        x1, y1 = a[:, 0], a[:, 1]
        x2, y2 = np.roll(x1, -1), np.roll(y1, -1)
        for k in range(len(a)):
            cond = (y1[k] > py) != (y2[k] > py)
            if not cond.any():
                continue
            xc = x1[k] + (py - y1[k]) * (x2[k] - x1[k]) / ((y2[k] - y1[k]) or 1e-30)
            inside ^= cond & (px < xc)
    return inside


def polyline_length(pts):
    a = np.asarray(pts)
    return float(np.sum(np.hypot(np.diff(a[:, 0]), np.diff(a[:, 1]))))


def resample_ring(ring, step):
    """Points every `step` mm along a closed ring, with their arc-length position."""
    a = np.asarray(ring, dtype=float)
    a = np.vstack([a, a[:1]])
    d = np.hypot(np.diff(a[:, 0]), np.diff(a[:, 1]))
    s = np.r_[0, np.cumsum(d)]
    L = s[-1]
    n = max(4, int(L / step))
    ss = np.linspace(0, L, n, endpoint=False)
    x = np.interp(ss, s, a[:, 0])
    y = np.interp(ss, s, a[:, 1])
    return x, y, ss, L


def cluster(xs, ys, tol=0.005):
    """Labels that join points closer than tol (union-find over a KD-tree)."""
    from scipy.spatial import cKDTree
    pts = np.c_[np.asarray(xs, float), np.asarray(ys, float)]
    n = len(pts)
    parent = np.arange(n)

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a
    if n:
        for i, j in cKDTree(pts).query_pairs(tol):
            ri, rj = find(i), find(j)
            if ri != rj:
                parent[ri] = rj
    return np.array([find(i) for i in range(n)])
