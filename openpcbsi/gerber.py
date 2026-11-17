"""Readers for the files a fabricator receives: Gerber (RS-274X / X2) images and Excellon drill files.

Used by the `output` check of `openpcbsi dfm` (MF20) to compare the files that are sent with the
geometry the other checks worked on. The Gerber reader covers what KiCad plots: standard apertures
C, R, O and P, aperture macros with primitives 1, 4, 5, 20 and 21 and arithmetic in their
parameters, linear and circular interpolation (multi-quadrant), regions, flashes and both polarities.
Anything else (apertures with holes, step-and-repeat, single-quadrant arcs, macro primitives it does
not know) is listed in `unsupported` and the image is then not trusted: the caller reports it instead
of comparing a partial image.
"""
import ast
import math
import operator
import os
import re

import numpy as np

from . import geom

_OPS = {ast.Add: operator.add, ast.Sub: operator.sub, ast.Mult: operator.mul, ast.Div: operator.truediv}


def _calc(expr, params):
    """A macro expression: numbers, $n, + - x / and brackets."""
    s = re.sub(r"\$(\d+)", lambda m: repr(float(params.get(int(m.group(1)), 0.0))), expr.strip())
    s = s.replace("x", "*").replace("X", "*")

    def ev(n):
        if isinstance(n, ast.Constant) and isinstance(n.value, (int, float)):
            return float(n.value)
        if isinstance(n, ast.BinOp) and type(n.op) in _OPS:
            return _OPS[type(n.op)](ev(n.left), ev(n.right))
        if isinstance(n, ast.UnaryOp) and isinstance(n.op, (ast.USub, ast.UAdd)):
            v = ev(n.operand)
            return -v if isinstance(n.op, ast.USub) else v
        raise ValueError(f"expression {expr!r}")
    return ev(ast.parse(s, mode="eval").body)


def _rot(pts, deg):
    if not deg:
        return pts
    c, s = math.cos(math.radians(deg)), math.sin(math.radians(deg))
    return [(x * c - y * s, x * s + y * c) for x, y in pts]


def _circle(cx, cy, r, n=48):
    return [(cx + r * math.cos(2 * math.pi * k / n), cy + r * math.sin(2 * math.pi * k / n)) for k in range(n)]


class Aperture:
    """Shapes of one aperture in its own coordinates: polygons (point lists) and circles (cx, cy, r)."""

    def __init__(self):
        self.polys, self.circles, self.diameter = [], [], None


def _standard(kind, vals, bad):
    a = Aperture()
    if kind == "C":
        if len(vals) > 1:
            bad.add("aperture with a hole")
        a.circles.append((0.0, 0.0, vals[0] / 2))
        a.diameter = vals[0]
    elif kind == "R":
        if len(vals) > 2:
            bad.add("aperture with a hole")
        w, h = vals[0] / 2, vals[1] / 2
        a.polys.append([(-w, -h), (w, -h), (w, h), (-w, h)])
    elif kind == "O":
        if len(vals) > 2:
            bad.add("aperture with a hole")
        w, h = vals[0], vals[1]
        r = min(w, h) / 2
        d = (max(w, h) - min(w, h)) / 2
        ends = ((-d, 0.0), (d, 0.0)) if w >= h else ((0.0, -d), (0.0, d))
        a.circles += [(ends[0][0], ends[0][1], r), (ends[1][0], ends[1][1], r)]
        a.polys.append([(-w / 2 + (r if w >= h else 0), -h / 2 + (0 if w >= h else r)),
                        (w / 2 - (r if w >= h else 0), -h / 2 + (0 if w >= h else r)),
                        (w / 2 - (r if w >= h else 0), h / 2 - (0 if w >= h else r)),
                        (-w / 2 + (r if w >= h else 0), h / 2 - (0 if w >= h else r))])
    elif kind == "P":
        if len(vals) > 3:
            bad.add("aperture with a hole")
        n = int(vals[1])
        rot = vals[2] if len(vals) > 2 else 0.0
        a.polys.append([(vals[0] / 2 * math.cos(math.radians(rot) + 2 * math.pi * k / n),
                         vals[0] / 2 * math.sin(math.radians(rot) + 2 * math.pi * k / n)) for k in range(n)])
    else:
        bad.add(f"aperture type {kind}")
    return a


def _macro(body, vals, bad):
    a = Aperture()
    params = {k + 1: v for k, v in enumerate(vals)}
    for st in body:
        st = st.strip()
        if not st or st.startswith("0 ") or st == "0":
            continue
        if st.startswith("$") and "=" in st:
            k, e = st.split("=", 1)
            params[int(k[1:])] = _calc(e, params)
            continue
        try:
            f = [_calc(x, params) for x in st.split(",")]
        except Exception:
            bad.add("macro expression")
            continue
        code = int(f[0])
        if code in (1, 4, 5, 20, 21) and f[1] != 1:
            bad.add("macro primitive with clear exposure")
            continue
        if code == 1:
            rot = f[5] if len(f) > 5 else 0.0
            (cx, cy), = _rot([(f[3], f[4])], rot)
            a.circles.append((cx, cy, f[2] / 2))
        elif code == 4:
            n = int(f[2])
            pts = [(f[3 + 2 * k], f[4 + 2 * k]) for k in range(n + 1)]
            a.polys.append(_rot(pts[:-1] if pts[0] == pts[-1] else pts, f[5 + 2 * n] if len(f) > 5 + 2 * n else 0.0))
        elif code == 5:
            n = int(f[2])
            a.polys.append(_rot([(f[3] + f[5] / 2 * math.cos(2 * math.pi * k / n), f[4] + f[5] / 2 * math.sin(2 * math.pi * k / n))
                                 for k in range(n)], f[6] if len(f) > 6 else 0.0))
        elif code == 20:
            w, (xs, ys, xe, ye) = f[2], f[3:7]
            L = math.hypot(xe - xs, ye - ys) or 1e-12
            nx, ny = -(ye - ys) / L * w / 2, (xe - xs) / L * w / 2
            a.polys.append(_rot([(xs + nx, ys + ny), (xe + nx, ye + ny), (xe - nx, ye - ny), (xs - nx, ys - ny)],
                                f[7] if len(f) > 7 else 0.0))
        elif code == 21:
            w, h, cx, cy = f[2] / 2, f[3] / 2, f[4], f[5]
            a.polys.append(_rot([(cx - w, cy - h), (cx + w, cy - h), (cx + w, cy + h), (cx - w, cy + h)],
                                f[6] if len(f) > 6 else 0.0))
        else:
            bad.add(f"macro primitive {code}")
    return a


def read(path):
    """A Gerber file as drawing operations in mm, Gerber axes (y up).

    Returns {"function", "polarity", "ops", "unsupported"}; ops are (dark, "poly", points) and
    (dark, "capsule", x1, y1, x2, y2, r) in drawing order."""
    text = open(path, encoding="latin-1").read()
    out = {"function": None, "polarity": "Positive", "ops": [], "unsupported": set(), "path": path,
           "flashes": [], "draws": []}     # flashes: (x, y, w, h, area, dark); draws: (x0, y0, x1, y1, r, dark), Gerber axes
    bad = out["unsupported"]
    macros, aps = {}, {}
    scale, dec = 1.0, 6
    x = y = 0.0
    mode, region, dark, cur = 1, False, True, None
    contour = []

    def num(s):
        return int(s) / 10 ** dec * scale

    def close():
        if len(contour) >= 3:
            out["ops"].append((dark, "poly", list(contour)))
        contour.clear()

    def flash(px, py):
        a = aps.get(cur)
        if a is None:
            bad.add("flash without an aperture")
            return
        xs, ys, area = [], [], 0.0
        for p in a.polys:
            out["ops"].append((dark, "poly", [(px + u, py + v) for u, v in p]))
            xs += [u for u, v in p]
            ys += [v for u, v in p]
            area += abs(0.5 * sum(p[k][0] * p[k - 1][1] - p[k - 1][0] * p[k][1] for k in range(len(p))))
        for cx, cy, r in a.circles:
            out["ops"].append((dark, "capsule", px + cx, py + cy, px + cx, py + cy, r))
            xs += [cx - r, cx + r]
            ys += [cy - r, cy + r]
            area += math.pi * r * r
        if xs:
            out["flashes"].append((px, py, max(xs) - min(xs), max(ys) - min(ys), area, dark))

    def draw(x0, y0, x1, y1):
        if region:
            if not contour:
                contour.append((x0, y0))
            contour.append((x1, y1))
            return
        a = aps.get(cur)
        if a is None or a.diameter is None:
            bad.add("draw with a non-circular aperture")
            return
        out["ops"].append((dark, "capsule", x0, y0, x1, y1, a.diameter / 2))
        out["draws"].append((x0, y0, x1, y1, a.diameter / 2, dark))

    i, n = 0, len(text)
    while i < n:
        c = text[i]
        if c in " \r\n\t":
            i += 1
            continue
        if c == "%":
            j = text.find("%", i + 1)
            if j < 0:
                break
            sts = [s.strip() for s in text[i + 1:j].split("*") if s.strip()]
            i = j + 1
            if not sts:
                continue
            h = sts[0]
            if h.startswith("FS"):
                m = re.match(r"FS([LT])([AI])X(\d)(\d)Y\d\d", h)
                if not m or m.group(1) != "L" or m.group(2) != "A":
                    bad.add("coordinate format " + h)
                else:
                    dec = int(m.group(4))
            elif h.startswith("MO"):
                scale = 25.4 if h[2:4] == "IN" else 1.0
            elif h.startswith("AM"):
                macros[h[2:]] = sts[1:]
            elif h.startswith("AD"):
                m = re.match(r"ADD(\d+)([A-Za-z_.$][\w.$]*)(?:,(.*))?$", h)
                if not m:
                    bad.add("aperture definition")
                    continue
                vals = [float(v) for v in m.group(3).split("X")] if m.group(3) else []
                name = m.group(2)
                if name in macros:
                    ap = _macro(macros[name], vals, bad)
                else:
                    ap = _standard(name, vals, bad)
                for shape in (ap.polys, ap.circles):
                    pass
                if scale != 1.0:
                    ap.polys = [[(u * scale, v * scale) for u, v in p] for p in ap.polys]
                    ap.circles = [(u * scale, v * scale, r * scale) for u, v, r in ap.circles]
                    ap.diameter = ap.diameter * scale if ap.diameter else None
                aps[int(m.group(1))] = ap
            elif h.startswith("LP"):
                dark = h[2:3] != "C"
            elif h.startswith("TF.FileFunction"):
                out["function"] = h.split(",", 1)[1] if "," in h else ""
            elif h.startswith("TF.FilePolarity"):
                out["polarity"] = h.split(",", 1)[1]
            elif h.startswith("SR"):
                if h not in ("SR", "SRX1Y1I0J0", "SRX1Y1I0.0J0.0"):
                    bad.add("step and repeat")
            elif h[:2] in ("TF", "TA", "TO", "TD", "IN", "LN", "IP", "OF", "AS", "MI", "SF", "IR"):
                if h[:2] in ("MI", "SF", "IR", "OF", "AS") and h not in ("IPPOS", "OFA0B0", "MIA0B0", "SFA1.0B1.0", "ASAXBY", "IR0"):
                    bad.add("deprecated image transform " + h[:2])
                if h == "IPNEG":
                    out["polarity"] = "Negative"
            else:
                bad.add("command %" + h[:2])
            continue
        j = text.find("*", i)
        if j < 0:
            break
        st = text[i:j].strip()
        i = j + 1
        if not st or st.startswith("G04"):
            continue
        while True:
            m = re.match(r"G(\d+)", st)
            if not m:
                break
            g = int(m.group(1))
            st = st[m.end():]
            if g in (1, 2, 3):
                mode = g
            elif g == 36:
                region = True
                contour.clear()
            elif g == 37:
                close()
                region = False
            elif g == 74:
                bad.add("single-quadrant arcs")
            elif g in (75, 90, 71, 54, 55):
                pass
            elif g == 70:
                scale = 25.4
            elif g == 91:
                bad.add("incremental coordinates")
            else:
                bad.add(f"G{g:02d}")
        if not st:
            continue
        if st == "M02" or st == "M00":
            break
        m = re.match(r"D(\d+)$", st)
        if m and int(m.group(1)) >= 10:
            cur = int(m.group(1))
            continue
        m = re.match(r"(?:X(-?\d+))?(?:Y(-?\d+))?(?:I(-?\d+))?(?:J(-?\d+))?(?:D0?([123]))?$", st)
        if not m:
            bad.add("statement " + st[:12])
            continue
        nx = num(m.group(1)) if m.group(1) is not None else x
        ny = num(m.group(2)) if m.group(2) is not None else y
        op = int(m.group(5)) if m.group(5) else 1
        if op == 2:
            if region:
                close()
        elif op == 3:
            flash(nx, ny)
        elif mode == 1:
            draw(x, y, nx, ny)
        else:
            ci = num(m.group(3)) if m.group(3) is not None else 0.0
            cj = num(m.group(4)) if m.group(4) is not None else 0.0
            cx, cy = x + ci, y + cj
            r = math.hypot(ci, cj)
            a0, a1 = math.atan2(y - cy, x - cx), math.atan2(ny - cy, nx - cx)
            if mode == 2:                       # clockwise (y up)
                sweep = -((a0 - a1) % (2 * math.pi)) or -2 * math.pi
            else:
                sweep = ((a1 - a0) % (2 * math.pi)) or 2 * math.pi
            k = max(2, int(abs(sweep) / math.radians(5)) + 1)
            px, py = x, y
            for s in range(1, k + 1):
                qx, qy = (nx, ny) if s == k else (cx + r * math.cos(a0 + sweep * s / k), cy + r * math.sin(a0 + sweep * s / k))
                draw(px, py, qx, qy)
                px, py = qx, qy
        x, y = nx, ny
    out["unsupported"] = sorted(bad)
    return out


def bbox(g):
    xs, ys = [], []
    for o in g["ops"]:
        if o[1] == "poly":
            xs += [p[0] for p in o[2]]
            ys += [p[1] for p in o[2]]
        else:
            xs += [o[2] - o[6], o[4] + o[6], o[2] + o[6], o[4] - o[6]]
            ys += [o[3] - o[6], o[5] + o[6], o[3] + o[6], o[5] - o[6]]
    return (min(xs), min(ys), max(xs), max(ys)) if xs else None


def raster(g, grid, dx=0.0, dy=0.0):
    """Image of a Gerber file on a board grid (board y = -Gerber y, then the offset)."""
    out = np.zeros(grid.shape, dtype=bool)
    for o in g["ops"]:
        if o[1] == "poly":
            ring = [(px + dx, -py + dy) for px, py in o[2]]
            if o[0]:
                geom.raster_rings(grid, [ring], out)
            else:
                out &= ~geom.raster_rings(grid, [ring])
        else:
            _, _, x1, y1, x2, y2, r = o
            if o[0]:
                geom.raster_capsule(grid, out, x1 + dx, -y1 + dy, x2 + dx, -y2 + dy, r)
            else:
                tmp = np.zeros(grid.shape, dtype=bool)
                geom.raster_capsule(grid, tmp, x1 + dx, -y1 + dy, x2 + dx, -y2 + dy, r)
                out &= ~tmp
    return out


# ------------------------------------------------------------------------------------------------
LAYER_HINTS = (("F_Cu", "F.Cu"), ("B_Cu", "B.Cu"), ("F_Mask", "F.Mask"), ("B_Mask", "B.Mask"),
               ("F_Paste", "F.Paste"), ("B_Paste", "B.Paste"), ("Edge_Cuts", "Edge.Cuts"))
PROTEL = {".gtl": "F.Cu", ".gbl": "B.Cu", ".gts": "F.Mask", ".gbs": "B.Mask", ".gtp": "F.Paste", ".gbp": "B.Paste",
          ".gm1": "Edge.Cuts", ".gko": "Edge.Cuts"}


def layer_of(path, g, copper_names):
    """The board layer a Gerber file is for: from its X2 FileFunction, else from its name."""
    f = (g.get("function") or "").split(",")
    if f and f[0] == "Copper" and len(f) > 1 and re.match(r"L\d+$", f[1]):
        k = int(f[1][1:]) - 1
        return copper_names[k] if 0 <= k < len(copper_names) else None
    if f and f[0] in ("Soldermask", "Paste") and len(f) > 1:
        return ("F." if f[1] == "Top" else "B.") + ("Mask" if f[0] == "Soldermask" else "Paste")
    if f and f[0] == "Profile":
        return "Edge.Cuts"
    name = os.path.basename(path)
    m = re.search(r"In(\d+)[_.]Cu", name)
    if m:
        return f"In{m.group(1)}.Cu"
    for hint, layer in LAYER_HINTS:
        if hint in name or hint.replace("_", ".") in name:
            return layer
    ext = os.path.splitext(name)[1].lower()
    m = re.match(r"\.g(\d+)$", ext)
    if m and int(m.group(1)) < len(copper_names) - 1:
        return copper_names[int(m.group(1)) - 1 + 1] if 1 <= int(m.group(1)) <= len(copper_names) - 2 else None
    return PROTEL.get(ext)


def read_drill(path):
    """Holes of an Excellon file: [{x, y, d, x2, y2, plated}] in mm, file axes (y up). A slot has x2, y2."""
    tools, holes, bad = {}, [], set()
    scale, cur, header = 1.0, None, True
    plated = None
    lz = None            # zero suppression of integer coordinates: None = decimal points expected
    route, down, rx, ry = False, False, 0.0, 0.0
    for line in open(path, encoding="latin-1"):
        s = line.strip()
        if not s:
            continue
        if s.startswith(";"):
            if "NonPlated" in s or "NPTH" in s:
                plated = False
            elif "Plated" in s and plated is None:
                plated = True
            continue
        if header:
            if s in ("%", "M95"):
                header = False
            elif s.startswith("METRIC") or s.startswith("M71"):
                scale = 1.0
                lz = "LZ" if ",LZ" in s else "TZ" if ",TZ" in s else lz
            elif s.startswith("INCH") or s.startswith("M72"):
                scale = 25.4
                lz = "LZ" if ",LZ" in s else "TZ" if ",TZ" in s else lz
            else:
                m = re.match(r"T(\d+).*?C([\d.]+)", s)
                if m:
                    tools[int(m.group(1))] = float(m.group(2)) * scale
            continue
        m = re.match(r"T(\d+)(?:C([\d.]+))?$", s)
        if m:
            cur = int(m.group(1))
            if m.group(2):
                tools[cur] = float(m.group(2)) * scale
            continue
        if s in ("G90", "G05", "M30", "M16", "M17", "G40", "G41", "G42") or s.startswith("M48") or s.startswith("FMAT"):
            if s in ("M16", "M17"):
                down = False
            continue
        if s == "M15":
            down = True
            continue

        def xy(t):
            mm = re.match(r"(?:X(-?[\d.]+))?(?:Y(-?[\d.]+))?", t)
            vals = []
            for v in (mm.group(1), mm.group(2)):
                if v is None:
                    vals.append(None)
                elif "." in v or lz is None:
                    vals.append(float(v) * scale if "." in v else float(v) * scale)
                else:
                    bad.add("coordinates without a decimal point")
                    vals.append(None)
            return vals
        if s.startswith("G00"):
            px, py = xy(s[3:])
            rx, ry = px if px is not None else rx, py if py is not None else ry
            route = True
            continue
        if s.startswith("G01") and route:
            px, py = xy(s[3:])
            px, py = px if px is not None else rx, py if py is not None else ry
            if down and cur in tools:
                holes.append({"x": rx, "y": ry, "x2": px, "y2": py, "d": tools[cur], "plated": plated})
            rx, ry = px, py
            continue
        if "G85" in s:
            a, b_ = s.split("G85")
            x1, y1 = xy(a)
            x2, y2 = xy(b_)
            if None in (x1, y1, x2, y2) or cur not in tools:
                bad.add("slot statement")
                continue
            holes.append({"x": x1, "y": y1, "x2": x2, "y2": y2, "d": tools[cur], "plated": plated})
            continue
        if s[0] in "XY":
            px, py = xy(s)
            if cur not in tools or px is None or py is None:
                bad.add("hole statement")
                continue
            holes.append({"x": px, "y": py, "d": tools[cur], "plated": plated})
            continue
        bad.add("statement " + s[:10])
    return {"holes": holes, "unsupported": sorted(bad), "path": path}
