"""Pads as KiCad's IPC-2581 export describes them: the independent reference for pad geometry.

    pads(xml_path) -> [{layer, x, y, ref, pin, w, h, area}]

Coordinates are KiCad's (the export's y is negated back). w and h are the pad's extent in its own
frame (the primitive before the pad's Xform rotation), area its area, both in mm.
"""
import math
import xml.etree.ElementTree as ET

NS = "{http://webstds.ipc.org/2581}"


def _f(e, k, d=0.0):
    v = e.get(k)
    return float(v) if v is not None else d


def _arc_pts(x0, y0, x1, y1, cx, cy, cw, max_deg=5.0):
    a0, a1 = math.atan2(y0 - cy, x0 - cx), math.atan2(y1 - cy, x1 - cx)
    r = math.hypot(x0 - cx, y0 - cy)
    sweep = a1 - a0
    if cw and sweep > 0:
        sweep -= 2 * math.pi
    if not cw and sweep < 0:
        sweep += 2 * math.pi
    if abs(sweep) < 1e-12:
        sweep = -2 * math.pi if cw else 2 * math.pi
    n = max(2, int(abs(math.degrees(sweep)) / max_deg) + 1)
    return [(cx + r * math.cos(a0 + sweep * k / n), cy + r * math.sin(a0 + sweep * k / n)) for k in range(1, n + 1)]


def _poly(pe):
    pts = []
    for c in pe:
        tag = c.tag.replace(NS, "")
        x, y = _f(c, "x"), _f(c, "y")
        if tag == "PolyBegin":
            pts = [(x, y)]
        elif tag == "PolyStepSegment":
            pts.append((x, y))
        elif tag == "PolyStepCurve":
            x0, y0 = pts[-1]
            pts += _arc_pts(x0, y0, x, y, _f(c, "centerX"), _f(c, "centerY"), c.get("clockwise", "false").lower() == "true")
    return pts


def _area(pts):
    return 0.5 * abs(sum(x0 * y1 - x1 * y0 for (x0, y0), (x1, y1) in zip(pts, pts[1:] + pts[:1])))


def _prim(e):
    """(w, h, area) of a standard primitive."""
    tag = e.tag.replace(NS, "")
    if tag == "RectCenter":
        w, h = _f(e, "width"), _f(e, "height")
        return w, h, w * h
    if tag == "Circle":
        d = _f(e, "diameter")
        return d, d, math.pi * d * d / 4
    if tag == "Oval":
        w, h = _f(e, "width"), _f(e, "height")
        s = min(w, h)
        return w, h, (max(w, h) - s) * s + math.pi * s * s / 4
    if tag == "RectRound":
        w, h, r = _f(e, "width"), _f(e, "height"), _f(e, "radius")
        return w, h, w * h - (4 - math.pi) * r * r
    if tag == "RectCham":
        w, h, c = _f(e, "width"), _f(e, "height"), _f(e, "chamfer")
        n = sum(e.get(k, "false").lower() == "true" for k in ("upperRight", "upperLeft", "lowerRight", "lowerLeft"))
        return w, h, w * h - n * c * c / 2
    if tag == "Contour":
        outline = e.find(NS + "Polygon")
        pts = _poly(outline) if outline is not None else []
        if not pts:
            return None
        xs, ys = [p[0] for p in pts], [p[1] for p in pts]
        a = _area(pts) - sum(_area(_poly(c)) for c in e.findall(NS + "Cutout"))
        return max(xs) - min(xs), max(ys) - min(ys), a
    return None


def pads(path):
    root = ET.parse(path).getroot()
    std = {}
    for es in root.iter(NS + "EntryStandard"):
        for c in es:
            v = _prim(c)
            if v:
                std[es.get("id")] = v
                break
    stacks = {}
    for psd in root.iter(NS + "PadStackDef"):
        per = {}
        for ppd in psd.findall(NS + "PadstackPadDef"):
            ref = ppd.find(NS + "StandardPrimitiveRef")
            if ref is not None and ref.get("id") in std:
                per[ppd.get("layerRef")] = std[ref.get("id")]
        stacks[psd.get("name")] = per
    out = []
    for lf in root.iter(NS + "LayerFeature"):
        layer = lf.get("layerRef")
        for pad in lf.iter(NS + "Pad"):
            loc = pad.find(NS + "Location")
            pin = pad.find(NS + "PinRef")
            prim = None
            sref = pad.find(NS + "StandardPrimitiveRef")
            if sref is not None:
                prim = std.get(sref.get("id"))
            if prim is None:
                prim = stacks.get(pad.get("padstackDefRef"), {}).get(layer)
            if loc is None or prim is None:
                continue
            out.append({"layer": layer, "x": _f(loc, "x"), "y": -_f(loc, "y"),
                        "ref": pin.get("componentRef") if pin is not None else None,
                        "pin": pin.get("pin") if pin is not None else None,
                        "w": prim[0], "h": prim[1], "area": prim[2]})
    return out
