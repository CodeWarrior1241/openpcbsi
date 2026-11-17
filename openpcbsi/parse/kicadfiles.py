"""The board file read with KiCadFiles: the independent cross-check of the primary reader (plan phase 2).

KiCadFiles 1.0.5 reads in lenient mode only and leaves several fields empty or at defaults without
saying so (pad and via layer lists, board graphics, custom-pad primitives, multi-layer zone layers,
DNP); it was released before KiCad 10. The cross-check therefore compares only what it does read:

    CROSS_CHECKED = nets, track and arc geometry, via position / size / drill / net, footprint reference,
    value, library id, position, rotation and side, pad number / type / net / size (the size written in
    the file: a custom pad's anchor) / drill / position, and zone fill area per net and layer.

Everything else is listed as not cross-checked in the review.
"""
import math

from . import use

CROSS_CHECKED = ("nets", "tracks", "vias", "zones", "footprints", "pads")
NOT_CROSS_CHECKED = ("pad and via copper layers", "custom-pad primitives", "mask and paste openings", "board outline",
                     "stackup (KiCadFiles drops dielectric plies)", "footprint attributes (DNP, BOM, position files)",
                     "rule areas")


NET_RE = None


def _blocks(text, start, end, head):
    """Top-level `(head ...)` blocks between start and end, in file order: (from, to) offsets."""
    out, i = [], start
    key = "(" + head
    while True:
        i = text.find(key, i, end)
        if i < 0:
            return out
        if text[i + len(key)] not in " \n\t)":
            i += 1
            continue
        depth, j = 0, i
        while j < end:
            c = text[j]
            if c == '"':
                j = text.index('"', j + 1)
                while text[j - 1] == "\\":
                    j = text.index('"', j + 1)
            elif c == "(":
                depth += 1
            elif c == ")":
                depth -= 1
                if depth == 0:
                    break
            j += 1
        out.append((i, j + 1))
        i = j + 1


def raw_nets(text):
    """Net names of tracks, arcs, vias and pads in file order, read from the file text by openpcbsi:
    KiCadFiles turns KiCad 10's name-only nets into net 0. Independent of the primary reader."""
    import re
    rx = re.compile(r'\(net (?:\d+ )?"((?:[^"\\]|\\.)*)"\)')

    def net_of(a, b):
        # the item's own (net ...), not one of a nested item
        depth, k = 0, a + 1
        while k < b:
            c = text[k]
            if c == '"':
                k = text.index('"', k + 1)
            elif c == "(":
                depth += 1
                if depth == 1:
                    m = rx.match(text, k)
                    if m:
                        return m.group(1)
            elif c == ")":
                depth -= 1
            k += 1
        return ""
    root = _blocks(text, 0, len(text), "kicad_pcb")[0]
    tops = {h: _blocks(text, root[0] + 1, root[1] - 1, h) for h in ("segment", "arc", "via", "footprint", "zone")}
    # zones defined inside footprints are not the board's (KiCadFiles reads only the board's)
    import bisect
    starts = [a for a, _ in tops["footprint"]]
    def in_fp(pos):
        k = bisect.bisect_right(starts, pos) - 1
        return k >= 0 and tops["footprint"][k][0] < pos < tops["footprint"][k][1]
    tops["zone"] = [z for z in tops["zone"] if not in_fp(z[0])]
    res = {h: [net_of(a, b) for a, b in tops[h]] for h in ("segment", "arc", "via", "zone")}
    res["pads"] = [[net_of(a, b) for a, b in _blocks(text, fa + 1, fb - 1, "pad")] for fa, fb in tops["footprint"]]
    return res


def _v(x, default=None):
    return getattr(x, "value", x) if x is not None else default


def extract(board_path):
    use()
    from kicadfiles import KicadPcb, ParseStrictness
    from ..snapshot import _arc_points
    pcb = KicadPcb.from_file(board_path, ParseStrictness.SILENT)
    names = {n.number: n.name for n in (pcb.nets or [])}
    rn = raw_nets(open(board_path, encoding="utf-8").read()) if not names else None

    def net(x):
        v = _v(x)
        if isinstance(v, str):
            return v
        return names.get(v, "") if v is not None else ""
    tracks = []
    pick = lambda kind, k, x: rn[kind][k] if rn and k < len(rn[kind]) else net(x)
    for k, s in enumerate(pcb.segments or []):
        tracks.append({"x1": s.start.x, "y1": s.start.y, "x2": s.end.x, "y2": s.end.y, "w": _v(s.width), "layer": s.layer.name,
                       "net": pick("segment", k, s.net)})
    for ka, a in enumerate(pcb.arcs or []):
        pts = _arc_points((a.start.x, a.start.y), (a.mid.x, a.mid.y), (a.end.x, a.end.y))
        for (x1, y1), (x2, y2) in zip(pts[:-1], pts[1:]):
            tracks.append({"x1": x1, "y1": y1, "x2": x2, "y2": y2, "w": _v(a.width), "layer": _v(a.layer), "net": pick("arc", ka, a.net), "arc": True})
    vias = [{"x": v.at.x, "y": v.at.y, "d": _v(v.size), "drill": _v(v.drill), "net": pick("via", k, v.net)} for k, v in enumerate(pcb.vias or [])]
    zones = []
    for kz, z in enumerate(pcb.zones or []):
        if "keepout" in getattr(z, "_parsed_fields", ()):
            continue
        zn = rn["zone"][kz] if rn and kz < len(rn["zone"]) else (z.net_name.value if z.net_name else "")
        for f in z.filled_polygons or []:
            pts = [(p.x, p.y) for p in f.pts.points]
            if len(pts) >= 3 and str(_v(f.layer)).endswith(".Cu"):          # copper fills only, as the other readers
                zones.append({"net": zn, "layer": _v(f.layer), "polys": [[pts]]})
    fps, pads = [], []
    for kf, fp in enumerate(pcb.footprints or []):
        props = {p.key: p.value for p in fp.properties or []}
        pnets = rn["pads"][kf] if rn and kf < len(rn["pads"]) else None
        ref = props.get("Reference", "")
        a = math.radians(-(fp.at.angle or 0.0))
        c, s_ = math.cos(a), math.sin(a)
        fps.append({"ref": ref, "value": props.get("Value", ""), "fpid": fp.library_link or "", "x": fp.at.x, "y": fp.at.y,
                    "rot": float(fp.at.angle or 0.0), "side": "B" if fp.layer.name == "B.Cu" else "F"})
        for kp, p in enumerate(fp.pads or []):
            kind = {"thru_hole": "pth", "smd": "smd", "connect": "conn", "np_thru_hole": "npth"}.get(str(_v(p.type)), "other")
            x, y = fp.at.x + p.at.x * c - p.at.y * s_, fp.at.y + p.at.x * s_ + p.at.y * c
            drill = dy = 0.0
            if kind in ("pth", "npth") and "drill" in getattr(p, "_parsed_fields", ()):
                drill = float(p.drill.diameter or 0.0)
                dy = drill
                if p.drill.oval is not None and p.drill.width:
                    dy = float(p.drill.width)
            sz = (float(p.size.width), float(p.size.height))
            pnet = pnets[kp] if pnets is not None and kp < len(pnets) else (p.net.name if p.net else "")
            pads.append({"ref": ref, "num": str(p.number or ""), "net": pnet, "kind": kind, "x": x, "y": y,
                         "rot": float(p.at.angle or 0.0), "shape": str(_v(p.shape)), "anchor": list(sz), "sx": sz[0], "sy": sz[1],
                         "drill": drill, "drill_y": dy})
    allnets = {n.name for n in pcb.nets or [] if n.name} | {t["net"] for t in tracks if t["net"]} | {v["net"] for v in vias if v["net"]} \
        | {p["net"] for p in pads if p["net"]} | {z["net"] for z in zones if z["net"]}
    return {"reader": "kicadfiles", "net_names_from": "file text (openpcbsi)" if rn else "kicadfiles",
            "nets": {n: {} for n in allnets}, "tracks": tracks, "vias": vias,
            "zones": zones, "footprints": fps, "pads": pads}


def cross_check(primary, board_path):
    """The primary reading against KiCadFiles' on the fields KiCadFiles reads: compare.differences."""
    from . import compare
    other = extract(board_path)
    # a custom pad's size in the file is its anchor: compare that
    prim = dict(primary, pads=[dict(p, sx=p.get("anchor", [p["sx"], p["sy"]])[0], sy=p.get("anchor", [p["sx"], p["sy"]])[1])
                                for p in primary["pads"]])
    d = compare.differences(prim, other, fields=CROSS_CHECKED)
    d["not_cross_checked"] = list(NOT_CROSS_CHECKED)
    return d
