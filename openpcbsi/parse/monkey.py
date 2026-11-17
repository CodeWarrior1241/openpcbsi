"""The board file read with kicad-monkey (the primary file reader), as openpcbsi's snapshot.

kicad-monkey parses the s-expressions into typed objects; what it does not model (dielectric plies,
custom-pad primitives other than polygons, the board outline as polygons) is taken from the raw
s-expressions it keeps, and every polygon is computed by openpcbsi (shapes.py). The result has the
same fields as the pcbnew snapshot (snapshot.extract), so every check reads either one unchanged.

Zones: the fills saved in the file are used. `snapshot.load_or_extract` first refills them with
`kicad-cli pcb drc --refill-zones --save-board` on a copy (the user's file is never written).
"""
import math
import os

from . import shapes, use

READER = "kicad-monkey"


def _layers(tokens, copper):
    """Layer tokens of a pad ("*.Cu", "F&B.Cu", "*.Mask", ...) as layer names."""
    out = []
    for t in tokens or ():
        t = str(t)
        if t == "*.Cu":
            out += list(copper)
        elif t == "F&B.Cu":
            out += [copper[0], copper[-1]]
        elif t.startswith("*."):
            out += ["F." + t[2:], "B." + t[2:]]
        elif t.startswith("F&B."):
            out += ["F." + t[4:], "B." + t[4:]]
        else:
            out.append(t)
    return list(dict.fromkeys(out))


def _xform(fp):
    a = math.radians(-(fp.at_angle or 0.0))
    c, s = math.cos(a), math.sin(a)
    return lambda x, y: (fp.at_x + x * c - y * s, fp.at_y + x * s + y * c)


def _enum(v):
    s = str(getattr(v, "value", v) or "")
    return s.split(".")[-1].lower()


def _raw_items(sx, head):
    return [it for it in (sx or [])[1:] if isinstance(it, list) and it and str(it[0]) == head]


def _graphic_chains(items, tf=None):
    """Point chains of graphic items (lines, arcs, rects, circles, polys) in board coordinates."""
    out = []
    t = tf or (lambda x, y: (x, y))
    for g in items:
        kind = type(g).__name__.lower()
        try:
            if "line" in kind:
                out.append([t(g.start_x, g.start_y), t(g.end_x, g.end_y)])
            elif "arc" in kind:
                out.append([t(x, y) for x, y in shapes._arc_pts((g.start_x, g.start_y), (g.mid_x, g.mid_y), (g.end_x, g.end_y))])
            elif "rect" in kind:
                x0, y0, x1, y1 = g.start_x, g.start_y, g.end_x, g.end_y
                out.append([t(x0, y0), t(x1, y0), t(x1, y1), t(x0, y1), t(x0, y0)])
            elif "circle" in kind:
                r = math.hypot(g.end_x - g.center_x, g.end_y - g.center_y)
                n = max(16, int(2 * math.pi / (2 * math.acos(max(-1, 1 - shapes.ARC_ERR / max(r, 1e-6))))))
                pts = [t(g.center_x + r * math.cos(2 * math.pi * k / n), g.center_y + r * math.sin(2 * math.pi * k / n)) for k in range(n)]
                out.append(pts + [pts[0]])
            elif "poly" in kind:
                pts = [t(x, y) for x, y in (g.points or [])]
                if pts:
                    out.append(pts + [pts[0]])
        except (AttributeError, TypeError):
            continue
    return out


def _filled(g):
    f = getattr(g, "fill", None)
    return str(getattr(f, "value", f) or "").lower() in ("yes", "solid", "true")


def refilled_copy(board_path, workdir):
    """A copy of the board (with its project files) whose zone fills KiCad has just recomputed:
    `kicad-cli pcb drc --refill-zones --save-board` on the copy. The user's file is never written.
    Returns (path of the copy, whether the refill ran)."""
    import shutil
    import subprocess
    os.makedirs(workdir, exist_ok=True)
    stem = os.path.splitext(os.path.basename(board_path))[0]
    for ext in (".kicad_pcb", ".kicad_pro", ".kicad_dru"):
        src = os.path.splitext(board_path)[0] + ext
        if os.path.exists(src):
            shutil.copy(src, os.path.join(workdir, stem + ext))
    cp = os.path.join(workdir, stem + ".kicad_pcb")
    try:
        r = subprocess.run(["kicad-cli", "pcb", "drc", "--refill-zones", "--save-board", "--format", "json",
                            "-o", os.path.join(workdir, "refill_drc.json"), cp], capture_output=True, text=True, timeout=1800)
        return cp, r.returncode in (0, 5)              # 5: DRC found violations; the board is saved either way
    except (OSError, subprocess.TimeoutExpired):
        return cp, False


def extract(board_path, refill_dir=None):
    """The snapshot of a board file, read with kicad-monkey (same fields as snapshot.extract). With
    `refill_dir`, the zones are refilled first by KiCad on a copy there."""
    lk = use()
    refilled = False
    if refill_dir:
        src_path = board_path
        board_path, refilled = refilled_copy(board_path, refill_dir)
    from kicad_monkey import KiCadPcb
    from kicad_monkey.kicad_pcb_layer_policy import PcbLayerFlashResolver, copper_layer_stack
    from shapely import geometry, ops
    from .. import padgeom
    from ..snapshot import file_sha256, parse_stackup, _arc_points, _dnp_field
    text = open(board_path, encoding="utf-8").read()
    pcb = KiCadPcb(board_path)
    copper = list(copper_layer_stack(pcb))
    flash = PcbLayerFlashResolver.from_board(pcb)
    snap = {"version": None, "board": os.path.abspath(src_path if refill_dir else board_path),
            "sha256": file_sha256(src_path if refill_dir else board_path),
            "kicad": f"file read by {READER} {lk['kicad_monkey']['version']} (KiCad file {pcb.version})", "reader": READER,
            "copper_layers": copper, "stackup": parse_stackup(text), "board_thickness": float(pcb.thickness or 0.0),
            "zones_refilled": refilled}
    from kicad_monkey.kicad_base import get_value
    smw = get_value(pcb.setup_sexp, "solder_mask_min_width") if pcb.setup_sexp else None
    snap["design"] = {"mask_expansion": float(pcb.pad_to_mask_clearance or 0.0), "mask_min_width": float(smw or 0.0)}
    # outline: Edge.Cuts of the board and of footprints
    chains = _graphic_chains(pcb.top_level_outline_items("Edge.Cuts"))
    for fp in pcb.footprints:
        tf = _xform(fp)
        chains += _graphic_chains([g for g in (fp.fp_lines + fp.fp_arcs + fp.fp_rects + fp.fp_circles + fp.fp_polys)
                                   if getattr(g, "layer", None) == "Edge.Cuts"], tf)
    snap["outline"] = shapes.outline(chains)
    # nets: KiCad 10 files carry names only; classes come from the project
    names = set()
    for coll in (pcb.segments, pcb.arcs, pcb.vias, pcb.zones):
        names |= {it.net.name for it in coll if getattr(it, "net", None) and it.net.name}
    for fp in pcb.footprints:
        names |= {p.net.name for p in fp.pads if getattr(p, "net", None) and p.net.name}
    for n in pcb.nets or []:
        if n.name:
            names.add(n.name)
    try:
        from ..project import Project
        pj = Project(board_path)
        cls = lambda n: pj.netclass(n)
    except Exception:
        cls = lambda n: "Default"
    snap["nets"] = {n: {"code": k + 1, "netclass": cls(n)} for k, n in enumerate(sorted(names))}
    tracks = []
    for t in pcb.segments:
        tracks.append({"x1": t.start_x, "y1": t.start_y, "x2": t.end_x, "y2": t.end_y, "w": t.width, "layer": t.layer, "net": t.net.name or ""})
    for a in pcb.arcs:
        pts = _arc_points((a.start_x, a.start_y), (a.mid_x, a.mid_y), (a.end_x, a.end_y))
        for (x1, y1), (x2, y2) in zip(pts[:-1], pts[1:]):
            tracks.append({"x1": x1, "y1": y1, "x2": x2, "y2": y2, "w": a.width, "layer": a.layer, "net": a.net.name or "", "arc": True})
    snap["tracks"] = tracks
    vias = []
    vtype = {None: 3, "through": 3, "blind": 2, "buried": 2, "micro": 1}    # pcbnew's VIATYPE values
    for v in pcb.vias:
        fl = [l for l in flash.via_flash_layers(v) if l in copper]
        span = list(v.layers or [copper[0], copper[-1]])
        ten = getattr(v, "tenting", None)
        # back-drilling (KiCad 10: a secondary or tertiary drill with the two layers it spans): from the
        # outer layer it starts at, up to the other one
        backdrill = None
        for attr in ("backdrill", "tertiary_drill"):
            bd = getattr(v, attr, None)
            ly = getattr(bd, "layers", None) if bd is not None else None
            if ly is None:
                continue
            pair = [getattr(ly, "start", None), getattr(ly, "end", None)]
            backdrill = backdrill or {"mode": None, "bottom_to": None, "top_to": None}
            if copper[-1] in pair:
                backdrill["bottom_to"] = next((l for l in pair if l != copper[-1]), None)
            elif copper[0] in pair:
                backdrill["top_to"] = next((l for l in pair if l != copper[0]), None)
        vias.append({"x": v.at_x, "y": v.at_y, "d": v.size, "drill": v.drill, "net": v.net.name or "", "span": span[:2],
                     "flashed": fl, "type": vtype.get(str(v.via_type).lower() if v.via_type else None, 3), "backdrill": backdrill,
                     "tented": [getattr(ten, "front", None), getattr(ten, "back", None)] if ten is not None else [None, None],
                     "mask_exp": None})
    snap["vias"] = vias
    zones = []
    for z in pcb.zones:
        if z.keepout is not None:
            outline = [[[(x, y) for x, y in poly.points]] for poly in (z.polygons or [])]
            zones.append({"rule_area": True, "name": z.name, "layers": [l for l in (z.layers or []) if l in copper], "outline": outline})
            continue
        for fpoly in z.filled_polygons or []:
            if fpoly.layer not in copper or len(fpoly.points) < 3:
                continue
            zones.append({"net": z.net.name if z.net else "", "layer": fpoly.layer, "priority": int(z.priority or 0), "name": z.name,
                          "teardrop": "teardrop" in str(getattr(z, "_raw_sexp", "")),
                          "polys": [[[(round(x, 5), round(y, 5)) for x, y in fpoly.points]]]})
    snap["zones"] = zones
    fps, pads, mask_shapes = [], [], []
    board_mask = float(pcb.pad_to_mask_clearance or 0.0)
    board_paste = float(getattr(pcb, "pad_to_paste_clearance", 0.0) or 0.0)
    board_ratio = float(getattr(pcb, "pad_to_paste_clearance_ratio", 0.0) or 0.0)
    for fp in pcb.footprints:
        tf = _xform(fp)
        ref = fp.get_property_value("Reference") or ""
        fields = {p.name: p.value for p in fp.properties if p.name not in ("Reference", "Value", "Footprint") and p.value}
        attrs = set(fp.attr or [])
        side = "B" if fp.layer == "B.Cu" else "F"
        crt = []
        for lay in ("F.CrtYd", "B.CrtYd"):
            ch = _graphic_chains([g for g in (fp.fp_lines + fp.fp_arcs + fp.fp_rects + fp.fp_circles + fp.fp_polys)
                                  if getattr(g, "layer", None) == lay], tf)
            crt += shapes.outline(ch)
        fab = _graphic_chains([g for g in (fp.fp_lines + fp.fp_arcs + fp.fp_rects + fp.fp_circles + fp.fp_polys)
                               if getattr(g, "layer", None) in ("F.Fab", "B.Fab")], tf)
        body = None
        pts = [p for c in fab for p in c]
        if len(pts) >= 3:
            hull = geometry.MultiPoint(pts).convex_hull
            if hull.geom_type == "Polygon":
                body = [[round(x, 5), round(y, 5)] for x, y in list(hull.exterior.coords)[:-1]]
        for g in fp.fp_lines + fp.fp_arcs + fp.fp_rects + fp.fp_circles + fp.fp_polys:      # graphics on the mask layers
            if getattr(g, "layer", None) in ("F.Mask", "B.Mask"):
                w = float(getattr(getattr(g, "stroke", None), "width", 0.0) or 0.0)
                for ch in _graphic_chains([g], tf):
                    gg = geometry.Polygon(ch).buffer(w / 2) if (_filled(g) and len(ch) >= 4) else geometry.LineString(ch).buffer(max(w, 1e-6) / 2)
                    mask_shapes.append({"layer": g.layer, "polys": shapes.rings(gg)})
        fpads = []
        for p in fp.pads:
            kind = {"thru_hole": "pth", "smd": "smd", "connect": "conn", "np_thru_hole": "npth"}.get(_enum(p.pad_type), "other")
            lays = _layers(p.layers, copper)
            cu_on = [l for l in copper if l in lays]
            if kind == "pth":
                cu_on = [l for l in flash.pad_item_layers(p, fp) if l in copper] or cu_on
            elif kind == "npth" and max(p.size_x, p.size_y) <= float(p.drill or 0.0) + 1e-6:
                cu_on = []                                 # a bare hole; a larger pad keeps its copper ring
            px, py = tf(p.at_x, p.at_y)
            rot = float(p.at_angle or 0.0)
            shp = _enum(p.shape)
            prims = [prim._raw_sexp for prim in (p.custom_primitives or []) if getattr(prim, "_raw_sexp", None)]
            anchor = _enum(getattr(getattr(p, "custom_options", None), "anchor", None) or "rect") or "rect"
            corners = [str(c).lower() for c in (p.chamfer_corners or [])]
            loc = shapes.pad_local(shp, p.size_x, p.size_y, p.roundrect_rratio or 0.0, p.chamfer_ratio or 0.0, corners,
                                   (p.rect_delta_x or 0.0, p.rect_delta_y or 0.0), anchor, prims)
            g = shapes.to_board(loc, px, py, rot)
            polys = shapes.rings(g) if cu_on else []
            mask, paste, paste_poly = {}, {}, {}
            # KiCad applies mask and paste margins only to pads on a copper layer: an aperture pad (paste-
            # or mask-only) is its drawn shape
            on_cu = any(l.endswith(".Cu") for l in lays)
            for mname, pname in (("F.Mask", "F.Paste"), ("B.Mask", "B.Paste")):
                if mname in lays:
                    exp = next((v for v in (p.solder_mask_margin, fp.solder_mask_margin) if v is not None), board_mask) if on_cu else 0.0
                    mask[mname] = {"exp": float(exp), "polys": shapes.rings(shapes.grow(g, float(exp)))}
                if pname in lays:
                    m = next((v for v in (p.solder_paste_margin, fp.solder_paste_margin) if v is not None), board_paste)
                    r = next((v for v in (p.solder_paste_margin_ratio, getattr(fp, "solder_paste_margin_ratio", None)) if v is not None), board_ratio)
                    mx, my = (float(m) + float(r) * p.size_x, float(m) + float(r) * p.size_y) if on_cu else (0.0, 0.0)
                    paste[pname] = [mx, my]
                    if shp == "custom":
                        pg = shapes.grow(g, mx)
                    elif abs(mx) > 1e-12 or abs(my) > 1e-12:
                        pg = shapes.to_board(shapes.pad_local(shp, max(p.size_x + 2 * mx, 1e-6), max(p.size_y + 2 * my, 1e-6),
                                                              p.roundrect_rratio or 0.0, p.chamfer_ratio or 0.0, corners,
                                                              (p.rect_delta_x or 0.0, p.rect_delta_y or 0.0), anchor, prims), px, py, rot)
                    else:
                        pg = g
                    paste_poly[pname] = shapes.rings(pg)
            own = polys or next((v for v in paste_poly.values() if v), None) or next((v["polys"] for v in mask.values() if v["polys"]), None) or []
            dims = padgeom.measure(own, px, py, rot)
            custom = shp == "custom"
            dx_ = float(p.drill or 0.0) if kind in ("pth", "npth") else 0.0
            dy_ = dx_
            if kind in ("pth", "npth") and getattr(p, "drill_oval", False):
                dx_, dy_ = float(p.drill_width or dx_), float(p.drill_height or dx_)
            rec = {"ref": ref, "num": p.number or "", "net": (p.net.name if getattr(p, "net", None) else "") or "", "kind": kind,
                   "x": px, "y": py, "rot": rot, "shape": shp,
                   "sx": dims["w"] if (dims and custom) else p.size_x, "sy": dims["h"] if (dims and custom) else p.size_y,
                   "anchor": [p.size_x, p.size_y], "dims": dims, "drill": dx_, "drill_y": dy_,
                   "layers": cu_on, "poly": polys, "pinfunc": p.pinfunction or "", "pintype": p.pintype or "",
                   "mask": mask, "paste": paste, "paste_poly": paste_poly}
            pads.append(rec)
            fpads.append(rec)
        allpts = [q for r_ in fpads for poly in r_["poly"] for q in poly[0]] + [q for c in fab for q in c]
        bb = ([min(q[0] for q in allpts), min(q[1] for q in allpts), max(q[0] for q in allpts), max(q[1] for q in allpts)]
              if allpts else [fp.at_x, fp.at_y, fp.at_x, fp.at_y])
        fps.append({"ref": ref, "value": fp.get_property_value("Value") or "", "fpid": fp.library_link or "",
                    "x": fp.at_x, "y": fp.at_y, "rot": float(fp.at_angle or 0.0), "side": side,
                    "dnp": "dnp" in attrs or _dnp_field(fields), "fields": fields, "courtyard": crt, "smd": "smd" in attrs,
                    "placed": not ({"exclude_from_pos_files", "board_only"} & attrs), "in_pos": "exclude_from_pos_files" not in attrs,
                    "in_bom": "exclude_from_bom" not in attrs, "board_only": "board_only" in attrs, "bbox": bb, "body": body})
    for g in pcb.top_level_outline_items("F.Mask") + pcb.top_level_outline_items("B.Mask") if hasattr(pcb, "top_level_outline_items") else []:
        w = float(getattr(getattr(g, "stroke", None), "width", 0.0) or 0.0)
        for ch in _graphic_chains([g]):
            gg = geometry.Polygon(ch).buffer(w / 2) if (_filled(g) and len(ch) >= 4) else geometry.LineString(ch).buffer(max(w, 1e-6) / 2)
            mask_shapes.append({"layer": g.layer, "polys": shapes.rings(gg)})
    snap["footprints"] = fps
    snap["pads"] = pads
    snap["mask_shapes"] = mask_shapes
    return snap
