"""Board snapshot: everything openpcbsi needs from a KiCad 10 board, read once.

Two readers give the same fields: the file reader (openpcbsi.parse, on kicad-monkey; the default) and
KiCad's SWIG `pcbnew` module (OPENPCBSI_READER=pcbnew; KiCad has deprecated it).

The snapshot is plain JSON (gzip), cached by the SHA-256 of the board file, so that every later
step (and the openEMS venv, which has no pcbnew) works from the same frozen data. Units: mm,
KiCad board coordinates (y grows downwards).

Only this module imports pcbnew, and only for the pcbnew reader.
"""
import gzip
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys

from . import padgeom

SNAPSHOT_VERSION = 12  # 12: copper graphics and text (copper_gfx); no-net copper is copper
                       # 11: coordinates kept to 1 nm (KiCad's unit), each standard pad's exact core and radius
                       #     (padgeom.core), each arc chord with its arc's centre and radius, padstacks per layer
                       # 10: footprint body outline (convex hull of its fabrication-layer shapes)
                       # 9: pad shape, anchor size kept apart, sizes measured from the pad's polygon, paste
                       #    apertures as polygons (margins applied), polygons of pads with no copper
                       # 4: zones refilled in memory, every footprint field, pin types, file version; 5: back-drill;
                       # 6: DNP also from a "DNP" field (KiCad 6 designs mark it that way)
                       # 7: solder-mask openings and paste margins of pads, mask shapes, via tenting (MF checks)
                       # 8: slot length of pad holes, via mask opening, footprint attributes, board mask settings
NM = 1e-6  # nm -> mm


def file_sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# ----------------------------------------------------------------------------------------------
# s-expression helpers (the stackup is not exposed through the SWIG API of KiCad 10)

def _sexpr_block(text, head):
    i = text.find("(" + head)
    if i < 0:
        return None
    depth = 0
    for j in range(i, len(text)):
        c = text[j]
        if c == "(":
            depth += 1
        elif c == ")":
            depth -= 1
            if depth == 0:
                return text[i:j + 1]
    return None


def parse_stackup(board_text):
    """Return the ordered stackup list [{name, type, thickness, epsilon_r, loss_tangent, material}]."""
    blk = _sexpr_block(board_text, "stackup")
    out = []
    if blk is None:
        return out
    for m in re.finditer(r'\(layer "([^"]+)"(.*?)\)\s*(?=\(layer|\(copper_finish|\(dielectric_constraints|\(edge|\)\s*$)',
                         blk, re.S):
        name, body = m.group(1), m.group(2)

        def val(key, cast=float, text=body):
            mm = re.search(r'\(%s "?([^")]+)"?\)' % key, text)
            return cast(mm.group(1)) if mm else None

        # a dielectric of several plies ("addsublayer"): the thicknesses add, permittivity and loss tangent
        # are weighted by thickness
        plies = re.split(r"\baddsublayer\b", body)
        ts = [val("thickness", text=p) for p in plies]
        er = [val("epsilon_r", text=p) for p in plies]
        tl = [val("loss_tangent", text=p) for p in plies]
        thick = sum(t for t in ts if t) if any(ts) else None

        def weighted(vals):
            pairs = [(t, v) for t, v in zip(ts, vals) if t and v is not None]
            if not pairs:
                return vals[0] if vals else None
            return sum(t * v for t, v in pairs) / sum(t for t, _ in pairs)

        out.append({
            "name": name,
            "type": val("type", str),
            "thickness": thick,
            "epsilon_r": weighted(er),
            "loss_tangent": weighted(tl),
            "material": val("material", str),
            "plies": len(plies),
        })
    fin = re.search(r'\(copper_finish "([^"]+)"\)', blk)
    if fin:
        out.append({"name": "_finish", "type": "finish", "material": fin.group(1)})
    return out


# ----------------------------------------------------------------------------------------------

def _chain_pts(chain):
    n = chain.PointCount()
    pts = []
    for k in range(n):
        p = chain.CPoint(k)
        pts.append((round(p.x * NM, 6), round(p.y * NM, 6)))
    return pts


def _polyset(ps):
    """SHAPE_POLY_SET -> list of polygons, each [outline, hole, hole...] as point lists."""
    polys = []
    for i in range(ps.OutlineCount()):
        rings = [_chain_pts(ps.Outline(i))]
        for h in range(ps.HoleCount(i)):
            rings.append(_chain_pts(ps.Hole(i, h)))
        polys.append(rings)
    return polys


def extract(board_path):
    import pcbnew  # noqa: only here

    with open(board_path, encoding="utf-8") as f:
        text = f.read()
    b = pcbnew.LoadBoard(board_path)
    # #3: the fills saved in the file may be stale (or missing) after an edit; refill in memory, as
    # KiCad's DRC does with --refill-zones. The file itself is never written.
    refilled = False
    if os.environ.get("OPENPCBSI_NO_REFILL") != "1":
        zones = [z for z in b.Zones() if not z.GetIsRuleArea()]
        if zones:
            refilled = bool(pcbnew.ZONE_FILLER(b).Fill(b.Zones())) or True
    cu_ids = list(b.GetEnabledLayers().CuStack())
    cu_names = [pcbnew.LayerName(l) for l in cu_ids]
    lid2name = dict(zip(cu_ids, cu_names))

    snap = {
        "version": SNAPSHOT_VERSION,
        "board": os.path.abspath(board_path),
        "sha256": file_sha256(board_path),
        "kicad": pcbnew.Version(),
        "copper_layers": cu_names,
        "stackup": parse_stackup(text),
        "board_thickness": b.GetDesignSettings().GetBoardThickness() * NM,
        "zones_refilled": refilled,
    }
    ds = b.GetDesignSettings()
    snap["design"] = {"mask_expansion": ds.m_SolderMaskExpansion * NM, "mask_min_width": ds.m_SolderMaskMinWidth * NM}

    ps = pcbnew.SHAPE_POLY_SET()
    b.GetBoardPolygonOutlines(ps, True)
    snap["outline"] = _polyset(ps)
    # the circles of the outline's arcs: a distance to an outline chord on one of them carries the
    # chord's sagitta as its error bound (exact.outline_index); straight edges are exact
    arcs = []
    for it in list(b.GetDrawings()) + [g for f in b.GetFootprints() for g in f.GraphicalItems()]:
        try:
            if it.GetLayer() == pcbnew.Edge_Cuts and it.GetShape() in (pcbnew.SHAPE_T_ARC, pcbnew.SHAPE_T_CIRCLE):
                c = it.GetCenter()
                arcs.append([c.x * NM, c.y * NM, it.GetRadius() * NM])
        except AttributeError:
            continue
    snap["outline_arcs"] = arcs

    # nets and classes
    nets = {}
    for name, n in b.GetNetsByName().items():
        name = str(name)
        if not name:
            continue
        nets[name] = {"code": n.GetNetCode(), "netclass": str(n.GetNetClassName())}
    snap["nets"] = nets

    tracks, vias = [], []
    for t in b.GetTracks():
        ty = t.Type()
        net = t.GetNetname()
        if ty == pcbnew.PCB_VIA_T:
            flashed = [lid2name[l] for l in cu_ids if t.FlashLayer(l)]
            span = [lid2name.get(t.TopLayer()), lid2name.get(t.BottomLayer())]
            vias.append({
                "x": t.GetPosition().x * NM, "y": t.GetPosition().y * NM,
                "d": t.GetWidth(cu_ids[0]) * NM, "drill": t.GetDrillValue() * NM,
                "net": net, "span": span, "flashed": flashed, "type": int(t.GetViaType()),
                "backdrill": _backdrill(pcbnew, t, lid2name), "tented": _tented(pcbnew, t),
                "mask_exp": _via_mask_exp(pcbnew, t), "protect": _protect(pcbnew, t, b.GetDesignSettings()),
            })
        elif ty == pcbnew.PCB_ARC_T:
            # approximate arcs by chords of at most 10 degrees
            s, m, e = t.GetStart(), t.GetMid(), t.GetEnd()
            pts = _arc_points((s.x * NM, s.y * NM), (m.x * NM, m.y * NM), (e.x * NM, e.y * NM))
            arc = _arc_circle((s.x * NM, s.y * NM), (m.x * NM, m.y * NM), (e.x * NM, e.y * NM)) or True
            for (x1, y1), (x2, y2) in zip(pts[:-1], pts[1:]):
                tracks.append({"x1": x1, "y1": y1, "x2": x2, "y2": y2, "w": t.GetWidth() * NM,
                               "layer": lid2name.get(t.GetLayer()), "net": net, "arc": arc})
        else:
            s, e = t.GetStart(), t.GetEnd()
            tracks.append({"x1": s.x * NM, "y1": s.y * NM, "x2": e.x * NM, "y2": e.y * NM,
                           "w": t.GetWidth() * NM, "layer": lid2name.get(t.GetLayer()), "net": net})
    snap["tracks"] = tracks
    snap["vias"] = vias

    zones = []
    for z in b.Zones():
        if z.GetIsRuleArea():
            zones.append({"rule_area": True, "name": z.GetZoneName(),
                          "layers": [lid2name[l] for l in z.GetLayerSet().Seq() if l in lid2name],
                          "outline": _polyset(z.Outline())})
            continue
        for l in z.GetLayerSet().Seq():
            if l not in lid2name:
                continue
            if not z.IsFilled():
                continue
            fp = z.GetFilledPolysList(l)
            zones.append({"net": z.GetNetname(), "layer": lid2name[l], "priority": z.GetAssignedPriority(),
                          "name": z.GetZoneName(), "teardrop": bool(getattr(z, "IsTeardropArea", lambda: False)()), "polys": _polyset(fp)})
    snap["zones"] = zones

    fps, pads = [], []
    for f in b.GetFootprints():
        ref = f.GetReference()
        fields = {}
        for fld in f.GetFields():           # every user field: design data lives here (#4)
            k, v = fld.GetName(), fld.GetText()
            if k not in ("Reference", "Value", "Footprint") and v:
                fields[k] = v
        crt = []
        for side in (pcbnew.F_CrtYd, pcbnew.B_CrtYd):
            try:
                c = f.GetCourtyard(side)
                crt += _polyset(c)
            except Exception:
                pass
        bb = f.GetBoundingBox(False)
        body = _fab_body(pcbnew, f)
        fps.append({
            "ref": ref, "value": f.GetValue(), "fpid": f.GetFPIDAsString(),
            "x": f.GetPosition().x * NM, "y": f.GetPosition().y * NM,
            "rot": f.GetOrientationDegrees(), "side": "B" if f.IsFlipped() else "F",
            "dnp": bool(f.IsDNP()) or _dnp_field(fields), "fields": fields, "courtyard": crt,
            "smd": bool(f.GetAttributes() & pcbnew.FP_SMD),
            "placed": not (f.GetAttributes() & (pcbnew.FP_EXCLUDE_FROM_POS_FILES | pcbnew.FP_BOARD_ONLY)),
            "in_pos": not (f.GetAttributes() & pcbnew.FP_EXCLUDE_FROM_POS_FILES),
            "in_bom": not (f.GetAttributes() & pcbnew.FP_EXCLUDE_FROM_BOM),
            "board_only": bool(f.GetAttributes() & pcbnew.FP_BOARD_ONLY),
            "bbox": [bb.GetX() * NM, bb.GetY() * NM, bb.GetRight() * NM, bb.GetBottom() * NM], "body": body,
        })
        for p in f.Pads():
            attr = p.GetAttribute()
            kind = {pcbnew.PAD_ATTRIB_PTH: "pth", pcbnew.PAD_ATTRIB_SMD: "smd",
                    pcbnew.PAD_ATTRIB_CONN: "conn", pcbnew.PAD_ATTRIB_NPTH: "npth"}.get(attr, "other")
            layers = [lid2name[l] for l in cu_ids if p.IsOnLayer(l) and (kind != "pth" or p.FlashLayer(l))]
            # one polygon per pad (from its first copper layer; pad stacks are not distinguished)
            polys = []
            lid0 = b.GetLayerID(layers[0]) if layers else cu_ids[0]
            if layers:
                s = pcbnew.SHAPE_POLY_SET()
                p.TransformShapeToPolygon(s, lid0, 0, 2000, pcbnew.ERROR_INSIDE)
                polys = _polyset(s)
            if kind == "npth" and layers and max(p.GetSize(lid0).x, p.GetSize(lid0).y) <= p.GetDrillSize().x + 1000:
                layers, polys = [], []                     # a bare hole: no copper ring (KiCad plots none)
            mp = _pad_mask_paste(pcbnew, p)
            # the layer that defines the pad's shape: its copper, else its paste or mask opening (pads
            # that exist only on a technical layer: paste-only apertures, mask-only openings)
            own = polys or next((v for v in mp["paste_poly"].values() if v), None) or \
                next((v["polys"] for v in mp["mask"].values() if v.get("polys")), None) or []
            slay = lid0 if layers else next((l for l in (pcbnew.F_Paste, pcbnew.B_Paste, pcbnew.F_Mask, pcbnew.B_Mask)
                                             if p.IsOnLayer(l)), lid0)
            sz = p.GetSize(slay)
            px, py, rot = p.GetPosition().x * NM, p.GetPosition().y * NM, p.GetOrientationDegrees()
            dims = padgeom.measure(own, px, py, rot)
            def core_of(lid):
                shp_ = _shape_name(pcbnew, p.GetShape(lid))
                sz_ = p.GetSize(lid)
                dl = p.GetDelta(lid) if shp_ == "trapezoid" else None
                return padgeom.core(shp_, sz_.x * NM, sz_.y * NM, px, py, rot, p.GetRoundRectCornerRadius(lid) * NM if shp_ in
                                    ("roundrect", "chamfered_rect") else 0.0, (dl.x * NM, dl.y * NM) if dl is not None else (0.0, 0.0),
                                    shp_ == "chamfered_rect" and p.GetChamferRectRatio(lid) > 0 and p.GetChamferPositions(lid) != 0)
            core = core_of(slay)
            # padstacks (KiCad 9+): a pad whose shape differs per copper layer gets every layer's polygon
            # and core; `poly` stays the shape of its first copper layer
            stack = {}
            try:
                if p.Padstack().Mode() != pcbnew.PADSTACK.MODE_NORMAL and len(layers) > 1:
                    for ln in layers:
                        lid = b.GetLayerID(ln)
                        s_ = pcbnew.SHAPE_POLY_SET()
                        p.TransformShapeToPolygon(s_, lid, 0, 2000, pcbnew.ERROR_INSIDE)
                        sl = p.GetSize(lid)
                        stack.setdefault("layer_poly", {})[ln] = _polyset(s_)
                        stack.setdefault("layer_core", {})[ln] = core_of(lid)
                        stack.setdefault("layer_size", {})[ln] = [sl.x * NM, sl.y * NM]
                    stack["padstack"] = {pcbnew.PADSTACK.MODE_FRONT_INNER_BACK: "front_inner_back",
                                         pcbnew.PADSTACK.MODE_CUSTOM: "custom"}.get(p.Padstack().Mode(), "other")
            except AttributeError:
                pass
            pads.append({
                "ref": ref, "num": p.GetNumber(), "net": p.GetNetname(), "kind": kind, "core": core,
                "x": px, "y": py, "rot": rot, "shape": _shape_name(pcbnew, p.GetShape(slay)),
                # sx, sy: the pad's extent in its own frame. A standard shape's size is exact; a custom
                # pad's is measured from its polygon (its anchor, often 0.2 mm, is only `anchor` and no
                # check reads it)
                "sx": dims["w"] if (dims and p.GetShape(slay) == pcbnew.PAD_SHAPE_CUSTOM) else sz.x * NM,
                "sy": dims["h"] if (dims and p.GetShape(slay) == pcbnew.PAD_SHAPE_CUSTOM) else sz.y * NM,
                "anchor": [sz.x * NM, sz.y * NM], "dims": dims,
                "drill": p.GetDrillSize().x * NM if kind in ("pth", "npth") else 0.0,
                "drill_y": p.GetDrillSize().y * NM if kind in ("pth", "npth") else 0.0,
                "layers": layers, "poly": polys, "pinfunc": p.GetPinFunction(), "pintype": p.GetPinType(),
                **mp, **stack,
            })
    snap["footprints"] = fps
    snap["pads"] = pads
    snap["mask_shapes"] = _mask_shapes(pcbnew, b)
    snap["copper_gfx"] = _copper_gfx(pcbnew, b, lid2name)
    return snap


MASK_PASTE = (("F.Mask", "F_Mask", "F.Paste", "F_Paste"), ("B.Mask", "B_Mask", "B.Paste", "B_Paste"))


def _shape_name(pcbnew, shape):
    for k in dir(pcbnew):
        if k.startswith("PAD_SHAPE_") and getattr(pcbnew, k) == shape:
            return k[len("PAD_SHAPE_"):].lower()
    return str(shape)


def _paste_polys(pcbnew, p, pl, mx, my):
    """The paste aperture of a pad on one paste layer, as KiCad plots it: the pad's shape on that
    layer with the paste margin applied (a standard shape is resized by the margin on each axis; a
    custom shape is grown or shrunk by it)."""
    s = pcbnew.SHAPE_POLY_SET()
    if abs(mx) < 1e-9 and abs(my) < 1e-9:
        p.TransformShapeToPolygon(s, pl, 0, 2000, pcbnew.ERROR_INSIDE)
        return _polyset(s)
    if p.GetShape(pl) == pcbnew.PAD_SHAPE_CUSTOM:
        p.TransformShapeToPolygon(s, pl, 0, 2000, pcbnew.ERROR_INSIDE)
        s.Inflate(int(round(mx / NM)), pcbnew.CORNER_STRATEGY_ROUND_ALL_CORNERS, 2000)
        return _polyset(s)
    q = pcbnew.PAD(p)
    sz = p.GetSize(pl)
    q.SetSize(pl, pcbnew.VECTOR2I(max(int(sz.x + 2 * mx / NM), 1), max(int(sz.y + 2 * my / NM), 1)))
    q.TransformShapeToPolygon(s, pl, 0, 2000, pcbnew.ERROR_INSIDE)
    return _polyset(s)


def _pad_mask_paste(pcbnew, p):
    """Solder-mask opening (the pad shape grown by its mask expansion) and paste margin and aperture
    of a pad, per side. `mask` is {mask layer: polygons}; a pad that is not on a mask layer is covered
    by mask. `paste_poly` is {paste layer: polygons}, the aperture as KiCad plots it."""
    mask, paste, paste_poly = {}, {}, {}
    for mname, mattr, pname, pattr in MASK_PASTE:
        ml, pl = getattr(pcbnew, mattr), getattr(pcbnew, pattr)
        if p.IsOnLayer(ml):
            try:
                exp = int(p.GetSolderMaskExpansion(ml))
                s = pcbnew.SHAPE_POLY_SET()
                p.TransformShapeToPolygon(s, ml, exp, 2000, pcbnew.ERROR_OUTSIDE)
                mask[mname] = {"exp": exp * NM, "polys": _polyset(s)}
            except Exception:
                pass
        if p.IsOnLayer(pl):
            try:
                m = p.GetSolderPasteMargin(pl)
                paste[pname] = [m.x * NM, m.y * NM]
            except Exception:
                paste[pname] = [0.0, 0.0]
            try:
                paste_poly[pname] = _paste_polys(pcbnew, p, pl, *paste[pname])
            except Exception:
                paste_poly[pname] = []
    return {"mask": mask, "paste": paste, "paste_poly": paste_poly}


def _via_mask_exp(pcbnew, v):
    """Mask expansion of a via's opening (mm), where it is not tented; None if unknown."""
    try:
        return v.GetSolderMaskExpansion() * NM
    except Exception:
        return None


def _tented(pcbnew, v):
    """[front, back]: whether the via is covered by solder mask on each side (None if unknown)."""
    out = []
    for attr in ("F_Mask", "B_Mask"):
        try:
            out.append(bool(v.IsTented(getattr(pcbnew, attr))))
        except Exception:
            out.append(None)
    return out


def _protect(pcbnew, v, ds):
    """Via protection (KiCad 10, IPC-4761): {filled, capped, plugged: [front, back]}, each True / False
    as the via or the board setup states it, None where nothing is stated (the setup's "no" is KiCad's
    default, not a statement)."""
    def one(mode, yes, no, board):
        if mode == yes:
            return True
        if mode == no:
            return False
        return True if board else None
    try:
        return {"filled": one(v.GetFillingMode(), pcbnew.FILLING_MODE_FILLED, pcbnew.FILLING_MODE_NOT_FILLED, ds.m_FillVias),
                "capped": one(v.GetCappingMode(), pcbnew.CAPPING_MODE_CAPPED, pcbnew.CAPPING_MODE_NOT_CAPPED, ds.m_CapVias),
                "plugged": [one(v.GetFrontPluggingMode(), pcbnew.PLUGGING_MODE_PLUGGED, pcbnew.PLUGGING_MODE_NOT_PLUGGED, ds.m_PlugViasFront),
                            one(v.GetBackPluggingMode(), pcbnew.PLUGGING_MODE_PLUGGED, pcbnew.PLUGGING_MODE_NOT_PLUGGED, ds.m_PlugViasBack)]}
    except AttributeError:                    # KiCad before 9
        return {"filled": None, "capped": None, "plugged": [None, None]}


def _copper_gfx(pcbnew, b, lid2name):
    """Copper that is not a track, pad, via or zone: graphic shapes and text on copper layers, board
    and footprint level. Each item: {layer, net, src, segs: [[x1, y1, x2, y2, width]], arcs: [[x1, y1,
    x2, y2, width, cx, cy, r]], polys, exact}. Lines, arcs and circles are exact strokes (a filled
    circle is a zero-length stroke as wide as the disc); filled rectangles and polygons are their exact
    outline plus the stroke; text is KiCad's own polygons (2 um inside the true strokes), not exact."""
    out = []
    items = [(d, "") for d in b.GetDrawings()]
    for f in b.GetFootprints():
        ref = f.GetReference()
        items += [(d, ref) for d in f.GraphicalItems()]
        items += [(fl, ref) for fl in f.GetFields() if fl.IsVisible()]
    for d, ref in items:
        lid = d.GetLayer()
        if lid not in lid2name:
            continue
        cls = d.GetClass()
        net = d.GetNetname() if hasattr(d, "GetNetname") else ""
        it = {"layer": lid2name[lid], "net": net or "", "src": (f"{ref} " if ref else "") + cls, "segs": [], "arcs": [], "polys": [],
              "exact": True}
        if cls in ("PCB_SHAPE", "FP_SHAPE"):
            w = d.GetWidth() * NM
            st = d.GetShape()
            P = lambda v: (v.x * NM, v.y * NM)
            if st == pcbnew.SHAPE_T_SEGMENT:
                (x1, y1), (x2, y2) = P(d.GetStart()), P(d.GetEnd())
                it["segs"].append([x1, y1, x2, y2, w])
            elif st == pcbnew.SHAPE_T_ARC:
                a, m, e = P(d.GetStart()), P(d.GetArcMid()), P(d.GetEnd())
                circ = _arc_circle(a, m, e)
                pts = _arc_points(a, m, e)
                for (x1, y1), (x2, y2) in zip(pts[:-1], pts[1:]):
                    it["arcs" if circ else "segs"].append([x1, y1, x2, y2, w] + (list(circ) if circ else []))
            elif st == pcbnew.SHAPE_T_CIRCLE:
                cx, cy = P(d.GetCenter())
                r = d.GetRadius() * NM
                if d.IsAnyFill():
                    it["segs"].append([cx, cy, cx, cy, 2 * r + w])
                else:
                    pts = [(cx + r * math.cos(2 * math.pi * k / 72), cy + r * math.sin(2 * math.pi * k / 72)) for k in range(73)]
                    for (x1, y1), (x2, y2) in zip(pts[:-1], pts[1:]):
                        it["arcs"].append([x1, y1, x2, y2, w, cx, cy, r])
            else:
                # rectangle, polygon, curve: the outline as strokes, and the filled area as a polygon
                if st == pcbnew.SHAPE_T_RECTANGLE:
                    ring = [P(v) for v in d.GetRectCorners()]
                elif st == pcbnew.SHAPE_T_POLY:
                    ring = _chain_pts(d.GetPolyShape().Outline(0))
                else:
                    ring = [P(v) for v in d.GetBezierPoints()]
                    it["exact"] = False
                closed = st != pcbnew.SHAPE_T_BEZIER
                edges = list(zip(ring, ring[1:] + (ring[:1] if closed else [])))
                if w > 0:
                    it["segs"] += [[a[0], a[1], c[0], c[1], w] for a, c in edges]
                if d.IsAnyFill() and closed:
                    it["polys"].append([[(round(x, 6), round(y, 6)) for x, y in ring]])
        elif cls in ("PCB_TEXT", "PCB_FIELD", "FP_TEXT", "PCB_TEXTBOX"):
            if hasattr(d, "IsVisible") and not d.IsVisible():
                continue
            s = pcbnew.SHAPE_POLY_SET()
            if getattr(d, "IsKnockout", lambda: False)():
                d.TransformShapeToPolygon(s, lid, 0, 1000, pcbnew.ERROR_INSIDE)
            else:
                d.TransformTextToPolySet(s, 0, 1000, pcbnew.ERROR_INSIDE)
            it["polys"] = _polyset(s)
            it["exact"] = False
        else:
            continue
        if it["segs"] or it["arcs"] or it["polys"]:
            out.append(it)
    return out


def _mask_shapes(pcbnew, b):
    """Graphic shapes on the mask layers (board and footprint level): further openings in the mask."""
    out = []
    ids = {pcbnew.F_Mask: "F.Mask", pcbnew.B_Mask: "B.Mask"}
    items = list(b.GetDrawings())
    for f in b.GetFootprints():
        items += list(f.GraphicalItems())
    for d in items:
        l = d.GetLayer()
        if l not in ids or not hasattr(d, "TransformShapeToPolygon"):
            continue
        try:
            s = pcbnew.SHAPE_POLY_SET()
            d.TransformShapeToPolygon(s, l, 0, 2000, pcbnew.ERROR_OUTSIDE)
            polys = _polyset(s)
        except Exception:
            continue
        if polys:
            out.append({"layer": ids[l], "polys": polys})
    return out


def _fab_body(pcbnew, f):
    """The part's body as drawn on its fabrication layer: the convex hull of the footprint's F.Fab / B.Fab
    shapes (text excluded), in board coordinates; None when the footprint draws none."""
    import numpy as np
    pts = []
    for d in f.GraphicalItems():
        if d.Type() != pcbnew.PCB_SHAPE_T or d.GetLayer() not in (pcbnew.F_Fab, pcbnew.B_Fab):
            continue
        s = pcbnew.SHAPE_POLY_SET()
        try:
            d.TransformShapeToPolygon(s, d.GetLayer(), 0, 2000, pcbnew.ERROR_INSIDE)
        except Exception:
            continue
        for poly in _polyset(s):
            for ring in poly:
                pts += ring
    if len(pts) < 3:
        return None
    from scipy.spatial import ConvexHull
    P = np.array(pts)
    try:
        h = ConvexHull(P)
    except Exception:
        return None
    return [[round(float(x), 5), round(float(y), 5)] for x, y in P[h.vertices]]


def _dnp_field(fields):
    """A 'DNP' / 'Do not populate' field with a true-ish value (the KiCad 6 convention, before the
    attribute existed); 'Populate = No' counts too."""
    for k, v in fields.items():
        kk, vv = k.strip().lower(), str(v).strip().lower()
        if kk in ("dnp", "do not populate", "do_not_populate", "nopop", "no_populate") and vv not in ("", "0", "no", "false", "~"):
            return True
        if kk in ("populate", "fitted", "mount") and vv in ("no", "0", "false", "n"):
            return True
    return False


def _backdrill(pcbnew, v, lid2name):
    """#26: back-drill of a via (KiCad 10): the layers up to which the barrel is drilled out, or None."""
    try:
        mode = int(v.GetBackdrillMode())
    except Exception:
        return None
    if mode == int(getattr(pcbnew, "BACKDRILL_MODE_NO_BACKDRILL", 0)):
        return None
    bot = lid2name.get(v.GetBottomBackdrillLayer()) if mode in (int(pcbnew.BACKDRILL_MODE_BACKDRILL_BOTTOM),
                                                                int(pcbnew.BACKDRILL_MODE_BACKDRILL_BOTH)) else None
    top = lid2name.get(v.GetTopBackdrillLayer()) if mode in (int(pcbnew.BACKDRILL_MODE_BACKDRILL_TOP),
                                                             int(pcbnew.BACKDRILL_MODE_BACKDRILL_BOTH)) else None
    return {"mode": mode, "bottom_to": bot, "top_to": top}


def _lid(pcbnew, board, name):
    return board.GetLayerID(name)


def _arc_circle(s, m, e):
    """[cx, cy, r] of the circle through an arc's start, mid and end; None for a straight one. Every
    chord of the arc carries it, so that distances to the arc itself can be exact (geom.arc_dist)."""
    (x1, y1), (x2, y2), (x3, y3) = s, m, e
    d = 2 * (x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2))
    if abs(d) < 1e-12:
        return None
    ux = ((x1 ** 2 + y1 ** 2) * (y2 - y3) + (x2 ** 2 + y2 ** 2) * (y3 - y1) + (x3 ** 2 + y3 ** 2) * (y1 - y2)) / d
    uy = ((x1 ** 2 + y1 ** 2) * (x3 - x2) + (x2 ** 2 + y2 ** 2) * (x1 - x3) + (x3 ** 2 + y3 ** 2) * (x2 - x1)) / d
    return [ux, uy, math.hypot(x1 - ux, y1 - uy)]


def _arc_points(s, m, e, max_deg=10.0):
    (x1, y1), (x2, y2), (x3, y3) = s, m, e
    d = 2 * (x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2))
    if abs(d) < 1e-12:
        return [s, e]
    ux = ((x1 ** 2 + y1 ** 2) * (y2 - y3) + (x2 ** 2 + y2 ** 2) * (y3 - y1) + (x3 ** 2 + y3 ** 2) * (y1 - y2)) / d
    uy = ((x1 ** 2 + y1 ** 2) * (x3 - x2) + (x2 ** 2 + y2 ** 2) * (x1 - x3) + (x3 ** 2 + y3 ** 2) * (x2 - x1)) / d
    r = math.hypot(x1 - ux, y1 - uy)
    a1, a2, a3 = (math.atan2(p[1] - uy, p[0] - ux) for p in (s, m, e))

    def norm(a):
        return (a + 2 * math.pi) % (2 * math.pi)
    sweep = norm(a3 - a1)
    if not (norm(a2 - a1) < sweep):
        sweep -= 2 * math.pi
    n = max(2, int(abs(math.degrees(sweep)) / max_deg) + 1)
    return [(ux + r * math.cos(a1 + sweep * k / n), uy + r * math.sin(a1 + sweep * k / n)) for k in range(n + 1)]


# ----------------------------------------------------------------------------------------------

class InputError(Exception):
    """The input is not something openpcbsi reads (reported without a traceback, exit code 3)."""


def file_version(board_path):
    """The `(version N)` of a .kicad_pcb file header, or None."""
    with open(board_path, encoding="utf-8", errors="replace") as f:
        head = f.read(400)
    m = re.search(r"\(version (\d+)\)", head)
    return int(m.group(1)) if m else None


def _cli_has_downgrade():
    try:
        p = subprocess.run(["kicad-cli", "pcb", "--help"], capture_output=True, text=True, timeout=60)
        return "downgrade" in (p.stdout + p.stderr)
    except (OSError, subprocess.TimeoutExpired):
        return False


def maybe_downgrade(board_path, workdir):
    """Input contract (#1): a KiCad 10.0 board file. A file from a newer KiCad is converted only when the
    `kicad-cli` on the PATH can downgrade it; otherwise a clear InputError says what to do."""
    import pcbnew
    if not os.path.exists(board_path):
        raise InputError(f"{board_path}: no such file")
    ver = file_version(board_path)
    if ver is None:
        raise InputError(f"{board_path} is not a KiCad board file (no '(version N)' header)")
    mine = int(getattr(pcbnew, "SEXPR_BOARD_FILE_VERSION", 0)) or None
    if mine is None or ver <= mine:
        try:
            if pcbnew.LoadBoard(board_path) is not None:
                return board_path
        except Exception as e:
            raise InputError(f"KiCad {pcbnew.Version()} could not read {board_path} (file version {ver}): {e}")
    if _cli_has_downgrade():
        os.makedirs(workdir, exist_ok=True)
        out = os.path.join(workdir, "downgraded_" + os.path.basename(board_path))
        p = subprocess.run(["kicad-cli", "pcb", "downgrade", "--target", "10.0", "-o", out, board_path],
                           capture_output=True, text=True)
        if p.returncode == 0 and os.path.exists(out):
            if (file_version(out) or 0) > FILE_VERSION_MAX:
                raise InputError(f"kicad-cli's downgrade of {os.path.basename(board_path)} is still file version {file_version(out)}")
            _check_constructs(out)
            print(f"openpcbsi: board file version {ver} downgraded to KiCad 10.0 -> {out}", file=sys.stderr)
            return out
    raise InputError(
        f"{board_path} has file version {ver}, newer than KiCad {pcbnew.Version()} reads (file version {mine}). "
        f"openpcbsi reads KiCad 10.0 boards: save the board as KiCad 10.0 in the tool that made it "
        f"(for example `kicad-cli pcb downgrade --target 10.0` of that KiCad) and run again.")


# The newest board-file version this release reads (KiCad 10.0). A newer file is downgraded with the
# installed kicad-cli when it can, else refused with a clear error (#1).
FILE_VERSION_MAX = 20260206
READERS = ("file", "pcbnew")


def reader_name(reader=None):
    """The board reader: "file" (openpcbsi's file reader on kicad-monkey; the default) or "pcbnew"
    (KiCad's SWIG module, deprecated by KiCad; kept while KiCad ships it). OPENPCBSI_READER chooses."""
    r = (reader or os.environ.get("OPENPCBSI_READER") or "file").lower()
    if r not in READERS:
        raise InputError(f"board reader {r!r}: one of {', '.join(READERS)}")
    return r


# constructs of KiCad files newer than 10.0 that the file reader does not read (a file holding one is
# refused, also after a downgrade that left it in: KiCad 10.99 places footprints with `(transform ...)`,
# which read as missing would put every footprint at (0, 0))
NEWER_CONSTRUCTS = (r"\(transform\b",)


def _check_constructs(path):
    with open(path, encoding="utf-8", errors="replace") as f:
        text = f.read()
    hit = [re.sub(r"[^a-z_]", "", c.replace("\\b", "")) for c in NEWER_CONSTRUCTS if re.search(c, text)]
    if hit:
        gen = re.search(r'\(generator_version "([^"]+)"\)', text[:2000])
        raise InputError(f"{os.path.basename(path)} uses {', '.join(hit)} (a KiCad newer than 10.0"
                         + (f"; generator {gen.group(1)}" if gen else "") + "), which this openpcbsi does not read. "
                         f"Save the board as KiCad 10.0 with KiCad 10.0 itself, and run again")


def _file_downgrade(board_path, workdir):
    """Input contract (#1) without pcbnew: a file newer than KiCad 10.0 is downgraded by kicad-cli when
    it can, else refused; a file (or a downgraded copy) that still holds a newer construct is refused."""
    ver = file_version(board_path)
    if ver is None:
        raise InputError(f"{board_path} is not a KiCad board file (no '(version N)' header)")
    if ver <= FILE_VERSION_MAX:
        _check_constructs(board_path)
        return board_path
    if _cli_has_downgrade():
        os.makedirs(workdir, exist_ok=True)
        out = os.path.join(workdir, "downgraded_" + os.path.basename(board_path))
        p = subprocess.run(["kicad-cli", "pcb", "downgrade", "--target", "10.0", "-o", out, board_path], capture_output=True, text=True)
        if p.returncode == 0 and os.path.exists(out):
            if (file_version(out) or 0) > FILE_VERSION_MAX:
                raise InputError(f"kicad-cli's downgrade of {os.path.basename(board_path)} is still file version {file_version(out)}")
            _check_constructs(out)
            print(f"openpcbsi: board file version {ver} downgraded to KiCad 10.0 -> {out}", file=sys.stderr)
            return out
    raise InputError(
        f"{board_path} has file version {ver}, newer than this openpcbsi reads (file version {FILE_VERSION_MAX}, KiCad 10.0). "
        f"Save the board as KiCad 10.0 in the tool that made it (for example `kicad-cli pcb downgrade --target 10.0` of that "
        f"KiCad) and run again.")


def extract_file(board_path, cache_dir):
    """The snapshot through openpcbsi's file reader: zones refilled by kicad-cli on a copy in the cache."""
    from .parse import monkey
    sha = file_sha256(board_path)
    refill = None if os.environ.get("OPENPCBSI_NO_REFILL") == "1" else os.path.join(cache_dir, "refill", sha[:16])
    snap = monkey.extract(board_path, refill_dir=refill)
    snap["version"] = SNAPSHOT_VERSION
    return snap


_CLI_VERSION = None


def kicad_cli_version():
    """The version of the kicad-cli on the PATH (it refills zones and downgrades files), or "none"."""
    global _CLI_VERSION
    if _CLI_VERSION is None:
        try:
            p = subprocess.run(["kicad-cli", "version"], capture_output=True, text=True, timeout=60)
            _CLI_VERSION = (p.stdout or p.stderr).strip().splitlines()[0] if p.returncode == 0 and (p.stdout or p.stderr).strip() else "unknown"
        except (OSError, subprocess.TimeoutExpired):
            _CLI_VERSION = "none"
    return _CLI_VERSION


def cache_tag(rd):
    """What a cached reading depends on besides the board's bytes: the reader, whether zones were
    refilled, and the kicad-cli that refilled or downgraded them. A run with other settings reads the
    board again instead of reusing (or poisoning) another's reading."""
    refill = "norefill" if os.environ.get("OPENPCBSI_NO_REFILL") == "1" else "refill"
    cli = re.sub(r"[^0-9A-Za-z.]+", "-", kicad_cli_version())[:40]
    return f"{rd}_{refill}_cli{cli}"


def _check_placement(snap):
    """A reading in which most footprints sit far from their own pads (or all at the origin) is wrong:
    refused, not reviewed. (A construct the reader does not know can drop positions silently.)"""
    fps = snap.get("footprints") or []
    pads = {}
    for p in snap.get("pads") or []:
        pads.setdefault(p["ref"], []).append(p)
    far = origin = n = 0
    for f in fps:
        ps = pads.get(f["ref"])
        if not ps:
            continue
        n += 1
        cx, cy = sum(p["x"] for p in ps) / len(ps), sum(p["y"] for p in ps) / len(ps)
        far += math.hypot(cx - f["x"], cy - f["y"]) > 50.0
        origin += abs(f["x"]) < 1e-9 and abs(f["y"]) < 1e-9
    if n >= 4 and (far > n / 2 or origin > n / 2):
        raise InputError(f"the reading places {max(far, origin)} of {n} footprints far from their own pads or at the origin: "
                         f"the board uses something the reader does not read (a newer KiCad?); not reviewed")


def load_or_extract(board_path, cache_dir, refresh=False, reader=None):
    """Return the snapshot dict, reading the board only if the cache is stale. Parallel runs on one
    board take a per-board lock, so the board is read once and the others reuse it (#33)."""
    import fcntl
    board_path = os.path.abspath(board_path)
    if not os.path.exists(board_path):
        raise InputError(f"{board_path}: no such file")
    rd = reader_name(reader)
    sha = file_sha256(board_path)
    os.makedirs(cache_dir, exist_ok=True)
    cpath = os.path.join(cache_dir, f"snapshot_{sha[:16]}_{cache_tag(rd)}.json.gz")
    with open(cpath + ".lock", "a") as lk:
        fcntl.flock(lk, fcntl.LOCK_EX)
        try:
            if os.path.exists(cpath) and not refresh:
                with gzip.open(cpath, "rt") as f:
                    snap = json.load(f)
                if snap.get("version") == SNAPSHOT_VERSION and snap.get("cache_tag") == cache_tag(rd):
                    snap["board"] = board_path
                    return snap
            if rd == "pcbnew":
                src = maybe_downgrade(board_path, os.path.join(cache_dir, "downgrade"))
                snap = extract(src)
                snap["reader"] = "pcbnew"
            else:
                src = _file_downgrade(board_path, os.path.join(cache_dir, "downgrade"))
                snap = extract_file(src, cache_dir)
            snap["board"] = board_path
            snap["sha256"] = sha
            snap["downgraded"] = src if src != board_path else None   # G22: the report records it
            snap["kicad_file_version"] = file_version(board_path)
            snap["cache_tag"] = cache_tag(rd)
            _check_placement(snap)
            tmp = f"{cpath}.{os.getpid()}.tmp"
            with gzip.open(tmp, "wt") as f:
                json.dump(snap, f)
            os.replace(tmp, cpath)
            return snap
        finally:
            fcntl.flock(lk, fcntl.LOCK_UN)


if __name__ == "__main__":
    s = load_or_extract(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else "cache", refresh=True)
    print(len(s["tracks"]), "tracks", len(s["vias"]), "vias", len(s["zones"]), "zones", len(s["pads"]), "pads")


def cross_check(snap, cache_dir):
    """The board read again by an independent parser (KiCadFiles) and compared with `snap` on the
    fields that parser reads (openpcbsi.parse.kicadfiles). Cached by the board's SHA-256.
    Returns {"reader", "compared", "problems", "not_cross_checked"}."""
    sha = snap["sha256"]
    cp = os.path.join(cache_dir, f"crosscheck_{sha[:16]}_v{SNAPSHOT_VERSION}_{snap.get('cache_tag', 'x')}.json")
    if os.path.exists(cp):
        try:
            return json.load(open(cp))
        except ValueError:
            pass
    from .parse import kicadfiles
    board = snap["board"]
    refilled = os.path.join(cache_dir, "refill", sha[:16], os.path.basename(board))
    d = kicadfiles.cross_check(snap, refilled if os.path.exists(refilled) else board)
    out = {"reader": "kicadfiles (KiCadFiles " + __import__("openpcbsi.parse", fromlist=["x"]).lock()["kicadfiles"]["version"] + ")",
           "compared": d["compared"], "counts": {k: list(v) for k, v in d["counts"].items()}, "problems": d["problems"][:50],
           "n_problems": len(d["problems"]), "not_cross_checked": d["not_cross_checked"]}
    tmp = f"{cp}.{os.getpid()}.tmp"
    json.dump(out, open(tmp, "w"))
    os.replace(tmp, cp)
    return out
