"""Board snapshot: everything openpcbsi needs from a KiCad 10 board, read once through pcbnew.

The snapshot is plain JSON (gzip), cached by the SHA-256 of the board file, so that every later
step (and the openEMS venv, which has no pcbnew) works from the same frozen data. Units: mm,
KiCad board coordinates (y grows downwards).

Only this module imports pcbnew. Run it with the system python3 of stock KiCad 10.
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

SNAPSHOT_VERSION = 8   # 4: zones refilled in memory, every footprint field, pin types, file version; 5: back-drill;
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

        def val(key, cast=float):
            mm = re.search(r'\(%s "?([^")]+)"?\)' % key, body)
            return cast(mm.group(1)) if mm else None

        out.append({
            "name": name,
            "type": val("type", str),
            "thickness": val("thickness"),
            "epsilon_r": val("epsilon_r"),
            "loss_tangent": val("loss_tangent"),
            "material": val("material", str),
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
        pts.append((round(p.x * NM, 5), round(p.y * NM, 5)))
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
                "mask_exp": _via_mask_exp(pcbnew, t),
            })
        elif ty == pcbnew.PCB_ARC_T:
            # approximate arcs by chords of at most 10 degrees
            s, m, e = t.GetStart(), t.GetMid(), t.GetEnd()
            pts = _arc_points((s.x * NM, s.y * NM), (m.x * NM, m.y * NM), (e.x * NM, e.y * NM))
            for (x1, y1), (x2, y2) in zip(pts[:-1], pts[1:]):
                tracks.append({"x1": x1, "y1": y1, "x2": x2, "y2": y2, "w": t.GetWidth() * NM,
                               "layer": lid2name.get(t.GetLayer()), "net": net, "arc": True})
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
                          "name": z.GetZoneName(), "polys": _polyset(fp)})
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
        fps.append({
            "ref": ref, "value": f.GetValue(), "fpid": f.GetFPIDAsString(),
            "x": f.GetPosition().x * NM, "y": f.GetPosition().y * NM,
            "rot": f.GetOrientationDegrees(), "side": "B" if f.IsFlipped() else "F",
            "dnp": bool(f.IsDNP()) or _dnp_field(fields), "fields": fields, "courtyard": crt,
            "smd": bool(f.GetAttributes() & pcbnew.FP_SMD),
            "placed": not (f.GetAttributes() & (pcbnew.FP_EXCLUDE_FROM_POS_FILES | pcbnew.FP_BOARD_ONLY)),
            "bbox": [bb.GetX() * NM, bb.GetY() * NM, bb.GetRight() * NM, bb.GetBottom() * NM],
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
            sz = p.GetSize(lid0)
            pads.append({
                "ref": ref, "num": p.GetNumber(), "net": p.GetNetname(), "kind": kind,
                "x": p.GetPosition().x * NM, "y": p.GetPosition().y * NM,
                "sx": sz.x * NM, "sy": sz.y * NM, "rot": p.GetOrientationDegrees(),
                "drill": p.GetDrillSize().x * NM if kind in ("pth", "npth") else 0.0,
                "drill_y": p.GetDrillSize().y * NM if kind in ("pth", "npth") else 0.0,
                "layers": layers, "poly": polys, "pinfunc": p.GetPinFunction(), "pintype": p.GetPinType(),
                **_pad_mask_paste(pcbnew, p),
            })
    snap["footprints"] = fps
    snap["pads"] = pads
    snap["mask_shapes"] = _mask_shapes(pcbnew, b)
    return snap


MASK_PASTE = (("F.Mask", "F_Mask", "F.Paste", "F_Paste"), ("B.Mask", "B_Mask", "B.Paste", "B_Paste"))


def _pad_mask_paste(pcbnew, p):
    """Solder-mask opening (the pad shape grown by its mask expansion) and paste margin of a pad, per
    side. `mask` is {mask layer: polygons}; a pad that is not on a mask layer is covered by mask."""
    mask, paste = {}, {}
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
    return {"mask": mask, "paste": paste}


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
            print(f"openpcbsi: board file version {ver} downgraded to KiCad 10.0 -> {out}", file=sys.stderr)
            return out
    raise InputError(
        f"{board_path} has file version {ver}, newer than KiCad {pcbnew.Version()} reads (file version {mine}). "
        f"openpcbsi reads KiCad 10.0 boards: save the board as KiCad 10.0 in the tool that made it "
        f"(for example `kicad-cli pcb downgrade --target 10.0` of that KiCad) and run again.")


def load_or_extract(board_path, cache_dir, refresh=False):
    """Return the snapshot dict, extracting through pcbnew only if the cache is stale. Parallel runs on
    one board take a per-board lock, so the board is extracted once and the others reuse it (#33)."""
    import fcntl
    board_path = os.path.abspath(board_path)
    if not os.path.exists(board_path):
        raise InputError(f"{board_path}: no such file")
    sha = file_sha256(board_path)
    os.makedirs(cache_dir, exist_ok=True)
    cpath = os.path.join(cache_dir, f"snapshot_{sha[:16]}.json.gz")
    with open(cpath + ".lock", "a") as lk:
        fcntl.flock(lk, fcntl.LOCK_EX)
        try:
            if os.path.exists(cpath) and not refresh:
                with gzip.open(cpath, "rt") as f:
                    snap = json.load(f)
                if snap.get("version") == SNAPSHOT_VERSION:
                    snap["board"] = board_path
                    return snap
            src = maybe_downgrade(board_path, os.path.join(cache_dir, "downgrade"))
            snap = extract(src)
            snap["board"] = board_path
            snap["sha256"] = sha
            snap["downgraded"] = src if src != board_path else None   # G22: the report records it
            snap["kicad_file_version"] = file_version(board_path)
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
