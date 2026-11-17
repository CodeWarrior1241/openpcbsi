"""Body size, heat capacity and conductance of parts, for the reflow checks (MF7, MF8).

KiCad carries no body height or mass. Chip parts (two terminals, EIA size in the footprint name) come
from the tables below; every other part is a block over its courtyard with an assumed height, unless
the board config (`[reflow.parts.REF]`) or a footprint field (`Mass`, `Height`) says otherwise. All
of these are approximate and the report states which source each part used.
"""
import math
import re

import numpy as np

# EIA code -> body length, width, height in mm (height: multilayer ceramic capacitor, nominal)
CHIP = {"01005": (0.4, 0.2, 0.2), "0201": (0.6, 0.3, 0.3), "0402": (1.0, 0.5, 0.5), "0603": (1.6, 0.8, 0.8),
        "0805": (2.0, 1.25, 0.85), "1206": (3.2, 1.6, 1.0), "1210": (3.2, 2.5, 1.5), "1812": (4.5, 3.2, 1.5),
        "2010": (5.0, 2.5, 0.55), "2512": (6.3, 3.2, 0.55)}
# thick-film chip resistors are thinner than capacitors of the same size
RES_H = {"01005": 0.13, "0201": 0.23, "0402": 0.35, "0603": 0.45, "0805": 0.5, "1206": 0.55, "1210": 0.55,
         "1812": 0.55, "2010": 0.55, "2512": 0.55}
ORDER = list(CHIP)
# density kg/m3, specific heat J/kgK, conductivity W/mK (along the body, terminal to terminal)
MATERIAL = {
    "mlcc": (6000.0, 500.0, 10.0),       # barium titanate with nickel electrodes
    "resistor": (3900.0, 800.0, 25.0),   # alumina substrate
    "ferrite": (5000.0, 700.0, 4.0),     # chip inductors and beads
    "moulded": (1900.0, 900.0, 0.8),     # epoxy-moulded parts (diodes, LEDs, ICs)
}
KIND_BY_NAME = (("C_", "mlcc"), ("CP_", "mlcc"), ("R_", "resistor"), ("L_", "ferrite"), ("FB_", "ferrite"),
                ("Fuse_", "resistor"), ("LED_", "moulded"), ("D_", "moulded"))
KIND_BY_REF = (("C", "mlcc"), ("R", "resistor"), ("L", "ferrite"), ("FB", "ferrite"), ("F", "resistor"),
               ("D", "moulded"), ("LED", "moulded"))
_SIZE = re.compile(r"(?:^|_)(01005|0201|0402|0603|0805|1206|1210|1812|2010|2512)(?:_|$)")


def _kind(fp):
    name = (fp.get("fpid") or "").split(":")[-1]
    for pre, k in KIND_BY_NAME:
        if name.startswith(pre):
            return k
    ref = re.match(r"[A-Za-z]+", fp.get("ref") or "")
    for pre, k in sorted(KIND_BY_REF, key=lambda x: -len(x[0])):
        if ref and ref.group(0).upper() == pre:
            return k
    return "moulded"


def chip_of(fp):
    """The chip package of a footprint (by the EIA size in its name), or None."""
    name = (fp.get("fpid") or "").split(":")[-1]
    m = _SIZE.search(name)
    if not m:
        return None
    code = m.group(1)
    L, W, H = CHIP[code]
    kind = _kind(fp)
    if kind == "resistor":
        H = RES_H[code]
    rho, cp, k = MATERIAL[kind]
    vol = L * W * H * 1e-9
    return {"code": code, "kind": kind, "L": L, "W": W, "H": H, "mass_mg": round(rho * vol * 1e6, 3),
            "C": rho * cp * vol, "k": k, "source": f"chip table ({code} {kind})"}


def size_rank(code):
    return ORDER.index(code) if code in ORDER else len(ORDER)


def _num(v):
    m = re.search(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", str(v))
    return float(m.group(0)) if m else None


def body_of(fp, cfg):
    """Heat capacity (J/K) and height of any part: chip table, else the board config, else footprint
    fields, else a block over the courtyard with the assumed height and moulded-epoxy properties."""
    r = cfg["reflow"]
    spec = r.get("parts", {}).get(fp["ref"], {})
    chip = chip_of(fp)
    rho, cp, _ = MATERIAL["moulded"]
    x0, y0, x1, y1 = fp["bbox"]
    area = abs(x1 - x0) * abs(y1 - y0)
    rings = [ring for poly in fp.get("courtyard", []) for ring in poly]
    if rings:
        a = sum(abs(_area(rg)) for rg in rings)
        if a > 0:
            area = min(area, a) if area > 0 else a
    height, hsrc = spec.get("height_mm"), "config"
    if height is None:
        fv = {k.lower(): v for k, v in fp.get("fields", {}).items()}
        height, hsrc = (_num(fv["height"]), "field Height") if "height" in fv and _num(fv["height"]) else (None, None)
    mass_g, msrc = spec.get("mass_g"), "config"
    if mass_g is None:
        fv = {k.lower(): v for k, v in fp.get("fields", {}).items()}
        for key in ("mass", "weight"):
            if key in fv and _num(fv[key]):
                mass_g, msrc = _num(fv[key]), f"field {key.capitalize()} (g)"
                break
    if mass_g is not None:
        return {"C": mass_g * 1e-3 * spec.get("cp", cp), "height": height or r["default_height_mm"],
                "area": area, "source": f"mass {mass_g:g} g ({msrc})", "assumed": False}
    if chip and height is None:
        return {"C": chip["C"], "height": chip["H"], "area": chip["L"] * chip["W"], "source": chip["source"],
                "assumed": False}
    if height is None:
        height, hsrc = r["default_height_mm"], "assumed height"
    return {"C": rho * cp * area * height * 1e-9, "height": height, "area": area,
            "source": f"block {area:.1f} mm2 x {height:g} mm ({hsrc}), moulded epoxy", "assumed": hsrc == "assumed height"}


def _area(ring):
    a = np.asarray(ring, dtype=float)
    x, y = a[:, 0], a[:, 1]
    return 0.5 * float(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))


PART_CLASSES = (("BGA", "BGA"), ("CSP", "BGA"), ("LGA", "BGA"), ("QFN", "QFN"), ("DFN", "QFN"), ("QFP", "QFP"), ("SOIC", "SOP"),
                ("SSOP", "SOP"), ("TSSOP", "SOP"), ("MSOP", "SOP"), ("SOP", "SOP"), ("SOT", "SOT"), ("SC-70", "SOT"), ("TSOT", "SOT"))


# KiCad's standard library nicknames by package class
LIB_CLASSES = (("PACKAGE_BGA", "BGA"), ("PACKAGE_CSP", "BGA"), ("PACKAGE_DFN_QFN", "QFN"), ("PACKAGE_QFP", "QFP"),
               ("PACKAGE_SO", "SOP"), ("PACKAGE_TO_SOT_SMD", "SOT"))


def part_class_how(fp, pads=None):
    """(class, how it was decided) for the fabricator's part-spacing table: an EIA chip size, or QFN / QFP /
    SOP / SOT / BGA. Decided from the footprint's EIA size, its KiCad library, its name, then its pads;
    (None, why) when none of these places it."""
    chip = chip_of(fp)
    if chip:
        return chip["code"], "EIA size in the footprint name"
    lib, _, name = (fp.get("fpid") or "").partition(":")
    if not name:
        lib, name = "", lib
    for key, cls in LIB_CLASSES:
        if lib.upper().startswith(key):
            return cls, f"KiCad library {lib}"
    if lib.upper().startswith("CONNECTOR") or (re.match(r"(J|P|CN|X)\d", fp.get("ref") or "") and not chip_of(fp)):
        return None, "a connector: the assembler's table has no class for connectors"
    text = name.upper() + " " + str(fp.get("value") or "").upper()
    for key, cls in PART_CLASSES:
        if key in text:
            return cls, f"'{key}' in the footprint name"
    if pads:
        smd = [p for p in pads if p["kind"] == "smd" and p.get("layers") and str(p["num"]).strip()]
        nums = {str(p["num"]) for p in smd}
        if len(nums) >= 9 and all(p.get("shape") == "circle" for p in smd):
            return "BGA", f"{len(nums)} round pads"
        if len(nums) == 3:
            return "SOT", "three pads"
        if len(nums) >= 4:
            xy = np.array([(p["x"], p["y"]) for p in smd])
            cx, cy = xy.mean(axis=0)
            rot = math.radians(fp.get("rot") or 0.0)
            u = (xy[:, 0] - cx) * math.cos(rot) - (xy[:, 1] - cy) * math.sin(rot)
            v = (xy[:, 0] - cx) * math.sin(rot) + (xy[:, 1] - cy) * math.cos(rot)
            rows = len(set(np.round(v, 2))) <= 3 and len(set(np.round(u, 2))) > 3
            cols = len(set(np.round(u, 2))) <= 3 and len(set(np.round(v, 2))) > 3
            if rows or cols:
                return "SOP", "pins on two sides"
            return None, f"{len(nums)} pins on more than two sides: QFN or QFP cannot be told from the pads alone"
        return None, f"{len(nums)} numbered SMD pads and no package class in its name or library"
    return None, "no package class in its name or library"


def part_class(fp, pads=None):
    return part_class_how(fp, pads)[0]


def land_pattern(chip, pads):
    """Geometry of a two-pad land pattern against the nominal body of the chip it is named for:
    gap between the pads, their outer span along the part, toe (pad beyond the body end), pad width
    over body width, and the larger pad's area over the smaller's."""
    from . import geom
    a = [np.vstack([np.asarray(r, dtype=float) for poly in p["poly"] for r in poly]) for p in pads]
    cx = [(p["x"], p["y"]) for p in pads]
    ux, uy = cx[1][0] - cx[0][0], cx[1][1] - cx[0][1]
    n = (ux * ux + uy * uy) ** 0.5 or 1.0
    ux, uy = ux / n, uy / n
    along = [q[:, 0] * ux + q[:, 1] * uy for q in a]
    across = [-q[:, 0] * uy + q[:, 1] * ux for q in a]
    span = max(al.max() for al in along) - min(al.min() for al in along)
    width = min(ac.max() - ac.min() for ac in across)
    gap = geom.rings_distance([r for poly in pads[0]["poly"] for r in poly], [r for poly in pads[1]["poly"] for r in poly])
    areas = [sum(abs(_area(r)) for poly in p["poly"] for r in poly) for p in pads]
    return {"gap_mm": round(gap, 3), "span_mm": round(span, 3), "toe_mm": round((span - chip["L"]) / 2, 3),
            "width_ratio": round(width / chip["W"], 2), "area_ratio": round(max(areas) / max(min(areas), 1e-9), 2),
            "body": f"{chip['L']:g} x {chip['W']:g} mm"}


def pin_pitch(pads, min_pins=4):
    """The smallest distance between two distinct pin positions of a part's SMD pads (a pin drawn as
    several pads counts once), or None for fewer than `min_pins` pins."""
    pins = {}
    for p in pads:
        if p["kind"] == "smd" and p.get("layers") and str(p["num"]).strip():
            pins.setdefault(str(p["num"]), (p["x"], p["y"]))
    if len(pins) < min_pins:
        return None
    from scipy.spatial import cKDTree
    xy = np.array(list(pins.values()))
    d = cKDTree(xy).query(xy, k=2)[0][:, 1]
    d = d[d > 1e-6]
    return float(d.min()) if len(d) else None


def stencil_prediction(b, guide):
    """The stencil thickness a fabricator's guide gives for this board: per placed SMD part, the row for
    its chip size or its pin / ball pitch; the board's finest part decides. Returns {thickness_mm, part,
    row, uncovered: [refs the guide does not cover]}."""
    rows = guide.get("rows", [])
    larger = guide.get("larger_mm")
    best, uncovered = None, []
    for ref, fp in b.fps.items():
        pads = [p for p in b.pads_by_ref.get(ref, []) if p["kind"] == "smd" and p.get("layers")]
        if fp.get("dnp") or not fp.get("placed", True) or not pads:
            continue
        chip = chip_of(fp)
        t, row = None, None
        if chip:
            r = next((r for r in rows if r.get("kind") == "chip" and r.get("size") == chip["code"]), None)
            if r:
                t, row = r["thickness_mm"], f"{chip['code']} chips"
            elif size_rank(chip["code"]) < min((size_rank(r["size"]) for r in rows if r.get("kind") == "chip"), default=0):
                uncovered.append(f"{ref} ({chip['code']})")
                continue
            else:
                t, row = larger, f"{chip['code']} chips (larger parts)"
        else:
            pitch = pin_pitch(pads)
            if pitch is None:
                t, row = larger, "parts with fewer than 4 pins (larger parts)"
            else:
                kind = "bga" if part_class(fp) == "BGA" else "ic"
                cand = sorted((r for r in rows if r.get("kind") == kind), key=lambda r: r["pitch_mm"])
                fit = [r for r in cand if r["pitch_mm"] <= pitch + 1e-3]
                if not fit:
                    uncovered.append(f"{ref} ({kind.upper()} pitch {pitch:.2f} mm)")
                    continue
                r = fit[-1]
                t, row = r["thickness_mm"], f"{kind.upper()} pitch {r['pitch_mm']:g} mm (this part: {pitch:.2f} mm)"
        if t is not None and (best is None or t < best[0]):
            best = (t, ref, row)
    return {"thickness_mm": None if (uncovered or best is None) else best[0], "part": best[1] if best else None,
            "row": best[2] if best else None, "uncovered": uncovered}


NOT_PARTS = ("FIDUCIAL", "TESTPOINT", "TEST_POINT", "MOUNTINGHOLE", "MOUNTING_HOLE")
NOT_PART_REFS = ("FID", "TP", "MH", "H")


def not_a_part(fp, pads):
    """Why a footprint is not a part that is soldered in reflow (a fiducial, a test point, a mounting hole,
    a board-only or BOM-excluded footprint, or one with no paste), else None."""
    if fp.get("board_only"):
        return "board-only footprint"
    if fp.get("in_bom") is False:
        return "excluded from the BOM"
    name = (fp.get("fpid") or "").upper()
    if any(k in name for k in NOT_PARTS):
        return "fiducial, test point or mounting hole (footprint)"
    m = re.match(r"[A-Za-z]+", fp.get("ref") or "")
    if m and m.group(0).upper() in NOT_PART_REFS:
        return f"reference {m.group(0).upper()}"
    if not any(p.get("paste") for p in pads):
        return "no paste on any pad"
    return None
