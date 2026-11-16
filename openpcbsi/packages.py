"""Body size, heat capacity and conductance of parts, for the reflow checks (MF7, MF8).

KiCad carries no body height or mass. Chip parts (two terminals, EIA size in the footprint name) come
from the tables below; every other part is a block over its courtyard with an assumed height, unless
the board config (`[reflow.parts.REF]`) or a footprint field (`Mass`, `Height`) says otherwise. All
of these are approximate and the report states which source each part used.
"""
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
