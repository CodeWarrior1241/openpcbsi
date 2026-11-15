"""G9: via model library. Every distinct via geometry used by nets in scope (drill, pad, antipads
measured on each plane it passes, connected layers, stub, nearest return via) gets a lumped model,
computed once and cached by geometry:

  C = 2 pi eps l / ln(r_antipad / r_barrel)      (barrel inside the plane antipads, inner pads removed)
  L = mu0 l / (2 pi) ln(s / r_barrel)            (return through the nearest ground via at distance s)
  Z = sqrt(L / C), delay = sqrt(L C)

Closed forms, not field-solved; G8 (openEMS) is the route to S-parameter models of critical vias."""
import hashlib
import json
import math
import os

import numpy as np

MU0 = 4e-7 * math.pi
IN = 25.4


def antipad_radius(ctx, li, x, y, r0, rmax=1.5):
    """Smallest clearance radius around (x, y) on plane layer li (rays in 16 directions)."""
    g = ctx.g
    ref = ctx.refmask(li)
    ang = np.linspace(0, 2 * np.pi, 16, endpoint=False)
    rs = np.arange(r0, rmax, g.res / 2)
    X = x + np.outer(rs, np.cos(ang))
    Y = y + np.outer(rs, np.sin(ang))
    hit = g.sample(ref.view(np.uint8), X, Y).astype(bool)
    first = np.where(hit.any(axis=0), np.argmax(hit, axis=0), len(rs) - 1)
    rr = rs[first]
    return float(rr.min()), float(np.median(rr))


def model(geom_, er):
    """Coaxial barrel model: the barrel (drill + plating; inner pads removed) inside the antipads of
    the planes it passes, return through the nearest ground via. Outer pads add a parallel-plate
    term over the first dielectric only where they are larger than the barrel."""
    l = geom_["barrel_mm"] * 1e-3
    rb = geom_["drill_mm"] / 2 * 1e-3 + 20e-6
    ra = max(geom_["antipad_mm"] / 2 * 1e-3, rb * 1.2)
    C = 2 * math.pi * 8.8541878128e-12 * er * l / math.log(ra / rb)
    s = geom_["return_mm"] * 1e-3 if geom_["return_mm"] else 2 * ra
    L = MU0 * l / (2 * math.pi) * math.log(max(s, 2 * rb) / rb)
    out = {"C_pF": round(C * 1e12, 4), "L_nH": round(L * 1e9, 4), "Lloop_nH": round(2 * L * 1e9, 4),
           "Z_ohm": round(math.sqrt(L / C), 1) if C > 0 else None,
           "delay_ps": round(math.sqrt(L * C) * 1e12, 3) if C > 0 else None,
           "model": "coaxial barrel (Johnson & Graham style), closed form"}
    return out


def run(ctx, cache_dir, log=print):
    b, cfg = ctx.b, ctx.cfg
    st = b.stack
    vtl, vzl = ctx.via_track_layers(), ctx.via_zone_layers()
    lib_path = os.path.join(cache_dir, "via_library.json")
    lib = json.load(open(lib_path)) if os.path.exists(lib_path) else {}
    V = b.vias
    gv = np.flatnonzero(np.isin(V["net"], ctx.ref_ids))
    planes = [li for li in range(len(b.layers)) if ctx.refmask(li).mean() > 0.3]
    rows, inst = [], {}
    nets = set(ctx.nets("via_stub"))
    for k in range(len(V["x"])):
        net = b.net_names[V["net"][k]]
        if net not in nets:
            continue
        conn = sorted(vtl[k] | vzl[k])
        if len(conn) < 2:
            continue
        x, y = V["x"][k], V["y"][k]
        aps = []
        for li in planes:
            if li in conn:
                continue
            aps.append(antipad_radius(ctx, li, x, y, V["d"][k] / 2)[0] * 2)
        ap = round(min(aps), 2) if aps else round(V["d"][k] * 1.5, 2)
        d = np.hypot(V["x"][gv] - x, V["y"][gv] - y) if len(gv) else np.array([np.inf])
        s = round(float(d.min()), 1)
        lo, hi = min(conn), max(conn)
        barrel = st.barrel(b.layers[lo], b.layers[hi])
        stub = max(st.barrel(b.layers[int(V["top"][k])], b.layers[lo]), st.barrel(b.layers[hi], b.layers[int(V["bot"][k])]))
        g = {"drill_mm": round(V["drill"][k], 3), "pad_mm": round(V["d"][k], 3), "antipad_mm": ap,
             "barrel_mm": round(barrel, 4), "stub_mm": round(stub, 4), "return_mm": s,
             "layers": f"{b.layers[lo]}-{b.layers[hi]}"}
        key = hashlib.sha1(json.dumps(g, sort_keys=True).encode()).hexdigest()[:12]
        if key not in lib:
            lib[key] = dict(g, **model(g, st.mean_er()))
        inst.setdefault(key, []).append(net)
    for key, nets_ in inst.items():
        r = dict(lib[key], key=key, instances=len(nets_), nets=sorted(set(nets_))[:6])
        rows.append(r)
    rows.sort(key=lambda r: -r["instances"])
    with open(lib_path, "w") as f:
        json.dump(lib, f, indent=1)
    for r in rows[:8]:
        log(f"  {r['layers']:15s} drill {r['drill_mm']} pad {r['pad_mm']} antipad {r['antipad_mm']} return {r['return_mm']} mm"
            f"  C {r['C_pF']} pF  L {r['L_nH']} nH  Z {r.get('Z_ohm')} ohm  x{r['instances']}")
    return {"rule": "vias", "spec": "G9", "title": "Via model library (closed form)", "validation": "closed-form",
            "rows": rows, "metrics": {f"{r['layers']} d{r['drill_mm']} ap{r['antipad_mm']}": {"C_pF": r["C_pF"], "L_nH": r["L_nH"]} for r in rows},
            "findings": [], "summary": {"geometries": len(rows)},
            "note": "Johnson & Graham closed forms; G8 (openEMS) is the route to S-parameter via models."}


def lookup(vias_report, layers):
    for r in (vias_report or {}).get("rows", []):
        if r["layers"] == layers:
            return r
    return None
