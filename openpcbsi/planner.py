"""G6: stackup and impedance planner. Impedance of every track geometry used on every layer
(single-ended widths and pair width/gap), compared with class targets and fabricator figures."""
import math
from collections import defaultdict

import numpy as np

from . import xsec
from .validate import status


def pair_geometries(ctx):
    """(layer, w, gap) -> length, from parallel P/N segment overlaps (exact vector geometry)."""
    b = ctx.b
    T = b.tracks
    out = defaultdict(float)
    seen = set()
    for pn, nn in ctx.partner.items():
        if pn in seen:
            continue
        seen |= {pn, nn}
        kp, kn = b.tracks_of(pn), b.tracks_of(nn)
        for i in kp:
            ax, ay = T["x2"][i] - T["x1"][i], T["y2"][i] - T["y1"][i]
            La = math.hypot(ax, ay)
            if La < 0.2:
                continue
            ux, uy = ax / La, ay / La
            for j in kn:
                if T["layer"][j] != T["layer"][i]:
                    continue
                bx, by = T["x2"][j] - T["x1"][j], T["y2"][j] - T["y1"][j]
                Lb = math.hypot(bx, by)
                if Lb < 0.2 or abs(ux * by - uy * bx) / Lb > 0.01:
                    continue
                # perpendicular distance and overlap along the direction
                d = abs((T["x1"][j] - T["x1"][i]) * uy - (T["y1"][j] - T["y1"][i]) * ux)
                s1 = (T["x1"][j] - T["x1"][i]) * ux + (T["y1"][j] - T["y1"][i]) * uy
                s2 = (T["x2"][j] - T["x1"][i]) * ux + (T["y2"][j] - T["y1"][i]) * uy
                ov = min(La, max(s1, s2)) - max(0.0, min(s1, s2))
                gap = d - (T["w"][i] + T["w"][j]) / 2
                if ov > 0.1 and 0 < gap < 0.5:
                    out[(int(T["layer"][i]), round(float(T["w"][i]), 4), round(gap, 3))] += ov
    return out


def run(ctx, cache_dir, log=print):
    b, cfg = ctx.b, ctx.cfg
    xc = cfg["xsec"]
    freqs = [f * 1e9 for f in xc["freqs_ghz"]]
    T = b.tracks
    rows = []
    # single-ended widths per layer (signal nets only; >= 1 mm in total)
    se = defaultdict(float)
    for k in range(len(T["x1"])):
        net = b.net_names[T["net"][k]]
        if net in cfg.reference_nets or b.kicad_class.get(net) == "PWR":
            continue
        se[(int(T["layer"][k]), round(float(T["w"][k]), 4))] += float(T["len"][k])
    pg = pair_geometries(ctx)
    # merge pair geometries closer than 2 um in gap
    pairs = defaultdict(float)
    for (li, w, g), L in pg.items():
        pairs[(li, w, round(g / 0.002) * 0.002)] += L
    jobs = [("se", li, w, None, L) for (li, w), L in se.items() if L >= 1.0] + \
           [("pair", li, w, g, L) for (li, w, g), L in pairs.items() if L >= 1.0]
    jobs.sort(key=lambda j: (j[1], j[0], j[2], j[3] or 0))
    fab = {(f["layer"], round(f["w"], 4), round(f.get("gap", 0) or 0, 3)): f for f in xc["fab"]}
    for kind, li, w, g, L in jobs:
        lname = b.layers[li]
        refs = [b.layers[r] for r in ctx.ref_layers(li)]
        try:
            sec = xsec.build_section(b.stack, lname, refs, [w] if kind == "se" else [w, w],
                                     [] if kind == "se" else [g], xc, xc["mask_er"])
            r = xsec.solve_lines(sec, xc["res_um"], freqs, xc["roughness_um"])
        except Exception as e:
            log(f"  {lname} {kind} w={w}: {e}")
            continue
        f0 = r["freqs"][0]
        row = {"layer": lname, "kind": kind, "w": w, "gap": g, "length_mm": round(L, 1), "refs": refs,
               "nodes": r["nodes"]}
        if kind == "se":
            row.update({"Z0": round(f0["Z0"], 2), "eeff": round(f0["eeff"], 3), "delay_ps_mm": round(f0["delay_ps_mm"], 3)})
            row["loss_db_mm"] = {f"{fr['f_hz']/1e9:g}GHz": round(fr["loss_db_mm"], 5) for fr in r["freqs"]}
        else:
            row.update({"Zdiff": round(f0["Zdiff"], 2), "Zodd": round(f0["Zodd"], 2), "Zeven": round(f0["Zeven"], 2),
                        "Zcm": round(f0["Zcm"], 2), "eeff_odd": round(f0["eeff_odd"], 3),
                        "delay_ps_mm": round(f0["delay_diff_ps_mm"], 3), "Kb": round(f0["Kb"], 4)})
            row["loss_db_mm"] = {f"{fr['f_hz']/1e9:g}GHz": round(fr["loss_diff_db_mm"], 5) for fr in r["freqs"]}
            row["Zdiff_by_f"] = {f"{fr['f_hz']/1e9:g}GHz": round(fr["Zdiff"], 2) for fr in r["freqs"]}
            tgt = None
            for cname, c in cfg["classes"].items():
                if c.get("z_diff"):
                    tgt = (c["z_diff"], c.get("z_tol", 10.0))
            if tgt:
                row["target"] = tgt[0]
                row["delta_target"] = round(row["Zdiff"] - tgt[0], 2)
                row["within_target"] = abs(row["Zdiff"] - tgt[0]) <= tgt[1]
        fk = (lname, round(w, 4), round(g or 0, 3))
        if fk in fab:
            row["fab_z"] = fab[fk]["z"]
            row["delta_fab"] = round((row.get("Zdiff") or row.get("Z0")) - fab[fk]["z"], 2)
        rows.append(row)
        log(f"  {lname:7s} {kind:4s} w={w:.3f}" + (f" gap={g:.3f}" if g else "          ")
            + (f"  Zdiff={row['Zdiff']:6.2f}" if kind == "pair" else f"  Z0={row['Z0']:6.2f}   ")
            + f"  {row['delay_ps_mm']:.2f} ps/mm  loss@{xc['freqs_ghz'][-1]}GHz {list(row['loss_db_mm'].values())[-1]:.4f} dB/mm"
            + f"  ({L:.0f} mm on board)")
    st = status(cache_dir, "xsec2d")
    return {"rule": "impedance", "spec": "G6", "title": "Stackup and impedance planner", "validation": st,
            "params": xc, "rows": rows,
            "metrics": {f"{r['layer']} {r['kind']} w={r['w']}" + (f" g={r['gap']}" if r['gap'] else ""):
                        {k: r.get(k) for k in ("Z0", "Zdiff", "delay_ps_mm")} for r in rows},
            "findings": [], "summary": {"geometries": len(rows), "validation": st}}
