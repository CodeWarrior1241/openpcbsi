"""G4: board-wide crosstalk scan. Every victim net in scope against every aggressor within
max_spacing_mm on the same layer: coupled length and spacing from the copper raster, backward and
forward coupling coefficients from the 2D solver (G6) for that layer, widths and spacing.

  NEXT = Kb * min(1, 2 Lc td / tr)          (saturates once the coupled length exceeds tr / 2td)
  FEXT = |Kf| Lc / tr                        (Kf in ps/mm; ~0 in homogeneous stripline)

tr is the aggressor's rise time (its class), default_rise_ps for nets without a class."""
import math
from collections import defaultdict

import numpy as np

from . import drc, xsec
from .context import runs
from .validate import status


def run(ctx, cache_dir, log=print):
    b, cfg = ctx.b, ctx.cfg
    p = cfg["crosstalk"]
    xc = cfg["xsec"]
    cls = p.get("classes") or list(cfg["classes"].keys())
    victims = [n for n in b.net_names[1:] if cfg.class_of(n) in cls]
    default_tr = p.get("default_rise_ps", 1000.0)
    table = {}

    def coeff(li, wv, wa, s):
        key = (li, round(wv, 3), round(wa, 3), round(max(s, 0.05) / 0.01) * 0.01)
        if key not in table:
            refs = [b.layers[r] for r in ctx.ref_layers(li)]
            sec = xsec.build_section(b.stack, b.layers[li], refs, [key[1], key[2]], [key[3]], xc, xc["mask_er"])
            r = xsec.solve_lines(sec, max(xc["res_um"], 5.0), (1e9,), 0.0)["freqs"][0]
            # delay of the victim alone (conductor 0): sqrt(L11 C11)
            table[key] = {"Kb": abs(r["Kb"]), "Kf_ps_mm": abs(r["Kf_ps_mm"]), "td_ps_mm": r["delay_diff_ps_mm"]}
        return table[key]

    rows = []
    T = b.tracks
    for v in victims:
        partner = ctx.partner.get(v)
        prof = drc.coupling_profile(ctx, v, p["max_spacing_mm"], step=0.05, exclude=[partner] if partner else [])
        for pr in prof:
            near = (pr["spacing"] <= p["max_spacing_mm"])
            for a in np.unique(pr["agg"][near]):
                an = b.net_names[a]
                m = near & (pr["agg"] == a)
                for a0, a1 in runs(m):
                    Lc = float(pr["dl"][a0:a1].sum())
                    if Lc < 0.5:
                        continue
                    s = float(np.median(pr["spacing"][a0:a1]))
                    wv = float(np.median(pr["w"][a0:a1]))
                    ka = b.tracks_of(an, pr["layer"])
                    wa = float(np.median(T["w"][ka])) if len(ka) else wv
                    c = coeff(pr["layer"], wv, wa, s)
                    tr = cfg.class_param(cfg.class_of(an), "rise_ps", default_tr) if cfg.class_of(an) else default_tr
                    next_ = c["Kb"] * min(1.0, 2 * Lc * c["td_ps_mm"] / tr)
                    fext = c["Kf_ps_mm"] * Lc / tr
                    i = (a0 + a1) // 2
                    rows.append({"victim": v, "aggressor": an, "layer": b.layers[pr["layer"]], "coupled_mm": round(Lc, 2),
                                 "spacing_mm": round(s, 3), "rise_ps": tr, "Kb": round(c["Kb"], 4),
                                 "NEXT_pct": round(100 * next_, 3), "FEXT_pct": round(100 * fext, 3),
                                 "x": round(float(pr["x"][i]), 3), "y": round(float(pr["y"][i]), 3),
                                 "area": ctx.area_name(pr["area"][i])})
    rows.sort(key=lambda r: -max(r["NEXT_pct"], r["FEXT_pct"]))
    top = rows[:p["top"]]
    for r in top[:10]:
        log(f"  {r['victim']:18s} <- {r['aggressor']:18s} {r['layer']:7s} {r['coupled_mm']:6.2f} mm at {r['spacing_mm']:.3f} mm"
            f"  NEXT {r['NEXT_pct']:5.2f} %  FEXT {r['FEXT_pct']:5.2f} %" + (f"  [{r['area']}]" if r["area"] else ""))
    per_victim = defaultdict(lambda: {"worst_NEXT_pct": 0.0, "worst_FEXT_pct": 0.0, "aggressors": 0})
    for r in rows:
        m = per_victim[r["victim"]]
        m["worst_NEXT_pct"] = max(m["worst_NEXT_pct"], r["NEXT_pct"])
        m["worst_FEXT_pct"] = max(m["worst_FEXT_pct"], r["FEXT_pct"])
        m["aggressors"] += 1
    findings = [drc.finding("crosstalk", "info" if r["area"] else "warning",
                            f"{r['aggressor']} -> NEXT {r['NEXT_pct']:.2f} %, FEXT {r['FEXT_pct']:.2f} % over {r['coupled_mm']} mm at {r['spacing_mm']} mm",
                            r["victim"], r["layer"], r["x"], r["y"], value=max(r["NEXT_pct"], r["FEXT_pct"]), unit="%",
                            area=r["area"], aggressor=r["aggressor"])
                for r in top if max(r["NEXT_pct"], r["FEXT_pct"]) >= p.get("report_pct", 1.0)]
    return {"rule": "crosstalk", "spec": "G4", "title": "Crosstalk scan (coupled-line estimates)",
            "validation": status(cache_dir, "xsec2d"), "params": p, "pairs": len(rows), "top": top,
            "to_region_solver": [{"victim": r["victim"], "aggressor": r["aggressor"], "x": r["x"], "y": r["y"],
                                  "layer": r["layer"]} for r in top[:5]],
            "metrics": dict(per_victim), "findings": findings,
            "summary": {"pairs": len(rows), "worst_pct": max([max(r["NEXT_pct"], r["FEXT_pct"]) for r in rows], default=0)}}
