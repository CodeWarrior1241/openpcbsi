"""G6: stackup and impedance planner. Impedance of every track geometry used on every layer
(single-ended widths and pair width/gap), compared with class targets and fabricator figures."""
import math
from collections import defaultdict

import numpy as np

from . import xsec
from .validate import status


def scaled_stack(stack, fh, fe):
    """A copy of the stackup with every dielectric `fh` times as thick and `fe` times the er (depths redone)."""
    import copy
    st = copy.deepcopy(stack)
    items = sorted([("cu", c) for c in st.copper] + [("d", d) for d in st.dielectrics], key=lambda kv: kv[1]["z_top"])
    z = 0.0
    for kind, it in items:
        t = it["t"] * (fh if kind == "d" else 1.0)
        it["t"], it["z_top"], it["z_bot"] = t, z, z + t
        if kind == "d":
            it["er"] = it["er"] * fe
        z += t
    st.thickness = z
    return st


def _corners(stack, lname, refs, kind, w, g, xc, freqs):
    c = xc.get("corners") or {}
    de, fh, fe = c.get("etch_um", 20.0) * 2e-3, c.get("thickness_pct", 10.0) / 100, c.get("er_pct", 5.0) / 100
    out = []
    for sw, sh, se in ((+1, -1, +1), (-1, +1, -1)):            # low-Z corner, high-Z corner
        st = scaled_stack(stack, 1 + sh * fh, 1 + se * fe)
        ww = max(0.01, w + sw * de)
        if kind == "se":
            sec = xsec.build_section(st, lname, refs, [ww], [], xc, xc["mask_er"])
            out.append(xsec.solve_lines(sec, xc["res_um"], freqs[:1], xc["roughness_um"])["freqs"][0]["Z0"])
        else:
            gg = max(0.01, g - sw * de)
            sec = xsec.build_section(st, lname, refs, [ww, ww], [gg], xc, xc["mask_er"])
            out.append(xsec.solve_lines(sec, xc["res_um"], freqs[:1], xc["roughness_um"])["freqs"][0]["Zdiff"])
    return min(out), max(out)


def pair_geometries(ctx, where=None):
    """(layer, w, gap) -> length, from parallel P/N segment overlaps (exact vector geometry).
    `where`, if given, collects (layer, w, gap) -> {nets, longest overlap and its midpoint}."""
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
                    key = (int(T["layer"][i]), round(float(T["w"][i]), 4), round(gap, 3))
                    out[key] += ov
                    if where is not None:
                        w_ = where.setdefault(key, {"nets": set(), "best": 0.0})
                        w_["nets"] |= {pn, nn}
                        if ov > w_["best"]:
                            sm = (max(0.0, min(s1, s2)) + min(La, max(s1, s2))) / 2
                            w_.update(best=ov, x=float(T["x1"][i] + ux * sm), y=float(T["y1"][i] + uy * sm))
    return out


def run(ctx, cache_dir, log=print):
    b, cfg = ctx.b, ctx.cfg
    xc = cfg["xsec"]
    freqs = [f * 1e9 for f in xc["freqs_ghz"]]
    T = b.tracks
    rows = []
    # single-ended widths per layer (signal nets only; >= 1 mm in total)
    se = defaultdict(float)
    where = {}
    for k in range(len(T["x1"])):
        net = b.net_names[T["net"][k]]
        if net in cfg.reference_nets or b.kicad_class.get(net) == "PWR":
            continue
        key = (int(T["layer"][k]), round(float(T["w"][k]), 4))
        se[key] += float(T["len"][k])
        if net in ctx.partner:            # pair tracks are judged as pairs, not as single lines
            continue
        w_ = where.setdefault(("se",) + key, {"nets": set(), "best": 0.0})
        w_["nets"].add(net)
        if float(T["len"][k]) > w_["best"]:
            w_.update(best=float(T["len"][k]), x=float((T["x1"][k] + T["x2"][k]) / 2), y=float((T["y1"][k] + T["y2"][k]) / 2))
    pwhere = {}
    pg = pair_geometries(ctx, pwhere)
    # merge pair geometries closer than 2 um in gap
    pairs = defaultdict(float)
    for (li, w, g), L in pg.items():
        k2 = (li, w, round(g / 0.002) * 0.002)
        pairs[k2] += L
        src, dst = pwhere.get((li, w, g)), where.setdefault(("pair",) + k2, {"nets": set(), "best": 0.0})
        if src:
            dst["nets"] |= src["nets"]
            if src["best"] > dst["best"]:
                dst.update(best=src["best"], x=src["x"], y=src["y"])
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
        # targets of the classes of the nets that use this geometry (#5)
        wh = where.get((kind, li, w) if kind == "se" else (kind, li, w, g), {})
        nets = sorted(wh.get("nets", ()))
        row["nets"] = nets[:12]
        row["net_example"] = nets[0] if nets else ""
        if "x" in wh:
            row["x"], row["y"] = round(wh["x"], 3), round(wh["y"], 3)
            row["area"] = ctx.area_name(ctx.area_at([wh["x"]], [wh["y"]])[0])
        targets = []
        for cname in sorted({cfg.class_of(n) for n in nets} - {None}):
            c = cfg["classes"][cname]
            zt = c.get("z_diff") if kind == "pair" else c.get("z_se")
            if zt:
                targets.append({"class": cname, "target": float(zt), "tol": float(c.get("z_tol", 10.0))})
        if targets:
            row["targets"] = targets
            try:
                lo_, hi_ = _corners(b.stack, lname, refs, kind, w, g, xc, freqs)
                k = "Z0" if kind == "se" else "Zdiff"
                row[k + "_min"], row[k + "_max"] = round(lo_, 2), round(hi_, 2)
            except Exception as e:
                log(f"  corners {lname} {kind} w={w}: {e}")
            t0 = targets[0]
            zv = row["Zdiff"] if kind == "pair" else row["Z0"]
            row["target"] = t0["target"]
            row["delta_target"] = round(zv - t0["target"], 2)
            row["within_target"] = all(abs(zv - t["target"]) <= t["tol"] for t in targets)
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
    from . import judge
    finds = judge.impedance(rows, cfg)
    judged = sum(1 for r in rows if r.get("targets"))
    return {"rule": "impedance", "spec": "G6", "title": "Stackup and impedance planner", "validation": st,
            "params": xc, "rows": rows,
            "metrics": {f"{r['layer']} {r['kind']} w={r['w']}" + (f" g={r['gap']}" if r['gap'] else ""):
                        {k: r.get(k) for k in ("Z0", "Zdiff", "delay_ps_mm")} for r in rows},
            "findings": finds, "judged": judged,
            "summary": {"geometries": len(rows), "judged": judged, "validation": st}}
