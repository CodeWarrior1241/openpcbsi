"""Automatic 3D regions (#19, FS9): one openEMS region job per discontinuity on a fast class, each with
the standard fix as a variant, so the 3D result says whether the fix is worth it.

    series_cap    two-pad part between nets of a fast class    variant: reference voided under the pads (+0.15 mm)
    layer_change  via of a fast-class net                      variant: two return vias added within 0.7 mm
    split         split / gap crossing (from G1)               variant: the gap filled with reference copper
    connector     connector pads of a fast-class net           (no variant)
    breakout      the longest breakout neck (from #24)          (no variant)

Fast class: Nyquist at or above `[em] fast_ghz` (default 3 GHz). A job is a region spec for region3d
(nets + box) plus `variant`; `evaluate --em plan` lists them, `--em run` runs them (minutes each).
"""
import math
import re

HALF = 3.0          # mm around a point discontinuity
CLEAR = 2.5         # mm from the features to the box edge: the PML plus 8 cells (see region3d pre-flight)
MAX_JOBS = 12


def _fast_nets(ctx, fast_ghz):
    out = set()
    for n in ctx.b.net_names[1:]:
        c = ctx.cfg.class_of(n)
        if c and (ctx.cfg.class_param(c, "nyquist_ghz", 0) or 0) >= fast_ghz:
            out.add(n)
    return out


def _box(x, y, h=HALF):
    return [round(x - h, 3), round(y - h, 3), round(x + h, 3), round(y + h, 3)]


def _bbox(pts, clear=CLEAR):
    """The features' bounding box grown by `clear` (pads: their size counted)."""
    x0 = min(p[0] - p[2] / 2 for p in pts) - clear
    y0 = min(p[1] - p[3] / 2 for p in pts) - clear
    x1 = max(p[0] + p[2] / 2 for p in pts) + clear
    y1 = max(p[1] + p[3] / 2 for p in pts) + clear
    return [round(v, 3) for v in (x0, y0, x1, y1)]


def plan(ctx, report, fast_ghz=3.0):
    b = ctx.b
    fast = _fast_nets(ctx, fast_ghz)
    jobs = []
    if not fast:
        return jobs
    pair = lambda n: [x for x in (n, ctx.partner.get(n)) if x]
    # series capacitors (AC coupling) between fast nets
    seen = set()
    for ref, pads in b.pads_by_ref.items():
        pads = [p for p in pads if p["net"]]
        if len(pads) != 2 or not re.match(r"C\d", ref) or b.fps[ref].get("dnp"):
            continue
        a, c = pads[0]["net"], pads[1]["net"]
        if a not in fast and c not in fast:
            continue
        key = tuple(sorted(pair(a) + pair(c)))
        if key in seen:
            continue
        seen.add(key)
        others = [r for r, ps in b.pads_by_ref.items() if r != ref and len(ps) == 2 and {q["net"] for q in ps} <= set(key)]
        refs = [ref] + others[:1]
        pts = [(p["x"], p["y"], p["sx"], p["sy"]) for r in refs for p in b.pads_by_ref[r]]
        jobs.append({"kind": "series_cap", "name": f"cap_{ref}", "nets": list(key), "parts": refs,
                     "box": _bbox(pts),
                     "variant": {"name": "cut-out", "void_under": refs, "margin_mm": 0.15}})
    # layer changes on fast nets
    seen = set()
    for k in range(len(b.vias["x"])):
        net = b.net_names[int(b.vias["net"][k])]
        if net not in fast:
            continue
        key = tuple(sorted(pair(net)))
        if key in seen:
            continue
        seen.add(key)
        x, y = float(b.vias["x"][k]), float(b.vias["y"][k])
        vs = [(float(b.vias["x"][q]), float(b.vias["y"][q]), float(b.vias["d"][q]), float(b.vias["d"][q]))
              for q in range(len(b.vias["x"])) if b.net_names[int(b.vias["net"][q])] in key]
        jobs.append({"kind": "layer_change", "name": f"via_{net}", "nets": list(key), "box": _bbox(vs + [(x + 0.7, y, 0.3, 0.3), (x - 0.7, y, 0.3, 0.3)]),
                     "variant": {"name": "return vias", "add_vias": [[round(x + 0.7, 3), round(y, 3)], [round(x - 0.7, 3), round(y, 3)]]}})
    # split crossings
    for rule in report["checks"].get("drc", []):
        if rule["rule"] != "split_crossing":
            continue
        done = set()
        for f in rule["findings"]:
            if f["severity"] != "error" or f["net"] in done:
                continue
            done |= set(pair(f["net"]))
            jobs.append({"kind": "split", "name": f"split_{f['net']}", "nets": pair(f["net"]), "box": _box(f["x"], f["y"]),
                         "variant": {"name": "solid reference", "fill_layer": f.get("ref_layer"), "fill_box": _box(f["x"], f["y"], 1.0)}})
            break
    # connector launches
    for ref, f in b.fps.items():
        if not re.match(r"(J|P|CN|CON)\d", ref):
            continue
        hits = [p for p in b.pads_by_ref[ref] if p["net"] in fast]
        if hits:
            p = hits[0]
            jobs.append({"kind": "connector", "name": f"launch_{ref}_{p['net']}", "nets": pair(p["net"]), "box": _box(p["x"], p["y"])})
    # the longest breakout neck
    bo = [r for r in report["checks"].get("drc", []) if r["rule"] == "breakout"]
    if bo and bo[0]["findings"]:
        f = max(bo[0]["findings"], key=lambda q: q.get("neck_mm", 0))
        jobs.append({"kind": "breakout", "name": f"breakout_{f['net']}", "nets": pair(f["net"]), "box": _box(f["x"], f["y"])})
    for j in jobs:
        j.setdefault("fmax_ghz", 20)
    return jobs[:MAX_JOBS]


def judge(job, res, base=None, cls_tol_pct=10.0):
    """FS9 findings from a region result: TDR excursion, Sdd11 at Nyquist, Scd21; the better variant."""
    out = []
    if not res.get("ran"):
        return [("warn", f"{job['name']}: not run ({'; '.join(res.get('preflight', {}).get('problems', [])[:2]) or res.get('reason', '')})")]
    conv = res.get("convergence", {})
    if not conv.get("converged"):
        out.append(("warn", f"{job['name']}: 3D result not accepted: {conv.get('reason')}"))
    ex = res.get("tdr_excursion")
    if ex and max(abs(ex["min_pct"]), abs(ex["max_pct"])) > cls_tol_pct:
        out.append(("warn", f"{job['name']}: TDR {ex['min_pct']:+.0f} / {ex['max_pct']:+.0f} % of the line (limit {cls_tol_pct:g} %)"))
    return out
