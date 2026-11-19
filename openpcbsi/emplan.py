"""Automatic 3D regions (#19, FS9): one openEMS region job per discontinuity on a fast class, each with
the standard fix as a variant, so the 3D result says whether the fix is worth it.

    series_cap    two-pad part between nets of a fast class    variant: reference voided under the pads (+0.15 mm)
    layer_change  via of a fast-class net                      variant: two return vias added within 0.7 mm
    split         split / gap crossing (from G1)               variant: the gap filled with reference copper
    connector     connector pads of a fast-class net           (no variant)
    breakout      the longest breakout neck (from #24)          (no variant)

Fast class: Nyquist at or above `[em] fast_ghz` (default 3 GHz). A job is a region spec for region3d
(nets + box) plus `variant`; `evaluate --em plan` lists them, `--em run` runs them (minutes each).

A job whose box sides cut its nets where they run at 45 degrees is planned in a frame turned by -45
degrees (orient), so a side cuts the route square. Each planned box is then fitted to region3d's own
pre-flight (fit_box): on a real board another net's
track, pad or via often ends near a box edge, which the pre-flight refuses. The fit moves each side that
has a problem, a step at a time, outwards up to FIT_GROW or inwards down to the features the job is about
(they stay FIT_CORE clear of the edge), keeping the same ports, and keeps the first box that passes. A job
whose box cannot be fitted says so (`fit.ok` false, `fit.why` the remaining problems).
"""
import math
import re

from .si import safe_name

HALF = 3.0          # mm around a point discontinuity
CLEAR = 2.5         # mm from the features to the box edge: the PML plus 8 cells (see region3d pre-flight)
MAX_JOBS = 12
FIT_STEP = 0.25     # mm: one move of a box side
FIT_GROW = 3.0      # mm: a side moves out by at most this
FIT_CORE = 0.5      # mm: the job's own features stay at least this far inside the box
FIT_POINT = 1.0     # mm: half-size of the core of a point discontinuity (split, connector, breakout)
FIT_TRIES = 48      # pre-flight builds per job
PLANNER_KEYS = ("kind", "parts", "core", "fit", "planned_box", "result", "variant_result")


def _fast_nets(ctx, fast_ghz):
    out = set()
    for n in ctx.b.net_names[1:]:
        c = ctx.cfg.class_of(n)
        if c and (ctx.cfg.class_param(c, "nyquist_ghz", 0) or 0) >= fast_ghz:
            out.add(n)
    return out


def _box(x, y, h=HALF):
    return [round(x - h, 3), round(y - h, 3), round(x + h, 3), round(y + h, 3)]


def _name(*parts):
    """A job name from part and net names: also a file and folder name (hierarchical nets are '/sheet/NET')."""
    return "_".join(safe_name(p) for p in parts)


def _bbox(pts, clear=CLEAR):
    """The features' bounding box grown by `clear` (pads: their size counted)."""
    x0 = min(p[0] - p[2] / 2 for p in pts) - clear
    y0 = min(p[1] - p[3] / 2 for p in pts) - clear
    x1 = max(p[0] + p[2] / 2 for p in pts) + clear
    y1 = max(p[1] + p[3] / 2 for p in pts) + clear
    return [round(v, 3) for v in (x0, y0, x1, y1)]


def plan(ctx, report, fast_ghz=3.0, fit=True, log=None):
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
        jobs.append({"kind": "series_cap", "name": _name("cap", ref), "nets": list(key), "parts": refs,
                     "box": _bbox(pts), "core": _bbox(pts, FIT_CORE),
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
        feat = vs + [(x + 0.7, y, 0.3, 0.3), (x - 0.7, y, 0.3, 0.3)]
        jobs.append({"kind": "layer_change", "name": _name("via", net), "nets": list(key), "box": _bbox(feat), "core": _bbox(feat, FIT_CORE),
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
            jobs.append({"kind": "split", "name": _name("split", f["net"]), "nets": pair(f["net"]), "box": _box(f["x"], f["y"]),
                         "core": _box(f["x"], f["y"], FIT_POINT),
                         "variant": {"name": "solid reference", "fill_layer": f.get("ref_layer"), "fill_box": _box(f["x"], f["y"], 1.0)}})
            break
    # connector launches
    for ref, f in b.fps.items():
        if not re.match(r"(J|P|CN|CON)\d", ref):
            continue
        hits = [p for p in b.pads_by_ref[ref] if p["net"] in fast]
        if hits:
            p = hits[0]
            jobs.append({"kind": "connector", "name": _name("launch", ref, p["net"]), "nets": pair(p["net"]), "box": _box(p["x"], p["y"]),
                         "core": _box(p["x"], p["y"], FIT_POINT)})
    # the longest breakout neck
    bo = [r for r in report["checks"].get("drc", []) if r["rule"] == "breakout"]
    if bo and bo[0]["findings"]:
        f = max(bo[0]["findings"], key=lambda q: q.get("neck_mm", 0))
        jobs.append({"kind": "breakout", "name": _name("breakout", f["net"]), "nets": pair(f["net"]), "box": _box(f["x"], f["y"]),
                     "core": _box(f["x"], f["y"], FIT_POINT)})
    used = set()
    for j in jobs:
        j.setdefault("fmax_ghz", 20)
        n, k = j["name"], 2
        while j["name"] in used:              # two names that differ only in characters safe_name drops
            j["name"], k = f"{n}_{k}", k + 1
        used.add(j["name"])
    jobs = jobs[:MAX_JOBS]
    inc = ctx.cfg.data.get("em", {}).get("include_mm")
    for j in jobs:
        if inc is not None:
            j["include_mm"] = inc
        orient(ctx, j)
    if fit:
        for j in jobs:
            fit_box(ctx, j, log=log)
    return jobs


def orient(ctx, job):
    """Turn the job's frame by -45 degrees when the box's sides cut its nets mostly where they run at 45
    or 135 degrees: an axis-aligned box cuts a diagonal pair at an angle, which the pre-flight refuses
    (see region3d.Rotated). The box and core are kept about their centre."""
    b = ctx.b
    x0, y0, x1, y1 = job["box"]
    diag = ortho = 0
    T = b.tracks
    # where the box's sides cut the job's nets: the ports' tracks (a 45-degree fan-out inside the box
    # does not matter; the route the sides cut does)
    for n in job["nets"]:
        if n not in b.net_id:
            continue
        for k in b.tracks_of(n):
            xa, ya, xb, yb = (float(T[c][k]) for c in ("x1", "y1", "x2", "y2"))
            L = math.hypot(xb - xa, yb - ya)
            if L <= 0:
                continue
            cuts = 0
            for c, (pa, pb, qa, qb, lo, hi) in ((x0, (xa, xb, ya, yb, y0, y1)), (x1, (xa, xb, ya, yb, y0, y1)),
                                               (y0, (ya, yb, xa, xb, x0, x1)), (y1, (ya, yb, xa, xb, x0, x1))):
                if (pa - c) * (pb - c) < 0:
                    q = qa + (c - pa) / (pb - pa) * (qb - qa)
                    cuts += lo < q < hi
            if not cuts:
                continue
            a = math.degrees(math.atan2(yb - ya, xb - xa)) % 90.0
            if abs(a - 45.0) < 10.0:
                diag += cuts
            elif a < 10.0 or a > 80.0:
                ortho += cuts
    if diag > ortho:
        bx = job["box"]
        job["rotate_deg"] = -45.0
        job["rotate_about"] = [round((bx[0] + bx[2]) / 2, 3), round((bx[1] + bx[3]) / 2, 3)]
    return job


def region_spec(job, variant=False):
    """The region3d spec of a job: the planner's own keys left out, the fix variant only when asked."""
    return {k: v for k, v in job.items() if k not in PLANNER_KEYS and (variant or k != "variant")}


def _check(ctx, job, box):
    """region3d's pre-flight of the job's base model in `box`: (pre-flight dict, its set of ports)."""
    from . import region3d
    spec = dict(region_spec(job), box=[round(v, 3) for v in box], preflight_only=True)
    reg = region3d.Region(ctx, spec)
    reg.build()
    pf = reg.preflight()
    return pf, sorted({p["net"] for p in reg.model["ports"]})


SIDE_AX = {"xmin": (0, -1), "ymin": (1, -1), "xmax": (2, +1), "ymax": (3, +1)}   # box index, outward sign


def fit_box(ctx, job, log=None, tries=FIT_TRIES):
    """Fit the job's box to region3d's pre-flight (see the module text); records `fit` in the job and
    updates `box`. The search is per side: a side whose pre-flight problems name it moves out or in, a
    step at a time, nearest first, until its problems are gone without adding any elsewhere; the nets with
    ports must stay those of the planned box (a side may move off a corner that cut a pair, so a port can
    change side), and the job's own features stay inside (`core`), so it still measures what it was
    planned for."""
    box0 = list(job["box"])
    core = job.get("core") or box0
    job["planned_box"] = box0
    try:
        pf, ports0 = _check(ctx, job, box0)
    except Exception as ex:                       # a model that cannot be built is reported, not raised
        job["fit"] = {"ok": False, "tries": 1, "why": f"model not built: {type(ex).__name__}: {ex}"}
        return job
    n, box = 1, list(box0)
    lim = {"xmin": (box0[0] - FIT_GROW, core[0]), "ymin": (box0[1] - FIT_GROW, core[1]),
           "xmax": (core[2], box0[2] + FIT_GROW), "ymax": (core[3], box0[3] + FIT_GROW)}
    for _ in range(3):                            # rounds: fixing one side can move another's problems
        if pf["ok"] or n >= tries:
            break
        bad = [s for s in SIDE_AX if pf.get("by_side", {}).get(s)]
        if not bad:
            break                                 # nothing a side move can fix (no reference plane, no port)
        progress = False
        for side in bad:
            if not pf.get("by_side", {}).get(side):
                continue                          # fixed by an earlier side's move
            i, sg = SIDE_AX[side]
            lo, hi = lim[side]
            steps = []
            for k in range(1, int(FIT_GROW / FIT_STEP) + 1):
                for d in (+k, -k):                # out, then in, nearest first
                    v = round(box[i] + sg * d * FIT_STEP, 3)
                    if lo - 1e-9 <= v <= hi + 1e-9:
                        steps.append(v)
            for v in steps:
                if n >= tries:
                    break
                cand = list(box)
                cand[i] = v
                pf2, ports2 = _check(ctx, job, cand)
                n += 1
                if ports2 != ports0:
                    continue
                if not pf2.get("by_side", {}).get(side) and pf2["n_problems"] <= pf["n_problems"] and \
                        (pf2["n_problems"] < pf["n_problems"] or len(pf2.get("by_side", {}).get("other", [])) < len(pf.get("by_side", {}).get("other", []))):
                    box, pf, progress = cand, pf2, True
                    break
            if pf["ok"]:
                break
        if not progress:
            break
    fit = {"ok": bool(pf["ok"]), "tries": n, "moved_mm": [round(b - a, 3) for a, b in zip(box0, box)]}
    if pf["ok"]:
        job["box"] = [round(v, 3) for v in box]
        if job.get("variant"):                   # the fix variant adds copper or vias: check it in the same box
            from . import region3d
            try:
                reg = region3d.Region(ctx, dict(region_spec(job, variant=True), preflight_only=True))
                reg.build()
                pv = reg.preflight()
                fit["variant_ok"] = pv["ok"]
                if not pv["ok"]:
                    fit["variant_why"] = "; ".join(pv["problems"][:3])
            except Exception as ex:
                fit["variant_ok"], fit["variant_why"] = False, f"model not built: {type(ex).__name__}: {ex}"
    else:
        fit["why"] = ("; ".join(pf["problems"][:3]) + (f" (and {pf['n_problems'] - 3} more)" if pf["n_problems"] > 3 else "")
                      + f"; searched {n} boxes, each side moved up to {FIT_GROW:g} mm out or in to the features")
    job["fit"] = fit
    if log:
        log(f"  3D region {job['name']}: " + (f"box {job['box']} passes pre-flight" + (f" (moved {fit['moved_mm']} mm)" if any(fit['moved_mm']) else "")
                                            if fit["ok"] else f"no box passes pre-flight: {fit['why']}"))
    return job


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
