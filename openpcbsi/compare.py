"""G15: A/B comparison of two layouts of the same design.

Both boards run with the same configuration (rules, areas, solver settings); a run is refused if
the solver sections differ. Nets are matched by name, then by their pad set (ref.pad) for nets that
were renamed. Each numeric metric gets a direction (lower or higher is better) and the report
lists what got better or worse per check and item."""
import html
import json
import math
import os

from . import report as rep
from .run import RUNS, _jd, run_checks

LOWER = ("overshoot", "ringback", "rl_nyq", "uncovered", "crossings", "worst_return", "stub_mm", "uncoupled", "violation", "skew", "drop_mv",
         "drop_pct", "max_ratio", "max_z", "next", "fext", "max_stitch_gap", "max_edge_gap", "dangling",
         "unreferenced", "via_mm", "load_mm", "area_mm2", "tj_est", "t_board", "j_a_mm2", "dt_c", "worst_")
HIGHER = ("il_nyq", "eye_height", "eye_width", "margin", "min_spacing", "min_f_res", "min_board_edge", "min_plane_edge", "vias_in_pad", "solid", "stitching_vias")


def direction(metric):
    """+1: higher is better, -1: lower is better, 0: no judgement. IL in dB is negative (S21),
    so a higher IL_nyq_db is better; RL in dB (S11) is better when lower."""
    m = metric.lower()
    if any(k in m for k in HIGHER):
        return +1
    if any(k in m for k in LOWER):
        return -1
    return 0


def pad_sets(board):
    out = {}
    for p in board.pads:
        if p["net"]:
            out.setdefault(p["net"], set()).add(f"{p['ref']}.{p['num']}")
    return out


def match_nets(ba, bb):
    na, nb = set(ba.net_names[1:]), set(bb.net_names[1:])
    common = na & nb
    only_a, only_b = sorted(na - nb), sorted(nb - na)
    pa, pb = pad_sets(ba), pad_sets(bb)
    renamed = {}
    for a in only_a:
        sa = pa.get(a, set())
        if not sa:
            continue
        best, score = None, 0.0
        for b_ in only_b:
            sb = pb.get(b_, set())
            if not sb:
                continue
            j = len(sa & sb) / len(sa | sb)
            if j > score:
                best, score = b_, j
        if best and score >= 0.6:
            renamed[a] = (best, round(score, 2))
    only_a = [a for a in only_a if a not in renamed]
    only_b = [b_ for b_ in only_b if b_ not in {v[0] for v in renamed.values()}]
    return common, renamed, only_a, only_b


def flat(metrics):
    """{item: {metric: number}} with nested dicts flattened one level and lists dropped."""
    out = {}
    for item, m in metrics.items():
        row = {}
        for k, v in m.items():
            if isinstance(v, bool):
                continue
            if isinstance(v, (int, float)) and v is not None and math.isfinite(v):
                row[k] = v
        out[item] = row
    return out


def diff_metrics(ma, mb, rename):
    fa, fb = flat(ma), flat(mb)
    rows = []
    for item, ra in fa.items():
        key_b = item
        for old, (new, _) in rename.items():
            if old in item:
                key_b = item.replace(old, new)
        rb = fb.get(key_b)
        if rb is None:
            continue
        for k, va in ra.items():
            vb = rb.get(k)
            if vb is None or va == vb:
                continue
            d = direction(k)
            delta = vb - va
            verdict = "changed" if d == 0 else ("better" if delta * d > 0 else "worse")
            rows.append({"item": item, "metric": k, "A": va, "B": vb, "delta": round(delta, 6), "verdict": verdict})
    return rows


def checks_of(report):
    out = {}
    for name, chk in report["checks"].items():
        if isinstance(chk, list):
            for r in chk:
                out[r["rule"]] = r
        elif name in ("dc", "pdn"):
            out[name] = {"metrics": {r["rail"]: {k: v for k, v in r.items() if isinstance(v, (int, float))}
                                     for r in chk["rails"]}, "findings": []}
        else:
            out[name] = chk
    return out


def run_compare(board_a, board_b, cfg_path, cfg_b=None, checks=None, label=None, log=print):
    if cfg_b and cfg_b != cfg_path:
        from .config import Config
        ca, cb = Config.load(cfg_path), Config.load(cfg_b)
        for sec in ("xsec", "pdn", "dc", "thermal", "areas", "rules", "classes"):
            if ca.hash(sec) != cb.hash(sec):
                raise SystemExit(f"G15: section [{sec}] differs between the two configs; A/B runs must use the same settings")
    log("=== board A")
    ra, ta, da = run_checks(board_a, cfg_path, checks, log=log)
    rep.write_all(ra, ta, da, markers=True, log=log)
    log("=== board B")
    rb, tb, db = run_checks(board_b, cfg_b or cfg_path, checks, log=log)
    rep.write_all(rb, tb, db, markers=True, log=log)
    common, renamed, only_a, only_b = match_nets(ta[0], tb[0])
    ca, cb = checks_of(ra), checks_of(rb)
    out = {"A": {"board": ra["board"], "sha": ra["board_sha256"], "run": da},
           "B": {"board": rb["board"], "sha": rb["board_sha256"], "run": db},
           "config": cfg_path, "config_hash": ra["config_hash"],
           "nets": {"common": len(common), "renamed": renamed, "only_A": only_a, "only_B": only_b},
           "checks": {}}
    for name in ca:
        if name not in cb:
            continue
        rows = diff_metrics(ca[name].get("metrics", {}), cb[name].get("metrics", {}), renamed)
        fa = [f for f in ca[name].get("findings", []) if f["severity"] in ("error", "warning")]
        fb = [f for f in cb[name].get("findings", []) if f["severity"] in ("error", "warning")]
        out["checks"][name] = {"findings_A": len(fa), "findings_B": len(fb),
                               "better": sum(r["verdict"] == "better" for r in rows),
                               "worse": sum(r["verdict"] == "worse" for r in rows),
                               "changed": sum(r["verdict"] == "changed" for r in rows), "rows": rows}
    label = label or f"compare_{os.path.basename(da)}_vs_{os.path.basename(db)}"
    od = os.path.join(RUNS, label)
    os.makedirs(od, exist_ok=True)
    charts = {}
    pa = {r["rail"]: r for r in ra["checks"].get("pdn", {}).get("rails", []) if r.get("Z_ohm")}
    pb = {r["rail"]: r for r in rb["checks"].get("pdn", {}).get("rails", []) if r.get("Z_ohm")}
    for k in pa:
        if k in pb:
            charts[k] = rep.z_chart(pa[k], od, pb[k], ("A", "B"))
    with open(os.path.join(od, "compare.json"), "w") as f:
        json.dump(out, f, indent=1, default=_jd)
    write_html(out, od, charts)
    log("=== A/B summary (B relative to A)")
    log(f"  nets: {len(common)} common, {len(renamed)} renamed, {len(only_a)} only in A, {len(only_b)} only in B")
    for name, c in out["checks"].items():
        log(f"  {name:15s} errors+warnings A {c['findings_A']:4d}  B {c['findings_B']:4d}   metrics better {c['better']:4d}  worse {c['worse']:4d}  changed {c['changed']:4d}")
    log(f"compare: {os.path.join(od, 'compare.html')}")
    return out, od


def write_html(out, od, charts):
    H = [f"<!doctype html><html><head><meta charset=utf-8><meta name=viewport content='width=device-width,initial-scale=1'>"
         f"<title>openpcbsi A/B</title><style>{rep.CSS}</style></head><body>"]
    H.append("<h1>openpcbsi A/B comparison</h1>")
    H.append(f"<p class=muted>A: {html.escape(out['A']['board'])} ({out['A']['sha'][:12]})<br>"
             f"B: {html.escape(out['B']['board'])} ({out['B']['sha'][:12]})<br>config {html.escape(str(out['config']))} "
             f"({out['config_hash']}), identical for both runs</p>")
    n = out["nets"]
    H.append(f"<h2>Nets</h2><p>{n['common']} matched by name, {len(n['renamed'])} matched by pads (renamed).</p>")
    if n["renamed"]:
        H.append(rep._t([{"A": a, "B": b_[0], "pad overlap": b_[1]} for a, b_ in n["renamed"].items()], ["A", "B", "pad overlap"]))
    H.append(f"<p>Only in A ({len(n['only_A'])}): {html.escape(', '.join(n['only_A'][:80]))}</p>")
    H.append(f"<p>Only in B ({len(n['only_B'])}): {html.escape(', '.join(n['only_B'][:80]))}</p>")
    H.append("<h2>Summary per check</h2>")
    H.append(rep._t([dict({"check": k}, **{kk: v for kk, v in c.items() if kk != "rows"}) for k, c in out["checks"].items()],
                    ["check", "findings_A", "findings_B", "better", "worse", "changed"],
                    num=("findings_A", "findings_B", "better", "worse", "changed")))
    if charts:
        H.append("<h2>PDN impedance, A against B</h2><div class=grid>" + "".join(
            f"<img src='{p}' alt='impedance of {k}, A and B' loading=lazy>" for k, p in charts.items()) + "</div>")
    for k, c in out["checks"].items():
        rows = sorted([r for r in c["rows"] if r["verdict"] != "changed"], key=lambda r: (r["verdict"] != "worse", -abs(r["delta"])))
        if not rows:
            continue
        H.append(f"<h3>{html.escape(k)}: {c['worse']} worse, {c['better']} better</h3>")
        H.append(rep._t([dict(r, verdict=r["verdict"]) for r in rows[:150]], ["verdict", "item", "metric", "A", "B", "delta"],
                        num=("A", "B", "delta")))
    H.append("</body></html>")
    with open(os.path.join(od, "compare.html"), "w") as f:
        f.write("".join(H))


def rediff(run_a, run_b, label=None, log=print):
    """Recompute the A/B comparison from two finished runs (report.json in each run directory)."""
    from .run import load
    ra = json.load(open(os.path.join(run_a, "report.json")))
    rb = json.load(open(os.path.join(run_b, "report.json")))
    ba, bb = load(ra["board"], ra["config"])[0], load(rb["board"], rb["config"])[0]
    if ra["config_hash"] != rb["config_hash"]:
        raise SystemExit("G15: the two runs used different configurations")
    common, renamed, only_a, only_b = match_nets(ba, bb)
    ca, cb = checks_of(ra), checks_of(rb)
    out = {"A": {"board": ra["board"], "sha": ra["board_sha256"], "run": run_a},
           "B": {"board": rb["board"], "sha": rb["board_sha256"], "run": run_b},
           "config": ra["config"], "config_hash": ra["config_hash"],
           "nets": {"common": len(common), "renamed": renamed, "only_A": only_a, "only_B": only_b}, "checks": {}}
    for name in ca:
        if name in cb:
            rows = diff_metrics(ca[name].get("metrics", {}), cb[name].get("metrics", {}), renamed)
            fa = [f for f in ca[name].get("findings", []) if f["severity"] in ("error", "warning")]
            fb = [f for f in cb[name].get("findings", []) if f["severity"] in ("error", "warning")]
            out["checks"][name] = {"findings_A": len(fa), "findings_B": len(fb),
                                   "better": sum(r["verdict"] == "better" for r in rows),
                                   "worse": sum(r["verdict"] == "worse" for r in rows),
                                   "changed": sum(r["verdict"] == "changed" for r in rows), "rows": rows}
    od = os.path.join(RUNS, label or f"compare_{os.path.basename(run_a.rstrip('/'))}_vs_{os.path.basename(run_b.rstrip('/'))}")
    os.makedirs(od, exist_ok=True)
    charts = {}
    pa = {r["rail"]: r for r in ra["checks"].get("pdn", {}).get("rails", []) if r.get("Z_ohm")}
    pb = {r["rail"]: r for r in rb["checks"].get("pdn", {}).get("rails", []) if r.get("Z_ohm")}
    for k in pa:
        if k in pb:
            charts[k] = rep.z_chart(pa[k], od, pb[k], ("A", "B"))
    with open(os.path.join(od, "compare.json"), "w") as f:
        json.dump(out, f, indent=1, default=_jd)
    write_html(out, od, charts)
    for name, c in out["checks"].items():
        log(f"  {name:15s} errors+warnings A {c['findings_A']:4d}  B {c['findings_B']:4d}   better {c['better']:4d}  worse {c['worse']:4d}  changed {c['changed']:4d}")
    return out, od
