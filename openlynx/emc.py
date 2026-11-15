"""G14: EMC estimate, a ranked list of likely radiators built from the other checks:

  - plane-cavity resonances (G12) whose largest response sits near a board edge, weighted by the
    longest unstitched stretch of that edge (PC9),
  - pairs with skew or uncoupled length (PC6): common-mode conversion ~ skew delay / rise time,
  - signal tracks that cross a split or gap in their reference plane (PC2) or run unreferenced:
    an interrupted return path radiates at the knee frequency 0.35 / tr.

Rule-of-thumb ranking for review, not a field computation (G8 / G12 are the routes to numbers)."""
import math


def run(ctx, report, log=print):
    b, cfg = ctx.b, ctx.cfg
    rows = []
    chk = {r["rule"]: r for r in report["checks"].get("drc", [])}
    edge_gap = (chk.get("stitching", {}).get("metrics", {}).get("_board_edge", {}) or {}).get("max_edge_gap_mm")
    x0, y0, x1, y1 = b.bbox
    for r in report["checks"].get("pdn", {}).get("rails", []):
        for res in r.get("resonances_openlynx", []) or r.get("resonances", []):
            if res.get("x") is None:
                continue
            if res["f_hz"] < 50e6:      # regulator / bulk-capacitor anti-resonances do not radiate
                continue
            d_edge = min(res["x"] - x0, x1 - res["x"], res["y"] - y0, y1 - res["y"])
            score = res["Z_ohm"] / max(r["target_ohm"], 1e-6) * (2.0 if d_edge < 5 else 1.0) * (1.5 if (edge_gap or 0) > 10 else 1.0)
            rows.append({"source": f"{r['rail']} plane resonance", "f_mhz": round(res["f_hz"] / 1e6, 1),
                         "where": f"{res['layer']} ({res['x']}, {res['y']}), {d_edge:.1f} mm from the edge",
                         "score": round(score, 2)})
    for key, m in chk.get("pair_coupling", {}).get("metrics", {}).items():
        pos = key.split("/")[0]
        tr = cfg.class_param(cfg.class_of(pos), "rise_ps", 300.0)
        td = 6.5  # ps/mm
        cm = (m["skew_mm"] * td + 0.5 * m["uncoupled_mm"] * td * 0.2) / tr
        if cm > 0.005:
            rows.append({"source": f"{key} common mode (skew {m['skew_mm']} mm, uncoupled {m['uncoupled_mm']} mm)",
                         "f_mhz": round(0.35 / (tr * 1e-12) / 1e6, 0), "where": "pair", "score": round(100 * cm, 2)})
    for rule in ("split_crossing", "unreferenced"):
        for f in chk.get(rule, {}).get("findings", []):
            if f["severity"] == "info":
                continue
            tr = cfg.class_param(cfg.class_of(f["net"]), "rise_ps", 1000.0) if cfg.class_of(f["net"]) else 1000.0
            rows.append({"source": f"{f['net']}: {f['msg']}", "f_mhz": round(0.35 / (tr * 1e-12) / 1e6, 0),
                         "where": f"{f['layer']} ({f['x']}, {f['y']})", "score": round((f.get("value") or 0.5) * 300.0 / tr, 2)})
    rows.sort(key=lambda r: -r["score"])
    for r in rows[:8]:
        log(f"  {r['score']:7.2f}  {r['f_mhz']:8.1f} MHz  {r['source'][:80]}  @ {r['where']}")
    if edge_gap:
        log(f"  board edge: longest stretch without a stitching via {edge_gap} mm")
    return {"rule": "emc", "spec": "G14", "title": "EMC estimate (ranked radiators, rule of thumb)", "validation": "rule-based",
            "rows": rows[:60], "metrics": {f"{i:02d} {r['source'][:40]}": {"score": r["score"]} for i, r in enumerate(rows[:30])},
            "findings": [], "summary": {"radiators": len(rows)}, "note": __doc__.split("\n\n")[-1].strip()}
