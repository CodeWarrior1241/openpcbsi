"""openpcbsi -> KiPI driver. Runs inside KiPI's own environment, never in openpcbsi's:

    deps/kipi-venv/bin/python -I deps/kipi_driver.py JOB.json

JOB.json (written by openpcbsi/kipi_bridge.py):
    {"board": ".../board.kicad_pcb", "board_sha256": "...", "setup": ".../x.kipi-setup.json",
     "mode": "dc" | "ac" | "dc+ac", "ac_rails": ["U11", ...], "ground": "GND",
     "cache_dir": ".../cache/kipi", "out": ".../kipi_result.json"}

Writes plain JSON (numbers, lists, strings only) so openpcbsi never imports KiPI. KiPI is
pinned in kipi.lock; its headless entry points (solve_tree, solve_ac) are internal API, so a
version change must be re-checked against this file.
"""
import json
import math
import sys
import time
import traceback
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "KiPI"))

import numpy as np  # noqa: E402

from kipi.core.model import snapshot_io  # noqa: E402
from kipi.dc.pi.model import PISetup  # noqa: E402
from kipi.dc.pi.tree import solve_tree  # noqa: E402
from kipi.sources.kicad_pcb.source import KiCadFileSource  # noqa: E402


def _f(x):
    """JSON-safe float."""
    if x is None:
        return None
    try:
        x = float(x)
    except (TypeError, ValueError):
        return None
    return x if math.isfinite(x) else None


def board_snapshot(job, log):
    cache = Path(job["cache_dir"])
    cache.mkdir(parents=True, exist_ok=True)
    p = cache / f"kipi_snapshot_{job['board_sha256'][:16]}.kipi.json.gz"
    if p.exists():
        return snapshot_io.load(p), True
    t = time.time()
    snap = KiCadFileSource(job["board"]).load()
    snapshot_io.save(snap, p)
    log(f"KiPI read the board in {time.time() - t:.1f} s")
    return snap, False


def dc(snap, setup):
    r = solve_tree(snap, setup)
    out = {"solve_s": _f(r.solve_s), "warnings": list(r.warnings), "loads": [], "vrms": [], "groups": []}
    for L in r.loads:
        out["loads"].append({
            "name": L.name, "group": L.group, "current": _f(L.current), "nominal": _f(L.nominal),
            "v_supply": _f(L.v_supply), "v_supply_min": _f(L.v_supply_min), "v_ground": _f(L.v_ground),
            "v_effective": _f(L.v_effective), "drop": _f(L.drop), "v_low": _f(L.v_low), "v_high": _f(L.v_high),
            "v_worst": _f(L.v_worst), "tolerance": L.tolerance, "passed": L.passed, "status": L.status,
            "limit_low": _f(L.limit_low), "limit_high": _f(L.limit_high),
            "pin_groups": [[g[0]] + [_f(v) for v in g[1:]] for g in L.pin_groups]})
    for v in r.vrms:
        out["vrms"].append({"name": v.name, "kind": v.kind, "group": v.group, "v_set": _f(v.v_set),
                            "v_out": _f(v.v_out), "v_sense": _f(v.v_sense), "i_out": _f(v.i_out),
                            "v_in": _f(v.v_in), "i_in": _f(v.i_in), "regulation": v.regulation,
                            "passed": v.passed})
    for name, g in r.groups.items():
        f = g.result.fields
        vias = sorted(f.vias, key=lambda x: -abs(x.current))[:10]
        out["groups"].append({
            "name": name, "kind": g.kind, "nets": list(g.nets), "driver": g.driver, "converged": g.converged,
            "summary": {k: _f(v) for k, v in f.summary.items()},
            "top_vias": [{"x": _f(x.position[0]), "y": _f(x.position[1]), "current": _f(x.current),
                          "drill": _f(x.drill_mm), "kind": x.kind} for x in vias],
            "series": [{"name": s.name, "resistance": _f(s.resistance), "current": _f(s.current),
                        "v_a": _f(s.v_a), "v_b": _f(s.v_b), "power": _f(s.power)} for s in f.series],
            "mesh_nodes": int(len(g.result.solution.node_v))})
    return out, r


def ac(snap, setup, rail, ground, dc_result):
    from kipi.ac.solver.sweep import solve_ac
    r = solve_ac(snap, setup, rail, ground, dc_result=dc_result)
    return {"rail": r.rail, "rail_nets": list(r.rail_nets), "ground": r.ground, "solve_s": _f(r.solve_s),
            "freqs_hz": [float(x) for x in r.freqs], "caps": list(r.caps), "missing_caps": list(r.missing_caps),
            "links": list(r.links), "plane_capacitance_f": _f(r.plane_capacitance_f),
            "bandwidth_hz": _f(r.bandwidth_hz), "warnings": list(r.warnings), "bias_note": r.bias_note,
            "probes": [{"name": p.name, "z_mag": [float(x) for x in np.abs(p.z)],
                        "z_phase_deg": [float(x) for x in np.degrees(np.angle(p.z))],
                        "target_ohm": _f(p.target_ohm), "target_hz": _f(p.target_hz), "above_hz": _f(p.above_hz)}
                       for p in r.probes]}


def main(job_path):
    job = json.loads(Path(job_path).read_text())
    msgs = []
    log = lambda m: (msgs.append(m), print(m, flush=True))
    res = {"kipi_commit": (HERE / "kipi.lock").read_text() if (HERE / "kipi.lock").exists() else None,
           "log": msgs, "dc": None, "ac": {}, "errors": []}
    t0 = time.time()
    snap, cached = board_snapshot(job, log)
    setup = PISetup.load(job["setup"])
    dc_result = None
    if "dc" in job["mode"]:
        try:
            res["dc"], dc_result = dc(snap, setup)
            log(f"KiPI DC: {len(res['dc']['loads'])} loads, {len(res['dc']['groups'])} groups in {res['dc']['solve_s']:.1f} s")
        except Exception as e:
            res["errors"].append(f"dc: {e}\n{traceback.format_exc()}")
    if "ac" in job["mode"]:
        for rail in job.get("ac_rails", []):
            try:
                t = time.time()
                res["ac"][rail] = ac(snap, setup, rail, job.get("ground", "GND"), dc_result)
                log(f"KiPI AC {rail}: {time.time() - t:.1f} s")
            except Exception as e:
                res["errors"].append(f"ac {rail}: {e}\n{traceback.format_exc()}")
    res["seconds"] = round(time.time() - t0, 1)
    Path(job["out"]).write_text(json.dumps(res))


if __name__ == "__main__":
    main(sys.argv[1])
