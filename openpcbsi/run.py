"""Run a set of checks on one board and write <out>/report.json (G16 adds HTML, images, markers)."""
import datetime
import json
import os
import time

from . import config as config_mod
from . import context, drc, model, snapshot, validate

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CACHE = os.path.join(ROOT, "cache")
RUNS = os.path.join(ROOT, "runs")

ALL_CHECKS = ["drc", "impedance", "dc", "pdn", "crosstalk", "thermal", "vias", "si", "emc"]


def load(board_path, cfg_path=None, refresh=False):
    """Board, config and context. Without a config file the config is detected from the board
    (EV5, detect.py) and `cfg.detected` holds what was assumed and what was left out."""
    snap = snapshot.load_or_extract(board_path, CACHE, refresh=refresh)
    if cfg_path:
        cfg = config_mod.Config.load(cfg_path)
    else:
        from . import detect
        data, notes = detect.detect(snap)
        cfg = config_mod.Config(data)
        cfg.detected = notes
    b = model.Board(snap, cfg, cache_dir=CACHE)
    return b, cfg, context.Context(b, cfg)


def run_checks(board_path, cfg_path=None, checks=None, out_dir=None, label=None, log=print, rules=None):
    t0 = time.time()
    b, cfg, ctx = load(board_path, cfg_path)
    checks = ALL_CHECKS if checks is None else checks
    label = label or f"{b.name}_{b.sha[:8]}"
    out_dir = out_dir or os.path.join(RUNS, label)
    os.makedirs(out_dir, exist_ok=True)
    log(f"openpcbsi: {b.path}")
    log(f"  sha256 {b.sha[:16]}  {len(b.layers)} layers  {len(b.net_names) - 1} nets  "
        f"{len(b.tracks['x1'])} tracks  {len(b.vias['x'])} vias  {len(b.zones)} zone fills")
    _ensure_validation(log, ctx)
    report = {
        "tool": "openpcbsi", "version": "0.1.0",
        "date": datetime.datetime.now().isoformat(timespec="seconds"),
        "board": b.path, "board_sha256": b.sha, "label": label, "config": cfg.path or "detected from the board",
        "detected": getattr(cfg, "detected", None), "downgraded": b.snap.get("downgraded"),
        "config_hash": cfg.hash(), "kicad": b.snap.get("kicad"),
        "stackup": b.stack.describe(), "layers": b.layers,
        "areas": {k: v for k, v in cfg["areas"].items()},
        "validation": {s: validate.status(CACHE, s) for s in ("xsec2d", "dc", "pdn", "thermal", "si")},
        "checks": {},
    }
    if "drc" in checks:
        log("G1 electrical rules:")
        report["checks"]["drc"] = drc.run(ctx, only=rules, log=log)
    if "impedance" in checks:
        from . import planner
        log("G6 impedance planner (2D field solver):")
        report["checks"]["impedance"] = planner.run(ctx, CACHE, log=log)
    eng = cfg.data.get("engines", {})
    use_kipi_dc = "dc" in checks and cfg["rails"] and eng.get("dc") == "kipi"
    use_kipi_pdn = "pdn" in checks and cfg["rails"] and eng.get("pdn") == "kipi"
    kres = None
    if use_kipi_dc or use_kipi_pdn:
        from . import kipi_bridge
        if not kipi_bridge.available():
            log("KiPI is not set up (deps/setup_kipi.sh); falling back to openpcbsi's own solvers")
            use_kipi_dc = use_kipi_pdn = False
        else:
            mode = "+".join(m for m, on in (("dc", True), ("ac", use_kipi_pdn)) if on)
            try:
                kres, ktree, krail_vrm, ksetup = kipi_bridge.run(ctx, mode, os.path.join(out_dir, "kipi"), CACHE, log=log)
            except Exception as e:
                log(f"KiPI run failed: {e}; falling back to openpcbsi's own solvers")
                use_kipi_dc = use_kipi_pdn = False
    cross = eng.get("cross_check", False)
    if "dc" in checks and cfg["rails"]:
        from . import pdn
        log("G10 DC drop" + (" (KiPI)" if use_kipi_dc else "") + ":")
        own = []
        if not use_kipi_dc or cross:
            for name, rail in cfg["rails"].items():
                r, _ = pdn.dc_drop(ctx, name, rail)
                r["engine"] = "openpcbsi"
                own.append(r)
        if use_kipi_dc:
            rows = kipi_bridge.dc_rows(ctx, kres, ktree)
            for r in rows:
                log(f"  {r['rail']:10s} {r.get('regulation', '-'):12s} regulator->ball {r.get('drop_mv', float('nan')):7.2f} mV = "
                    f"{r.get('drop_pct_of_window', float('nan')):6.1f} % of window ({'PASS' if r.get('passed') else 'FAIL'})"
                    f"  from shunt {r.get('drop_from_source_mv', float('nan')):7.2f} mV  ground {r.get('v_ground_mv', float('nan')):.2f} mV"
                    f"  max via {r.get('max_via_current_a', 0):.2f} A")
            xc = []
            if cross:
                ob = {r["rail"]: r for r in own}
                for r in rows:
                    o = ob.get(r["rail"], {})
                    a_, b2 = r.get("drop_from_source_mv"), o.get("drop_mv")
                    xc.append({"rail": r["rail"], "kipi_mv": a_, "openpcbsi_mv": b2,
                               "diff_pct": round(100 * (b2 - a_) / a_, 1) if (a_ and b2 is not None) else None,
                               "kipi_via_a": r.get("max_via_current_a"), "openpcbsi_via_a": o.get("max_via_current_a")})
                log("  cross-check, drop from the shunt (or the sense point) to the worst ball, supply side:")
                for x in xc:
                    log(f"    {x['rail']:10s} KiPI {x['kipi_mv'] if x['kipi_mv'] is not None else float('nan'):7.2f} mV"
                        f"   openpcbsi {x['openpcbsi_mv'] if x['openpcbsi_mv'] is not None else float('nan'):7.2f} mV"
                        f"   diff {x['diff_pct'] if x['diff_pct'] is not None else float('nan'):+6.1f} %")
            report["checks"]["dc"] = {"spec": "G10", "engine": "kipi", "rails": rows, "cross_check": xc,
                                      "openpcbsi_rails": own, "kipi_warnings": (kres.get("dc") or {}).get("warnings", []),
                                      "kipi_vrms": (kres.get("dc") or {}).get("vrms", []),
                                      "validation": "KiPI accuracy suite (upstream)"}
        else:
            for r in own:
                rail = cfg["rails"][r["rail"]]
                log(f"  {r['rail']:10s} {r.get('regulation', '-'):12s} drop {r.get('drop_mv', float('nan')):7.2f} mV = "
                    f"{r.get('drop_pct_of_window', float('nan')):6.1f} % of +/-{rail.get('tol_pct', 3)} % window at {rail['current_a']} A"
                    f"  worst {r.get('worst_ball', '-')}  max via {r.get('max_via_current_a', 0):.2f} A"
                    + (f"  PROBLEM: {r['problems']}" if r["problems"] else ""))
            report["checks"]["dc"] = {"spec": "G10", "engine": "openpcbsi", "rails": own, "validation": validate.status(CACHE, "dc")}
    if "pdn" in checks and cfg["rails"]:
        from . import pdn
        log("G11/G12 PDN impedance" + (" (KiPI)" if use_kipi_pdn else "") + ", resonances, decoupling optimiser:")
        own = []
        if not use_kipi_pdn or cross:
            for name, rail in cfg["rails"].items():
                if "ripple_pct" not in rail:
                    continue
                r = pdn.pdn_impedance(ctx, name, rail)
                r["engine"] = "openpcbsi"
                own.append(r)
        ob = {r["rail"]: r for r in own}
        if use_kipi_pdn:
            rows = kipi_bridge.pdn_rows(ctx, kres, krail_vrm)
            for r in rows:
                o = ob.get(r["rail"], {})
                if o.get("optimiser"):
                    r["optimiser"] = dict(o["optimiser"], engine="openpcbsi network")
                if o.get("resonances"):
                    r["resonances_openpcbsi"] = o["resonances"]
                log(f"  {r['rail']:10s} Zt {r.get('target_ohm', 0)*1000:7.2f} mOhm  KiPI max |Z| {r.get('max_Z_ohm', 0)*1000:8.2f} mOhm "
                    f"({r.get('max_ratio', 0):.2f} x) at {r.get('f_at_max_hz', 0)/1e6:6.1f} MHz, {r.get('caps')} caps, "
                    f"plane C {r.get('plane_capacitance_nf')} nF"
                    + (f";  openpcbsi {o['max_Z_ohm']*1000:.2f} mOhm ({o['max_ratio']:.2f} x)" if o.get("max_Z_ohm") else ""))
            report["checks"]["pdn"] = {"spec": "G11/G12", "engine": "kipi", "rails": rows, "openpcbsi_rails": own,
                                       "validation": "KiPI accuracy suite (upstream)"}
        else:
            for r in own:
                opt = r.get("optimiser", {})
                log(f"  {r['rail']:10s} Zt {r.get('target_ohm', 0)*1000:7.2f} mOhm  max |Z| {r.get('max_Z_ohm', 0)*1000:8.2f} mOhm "
                    f"({r.get('max_ratio', 0):.2f} x) at {r.get('f_at_max_hz', 0)/1e6:6.1f} MHz;  optimiser: "
                    f"{'meets' if opt.get('met') else 'does not meet'} target with +{opt.get('add')}"
                    + (f", could drop {opt.get('removable')}" if opt.get("removable") else ""))
            report["checks"]["pdn"] = {"spec": "G11/G12", "engine": "openpcbsi", "rails": own, "validation": validate.status(CACHE, "pdn")}
    if kres is not None:
        report["kipi"] = {"lock": kipi_bridge.lock(), "setup": ksetup, "seconds": kres.get("seconds"),
                          "errors": [e.splitlines()[0] for e in kres.get("errors", [])],
                          "tree_notes": ktree.notes}
    if "crosstalk" in checks:
        from . import crosstalk
        log("G4 crosstalk scan:")
        report["checks"]["crosstalk"] = crosstalk.run(ctx, CACHE, log=log)
    if "thermal" in checks and cfg["thermal"]["parts"]:
        from . import thermal
        log("G13 board thermal:")
        report["checks"]["thermal"] = thermal.run(ctx, CACHE, out_dir, log=log, dc=report["checks"].get("dc"))
    if "vias" in checks:
        from . import vias
        log("G9 via models:")
        report["checks"]["vias"] = vias.run(ctx, CACHE, log=log)
    if "si" in checks:
        from . import si
        log("G2/G3 net extraction and linear SI:")
        report["checks"]["si"] = si.run(ctx, CACHE, out_dir, log=log, vias=report["checks"].get("vias"))
    if "emc" in checks:
        from . import emc
        log("G14 EMC estimate:")
        report["checks"]["emc"] = emc.run(ctx, report, log=log)
    report["seconds"] = round(time.time() - t0, 1)
    path = os.path.join(out_dir, "report.json")
    with open(path + ".tmp", "w") as f:
        json.dump(report, f, indent=1, default=_jd)
    os.replace(path + ".tmp", path)
    log(f"report: {path}  ({report['seconds']} s)")
    return report, (b, cfg, ctx), out_dir


def _ensure_validation(log, ctx=None):
    """Run the closed-form cases once per solver version (G17)."""
    from . import pdn
    if ctx is not None and validate.status(CACHE, "thermal") != "validated":
        from . import thermal
        validate.register(CACHE, "thermal", thermal.validation_cases(ctx))
        log(f"  G17 validation thermal: {validate.status(CACHE, 'thermal')}")
    for solver in ("xsec2d", "dc", "pdn"):
        if validate.status(CACHE, solver) == "validated":
            continue
        if solver == "xsec2d":
            validate.run(CACHE, ("xsec2d",), log=lambda *a: None)
        else:
            cases = [c for c in pdn.validation_cases() if c["solver"] == solver]
            validate.register(CACHE, solver, cases)
        log(f"  G17 validation {solver}: {validate.status(CACHE, solver)}")


def _jd(o):
    try:
        import numpy as np
        if isinstance(o, (np.integer,)):
            return int(o)
        if isinstance(o, (np.floating,)):
            return float(o)
        if isinstance(o, np.ndarray):
            return o.tolist()
        if isinstance(o, (np.bool_,)):
            return bool(o)
    except Exception:
        pass
    if isinstance(o, set):
        return sorted(o)
    return str(o)
