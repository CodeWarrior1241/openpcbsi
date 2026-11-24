"""Run a set of checks on one board and write <out>/report.json (G16 adds HTML, images, markers)."""
import datetime
import json
import os
import time

from . import config as config_mod
from . import context, drc, model, snapshot, validate

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# #36: a checkout keeps cache/ and runs/ beside the package; an installed copy uses the working
# directory. OPENPCBSI_CACHE / OPENPCBSI_RUNS override both.
_CHECKOUT = os.path.isdir(os.path.join(ROOT, "boards")) or os.path.isdir(os.path.join(ROOT, ".git")) or \
    os.path.exists(os.path.join(ROOT, ".git"))
_BASE = ROOT if _CHECKOUT else os.path.join(os.getcwd(), "openpcbsi")
DETECTED = "detected from the board"   # what a report records as its config when none was given
CACHE = os.environ.get("OPENPCBSI_CACHE") or os.path.join(_BASE, "cache")
RUNS = os.environ.get("OPENPCBSI_RUNS") or os.path.join(_BASE, "runs")

ALL_CHECKS = ["drc", "impedance", "dc", "pdn", "crosstalk", "thermal", "vias", "si", "timing", "emc"]


def _standards(b):
    from . import standards
    out = {}
    for n in b.net_names[1:]:
        s = standards.of_net(n)
        if s:
            out[n] = s[0]
    return out


def configure(snap, cfg=None):
    """Board, config and context from a snapshot (the library API works from here, #32). `cfg` is a
    config file path, a dict, a Config, or None (detected from the board, EV5). The KiCad project next to
    the snapshot's board, if any, gives netclasses and custom rules (#2)."""
    from . import project as project_mod
    proj = project_mod.load(snap["board"])
    proj.set_nets(snap["nets"])
    for name, v in snap["nets"].items():
        v["netclass"] = proj.netclass(name, v.get("netclass"))
    if isinstance(cfg, config_mod.Config):
        pass
    elif isinstance(cfg, dict):
        cfg = config_mod.Config(cfg)
    elif cfg:
        cfg = config_mod.Config.load(cfg)
    else:
        from . import detect
        data, notes = detect.detect(snap)
        cfg = config_mod.Config(data)
        cfg.detected = notes
    cfg.project = proj
    b = model.Board(snap, cfg, cache_dir=CACHE)
    return b, cfg, context.Context(b, cfg)


def load(board_path, cfg_path=None, refresh=False):
    """Board, config and context of a board file. Without a config file the config is detected from the
    board (EV5, detect.py) and `cfg.detected` holds what was assumed and what was left out."""
    snap = snapshot.load_or_extract(board_path, CACHE, refresh=refresh)
    return configure(snap, cfg_path)


def _progress(out_dir, report, t0, starting):
    """#30: what is done, what runs now, and the partial report so far (written before each check)."""
    done = list(report["checks"])
    st = {"done": done, "running": starting, "seconds": round(time.time() - t0, 1),
          "updated": datetime.datetime.now().isoformat(timespec="seconds")}
    try:
        with open(os.path.join(out_dir, "progress.json"), "w") as f:
            json.dump(st, f)
        if done:
            from .paths import portable
            with open(os.path.join(out_dir, "report.partial.json.tmp"), "w") as f:
                json.dump(portable(json.loads(json.dumps(report, default=str)), out_dir), f)
            os.replace(os.path.join(out_dir, "report.partial.json.tmp"), os.path.join(out_dir, "report.partial.json"))
    except OSError:
        pass


def run_checks(board_path, cfg_path=None, checks=None, out_dir=None, label=None, log=print, rules=None):
    b, cfg, ctx = load(board_path, cfg_path)
    return run_loaded(b, cfg, ctx, checks, out_dir, label, log, rules)


def run_loaded(b, cfg, ctx, checks=None, out_dir=None, label=None, log=print, rules=None):
    t0 = time.time()
    checks = ALL_CHECKS if checks is None else checks
    label = label or f"{b.name}_{b.sha[:8]}"
    out_dir = out_dir or os.path.join(RUNS, label)
    os.makedirs(out_dir, exist_ok=True)
    log(f"openpcbsi: {b.path}")
    log(f"  sha256 {b.sha[:16]}  {len(b.layers)} layers  {len(b.net_names) - 1} nets  "
        f"{len(b.tracks['x1'])} tracks  {len(b.vias['x'])} vias  {len(b.zones)} zone fills")
    _ensure_validation(log, ctx)
    report = {
        "tool": "openpcbsi", "version": __import__("openpcbsi").__version__,
        "date": datetime.datetime.now().isoformat(timespec="seconds"),
        "board": b.path, "board_sha256": b.sha, "label": label, "config": cfg.path or DETECTED,
        "detected": getattr(cfg, "detected", None), "downgraded": b.snap.get("downgraded"),
        "kicad_file_version": b.snap.get("kicad_file_version") or snapshot.file_version(b.path),
        "zones_refilled": b.snap.get("zones_refilled"),
        "netclasses": {n: v["netclass"] for n, v in b.snap["nets"].items() if v.get("netclass") not in (None, "", "Default")},
        "project": cfg.project.describe() if getattr(cfg, "project", None) else None,
        "config_used": {k: cfg.data.get(k) for k in ("board", "classes", "areas", "rails", "engines")},
        "standards": _standards(b),
        "config_hash": cfg.hash(), "kicad": b.snap.get("kicad"),
        "stackup": b.stack.describe(), "layers": b.layers,
        "areas": {k: v for k, v in cfg["areas"].items()},
        "validation": {s: validate.status(CACHE, s) for s in ("xsec2d", "dc", "pdn", "thermal", "si")},
        "checks": {},
    }
    _progress(out_dir, report, t0, 'drc')
    if "drc" in checks:
        log("G1 electrical rules:")
        report["checks"]["drc"] = drc.run(ctx, only=rules, log=log)
    _progress(out_dir, report, t0, 'impedance')
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
    _progress(out_dir, report, t0, 'dc')
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
    _progress(out_dir, report, t0, 'pdn')
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
    _progress(out_dir, report, t0, 'crosstalk')
    if "crosstalk" in checks:
        from . import crosstalk
        log("G4 crosstalk scan:")
        report["checks"]["crosstalk"] = crosstalk.run(ctx, CACHE, log=log)
    _progress(out_dir, report, t0, 'thermal')
    if "thermal" in checks and cfg["thermal"]["parts"]:
        from . import thermal
        log("G13 board thermal:")
        report["checks"]["thermal"] = thermal.run(ctx, CACHE, out_dir, log=log, dc=report["checks"].get("dc"))
    _progress(out_dir, report, t0, 'vias')
    if "vias" in checks:
        from . import vias
        log("G9 via models:")
        report["checks"]["vias"] = vias.run(ctx, CACHE, log=log)
    _progress(out_dir, report, t0, 'si')
    if "si" in checks:
        from . import si
        log("G2/G3 net extraction and linear SI:")
        report["checks"]["si"] = si.run(ctx, CACHE, out_dir, log=log, vias=report["checks"].get("vias"))
    _progress(out_dir, report, t0, 'timing')
    if "timing" in checks:
        from . import timing
        log("Timing (delay and skew):")
        report["checks"]["timing"] = timing.run(ctx, CACHE, log=log)
    _progress(out_dir, report, t0, 'emc')
    if "emc" in checks:
        from . import emc
        log("G14 EMC estimate:")
        report["checks"]["emc"] = emc.run(ctx, report, log=log)
    _progress(out_dir, report, t0, None)
    report["seconds"] = round(time.time() - t0, 1)
    path = os.path.join(out_dir, "report.json")
    from .paths import portable
    with open(path + ".tmp", "w") as f:
        json.dump(portable(json.loads(json.dumps(report, default=_jd)), out_dir), f, indent=1)
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
