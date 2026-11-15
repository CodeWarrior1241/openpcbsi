"""openpcbsi command line.

  python3 -m openpcbsi evaluate BOARD [-c CONFIG] [--waivers FILE] [--families g1,rp,sd,si,th,emc,ev]
                                     [--checks ...] [--json OUT] [-o OUT] [--no-drc] [--no-markers]
                                     (G18/G21: scorecard and exit code 0 ok / 1 warn / 2 fail / 3 error)
  python3 -m openpcbsi init    BOARD [-o boards/NAME.toml] [--force]   (EV5: draft a config from the board)
  python3 -m openpcbsi check   BOARD [-c CONFIG] [--checks drc,impedance,...] [--rules ref_plane,...] [-o OUT]
  python3 -m openpcbsi compare BOARD_A BOARD_B -c CONFIG [--checks ...]
  python3 -m openpcbsi rediff  RUN_DIR_A RUN_DIR_B            (A/B from two finished runs)
  python3 -m openpcbsi validate
  python3 -m openpcbsi whatif  BOARD -c CONFIG SPEC.toml
  python3 -m openpcbsi region  BOARD -c CONFIG REGION.toml       (G8, openEMS)
  python3 -m openpcbsi snapshot BOARD

Checks: drc (G1), impedance (G6), dc (G10), pdn (G11/G12), crosstalk (G4), thermal (G13),
vias (G9), si (G2/G3/G7), emc (G14). Outputs go to runs/<board>_<sha8>/ under the openpcbsi folder.
Without -c, the config is detected from the board (detect.py) and the report says what was assumed.
"""
import argparse
import os
import sys

# keep solver threads modest (other jobs share this machine) and caches off the root disk
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("OMP_NUM_THREADS", "4")
_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.environ.setdefault("MPLCONFIGDIR", os.path.join(_root, "cache", "mpl"))


def main(argv=None):
    ap = argparse.ArgumentParser(prog="openpcbsi", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    e = sub.add_parser("evaluate", help="G18 design evaluation: every built check, KiCad DRC, scorecard, exit code")
    e.add_argument("board")
    e.add_argument("-c", "--config")
    e.add_argument("--waivers")
    e.add_argument("--families", default=None)
    e.add_argument("--checks", default=None)
    e.add_argument("--json", dest="json_out")
    e.add_argument("-o", "--out")
    e.add_argument("--no-drc", action="store_true")
    e.add_argument("--no-markers", action="store_true")
    i = sub.add_parser("init", help="draft a board config from the board (EV5 detection)")
    i.add_argument("board")
    i.add_argument("-o", "--out")
    i.add_argument("--force", action="store_true")
    c = sub.add_parser("check")
    c.add_argument("board")
    c.add_argument("-c", "--config")
    c.add_argument("--checks", default=None)
    c.add_argument("--rules", default=None)
    c.add_argument("-o", "--out")
    c.add_argument("--no-markers", action="store_true")
    k = sub.add_parser("compare")
    k.add_argument("board_a")
    k.add_argument("board_b")
    k.add_argument("-c", "--config", required=True)
    k.add_argument("--config-b")
    k.add_argument("--checks", default=None)
    rd = sub.add_parser("rediff", help="recompute an A/B comparison from two finished run directories")
    rd.add_argument("run_a")
    rd.add_argument("run_b")
    sub.add_parser("validate")
    w = sub.add_parser("whatif")
    w.add_argument("board")
    w.add_argument("-c", "--config")
    w.add_argument("spec")
    r = sub.add_parser("region")
    r.add_argument("board")
    r.add_argument("-c", "--config")
    r.add_argument("spec")
    s = sub.add_parser("snapshot")
    s.add_argument("board")
    a = ap.parse_args(argv)
    from . import run as runmod
    if a.cmd == "evaluate":
        from . import evaluate
        try:
            code, _ = evaluate.run(a.board, a.config, a.waivers, a.families,
                                   a.checks.split(",") if a.checks else None, a.out, a.json_out,
                                   markers=not a.no_markers, drc=not a.no_drc)
        except Exception as ex:  # G21: 3 = could not evaluate
            import traceback
            traceback.print_exc()
            print(f"openpcbsi evaluate: could not evaluate {a.board}: {ex}", file=sys.stderr)
            return 3
        return code
    if a.cmd == "init":
        from . import detect, snapshot
        snap = snapshot.load_or_extract(a.board, runmod.CACHE)
        cfg, notes = detect.detect(snap)
        out = a.out or os.path.join(runmod.ROOT, "boards", os.path.splitext(os.path.basename(a.board))[0] + ".toml")
        if os.path.exists(out) and not a.force:
            print(f"init: {out} exists; use --force or -o", file=sys.stderr)
            return 1
        with open(out, "w") as fh:
            fh.write(detect.to_toml(cfg, notes, os.path.abspath(a.board)))
        for x in notes["assumptions"]:
            print("  assumed:  " + x)
        for x in notes["left_out"]:
            print("  left out: " + x)
        print(f"init: {out}")
        return 0
    if a.cmd == "check":
        from . import report
        rep, t, od = runmod.run_checks(a.board, a.config, a.checks.split(",") if a.checks else None, a.out,
                                       rules=a.rules.split(",") if a.rules else None)
        report.write_all(rep, t, od, markers=not a.no_markers)
    elif a.cmd == "compare":
        from . import compare
        compare.run_compare(a.board_a, a.board_b, a.config, a.config_b, a.checks.split(",") if a.checks else None)
    elif a.cmd == "rediff":
        from . import compare
        out, od = compare.rediff(a.run_a, a.run_b)
        print(f"compare: {od}/compare.html")
    elif a.cmd == "validate":
        from . import validate, pdn, thermal
        validate.run(runmod.CACHE, ("xsec2d",))
        for solver in ("dc", "pdn"):
            cs = [x for x in pdn.validation_cases() if x["solver"] == solver]
            for x in cs:
                print(f"  [{'PASS' if x['pass'] else 'FAIL'}] {x['case']:55s} {x['value']:.6g} vs {x['reference']:.6g} ({x['error_pct']:+.3f} %)")
            validate.register(runmod.CACHE, solver, cs)
    elif a.cmd == "whatif":
        import tomllib
        from . import si
        b, cfg, ctx = runmod.load(a.board, a.config)
        spec = tomllib.load(open(a.spec, "rb"))["whatif"]
        od = os.path.join(runmod.RUNS, f"whatif_{os.path.splitext(os.path.basename(a.spec))[0]}")
        os.makedirs(od, exist_ok=True)
        si.whatif(ctx, spec, od)
        print(f"whatif: {od}/whatif.json")
    elif a.cmd == "region":
        import tomllib
        from . import region
        b, cfg, ctx = runmod.load(a.board, a.config)
        spec = tomllib.load(open(a.spec, "rb"))["region"]
        region.run(ctx, spec, runmod.CACHE, runmod.RUNS)
    elif a.cmd == "snapshot":
        from . import snapshot
        s_ = snapshot.load_or_extract(a.board, runmod.CACHE, refresh=True)
        print(f"{len(s_['tracks'])} tracks, {len(s_['vias'])} vias, {len(s_['zones'])} zone fills, {len(s_['pads'])} pads")


if __name__ == "__main__":
    sys.exit(main())
