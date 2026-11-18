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
  python3 -m openpcbsi dfm     BOARD [-c CONFIG] [--fab jlcpcb|FILE] [--gerbers DIR] [--observed FILE]
                                     [--checks audit,output,copper,drill,mask,paste,parts,reflow,chips]
                                     [--waivers FILE] [--json OUT] [-o OUT] [--no-markers]
                                     (5.17: manufacturing review against a fabricator profile. Rules give the
                                     verdict and the exit code; rankings order places and judge nothing)

Checks: drc (G1), impedance (G6), dc (G10), pdn (G11/G12), crosstalk (G4), thermal (G13),
vias (G9), si (G2/G3/G7), emc (G14). Outputs go to runs/<board>_<sha8>/ under the openpcbsi folder.
Without -c, the config is detected from the board (detect.py) and the report says what was assumed.
"""
import argparse
import os
import sys

# keep solver threads modest (other jobs share this machine) and caches off the root disk
_T = "1" if "dfm" in sys.argv[1:2] else "4"     # dfm: many small sparse solves, in worker processes
os.environ.setdefault("OPENBLAS_NUM_THREADS", _T)
os.environ.setdefault("OMP_NUM_THREADS", _T)
from .run import CACHE as _CACHE  # noqa: E402
os.environ.setdefault("MPLCONFIGDIR", os.path.join(_CACHE, "mpl"))


def main(argv=None):
    ap = argparse.ArgumentParser(prog="openpcbsi", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    e = sub.add_parser("evaluate", help="G18 design evaluation: every built check, KiCad DRC, scorecard, exit code")
    e.add_argument("board", nargs="?", help="board file (omit with --live)")
    e.add_argument("--live", action="store_true", help="the board open in KiCad, through its IPC API (#32)")
    e.add_argument("--socket", help="KiCad API socket (default: kipy's)")
    e.add_argument("-c", "--config")
    e.add_argument("--waivers")
    e.add_argument("--families", default=None)
    e.add_argument("--checks", default=None)
    e.add_argument("--json", dest="json_out")
    e.add_argument("-o", "--out")
    e.add_argument("--no-drc", action="store_true")
    e.add_argument("--no-markers", action="store_true")
    e.add_argument("--stage", choices=["placement", "routing", "final"], help="layout stage (default: detected)")
    e.add_argument("--since", help="evaluation.json of an earlier run: report new, resolved and changed findings")
    e.add_argument("--em", choices=["none", "plan", "run"], default="none",
                   help="automatic 3D regions at discontinuities: list them (plan) or run them (minutes each)")
    so = sub.add_parser("signoff", help="evaluate and write a sign-off record tied to the exact inputs")
    for x in e._actions[1:]:
        if x.dest not in ("help",):
            so._add_action(x)
    so.add_argument("--reviewer")
    so.add_argument("--sign-key", help="SSH private key: sign the record (ssh-keygen -Y sign), authenticating it and the reviewer")
    vs = sub.add_parser("verify-signoff", help="check that a sign-off record still matches its inputs")
    vs.add_argument("record")
    vs.add_argument("--allowed-signers", help="ssh allowed_signers file: check the record's signature against the reviewer")
    vs.add_argument("--rerun", action="store_true", help="also evaluate the board again and compare the verdict")
    fc = sub.add_parser("fab-check", help="re-read a fabricator profile's pages and check that every quoted figure is still there")
    fc.add_argument("profile", nargs="?", default="jlcpcb")
    fc.add_argument("--json", dest="json_out")
    i = sub.add_parser("init", help="draft a board config from the board (EV5 detection)")
    i.add_argument("board")
    i.add_argument("-o", "--out")
    i.add_argument("--force", action="store_true")
    c = sub.add_parser("check")
    c.add_argument("board", nargs="?")
    c.add_argument("--live", action="store_true")
    c.add_argument("--socket")
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
    va = sub.add_parser("validate")
    va.add_argument("--3d", dest="three_d", action="store_true", help="also the openEMS cases (minutes each)")
    w = sub.add_parser("whatif")
    w.add_argument("board")
    w.add_argument("-c", "--config")
    w.add_argument("spec")
    r = sub.add_parser("region")
    r.add_argument("board")
    r.add_argument("-c", "--config")
    r.add_argument("spec")
    pl = sub.add_parser("plan", help="width (and gap) for a target impedance on a layer of this board's stackup")
    pl.add_argument("board")
    pl.add_argument("-c", "--config")
    pl.add_argument("--layer", required=True)
    pl.add_argument("--zdiff", type=float)
    pl.add_argument("--z0", type=float)
    pl.add_argument("--gap", type=float, default=None)
    co = sub.add_parser("correlate", help="compare with characterised boards (published measurements, deps/references.lock)")
    co.add_argument("--fetch", action="store_true", help="download the reference files at the pinned commits")
    co.add_argument("--reference", default="antmicro_si")
    s = sub.add_parser("snapshot")
    s.add_argument("board")
    dm = sub.add_parser("dfm", help="manufacturing review: rules against a fabricator profile, and rankings")
    dm.add_argument("board")
    dm.add_argument("-c", "--config")
    dm.add_argument("--fab", help="fabricator profile: a name (jlcpcb) or a TOML file; default from [dfm] fab")
    dm.add_argument("--checks", default=None, help="audit,output,copper,drill,mask,paste,parts,reflow,chips (default: all)")
    dm.add_argument("--gerbers", help="folder with the Gerber and drill files that will be sent (default: plot them now)")
    dm.add_argument("--observed", help="TOML list of observed defects ([[defect]] ref, kind): score the rankings against it")
    dm.add_argument("--waivers")
    dm.add_argument("--json", dest="json_out")
    dm.add_argument("-o", "--out")
    dm.add_argument("--no-markers", action="store_true")
    dm.add_argument("--no-cross-check", action="store_true", help="do not read the board a second time with KiCadFiles")
    a = ap.parse_args(argv)
    from . import run as runmod
    if getattr(a, "live", False):                 # #32: the board open in KiCad
        from . import live
        from .snapshot import InputError
        try:
            a.board = live.fetch_board(runmod.CACHE, a.socket)
        except InputError as ex:
            print(f"openpcbsi: {ex}", file=sys.stderr)
            return 3
        print(f"live board from KiCad -> {a.board}", file=sys.stderr)
    elif getattr(a, "cmd", None) in ("evaluate", "signoff", "check") and not a.board:
        print("openpcbsi: give a board file or --live", file=sys.stderr)
        return 2
    if a.cmd == "fab-check":
        from . import fab
        r = fab.check_pages(a.profile)
        print(f"{r['name']} profile, read {r['retrieved']} ({fab.age_days(r['retrieved'])} days ago)")
        for k, v in r["pages"].items():
            print(f"  page {k}: {'could not be read: ' + v['error'] if v['error'] else 'read'}  {v['url']}")
        n = len(r["quotes"])
        print(f"  {n - len(r['missing'])} of {n} quoted figures found on their pages")
        for m in r["missing"]:
            print(f"  NOT FOUND [{m['page']}] {m['where']}: \"{m['quote'][:110]}\"" + (f" ({m['note']})" if m.get("note") else ""))
        if a.json_out:
            import json
            json.dump(r, open(a.json_out, "w"), indent=1)
        print("every quote is still on its page: the profile can be dated today" if not r["missing"] else
              "the quotes not found above must be checked by hand and the profile updated (figures, quotes, retrieved date)")
        return 0 if not r["missing"] else 1
    if a.cmd == "verify-signoff":
        from . import signoff
        good, bad, rec = signoff.verify(a.record, a.allowed_signers, a.rerun)
        print(f"sign-off of {rec.get('board')} ({rec.get('date')}, verdict: {rec.get('verdict')})")
        for x in bad:
            print("  " + x)
        for x in rec.get("_notes", []):
            print("  note: " + x)
        print("record matches its inputs and its evaluation" if good else "record does NOT match its inputs or its evaluation")
        return 0 if good else 1
    if a.cmd in ("evaluate", "signoff"):
        from . import evaluate
        try:
            code, _ = evaluate.run(a.board, a.config, a.waivers, a.families,
                                   a.checks.split(",") if a.checks else None, a.out, a.json_out,
                                   markers=not a.no_markers, drc=not a.no_drc, stage=a.stage, since=a.since, em=a.em)
            if a.cmd == "signoff":
                from . import signoff
                od = os.path.dirname(_["report"])
                p, rec = signoff.write(_, od, a.config, a.waivers, a.reviewer, a.sign_key,
                                       {"families": a.families, "checks": a.checks, "stage": a.stage, "no_drc": a.no_drc})
                print(f"sign-off: {p}\nverdict:  {rec['verdict']}")
        except Exception as ex:  # G21: 3 = could not evaluate
            from .snapshot import InputError
            if not isinstance(ex, InputError):
                import traceback
                traceback.print_exc()
            print(f"openpcbsi evaluate: could not evaluate {a.board}: {ex}", file=sys.stderr)
            return 3
        return code
    if a.cmd == "dfm":
        from . import dfm
        from .snapshot import InputError
        try:
            code, _ = dfm.run(a.board, a.config, a.fab, a.checks.split(",") if a.checks else None, a.out, a.waivers,
                              a.json_out, markers=not a.no_markers, gerbers=a.gerbers, observed=a.observed,
                              cross=not a.no_cross_check)
        except Exception as ex:  # 3 = could not evaluate
            if not isinstance(ex, (InputError, FileNotFoundError, ValueError)):
                import traceback
                traceback.print_exc()
            print(f"openpcbsi dfm: could not evaluate {a.board}: {ex}", file=sys.stderr)
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
    if a.cmd == "correlate":
        from . import references
        out = references.correlate(runmod.CACHE, runmod.RUNS, a.reference, do_fetch=a.fetch)
        print(f"correlation: {os.path.join(runmod.RUNS, 'correlation', a.reference + '.json')}  "
              f"({'PASS' if out['pass'] else 'FAIL'}, tolerance {out['tolerance_pct']:g} % on resonance frequencies)")
        return 0 if out["pass"] else 1
    if a.cmd == "plan":
        from . import xsec
        b, cfg, ctx = runmod.load(a.board, a.config)
        li = b.L[a.layer]
        refs = [b.layers[r] for r in ctx.ref_layers(li)]
        xc = cfg["xsec"]
        target = a.zdiff or a.z0
        if not target:
            print("plan: give --zdiff or --z0", file=sys.stderr)
            return 2
        def z_of(w):
            if a.zdiff:
                sec = xsec.build_section(b.stack, a.layer, refs, [w, w], [a.gap or w], xc, xc["mask_er"])
                return xsec.solve_lines(sec, xc["res_um"], (1e9,))["freqs"][0]["Zdiff"]
            sec = xsec.build_section(b.stack, a.layer, refs, [w], [], xc, xc["mask_er"])
            return xsec.solve_lines(sec, xc["res_um"], (1e9,))["freqs"][0]["Z0"]
        lo, hi = 0.03, 2.0
        for _ in range(22):                    # Z falls with width: bisection
            mid = (lo + hi) / 2
            if z_of(mid) > target:
                lo = mid
            else:
                hi = mid
        w = round((lo + hi) / 2, 4)
        z = z_of(w)
        print(f"{a.layer} (references {', '.join(refs)}): w = {w:.4f} mm"
              + (f"  gap = {a.gap or w:.4f} mm  Zdiff = {z:.2f} ohm" if a.zdiff else f"  Z0 = {z:.2f} ohm"))
        return 0
    if a.cmd == "check":
        from . import report
        from .snapshot import InputError
        try:
            runmod.load(a.board, a.config)
        except InputError as ex:
            print(f"openpcbsi check: {ex}", file=sys.stderr)
            return 3
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
        from . import etch, reflow
        for solver, mod in (("etch", etch), ("reflow", reflow)):
            cs = mod.validation_cases()
            for x in cs:
                print(f"  [{'PASS' if x['pass'] else 'FAIL'}] {x['case']:55s} {x['value']:.6g} vs {x['reference']:.6g} ({x['error_pct']:+.3f} %)")
            validate.register(runmod.CACHE, solver, cs)
        if a.three_d:
            cs = validate.run3d(runmod.CACHE, runmod.RUNS)
            return 0 if all(c["pass"] for c in cs) else 1
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
        if "nets" in spec:                       # 3D regions with line ports (region3d.py)
            from . import region3d
            from .snapshot import InputError
            try:
                res = region3d.run(ctx, spec, runmod.CACHE, runmod.RUNS)
            except InputError as ex:             # e.g. OPENPCBSI_OPENEMS_PYTHON not set
                print(f"openpcbsi region: {ex}", file=sys.stderr)
                return 3
            path = os.path.join(res["dir"], "result.json")
            pf = res["preflight"]
            if not pf["ok"]:
                print("region pre-flight refused the model:\n  " + "\n  ".join(pf["problems"]), file=sys.stderr)
            print(f"region result: {path}")
            return 0 if pf["ok"] else 2
        region.run(ctx, spec, runmod.CACHE, runmod.RUNS)
    elif a.cmd == "snapshot":
        from . import snapshot
        s_ = snapshot.load_or_extract(a.board, runmod.CACHE, refresh=True)
        print(f"{len(s_['tracks'])} tracks, {len(s_['vias'])} vias, {len(s_['zones'])} zone fills, {len(s_['pads'])} pads")


if __name__ == "__main__":
    sys.exit(main())
