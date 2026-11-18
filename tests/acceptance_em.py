"""openpcbsi acceptance tests for the 3D-region planner and runner (#9-#20) and the pad-void conflict
between pad_plane and ref_plane. Each test builds its own boards (tests/harness/boards_em.py), runs
openpcbsi as a user would (the CLI, in a subprocess) and compares with an oracle: the planted defect, the
pre-flight itself, or a count of solver runs. Each test function returns a short note and raises
AssertionError on failure; their docstrings' first line is the title.

    python3 tests/acceptance_em.py                # all tests here except the openEMS one
    python3 tests/acceptance_em.py --slow         # also the openEMS one (minutes; needs OPENPCBSI_OPENEMS_PYTHON)
    python3 tests/acceptance_em.py t_em_fit ...   # some tests, by name

Run with KiCad's python3 (the boards are built through pcbnew).
"""
import json
import os
import re
import shutil
import stat
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(HERE, "out", "em")
sys.path.insert(0, os.path.join(HERE, "harness"))
sys.path.insert(0, ROOT)

import boards_em  # noqa: E402  (needs pcbnew)

LANE = ["LANE0_RX_P", "LANE0_RX_N", "LANE0_RXC_P", "LANE0_RXC_N"]
CAP_BOX = [16.75, 6.25, 23.25, 13.75]      # the series-cap box the planner draws around C1 / C2 (boards_em)


# ---- helpers -------------------------------------------------------------------------------------
def need(cond, msg):
    if not cond:
        raise AssertionError(msg)


def out_dir(name):
    d = os.path.join(OUT, name)
    if os.path.isdir(d):
        shutil.rmtree(d)
    os.makedirs(d)
    return d


def cli(args, env=None, unset=(), timeout=1800):
    e = dict(os.environ, **(env or {}))
    for k in unset:
        e.pop(k, None)
    p = subprocess.run([sys.executable, "-m", "openpcbsi"] + args, cwd=ROOT, capture_output=True, text=True,
                       timeout=timeout, env=e)
    return p.returncode, p.stdout, p.stderr


def _load(p):
    return json.load(open(p)) if p and os.path.exists(p) else None


def evaluate(board, *args, env=None):
    out = os.path.join(os.path.dirname(board), "run")
    if os.path.isdir(out):
        shutil.rmtree(out)
    rc, so, se = cli(["evaluate", board, "-o", out, "--no-markers"] + list(args), env=env)
    return {"rc": rc, "out": so, "err": se, "ev": _load(os.path.join(out, "evaluation.json")), "dir": out}


def region(board, spec, env=None, unset=(), fname="region.toml"):
    """`openpcbsi region BOARD SPEC.toml` with `spec` (a dict) as the [region] table."""
    def val(v):
        return json.dumps(v)
    path = os.path.join(os.path.dirname(board), fname)
    with open(path, "w") as fh:
        fh.write("[region]\n" + "".join(f"{k} = {val(v)}\n" for k, v in spec.items() if not isinstance(v, dict)))
    rc, so, se = cli(["region", board, path], env=env, unset=unset)
    m = re.search(r"region result: (\S+)", so)
    return rc, so, se, (_load(m.group(1)) if m else None), (m.group(1) if m else None)


def check_rules(board, rules, cfg=None):
    out = os.path.join(os.path.dirname(board), "check_" + (os.path.splitext(os.path.basename(cfg))[0] if cfg else "detected"))
    if os.path.isdir(out):
        shutil.rmtree(out)
    rc, so, se = cli(["check", board, "--checks", "drc", "--rules", rules, "-o", out, "--no-markers"] + (["-c", cfg] if cfg else []))
    rep = _load(os.path.join(out, "report.json"))
    need(rep, f"no report (rc {rc}): {(se or so)[-300:]}")
    return {r["rule"]: r for r in rep["checks"]["drc"]}


def job(ev, kind):
    js = [j for j in (ev or {}).get("em_jobs") or [] if j["kind"] == kind]
    need(js, f"no {kind} job planned: {[j['kind'] for j in (ev or {}).get('em_jobs') or []]}")
    return js[0]


# ---- 1. names ------------------------------------------------------------------------------------
def t_em_names():
    """Automatic 3D regions on a multi-sheet design: hierarchical and '../' net names give safe job files and run folders"""
    from openpcbsi.si import inside, safe_name
    d = out_dir("names")
    notes = []
    for name, prefix in (("hier", "/serdes/"), ("trav", "/../../../../ESC3/")):
        b = boards_em.hier_lane(d, name, prefix=prefix)
        r = evaluate(b, "--checks", "drc", "--no-drc", "--em", "plan")
        need(r["ev"] and r["rc"] != 3, f"{name}: evaluate --em plan failed (rc {r['rc']}): {(r['err'] or r['out'])[-400:]}")
        jobs = r["ev"].get("em_jobs") or []
        need({"series_cap", "layer_change"} <= {j["kind"] for j in jobs}, f"{name}: planned {[j['kind'] for j in jobs]}")
        jd = os.path.join(r["dir"], "em_jobs")
        for j in jobs:
            need(re.fullmatch(r"[A-Za-z0-9_+\-][A-Za-z0-9_.+\-]*", j["name"]) and ".." not in j["name"].split("_"),
                 f"{name}: job name {j['name']!r} is not a safe file name")
            need(os.path.exists(os.path.join(jd, j["name"] + ".json")), f"{name}: no job file for {j['name']}")
        for f in os.listdir(jd):
            need(os.path.realpath(os.path.join(jd, f)).startswith(os.path.realpath(jd) + os.sep), f"{name}: {f} outside em_jobs")
        notes.append(f"{prefix} -> {sorted(j['name'] for j in jobs)}")
        # the region command with the raw names, and a region name that tries to climb out of runs/em
        runs = os.path.join(d, name, "runs")
        rc, so, se, res, path = region(b, {"name": "../../../ESC3/../x", "nets": [prefix + n for n in LANE], "box": CAP_BOX,
                                           "preflight_only": True}, env={"OPENPCBSI_RUNS": runs})
        need(res is not None, f"{name}: region failed (rc {rc}): {(se or so)[-400:]}")
        need(os.path.realpath(path).startswith(os.path.realpath(os.path.join(runs, "em")) + os.sep),
             f"{name}: region result {path} outside {runs}/em")
    # nothing was written where a '../' would have led
    for k in range(1, 7):
        up = os.path.join(d, *([".."] * k), "ESC3")
        need(not os.path.exists(up), f"{os.path.realpath(up)} exists: a name reached outside its folder")
    # the guard itself
    need(safe_name("..") == "net" and safe_name("/../../x") == "x", f"safe_name('..') = {safe_name('..')!r}")
    try:
        inside(d, os.path.join(d, "..", "escape.json"))
        raise AssertionError("inside() let a path leave its folder")
    except ValueError:
        pass
    return "; ".join(notes) + "; region run folders inside runs/em; nothing written outside"


# ---- 2. pre-flight -------------------------------------------------------------------------------
def t_em_preflight_crossings():
    """3D pre-flight: another net's track crossing the absorbing side cleanly (45 degrees) passes; one ending near it or running along it is refused"""
    d = out_dir("preflight")
    spec = {"name": "cap", "nets": LANE, "box": CAP_BOX, "preflight_only": True}
    b = boards_em.diagonal_lane(d)
    rc, so, se, res, _ = region(b, spec)
    need(res and rc == 0 and res["preflight"]["ok"], f"45-degree crossing refused: {(res or {}).get('preflight', {}).get('problems') or (se or so)[-300:]}")
    b = boards_em.stub_lane(d)
    rc, so, se, res, _ = region(b, spec)
    pr = (res or {}).get("preflight", {}).get("problems", [])
    need(rc == 2 and any("SLOW_B" in p and "ends" in p and "ymax" in p for p in pr), f"track ending near ymax not refused: rc {rc}, {pr}")
    b = boards_em.parallel_lane(d)
    rc, so, se, res, _ = region(b, spec)
    pr = (res or {}).get("preflight", {}).get("problems", [])
    need(rc == 2 and any("SLOW_C" in p and "runs along" in p for p in pr), f"track along ymax not refused: rc {rc}, {pr}")
    # the clearance does not shrink as the box edge comes closer (the mesh refines there): the same track
    # 0.49 mm inside the via box's ymax side, where PML + 8 cells is only ~0.36 mm, is refused
    rc, so, se, res, _ = region(b, {"name": "via", "nets": ["LANE0_RXC_P", "LANE0_RXC_N"], "box": [30.65, 6.665, 37.35, 13.085],
                                    "preflight_only": True})
    pf = (res or {}).get("preflight", {})
    need(rc == 2 and any("SLOW_C" in p for p in pf.get("problems", [])), f"track 0.49 mm from the boundary passed: {pf}")
    return f"diagonal crossing passes; end and parallel run refused; clearance floor {pf['keep_um']['ymax']:.0f} um"


# ---- 3. the planner's box search ---------------------------------------------------------------
def t_em_fit():
    """Automatic 3D regions: a planned box the pre-flight refuses is moved until it passes; when no box passes the job says why"""
    d = out_dir("fit")
    b = boards_em.stub_lane(d)
    r = evaluate(b, "--checks", "drc", "--no-drc", "--em", "plan")
    need(r["ev"], f"no evaluation: {(r['err'] or r['out'])[-300:]}")
    j = job(r["ev"], "series_cap")
    fit = j.get("fit") or {}
    need(fit.get("ok") is True and j.get("planned_box") and j["box"] != j["planned_box"], f"box not fitted: {fit}, box {j.get('box')}")
    # oracle: the pre-flight itself, through the region command, on both boxes
    spec = {k: v for k, v in j.items() if k in ("name", "nets", "fmax_ghz")}
    rc0, so, se, r0, _ = region(b, dict(spec, box=j["planned_box"], preflight_only=True))
    rc1, so, se, r1, _ = region(b, dict(spec, box=j["box"], preflight_only=True))
    need(r0 and not r0["preflight"]["ok"], f"the planned box {j['planned_box']} passes pre-flight: nothing to fit")
    need(r1 and r1["preflight"]["ok"], f"the fitted box {j['box']} fails pre-flight: {(r1 or {}).get('preflight')}")
    # no box possible: a bus of tracks along every reachable position of the ymax side
    b2 = boards_em.fenced_lane(d)
    r2 = evaluate(b2, "--checks", "drc", "--no-drc", "--em", "plan")
    j2 = job(r2["ev"], "series_cap")
    f2 = j2.get("fit") or {}
    need(f2.get("ok") is False and "BUS" in f2.get("why", "") and "searched" in f2.get("why", ""), f"fenced box: {f2}")
    fi = [f for f in r2["ev"]["findings"] if f["check"] == "em_region" and "no box passes" in f["reason"]]
    need(fi, "no finding says that no box passes")
    need("no box passes" in r2["out"] + r2["err"], "the run's log does not say that no box passes")
    return f"moved {fit['moved_mm']} mm in {fit['tries']} pre-flights; fenced: '{f2['why'][:80]}...'"


# ---- 4. grid reuse -------------------------------------------------------------------------------
FAKE_RUNNER = r'''
import json, os, sys
import numpy as np
mp, gdir = sys.argv[1], sys.argv[2]
with open(os.environ["FAKE_EM_CALLS"], "a") as fh:
    fh.write(os.path.basename(gdir) + "\n")
if os.environ.get("FAKE_EM_FAIL") == os.path.basename(gdir):
    print("fake openEMS: interrupted", file=sys.stderr)
    sys.exit(1)
m = json.load(open(mp))
f = np.linspace(0, m["fmax"], 401); f[0] = 1e6
out = {"freq_hz": f.tolist(), "cells": len(m["mesh"]["x"]) * len(m["mesh"]["y"]) * len(m["mesh"]["z"]), "ports": {}}
for p in m["ports"]:
    a = np.full(len(f), 1.0 if p["excite"] else 1e-6, complex)
    bb = (0.05 if p["excite"] else 0.9) * a[0] * np.exp(-2j * np.pi * f * 1e-10)
    out["ports"][str(p["n"])] = {"a_re": a.real.tolist(), "a_im": a.imag.tolist(), "b_re": bb.real.tolist(),
                                 "b_im": bb.imag.tolist(), "zref_re": [50.0] * len(f)}
json.dump(out, open(os.path.join(gdir, "ports.json"), "w"))
print("fake openEMS done")
'''


def t_em_grid_reuse():
    """3D runner: a rerun reuses a finished grid whose inputs hash the same; an interrupted, damaged or changed grid is solved again"""
    import boards
    d = out_dir("reuse")
    b = boards.line(d, "reuse", w=0.22)
    runs = os.path.join(d, "runs")
    calls = os.path.join(d, "calls.txt")
    fake = os.path.join(d, "fake_runner.py")
    open(fake, "w").write(FAKE_RUNNER)
    py = os.path.join(d, "fake_python")            # stands in for the openEMS Python: runs the fake runner
    open(py, "w").write(f'#!/bin/sh\n# args: -I em3d.py MODEL GRID_DIR THREADS\nexec "{sys.executable}" -I "{fake}" "$3" "$4"\n')
    os.chmod(py, os.stat(py).st_mode | stat.S_IXUSR)
    spec = {"name": "reuse", "nets": ["LINE0"], "box": [8, 6, 32, 14], "fmax_ghz": 10}
    env = {"OPENPCBSI_RUNS": runs, "OPENPCBSI_OPENEMS_PYTHON": py, "FAKE_EM_CALLS": calls}

    def go(**extra):
        n0 = len(open(calls).read().split()) if os.path.exists(calls) else 0
        rc, so, se, res, path = region(b, dict(spec, **extra.pop("spec", {})), env=dict(env, **extra))
        new = open(calls).read().split()[n0:] if os.path.exists(calls) else []
        return rc, so + se, res, new
    # without the variable: a clear error, no solve
    rc, so, se, res, _ = region(b, spec, env={"OPENPCBSI_RUNS": runs}, unset=("OPENPCBSI_OPENEMS_PYTHON",))
    need(rc == 3 and "OPENPCBSI_OPENEMS_PYTHON" in se, f"no clear error without OPENPCBSI_OPENEMS_PYTHON (rc {rc}): {se[-300:]}")
    rc, log, res, new = go(FAKE_EM_FAIL="grid1")                      # interrupted in the second grid
    need(rc != 0 and new == ["grid0", "grid1"], f"first run: rc {rc}, solved {new}")
    od = next(os.path.join(runs, "em", x) for x in os.listdir(os.path.join(runs, "em")))
    need(os.path.exists(os.path.join(od, "grid0", "done.json")) and not os.path.exists(os.path.join(od, "grid1", "done.json")),
         "grid records after the interruption are wrong")
    rc, log, res, new = go()
    need(rc == 0 and res and res.get("ran") and new == ["grid1"], f"rerun after the interruption solved {new} (rc {rc}): {log[-300:]}")
    need("reused" in log, "the rerun does not say that it reused grid 0")
    rc, log, res2, new = go()
    need(rc == 0 and new == [], f"a complete rerun solved {new}")
    need(res2.get("S11_db") == res.get("S11_db"), "a reused grid gives a different result")
    with open(os.path.join(od, "grid0", "ports.json"), "r+") as fh:   # a damaged grid: solved again
        fh.truncate(100)
    rc, log, res, new = go()
    need(rc == 0 and new == ["grid0"], f"after damaging grid 0's ports, solved {new}")
    dj = os.path.join(od, "grid1", "done.json")                       # a grid recorded for other inputs
    rec = json.load(open(dj))
    rec["key"] = "0" * 64
    json.dump(rec, open(dj, "w"))
    rc, log, res, new = go()
    need(rc == 0 and new == ["grid1"], f"after changing grid 1's input hash, solved {new}")
    rc, log, res, new = go(spec={"end_criteria": 1e-5})              # a solver setting changed
    need(rc == 0 and new == ["grid0", "grid1"], f"after changing end_criteria, solved {new}")
    return "interrupted run: grid 0 reused, grid 1 solved; complete rerun solves nothing; damaged/changed grids solved again"


def t_em_grid_reuse_openems():
    """3D runner with openEMS: the second run of a solved region reuses both grids and gives the same result (slow)"""
    if not os.environ.get("OPENPCBSI_OPENEMS_PYTHON"):
        return "blocked: OPENPCBSI_OPENEMS_PYTHON is not set"
    import boards
    d = out_dir("reuse_openems")
    b = boards.line(d, "line", w=0.22)
    env = {"OPENPCBSI_RUNS": os.path.join(d, "runs")}
    spec = {"name": "line", "nets": ["LINE0"], "box": [8, 6, 32, 14], "fmax_ghz": 10}
    t = time.time()
    rc, so, se, r1, _ = region(b, spec, env=env)
    t1 = time.time() - t
    need(rc == 0 and r1 and r1.get("ran"), f"first run failed (rc {rc}): {(se or so)[-400:]}")
    t = time.time()
    rc, so, se, r2, _ = region(b, spec, env=env)
    t2 = time.time() - t
    need(rc == 0 and r2 and r2.get("ran"), f"second run failed (rc {rc}): {(se or so)[-400:]}")
    need(so.count("reused the finished solve") == 2, f"second run did not reuse both grids: {so[-400:]}")
    need(r2.get("S21_db") == r1.get("S21_db") and r2.get("z_line_ohm") == r1.get("z_line_ohm"), "reused grids give another result")
    return f"first run {t1 / 60:.1f} min, rerun {t2:.0f} s with both grids reused"


# ---- 5. pad_plane vs ref_plane -----------------------------------------------------------------
def _near(fs, net, x, tol=1.0, sev=("error", "warning")):
    return [f for f in fs if f["net"] == net and f["severity"] in sev
            and min(abs((f["x"] or 0) - x), abs((f.get("x2") if f.get("x2") is not None else f["x"] or 0) - x)) <= tol]


def t_em_pad_void():
    """ref_plane does not fail a track over the void pad_plane asks for under its own pad; any other void, or more void than asked, still fails"""
    d = out_dir("padvoid")
    rules = "ref_plane,split_crossing,pad_plane"
    cls = '[classes.FAST]\npatterns = ["LANE*"]\nrise_ps = 35\nnyquist_ghz = 6.25\nz_diff = 100.0\n'
    on, off = os.path.join(d, "pad_plane_on.toml"), os.path.join(d, "pad_plane_off.toml")
    open(on, "w").write(cls)
    open(off, "w").write(cls + "[rules.pad_plane]\nenabled = false\n")
    b = boards_em.pad_void(d)
    R = check_rules(b, rules, on)
    pp = {f["pad"]: f["state"] for f in R["pad_plane"]["findings"]}
    need(pp.get("C1.1") == "cut-out" and pp.get("C2.1") == "cut-out", f"pad_plane does not see the voids: {pp}")
    need(not any(f["severity"] == "error" for f in R["pad_plane"]["findings"]), "pad_plane fails a voided pad")
    rp, sc = R["ref_plane"]["findings"], R["split_crossing"]["findings"]
    for net, x in (("LANE0_RX_P", 19.5), ("LANE0_RX_N", 19.5), ("LANE0_RXC_P", 20.5), ("LANE0_RXC_N", 20.5)):
        need(not _near(rp, net, x), f"ref_plane fails {net} at its own pad's void: {[f['msg'] for f in _near(rp, net, x)]}")
        need(not _near(sc, net, x, sev=("error",)), f"split_crossing at {net}'s own pad's void")
    need(_near(rp, "LANE0_RX_P", 10.0, sev=("error",)), "the slot under LANE0_RX_P at x = 10 is not failed by ref_plane")
    need(_near(sc, "LANE0_RX_P", 10.0, sev=("error",)), "the slot under LANE0_RX_P at x = 10 is not a split crossing")
    ex = R["ref_plane"]["metrics"]["LANE0_RX_P"].get("excused_own_pad_void_mm", 0)
    need(ex > 0, "no excused length recorded")
    # oracle: the same board and class with pad_plane off (nothing asks for the voids): they fail, and
    # the slot's findings are the same
    R0 = check_rules(b, "ref_plane,split_crossing", off)
    for net, x in (("LANE0_RX_P", 19.5), ("LANE0_RXC_N", 20.5)):
        need(_near(R0["ref_plane"]["findings"], net, x), f"with pad_plane off, the void under {net}'s pad is not failed")
    need(len(_near(R0["ref_plane"]["findings"], "LANE0_RX_P", 10.0)) == len(_near(rp, "LANE0_RX_P", 10.0)), "the slot is judged differently")
    # more void than pad_plane asks for (pad + 0.6 mm under C2): the excess fails, C1's does not
    b2 = boards_em.pad_void(d, "padvoid_big", c2_margin=0.6)
    R2 = check_rules(b2, rules, on)
    rp2 = R2["ref_plane"]["findings"]
    need(_near(rp2, "LANE0_RX_N", 19.0, tol=0.5, sev=("error",)), f"the excess void under C2 is not failed: {[f['msg'] for f in rp2]}")
    need(not _near(rp2, "LANE0_RX_P", 19.5), "C1's exact void is failed on the second board")
    need(not _near(R2["split_crossing"]["findings"], "LANE0_RX_N", 19.0, sev=("error",)),
         "a track ending over its own pad's (oversized) void is called a crossing")
    return f"own-pad voids excused ({ex:.2f} mm on LANE0_RX_P); the slot and the oversized void still fail; with pad_plane off they fail"


# ---- 6. no machine-specific paths ------------------------------------------------------------------
def t_no_user_paths():
    """No interpreter or folder of the author's machine in the shipped files; an unset interpreter variable is a clear error"""
    pat = re.compile(r"/(?:media|home|Users|mnt)/[A-Za-z0-9_.\-]+/|[A-Za-z]:\\\\Users\\\\")
    skip_dirs = {"kicad_monkey", "KiCadFiles", "parsers-site", "KiPI", "kipi-venv", "kipi-home", "corpus", "runs", "cache", "out",
                 "__pycache__", ".git", "build", "openpcbsi.egg-info"}
    exts = (".py", ".md", ".toml", ".sh", ".txt", ".in", ".lock", ".json", ".yml", ".cfg")
    # the files that ship: git's list (tracked, and new files not ignored) in a checkout, else a walk
    g = subprocess.run(["git", "ls-files", "--cached", "--others", "--exclude-standard"], cwd=ROOT, capture_output=True, text=True)
    if g.returncode == 0 and g.stdout.strip():
        files = [f for f in g.stdout.splitlines() if not set(f.split("/")[:-1]) & skip_dirs]
    else:
        files = []
        for top in ("openpcbsi", "tests", "docs", "boards", "deps"):
            for dp, dns, fns in os.walk(os.path.join(ROOT, top)):
                dns[:] = [x for x in dns if x not in skip_dirs]
                files += [os.path.relpath(os.path.join(dp, fn), ROOT) for fn in fns]
        files += ["README.md", "pyproject.toml"]
    hits = []
    for f in files:
        p = os.path.join(ROOT, f)
        if not (f.endswith(exts) or "." not in os.path.basename(f)) or not os.path.isfile(p):
            continue
        for i, line in enumerate(open(p, errors="replace"), 1):
            if pat.search(line) and "re.compile" not in line:
                hits.append(f"{f}:{i}")
    need(not hits, f"machine-specific paths in: {hits[:10]}")
    from openpcbsi import live, region3d
    from openpcbsi.snapshot import InputError
    for var, fn in (("OPENPCBSI_OPENEMS_PYTHON", region3d.openems_python), ("OPENPCBSI_KIPY_PYTHON", live.kipy_python)):
        old = os.environ.pop(var, None)
        try:
            fn()
            raise AssertionError(f"{var} unset, but an interpreter was found")
        except InputError as ex:
            need(var in str(ex), f"the error does not name {var}: {ex}")
        finally:
            if old is not None:
                os.environ[var] = old
    return f"no absolute user paths in {len(files)} shipped files; both interpreter variables required by name"


TESTS = [t_em_names, t_em_preflight_crossings, t_em_fit, t_em_grid_reuse, t_em_pad_void, t_no_user_paths]
SLOW = [t_em_grid_reuse_openems]


def main(argv):
    slow = "--slow" in argv
    names = [a for a in argv if not a.startswith("-")]
    sel = [t for t in TESTS + SLOW if (t.__name__ in names) or (not names and (t in TESTS or slow))]
    bad = 0
    for t in sel:
        t0 = time.time()
        try:
            note, st = t(), "PASS"
        except AssertionError as e:
            note, st, bad = str(e), "FAIL", bad + 1
        except Exception as e:
            note, st, bad = f"{type(e).__name__}: {e}", "ERROR", bad + 1
        title = (t.__doc__ or "").strip().splitlines()[0]
        print(f"{st:6s} {time.time() - t0:7.1f}s  {t.__name__}: {title}\n               {note[:400]}", flush=True)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
