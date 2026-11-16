"""openpcbsi acceptance tests: one or more per gap item (#1-#36). Each test builds its own boards with a
known answer (tests/harness/boards.py), runs openpcbsi as a user would, and compares with an oracle:
the planted defect, a closed-form result, an independent solver, or a physics invariant.

    python3 tests/acceptance.py                 # all fast tests
    python3 tests/acceptance.py --slow          # also the 3D (openEMS) tests, which take minutes each
    python3 tests/acceptance.py --items 1,3,5   # some items
    python3 tests/acceptance.py --list

Writes tests/out/acceptance.json and docs/acceptance.md (the status table). Exit code 1 if any test fails.
Status: pass, fail, blocked (needs an input this machine does not have), error (the test itself broke).
"""
import argparse
import glob
import json
import math
import os
import re
import shutil
import subprocess
import sys
import time
import traceback

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(HERE, "harness"))
OUT = os.path.join(HERE, "out")
CASES = []


def case(item, title, slow=False):
    def deco(fn):
        CASES.append({"item": item, "title": title, "slow": slow, "fn": fn, "name": fn.__name__})
        return fn
    return deco


class Fail(Exception):
    pass


def need(cond, msg):
    if not cond:
        raise Fail(msg)


def blocked(msg):
    return {"status": "blocked", "note": msg}


def ok(note=""):
    return {"status": "pass", "note": note}


# ---- running openpcbsi ---------------------------------------------------------------------------
def run(args, timeout=1800, env=None):
    t = time.time()
    p = subprocess.run([sys.executable, "-m", "openpcbsi"] + args, cwd=ROOT, capture_output=True, text=True,
                       timeout=timeout, env=dict(os.environ, **(env or {})))
    return p.returncode, p.stdout, p.stderr, time.time() - t


def evaluate(board, *args, out=None, timeout=1800):
    out = out or os.path.join(os.path.dirname(board), "run")
    if os.path.isdir(out):
        shutil.rmtree(out)
    rc, so, se, dt = run(["evaluate", board, "-o", out, "--no-markers"] + list(args), timeout=timeout)
    ev = _load(os.path.join(out, "evaluation.json"))
    rep = _load(os.path.join(out, "report.json"))
    return {"rc": rc, "out": so, "err": se, "ev": ev, "rep": rep, "dir": out, "seconds": dt}


def check(board, checks, *args, out=None, timeout=1800):
    out = out or os.path.join(os.path.dirname(board), "check_" + checks.replace(",", "_"))
    if os.path.isdir(out):
        shutil.rmtree(out)
    rc, so, se, dt = run(["check", board, "--checks", checks, "-o", out, "--no-markers"] + list(args), timeout=timeout)
    return {"rc": rc, "out": so, "err": se, "rep": _load(os.path.join(out, "report.json")), "dir": out, "seconds": dt}


def _load(p):
    return json.load(open(p)) if os.path.exists(p) else None


def toml(path, text):
    open(path, "w").write(text)
    return path


def findings(rep, rule=None, sev=None):
    out = []
    for r in (rep or {}).get("checks", {}).get("drc", []):
        if rule and r["rule"] != rule:
            continue
        out += [f for f in r["findings"] if sev is None or f["severity"] in sev]
    return out


def ev_findings(ev, check=None, sev=None):
    return [f for f in (ev or {}).get("findings", []) if (check is None or f["check"] == check)
            and (sev is None or f["severity"] in sev)]


def out_dir(name):
    d = os.path.join(OUT, name)
    if os.path.isdir(d):
        shutil.rmtree(d)
    os.makedirs(d)
    return d


import boards  # noqa: E402  (needs pcbnew: run with KiCad's python3)


# =================================================================================================
# A. KiCad input
# =================================================================================================
@case(1, "Newer KiCad file: clear error naming the version; KiCad version recorded")
def t01_version():
    d = out_dir("t01")
    new = boards.pair_over_plane(d, "newer", version=20260901)
    r = evaluate(new, "--checks", "drc", "--no-drc")
    need(r["rc"] == 3, f"exit code {r['rc']}, want 3")
    msg = r["err"] + r["out"]
    need("20260901" in msg, "the error does not name the file version")
    need(re.search(r"KiCad 10\.0", msg) and re.search(r"save", msg, re.I), "the error does not say to save as KiCad 10.0")
    need("Traceback" not in r["err"], "a traceback instead of a clear error")
    ok_b = boards.pair_over_plane(d, "current")
    r2 = evaluate(ok_b, "--checks", "drc", "--no-drc")
    need(r2["rep"] and r2["rep"].get("kicad_file_version"), "report has no kicad_file_version")
    need(r2["ev"] and r2["ev"].get("kicad_file_version"), "evaluation has no kicad_file_version")
    return ok(f"exit 3 with version named; file version {r2['rep']['kicad_file_version']} recorded")


@case(2, "Netclasses (.kicad_pro patterns) and .kicad_dru limits are used without a config")
def t02_project():
    d = out_dir("t02")
    b = boards.pair_uncoupled(d, "proj", "HSL0_TX", split_len=3.0,
                              netclass=("HS", ["HSL*"], {"diff_pair_width": 0.12, "diff_pair_gap": 0.15}),
                              dru=['(rule hs_unc (condition "A.NetClass == \'HS\'") (constraint diff_pair_uncoupled (max 2mm)))'])
    r = evaluate(b, "--checks", "drc", "--no-drc")
    need(r["rep"], f"no report: {r['err'][-300:]}")
    nets = (r["rep"].get("netclasses") or {})
    need(nets.get("HSL0_TX_P") == "HS", f"netclass of HSL0_TX_P is {nets.get('HSL0_TX_P')}, want HS (from .kicad_pro patterns)")
    det = json.dumps(r["ev"].get("detected") or {})
    need("HS" in det, "detection does not name the HS netclass as a class")
    pc = [f for f in findings(r["rep"], "pair_coupling") if f["net"].startswith("HSL0")]
    need(pc, "no pair_coupling finding for the 3 mm split")
    lim = {f.get("limit") for f in pc}
    need(2.0 in lim, f"pair_coupling limit {lim}, want 2.0 from .kicad_dru")
    need(any(f["severity"] == "error" for f in pc), "3 mm uncoupled against a 2 mm limit is not an error")
    return ok("class HS from .kicad_pro, uncoupled limit 2.0 mm from .kicad_dru")


@case(3, "Zone fills are refilled before extraction (a board saved unfilled gives the right answer)")
def t03_refill():
    d = out_dir("t03")
    clean = boards.pair_over_plane(d, "clean_unfilled", fill=False)
    slot = boards.pair_over_plane(d, "slot_unfilled", slot=True, fill=False)
    rc = check(clean, "drc", "--rules", "ref_plane,split_crossing")
    rs = check(slot, "drc", "--rules", "ref_plane,split_crossing")
    need(rc["rep"] and rs["rep"], "no report")
    ec = findings(rc["rep"], sev=("error",))
    es = findings(rs["rep"], sev=("error",))
    need(not ec, f"clean board saved unfilled gives {len(ec)} errors (fills not refilled)")
    need(len(es) == 4, f"slotted board gives {len(es)} errors, want 4")
    need(rc["rep"].get("zones_refilled") is True, "report does not record zones_refilled")
    return ok("0 errors clean, 4 at the slot, zones_refilled recorded")


@case(4, "Design data from symbol fields is used and marked sourced; guesses are marked assumed")
def t04_fields():
    d = out_dir("t04")
    b = boards.rail_board(d, "rail", current_field="1.5A")
    r = evaluate(b, "--checks", "dc", "--no-drc")
    need(r["rep"], f"no report: {r['err'][-300:]}")
    rails = (r["rep"].get("config_used") or {}).get("rails") or {}
    r18 = next((v for k, v in rails.items() if v.get("net") == "+1V8"), None)
    need(r18, f"+1V8 rail not detected (rails {list(rails)})")
    need(abs(r18.get("current_a", 0) - 1.5) < 1e-6, f"+1V8 current {r18.get('current_a')}, want 1.5 from U2's Current field")
    need("U2" in str(r18.get("current_source", "")), f"current_source {r18.get('current_source')}, want U2.Current")
    r33 = next((v for k, v in rails.items() if v.get("net") == "+3V3"), None)
    need(r33 and r33.get("current_assumed") is True, "+3V3 (no field) is not marked current_assumed")
    return ok("1.5 A from U2.Current (sourced); +3V3 assumed")


# =================================================================================================
# B. Verdicts
# =================================================================================================
def _cfg_class(d, name, z=100.0, tol=10.0, extra=""):
    return toml(os.path.join(d, name + ".toml"), f"""
[classes.HS]
patterns = ["LANE*"]
rise_ps = 35
nyquist_ghz = 6.25
z_diff = {z}
z_tol = {tol}
{extra}
""")


@case(5, "Every computed judgement that is out of limit has a finding; 'nothing to check' is not 'pass'")
def t05_judged():
    d = out_dir("t05")
    notes = []
    # (a) impedance: a 0.30 / 0.10 pair is far below 100 ohm
    b = boards.pair_over_plane(d, "wide", w=0.30, gap=0.10)
    r = evaluate(b, "-c", _cfg_class(d, "wide"), "--checks", "impedance", "--no-drc")
    need(r["ev"], f"no evaluation: {r['err'][-300:]}")
    rows = r["rep"]["checks"]["impedance"]["rows"]
    off = [x for x in rows if x.get("within_target") is False]
    need(off, "the seeded off-target geometry was not computed as off target")
    fi = ev_findings(r["ev"], "impedance", ("fail", "warn"))
    need(len(fi) >= len(off), f"{len(off)} geometries off target, {len(fi)} impedance findings")
    notes.append(f"impedance {len(fi)} findings")
    # (b) thermal: 3 W at theta_jb 20 on a part with Tj_max = 85
    b = boards.thermal_part(d)
    cfg = toml(os.path.join(d, "th.toml"), "[thermal.parts]\nU1 = { power_w = 3.0, theta_jb = 20.0 }\n")
    r = evaluate(b, "-c", cfg, "--checks", "thermal", "--no-drc")
    need(r["ev"], f"no evaluation: {r['err'][-300:]}")
    tj = r["rep"]["checks"]["thermal"]["parts"][0]["Tj_est_c"]
    ft = ev_findings(r["ev"], "thermal", ("fail", "warn"))
    need(tj > 85 and ft, f"Tj {tj} C against Tj_max 85 gives {len(ft)} findings")
    notes.append(f"thermal Tj {tj:.0f} C -> {ft[0]['severity']}")
    # (c) a board none of whose nets is in a signal class: class-scoped rules checked nothing
    from gen import Board
    p = os.path.join(d, "noclass", "noclass.kicad_pcb")
    os.makedirs(os.path.dirname(p), exist_ok=True)
    bb = Board(p, layers=2, size=(40, 20))
    bb.track("F.Cu", [(5, 8), (35, 8)], 0.3, "ENABLE")
    bb.track("F.Cu", [(5, 12), (35, 12)], 0.3, "RESET")
    bb.save()
    r = evaluate(p, "--checks", "drc", "--no-drc")
    need(r["ev"], f"no evaluation: {r['err'][-300:]}")
    st = {c["check"]: c["status"] for c in r["ev"]["checks"]}
    need(st.get("ref_plane") == "nothing to check", f"ref_plane with no net in scope: {st.get('ref_plane')}")
    fam = {F["family"]: F["status"] for F in r["ev"]["families"]}
    need(fam.get("RP") != "pass", f"RP family with nothing in scope: {fam.get('RP')}")
    notes.append("no net in scope: nothing to check")
    # (d) a channel mask the channel cannot meet
    b = boards.pair_over_plane(d, "chan")
    cfg = _cfg_class(d, "chan", extra='[channels.C0]\npath = ["LANE0_TX_P"]\nbitrate_gbps = 12.5\nrise_ps = 35\n'
                                     'masks = { sdd11_db = -60.0, sdd11_fmax_ghz = 6.25 }\n')
    r = evaluate(b, "-c", cfg, "--checks", "si", "--no-drc")
    need(r["ev"], f"no evaluation: {r['err'][-300:]}")
    fs = ev_findings(r["ev"], "si", ("fail",))
    need(fs, "a channel failing its SDD11 mask gives no si finding")
    notes.append("channel mask fail -> finding")
    return ok("; ".join(notes))


@case(6, "Interface standards: the same layout is judged by the standard its nets belong to")
def t06_standards():
    d = out_dir("t06")
    res = {}
    for nm in ("PCIE_TX0", "LVDS_D0"):
        b = boards.pair_uncoupled(d, nm.lower(), nm, split_len=3.0)
        r = evaluate(b, "--checks", "drc", "--no-drc")
        need(r["rep"], f"no report for {nm}: {r['err'][-200:]}")
        cls = (r["rep"].get("standards") or {}).get(nm + "_P")
        pc = [f for f in findings(r["rep"], "pair_coupling") if f["net"].startswith(nm)]
        res[nm] = (cls, {f["severity"] for f in pc}, {f.get("limit") for f in pc})
    need(res["PCIE_TX0"][0] and "pcie" in str(res["PCIE_TX0"][0]).lower(), f"PCIE nets not mapped to a PCIe standard: {res['PCIE_TX0'][0]}")
    need(res["LVDS_D0"][0] and "lvds" in str(res["LVDS_D0"][0]).lower(), f"LVDS nets not mapped to LVDS: {res['LVDS_D0'][0]}")
    need(res["PCIE_TX0"][2] != res["LVDS_D0"][2], f"same uncoupled limits for both standards: {res}")
    return ok(f"PCIe limit {res['PCIE_TX0'][2]}, LVDS limit {res['LVDS_D0'][2]}")


def _hj(w, h, t, er):
    """Hammerstad-Jensen microstrip Z0 (thickness-corrected width), the closed-form oracle."""
    if t > 0:
        w = w + t / math.pi * (1 + math.log(4 * math.e / math.sqrt((t / h) ** 2 + (1 / math.pi / (w / t + 1.1)) ** 2)))
    u = w / h
    a = 1 + math.log((u ** 4 + (u / 52) ** 2) / (u ** 4 + 0.432)) / 49 + math.log(1 + (u / 18.1) ** 3) / 18.7
    b = 0.564 * ((er - 0.9) / (er + 3)) ** 0.053
    ee = (er + 1) / 2 + (er - 1) / 2 * (1 + 10 / u) ** (-a * b)
    f = 6 + (2 * math.pi - 6) * math.exp(-(30.666 / u) ** 0.7528)
    return 60 / math.sqrt(ee) * math.log(f / u + math.sqrt(1 + (2 / u) ** 2))


@case(7, "Tolerance corners bracket the nominal impedance and match closed-form sensitivity")
def t07_corners():
    d = out_dir("t07")
    b = boards.line(d, "se", w=0.22)
    cfg = toml(os.path.join(d, "c.toml"), '[classes.SE]\npatterns = ["LINE*"]\nz_se = 50.0\nz_tol = 2.0\n'
                                          '[xsec]\nmask_er = 1.0\n[xsec.corners]\netch_um = 20\nthickness_pct = 10\ner_pct = 5\n')
    r = evaluate(b, "-c", cfg, "--checks", "impedance", "--no-drc")
    need(r["rep"], f"no report: {r['err'][-300:]}")
    rows = [x for x in r["rep"]["checks"]["impedance"]["rows"] if x.get("layer") == "F.Cu" and x.get("kind", "se") == "se"]
    need(rows and "Z0_min" in rows[0], "impedance rows carry no Z0_min / Z0_max corners")
    x = rows[0]
    need(x["Z0_min"] < x["Z0"] < x["Z0_max"], f"corners do not bracket: {x['Z0_min']} {x['Z0']} {x['Z0_max']}")
    # closed-form spread at the same corners (no mask)
    h, t, er = 0.12, 0.035, 4.2
    hj = [_hj(0.22 + dw, h * fh, t, er * fe) for dw, fh, fe in ((0.04, 0.9, 1.05), (-0.04, 1.1, 0.95))]
    spread_hj, spread = hj[1] - hj[0], x["Z0_max"] - x["Z0_min"]
    need(abs(spread - spread_hj) / spread_hj < 0.3, f"corner spread {spread:.2f} ohm vs closed form {spread_hj:.2f} ohm")
    need(ev_findings(r["ev"], "impedance", ("warn", "fail")), "corners leave the 2 ohm tolerance but no finding")
    return ok(f"Z0 {x['Z0_min']:.1f}-{x['Z0_max']:.1f} ohm, closed form spread {spread_hj:.1f} ohm")


@case(8, "Correlation with a measured board: Antmicro's SI test board (published VNA data)")
def t08_measure():
    rc, so, se, dt = run(["correlate", "--fetch"], timeout=1800)
    out = _load(os.path.join(ROOT, "runs", "correlation", "antmicro_si.json"))
    need(out, f"no correlation result (rc {rc}): {(se or so)[-400:]}")
    lines = []
    for c in out["cases"]:
        need(c["stubs"], f"{c['case']}: the model has no stub (path extraction)")
        for k in c["checks"]:
            lines.append(f"{c['case']} {k['feature']} {k['error_pct']:+.1f} %")
            need(k["pass"], f"{c['case']} {k['feature']}: model {k['model_ghz']} vs measured {k['measured_ghz']} GHz "
                            f"({k['error_pct']:+.1f} %, limit {out['tolerance_pct']:g} %)")
    return ok("; ".join(lines) + " (one board, 0.2-6 GHz, outer-layer microstrip; see README)")


# =================================================================================================
# C. 3D regions (openEMS)
# =================================================================================================
def _region(board, spec_text, cfg=None, timeout=5400):
    d = os.path.dirname(board)
    spec = toml(os.path.join(d, "region.toml"), spec_text)
    args = ["region", board, spec] + (["-c", cfg] if cfg else [])
    rc, so, se, dt = run(args, timeout=timeout)
    m = re.search(r"region result: (\S+)", so)
    res = _load(m.group(1)) if m else None
    return rc, so, se, res, dt


@case(9, "3D line ports: a straight microstrip reads its 2D impedance within 2 % and is passive", slow=True)
def t09_ports():
    d = out_dir("t09")
    b = boards.line(d, "line", w=0.22)
    rc, so, se, res, dt = _region(b, '[region]\nname = "line"\nnets = ["LINE0"]\nbox = [8, 6, 32, 14]\nfmax_ghz = 20\n')
    need(res, f"no region result (rc {rc}): {(se or so)[-400:]}")
    need(res.get("ports_kind") == "line", f"ports are {res.get('ports_kind')}, want transmission-line ports")
    need(res.get("passive") is True, "result is not passive")
    z, z2 = res.get("z_line_ohm"), res.get("z_2d_ohm")
    need(z and z2 and abs(z - z2) / z2 < 0.02, f"3D {z} ohm vs 2D {z2} ohm")
    return ok(f"3D {z:.1f} vs 2D {z2:.1f} ohm in {dt / 60:.0f} min")


@case(10, "3D differential drive: Sdd and Scd from one run; a straight pair has Scd21 < -40 dB", slow=True)
def t10_mixed():
    d = out_dir("t10")
    b = boards.pair_over_plane(d, "pair")
    rc, so, se, res, dt = _region(b, '[region]\nname = "pair"\nnets = ["LANE0_TX_P", "LANE0_TX_N"]\nbox = [8, 7, 32, 13]\n')
    need(res, f"no region result (rc {rc}): {(se or so)[-400:]}")
    need(all(k in res for k in ("Sdd11_db", "Sdd21_db", "Scd21_db")), "no mixed-mode parameters")
    need(max(res["Scd21_db"]) < -40, f"Scd21 max {max(res['Scd21_db']):.1f} dB on a symmetric pair")
    need(res.get("runs_per_grid") == 1, "more than one run per grid for a differential pair")
    need(res.get("passive") is True and max(res["Sdd21_db"]) <= 0.21, "result is not passive")
    return ok(f"Scd21 max {max(res['Scd21_db']):.1f} dB")


@case(11, "3D series parts: an AC-coupled pair transmits (the capacitor is modelled)", slow=True)
def t11_parts():
    d = out_dir("t11")
    from gen import Board
    p = os.path.join(d, "accap", "accap.kicad_pcb")
    os.makedirs(os.path.dirname(p), exist_ok=True)
    b = Board(p, layers=4, size=(40, 20))
    b.plane("In1.Cu", "GND"); b.plane("In2.Cu", "GND")
    for sfx, y in (("_P", 9.0), ("_N", 11.0)):
        b.track("F.Cu", [(5, y), (19.5, y)], 0.12, "A" + sfx)
        b.two_pin("C1" if sfx == "_P" else "C2", "100n", 20, y, "A" + sfx, "B" + sfx, pitch=1.0, pad=(0.6, 0.55))
        b.track("F.Cu", [(20.5, y), (35, y)], 0.12, "B" + sfx)
    b.save()
    rc, so, se, res, dt = _region(p, '[region]\nname = "accap"\nnets = ["A_P", "A_N", "B_P", "B_N"]\nbox = [12, 5, 28, 15]\nfmax_ghz = 10\n')
    need(res, f"no region result (rc {rc}): {(se or so)[-400:]}")
    need(res.get("series_parts"), "no series parts in the model")
    sdd21 = dict(zip(res["freq_ghz"], res["Sdd21_db"]))
    f5 = min(sdd21, key=lambda f: abs(f - 5))
    need(res.get("passive") is True, f"result is not passive (max Sdd21 {max(res['Sdd21_db']):.1f} dB)")
    need(res.get("convergence", {}).get("converged") is True, f"not converged: {res.get('convergence', {}).get('reason')}")
    need(-3 < sdd21[f5] <= 0.21, f"Sdd21 at 5 GHz {sdd21[f5]:.2f} dB (open, or not passive?)")
    return ok(f"Sdd21(5 GHz) {sdd21[f5]:.2f} dB; passive, converged")


@case(12, "3D clipping keeps the signal path: a layer change through a via still transmits", slow=True)
def t12_clip():
    d = out_dir("t12")
    from gen import Board
    p = os.path.join(d, "lc", "lc.kicad_pcb")
    os.makedirs(os.path.dirname(p), exist_ok=True)
    b = Board(p, layers=4, size=(40, 20))
    b.plane("In1.Cu", "GND"); b.plane("In2.Cu", "GND")
    b.track("F.Cu", [(5, 10), (20, 10)], 0.22, "S")
    b.via(20, 10, "S")
    b.track("B.Cu", [(20, 10), (35, 10)], 0.22, "S")
    for x, y in ((19, 9), (21, 11)):
        b.via(x, y, "GND")
    b.save()
    rc, so, se, res, dt = _region(p, '[region]\nname = "lc"\nnets = ["S"]\nbox = [12, 6, 28, 14]\nfull_stack = true\nfmax_ghz = 10\n')
    need(res, f"no region result (rc {rc}): {(se or so)[-400:]}")
    need(res.get("path_connected") is True, "the clipped model is not connected port to port")
    s21 = dict(zip(res["freq_ghz"], res["S21_db"]))
    f5 = min(s21, key=lambda f: abs(f - 5))
    need(res.get("passive") is True and -10 < s21[f5] <= 0.21, f"S21 at 5 GHz {s21[f5]:.1f} dB, passive {res.get('passive')}")
    return ok(f"connected, S21(5 GHz) {s21[f5]:.2f} dB")


@case(13, "3D pre-flight refuses a bad setup; full stack with absorbing boundaries runs")
def t13_preflight():
    d = out_dir("t13")
    b = boards.line(d, "pf", w=0.22)
    rc, so, se, res, dt = _region(b, '[region]\nname = "pf"\nnets = ["LINE0"]\nbox = [8, 9.6, 32, 10.4]\npreflight_only = true\n')
    need(rc != 0 and re.search(r"pre-?flight", so + se, re.I), f"a 0.6 mm wide box with copper at the boundary was not refused (rc {rc})")
    rc, so, se, res, dt = _region(b, '[region]\nname = "pf_ok"\nnets = ["LINE0"]\nbox = [8, 6, 32, 14]\npreflight_only = true\n')
    need(rc == 0 and res and res.get("preflight", {}).get("ok") is True, f"a good box failed pre-flight: {(res or {}).get('preflight')}")
    need(res.get("boundaries") and all(x.startswith("PML") for x in res["boundaries"][:4]), f"boundaries {res.get('boundaries')}")
    return ok("bad box refused, good box passes with PML sides")


@case(14, "3D mesh: >= 4 cells across tracks and gaps, no lines closer than 5 um")
def t14_mesh():
    d = out_dir("t14")
    b = boards.pair_over_plane(d, "mesh")
    rc, so, se, res, dt = _region(b, '[region]\nname = "mesh"\nnets = ["LANE0_TX_P", "LANE0_TX_N"]\nbox = [8, 7, 32, 13]\npreflight_only = true\n')
    need(res and res.get("mesh"), f"no mesh report (rc {rc}): {(se or so)[-300:]}")
    m = res["mesh"]
    need(m["min_cells_across_track"] >= 4 and m["min_cells_across_gap"] >= 4, f"cells across track/gap {m['min_cells_across_track']}/{m['min_cells_across_gap']}")
    need(m["min_spacing_um"] >= 5.0, f"mesh lines {m['min_spacing_um']} um apart")
    return ok(f"{m['min_cells_across_track']} cells across track, {m['min_cells_across_gap']} across gap")


@case(15, "3D convergence: each run records how it ended; a run cut short is not accepted", slow=True)
def t15_conv():
    d = out_dir("t15")
    b = boards.line(d, "conv", w=0.22)
    rc, so, se, res, dt = _region(b, '[region]\nname = "short"\nnets = ["LINE0"]\nbox = [8, 6, 32, 14]\nmax_timesteps = 2000\n')
    need(res, f"no region result (rc {rc}): {(se or so)[-300:]}")
    ends = [x.get("ended_by") for x in res.get("convergence", {}).get("runs", [])]
    need(ends and all(e == "timesteps" for e in ends), f"ended_by {ends}")
    need(res["convergence"].get("converged") is False and "step" in res["convergence"].get("reason", ""), f"{res['convergence']}")
    need("z_line_diff_pct" in res["convergence"] or "tdr_diff_pct" in res["convergence"], "convergence does not compare Z/TDR between grids")
    return ok("runs cut short are flagged and not accepted")


@case(16, "3D calibration: each class/layer line is solved in the same mesh and compared with 2D", slow=True)
def t16_cal():
    d = out_dir("t16")
    b = boards.pair_over_plane(d, "cal")
    rc, so, se, res, dt = _region(b, '[region]\nname = "cal"\nnets = ["LANE0_TX_P", "LANE0_TX_N"]\nbox = [8, 7, 32, 13]\n')
    need(res and res.get("calibration"), f"no calibration record (rc {rc})")
    c = res["calibration"]
    need("z_model_ohm" in c and "z_2d_ohm" in c and "agree_pct" in c, f"{c}")
    need(c.get("renormalised") == (abs(c["agree_pct"]) > 2.0), "renormalisation does not follow the 2 % gate")
    return ok(f"model line {c['z_model_ohm']:.1f} vs 2D {c['z_2d_ohm']:.1f} ohm")


@case(17, "TDR recovers a known impedance step within 2 %")
def t17_tdr():
    sys.path.insert(0, ROOT)
    from openpcbsi import tdr as T
    import numpy as np
    f = np.linspace(0, 20e9, 401)
    z0, z1, tau = 100.0, 80.0, 200e-12             # a 100 ohm reference into an 80 ohm line, then matched
    g = (z1 - z0) / (z1 + z0)
    S11 = g + 0 * f                                  # a step to 80 ohm (infinite line)
    t, z = T.tdr(f, S11, Z0=z0, rise=35e-12)
    zr = z[(t > 100e-12) & (t < 300e-12)].mean()
    need(abs(zr - z1) / z1 < 0.02, f"TDR {zr:.2f} ohm, want {z1} (step at the port)")
    need(len(t) > 16 * 400, "TDR is not zero-padded")
    # a 100 ohm line of 150 ps round trip, then the 80 ohm step
    S11d = g * np.exp(-2j * np.pi * f * 150e-12)
    t, z = T.tdr(f, S11d, Z0=z0, rise=35e-12)
    before, after = z[(t > 20e-12) & (t < 110e-12)].mean(), z[(t > 220e-12) & (t < 400e-12)].mean()
    need(abs(before - z0) / z0 < 0.02 and abs(after - z1) / z1 < 0.02, f"delayed step: {before:.1f} then {after:.1f} ohm")
    return ok(f"step at the port {zr:.2f} ohm; delayed step {before:.1f} -> {after:.1f} ohm")


@case(18, "3D lossy materials: microstrip loss within 25 % of the closed form", slow=True)
def t18_loss():
    d = out_dir("t18")
    b = boards.line(d, "loss", w=0.22)
    rc, so, se, res, dt = _region(b, '[region]\nname = "loss"\nnets = ["LINE0"]\nbox = [6, 6, 34, 14]\nlossy = true\nfmax_ghz = 20\n')
    need(res and res.get("lossy") is True, f"no lossy result (rc {rc})")
    s21 = dict(zip(res["freq_ghz"], res["S21_db"]))
    f10 = min(s21, key=lambda f: abs(f - 10))
    L = res.get("line_length_mm")
    ee = res.get("eps_eff") or 3.3
    need(L, "no line length")
    need(res.get("passive") is True, "result is not passive")
    # closed form at 10 GHz: dielectric (filling-factor) loss plus parallel-plate conductor loss (strip and
    # plane), w = 0.22 mm, Z0 from the 2D solver; edge crowding adds some conductor loss (within the 25 %)
    er, tand, f = 4.2, 0.02, 10e9
    alpha_d = 27.3 * (er / math.sqrt(ee)) * ((ee - 1) / (er - 1)) * tand * f / 2.998e8 / 1000      # dB/mm
    rs = math.sqrt(math.pi * f * 4e-7 * math.pi / 5.8e7)
    z0 = res.get("z_2d_ohm") or 47.0
    alpha_c = 8.686 * rs / (z0 * 0.22e-3) / 1000          # dB/mm: R' = 2 Rs / w (strip and plane), alpha = R' / 2 Z0
    want = -(alpha_d + alpha_c) * L
    need(abs(s21[f10] - want) / abs(want) < 0.25, f"S21(10 GHz) {s21[f10]:.3f} dB vs closed form {want:.3f} dB "
                                                   f"(dielectric {-alpha_d * L:.3f}, conductor {-alpha_c * L:.3f})")
    return ok(f"S21(10 GHz) {s21[f10]:.3f} dB vs closed form {want:.3f} dB over {L} mm")


@case(19, "Automatic 3D regions with fix variants at each discontinuity type")
def t19_auto():
    d = out_dir("t19")
    from gen import Board
    p = os.path.join(d, "auto", "auto.kicad_pcb")
    os.makedirs(os.path.dirname(p), exist_ok=True)
    bb = Board(p, layers=4, size=(50, 20))
    bb.plane("In1.Cu", "GND"); bb.plane("In2.Cu", "GND")
    for sfx, y in (("_P", 9.865), ("_N", 10.135)):
        bb.track("F.Cu", [(5, y), (19.5, y)], 0.12, "LANE0_RX" + sfx)
    bb.two_pin("C1", "100n", 20, 9.0, "LANE0_RX_P", "LANE0_RXC_P", pitch=1.0, pad=(0.5, 0.5))
    bb.two_pin("C2", "100n", 20, 11.0, "LANE0_RX_N", "LANE0_RXC_N", pitch=1.0, pad=(0.5, 0.5))
    for sfx, y in (("_P", 9.865), ("_N", 10.135)):
        bb.track("F.Cu", [(20.5, y), (34, y)], 0.12, "LANE0_RXC" + sfx)
        bb.via(34, y + (-0.4 if sfx == "_P" else 0.4), "LANE0_RXC" + sfx)
        bb.track("B.Cu", [(34, y + (-0.4 if sfx == "_P" else 0.4)), (45, y)], 0.12, "LANE0_RXC" + sfx)
    bb.save()
    r = evaluate(p, "--checks", "drc", "--no-drc", "--em", "plan")
    need(r["ev"], f"no evaluation: {r['err'][-300:]}")
    jobs = r["ev"].get("em_jobs")
    need(jobs is not None, "evaluation has no em_jobs (automatic regions not built)")
    kinds = {j["kind"] for j in jobs}
    need({"series_cap", "layer_change"} <= kinds, f"planned kinds {sorted(kinds)}")
    need(all(j.get("variant") for j in jobs if j["kind"] in ("series_cap", "layer_change")), "a job has no fix variant")
    # each planned job builds a model that passes pre-flight
    for j in jobs:
        spec = {k: v for k, v in j.items() if k not in ("kind", "parts")}
        spec["preflight_only"] = True
        rc, so, se, res, dt = _region(p, "[region]\n" + "\n".join(f"{k} = {json.dumps(v)}" for k, v in spec.items()
                                                                 if not isinstance(v, dict)) + "\n")
        need(res and res.get("preflight", {}).get("ok"), f"job {j['name']}: pre-flight {(res or {}).get('preflight')}")
    return ok(f"{len(jobs)} regions planned ({', '.join(sorted(kinds))}), each with its fix variant, all pass pre-flight")


@case(20, "3D regression cases pass: calibration line, coupled pair, symmetric via", slow=True)
def t20_reg():
    rc, so, se, dt = run(["validate", "--3d"], timeout=7200)
    fails = [l.strip() for l in so.splitlines() if "[FAIL]" in l]
    need(rc == 0, f"validate --3d exit {rc}: " + ("; ".join(fails) if fails else (se or so)[-400:]))
    for name in ("calibration line", "coupled pair", "symmetric via"):
        need(re.search(rf"\[PASS\].*{name}", so, re.I), f"3D case '{name}' did not pass")
    return ok(f"three 3D cases pass in {dt / 60:.0f} min")


# =================================================================================================
# D. Coverage
# =================================================================================================
@case(21, "Multi-terminal nets (T and fly-by) are extracted, one result per receiver")
def t21_topo():
    d = out_dir("t21")
    b = boards.topology(d)
    cfg = toml(os.path.join(d, "t.toml"), '[classes.S]\npatterns = ["TNET", "FLYBY"]\nrise_ps = 200\n[si]\nnets = ["TNET", "FLYBY"]\n')
    r = evaluate(b, "-c", cfg, "--checks", "si", "--no-drc")
    need(r["rep"], f"no report: {r['err'][-300:]}")
    rows = r["rep"]["checks"]["si"]["rows"]
    probs = [x for x in rows if "problem" in x]
    need(not probs, f"nets skipped: {[x.get('problem') for x in probs]}")
    rx = {(x["net"].split("/")[0], x.get("to")) for x in rows}
    need(len([x for x in rx if x[0] == "TNET"]) == 2 and len([x for x in rx if x[0] == "FLYBY"]) == 3,
         f"receivers per net: {sorted(rx)}")
    return ok("T: 2 receivers, fly-by: 3 receivers")


@case(22, "PDN: target judged only inside the rail's band; package/die model and load step reported")
def t22_pdn():
    d = out_dir("t22")
    b = boards.rail_board(d, "pdn", second_rail=False)
    cfg = toml(os.path.join(d, "p.toml"), """
[rails.V18]
net = "+1V8"
v = 1.8
tol_pct = 3.0
current_a = 1.0
source = "U1"
loads = ["U2"]
ripple_pct = 3.0
step_a = 0.5
band_hz = [1e4, 2e7]
package = { l_nh = 0.1, c_die_nf = 50.0 }
[engines]
dc = "openpcbsi"
pdn = "openpcbsi"
""")
    r = evaluate(b, "-c", cfg, "--checks", "pdn", "--no-drc")
    need(r["rep"], f"no report: {r['err'][-300:]}")
    row = r["rep"]["checks"]["pdn"]["rails"][0]
    need(row.get("band_hz") == [1e4, 2e7], f"band not used: {row.get('band_hz')}")
    need(row.get("f_at_max_hz", 0) <= 2e7 * 1.01, f"worst point {row.get('f_at_max_hz')} outside the band")
    need("package" in row and "load_step_mv" in row, "no package model or load-step result")
    # closed-form load steps: a resistance droops step x R; an inductance L x dI / t_rise
    sys.path.insert(0, ROOT)
    import numpy as np
    from openpcbsi.pdn import load_step
    fq = np.logspace(5, 9, 200)
    r1 = load_step(fq, np.full(len(fq), 0.01 + 0j), 2.0, 1.0)["load_step_mv"]
    need(abs(r1 - 20.0) / 20.0 < 0.02, f"10 mohm x 2 A: {r1:.2f} mV, want 20")
    r2 = load_step(fq, 1j * 2 * np.pi * fq * 1e-9, 1.0, 10.0)["load_step_mv"]
    need(abs(r2 - 100.0) / 100.0 < 0.1, f"1 nH x 1 A / 10 ns: {r2:.1f} mV, want 100")
    return ok(f"judged in band; load step R {r1:.1f} mV (20), L {r2:.0f} mV (100); fixture {row['load_step_mv']:.0f} mV (no capacitors)")


@case(23, "Board-wide stitching map: percentiles and the largest via-free area match the construction")
def t23_stitch():
    d = out_dir("t23")
    b = boards.via_grid(d, pitch=3.0, void=(15, 6, 27, 14))
    r = evaluate(b, "--checks", "drc", "--no-drc")
    need(r["rep"], f"no report: {r['err'][-300:]}")
    sm = [x for x in r["rep"]["checks"]["drc"] if x["rule"] == "stitching_map"]
    need(sm and "_map" in sm[0]["metrics"], f"no stitching map: {sm[0]['summary'] if sm else 'rule missing'}")
    m = dict(sm[0]["metrics"]["_map"], void_area_mm2=sm[0]["metrics"]["_map"]["void"]["area_mm2"])
    need(abs(m["void_area_mm2"] - 12 * 8) / 96 < 0.35, f"largest via-free area {m.get('void_area_mm2')} mm2, want ~96")
    need(m.get("p50_mm") and 0.25 * 3 < m["p50_mm"] < 0.75 * 3, f"median distance {m.get('p50_mm')} mm for a 3 mm grid")
    return ok(f"median {m['p50_mm']:.2f} mm, void {m['void_area_mm2']:.0f} mm2")


@case(24, "Breakout necks inside a BGA field are measured and judged")
def t24_breakout():
    d = out_dir("t24")
    from gen import Board
    p = os.path.join(d, "bga", "bga.kicad_pcb")
    os.makedirs(os.path.dirname(p), exist_ok=True)
    b = Board(p, layers=4, size=(40, 30))
    b.plane("In1.Cu", "GND"); b.plane("In2.Cu", "GND")
    pads = []
    for i in range(16):
        for j in range(16):
            pads.append({"num": f"{chr(65 + i)}{j + 1}", "net": "GND" if (i + j) % 3 == 0 else None, "dx": (j - 7.5) * 0.8,
                         "dy": (i - 7.5) * 0.8, "sx": 0.4, "shape": "circle"})
    pads[0]["net"] = "LANE0_TX_P"; pads[1]["net"] = "LANE0_TX_N"
    b.footprint("U1", "BGA256", 15, 15, pads, courtyard=(-7, -7, 7, 7))
    x0, y0 = 15 - 6.0, 15 - 6.0
    # 0.10 mm necks from the balls to the field edge (x = 8), then the 0.12 mm pair outside
    b.track("F.Cu", [(x0, y0), (8.0, y0)], 0.10, "LANE0_TX_P")
    b.track("F.Cu", [(x0 + 0.8, y0), (x0 + 0.8, y0 - 0.55), (8.0, y0 - 0.55)], 0.10, "LANE0_TX_N")
    b.track("F.Cu", [(8.0, y0), (8.0, y0 - 0.135), (2.0, y0 - 0.135)], 0.12, "LANE0_TX_P")
    b.track("F.Cu", [(8.0, y0 - 0.55), (8.0, y0 - 0.405), (2.0, y0 - 0.405)], 0.12, "LANE0_TX_N")
    b.save()
    r = evaluate(p, "-c", _cfg_class(d, "bga", extra='[areas.bga]\nfootprints = ["U1"]\nmargin = 0.0\n'), "--checks", "drc", "--no-drc")
    need(r["rep"], f"no report: {r['err'][-300:]}")
    nk = [f for f in ev_findings(r["ev"]) if f["check"] == "breakout"]
    need(nk, "no breakout findings for 0.10 mm necks in the BGA field")
    return ok(f"{len(nk)} breakout findings")


@case(25, "Timing: delay = length x sqrt(eps_eff) / c within 2 %; skew per matching group")
def t25_timing():
    d = out_dir("t25")
    b = boards.pair_over_plane(d, "tm", length=50.0)
    cfg = _cfg_class(d, "tm", extra='[timing.groups]\nG0 = ["LANE0_TX_P", "LANE0_TX_N"]\n')
    r = evaluate(b, "-c", cfg, "--checks", "timing,impedance", "--no-drc")
    need(r["rep"] and "timing" in r["rep"]["checks"], f"no timing check: {(r['err'] or '')[-200:]}")
    t = r["rep"]["checks"]["timing"]
    row = next(x for x in t["nets"] if x["net"] == "LANE0_TX_P")
    want = row["length_mm"] * 1e-3 * math.sqrt(row["eps_eff"]) / 2.998e8 * 1e12
    need(abs(row["delay_ps"] - want) / want < 0.02, f"delay {row['delay_ps']} ps vs {want:.1f}")
    need(t["groups"]["G0"]["skew_ps"] < 1.0, "matched pair shows skew")
    return ok(f"{row['delay_ps']:.1f} ps over {row['length_mm']:.1f} mm")


@case(26, "Via technologies: a blind via has no stub; a back-drilled via keeps only its remaining stub")
def t26_vias():
    d = out_dir("t26")
    from gen import Board
    p = os.path.join(d, "bv", "bv.kicad_pcb")
    os.makedirs(os.path.dirname(p), exist_ok=True)
    b = Board(p, layers=4, size=(40, 20))
    b.plane("In2.Cu", "GND")
    b.track("F.Cu", [(5, 6), (20, 6)], 0.12, "LANE0_TX_P")
    b.via(20, 6, "LANE0_TX_P", kind="blind", top="F.Cu", bottom="In1.Cu", drill=0.15, size=0.35)
    b.track("In1.Cu", [(20, 6), (35, 6)], 0.12, "LANE0_TX_P")
    b.track("F.Cu", [(5, 14), (20, 14)], 0.12, "LANE1_TX_P")
    b.via(20, 14, "LANE1_TX_P")
    b.track("In1.Cu", [(20, 14), (35, 14)], 0.12, "LANE1_TX_P")
    b.track("F.Cu", [(5, 18), (20, 18)], 0.12, "LANE2_TX_P")              # through via, back-drilled up to In2
    b.via(20, 18, "LANE2_TX_P", backdrill_bottom="In2.Cu")
    b.track("In1.Cu", [(20, 18), (35, 18)], 0.12, "LANE2_TX_P")
    b.save()
    r = evaluate(p, "-c", _cfg_class(d, "bv"), "--checks", "drc", "--no-drc")
    need(r["rep"], f"no report: {r['err'][-300:]}")
    st = {f["net"]: f for f in findings(r["rep"], "via_stub")}
    m = {}
    for rr in r["rep"]["checks"]["drc"]:
        if rr["rule"] == "via_stub":
            m = rr["metrics"]
    need("LANE0_TX_P" not in st, "a blind via F.Cu-In1 is reported with a stub")
    need("LANE1_TX_P" in st or any("LANE1" in k for k in m), "the through via F.Cu-In1 has no stub result")
    s1 = m.get("LANE1_TX_P", {}).get("max_stub_mm")
    s2 = m.get("LANE2_TX_P", {}).get("max_stub_mm")
    need(s1 and s2 is not None and 0 < s2 < s1, f"back-drilled stub {s2} mm vs full through via {s1} mm")
    return ok(f"blind via: no stub; through via stub {s1} mm; back-drilled to In2: {s2} mm")


@case(27, "Planeless boards and split grounds: return path judged without planes; ground network solved")
def t27_ground():
    d = out_dir("t27")
    b = boards.planeless(d)
    r = evaluate(b, "--checks", "drc", "--no-drc")
    need(r["rep"], f"no report: {r['err'][-300:]}")
    rr = [x for x in r["rep"]["checks"]["drc"] if x["rule"] == "return_proximity"]
    need(rr and rr[0]["findings"], "no return_proximity rule for a board without planes")
    from gen import Board
    # (b) split grounds: AGND under the left half of In1, DGND under the right; the pair crosses the split
    p = os.path.join(d, "split", "split.kicad_pcb")
    os.makedirs(os.path.dirname(p), exist_ok=True)
    bb = Board(p, layers=4, size=(40, 20))
    bb.zone("In1.Cu", [(0.5, 0.5), (19.75, 0.5), (19.75, 19.5), (0.5, 19.5)], "AGND")
    bb.zone("In1.Cu", [(20.25, 0.5), (39.5, 0.5), (39.5, 19.5), (20.25, 19.5)], "DGND")
    bb.plane("In2.Cu", "DGND")
    bb.pair("F.Cu", "LANE0_TX", 5, 35, 10)
    bb.save()
    r = evaluate(p, "--checks", "drc", "--no-drc")
    sc = [f for f in findings(r["rep"], "split_crossing") if f["severity"] == "error"]
    need(sc, "a pair crossing from AGND to DGND gives no split crossing")
    # (c) ground network: 1 A returning through one 0.3 mm x 20 mm F.Cu ground track: rise = I R
    p = os.path.join(d, "gret", "gret.kicad_pcb")
    os.makedirs(os.path.dirname(p), exist_ok=True)
    bb = Board(p, layers=2, size=(40, 20))
    bb.footprint("U1", "TPS62133", 10, 10, [{"num": 1, "net": "+1V8", "dx": 0, "dy": 2, "sx": 0.6},
                                            {"num": 2, "net": "GND", "dx": 0, "dy": 0, "sx": 0.6}])
    bb.footprint("U2", "LOAD", 30, 10, [{"num": 1, "net": "+1V8", "dx": 0, "dy": 2, "sx": 0.6},
                                        {"num": 2, "net": "GND", "dx": 0, "dy": 0, "sx": 0.6}])
    bb.track("F.Cu", [(10, 12), (30, 12)], 1.0, "+1V8")
    bb.track("F.Cu", [(10, 10), (30, 10)], 0.3, "GND")
    bb.save()
    cfg = toml(os.path.join(d, "g.toml"), '[rails.V18]\nnet = "+1V8"\nv = 1.8\ntol_pct = 5.0\ncurrent_a = 1.0\n'
               'source = "U1"\nregulator = "U1"\nloads = ["U2"]\n[engines]\ndc = "openpcbsi"\n')
    r = evaluate(p, "-c", cfg, "--checks", "dc", "--no-drc")
    need(r["rep"], f"no report: {r['err'][-300:]}")
    row = r["rep"]["checks"]["dc"]["rails"][0]
    want = 1.724e-8 * 0.020 / (0.3e-3 * 35e-6) * 1000
    need("ground_rise_mv" in row, f"ground not solved: {row.get('ground')}")
    need(abs(row["ground_rise_mv"] - want) / want < 0.1, f"ground rise {row['ground_rise_mv']:.2f} mV vs I R {want:.2f} mV")
    return ok(f"planeless rule; AGND/DGND split crossing found; ground rise {row['ground_rise_mv']:.1f} mV (I R {want:.1f})")


@case(28, "Current and heat: IPC-2152 trace rise within 30 % of IPC-2221; via current judged")
def t28_current():
    d = out_dir("t28")
    b = boards.rail_board(d, "cur", current_field="4A", second_rail=False)
    r = evaluate(b, "--checks", "dc", "--no-drc")
    need(r["ev"], f"no evaluation: {r['err'][-300:]}")
    hot = ev_findings(r["ev"], "trace_heating")
    need(hot, "no trace_heating findings for 4 A in 0.3 mm tracks")
    f = hot[0]
    w_mil, t_mil = 0.3 / 0.0254, 0.035 / 0.0254
    dT = (4.0 / (0.048 * (w_mil * t_mil) ** 0.725)) ** (1 / 0.44)
    need(abs(f["value"] - dT) / dT < 0.3, f"rise {f['value']:.1f} C vs IPC-2221 {dT:.1f} C")
    need(ev_findings(r["ev"], "via_current"), "no via_current judgement")
    return ok(f"rise {f['value']:.0f} C (IPC-2221 {dT:.0f} C)")


# =================================================================================================
# E. Use during layout
# =================================================================================================
@case(29, "Layout stage: an unrouted board gives 'pending', not fails")
def t29_stage():
    d = out_dir("t29")
    b = boards.rail_board(d, "unrouted", unrouted=True)
    r = evaluate(b, "--checks", "drc,dc")
    need(r["ev"], f"no evaluation: {r['err'][-300:]}")
    need(r["ev"].get("stage") in ("placement", "routing"), f"stage {r['ev'].get('stage')}")
    fails = [f for f in ev_findings(r["ev"], sev=("fail",)) if "unconnected" in f["reason"] or "not connected" in f["reason"]]
    need(not fails, f"{len(fails)} fails for unrouted connections")
    need(ev_findings(r["ev"], sev=("pending",)), "no pending findings")
    return ok("unrouted connections are pending")


@case(30, "Speed: rule checks on a 100 x 100 mm 6-layer board in < 30 s; results written per check")
def t30_speed():
    d = out_dir("t30")
    b = boards.big_board(d)
    check(b, "drc")                                   # warm the snapshot/raster cache
    r = check(b, "drc")
    need(r["rep"], f"no report: {r['err'][-300:]}")
    need(r["seconds"] < 30, f"rule checks took {r['seconds']:.0f} s")
    need(os.path.exists(os.path.join(r["dir"], "progress.json")), "no progress.json written during the run")
    return ok(f"{r['seconds']:.1f} s")


@case(31, "Guidance: DC drop budget sums to the total; fix hints; width for a target impedance")
def t31_guidance():
    d = out_dir("t31")
    b = boards.rail_board(d, "g", current_field="1.5A", second_rail=False)
    r = evaluate(b, "--checks", "dc", "--no-drc")
    need(r["rep"], f"no report: {r['err'][-300:]}")
    row = r["rep"]["checks"]["dc"]["rails"][0]
    bud = row.get("budget")
    need(bud, "no drop budget per rail")
    need(abs(sum(x["mv"] for x in bud) - row["drop_mv"]) / row["drop_mv"] < 0.02, "budget does not sum to the drop")
    s = boards.pair_over_plane(d, "hint", slot=True)
    r = evaluate(s, "--checks", "drc", "--no-drc")
    need(any(f.get("hint") for f in ev_findings(r["ev"], "split_crossing")), "split crossing has no fix hint")
    rc, so, se, dt = run(["plan", s, "--layer", "F.Cu", "--zdiff", "100", "--gap", "0.15"])
    m = re.search(r"w\s*=\s*([\d.]+).*Zdiff\s*=\s*([\d.]+)", so)
    need(rc == 0 and m and abs(float(m.group(2)) - 100) < 1.0, f"plan: {(so or se)[-200:]}")
    return ok(f"budget sums; hint present; plan w={m.group(1)} mm")


@case(32, "Library API and live board: a snapshot, and the board open in KiCad, give the CLI's findings")
def t32_api():
    d = out_dir("t32")
    b = boards.pair_over_plane(d, "api", slot=True)
    r = evaluate(b, "--checks", "drc", "--no-drc")
    sys.path.insert(0, ROOT)
    try:
        from openpcbsi import api
    except ImportError:
        raise Fail("no openpcbsi.api module")
    from openpcbsi import snapshot, run as R
    snap = snapshot.load_or_extract(b, R.CACHE)
    res = api.evaluate(snap, None, families=["rp"])
    ids_api = {f["id"] for f in res["findings"] if f["severity"] in ("fail", "warn")}
    ids_cli = {f["id"] for f in ev_findings(r["ev"], sev=("fail", "warn")) if f["family"] == "RP"}
    need(ids_api == ids_cli, f"API {len(ids_api)} vs CLI {len(ids_cli)} findings")
    # the same board open in a (private, headless) KiCad, read through its IPC API
    from live_kicad import LiveKiCad
    with LiveKiCad(b, os.path.join(d, "kicad")) as k:
        rl = run(["evaluate", "--live", "--socket", k.socket, "--checks", "drc", "--no-drc", "--no-markers",
                  "-o", os.path.join(d, "live_run")], timeout=600)
        ls = api.live_snapshot(k.socket)
    ev_live = _load(os.path.join(d, "live_run", "evaluation.json"))
    need(ev_live, f"evaluate --live failed (rc {rl[0]}): {(rl[2] or rl[1])[-300:]}")
    ids_live = {f["id"] for f in ev_live["findings"] if f["severity"] in ("fail", "warn")}
    ids_file = {f["id"] for f in r["ev"]["findings"] if f["severity"] in ("fail", "warn")}
    need(ids_live == ids_file, f"live board {len(ids_live)} vs saved file {len(ids_file)} findings")
    need(len(ls["tracks"]) == len(snap["tracks"]) and len(ls["vias"]) == len(snap["vias"]), "live snapshot differs from the file's")
    return ok(f"API = CLI ({len(ids_cli)} RP findings); live board through KiCad's IPC API = saved file ({len(ids_live)} findings)")


@case(33, "Concurrent runs: a job lock and a shared thread budget; parallel runs agree")
def t33_conc():
    d = out_dir("t33")
    b = boards.pair_over_plane(d, "conc", slot=True)
    procs = []
    for i in range(3):
        o = os.path.join(d, f"run{i}")
        procs.append((o, subprocess.Popen([sys.executable, "-m", "openpcbsi", "check", b, "--checks", "drc", "-o", o, "--no-markers"],
                                          cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)))
    outs = []
    for o, p in procs:
        so, se = p.communicate(timeout=900)
        need(p.returncode == 0, f"parallel run failed: {se[-200:]}")
        outs.append(json.dumps([f["id"] if "id" in f else (f["rule"], f["net"], f["x"]) for f in findings(_load(os.path.join(o, "report.json")))]))
    need(len(set(outs)) == 1, "parallel runs disagree")
    need(os.path.exists(os.path.join(ROOT, "cache", "jobs.lock")) or any("lock" in x for x in os.listdir(os.path.join(ROOT, "cache"))),
         "no job lock in the cache")
    return ok("3 parallel runs agree, lock present")


# =================================================================================================
# F. Sign-off
# =================================================================================================
@case(34, "Sign-off record: hashes and versions; a changed board invalidates it")
def t34_signoff():
    d = out_dir("t34")
    b = boards.pair_over_plane(d, "so")
    o = os.path.join(d, "so_run")
    rc, so, se, dt = run(["signoff", b, "--checks", "drc", "--no-drc", "-o", o])
    rec = os.path.join(o, "signoff.json")
    need(os.path.exists(rec), f"no signoff.json (rc {rc}): {(se or so)[-300:]}")
    s = _load(rec)
    for k in ("board_sha256", "project_files", "openpcbsi_version", "kicad_file_version", "checks", "verdict"):
        need(k in s, f"signoff.json lacks {k}")
    rc, so, se, dt = run(["verify-signoff", rec])
    need(rc == 0, f"verify-signoff fails on an untouched board: {so + se}")
    open(b, "a").write("\n")
    rc, so, se, dt = run(["verify-signoff", rec])
    need(rc != 0, "verify-signoff passes after the board changed")
    return ok("record written, verified, invalidated by a change")


@case(35, "Regression since the last sign-off: new / resolved findings by stable ID; waivers carried")
def t35_since():
    d = out_dir("t35")
    a = boards.pair_over_plane(d, "rev", slot=True)
    ra = evaluate(a, "--checks", "drc", "--no-drc", out=os.path.join(d, "runA"))
    need(ra["ev"], "no first evaluation")
    boards.pair_over_plane(d, "rev", slot=False)
    rb = evaluate(a, "--checks", "drc", "--no-drc", "--since", os.path.join(d, "runA", "evaluation.json"), out=os.path.join(d, "runB"))
    need(rb["ev"] and rb["ev"].get("since"), f"no 'since' section: {(rb['err'] or '')[-200:]}")
    s = rb["ev"]["since"]
    need(not s.get("new"), f"{len(s.get('new', []))} new findings after removing the slot")
    prev = {f["id"]: f for f in ra["ev"]["findings"]}
    res = [prev[x["id"]] for x in s.get("resolved", [])]
    errs = [f for f in res if f["severity"] == "fail"]
    need(len(errs) == 4, f"{len(errs)} resolved fails, want the 4 slot errors")
    away = [f for f in res if f["x"] is not None and abs(f["x"] - 20.0) > 1.0]
    need(not away, f"{len(away)} resolved findings away from the slot")
    return ok(f"{len(res)} resolved, all at the slot (4 fails), 0 new")


@case(36, "Packaging: pip install gives an `openpcbsi` command; cache location configurable")
def t36_pkg():
    d = out_dir("t36")
    need(os.path.exists(os.path.join(ROOT, "pyproject.toml")), "no pyproject.toml")
    venv = os.path.join(d, "venv")
    subprocess.run([sys.executable, "-m", "venv", "--system-site-packages", venv], check=True)
    p = subprocess.run([os.path.join(venv, "bin", "pip"), "install", "--no-deps", "--no-build-isolation", "-q", ROOT],
                       capture_output=True, text=True, env=dict(os.environ, PIP_CACHE_DIR=os.path.join(d, "pipcache")))
    need(p.returncode == 0, f"pip install failed: {p.stderr[-300:]}")
    exe = os.path.join(venv, "bin", "openpcbsi")
    need(os.path.exists(exe), "no openpcbsi entry point")
    cache = os.path.join(d, "mycache")
    b = boards.pair_over_plane(d, "pkg")
    p = subprocess.run([exe, "check", b, "--checks", "drc", "-o", os.path.join(d, "o"), "--no-markers"], capture_output=True,
                       text=True, env=dict(os.environ, OPENPCBSI_CACHE=cache), cwd=d)
    need(p.returncode == 0, f"installed command failed: {p.stderr[-300:]}")
    need(glob.glob(os.path.join(cache, "snapshot_*")), "OPENPCBSI_CACHE not used")
    return blocked("Linux install works; Windows and macOS need CI runners")


# =================================================================================================
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--items")
    ap.add_argument("--slow", action="store_true")
    ap.add_argument("--list", action="store_true")
    a = ap.parse_args()
    sel = [c for c in CASES if (a.slow or not c["slow"])]
    if a.items:
        want = {int(x) for x in a.items.split(",")}
        sel = [c for c in CASES if c["item"] in want]
    if a.list:
        for c in CASES:
            print(f"#{c['item']:2d} {'(slow) ' if c['slow'] else ''}{c['title']}")
        return 0
    os.makedirs(OUT, exist_ok=True)
    resf = os.path.join(OUT, "acceptance.json")
    results = _load(resf) or {}
    for c in sel:
        t = time.time()
        try:
            r = c["fn"]()
        except Fail as e:
            r = {"status": "fail", "note": str(e)}
        except Exception as e:
            r = {"status": "error", "note": f"{type(e).__name__}: {e}", "trace": traceback.format_exc()[-1500:]}
        r.update(item=c["item"], title=c["title"], seconds=round(time.time() - t, 1), date=time.strftime("%Y-%m-%d %H:%M"))
        print(f"#{c['item']:2d} {r['status'].upper():8s} {r['seconds']:7.1f}s  {c['title']}\n            {r['note'][:300]}", flush=True)
        results = _record(resf, str(c["item"]), r)
    write_md(results)
    return 1 if any(results[str(c["item"])]["status"] in ("fail", "error") for c in sel) else 0


def _record(resf, key, r):
    """Merge one item's result into the results file under a lock: two runs (fast and slow items) may
    write at the same time, and each must keep the other's items."""
    import fcntl
    with open(resf + ".lock", "w") as lk:
        fcntl.flock(lk, fcntl.LOCK_EX)
        results = _load(resf) or {}
        results[key] = r
        tmp = f"{resf}.{os.getpid()}.tmp"
        json.dump(results, open(tmp, "w"), indent=1)
        os.replace(tmp, resf)
    return results


def write_md(results):
    L = ["# Acceptance status", "", "Generated by `python3 tests/acceptance.py`. Each item's test builds boards with a known",
         "answer and compares openpcbsi's result with an oracle (planted defect, closed form, independent solver or",
         "physics invariant). `blocked` means the test needs an input this machine does not have.", "",
         "| # | Test | Status | Note |", "|---|------|--------|------|"]
    for c in CASES:
        r = results.get(str(c["item"]))
        st, note = (r["status"], r["note"].replace("|", "/")[:160]) if r else ("not run", "")
        L.append(f"| {c['item']} | {c['title']} | {st} | {note} |")
    open(os.path.join(ROOT, "docs", "acceptance.md"), "w").write("\n".join(L) + "\n")


if __name__ == "__main__":
    sys.exit(main())
