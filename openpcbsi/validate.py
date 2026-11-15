"""G17: validation cases. A solver-based check reports its results as unvalidated unless the cases
for its solver have passed (results kept in <cache>/validation.json, keyed by solver version)."""
import json
import math
import os
import time

from . import xsec

SOLVER_VERSIONS = {"xsec2d": "xsec-1", "pdn": "pdn-1", "dc": "dc-1", "thermal": "th-1", "si": "si-1"}


def _case(name, solver, value, ref, tol_pct, unit, source):
    err = 100.0 * (value - ref) / ref
    return {"case": name, "solver": solver, "value": round(value, 4), "reference": round(ref, 4),
            "error_pct": round(err, 3), "tol_pct": tol_pct, "unit": unit, "source": source,
            "pass": bool(abs(err) <= tol_pct)}


def cases_xsec2d(res_um=5.0):
    out = []
    er = 4.2
    # 1-2: zero-thickness microstrip against Hammerstad-Jensen (stated accuracy 0.01 % for Z01, 0.2 % for eeff)
    for w, h in ((0.12, 0.1164), (0.30, 0.15)):
        r = xsec.solve_lines(xsec.simple_section("microstrip", w, h, er), res_um, (1e9,))["freqs"][0]
        z, ee = xsec.hammerstad_jensen(w, h, er)
        out.append(_case(f"microstrip w={w} h={h} er={er}: Z0", "xsec2d", r["Z0"], z, 2.0, "ohm", "Hammerstad-Jensen 1980"))
        out.append(_case(f"microstrip w={w} h={h} er={er}: eeff", "xsec2d", r["eeff"], ee, 1.5, "", "Hammerstad-Jensen 1980"))
    # 3: zero-thickness centred stripline against Cohn's exact solution
    for w, b in ((0.12, 0.45), (0.25, 0.6)):
        r = xsec.solve_lines(xsec.simple_section("stripline", w, None, er, b=b), res_um, (1e9,))["freqs"][0]
        out.append(_case(f"stripline w={w} b={b}: Z0", "xsec2d", r["Z0"], xsec.cohn_stripline(w, b, er), 1.5,
                         "ohm", "Cohn 1954 (exact, t=0)"))
    # 4: edge-coupled stripline, odd and even modes, exact
    w, s_, b = 0.12, 0.15, 0.45
    r = xsec.solve_lines(xsec.simple_section("stripline", w, None, er, gap=s_, b=b), res_um, (1e9,))["freqs"][0]
    z0e, z0o = xsec.cohn_coupled_stripline(w, s_, b, er)
    out.append(_case(f"coupled stripline w={w} s={s_} b={b}: Zodd", "xsec2d", r["Zodd"], z0o, 1.5, "ohm", "Cohn 1955 (exact, t=0)"))
    out.append(_case(f"coupled stripline w={w} s={s_} b={b}: Zeven", "xsec2d", r["Zeven"], z0e, 1.5, "ohm", "Cohn 1955 (exact, t=0)"))
    # 5: dielectric loss of a homogeneous line: alpha_d = pi f sqrt(er) tand / c, exact
    f = 5e9
    r = xsec.solve_lines(xsec.simple_section("stripline", 0.2, None, er, b=0.5, tand=0.02), res_um, (f,))["freqs"][0]
    ad = math.pi * f * math.sqrt(er) * 0.02 / xsec.C0 * 20 / math.log(10) * 1e-3
    out.append(_case("stripline dielectric loss at 5 GHz, tand 0.02", "xsec2d", r["loss_d_db_mm"], ad, 2.0, "dB/mm",
                     "exact for a homogeneous TEM line"))
    # 6: conductor loss in the wide-strip limit: two parallel-plate halves, R = Rs / w, exact as w/b -> inf
    w, b, t = 4.0, 0.2, 0.01
    r = xsec.solve_lines(xsec.simple_section("stripline", w, None, er, t=t, b=b), res_um, (f,))["freqs"][0]
    ref = (xsec.skin_rs(f) / (w * 1e-3)) / (2 * r["Z0"]) * 20 / math.log(10) * 1e-3
    out.append(_case("wide stripline conductor loss at 5 GHz (w/b=20)", "xsec2d", r["loss_c_db_mm"], ref, 3.0,
                     "dB/mm", "parallel-plate limit R = Rs/w (exact as w/b grows; fringing adds ~1 %)"))
    # 7: narrow stripline conductor loss against Pozar's approximate formula (information, wide tolerance)
    w, b, t = 0.2, 0.5, 0.035
    r = xsec.solve_lines(xsec.simple_section("stripline", w, None, er, t=t, b=b), res_um, (f,))["freqs"][0]
    ac = xsec.stripline_alpha_c(w, b, t, er, f) * 20 / math.log(10) * 1e-3
    out.append(_case("stripline conductor loss at 5 GHz, w=0.2 t=35 um", "xsec2d", r["loss_c_db_mm"], ac, 15.0, "dB/mm",
                     "Pozar eq. 3.198 (approximate formula; informative)"))
    return out


def run(cache_dir, which=("xsec2d",), log=print):
    path = os.path.join(cache_dir, "validation.json")
    db = {}
    if os.path.exists(path):
        with open(path) as f:
            db = json.load(f)
    for solver in which:
        t = time.time()
        if solver == "xsec2d":
            cases = cases_xsec2d()
        else:
            continue
        ok = bool(all(c["pass"] for c in cases))
        db[solver] = {"version": SOLVER_VERSIONS[solver], "pass": ok, "cases": cases,
                      "date": time.strftime("%Y-%m-%d %H:%M"), "seconds": round(time.time() - t, 1)}
        for c in cases:
            log(f"  [{'PASS' if c['pass'] else 'FAIL'}] {c['case']:55s} {c['value']:>10.4f} vs {c['reference']:>10.4f}"
                f"  ({c['error_pct']:+.2f} %, tol {c['tol_pct']} %)  {c['source']}")
    _save(path, db)
    return db


def status(cache_dir, solver):
    """'validated' if the cases for this solver version passed, else 'unvalidated'."""
    path = os.path.join(cache_dir, "validation.json")
    if not os.path.exists(path):
        return "unvalidated"
    with open(path) as f:
        db = json.load(f)
    e = db.get(solver)
    if e and e.get("pass") and e.get("version") == SOLVER_VERSIONS.get(solver):
        return "validated"
    return "unvalidated"


def register(cache_dir, solver, cases):
    """Record externally computed cases (PDN, DC, thermal closed forms)."""
    path = os.path.join(cache_dir, "validation.json")
    db = {}
    if os.path.exists(path):
        with open(path) as f:
            db = json.load(f)
    db[solver] = {"version": SOLVER_VERSIONS[solver], "pass": all(c["pass"] for c in cases), "cases": cases,
                  "date": time.strftime("%Y-%m-%d %H:%M")}
    _save(path, db)
    return db[solver]


def _save(path, db):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path + ".tmp", "w") as f:
        json.dump(db, f, indent=1)
    os.replace(path + ".tmp", path)
