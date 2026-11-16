"""Correlation with characterised boards (#8): fabricated KiCad boards whose measurements are published.

    python3 -m openpcbsi correlate [--fetch] [--reference antmicro_si]

The boards and measurements are listed in deps/references.lock with their repositories, pinned commits
and the SHA-256 of every file. `--fetch` downloads them into the cache; every run checks the hashes.

For each case openpcbsi models the measured structure from the KiCad board alone, as it would any
board: the net's path from the measured connector, with every branch as a stub (open where the part
at its end is not fitted), lines from the 2D solver of the board's stackup, the far connector
terminated in 50 ohm as the VNA's second port was. It then compares the model's S11 with the
measurement.

The gate is set on the features a fixture and its cables do not move: the frequency of the stub's
quarter-wave reflection peak and its half-wave null, within FREQ_TOL. The |S11| curve error over the
band is reported but not judged, because the measurement includes the SMA launches and cables, which
the board model does not.
"""
import hashlib
import json
import os
import subprocess
import tomllib

import numpy as np

FREQ_TOL = 0.05          # 5 % on resonance frequencies (~10 % on eps_eff): set before the first comparison
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LOCK = os.environ.get("OPENPCBSI_REFERENCES", os.path.join(HERE, "deps", "references.lock"))


def _sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        h.update(f.read())
    return h.hexdigest()


def load_lock():
    return tomllib.load(open(LOCK, "rb"))


def ref_dir(cache_dir, name, ref):
    return os.path.join(cache_dir, "references", f"{name}_{ref['board_commit'][:10]}_{ref['data_commit'][:10]}")


def fetch(cache_dir, name, ref, log=print):
    """Download the reference files at the pinned commits (raw GitHub URLs) and check their hashes."""
    d = ref_dir(cache_dir, name, ref)
    os.makedirs(d, exist_ok=True)
    for local, spec in ref["files"].items():
        dst = os.path.join(d, local)
        if os.path.exists(dst) and _sha(dst) == spec["sha256"]:
            continue
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        repo = ref["board_repo"] if spec["repo"] == "board" else ref["data_repo"]
        commit = ref["board_commit"] if spec["repo"] == "board" else ref["data_commit"]
        url = repo.replace("https://github.com/", "https://raw.githubusercontent.com/") + f"/{commit}/{spec['path']}"
        log(f"  fetching {url}")
        subprocess.run(["curl", "-sfL", "-o", dst, url], check=True)
    return verify(cache_dir, name, ref)


def verify(cache_dir, name, ref):
    d = ref_dir(cache_dir, name, ref)
    bad = []
    for local, spec in ref["files"].items():
        p = os.path.join(d, local)
        if not os.path.exists(p):
            bad.append(f"{local}: missing (run with --fetch)")
        elif _sha(p) != spec["sha256"]:
            bad.append(f"{local}: SHA-256 differs from deps/references.lock")
    return d, bad


def read_vna(path):
    a = np.loadtxt(path, delimiter=",", skiprows=1)
    return a[:, 0] * 1e6, a[:, 1] + 1j * a[:, 2]


SMOOTH_HZ = 300e6        # one period of the measurement's cable ripple (~260 MHz): averages it out


def _smooth(db, n=9, f=None):
    if f is not None and len(f) > 2:
        n = max(3, int(round(SMOOTH_HZ / (f[1] - f[0]))) | 1)
    k = np.ones(n) / n
    return np.convolve(np.pad(db, n // 2, mode="edge"), k, mode="valid")


def features(f, s11, fmin=0.5e9):
    """First quarter-wave reflection peak and the half-wave null after it (smoothed |S11|: the
    measurement's fine ripple comes from its cables)."""
    db = _smooth(20 * np.log10(np.abs(s11) + 1e-12), f=f)
    m = f >= fmin
    fi, di = f[m], db[m]
    ip = int(np.argmax(di))
    out = {"peak_hz": float(fi[ip]), "peak_db": round(float(di[ip]), 2)}
    after = (fi > fi[ip] * 1.3) & (fi < fi[ip] * 2.6)
    if after.any():
        j = np.flatnonzero(after)[int(np.argmin(di[after]))]
        if fi[j] < f[-1] * 0.98:
            out.update(null_hz=float(fi[j]), null_db=round(float(di[j]), 2))
    return out


def model_s11(ctx, case, f):
    from . import si
    b = ctx.b
    ctx.cfg.data["si"]["drivers"] = [case["driver"]]
    ctx.cfg.data["si"]["load_pf"] = 0.0                # an open stub end: its pad only
    paths, err = si.net_paths(ctx, case["net"])
    if not paths:
        raise RuntimeError(f"{case['net']}: {err}")
    path = next((p for p in paths if p["to"].split(".")[0] == case["far"]), None)
    if path is None:
        raise RuntimeError(f"{case['net']}: no path from {case['driver']} to {case['far']}")
    lib = si.LineLib(ctx)
    M, _, _ = si.net_abcd(ctx, lib, path, f, None, None)
    S = si.abcd_to_s(M, 50.0)
    stubs = [{"to": s[2]["ref"], "length_mm": round(sum(e[3] for e in s[1] if e[0] == "line"), 3),
              "end_fitted": not b.fps[s[2]["ref"]]["dnp"]} for s in path.get("stubs", [])]
    return S["S11"], {"path_mm": path["length_mm"], "stubs": stubs}


def correlate(cache_dir, runs_dir, name="antmicro_si", log=print, do_fetch=False):
    from . import run as runmod
    lock = load_lock()
    ref = lock[name]
    d, bad = fetch(cache_dir, name, ref, log) if do_fetch else verify(cache_dir, name, ref)
    if bad:
        raise RuntimeError("reference data: " + "; ".join(bad))
    board = os.path.join(d, "si-simulation-test-board.kicad_pcb")
    b, cfg, ctx = runmod.load(board)
    rows = []
    for case in ref["cases"]:
        fm, sm = read_vna(os.path.join(d, case["measured"]))
        s_mod, info = model_s11(ctx, case, fm)
        fe_m, fe_x = features(fm, sm), features(fm, s_mod)
        row = {"case": case["name"], "net": case["net"], "measured_from": case["driver"], **info,
               "measured": fe_m, "model": fe_x}
        checks = []
        for k in ("peak_hz", "null_hz"):
            if k in fe_m and k in fe_x:
                e = (fe_x[k] - fe_m[k]) / fe_m[k]
                checks.append({"feature": k.replace("_hz", ""), "measured_ghz": round(fe_m[k] / 1e9, 4),
                               "model_ghz": round(fe_x[k] / 1e9, 4), "error_pct": round(100 * e, 2),
                               "pass": abs(e) <= FREQ_TOL})
        dm, dx = 20 * np.log10(np.abs(sm) + 1e-12), 20 * np.log10(np.abs(s_mod) + 1e-12)
        band = fm <= 4e9
        row["s11_rms_db_to_4ghz"] = round(float(np.sqrt(np.mean((_smooth(dm, f=fm) - _smooth(dx, f=fm))[band] ** 2))), 2)
        row["checks"] = checks
        row["pass"] = bool(checks) and all(c["pass"] for c in checks)
        rows.append(row)
        log(f"  {case['name']:12s} " + "  ".join(f"{c['feature']} {c['model_ghz']:.3f} vs {c['measured_ghz']:.3f} GHz ({c['error_pct']:+.1f} %)"
                                                  for c in checks) + f"  |S11| rms {row['s11_rms_db_to_4ghz']} dB  -> {'PASS' if row['pass'] else 'FAIL'}")
    out = {"reference": name, "title": ref["title"], "license": ref["license"], "board_commit": ref["board_commit"],
           "data_commit": ref["data_commit"], "tolerance_pct": 100 * FREQ_TOL, "cases": rows,
           "pass": all(r["pass"] for r in rows)}
    od = os.path.join(runs_dir, "correlation")
    os.makedirs(od, exist_ok=True)
    json.dump(out, open(os.path.join(od, f"{name}.json"), "w"), indent=1)
    try:
        _plot(d, ref, ctx, od, name)
    except Exception as e:
        log(f"  (no plot: {e})")
    return out


def _plot(d, ref, ctx, od, name):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, len(ref["cases"]), figsize=(6 * len(ref["cases"]), 3.6))
    for ax, case in zip(np.atleast_1d(axs), ref["cases"]):
        fm, sm = read_vna(os.path.join(d, case["measured"]))
        s_mod, _ = model_s11(ctx, case, fm)
        ax.plot(fm / 1e9, 20 * np.log10(np.abs(sm)), lw=1, label="measured (Antmicro VNA)")
        ax.plot(fm / 1e9, 20 * np.log10(np.abs(s_mod)), lw=1.4, label="openpcbsi model")
        ax.set_title(case["name"])
        ax.set_xlabel("GHz")
        ax.set_ylabel("|S11| dB")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(od, f"{name}.png"), dpi=110)
