"""Sign-off record (#34): one evaluation tied to the exact inputs it judged.

    python3 -m openpcbsi signoff BOARD [-c CONFIG] [--waivers FILE] [--reviewer NAME] [--sign-key KEY] [evaluate options] [-o DIR]
    python3 -m openpcbsi verify-signoff DIR/signoff.json [--allowed-signers FILE] [--rerun]

signoff.json (schema openpcbsi-signoff/2) holds the SHA-256 of the board, of each KiCad project file
next to it, of the config, the waivers and the evaluation; the KiCad and openpcbsi versions; every
check's status; and the verdict. Paths are stored relative to the record's folder, so a record that
is passed on does not show the author's folders.

verify-signoff checks, in turn:
- every input's hash: a record cannot be carried over to a board it did not judge;
- every field that follows from the hashed evaluation (exit code, counts, checks, families, blockers,
  verdict): a record whose verdict, blockers or counts were edited does not verify;
- with --allowed-signers, the record's SSH signature (signoff.json.sig, written with --sign-key):
  only that authenticates the record as a whole and the reviewer's name. A record and an evaluation
  edited together, consistently, are caught only by the signature, or by --rerun;
- with --rerun, the board is evaluated again with the recorded config and waivers and the verdict
  compared.
"""
import datetime
import hashlib
import json
import os
import subprocess
import sys
import tempfile

from . import __version__

SCHEMA = "openpcbsi-signoff/2"
VERDICT = {0: "release", 1: "release with warnings", 2: "blocked: failures", 3: "blocked: could not evaluate"}
NAMESPACE = "openpcbsi-signoff"
NOT_JUDGED = ("not built", "not run", "nothing to check", "not judged", "unvalidated", "assumed inputs")


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def tool_version():
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    try:
        g = subprocess.run(["git", "-C", root, "describe", "--always", "--dirty"], capture_output=True, text=True, timeout=10)
        if g.returncode == 0:
            return f"{__version__} ({g.stdout.strip()})"
    except (OSError, subprocess.TimeoutExpired):
        pass
    return __version__


def _rel(path, base):
    return os.path.relpath(os.path.abspath(path), base) if path else None


def _abs(path, base):
    return os.path.normpath(os.path.join(base, path)) if path else None


def derive(ev, warnings_exit=True):
    """Everything in a record that follows from the evaluation alone."""
    live = [f for f in ev["findings"] if not f.get("waived") and not f.get("area")]
    code = 2 if any(f["severity"] == "fail" for f in live) else \
        1 if (warnings_exit and any(f["severity"] == "warn" for f in live)) else 0
    blockers = [{"id": f["id"], "check": f["check"], "reason": f["reason"]} for f in live if f["severity"] == "fail"]
    counts = {s: sum(1 for f in live if f["severity"] == s) for s in ("fail", "warn", "info", "pass", "pending")}
    counts["in_area"] = sum(1 for f in ev["findings"] if f.get("area") and not f.get("waived"))
    counts["waived"] = sum(1 for f in ev["findings"] if f.get("waived"))
    verdict = VERDICT.get(code, "blocked")
    if ev.get("stage") != "final":
        verdict = f"blocked: layout stage '{ev.get('stage')}'"
    return {"exit_code": code, "counts": counts, "blockers": blockers[:200], "families": ev["families"], "stage": ev.get("stage"),
            "checks": [{k: c.get(k) for k in ("family", "check", "requirement", "status", "fail", "warn", "pending")} for c in ev["checks"]],
            "not_judged": sorted(c["check"] for c in ev["checks"] if c["status"] in NOT_JUDGED), "verdict": verdict}


def _cross_check(board):
    from . import run as runmod, snapshot as snap_mod
    snap = snap_mod.load_or_extract(board, runmod.CACHE)
    if snap.get("reader") != "kicad-monkey":
        return None
    x = snap_mod.cross_check(snap, runmod.CACHE)
    return {"reader": x["reader"], "compared": x["compared"], "differences": x["n_problems"],
            "not_cross_checked": x["not_cross_checked"], "first": (x["problems"] or [None])[0]}


def write(ev, out_dir, cfg_path=None, waivers=None, reviewer=None, sign_key=None, args=None):
    board = os.path.abspath(ev["board"])
    base = os.path.abspath(out_dir)
    stem = os.path.splitext(board)[0]
    proj = {_rel(p, base): sha(p) for p in (stem + ext for ext in (".kicad_pro", ".kicad_dru", ".kicad_sch")) if os.path.exists(p)}
    evp = os.path.join(out_dir, "evaluation.json")
    rec = {
        "schema": SCHEMA, "date": datetime.datetime.now().isoformat(timespec="seconds"), "reviewer": reviewer,
        "board": _rel(board, base), "board_sha256": sha(board), "kicad_file_version": ev.get("kicad_file_version"),
        "kicad": ev.get("kicad"), "project_files": proj,
        "config": {"path": _rel(cfg_path, base), "sha256": sha(cfg_path)} if cfg_path else {"path": None, "detected": True},
        "waivers": {"path": _rel(waivers, base), "sha256": sha(waivers)} if waivers else None,
        "openpcbsi_version": tool_version(), "solver_versions": ev.get("solver_versions"), "kipi": ev.get("kipi"),
        "evaluation": {"path": _rel(evp, base), "sha256": sha(evp)}, "evaluate_args": args or {},
    }
    warnings_exit = True
    if cfg_path:
        from .config import Config
        warnings_exit = Config.load(cfg_path).data.get("evaluate", {}).get("warnings_exit", True)
    rec.update(derive(ev, warnings_exit))
    # the board read again by an independent parser: a sign-off needs both readings to agree
    try:
        x = _cross_check(board)
        if x is not None:
            rec["cross_check"] = x
            if x["differences"]:
                rec["verdict"] = f"blocked: the two readings of the board differ ({x['differences']}): {x['first']}"
    except Exception as e:
        rec["cross_check"] = {"error": str(e)}
        rec["verdict"] = f"blocked: the cross-check could not run ({e})"
    p = os.path.join(out_dir, "signoff.json")
    json.dump(rec, open(p, "w"), indent=1)
    if sign_key:
        s = subprocess.run(["ssh-keygen", "-Y", "sign", "-f", sign_key, "-n", NAMESPACE, p], capture_output=True, text=True)
        if s.returncode != 0:
            raise RuntimeError(f"signing the record failed: {(s.stderr or s.stdout).strip()[-300:]}")
        rec["_signed"] = p + ".sig"
    return p, rec


def verify(path, allowed_signers=None, rerun=False):
    rec = json.load(open(path))
    base = os.path.dirname(os.path.abspath(path))
    bad, notes = [], []

    def chk(label, p, want):
        if p is None:
            return
        if not os.path.exists(p):
            bad.append(f"{label} {p} is missing")
        elif sha(p) != want:
            bad.append(f"{label} {p} changed since the sign-off")
    if rec.get("schema") != SCHEMA:
        bad.append(f"schema {rec.get('schema')}: this version verifies {SCHEMA} records (sign off again)")
        return False, bad, rec
    board = _abs(rec["board"], base)
    chk("board", board, rec["board_sha256"])
    for p, h in rec.get("project_files", {}).items():
        chk("project file", _abs(p, base), h)
    stem = os.path.splitext(board)[0]
    for ext in (".kicad_pro", ".kicad_dru", ".kicad_sch"):
        if os.path.exists(stem + ext) and _rel(stem + ext, base) not in rec.get("project_files", {}):
            bad.append(f"project file {os.path.basename(stem + ext)} appeared after the sign-off")
    cfg = _abs((rec.get("config") or {}).get("path"), base)
    if cfg:
        chk("config", cfg, rec["config"]["sha256"])
    if rec.get("waivers"):
        chk("waivers", _abs(rec["waivers"]["path"], base), rec["waivers"]["sha256"])
    evp = _abs(rec["evaluation"]["path"], base)
    chk("evaluation", evp, rec["evaluation"]["sha256"])
    # the fields that follow from the evaluation: recomputed, never trusted
    if os.path.exists(evp) and not any("evaluation" in b for b in bad):
        ev = json.load(open(evp))
        warnings_exit = True
        if cfg and os.path.exists(cfg):
            from .config import Config
            warnings_exit = Config.load(cfg).data.get("evaluate", {}).get("warnings_exit", True)
        d = derive(ev, warnings_exit)
        if d["exit_code"] != ev.get("exit_code"):
            bad.append(f"the evaluation's exit code {ev.get('exit_code')} does not follow from its findings ({d['exit_code']})")
        x = rec.get("cross_check") or {}
        if not x.get("error") and os.path.exists(board):
            try:
                x2 = _cross_check(board)
                if x2 is not None and x2["differences"] != x.get("differences"):
                    bad.append(f"cross-check: the record says {x.get('differences')} differences, reading the board again gives "
                               f"{x2['differences']}")
                x = x2 or x
            except Exception as e:
                notes.append(f"the cross-check could not be repeated ({e})")
        if x.get("error"):
            d["verdict"] = f"blocked: the cross-check could not run ({x['error']})"
        elif x.get("differences"):
            d["verdict"] = f"blocked: the two readings of the board differ ({x['differences']}): {x.get('first')}"
        for k in ("exit_code", "counts", "blockers", "families", "stage", "checks", "not_judged", "verdict"):
            if rec.get(k) != d[k]:
                bad.append(f"{k}: the record says {json.dumps(rec.get(k))[:120]}, the evaluation gives {json.dumps(d[k])[:120]}")
    # authenticity: only a signature covers the record as a whole and the reviewer's name
    sig = path + ".sig"
    if os.path.exists(sig) and allowed_signers:
        ident = rec.get("reviewer") or ""
        with open(path, "rb") as f:
            s = subprocess.run(["ssh-keygen", "-Y", "verify", "-f", allowed_signers, "-I", ident, "-n", NAMESPACE, "-s", sig],
                               stdin=f, capture_output=True, text=True)
        if s.returncode != 0:
            bad.append(f"the signature does not verify for reviewer '{ident}': {(s.stderr or s.stdout).strip()[-200:]}")
        else:
            notes.append(f"signed by '{ident}' ({(s.stdout or '').strip()[-120:]})")
    elif os.path.exists(sig):
        notes.append("signed, but no --allowed-signers given: the signature was not checked")
    else:
        notes.append("unsigned: the reviewer's name and the record as a whole are not authenticated (sign with --sign-key)")
    if rerun and not bad:
        with tempfile.TemporaryDirectory() as td:
            a = rec.get("evaluate_args") or {}
            cmd = [sys.executable, "-m", "openpcbsi", "evaluate", board, "-o", td, "--no-markers"]
            if cfg:
                cmd += ["-c", cfg]
            if rec.get("waivers"):
                cmd += ["--waivers", _abs(rec["waivers"]["path"], base)]
            for k in ("families", "checks", "stage"):
                if a.get(k):
                    cmd += [f"--{k}", str(a[k])]
            if a.get("no_drc"):
                cmd.append("--no-drc")
            r = subprocess.run(cmd, capture_output=True, text=True)
            ev2p = os.path.join(td, "evaluation.json")
            if not os.path.exists(ev2p):
                bad.append(f"--rerun: the evaluation could not run (exit {r.returncode})")
            else:
                ev2 = json.load(open(ev2p))
                if ev2["exit_code"] != rec["exit_code"] or ev2["counts"].get("fail") != rec["counts"].get("fail"):
                    bad.append(f"--rerun: the board evaluates to exit {ev2['exit_code']} with {ev2['counts'].get('fail')} fails, "
                               f"the record says exit {rec['exit_code']} with {rec['counts'].get('fail')}")
                else:
                    notes.append(f"--rerun: the board evaluates again to exit {ev2['exit_code']}")
    rec["_notes"] = notes
    return not bad, bad, rec
