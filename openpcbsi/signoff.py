"""Sign-off record (#34): one evaluation tied to the exact inputs it judged.

    python3 -m openpcbsi signoff BOARD [-c CONFIG] [--waivers FILE] [--reviewer NAME] [evaluate options] [-o DIR]
    python3 -m openpcbsi verify-signoff DIR/signoff.json

signoff.json (schema openpcbsi-signoff/1) holds the SHA-256 of the board, of each KiCad project file
next to it, of the config, the waivers and the evaluation; the KiCad and openpcbsi versions; every
check's status; and the verdict. verify-signoff recomputes every hash and fails if any input changed,
so a record cannot be carried over to a board it did not judge.
"""
import datetime
import hashlib
import json
import os
import subprocess

from . import __version__

SCHEMA = "openpcbsi-signoff/1"
VERDICT = {0: "release", 1: "release with warnings", 2: "blocked: failures", 3: "blocked: could not evaluate"}


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


def write(ev, out_dir, cfg_path=None, waivers=None, reviewer=None):
    board = ev["board"]
    stem = os.path.splitext(board)[0]
    proj = {p: sha(p) for p in (stem + ext for ext in (".kicad_pro", ".kicad_dru", ".kicad_sch")) if os.path.exists(p)}
    evp = os.path.join(out_dir, "evaluation.json")
    blockers = [{"id": f["id"], "check": f["check"], "reason": f["reason"]} for f in ev["findings"]
                if f["severity"] == "fail" and not f.get("waived") and not f.get("area")]
    rec = {
        "schema": SCHEMA, "date": datetime.datetime.now().isoformat(timespec="seconds"), "reviewer": reviewer,
        "board": board, "board_sha256": sha(board), "kicad_file_version": ev.get("kicad_file_version"),
        "kicad": ev.get("kicad"), "project_files": proj,
        "config": {"path": os.path.abspath(cfg_path), "sha256": sha(cfg_path)} if cfg_path else {"path": None, "detected": True},
        "waivers": {"path": os.path.abspath(waivers), "sha256": sha(waivers)} if waivers else None,
        "openpcbsi_version": tool_version(), "solver_versions": ev.get("solver_versions"), "kipi": ev.get("kipi"),
        "evaluation": {"path": evp, "sha256": sha(evp)}, "stage": ev.get("stage"),
        "checks": [{k: c.get(k) for k in ("family", "check", "requirement", "status", "fail", "warn", "pending")} for c in ev["checks"]],
        "families": ev["families"], "counts": ev["counts"], "exit_code": ev["exit_code"],
        "verdict": VERDICT.get(ev["exit_code"], "blocked"),
        "not_judged": sorted(c["check"] for c in ev["checks"] if c["status"] in ("not built", "not run", "nothing to check", "not judged", "unvalidated", "assumed inputs")),
        "blockers": blockers[:200],
    }
    if ev.get("stage") != "final":
        rec["verdict"] = f"blocked: layout stage '{ev.get('stage')}'"
    p = os.path.join(out_dir, "signoff.json")
    json.dump(rec, open(p, "w"), indent=1)
    return p, rec


def verify(path):
    rec = json.load(open(path))
    bad = []
    def chk(label, p, want):
        if p is None:
            return
        if not os.path.exists(p):
            bad.append(f"{label} {p} is missing")
        elif sha(p) != want:
            bad.append(f"{label} {p} changed since the sign-off")
    chk("board", rec["board"], rec["board_sha256"])
    for p, h in rec.get("project_files", {}).items():
        chk("project file", p, h)
    stem = os.path.splitext(rec["board"])[0]
    for ext in (".kicad_pro", ".kicad_dru", ".kicad_sch"):
        if os.path.exists(stem + ext) and (stem + ext) not in rec.get("project_files", {}):
            bad.append(f"project file {stem + ext} appeared after the sign-off")
    if rec.get("config", {}).get("path"):
        chk("config", rec["config"]["path"], rec["config"]["sha256"])
    if rec.get("waivers"):
        chk("waivers", rec["waivers"]["path"], rec["waivers"]["sha256"])
    chk("evaluation", rec["evaluation"]["path"], rec["evaluation"]["sha256"])
    return not bad, bad, rec
