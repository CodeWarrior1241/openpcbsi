"""Reference corpus (plan phase 0): real boards recorded once from KiCad itself.

    python3 tests/harness/corpus.py build [NAME ...]     record the boards of tests/corpus.toml
    python3 tests/harness/corpus.py list                 what is recorded

The corpus lives in corpus/ at the repository root (OPENPCBSI_CORPUS overrides it). It is gitignored:
the boards belong to their authors. Recording needs KiCad's pcbnew (the golden snapshot) and
kicad-cli; it is done while both are installed, and later readers are compared with it.
"""
import gzip
import json
import os
import shutil
import subprocess
import sys
import time
import tomllib

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
CORPUS = os.environ.get("OPENPCBSI_CORPUS", os.path.join(ROOT, "corpus"))
MANIFEST = os.path.join(ROOT, "tests", "corpus.toml")


def manifest():
    env = dict(os.environ, OPENPCBSI_ROOT=ROOT, KICAD_DEMOS=os.environ.get("KICAD_DEMOS", "/usr/share/kicad/demos"))
    out = []
    for e in tomllib.load(open(MANIFEST, "rb"))["board"]:
        src = e["source"]
        for k, v in env.items():
            src = src.replace("${" + k + "}", v)
        out.append(dict(e, path=src if "${" not in src else None))
    return out


def board_dir(name):
    return os.path.join(CORPUS, name)


def recorded(name):
    """Paths of a recorded board, or None."""
    d = board_dir(name)
    meta = os.path.join(d, "meta.json")
    if not os.path.exists(meta):
        return None
    m = json.load(open(meta))
    return dict(m, dir=d, board=os.path.join(d, m["board_file"]), snapshot=os.path.join(d, "snapshot_pcbnew.json.gz"),
                ipc2581=os.path.join(d, "ipc2581.xml"), drc=os.path.join(d, "drc.json"), gerbers=os.path.join(d, "gerbers"),
                pos=os.path.join(d, "pos.csv"))


def load_snapshot(name):
    r = recorded(name)
    with gzip.open(r["snapshot"], "rt") as f:
        return json.load(f)


def _cli(args, log):
    p = subprocess.run(["kicad-cli"] + args, capture_output=True, text=True)
    if p.returncode != 0:
        log(f"    kicad-cli {' '.join(args[:3])}: exit {p.returncode}: {(p.stderr or p.stdout).strip().splitlines()[-1:]}")
    return p.returncode


def build(names=None, log=print):
    sys.path.insert(0, ROOT)
    from openpcbsi import snapshot
    done = []
    for e in manifest():
        if names and e["name"] not in names:
            continue
        if not e["path"] or not os.path.exists(e["path"]):
            log(f"  {e['name']}: source not found ({e['source']}), skipped")
            continue
        d = board_dir(e["name"])
        if os.path.isdir(d):
            shutil.rmtree(d)
        os.makedirs(d)
        t = time.time()
        stem = os.path.splitext(os.path.basename(e["path"]))[0]
        for ext in (".kicad_pcb", ".kicad_pro", ".kicad_dru"):
            s = os.path.splitext(e["path"])[0] + ext
            if os.path.exists(s):
                shutil.copy(s, os.path.join(d, stem + ext))
        brd = os.path.join(d, stem + ".kicad_pcb")
        # the zones are refilled by kicad-cli and the board saved (in the corpus copy), as openpcbsi's file
        # reader does, so that both readers see the same fills: a difference is then a reading difference,
        # not one of fill engines (pcbnew's in-memory filler does not regenerate teardrops, kicad-cli does)
        _cli(["pcb", "drc", "--refill-zones", "--save-board", "--format", "json", "-o", os.path.join(d, "refill_drc.json"), brd], log)
        os.environ["OPENPCBSI_NO_REFILL"] = "1"
        try:
            snap = snapshot.extract(brd)
        finally:
            os.environ.pop("OPENPCBSI_NO_REFILL", None)
        with gzip.open(os.path.join(d, "snapshot_pcbnew.json.gz"), "wt") as f:
            json.dump(snap, f)
        if e.get("ipc2581", True):
            _cli(["pcb", "export", "ipc2581", "-o", os.path.join(d, "ipc2581.xml"), brd], log)
        _cli(["pcb", "drc", "--format", "json", "--all-track-errors", "--refill-zones", "-o", os.path.join(d, "drc.json"), brd], log)
        g = os.path.join(d, "gerbers")
        os.makedirs(g)
        _cli(["pcb", "export", "gerbers", "-o", g + os.sep, brd], log)
        _cli(["pcb", "export", "drill", "-o", g + os.sep, brd], log)
        _cli(["pcb", "export", "pos", "--format", "csv", "--units", "mm", "-o", os.path.join(d, "pos.csv"), brd], log)
        meta = {"name": e["name"], "source": e["path"], "board_file": stem + ".kicad_pcb", "sha256": snap["sha256"],
                "kicad": snap["kicad"], "snapshot_version": snap["version"], "recorded": time.strftime("%Y-%m-%d %H:%M"),
                "pads": len(snap["pads"]), "seconds": round(time.time() - t, 1)}
        json.dump(meta, open(os.path.join(d, "meta.json"), "w"), indent=1)
        log(f"  {e['name']}: {meta['pads']} pads, {len(snap['tracks'])} track segments, {meta['seconds']} s")
        done.append(e["name"])
    return done


if __name__ == "__main__":
    a = sys.argv[1:]
    if not a or a[0] == "list":
        for e in manifest():
            r = recorded(e["name"])
            print(f"{e['name']:14s} {'recorded ' + r['recorded'] if r else 'not recorded'}  {e['why']}")
    elif a[0] == "build":
        build(a[1:] or None)
