"""The board open in KiCad, through KiCad's IPC API (#32): evaluate a live board without saving it.

    python3 -m openpcbsi evaluate --live [--socket PATH] [evaluate options]
    from openpcbsi import api; snap = api.live_snapshot()

The kipy client runs in its own Python (OPENPCBSI_KIPY_PYTHON, by default the kicad-python venv in
Dev_Tools/virtual_environments) through deps/kipy_driver.py, so openpcbsi never imports kipy. The
driver takes the open board as KiCad file text (unsaved edits included); openpcbsi writes it beside
copies of the project's .kicad_pro / .kicad_dru (netclasses and custom rules still apply) and reads it
like any board file, so a live board and the same board saved give the same findings.

Needs a running KiCad with its API server enabled (Preferences > Plugins > Enable IPC API server).
"""
import json
import os
import shutil
import subprocess

from .snapshot import InputError

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
KIPY_PY = os.environ.get("OPENPCBSI_KIPY_PYTHON",
                         "/media/fpgadev/Dev_Tools/virtual_environments/kicad-python/bin/python")
DRIVER = os.path.join(os.environ.get("OPENPCBSI_DEPS") or os.path.join(HERE, "deps"), "kipy_driver.py")


def fetch_board(cache_dir, socket=None, timeout=120):
    """Write the live board into the cache and return the path of the board file to read."""
    if not os.path.exists(KIPY_PY):
        raise InputError(f"no kicad-python (kipy) Python at {KIPY_PY}; set OPENPCBSI_KIPY_PYTHON")
    tmp = os.path.join(cache_dir, "live", "incoming")
    os.makedirs(tmp, exist_ok=True)
    cmd = [KIPY_PY, "-I", DRIVER, tmp] + (["--socket", socket] if socket else [])
    p = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    if p.returncode != 0:
        msg = (p.stderr or p.stdout).strip().splitlines()[-1:] or ["no output"]
        raise InputError("could not read the board from KiCad's IPC API (is KiCad running with its API server "
                         f"enabled, and a board open?): {msg[0]}")
    meta = json.load(open(os.path.join(tmp, "live.json")))
    import hashlib
    h = hashlib.sha256(open(os.path.join(tmp, "live.kicad_pcb"), "rb").read()).hexdigest()[:16]
    name = os.path.splitext(meta.get("name") or "live")[0] or "live"
    d = os.path.join(cache_dir, "live", h)
    os.makedirs(d, exist_ok=True)
    dst = os.path.join(d, name + ".kicad_pcb")
    shutil.copy(os.path.join(tmp, "live.kicad_pcb"), dst)
    proj = meta.get("project_path")
    if proj and os.path.isdir(proj):
        for ext in (".kicad_pro", ".kicad_dru"):
            src = os.path.join(proj, name + ext)
            if os.path.exists(src):
                shutil.copy(src, os.path.join(d, name + ext))
    json.dump(dict(meta, sha256_16=h), open(os.path.join(d, "live.json"), "w"))
    return dst
