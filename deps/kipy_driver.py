"""openpcbsi -> KiCad IPC API driver (#32). Runs in a Python that has kicad-python (kipy), never in
openpcbsi's own:

    <kipy python> -I deps/kipy_driver.py OUT_DIR [--socket PATH]

Connects to a running KiCad (its API server must be enabled), takes the board that is open in the PCB
editor as KiCad file text (unsaved edits included), and writes OUT_DIR/live.kicad_pcb and
OUT_DIR/live.json (the board's own file name, its project folder, KiCad's version).
"""
import json
import os
import sys


def main(out_dir, socket=None):
    from kipy import KiCad
    kw = {"socket_path": socket} if socket else {}
    kicad = KiCad(**kw)
    board = kicad.get_board()
    text = board.get_as_string()
    os.makedirs(out_dir, exist_ok=True)
    with open(os.path.join(out_dir, "live.kicad_pcb"), "w", encoding="utf-8") as f:
        f.write(text)
    meta = {"name": getattr(board, "name", None), "kicad_version": str(kicad.get_version())}
    try:
        proj = board.get_project()
        meta["project_path"] = getattr(proj, "path", None)
        meta["project_name"] = getattr(proj, "name", None)
    except Exception as e:
        meta["project_error"] = str(e)
    with open(os.path.join(out_dir, "live.json"), "w") as f:
        json.dump(meta, f)
    print(json.dumps(meta))


if __name__ == "__main__":
    a = sys.argv[1:]
    sock = None
    if "--socket" in a:
        i = a.index("--socket")
        sock = a[i + 1]
        del a[i:i + 2]
    main(a[0], sock)
