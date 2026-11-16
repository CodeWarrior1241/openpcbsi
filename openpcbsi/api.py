"""Library API (#32): evaluate a board snapshot without a file, for tools that call openpcbsi from
inside a layout (Backplane, an IPC plugin, a script).

    from openpcbsi import api, snapshot
    snap = snapshot.load_or_extract("board.kicad_pcb", cache_dir)      # or a snapshot built elsewhere
    res = api.evaluate(snap, config=None, families=["rp", "g1"], changed={"nets": ["LANE0_TX_P"]})
    for f in res["findings"]: ...

`snapshot` is the dict format of snapshot.py (documented there). `config` is a path, a dict or None
(detected from the board). `families` as for `evaluate --families` (default: the rule checks, which
run in seconds). `changed` limits the returned findings to the given nets and references, for an edit
loop. Finding IDs are the same as the command line's, so a caller can track findings across edits.

An adapter that builds the snapshot from a live board through KiCad's IPC API (kipy) is the remaining
part of #32; it needs kipy and a running KiCad.
"""
import copy
import os
import tempfile

from . import evaluate as ev_mod
from . import run as run_mod
from . import stage as stage_mod


def live_snapshot(socket=None):
    """The snapshot of the board open in KiCad (IPC API, unsaved edits included); see live.py."""
    from . import live, snapshot as snap_mod
    path = live.fetch_board(run_mod.CACHE, socket)
    return snap_mod.load_or_extract(path, run_mod.CACHE)


def evaluate(snapshot, config=None, families=None, changed=None, out_dir=None, log=None):
    snap = copy.deepcopy(snapshot)
    b, cfg, ctx = run_mod.configure(snap, config)
    fams = [f.lower() for f in families] if families else ["g1", "rp"]
    checks = sorted({c for f in fams for c in ev_mod.FAMILY_CHECKS.get(f, [])})
    od = out_dir or tempfile.mkdtemp(prefix="api_", dir=os.path.join(run_mod.CACHE))
    rep, _, od = run_mod.run_loaded(b, cfg, ctx, checks, od, log=log or (lambda *a, **k: None))
    findings, val = ev_mod.collect(rep, cfg, checks)
    st = stage_mod.analyse(b.snap)
    ev_mod.apply_stage(findings, st)
    findings = [f for f in findings if f["family"].lower() in fams]          # only the families asked for
    if changed:
        nets, refs = set(changed.get("nets", [])), set(changed.get("refs", []))
        findings = [f for f in findings if set(f["nets"]) & nets or set(f["refs"]) & refs]
    per, fam = ev_mod.scorecard(findings, val, {})
    return {"findings": findings, "checks": per, "families": fam, "stage": st["stage"],
            "incomplete_nets": st["incomplete_nets"], "report": rep, "dir": od}
