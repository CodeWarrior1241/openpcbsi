"""File readers for KiCad boards (plan phase 2): openpcbsi reads the .kicad_pcb file itself instead of
through KiCad's SWIG `pcbnew` module, which KiCad has deprecated.

    monkey.py       primary reader, on kicad-monkey (deps/parsers.lock pins it); produces the snapshot
    kicadfiles.py   independent cross-check on KiCadFiles: the fields it can read
    shapes.py       pad shapes, margins and outlines as polygons, computed by openpcbsi for every reader

The pinned parsers and their dependencies live in deps/ (deps/setup_parsers.sh); `use()` puts them on
the import path. Only this package imports a parser, so a parser's API changes stay here.
"""
import os
import sys
import tomllib

HERE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DEPS = os.environ.get("OPENPCBSI_DEPS") or os.path.join(HERE, "deps")


def lock():
    return tomllib.load(open(os.path.join(DEPS, "parsers.lock"), "rb"))


def use():
    """Put the pinned parsers (and kicad-monkey's dependencies) on sys.path; returns the lock."""
    lk = lock()
    for p in (lk["site"]["dir"], lk["kicad_monkey"]["python_path"], lk["kicadfiles"]["python_path"]):
        path = os.path.join(DEPS, p)
        if not os.path.isdir(path):
            from ..snapshot import InputError     # reported as a message, not a traceback
            raise InputError(f"the board file parsers are not set up ({path} is missing): run deps/setup_parsers.sh "
                             "once from the repository (see README, Quick start)")
        if path not in sys.path:
            sys.path.insert(0, path)
    return lk
