"""Paths in files meant to be shared (reports, JSON): relative to the folder they are written to, so
that a report sent to a fabricator or a colleague does not show the author's home or disk layout.

The objects in memory keep absolute paths; only what is written is rewritten. Set
OPENPCBSI_ABSOLUTE_PATHS=1 to write absolute paths (for debugging on one's own machine).
"""
import os
import re

NET_KEYS = {"net", "nets", "aggressor", "victim", "first_net", "last_net", "net_a", "net_b", "pair"}   # never paths
_PATH = re.compile(r"(?<![\w.~])(/(?:[^\s'\"`<>|,;()\[\]{}]+/)*[^\s'\"`<>|,;()\[\]{}]*)")


def _rel(p, base):
    try:
        return os.path.relpath(p, base)
    except ValueError:
        return os.path.basename(p)


def portable(obj, base):
    """A copy of `obj` (dicts, lists, strings) in which every absolute path, alone or inside a sentence,
    is made relative to `base`. A token counts as a path when it lies under `base` or the home folder,
    or when it has three or more parts and it, or its folder, exists here (so that hierarchical net
    names such as "/DIFF_N" stay as they are)."""
    if os.environ.get("OPENPCBSI_ABSOLUTE_PATHS"):
        return obj
    base = os.path.abspath(base)
    home = os.path.expanduser("~")
    seen = {}

    def fix_token(m):
        t = m.group(1)
        if t in seen:
            return seen[t]
        tt = t.rstrip(".:")
        tail = t[len(tt):]
        # a hierarchical net name ("/DIFF_N", "/serdes/LANE0_P") looks like a path: a token counts as one only
        # under the output or home folder, or with three or more parts that exist here (or whose folder does)
        parts = [q for q in tt.split("/") if q]
        ok = len(tt) > 1 and (tt.startswith(base + os.sep) or tt.startswith(home + os.sep)
                              or (len(parts) >= 3 and (os.path.exists(tt) or os.path.isdir(os.path.dirname(tt)))))
        seen[t] = (_rel(tt, base) + tail) if ok else t
        return seen[t]

    def walk(x):
        if isinstance(x, str):
            return _PATH.sub(fix_token, x) if "/" in x else x
        if isinstance(x, dict):
            return {k: (v if k in NET_KEYS else walk(v)) for k, v in x.items()}
        if isinstance(x, (list, tuple)):
            return [walk(v) for v in x]
        return x
    return walk(obj)
