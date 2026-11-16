"""The KiCad project around a board (#2): netclasses from the .kicad_pro (pattern and explicit
assignments, which pcbnew's LoadBoard does not apply) and custom rules from the .kicad_dru.

openpcbsi's rules ask `limit(rule_key, net)`; a .kicad_dru constraint that covers the net wins over
the openpcbsi config, which wins over openpcbsi's default. Conditions are understood in the forms
KiCad's rule editor writes for nets and classes:

    A.NetClass == 'HS'     A.hasNetclass('HS')     A.NetName == 'CLK*'     A.inDiffPair('/USB*')
    joined with && and ||; a leading ! negates. A rule whose condition is not understood is listed in the
    report as not applied, never applied to every net.
"""
import fnmatch
import hashlib
import json
import os
import re

# .kicad_dru constraint -> (openpcbsi rule, parameter, which bound)
MAP = {
    "diff_pair_uncoupled": ("pair_coupling", "limit_mm", "max"),
    "diff_pair_gap": ("pair_coupling", "max_gap_mm", "max"),
    "edge_clearance": ("edge_distance", "board_edge_mm", "min"),
    "skew": ("timing", "skew_mm", "max"),
    "length": ("timing", "length_mm", "max"),
}
UNIT = {"mm": 1.0, "um": 1e-3, "mil": 0.0254, "in": 25.4, "": 1.0}


def _sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        h.update(f.read())
    return h.hexdigest()


def _tokens(text):
    return re.findall(r'"(?:[^"\\]|\\.)*"|\(|\)|[^\s()]+', text)


def _parse(tokens):
    stack, cur = [], []
    for t in tokens:
        if t == "(":
            stack.append(cur)
            cur = []
        elif t == ")":
            done, cur = cur, stack.pop()
            cur.append(done)
        else:
            cur.append(t[1:-1].replace('\\"', '"') if t.startswith('"') else t)
    return cur


def _num(v):
    m = re.fullmatch(r"(-?[\d.]+)\s*(mm|um|mil|in)?", str(v))
    return float(m.group(1)) * UNIT[m.group(2) or ""] if m else None


class Project:
    def __init__(self, board_path):
        self.board = os.path.abspath(board_path)
        stem = os.path.splitext(self.board)[0]
        self.files = {}
        self.netclasses = {}          # name -> settings from the .kicad_pro
        self.patterns = []            # (pattern, class)
        self.assigned = {}            # net -> class (explicit)
        self.rules = []               # parsed .kicad_dru rules
        self.not_applied = []         # rules whose condition was not understood
        pro, dru = stem + ".kicad_pro", stem + ".kicad_dru"
        if os.path.exists(pro):
            self.files[pro] = _sha(pro)
            try:
                ns = json.load(open(pro)).get("net_settings", {})
            except (ValueError, OSError):
                ns = {}
            for c in ns.get("classes", []):
                self.netclasses[c.get("name")] = c
            for p in ns.get("netclass_patterns", []) or []:
                self.patterns.append((p.get("pattern", ""), p.get("netclass")))
            for net, cls in (ns.get("netclass_assignments") or {}).items():
                self.assigned[net] = cls[0] if isinstance(cls, list) and cls else cls
        if os.path.exists(dru):
            self.files[dru] = _sha(dru)
            self._read_dru(open(dru).read())
        sch = stem + ".kicad_sch"
        if os.path.exists(sch):
            self.files[sch] = _sha(sch)

    # ---- netclasses ---------------------------------------------------------------------------
    def netclass(self, net, board_class=None):
        """Effective netclass: explicit assignment, then the first matching pattern, then the class
        the board file gives."""
        if net in self.assigned:
            return self.assigned[net]
        for pat, cls in self.patterns:
            if fnmatch.fnmatchcase(net, pat) or fnmatch.fnmatchcase(net.lstrip("/"), pat.lstrip("/")):
                return cls
        return board_class or "Default"

    # ---- custom rules -------------------------------------------------------------------------
    def _read_dru(self, text):
        for item in _parse(_tokens(text)):
            if not isinstance(item, list) or not item or item[0] != "rule":
                continue
            name = item[1] if len(item) > 1 and isinstance(item[1], str) else "?"
            cond, cons, layer = None, [], None
            for x in item[2:]:
                if isinstance(x, list) and x:
                    if x[0] == "condition" and len(x) > 1:
                        cond = x[1]
                    elif x[0] == "layer" and len(x) > 1:
                        layer = x[1]
                    elif x[0] == "constraint" and len(x) > 1:
                        b = {"type": x[1]}
                        for y in x[2:]:
                            if isinstance(y, list) and len(y) == 2 and y[0] in ("min", "max", "opt"):
                                b[y[0]] = _num(y[1])
                        cons.append(b)
            pred = self._condition(cond) if cond is not None else (lambda net, cls: True)
            r = {"name": name, "condition": cond, "constraints": cons, "layer": layer, "pred": pred}
            if pred is None:
                self.not_applied.append({"rule": name, "condition": cond, "why": "condition not understood"})
                continue
            self.rules.append(r)

    def _condition(self, cond):
        """A predicate (net, netclass) -> bool, or None if the condition is not understood."""
        terms = []
        for disj in re.split(r"\|\|", cond):
            conj = []
            for t in re.split(r"&&", disj):
                t = t.strip().strip("()").strip()
                neg = t.startswith("!")
                t = t.lstrip("!").strip()
                m = (re.fullmatch(r"A\.NetClass\s*==\s*'([^']*)'", t) or re.fullmatch(r"A\.hasNetclass\('([^']*)'\)", t))
                if m:
                    conj.append((neg, "class", m.group(1)))
                    continue
                m = re.fullmatch(r"A\.NetName\s*==\s*'([^']*)'", t)
                if m:
                    conj.append((neg, "net", m.group(1)))
                    continue
                m = re.fullmatch(r"A\.inDiffPair\('([^']*)'\)", t)
                if m:
                    conj.append((neg, "pair", m.group(1)))
                    continue
                return None
            terms.append(conj)

        def pred(net, cls):
            for conj in terms:
                okc = True
                for neg, kind, pat in conj:
                    if kind == "class":
                        v = fnmatch.fnmatchcase(cls or "", pat)
                    elif kind == "net":
                        v = fnmatch.fnmatchcase(net, pat)
                    else:                     # inDiffPair('X*'): the pair's base name matches X*
                        base = re.sub(r"(_?[PN]|[+-])$", "", net)
                        v = fnmatch.fnmatchcase(base, re.sub(r"\*$", "", pat) + "*") or fnmatch.fnmatchcase(net, pat)
                    if v == neg:
                        okc = False
                        break
                if okc:
                    return True
            return False
        return pred

    def limit(self, rule_key, param, net, cls):
        """(value, source) from the last .kicad_dru rule covering `net` (KiCad: later rules win), or None."""
        hit = None
        for r in self.rules:
            for c in r["constraints"]:
                m = MAP.get(c["type"])
                if not m or m[0] != rule_key or m[1] != param or c.get(m[2]) is None:
                    continue
                if r["pred"](net, cls):
                    hit = (c[m[2]], f".kicad_dru rule '{r['name']}' {c['type']} {m[2]}")
        return hit

    def describe(self):
        return {"files": self.files, "netclasses": sorted(self.netclasses),
                "patterns": [{"pattern": p, "netclass": c} for p, c in self.patterns],
                "rules": [{"name": r["name"], "condition": r["condition"],
                           "constraints": [{k: v for k, v in c.items()} for c in r["constraints"]]} for r in self.rules],
                "rules_not_applied": self.not_applied}


def load(board_path):
    return Project(board_path)
