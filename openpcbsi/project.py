"""The KiCad project around a board (#2): netclasses from the .kicad_pro (pattern and explicit
assignments, which pcbnew's LoadBoard does not apply) and custom rules from the .kicad_dru.

openpcbsi's rules ask `limit(rule_key, net)`; a .kicad_dru constraint that covers the net wins over
the openpcbsi config, which wins over openpcbsi's default. Conditions are understood in the forms
KiCad's rule editor writes for nets and classes:

    A.NetClass == 'HS'     A.hasNetclass('HS')     A.NetName == 'CLK*'     A.inDiffPair('/USB*')
    joined with && and || (&& binds tighter), grouped with parentheses; ! negates; != for the two
    comparisons. inDiffPair follows KiCad (DRC_ENGINE::MatchDpSuffix, pcbexpr_functions.cpp): the net
    ends in P/N/+/- followed only by digits and '_', its complement net exists, and the base name, or the
    base name without its trailing '_', matches the pattern exactly ('*' and '?' wildcards, no implied
    '*'). A rule whose condition is not understood (B.*, insideCourtyard(), Type, ...) is listed in the
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


def wx_matches(text, pattern):
    """wxString::Matches: '*' any run, '?' one character, everything else literal, case-sensitive."""
    rx = "".join(".*" if c == "*" else "." if c == "?" else re.escape(c) for c in pattern)
    return re.fullmatch(rx, text, re.S) is not None


def dp_suffix(net):
    """KiCad's DRC_ENGINE::MatchDpSuffix: (polarity, base name, complement net) or (0, None, None).
    Trailing digits and '_' are skipped; the first other character must be P, N, + or -."""
    count = 0
    for ch in reversed(net):
        count += 1
        if ch.isdigit() and ch.isascii() or ch == "_":
            continue
        comp = {"+": "-", "-": "+", "N": "P", "P": "N"}.get(ch)
        if comp is None:
            break
        base = net[:len(net) - count]
        return (1 if ch in "+P" else -1), base, base + comp + net[len(net) - count + 1:]
    return 0, None, None


_ATOM = re.compile(r"""\s*(?:
      (?P<lp>\() | (?P<rp>\)) | (?P<and>&&) | (?P<or>\|\|) | (?P<not>!(?!=))
    | A\.(?P<prop>NetClass|NetName)\s*(?P<op>==|!=)\s*'(?P<val>[^']*)'
    | A\.(?P<fn>hasNetclass|inDiffPair)\(\s*'(?P<arg>[^']*)'\s*\)
    )""", re.X)


def _cond_tokens(cond):
    """Tokens of a condition, or None if any part of it is outside the understood forms."""
    out, i = [], 0
    while i < len(cond):
        if cond[i:].strip() == "":
            break
        m = _ATOM.match(cond, i)
        if not m or m.end() == i:
            return None
        if m.group("prop"):     # (kind, pattern, negated): NetClass/NetName == or !=
            out.append(("atom", ("classeq" if m.group("prop") == "NetClass" else "neteq", m.group("val"), m.group("op") == "!=")))
        elif m.group("fn"):
            out.append(("atom", ("hasclass" if m.group("fn") == "hasNetclass" else "pair", m.group("arg"), False)))
        else:
            out.append((next(k for k in ("lp", "rp", "and", "or", "not") if m.group(k)), None))
        i = m.end()
    return out


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
        self.via_height_in_length = True   # KiCad's rules.use_height_for_length_calcs (default true)
        pro, dru = stem + ".kicad_pro", stem + ".kicad_dru"
        if os.path.exists(pro):
            self.files[pro] = _sha(pro)
            try:
                pj = json.load(open(pro))
            except (ValueError, OSError):
                pj = {}
            ns = pj.get("net_settings", {}) or {}
            rules = ((pj.get("board") or {}).get("design_settings") or {}).get("rules") or {}
            if isinstance(rules.get("use_height_for_length_calcs"), bool):
                self.via_height_in_length = rules["use_height_for_length_calcs"]
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
                            elif isinstance(y, list) and len(y) == 1 and isinstance(y[0], str):
                                b.setdefault("flags", []).append(y[0])      # e.g. (within_diff_pairs)
                        cons.append(b)
            pred = self._condition(cond) if cond is not None else (lambda net, cls: True)
            r = {"name": name, "condition": cond, "constraints": cons, "layer": layer, "pred": pred}
            if pred is None:
                self.not_applied.append({"rule": name, "condition": cond, "why": "condition not understood"})
                continue
            self.rules.append(r)

    def set_nets(self, names):
        """The board's net names: inDiffPair() needs the complement net to exist, as in KiCad."""
        self.net_names = set(names)

    def _in_diff_pair(self, net, pat):
        pol, base, comp = dp_suffix(net)
        if not pol:
            return False
        names = getattr(self, "net_names", None)
        if names is not None and comp not in names:
            return False
        return wx_matches(base, pat) or (base.endswith("_") and wx_matches(base[:-1], pat))

    def _atom(self, kind, pat, net, cls):
        if kind == "hasclass":                         # NETCLASS::ContainsNetclassWithName: wxString::Matches
            return wx_matches(cls or "", pat)
        if kind == "pair":
            return self._in_diff_pair(net, pat)
        val = (cls or "") if kind == "classeq" else net   # VALUE::EqualTo: strings compared ignoring case
        return wx_matches(val.lower(), pat.lower())

    def _condition(self, cond):
        """A predicate (net, netclass) -> bool, or None if the condition is not understood."""
        toks = _cond_tokens(cond)
        if not toks:
            return None
        pos = 0

        def peek():
            return toks[pos][0] if pos < len(toks) else None

        def take():
            nonlocal pos
            pos += 1
            return toks[pos - 1]

        def factor():
            k = peek()
            if k == "not":
                take()
                f = factor()
                return None if f is None else (lambda n, c: not f(n, c))
            if k == "lp":
                take()
                e = expr()
                if e is None or peek() != "rp":
                    return None
                take()
                return e
            if k == "atom":
                kind, pat, ne = take()[1]
                return lambda n, c: self._atom(kind, pat, n, c) != ne
            return None

        def term():
            f = factor()
            while f is not None and peek() == "and":
                take()
                g = factor()
                if g is None:
                    return None
                f = (lambda a, b: lambda n, c: a(n, c) and b(n, c))(f, g)
            return f

        def expr():
            t = term()
            while t is not None and peek() == "or":
                take()
                u = term()
                if u is None:
                    return None
                t = (lambda a, b: lambda n, c: a(n, c) or b(n, c))(t, u)
            return t

        e = expr()
        if e is None or pos != len(toks):
            return None
        return e

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

    def skew_rules(self):
        """[(rule name, max mm, within_diff_pairs, pred)] of the .kicad_dru skew constraints. KiCad judges
        the skew across every net the rule covers, or per pair with (within_diff_pairs)."""
        return [(r["name"], c["max"], "within_diff_pairs" in c.get("flags", []), r["pred"])
                for r in self.rules for c in r["constraints"] if c["type"] == "skew" and c.get("max") is not None]

    def design_gap(self, net, cls):
        """(gap mm, source) the pair is designed to: the last .kicad_dru diff_pair_gap opt covering the
        net, else the diff_pair_gap of its netclass in the .kicad_pro; None if neither gives one. The
        Default class does not count: every project has one, with KiCad's own figure."""
        hit = None
        for r in self.rules:
            for c in r["constraints"]:
                if c["type"] == "diff_pair_gap" and c.get("opt") is not None and r["pred"](net, cls):
                    hit = (c["opt"], f".kicad_dru rule '{r['name']}' diff_pair_gap opt")
        if hit:
            return hit
        if not cls or cls == "Default":
            return None
        nc = self.netclasses.get(cls) or {}
        g = nc.get("diff_pair_gap")
        if isinstance(g, (int, float)) and g > 0:
            return float(g), f"netclass {cls or 'Default'} diff_pair_gap (.kicad_pro)"
        return None

    def describe(self):
        return {"files": self.files, "netclasses": sorted(self.netclasses),
                "patterns": [{"pattern": p, "netclass": c} for p, c in self.patterns],
                "rules": [{"name": r["name"], "condition": r["condition"],
                           "constraints": [{k: v for k, v in c.items()} for c in r["constraints"]]} for r in self.rules],
                "rules_not_applied": self.not_applied}


def load(board_path):
    return Project(board_path)
