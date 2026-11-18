"""Design evaluation run (REQUIREMENTS 5.10-5.16): G18 scorecard, G19 finding model, G20 waivers,
G21 exit codes, EV1 KiCad DRC in the run, SD1 DC drop judged.

    python3 -m openpcbsi evaluate BOARD [-c CONFIG] [--waivers FILE] [--families g1,rp,sd,si,th,emc,ev]
                                       [--checks ...] [--json OUT] [-o OUT_DIR] [--no-markers]

Writes <run>/evaluation.json (schema `openpcbsi-evaluation/1`, docs/evaluation.md) and
<run>/evaluation.html next to the usual report.json / report.html / marker board.

Exit codes: 0 nothing fails or warns; 1 warnings only (0 if `[evaluate] warnings_exit = false`);
2 any fail; 3 the board could not be evaluated.

Not built yet (reported as "not built" in the scorecard, never as a pass): the placement family
(PL), the field-solver gates and automatic regions (FS), SD2-SD9, RP3 percentiles, EV2-EV4
(version-difference DRC, claims file, A/B scorecard), G23 marker shapes, G24 library API.
"""
import datetime
import fnmatch
import hashlib
import html
import json
import os
import re
import subprocess
import tomllib

SCHEMA = "openpcbsi-evaluation/1"
SEV = {"error": "fail", "warning": "warn", "info": "info", "fail": "fail", "warn": "warn", "pass": "pass"}
ORDER = {"fail": 6, "warn": 5, "pending": 4.5, "unvalidated": 4, "not judged": 3, "assumed inputs": 2, "nothing to check": 1,
         "not run": 1, "not built": 1, "pass": 0}
NOT_REAL = ("not built", "not run", "nothing to check")
FAMILY_ORDER = ["G1", "RP", "SD", "SI", "TH", "EMC", "EV", "PL", "FS", "MF"]

# G1 rule key -> (family, requirement) as REQUIREMENTS 5.13 names them; the rest stay G1
RULE_FAMILY = {
    "pad_plane": ("RP", "RP1"), "return_via": ("RP", "RP2"), "stitching": ("RP", "RP4/RP5"),
    "ref_plane": ("RP", "RP6"), "split_crossing": ("RP", "RP6"), "pair_coupling": ("RP", "RP7"),
    "parallel": ("RP", "RP8"), "via_stub": ("RP", "RP9"), "stitching_map": ("RP", "RP3"),
    "return_proximity": ("RP", "RP (no plane)"), "breakout": ("SI", "breakout"),
    "unreferenced": ("G1", "G1"), "edge_distance": ("G1", "PC8"), "decap": ("G1", "PC11"),
    "thermal_vias": ("G1", "PC13"), "power_overlap": ("G1", "G1"), "stubs": ("G1", "G1"),
    # manufacturing family (5.17): run by `openpcbsi dfm`, which keeps its own scorecard
    "rule_audit": ("MF", "MF19"),
    "output_match": ("MF", "MF20"),
    "fab_track": ("MF", "MF3"),
    "fab_space": ("MF", "MF12"),
    "same_net_gap": ("MF", "MF2"),
    "fab_drill": ("MF", "MF13"),
    "mask_web": ("MF", "MF4"),
    "mask_exposed": ("MF", "MF5"),
    "via_pad": ("MF", "MF14"),
    "paste": ("MF", "MF15"),
    "footprint": ("MF", "MF16"),
    "part_spacing": ("MF", "MF17"),
    "assembly": ("MF", "MF18"),
    "acid_trap": ("MF", "MF1"),
    "copper": ("MF", "MF21"),
    "reflow": ("MF", "MF7"),
    "chips": ("MF", "MF8"),
}
CHECK_FAMILY = {
    "impedance": ("SI", "G6"), "crosstalk": ("SI", "G4"), "si": ("SI", "G2/G3/G7"), "vias": ("SI", "G9"),
    "dc": ("SD", "SD1"), "pdn": ("SD", "G11"), "thermal": ("TH", "G13"), "emc": ("EMC", "G14"),
    "kicad_drc": ("EV", "EV1"), "em_region": ("FS", "FS9"), "trace_heating": ("SD", "SD3"), "via_current": ("SD", "SD4"), "timing": ("SI", "timing"),
}
FAMILY_CHECKS = {"g1": ["drc"], "rp": ["drc"], "sd": ["dc", "pdn"], "si": ["impedance", "crosstalk", "si", "vias", "timing"],
                 "th": ["thermal"], "emc": ["pdn", "emc"], "ev": []}
NOT_BUILT = [
    ("PL", "PL1-PL10", "placement evaluation"),
    ("FS", "FS1-FS10", "field-solver gates and automatic regions"),
    ("SD", "SD2-SD9", "BGA supply entry, hot copper, current per via, shunt isolation, Kelvin sense, sense side, PDN band, necks"),
    ("EV", "EV2-EV4", "DRC on the downgraded copy, claims file, A/B scorecard"),
]
FAMILY_TITLE = {"G1": "Electrical rules (G1)", "RP": "Return path and grounding (RP)", "SD": "Supply delivery (SD)",
                "SI": "Signal integrity (G2-G9)", "TH": "Thermal (G13)", "EMC": "EMC (G14)", "EV": "Process (EV)",
                "PL": "Placement (PL)", "FS": "Field-solver gates (FS)", "MF": "Manufacturing (MF)"}


# ------------------------------------------------------------------------------------------------
def finding_id(check, nets=(), refs=(), layer=None, x=None, y=None, extra=None):
    """G19 stable ID: the check plus a hash of nets, references and location rounded to 0.1 mm."""
    key = [check, sorted(n for n in nets if n), sorted(r for r in refs if r), layer or "",
           None if x is None else round(float(x), 1), None if y is None else round(float(y), 1), extra or ""]
    return f"{check}-{hashlib.sha1(json.dumps(key).encode()).hexdigest()[:10]}"


def _f(check, family, req, sev, reason, nets=(), refs=(), layer=None, x=None, y=None, value=None, limit=None,
       unit=None, area="", evidence=None, best=None, extra=None):
    return {"id": finding_id(check, nets, refs, layer, x, y, extra), "severity": sev, "check": check,
            "family": family, "requirement": req, "value": value, "limit": limit, "unit": unit, "best": best,
            "nets": [n for n in nets if n], "refs": [r for r in refs if r], "layer": layer, "x": x, "y": y,
            "reason": reason, "evidence": evidence, "area": area or "", "waived": None}


HINTS = {   # #31: a first fix for each geometric finding (advice, not a constraint)
    "split_crossing": "route over unbroken reference copper on {ref_layer}, or bridge the gap with a stitching capacitor next to the crossing",
    "ref_plane": "restore solid {ref_layer} copper under the track, or move the track over solid copper",
    "return_via": "add a ground via within 1 mm of the layer change, joining the two reference planes",
    "pad_plane": "void the reference under the pad (pad + 0.15 mm) so the pad references the next ground layer",
    "pair_coupling": "keep P and N at the pair gap through the fan-out; match skew with small bumps close to the mismatch",
    "stitching": "add stitching vias at no more than {limit} mm spacing along this stretch",
    "parallel": "space the other net at least {limit} mm away, or shorten the parallel run",
    "edge_distance": "keep the track at least {limit} mm inside its reference copper and the board edge",
    "decap": "move the capacitor's via within {limit} mm of its pad",
    "via_stub": "back-drill the via, or change to the layer that uses the whole via length",
    "unreferenced": "route the segment over a plane, or add reference copper on the adjacent layer",
    "stubs": "remove the dangling end or complete the connection",
    "impedance": "adjust width/gap for {target} ohm: `openpcbsi plan BOARD --layer {layer} --zdiff {target}`",
    "thermal": "spread the heat: thermal vias and copper under the part, airflow or a heatsink; check the dissipation",
    "vias": "change the antipad or the return-via spacing; model it in 3D before relying on the closed form",
    "acid_trap": "meet the copper at 90 degrees or more: enter the pad square-on, or add a teardrop or fillet at the corner",
    "same_net_gap": "close the gap with copper or widen it to at least {limit} mm",
    "fab_track": "widen the track to at least {limit} mm, or choose a process of the fabricator that makes this width",
    "fab_space": "move the copper at least {limit} mm apart, or choose a process of the fabricator that makes this spacing",
    "fab_drill": "use a hole, pad or spacing the fabricator makes: at least {limit} {unit}",
    "mask_web": "increase the pad spacing or reduce the mask expansion so that a web of the minimum width remains",
    "mask_exposed": "move the neighbouring copper away from the mask opening, or reduce the mask expansion",
    "via_pad": "move the via off the pad and tent it, or order filled and capped vias and say so in [dfm] via_cover",
    "paste": "enlarge the aperture, or order a thinner stencil and set [dfm] stencil_mm",
    "footprint": "use the land pattern of the part's data sheet",
    "part_spacing": "move the parts at least {limit} mm apart",
    "assembly": "choose the assembly service that takes this part, or change the part",
    "rule_audit": "set the board's design rules to the fabricator's limits, so that KiCad's DRC checks them while routing",
    "output_match": "plot the files again from this board and send those; if copper has no net, give it one or check it by hand",
    "reflow": "check the part's mass in [reflow.parts]; a heavy part next to small ones needs a longer soak or a different position",
    "chips": "bring the two pads to the same copper: the same track width, the same via arrangement, a thermal relief in a pour",
}


class _Safe(dict):
    def __missing__(self, k):
        return "?"


def hint(rule, f):
    t = HINTS.get(rule)
    return t.format_map(_Safe({k: v for k, v in f.items() if v is not None})) if t else None


def _layer_png(layer):
    return f"layer_{layer.replace('.', '_')}.png" if layer and "/" not in layer else None


def from_rule_findings(check, family, req, fs):
    out = []
    for f in fs:
        nets = [f.get("net"), f.get("aggressor")]
        refs = [f.get("part"), (f.get("pad") or "").split(".")[0] or None]
        out.append(_f(check, family, req, SEV.get(f.get("severity"), "info"), f.get("msg", ""), nets, refs,
                      f.get("layer"), f.get("x"), f.get("y"), f.get("value"), f.get("limit"), f.get("unit"),
                      f.get("area", ""), evidence=_layer_png(f.get("layer")) or "marker board User.9",
                      extra=f.get("kind") or f.get("rule")))
        if out[-1]["severity"] in ("fail", "warn"):
            out[-1]["hint"] = hint(check, f)
    return out


def judge_dc(rep, cfg):
    """SD1: drop as a share of the window. fail > 100 %, warn > warn_pct; at an assumed current never
    worse than warn. A disagreement between the two engines above `engine_disagree_pct` is a warn."""
    ev = cfg.data.get("evaluate", {})
    warn_pct = ev.get("dc_warn_pct", 50.0)
    dis = ev.get("engine_disagree_pct", 20.0)
    d = rep["checks"].get("dc")
    if not d:
        return [], None
    out = []
    for r in d.get("rails", []):
        rail = cfg["rails"].get(r["rail"], {})
        assumed = rail.get("current_assumed", False)
        refs = [rail.get("source")] + list(rail.get("loads", []))
        if r.get("problems") or r.get("problem"):
            out.append(_f("dc", "SD", "SD1", "warn", f"{r['rail']}: not evaluated: {r.get('problems') or r.get('problem')}",
                          [r.get("net")], refs, extra="problem"))
            continue
        pct = r.get("drop_pct_of_window")
        if pct is None:
            continue
        cur = f"{r.get('current_a', rail.get('current_a'))} A" + (" (assumed: no current in the config)" if assumed else "")
        sev = "fail" if pct > 100 else "warn" if pct > warn_pct else "pass"
        if assumed and sev == "fail":
            sev = "warn"
        msg = (f"{r['rail']}: drop {r.get('drop_mv', float('nan')):.1f} mV = {pct:.0f} % of the "
               f"+/-{rail.get('tol_pct', 3)} % window at {cur}")
        if r.get("worst_ball"):
            msg += f", worst {r['worst_ball']}"
        out.append(_f("dc", "SD", "SD1", sev, msg, [r.get("net")], refs, value=round(pct, 1),
                      limit=100.0, unit="% of window", evidence="report.html#dc", extra="drop"))
    # SD3 / SD4: copper heating and via current from the solved currents (openpcbsi's own solution)
    th = ev.get("trace_warn_c", 10.0), ev.get("trace_fail_c", 30.0)
    va = ev.get("via_amps_per_0p2mm", 0.5)
    own = {r["rail"]: r for r in (d.get("openpcbsi_rails") or [])} or {r["rail"]: r for r in d.get("rails", [])}
    for name, r in own.items():
        rail = cfg["rails"].get(name, {})
        assumed = rail.get("current_assumed", False)
        for t in r.get("top_tracks", []):
            dT = t.get("dT_ipc2221_C", 0)
            if dT > th[0]:
                sev = "fail" if dT > th[1] and not assumed else "warn"
                out.append(_f("trace_heating", "SD", "SD3", sev, f"{name}: {t['I']:.2f} A in a {t['w']} mm track on {t['layer']} "
                              f"rises ~{dT:.0f} C (IPC-2221; warn {th[0]:g}, fail {th[1]:g} C)" + (" at an assumed current" if assumed else ""),
                              [rail.get("net")], layer=t["layer"], x=t.get("x"), y=t.get("y"), value=round(dT, 1),
                              limit=th[1], unit="C", evidence="report.html#dc", extra="heat"))
        for v in r.get("top_vias", []):
            lim = va * max(v.get("drill_mm", 0.2), 0.05) / 0.2
            if v["I"] > lim:
                out.append(_f("via_current", "SD", "SD4", "warn", f"{name}: {v['I']:.2f} A in one via of {v.get('drill_mm', '?')} mm drill "
                              f"(limit {lim:.2f} A: {va:g} A per 0.2 mm, proposed)", [rail.get("net")], x=v["x"], y=v["y"],
                              value=round(v["I"], 3), limit=round(lim, 3), unit="A", extra="via"))
    for x in d.get("cross_check", []):
        if x.get("diff_pct") is not None and abs(x["diff_pct"]) > dis:
            out.append(_f("dc", "SD", "SD1", "warn", f"{x['rail']}: KiPI {x['kipi_mv']:.2f} mV and openpcbsi "
                          f"{x['openpcbsi_mv']:.2f} mV differ by {x['diff_pct']:+.1f} % (> {dis:g} %)",
                          [cfg['rails'].get(x['rail'], {}).get('net')], value=x["diff_pct"], limit=dis, unit="%",
                          evidence="report.html#dc", extra="engines"))
    return out, d.get("validation")


def judge_pdn(rep, cfg):
    """G11: |Z| at the loads against V * ripple / step, up to [pdn] target_fmax_hz."""
    d = rep["checks"].get("pdn")
    if not d:
        return [], None
    out = []
    for r in d.get("rails", []):
        if r.get("problem") or r.get("max_ratio") is None:
            if r.get("problem"):
                out.append(_f("pdn", "SD", "G11", "warn", f"{r['rail']}: not evaluated: {r['problem']}", [r.get("net")],
                              extra="problem"))
            continue
        ratio = r["max_ratio"]
        sev = "fail" if ratio > 1.0 else "pass"
        rail = cfg["rails"].get(r["rail"], {})
        if rail.get("current_assumed") and sev == "fail":
            sev = "warn"
        f_at = r.get("f_at_max_hz")
        msg = (f"{r['rail']}: max |Z| {r.get('max_Z_ohm', 0) * 1000:.2f} mOhm = {ratio:.2f} x the target "
               f"{r.get('target_ohm', 0) * 1000:.2f} mOhm" + (f" at {f_at / 1e6:.1f} MHz" if f_at else ""))
        out.append(_f("pdn", "SD", "G11", sev, msg, [r.get("net")], [r.get("vrm")], value=round(ratio, 3), limit=1.0,
                      unit="x target", evidence=f"pdn_{r['rail']}.png", extra="z"))
    return out, d.get("validation")


# ------------------------------------------------------------------------------------------------
def kicad_drc(board, out_dir, cfg, log=print):
    """EV1: kicad-cli pcb drc --severity-all --schematic-parity --refill-zones (the board is not saved).
    KiCad errors, unconnected items and parity issues -> fail; KiCad warnings -> `[evaluate]
    drc_warning` (default warn); exclusions -> info."""
    ev = cfg.data.get("evaluate", {})
    stem = os.path.splitext(board)[0]
    parity = os.path.exists(stem + ".kicad_sch")
    out = os.path.join(out_dir, "kicad_drc.json")
    cmd = ["kicad-cli", "pcb", "drc", "--severity-all", "--refill-zones", "--format", "json", "-o", out]
    if parity:
        cmd.append("--schematic-parity")
    cmd.append(board)
    log("EV1 KiCad DRC: " + " ".join(cmd[3:-1]) + ("" if parity else "  (no schematic next to the board: parity not checked)"))
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=ev.get("drc_timeout_s", 900))
    except (OSError, subprocess.TimeoutExpired) as e:
        return [_f("kicad_drc", "EV", "EV1", "warn", f"KiCad DRC could not run: {e}", extra="run")], {"ran": False}
    if not os.path.exists(out):
        return [_f("kicad_drc", "EV", "EV1", "warn", f"KiCad DRC wrote no report: {p.stderr.strip()[-300:]}",
                   extra="run")], {"ran": False}
    d = json.load(open(out))
    warn_as = ev.get("drc_warning", "warn")
    fs = []
    for group, base in (("violations", None), ("unconnected_items", "fail"), ("schematic_parity", "fail")):
        for v in d.get(group, []):
            ks = v.get("severity", "error")
            sev = "info" if ks == "exclusion" else base or ("fail" if ks == "error" else warn_as)
            items = v.get("items", [])
            text = " ".join(i.get("description", "") for i in items)
            refs = re.findall(r"\b(?:of|Footprint)\s+([A-Za-z_]+\d+)", text)
            nets = re.findall(r"\[([^\]]+)\]", text)
            pos = items[0].get("pos", {}) if items else {}
            fs.append(_f("kicad_drc", "EV", "EV1", sev, f"{v.get('type')}: {v.get('description')}", nets, refs,
                         None, pos.get("x"), pos.get("y"), evidence="kicad_drc.json", extra=v.get("type")))
    meta = {"ran": True, "kicad": d.get("kicad_version"), "parity": parity, "command": " ".join(cmd),
            "violations": len(d.get("violations", [])), "unconnected": len(d.get("unconnected_items", [])),
            "parity_issues": len(d.get("schematic_parity", []))}
    log(f"  {meta['violations']} violations, {meta['unconnected']} unconnected, {meta['parity_issues']} parity issues")
    return fs, meta


# ------------------------------------------------------------------------------------------------
def load_waivers(path):
    """G20: [[waiver]] entries with `id`, or `check` with optional `nets` / `refs` / `area` globs, and
    `reason` and `reviewer` (both required)."""
    if not path:
        return [], None
    raw = open(path, "rb").read()
    ws = tomllib.loads(raw.decode()).get("waiver", [])
    for i, w in enumerate(ws):
        w["_index"] = i
        w["_problem"] = None
        if not w.get("reason") or not w.get("reviewer"):
            w["_problem"] = "needs reason and reviewer"
        elif not (w.get("id") or w.get("check")):
            w["_problem"] = "needs id or check"
    return ws, hashlib.sha256(raw).hexdigest()[:12]


def _waiver_matches(w, f):
    if w.get("id"):
        return w["id"] == f["id"]
    if w["check"] != f["check"]:
        return False
    for key, vals in (("nets", f["nets"]), ("refs", f["refs"])):
        if w.get(key) and not any(fnmatch.fnmatchcase(v, p) for v in vals for p in w[key]):
            return False
    if w.get("area") and w["area"] != f["area"]:
        return False
    return True


def apply_waivers(findings, waivers):
    used = set()
    for f in findings:
        for w in waivers:
            if w["_problem"] is None and _waiver_matches(w, f):
                f["waived"] = {"reason": w["reason"], "reviewer": w["reviewer"], "waiver": w.get("id") or w["check"]}
                used.add(w["_index"])
                break
    out = []
    for w in waivers:
        st = "invalid: " + w["_problem"] if w["_problem"] else "applied" if w["_index"] in used else "stale"
        out.append({k: v for k, v in w.items() if not k.startswith("_")} | {"status": st})
    return out


# ------------------------------------------------------------------------------------------------
def collect(rep, cfg, checks):
    """All findings of a finished check run in the G19 model, plus each check's validation state."""
    fs, val = [], {}
    for rule in rep["checks"].get("drc", []):
        fam, req = RULE_FAMILY.get(rule["rule"], ("G1", rule.get("spec", "G1")))
        fs += from_rule_findings(rule["rule"], fam, req, rule.get("findings", []))
        # a rule that found no item in scope checked nothing: never a pass (#5)
        val[rule["rule"]] = "geometry" if (rule.get("metrics") or rule.get("findings")) else "nothing to check"
    for c in ("impedance", "crosstalk", "si", "vias", "timing", "thermal", "emc"):
        d = rep["checks"].get(c)
        if d is None:
            continue
        fam, req = CHECK_FAMILY[c]
        fs += from_rule_findings(c, fam, req, d.get("findings", []))
        val[c] = d.get("validation") or ("rule of thumb" if c == "emc" else "unvalidated")
        n_items = len(d.get("rows") or d.get("parts") or []) + len(d.get("channels") or [])
        if not n_items and not d.get("findings"):
            val[c] = "nothing to check"
        elif c == "impedance" and d.get("judged") == 0:
            val[c] = "not judged"           # impedances computed, but no class gives a target
    if "dc" in rep["checks"]:
        f2, v = judge_dc(rep, cfg)
        fs += f2
        val["dc"] = v
    if "pdn" in rep["checks"]:
        f2, v = judge_pdn(rep, cfg)
        fs += f2
        val["pdn"] = v
    seen = {}
    for f in fs:  # same cause twice in one run (e.g. two segments at one spot): keep IDs unique
        n = seen.get(f["id"], 0)
        seen[f["id"]] = n + 1
        if n:
            f["id"] = f"{f['id']}.{n}"
    return fs, val


def _validated(v):
    """Geometry and rule checks need no solver case; solver checks need a passed G17 case (or KiPI's
    upstream suite)."""
    return isinstance(v, str) and (v in ("geometry", "rule-based", "closed-form") or v.startswith(("validated", "KiPI")))


def scorecard(findings, val, skipped, assumed=()):
    per = {}
    for f in findings:
        c = per.setdefault(f["check"], {"check": f["check"], "family": f["family"], "requirement": f["requirement"],
                                        "fail": 0, "warn": 0, "info": 0, "pass": 0, "in_area": 0, "waived": 0, "pending": 0})
        if f["waived"]:
            c["waived"] += 1
        elif f["area"]:
            c["in_area"] += 1
        else:
            c[f["severity"]] += 1
    for chk, v in val.items():
        fam, req = RULE_FAMILY.get(chk) or CHECK_FAMILY.get(chk, ("G1", chk))
        per.setdefault(chk, {"check": chk, "family": fam, "requirement": req, "fail": 0, "warn": 0, "info": 0,
                             "pass": 0, "in_area": 0, "waived": 0})
    for c in per.values():
        v = val.get(c["check"])
        c["validation"] = v
        c["status"] = ("fail" if c["fail"] else "warn" if c["warn"] else "pending" if c.get("pending") else
                       v if v in ("nothing to check", "not judged") else
                       "unvalidated" if not _validated(v) else
                       "assumed inputs" if c["check"] in assumed else "pass")
        if c["check"] in assumed:
            c["why"] = assumed[c["check"]]
    for chk, why in skipped.items():
        fam, req = RULE_FAMILY.get(chk) or CHECK_FAMILY.get(chk, ("G1", chk))
        per[chk] = {"check": chk, "family": fam, "requirement": req, "status": "not run", "why": why,
                    "fail": 0, "warn": 0, "info": 0, "pass": 0, "in_area": 0, "waived": 0}
    for fam, req, what in NOT_BUILT:
        per[f"{fam}:{req}"] = {"check": req, "family": fam, "requirement": req, "status": "not built", "why": what,
                               "fail": 0, "warn": 0, "info": 0, "pass": 0, "in_area": 0, "waived": 0}
    fams = {}
    for c in per.values():
        fams.setdefault(c["family"], []).append(c)
    card = []
    for fam, cs in fams.items():
        real = [c for c in cs if c["status"] not in NOT_REAL]
        if real:
            status = max((c["status"] for c in real), key=ORDER.get)
        else:
            sts = {c["status"] for c in cs}
            status = "not built" if sts == {"not built"} else "nothing to check" if "nothing to check" in sts else "not run"
        card.append({"family": fam, "title": FAMILY_TITLE.get(fam, fam), "status": status,
                     "partial": bool(real) and len(real) < len(cs), "checks": len(cs),
                     "fail": sum(c["fail"] for c in cs), "warn": sum(c["warn"] for c in cs)})
    fams = {F["family"]: F for F in card}
    fo = lambda fam: FAMILY_ORDER.index(fam) if fam in FAMILY_ORDER else 99
    return (sorted(per.values(), key=lambda c: (fo(c["family"]), c["check"])),
            sorted(fams.values(), key=lambda F: fo(F["family"])))


def exit_code(findings, cfg):
    live = [f for f in findings if not f["waived"] and not f["area"]]
    if any(f["severity"] == "fail" for f in live):
        return 2
    if any(f["severity"] == "warn" for f in live) and cfg.data.get("evaluate", {}).get("warnings_exit", True):
        return 1
    return 0


# ------------------------------------------------------------------------------------------------
def since_report(prev_path, findings):
    """#35: findings against an earlier evaluation, by stable ID (fail, warn and pending outside areas)."""
    prev = json.load(open(prev_path))
    live = lambda fs: {f["id"]: f for f in fs if f["severity"] in ("fail", "warn", "pending") and not f.get("area")
                       and not f.get("waived")}
    a, b = live(prev.get("findings", [])), live(findings)
    changed = [{"id": i, "was": a[i]["severity"], "now": b[i]["severity"], "value_was": a[i].get("value"),
                "value_now": b[i].get("value")} for i in sorted(set(a) & set(b))
               if a[i]["severity"] != b[i]["severity"] or a[i].get("value") != b[i].get("value")]
    pw = {f["id"]: f["waived"] for f in prev.get("findings", []) if f.get("waived")}
    now_ids = {f["id"] for f in findings}
    return {"against": os.path.abspath(prev_path), "against_board_sha256": prev.get("board_sha256"),
            "against_date": prev.get("date"),
            "new": [{"id": i, "check": b[i]["check"], "severity": b[i]["severity"], "reason": b[i]["reason"]} for i in sorted(set(b) - set(a))],
            "resolved": [{"id": i, "check": a[i]["check"], "severity": a[i]["severity"], "reason": a[i]["reason"]} for i in sorted(set(a) - set(b))],
            "changed": changed,
            "waived_before": {"still_present": sorted(i for i in pw if i in now_ids), "gone": sorted(i for i in pw if i not in now_ids)}}


def apply_stage(findings, st):
    """#29: before the final stage, 'not connected yet' is pending, not a failure."""
    if st["stage"] == "final":
        return 0
    inc = set(st["incomplete_nets"])
    n = 0
    for f in findings:
        if f["severity"] not in ("fail", "warn"):
            continue
        r = f["reason"]
        if (f["check"] == "kicad_drc" and r.startswith("unconnected_items")) or \
           (f["check"] in ("dc", "pdn") and ("not connected" in r or "no load pad" in r)) or \
           (f["check"] in ("stubs",) and set(f["nets"]) & inc):
            f["severity"], f["pending_reason"] = "pending", f"layout stage '{st['stage']}': not routed yet"
            n += 1
    return n


def run(board, cfg_path=None, waivers=None, families=None, checks=None, out_dir=None, json_out=None,
        markers=True, drc=True, log=print, stage=None, since=None, em="none"):
    from . import report as report_mod
    from . import run as runmod
    from . import validate
    fams = [x.strip().lower() for x in families.split(",")] if families else None
    if checks is None:
        checks = sorted({c for f in (fams or FAMILY_CHECKS) for c in FAMILY_CHECKS.get(f, [])}) if fams else None
    if fams and set(fams) & {"pl", "fs"}:
        log(f"evaluate: families {sorted(set(fams) & {'pl', 'fs'})} are not built yet; reported as 'not built'")
    rep, t, od = runmod.run_checks(board, cfg_path, checks, out_dir, log=log)
    b, cfg, ctx = t
    report_mod.write_all(rep, t, od, markers=markers, log=log)
    want = runmod.ALL_CHECKS if checks is None else checks
    skipped = {}
    for c in want:
        if c in ("dc", "pdn") and not cfg["rails"]:
            skipped[c] = "no rails (config or detection)"
        elif c == "pdn" and not any("ripple_pct" in r for r in cfg["rails"].values()):
            skipped[c] = "no rail has ripple_pct / step_a"
        elif c == "thermal" and not cfg["thermal"]["parts"]:
            skipped[c] = "no part dissipations in [thermal.parts]"
    assumed = {}
    am = [k for k, r in cfg["rails"].items() if r.get("current_assumed")]
    if am:
        assumed["dc"] = f"currents assumed ({len(am)} of {len(cfg['rails'])} rails): judged at most warn"
        assumed["pdn"] = assumed["dc"]
    findings, val = collect(rep, cfg, want)
    drc_meta = None
    if drc and (fams is None or "ev" in fams):
        f2, drc_meta = kicad_drc(os.path.abspath(board), od, cfg, log=log)
        findings += f2
        if drc_meta.get("ran"):
            val["kicad_drc"] = "geometry"
    elif not drc:
        skipped["kicad_drc"] = "--no-drc"
    em_jobs = None
    if em in ("plan", "run"):                    # #19: automatic 3D regions with fix variants
        from . import emplan, region3d
        from .si import inside, safe_name
        em_jobs = emplan.plan(ctx, rep, cfg.data.get("em", {}).get("fast_ghz", 3.0), log=log)
        jd = os.path.join(od, "em_jobs")
        os.makedirs(jd, exist_ok=True)
        for j in em_jobs:
            # job names come from net and part names ("/sheet/NET"): a file name only through safe_name
            with open(inside(jd, os.path.join(jd, safe_name(j["name"]) + ".json")), "w") as fh:
                json.dump(j, fh, indent=1)
        nfit = sum(1 for j in em_jobs if j.get("fit", {}).get("ok"))
        log(f"3D regions: {len(em_jobs)} planned, {nfit} with a box that passes pre-flight"
            + (" (running: minutes each)" if em == "run" else ""))
        for j in em_jobs:
            if not j.get("fit", {}).get("ok"):
                msg = f"{j['name']}: no box passes the 3D pre-flight: {j.get('fit', {}).get('why', '')}"
                log("  " + msg)
                if em == "plan":
                    findings.append(_f("em_region", "FS", "FS9", "info", msg, j["nets"], j.get("parts", []), extra=j["name"]))
        if em == "run":
            for j in em_jobs:
                spec = emplan.region_spec(j)
                base = region3d.run(ctx, spec, runmod.CACHE, runmod.RUNS, log)
                j["result"] = {k: base.get(k) for k in ("dir", "ran", "convergence", "calibration", "z_line_ohm", "preflight")}
                if j.get("variant"):
                    var = region3d.run(ctx, emplan.region_spec(j, variant=True), runmod.CACHE, runmod.RUNS, log)
                    j["variant_result"] = {k: var.get(k) for k in ("dir", "ran", "convergence", "z_line_ohm")}
                for sev, msg in emplan.judge(j, base):
                    findings.append(_f("em_region", "FS", "FS9", sev, msg, j["nets"], j.get("parts", []), extra=j["name"]))
    from . import stage as stage_mod
    st = stage_mod.analyse(b.snap)
    if stage:
        st["stage"], st["stage_set_by"] = stage, "--stage"
    st["pending"] = apply_stage(findings, st)
    ws, whash = load_waivers(waivers)
    wstat = apply_waivers(findings, ws)
    per, famcard = scorecard(findings, val, skipped, assumed)
    code = exit_code(findings, cfg)
    ev = {
        "schema": SCHEMA, "tool": "openpcbsi", "version": rep.get("version"),
        "date": datetime.datetime.now().isoformat(timespec="seconds"),
        "board": rep["board"], "board_sha256": rep["board_sha256"], "downgraded": rep.get("downgraded"),
        "kicad_file_version": rep.get("kicad_file_version"), "kicad": rep.get("kicad"),
        "config": rep["config"], "config_hash": rep["config_hash"], "detected": rep.get("detected"),
        "waivers_file": os.path.abspath(waivers) if waivers else None, "waivers_hash": whash,
        "solver_versions": dict(validate.SOLVER_VERSIONS), "kipi": (rep.get("kipi") or {}).get("lock"),
        "kicad_drc": drc_meta, "exit_code": code, "families": famcard, "checks": per, "waivers": wstat,
        "stage": st["stage"], "incomplete_nets": st["incomplete_nets"][:200], "unrouted_pads": st["unrouted_pads"],
        "since": since_report(since, findings) if since else None,
        "em_jobs": em_jobs,
        "counts": {s: sum(1 for f in findings if f["severity"] == s and not f["waived"] and not f["area"])
                   for s in ("fail", "warn", "info", "pass", "pending")} | {
                      "in_area": sum(1 for f in findings if f["area"] and not f["waived"]),
                      "waived": sum(1 for f in findings if f["waived"])},
        "findings": findings, "report": os.path.join(od, "report.json"),
    }
    path = os.path.join(od, "evaluation.json")
    from .paths import portable
    plain = json.loads(json.dumps(ev, default=str))
    with open(path, "w") as fh:
        json.dump(portable(plain, od), fh, indent=1)
    if json_out:
        with open(json_out, "w") as fh:
            json.dump(portable(plain, os.path.dirname(os.path.abspath(json_out))), fh, indent=1)
    write_html(portable(plain, od), os.path.join(od, "evaluation.html"))
    log("")
    log(f"scorecard ({'config ' + rep['config'] if cfg.path else 'config detected from the board'}):")
    for F in famcard:
        log(f"  {F['status'] + (' (partial)' if F['partial'] else ''):22s} {F['title']:38s} fail {F['fail']:4d}  warn {F['warn']:4d}")
    c = ev["counts"]
    log(f"  {c['fail']} fail, {c['warn']} warn, {c['info']} info; {c['in_area']} in exclusion areas, {c['waived']} waived"
        + (f"; waivers: {sum(1 for w in wstat if w['status'] == 'stale')} stale, "
           f"{sum(1 for w in wstat if w['status'].startswith('invalid'))} invalid" if wstat else ""))
    log(f"evaluation: {path}")
    log(f"html:       {os.path.join(od, 'evaluation.html')}")
    log(f"exit code {code}")
    return code, ev


# ------------------------------------------------------------------------------------------------
CSS = """
:root{--bg:#fff;--fg:#1d2330;--mut:#5d6675;--line:#d9dee7;--fail:#b3261e;--warn:#9a6700;--pass:#1e7a3d;--un:#5b5fc7}
@media (prefers-color-scheme:dark){:root{--bg:#14171c;--fg:#e4e7ec;--mut:#9aa3b2;--line:#2c323c;--fail:#ff8a80;--warn:#f2c14e;--pass:#7ad19a;--un:#a9adff}}
body{background:var(--bg);color:var(--fg);font:14px/1.45 system-ui,sans-serif;margin:24px auto;max-width:1200px;padding:0 16px}
h1{font-size:20px}h2{font-size:16px;margin-top:28px}table{border-collapse:collapse;width:100%;margin:8px 0}
td,th{border-bottom:1px solid var(--line);padding:4px 8px;text-align:left;vertical-align:top}th{color:var(--mut);font-weight:600}
.fail{color:var(--fail);font-weight:600}.warn{color:var(--warn);font-weight:600}.pass{color:var(--pass)}
.unvalidated,.notrun,.notbuilt{color:var(--un)}.mut{color:var(--mut)}td.n{text-align:right;font-variant-numeric:tabular-nums}
code{font-size:12px}
"""


def write_html(ev, path):
    e = html.escape
    st = lambda s: f'<span class="{s.replace(" ", "")}">{e(s)}</span>'
    H = [f"<!doctype html><html><head><meta charset=utf-8><meta name=viewport content='width=device-width,initial-scale=1'>"
         f"<title>openpcbsi evaluation</title><style>{CSS}</style></head><body>"]
    H.append(f"<h1>openpcbsi evaluation: {e(os.path.basename(ev['board']))}</h1>")
    H.append(f"<p class=mut>{e(ev['board'])}<br>sha256 {e(ev['board_sha256'][:16])} · config {e(str(ev['config']))} "
             f"(hash {e(ev['config_hash'])}) · {e(ev['date'])} · exit code {ev['exit_code']} · "
             f"<a href='report.html'>full report</a></p>")
    if ev.get("detected"):
        d = ev["detected"]
        H.append("<h2>Detected from the board (no config file)</h2><ul>")
        H += [f"<li>{e(a)}</li>" for a in d.get("assumptions", [])]
        H.append("</ul><p class=mut>Left out:</p><ul>")
        H += [f"<li class=mut>{e(a)}</li>" for a in d.get("left_out", [])]
        H.append("</ul>")
    H.append("<h2>Scorecard</h2><table><tr><th>Family</th><th>Status</th><th class=n>fail</th><th class=n>warn</th><th class=n>checks</th></tr>")
    for F in ev["families"]:
        H.append(f"<tr><td>{e(F['title'])}</td><td>{st(F['status'])}{' <span class=mut>(partial: some checks not built)</span>' if F['partial'] else ''}</td><td class=n>{F['fail']}</td><td class=n>{F['warn']}</td><td class=n>{F['checks']}</td></tr>")
    H.append("</table><h2>Checks</h2><table><tr><th>Family</th><th>Check</th><th>Req.</th><th>Status</th><th class=n>fail</th>"
             "<th class=n>warn</th><th class=n>info</th><th class=n>in areas</th><th class=n>waived</th><th>Validation / note</th></tr>")
    for c in ev["checks"]:
        H.append(f"<tr><td>{e(c['family'])}</td><td>{e(c['check'])}</td><td>{e(c['requirement'])}</td><td>{st(c['status'])}</td>"
                 f"<td class=n>{c['fail']}</td><td class=n>{c['warn']}</td><td class=n>{c['info']}</td><td class=n>{c['in_area']}</td>"
                 f"<td class=n>{c['waived']}</td><td class=mut>{e(str(c.get('why') or c.get('validation') or ''))}</td></tr>")
    H.append("</table>")
    if ev["waivers"]:
        H.append("<h2>Waivers</h2><table><tr><th>Waiver</th><th>Status</th><th>Reason</th><th>Reviewer</th></tr>")
        for w in ev["waivers"]:
            H.append(f"<tr><td><code>{e(str(w.get('id') or w.get('check')))}</code></td><td>{e(w['status'])}</td>"
                     f"<td>{e(str(w.get('reason', '')))}</td><td>{e(str(w.get('reviewer', '')))}</td></tr>")
        H.append("</table>")
    rank = {"fail": 0, "warn": 1}
    live = sorted((f for f in ev["findings"] if f["severity"] in rank),
                  key=lambda f: (bool(f["waived"]), bool(f["area"]), rank[f["severity"]], f["family"], f["check"]))
    H.append(f"<h2>Findings: fail and warn ({len(live)})</h2><table><tr><th>ID</th><th>Sev.</th><th>Check</th><th>Reason</th>"
             "<th>Nets / refs</th><th>Layer</th><th>x, y (mm)</th><th>Note</th></tr>")
    for f in live[:2000]:
        note = "waived: " + f["waived"]["reason"] if f["waived"] else ("in area " + f["area"] if f["area"] else "")
        xy = f"{f['x']:.2f}, {f['y']:.2f}" if f["x"] is not None and f["y"] is not None else ""
        H.append(f"<tr><td><code>{e(f['id'])}</code></td><td>{st(f['severity'])}</td><td>{e(f['check'])} <span class=mut>{e(f['requirement'])}</span></td>"
                 f"<td>{e(f['reason'])}</td><td>{e(', '.join(f['nets'] + f['refs']))}</td><td>{e(f['layer'] or '')}</td>"
                 f"<td>{xy}</td><td class=mut>{e(note)}</td></tr>")
    H.append("</table><p class=mut>Info findings are in evaluation.json and report.html.</p></body></html>")
    with open(path, "w") as fh:
        fh.write("\n".join(H))
