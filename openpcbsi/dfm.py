"""Manufacturing review (REQUIREMENTS 5.17): `openpcbsi dfm BOARD`.

    python3 -m openpcbsi dfm BOARD [-c CONFIG] [--fab jlcpcb | FILE] [--gerbers DIR] [--checks ...]
                                   [--waivers FILE] [--observed FILE] [--json OUT] [-o OUT_DIR] [--no-markers]

This is a review, not a sign-off: it cannot see what the fabricator's CAM does to the files, and it
lists what it does not check. It has two kinds of output and keeps them apart.

RULES judge geometry against a limit that the fabricator publishes (the profile), that an industry
standard gives, or that the board config sets. Only rules give a verdict and an exit code. They
apply everywhere on the board: there are no exempt areas, and where the fabricator publishes a
different figure for a BGA field that figure is used there.

  MF19 rule_audit     KiCad's own design rules against the fabricator's limits (is a clean DRC evidence?)
  MF20 output_match   the Gerber and drill files against the geometry the checks worked on
  MF3  fab_track      track width
  MF12 fab_space      copper spacing between nets, measured (track to track, pad to track)
  MF2  same_net_gap   gaps beside a track narrower than the same-net track spacing
  MF13 fab_drill      hole sizes, annular rings, hole to hole, hole to copper, copper to board edge
  MF4  mask_web       solder-mask webs and pad spacing, as drawn and at the fabricator's own expansion
  MF5  mask_exposed   copper of another net inside or beside a mask opening
  MF14 via_pad        vias in pads, untented vias beside pads
  MF15 paste          stencil apertures: area ratio and aspect ratio (IPC-7525), minimum aperture
  MF16 footprint      two-pad land patterns against the body of the chip they are named for
  MF17 part_spacing   spacing between parts against the assembler's table
  MF18 assembly       smallest package, pin pitch, sides populated, fiducials

RANKINGS order places on the board by a measurement or a model that no published limit covers and
that has not been compared with fabricated boards. They never pass, never fail and do not change
the exit code, unless `[dfm.judge]` gives a limit and says where it comes from.

  MF1  acid_trap      acute corners in copper of one net, thinnest track first
  MF21 copper         copper share of each layer and of each pair of layers about the centre
  MF7  reflow         how far each part's modelled peak lies from the board median, with its range
  MF8  chips          modelled time between the two joints of a chip melting, with the land geometry

Writes <run>/dfm/dfm.json (schema `openpcbsi-dfm/2`), dfm.html, one image per copper layer, the
files that were compared (fab_files/) and a marker board.
"""
import datetime
import html
import json
import hashlib
import math
import re
import os
import subprocess
import time
import tomllib

import numpy as np
from scipy import ndimage
from scipy.spatial import cKDTree

from . import padgeom, fab as fab_mod
from . import geom, packages
from .drc import finding, result

SCHEMA = "openpcbsi-dfm/2"
GROUPS = {"audit": ["rule_audit"], "output": ["output_match"], "copper": ["fab_track", "fab_space", "same_net_gap", "acid_trap", "copper"],
          "drill": ["fab_drill"], "mask": ["mask_web", "mask_exposed", "via_pad"], "paste": ["paste"],
          "parts": ["footprint", "part_spacing", "assembly"], "reflow": ["reflow"], "chips": ["chips"]}
REQ = {"rule_audit": "MF19", "output_match": "MF20", "fab_track": "MF3", "fab_space": "MF12", "same_net_gap": "MF2",
       "fab_drill": "MF13", "mask_web": "MF4", "mask_exposed": "MF5", "via_pad": "MF14", "paste": "MF15", "footprint": "MF16",
       "part_spacing": "MF17", "assembly": "MF18", "acid_trap": "MF1", "copper": "MF21", "reflow": "MF7", "chips": "MF8"}
RANKINGS = ("acid_trap", "copper", "reflow", "chips")
ORDER = ["rule_audit", "output_match", "fab_track", "fab_space", "same_net_gap", "fab_drill", "mask_web", "mask_exposed", "via_pad",
         "paste", "footprint", "part_spacing", "assembly", "acid_trap", "copper", "reflow", "chips"]
NOT_CHECKED = [
    "what the fabricator's CAM does to the files: etch compensation, changes to mask openings, panelisation",
    "the etched cross-section (etch factor) and finished track width (the published tolerance is recorded, not judged)",
    "drill aspect ratio, copper balance and warpage against a limit (no published figure; the copper shares are listed)",
    "part height, shadowing in reflow and the second pass of a double-sided board (KiCad carries no part heights)",
    "a land pattern against its part's data sheet (only two-pad chips, against the nominal body of their EIA size)",
    "stencil aperture shapes and paste volume beyond the IPC-7525 ratios",
    "silkscreen, test access, panel rails and tooling holes",
    "clearance, annular rings and holes as KiCad's DRC sees them: MF19 says whether its rules match the fabricator",
]
TOL = 0.003          # mm: polygons approximate arcs to 2 um
FINE_EDGE = 0.0025   # mm: the fine raster of the hole-to-copper and copper-to-edge measurements


class Agg:
    """One finding per kind (and layer or part): the worst value, with the count of places."""

    def __init__(self, rule):
        self.rule, self.d = rule, {}

    def add(self, key, sev, val, msg, x=None, y=None, n=1, **kw):
        e = self.d.setdefault(key, {"n": 0, "val": None})
        e["n"] += n
        if e["val"] is None or val < e["val"]:
            e.update(val=val, sev=sev, msg=msg, x=x, y=y, kw=kw)

    def findings(self, what="places"):
        out = []
        for key, e in self.d.items():
            kw = dict(e["kw"])
            more = f" ({e['n']} {what}; worst shown)" if e["n"] > 1 else ""
            out.append(finding(self.rule, e["sev"], e["msg"] + more, kw.pop("net", ""), kw.pop("layer", ""), e["x"], e["y"],
                               value=None if e["val"] is None else round(float(e["val"]), 4), count=e["n"],
                               kind=kw.pop("kind", str(key)), **kw))
        out.sort(key=lambda f: (f["severity"] != "error", f["severity"] != "warning"))
        return out


def _rings(polys):
    return [ring for poly in polys for ring in poly]


def _bbox(rings):
    a = np.vstack([np.asarray(r, dtype=float) for r in rings])
    return a[:, 0].min(), a[:, 1].min(), a[:, 0].max(), a[:, 1].max()


def _pt_rings(x, y, rings):
    """Distance from a point to a polygon (0 inside)."""
    if geom.points_in_rings([x], [y], rings)[0]:
        return 0.0
    s = geom.ring_segments(rings)
    return float(geom.point_seg_dist(x, y, s[:, 0], s[:, 1], s[:, 2], s[:, 3]).min())


def _mask_pads(b, mname):
    out = []
    for p in b.pads:
        m = p.get("mask", {}).get(mname)
        if not m or not p["poly"]:
            continue
        rings = _rings(m["polys"])
        if not rings:
            continue
        out.append({"pad": p, "mask": rings, "cu": _rings(p["poly"]), "bb": _bbox(rings), "exp": m["exp"]})
    return out


def _side(b, ref):
    """The side a part's SMD pads are on (the footprint's own flag can disagree with its pads)."""
    last = len(b.layers) - 1
    ls = {l for p in b.pads_by_ref.get(ref, []) if p["kind"] == "smd" for l in p["layer_ids"]}
    return "B" if ls == {last} and last > 0 else "F"


def bga_field(ctx):
    """Raster of the BGA fields (the pads of each BGA, grown by 0.3 mm) and their references."""
    b = ctx.b
    g = b.grid()
    refs = [ref for ref, fp in b.fps.items() if packages.part_class(fp) == "BGA" and len(b.pads_by_ref.get(ref, [])) >= 16]
    refs += [r for r in ctx.cfg["dfm"]["bga_refs"] if r in b.fps and r not in refs]
    m = np.zeros(g.shape, bool)
    for ref in refs:
        pts = [(p["x"], p["y"]) for p in b.pads_by_ref.get(ref, [])]
        if not pts:
            continue
        a = np.array(pts)
        x0, y0, x1, y1 = a[:, 0].min() - 0.6, a[:, 1].min() - 0.6, a[:, 0].max() + 0.6, a[:, 1].max() + 0.6
        geom.raster_rings(g, [[(x0, y0), (x1, y0), (x1, y1), (x0, y1)]], m)
    return m, refs


# ------------------------------------------------------------------------------------------------
# MF19: KiCad's own rules against the fabricator's

AUDIT = [   # board-setup rule, DRC test, every profile limit of the same measurand, what, the dfm check that measures it
    # one KiCad number cannot hold every figure: SMD pad to pad (0.15 mm at JLCPCB) would forbid finer tracks,
    # and the non-plated hole minimum would forbid vias; those are measured by mask_web and fab_drill only
    ("min_clearance", "clearance", ["track_space_mm", "pad_to_track_mm"], "copper clearance", "fab_space"),
    ("min_track_width", "track_width", ["track_width_mm"], "track width", "fab_track"),
    ("min_via_annular_width", "annular_width", ["via_ring_min_mm"], "via annular ring", "fab_drill"),
    ("min_via_diameter", None, ["via_diameter_min_mm"], "via diameter", "fab_drill"),
    ("min_through_hole_diameter", "drill_out_of_range", ["via_hole_min_mm"], "hole diameter", "fab_drill"),
    ("min_hole_to_hole", "hole_to_hole", ["via_hole_spacing_mm", "pad_hole_spacing_mm"], "hole to hole", "fab_drill"),
    ("min_hole_clearance", "hole_clearance", ["via_hole_to_copper_mm", "pth_hole_to_copper_mm", "inner_pth_hole_to_copper_mm",
                                              "npth_hole_to_copper_mm"], "hole to copper", "fab_drill"),
    ("min_copper_edge_clearance", "copper_edge_clearance", ["copper_to_edge_mm"], "copper to board edge", "fab_drill"),
    ("solder_mask_min_width", "solder_mask_bridge", ["mask_web_mm"], "solder-mask web", "mask_web"),
]
# .kicad_dru constraint type -> the profile limits of the same measurand
DRU = {"clearance": ["track_space_mm", "pad_to_track_mm", "smd_pad_to_pad_mm"], "track_width": ["track_width_mm"],
       "hole_clearance": ["via_hole_to_copper_mm", "pth_hole_to_copper_mm", "inner_pth_hole_to_copper_mm", "npth_hole_to_copper_mm"],
       "hole_to_hole": ["via_hole_spacing_mm", "pad_hole_spacing_mm"], "edge_clearance": ["copper_to_edge_mm"],
       "annular_width": ["via_ring_min_mm", "pth_ring_min_mm"], "hole_size": ["via_hole_min_mm"], "via_diameter": ["via_diameter_min_mm"]}
KEY_WHAT = {"track_space_mm": "track to track", "pad_to_track_mm": "pad to track", "smd_pad_to_pad_mm": "SMD pad to pad",
            "via_hole_to_copper_mm": "via holes", "pth_hole_to_copper_mm": "plated holes (outer layers)",
            "inner_pth_hole_to_copper_mm": "plated holes (inner layers)", "npth_hole_to_copper_mm": "non-plated holes",
            "via_hole_min_mm": "via holes", "npth_min_mm": "non-plated holes", "via_hole_spacing_mm": "via holes",
            "pad_hole_spacing_mm": "pad holes", "via_ring_min_mm": "vias", "pth_ring_min_mm": "plated pad holes"}


def _compare(fabr, val, keys):
    """KiCad's value against every profile limit of the same measurand: (as strict for [...], looser for
    [(what, limit)], measurand notes)."""
    strict, looser, notes = [], [], []
    for k in keys:
        if not fabr.judged(k):
            continue
        lim = fabr.get(k)
        what = KEY_WHAT.get(k, k)
        meas = fabr.limits[k].get("measurand", "")
        if meas == "hole_centre_to_centre" and val < lim - 1e-6:
            # KiCad measures wall to wall: a smaller wall-to-wall figure may still meet a centre-to-centre one,
            # depending on the hole sizes; it is not called looser
            notes.append(f"KiCad's {val:g} mm (wall to wall) against {fabr.name}'s {lim:g} mm for {what}, read centre to centre: "
                         f"it holds for holes of {lim - val:g} mm and larger; `dfm` measures it")
            continue
        (strict if val >= lim - 1e-6 else looser).append((what, lim, k))
    return strict, looser, notes


def _dru_items(path):
    from . import project as pj
    try:
        items = pj._parse(pj._tokens(open(path).read()))
    except (OSError, ValueError):
        return []
    out = []
    for it in items:
        if not (isinstance(it, list) and it and it[0] == "rule"):
            continue
        name = it[1] if len(it) > 1 and isinstance(it[1], str) else "?"
        cond = next((x[1] for x in it[2:] if isinstance(x, list) and x and x[0] == "condition" and len(x) > 1), None)
        for x in it[2:]:
            if isinstance(x, list) and x and x[0] == "constraint" and len(x) > 1:
                mn = next((pj._num(y[1]) for y in x[2:] if isinstance(y, list) and len(y) == 2 and y[0] == "min"), None)
                out.append({"rule": name, "condition": cond, "type": x[1], "min": mn})
    return out


def order_checklist(ctx, fabr, rules):
    """What the order must state for the verdict to hold: options the fabricator otherwise chooses."""
    cfg, out = ctx.cfg, []
    vp = rules.get("via_pad")
    if vp:
        n_in, cover = vp["metrics"]["vias"]["in_pads"], vp["params"]["via_cover"]
        note = fabr.limits.get("filled_order_note", {})
        if n_in:
            out.append({"item": "via covering", "state": f"{n_in} vias are in SMD pads: order them epoxy filled and capped"
                        + (f" and say which ({fabr.name}: \"{note.get('quote')}\")" if note.get("quote") else ""),
                        "set_in_config": bool(cfg["dfm"]["via_cover"]), "basis": vp["params"]["via_cover_source"]})
        else:
            out.append({"item": "via covering", "state": f"{cover}", "set_in_config": bool(cfg["dfm"]["via_cover"]),
                        "basis": vp["params"]["via_cover_source"]})
    pa = rules.get("paste")
    if pa and pa["params"].get("stencil_mm"):
        out.append({"item": "stencil thickness", "state": f"{pa['params']['stencil_mm']:g} mm"
                    + ("" if cfg["dfm"]["stencil_mm"] else f" (predicted; {fabr.name} chooses unless the order states it)"),
                    "set_in_config": bool(cfg["dfm"]["stencil_mm"]), "basis": pa["params"]["stencil_source"]})
    if fabr.get("rails_required"):
        out.append({"item": "edge rails", "state": f"required for the {cfg['dfm'].get('assembly_service', 'standard')} assembly service: "
                    f"add them to the panel or order them", "set_in_config": False, "basis": fabr.src("rails_required")})
    out.append({"item": "solder-mask colour", "state": str(cfg["dfm"].get("mask_color", "green")), "set_in_config": True,
                "basis": "board config (the pad-spacing figure depends on it)"})
    return out


def rule_audit(ctx, fabr):
    b = ctx.b
    stem = os.path.splitext(b.path)[0]
    pro = stem + ".kicad_pro"
    fs, rows = [], []
    if not os.path.exists(pro):
        fs.append(finding("rule_audit", "warning", f"no KiCad project file beside the board: KiCad's DRC judges clearance, annular rings "
                          f"and holes against KiCad's default rules, not against {fabr.name}'s limits", kind="no_project"))
        return dict(result("rule_audit", "MF19", "KiCad's design rules against the fabricator's limits", {}, {}, fs), rows=rows, checked=1)
    try:
        ds = json.load(open(pro)).get("board", {}).get("design_settings", {})
    except (ValueError, OSError) as e:
        fs.append(finding("rule_audit", "warning", f"the KiCad project file could not be read ({e})", kind="unreadable"))
        return dict(result("rule_audit", "MF19", "KiCad's design rules against the fabricator's limits", {}, {}, fs), rows=rows, checked=1)
    rules, sev = ds.get("rules", {}), ds.get("rule_severities", {})
    unset = []
    for key, test, keys, what, by in AUDIT:
        keys = [k for k in keys if fabr.judged(k)]
        if not keys:
            continue
        val = rules.get(key)
        s = sev.get(test) if test else None
        row = {"rule": key, "kicad_mm": val, "fab": {k: fabr.get(k) for k in keys}, "drc_test": test, "severity": s, "measured_by": by}
        if val is None:
            row["verdict"] = "not set in the project"
            unset.append(what)
        else:
            strict, looser, notes = _compare(fabr, val, keys)
            parts = []
            if strict:
                parts.append("as strict as " + fabr.name + " for " + ", ".join(f"{w} ({l:g})" for w, l, _ in strict))
            if looser:
                parts.append("looser for " + ", ".join(f"{w} ({val:g} < {l:g})" for w, l, _ in looser))
                fs.append(finding("rule_audit", "warning", f"KiCad's DRC accepts {what} down to {val:g} mm; {fabr.name} needs "
                                  + "; ".join(f"{l:g} mm for {w} ({fabr.src(k)})" for w, l, k in looser)
                                  + f". A clean DRC does not show these limits are met; `dfm` measures them in {by}",
                                  value=val, limit=max(l for _, l, _ in looser), unit="mm", kind=key, basis=fabr.src(looser[0][2])))
            row["verdict"] = "; ".join(parts + notes)
        if s == "ignore":
            row["verdict"] += "; test switched off"
            fs.append(finding("rule_audit", "warning", f"KiCad's DRC test '{test}' ({what}) is set to ignore in the project: "
                              f"DRC does not check it at all; `dfm` measures it in {by}", kind=key + ":ignore", basis="KiCad project"))
        rows.append(row)
    if unset:
        fs.append(finding("rule_audit", "warning", f"the KiCad project sets no rule for {', '.join(unset)}: KiCad's DRC uses its defaults "
                          f"there, not {fabr.name}'s limits; `dfm` measures them itself", kind="unset"))
    dru = stem + ".kicad_dru"
    if os.path.exists(dru):
        n = 0
        for c in _dru_items(dru):
            keys = [k for k in DRU.get(c["type"], []) if fabr.judged(k)]
            if c["min"] is None or not keys:
                continue
            n += 1
            _, looser, _ = _compare(fabr, c["min"], keys)
            row = {"rule": f"{c['rule']} ({c['type']})", "kicad_mm": c["min"], "fab": {k: fabr.get(k) for k in keys}, "drc_test": c["type"],
                   "severity": None, "measured_by": "custom rule", "condition": c["condition"],
                   "verdict": ("looser for " + ", ".join(f"{w} ({c['min']:g} < {l:g})" for w, l, _ in looser)) if looser else "as strict"}
            rows.append(row)
            if looser:
                fs.append(finding("rule_audit", "warning", f"custom rule '{c['rule']}' (.kicad_dru{', where ' + c['condition'] if c['condition'] else ''}) "
                                  f"lets KiCad's DRC accept {c['type'].replace('_', ' ')} down to {c['min']:g} mm; {fabr.name} needs "
                                  + "; ".join(f"{l:g} mm for {w}" for w, l, _ in looser), value=c["min"],
                                  limit=max(l for _, l, _ in looser), unit="mm", kind="dru:" + c["rule"], basis=fabr.src(looser[0][2])))
        if not n:
            fs.append(finding("rule_audit", "info", "the project has custom rules (.kicad_dru) with no minimum of a measurand the profile "
                              "covers", kind="custom_rules"))
    return dict(result("rule_audit", "MF19", "KiCad's design rules against the fabricator's limits", {"project": pro},
                       {r["rule"]: {"kicad_mm": r["kicad_mm"], "fab": r["fab"]} for r in rows}, fs), rows=rows, checked=len(rows))


# ------------------------------------------------------------------------------------------------
# MF3

def rule_fab_track(ctx, fabr, bga):
    b = ctx.b
    lim = fabr.judge("track_width_mm", "fab_track")
    local = fabr.judge("bga_track_mm", "fab_track")
    T = b.tracks
    agg, metrics = Agg("fab_track"), {}
    if lim is None or not len(T["w"]):
        return dict(result("fab_track", "MF3", "Track width against the fabricator's minimum", {"limit_mm": lim}, {}, []), checked=0)
    g = b.grid()
    mx, my = (T["x1"] + T["x2"]) / 2, (T["y1"] + T["y2"]) / 2
    field = g.sample(bga, mx, my)
    tol = fabr.get("track_tolerance_pct")
    for li, lname in enumerate(b.layers):
        on = T["layer"] == li
        if not on.any():
            continue
        wmin = float(T["w"][on].min())
        metrics[lname] = {"min_width_mm": round(wmin, 4)}
        if tol:
            metrics[lname]["finished_mm"] = [round(wmin * (1 - tol / 100), 4), round(wmin * (1 + tol / 100), 4)]
        kl = fabr.layer_key("track_width_mm", lname)
        lim_l = fabr.judge(kl, "fab_track")
        for w in np.unique(np.round(T["w"][on], 4)):
            if w >= lim_l - 1e-6 and (local is None or w >= local - 1e-6):
                continue
            sel = on & (np.abs(T["w"] - w) < 5e-5)
            for inb, limit, key in ((False, lim_l, kl), (True, local if local is not None else lim_l,
                                                         "bga_track_mm" if local is not None else kl)):
                s = sel & (field == inb)
                if not s.any() or w >= limit - 1e-6:
                    continue
                k = int(np.flatnonzero(s)[np.argmax(T["len"][s])])
                agg.add((lname, float(w), inb), "error", float(w),
                        f"{int(s.sum())} track segments of {w:g} mm ({float(T['len'][s].sum()):.1f} mm in all)"
                        f"{' in a BGA field' if inb else ''}, below the {limit:g} mm minimum ({fabr.src(key)}: {fabr.limits[key]['note']})",
                        mx[k], my[k], net=b.net_names[T["net"][k]], layer=lname, limit=limit, unit="mm", basis=fabr.src(key),
                        kind=f"{w:g}{'b' if inb else ''}")
    return dict(result("fab_track", "MF3", "Track width against the fabricator's minimum",
                       {"limit_mm": lim, "bga_limit_mm": local, "tolerance_pct": tol, "source": fabr.src("track_width_mm")},
                       metrics, agg.findings()), checked=int(len(T["w"])))


# ------------------------------------------------------------------------------------------------
# MF13

def rule_fab_drill(ctx, fabr):
    from . import thermal
    b = ctx.b
    g = b.grid()
    V = b.vias
    L = len(b.layers)
    lim = lambda k: fabr.judge(k, "fab_drill")
    agg = Agg("fab_drill")
    metrics = {}
    # ---- vias
    vh, vd, vr = lim("via_hole_min_mm"), lim("via_diameter_min_mm"), lim("via_ring_min_mm")
    if len(V["x"]):
        ring = (V["d"] - V["drill"]) / 2
        for name, arr, limit, key, what in (("via_hole", V["drill"], vh, "via_hole_min_mm", "via holes of {v:g} mm, below the {l:g} mm minimum"),
                                            ("via_diameter", V["d"], vd, "via_diameter_min_mm", "via pads of {v:g} mm, below the {l:g} mm minimum"),
                                            ("via_ring", ring, vr, "via_ring_min_mm", "via annular rings of {v:g} mm, below {l:g} mm")):
            if limit is None:
                continue
            for u in np.unique(np.round(arr, 4)):
                if u >= limit - 1e-6:
                    continue
                sel = np.flatnonzero(np.abs(arr - u) < 5e-5)
                agg.add((name, float(u)), "error", float(u), what.format(v=u, l=limit) + f" ({fabr.src(key)})", V["x"][sel[0]], V["y"][sel[0]],
                        n=len(sel), net=b.net_names[V["net"][sel[0]]], limit=limit, unit="mm", basis=fabr.src(key), kind=name)
        metrics["vias"] = {"count": int(len(V["x"])), "min_hole_mm": round(float(V["drill"].min()), 3),
                           "min_ring_mm": round(float(ring.min()), 3)}
    # ---- pad holes
    holes = [p for p in b.pads if p["kind"] in ("pth", "npth") and p.get("drill", 0) > 0]
    rmin, rrec = lim("pth_ring_min_mm"), lim("pth_ring_rec_mm")
    nmin, dmax, slot = lim("npth_min_mm"), lim("drill_max_mm"), lim("plated_slot_min_mm")
    for p in holes:
        dx, dy = p["drill"], p.get("drill_y") or p["drill"]
        small = min(dx, dy)
        if p["kind"] == "npth":
            if nmin is not None and small < nmin - 1e-6:
                agg.add(("npth", round(small, 3)), "error", small, f"non-plated hole of {small:g} mm, below the {nmin:g} mm minimum "
                        f"({fabr.src('npth_min_mm')})", p["x"], p["y"], part=p["ref"], limit=nmin, unit="mm", basis=fabr.src("npth_min_mm"), kind="npth")
        else:
            rg = min(p["sx"] - dx, p["sy"] - dy) / 2
            if rmin is not None and rg < rmin - TOL:
                agg.add(("pth_ring", p["ref"]), "error", rg, f"{p['ref']}.{p['num']}: annular ring of {rg:.3f} mm, below the absolute minimum "
                        f"of {rmin:g} mm ({fabr.src('pth_ring_min_mm')}: {fabr.limits['pth_ring_min_mm']['note']})", p["x"], p["y"],
                        part=p["ref"], net=p["net"], limit=rmin, unit="mm", basis=fabr.src("pth_ring_min_mm"), kind="pth_ring:" + p["ref"])
            elif rrec is not None and rg < rrec - TOL:
                agg.add(("pth_ring_rec", p["ref"]), "warning", rg, f"{p['ref']}.{p['num']}: annular ring of {rg:.3f} mm, below the recommended "
                        f"{rrec:g} mm ({fabr.src('pth_ring_rec_mm')})", p["x"], p["y"], part=p["ref"], net=p["net"], limit=rrec, unit="mm",
                        basis=fabr.src("pth_ring_rec_mm"), kind="pth_ring_rec:" + p["ref"])
            if slot is not None and abs(dx - dy) > 1e-6 and small < slot - 1e-6:
                agg.add(("slot", p["ref"]), "error", small, f"{p['ref']}.{p['num']}: plated slot {small:g} mm wide, below the {slot:g} mm "
                        f"minimum ({fabr.src('plated_slot_min_mm')})", p["x"], p["y"], part=p["ref"], limit=slot, unit="mm",
                        basis=fabr.src("plated_slot_min_mm"), kind="slot:" + p["ref"])
        if dmax is not None and max(dx, dy) > dmax + 1e-6 and abs(dx - dy) < 1e-6:
            agg.add(("drill_max", p["ref"]), "error", -max(dx, dy), f"{p['ref']}.{p['num']}: hole of {max(dx, dy):g} mm, above the largest drill of "
                    f"{dmax:g} mm ({fabr.src('drill_max_mm')})", p["x"], p["y"], part=p["ref"], limit=dmax, unit="mm",
                    basis=fabr.src("drill_max_mm"), kind="drill_max:" + p["ref"])
    # ---- holes as stadiums: a round hole is a point, a slot the segment between its two end centres
    def stadium(x, y, dx, dy, rot):
        if abs(dx - dy) < 1e-6:
            return x, y, x, y, dx / 2
        r = min(dx, dy) / 2
        half = (max(dx, dy) - min(dx, dy)) / 2
        a = math.radians(rot or 0.0) + (0.0 if dx >= dy else math.pi / 2)
        ux, uy = math.cos(a), -math.sin(a)            # KiCad: y grows downwards, angles turn counter-clockwise
        return x - half * ux, y - half * uy, x + half * ux, y + half * uy, r
    H = [stadium(V["x"][k], V["y"][k], V["drill"][k], V["drill"][k], 0) + ("via", int(V["net"][k]), int(V["top"][k]), int(V["bot"][k]), None)
         for k in range(len(V["x"]))]
    H += [stadium(p["x"], p["y"], p["drill"], p.get("drill_y") or p["drill"], p.get("rot", 0)) + (p["kind"], p["net_id"], 0, L - 1, p)
          for p in holes]
    hseg = np.array([h[:4] for h in H]) if H else np.zeros((0, 4))
    hr = np.array([h[4] for h in H]) if H else np.zeros(0)

    def seg_seg(a, b_):
        return geom.seg_seg_dist(tuple(a), tuple(b_))

    # ---- castellated holes: plated holes the board edge cuts through
    outline = b.outline_rings()
    cmin, cspace = lim("castellated_hole_min_mm"), lim("castellated_spacing_mm")
    cast = []
    for k, h in enumerate(H):
        if h[5] != "pth" or not outline:
            continue
        cx, cy = (h[0] + h[2]) / 2, (h[1] + h[3]) / 2
        s_ = geom.ring_segments(outline)
        dmin = float(geom.point_seg_dist(cx, cy, s_[:, 0], s_[:, 1], s_[:, 2], s_[:, 3]).min())
        if dmin < h[4] - 1e-3:
            cast.append(k)
            p = h[9]
            if cmin is not None and 2 * h[4] < cmin - 1e-6:
                agg.add(("castellated", p["ref"]), "error", 2 * h[4], f"{p['ref']}.{p['num']}: castellated hole of {2 * h[4]:g} mm, below the "
                        f"{cmin:g} mm minimum ({fabr.src('castellated_hole_min_mm')})", p["x"], p["y"], part=p["ref"], limit=cmin, unit="mm",
                        basis=fabr.src("castellated_hole_min_mm"), kind="castellated:" + p["ref"])
    if cspace is not None:
        for a_ in range(len(cast)):
            for c_ in range(a_ + 1, len(cast)):
                i_, j_ = cast[a_], cast[c_]
                wall = seg_seg(hseg[i_], hseg[j_]) - hr[i_] - hr[j_]
                if wall < cspace - 1e-3:
                    p = H[i_][9]
                    agg.add(("castellated_spacing", p["ref"]), "error", wall, f"castellated holes {max(wall, 0):.3f} mm apart hole edge to hole "
                            f"edge, below {cspace:g} mm ({fabr.src('castellated_spacing_mm')})", p["x"], p["y"], part=p["ref"], limit=cspace,
                            unit="mm", basis=fabr.src("castellated_spacing_mm"), kind="castellated_spacing")
    metrics["castellated"] = {"holes": len(cast)}

    # ---- hole to hole: holes that cut into each other fail as such (geometry); the fabricator's figure
    # is judged in the measurand the profile records (centre to centre, the owner's choice), and a pair
    # that meets it only in that reading is noted
    hv, hp = lim("via_hole_spacing_mm"), lim("pad_hole_spacing_mm")
    is_via = np.array([h[5] == "via" for h in H], bool)
    if len(H) > 1:
        cxy = np.c_[(hseg[:, 0] + hseg[:, 2]) / 2, (hseg[:, 1] + hseg[:, 3]) / 2]
        span = np.hypot(hseg[:, 2] - hseg[:, 0], hseg[:, 3] - hseg[:, 1]) / 2
        reach = 2 * float((hr + span).max()) + max([v for v in (hv, hp) if v is not None] + [0.0])
        for i_, j_ in cKDTree(cxy).query_pairs(reach):
            cc = seg_seg(hseg[i_], hseg[j_])
            wall = cc - hr[i_] - hr[j_]
            mx, my = (cxy[i_][0] + cxy[j_][0]) / 2, (cxy[i_][1] + cxy[j_][1]) / 2
            if wall < 1e-3:
                agg.add(("hole_overlap",), "error", wall, f"two holes {'touch' if wall > -1e-3 else f'overlap by {-wall:.3f} mm'}: "
                        f"drills that cut into each other make one ragged hole", mx, my, limit=0.0, unit="mm", basis="geometry",
                        kind="hole_overlap")
                continue
            limit, key = (hv, "via_hole_spacing_mm") if (is_via[i_] or is_via[j_]) else (hp, "pad_hole_spacing_mm")
            if limit is None:
                continue
            centre = fabr.limits[key].get("measurand", "hole_wall_to_wall") == "hole_centre_to_centre"
            val = cc if centre else wall
            how = "centre to centre" if centre else "wall to wall"
            if val < limit - 1e-3:
                agg.add(("hole_hole", key), "error", val, f"two holes {val:.3f} mm apart {how} ({wall:.3f} mm wall to wall), below "
                        f"{limit:g} mm ({fabr.src(key)})", mx, my, limit=limit, unit="mm", basis=fabr.src(key), kind="hole_hole:" + key)
            elif centre and wall < limit - 1e-3:
                agg.add(("hole_hole_wall", key), "info", wall, f"two holes {cc:.3f} mm apart centre to centre and {wall:.3f} mm wall to wall: "
                        f"they meet {fabr.name}'s {limit:g} mm only if it is measured centre to centre ({fabr.src(key)})", mx, my,
                        limit=limit, unit="mm", basis=fabr.src(key), kind="hole_hole_wall:" + key)
    # ---- hole to copper of another net, on every layer the hole passes through. Candidates come from
    # the board raster; each is measured again on a fine local raster against the exact hole shape
    # (stadium), so the measuring slack is FINE / 2 plus the 2 um arc approximation.
    lim_c = {"via": (lim("via_hole_to_copper_mm"), lim("inner_via_hole_to_copper_mm"), "via_hole_to_copper_mm", "inner_via_hole_to_copper_mm"),
             "pth": (lim("pth_hole_to_copper_mm"), lim("inner_pth_hole_to_copper_mm"), "pth_hole_to_copper_mm", "inner_pth_hole_to_copper_mm"),
             "npth": (lim("npth_hole_to_copper_mm"), lim("npth_hole_to_copper_mm"), "npth_hole_to_copper_mm", "npth_hole_to_copper_mm")}
    rast = [b.raster(l) for l in range(L)]
    slack = FINE_EDGE / 2 + 0.002
    cand = 0
    for (x1_, y1_, x2_, y2_, r, kind, nid, top, bot, _p) in H:
        outer, inner, ko, ki = lim_c[kind]
        big = max(v for v in (outer, inner, 0.0) if v is not None)
        if big <= 0:
            continue
        R = r + big + 2 * g.res
        i0, j0 = g.ij(min(x1_, x2_) - R, min(y1_, y2_) - R)
        i1, j1 = g.ij(max(x1_, x2_) + R, max(y1_, y2_) + R)
        i0, j0, i1, j1 = max(int(i0), 0), max(int(j0), 0), min(int(i1) + 1, g.ny), min(int(j1) + 1, g.nx)
        if i0 >= i1 or j0 >= j1:
            continue
        xs = g.x0 + (np.arange(j0, j1) + 0.5) * g.res
        ys = g.y0 + (np.arange(i0, i1) + 0.5) * g.res
        X, Y = np.meshgrid(xs, ys)
        dist = geom.point_seg_dist(X, Y, x1_, y1_, x2_, y2_) - r - g.res / 2
        for l in range(top, bot + 1):
            limit, key = (outer, ko) if l in (0, L - 1) else (inner, ki)
            if limit is None:
                continue
            sub = rast[l][i0:i1, j0:j1]
            other = (sub != 0) & (sub != nid) if nid else (sub != 0)
            if not (other & (dist < limit + slack + 1.5 * g.res)).any():
                continue
            cand += 1
            # fine measurement around the hole
            W = r + limit + 0.05
            fg, fr = b.raster_window(l, min(x1_, x2_) - W, min(y1_, y2_) - W, max(x1_, x2_) + W, max(y1_, y2_) + W, FINE_EDGE)
            fx = fg.x0 + (np.arange(fg.nx) + 0.5) * fg.res
            fy = fg.y0 + (np.arange(fg.ny) + 0.5) * fg.res
            FX, FY = np.meshgrid(fx, fy)
            fd = geom.point_seg_dist(FX, FY, x1_, y1_, x2_, y2_) - r - fg.res / 2
            fo = (fr != 0) & (fr != nid) if nid else (fr != 0)
            if not fo.any():
                continue
            d = float(fd[fo].min())
            if d >= limit + slack:
                continue
            k = int(np.argmin(np.where(fo, fd, 1e9)))
            net = b.net_names[int(fr.ravel()[k])]
            what = "via" if kind == "via" else kind.upper() + " hole"
            x, y = (x1_ + x2_) / 2, (y1_ + y2_) / 2
            who = _p["ref"] if _p else None
            if d < limit - slack:
                agg.add(("hole_copper", kind, b.layers[l], who), "error", d, f"{what} wall {max(d, 0):.3f} mm from copper of {net}, below "
                        f"{limit:g} mm ({fabr.src(key)})", x, y, net=net, layer=b.layers[l], limit=limit, unit="mm", basis=fabr.src(key),
                        kind=f"hole_copper:{kind}")
            else:
                agg.add(("hole_copper_marginal", kind, l in (0, L - 1)), "warning", d, f"{what} wall {max(d, 0):.3f} mm from copper of {net}: "
                        f"within the measuring slack ({slack * 1000:.1f} um) of the {limit:g} mm limit ({fabr.src(key)}), so neither a "
                        f"pass nor a fail can be shown here; KiCad's DRC with the hole clearance set to {limit:g} mm settles it exactly",
                        x, y, net=net, layer=b.layers[l], limit=limit, unit="mm", basis=fabr.src(key),
                        kind=f"hole_copper_marginal:{kind}")
    metrics["hole_copper"] = {"candidates_measured_fine": cand, "slack_mm": round(slack, 4)}
    # ---- copper to the board edge: candidates from the raster, then the fine raster against the
    # outline polygon itself
    edge = lim("copper_to_edge_mm")
    if edge is not None and outline:
        inside = thermal._inside(ctx)
        d_out = ndimage.distance_transform_edt(inside) * g.res
        s_ = geom.ring_segments(outline)
        # castellated pads reach the edge by design: their copper is not judged against the edge
        cast_rings = [ring for k in cast for poly in H[k][9]["poly"] for ring in poly]
        cast_mask = ndimage.binary_dilation(geom.raster_rings(g, cast_rings), iterations=2) if cast_rings else np.zeros(g.shape, bool)
        for l in range(L):
            bad = (rast[l] != 0) & ~cast_mask & (d_out < edge + slack + 2 * g.res)
            if not bad.any():
                continue
            lbl, n = ndimage.label(bad)
            # 1 mm tiles of candidate pixels, nearest the edge first; a tile whose coarse distance cannot
            # beat the best fine reading so far is not refined
            ii, jj = np.nonzero(bad)
            tpx = max(int(1.0 / g.res), 1)
            tiles = {}
            for a_, c_ in zip(ii // tpx, jj // tpx):
                tiles.setdefault((int(a_), int(c_)), None)
            coarse = sorted(((float(d_out[a_ * tpx:(a_ + 1) * tpx, c_ * tpx:(c_ + 1) * tpx][bad[a_ * tpx:(a_ + 1) * tpx, c_ * tpx:(c_ + 1) * tpx]].min()),
                              a_, c_) for a_, c_ in tiles))
            worst = None
            for cd, a_, c_ in coarse:
                if worst is not None and cd - 2 * g.res > worst[0]:
                    break
                x0_, y0_ = g.xy(a_ * tpx, c_ * tpx)
                x1w, y1w = g.xy((a_ + 1) * tpx, (c_ + 1) * tpx)
                fg, fr = b.raster_window(l, min(x0_, x1w) - g.res, min(y0_, y1w) - g.res, max(x0_, x1w) + g.res, max(y0_, y1w) + g.res,
                                         FINE_EDGE)
                if not (fr != 0).any():
                    continue
                fx = fg.x0 + (np.arange(fg.nx) + 0.5) * fg.res
                fy = fg.y0 + (np.arange(fg.ny) + 0.5) * fg.res
                FX, FY = np.meshgrid(fx, fy)
                pts = np.c_[FX[fr != 0], FY[fr != 0]]
                ids = fr[fr != 0]
                near = np.full(len(pts), 1e9)
                for a0, b0, a1, b1 in s_:
                    near = np.minimum(near, geom.point_seg_dist(pts[:, 0], pts[:, 1], a0, b0, a1, b1))
                # copper outside the outline (a fill drawn past the edge) is at distance 0
                outside = ~geom.points_in_rings(pts[:, 0], pts[:, 1], outline)
                near[outside] = 0.0
                if cast_rings:
                    near[geom.points_in_rings(pts[:, 0], pts[:, 1], cast_rings)] = 1e9
                k = int(np.argmin(near))
                d = float(near[k]) - fg.res / 2
                if worst is None or d < worst[0]:
                    worst = (d, pts[k][0], pts[k][1], int(ids[k]))
            if worst is None or worst[0] >= edge + slack:
                continue
            d, x, y, nid_ = worst
            sev, how = ("error", "below") if d < edge - slack else ("warning", "within the measuring slack of")
            agg.add(("edge" if sev == "error" else "edge_marginal", b.layers[l]), sev, d, f"copper of {b.net_names[nid_]} {max(d, 0):.3f} mm "
                    f"from the board edge, {how} {edge:g} mm ({fabr.src('copper_to_edge_mm')})", float(x), float(y), n=n,
                    net=b.net_names[nid_], layer=b.layers[l], limit=edge, unit="mm", basis=fabr.src("copper_to_edge_mm"),
                    kind="edge" if sev == "error" else "edge_marginal")
    # ---- board thickness
    t = b.stack.thickness
    tmin, tmax = lim("thickness_min_mm"), lim("thickness_max_mm")
    if tmin is not None and tmax is not None and not (tmin - 1e-3 <= t <= tmax + 1e-3):
        agg.add(("thickness",), "error", t, f"board thickness {t:.2f} mm, outside {tmin:g} to {tmax:g} mm ({fabr.src('thickness_min_mm')})",
                limit=tmin if t < tmin else tmax, unit="mm", basis=fabr.src("thickness_min_mm"), kind="thickness")
    plated = np.r_[V["drill"], [min(p["drill"], p.get("drill_y") or p["drill"]) for p in holes if p["kind"] == "pth"]]
    metrics["board"] = {"thickness_mm": round(t, 3), "holes": int(len(H)),
                        "aspect_ratio": round(t / float(plated.min()), 1) if len(plated) else None}
    return dict(result("fab_drill", "MF13", "Holes, annular rings and the clearances around them", {
        "measuring_slack_mm": round(slack, 4)}, metrics, agg.findings()), checked=int(len(H)))


# ------------------------------------------------------------------------------------------------
# MF4 / MF5 / MF14

def rule_mask_web(ctx, fabr):
    b = ctx.b
    lim = lambda k: fabr.judge(k, "mask_web")
    web, spacing, pp, fexp = lim("mask_web_mm"), lim("mask_pad_spacing_mm"), lim("smd_pad_to_pad_mm"), lim("mask_expansion_mm")
    pmin = lim("smd_pad_min_mm")
    reach = max(v for v in (web, spacing, pp, 0.05) if v is not None) + 0.02
    worst, metrics, npads = {}, {}, 0
    for mname, li in (("F.Mask", 0), ("B.Mask", len(b.layers) - 1)):
        P = sorted(_mask_pads(b, mname), key=lambda q: q["bb"][0])
        npads += len(P)
        pairs, min_web = 0, None
        for a_, A in enumerate(P):
            for B in P[a_ + 1:]:
                if B["bb"][0] > A["bb"][2] + reach:
                    break
                if B["bb"][1] > A["bb"][3] + reach or B["bb"][3] < A["bb"][1] - reach:
                    continue
                pa, pb = A["pad"], B["pad"]
                if pa["ref"] == pb["ref"] and pa["num"] == pb["num"]:
                    continue
                w = geom.rings_distance(A["mask"], B["mask"])
                if w > reach:
                    continue
                gap = geom.rings_distance(A["cu"], B["cu"])
                same = pa["net_id"] == pb["net_id"] and pa["net_id"] != 0
                if same and gap <= TOL:
                    continue                                   # one piece of copper
                pairs += 1
                min_web = w if min_web is None else min(min_web, w)
                x, y = (pa["x"] + pb["x"]) / 2, (pa["y"] + pb["y"]) / 2
                name = f"{pa['ref']}.{pa['num']} - {pb['ref']}.{pb['num']}"
                at_fab = "" if fexp is None else (f"; at the {fexp * 1000:.0f} um expansion {fabr.name} says it plots, "
                                                  f"{max(gap - 2 * fexp, 0):.3f} mm would remain")
                found = []
                if same:
                    if web is not None and w < web - TOL:
                        sev = "warning" if pa["ref"] != pb["ref"] else "info"
                        found.append(("same_net", sev, w, web, "mask_web_mm", f"pads of one net ({pa['net']}) {gap:.3f} mm apart with "
                                      f"{'no' if w <= TOL else f'a {w:.3f} mm'} solder-mask web between them: the joints share solder, "
                                      f"which can starve one or pull a part{at_fab}"))
                else:
                    if pp is not None and pa["kind"] == pb["kind"] == "smd" and gap < pp - TOL:
                        found.append(("pad_gap", "error", gap, pp, "smd_pad_to_pad_mm", f"pads {gap:.3f} mm apart, below the {pp:g} mm "
                                      f"minimum between SMD pads of different nets ({fabr.src('smd_pad_to_pad_mm')})"))
                    if w <= TOL:
                        if spacing is not None and gap < spacing - TOL:
                            found.append(("no_dam", "warning", gap, spacing, "mask_pad_spacing_mm", f"no solder-mask web: pads {gap:.3f} mm "
                                          f"apart, less than the {spacing:g} mm that keeps a web ({fabr.src('mask_pad_spacing_mm')}); "
                                          f"solder can bridge"))
                        else:
                            found.append(("no_dam", "warning", w, web, "mask_web_mm", f"the mask openings merge as drawn although the pads "
                                          f"are {gap:.3f} mm apart (expansion {max(A['exp'], B['exp']) * 1000:.0f} um){at_fab}"))
                    elif web is not None and w < web - TOL:
                        found.append(("web", "error", w, web, "mask_web_mm", f"solder-mask web of {w:.3f} mm as drawn, below the {web:g} mm "
                                      f"minimum ({fabr.src('mask_web_mm')}): it will be removed or lift{at_fab}"))
                for kind, sev, val, limit, key, msg in found:
                    k = (mname, tuple(sorted((pa["ref"], pb["ref"]))), kind)
                    e = worst.setdefault(k, {"n": 0, "val": 1e9})
                    e["n"] += 1
                    if val < e["val"]:
                        e.update(val=val, sev=sev, limit=limit, msg=msg, name=name, x=x, y=y, layer=b.layers[li],
                                 nets=(pa["net"], pb["net"]), part=pa["ref"], key=key)
        metrics[mname] = {"openings": len(P), "close_pairs": pairs, "min_web_mm": None if min_web is None else round(max(min_web, 0.0), 4)}
        if pmin is not None:
            for q in P:
                p = q["pad"]
                # a copper pad smaller than the smallest pad both ways; its size is measured from its
                # copper (a custom pad's anchor is not its size; mask- or paste-only openings are not pads)
                if p["kind"] == "smd" and p["layers"] and max(p["sx"], p["sy"]) < pmin - 1e-6:
                    e = worst.setdefault((mname, (p["ref"], p["ref"]), "pad_min"), {"n": 0, "val": 1e9})
                    e["n"] += 1
                    if max(p["sx"], p["sy"]) < e["val"]:
                        e.update(val=max(p["sx"], p["sy"]), sev="error", limit=pmin, name=f"{p['ref']}.{p['num']}", x=p["x"], y=p["y"],
                                 layer=b.layers[li], nets=(p["net"], p["net"]), part=p["ref"], key="smd_pad_min_mm",
                                 msg=f"SMD pad {p['sx']:g} x {p['sy']:g} mm, below the smallest pad of {pmin:g} mm ({fabr.src('smd_pad_min_mm')})")
    fs = []
    for (mname, refs, kind), e in sorted(worst.items(), key=lambda kv: (kv[1]["sev"] != "error", kv[0][2], kv[1]["val"])):
        more = (f" ({e['n']} {'pads' if kind == 'pad_min' else 'pad pairs'} of {' / '.join(dict.fromkeys(refs))}; worst shown)"
                if e["n"] > 1 else "")
        fs.append(finding("mask_web", e["sev"], f"{e['name']}: {e['msg']}{more}", e["nets"][0], e["layer"], e["x"], e["y"],
                          value=round(e["val"], 4), limit=e["limit"], unit="mm", aggressor=e["nets"][1] if e["nets"][1] != e["nets"][0] else "",
                          part=e["part"], kind=f"{kind}:{'/'.join(refs)}", pairs=e["n"], basis=fabr.src(e["key"])))
    return dict(result("mask_web", "MF4", "Solder-mask web and pad spacing", {"web_mm": web, "pad_spacing_mm": spacing,
                       "smd_pad_to_pad_mm": pp, "fab_expansion_mm": fexp}, metrics, fs), checked=npads)


def rule_mask_exposed(ctx, fabr, res):
    """Per pad with a mask opening: copper of another net inside the opening (fail) or closer to it than
    the fabricator's opening-to-track clearance (warn). Also the copper of the pad's own net that the
    opening uncovers beyond the pad (`exposure`, listed with the chips)."""
    b = ctx.b
    g = b.grid()
    lim = fabr.judge("mask_to_track_mm", "mask_exposed")
    reach = (lim or 0.1) + 0.1
    worst, exposure, checked, total = {}, {}, 0, 0
    for mname, li in (("F.Mask", 0), ("B.Mask", len(b.layers) - 1)):
        r = b.raster(li)
        P = _mask_pads(b, mname)
        total += len(P)
        bbs = np.array([q["bb"] for q in P]) if P else np.zeros((0, 4))
        for kq, q in enumerate(P):
            p = q["pad"]
            nid = p["net_id"]
            x0, y0, x1, y1 = q["bb"]
            i0, j0 = g.ij(x0 - reach, y0 - reach)
            i1, j1 = g.ij(x1 + reach, y1 + reach)
            sub = r[max(int(i0), 0):int(i1) + 1, max(int(j0), 0):int(j1) + 1]
            if not ((sub != 0) & (sub != nid)).any() and q["exp"] <= 0:
                continue
            checked += 1
            gw, lab = b.raster_window(li, x0 - reach, y0 - reach, x1 + reach, y1 + reach, res)
            O = geom.raster_rings(gw, q["mask"])
            Pm = geom.raster_rings(gw, q["cu"])
            own = (lab == nid) if nid else Pm
            if Pm.any():
                exposure[(p["ref"], p["num"])] = 1.0 + int((O & own & ~Pm).sum()) / float(Pm.sum())
            foreign = (lab != 0) & ~own
            # copper under another pad's own opening is that pad (MF4 judges the web between the two)
            for ko in np.flatnonzero((bbs[:, 0] <= x1 + reach) & (bbs[:, 2] >= x0 - reach) & (bbs[:, 1] <= y1 + reach) & (bbs[:, 3] >= y0 - reach)):
                if ko != kq:
                    foreign &= ~geom.raster_rings(gw, P[ko]["mask"])
            if not foreign.any():
                continue
            inside = O & foreign
            found = None
            if inside.sum() >= 6:
                net = b.net_names[int(np.bincount(lab[inside]).argmax())]
                area = float(inside.sum()) * res * res
                found = ("exposed", "error", -area, 0.0, net, "as drawn",
                         f"{area * 1e6:.0f} um2 of {net} lies inside the mask opening of the pad as drawn: solder can bridge to it")
            elif lim is not None:
                d = ndimage.distance_transform_edt(~O)
                dmin = (float(d[foreign].min()) - 1.0) * res
                if dmin < lim - 1.5 * res:
                    k = np.flatnonzero(foreign.ravel())[int(np.argmin(d[foreign]))]
                    net = b.net_names[int(lab.ravel()[k])]
                    found = ("near", "warning", dmin, lim, net, fabr.src("mask_to_track_mm"),
                             f"{net} is {max(dmin, 0) * 1000:.0f} um from the mask opening, closer than {lim * 1000:.0f} um "
                             f"({fabr.src('mask_to_track_mm')}): mask misregistration can uncover it")
            if found:
                kind, sev, val, limit, net, basis, msg = found
                e = worst.setdefault((p["ref"], kind), {"n": 0, "val": 1e9})
                e["n"] += 1
                if val < e["val"]:
                    e.update(val=val, sev=sev, limit=limit, msg=msg, pad=p, layer=b.layers[li], other=net, basis=basis)
    fs = []
    for (ref, kind), e in sorted(worst.items(), key=lambda kv: (kv[0][1], kv[1]["val"])):
        p = e["pad"]
        more = f" ({e['n']} pads of {ref}; worst shown)" if e["n"] > 1 else ""
        val = -e["val"] if kind == "exposed" else e["val"]
        fs.append(finding("mask_exposed", e["sev"], f"{ref}.{p['num']}: {e['msg']}{more}", p["net"], e["layer"], p["x"], p["y"],
                          value=round(val, 5), limit=e["limit"], unit="mm2" if kind == "exposed" else "mm", aggressor=e["other"],
                          part=ref, kind=kind, pads=e["n"], basis=e["basis"]))
    metrics = {"pads": {"openings": total, "measured": checked, "with_exposure": sum(1 for v in exposure.values() if v > 1.02)}}
    return dict(result("mask_exposed", "MF5", "Copper of another net inside or beside a mask opening",
                       {"opening_to_track_mm": lim, "res_mm": res}, metrics, fs), checked=total), exposure


def rule_via_pad(ctx, fabr):
    b, cfg = ctx.b, ctx.cfg
    V = b.vias
    L = len(b.layers)
    cover = cfg["dfm"]["via_cover"] or fabr.get("via_cover") or "tented"
    csrc = "board config" if cfg["dfm"]["via_cover"] else (fabr.src("via_cover") + ": " + fabr.limits["via_cover"]["note"]
                                                           if "via_cover" in fabr.limits else "assumed")
    web = fabr.judge("mask_web_mm", "via_pad")
    fmin, fmax = fabr.judge("filled_via_min_mm", "via_pad"), fabr.judge("filled_via_max_mm", "via_pad")
    plug = fabr.judge("plugged_to_opening_mm", "via_pad")
    pmax = fabr.judge("plugged_via_max_mm", "via_pad")
    agg = Agg("via_pad")
    n_in = n_pads = 0
    for mname, li, side in (("F.Mask", 0, 0), ("B.Mask", L - 1, 1)):
        sel = np.flatnonzero((V["top"] <= li) & (V["bot"] >= li)) if len(V["x"]) else np.zeros(0, int)
        P = [q for q in _mask_pads(b, mname) if q["pad"]["kind"] == "smd"]
        n_pads += len(P)
        if not len(sel) or not P:
            continue
        tree = cKDTree(np.c_[V["x"][sel], V["y"][sel]])
        vmax = float(V["d"][sel].max()) / 2
        for q in P:
            p = q["pad"]
            x0, y0, x1, y1 = q["bb"]
            reach = math.hypot(x1 - x0, y1 - y0) / 2 + vmax + max(web or 0.1, plug or 0.0) + 0.1
            for kk in tree.query_ball_point(((x0 + x1) / 2, (y0 + y1) / 2), reach):
                k = int(sel[kk])
                vx, vy, vdia, vdr = V["x"][k], V["y"][k], V["d"][k], V["drill"][k]
                same = int(V["net"][k]) == p["net_id"] and p["net_id"] != 0
                dcu = _pt_rings(vx, vy, q["cu"])
                if same and dcu < vdr / 2:
                    n_in += 1
                    if cover == "filled":
                        # the range is of drill sizes (the profile's measurand; the via-covering article
                        # gives its upper bound as a hole size)
                        if fmin is not None and not (fmin - 1e-6 <= vdr <= fmax + 1e-6):
                            agg.add(("fill_range", p["ref"]), "error", vdr, f"{p['ref']}.{p['num']}: via hole of {vdr:g} mm in the pad; the "
                                    f"filled and capped process takes holes of {fmin:g} to {fmax:g} mm ({fabr.src('filled_via_max_mm')})", vx, vy,
                                    part=p["ref"], net=p["net"], layer=b.layers[li], limit=fmax if vdr > fmax else fmin, unit="mm",
                                    basis=fabr.src("filled_via_max_mm"), kind="fill_range:" + p["ref"])
                    else:
                        agg.add(("in_pad", p["ref"]), "error", vdr, f"{p['ref']}.{p['num']}: via in an SMD pad, and the vias are {cover} "
                                f"({csrc}): solder drains into the hole. Order filled and capped vias or move the via off the pad",
                                vx, vy, part=p["ref"], net=p["net"], layer=b.layers[li], unit="mm", basis=csrc, kind="in_pad:" + p["ref"])
                    continue
                if cover == "filled" or dcu < vdr / 2:
                    continue
                if cover == "plugged" and pmax is not None and vdia > pmax + 1e-6:
                    agg.add(("plug_max", round(float(vdia), 3)), "error", -vdia, f"plugged via of {vdia:g} mm; a via filled with solder mask "
                            f"may be at most {pmax:g} mm ({fabr.src('plugged_via_max_mm')})", vx, vy, net=b.net_names[int(V['net'][k])],
                            layer=b.layers[li], limit=pmax, unit="mm", basis=fabr.src("plugged_via_max_mm"), kind="plug_max")
                tented = (b.snap["vias"][k].get("tented") or [None, None])[side]
                dm = _pt_rings(vx, vy, q["mask"])
                if cover == "plugged" and plug is not None and dm - vdia / 2 < plug - TOL:
                    agg.add(("plug", p["ref"]), "warning", dm - vdia / 2, f"{p['ref']}.{p['num']}: plugged via {max(dm - vdia / 2, 0):.2f} mm "
                            f"from the pad's mask opening, below {plug:g} mm ({fabr.src('plugged_to_opening_mm')})", vx, vy, part=p["ref"],
                            net=p["net"], layer=b.layers[li], limit=plug, unit="mm", basis=fabr.src("plugged_to_opening_mm"),
                            kind="plug:" + p["ref"])
                elif tented is False and web is not None:
                    exp = b.snap["vias"][k].get("mask_exp") or 0.0
                    w = dm - vdia / 2 - exp
                    if w < web - TOL:
                        agg.add(("open_via", p["ref"]), "warning", w, f"{p['ref']}.{p['num']}: untented via with "
                                f"{'no' if w <= TOL else f'a {w:.3f} mm'} solder-mask web to the pad opening (minimum {web:g} mm, "
                                f"{fabr.src('mask_web_mm')}): solder can run from the pad into the via", vx, vy, part=p["ref"], net=p["net"],
                                layer=b.layers[li], limit=web, unit="mm", basis=fabr.src("mask_web_mm"), kind="open_via:" + p["ref"])
    assumption = None
    if not cfg["dfm"]["via_cover"] and len(V["x"]):
        assumption = (f"vias {cover} ({csrc}); {n_in} vias are in SMD pads" if n_in else f"vias {cover} ({csrc})") + \
            "; state the via covering in the order or set [dfm] via_cover"
    return dict(result("via_pad", "MF14", "Vias in and beside SMD pads", {"via_cover": cover, "via_cover_source": csrc},
                       {"vias": {"in_pads": n_in, "cover": cover}}, agg.findings("vias")), checked=n_pads if len(V["x"]) else 0,
                assumption=assumption)


# ------------------------------------------------------------------------------------------------
# MF15 / MF16 / MF17 / MF18

def _area_perimeter(rings):
    a = max(abs(geom.ring_area(r)) for r in rings)
    r = max(rings, key=lambda r: abs(geom.ring_area(r)))
    q = np.asarray(r + r[:1], dtype=float)
    return a, float(np.hypot(np.diff(q[:, 0]), np.diff(q[:, 1])).sum())


def rule_paste(ctx, fabr):
    b, cfg = ctx.b, ctx.cfg
    ar_min, as_min = fabr.judge("stencil_area_ratio_min", "paste"), fabr.judge("stencil_aspect_ratio_min", "paste")
    ap_min = fabr.judge("stencil_aperture_min_mm", "paste")
    ts = sorted(fabr.get("stencil_thickness_mm") or [])
    agg = Agg("paste")
    assumption, pred = None, None
    if cfg["dfm"]["stencil_mm"]:
        t, tsrc = cfg["dfm"]["stencil_mm"], "board config"
    else:
        guide = fabr.get("stencil_guide")
        pred = packages.stencil_prediction(ctx.b, guide) if guide else {"thickness_mm": None, "uncovered": [], "part": None}
        has_paste = any(p["kind"] in ("smd", "conn", "other") and p.get("paste") for p in ctx.b.pads)
        if pred["thickness_mm"] is None and not has_paste:
            return dict(result("paste", "MF15", "Stencil apertures", {"stencil_mm": None}, {}, []), checked=0)
        if pred["thickness_mm"] is None:
            why = (f"{fabr.name}'s stencil guide does not cover {', '.join(pred['uncovered'][:6])}" if pred["uncovered"] else
                   f"the profile has no stencil guidance for {fabr.name}")
            agg.add(("stencil_unknown",), "warning", 0.0, f"the stencil thickness is not known: {why}. Set [dfm] stencil_mm in the board "
                    f"config to the thickness that will be ordered", basis="board config needed", kind="stencil_unknown")
            return dict(result("paste", "MF15", "Stencil apertures", {"stencil_mm": None, "stencil_source": why}, {"prediction": pred},
                               agg.findings("pads")), checked=0, unchecked=True)
        t = pred["thickness_mm"]
        tsrc = (f"predicted: {fabr.name}'s guide gives {t:g} mm for {pred['row']}, the board's finest ({pred['part']}); "
                f"the fabricator chooses unless the order states it")
        assumption = f"stencil {t:g} mm, as {fabr.name}'s guide gives for {pred['part']} ({pred['row']}); state it in the order or set [dfm] stencil_mm"
    n, worst_ar, need = 0, None, None
    with_paste = set()
    for p in b.pads:
        if p["kind"] not in ("smd", "conn", "other") or not p.get("paste"):
            continue
        fp = b.fps.get(p["ref"], {})
        if fp.get("dnp"):
            continue
        # the aperture is the paste polygon KiCad plots (margins applied, custom shapes as drawn);
        # its width for the aspect ratio is the largest disc that fits in it
        pl = next(iter(p["paste"]))
        ap = (p.get("paste_poly") or {}).get(pl)
        dm = padgeom.measure(ap, p["x"], p["y"], p.get("rot", 0.0)) if ap else None
        if dm is None or dm["area"] <= 0:
            continue
        w, h = dm["w"], dm["h"]
        with_paste.add(p["ref"])
        n += 1
        A, P = dm["area"], dm["perim"]
        small = dm["inscribed_d"]
        lay = b.layers[0] if pl == "F.Paste" else b.layers[-1]
        pname = f"{p['ref']}.{p['num']}" if str(p["num"]).strip() else f"{p['ref']} (aperture without a pin number)"
        ar, asp = A / (P * t), small / t
        worst_ar = ar if worst_ar is None else min(worst_ar, ar)
        if ar_min and as_min:                                  # thickest stencil at which this aperture still prints
            tmax = min(A / (P * ar_min), small / as_min)
            need = tmax if need is None else min(need, tmax)
        if ap_min is not None and small <= ap_min + 1e-6:
            agg.add(("aperture", p["ref"]), "error", small, f"{pname}: paste aperture {small:.3f} mm wide; the stencil's "
                    f"minimum aperture is above {ap_min:g} mm ({fabr.src('stencil_aperture_min_mm')})", p["x"], p["y"], part=p["ref"],
                    net=p["net"], layer=lay, limit=ap_min, unit="mm", basis=fabr.src("stencil_aperture_min_mm"), kind="aperture:" + p["ref"])
        elif (ar_min and ar < ar_min - 1e-3) or (as_min and asp < as_min - 1e-3):
            agg.add(("ratio", p["ref"]), "warning", ar, f"{pname}: paste aperture {w:.2f} x {h:.2f} mm has an area ratio of "
                    f"{ar:.2f} and an aspect ratio of {asp:.1f} with a {t:g} mm stencil ({tsrc}); IPC-7525 asks for at least {ar_min:g} and "
                    f"{as_min:g}: the paste may not release", p["x"], p["y"], part=p["ref"], net=p["net"], layer=lay,
                    limit=ar_min, unit="", basis="industry standard (IPC-7525)", kind="ratio:" + p["ref"])
    for ref, fp in b.fps.items():                              # a placed SMD part with pads but no paste at all
        pads = [p for p in b.pads_by_ref.get(ref, []) if p["kind"] == "smd"]
        if len(pads) >= 2 and ref not in with_paste and not fp.get("dnp") and fp.get("placed", True):
            agg.add(("no_paste", ref), "warning", 0.0, f"{ref}: {len(pads)} SMD pads and no paste aperture: the part cannot be reflowed",
                    fp["x"], fp["y"], part=ref, layer=pads[0]["layers"][0] if pads[0]["layers"] else "", basis="geometry", kind="no_paste:" + ref)
    fit = [x for x in ts if need is None or x <= need + 1e-9]
    metrics = {"stencil": {"apertures": n, "judged_at_mm": t, "min_area_ratio": None if worst_ar is None else round(worst_ar, 2),
                           "thickest_stencil_that_prints_all_mm": None if need is None else round(need, 3),
                           "standard_thicknesses_that_fit_mm": fit}}
    if pred:
        metrics["prediction"] = pred
    return dict(result("paste", "MF15", "Stencil apertures", {"stencil_mm": t, "stencil_source": tsrc, "area_ratio_min": ar_min,
                       "aspect_ratio_min": as_min, "aperture_min_mm": ap_min}, metrics, agg.findings("pads")), checked=n,
                assumption=assumption)


def rule_footprint(ctx):
    """A two-pad chip land pattern whose pads are further apart than the chip's body: the part cannot
    reach both pads. The body length is the EIA nominal (the chip table holds no sourced tolerance, so
    the nominal length is the limit and the finding says so). Nothing else about the land pattern is
    judged; the toe and the ratio of the pad areas are listed in the chip ranking."""
    b = ctx.b
    agg = Agg("footprint")
    n, metrics = 0, {}
    for ref, fp in b.fps.items():
        pads = b.pads_by_ref.get(ref, [])
        chip = packages.chip_of(fp)
        if chip is None or len(pads) != 2 or fp.get("dnp") or any(p["kind"] != "smd" or not p["poly"] for p in pads):
            continue
        n += 1
        lp = packages.land_pattern(chip, pads)
        metrics[ref] = lp
        layer = pads[0]["layers"][0] if pads[0]["layers"] else ""
        lmax = chip.get("L_max") or chip["L"]
        if lp["gap_mm"] >= lmax - 1e-6:
            agg.add(("gap", ref), "error", chip["L"] - lp["gap_mm"], f"{ref}: the pads are {lp['gap_mm']:.2f} mm apart and the footprint is "
                    f"named for a {chip['code']} body {lmax:g} mm long ({lp['body']}): its terminations cannot reach both pads", fp["x"], fp["y"],
                    part=ref, layer=layer, limit=round(lmax, 3), unit="mm",
                    basis="EIA body length" + (" with tolerance" if chip.get("L_max") else " (nominal; no sourced tolerance)"), kind="gap:" + ref)
    return dict(result("footprint", "MF16", "Two-pad land patterns against the chip they are named for", {}, metrics, agg.findings()), checked=n)


def rule_part_spacing(ctx, fabr):
    b = ctx.b
    tab = fabr.judge("part_spacing", "part_spacing")
    if not tab:
        return dict(result("part_spacing", "MF17", "Spacing between parts", {}, {}, []), checked=0)
    cls_i = {c: k for k, c in enumerate(tab["classes"])}
    reach = max(max(r) for r in tab["matrix"]) + 0.05
    agg = Agg("part_spacing")
    parts, unknown = [], []
    no_body = 0
    for ref, fp in b.fps.items():
        pads = [p for p in b.pads_by_ref.get(ref, []) if p["poly"] and p["kind"] == "smd"]
        if not pads or fp.get("dnp") or not fp.get("placed", True):
            continue
        c, how = packages.part_class_how(fp, b.pads_by_ref.get(ref, []))
        if c not in cls_i:
            unknown.append({"ref": ref, "footprint": fp.get("fpid"), "why": how if c is None else f"class {c} is not in the table"})
            continue
        # the part's outline as the assembler measures it: its pads and its body (fabrication-layer outline)
        rings = [_rings(p["poly"]) for p in pads]
        if fp.get("body"):
            rings.append([fp["body"]])
        else:
            no_body += 1
        bbs = np.array([_bbox(r) for r in rings])
        parts.append({"ref": ref, "cls": c, "side": _side(b, ref), "rings": rings, "bbs": bbs,
                      "bb": (bbs[:, 0].min(), bbs[:, 1].min(), bbs[:, 2].max(), bbs[:, 3].max()), "x": fp["x"], "y": fp["y"]})
    parts.sort(key=lambda q: q["bb"][0])
    judged = 0
    for a_, A in enumerate(parts):
        for B in parts[a_ + 1:]:
            if B["bb"][0] > A["bb"][2] + reach:
                break
            if A["side"] != B["side"] or B["bb"][1] > A["bb"][3] + reach or B["bb"][3] < A["bb"][1] - reach:
                continue
            limit = tab["matrix"][cls_i[A["cls"]]][cls_i[B["cls"]]]
            d = 1e9
            for i, ra in enumerate(A["rings"]):
                near = np.flatnonzero((B["bbs"][:, 0] <= A["bbs"][i, 2] + limit) & (B["bbs"][:, 2] >= A["bbs"][i, 0] - limit) &
                                      (B["bbs"][:, 1] <= A["bbs"][i, 3] + limit) & (B["bbs"][:, 3] >= A["bbs"][i, 1] - limit))
                for j in near:
                    d = min(d, geom.rings_distance(ra, B["rings"][j]))
            judged += 1
            if d < limit - TOL:
                big, small = (A, B) if cls_i[A["cls"]] >= cls_i[B["cls"]] else (B, A)
                chips = cls_i[big["cls"]] < 5                      # two chips: one finding per pair; else one per larger part
                agg.add((big["ref"], small["ref"]) if chips else (big["ref"],), "warning", d - limit,
                        f"{A['ref']} ({A['cls']}) and {B['ref']} ({B['cls']}): {d:.2f} mm apart outline to outline (pads and body), "
                        f"below the {limit:g} mm recommended between these packages ({fabr.src('part_spacing')})", (A["x"] + B["x"]) / 2,
                        (A["y"] + B["y"]) / 2, part=A["ref"], layer=b.layers[0 if A["side"] == "F" else -1], limit=limit, unit="mm",
                        basis=fabr.src("part_spacing") + ", recommended", kind=f"{big['ref']}/{small['ref']}" if chips else big["ref"])
    fs = agg.findings("parts closer than recommended")
    if unknown:
        fs.append(finding("part_spacing", "info", f"{len(unknown)} parts have no package class of the table and were not judged: "
                          + "; ".join(f"{u['ref']} ({u['why']})" for u in unknown[:12]) + (" ..." if len(unknown) > 12 else "")
                          + ". Name the class in the footprint or its library to have them judged", kind="not_classified"))
    return dict(result("part_spacing", "MF17", "Spacing between parts against the assembler's table",
                       {"measured": "outline to outline (pads and fabrication-layer body), as the assembler's diagram shows"},
                       {"parts": {"classified": len(parts), "not_classified": len(unknown), "without_body_outline": no_body,
                                  "pairs_judged": judged}, "not_classified": unknown}, fs),
                checked=judged, note=f"{len(unknown)} parts not classified (listed)" if unknown else "")


def rule_assembly(ctx, fabr):
    b = ctx.b
    agg = Agg("assembly")
    placed = {ref: fp for ref, fp in b.fps.items() if not fp.get("dnp") and fp.get("placed", True) and
              any(p["kind"] == "smd" for p in b.pads_by_ref.get(ref, []))}
    sides = {s: [r for r in placed if _side(b, r) == s] for s in ("F", "B")}
    svc = fabr.limits.get("min_package", {}).get("note", "")
    if fabr.judge("double_sided", "assembly") is False and sides["F"] and sides["B"]:
        agg.add(("sides",), "error", 0.0, f"SMD parts on both sides ({len(sides['F'])} front, {len(sides['B'])} back); the {svc} places "
                f"one side only ({fabr.src('double_sided')})", basis=fabr.src("double_sided"), kind="sides")
    minpkg = fabr.judge("min_package", "assembly")
    pitch, bpitch = (fabr.judge(k, "assembly") for k in ("pin_pitch_mm", "bga_pitch_mm"))
    tight = None
    for ref, fp in placed.items():
        chip = packages.chip_of(fp)
        if chip and minpkg and packages.size_rank(chip["code"]) < packages.size_rank(minpkg):
            agg.add(("package", chip["code"]), "error", packages.size_rank(chip["code"]), f"{ref}: {chip['code']} is smaller than the smallest "
                    f"package of the {svc}, {minpkg} ({fabr.src('min_package')})", fp["x"], fp["y"], part=ref, basis=fabr.src("min_package"),
                    kind="package:" + chip["code"])
        pins = {}
        for p in b.pads_by_ref.get(ref, []):                   # one position per pin: a pin drawn as several pads is one pin
            if p["kind"] == "smd" and p["poly"] and str(p["num"]).strip():
                pins.setdefault(str(p["num"]), (p["x"], p["y"]))
        if len(pins) < 6:
            continue
        xy = np.array(list(pins.values()))
        d = cKDTree(xy).query(xy, k=2)[0][:, 1]
        d = float(d[d > 1e-6].min()) if (d > 1e-6).any() else None
        if d is None:
            continue
        tight = d if tight is None else min(tight, d)
        bga = packages.part_class(fp) == "BGA"
        limit, key = (bpitch, "bga_pitch_mm") if bga else (pitch, "pin_pitch_mm")
        if limit is not None and d < limit - 1e-3:
            agg.add(("pitch", ref), "error", d, f"{ref}: {'ball' if bga else 'pin'} pitch {d:.2f} mm, below the {limit:g} mm of the {svc} "
                    f"({fabr.src(key)})", fp["x"], fp["y"], part=ref, limit=limit, unit="mm", basis=fabr.src(key), kind="pitch:" + ref)
    fid = {s: [r for r, fp in b.fps.items() if _side(b, r) == s and ("FIDUCIAL" in (fp.get("fpid") or "").upper()
                                                                   or r.upper().startswith(("FID", "FD")))] for s in ("F", "B")}
    if fabr.judge("fiducials_required", "assembly"):
        for s in ("F", "B"):
            if sides[s] and not fid[s]:
                agg.add(("fiducials", s), "warning", 0.0, f"{len(sides[s])} SMD parts on the {'front' if s == 'F' else 'back'} and no fiducial "
                        f"there; the {svc} lists fiducials as necessary ({fabr.src('fiducials_required')}). They may instead be on the panel "
                        f"rails, which this review does not see", layer=b.layers[0 if s == "F" else -1], basis=fabr.src("fiducials_required"),
                        kind="fiducials:" + s)
    metrics = {"board": {"smd_parts_front": len(sides["F"]), "smd_parts_back": len(sides["B"]), "fiducials_front": len(fid["F"]),
                         "fiducials_back": len(fid["B"]), "tightest_pitch_mm": None if tight is None else round(tight, 3),
                         "rails_required": fabr.get("rails_required")}}
    return dict(result("assembly", "MF18", "What the assembly service accepts", {"service": svc}, metrics, agg.findings("parts")),
                checked=len(placed))


# ------------------------------------------------------------------------------------------------
# MF20: the files that are sent

def _plot(board_path, layers, od, log):
    os.makedirs(od, exist_ok=True)
    for f in os.listdir(od):
        os.remove(os.path.join(od, f))
    cmds = [["kicad-cli", "pcb", "export", "gerbers", "-o", od + os.sep, "-l", ",".join(layers), "--no-protel-ext", "--check-zones", board_path],
            ["kicad-cli", "pcb", "export", "drill", "-o", od + os.sep, "--format", "excellon", "--excellon-units", "mm",
             "--excellon-zeros-format", "decimal", "--drill-origin", "absolute", "--excellon-separate-th", board_path],
            ["kicad-cli", "pcb", "export", "pos", "--format", "csv", "--units", "mm", "--side", "both", "-o",
             os.path.join(od, os.path.splitext(os.path.basename(board_path))[0] + "-pos.csv"), board_path]]
    for c in cmds:
        p = subprocess.run(c, capture_output=True, text=True, timeout=900)
        if p.returncode != 0:
            raise RuntimeError(f"{' '.join(c[:4])} failed: {(p.stderr or p.stdout).strip()[-300:]}")
    return " ; ".join(" ".join(c[:-1]) for c in cmds)


def _blobs(diff, g, cap=4):
    """Places where two images differ by more than a pixel-wide rim: [(area mm2, x, y)], largest first."""
    thick = ndimage.binary_opening(diff, structure=np.ones((2, 2), bool))
    if not thick.any():
        return [], 0
    lbl, n = ndimage.label(thick)
    sizes = ndimage.sum(thick, lbl, np.arange(1, n + 1))
    out = []
    for k in np.argsort(-sizes)[:cap]:
        sl = ndimage.find_objects(lbl == k + 1)[0]
        ii, jj = np.nonzero(lbl[sl] == k + 1)
        x, y = g.xy(sl[0].start + int(ii[0]), sl[1].start + int(jj[0]))
        out.append((float(sizes[k]) * g.res * g.res, float(x), float(y)))
    return out, n


def rule_output(ctx, fabr, od, gerber_dir, log):
    """The Gerber and drill files against the geometry that was checked. With `gerber_dir` the files
    there are read (the ones that will be sent); without it they are plotted now with kicad-cli and
    kept in <run>/dfm/fab_files, and the check then only shows that the reader and the plot agree."""
    from . import gerber
    b = ctx.b
    g = b.grid()
    L = len(b.layers)
    fs = []
    supplied = bool(gerber_dir)
    src = os.path.abspath(gerber_dir) if supplied else os.path.join(od, "fab_files")
    cmd = None
    if not supplied:
        try:
            cmd = _plot(b.path, list(b.layers) + ["F.Mask", "B.Mask", "F.Paste", "B.Paste", "Edge.Cuts"], src, log)
        except (OSError, RuntimeError, subprocess.TimeoutExpired) as e:
            fs.append(finding("output_match", "warning", f"the fabrication files could not be plotted ({e}); the checks ran on the board "
                              f"file only", kind="plot"))
            return dict(result("output_match", "MF20", "The files that are sent against the geometry that was checked", {}, {}, fs),
                        checked=0, unchecked="the files could not be plotted")
    if not os.path.isdir(src):
        raise FileNotFoundError(f"--gerbers {gerber_dir}: no such folder")
    files, drills, unread = {}, [], []
    for name in sorted(os.listdir(src)):
        path = os.path.join(src, name)
        ext = os.path.splitext(name)[1].lower()
        if ext in (".drl", ".xln", ".exc", ".nc") or (ext == ".txt" and "drl" in name.lower()):
            drills.append(gerber.read_drill(path))
            continue
        if not (ext == ".gbr" or ext in gerber.PROTEL or (len(ext) >= 3 and ext[1] == "g" and ext[2:].isdigit())):
            continue
        try:
            G = gerber.read(path)
        except Exception as e:
            unread.append(f"{name}: {e}")
            continue
        layer = gerber.layer_of(path, G, b.layers)
        if layer:
            files[layer] = G
    for msg in unread:
        fs.append(finding("output_match", "warning", f"could not be read: {msg}", kind="unread"))
    # position of the files against the board (an auxiliary origin moves them)
    dx = dy = 0.0
    ref_layer = next((l for l in b.layers if l in files and gerber.bbox(files[l])), None)
    if ref_layer:
        r = b.raster(ref_layer) != 0
        ii, jj = np.nonzero(r)
        if len(ii):
            gx0, gy0, gx1, gy1 = gerber.bbox(files[ref_layer])
            ox, oy = g.x0 + jj.min() * g.res - gx0, g.y0 + ii.min() * g.res + gy1
            if abs(ox) > 0.2 or abs(oy) > 0.2:
                dx, dy = round(ox, 3), round(oy, 3)
    metrics = {}
    compared = 0

    def compare(layer, seen, what, sev):
        nonlocal compared
        G = files.get(layer)
        if G is None:
            if supplied:
                fs.append(finding("output_match", "error", f"no Gerber file for {layer} in {src}", layer=layer, kind="missing:" + layer))
            return
        if G["unsupported"]:
            fs.append(finding("output_match", "warning", f"{os.path.basename(G['path'])} uses {', '.join(G['unsupported'])}, which the reader "
                              f"does not draw: {layer} was not compared", layer=layer, kind="unsupported:" + layer))
            return
        img = gerber.raster(G, g, dx, dy)
        compared += 1
        blobs, n = _blobs(img ^ seen, g)
        metrics[layer] = {"file": os.path.basename(G["path"]), "file_mm2": round(float(img.sum()) * g.res ** 2, 2),
                          "checked_mm2": round(float(seen.sum()) * g.res ** 2, 2), "differences": n}
        for area, x, y in blobs:
            i, j = g.ij(x, y)
            side = "is in the file and was not checked" if img[int(i), int(j)] else "was checked and is not in the file"
            fs.append(finding("output_match", sev, f"{layer}: {what} of {area:.3f} mm2 {side}"
                              + (f" ({n} such places on this layer; the largest are listed)" if n > len(blobs) else ""),
                              layer=layer if layer in b.layers else b.layers[0 if layer.startswith("F.") else -1], x=x, y=y,
                              value=round(area, 4), unit="mm2", kind=f"diff:{layer}", basis="the fabrication files"))

    for li, lname in enumerate(b.layers):
        compare(lname, b.raster(li) != 0, "copper", "error")
    for mname, li, side in (("F.Mask", 0, 0), ("B.Mask", L - 1, 1)):
        seen = np.zeros(g.shape, bool)
        for p in b.pads:                                      # every opening, also apertures without copper of their own
            m = p.get("mask", {}).get(mname)
            if m and m["polys"]:
                geom.raster_rings(g, _rings(m["polys"]), seen)
        for s in b.snap.get("mask_shapes", []):
            if s["layer"] == mname:
                geom.raster_rings(g, _rings(s["polys"]), seen)
        dflt = b.snap.get("design", {}).get("mask_expansion", 0.0)
        for k, v in enumerate(b.snap["vias"]):
            if (v.get("tented") or [True, True])[side] is False and li in b.via_flashed[k]:
                rr = v["d"] / 2 + (v.get("mask_exp") if v.get("mask_exp") is not None else dflt)
                geom.raster_capsule(g, seen, v["x"], v["y"], v["x"], v["y"], rr)
        compare(mname, seen, "mask opening", "error")
    # holes
    have = [(v["x"], v["y"], v["drill"]) for v in b.snap["vias"]]
    have += [(p["x"], p["y"], min(p["drill"], p.get("drill_y") or p["drill"])) for p in b.pads if p["kind"] in ("pth", "npth") and p.get("drill", 0) > 0]
    if drills:
        bad = sorted({u for d in drills for u in d["unsupported"]})
        got = [((h["x"] + h.get("x2", h["x"])) / 2 + dx, -(h["y"] + h.get("y2", h["y"])) / 2 + dy, h["d"]) for d in drills for h in d["holes"]]
        if bad:
            fs.append(finding("output_match", "warning", f"the drill files use {', '.join(bad)}, which the reader does not follow: holes "
                              f"were not compared", kind="unsupported:drill"))
        elif have or got:
            compared += 1
            A, B = np.array(have).reshape(-1, 3), np.array(got).reshape(-1, 3)
            miss_f = miss_b = wrong = 0
            where = None
            if len(A) and len(B):
                d, k = cKDTree(B[:, :2]).query(A[:, :2])
                lost = d > 0.02
                miss_f = int(lost.sum())
                size = ~lost & (np.abs(B[k, 2] - A[:, 2]) > 0.006)
                wrong = int(size.sum())
                d2, _ = cKDTree(A[:, :2]).query(B[:, :2])
                miss_b = int((d2 > 0.02).sum())
                for arr, m in ((A, lost), (A, size), (B, d2 > 0.02)):
                    if where is None and m.any():
                        where = arr[np.flatnonzero(m)[0]]
            else:
                miss_f, miss_b = len(A), len(B)
            metrics["holes"] = {"in_files": len(B), "on_board": len(A), "not_in_files": miss_f, "not_on_board": miss_b, "other_size": wrong}
            if miss_f or miss_b or wrong:
                fs.append(finding("output_match", "error", f"the drill files differ from the board: {miss_f} holes of the board are not in "
                                  f"the files, {miss_b} holes in the files are not on the board, {wrong} have another diameter",
                                  x=None if where is None else float(where[0]), y=None if where is None else float(where[1]),
                                  value=miss_f + miss_b + wrong, kind="diff:drill", basis="the fabrication files"))
    elif supplied:
        fs.append(finding("output_match", "error", f"no drill file in {src}", kind="missing:drill"))
    # ---- feature by feature: every track, pad and via of the board as a draw or flash in the file, to
    # 1 um in position and size (the images above see only differences wider than two pixels)
    vec = _vector_match(b, files, dx, dy, fs)
    metrics["features"] = vec
    compared += bool(vec)
    # paste layers: images and apertures
    for pname, side in (("F.Paste", "F"), ("B.Paste", "B")):
        seen = np.zeros(g.shape, bool)
        for p in b.pads:
            ap = (p.get("paste_poly") or {}).get(pname)
            if ap:
                geom.raster_rings(g, _rings(ap), seen)
        if seen.any() or pname in files:
            compare(pname, seen, "paste aperture", "error")
    # outline
    _outline_match(b, files.get("Edge.Cuts"), dx, dy, fs, metrics, supplied)
    # placement and BOM files, when they are among the files
    _placement_match(b, src, dx, dy, fs, metrics, supplied)
    if not supplied and compared:
        fs.append(finding("output_match", "info", f"the files in {src} were plotted from the board just now. They agree with what was checked"
                          f"{'' if not any(f['severity'] == 'error' for f in fs) else ' except as listed'}, but that says nothing about files "
                          f"plotted earlier: send these files, or check the ones you will send with --gerbers DIR", kind="plotted"))
    params = {"files": src, "source": "supplied with --gerbers" if supplied else "plotted now by kicad-cli", "command": cmd,
              "offset_mm": [dx, dy], "resolution_mm": g.res,
              "detects": "tracks, pads and vias to 1 um in position and size; other copper, mask and paste differences wider than two pixels"}
    unchecked = "" if compared else "no file could be compared"
    if not supplied:
        unchecked = ("the files that will be sent were not given (--gerbers DIR): only the plotter and the reader were compared, "
                     "a self-test of this check")
    return dict(result("output_match", "MF20", "The files that are sent against the geometry that was checked", params, metrics, fs),
                checked=compared, unchecked=unchecked)


VEC_TOL = 0.001      # mm: position and size of a feature in the files against the board
POLY_ERR = 0.002     # mm: KiCad's polygons approximate arcs to 2 um (inside the shape)


def _expected_wh(p, poly_w, poly_h, margin=(0.0, 0.0)):
    """The extent a standard pad shape plots with, in board axes: exact from its size when the pad is at
    a multiple of 90 degrees, else its polygon's extent (then with the polygon's arc error allowed).
    Returns (w, h, tolerance)."""
    rot = (p.get("rot") or 0.0) % 180
    if p.get("shape") in ("circle", "rect", "oval", "roundrect", "trapezoid", "chamfered_rect") and min(rot, 180 - rot) < 1e-6 or abs(rot - 90) < 1e-6:
        w, h = p["sx"] + 2 * margin[0], p["sy"] + 2 * margin[1]
        if abs(rot - 90) < 1e-6:
            w, h = h, w
        return w, h, 3 * VEC_TOL
    return poly_w, poly_h, 3 * VEC_TOL + 2 * POLY_ERR


def _vector_match(b, files, dx, dy, fs):
    """Tracks, pads and vias of each copper layer against the draws and flashes of its Gerber file."""
    out = {}
    T, V = b.tracks, b.vias
    for li, lname in enumerate(b.layers):
        G = files.get(lname)
        if G is None or G["unsupported"]:
            continue
        D = np.array([(x0 + dx, -y0 + dy, x1 + dx, -y1 + dy, r) for x0, y0, x1, y1, r, dark in G["draws"] if dark]).reshape(-1, 5)
        F = np.array([(x + dx, -y + dy, w, h, a) for x, y, w, h, a, dark in G["flashes"] if dark]).reshape(-1, 5)
        st = {"tracks": 0, "tracks_missing": 0, "tracks_other_width": 0, "pads": 0, "pads_missing": 0, "pads_other_size": 0,
              "vias": 0, "vias_missing": 0}
        ex = {}
        # tracks (straight segments; arcs are cut into other chords by the plotter and the reader)
        on = np.flatnonzero((T["layer"] == li) & ~T["arc"]) if len(T["w"]) else []
        if len(on):
            if len(D):
                ends = cKDTree(np.r_[D[:, :2], D[:, 2:4]])
                n_d = len(D)
            for k in on:
                st["tracks"] += 1
                a = (T["x1"][k], T["y1"][k])
                c = (T["x2"][k], T["y2"][k])
                hit = None
                if len(D):
                    for j in ends.query_ball_point(a, VEC_TOL * 2):
                        jj = j % n_d
                        other = D[jj, 2:4] if j < n_d else D[jj, :2]
                        if math.hypot(other[0] - c[0], other[1] - c[1]) <= VEC_TOL * 2:
                            if hit is None or abs(D[jj, 4] - T["w"][k] / 2) < abs(D[hit, 4] - T["w"][k] / 2):
                                hit = jj
                if hit is None:
                    st["tracks_missing"] += 1
                    ex.setdefault("track_missing", (T["x1"][k], T["y1"][k], f"a {T['w'][k]:g} mm track"))
                elif abs(2 * D[hit, 4] - T["w"][k]) > VEC_TOL:
                    st["tracks_other_width"] += 1
                    ex.setdefault("track_width", (T["x1"][k], T["y1"][k], f"a {T['w'][k]:g} mm track plotted {2 * D[hit, 4]:g} mm wide"))
        # pads (standard shapes; custom pads may be plotted as regions and are left to the image comparison)
        if len(F):
            ft = cKDTree(F[:, :2])
        for p in b.pads:
            if lname not in p["layers"] or p.get("shape") == "custom" or not p["poly"]:
                continue
            st["pads"] += 1
            rings = _rings(p["poly"])
            P = np.vstack([np.asarray(r) for r in rings])
            w, h = P[:, 0].max() - P[:, 0].min(), P[:, 1].max() - P[:, 1].min()
            cands = ft.query_ball_point((p["x"], p["y"]), VEC_TOL * 2) if len(F) else []
            if not cands:
                st["pads_missing"] += 1
                ex.setdefault("pad_missing", (p["x"], p["y"], f"pad {p['ref']}.{p['num']}"))
                continue
            w, h, tol = _expected_wh(p, w, h)
            best = min(cands, key=lambda j: abs(F[j, 2] - w) + abs(F[j, 3] - h))
            if abs(F[best, 2] - w) > tol or abs(F[best, 3] - h) > tol:
                st["pads_other_size"] += 1
                ex.setdefault("pad_size", (p["x"], p["y"], f"pad {p['ref']}.{p['num']} ({w:.3f} x {h:.3f} mm) plotted "
                                                           f"{F[best, 2]:.3f} x {F[best, 3]:.3f} mm"))
        # vias on the layers they are flashed on
        for k in range(len(V["x"])):
            if li not in b.via_flashed[k]:
                continue
            st["vias"] += 1
            cands = ft.query_ball_point((V["x"][k], V["y"][k]), VEC_TOL * 2) if len(F) else []
            if not any(abs(F[j, 2] - V["d"][k]) <= 3 * VEC_TOL for j in cands):
                st["vias_missing"] += 1
                ex.setdefault("via", (V["x"][k], V["y"][k], f"a {V['d'][k]:g} mm via"))
        out[lname] = st
        bad = st["tracks_missing"] + st["tracks_other_width"] + st["pads_missing"] + st["pads_other_size"] + st["vias_missing"]
        if bad:
            x, y, what = next(iter(ex.values()))
            fs.append(finding("output_match", "error", f"{lname}: {bad} features of the board differ in the file by more than "
                              f"{VEC_TOL * 1000:g} um ({st['tracks_missing']} tracks missing, {st['tracks_other_width']} tracks of another width, "
                              f"{st['pads_missing']} pads missing, {st['pads_other_size']} pads of another size, {st['vias_missing']} vias "
                              f"missing); for example {what}", layer=lname, x=float(x), y=float(y), value=bad, kind="features:" + lname,
                              basis="the fabrication files"))
    # paste apertures (standard shapes; custom apertures are left to the image comparison)
    for pname in ("F.Paste", "B.Paste"):
        G = files.get(pname)
        if G is None or G["unsupported"]:
            continue
        F = np.array([(x + dx, -y + dy, w, h, a) for x, y, w, h, a, dark in G["flashes"] if dark]).reshape(-1, 5)
        ft = cKDTree(F[:, :2]) if len(F) else None
        n = miss = size = 0
        ex = None
        for p in b.pads:
            ap = (p.get("paste_poly") or {}).get(pname)
            if not ap or p.get("shape") == "custom":
                continue
            n += 1
            P = np.vstack([np.asarray(r) for r in _rings(ap)])
            w, h = P[:, 0].max() - P[:, 0].min(), P[:, 1].max() - P[:, 1].min()
            w, h, tol = _expected_wh(p, w, h, tuple((p.get("paste") or {}).get(pname) or (0.0, 0.0)))
            cands = ft.query_ball_point((p["x"], p["y"]), VEC_TOL * 2) if ft is not None else []
            if not cands:
                miss += 1
                ex = ex or (p["x"], p["y"], f"the aperture of {p['ref']}.{p['num']}")
            elif min(max(abs(F[j, 2] - w), abs(F[j, 3] - h)) for j in cands) > tol:
                size += 1
                ex = ex or (p["x"], p["y"], f"the aperture of {p['ref']}.{p['num']} ({w:.3f} x {h:.3f} mm)")
        out[pname] = {"apertures": n, "missing": miss, "other_size": size}
        if miss or size:
            fs.append(finding("output_match", "error", f"{pname}: {miss} paste apertures missing and {size} of another size in the file; "
                              f"for example {ex[2]}", layer=b.layers[0 if pname == "F.Paste" else -1], x=float(ex[0]), y=float(ex[1]),
                              value=miss + size, kind="features:" + pname, basis="the fabrication files"))
    return out


def _outline_match(b, G, dx, dy, fs, metrics, supplied):
    """The board outline in the Edge.Cuts file: every draw on the outline, and its length that of the outline."""
    rings = b.outline_rings()
    if not rings:
        return
    if G is None:
        if supplied:
            fs.append(finding("output_match", "error", "no Edge.Cuts (outline) file among the files", kind="missing:Edge.Cuts"))
        return
    segs = geom.ring_segments(rings)
    per = float(np.hypot(segs[:, 2] - segs[:, 0], segs[:, 3] - segs[:, 1]).sum())
    D = [(x0 + dx, -y0 + dy, x1 + dx, -y1 + dy) for x0, y0, x1, y1, r, dark in G["draws"] if dark]
    length, off = 0.0, 0
    for x0, y0, x1, y1 in D:
        length += math.hypot(x1 - x0, y1 - y0)
        for x, y in ((x0, y0), (x1, y1), ((x0 + x1) / 2, (y0 + y1) / 2)):
            if float(geom.point_seg_dist(x, y, segs[:, 0], segs[:, 1], segs[:, 2], segs[:, 3]).min()) > 0.005:
                off += 1
                break
    metrics["outline"] = {"board_mm": round(per, 3), "file_mm": round(length, 3), "draws_off_outline": off}
    if off or abs(length - per) > max(0.05, 0.002 * per):
        fs.append(finding("output_match", "error", f"the outline in the file differs from the board's: {length:.2f} mm drawn against "
                          f"{per:.2f} mm, {off} draws off the board outline", kind="outline", value=round(length - per, 3),
                          basis="the fabrication files"))


def _placement_match(b, src, dx, dy, fs, metrics, supplied):
    """A placement (pos) file and a BOM among the files: references, sides and positions of the fitted parts."""
    import csv
    pos = bom = None
    for name in sorted(os.listdir(src)):
        low = name.lower()
        path = os.path.join(src, name)
        if not low.endswith((".csv", ".pos")):
            continue
        try:
            head = open(path, encoding="utf-8", errors="replace").readline().lower()
        except OSError:
            continue
        if "posx" in head or "mid x" in head:
            pos = pos or path
        elif "ref" in head or "designator" in head:
            bom = bom or path
    fitted = {r for r, fp in b.fps.items() if not fp.get("dnp") and fp.get("in_pos", fp.get("placed", True))
              and any(p["kind"] in ("smd", "pth") for p in b.pads_by_ref.get(r, []))}
    if pos:
        rows = list(csv.DictReader(open(pos, encoding="utf-8", errors="replace")))
        key = {k.lower(): k for k in (rows[0].keys() if rows else [])}
        got = {}
        for r in rows:
            try:
                got[r[key.get("ref", "Ref")]] = (float(r[key.get("posx", "PosX")]), float(r[key.get("posy", "PosY")]),
                                                 float(r[key.get("rot", "Rot")]), r[key.get("side", "Side")].strip().lower())
            except (KeyError, ValueError):
                continue
        common = [r for r in got if r in b.fps]
        # the file's origin and y axis: from the median offset of the common parts
        if common:
            ox = np.median([got[r][0] - b.fps[r]["x"] for r in common])
            oy = np.median([got[r][1] + b.fps[r]["y"] for r in common])
            moved = [r for r in common if math.hypot(got[r][0] - ox - b.fps[r]["x"], -(got[r][1] - oy) - b.fps[r]["y"]) > 0.002]
            side = [r for r in common if (got[r][3].startswith("b")) != (b.fps[r]["side"] == "B")]
        else:
            moved, side = [], []
        missing = sorted(fitted - set(got))
        extra = sorted(set(got) - set(b.fps))
        metrics["placement"] = {"file": os.path.basename(pos), "parts": len(got), "missing": missing[:20], "not_on_board": extra[:20],
                                "moved": moved[:20], "other_side": side[:20]}
        if missing or extra or moved or side:
            fs.append(finding("output_match", "error", f"the placement file {os.path.basename(pos)} differs from the board: "
                              f"{len(missing)} fitted parts missing, {len(extra)} not on the board, {len(moved)} at another position, "
                              f"{len(side)} on the other side" + (f" (e.g. {(missing + extra + moved + side)[0]})"),
                              kind="placement", value=len(missing) + len(extra) + len(moved) + len(side), basis="the fabrication files"))
    elif supplied:
        fs.append(finding("output_match", "info", "no placement (pos) file among the files: placement was not compared", kind="no_pos"))
    if bom:
        rows = list(csv.reader(open(bom, encoding="utf-8", errors="replace")))
        head = [h.lower() for h in rows[0]] if rows else []
        col = next((i for i, h in enumerate(head) if h in ("reference", "references", "designator", "ref", "refs")), None)
        refs = set()
        for r in rows[1:]:
            if col is not None and col < len(r):
                refs |= {x.strip() for x in re.split(r"[,\s]+", r[col]) if x.strip()}
        want = {r for r, fp in b.fps.items() if not fp.get("dnp") and fp.get("in_bom", True) and b.pads_by_ref.get(r)}
        metrics["bom"] = {"file": os.path.basename(bom), "refs": len(refs), "missing": sorted(want - refs)[:20],
                          "not_fitted_or_absent": sorted(refs - want)[:20]}
        if want - refs or refs - want:
            fs.append(finding("output_match", "error", f"the BOM {os.path.basename(bom)} differs from the board's fitted parts: "
                              f"{len(want - refs)} missing, {len(refs - want)} not fitted or not on the board", kind="bom",
                              value=len(want - refs) + len(refs - want), basis="the fabrication files"))
    elif supplied:
        fs.append(finding("output_match", "info", "no BOM among the files: the BOM was not compared", kind="no_bom"))


# ------------------------------------------------------------------------------------------------
# MF21 (ranking): copper share per layer

def rank_copper(ctx):
    from . import thermal
    b = ctx.b
    inside = thermal._inside(ctx)
    area = float(inside.sum()) or 1.0
    share = [float(((b.raster(l) != 0) & inside).sum()) / area for l in range(len(b.layers))]
    L = len(share)
    rows = [{"layer": b.layers[l], "copper_pct": round(100 * share[l], 1), "thickness_um": round(b.stack.copper[l]["t"] * 1000, 1),
             "mirror": b.layers[L - 1 - l], "difference_pct": round(100 * abs(share[l] - share[L - 1 - l]), 1)} for l in range(L)]
    return {"rule": "copper", "spec": "MF21", "title": "Copper share of each layer and of each pair of layers about the centre",
            "params": {}, "findings": [], "rows": rows,
            "summary": {"largest_pair_difference_pct": max((r["difference_pct"] for r in rows), default=0.0)}}


# ------------------------------------------------------------------------------------------------
JUDGE_KINDS = {
    "observed": ("file", "sha256"),        # defects counted on built boards: the --observed file and its SHA-256
    "assembler": ("document", "date"),     # the assembler's written defect data or limit
    "standard": ("document",),             # a published standard
}


def judge_provenance(cfg, base):
    """`[dfm.judge]` gives a ranking a limit only with a source that can be checked: a table naming its
    kind (observed, assembler or standard) and what the kind requires. A free-text source is refused.
    The source is replaced by a description that findings carry as their basis."""
    from .snapshot import InputError
    j = cfg["dfm"]["judge"]
    limits = [k for k in j if k != "source"]
    if not limits:
        return
    src = j.get("source")
    if not isinstance(src, dict):
        raise InputError("[dfm.judge] needs `source` as a table that can be checked, e.g. "
                         "source = { kind = \"observed\", file = \"defects.toml\", sha256 = \"...\" } "
                         "(kinds: observed: file, sha256; assembler: document, date; standard: document); "
                         f"a free-text source ({src!r}) is not enough to give a ranking a limit")
    kind = src.get("kind")
    if kind not in JUDGE_KINDS:
        raise InputError(f"[dfm.judge] source kind {kind!r}: one of {', '.join(JUDGE_KINDS)}")
    miss = [k for k in JUDGE_KINDS[kind] if not src.get(k)]
    if miss:
        raise InputError(f"[dfm.judge] source of kind {kind} needs {', '.join(miss)}")
    if kind == "observed":
        path = src["file"] if os.path.isabs(src["file"]) else os.path.join(base, src["file"])
        if not os.path.exists(path):
            raise InputError(f"[dfm.judge] source file {path} not found")
        h = hashlib.sha256(open(path, "rb").read()).hexdigest()
        if not h.startswith(str(src["sha256"]).lower()):
            raise InputError(f"[dfm.judge] source file {path} does not match the SHA-256 given (it is {h[:16]}...)")
        desc = f"observed defects, {os.path.basename(path)} (sha256 {h[:12]})"
    elif kind == "assembler":
        desc = f"assembler: {src['document']}, {src['date']}"
    else:
        desc = f"standard: {src['document']}"
    j["source_table"] = dict(src)
    j["source"] = desc


MIN_DEFECTIVE = 5          # fewer defective parts than this: no order score is given


def score_observed(path, rankings, boot=1000):
    """How well each ranking orders the parts that showed a defect: the share of (defective, sound) pairs
    it puts in the right order (0.5 = chance, 1 = every defective part above every sound one), with a
    bootstrap 95 % interval and the counts it rests on. The file must say how many boards were built;
    below MIN_DEFECTIVE defective parts a ranking gets no score ("insufficient data")."""
    data = tomllib.loads(open(path, "rb").read().decode())
    obs = data.get("defect", [])
    built = data.get("boards_built")
    out = {"file": os.path.abspath(path), "defects": len(obs), "boards_built": built, "scores": []}
    if not isinstance(built, int) or built < 1:
        out["error"] = "the file must give boards_built (the number of boards assembled) as a positive integer"
        return out
    kinds = {}
    for d in obs:
        kinds.setdefault(str(d.get("kind", "any")).lower(), {})[d.get("ref")] = int(d.get("count", 1))
    out["defect_count"] = sum(sum(v.values()) for v in kinds.values())
    cols = []
    if "chips" in rankings:
        base = {}
        for x in rankings["chips"]["rows"]:
            for ref in x["refs"]:
                base[ref] = x
        cols.append(("chips: time between the joints", {r: x["dt_s"] for r, x in base.items()}))
        cols.append(("chips: copper leaving the pads, ratio", {r: max(q["entry_mm"] for q in x["pads"]) / max(min(q["entry_mm"] for q in x["pads"]), 0.1)
                                                              for r, x in base.items()}))
        cols.append(("chips: via into a plane at one pad only", {r: float(bool(x["pads"][0]["plane_vias"]) != bool(x["pads"][1]["plane_vias"]))
                                                                for r, x in base.items()}))
        cols.append(("chips: toe of the pads", {r: x["land"]["toe_mm"] for r, x in base.items()}))
    if "reflow" in rankings:
        cols.append(("reflow: below the median peak", {x["ref"]: -x["d_peak_c"] for x in rankings["reflow"]["rows"]}))
    rng = np.random.default_rng(0)

    def auc(pos, neg):
        return float(((pos[:, None] > neg[None, :]).sum() + 0.5 * (pos[:, None] == neg[None, :]).sum()) / (len(pos) * len(neg)))
    for kind, refs in kinds.items():
        for name, score in cols:
            pos = np.array([score[r] for r in refs if r in score])
            neg = np.array([v for r, v in score.items() if r not in refs])
            if not len(pos) or not len(neg):
                continue
            row = {"kind": kind, "ranking": name, "defective_ranked": len(pos), "sound": len(neg), "boards_built": built,
                   "defects": sum(refs[r] for r in refs if r in score)}
            if len(pos) < MIN_DEFECTIVE:
                row.update(order=None, note=f"insufficient data: {len(pos)} defective parts, at least {MIN_DEFECTIVE} needed")
            else:
                bs = [auc(rng.choice(pos, len(pos)), rng.choice(neg, len(neg))) for _ in range(boot)]
                row.update(order=round(auc(pos, neg), 2), ci95=[round(float(np.percentile(bs, 2.5)), 2), round(float(np.percentile(bs, 97.5)), 2)])
            out["scores"].append(row)
    return out


def _ensure_validation(cache, log):
    from . import etch, reflow, validate
    for solver, mod in (("etch", etch), ("reflow", reflow)):
        if validate.status(cache, solver) != "validated":
            validate.register(cache, solver, mod.validation_cases())
            log(f"  measurement and solver cases, {solver}: {'pass' if validate.status(cache, solver) == 'validated' else 'FAIL'}")


def run(board, cfg_path=None, fab=None, checks=None, out_dir=None, waivers=None, json_out=None, markers=True, log=print, cross=True,
        gerbers=None, observed=None):
    from . import etch, evaluate, reflow, report as report_mod, run as runmod, validate
    t0 = time.time()
    b, cfg, ctx = runmod.load(board, cfg_path)
    b.label_nonet_pads()                 # copper without a net (fiducials, mounting pads) is copper too
    ctx.__init__(b, cfg)
    if fab:
        cfg.data["dfm"]["fab"] = fab
    gerbers = gerbers or cfg["dfm"]["gerbers"] or None
    fabr = fab_mod.resolve(fab_mod.load(cfg["dfm"]["fab"]), b, cfg)
    judge_provenance(cfg, os.path.dirname(cfg.path) if cfg.path else os.getcwd())
    # the board read a second time by an independent parser: a disagreement means the input is not known
    xchk = None
    if cross and b.snap.get("reader") == "kicad-monkey":
        from . import snapshot as snap_mod
        xchk = snap_mod.cross_check(b.snap, runmod.CACHE)
        log(f"cross-check ({xchk['reader']}): {sum(xchk['compared'].values())} items compared, {xchk['n_problems']} differences")
    want = list(GROUPS) if not checks else [c.strip().lower() for c in checks]
    unknown = [c for c in want if c not in GROUPS]
    if unknown:
        raise ValueError(f"unknown dfm checks {unknown}; choose from {', '.join(GROUPS)}")
    od = out_dir or os.path.join(runmod.RUNS, f"{b.name}_{b.sha[:8]}", "dfm")
    os.makedirs(od, exist_ok=True)
    log(f"openpcbsi dfm: {b.path}")
    log(f"  sha256 {b.sha[:16]}  {len(b.layers)} layers  fabricator profile {fabr.name} ({os.path.basename(fabr.path)}, read {fabr.retrieved})")
    _ensure_validation(runmod.CACHE, log)
    verified = {s: validate.status(runmod.CACHE, s) == "validated" for s in ("etch", "reflow")}
    res_hi = b.grid().res / max(1, int(cfg["dfm"]["acid"]["hires_div"]))
    bga, bga_refs = bga_field(ctx)
    g = b.grid()
    in_bga = lambda x, y: bool(g.sample(bga, [x], [y])[0])
    rules, ranks, took = {}, {}, {}
    exposure = {}

    def timed(name, fn):
        t = time.time()
        out = fn()
        took[name] = round(time.time() - t, 1)
        return out

    if "audit" in want:
        log("MF19 KiCad's rules against the fabricator's:")
        rules["rule_audit"] = timed("audit", lambda: rule_audit(ctx, fabr))
        for r in rules["rule_audit"]["rows"]:
            log(f"  {r['rule']:28s} KiCad {r['kicad_mm']}  {fabr.name} {' / '.join(f'{v:g}' for v in r['fab'].values())}  {r['verdict']}")
    if "copper" in want:
        log(f"MF3/MF12/MF2/MF1 copper: widths, spacing, same-net gaps, acute corners (measured at {res_hi * 1000:g} um):")
        rules["fab_track"] = timed("track", lambda: rule_fab_track(ctx, fabr, bga))
        o = timed("copper", lambda: etch.run(ctx, fabr, in_bga, log=log))
        for k in ("fab_space", "same_net_gap"):
            rules[k] = dict(o[k], checked=len(b.tracks["w"]) + len(b.pads))
        ranks["acid_trap"] = o["acid_trap"]
        j = cfg["dfm"]["judge"]
        if j.get("source") and j.get("acute_deg") is not None:
            for x in o["acid_trap"]["rows"]:
                if x["track_mm"] is not None and x["track_mm"] <= j.get("acute_track_mm", 1e9) and x["angle_deg"] < j["acute_deg"]:
                    o["acid_trap"]["findings"].append(finding(
                        "acid_trap", "warning", f"corner of {x['angle_deg']:.0f} deg on a {x['track_mm']:g} mm track, sharper than "
                        f"{j['acute_deg']:g} deg (board config: {j['source']})", x["net"], x["layer"], x["x"], x["y"], value=x["angle_deg"],
                        limit=j["acute_deg"], unit="deg", basis=f"board config: {j['source']}"))
        ranks["copper"] = rank_copper(ctx)
    if "drill" in want:
        log("MF13 holes:")
        rules["fab_drill"] = timed("drill", lambda: rule_fab_drill(ctx, fabr))
        log(f"  {rules['fab_drill']['checked']} holes, {len(rules['fab_drill']['findings'])} findings; {rules['fab_drill']['metrics'].get('board')}")
    if "mask" in want:
        log("MF4/MF5/MF14 solder mask and vias at pads:")
        rules["mask_web"] = timed("mask_web", lambda: rule_mask_web(ctx, fabr))
        rules["mask_exposed"], exposure = timed("mask_exposed", lambda: rule_mask_exposed(ctx, fabr, res_hi))
        rules["via_pad"] = timed("via_pad", lambda: rule_via_pad(ctx, fabr))
        log(f"  {rules['mask_web']['checked']} openings: {len(rules['mask_web']['findings'])} web, {len(rules['mask_exposed']['findings'])} exposure, "
            f"{len(rules['via_pad']['findings'])} via findings (vias {rules['via_pad']['params']['via_cover']})")
    if "paste" in want:
        log("MF15 stencil apertures:")
        rules["paste"] = timed("paste", lambda: rule_paste(ctx, fabr))
        log(f"  {rules['paste']['metrics'].get('stencil') or rules['paste']['params'].get('stencil_source') or 'no paste apertures'}")
    if "parts" in want:
        log("MF16/MF17/MF18 parts:")
        rules["footprint"] = timed("footprint", lambda: rule_footprint(ctx))
        rules["part_spacing"] = timed("part_spacing", lambda: rule_part_spacing(ctx, fabr))
        rules["assembly"] = timed("assembly", lambda: rule_assembly(ctx, fabr))
        log(f"  {rules['footprint']['checked']} chip land patterns, {rules['part_spacing']['checked']} part pairs, "
            f"{rules['assembly']['metrics']['board']}")
    if "output" in want:
        log(f"MF20 the files that are sent ({'supplied: ' + gerbers if gerbers else 'plotted now'}):")
        rules["output_match"] = timed("output", lambda: rule_output(ctx, fabr, od, gerbers, log))
        m = rules["output_match"]["metrics"]
        log(f"  {rules['output_match']['checked']} files compared; differences: "
            + (", ".join(f"{k} {v['differences']}" for k, v in m.items() if "differences" in v and v["differences"]) or "none")
            + (f"; holes {m['holes']}" if "holes" in m else ""))
    internal = None
    if "reflow" in want or "chips" in want:
        log("MF7 reflow ranking (model):")
        rb = timed("reflow", lambda: reflow.board(ctx, fabr, od, log=log))
        internal = rb.pop("_internal")
        if "reflow" in want:
            ranks["reflow"] = rb
    if "chips" in want:
        log("MF8 chip ranking (model):")
        ranks["chips"] = timed("chips", lambda: reflow.chips(ctx, fabr, internal, exposure, log=log))
    internal = None
    # ---- findings in the G19 model
    findings = []
    for name in ORDER:
        r = rules.get(name) or ranks.get(name)
        if not r:
            continue
        tier = "ranking" if name in RANKINGS else "rule"
        for f in r.get("findings", []):
            e = evaluate._f(name, "MF", REQ[name], evaluate.SEV.get(f.get("severity"), "info"), f.get("msg", ""),
                            [f.get("net"), f.get("aggressor")], [f.get("part")], f.get("layer") or None, f.get("x"), f.get("y"),
                            f.get("value"), f.get("limit"), f.get("unit"), extra=f.get("kind"))
            e.update(tier=tier, basis=f.get("basis") or ("geometry" if tier == "rule" else ""), field=f.get("field") or "")
            # a limit whose reading the profile had to interpret (the fabricator does not say what the
            # figure measures) can warn, never fail
            if e["severity"] == "fail" and "reading interpreted" in e["basis"]:
                e["severity"] = "warn"
                e["capped"] = "a fail on a limit whose reading is interpreted, not stated, is a warning"
            if e["severity"] in ("fail", "warn"):
                e["hint"] = evaluate.hint(name, f)
            findings.append(e)
    seen = {}
    for f in findings:
        n = seen.get(f["id"], 0)
        seen[f["id"]] = n + 1
        if n:
            f["id"] = f"{f['id']}.{n}"
    ws, whash = evaluate.load_waivers(waivers)
    wstat = evaluate.apply_waivers(findings, ws)
    # ---- scorecard
    judge = cfg["dfm"]["judge"]
    per = []
    for name in ORDER:
        grp = next(k for k, v in GROUPS.items() if name in v)
        r = rules.get(name) or ranks.get(name)
        fs = [f for f in findings if f["check"] == name]
        c = {"check": name, "requirement": REQ[name], "tier": "ranking" if name in RANKINGS else "rule",
             "fail": sum(1 for f in fs if f["severity"] == "fail" and not f["waived"]),
             "warn": sum(1 for f in fs if f["severity"] == "warn" and not f["waived"]),
             "info": sum(1 for f in fs if f["severity"] == "info" and not f["waived"]),
             "waived": sum(1 for f in fs if f["waived"]), "seconds": took.get(grp)}
        if r is None:
            c.update(status="not run", note="not selected (--checks)" if grp not in want else "")
        elif name in RANKINGS:
            n = len(r.get("rows", []))
            judged = bool(judge.get("source")) and name in ("reflow", "chips", "acid_trap")
            c["status"] = "warn" if c["warn"] else ("ranked" if n else "nothing to rank")
            c["note"] = (f"{n} rows; " if n else "") + (f"judged with limits from the board config ({judge['source']})" if judged and c["warn"]
                                                       else "not judged: no published limit, not compared with fabricated boards")
            if name in ("reflow", "chips") and not verified["reflow"]:
                c["note"] += "; SOLVER CASES FAIL"
        else:
            c["status"] = ("fail" if c["fail"] else "could not check" if r.get("unchecked") else "warn" if c["warn"] else
                           "pass" if r.get("checked") else "nothing to check")
            c["note"] = r.get("unchecked") or r.get("note") or (f"{r.get('checked')} items" if r.get("checked") else "")
            if name in ("fab_space", "same_net_gap") and not verified["etch"]:
                c["status"], c["note"] = "could not check", "the measurement cases fail"
        per.append(c)
    live = [f for f in findings if not f["waived"]]
    counts = {s: sum(1 for f in live if f["severity"] == s) for s in ("fail", "warn", "info")}
    counts["waived"] = sum(1 for f in findings if f["waived"])
    code = 2 if counts["fail"] else (1 if counts["warn"] and cfg.data.get("evaluate", {}).get("warnings_exit", True) else 0)
    status = "fail" if counts["fail"] else "warn" if counts["warn"] else "pass"
    # order options the verdict relies on without the board config stating them
    ran = {c["check"] for c in per if c["status"] not in ("not run",)}
    assumptions = [{"check": n, "assumption": r["assumption"]} for n, r in rules.items() if r.get("assumption") and n in ran]
    if status == "pass" and assumptions:
        status = "pass with assumptions"
        code = max(code, 1)
    checklist = order_checklist(ctx, fabr, rules)
    if xchk and xchk["n_problems"]:
        status, code = "could not check", 3
        log(f"the two readings of the board differ ({xchk['n_problems']}): the review cannot rely on its input; e.g. {xchk['problems'][0]}")
    obs = None
    if observed:
        obs = score_observed(observed, ranks)
        if obs.get("error"):
            log(f"observed defects: {os.path.basename(observed)}: {obs['error']}")
        else:
            log(f"observed defects ({obs['defects']} parts, {obs['boards_built']} boards built, in {os.path.basename(observed)}): how each "
                f"ranking orders them (0.5 = chance, 1 = perfect)")
        for sc in obs["scores"]:
            val = f"{sc['order']:.2f} (95 % {sc['ci95'][0]:.2f} to {sc['ci95'][1]:.2f})" if sc["order"] is not None else sc["note"]
            log(f"  {sc['kind']:12s} {sc['ranking']:42s} {val}  ({sc['defective_ranked']} defective, {sc['sound']} sound)")
    rep = {"checks": {"drc": [r for r in list(rules.values()) + list(ranks.values())]}}
    images = {}
    try:
        images = report_mod.layer_images(rep, b, ctx, od)
    except Exception as e:
        log(f"layer images not written: {e}")
    marker = None
    if markers:
        try:
            marker, n = report_mod.marker_board(rep, b.path, od)
            log(f"marker board: {os.path.join(od, marker)} ({n} markers on User.9)")
        except Exception as e:
            log(f"marker board not written: {e}")
    unjudged = fabr.not_judged()
    ev = {"schema": SCHEMA, "tool": "openpcbsi", "version": __import__("openpcbsi").__version__,
          "date": datetime.datetime.now().isoformat(timespec="seconds"), "board": b.path, "board_sha256": b.sha,
          "kicad": b.snap.get("kicad"), "config": cfg.path or "detected from the board", "config_hash": cfg.hash(),
          "fab": fabr.describe(), "limits_not_judged": unjudged if not checks else None,
          "waivers_file": os.path.abspath(waivers) if waivers else None, "waivers_hash": whash,
          "solver_versions": {k: validate.SOLVER_VERSIONS[k] for k in ("etch", "reflow", "thermal")}, "cases_pass": verified,
          "reader": {"primary": b.snap.get("reader") or "pcbnew", "kicad": b.snap.get("kicad"), "cross_check": xchk},
          "exit_code": code, "status": status, "assumptions": assumptions, "order_checklist": checklist, "checks": per, "waivers": wstat, "counts": counts, "findings": findings,
          "rules": list(rules.values()), "rankings": list(ranks.values()), "bga_fields": bga_refs, "not_checked": NOT_CHECKED,
          "observed": obs, "judge": judge or None, "images": images, "marker_board": marker, "seconds": round(time.time() - t0, 1)}
    path = os.path.join(od, "dfm.json")
    for pth in [path] + ([json_out] if json_out else []):
        with open(pth, "w") as fh:
            json.dump(ev, fh, indent=1, default=runmod._jd)
    write_html(ev, os.path.join(od, "dfm.html"))
    log("")
    log(f"rules ({fabr.name} profile; these give the verdict):")
    for c in per:
        if c["tier"] == "rule":
            log(f"  {c['status']:16s} {c['requirement']:5s} {c['check']:13s} fail {c['fail']:3d}  warn {c['warn']:3d}  info {c['info']:3d}"
                + (f"  waived {c['waived']}" if c["waived"] else "") + (f"   {c['note']}" if c["note"] else ""))
    log("rankings (models and measurements without a published limit; not part of the verdict):")
    for c in per:
        if c["tier"] == "ranking":
            log(f"  {c['status']:16s} {c['requirement']:5s} {c['check']:13s} {c['note']}")
    for a in assumptions:
        log(f"assumed ({a['check']}): {a['assumption']}")
    if checklist:
        log("order checklist:")
        for c in checklist:
            log(f"  {c['item']}: {c['state']}")
    log(f"verdict of the rules: {status}  ({counts['fail']} fail, {counts['warn']} warn, {counts['info']} info, {counts['waived']} waived)")
    log(f"dfm:  {path}")
    log(f"html: {os.path.join(od, 'dfm.html')}")
    log(f"exit code {code}  ({ev['seconds']} s)")
    return code, ev


# ------------------------------------------------------------------------------------------------
def write_html(ev, path):
    from .evaluate import CSS
    e = html.escape
    st = lambda s: f'<span class="{ {"could not check": "unvalidated", "ranked": "mut", "nothing to rank": "mut"}.get(s, s.replace(" ", ""))}">{e(s)}</span>'
    num = lambda v, d=2: "" if v is None else f"{v:.{d}f}"
    fabd = ev["fab"]
    rules = {r["rule"]: r for r in ev["rules"]}
    ranks = {r["rule"]: r for r in ev["rankings"]}
    H = [f"<!doctype html><html><head><meta charset=utf-8><meta name=viewport content='width=device-width,initial-scale=1'>"
         f"<title>openpcbsi dfm</title><style>{CSS} img{{max-width:100%;height:auto}} .wrap{{overflow-x:auto}}</style></head><body>"]
    H.append(f"<h1>Manufacturing review: {e(os.path.basename(ev['board']))} &mdash; rules {st(ev['status'])}</h1>")
    H.append(f"<p class=mut>{e(ev['board'])}<br>sha256 {e(ev['board_sha256'][:16])} · config {e(str(ev['config']))} · "
             f"fabricator profile {e(fabd['name'])} (read {e(str(fabd['retrieved']))}, hash {e(str(fabd['profile_hash']))}) · "
             f"{e(ev['date'])} · exit code {ev['exit_code']} · {ev['seconds']} s</p>")
    rd = ev.get("reader") or {}
    xc = rd.get("cross_check")
    H.append(f"<p class=mut>Board read by {e(str(rd.get('kicad') or rd.get('primary')))}"
             + (f"; read again by {e(xc['reader'])}: {sum(xc['compared'].values())} items compared, "
                f"{'no difference' if not xc['n_problems'] else str(xc['n_problems']) + ' differences (the review cannot rely on its input)'}; "
                f"not cross-checked: {e(', '.join(xc['not_cross_checked']))}" if xc else "; not read a second time") + "</p>")
    H.append("<p>This is a review of the layout against published limits, not a sign-off. The <b>rules</b> judge drawn geometry against a "
             "limit that the fabricator publishes, that an industry standard gives or that the board config sets; they alone give the verdict "
             "and the exit code. The <b>rankings</b> order places on this board by a measurement or a model that no published limit covers and "
             "that has not been compared with fabricated boards; they never pass or fail.</p>")
    H.append("<h2>Rules</h2><div class=wrap><table><tr><th>Req.</th><th>Check</th><th>Status</th><th class=n>fail</th><th class=n>warn</th>"
             "<th class=n>info</th><th class=n>waived</th><th>Note</th></tr>")
    for c in ev["checks"]:
        if c["tier"] == "rule":
            H.append(f"<tr><td>{e(c['requirement'])}</td><td>{e(c['check'])}</td><td>{st(c['status'])}</td><td class=n>{c['fail']}</td>"
                     f"<td class=n>{c['warn']}</td><td class=n>{c['info']}</td><td class=n>{c['waived']}</td><td class=mut>{e(str(c.get('note') or ''))}</td></tr>")
    H.append("</table></div><h2>Rankings</h2><div class=wrap><table><tr><th>Req.</th><th>Ranking</th><th>Status</th><th>Note</th></tr>")
    for c in ev["checks"]:
        if c["tier"] == "ranking":
            H.append(f"<tr><td>{e(c['requirement'])}</td><td>{e(c['check'])}</td><td>{st(c['status'])}</td><td class=mut>{e(str(c.get('note') or ''))}</td></tr>")
    H.append("</table></div>")
    if ev.get("assumptions"):
        H.append("<h2>Assumed, not stated</h2><p>The verdict relies on these order options, which the board config does not state:</p><ul>"
                 + "".join(f"<li>{e(a['assumption'])} <span class=mut>({e(a['check'])})</span></li>" for a in ev["assumptions"]) + "</ul>")
    if ev.get("order_checklist"):
        H.append("<h2>Order checklist</h2><p>State these in the order; otherwise the fabricator chooses.</p><div class=wrap><table>"
                 "<tr><th>Item</th><th>What the order needs</th><th>Source</th></tr>"
                 + "".join(f"<tr><td>{e(c['item'])}</td><td>{e(c['state'])}</td><td class=mut>{e(c['basis'])}</td></tr>" for c in ev["order_checklist"])
                 + "</table></div>")
    H.append("<h2>Not checked</h2><ul>" + "".join(f"<li>{e(x)}</li>" for x in ev["not_checked"]) + "</ul>")
    rank = {"fail": 0, "warn": 1}
    live = sorted((f for f in ev["findings"] if f["severity"] in rank),
                  key=lambda f: (bool(f["waived"]), rank[f["severity"]], f["requirement"], f["check"]))
    H.append(f"<h2>Findings: fail and warn ({len(live)})</h2><div class=wrap><table><tr><th>ID</th><th>Sev.</th><th>Check</th><th>Reason</th>"
             "<th>Limit from</th><th>Nets / refs</th><th>Layer</th><th>x, y (mm)</th><th>Note</th></tr>")
    for f in live[:2000]:
        note = "waived: " + f["waived"]["reason"] if f["waived"] else (f.get("field") and "in a " + f["field"] + " field" or "")
        xy = f"{f['x']:.2f}, {f['y']:.2f}" if f["x"] is not None and f["y"] is not None else ""
        H.append(f"<tr><td><code>{e(f['id'])}</code></td><td>{st(f['severity'])}</td><td>{e(f['check'])} <span class=mut>{e(f['requirement'])}</span></td>"
                 f"<td>{e(f['reason'])}</td><td class=mut>{e(f.get('basis') or '')}</td><td>{e(', '.join(f['nets'] + f['refs']))}</td>"
                 f"<td>{e(f['layer'] or '')}</td><td>{xy}</td><td class=mut>{e(note)}</td></tr>")
    H.append("</table></div>")
    r = rules.get("rule_audit")
    if r and r.get("rows"):
        H.append("<h2>KiCad's design rules against the fabricator's (MF19)</h2><p class=mut>A clean KiCad DRC shows that the board meets "
                 "the rules below, in the left column. Where they are looser than the fabricator's, it shows nothing about the fabricator's "
                 "limit; the last column names the check here that measures it.</p><div class=wrap><table><tr><th>KiCad rule</th>"
                 "<th class=n>KiCad (mm)</th><th class=n>fabricator (mm)</th><th>DRC test</th><th>Verdict</th><th>Measured by</th></tr>")
        for x in r["rows"]:
            H.append(f"<tr><td><code>{e(x['rule'])}</code></td><td class=n>{e(str(x['kicad_mm']))}</td><td class=n>{e(' / '.join(f'{v:g}' for v in x['fab'].values()))}</td>"
                     f"<td>{e(str(x['drc_test'] or ''))} {e(str(x['severity'] or ''))}</td><td>{e(x['verdict'])}</td><td>{e(x['measured_by'])}</td></tr>")
        H.append("</table></div>")
    r = rules.get("output_match")
    if r:
        p = r["params"]
        H.append(f"<h2>The files that are sent (MF20)</h2><p class=mut>Files: {e(str(p.get('files')))} ({e(str(p.get('source')))}). "
                 f"Compared at {p.get('resolution_mm')} mm; detects {e(str(p.get('detects')))}. What the fabricator's CAM does with the files "
                 f"is not visible here.</p><div class=wrap><table><tr><th>Layer</th><th>File</th><th class=n>in the file (mm2)</th>"
                 f"<th class=n>checked (mm2)</th><th class=n>differences</th></tr>")
        for k, m in r["metrics"].items():
            if "file_mm2" in m:
                H.append(f"<tr><td>{e(k)}</td><td>{e(m['file'])}</td><td class=n>{m['file_mm2']}</td><td class=n>{m['checked_mm2']}</td>"
                         f"<td class=n>{m['differences']}</td></tr>")
        H.append("</table></div>")
        if "holes" in r["metrics"]:
            H.append(f"<p class=mut>Holes: {e(json.dumps(r['metrics']['holes']))}</p>")
        feat = r["metrics"].get("features") or {}
        if feat:
            H.append("<p>Feature by feature (1 um):</p><div class=wrap><table><tr><th>Layer</th><th>Compared</th></tr>"
                     + "".join(f"<tr><td>{e(k)}</td><td class=mut>{e(json.dumps(v))}</td></tr>" for k, v in feat.items()) + "</table></div>")
        for k in ("outline", "placement", "bom"):
            if k in r["metrics"]:
                H.append(f"<p class=mut>{k.title()}: {e(json.dumps(r['metrics'][k]))}</p>")
    H.append("<h2>Limits used</h2><div class=wrap><table><tr><th>Limit</th><th>Value</th><th>Origin</th><th>Measures</th><th>Certainty</th>"
             "<th>Judged by</th><th>The page's words</th><th>Note</th></tr>")
    for l in fabd["limits"]:
        v = l["value"] if not isinstance(l["value"], dict) else "table"
        src_ = e(l["source"]) if not l.get("url") else f"<a href='{e(l['url'])}'>{e(l['source'])}</a>"
        judged = ", ".join(l.get("judged_by") or []) or ({"order": "order checklist", "input": "input", "record": "recorded"}.get(l.get("use"), "no"))
        H.append(f"<tr><td><code>{e(l['key'])}</code></td><td>{e(str(v))}</td><td>{src_}</td><td>{e((l.get('measurand') or '').replace('_', ' '))}</td>"
                 f"<td>{e(l.get('certainty') or '')}</td><td>{e(judged)}</td><td class=mut>{e(l.get('quote') or '')}</td>"
                 f"<td class=mut>{e((l.get('note') or '') + ('; pages disagree: ' + l['conflict'] if l.get('conflict') else ''))}</td></tr>")
    H.append("</table></div><p class=mut>published: on the fabricator's pages ("
             + ", ".join(f"<a href='{e(u)}'>{e(u)}</a>" for u in fabd["sources"])
             + "); standard: an industry standard; assumed: neither, used by the rankings only; config: set in the board config.</p>")
    for n in fabd["notes"]:
        H.append(f"<p class=mut>{e(n)}</p>")
    r = ranks.get("acid_trap")
    if r:
        p, s = r["params"], r["summary"]
        H.append("<h2>Ranking: acute corners in copper of one net (MF1)</h2>"
                 f"<p class=mut>Corners with an opening angle below {p['listed_below_deg']:g} deg ({p['acute_deg']:g} deg less the measurement "
                 f"margin). This is the drafting rule of thumb about acute angles; no fabricator in the profile publishes a limit for it and "
                 f"nothing here models etching. {s['on_track']} corners are formed by a track and are listed thinnest track first, since a "
                 f"thin track has the least copper to lose; {s['wide_copper']} lie between pads and pour copper only and are not listed.</p>"
                 "<div class=wrap><table><tr><th>Layer</th><th>x, y (mm)</th><th>Net</th><th class=n>angle (deg)</th><th class=n>track (mm)</th><th>BGA field</th></tr>")
        for x in [x for x in r["rows"] if x["track_mm"] is not None][:40]:
            H.append(f"<tr><td>{e(x['layer'])}</td><td>{x['x']:.2f}, {x['y']:.2f}</td><td>{e(x['net'])}</td><td class=n>{x['angle_deg']:.0f}</td>"
                     f"<td class=n>{x['track_mm']:g}</td><td>{'yes' if x['bga'] else ''}</td></tr>")
        H.append("</table></div>")
    r = ranks.get("copper")
    if r:
        H.append("<h2>Ranking: copper share by layer (MF21)</h2><p class=mut>No limit is published for copper balance; an uneven pair of "
                 "layers about the centre of the stack is a known cause of bow and twist.</p><div class=wrap><table><tr><th>Layer</th>"
                 "<th class=n>copper (%)</th><th class=n>thickness (um)</th><th>mirror layer</th><th class=n>difference (%)</th></tr>")
        for x in r["rows"]:
            H.append(f"<tr><td>{e(x['layer'])}</td><td class=n>{x['copper_pct']}</td><td class=n>{x['thickness_um']}</td><td>{e(x['mirror'])}</td>"
                     f"<td class=n>{x['difference_pct']}</td></tr>")
        H.append("</table></div>")
    r = ranks.get("reflow")
    if r and r["rows"]:
        p, s = r["params"], r["summary"]
        H.append("<h2>Ranking: reflow, how far each part's peak lies from the board median (MF7, model)</h2>"
                 f"<p class=mut>Generic oven profile (air by zone: {e(', '.join(f'{z:g}' for z in p['zones_c']))} C for {p['zone_s']:g} s each), "
                 f"spike zones shifted {p['shift_c']:+g} C so that the median peak is {s['median_peak_c']:g} C ({e(p['target_source'])}). The "
                 f"absolute peak is therefore an input; only the differences are a result. Film coefficient {p['h_top']:g} W/m2K "
                 f"({e(p['h_source'])}), belt {p['belt_mm_s']:g} mm/s along {e(p['travel'])}. The range is over {len(p['variants'])} model "
                 f"variants ({e(', '.join(p['variants'][1:]))}). {s['assumed_bodies']} of {s['parts']} part bodies are an assumed block. "
                 f"Spread over the parts: {s['spread_c']:g} C (range {s['spread_range_c'][0]:g} to {s['spread_range_c'][1]:g} C). "
                 f"A part is as warm as the board under it: a temperature difference inside a part is not modelled.</p>")
        if r.get("image"):
            H.append(f"<p><img src='{e(r['image'])}' alt='modelled peak temperature'></p>")
        rows = r["rows"]
        show = rows[:12] + ([None] if len(rows) > 18 else []) + (rows[-6:] if len(rows) > 18 else rows[12:])
        H.append("<div class=wrap><table><tr><th>Part</th><th>Side</th><th class=n>from the median (C)</th><th>sensitivity to parameters (C)</th>"
                 "<th class=n>lag (s)</th><th class=n>under the part (C)</th><th>Body</th></tr>")
        for x in show:
            if x is None:
                H.append(f"<tr><td colspan=7 class=mut>... {len(rows) - 18} parts between ...</td></tr>")
                continue
            H.append(f"<tr><td>{e(x['ref'])}</td><td>{x['side']}</td><td class=n>{x['d_peak_c']:+.1f}</td>"
                     f"<td>{x['d_range_c'][0]:+.1f} to {x['d_range_c'][1]:+.1f}</td>"
                     f"<td class=n>{x['lag_s']:.0f}</td><td class=n>{num(x.get('within_c'), 1)}</td><td class=mut>{e(x['body'])}</td></tr>")
        H.append("</table></div><p class=mut>The 12 coldest and 6 hottest parts; all parts are in dfm.json.</p>")
    r = ranks.get("chips")
    if r and r["rows"]:
        p = r["params"]
        H.append("<h2>Ranking: chips, time between the two joints melting (MF8, model)</h2>"
                 "<p class=mut>Timing is one input to tombstoning; the wetting forces (pad geometry, paste volume and offset, placement, "
                 "alloy) are not modelled, so the land geometry is listed beside it. <i>dt</i> is the difference of the two joints' lag "
                 "times; <i>range</i> is over the assumed paste conductance (x 0.5, x 2), the film coefficient (x 0.6, x 1.6) and a window "
                 "twice as large; <i>in this oven</i> is the same window marched through the oven profile of the board run, which also "
                 "carries one pad into each zone before the other (<i>belt</i>: the first pad's lead). Each figure's sensitivity to the "
                 "model's parameters (paste, window, film coefficient; not an error bar) and the oven "
                 f"figure are given for the first {p['detail']} layouts. Parts with the same layout and result share a row. The window's "
                 f"edge follows the {e(p['edge'])}.</p><div class=wrap><table><tr><th>Parts</th><th>Package</th><th class=n>dt (s)</th>"
                 "<th>dt sensitivity (s)</th><th class=n>in this oven (s)</th><th>oven sensitivity (s)</th><th class=n>belt (s)</th><th>melts first</th><th>melts last</th><th>copper leaving the pads (mm)</th>"
                 "<th class=n>toe (mm)</th><th class=n>pad areas</th></tr>")
        for x in r["rows"][:30]:
            refs = ", ".join(x["refs"][:5]) + (f" +{x['n'] - 5}" if x["n"] > 5 else "")
            ent = "; ".join(f"{q['pad']}: {q['entry_mm']:.2f}" + (" in a pour" if q["share"] >= 0.6 else "")
                            + (" +via in pad" if q["vias_in_pad"] else "") + (f" +{q['plane_vias']} plane via" if q["plane_vias"] else "")
                            for q in x["pads"])
            rng = f"{x['dt_range_s'][0]:g} to {x['dt_range_s'][1]:g}" if x.get("dt_range_s") else ""
            H.append(f"<tr><td>{e(refs)}</td><td>{e(x['package'])}</td><td class=n>{x['dt_s']:g}</td><td>{rng}</td>"
                     f"<td class=n>{'' if x.get('dt_oven_s') is None else format(x['dt_oven_s'], 'g')}</td>"
                     f"<td>{'' if not x.get('dt_oven_range_s') else format(x['dt_oven_range_s'][0], 'g') + ' to ' + format(x['dt_oven_range_s'][1], 'g')}</td>"
                     f"<td class=n>{x.get('belt_s', 0):+.2f}</td>"
                     f"<td>{e(str(x['first_pad']))} {e(x['first_net'])}</td><td>{e(str(x['last_pad']))} {e(x['last_net'])}</td>"
                     f"<td class=mut>{e(ent)}</td><td class=n>{x['land']['toe_mm']:.2f}</td><td class=n>{x['land']['area_ratio']:g}</td></tr>")
        H.append(f"</table></div><p class=mut>The {min(30, len(r['rows']))} layouts with the largest difference of {r['summary']['layouts']} "
                 f"({r['summary']['chips']} chips); all are in dfm.json.</p>")
    if ev.get("observed"):
        o = ev["observed"]
        H.append(f"<h2>Observed defects</h2><p class=mut>{o['defects']} defective parts on {o.get('boards_built')} boards built, from "
                 f"{e(o['file'])}{' &mdash; ' + e(o['error']) if o.get('error') else ''}. For each ranking, the share of (defective, sound) "
                 f"pairs it puts in the right order: 0.5 is chance, 1 is perfect; with its bootstrap 95 % interval. Below {MIN_DEFECTIVE} "
                 "defective parts no score is given.</p><div class=wrap><table><tr><th>Defect</th><th>Ranking</th><th class=n>order</th>"
                 "<th>95 % interval</th><th class=n>defective</th><th class=n>sound</th></tr>")
        for sc in o["scores"]:
            H.append(f"<tr><td>{e(sc['kind'])}</td><td>{e(sc['ranking'])}</td><td class=n>{'' if sc['order'] is None else format(sc['order'], '.2f')}</td>"
                     f"<td>{e(sc['note']) if sc['order'] is None else format(sc['ci95'][0], '.2f') + ' to ' + format(sc['ci95'][1], '.2f')}</td>"
                     f"<td class=n>{sc['defective_ranked']}</td><td class=n>{sc['sound']}</td></tr>")
        H.append("</table></div>")
    if ev["images"]:
        H.append("<h2>Findings on the copper</h2>")
        for lname, img in ev["images"].items():
            H.append(f"<h3>{e(lname)}</h3><p><img src='{e(img)}' alt='{e(lname)}'></p>")
    if ev["waivers"]:
        H.append("<h2>Waivers</h2><table><tr><th>Waiver</th><th>Status</th><th>Reason</th><th>Reviewer</th></tr>")
        for w in ev["waivers"]:
            H.append(f"<tr><td><code>{e(str(w.get('id') or w.get('check')))}</code></td><td>{e(w['status'])}</td>"
                     f"<td>{e(str(w.get('reason', '')))}</td><td>{e(str(w.get('reviewer', '')))}</td></tr>")
        H.append("</table>")
    H.append("<p class=mut>Information findings and every row of every ranking are in dfm.json.</p></body></html>")
    with open(path, "w") as fh:
        fh.write("\n".join(H))
