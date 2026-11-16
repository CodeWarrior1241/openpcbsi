"""Manufacturing sign-off (REQUIREMENTS 5.17): `openpcbsi dfm BOARD`.

    python3 -m openpcbsi dfm BOARD [-c CONFIG] [--fab jlcpcb | FILE] [--checks track,etch,mask,balance,reflow,tombstone]
                                   [--waivers FILE] [--json OUT] [-o OUT_DIR] [--no-markers]

Checks the layout against what a fabricator and an assembler can make, before the files are sent:

  MF0  fabricator profile (fab.py): the limits, each marked published / assumed / board config
  MF1  acid_trap      acute corners in copper of one net (etch.py)
  MF2  same_net_gap   gaps narrower than the fabricator's same-net spacing (etch.py)
  MF3  fab_track      track width against the fabricator's minimum for the layer count and copper weight
  MF4  mask_web       solder-mask web between pads, pad spacing that leaves no room for a web
  MF5  mask_exposed   copper of another net inside or next to a mask opening
  MF6  pad_balance    copper leaving the two pads of a chip (tombstoning, geometric)
  MF7  reflow         peak temperature, time above liquidus and lag of every part (reflow.py)
  MF8  tombstone      time between the two joints of a chip reaching liquidus (reflow.py)

Track-to-track and pad-to-track clearance between different nets, annular rings, copper slivers and
starved thermals stay with KiCad's DRC (EV1, part of `evaluate`); they are not repeated here.

Writes <run>/dfm/dfm.json (schema `openpcbsi-dfm/1`), dfm.html, one image per copper layer, the peak
temperature map and a marker board. Findings use the G19 model, waivers the G20 format, and the exit
code is that of `evaluate` (0 / 1 warnings / 2 fails / 3 could not run).
"""
import datetime
import html
import json
import os

import numpy as np
from scipy import ndimage

from . import fab as fab_mod
from . import geom, packages
from .drc import finding, result

SCHEMA = "openpcbsi-dfm/1"
GROUPS = {"track": ["fab_track"], "etch": ["acid_trap", "same_net_gap"], "mask": ["mask_web", "mask_exposed"],
          "balance": ["pad_balance"], "reflow": ["reflow"], "tombstone": ["tombstone"]}
REQ = {"acid_trap": "MF1", "same_net_gap": "MF2", "fab_track": "MF3", "mask_web": "MF4", "mask_exposed": "MF5",
       "pad_balance": "MF6", "reflow": "MF7", "tombstone": "MF8"}
TOL = 0.003          # mm: polygons approximate arcs to 2 um


# ------------------------------------------------------------------------------------------------
# MF3

def rule_fab_track(ctx, fabr):
    b = ctx.b
    lim, bga = fabr.get("track_width_mm"), fabr.get("bga_track_mm")
    T = b.tracks
    fs, metrics = [], {}
    if lim is None or not len(T["w"]):
        return result("fab_track", "MF3", "Track width against the fabricator's minimum", {"limit_mm": lim}, {}, [])
    g = b.grid()
    mx, my = (T["x1"] + T["x2"]) / 2, (T["y1"] + T["y2"]) / 2
    in_area = g.sample(ctx.areas(), mx, my)
    for li, lname in enumerate(b.layers):
        on = T["layer"] == li
        if not on.any():
            continue
        metrics[lname] = {"min_width_mm": round(float(T["w"][on].min()), 4)}
        for w in np.unique(np.round(T["w"][on], 4)):
            if w >= lim - 1e-6:
                continue
            sel = on & (np.abs(T["w"] - w) < 5e-5)
            for field, local in ((False, lim), (True, bga if bga is not None else lim)):
                s = sel & ((in_area != 0) == field)
                if not s.any() or w >= local - 1e-6:
                    continue
                k = int(np.flatnonzero(s)[np.argmax(T["len"][s])])
                where = " inside a BGA or connector field" if field else ""
                fs.append(finding("fab_track", "error", f"{int(s.sum())} track segments of {w:g} mm ({float(T['len'][s].sum()):.1f} mm in all)"
                                  f"{where}, below the {local:g} mm minimum ({fabr.src('bga_track_mm' if field and bga is not None else 'track_width_mm')}: "
                                  f"{fabr.limits['track_width_mm']['note']})", b.net_names[T["net"][k]], lname, mx[k], my[k],
                                  value=float(w), limit=local, unit="mm", kind=f"{w:g}{'f' if field else ''}"))
    return result("fab_track", "MF3", "Track width against the fabricator's minimum",
                  {"limit_mm": lim, "bga_limit_mm": bga, "source": fabr.src("track_width_mm")}, metrics, fs)


# ------------------------------------------------------------------------------------------------
# MF4 / MF5

def _rings(polys):
    return [ring for poly in polys for ring in poly]


def _bbox(rings):
    a = np.vstack([np.asarray(r, dtype=float) for r in rings])
    return a[:, 0].min(), a[:, 1].min(), a[:, 0].max(), a[:, 1].max()


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


def rule_mask_web(ctx, fabr):
    b = ctx.b
    web, spacing, pp = fabr.get("mask_web_mm"), fabr.get("mask_pad_spacing_mm"), fabr.get("smd_pad_to_pad_mm")
    reach = max(v for v in (web, spacing, pp, 0.05) if v is not None) + 0.02
    worst, metrics = {}, {}
    for mname, li in (("F.Mask", 0), ("B.Mask", len(b.layers) - 1)):
        P = sorted(_mask_pads(b, mname), key=lambda q: q["bb"][0])
        pairs = 0
        min_web = None
        for a_, A in enumerate(P):
            for B in P[a_ + 1:]:
                if B["bb"][0] > A["bb"][2] + reach:
                    break
                if B["bb"][1] > A["bb"][3] + reach or B["bb"][3] < A["bb"][1] - reach:
                    continue
                pa, pb = A["pad"], B["pad"]
                if pa["net"] == pb["net"] and pa["net"]:
                    continue
                if pa["ref"] == pb["ref"] and pa["num"] == pb["num"]:
                    continue
                w = geom.rings_distance(A["mask"], B["mask"])
                if w > reach:
                    continue
                gap = geom.rings_distance(A["cu"], B["cu"])
                pairs += 1
                min_web = w if min_web is None else min(min_web, w)
                x, y = (pa["x"] + pb["x"]) / 2, (pa["y"] + pb["y"]) / 2
                name = f"{pa['ref']}.{pa['num']} - {pb['ref']}.{pb['num']}"
                found = []
                if pp is not None and pa["kind"] == pb["kind"] == "smd" and gap < pp - TOL:
                    found.append(("pad_gap", "error", gap, pp, f"pads {gap:.3f} mm apart, below the {pp:g} mm minimum between "
                                  f"SMD pads of different nets ({fabr.src('smd_pad_to_pad_mm')})"))
                if w <= TOL:
                    if spacing is not None and gap < spacing - TOL:
                        found.append(("no_dam", "warning", gap, spacing, f"no solder-mask web: pads {gap:.3f} mm apart, less than the "
                                      f"{spacing:g} mm that keeps a web ({fabr.src('mask_pad_spacing_mm')}); solder can bridge"))
                    else:
                        found.append(("no_dam", "warning", w, web, f"the mask openings merge although the pads are {gap:.3f} mm apart: "
                                      f"reduce the mask expansion ({max(A['exp'], B['exp']) * 1000:.0f} um) to leave a web"))
                elif web is not None and w < web - TOL:
                    found.append(("web", "error", w, web, f"solder-mask web of {w:.3f} mm, below the {web:g} mm minimum "
                                  f"({fabr.src('mask_web_mm')}): it will be removed or lift"))
                for kind, sev, val, limit, msg in found:
                    key = (mname, tuple(sorted((pa["ref"], pb["ref"]))), kind)
                    e = worst.setdefault(key, {"n": 0, "val": 1e9})
                    e["n"] += 1
                    if val < e["val"]:
                        e.update(val=val, sev=sev, limit=limit, msg=msg, name=name, x=x, y=y, layer=b.layers[li],
                                 nets=(pa["net"], pb["net"]), part=pa["ref"])
        metrics[mname] = {"openings": len(P), "close_pairs": pairs, "min_web_mm": None if min_web is None else round(max(min_web, 0.0), 4)}
    fs = []
    for (mname, refs, kind), e in sorted(worst.items(), key=lambda kv: (kv[0][2], kv[1]["val"])):
        more = f" ({e['n']} pad pairs of {' / '.join(dict.fromkeys(refs))}; worst shown)" if e["n"] > 1 else ""
        fs.append(finding("mask_web", e["sev"], f"{e['name']}: {e['msg']}{more}", e["nets"][0], e["layer"], e["x"], e["y"],
                          value=round(e["val"], 4), limit=e["limit"], unit="mm", aggressor=e["nets"][1], part=e["part"],
                          kind=f"{kind}:{'/'.join(refs)}", pairs=e["n"]))
    return result("mask_web", "MF4", "Solder-mask web and pad spacing", {"web_mm": web, "pad_spacing_mm": spacing,
                  "smd_pad_to_pad_mm": pp}, metrics, fs)


def rule_mask_exposed(ctx, fabr, res):
    """Per pad with a mask opening: copper of another net inside the opening (fail) or closer to it than
    the fabricator's opening-to-track clearance (warn). Also the copper of the pad's own net that the
    opening uncovers beyond the pad (`exposure`, used by MF6)."""
    b = ctx.b
    g = b.grid()
    lim = fabr.get("mask_to_track_mm")
    reach = (lim or 0.1) + 0.1
    worst, exposure, checked = {}, {}, 0
    for mname, li in (("F.Mask", 0), ("B.Mask", len(b.layers) - 1)):
        r = b.raster(li)
        P = _mask_pads(b, mname)
        bbs = np.array([q["bb"] for q in P]) if P else np.zeros((0, 4))
        for kq, q in enumerate(P):
            p = q["pad"]
            nid = p["net_id"]
            x0, y0, x1, y1 = q["bb"]
            i0, j0 = g.ij(x0 - reach, y0 - reach)
            i1, j1 = g.ij(x1 + reach, y1 + reach)
            sub = r[max(int(i0), 0):int(i1) + 1, max(int(j0), 0):int(j1) + 1]
            foreign_near = ((sub != 0) & (sub != nid)).any()
            if not foreign_near and q["exp"] <= 0:
                continue
            checked += 1
            gw, lab = b.raster_window(li, x0 - reach, y0 - reach, x1 + reach, y1 + reach, res)
            O = geom.raster_rings(gw, q["mask"])
            Pm = geom.raster_rings(gw, q["cu"])
            own = (lab == nid) if nid else Pm
            extra = int((O & own & ~Pm).sum())
            if Pm.any():
                exposure[(p["ref"], p["num"])] = 1.0 + extra / float(Pm.sum())
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
                found = ("exposed", "error", -area, 0.0, net,
                         f"{area * 1e6:.0f} um2 of {net} lies inside the mask opening of the pad: solder can bridge to it")
            elif lim is not None:
                d = ndimage.distance_transform_edt(~O)
                dmin = (float(d[foreign].min()) - 1.0) * res
                if dmin < lim - 1.5 * res:
                    k = np.flatnonzero(foreign.ravel())[int(np.argmin(d[foreign]))]
                    net = b.net_names[int(lab.ravel()[k])]
                    found = ("near", "warning", dmin, lim, net, f"{net} is {max(dmin, 0) * 1000:.0f} um from the mask opening, closer than "
                             f"{lim * 1000:.0f} um ({fabr.src('mask_to_track_mm')}): mask misregistration can uncover it")
            if found:
                kind, sev, val, limit, net, msg = found
                key = (p["ref"], kind)
                e = worst.setdefault(key, {"n": 0, "val": 1e9})
                e["n"] += 1
                if val < e["val"]:
                    e.update(val=val, sev=sev, limit=limit, msg=msg, pad=p, layer=b.layers[li], other=net)
    fs = []
    for (ref, kind), e in sorted(worst.items(), key=lambda kv: (kv[0][1], kv[1]["val"])):
        p = e["pad"]
        more = f" ({e['n']} pads of {ref}; worst shown)" if e["n"] > 1 else ""
        val = -e["val"] if kind == "exposed" else e["val"]
        fs.append(finding("mask_exposed", e["sev"], f"{ref}.{p['num']}: {e['msg']}{more}", p["net"], e["layer"], p["x"], p["y"],
                          value=round(val, 5), limit=e["limit"], unit="mm2" if kind == "exposed" else "mm", aggressor=e["other"],
                          part=ref, kind=kind, pads=e["n"]))
    metrics = {"pads": {"checked": checked, "with_exposure": sum(1 for v in exposure.values() if v > 1.02)}}
    return result("mask_exposed", "MF5", "Copper of another net inside or next to a mask opening",
                  {"opening_to_track_mm": lim, "res_mm": res}, metrics, fs), exposure


# ------------------------------------------------------------------------------------------------
# MF6

def _pad_entry(b, pad, ring_px, g):
    """Copper of the pad's net crossing a ring just outside the pad: (width mm, share of the ring)."""
    rings = _rings(pad["poly"])
    x0, y0, x1, y1 = _bbox(rings)
    m = (ring_px + 3) * g.res
    i0, j0 = g.ij(x0 - m, y0 - m)
    i1, j1 = g.ij(x1 + m, y1 + m)
    i0, j0 = max(int(i0), 0), max(int(j0), 0)
    i1, j1 = min(int(i1) + 1, g.ny), min(int(j1) + 1, g.nx)
    sub = geom.Grid.__new__(geom.Grid)
    sub.res, sub.x0, sub.y0, sub.nx, sub.ny = g.res, g.x0 + j0 * g.res, g.y0 + i0 * g.res, j1 - j0, i1 - i0
    M = geom.raster_rings(sub, rings)
    d = ndimage.distance_transform_edt(~M)
    ring = (d > ring_px - 0.5) & (d <= ring_px + 0.5)
    lab = b.raster(pad["layer_ids"][0])[i0:i1, j0:j1]
    n = int(ring.sum())
    cu = int((ring & (lab == pad["net_id"])).sum()) if pad["net_id"] else 0
    return cu * g.res, (cu / n if n else 0.0)


def rule_pad_balance(ctx, exposure):
    from . import reflow
    b, cfg = ctx.b, ctx.cfg
    p = cfg["dfm"]["pad_balance"]
    g = b.grid()
    ring_px = max(2, int(round(p["ring_mm"] / g.res)))
    small = packages.size_rank(p["small_package"])
    V = b.vias
    fs, metrics, rows = [], {}, []
    reach = p["via_reach_mm"]
    plane = {}                                  # via index -> tied to a zone on some layer

    def to_plane(k):
        if k not in plane:
            plane[k] = bool(b.zone_connected_layers(k))
        return plane[k]
    for fp, chip, pads in reflow.chips_in_scope(ctx):
        info = []
        for pad in pads:
            w, share = _pad_entry(b, pad, ring_px, g)
            rings = _rings(pad["poly"])
            x0, y0, x1, y1 = _bbox(rings)
            near = np.flatnonzero((V["net"] == pad["net_id"]) & (V["x"] >= x0) & (V["x"] <= x1) & (V["y"] >= y0) & (V["y"] <= y1))
            vias = int(geom.points_in_rings(V["x"][near], V["y"][near], rings).sum()) if len(near) else 0
            close = np.flatnonzero((V["net"] == pad["net_id"]) & (V["x"] >= x0 - reach) & (V["x"] <= x1 + reach)
                                   & (V["y"] >= y0 - reach) & (V["y"] <= y1 + reach)) if pad["net_id"] else []
            info.append({"pad": pad["num"], "net": pad["net"], "entry_mm": round(w, 3), "share": round(share, 3),
                         "vias_in_pad": vias, "plane_vias": sum(1 for k in close if to_plane(int(k))),
                         "exposure": round(exposure.get((pad["ref"], pad["num"]), 1.0), 3)})
        a, c = info
        hi, lo = (a, c) if a["entry_mm"] >= c["entry_mm"] else (c, a)
        ratio = hi["entry_mm"] / max(lo["entry_mm"], p["min_entry_mm"])
        why = []
        if ratio >= p["ratio_warn"]:
            why.append(f"copper leaving pad {hi['pad']} ({hi['net']}) is {hi['entry_mm']:.2f} mm wide against {lo['entry_mm']:.2f} mm "
                       f"at pad {lo['pad']} ({lo['net']}): {ratio:.1f} x (limit {p['ratio_warn']:g} x, proposed)")
        solid = [x["share"] >= p["solid_fraction"] for x in info]
        if solid[0] != solid[1]:
            s, o = (a, c) if solid[0] else (c, a)
            why.append(f"pad {s['pad']} ({s['net']}) sits in a pour ({s['share']:.0%} of its edge in copper), pad {o['pad']} does not "
                       f"({o['share']:.0%}): give the poured pad a thermal relief")
        if bool(a["vias_in_pad"]) != bool(c["vias_in_pad"]):
            s = a if a["vias_in_pad"] else c
            why.append(f"via in pad {s['pad']} ({s['net']}) only")
        note = []                               # information only: MF8 puts a number on it
        if bool(a["plane_vias"]) != bool(c["plane_vias"]):
            s, o = (a, c) if a["plane_vias"] else (c, a)
            note.append(f"pad {s['pad']} ({s['net']}) has {s['plane_vias']} via(s) into a plane within {reach:g} mm, pad {o['pad']} "
                        f"({o['net']}) has none (see tombstone, MF8, for the estimated effect)")
        if abs(a["exposure"] - c["exposure"]) >= 0.25:
            s = a if a["exposure"] > c["exposure"] else c
            why.append(f"the mask opening of pad {s['pad']} uncovers {100 * (s['exposure'] - 1):.0f} % more copper than the pad")
        metrics[fp["ref"]] = {"ratio": round(ratio, 2), "entry_mm": [a["entry_mm"], c["entry_mm"]]}
        rows.append({"ref": fp["ref"], "package": chip["code"], "ratio": round(ratio, 2), "pads": info, "reasons": len(why)})
        if why or note:
            sev = "warning" if why and packages.size_rank(chip["code"]) <= small else "info"
            fs.append(finding("pad_balance", sev, f"{fp['ref']} ({chip['code']}): " + "; ".join(why + note), hi["net"],
                              b.layers[pads[0]["layer_ids"][0]], fp["x"], fp["y"], value=round(ratio, 2), limit=p["ratio_warn"],
                              unit="x", part=fp["ref"], aggressor=lo["net"]))
    rows.sort(key=lambda r: -r["ratio"])
    return dict(result("pad_balance", "MF6", "Copper leaving the two pads of a chip (tombstoning)", p, metrics, fs), parts=rows)


# ------------------------------------------------------------------------------------------------
def _ensure_validation(cache, log):
    from . import etch, reflow, validate
    for solver, mod in (("etch", etch), ("reflow", reflow)):
        if validate.status(cache, solver) != "validated":
            validate.register(cache, solver, mod.validation_cases())
            log(f"  G17 validation {solver}: {validate.status(cache, solver)}")


def run(board, cfg_path=None, fab=None, checks=None, out_dir=None, waivers=None, json_out=None, markers=True, log=print):
    from . import etch, evaluate, reflow, report as report_mod, run as runmod, validate
    import time
    t0 = time.time()
    b, cfg, ctx = runmod.load(board, cfg_path)
    if fab:
        cfg.data["dfm"]["fab"] = fab
    fabr = fab_mod.resolve(fab_mod.load(cfg["dfm"]["fab"]), b, cfg)
    want = list(GROUPS) if not checks else [c.strip().lower() for c in checks]
    unknown = [c for c in want if c not in GROUPS]
    if unknown:
        raise ValueError(f"unknown dfm checks {unknown}; choose from {', '.join(GROUPS)}")
    od = out_dir or os.path.join(runmod.RUNS, f"{b.name}_{b.sha[:8]}", "dfm")
    os.makedirs(od, exist_ok=True)
    log(f"openpcbsi dfm: {b.path}")
    log(f"  sha256 {b.sha[:16]}  {len(b.layers)} layers  fabricator profile {fabr.name} ({os.path.basename(fabr.path)}, read {fabr.retrieved})")
    for n in fabr.notes:
        log(f"  note: {n}")
    _ensure_validation(runmod.CACHE, log)
    res_hi = b.grid().res / max(1, int(cfg["dfm"]["acid"]["hires_div"]))
    rules, val, assumed, skipped = [], {}, {}, {}
    exposure = {}
    geo = lambda r: "geometry" if (r.get("metrics") or r.get("findings")) else "nothing to check"
    if "track" in want:
        log("MF3 track width:")
        rules.append(rule_fab_track(ctx, fabr))
        val["fab_track"] = geo(rules[-1])
        log(f"  {len(rules[-1]['findings'])} widths below the minimum of {fabr.get('track_width_mm')} mm")
    if "etch" in want:
        log(f"MF1/MF2 acid traps and same-net gaps (ball {cfg['dfm']['acid']['ball_mm'] or (fabr.get('same_net_spacing_mm') or 0.25) / 2:g} mm, "
            f"measured at {res_hi * 1000:g} um):")
        for r in etch.run(ctx, fabr, log=log):
            rules.append(r)
            val[r["rule"]] = validate.status(runmod.CACHE, "etch")
        assumed["acid_trap"] = (f"contact angle {cfg['dfm']['acid']['contact_angle_deg']:g} deg assumed "
                                f"(limit {rules[-2]['params']['limit_deg']:g} deg)")
    if "mask" in want:
        log("MF4/MF5 solder mask:")
        rules.append(rule_mask_web(ctx, fabr))
        val["mask_web"] = geo(rules[-1])
        r, exposure = rule_mask_exposed(ctx, fabr, res_hi)
        rules.append(r)
        val["mask_exposed"] = "geometry" if r["metrics"]["pads"]["checked"] or any(m["openings"] for m in rules[-2]["metrics"].values()) else "nothing to check"
        log(f"  {sum(m['openings'] for m in rules[-2]['metrics'].values())} openings, {len(rules[-2]['findings'])} web findings, "
            f"{len(r['findings'])} exposure findings")
    order = None
    if "balance" in want:
        log("MF6 pad balance of chips:")
        r = rule_pad_balance(ctx, exposure)
        rules.append(r)
        val["pad_balance"] = "geometry" if r["parts"] else "nothing to check"
        order = [x["ref"] for x in r["parts"]]
        log(f"  {len(r['parts'])} chips, {len(r['findings'])} unbalanced")
    if "reflow" in want:
        log("MF7 reflow, board level:")
        r = reflow.board(ctx, fabr, od, log=log)
        rules.append(r)
        val["reflow"] = validate.status(runmod.CACHE, "reflow") if r["parts"] else "nothing to check"
        assumed["reflow"] = (f"generic oven profile, h {r['params']['h_top']:g} W/m2K ({r['params']['h_source']}), "
                             f"{r['summary']['assumed_bodies']} of {r['summary']['parts']} part bodies assumed")
    if "tombstone" in want:
        log("MF8 tombstoning, pad level:")
        r = reflow.tombstone(ctx, fabr, order=order, log=log)
        rules.append(r)
        val["tombstone"] = validate.status(runmod.CACHE, "reflow") if r["parts"] else "nothing to check"
        assumed["tombstone"] = "chip bodies from the package table, paste conductance and film coefficient assumed"
    for grp, names in GROUPS.items():
        if grp not in want:
            for n in names:
                skipped[n] = "not selected (--checks)"
    findings = []
    for r in rules:
        findings += evaluate.from_rule_findings(r["rule"], "MF", REQ[r["rule"]], r.get("findings", []))
    seen = {}
    for f in findings:
        n = seen.get(f["id"], 0)
        seen[f["id"]] = n + 1
        if n:
            f["id"] = f"{f['id']}.{n}"
    ws, whash = evaluate.load_waivers(waivers)
    wstat = evaluate.apply_waivers(findings, ws)
    per, _ = evaluate.scorecard(findings, val, skipped, assumed)
    per = sorted((c for c in per if c["family"] == "MF"), key=lambda c: c["requirement"])
    code = evaluate.exit_code(findings, cfg)
    rep = {"checks": {"drc": rules}}
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
    counts = {s: sum(1 for f in findings if f["severity"] == s and not f["waived"] and not f["area"]) for s in ("fail", "warn", "info")}
    counts |= {"in_area": sum(1 for f in findings if f["area"] and not f["waived"]), "waived": sum(1 for f in findings if f["waived"])}
    status = "fail" if counts["fail"] else "warn" if counts["warn"] else "pass"
    ev = {"schema": SCHEMA, "tool": "openpcbsi", "version": __import__("openpcbsi").__version__,
          "date": datetime.datetime.now().isoformat(timespec="seconds"), "board": b.path, "board_sha256": b.sha,
          "kicad": b.snap.get("kicad"), "config": cfg.path or "detected from the board", "config_hash": cfg.hash(),
          "fab": fabr.describe(), "waivers_file": os.path.abspath(waivers) if waivers else None, "waivers_hash": whash,
          "solver_versions": {k: validate.SOLVER_VERSIONS[k] for k in ("etch", "reflow", "thermal")},
          "exit_code": code, "status": status, "checks": per, "waivers": wstat, "counts": counts, "findings": findings,
          "rules": rules, "images": images, "marker_board": marker, "seconds": round(time.time() - t0, 1)}
    path = os.path.join(od, "dfm.json")
    for pth in [path] + ([json_out] if json_out else []):
        with open(pth, "w") as fh:
            json.dump(ev, fh, indent=1, default=runmod._jd)
    write_html(ev, os.path.join(od, "dfm.html"))
    log("")
    log(f"manufacturing scorecard ({fabr.name} profile):")
    for c in per:
        log(f"  {c['status']:18s} {c['requirement']:4s} {c['check']:14s} fail {c['fail']:4d}  warn {c['warn']:4d}  info {c['info']:4d}"
            + (f"  in areas {c['in_area']}" if c["in_area"] else ""))
    log(f"  {counts['fail']} fail, {counts['warn']} warn, {counts['info']} info; {counts['in_area']} in exclusion areas, {counts['waived']} waived")
    log(f"dfm:  {path}")
    log(f"html: {os.path.join(od, 'dfm.html')}")
    log(f"exit code {code}  ({ev['seconds']} s)")
    return code, ev


# ------------------------------------------------------------------------------------------------
def write_html(ev, path):
    from .evaluate import CSS
    e = html.escape
    st = lambda s: f'<span class="{s.replace(" ", "")}">{e(s)}</span>'
    num = lambda v, d=2: "" if v is None else f"{v:.{d}f}"
    fabd = ev["fab"]
    rules = {r["rule"]: r for r in ev["rules"]}
    H = [f"<!doctype html><html><head><meta charset=utf-8><meta name=viewport content='width=device-width,initial-scale=1'>"
         f"<title>openpcbsi dfm</title><style>{CSS} img{{max-width:100%;height:auto}} .wrap{{overflow-x:auto}}</style></head><body>"]
    H.append(f"<h1>Manufacturing sign-off: {e(os.path.basename(ev['board']))} &mdash; {st(ev['status'])}</h1>")
    H.append(f"<p class=mut>{e(ev['board'])}<br>sha256 {e(ev['board_sha256'][:16])} · config {e(str(ev['config']))} · "
             f"fabricator profile {e(fabd['name'])} (read {e(str(fabd['retrieved']))}, hash {e(str(fabd['profile_hash']))}) · "
             f"{e(ev['date'])} · exit code {ev['exit_code']}</p>")
    H.append("<h2>Scorecard</h2><div class=wrap><table><tr><th>Req.</th><th>Check</th><th>Status</th><th class=n>fail</th><th class=n>warn</th>"
             "<th class=n>info</th><th class=n>in areas</th><th class=n>waived</th><th>Validation / inputs</th></tr>")
    for c in ev["checks"]:
        H.append(f"<tr><td>{e(c['requirement'])}</td><td>{e(c['check'])}</td><td>{st(c['status'])}</td><td class=n>{c['fail']}</td>"
                 f"<td class=n>{c['warn']}</td><td class=n>{c['info']}</td><td class=n>{c['in_area']}</td><td class=n>{c['waived']}</td>"
                 f"<td class=mut>{e(str(c.get('why') or c.get('validation') or ''))}</td></tr>")
    H.append("</table></div>")
    H.append("<h2>Fabricator limits used</h2><div class=wrap><table><tr><th>Limit</th><th>Value</th><th>Origin</th><th>Applies to</th></tr>")
    for l in fabd["limits"]:
        H.append(f"<tr><td><code>{e(l['key'])}</code></td><td>{e(str(l['value']))}</td>"
                 f"<td>{e(l['source'])}</td><td class=mut>{e(l.get('note') or '')}</td></tr>")
    H.append("</table></div><p class=mut>published: on the fabricator's capability pages ("
             + ", ".join(f"<a href='{e(u)}'>{e(u)}</a>" for u in fabd["sources"])
             + "); assumed: a typical figure the fabricator does not state; config: set in the board config.</p>")
    for n in fabd["notes"]:
        H.append(f"<p class=warn>{e(n)}</p>")
    rank = {"fail": 0, "warn": 1}
    live = sorted((f for f in ev["findings"] if f["severity"] in rank),
                  key=lambda f: (bool(f["waived"]), bool(f["area"]), rank[f["severity"]], f["requirement"], f["check"]))
    H.append(f"<h2>Findings: fail and warn ({len(live)})</h2><div class=wrap><table><tr><th>ID</th><th>Sev.</th><th>Check</th><th>Reason</th>"
             "<th>Nets / refs</th><th>Layer</th><th>x, y (mm)</th><th>Note</th></tr>")
    for f in live[:2000]:
        note = "waived: " + f["waived"]["reason"] if f["waived"] else ("in area " + f["area"] if f["area"] else "")
        xy = f"{f['x']:.2f}, {f['y']:.2f}" if f["x"] is not None and f["y"] is not None else ""
        H.append(f"<tr><td><code>{e(f['id'])}</code></td><td>{st(f['severity'])}</td><td>{e(f['check'])} <span class=mut>{e(f['requirement'])}</span></td>"
                 f"<td>{e(f['reason'])}</td><td>{e(', '.join(f['nets'] + f['refs']))}</td><td>{e(f['layer'] or '')}</td>"
                 f"<td>{xy}</td><td class=mut>{e(note)}</td></tr>")
    H.append("</table></div>")
    r = rules.get("acid_trap")
    if r:
        p = r["params"]
        H.append("<h2>Etching: acute corners and same-net gaps (MF1, MF2)</h2>"
                 f"<p class=mut>Ball radius {p['ball_mm']:g} mm, measured at {p['res_mm'] * 1000:g} um. A corner is reported below "
                 f"{p['judged_below_deg']:g} deg (limit {p['limit_deg']:g} deg from a {p['contact_angle_deg']:g} deg contact angle, less the "
                 f"measurement margin) and warns from a flush time of {p['flush_warn_s']:g} s. Flush time: diffusion out of the pocket at "
                 f"D = {p['diffusivity_m2s']:g} m2/s, a rank, not a process prediction.</p>"
                 "<div class=wrap><table><tr><th>Layer</th><th class=n>candidates</th><th class=n>acute corners</th><th class=n>slow to flush</th><th class=n>sharpest (deg)</th>"
                 "<th class=n>longest flush (s)</th><th class=n>narrow gaps</th><th class=n>narrowest (mm)</th></tr>")
        g2 = rules.get("same_net_gap", {}).get("metrics", {})
        for lname, m in r["metrics"].items():
            mg = g2.get(lname, {})
            H.append(f"<tr><td>{e(lname)}</td><td class=n>{m['candidates']}</td><td class=n>{m['wedges']}</td><td class=n>{m['slow']}</td><td class=n>{num(m['min_angle_deg'], 0)}</td>"
                     f"<td class=n>{num(m['max_flush_s'])}</td><td class=n>{mg.get('gaps', '')}</td><td class=n>{num(mg.get('min_width_mm'), 3)}</td></tr>")
        H.append("</table></div>")
    r = rules.get("reflow")
    if r:
        p = r["params"]
        H.append("<h2>Reflow, board level (MF7)</h2>"
                 f"<p class=mut>Oven air by zone: {e(', '.join(f'{z:g}' for z in p['zones_c']))} C for {p['zone_s']:g} s each, then "
                 f"{e(', '.join(f'{z:g}' for z in p['cool_c']))} C; spike zones shifted {p['shift_c']:+g} C so the median peak is "
                 f"{r['summary']['median_peak_c']:g} C (target {e(str(p['target_c']))} C, {e(p['target_source'])}). Film coefficient "
                 f"{p['h_top']:g} / {p['h_bot']:g} W/m2K ({e(p['h_source'])}). Limits: {e(json.dumps(p['limits']))} ({e(p['limits_source'])}). "
                 f"Peak spread over the parts: {e(str(r['summary']['spread_c']))} C.</p>")
        if r.get("image"):
            H.append(f"<p><img src='{e(r['image'])}' alt='peak temperature map'></p>")
        rows = sorted(r["parts"], key=lambda x: x["peak_c"])
        show = rows[:15] + ([None] if len(rows) > 30 else []) + (rows[-15:] if len(rows) > 30 else rows[15:])
        H.append("<div class=wrap><table><tr><th>Part</th><th>Side</th><th class=n>peak (C)</th><th class=n>above liquidus (s)</th>"
                 "<th class=n>ramp (C/s)</th><th class=n>lag (s)</th><th>Body</th></tr>")
        for x in show:
            if x is None:
                H.append(f"<tr><td colspan=7 class=mut>... {len(rows) - 30} parts between ...</td></tr>")
                continue
            H.append(f"<tr><td>{e(x['ref'])}</td><td>{x['side']}</td><td class=n>{x['peak_c']:.1f}</td><td class=n>{x['tal_s']:.0f}</td>"
                     f"<td class=n>{x['ramp_c_s']:.2f}</td><td class=n>{x['lag_s']:.0f}</td><td class=mut>{e(x['body'])}</td></tr>")
        H.append("</table></div><p class=mut>The 15 coldest and 15 hottest parts; all parts are in dfm.json.</p>")
    r, rb = rules.get("tombstone"), rules.get("pad_balance")
    if r or rb:
        H.append("<h2>Chips: pad balance and tombstoning (MF6, MF8)</h2>")
        bal = {x["ref"]: x for x in (rb or {}).get("parts", [])}
        if r:
            H.append(f"<p class=mut>Time between the two joints reaching liquidus in a ramp (independent of the ramp rate). "
                     f"Reported from {r['params']['dt_warn_s']:g} s (proposed limit; calibrate against the assembler's defect data).</p>")
        rows = (r or {}).get("parts") or [{"ref": k, "package": v["package"], "dt_s": None} for k, v in bal.items()]
        H.append("<div class=wrap><table><tr><th>Part</th><th>Package</th><th class=n>dt (s)</th><th>melts first</th><th>melts last</th>"
                 "<th class=n>copper ratio</th><th>copper leaving the pads (mm)</th></tr>")
        for x in rows[:40]:
            bx = bal.get(x["ref"], {})
            ent = ", ".join(f"{q['pad']}: {q['entry_mm']:.2f}" + (" +via in pad" if q["vias_in_pad"] else "")
                            + (f" +{q['plane_vias']} plane via" if q["plane_vias"] else "") for q in bx.get("pads", []))
            H.append(f"<tr><td>{e(x['ref'])}</td><td>{e(x['package'])}</td><td class=n>{num(x.get('dt_s'), 3)}</td>"
                     f"<td>{e(str(x.get('first_pad', '')))} {e(str(x.get('first_net', '')))}</td>"
                     f"<td>{e(str(x.get('last_pad', '')))} {e(str(x.get('last_net', '')))}</td>"
                     f"<td class=n>{num(bx.get('ratio'), 1)}</td><td class=mut>{e(ent)}</td></tr>")
        H.append(f"</table></div><p class=mut>The {min(40, len(rows))} chips with the largest difference; all are in dfm.json.</p>")
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
    H.append("<p class=mut>Clearance between different nets, annular rings, copper slivers and starved thermals are checked by KiCad's "
             "DRC (<code>openpcbsi evaluate</code>). Information findings are in dfm.json.</p></body></html>")
    with open(path, "w") as fh:
        fh.write("\n".join(H))
