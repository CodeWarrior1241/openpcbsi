"""Timing (#25): propagation delay per net, skew per matching group and per pair.

delay = sum over the net's tracks of length x sqrt(eps_eff) / c, with eps_eff from the 2D solver for
the track's layer and width (odd-mode for the tracks of a pair, where the pair geometry is known).
Vias are not included in the delay (a via adds a few ps: see the G9 via table).
A .kicad_dru skew rule is judged as KiCad judges it, in length: tracks plus via spans (with the project's
"use stackup height" setting), over the nets the rule covers or per pair with (within_diff_pairs). Pairs it
covers are not judged again in ps.

    [timing.groups]
    DQ0 = ["DDR_DQ0", "DDR_DQ1", "DDR_DQS0_P"]   # nets (fnmatch patterns allowed)
    [timing]
    skew_ps = 10.0                                  # limit for group skew (default: none, reported only)
    pair_skew_ps = 1.0                              # limit for P/N skew (default: the standard's, else 2 ps)
"""
import fnmatch
import math
from collections import defaultdict

from . import xsec
from .drc import finding

C = 299792458.0


def run(ctx, cache_dir, log=print):
    b, cfg = ctx.b, ctx.cfg
    tc = cfg.data.get("timing", {}) or {}
    groups = tc.get("groups", {}) or {}
    names = b.net_names[1:]
    gnets = {g: sorted({n for p in pats for n in names if fnmatch.fnmatchcase(n, p)}) for g, pats in groups.items()}
    nets = sorted({n for ns in gnets.values() for n in ns} | {n for n in names if cfg.class_of(n)})
    xc = cfg["xsec"]
    ee_cache = {}

    def eeff(li, w, pair_gap=None):
        key = (li, round(w, 4), pair_gap)
        if key not in ee_cache:
            lname = b.layers[li]
            refs = [b.layers[r] for r in ctx.ref_layers(li)]
            try:
                if pair_gap:
                    sec = xsec.build_section(b.stack, lname, refs, [w, w], [pair_gap], xc, xc["mask_er"])
                    ee_cache[key] = xsec.solve_lines(sec, xc["res_um"], (1e9,))["freqs"][0]["eeff_odd"]
                else:
                    sec = xsec.build_section(b.stack, lname, refs, [w], [], xc, xc["mask_er"])
                    ee_cache[key] = xsec.solve_lines(sec, xc["res_um"], (1e9,))["freqs"][0]["eeff"]
            except Exception as e:
                log(f"  timing: no eps_eff for {lname} w={w}: {e}")
                ee_cache[key] = None
        return ee_cache[key]

    # pair gap per pair net and layer (median of the planner's pair geometries)
    from .planner import pair_geometries
    where = {}
    pair_geometries(ctx, where)
    gap_of = {}
    for (li, w, g), wh in where.items():
        for n in wh["nets"]:
            gap_of.setdefault((n, li, round(w, 4)), g)
    T = b.tracks
    rows, finds = [], []
    for n in nets:
        if n not in b.net_id:
            continue
        L, tsum, eesum = 0.0, 0.0, 0.0
        for t in b.tracks_of(n):
            li, w, l = int(T["layer"][t]), float(T["w"][t]), float(T["len"][t])
            e = eeff(li, w, gap_of.get((n, li, round(w, 4))))
            if e is None:
                continue
            L += l
            eesum += e * l
            tsum += l * 1e-3 * math.sqrt(e) / C * 1e12
        if L <= 0:
            continue
        rows.append({"net": n, "length_mm": round(L, 3), "eps_eff": round(eesum / L, 4), "delay_ps": round(tsum, 3)})
    by = {r["net"]: r for r in rows}
    gout = {}
    lim = tc.get("skew_ps")
    for g, ns in gnets.items():
        ds = [(by[n]["delay_ps"], n) for n in ns if n in by]
        if len(ds) < 2:
            continue
        sk = max(ds)[0] - min(ds)[0]
        gout[g] = {"nets": [n for _, n in ds], "skew_ps": round(sk, 3), "slowest": max(ds)[1], "fastest": min(ds)[1]}
        if lim is not None and sk > lim:
            finds.append(finding("timing", "error", f"group {g}: skew {sk:.1f} ps above {lim:g} ps ({max(ds)[1]} vs {min(ds)[1]})",
                                 max(ds)[1], "", value=round(sk, 3), limit=lim, unit="ps", group=g))
    # .kicad_dru skew rules (#2), judged on KiCad's own measurand
    dru = _dru_skew(ctx, names)
    dru_nets = {n for r in dru for n in r["nets"]}
    for r in dru:
        if r["skew_mm"] > r["limit_mm"]:
            finds.append(finding("timing", "warning", f"{r['what']}: length skew {r['skew_mm']:.3f} mm above {r['limit_mm']:g} mm "
                                 f"({r['source']}; {r['measurand']})", r["longest"], "", value=r["skew_mm"],
                                 limit=r["limit_mm"], unit="mm", limit_source=r["source"]))
    # intra-pair skew
    from . import standards
    done = set()
    for pn, nn in ctx.partner.items():
        if pn in done or pn not in by or nn not in by:
            continue
        done |= {pn, nn}
        if pn in dru_nets and nn in dru_nets:          # a .kicad_dru skew rule covers the pair: it wins
            continue
        sk = abs(by[pn]["delay_ps"] - by[nn]["delay_ps"])
        plim = tc.get("pair_skew_ps")
        src = "[timing] pair_skew_ps"
        if plim is None:
            st = standards.of_net(pn)
            if st:
                per_mm = by[pn]["delay_ps"] / max(by[pn]["length_mm"], 1e-6)
                plim, src = st[1]["skew_max_mm"] * per_mm, f"standard {st[0]} ({st[1]['skew_max_mm']} mm)"
            else:
                plim, src = 2.0, "default 2 ps"
        if sk > plim:
            finds.append(finding("timing", "warning", f"pair {pn}/{nn}: skew {sk:.2f} ps above {plim:.2f} ps ({src})",
                                 pn, "", value=round(sk, 3), limit=round(plim, 3), unit="ps"))
    log(f"  {len(rows)} nets timed, {len(gout)} groups")
    from .validate import status
    return {"rule": "timing", "spec": "timing", "title": "Propagation delay and skew",
            "validation": status(cache_dir, "xsec2d") + " (2D eps_eff)",
            "nets": rows, "groups": gout, "rows": rows, "dru_skew": dru, "metrics": {r["net"]: {"delay_ps": r["delay_ps"]} for r in rows},
            "findings": finds, "summary": {"nets": len(rows), "groups": len(gout)}}


def kicad_layer_distance(stack, a, b, layers):
    """BOARD_STACKUP::GetLayerDistance: the stack between copper layers a and b, with an outer layer
    counted whole and an inner end layer counted half (centre to centre plus half of each outer end)."""
    if a == b:
        return 0.0
    d = stack.barrel(a, b)
    for name in (a, b):
        if name in (layers[0], layers[-1]):
            d += stack.cu(name)["t"] / 2
    return d


def _dru_skew(ctx, names):
    """The .kicad_dru skew rules as KiCad judges them (drc_test_provider_matched_length): over every net
    the rule covers, or per pair with (within_diff_pairs); a net's length is its tracks plus, when the
    project uses stackup height for lengths (the default), each via's span between the outermost layers
    it joins tracks or pads on."""
    b, cfg = ctx.b, ctx.cfg
    proj = getattr(cfg, "project", None)
    rules = proj.skew_rules() if proj is not None else []
    if not rules:
        return []
    T, V = b.tracks, b.vias
    vtl = ctx.via_track_layers()
    use_h = proj.via_height_in_length
    cache = {}

    def length(n):
        if n not in cache:
            tl = float(T["len"][b.tracks_of(n)].sum())
            vl = 0.0
            if use_h:
                for v in b.vias_of(n):
                    ls = sorted(vtl[v])
                    if len(ls) >= 2:
                        vl += kicad_layer_distance(b.stack, b.layers[ls[0]], b.layers[ls[-1]], b.layers)
            cache[n] = (tl, vl)
        return cache[n]

    measurand = "tracks plus via spans, as KiCad measures it" if use_h else "tracks only (the project does not use stackup height)"
    out = []
    for name, mx, per_pair, pred in rules:
        covered = sorted(n for n in names if pred(n, b.kicad_class.get(n)) and len(b.tracks_of(n)))
        if per_pair:
            groups, seen = [], set()
            for n in covered:
                q = ctx.partner.get(n)
                if q in covered and n not in seen:
                    seen |= {n, q}
                    groups.append(sorted([n, q]))
        else:
            groups = [covered] if len(covered) >= 2 else []
        for grp in groups:
            ls = {n: sum(length(n)) for n in grp}
            hi, lo = max(ls, key=ls.get), min(ls, key=ls.get)
            out.append({"rule": name, "nets": grp, "what": f"rule '{name}' ({', '.join(grp) if len(grp) <= 4 else f'{len(grp)} nets'})",
                        "skew_mm": round(ls[hi] - ls[lo], 4), "limit_mm": mx, "longest": hi, "shortest": lo,
                        "lengths_mm": {n: {"tracks": round(length(n)[0], 4), "vias": round(length(n)[1], 4)} for n in grp},
                        "source": f".kicad_dru rule '{name}' skew max", "measurand": measurand})
    return out
