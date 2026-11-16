"""Judgements for checks that compute a result against a limit (#5): every out-of-limit result becomes a
finding in the same form as the G1 rules (drc.finding), so evaluate.py scores it.

Limits come from the signal class (config or detection), then the check's config section, then the
defaults below; each finding states which.
"""
from .drc import finding

DEFAULTS = {
    "via_z_tol_pct": 20.0,        # differential via impedance (2 x Z_se) against the class Zdiff
    "overshoot_pct": 15.0,
    "ringback_pct": 10.0,
    "rl_nyq_db": -10.0,           # return loss at Nyquist worse than this is a warning
    "tj_max_c": 100.0,            # junction limit when neither the config nor a footprint field gives one
    "emc_score": 5.0,             # ranked radiator score above which the EMC estimate warns
    "edge_gap_mm": 20.0,          # unstitched board edge above which the EMC estimate warns
}


def _lim(cfg, cls, key, section=None):
    if cls:
        v = cfg.class_param(cls, key)
        if v is not None:
            return v, f"class {cls}"
    if section and cfg.get(section, key) is not None:
        return cfg.get(section, key), f"[{section}] {key}"
    return DEFAULTS[key], "default"


def impedance(rows, cfg):
    """Rows of the planner carry `targets` [(class, target, tol)] and a representative position."""
    out = []
    for r in rows:
        z = r.get("Zdiff") if r["kind"] == "pair" else r.get("Z0")
        for t in r.get("targets", []):
            dz = z - t["target"]
            r.setdefault("checks", []).append({"class": t["class"], "target": t["target"], "tol": t["tol"],
                                                "within": abs(dz) <= t["tol"]})
            lo, hi = r.get("Zdiff_min" if r["kind"] == "pair" else "Z0_min"), r.get("Zdiff_max" if r["kind"] == "pair" else "Z0_max")
            if abs(dz) > t["tol"]:
                sev = "error"
                msg = (f"{'Zdiff' if r['kind'] == 'pair' else 'Z0'} {z:.1f} ohm against {t['target']:g} +/- {t['tol']:g} ohm "
                       f"(class {t['class']}): w={r['w']}" + (f" gap={r['gap']}" if r.get("gap") else "") + f", {r['length_mm']} mm")
            elif lo is not None and (lo < t["target"] - t["tol"] or hi > t["target"] + t["tol"]):
                sev = "warning"
                msg = (f"nominal {z:.1f} ohm is inside {t['target']:g} +/- {t['tol']:g} ohm but the tolerance corners reach "
                       f"{lo:.1f}-{hi:.1f} ohm (class {t['class']})")
            else:
                continue
            out.append(finding("impedance", sev, msg, net=r.get("net_example", ""), layer=r["layer"], x=r.get("x"), y=r.get("y"),
                               value=round(z, 2), limit=t["tol"], unit="ohm", target=t["target"], cls=t["class"],
                               area=r.get("area", "")))
    return out


def thermal(parts, cfg, fields_of):
    out = []
    for p in parts:
        spec = cfg["thermal"]["parts"].get(p["ref"], {})
        lim, src = spec.get("tj_max_c"), "config"
        if lim is None:
            fv = fields_of(p["ref"])
            for k in ("Tj_max", "Tjmax", "TJ_MAX", "Tj max"):
                if k in fv:
                    try:
                        lim, src = float(str(fv[k]).replace("C", "").replace("°", "").strip()), f"field {p['ref']}.{k}"
                    except ValueError:
                        pass
                    break
        if lim is None:
            lim, src = cfg.get("thermal", "tj_max_c") or DEFAULTS["tj_max_c"], "default (assumed)"
        p["tj_max_c"], p["tj_limit_source"] = lim, src
        if p["Tj_est_c"] > lim:
            out.append(finding("thermal", "error", f"{p['ref']} Tj ~ {p['Tj_est_c']:.0f} C above {lim:g} C ({src}) at "
                               f"{p['power_w']} W, theta_jb {p.get('theta_jb')}", value=p["Tj_est_c"], limit=lim, unit="C",
                               part=p["ref"], x=p.get("x"), y=p.get("y")))
    return out


def vias(rows, cfg):
    out = []
    for r in rows:
        cls = {cfg.class_of(n) for n in r.get("nets", [])} - {None}
        for c in sorted(cls):
            zd = cfg.class_param(c, "z_diff")
            if not zd:
                continue
            tol, src = _lim(cfg, c, "via_z_tol_pct", "vias")
            zv = 2 * r["Z_ohm"]
            if abs(zv - zd) / zd * 100 > tol:
                out.append(finding("vias", "warning", f"via {r['layers']} drill {r['drill_mm']} pad {r['pad_mm']} antipad "
                                   f"{r['antipad_mm']}: differential ~ {zv:.0f} ohm against {zd:g} ohm +/- {tol:g} % ({src}), "
                                   f"{r['instances']} instances", net=(r.get("nets") or [""])[0], layer=r["layers"],
                                   value=round(zv, 1), limit=tol, unit="ohm", cls=c))
    return out


def si(rep, cfg):
    out = []
    for r in rep.get("rows", []):
        if "problem" in r:
            continue
        net = r["net"].split("/")[0]
        c = cfg.class_of(net)
        for key, label in (("overshoot_pct", "overshoot"), ("ringback_pct", "ringback")):
            lim, src = _lim(cfg, c, key, "si")
            if r.get(key, 0) > lim:
                out.append(finding("si", "warning", f"{r['net']} {r['from']} -> {r['to']}: {label} {r[key]:.1f} % above {lim:g} % ({src})",
                                   net=net, value=r[key], limit=lim, unit="%"))
        if r.get("monotonic") is False:
            out.append(finding("si", "warning", f"{r['net']}: non-monotonic edge at the receiver {r['to']}", net=net))
        lim, src = _lim(cfg, c, "rl_nyq_db", "si")
        if r.get("RL_nyq_db") is not None and r["RL_nyq_db"] > lim:
            out.append(finding("si", "warning", f"{r['net']}: return loss at Nyquist {r['RL_nyq_db']:.1f} dB above {lim:g} dB ({src})",
                               net=net, value=r["RL_nyq_db"], limit=lim, unit="dB"))
    for ch in rep.get("channels", []):
        if ch.get("problem"):
            out.append(finding("si", "warning", f"channel {ch['channel']} not evaluated: {ch['problem']}", channel=ch["channel"]))
        for m in ch.get("masks", []):
            if not m.get("pass"):
                out.append(finding("si", "error", f"channel {ch['channel']}: {m['mask']} fails, worst {m['worst']} "
                                   f"(margin {m['margin_db']} dB)", net=(ch.get("nets") or [""])[0], value=m["worst"],
                                   limit=m.get("limit"), unit="dB", channel=ch["channel"]))
    return out


def emc(rep, cfg):
    out = []
    lim, src = cfg.get("emc", "score") or DEFAULTS["emc_score"], "[emc] score" if cfg.get("emc", "score") else "default"
    for r in rep.get("rows", []):
        if r["score"] > lim:
            out.append(finding("emc", "warning", f"likely radiator: {r['source']} near {r['f_mhz']} MHz (score {r['score']} > {lim:g}, {src})",
                               value=r["score"], limit=lim, unit="score"))
    eg = rep.get("edge_gap_mm")
    lim, src = cfg.get("emc", "edge_gap_mm") or DEFAULTS["edge_gap_mm"], "default"
    if eg and eg > lim:
        out.append(finding("emc", "warning", f"{eg:.0f} mm of board edge without a stitching via (> {lim:g} mm, {src})",
                           value=eg, limit=lim, unit="mm"))
    return out
