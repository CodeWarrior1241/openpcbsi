"""G2 net extraction, G3 batch simulation (linear driver), G5 what-if sweeps, G7 channel masks.

G2: each point-to-point net in scope is walked from one pad to the other through its tracks and
vias. Every track run becomes a transmission line with the 2D solver's (G6) impedance, delay and
loss for its layer and width (differential pairs: the odd-mode parameters of the pair geometry).
Vias are lumped from G9, wide pads are shunt capacitance to the plane under them. Output per net:
a SPICE subcircuit (lossless T lines + lumped vias, for ngspice) and a Touchstone file from the
ABCD cascade of the lossy lines.

G3: drivers here are linear (Thevenin source with a rise time), receivers a termination or a
capacitance. IBIS buffers need ngspice with KIBIS, which is not installed on this machine; the
results carry that limitation. Metrics: overshoot, ringback, monotonic edge, delay to 50 %, and a
suggested series termination.

G7: nets joined through series capacitors are cascaded into a channel and checked against masks
(SDD11, insertion-loss deviation, IL at Nyquist); the eye comes from a peak-distortion analysis of
the pulse response."""
import cmath
import heapq
import json
import math
import os
import re
from collections import defaultdict

import numpy as np

from . import xsec
from .validate import status

EPS0 = 8.8541878128e-12


class LineLib:
    """Line parameters per (layer, w, gap) from the 2D solver, cached for the run."""

    def __init__(self, ctx):
        self.ctx = ctx
        self.cache = {}

    def get(self, li, w, gap=None):
        key = (li, round(w, 3), round(gap, 3) if gap else None)
        if key in self.cache:
            return self.cache[key]
        b, xc = self.ctx.b, self.ctx.cfg["xsec"]
        refs = [b.layers[r] for r in self.ctx.ref_layers(li)]
        freqs = [0.1e9, 1e9, 5e9, 12.5e9, 25e9]
        sec = xsec.build_section(b.stack, b.layers[li], refs, [w] if not gap else [w, w], [] if not gap else [gap], xc, xc["mask_er"])
        r = xsec.solve_lines(sec, max(xc["res_um"], 5.0), freqs, xc["roughness_um"])
        if gap:
            z = r["freqs"][1]["Zdiff"]
            td = r["freqs"][1]["delay_diff_ps_mm"]
            loss = [fr["loss_diff_db_mm"] for fr in r["freqs"]]
        else:
            z = r["freqs"][1]["Z0"]
            td = r["freqs"][1]["delay_ps_mm"]
            loss = [fr["loss_db_mm"] for fr in r["freqs"]]
        self.cache[key] = {"Z": z, "td_ps_mm": td, "f": np.array(freqs), "loss_db_mm": np.array(loss)}
        return self.cache[key]


def line_abcd(lp, length_mm, f):
    a = np.interp(np.log(np.maximum(f, 1e6)), np.log(lp["f"]), np.log(np.maximum(lp["loss_db_mm"], 1e-9)))
    alpha = np.exp(a) * length_mm / 8.686          # Np
    beta = 2 * np.pi * f * lp["td_ps_mm"] * 1e-12 * length_mm
    g = alpha + 1j * beta
    Z = lp["Z"]
    return np.array([[np.cosh(g), Z * np.sinh(g)], [np.sinh(g) / Z, np.cosh(g)]])


def series_abcd(Z):
    one = np.ones_like(Z)
    return np.array([[one, Z], [0 * one, one]])


def shunt_abcd(Y):
    one = np.ones_like(Y)
    return np.array([[one, 0 * one], [Y, one]])


def cascade(mats):
    out = None
    for m in mats:
        out = m if out is None else np.einsum("ijf,jkf->ikf", out, m)
    return out


def abcd_to_s(M, z0):
    A, B, C, D = M[0, 0], M[0, 1], M[1, 0], M[1, 1]
    den = A + B / z0 + C * z0 + D
    return {"S11": (A + B / z0 - C * z0 - D) / den, "S21": 2 / den,
            "S12": 2 * (A * D - B * C) / den, "S22": (-A + B / z0 - C * z0 + D) / den}


# ----------------------------------------------------------------------------------------------
# G2: path extraction

def net_path(ctx, net):
    """Ordered elements from one end pad to the other: ('line', li, w, L, k) / ('via', li_a, li_b, v) / ('pad', ...)."""
    b = ctx.b
    T, V = b.tracks, b.vias
    ks = b.tracks_of(net)
    pads = [p for p in b.pads_of(net) if p["kind"] in ("smd", "conn", "pth")]
    if len(pads) != 2 or not len(ks):
        return None, "needs exactly two pads and tracks" if len(pads) != 2 else "no tracks"
    tol = 0.002
    from .geom import cluster
    n_ = len(ks)
    lab = cluster(np.r_[T["x1"][ks], T["x2"][ks]] + 1000.0 * np.r_[T["layer"][ks], T["layer"][ks]],
                  np.r_[T["y1"][ks], T["y2"][ks]], 0.005)
    rep = {}

    def key_of(i, e, li, x, y):
        c = lab[i + e * n_]
        if c not in rep:
            rep[c] = (li, round(x / tol), round(y / tol))
        return rep[c]
    adj = defaultdict(list)
    for i, k in enumerate(ks):
        li = int(T["layer"][k])
        a, c = key_of(i, 0, li, T["x1"][k], T["y1"][k]), key_of(i, 1, li, T["x2"][k], T["y2"][k])
        adj[a].append((c, float(T["len"][k]), ("line", li, float(T["w"][k]), float(T["len"][k]), int(k))))
        adj[c].append((a, float(T["len"][k]), ("line", li, float(T["w"][k]), float(T["len"][k]), int(k))))
    nodes = list(adj.keys())
    vtl = ctx.via_track_layers()
    for v in b.vias_of(net):
        x, y = V["x"][v], V["y"][v]
        at = [n for n in nodes if math.hypot(n[1] * tol - x, n[2] * tol - y) <= V["d"][v] / 2 + 0.06]
        for n1 in at:
            for n2 in at:
                if n1[0] != n2[0]:
                    adj[n1].append((n2, 0.0, ("via", n1[0], n2[0], int(v))))
    from .pdn import _pad_dist
    term = []
    for p in pads:
        tn = [n for n in nodes if n[0] in p["layer_ids"] and _pad_dist(p, n[1] * tol, n[2] * tol) <= 0.07]
        if not tn:
            # through-hole / via-in-pad: any node at the pad on any layer
            tn = [n for n in nodes if _pad_dist(p, n[1] * tol, n[2] * tol) <= 0.07]
        term.append((p, tn))
    if not term[0][1] or not term[1][1]:
        return None, "pad not reached by a track"
    src = {n: 0.0 for n in term[0][1]}
    dist, prev = {}, {}
    pq = [(0.0, n) for n in src]
    heapq.heapify(pq)
    while pq:
        d, n = heapq.heappop(pq)
        if n in dist:
            continue
        dist[n] = d
        for m, w, el in adj[n]:
            if m not in dist:
                nd = d + w + (1e-6 if el[0] == "via" else 0)
                if nd < dist.get(m, math.inf) and (m not in prev or nd < prev[m][2]):
                    prev[m] = (n, el, nd)
                heapq.heappush(pq, (nd, m))
    ends = [n for n in term[1][1] if n in dist]
    if not ends:
        return None, "pads not connected through tracks"
    n = min(ends, key=lambda e: dist[e])
    els = []
    while n in prev and n not in src:
        p_, el, _ = prev[n]
        els.append(el)
        n = p_
    els.reverse()
    return {"from": f"{term[0][0]['ref']}.{term[0][0]['num']}", "to": f"{term[1][0]['ref']}.{term[1][0]['num']}",
            "pads": [term[0][0], term[1][0]], "elements": els, "length_mm": round(dist[min(ends, key=lambda e: dist[e])], 3)}, None


def merge_lines(els):
    out = []
    for e in els:
        if e[0] == "line" and out and out[-1][0] == "line" and out[-1][1] == e[1] and abs(out[-1][2] - e[2]) < 1e-4:
            out[-1] = ("line", e[1], e[2], out[-1][3] + e[3], out[-1][4])
        else:
            out.append(e)
    return out


def pad_cap(ctx, pad, w_track):
    """Excess capacitance of a pad wider than its track over the plane below it."""
    b = ctx.b
    if not pad["layer_ids"]:
        return 0.0
    li = pad["layer_ids"][0]
    refs = ctx.ref_layers(li)
    if not refs:
        return 0.0
    h = min(b.stack.gap(b.layers[li], b.layers[r]) for r in refs) * 1e-3
    er = b.stack.er_between(b.layers[li], b.layers[refs[0]])
    A = max(pad["sx"] * pad["sy"] - w_track * min(pad["sx"], pad["sy"]), 0) * 1e-6
    return EPS0 * er * A / h


def net_abcd(ctx, lib, path, f, gap=None, vias_rep=None):
    from . import vias as vmod
    b = ctx.b
    mats, spice = [], []
    els = merge_lines(path["elements"])
    w0 = next((e[2] for e in els if e[0] == "line"), 0.12)
    cp = [pad_cap(ctx, p, w0) for p in path["pads"]]
    diff = gap is not None
    k = 0.5 if diff else 1.0       # differential: two shunt caps in series
    mats.append(shunt_abcd(1j * 2 * np.pi * f * cp[0] * k))
    spice.append(f"C_pad1 n0 0 {cp[0] * k:.4e}")
    node = 0
    for e in els:
        if e[0] == "line":
            lp = lib.get(e[1], e[2], gap)
            mats.append(line_abcd(lp, e[3], f))
            spice.append(f"T{node} n{node} 0 n{node + 1} 0 Z0={lp['Z']:.2f} TD={lp['td_ps_mm'] * e[3]:.2f}p  * {b.layers[e[1]]} w={e[2]} L={e[3]:.3f}mm")
            node += 1
        else:
            layers = f"{b.layers[min(e[1], e[2])]}-{b.layers[max(e[1], e[2])]}"
            vm = vmod.lookup(vias_rep, layers) or {"C_pF": 0.3, "Lloop_nH": 0.4, "L_nH": 0.6}
            L = vm["L_nH"] * 1e-9 * (2.0 if diff else 1.0)    # differential: both barrels in series
            C = vm["C_pF"] * 1e-12 * k
            mats += [shunt_abcd(1j * 2 * np.pi * f * C / 2), series_abcd(1j * 2 * np.pi * f * L), shunt_abcd(1j * 2 * np.pi * f * C / 2)]
            spice += [f"C_v{node}a n{node} 0 {C / 2:.4e}", f"L_v{node} n{node} n{node + 1} {L:.4e}", f"C_v{node}b n{node + 1} 0 {C / 2:.4e}"]
            node += 1
    mats.append(shunt_abcd(1j * 2 * np.pi * f * cp[1] * k))
    spice.append(f"C_pad2 n{node} 0 {cp[1] * k:.4e}")
    return cascade(mats), spice, node


# ----------------------------------------------------------------------------------------------
# G3: linear transient from the frequency response

def safe_name(net):
    """A file and SPICE name for a net: KiCad's hierarchical names contain '/', and others '(' or '{'."""
    return re.sub(r"[^A-Za-z0-9_.+\-]+", "_", net).strip("_") or "net"


def step_response(M, f, rs, zl, rise_ps, tmax_ns):
    """Receiver voltage for a 0 -> 1 V step behind rs with a linear ramp of rise_ps (10-90 % ~ 0.8 tr)."""
    A, B, C, D = M[0, 0], M[0, 1], M[1, 0], M[1, 1]
    H = 1.0 / (A + B / zl + rs * (C + D / zl))
    H = np.where(np.isfinite(H), H, 0)
    n = 2 * (len(f) - 1)
    df = f[1] - f[0]
    dt = 1.0 / (n * df)
    h = np.fft.irfft(H, n) / dt       # impulse response
    s = np.cumsum(h) * dt             # step response
    m = max(1, int(round(rise_ps * 1e-12 / dt)))
    s = np.convolve(s, np.ones(m) / m)[: len(s)]
    t = np.arange(len(s)) * dt
    keep = t <= tmax_ns * 1e-9
    return t[keep], s[keep], dt


def metrics(t, v, rise_ps, vfinal):
    vf = float(np.mean(v[-max(5, len(v) // 20):]))     # settled value at the end of the window
    i50 = int(np.argmax(v >= 0.5 * vf)) if (v >= 0.5 * vf).any() else None
    delay = (t[i50] - 0.5 * rise_ps * 1e-12) * 1e12 if i50 is not None else None
    pk = int(np.argmax(v))
    over = (v[pk] - vf) / vf * 100
    i90 = int(np.argmax(v >= 0.9 * vf)) if (v >= 0.9 * vf).any() else None
    # ringback: how far the signal falls back below the final value after its first peak
    ring = max(0.0, (vf - v[pk:].min()) / vf * 100) if (i90 is not None and pk >= i90) else None
    i10 = int(np.argmax(v >= 0.1 * vf)) if (v >= 0.1 * vf).any() else None
    mono = bool(np.all(np.diff(v[i10:i90 + 1]) >= -1e-3 * vf)) if (i10 is not None and i90 is not None and i90 > i10) else False
    return {"delay_ps": round(delay, 1) if delay is not None else None, "overshoot_pct": round(over, 2),
            "ringback_pct": round(ring, 2) if ring is not None else None, "monotonic": mono}


def run(ctx, cache_dir, out_dir, log=print, vias=None):
    b, cfg = ctx.b, ctx.cfg
    p = cfg["si"]
    lib = LineLib(ctx)
    nets = p["nets"] or [n for n in b.net_names[1:] if cfg.class_of(n)]
    sdir = os.path.join(out_dir, "nets")
    os.makedirs(sdir, exist_ok=True)
    f = np.linspace(0, p.get("fmax_ghz", 40.0) * 1e9, 4097)
    f[0] = 1e6
    rows, done = [], set()
    for net in sorted(nets):
        if net in done:
            continue
        partner = ctx.partner.get(net)
        diff = partner is not None
        if diff:
            pos, neg = sorted([net, partner])
            done |= {pos, neg}
            gap = (ctx_gap(ctx, pos) or 0.15)
            name = pos
        else:
            gap, name = None, net
            done.add(net)
        path, err = net_path(ctx, name)
        if path is None:
            rows.append({"net": name, "problem": err})
            continue
        M, spice, nn = net_abcd(ctx, lib, path, f, gap, vias)
        cls = cfg.class_of(name)
        rise = cfg.class_param(cls, "rise_ps", p["driver"]["rise_ps"])
        z_ref = 100.0 if diff else 50.0
        term = cfg.class_param(cls, "term_ohm", z_ref)       # LVDS / GTH: internal 100 ohm diff / 50 ohm
        rs = cfg.class_param(cls, "driver_r_ohm", z_ref if diff else p["driver"]["r_ohm"])
        zl = term if term else 1.0 / (1j * 2 * np.pi * f * p["load_pf"] * 1e-12)
        t, v, dt = step_response(M, f, rs, zl, rise, 5.0)
        vfin = term / (rs + term) if term else 1.0
        m = metrics(t, v, rise, vfin)
        S = abcd_to_s(M, z_ref)
        nyq = cfg.class_param(cls, "nyquist_ghz", 1.0) * 1e9
        iq = int(np.argmin(abs(f - nyq)))
        row = {"net": name + ("/" + neg if diff else ""), "from": path["from"], "to": path["to"], "length_mm": path["length_mm"],
               "elements": len(merge_lines(path["elements"])), "vias": sum(1 for e in path["elements"] if e[0] == "via"),
               "rise_ps": rise, "driver_ohm": rs, "term_ohm": term, **m,
               "IL_nyq_db": round(20 * math.log10(abs(S["S21"][iq])), 3), "RL_nyq_db": round(20 * math.log10(abs(S["S11"][iq]) + 1e-12), 2)}
        if m["overshoot_pct"] > 10 and not term:
            row["suggest"] = f"series R ~ {max(0, lib.get(*first_line(path), gap)['Z'] - rs):.0f} ohm at the driver"
        rows.append(row)
        fn = safe_name(name)  # hierarchical KiCad nets: "/sheet/NET"
        row["files"] = f"nets/{fn}.sp, nets/{fn}.s2p"
        with open(os.path.join(sdir, f"{fn}.sp"), "w") as fh:
            fh.write(f"* openpcbsi G2: {name}{' (differential, odd-mode equivalent)' if diff else ''} {path['from']} -> {path['to']}\n")
            fh.write(f".subckt {fn} n0 n{nn}\n" + "\n".join(spice) + "\n.ends\n")
        write_touchstone(os.path.join(sdir, f"{fn}.s2p"), f, S, z_ref)
    for r in rows[:6]:
        if "problem" in r:
            continue
        log(f"  {r['net']:32s} {r['length_mm']:6.1f} mm {r['vias']} via  delay {r['delay_ps']} ps  overshoot {r['overshoot_pct']} %"
            f"  IL@Nyq {r['IL_nyq_db']} dB  RL@Nyq {r['RL_nyq_db']} dB")
    probs = [r for r in rows if "problem" in r]
    log(f"  {len(rows) - len(probs)} nets/pairs extracted, {len(probs)} not point-to-point")
    out = {"rule": "si", "spec": "G2/G3", "title": "Net extraction and linear SI", "validation": status(cache_dir, "xsec2d") + " (lines); linear driver",
           "rows": [r for r in rows if "problem" not in r], "problems": probs,
           "metrics": {r["net"]: {k: r[k] for k in ("delay_ps", "overshoot_pct", "IL_nyq_db", "RL_nyq_db") if k in r} for r in rows if "problem" not in r},
           "findings": [], "summary": {"nets": len(rows) - len(probs)},
           "note": "Linear drivers (no IBIS: ngspice is not installed). Lines from the validated 2D solver; vias from G9 closed forms; package and connector not included."}
    # G7 channels
    if cfg["channels"]:
        out["channels"] = [channel(ctx, lib, name, ch, vias, log) for name, ch in cfg["channels"].items()]
    return out


def first_line(path):
    for e in path["elements"]:
        if e[0] == "line":
            return e[1], e[2]
    return 0, 0.12


def ctx_gap(ctx, net):
    """Median P-to-N gap of a pair from the vector geometry (planner)."""
    from .planner import pair_geometries
    if not hasattr(ctx, "_pg"):
        ctx._pg = pair_geometries(ctx)
    b = ctx.b
    ks = b.tracks_of(net)
    if not len(ks):
        return None
    li = int(np.bincount(b.tracks["layer"][ks]).argmax())
    w = round(float(np.median(b.tracks["w"][ks])), 4)
    cand = [(L, g) for (l, ww, g), L in ctx._pg.items() if l == li and abs(ww - w) < 1e-3]
    return max(cand)[1] if cand else None


def write_touchstone(path, f, S, z0):
    with open(path, "w") as fh:
        fh.write(f"! openpcbsi G2\n# HZ S DB R {z0:g}\n")
        for i in range(0, len(f), 8):
            vals = []
            for k in ("S11", "S21", "S12", "S22"):
                s = S[k][i]
                vals += [f"{20 * math.log10(abs(s) + 1e-15):.4f}", f"{math.degrees(cmath.phase(s)):.3f}"]
            fh.write(f"{f[i]:.6e} " + " ".join(vals) + "\n")


# ----------------------------------------------------------------------------------------------
# G7

def channel(ctx, lib, name, ch, vias_rep, log=print):
    """Cascade of pair nets (and series capacitors) with masks and a peak-distortion eye."""
    from .pdn import parse_value
    f = np.linspace(0, ch.get("fmax_ghz", 40.0) * 1e9, 4097)
    f[0] = 1e6
    mats = []
    descr = []
    for item in ch["path"]:
        if item.startswith("cap:"):
            c = parse_value(item[4:]) or 1e-7
            esl = ch.get("cap_esl_nh", 0.4) * 1e-9
            # differential: the two series caps in series around the loop
            Z = 2 * (1 / (1j * 2 * np.pi * f * c) + 1j * 2 * np.pi * f * esl + 0.02)
            mats.append(series_abcd(Z))
            descr.append(item)
            continue
        path, err = net_path(ctx, item)
        if path is None:
            return {"channel": name, "problem": f"{item}: {err}"}
        gap = ctx_gap(ctx, item) or 0.15
        M, _, _ = net_abcd(ctx, lib, path, f, gap, vias_rep)
        mats.append(M)
        descr.append(f"{item} {path['from']}->{path['to']} {path['length_mm']} mm")
    M = cascade(mats)
    S = abcd_to_s(M, 100.0)
    il = -20 * np.log10(np.abs(S["S21"]) + 1e-15)
    rl = 20 * np.log10(np.abs(S["S11"]) + 1e-15)
    br = ch["bitrate_gbps"] * 1e9
    nyq = br / 2
    res = {"channel": name, "path": descr, "bitrate_gbps": ch["bitrate_gbps"],
           "IL_nyq_db": round(float(np.interp(nyq, f, il)), 3), "masks": []}
    m = ch.get("masks", {})
    if "sdd11_db" in m:
        band = f <= m.get("sdd11_fmax_ghz", nyq / 1e9) * 1e9
        worst = float(rl[band].max())
        res["masks"].append({"mask": f"SDD11 < {m['sdd11_db']} dB to {m.get('sdd11_fmax_ghz', nyq / 1e9)} GHz",
                             "worst": round(worst, 2), "margin_db": round(m["sdd11_db"] - worst, 2), "pass": worst < m["sdd11_db"]})
    if "ild_db" in m:
        fm = m.get("ild_fmax_ghz", br / 1e9) * 1e9
        band = (f >= 50e6) & (f <= fm)
        x = f[band] / 1e9
        X = np.c_[np.ones_like(x), np.sqrt(x), x, x * x]
        coef, *_ = np.linalg.lstsq(X, il[band], rcond=None)
        ild = float(np.abs(il[band] - X @ coef).max())
        res["masks"].append({"mask": f"ILD < {m['ild_db']} dB to {fm / 1e9:g} GHz", "worst": round(ild, 3),
                             "margin_db": round(m["ild_db"] - ild, 3), "pass": ild < m["ild_db"]})
    if "il_nyq_max_db" in m:
        res["masks"].append({"mask": f"IL at Nyquist < {m['il_nyq_max_db']} dB", "worst": res["IL_nyq_db"],
                             "margin_db": round(m["il_nyq_max_db"] - res["IL_nyq_db"], 3), "pass": res["IL_nyq_db"] < m["il_nyq_max_db"]})
    # pulse response and peak-distortion eye (linear TX, 50 ohm per leg = 100 ohm diff, matched RX)
    rise = ch.get("rise_ps", 35.0)
    t, s, dt = step_response(M, f, 100.0, 100.0, rise, 20.0)
    ui = 1.0 / br
    n_ui = int(round(ui / dt))
    pulse = s - np.r_[np.zeros(n_ui), s[:-n_ui]]
    pk = int(np.argmax(pulse))
    best = (-1, None)
    heights = []
    for off in range(-n_ui // 2, n_ui // 2):
        i0 = pk + off
        cur = pulse[i0]
        isi = sum(abs(pulse[i0 + k * n_ui]) for k in range(-3, 40) if k != 0 and 0 <= i0 + k * n_ui < len(pulse))
        h = cur - isi
        heights.append(h)
        if h > best[0]:
            best = (h, off)
    hh = np.array(heights)
    width = min(float((hh > 0).sum()) * dt, ui)
    res.update({"eye_height_pct": round(100 * best[0] / 0.5, 1), "eye_height_v": round(2 * best[0], 4), "eye_width_ps": round(width * 1e12, 1),
                "eye_width_ui": round(width / ui, 3), "pulse_peak_v": round(float(pulse[pk]), 4),
                "eye_method": "peak distortion, linear driver (no IBIS-AMI); height in % of a lossless channel's"})
    log(f"  channel {name}: IL@{nyq / 1e9:g} GHz {res['IL_nyq_db']} dB, eye {res['eye_height_pct']} % x {res['eye_width_ui']} UI; "
        + "; ".join(f"{mk['mask']}: {'PASS' if mk['pass'] else 'FAIL'} ({mk['margin_db']:+} dB)" for mk in res["masks"]))
    return res


# ----------------------------------------------------------------------------------------------
# G5: what-if from a topology file

def whatif(ctx, spec, out_dir, log=print):
    """spec: {layer, w, gap?, length_mm, driver_r_ohm, rise_ps, term_ohm|load_pf, sweep: {key: [values]}}"""
    import itertools
    lib = LineLib(ctx)
    b = ctx.b
    f = np.linspace(0, 40e9, 4097)
    f[0] = 1e6
    keys = list(spec.get("sweep", {}).keys())
    vals = [spec["sweep"][k] for k in keys]
    rows = []
    for combo in itertools.product(*vals) if keys else [()]:
        c = dict(spec)
        c.update(dict(zip(keys, combo)))
        li = b.L[c["layer"]]
        lp = lib.get(li, c["w"], c.get("gap"))
        M = line_abcd(lp, c["length_mm"], f)
        rs = c.get("driver_r_ohm", 40.0)
        if c.get("series_r_ohm"):
            M = cascade([series_abcd(np.full_like(f, c["series_r_ohm"], dtype=complex)), M])
        term = c.get("term_ohm")
        zl = term if term else 1.0 / (1j * 2 * np.pi * f * c.get("load_pf", 2.0) * 1e-12)
        t, v, dt = step_response(M, f, rs, zl, c.get("rise_ps", 100.0), 10.0)
        vfin = term / (rs + term) if term else 1.0
        m = metrics(t, v, c.get("rise_ps", 100.0), vfin)
        rows.append(dict({k: c.get(k) for k in keys}, Z=round(lp["Z"], 2), **m))
    with open(os.path.join(out_dir, "whatif.json"), "w") as fh:
        json.dump(rows, fh, indent=1)
    for r in rows:
        log("  " + "  ".join(f"{k}={r[k]}" for k in keys) + f"  Z={r['Z']}  overshoot {r['overshoot_pct']} %  ringback {r['ringback_pct']} %  delay {r['delay_ps']} ps")
    return rows
