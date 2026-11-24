"""EV5: a first configuration detected from the board alone (netlist, footprints, zones, netclasses).

Used when `check` or `evaluate` runs without `-c`, and by `init`, which writes the same result as
a TOML file to review and refine. Everything detected here is a guess with a stated default; a
reviewed board config replaces it (with `-c`, detection is not used at all).

Detected:
  reference nets   GND-like names (GND, AGND, DGND, VSS, ...), largest first
  pair classes     P/N nets by `[pairs] suffixes`. Serial-lane names (DP, GT, MGT, PCIe, SATA, USB3,
                   HDMI, TMDS, DPHY, CSI, DSI, REFCLK, TX/RX) -> class `serdes`; other pairs -> their
                   KiCad netclass if it is not Default/Power, else `diff`. Zdiff from a netclass named
                   like `90ohm`, 90 ohm for USB names, else 100 ohm. Edge rates are the defaults below.
  clock nets       single-ended nets with CLK/CK in the name -> class `clk`
  fields           BGAs (footprint name, or >= min_balls pads) and connectors carrying a pair
  rails            power nets (KiCad netclass Power/PWR, or names like +3V3, 1V8, VCC*, VDD*, VIN*,
                   VBUS*, VADJ*) with a voltage read from the name. Source: a shunt (resistor of
                   1 ohm or less to another power net), else a regulator (part value like TPS..,
                   LTC.., AP.., ...), else a ferrite or inductor, else a connector. Loads: the ICs
                   with the most pins on the net. The current is not on the board: rails get
                   `current_a = DEFAULT_CURRENT_A` and `current_assumed = true`, and a judgement made
                   at an assumed current is never worse than `warn` (see evaluate.py).
"""
import fnmatch
import os
import re
from collections import defaultdict

from .pdn import parse_value

DEFAULT_CURRENT_A = 0.5
MIN_BALLS = 196

CLASS_DEFAULTS = {
    "serdes": {"rise_ps": 35.0, "nyquist_ghz": 6.25, "z_diff": 100.0, "z_tol": 10.0,
               "note": "serial-lane names; 12.5 Gb/s edge rate assumed"},
    "diff": {"rise_ps": 300.0, "nyquist_ghz": 1.0, "z_diff": 100.0, "z_tol": 10.0,
             "note": "other differential pairs; LVDS-like edge rate assumed"},
    "clk": {"rise_ps": 500.0, "nyquist_ghz": 0.7, "note": "single-ended clock names; 500 ps edge assumed"},
}

GND_RE = re.compile(r"^(?:.*/)?(?:[ADPSC]?GND\w*|GND_?\w*|VSS\w*|0V)$", re.I)
SERDES_RE = re.compile(r"(?:^|[/_.\-])(?:DP\d|GT[HXYP]?|MGT|PCIE|SATA|USB3|SS(?:TX|RX)|HDMI|TMDS|DPHY|CSI|DSI|"
                       r"REFCLK|GBTCLK|TX\d*|RX\d*|TRD\d)", re.I)
USB_RE = re.compile(r"USB|(?:^|[/_.])D[_.]?[PN]$", re.I)
CLK_RE = re.compile(r"CLK|(?:^|[/_.])CK(?:$|[/_.])", re.I)
CONTROL_RE = re.compile(r"(?:^|[/_.\-])(?:n?EN|PG|POK|PGOOD|PWR_?(?:EN|BUT|ON|OK|GOOD)|LED|SNS|SENSE|FB|SD_PWR|"
                        r"ON|OK|BUT|nLED|ADJ_SEL|SEL|DET|ALERT|INT)(?:$|[/_.\-\d])|LED", re.I)
POWER_RE = re.compile(r"^(?:.*/)?(?:[+\-]\d|\d+V\d*|\d+P\d+V|V(?:CC|DD|IN|BUS|BAT|ADJ|TT|DDQ|PP)|"
                      r"\w*VCC|\w*VDD|\w*AVTT|\+\w+|\w*_\d+V\d*|\w*\d+V\d+\w*|\w*BATT?)", re.I)
REG_RE = re.compile(r"^(TPSM?|LTC|LTM|LT\d|LM\d|ADP|MPM?\d|RT\d|TLV|NCP|NCV|AP\d|SY\d|XC6|MIC|ISL|MAX\d|LP\d|"
                    r"BD\d|IR\d|AMS1|AZ1|LD\d|RP\d|TPS)", re.I)
PASSIVE = ("R", "C", "L", "FB", "D", "TP", "H", "MH", "JP", "SW", "Y", "X", "F", "Q", "LED", "NT", "G", "FID", "SJ")


def _prefix(ref):
    m = re.match(r"[A-Za-z_]+", ref)
    return m.group(0).upper() if m else ""


def _glob_escape(s):
    return re.sub(r"([\[\]*?])", r"[\1]", s)


def voltage_from_name(name):
    """1V8 -> 1.8, +3V3 -> 3.3, 12P0V -> 12.0, +5V -> 5.0, VCC0V85 -> 0.85; None if not in the name."""
    base = name.split("/")[-1]
    m = re.search(r"(\d+)P(\d+)V", base, re.I)
    if m:
        return float(f"{m.group(1)}.{m.group(2)}")
    m = re.search(r"(\d+)V(\d+)", base, re.I)
    if m:
        return float(f"{m.group(1)}.{m.group(2)}")
    m = re.search(r"(?:^|[+_\-])(\d+(?:\.\d+)?)V(?:$|[_\-])", base, re.I)
    if m:
        return float(m.group(1))
    return None


class _Board:
    """Light view of a snapshot for detection (no rasters)."""

    def __init__(self, snap):
        self.snap = snap
        self.fps = {f["ref"]: f for f in snap["footprints"]}
        self.pads_by_ref = defaultdict(list)
        self.pads_by_net = defaultdict(list)
        for p in snap["pads"]:
            self.pads_by_ref[p["ref"]].append(p)
            if p["net"]:
                self.pads_by_net[p["net"]].append(p)
        self.netclass = {n: v.get("netclass", "Default") for n, v in snap["nets"].items()}
        self.zone_nets = defaultdict(int)
        for z in snap["zones"]:
            if z.get("net"):
                self.zone_nets[z["net"]] += 1

    def parts_on(self, net):
        return sorted({p["ref"] for p in self.pads_by_net.get(net, []) if not self.fps.get(p["ref"], {}).get("dnp")})

    def nets_of(self, ref):
        return {p["net"] for p in self.pads_by_ref.get(ref, []) if p["net"]}


def _partner(net, names, suffixes):
    for pat, rep in suffixes:
        if re.search(pat, net):
            cand = re.sub(pat, rep, net)
            if cand in names and cand != net:
                return cand
    return None


def ground_nets(B):
    vias = defaultdict(int)
    for v in B.snap["vias"]:
        vias[v["net"]] += 1
    # a ground net may be zones and vias only (no pads): count every kind of copper
    g = [n for n in B.snap["nets"] if GND_RE.match(n)
         and len(B.pads_by_net.get(n, [])) + 2 * B.zone_nets.get(n, 0) + vias[n] >= 2]
    return sorted(g, key=lambda n: (-B.zone_nets.get(n, 0), -len(B.pads_by_net.get(n, []))))


def is_power(B, n, gnd):
    if n in gnd or n.startswith(("unconnected-", "Net-(")) and B.netclass.get(n, "").upper() not in ("POWER", "PWR"):
        return False
    if B.netclass.get(n, "").upper() in ("POWER", "PWR"):
        return True
    base = n.split("/")[-1]
    if CONTROL_RE.search(base) or base.upper().startswith("VREF"):
        return False
    return bool(POWER_RE.match(n)) and len(B.pads_by_net.get(n, [])) >= 2


TIER = {}          # class -> edge-rate tier ("serdes" or "diff") of the nets in it


def pair_classes(B, suffixes, gnd, power):
    names = set(B.snap["nets"])
    classes = defaultdict(list)
    seen = set()
    zdiff = {}
    for n in sorted(names):
        if n in seen or n in gnd or n in power or n.startswith("unconnected-"):
            continue
        q = _partner(n, names, suffixes)
        if not q:
            continue
        seen |= {n, q}
        nc = B.netclass.get(n, "Default")
        zm = re.search(r"(\d+(?:\.\d+)?)\s*(?:ohm|r|Ω)", nc, re.I)
        serdes = bool(SERDES_RE.search(n.split("/")[-1]) or SERDES_RE.search(n))
        if nc not in ("Default", "") and nc.upper() not in ("POWER", "PWR"):
            cname = re.sub(r"\W+", "_", nc)            # the designer's own netclass comes first (#2)
        elif serdes:
            cname = "serdes" + (f"_{zm.group(1)}ohm" if zm else "")
        else:
            cname = "diff"
        TIER[cname] = "serdes" if (serdes or TIER.get(cname) == "serdes") else "diff"
        classes[cname] += [n, q]
        if zm:                                         # the designer's figure, in the netclass name
            zdiff.setdefault(cname, (float(zm.group(1)), "netclass"))
        elif USB_RE.search(n) and not SERDES_RE.search(n.split("/")[-1]):
            zdiff.setdefault(cname, (90.0, "name"))       # a guess from the net name
    return classes, zdiff, seen


def fields(B, pair_nets):
    out = {}
    for ref, f in B.fps.items():
        pads = B.pads_by_ref.get(ref, [])
        name = (f.get("fpid") or "").split(":")[-1]
        if re.search(r"BGA", name, re.I) or len(pads) >= MIN_BALLS:
            out[f"bga_{ref}"] = {"footprints": [ref], "margin": 0.0, "kind": "bga"}
        elif _prefix(ref) in ("J", "P", "CN", "CON", "X", "MODULE") and len(pads) >= 20 and B.nets_of(ref) & pair_nets:
            out[f"conn_{ref}"] = {"footprints": [ref], "margin": 0.0, "kind": "connector"}
    return out


def _pad_area(p):
    return (p.get("sx") or 0) * (p.get("sy") or 0)


OUT_FUNC = re.compile(r"^(?:V?OUT\w*|VO\d*|SW\d*|LX\d*|PH|VSW)$", re.I)
IN_FUNC = re.compile(r"^(?:V?IN\w*|PVIN\w*)$", re.I)
SWITCH_FUNC = re.compile(r"^(?:SW\d*|LX\d*|PH|VSW)$", re.I)
MON_RE = re.compile(r"^(INA\d|PAC19|LTC29\d\d|LTC41\d\d|MAX9\d{3}|MAX4\d{3}|ADM1\d{3}|ACS7|LM5056|ADS7|TMP\d)", re.I)
CONN = ("J", "P", "CN", "CON", "USB", "X")


def is_regulator(B, ref, gset):
    f = B.fps[ref]
    if f.get("dnp") or _prefix(ref) not in ("U", "IC", "PS", "VR") or MON_RE.search(f.get("value", "")):
        return False
    funcs = {(p.get("pinfunc") or "").upper() for p in B.pads_by_ref[ref]}
    by_func = any(OUT_FUNC.match(x) for x in funcs) and any(IN_FUNC.match(x) for x in funcs)
    return (by_func or bool(REG_RE.search(f.get("value", "")))) and bool(B.nets_of(ref) & gset)


def regulator_outputs(B, ref, power, gset, regs):
    """Nets a regulator drives: from pin functions when the footprint has them, else the power net
    of the part that no other regulator shares (inputs and bias supplies usually are shared) and
    that is not its highest-voltage net."""
    outs = {p["net"] for p in B.pads_by_ref[ref] if p["net"] and OUT_FUNC.match(p.get("pinfunc") or "")}
    if outs:
        return outs - gset
    cands = (B.nets_of(ref) & power) - gset
    shared = {n for n in cands if any(r != ref and n in B.nets_of(r) for r in regs)}
    if cands - shared:
        cands = cands - shared
    known = {n: voltage_from_name(n) for n in cands if voltage_from_name(n) is not None}
    if len(cands) > 1 and known:
        cands = cands - {max(known, key=known.get)}
    return cands


def _series_parts(B, supply, gset, max_ohm=1.0):
    """Series parts between supply nets: shunts (R <= 1 ohm), ferrites, inductors, fuses, and pass
    FETs (Q with two supply nets and a gate)."""
    out = []
    for ref, f in B.fps.items():
        nets = B.nets_of(ref)
        px = _prefix(ref)
        if f.get("dnp") or nets & gset:
            continue
        if px == "Q" and len(nets) == 3 and len(nets & supply) == 2:
            a, b = sorted(nets & supply)
            out.append((ref, "switch", a, b))
            continue
        if len(nets) != 2 or not nets <= supply:
            continue
        if px == "R":
            v = parse_value(f.get("value", ""))
            if v is None or v > max_ohm:
                continue
            kind = "shunt"
        elif px in ("FB", "L"):
            kind = "ferrite" if px == "FB" else "inductor"
        elif px == "F":
            kind = "fuse"
        else:
            continue
        a, b = sorted(nets)
        out.append((ref, kind, a, b))
    return out


def entry_connectors(B, power, gset, reg_out):
    """Nets where power enters the board: a connector whose other pads are all supply or ground
    (jack, terminal, power header), or a USB receptacle's VBUS that no regulator drives."""
    out = {}
    for ref in B.fps:
        if _prefix(ref) not in CONN or B.fps[ref].get("dnp"):
            continue
        nets = B.nets_of(ref) - gset
        pw = nets & power
        if not pw:
            continue
        name = (B.fps[ref].get("value", "") + " " + B.fps[ref].get("fpid", "")).upper()
        sig = {n for n in nets - pw if not n.startswith(("unconnected-", "Net-("))}
        if not sig or ("USB" in name or "TYPE-C" in name or "TYPEC" in name):
            for n in pw:
                if n not in reg_out:
                    out.setdefault(n, ref)
    return out


def _has_conn(B, net):
    return any(_prefix(r) in CONN for r in B.parts_on(net))


def rails(B, power, gnd):
    """Each power net with a voltage in its name -> rail dict; the others are listed apart.

    The supply tree: regulator outputs (and their switch nodes) are the roots, series parts
    (shunts, ferrites, inductors, fuses) the edges. A net's source is the series part towards the
    nearest root, else its regulator, else a connector (assumed power entry)."""
    gset = set(gnd)
    regs = [r for r in B.fps if is_regulator(B, r, gset)]
    reg_out, switch = {}, set()
    for r in regs:
        for n in regulator_outputs(B, r, power, gset, regs):
            is_sw = any(SWITCH_FUNC.match(p.get("pinfunc") or "") for p in B.pads_by_ref[r] if p["net"] == n)
            if n in power or is_sw:
                reg_out.setdefault(n, r)
                if n not in power:
                    switch.add(n)
    edges = _series_parts(B, power | switch, gset)
    adj = defaultdict(list)
    for ref, kind, a, b in edges:
        adj[a].append((b, ref, kind))
        adj[b].append((a, ref, kind))
    dist, parent = {}, {}

    def walk(roots):
        frontier = [n for n in roots if n not in dist]
        for n in frontier:
            dist[n] = 0
        while frontier:
            nxt = []
            for n in frontier:
                for m, ref, kind in adj[n]:
                    if m not in dist:
                        dist[m], parent[m] = dist[n] + 1, (ref, kind)
                        nxt.append(m)
            frontier = nxt

    walk(reg_out)
    entries = entry_connectors(B, power, gset, reg_out)
    walk(entries)
    # still unreached: the series part's supply side is the side with fewer parts
    for ref, kind, a, b in edges:
        for n, m in ((a, b), (b, a)):
            if n not in dist and n not in parent and len(B.parts_on(m)) < len(B.parts_on(n)):
                parent[n] = (ref, kind)
    # a net fed through series parts carries its supply's voltage
    volt = {}

    def v_of(n, depth=0):
        if n in volt:
            return volt[n]
        v = voltage_from_name(n)
        if v is None and n in parent and depth < 20:
            ref = parent[n][0]
            up = next((x for r, k, a, b in edges if r == ref for x in (a, b) if x != n), None)
            v = v_of(up, depth + 1) if up else None
        volt[n] = v
        return v

    out, left = {}, []
    for net in sorted(power):
        if net in parent:
            src, kind = parent[net]
        elif net in reg_out:
            src, kind = reg_out[net], "regulator"
        elif net in entries:
            src, kind = entries[net], "connector"
        else:
            src = next((r for r in B.parts_on(net) if _prefix(r) in CONN), None)
            kind = "connector" if src else None
        loads = _loads(B, net, src, edges, parent)
        entry = {"net": net, "source": src, "source_kind": kind, "loads": loads}
        v = v_of(net)
        if v is None or not src or not loads:
            entry["skip"] = ("no voltage in the net name" if v is None else
                             "no source found" if not src else "no load found")
            left.append(entry)
            continue
        key = re.sub(r"\W+", "_", net.split("/")[-1]).strip("_") or net
        out[key] = {"net": net, "v": v, "tol_pct": 3.0 if v <= 1.2 else 5.0, "current_a": DEFAULT_CURRENT_A,
                    "current_assumed": True, "source": src, "source_kind": kind, "loads": loads}
        # #4: currents from the loads' footprint fields (Current = 1.5A / 1500mA) are sourced, not assumed
        cur, srcs = 0.0, []
        for ref in loads:
            fv = B.fps.get(ref, {}).get("fields", {})
            for k in CURRENT_FIELDS:
                if k in fv and current_value(fv[k]) is not None:
                    cur += current_value(fv[k])
                    srcs.append(f"{ref}.{k}")
                    break
        if srcs:
            out[key].update(current_a=round(cur, 6), current_assumed=False, current_source=", ".join(srcs))
        if voltage_from_name(net) is None:
            out[key]["v_from"] = "supply side of " + src
    return out, left


CURRENT_FIELDS = ("Current", "current", "I", "Iload", "I_load", "Current_A", "Imax")


def current_value(s):
    """'1.5A', '1500mA', '1.5', '2 A' -> amperes; None if not a current."""
    m = re.fullmatch(r"\s*([\d.]+)\s*(mA|uA|A)?\s*", str(s), re.I)
    if not m:
        return None
    v = float(m.group(1))
    u = (m.group(2) or "A").lower()
    return v * {"a": 1.0, "ma": 1e-3, "ua": 1e-6}[u]


def _loads(B, net, src, edges, parent):
    """The parts that draw the rail's current: ICs and connectors with the most pins on the net
    (current monitors and the source excluded); for a net that only feeds other rails, the series
    parts towards them."""
    count = {}
    for ref in B.parts_on(net):
        if ref == src or _prefix(ref) in PASSIVE or MON_RE.search(B.fps[ref].get("value", "")):
            continue
        count[ref] = sum(1 for p in B.pads_by_ref[ref] if p["net"] == net)
    if count:
        top = max(count.values())
        return sorted(r for r, c in count.items() if c >= max(1, top // 4))
    return sorted(ref for ref, kind, a, b in edges
                  if net in (a, b) and ref != src and parent.get(b if a == net else a, (None,))[0] == ref)


def detect(snap, suffixes=None):
    """-> (config dict, notes dict). notes lists what was assumed and what was left out."""
    from .config import DEFAULTS
    suffixes = suffixes or DEFAULTS["pairs"]["suffixes"]
    B = _Board(snap)
    TIER.clear()
    gnd = ground_nets(B)
    power = {n for n in B.snap["nets"] if is_power(B, n, gnd)}
    classes, zdiff, pair_nets = pair_classes(B, suffixes, gnd, power)
    cfg = {"board": {"reference_nets": gnd[:3] or ["GND"]}, "classes": {}, "rules": {}, "rails": {}, "areas": {}}
    notes = {"detected": True, "assumptions": [], "left_out": []}
    for cname, nets in sorted(classes.items()):
        d = dict(CLASS_DEFAULTS["serdes" if TIER.get(cname) == "serdes" or cname.startswith("serdes") else
                               cname if cname in CLASS_DEFAULTS else "diff"])
        note = d.pop("note")
        from . import standards
        std = standards.of_net(nets[0])
        z, zsrc = zdiff.get(cname, (None, None))
        if std:                                # #6: the interface standard of the class's nets
            d.update({k: std[1][k] for k in ("z_diff", "z_tol", "rise_ps", "nyquist_ghz")})
            note = f"standard {std[0]}: {std[1]['source']}"
        if z is not None and (zsrc == "netclass" or not std):
            if std and z != d["z_diff"]:       # the netclass's impedance, at the standard's relative tolerance
                d["z_tol"] = round(d["z_tol"] * z / d["z_diff"], 3)
                note += f"; Zdiff {z:g} ohm from the netclass name"
            d["z_diff"] = z
        d["patterns"] = sorted(_glob_escape(n) for n in nets)
        cfg["classes"][cname] = d
        notes["assumptions"].append(f"class {cname}: {len(nets) // 2} pairs, rise {d['rise_ps']:g} ps, "
                                    f"Nyquist {d['nyquist_ghz']:g} GHz, Zdiff {d['z_diff']:g} ohm ({note})")
    clk = [n for n in B.snap["nets"] if n not in pair_nets and n not in power and n not in gnd
           and not n.startswith("unconnected-") and CLK_RE.search(n.split("/")[-1])
           and len(B.pads_by_net.get(n, [])) >= 2]
    if clk:
        d = dict(CLASS_DEFAULTS["clk"])
        note = d.pop("note")
        d["patterns"] = sorted(_glob_escape(n) for n in clk)
        cfg["classes"]["clk"] = d
        notes["assumptions"].append(f"class clk: {len(clk)} nets, rise {d['rise_ps']:g} ps ({note})")
    sig = [c for c in cfg["classes"]]
    fast = [c for c in sig if c != "clk"]
    for rule in ("ref_plane", "split_crossing", "return_via", "via_stub", "edge_distance", "pair_coupling",
                 "parallel", "pad_plane", "stubs"):
        cfg["rules"][rule] = {"classes": sig if rule in ("ref_plane", "split_crossing", "return_via", "stubs") else fast}
    cfg["rules"]["stitching"] = {"classes": [c for c in fast if TIER.get(c) == "serdes" or c.startswith("serdes")] or fast}
    cfg["crosstalk"] = {"classes": sig}
    ar = fields(B, pair_nets)
    cfg["areas"] = {k: {"footprints": v["footprints"], "margin": v["margin"]} for k, v in ar.items()}
    bgas = [v["footprints"][0] for v in ar.values() if v["kind"] == "bga"]
    big_ics = sorted((r for r in B.fps if _prefix(r) in ("U", "IC") and len(B.pads_by_ref[r]) >= 32),
                     key=lambda r: -len(B.pads_by_ref[r]))
    cfg["rules"]["decap"] = {"load_refs": bgas or big_ics[:3]}
    rl, no_v = rails(B, power, gnd)
    cfg["rails"] = {k: {kk: vv for kk, vv in v.items() if kk != "source_kind"} for k, v in rl.items()}
    regs = sorted(r for r in B.fps if is_regulator(B, r, set(gnd)))
    if regs:
        cfg["rules"]["thermal_vias"] = {"parts": regs}
    cfg["engines"] = {"dc": "openpcbsi", "pdn": "openpcbsi", "cross_check": False}
    for k, v in rl.items():
        cur = (f"current {v['current_a']:g} A from {v['current_source']}" if not v["current_assumed"]
               else f"current assumed {DEFAULT_CURRENT_A} A")
        notes["assumptions"].append(f"rail {k}: {v['v']:g} V +/-{v['tol_pct']:g} % from {v['source_kind']} {v['source']} "
                                    f"to {', '.join(v['loads'])}, {cur}")
    for e in no_v:
        why = e.get("skip") or "no voltage in the net name"
        notes["left_out"].append(f"power net {e['net']}: {why} (source {e['source'] or '-'}, loads {', '.join(e['loads']) or '-'})")
    if not any("ripple_pct" in r for r in cfg["rails"].values()):
        notes["left_out"].append("PDN impedance (G11): needs ripple_pct and step_a per rail, not on the board")
    notes["left_out"].append("thermal dissipation per part (G13): only I2R of the assumed rail currents")
    notes["left_out"].append("G7 channels and G3 drivers/terminations: none configured")
    notes["fields"] = {k: v["kind"] for k, v in ar.items()}
    notes["reference_nets"] = cfg["board"]["reference_nets"]
    return cfg, notes


def matches(cfg, net):
    for cname, c in cfg.get("classes", {}).items():
        if any(fnmatch.fnmatchcase(net, p) for p in c.get("patterns", [])):
            return cname
    return None


# ------------------------------------------------------------------------------------------------
def _tv(v):
    """TOML value."""
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, (int, float)):
        return repr(float(v)) if isinstance(v, float) else str(v)
    if isinstance(v, str):
        return '"' + v.replace("\\", "\\\\").replace('"', '\\"') + '"'
    if isinstance(v, list):
        items = [_tv(x) for x in v]
        one = "[" + ", ".join(items) + "]"
        return one if len(one) <= 96 else "[\n" + "".join(f"    {x},\n" for x in items) + "]"
    if isinstance(v, dict):
        return "{ " + ", ".join(f"{_tk(k)} = {_tv(x)}" for k, x in v.items()) + " }"
    raise TypeError(v)


def _tk(k):
    return k if re.fullmatch(r"[A-Za-z0-9_\-]+", k) else _tv(k)


def to_toml(cfg, notes, board_path):
    """The detected config as a TOML file to review: every assumption is commented."""
    L = [f"# openpcbsi board config drafted by `python3 -m openpcbsi init` from",
         f"#   {os.path.basename(board_path)}   (file name only: a config in boards/ is shipped, #63)",
         "# Everything below was detected from the board (detect.py). Review it before trusting a",
         "# result: edge rates, rail currents, tolerances and ripple targets are not on the board.",
         "# The pilot's boards/au7p-fmc.toml shows every section with comments.", ""]
    L += ["# Assumed:"] + [f"#   {a}" for a in notes["assumptions"]]
    L += ["# Left out (add by hand where needed):"] + [f"#   {a}" for a in notes["left_out"]] + [""]
    L += ["[board]", f"reference_nets = {_tv(cfg['board']['reference_nets'])}", ""]
    L.append("# --- signal classes: edge rates are defaults, set them from the interface standards ---------")
    for cname, c in cfg["classes"].items():
        L.append(f"[classes.{_tk(cname)}]")
        for k in ("rise_ps", "nyquist_ghz", "z_diff", "z_tol"):
            if k in c:
                L.append(f"{k} = {_tv(c[k])}")
        L.append(f"patterns = {_tv(c['patterns'])}")
        L.append("")
    L.append("# --- exclusion areas (BGA and connector fields) ----------------------------------------------")
    for k, a in cfg["areas"].items():
        L += [f"[areas.{_tk(k)}]", f"footprints = {_tv(a['footprints'])}", f"margin = {_tv(a['margin'])}", ""]
    L.append("# --- rules: which classes each G1 rule applies to (limits stay at their defaults) ------------")
    for r, d in cfg["rules"].items():
        L.append(f"[rules.{r}]")
        for k, v in d.items():
            L.append(f"{k} = {_tv(v)}")
        L.append("")
    L.append("# --- rails: CURRENTS ARE ASSUMED. Set current_a (and ripple_pct / step_a for G11) per rail; --")
    L.append("# --- then delete current_assumed. While it is true a DC drop is never judged worse than warn. -")
    for k, r in cfg["rails"].items():
        L.append(f"[rails.{_tk(k)}]")
        for kk in ("net", "v", "tol_pct", "current_a", "current_assumed", "source", "loads"):
            if kk in r:
                L.append(f"{kk} = {_tv(r[kk])}")
        if r.get("v_from"):
            L.append(f"# v taken from the {r['v_from']}")
        L.append("# ripple_pct = 3.0   # PDN target: V * ripple_pct / step_a")
        L.append("# step_a = 0.5")
        L.append("")
    L += ["[engines]", 'dc = "openpcbsi"     # "kipi" needs [kipi] and [kipi.input] (see boards/au7p-fmc.toml)',
          'pdn = "openpcbsi"', "cross_check = false", "",
          "[crosstalk]", f"classes = {_tv(cfg['crosstalk']['classes'])}", "",
          "[evaluate]", "warnings_exit = true      # exit code 1 when there are warnings and no fails",
          'drc_warning = "warn"      # KiCad DRC warnings as "warn" or "fail"', "dc_warn_pct = 50.0", ""]
    return "\n".join(L)
