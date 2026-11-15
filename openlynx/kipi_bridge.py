"""KiPI as the G10 (DC drop) and G11/G12 (PDN impedance) engine.

KiPI lives in deps/ (pinned in deps/kipi.lock, set up by deps/setup_kipi.sh) and runs in its own
environment; openlynx never imports it. This module:

  1. derives KiPI's power tree from the openlynx rail config and the board: each rail's shunt
     (`source`) is a series link, and the walk continues upstream through low-value series parts
     (resistors <= 1 ohm, ferrites, inductors) to the regulator that drives the net (`regulator`
     in the rail, or the U part with the most pad area on that net),
  2. writes the setup as JSON (KiPI reads JSON setups), plus a job file,
  3. runs deps/kipi_driver.py with deps/kipi-venv/bin/python, and
  4. turns KiPI's results into openlynx's dc / pdn rows (same keys as pdn.py, plus KiPI's own).

Differences from openlynx's own solver, by design: KiPI solves the ground net as well (loads and
regulators have ground pins), and its meshes are FEM triangles rather than raster cells.
"""
import json
import math
import os
import re
import subprocess
from collections import defaultdict

import numpy as np

from .pdn import esr_default, package_of, parse_value

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEPS = os.path.join(ROOT, "deps")
VENV_PY = os.path.join(DEPS, "kipi-venv", "bin", "python")
DRIVER = os.path.join(DEPS, "kipi_driver.py")

LDO_PAT = re.compile(r"TPS7|TLV7|LP5|LT30|ADP1|NCP1|LDO|AP2", re.I)
REG_PAT = re.compile(r"^(TPSM?|LTC|LTM|LT\d|LM\d|ADP|MPM?\d|RT\d|TLV|NCP|NCV|AP\d|SY\d|XC6|MIC|ISL|MAX\d|LP\d|BD\d|IR\d)", re.I)


def available():
    return os.path.exists(VENV_PY) and os.path.exists(DRIVER) and os.path.isdir(os.path.join(DEPS, "KiPI", "kipi"))


def lock():
    p = os.path.join(DEPS, "kipi.lock")
    out = {}
    if os.path.exists(p):
        for line in open(p):
            m = re.match(r'\s*(\w+)\s*=\s*"([^"]*)"', line)
            if m:
                out[m.group(1)] = m.group(2)
    return out


# ----------------------------------------------------------------------------------------------
# power tree

def _pin(ref, net, pads):
    return {"ref": ref, "net": net, "pins": sorted({p["num"] for p in pads if p["net"] == net and p["num"]})}


def _pad_area(p):
    return p["sx"] * p["sy"]


class Tree:
    def __init__(self, ctx):
        self.ctx, self.b, self.cfg = ctx, ctx.b, ctx.cfg
        self.kc = self.cfg.data.get("kipi", {})
        self.gnd = self.cfg.reference_nets[0]
        self.vrms = {}          # output net -> vrm dict
        self.links = {}         # ref -> link dict
        self.loads = []
        self.probes = []
        self.nets = {self.gnd}
        self.notes = []

    def parts_on(self, net):
        return {p["ref"] for p in self.b.pads if p["net"] == net}

    def series_value(self, ref):
        """DC resistance of a 2-terminal series part, or None if it is not one."""
        f = self.b.fps[ref]
        pads = [p for p in self.b.pads_by_ref[ref] if p["net"]]
        nets = {p["net"] for p in pads}
        if f["dnp"] or len(nets) != 2 or self.gnd in nets:
            return None
        if ref.startswith(("FB", "L")):
            return self.kc.get("ferrite_dcr_ohm", 0.1)
        if ref.startswith("R"):
            v = parse_value(f["value"])
            if v is not None and v <= self.kc.get("max_series_ohm", 1.0):
                return max(v, 1e-6)
        return None

    def add_link(self, ref):
        if ref in self.links:
            return
        pads = self.b.pads_by_ref[ref]
        a, c = sorted({p["net"] for p in pads if p["net"]})
        self.links[ref] = {"name": ref, "side_a": _pin(ref, a, pads), "side_b": _pin(ref, c, pads),
                           "resistance": self.series_value(ref)}
        self.nets |= {a, c}

    def regulator_on(self, net, hint=None):
        """The regulator driving `net`: the configured one if it has a pad on the net, else a part
        whose value looks like a regulator (never a load of any rail), most pad area on the net."""
        if hint and hint in self.b.fps and any(p["net"] == net for p in self.b.pads_by_ref[hint]):
            return hint
        load_refs = {r for rail in self.cfg["rails"].values() for r in rail.get("loads", [])}
        best = None
        for ref in self.parts_on(net):
            f = self.b.fps[ref]
            if not ref.startswith("U") or f["dnp"] or ref in load_refs or not REG_PAT.search(f["value"]):
                continue
            pads = self.b.pads_by_ref[ref]
            if not any(p["net"] == self.gnd for p in pads):
                continue
            a = sum(_pad_area(p) for p in pads if p["net"] == net)
            if best is None or a > best[0]:
                best = (a, ref)
        return best[1] if best else None

    def add_rail(self, name, rail):
        b = self.b
        net = rail["net"]
        joined_before = set(self.nets)
        self.nets.add(net)
        cur, used, hops = net, set(), []
        hint = rail.get("regulator")
        for depth in range(8):
            if cur in self.vrms or (depth > 0 and cur in joined_before):
                break                                   # joined the tree built by earlier rails
            if depth == 0 and rail.get("source") and self.series_value(rail["source"]) is not None:
                cand = [rail["source"]]                 # the rail's shunt comes first
            else:
                reg = self.regulator_on(cur, hint)
                if reg:
                    self.add_vrm(cur, reg, rail)
                    break
                cand = sorted(r for r in self.parts_on(cur) if r not in used and self.series_value(r) is not None)
            if not cand:
                self.notes.append(f"{name}: no regulator or series part found upstream of {cur}")
                break
            ref = cand[0]
            used.add(ref)
            self.add_link(ref)
            hops.append(ref)
            cur = ({p["net"] for p in b.pads_by_ref[ref] if p["net"]} - {cur}).pop()
        # the rail's loads
        for lref in rail.get("loads", []):
            pads = b.pads_by_ref.get(lref, [])
            if not any(p["net"] == net for p in pads):
                continue
            tol = rail.get("tol_pct", 3.0)
            self.loads.append({"name": f"{name}:{lref}", "supply": _pin(lref, net, pads), "current": rail["current_a"],
                               "ground": _pin(lref, self.gnd, pads), "tol_min": -tol, "tol_max": tol, "tol_pct": True})
            if "ripple_pct" in rail:
                zt = rail["v"] * rail["ripple_pct"] / 100 / rail.get("step_a", rail["current_a"] / 2)
                self.probes.append({"name": name, "rail": _pin(lref, net, pads), "ground": _pin(lref, self.gnd, pads),
                                    "target_ohm": zt, "target_hz": rail.get("target_fmax_hz",
                                                                            self.cfg["pdn"].get("target_fmax_hz", 1e8)),
                                    "step_a": rail.get("step_a")})
        return hops

    def add_vrm(self, net, reg, rail):
        if net in self.vrms:
            return
        b = self.b
        pads = b.pads_by_ref[reg]
        f = b.fps[reg]
        kind = "ldo" if LDO_PAT.search(f["value"]) else "smps"
        others = {p["net"] for p in pads if p["net"] and p["net"] not in (net, self.gnd)}
        vrm = {"name": reg, "output": _pin(reg, net, pads), "voltage": rail["v"], "kind": kind,
               "ground": _pin(reg, self.gnd, pads), "bandwidth_khz": rail.get("bandwidth_khz",
                                                                              self.kc.get("bandwidth_khz", 50.0))}
        self.vrms[net] = vrm
        self.nets.add(net)
        vrm["_inputs"] = others            # resolved once every rail is in (inputs inside the tree)
        sense = self.sense_pin(rail)
        if sense:
            vrm["sense"] = sense

    def sense_pin(self, rail):
        """The pad on the rail where the regulator senses: a part joining the sense net to the rail,
        or one series part away (FB -> divider -> CORE_SNS -> 0 ohm -> VCCINT)."""
        sn = rail.get("sense")
        if not sn:
            return None
        b, net = self.b, rail["net"]
        nets_of = {r: {p["net"] for p in ps} for r, ps in b.pads_by_ref.items() if not b.fps[r]["dnp"]}
        for r, ns in nets_of.items():
            if sn in ns and net in ns:
                return _pin(r, net, b.pads_by_ref[r])
        for r, ns in nets_of.items():
            if sn in ns and len(ns) == 2:
                other = (ns - {sn}).pop()
                for r2, ns2 in nets_of.items():
                    if r2 != r and other in ns2 and net in ns2:
                        return _pin(r2, net, b.pads_by_ref[r2])
        return None

    def add_input(self):
        """[kipi.input]: the board's supply input as an ideal source at `source` (its pins on the
        source net; its ground pins are the ground reference), joined to `net` (where the
        switchers draw from) through `chain`, the series parts in order from `net` upstream.
        Parts that are not resistors take their resistance from `resistances`."""
        inp = self.kc.get("input")
        if not inp:
            self.notes.append("no [kipi.input]: switchers have no input and ground has no reference")
            return
        b = self.b
        cur = inp["net"]
        chain = inp.get("chain", [])
        rmap = inp.get("resistances", {})
        for k, ref in enumerate(chain):
            pads = b.pads_by_ref[ref]
            nets = {p["net"] for p in pads if p["net"] and p["net"] != self.gnd}
            if cur not in nets:
                raise ValueError(f"[kipi.input] chain: {ref} has no pin on {cur}")
            nxt_nets = ({p["net"] for p in b.pads_by_ref[chain[k + 1]]} if k + 1 < len(chain)
                        else {p["net"] for p in b.pads_by_ref[inp["source"]]})
            up = [n for n in nets - {cur} if n in nxt_nets]
            if not up:
                raise ValueError(f"[kipi.input] chain: cannot find the net between {ref} and the next part")
            r = rmap.get(ref, parse_value(b.fps[ref]["value"]) if ref.startswith("R") else None)
            if r is None:
                raise ValueError(f"[kipi.input] resistances: give a DC resistance for {ref}")
            self.links[ref] = {"name": ref, "side_a": _pin(ref, cur, pads), "side_b": _pin(ref, up[0], pads),
                               "resistance": max(float(r), 1e-6)}
            self.nets |= {cur, up[0]}
            cur = up[0]
        src = inp["source"]
        pads = b.pads_by_ref[src]
        self.vrms[cur] = {"name": src, "output": _pin(src, cur, pads), "voltage": inp["v"], "kind": "source",
                          "ground": _pin(src, self.gnd, pads), "_inputs": set()}
        self.nets.add(cur)

    def caps(self):
        """Capacitor models grouped by value and package, for every capacitor between a tree net and ground."""
        b = self.b
        esl_tab = self.cfg["pdn"]["esl_nh"]
        groups = defaultdict(list)
        for ref, pads in b.pads_by_ref.items():
            if not ref.startswith("C") or b.fps[ref]["dnp"]:
                continue
            nets = {p["net"] for p in pads}
            if self.gnd not in nets or not (nets & self.nets - {self.gnd}):
                continue
            f = b.fps[ref]
            c = parse_value(f["value"])
            if not c:
                continue
            pkg = package_of(f["fpid"])
            groups[(f["value"].split()[0], pkg)].append((ref, c))
        out = []
        for (val, pkg), items in sorted(groups.items()):
            c = items[0][1]
            out.append({"name": f"{val} {pkg}", "capacitance_f": c, "esr_ohm": esr_default(c, pkg),
                        "esl_h": esl_tab.get(pkg, 1.5 if pkg == "elec" else 0.6) * 1e-9,
                        "refs": [r for r, _ in items]})
        return out

    def setup(self):
        vrm_nets = set(self.vrms)
        vrms = []
        for net, v in self.vrms.items():
            ins = [n for n in v.pop("_inputs") if n in vrm_nets or any(n in (l["side_a"]["net"], l["side_b"]["net"])
                                                                       for l in self.links.values())]
            if ins:
                v["input"] = _pin(v["name"], ins[0], self.b.pads_by_ref[v["name"]])
            vrms.append(v)
        dc = self.cfg["dc"]
        kc = self.kc
        pdn = self.cfg["pdn"]
        return {
            "format_version": 2, "vrms": vrms, "loads": self.loads, "links": list(self.links.values()),
            "nets": sorted(self.nets), "caps": self.caps(), "probes": self.probes,
            "settings": {"max_edge_mm": kc.get("max_edge_mm", 0.5),
                         "via_rules": [{"kind": "any", "plating_mm": dc["plating_um"] / 1000.0, "fill": "nonconductive"}],
                         "rho_ohm_mm": 1.0 / dc["sigma"] * 1000.0},
            "ac": {"f_start_hz": pdn["fmin_hz"], "f_stop_hz": pdn["fmax_hz"],
                   "points_per_decade": kc.get("ac_points_per_decade", 25),
                   "field_points_per_decade": kc.get("ac_field_points_per_decade", 8),
                   "max_edge_mm": kc.get("ac_max_edge_mm", 1.0), "ground": self.gnd, "dc_bias": False},
        }


def build(ctx, log=print):
    t = Tree(ctx)
    for name, rail in ctx.cfg["rails"].items():
        hops = t.add_rail(name, rail)
        log(f"  {name:10s} {rail['net']:10s} <- " + " <- ".join(hops + [t._driver_of(rail["net"], hops)]))
    t.add_input()
    inp = ctx.cfg.data.get("kipi", {}).get("input")
    if inp:
        log(f"  input      {inp['net']:10s} <- " + " <- ".join(inp.get("chain", []) + [f"{inp['source']} (source {inp['v']} V)"]))
    for n in t.notes:
        log(f"  note: {n}")
    return t


def _driver_of(self, net, hops):
    # the regulator at the top of a rail's chain, for the log
    cur = net
    for ref in hops:
        cur = ({p["net"] for p in self.b.pads_by_ref[ref] if p["net"]} - {cur}).pop()
    v = self.vrms.get(cur)
    return f"{v['name']} ({v['kind']}{', sense ' + v['sense']['ref'] if v.get('sense') else ''}) on {cur}" if v else f"(joins {cur})"


Tree._driver_of = _driver_of


# ----------------------------------------------------------------------------------------------

def run(ctx, mode, out_dir, cache_dir, log=print, ac_rails=None):
    if not available():
        raise RuntimeError("KiPI is not set up: run deps/setup_kipi.sh")
    log(f"KiPI {lock().get('tag')} ({lock().get('commit', '')[:10]}) power tree:")
    tree = build(ctx, log)
    setup = tree.setup()
    os.makedirs(out_dir, exist_ok=True)
    sp = os.path.join(out_dir, f"{ctx.b.name}.kipi-setup.json")
    with open(sp, "w") as f:
        json.dump(setup, f, indent=1)
    want = [r for r in ctx.cfg["rails"] if "ripple_pct" in ctx.cfg["rails"][r]] if ac_rails is None else ac_rails
    rail_vrm = {}
    for name in want:
        net = ctx.cfg["rails"][name]["net"]
        vrm = _vrm_for(tree, net)
        if vrm:
            rail_vrm[name] = vrm
    job = {"board": ctx.b.path, "board_sha256": ctx.b.sha, "setup": sp, "mode": mode,
           "ac_rails": sorted(set(rail_vrm.values())), "ground": tree.gnd,
           "cache_dir": os.path.join(cache_dir, "kipi"), "out": os.path.join(out_dir, "kipi_result.json")}
    jp = os.path.join(out_dir, "kipi_job.json")
    with open(jp, "w") as f:
        json.dump(job, f, indent=1)
    env = dict(os.environ, KIPI_HOME=os.path.join(DEPS, "kipi-home"), MPLCONFIGDIR=os.path.join(DEPS, "kipi-home", "mpl"),
               OPENBLAS_NUM_THREADS=os.environ.get("OPENBLAS_NUM_THREADS", "4"), QT_QPA_PLATFORM="offscreen")
    p = subprocess.run([VENV_PY, "-I", DRIVER, jp], env=env, capture_output=True, text=True)
    with open(os.path.join(out_dir, "kipi_driver.log"), "w") as f:
        f.write(p.stdout + "\n" + p.stderr)
    if p.returncode != 0 or not os.path.exists(job["out"]):
        raise RuntimeError(f"KiPI driver failed (exit {p.returncode}); see {out_dir}/kipi_driver.log\n{p.stderr[-2000:]}")
    res = json.load(open(job["out"]))
    for line in res["log"]:
        log(f"  {line}")
    for e in res["errors"]:
        log(f"  KiPI error: {e.splitlines()[0]}")
    return res, tree, rail_vrm, sp


def _vrm_for(tree, net):
    """The VRM whose group contains `net` (walking series links upstream)."""
    seen = {net}
    frontier = [net]
    while frontier:
        n = frontier.pop()
        if n in tree.vrms:
            return tree.vrms[n]["name"]
        for l in tree.links.values():
            a, c = l["side_a"]["net"], l["side_b"]["net"]
            for x, y in ((a, c), (c, a)):
                if x == n and y not in seen:
                    seen.add(y)
                    frontier.append(y)
    return None


# ----------------------------------------------------------------------------------------------
# results -> openlynx rows

def dc_rows(ctx, res, tree=None):
    """G10 rows from KiPI. drop_mv: regulator target to the worst ball, supply minus ground (what
    the part sees: shunts, shared upstream copper and the ground return included).
    drop_from_source_mv: from the rail's own shunt (local sense) or the target (remote sense) to the
    worst ball, supply side only, the same definition as openlynx's own solver (the cross-check)."""
    rows = []
    if not res.get("dc"):
        return rows
    loads = {l["name"].split(":")[0]: l for l in res["dc"]["loads"]}
    groups = {g["name"]: g for g in res["dc"]["groups"]}
    for name, rail in ctx.cfg["rails"].items():
        L = loads.get(name)
        if not L:
            rows.append({"rail": name, "net": rail["net"], "engine": "kipi", "problems": ["not in KiPI's results"]})
            continue
        v = rail["v"]
        win = v * rail.get("tol_pct", 3.0) / 100
        v_low = L["v_low"] if L["v_low"] is not None else L["v_effective"]
        drop_eff = v - v_low if v_low is not None else None
        drop_sup = v - L["v_supply_min"] if L["v_supply_min"] is not None else None
        g = groups.get(L["group"], {})
        s = g.get("summary", {})
        # voltage at the rail side of its shunt, from KiPI's series-link results
        v_src = None
        src = rail.get("source")
        if tree is not None and src in tree.links:
            lk = tree.links[src]
            side = "v_a" if lk["side_a"]["net"] == rail["net"] else "v_b"
            for se in g.get("series", []):
                if se["name"] == src:
                    v_src = se[side]
        remote = rail.get("sense") is not None
        if remote:
            # from the sense pin as KiPI solved it (absolute, like the supply pins), not from
            # nominal: the regulator holds sense relative to its own ground, which sits ~1 mV up
            vs = next((x["v_sense"] for x in res["dc"]["vrms"] if x["group"] == L["group"] and x["v_sense"] is not None), None)
            ref = vs if vs is not None else v
            d_src = ref - L["v_supply_min"] if L["v_supply_min"] is not None else None
        else:
            d_src = (v_src - L["v_supply_min"]) if (v_src is not None and L["v_supply_min"] is not None) else None
        gnd = [x for x in res["dc"]["groups"] if x["kind"] == "ground"]
        rows.append({
            "rail": name, "net": rail["net"], "engine": "kipi", "current_a": rail["current_a"], "v_nom": v,
            "drop_mv": round(1000 * drop_eff, 3) if drop_eff is not None else None,
            "drop_pct_of_window": round(100 * drop_eff / win, 2) if drop_eff is not None else None,
            "drop_supply_only_mv": round(1000 * drop_sup, 3) if drop_sup is not None else None,
            "drop_from_source_mv": round(1000 * d_src, 3) if d_src is not None else None,
            "v_ground_mv": round(1000 * L["v_ground"], 3) if L["v_ground"] is not None else None,
            "passed": L["passed"], "status": L["status"], "group": L["group"],
            "max_via_current_a": round(s.get("via_current_max") or 0.0, 4),
            "max_J_A_mm2": round(s.get("j_max") or 0.0, 2), "loss_w": round(s.get("power_loss") or 0.0, 5),
            "mesh_nodes": g.get("mesh_nodes"),
            "top_vias": g.get("top_vias", [])[:5],
            "regulation": "remote sense" if any(v_.get("v_sense") is not None for v_ in res["dc"]["vrms"]
                                                if v_["group"] == L["group"]) else "local",
            "problems": [w for w in res["dc"]["warnings"] if rail["net"] in w][:3]})
    return rows


def pdn_rows(ctx, res, rail_vrm):
    rows = []
    for name, vrm in rail_vrm.items():
        r = res.get("ac", {}).get(vrm)
        if not r:
            rows.append({"rail": name, "engine": "kipi", "problem": "no AC result"})
            continue
        pr = next((p for p in r["probes"] if p["name"] == name), None)
        if pr is None:
            continue
        f = np.array(r["freqs_hz"])
        mag = np.array(pr["z_mag"])
        zt = pr["target_ohm"]
        fmax_t = pr["target_hz"] or f.max()
        band = f <= fmax_t
        viol = (mag > zt) & band
        res_ = []
        for i in range(1, len(f) - 1):
            if mag[i] > mag[i - 1] and mag[i] > mag[i + 1] and mag[i] > 0.3 * zt:
                res_.append({"f_hz": float(f[i]), "Z_ohm": round(float(mag[i]), 5), "x": None, "y": None, "layer": ""})
        rows.append({"rail": name, "net": ctx.cfg["rails"][name]["net"], "engine": "kipi", "vrm": vrm,
                     "target_ohm": round(zt, 5), "target_fmax_hz": fmax_t, "caps": len(r["caps"]),
                     "missing_caps": r["missing_caps"], "freqs_hz": f.tolist(), "Z_ohm": mag.tolist(),
                     "max_Z_ohm": round(float(mag[band].max()), 5), "f_at_max_hz": float(f[band][int(np.argmax(mag[band]))]),
                     "max_ratio": round(float((mag[band] / zt).max()), 3),
                     "violation_band_hz": [float(f[viol].min()), float(f[viol].max())] if viol.any() else None,
                     "plane_capacitance_nf": round((r["plane_capacitance_f"] or 0) * 1e9, 3),
                     "resonances": res_, "solve_s": r["solve_s"], "warnings": r["warnings"][:5]})
    return rows
