"""Per-board configuration (TOML). Every key has a default (DEFAULTS below).

A run without a config file does not use DEFAULTS alone: run.load() detects classes, fields and
rails from the board (detect.py, EV5) and `python3 -m openpcbsi init BOARD` writes that detection out
as a starting config. boards/au7p-fmc.toml is the reviewed, commented example.
"""
import copy
import fnmatch
import hashlib
import json
import re
import tomllib

DEFAULTS = {
    "board": {
        "raster_res": 0.025,           # mm, copper raster for plane lookups
        "reference_nets": ["GND"],     # nets that count as return-path copper
        "power_as_reference": False,   # whether power areas count as a reference (PC1)
    },
    # signal classes by edge rate (separate from KiCad net classes)
    "classes": {},
    # exclusion areas: footprints (courtyard + margin) or rectangles [x0, y0, x1, y1]
    "areas": {},
    # override of reference layers per signal layer, e.g. {"In5.Cu": ["In4.Cu", "In6.Cu"]}
    "reference_layers": {},
    "pairs": {
        # regular expressions that turn a P net name into its N partner
        "suffixes": [["_P$", "_N"], ["P$", "N"], ["\\+$", "-"]],
    },
    "rules": {
        "ref_plane": {"enabled": True, "classes": [], "margin_h": 3.0, "step": 0.05, "limit_mm": 0.0},
        "split_crossing": {"enabled": True, "classes": [], "min_gap_mm": 0.0},
        "unreferenced": {"enabled": True, "classes": ["*"], "min_len_mm": 0.5},
        "return_via": {"enabled": True, "classes": [], "limit_mm": 1.0, "count_pads_in_areas": True},
        "via_stub": {"enabled": True, "classes": [], "factor": 3.0},
        "edge_distance": {"enabled": True, "classes": [], "board_edge_mm": 1.0, "plane_edge_mm": 1.0},
        "stitching": {"enabled": True, "classes": [], "max_gap_mm": 2.3, "band_mm": 1.5,
                      "board_edge_gap_mm": 5.0, "edge_band_mm": 2.0},
        "decap": {"enabled": True, "via_mm": 0.3, "prefixes": ["C"], "load_refs": [], "load_mm": 10.0},
        "parallel": {"enabled": True, "classes": [], "max_spacing_mm": 1.0, "limits": [[0.5, 0.0]],
                     "nonpair_min_mm": 0.36},
        "pair_coupling": {"enabled": True, "classes": [], "max_gap_mm": 0.18, "limit_mm": 1.0},
        "pad_plane": {"enabled": True, "classes": [], "ratio": 1.5, "fast_ghz": 3.0},
        "power_overlap": {"enabled": True, "max_dielectric_mm": 0.2, "min_area_mm2": 1.0},
        "stubs": {"enabled": True, "classes": [], "max_stub_mm": 2.0, "dangling_min_mm": 0.05},
        "thermal_vias": {"enabled": True, "parts": []},
        # #23 board-wide stitching map; #27 return path where there is no plane
        "stitching_map": {"enabled": True, "cell_mm": 0.25, "p90_mm": 5.0, "void_mm": 2.0, "void_warn_mm2": 100.0},
        "return_proximity": {"enabled": True, "classes": [], "max_mm": 1.0, "min_len_mm": 1.0},
        "breakout": {"enabled": True, "classes": [], "neck_min_mm": 0.5},       # #24
    },
    "xsec": {
        "res_um": 5.0,
        "mask_er": 3.8,
        "roughness_um": 0.0,          # RMS; 0 = smooth
        "etch_factor": 0.0,           # (top width reduction per side) / thickness
        "freqs_ghz": [0.1, 1.0, 5.0, 10.0],
        "er_ref_ghz": 1.0,            # frequency at which the stackup er / tand are quoted
        "fab": [],                    # [{layer, w, gap, z}] fabricator figures for comparison
        # #7: manufacturing corners for judged geometries: width +/- 2 x etch (gap the other way),
        # dielectric thickness +/- %, er +/- %
        "corners": {"etch_um": 20.0, "thickness_pct": 10.0, "er_pct": 5.0},
    },
    "rails": {},
    "pdn": {"cell_mm": 0.5, "fmin_hz": 1e5, "fmax_hz": 1e9, "points": 120,
            "esl_nh": {"0201": 0.25, "0402": 0.35, "0603": 0.5, "0805": 0.6, "1206": 0.8, "1210": 0.9},
            "candidates": ["100n 0402", "1u 0402", "4.7u 0603", "10u 0603", "47u 0805"]},
    "dc": {"cell_mm": 0.1, "sigma": 5.8e7, "plating_um": 20.0},
    "thermal": {"cell_mm": 1.0, "ambient_c": 25.0, "h_top": 10.0, "h_bot": 10.0, "parts": {}},
    "crosstalk": {"classes": [], "victims": [], "max_spacing_mm": 0.6, "top": 30},
    "si": {"nets": [], "driver": {"r_ohm": 40.0, "rise_ps": 100.0, "v": 1.0}, "load_pf": 2.0},
    "channels": {},
    # manufacturing review (`openpcbsi dfm`, REQUIREMENTS 5.17). A rule check judges geometry against a
    # limit the fabricator publishes (the profile, [dfm] fab), an industry standard, or [dfm.fab_limits].
    # A ranking orders places on the board and judges nothing unless [dfm.judge] gives it a limit.
    "dfm": {
        "fab": "jlcpcb",               # profile name (openpcbsi/fabs/*.toml) or a file path
        "mask_color": "green",
        "assembly_service": "standard",
        "fab_limits": {},              # overrides of resolved limits, e.g. mask_web_mm = 0.12
        "via_cover": "",               # how the vias are ordered: tented, plugged or filled; "" = the profile's default
        "stencil_mm": 0.0,             # stencil thickness; 0 = the thinnest standard stencil of the profile
        "bga_refs": [],                # parts treated as BGA fields besides those whose footprint name says BGA
        "gerbers": "",                 # folder of the Gerber and drill files that will be sent (also --gerbers)
        "workers": 0,                  # worker processes; 0 = a quarter of the CPUs, at most 8
        "acid": {
            "acute_deg": 90.0,            # corners below this opening angle are listed (the drafting rule of thumb)
            "angle_tol_deg": 7.5,         # measurement margin
            "ball_mm": 0.125,             # measuring scale of the pocket search; not a property of a process
            "hires_div": 4,               # the fine raster is raster_res / hires_div
            "tile_mm": 2.0,
            "min_gap_len_mm": 0.3,        # a narrow same-net gap shorter than this is not reported
            "list": 40,                   # corners shown in the report
        },
        "space": {"max_tiles": 600},      # per layer: tiles measured at the fine resolution, closest first
        "chips": {
            "max_package": "1206",        # two-terminal chips up to this size are ranked
            "ring_mm": 0.075,             # where the copper leaving a pad is measured, outside the pad edge
            "via_reach_mm": 1.0,          # a via of the pad's net into a plane, this close, is listed
            "cell_mm": 0.05, "coarse": 4, "margin_mm": 2.0,
            "paste_k": 2.0, "paste_um": 100.0,   # unmelted paste between pad and termination (assumed)
            "max_parts": 0,               # 0 = every chip in scope
            "detail": 8,                  # layouts that get an error bar and a transient run
            "list": 30,
        },
        # Limits for the rankings, from the user's own data. Without them a ranking never warns.
        #   source = "tombstone rate of lot 12, assembler X"   (required: where the limits come from)
        #   joint_dt_s = 0.5        chips whose joints melt further apart warn
        #   reflow_cold_c = 8.0     parts peaking more than this below the board median warn
        #   acute_deg = 60.0, acute_track_mm = 0.15   corners sharper than this on a track this thin warn
        "judge": {},
    },
    "reflow": {
        "cell_mm": 2.0, "dt_s": 1.0,
        "h_top": 0.0, "h_bot": 0.0,       # W/m2K; 0 = the profile's figure (assumed)
        "start_c": 25.0,
        # air temperature of each oven zone (C) and the time the board spends in each (s): a generic
        # lead-free ramp-soak-spike profile, not any assembler's. With `calibrate`, the last
        # `spike_zones` zones are shifted until the median peak over the parts is the target.
        "zones_c": [110.0, 140.0, 160.0, 175.0, 185.0, 200.0, 260.0, 270.0],
        "zone_s": 30.0,
        "cool_c": [120.0, 60.0], "cool_s": 30.0,
        "zone_mm": 350.0,                 # length of a zone: with zone_s it gives the belt speed
        "travel": "x",                    # board axis along the belt: x, -x, y or -y
        "calibrate": True, "spike_zones": 2,
        "peak_target_c": 0.0,             # 0 = the fabricator's stated reflow temperature
        "liquidus_c": 217.0,              # SAC305
        "sensitivity": True,              # repeat with the film coefficient and the assumed bodies varied
        "default_height_mm": 1.0,         # body height of a part with no data
        "list": 12,
        "parts": {},                      # REF = {mass_g = .., height_mm = .., cp = ..}
    },
}


def deep_merge(a, b):
    out = copy.deepcopy(a)
    for k, v in b.items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = deep_merge(out[k], v)
        else:
            out[k] = copy.deepcopy(v)
    return out


class Config:
    def __init__(self, data=None, path=None):
        self.path = path
        self.detected = None
        self.project = None
        self.data = deep_merge(DEFAULTS, data or {})

    @classmethod
    def load(cls, path):
        if not path:
            return cls()
        with open(path, "rb") as f:
            return cls(tomllib.load(f), path)

    def __getitem__(self, k):
        return self.data[k]

    def get(self, *keys, default=None):
        d = self.data
        for k in keys:
            if not isinstance(d, dict) or k not in d:
                return default
            d = d[k]
        return d

    @property
    def raster_res(self):
        return self.data["board"]["raster_res"]

    @property
    def reference_nets(self):
        return self.data["board"]["reference_nets"]

    def hash(self, *sections):
        d = {s: self.data.get(s) for s in sections} if sections else self.data
        return hashlib.sha256(json.dumps(d, sort_keys=True, default=str).encode()).hexdigest()[:12]

    # ---------------------------------------------------------------------------------------
    def class_of(self, net):
        for cname, c in self.data["classes"].items():
            for p in c.get("patterns", []):
                if fnmatch.fnmatchcase(net, p):
                    return cname
        return None

    def class_param(self, cname, key, default=None):
        return self.data["classes"].get(cname, {}).get(key, default)

    def rule_classes(self, rule):
        cl = self.data["rules"][rule].get("classes") or []
        if not cl:
            cl = list(self.data["classes"].keys())
        return cl

    def nets_in_rule(self, board, rule):
        cl = self.rule_classes(rule)
        out = []
        for n in board.net_names[1:]:
            c = self.class_of(n)
            if "*" in cl:
                if n in self.reference_nets or board.kicad_class.get(n) == "PWR":
                    continue
                if n.startswith("unconnected-"):
                    continue
                out.append(n)
            elif c in cl:
                out.append(n)
        return out

    def partner(self, net, names):
        for pat, rep in self.data["pairs"]["suffixes"]:
            if re.search(pat, net):
                cand = re.sub(pat, rep, net)
                if cand in names and cand != net:
                    return cand
        # reverse: N -> P
        for pat, rep in self.data["pairs"]["suffixes"]:
            rpat = re.escape(rep.replace("\\", "")) + "$"
            if re.search(rpat, net):
                p = pat.rstrip("$").replace("\\", "")
                cand = re.sub(rpat, p, net)
                if cand in names and cand != net:
                    return cand
        return None
