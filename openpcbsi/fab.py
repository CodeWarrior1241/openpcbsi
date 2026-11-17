"""MF0: fabricator capability profiles for the manufacturing review (`openpcbsi dfm`).

A profile is a TOML file with three kinds of figure (openpcbsi/fabs/jlcpcb.toml is the exemplar):

  [published]  the fabricator's own limits, from its capability pages;
  [standard]   figures of a published industry standard (IPC), not of this fabricator;
  [assumed]    typical figures that neither states. Only the ranking models use these; no check
               that can fail a board does.

`resolve()` picks the limits that apply to one board (layer count, copper weight, mask colour,
assembly service) and keeps where each came from, so that every finding names its limit's origin.
"""
import hashlib
import os
import tomllib

HERE = os.path.dirname(os.path.abspath(__file__))
OZ_MM = 0.035          # finished copper thickness of 1 oz/ft2


def load(name_or_path):
    """A profile by file path, or by name from openpcbsi/fabs/<name>.toml."""
    path = name_or_path if os.path.exists(name_or_path) else os.path.join(HERE, "fabs", name_or_path + ".toml")
    if not os.path.exists(path):
        have = sorted(f[:-5] for f in os.listdir(os.path.join(HERE, "fabs")) if f.endswith(".toml"))
        raise FileNotFoundError(f"no fabricator profile '{name_or_path}' (a file path, or one of: {', '.join(have)})")
    raw = open(path, "rb").read()
    prof = tomllib.loads(raw.decode())
    prof["_path"], prof["_hash"] = path, hashlib.sha256(raw).hexdigest()[:12]
    return prof


class Fab:
    """Limits for one board: `fab["mask_web_mm"]` is the value, `fab.src("mask_web_mm")` its origin."""

    def __init__(self, prof):
        self.name = prof.get("fab", {}).get("name", "?")
        self.retrieved = prof.get("fab", {}).get("retrieved")
        self.sources = prof.get("fab", {}).get("sources", [])
        self.path, self.hash = prof.get("_path"), prof.get("_hash")
        self.limits = {}       # key -> {"value", "source": published | standard | assumed | config, "note"}
        self.notes = []
        self.used = set()      # keys a check has read (MF0: a limit that nothing reads is reported)

    def set(self, key, value, source, note=""):
        if value is not None:
            self.limits[key] = {"value": value, "source": source, "note": note}

    def __getitem__(self, key):
        self.used.add(key)
        return self.limits[key]["value"]

    def get(self, key, default=None):
        if key in self.limits:
            self.used.add(key)
            return self.limits[key]["value"]
        return default

    def judged(self, key):
        """A limit a rule check may judge with: published, standard or set in the board config."""
        return key in self.limits and self.limits[key]["source"] in ("published", "standard", "config")

    def src(self, key):
        """How a finding names its limit."""
        s = self.limits.get(key, {}).get("source")
        return {"published": f"{self.name} published", "standard": "industry standard",
                "assumed": f"assumed, not published by {self.name}", "config": "board config"}.get(s, "default")

    def describe(self):
        return {"name": self.name, "retrieved": self.retrieved, "sources": self.sources, "profile": self.path,
                "profile_hash": self.hash, "notes": self.notes,
                "limits": [dict(key=k, used=k in self.used, **v) for k, v in self.limits.items()]}


def _row(rows, n, oz):
    """The row of a by-layer-count, by-copper-weight table that applies; (row, note)."""
    rows = [r for r in rows if r["layers"][0] <= n <= r["layers"][1]]
    if not rows:
        return None, None
    fit = sorted((r for r in rows if r["oz"] >= oz - 0.15), key=lambda r: r["oz"])
    if fit:
        return fit[0], None
    row = max(rows, key=lambda r: r["oz"])
    return row, f"outer copper {oz:.1f} oz is heavier than any row for {n} layers; the {row['oz']:g} oz row is used"


def resolve(prof, board, cfg):
    """The limits that apply to this board. `[dfm.fab_limits]` in the board config overrides any key."""
    fab = Fab(prof)
    pub, std, asm = prof.get("published", {}), prof.get("standard", {}), prof.get("assumed", {})
    d = cfg["dfm"]
    n = len(board.layers)
    oz = board.stack.copper[0]["t"] / OZ_MM if board.stack.copper else 1.0
    what = f"{n} layers, {board.stack.copper[0]['t'] * 1000:.0f} um outer copper" if board.stack.copper else f"{n} layers"
    row, note = _row(pub.get("track", []), n, oz)
    if note:
        fab.notes.append("track: " + note)
    if row:
        fab.set("track_width_mm", row["width_mm"], "published", what)
        fab.set("track_space_mm", row["space_mm"], "published", what)
        fab.set("bga_track_mm", row.get("bga_mm"), "published", "width and spacing accepted locally in a BGA fan-out")
    else:
        fab.notes.append(f"the profile has no track row for {n} layers: track width and spacing are not judged")
    e = pub.get("etch", {})
    fab.set("same_net_spacing_mm", e.get("same_net_spacing_mm"), "published", "between tracks of one net")
    fab.set("track_tolerance_pct", e.get("track_tolerance_pct"), "published", "finished width; recorded, not judged")
    dr = pub.get("drill", {})
    for k in ("via_hole_min_mm", "via_diameter_min_mm", "npth_min_mm", "drill_max_mm", "via_hole_spacing_mm", "pad_hole_spacing_mm"):
        fab.set(k, dr.get(k), "published")
    if dr.get("via_diameter_over_hole_mm") is not None:
        fab.set("via_ring_min_mm", dr["via_diameter_over_hole_mm"] / 2, "published",
                f"via diameter at least {dr['via_diameter_over_hole_mm']:g} mm larger than its hole")
    fab.set("plated_slot_min_mm", dr.get("plated_slot_min_2l_mm" if n <= 2 else "plated_slot_min_ml_mm"), "published", what)
    ring, note = _row(pub.get("pth_ring", []), n, oz)
    if note:
        fab.notes.append("annular ring: " + note)
    if ring:
        fab.set("pth_ring_min_mm", ring["min_mm"], "published", what + ", absolute minimum")
        fab.set("pth_ring_rec_mm", ring.get("recommended_mm"), "published", what + ", recommended")
    c = pub.get("clearance", {})
    for k in ("via_hole_to_copper_mm", "pth_hole_to_copper_mm", "npth_hole_to_copper_mm", "inner_via_hole_to_copper_mm",
              "inner_pth_hole_to_copper_mm", "pad_to_track_mm", "bga_pad_to_track_mm", "copper_to_edge_mm"):
        fab.set(k, c.get(k), "published", "from the hole wall" if "hole" in k else "")
    fab.set("copper_to_vcut_mm", c.get("copper_to_vcut_mm"), "published", "recorded only: the board file does not say which edges are V-cut")
    bd = pub.get("board", {})
    fab.set("thickness_min_mm", bd.get("thickness_min_mm"), "published")
    fab.set("thickness_max_mm", bd.get("thickness_max_mm"), "published")
    vc = pub.get("via_cover", {})
    if vc.get("filled_default_from_layers") is not None:
        fab.set("via_cover", "filled" if n >= vc["filled_default_from_layers"] else "tented", "published",
                f"epoxy filled and capped is the default from {vc['filled_default_from_layers']} layers; this board has {n}")
    for k in ("filled_via_min_mm", "filled_via_max_mm", "plugged_via_max_mm", "plugged_to_opening_mm"):
        fab.set(k, vc.get(k), "published")
    m = pub.get("mask", {})
    colour = str(d.get("mask_color", "green")).lower()
    dark = colour in [x.lower() for x in m.get("dark_colors", [])]
    fab.set("mask_web_mm", m.get("web_mm"), "published")
    if oz >= 1.5 and m.get("pad_spacing_2oz_mm") is not None:
        fab.set("mask_pad_spacing_mm", m["pad_spacing_2oz_mm"], "published", "2 oz copper, any colour")
    else:
        fab.set("mask_pad_spacing_mm", m.get("pad_spacing_1oz_dark_mm" if dark else "pad_spacing_1oz_mm"), "published",
                f"1 oz copper, {colour} mask")
    fab.set("mask_expansion_mm", m.get("expansion_mm"), "published", "the expansion the fabricator says it can plot")
    fab.set("mask_to_track_mm", m.get("opening_to_track_mm"), "published")
    p = pub.get("pads", {})
    fab.set("smd_pad_to_pad_mm", p.get("smd_pad_to_pad_mm"), "published")
    fab.set("smd_pad_min_mm", p.get("smd_pad_min_mm"), "published")
    a = pub.get("assembly", {})
    service = str(d.get("assembly_service", "standard")).lower()
    sv = f"{service} assembly service"
    fab.set("reflow_peak_c", a.get(f"reflow_peak_{service}_c"), "published", sv)
    fab.set("reflow_peak_tol_c", a.get("reflow_peak_tol_c"), "published", "recorded only")
    fab.set("min_package", a.get(f"min_package_{service}"), "published", sv)
    fab.set("pin_pitch_mm", a.get(f"pin_pitch_{service}_mm"), "published", sv)
    fab.set("bga_pitch_mm", a.get(f"bga_pitch_{service}_mm"), "published", sv)
    fab.set("double_sided", a.get(f"double_sided_{service}"), "published", sv)
    fab.set("fiducials_required", a.get(f"fiducials_{service}"), "published", sv)
    fab.set("rails_required", a.get(f"rails_{service}"), "published", sv)
    st = pub.get("stencil", {})
    fab.set("stencil_thickness_mm", st.get("thickness_mm"), "published", "standard thicknesses; the fabricator chooses")
    fab.set("stencil_aperture_min_mm", st.get("aperture_min_mm"), "published")
    ps = pub.get("part_spacing", {})
    if ps.get("matrix"):
        fab.set("part_spacing", {"classes": ps["classes"], "matrix": ps["matrix"]}, "published", "recommended minimum")
    fab.set("stencil_area_ratio_min", std.get("stencil", {}).get("area_ratio_min"), "standard", "IPC-7525")
    fab.set("stencil_aspect_ratio_min", std.get("stencil", {}).get("aspect_ratio_min"), "standard", "IPC-7525")
    fab.set("etch_factor", asm.get("etch", {}).get("etch_factor"), "assumed", "recorded only")
    fab.set("etchant", asm.get("etch", {}).get("etchant"), "assumed", "recorded only")
    fab.set("oven", asm.get("assembly", {}).get("oven"), "assumed", "recorded only")
    fab.set("h_w_m2k", asm.get("assembly", {}).get("h_w_m2k"), "assumed", "ranking models only")
    for k, v in prof.get("notes", {}).items():
        fab.notes.append(f"{k}: {v}")
    for k, v in (d.get("fab_limits") or {}).items():
        fab.set(k, v, "config", "set in [dfm.fab_limits] of the board config")
    return fab
