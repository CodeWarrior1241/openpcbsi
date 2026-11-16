"""MF0: fabricator capability profiles for the manufacturing checks (`openpcbsi dfm`).

A profile is a TOML file with the fabricator's published limits under `[published]` and typical
figures the fabricator does not state under `[assumed]` (openpcbsi/fabs/jlcpcb.toml is the exemplar).
`resolve()` picks the limits that apply to one board (layer count, copper weight, mask colour,
assembly service) and keeps where each came from, so a finding can say whether its limit is the
fabricator's published figure, an assumption, or a value set in the board config.
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
        self.limits = {}       # key -> {"value", "source": published | assumed | config, "note"}
        self.notes = []

    def set(self, key, value, source, note=""):
        if value is not None:
            self.limits[key] = {"value": value, "source": source, "note": note}

    def __getitem__(self, key):
        return self.limits[key]["value"]

    def get(self, key, default=None):
        return self.limits[key]["value"] if key in self.limits else default

    def src(self, key):
        """How a finding names its limit: 'JLCPCB published', 'assumed (not published by JLCPCB)', 'board config'."""
        s = self.limits.get(key, {}).get("source")
        return {"published": f"{self.name} published", "assumed": f"assumed, not published by {self.name}",
                "config": "board config"}.get(s, "default")

    def describe(self):
        return {"name": self.name, "retrieved": self.retrieved, "sources": self.sources, "profile": self.path,
                "profile_hash": self.hash, "notes": self.notes,
                "limits": [dict(key=k, **v) for k, v in self.limits.items()]}


def resolve(prof, board, cfg):
    """The limits that apply to this board. `[dfm.fab_limits]` in the board config overrides any key."""
    fab = Fab(prof)
    pub, asm = prof.get("published", {}), prof.get("assumed", {})
    d = cfg["dfm"]
    n = len(board.layers)
    oz = board.stack.copper[0]["t"] / OZ_MM if board.stack.copper else 1.0
    rows = [r for r in pub.get("track", []) if r["layers"][0] <= n <= r["layers"][1]]
    row = None
    if rows:
        fit = sorted((r for r in rows if r["oz"] >= oz - 0.15), key=lambda r: r["oz"])
        row = fit[0] if fit else max(rows, key=lambda r: r["oz"])
        if not fit:
            fab.notes.append(f"outer copper {oz:.1f} oz is heavier than any row of the profile for {n} layers; "
                             f"the {row['oz']:g} oz row is used")
    if row:
        what = f"{n} layers, {row['oz']:g} oz outer copper (board: {board.stack.copper[0]['t'] * 1000:.0f} um)"
        fab.set("track_width_mm", row["width_mm"], "published", what)
        fab.set("track_space_mm", row["space_mm"], "published", what)
        fab.set("bga_track_mm", row.get("bga_mm"), "published", "accepted locally in a BGA fan-out")
    else:
        fab.notes.append(f"the profile has no track row for {n} layers: track width is not judged")
    e = pub.get("etch", {})
    fab.set("same_net_spacing_mm", e.get("same_net_spacing_mm"), "published")
    fab.set("track_tolerance_pct", e.get("track_tolerance_pct"), "published")
    m = pub.get("mask", {})
    colour = str(d.get("mask_color", "green")).lower()
    dark = colour in [c.lower() for c in m.get("dark_colors", [])]
    fab.set("mask_web_mm", m.get("web_mm"), "published")
    if oz >= 1.5 and m.get("pad_spacing_2oz_mm") is not None:
        fab.set("mask_pad_spacing_mm", m["pad_spacing_2oz_mm"], "published", "2 oz copper, any colour")
    else:
        key = "pad_spacing_1oz_dark_mm" if dark else "pad_spacing_1oz_mm"
        fab.set("mask_pad_spacing_mm", m.get(key), "published", f"1 oz copper, {colour} mask")
    fab.set("mask_expansion_mm", m.get("expansion_mm"), "published")
    fab.set("mask_to_track_mm", m.get("opening_to_track_mm"), "published")
    p = pub.get("pads", {})
    fab.set("smd_pad_to_pad_mm", p.get("smd_pad_to_pad_mm"), "published")
    a = pub.get("assembly", {})
    service = str(d.get("assembly_service", "standard")).lower()
    fab.set("reflow_peak_c", a.get(f"reflow_peak_{service}_c"), "published", f"{service} assembly service")
    fab.set("reflow_peak_tol_c", a.get("reflow_peak_tol_c"), "published")
    fab.set("etch_factor", asm.get("etch", {}).get("etch_factor"), "assumed")
    fab.set("etchant", asm.get("etch", {}).get("etchant"), "assumed")
    fab.set("oven", asm.get("assembly", {}).get("oven"), "assumed")
    fab.set("h_w_m2k", asm.get("assembly", {}).get("h_w_m2k"), "assumed")
    for k, v in (d.get("fab_limits") or {}).items():
        fab.set(k, v, "config", "set in [dfm.fab_limits] of the board config")
    return fab
