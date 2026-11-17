"""MF0: fabricator capability profiles for the manufacturing review (`openpcbsi dfm`).

A profile is a TOML file (openpcbsi/fabs/jlcpcb.toml is the exemplar) with:

  [published]  the fabricator's own limits, from its capability pages. Each value carries the page,
               the page's own words, what it measures, and how certain that reading is:
                 stated       the page says what it measures (a rule may fail a board with it)
                 interpreted  the profile had to choose a reading; a rule may warn, never fail
                 chosen       the board owner chose the reading (named in the note); a rule may fail
                 confirmed    the fabricator answered the question (note says who and when)
               and how it is used: limit (judged), order (the order checklist), input, record.
  [standard]   figures of a published industry standard (IPC), not of this fabricator;
  [guidance]   the fabricator's design articles (not capability pages): predictions, never limits;
  [assumed]    typical figures that neither states. Only the ranking models use these.

`resolve()` picks the limits that apply to one board (layer count, copper weight per layer, mask
colour, assembly service) and keeps where each came from, so every finding names its limit's origin.
"""
import hashlib
import os
import tomllib

HERE = os.path.dirname(os.path.abspath(__file__))
OZ_MM = 0.035          # finished copper thickness of 1 oz/ft2
META = ("page", "quote", "measurand", "certainty", "use", "note", "conflict")


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


def split(x):
    """(value, meta) of a profile entry: an inline table {value, page, quote, ...} or a bare value."""
    if isinstance(x, dict) and "value" in x:
        return x["value"], {k: x[k] for k in META if k in x}
    return x, {}


def lint(prof):
    """Problems of a profile: every published value needs a page, the page's words and a measurand."""
    out = []
    pages = prof.get("fab", {}).get("pages", {})
    for sec, body in prof.get("published", {}).items():
        items = body if isinstance(body, list) else [body]
        for i, it in enumerate(items):
            if isinstance(body, list) or sec == "part_spacing":          # tables: one page and quote per row
                for k in ("page", "quote"):
                    if not it.get(k):
                        out.append(f"published.{sec}[{i}]: no {k}")
                if it.get("page") and it["page"] not in pages:
                    out.append(f"published.{sec}[{i}]: page '{it['page']}' is not in [fab.pages]")
                continue
            for key, v in it.items():
                val, meta = split(v)
                if not meta:
                    out.append(f"published.{sec}.{key}: a bare value (no page, quote, measurand)")
                    continue
                for k in ("page", "quote", "measurand"):
                    if not meta.get(k):
                        out.append(f"published.{sec}.{key}: no {k}")
                if meta.get("page") and meta["page"] not in pages:
                    out.append(f"published.{sec}.{key}: page '{meta['page']}' is not in [fab.pages]")
                if meta.get("certainty", "stated") not in ("stated", "interpreted", "chosen", "confirmed"):
                    out.append(f"published.{sec}.{key}: certainty '{meta['certainty']}'")
                if meta.get("certainty") in ("interpreted", "chosen", "confirmed") and not meta.get("note"):
                    out.append(f"published.{sec}.{key}: {meta['certainty']} without a note saying how")
    return out


class Fab:
    """Limits for one board: `fab["mask_web_mm"]` is the value, `fab.src("mask_web_mm")` its origin."""

    def __init__(self, prof):
        self.name = prof.get("fab", {}).get("name", "?")
        self.retrieved = prof.get("fab", {}).get("retrieved")
        self.pages = prof.get("fab", {}).get("pages", {})
        self.sources = sorted(set(self.pages.values())) or prof.get("fab", {}).get("sources", [])
        self.path, self.hash = prof.get("_path"), prof.get("_hash")
        self.limits = {}       # key -> {"value", "source": published | standard | guidance | assumed | config, meta...}
        self.notes = []
        self.read = set()      # keys a check has looked up
        self.judged_by = {}    # key -> checks that compared geometry with it

    def set(self, key, value, source, note="", meta=None):
        if value is not None:
            e = {"value": value, "source": source, "note": note}
            m = dict(meta or {})
            if m.get("note"):
                e["note"] = (note + "; " if note else "") + m.pop("note")
            if m.get("page") in self.pages:
                m["url"] = self.pages[m["page"]]
            e.update(m)
            e.setdefault("certainty", "stated" if source in ("published", "standard", "config") else None)
            e.setdefault("use", "limit")
            self.limits[key] = e

    def __getitem__(self, key):
        self.read.add(key)
        return self.limits[key]["value"]

    def get(self, key, default=None):
        if key in self.limits:
            self.read.add(key)
            return self.limits[key]["value"]
        return default

    def judge(self, key, check):
        """The value of a limit a rule compares geometry with (None if the profile has none); records
        that `check` judged it."""
        if key not in self.limits or not self.judged(key):
            return None
        self.read.add(key)
        self.judged_by.setdefault(key, set()).add(check)
        return self.limits[key]["value"]

    def layer_key(self, key, layer):
        """`key@layer` where that layer's copper weight gives it its own limit, else `key`."""
        k = f"{key}@{layer}"
        return k if k in self.limits else key

    def judged(self, key):
        """A limit a rule check may judge with: published, standard or set in the board config."""
        return key in self.limits and self.limits[key]["source"] in ("published", "standard", "config")

    def certainty(self, key):
        return self.limits.get(key, {}).get("certainty") or "stated"

    def severity(self, key, sev):
        """A fail on a limit whose reading the profile had to interpret is a warning: the fabricator
        has not said what the figure measures."""
        if sev == "error" and self.certainty(key) == "interpreted":
            return "warning"
        return sev

    def src(self, key):
        """How a finding names its limit: its origin, and how its reading was settled when the page
        does not say."""
        e = self.limits.get(key, {})
        s = e.get("source")
        base = {"published": f"{self.name} published", "standard": "industry standard", "guidance": f"{self.name} guidance",
                "assumed": f"assumed, not published by {self.name}", "config": "board config"}.get(s, "default")
        c = e.get("certainty")
        if s == "published" and c == "interpreted":
            base += f" (reading interpreted, not stated: {e.get('note', '')})"
        elif s == "published" and c == "chosen":
            base += f" (reading chosen by the board owner: {e.get('measurand', '').replace('_', ' ')})"
        return base

    def describe(self):
        out = []
        for k, v in self.limits.items():
            d = dict(key=k, used=k in self.read, judged_by=sorted(self.judged_by.get(k, ())), **v)
            out.append(d)
        return {"name": self.name, "retrieved": self.retrieved, "sources": self.sources, "pages": self.pages,
                "profile": self.path, "profile_hash": self.hash, "notes": self.notes, "limits": out}

    def not_judged(self):
        """Published and standard limits (use = limit) that no check compared geometry with."""
        return sorted(k for k, v in self.limits.items() if v["source"] in ("published", "standard")
                      and v.get("use", "limit") == "limit" and k not in self.judged_by)


def _row(rows, n, oz):
    """The row of a by-layer-count, by-copper-weight table that applies; (row, note)."""
    rows = [r for r in rows if r["layers"][0] <= n <= r["layers"][1]]
    if not rows:
        return None, None
    fit = sorted((r for r in rows if r["oz"] >= oz - 0.15), key=lambda r: r["oz"])
    if fit:
        note = None
        if fit[0]["oz"] > oz + 0.15:
            note = (f"{oz:.1f} oz copper has no row for {n} layers; the {fit[0]['oz']:g} oz row is used (heavier copper needs "
                    f"the wider figure, so this is the strict side)")
        return fit[0], note
    row = max(rows, key=lambda r: r["oz"])
    return row, f"copper {oz:.1f} oz is heavier than any row for {n} layers; the {row['oz']:g} oz row is used"


def _rowmeta(row, measurand):
    return {"page": row.get("page"), "quote": row.get("quote"), "measurand": measurand} if row.get("quote") else {}


def layer_rows(prof, board):
    """Per copper layer: its copper weight (oz) and the track and annular-ring rows that apply to it."""
    pub = prof.get("published", {})
    n = len(board.layers)
    out = {}
    for c in board.stack.copper:
        oz = c["t"] / OZ_MM if c.get("t") else 1.0
        tr, tnote = _row(pub.get("track", []), n, oz)
        ring, rnote = _row(pub.get("pth_ring", []), n, oz)
        out[c["name"]] = {"oz": round(oz, 2), "t_um": round(c["t"] * 1000, 1), "track": tr, "track_note": tnote,
                          "ring": ring, "ring_note": rnote}
    return out


def resolve(prof, board, cfg):
    """The limits that apply to this board. `[dfm.fab_limits]` in the board config overrides any key.
    Copper-weight limits are kept per layer (`fab.per_layer`); the board-wide key holds the outer
    layer's value."""
    fab = Fab(prof)
    pub, std, asm = prof.get("published", {}), prof.get("standard", {}), prof.get("assumed", {})
    d = cfg["dfm"]
    n = len(board.layers)

    def put(sec, k, key=None, note="", source="published"):
        if k in pub.get(sec, {}):
            v, meta = split(pub[sec][k])
            fab.set(key or k, v, source, note, meta)

    oz = board.stack.copper[0]["t"] / OZ_MM if board.stack.copper else 1.0
    what = f"{n} layers, {board.stack.copper[0]['t'] * 1000:.0f} um outer copper" if board.stack.copper else f"{n} layers"
    # per layer: the copper weight of each layer picks its own row (inner layers are often 0.5 oz)
    fab.per_layer = layer_rows(prof, board) if board.stack.copper else {}
    for ln, r in fab.per_layer.items():
        for k in ("track_note", "ring_note"):
            if r[k]:
                fab.notes.append(f"{ln} ({r['t_um']:g} um): {r[k]}")
    row, note = _row(pub.get("track", []), n, oz)
    if row:
        fab.set("track_width_mm", row["width_mm"], "published", what, _rowmeta(row, "track_width"))
        fab.set("track_space_mm", row["space_mm"], "published", what, _rowmeta(row, "copper_gap"))
        fab.set("bga_track_mm", row.get("bga_mm"), "published", "width and spacing accepted locally in a BGA fan-out",
                _rowmeta(row, "track_width"))
        # a layer whose own copper weight picks another row gets its own keys (key@layer)
        for ln, r in fab.per_layer.items():
            tr = r["track"]
            if tr and tr is not row:
                lw = f"{ln}, {r['t_um']:g} um copper" + (f" ({r['track_note']})" if r["track_note"] else "")
                fab.set(f"track_width_mm@{ln}", tr["width_mm"], "published", lw, _rowmeta(tr, "track_width"))
                fab.set(f"track_space_mm@{ln}", tr["space_mm"], "published", lw, _rowmeta(tr, "copper_gap"))
    else:
        fab.notes.append(f"the profile has no track row for {n} layers: track width and spacing are not judged")
    put("etch", "same_net_spacing_mm", note="between tracks of one net")
    put("etch", "track_tolerance_pct", note="finished width; recorded, not judged")
    for k in ("via_hole_min_mm", "via_diameter_min_mm", "npth_min_mm", "drill_max_mm", "via_hole_spacing_mm", "pad_hole_spacing_mm",
              "castellated_hole_min_mm", "castellated_spacing_mm"):
        put("drill", k)
    if "via_diameter_over_hole_mm" in pub.get("drill", {}):
        v, meta = split(pub["drill"]["via_diameter_over_hole_mm"])
        fab.set("via_ring_min_mm", v / 2, "published", f"via diameter at least {v:g} mm larger than its hole", dict(meta, measurand="via_ring"))
    put("drill", "plated_slot_min_2l_mm" if n <= 2 else "plated_slot_min_ml_mm", "plated_slot_min_mm", what)
    ring, note = _row(pub.get("pth_ring", []), n, oz)
    if ring:
        fab.set("pth_ring_min_mm", ring["min_mm"], "published", what + ", absolute minimum", _rowmeta(ring, "annular_ring"))
        fab.set("pth_ring_rec_mm", ring.get("recommended_mm"), "published", what + ", recommended", _rowmeta(ring, "annular_ring"))
    for k in ("via_hole_to_copper_mm", "pth_hole_to_copper_mm", "npth_hole_to_copper_mm", "inner_via_hole_to_copper_mm",
              "inner_pth_hole_to_copper_mm", "pad_to_track_mm", "copper_to_edge_mm", "copper_to_vcut_mm"):
        put("clearance", k)
    if n >= 3:
        put("clearance", "bga_pad_to_track_mm")
    put("board", "thickness_min_mm")
    put("board", "thickness_max_mm")
    vc = pub.get("via_cover", {})
    if "filled_default_from_layers" in vc:
        v, meta = split(vc["filled_default_from_layers"])
        fab.set("via_cover", "filled" if n >= v else "tented", "published",
                f"epoxy filled and capped is the default from {v} layers; this board has {n}", dict(meta, use="order"))
    for k in ("filled_order_note", "filled_via_min_mm", "filled_via_max_mm", "plugged_via_max_mm", "plugged_to_opening_mm"):
        put("via_cover", k)
    m = pub.get("mask", {})
    colour = str(d.get("mask_color", "green")).lower()
    dark_v = split(m.get("dark_colors", []))[0] or []
    dark = colour in [x.lower() for x in dark_v]
    put("mask", "web_mm", "mask_web_mm")
    if oz >= 1.5 and "pad_spacing_2oz_mm" in m:
        put("mask", "pad_spacing_2oz_mm", "mask_pad_spacing_mm", "2 oz copper, any colour")
    else:
        put("mask", "pad_spacing_1oz_dark_mm" if dark else "pad_spacing_1oz_mm", "mask_pad_spacing_mm", f"1 oz copper, {colour} mask")
    put("mask", "expansion_mm", "mask_expansion_mm", "the expansion the fabricator says it can plot")
    put("mask", "opening_to_track_mm", "mask_to_track_mm")
    put("pads", "smd_pad_to_pad_mm")
    put("pads", "smd_pad_min_mm")
    a = pub.get("assembly", {})
    service = str(d.get("assembly_service", "standard")).lower()
    sv = f"{service} assembly service"
    put("assembly", f"reflow_peak_{service}_c", "reflow_peak_c", sv)
    put("assembly", "reflow_peak_tol_c", note="recorded only")
    put("assembly", f"min_package_{service}", "min_package", sv)
    put("assembly", f"pin_pitch_{service}_mm", "pin_pitch_mm", sv)
    put("assembly", f"bga_pitch_{service}_mm", "bga_pitch_mm", sv)
    put("assembly", f"double_sided_{service}", "double_sided", sv)
    put("assembly", f"fiducials_{service}", "fiducials_required", sv)
    put("assembly", f"rails_{service}", "rails_required", sv)
    put("stencil", "thickness_mm", "stencil_thickness_mm", "standard thicknesses; the fabricator chooses")
    put("stencil", "engineer_chooses", "stencil_engineer_chooses")
    put("stencil", "aperture_min_mm", "stencil_aperture_min_mm")
    ps = pub.get("part_spacing", {})
    if ps.get("matrix"):
        fab.set("part_spacing", {"classes": ps["classes"], "matrix": ps["matrix"]}, "published", "recommended minimum",
                {"page": ps.get("page"), "quote": ps.get("quote"), "measurand": ps.get("measurand", "part_outline_to_outline")})
    fab.set("stencil_area_ratio_min", std.get("stencil", {}).get("area_ratio_min"), "standard", "IPC-7525")
    fab.set("stencil_aspect_ratio_min", std.get("stencil", {}).get("aspect_ratio_min"), "standard", "IPC-7525")
    g = prof.get("guidance", {}).get("stencil")
    if g:
        fab.set("stencil_guide", g, "guidance", "thickness the fabricator's article gives per package and pitch",
                {"page": g.get("page"), "quote": g.get("quote"), "measurand": "stencil_thickness", "use": "input",
                 "conflict": g.get("conflict")})
    fab.set("etch_factor", asm.get("etch", {}).get("etch_factor"), "assumed", "recorded only")
    fab.set("etchant", asm.get("etch", {}).get("etchant"), "assumed", "recorded only")
    fab.set("oven", asm.get("assembly", {}).get("oven"), "assumed", "recorded only")
    fab.set("h_w_m2k", asm.get("assembly", {}).get("h_w_m2k"), "assumed", "ranking models only")
    for k, v in prof.get("notes", {}).items():
        fab.notes.append(f"{k}: {v}")
    for k, e in fab.limits.items():
        if e.get("conflict"):
            fab.notes.append(f"{k}: {fab.name}'s pages disagree: {e['conflict']}")
    for k, v in (d.get("fab_limits") or {}).items():
        fab.set(k, v, "config", "set in [dfm.fab_limits] of the board config")
    for e in fab.limits.values():                  # assumed and guidance figures are never judged limits
        if e["source"] in ("assumed", "guidance") and e.get("use", "limit") == "limit":
            e["use"] = "input"
    return fab
