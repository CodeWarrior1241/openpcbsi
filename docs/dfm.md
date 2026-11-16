# Manufacturing sign-off: checks, models and files

`python3 -m openpcbsi dfm BOARD` checks a layout against what a fabricator and an assembler can
make, before the files are sent. This page describes each check, the model behind it, where its
limit comes from, the fabricator profile, the configuration keys and `dfm.json`.

```sh
python3 -m openpcbsi dfm board.kicad_pcb                       # JLCPCB profile, every check
python3 -m openpcbsi dfm board.kicad_pcb --fab myfab.toml -c boards/board.toml --waivers waivers.toml
python3 -m openpcbsi dfm board.kicad_pcb --checks etch,mask    # the geometric checks only
```

Exit codes are those of `evaluate`: 0, 1 (warnings only), 2 (any fail), 3 (could not run).
Findings use the same model and the same waiver file format ([evaluation.md](evaluation.md)).

## What is judged, and against what

| Check | Group (`--checks`) | Finds | Limit | Origin of the limit |
|-------|--------------------|-------|-------|---------------------|
| `fab_track` (MF3) | `track` | Tracks narrower than the fabricator's minimum for the board's layer count and outer copper weight; a lower figure applies inside BGA and connector fields when the profile gives one | fail below | profile, published |
| `acid_trap` (MF1) | `etch` | Corners in copper of one net with an opening angle below the capillary limit | warn when the pocket's flush time is at least `flush_warn_s`; information below | contact angle assumed; flush limit proposed |
| `same_net_gap` (MF2) | `etch` | Slots, channels and enclosed openings in copper of one net narrower than the fabricator's same-net spacing | warn | profile, published |
| `mask_web` (MF4) | `mask` | SMD pads of different nets closer than the fabricator's minimum (fail); a mask web narrower than the minimum (fail); openings that merge (warn) | | profile, published |
| `mask_exposed` (MF5) | `mask` | Copper of another net inside a pad's mask opening (fail), or closer to it than the opening-to-track clearance (warn) | | profile, published |
| `pad_balance` (MF6) | `balance` | Two-terminal chips whose pads differ in the copper leaving them, in sitting in a pour, or in a via in the pad (warn up to `small_package`, information above); a via into a plane near one pad only (information) | `ratio_warn` | proposed |
| `reflow` (MF7) | `reflow` | Parts whose peak temperature, time above liquidus or ramp falls outside the limits; peak spread over the board | warn | proposed (`[reflow]`) |
| `tombstone` (MF8) | `tombstone` | Chips whose two joints reach liquidus more than `dt_warn_s` apart | warn up to `small_package` | proposed |

A check status is one of `fail`, `warn`, `pass`, `assumed inputs` (nothing found, but the judgement
rests on an assumed input, named in the scorecard), `nothing to check`, `unvalidated` and `not run`.

Not repeated here, because KiCad's DRC does them and `evaluate` runs it: clearance between different
nets, annular rings, copper slivers, starved thermal reliefs, hole sizes, silkscreen.

Not built: paste aperture and stencil checks, part-to-part spacing, via tenting and via-in-pad rules,
the etched cross-section (etch factor), shadowing by tall parts.

## Models

### Etching (MF1, MF2)

A ball of radius r is rolled through the space between copper; what it cannot reach is a pocket.
r stands for the reach of the etchant flow and defaults to half the fabricator's same-net spacing
(0.125 mm for JLCPCB), so that every same-net gap below that spacing is a pocket. A pocket counts
only where all its walls are copper of one net: a gap between different nets is clearance.

- **Wedge.** With E the distance from the wedge tip to the nearest ball centre, the opening angle is
  2 asin(r / E). A rounded or filleted corner gives a larger angle, as it should.
- **Capillary limit.** Liquid with contact angle θ is held in a wedge of half-angle α when
  α + θ < 90° (Concus and Finn, 1969). The limit on the opening angle is therefore 180° − 2θ. With the
  assumed θ = 45° that is 90°, the usual rule of thumb. θ for an etchant on resist and laminate is
  not a published figure: set `contact_angle_deg` if the fabricator gives one. A corner is reported
  below the limit less `angle_tol_deg` (measurement margin).
- **Flush time.** D ∇²τ = −1 on the pocket, τ = 0 at its mouth, no flux through the copper: the mean
  time for a dissolved ion to diffuse out. A dead-end slot of length L gives L²/2D, a sector of
  radius R gives R²/4D. Exchange through the top of the resist channel is ignored, so τ ranks
  pockets on one board; it does not predict a fabricator's process. A full free-surface simulation
  of spray etching is not attempted: its inputs (spray pressure, chemistry, resist profile) are not
  known to the designer.
- **Resolution.** Candidates come from the board raster (25 µm); each is then measured on a raster
  `hires_div` times finer (6.25 µm), drawn from the vector data in `tile_mm` tiles. A slot that
  crosses a tile border is reported with "at least" its measured length and without a flush time.

### Solder mask (MF4, MF5)

The snapshot carries each pad's mask opening (pad shape grown by its mask expansion). MF4 measures
the distance between openings and between the pads' copper as polygons. MF5 draws each opening and
the copper around it at the fine resolution; copper under another pad's own opening belongs to MF4.

### Reflow, board level (MF7)

The 2.5-D finite-volume model of the thermal check (one node per copper layer per cell, copper
fraction from the raster, vias through the stack) with heat capacity: laminate, copper, and each
part as a lump on the cells under it. Both faces, and the side faces of each part, exchange heat with
the oven air through a film coefficient h. The air temperature follows the oven's zones in time
(backward Euler). With `calibrate` the spike zones are shifted until the median peak over the board
is the assembler's stated reflow temperature.

Part bodies: chips from a package table (size from the footprint name; ceramic capacitor, resistor,
ferrite or moulded by the reference prefix); other parts from `[reflow.parts.REF]`, from footprint
fields `Mass` (g) and `Height` (mm), or else a block over the courtyard with `default_height_mm`.
The report says which, and counts the assumed ones.

### Tombstoning, pad level (MF8)

In a ramp every node rises at the ramp rate and lags the air by its lag time
τ = (K + H)⁻¹ C 1, with K the conduction matrix, H the film conductances and C the heat capacities.
Two joints of one chip therefore reach liquidus τ₁ − τ₂ seconds apart at any ramp rate. τ is solved
in a window around the chip: its own copper layer at 50 µm, the other layers at 200 µm, the chip as
two nodes that take heat from the air over their surface, pass it to the pads through the unmelted
paste and are joined through the body.

In this model a small chip runs warmer than the board, so heat flows from the chip into the pads.
The pad with the lower resistance into the board stays cooler and melts last. On the pilot board the
dominant path is a via into a plane near one pad, not the width of the copper on the outer layer:
the pad with such a via melts last in every chip that has one on one side only, while a pad in an
outer-layer pour differs little from a pad on a track. MF6 therefore reports both, and MF8 puts a
number on them.

### Validation

`python3 -m openpcbsi validate` runs the closed-form cases; a check whose cases have not passed is
`unvalidated` and never a pass.

| Solver | Cases |
|--------|-------|
| `etch` | opening angle of drawn wedges of 30° to 90° within 6 %; width of a drawn slot; flush time of a dead-end slot (L²/2D) and of a sector (R²/4D) |
| `reflow` | step response and ramp lag of a uniform board (lumped body); lag time from the steady solve; the two-resolution window model against the lumped value |

None of this is correlated with fabricated boards or with an assembler's defect data. The limits
marked *proposed* are starting points to be calibrated.

## Fabricator profile

A TOML file: `[published]` holds the fabricator's own figures, `[assumed]` typical figures the
fabricator does not state. [openpcbsi/fabs/jlcpcb.toml](../openpcbsi/fabs/jlcpcb.toml) is the
exemplar, from JLCPCB's capability pages on the date in the file. `--fab NAME` reads
`openpcbsi/fabs/NAME.toml`; `--fab FILE` reads a file. The report lists every limit with its origin.

| Limit key | From the profile | Used by |
|-----------|------------------|---------|
| `track_width_mm`, `track_space_mm`, `bga_track_mm` | `[[published.track]]` row for the board's layer count and outer copper weight | MF3 |
| `same_net_spacing_mm`, `track_tolerance_pct` | `[published.etch]` | MF1, MF2 |
| `mask_web_mm`, `mask_pad_spacing_mm`, `mask_expansion_mm`, `mask_to_track_mm` | `[published.mask]` (pad spacing by copper weight and mask colour) | MF4, MF5 |
| `smd_pad_to_pad_mm` | `[published.pads]` | MF4 |
| `reflow_peak_c`, `reflow_peak_tol_c` | `[published.assembly]`, by `assembly_service` | MF7 |
| `etch_factor`, `etchant`, `oven`, `h_w_m2k` | `[assumed]` | MF7, MF8 (`h_w_m2k`); the others are recorded only |

## Configuration

Every key has a default.

```toml
[dfm]
fab = "jlcpcb"                 # profile name or file
mask_color = "green"           # black and white need a wider pad spacing at some fabricators
assembly_service = "standard"  # which reflow temperature of the profile applies

[dfm.fab_limits]               # override any resolved limit; marked "config" in the report
mask_web_mm = 0.12

[dfm.acid]
contact_angle_deg = 45.0
flush_warn_s = 5.0

[dfm.pad_balance]
max_package = "1206"           # chips up to this size are checked
small_package = "0603"         # up to this size a finding warns

[dfm.tombstone]
dt_warn_s = 0.3

[reflow]
zones_c = [110, 140, 160, 175, 185, 200, 260, 270]   # oven air per zone
zone_s = 30
liquidus_c = 217
peak_min_c = 230
peak_max_c = 250

[reflow.parts.L1]
mass_g = 2.4
height_mm = 4.0
```

## Files

Written to `runs/<board>_<hash>/dfm/` (or `-o DIR`):

| File | Content |
|------|---------|
| `dfm.html` | Scorecard, the limits used and their origin, findings, per-layer etch table, reflow table and map, chip table. |
| `dfm.json` | Schema `openpcbsi-dfm/1`, below. |
| `layer_<layer>.png` | Findings on each copper layer. |
| `reflow_peak_top.png` | Peak temperature of the top layer. |
| `<board>_openpcbsi.kicad_pcb` | The board with a marker per warning and failure on User.9. |

### `dfm.json`

| Field | Content |
|-------|---------|
| `schema`, `tool`, `version`, `date`, `seconds` | Schema id, `openpcbsi`, its version, time and duration of the run. |
| `board`, `board_sha256`, `kicad`, `config`, `config_hash` | The inputs. |
| `fab` | `name`, `retrieved`, `sources`, `profile`, `profile_hash`, `notes`, and `limits`: a list of `{key, value, source: published / assumed / config, note}`. |
| `waivers_file`, `waivers_hash`, `waivers` | As in `evaluation.json`. |
| `solver_versions` | Versions of the `etch`, `reflow` and `thermal` solvers. |
| `exit_code`, `status`, `counts` | Verdict (`fail`, `warn`, `pass`) and totals outside areas and waivers. |
| `checks` | Per check: `check`, `requirement`, `status`, counts, `validation`, `why`. |
| `findings` | Every finding in the model of `evaluation.json` (stable `id`, `severity`, `check`, `requirement`, `value`, `limit`, `unit`, `nets`, `refs`, `layer`, `x`, `y`, `reason`, `hint`, `waived`). |
| `rules` | Per check: `params`, `metrics`, `summary`, its own `findings`, and for `reflow`, `tombstone` and `pad_balance` a `parts` list (per part: peak, time above liquidus, ramp, lag; time between the joints and which melts first; copper leaving each pad). |
| `images`, `marker_board` | File names of the layer images and the marker board. |
