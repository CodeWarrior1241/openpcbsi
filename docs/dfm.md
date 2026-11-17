# Manufacturing review: rules, rankings, limits and files

`python3 -m openpcbsi dfm BOARD` reviews a layout against what a fabricator and an assembler
publish about their processes, before the files are sent.

```sh
python3 -m openpcbsi dfm board.kicad_pcb                              # JLCPCB profile, everything
python3 -m openpcbsi dfm board.kicad_pcb --gerbers out/gerbers        # check the files that will be sent
python3 -m openpcbsi dfm board.kicad_pcb --fab myfab.toml -c boards/board.toml --waivers waivers.toml
python3 -m openpcbsi dfm board.kicad_pcb --checks copper,drill,mask   # some groups only
```

## What it is, and what it is not

It is a review, not a sign-off. It reads the board and the files that will be sent; it cannot see
what the fabricator's CAM does to them, and the report lists what it does not check.

It has two kinds of output and keeps them apart.

- **Rules** judge drawn geometry against a limit that the fabricator publishes, that an industry
  standard gives, or that the board config sets. Only rules give the verdict and the exit code.
- **Rankings** order places on one board by a measurement or a model that no published limit covers
  and that has not been compared with fabricated boards or with an assembler's defect data. A
  ranking never passes, never fails and does not change the exit code. It gets a limit only when the
  board config supplies one with a source that can be checked (`[dfm.judge]`).

Rules apply everywhere on the board. There are no exempt areas: inside a BGA field the figure the
fabricator publishes for BGA fields is used, where it publishes one.

### Verdict and exit code

| Status | Exit | Means |
|--------|------|-------|
| `pass` | 0 | Every rule that ran passed, and nothing it relied on was assumed. |
| `pass with assumptions` | 1 | Every rule passed, but the verdict rests on an order option the board config does not state (the via covering, the stencil thickness). The report lists them under *Assumed, not stated*; set them in `[dfm]` or state them in the order. |
| `warn` | 1 | At least one warning, no fail. |
| `fail` | 2 | At least one fail. |
| `could not check` | 3 | The review could not rely on its input: the board could not be read, or the two independent readings of the board differ (see [Reading the board](#reading-the-board)). |

A rule's own status is `fail`, `could not check`, `warn`, `pass`, `nothing to check` or `not run`.

## What a limit measures

Every published value in a profile carries the page it is from, the page's own words, what it
measures, and how certain that reading is. Findings name the limit's origin, and the report's
*Limits used* table shows all of it.

| Certainty | Means | In a rule |
|-----------|-------|-----------|
| `stated` | The page says what the figure measures. | May fail a board. |
| `interpreted` | The page does not say; the profile had to choose a reading (named in its note). | Warns at most: a fail on it is shown as a warning, marked "reading interpreted". |
| `chosen` | The board owner chose the reading (named in the note). | May fail a board. |
| `confirmed` | The fabricator answered the question (the note says who and when). | May fail a board. |

A value's *use* is `limit` (judged by a rule), `order` (an order option: the order checklist),
`input` (an input to a model or a choice, such as the standard stencil thicknesses) or `record`
(recorded, not judged). A full run lists published limits that no check compared geometry with in
`limits_not_judged` (acceptance test 39 requires the list to be empty).

For JLCPCB, read from its pages on 2026-11-16:

- **Hole-to-hole spacing** (0.2 mm vias, 0.45 mm pad holes): no page says whether it is wall to wall
  or centre to centre. It is judged centre to centre, the board owner's choice. Holes that touch or
  overlap fail as such (geometry), and a pair that meets the figure centre to centre but not wall to
  wall is noted (`info`).
- **The smallest SMD pad** ("Minimum SMD pad: 0.25mm × 0.25mm") is read as too small in both
  directions: interpreted.
- **Filled vias**: the capability page says "Compatible with via diameters from 0.15 to 0.55 mm"; the
  via-covering article says "Plugged and epoxy-filled vias holes should not be larger than 0.5 mm".
  The drill is judged, 0.15 to 0.5 mm, and the conflict is recorded.
- **Castellated holes**: 0.5 mm diameter and 0.5 mm edge to edge from the capability page; the
  castellated-PCB guide's 0.6 mm is recorded as conflicting.
- **The reflow temperature** ("Reflow temperature: 240+/-5 °C") is read as the peak: interpreted.
- **Copper weight**: limits that depend on it are chosen per layer from that layer's own copper in
  the stackup. JLCPCB has no row for inner layers or for 0.5 oz (its default inner copper); the 1 oz
  row is used there, the strict side, and the report says so.

## Rules

| Check | Group (`--checks`) | Judges | Limit from |
|-------|--------------------|--------|------------|
| `rule_audit` (MF19) | `audit` | Each of KiCad's board-setup rules against every profile figure of the same measurand (KiCad's one hole clearance against the via, plated-hole, inner-layer and NPTH figures; its clearance against track-to-track and pad-to-track), the custom rules of the `.kicad_dru` that relax one of them, DRC tests set to ignore, rules not set. A rule is "as strict" only for the figures it meets, and the verdict names the rest. SMD pad-to-pad and the non-plated hole minimum have no KiCad board-setup rule that could hold them without forbidding other features; `mask_web` and `fab_drill` measure them. | profile |
| `output_match` (MF20) | `output` | The files that are sent against the geometry the other checks worked on (see [The files](#the-files)). | the files |
| `fab_track` (MF3) | `copper` | Track width per layer, at the row of that layer's copper weight; the BGA figure inside BGA fields. | published |
| `fab_space` (MF12) | `copper` | The gap between copper of different nets, measured: track to track and pad to track, per layer's copper weight, with the BGA figures inside BGA fields. A gap within the measuring slack of the limit is `marginal` (a warning). | published |
| `same_net_gap` (MF2) | `copper` | Gaps in copper of one net narrower than the same-net track spacing. Warns where a track bounds the gap; a gap between pads or pour copper only is information, since the published figure is for tracks. | published |
| `fab_drill` (MF13) | `drill` | Via hole, via pad and via ring; through-hole annular ring (fail below the absolute minimum, warn below the recommended); non-plated hole, largest drill, plated slot width; castellated holes (diameter, spacing); hole to hole and holes that cut into each other; hole to copper of another net on every layer the hole passes; copper to the board edge (castellated pads excepted); board thickness. Slots are slots: an oval hole is a stadium, so its ends are judged. | published |
| `mask_web` (MF4) | `mask` | Pad spacing; the mask web as drawn (fail below the minimum, warn where openings merge); the web that would remain at the expansion the fabricator says it plots; pads of one net with no web between them (warn between parts); the smallest SMD pad, measured from its copper (a custom pad's anchor is not its size). | published |
| `mask_exposed` (MF5) | `mask` | Copper of another net inside a pad's opening as drawn (fail), or closer to it than the opening-to-track clearance (warn). | published |
| `via_pad` (MF14) | `mask` | A via in an SMD pad: with filled and capped vias, its drill against the process range; otherwise a fail. An untented via with no mask web to a pad; plugged-via rules. The covering is `[dfm] via_cover`; without it the fabricator's default for the layer count is assumed and the verdict is at best `pass with assumptions`. | published |
| `paste` (MF15) | `paste` | Every paste aperture, from its own polygon (margins applied, custom apertures as drawn): area ratio and aspect ratio at the stencil thickness (IPC-7525: at least 0.66 and 1.5; the aperture's width is the largest disc that fits in it); the aperture against the stencil's minimum; an SMD part with no paste. See [Stencil thickness](#stencil-thickness). Also reports which standard thicknesses print every aperture. | standard, published |
| `footprint` (MF16) | `parts` | Two-pad chip land patterns: fails when the pads are further apart than the chip's EIA body is long, so the part cannot reach both pads. The chip table holds nominal lengths with no sourced tolerance, so the nominal length is the limit and the finding says so. Nothing else about the land pattern is judged; the toe of the pads beyond the body and the ratio of the pad areas are listed in the `chips` ranking. | EIA body size |
| `part_spacing` (MF17) | `parts` | The distance between neighbouring parts, outline to outline (pads and the body drawn on the fabrication layer), as the assembler's measurement diagram shows, against its table by package class. A part's class comes from its EIA size, KiCad's library (`Package_DFN_QFN`, `Package_QFP`, ...), its name, then its pads; a part of no class is listed with the reason, not judged. | published (recommended) |
| `assembly` (MF18) | `parts` | Smallest package, pin and ball pitch, parts on both sides, fiducials, for the assembly service chosen. | published |

### Order checklist

The report lists what the order must state for the verdict to hold, since the fabricator otherwise
chooses: the via covering (with the number of vias in pads and the fabricator's own words on saying
which to fill), the stencil thickness, edge rails for the assembly service, the solder-mask colour.

### Stencil thickness

The paste rules depend on the stencil thickness, and the fabricator chooses it unless the order
states one. JLCPCB's ordering page: "You can either select a thickness by yourself, or you can just
let JLCPCB engineer to select the appropriate thickness for you based on your Gerber file features."

- With `[dfm] stencil_mm` set, the paste rules judge at that thickness.
- Without it, the review predicts the fabricator's choice from the profile's `[guidance.stencil]`
  table. For JLCPCB that is the "commonly used" column of its article "How to choose a solder paste
  stencil": a thickness per package type and pitch. The board's finest-pitch part decides, as
  JLCPCB's stencil design guide says ("Give most of your attention to the board's smallest pitch
  component"). The report names the part and the table row.
- A predicted thickness is an assumption: the verdict is at best `pass with assumptions`, and the
  order checklist asks for that thickness. Setting `stencil_mm` removes the assumption.
- A part the table does not cover (01005 chips, pitches finer than its rows) leaves the thickness
  unknown: the paste check is `could not check` until `stencil_mm` is set.

The guide is a JLCPCB article, not its capability page, so its figures are marked as guidance. Its
pages disagree in places: the article counts 0.13 mm as a free thickness, which the capability and
ordering pages do not, and its quick answer gives 0.10 mm for every pitch up to 0.5 mm while its table
gives 0.12 mm for 0.5 mm-pitch ICs. The table is used.

### How the measurements are made, and how exact they are

- **Pads** are measured from their polygons. A custom pad's anchor (often 0.2 mm) is kept as
  `anchor` and no check reads it; paste-only and mask-only pads are measured from their own openings.
- **Rasters.** Copper is drawn per layer at 25 µm with one label per net. Pads without a net
  (fiducials, mounting pads) get a label of their own, so they are copper for every check. Candidate
  places are found there and measured again on a finer raster drawn from the vector data.
- **Spacing (MF12)** is measured at 6.25 µm. A gap more than one fine pixel below the limit fails; a
  gap within one pixel of the limit either side is `marginal`, a warning: neither a pass nor a fail
  can be shown.
- **Hole to copper and copper to the edge** are measured at 2.5 µm against the exact hole (a circle
  or a slot's stadium) and the board outline polygon. The slack is 3.3 µm (half a pixel and KiCad's
  2 µm arc approximation). A reading within the slack of the limit is `marginal`; KiCad's DRC with the
  clearance set to the limit settles such a place exactly.
- **Same-net gaps (MF2).** A ball of radius 0.125 mm is rolled through the space between copper;
  what it cannot reach, where every wall is copper of one net, is a pocket. A pocket of constant
  width is a gap. The radius is a measuring scale, not a property of a process.

A rule that passes means: the geometry meets the limit plus the stated slack.

### The files

With `--gerbers DIR` the files in that folder are read: these should be the files that are sent.

- **Feature by feature**, to 1 µm: every straight track as a draw with the same ends and width,
  every pad and via as a flash at its centre with its size, every paste aperture as a flash with the
  size of the board's aperture. Arc tracks and custom pads, which plotter and reader cut into other
  chords or regions, are left to the image comparison.
- **Images**, at 25 µm, for everything else on copper, mask and paste: differences wider than two
  pixels.
- **Holes** from the drill files, by position and size.
- **Outline**: every draw of the Edge.Cuts file on the board outline, and its length the outline's.
- **Placement** (a `pos` file among the files): fitted parts, positions to 2 µm, sides.
- **BOM** (a CSV with a Reference or Designator column among the files): the board's fitted parts,
  DNP and BOM-excluded footprints respected.

Without `--gerbers`, the files are plotted with `kicad-cli` into `dfm/fab_files/`. Comparing those
only shows that the plotter and the reader agree: the check's status is `could not check` (a
self-test), never a pass. A construct the reader does not draw (an aperture with a hole,
step-and-repeat) is named and that file is not compared in part.

### What has been compared with something independent

- **The board reading** against KiCad's own (acceptance tests 50 and 51): see
  [Reading the board](#reading-the-board).
- **Spacing against KiCad's own DRC** (test 46). On Antmicro's SI test board and, when it is given,
  any other board, every clearance is set to one figure and KiCad's DRC is run with all track errors
  reported. Every net pair KiCad reports clearly below the limit must be measured below it, nothing
  may be reported here that KiCad lacks, and the gaps must agree. On an 8-layer FPGA board at 0.13 mm:
  107 net pairs reported by both, gaps within 4 µm.
- **The files against the board** compare two independent paths: KiCad's plotter and this tool's own
  reading of the board. Test 45 plants file faults the images cannot see: tracks plotted 40 µm narrow,
  a layer moved by 10 µm, a missing paste aperture, a cut outline, a moved part, a BOM missing a part.
- **The pilot board's known answers** (test 51): a custom pad is judged by its copper, paste
  apertures by their own outline, 0.3 mm vias in pads against the drill range, part spacing outline to
  outline.
- The other rules are tested on generated boards with one planted defect each and the same board
  without it (tests 37 to 49). That shows each check finds what it was written to find, and no more.

## Reading the board

The board is read from the file by openpcbsi's file reader, on kicad-monkey (pinned in
`deps/parsers.lock`). Zones are refilled first by `kicad-cli pcb drc --refill-zones --save-board` on a
copy in the cache; the board file itself is never written. Every polygon (pads of every shape,
custom pads with all their primitives, mask openings and paste apertures with their margins, the
board outline) is computed by openpcbsi.

For the review, the board is read a second time by an independent parser, KiCadFiles, and the two
readings are compared on the fields KiCadFiles reads: nets, tracks and arcs, via position, size, drill
and net, footprint placement and value, pad number, type, net, size, drill and position, and zone fill
area per net and layer. KiCadFiles turns KiCad 10's name-only nets into net 0; for KiCad 10 files the
net names of its items are read from the file text by openpcbsi instead (`net_names_from` in the
report). It does not read pad and via copper layers, custom-pad primitives, mask and paste openings,
the outline, dielectric plies or footprint attributes; the report lists these as not cross-checked. A
difference makes the review `could not check` (exit 3) and blocks `signoff`. `--no-cross-check` skips
the second reading.

Both readers are compared with KiCad's own reading (the pcbnew snapshot) and KiCad's IPC-2581 export
on a reference corpus of real boards (acceptance test 50): the pilot board, Antmicro's SI test board
and five of KiCad's demo boards, among them boards with 71, 39 and 16 custom pads. KiCad's pcbnew
reader remains available (`OPENPCBSI_READER=pcbnew`) while KiCad ships it.

## Rankings

| Ranking | Group | Orders | Basis |
|---------|-------|--------|-------|
| `acid_trap` (MF1) | `copper` | Corners in copper of one net with an opening angle below 90° (less a measurement margin), thinnest track first. | The drafting rule of thumb about acute angles. No fabricator in the profile publishes a limit for it, and nothing here models etching. |
| `copper` (MF21) | `copper` | Copper share of each layer and the difference between each pair of layers about the centre of the stack. | No published limit. |
| `reflow` (MF7) | `reflow` | Parts by how far their modelled peak temperature lies from the board's median, with its sensitivity to the model's parameters. | Thermal model, below. |
| `chips` (MF8) | `chips` | Two-terminal chips by the modelled time between their two joints reaching liquidus, with its sensitivity, and the land geometry. | Thermal model, below. |

**Acute corners.** With E the distance from a wedge's tip to the nearest ball centre, the opening
angle is 2 asin(r / E), which does not depend on r. A corner formed by a track is listed with the
track's width, since a thin track has the least copper to lose; corners between pads and pour copper
only are counted and not listed. There is no flush-time or capillary model.

**Reflow (MF7).** The 2.5-D finite-volume model of the thermal check with heat capacity: laminate,
copper, and each part as a lump on the cells under it. Both faces and the sides of each part exchange
heat with the oven air through a film coefficient; each point of the board enters a zone when the
belt carries it there. The spike zones are shifted until the median peak over the parts is the
assembler's stated reflow temperature, so the absolute peak is an input and only the differences
between parts are a result. The run is repeated with the film coefficient at 0.6 and 1.6 times; each
part's difference from the median is given with its spread over those runs. That spread is the
result's sensitivity to that parameter, not an error bar: the runs share one model form.

Only soldered parts are modelled. Fiducials, test points, mounting holes, board-only and BOM-excluded
footprints and footprints with no paste are left out and listed. A part's body is its `[reflow.parts]`
entry in the board config, a Mass or Height field, or the chip table; a part with none of these is
left out of the model and listed (`left_out`), with what to give for it, instead of being given an
assumed body. Its thermal mass is then missing from the model, which the report says.

**Chips (MF8).** A window around each chip: its own copper layer at 50 µm, the others at 200 µm, the
chip as two nodes heated over their own surface and joined to the pads through the unmelted paste and
to each other through the body. The window's edge takes its values from the board model, so planes
do not end at the window. Two numbers, each with its own sensitivity:

- *dt*: the difference of the two joints' lag times, τ = (K + H)⁻¹ C 1. In a slow ramp the two
  joints reach liquidus that far apart. Solved for every chip; its sensitivity is over the paste
  conductance (0.5 and 2 times), the film coefficient (0.6 and 1.6 times) and a window twice as large.
- *in this oven*: the same window marched through the oven profile, giving the time between the two
  joints crossing liquidus directly. It also includes the belt carrying one pad into each zone before
  the other, so it is another quantity than *dt*. Its sensitivity is over the paste conductance and
  the window, each run through the oven. Run for the first `detail` layouts.

A self-check asserts that each nominal figure lies inside its own sensitivity range. Chips with the
same layout and nearly the same result share a row. Timing is one input to tombstoning. The wetting
forces (pad geometry, paste volume and offset, placement, alloy) are not modelled; the toe of the pads
beyond the body and the ratio of the pad areas are listed beside the timing, with the copper leaving
each pad, whether it sits in a pour, and vias in or near it.

**What the models leave out.** Latent heat, flux, shadowing by tall parts, the non-uniform film
coefficient of real ovens, infrared heating, the temperature inside a part (a part is as warm as the
board under it, so the centre of a BGA cannot lag its corners), panels, and the second pass of a
double-sided board. The oven profile is generic.

**Solver cases.** `python3 -m openpcbsi validate` runs closed-form cases: step response and ramp lag
of a uniform board, the lag time, the two-resolution window against the lumped value; for the
measurements, the angle of drawn wedges, the width of a drawn slot and the gap between two nets
drawn at five angles. They show that the code solves its equations, not that the equations describe
an oven.

### Giving a ranking a limit, and scoring it

```toml
[dfm.judge]
source = { kind = "observed", file = "defects.toml", sha256 = "3f2a..." }   # required, and checked
joint_dt_s = 0.5          # chips whose joints melt further apart warn
reflow_cold_c = 8.0       # parts peaking more than this below the median warn
acute_deg = 60.0          # corners sharper than this ...
acute_track_mm = 0.15     # ... on a track this thin or thinner warn
```

The source must be a table that can be checked: `kind = "observed"` names the `--observed` file and
its SHA-256 (verified), `kind = "assembler"` needs a `document` and a `date`, `kind = "standard"` a
`document`. A free-text source is refused (exit 3). Findings from these limits carry the source as
their basis and count in the exit code.

`--observed defects.toml` scores the rankings against defects seen on built boards:

```toml
boards_built = 40          # required: how many boards were assembled

[[defect]]
ref = "C12"
kind = "tombstone"
count = 3                  # on how many of the boards (default 1)
```

For each ranking it reports the share of (defective, sound) pairs it puts in the right order (0.5 is
chance, 1 is perfect), with its bootstrap 95 % interval and the counts it rests on. With fewer than 5
defective parts it reports "insufficient data" instead of a score. That is the evidence a limit in
`[dfm.judge]` should rest on.

## Not checked

- What the fabricator's CAM does to the files: etch compensation, changes to mask openings,
  panelisation.
- The etched cross-section and the finished track width (the published tolerance is recorded).
- Drill aspect ratio, copper balance and warpage against a limit (none is published; the figures are
  listed).
- Part height, shadowing in reflow, the second pass of a double-sided board.
- A land pattern against its part's data sheet. Only two-pad chips are judged, and only on whether
  their EIA body can reach both pads.
- Stencil aperture shapes and paste volume beyond the IPC-7525 ratios.
- Silkscreen, test access, panel rails and tooling holes.

## Fabricator profile

A TOML file. [openpcbsi/fabs/jlcpcb.toml](../openpcbsi/fabs/jlcpcb.toml) is the exemplar, from
JLCPCB's capability, assembly, stencil, via-covering and part-spacing pages on the date in the file.
`--fab NAME` reads `openpcbsi/fabs/NAME.toml`; `--fab FILE` reads a file.

| Section | Holds | Used by |
|---------|-------|---------|
| `[fab.pages]` | The pages the figures are from, by a short key. | Every published value names its page. |
| `[published]` | The fabricator's own limits, each as `{ value, page, quote, measurand, certainty, use, note, conflict }`. | Rules, the order checklist. |
| `[standard]` | Figures of an industry standard (IPC-7525 stencil ratios). | Rules. |
| `[guidance]` | The fabricator's design articles: the stencil thickness per package and pitch. | Predictions, never limits. |
| `[assumed]` | Typical figures neither states: film coefficient, etch factor, oven type. | Rankings only. No rule reads them. |
| `[notes]` | What the pages do not publish. | The report. |

`openpcbsi.fab.lint(profile)` lists published values that lack a page, the page's words or a
measurand, and readings marked interpreted or chosen without a note; acceptance test 39 requires
the JLCPCB profile to have none.

## Configuration

Every key has a default.

```toml
[dfm]
fab = "jlcpcb"                 # profile name or file
mask_color = "green"           # black and white need a wider pad spacing at some fabricators
assembly_service = "standard"  # which column of the assembler's table applies
via_cover = ""                 # tented, plugged or filled; "" = the fabricator's default (an assumption)
stencil_mm = 0.0               # the stencil thickness that will be ordered; 0 = the fabricator guide's choice (an assumption)
bga_refs = []                  # parts treated as BGA fields besides those whose footprint name says BGA
gerbers = ""                   # folder of the files that will be sent (also --gerbers)
workers = 0                    # worker processes; 0 = a quarter of the CPUs, at most 8

[dfm.fab_limits]               # override any resolved limit; marked "config" in the report
mask_web_mm = 0.12

[reflow]
zones_c = [110, 140, 160, 175, 185, 200, 260, 270]   # oven air per zone
zone_s = 30
zone_mm = 350                  # with zone_s: the belt speed
travel = "x"                   # board axis along the belt

[reflow.parts.L1]
mass_g = 2.4
height_mm = 4.0
```

## Files

Written to `runs/<board>_<hash>/dfm/` (or `-o DIR`):

| File | Content |
|------|---------|
| `dfm.html` | Verdict of the rules, assumptions and the order checklist, the rule and ranking tables, what is not checked, findings with the origin of each limit, the audit of KiCad's rules, the file comparison, the limits used (with the page's words, measurand and certainty), and each ranking. |
| `dfm.json` | Schema `openpcbsi-dfm/2`, below. |
| `fab_files/` | The Gerber, drill and placement files that were plotted and compared (without `--gerbers`). |
| `layer_<layer>.png` | Findings on each copper layer. |
| `reflow_peak_top.png` | Modelled peak temperature of the top layer. |
| `<board>_openpcbsi.kicad_pcb` | A copy of the board with a marker per warning and failure on User.9 (Cmts.User where the board has no User.9), written as text into the copy. |

### `dfm.json`

| Field | Content |
|-------|---------|
| `schema`, `tool`, `version`, `date`, `seconds` | Schema id, `openpcbsi`, its version, time and duration of the run. |
| `board`, `board_sha256`, `kicad`, `config`, `config_hash` | The inputs. |
| `reader` | `primary` (kicad-monkey or pcbnew), `kicad`, and `cross_check`: the second reader, items compared, differences, what is not cross-checked. |
| `fab` | `name`, `retrieved`, `pages`, `profile`, `profile_hash`, `notes`, and `limits`: a list of `{key, value, source, page, url, quote, measurand, certainty, use, note, conflict, judged_by}`. |
| `limits_not_judged` | Published and standard limits no check compared geometry with (a full run only). |
| `exit_code`, `status`, `counts` | Verdict of the rules (see the table above) and totals outside waivers. |
| `assumptions`, `order_checklist` | What the verdict assumed, and what the order must state. |
| `checks` | Per check: `check`, `requirement`, `tier` (`rule` or `ranking`), `status`, counts, `note`, `seconds`. |
| `findings` | Every finding in the model of `evaluation.json`, with `tier`, `basis` (the origin of its limit), `field` (BGA) and, for a fail turned into a warning because its limit's reading is interpreted, `capped`. |
| `rules` | Per rule: `params`, `metrics`, its own `findings`; `rule_audit` has `rows`. |
| `rankings` | Per ranking: `params`, `summary`, `rows` (every row), and `findings` when `[dfm.judge]` applies. |
| `not_checked`, `bga_fields`, `cases_pass`, `solver_versions` | What the review does not cover; the BGA parts; whether the closed-form cases pass. |
| `observed`, `judge` | With `--observed`: the score of each ranking; the `[dfm.judge]` table in force. |
| `waivers_file`, `waivers_hash`, `waivers`, `images`, `marker_board` | As in `evaluation.json`. |
