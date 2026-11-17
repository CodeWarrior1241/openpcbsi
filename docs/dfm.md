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

It is a review, not a sign-off. It reads the board and the plotted files; it cannot see what the
fabricator's CAM does to them, and the report lists what it does not check.

It has two kinds of output and keeps them apart.

- **Rules** judge drawn geometry against a limit that the fabricator publishes, that an industry
  standard gives, or that the board config sets. Only rules give the verdict and the exit code
  (0, 1 for warnings only, 2 for any fail, 3 if the review could not run). A rule that passes means
  one thing: the geometry meets that published number.
- **Rankings** order places on one board by a measurement or a model that no published limit covers
  and that has not been compared with fabricated boards or with an assembler's defect data. A
  ranking never passes, never fails and does not change the exit code. It gets a limit only when the
  board config supplies one together with its source (`[dfm.judge]`).

Rules apply everywhere on the board. There are no exempt areas: inside a BGA field the figure the
fabricator publishes for BGA fields is used, where it publishes one.

## Rules

| Check | Group (`--checks`) | Judges | Limit from |
|-------|--------------------|--------|------------|
| `rule_audit` (MF19) | `audit` | Each of KiCad's board-setup rules against the fabricator's limit; DRC tests set to ignore; rules not set. Warns where KiCad's DRC is looser than the fabricator, and names the check here that measures the limit itself. | profile |
| `output_match` (MF20) | `output` | The Gerber and drill files against the geometry the other checks worked on: every copper layer, both mask layers, every hole. Fails on a difference wider than two pixels (50 µm), a missing layer file, a hole that is missing, extra or of another size. | the files |
| `fab_track` (MF3) | `copper` | Track width, with the BGA figure inside BGA fields. | published |
| `fab_space` (MF12) | `copper` | The gap between copper of different nets, measured: track to track and pad to track, with the BGA figures inside BGA fields. | published |
| `same_net_gap` (MF2) | `copper` | Gaps in copper of one net narrower than the same-net track spacing. Warns where a track bounds the gap; a gap between pads or pour copper only is information, since the published figure is for tracks. | published |
| `fab_drill` (MF13) | `drill` | Via hole, via pad and via ring; through-hole annular ring (fail below the absolute minimum, warn below the recommended); non-plated hole, largest drill, plated slot; hole to hole; hole to copper of another net on every layer the hole passes; copper to the board edge; board thickness. | published |
| `mask_web` (MF4) | `mask` | Pad spacing; the mask web as drawn (fail below the minimum, warn where openings merge); the web that would remain at the expansion the fabricator says it plots; pads of one net with no web between them (warn between parts); the smallest SMD pad. | published |
| `mask_exposed` (MF5) | `mask` | Copper of another net inside a pad's opening as drawn (fail), or closer to it than the opening-to-track clearance (warn). | published |
| `via_pad` (MF14) | `mask` | A via in an SMD pad (fail unless the vias are filled and capped, then its diameter against the process range); an untented via with no mask web to a pad; plugged-via rules. How the vias are covered comes from `[dfm] via_cover`, else the fabricator's default for the layer count. | published |
| `paste` (MF15) | `paste` | Every paste aperture: area ratio and aspect ratio at the stencil thickness (IPC-7525: at least 0.66 and 1.5); the aperture against the stencil's minimum; an SMD part with no paste. Reports the thickest standard stencil at which every aperture prints. | standard, published |
| `footprint` (MF16) | `parts` | Two-pad land patterns against the nominal body of the chip their footprint is named for: pads further apart than the body can reach (fail), pads ending inside the body, unequal pads. | EIA body size |
| `part_spacing` (MF17) | `parts` | The distance between the pads of neighbouring parts against the assembler's table by package class. | published (recommended) |
| `assembly` (MF18) | `parts` | Smallest package, pin and ball pitch, parts on both sides, fiducials, for the assembly service chosen. | published |

A rule's status is `fail`, `warn`, `pass`, `nothing to check`, `could not check` or `not run`.

### How the measurements are made

- **Rasters.** Copper is drawn per layer at 25 µm with one label per net. Pads without a net
  (fiducials, mounting pads) get a label of their own, so they are copper for every check. Candidate
  places are found there and measured on a raster four times finer (6.25 µm), drawn from the vector
  data in 2 mm tiles.
- **Spacing (MF12).** At two adjacent pixels whose nearest copper belongs to different nets, the gap
  is the distance between those two copper pixels less one pixel's extent, averaged along the
  border. A gap is reported when it reads more than one fine pixel below the limit. One finding per
  pair of nets and layer, with the smallest gap and how much of the border is below the limit. Pad to
  pad is MF4's.
- **Same-net gaps (MF2).** A ball of radius 0.125 mm is rolled through the space between copper;
  what it cannot reach, where every wall is copper of one net, is a pocket. A pocket of constant
  width is a gap. The radius is a measuring scale, not a property of a process.
- **Files (MF20).** A Gerber reader (standard apertures, aperture macros, regions, arcs, both
  polarities) and an Excellon reader. With `--gerbers DIR` the files in that folder are read: these
  should be the files that are sent. Without it the files are plotted with `kicad-cli` into
  `dfm/fab_files/`, and the check then only shows that the plot and the reader agree; the report
  says so. A construct the reader does not draw (an aperture with a hole, step-and-repeat) is named
  and that file is `could not check`, never compared in part.

### What has been compared with something independent

- **Spacing against KiCad's own DRC** (acceptance test 46). On Antmicro's SI test board and, when
  it is given, any other board, every clearance is set to one figure and KiCad's DRC is run with all
  track errors reported. Every net pair KiCad reports clearly below the limit must be measured below
  it, nothing may be reported here that KiCad lacks, and the gaps must agree. On an 8-layer FPGA
  board at 0.13 mm: 107 net pairs reported by both, gaps within 4 µm.
- **The files against the board** is itself a comparison of two independent paths: KiCad's plotter
  and this tool's own drawing of the copper.
- The other rules are tested on generated boards with one planted defect each and the same board
  without it (tests 37 to 45). That shows each check finds what it was written to find, and no more.

## Rankings

| Ranking | Group | Orders | Basis |
|---------|-------|--------|-------|
| `acid_trap` (MF1) | `copper` | Corners in copper of one net with an opening angle below 90° (less a measurement margin), thinnest track first. | The drafting rule of thumb about acute angles. No fabricator in the profile publishes a limit for it, and nothing here models etching. |
| `copper` (MF21) | `copper` | Copper share of each layer and the difference between each pair of layers about the centre of the stack. | No published limit. |
| `reflow` (MF7) | `reflow` | Parts by how far their modelled peak temperature lies from the board's median, with a range. | Thermal model, below. |
| `chips` (MF8) | `chips` | Two-terminal chips by the modelled time between their two joints reaching liquidus, with a range and the land geometry. | Thermal model, below. |

**Acute corners.** With E the distance from a wedge's tip to the nearest ball centre, the opening
angle is 2 asin(r / E), which does not depend on r. A corner formed by a track is listed with the
track's width, since a thin track has the least copper to lose; corners between pads and pour copper
only are counted and not listed. There is no flush-time or capillary model: an earlier one only
restated the angle.

**Reflow (MF7).** The 2.5-D finite-volume model of the thermal check with heat capacity: laminate,
copper, and each part as a lump on the cells under it. Both faces and the sides of each part exchange
heat with the oven air through a film coefficient; each point of the board enters a zone when the
belt carries it there. The spike zones are shifted until the median peak over the parts is the
assembler's stated reflow temperature, so the absolute peak is an input and only the differences
between parts are a result. The run is repeated with the film coefficient at 0.6 and 1.6 times and
the assumed part bodies at 0.5 and 2 times; each part's difference from the median is given with its
range over those variants, and the parts in the coldest tenth of every variant are marked. Where a
part covers several cells, the spread of the board's peak under it is listed.

**Chips (MF8).** A window around each chip: its own copper layer at 50 µm, the others at 200 µm, the
chip as two nodes heated over their own surface and joined to the pads through the unmelted paste and
to each other through the body. The window's edge takes its values from the board model, so planes
do not end at the window. Two numbers:

- *dt*: the difference of the two joints' lag times, τ = (K + H)⁻¹ C 1. In a slow ramp the two
  joints reach liquidus that far apart. Solved for every chip.
- *in this oven*: the same window marched through the oven profile, giving the time between the two
  joints crossing liquidus directly. It also includes the belt carrying one pad into each zone before
  the other. Run for the first `detail` layouts, with a range of *dt* over the paste conductance (0.5
  and 2 times), the film coefficient (0.6 and 1.6 times) and a window twice as large.

On the pilot board the oven figure is within 0.16 s of *dt* for the eight layouts at the top of the
list (0.86 to 1.4 times *dt*), and the range is within 20 % of *dt*. Chips with the same layout and nearly the same result share a row. Timing is one input to
tombstoning. The wetting forces (pad geometry, paste volume and offset, placement, alloy) are not
modelled; the toe of the pads beyond the body and the ratio of the pad areas are listed beside the
timing, with the copper leaving each pad, whether it sits in a pour, and vias in or near it.

**What the models leave out.** Latent heat, flux, shadowing by tall parts, the non-uniform film
coefficient of real ovens, infrared heating, the temperature inside a part (a part is as warm as the
board under it, so the centre of a BGA cannot lag its corners), panels, and the second pass of a
double-sided board. The oven profile is generic. A quarter of the bodies on a typical board are an
assumed block over the courtyard; the report counts them.

**Solver cases.** `python3 -m openpcbsi validate` runs closed-form cases: step response and ramp lag
of a uniform board, the lag time, the two-resolution window against the lumped value; for the
measurements, the angle of drawn wedges, the width of a drawn slot and the gap between two nets
drawn at five angles. They show that the code solves its equations, not that the equations describe
an oven.

### Giving a ranking a limit, and scoring it

```toml
[dfm.judge]
source = "tombstone count of lot 12, assembler X"   # required: where the limits come from
joint_dt_s = 0.5          # chips whose joints melt further apart warn
reflow_cold_c = 8.0       # parts peaking more than this below the median warn
acute_deg = 60.0          # corners sharper than this ...
acute_track_mm = 0.15     # ... on a track this thin or thinner warn
```

Findings from these carry "board config: <source>" as their basis and then count in the exit code.

`--observed defects.toml` scores the rankings against defects seen on built boards:

```toml
[[defect]]
ref = "C12"
kind = "tombstone"
```

For each ranking it reports the share of (defective, sound) pairs it puts in the right order: 0.5 is
chance, 1 is perfect. The columns scored are the time between the joints, the ratio of copper
leaving the pads, a plane via at one pad only, the toe of the pads, and the reflow ranking. That is
the evidence a limit in `[dfm.judge]` should rest on.

## Not checked

- What the fabricator's CAM does to the files: etch compensation, changes to mask openings,
  panelisation.
- The etched cross-section and the finished track width (the published tolerance is recorded).
- Drill aspect ratio, copper balance and warpage against a limit (none is published; the figures are
  listed).
- Part height, shadowing in reflow, the second pass of a double-sided board.
- A land pattern against its part's data sheet (only two-pad chips, against their EIA body).
- Stencil aperture shapes and paste volume beyond the IPC-7525 ratios.
- Silkscreen, test access, panel rails and tooling holes.

KiCad's DRC still checks clearance, annular rings and holes against the board's own rules;
`rule_audit` says whether those rules are as strict as the fabricator's.

## Fabricator profile

A TOML file with three kinds of figure. [openpcbsi/fabs/jlcpcb.toml](../openpcbsi/fabs/jlcpcb.toml)
is the exemplar, from JLCPCB's capability, assembly, stencil and part-spacing pages on the date in
the file. `--fab NAME` reads `openpcbsi/fabs/NAME.toml`; `--fab FILE` reads a file.

| Section | Holds | Used by |
|---------|-------|---------|
| `[published]` | The fabricator's own limits. | Rules. |
| `[standard]` | Figures of an industry standard (IPC-7525 stencil ratios). | Rules. |
| `[assumed]` | Typical figures neither states: film coefficient, etch factor, oven type. | Rankings only. No rule reads them. |
| `[notes]` | Where the pages are ambiguous and what was assumed about them. | The report. |

The report lists every limit with its origin and whether a check read it; a published limit that
no check reads is listed in `limits_not_read_by_a_check_run` (acceptance test 39 requires that list
to be empty). Two published figures are recorded and not judged, and marked so: the finished track
width tolerance, and the copper clearance from V-cut edges (the board file does not say which edges
are V-cut).

## Configuration

Every key has a default.

```toml
[dfm]
fab = "jlcpcb"                 # profile name or file
mask_color = "green"           # black and white need a wider pad spacing at some fabricators
assembly_service = "standard"  # which column of the assembler's table applies
via_cover = ""                 # tented, plugged or filled; "" = the fabricator's default for the layer count
stencil_mm = 0.0               # 0 = the thinnest standard stencil of the profile
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
| `dfm.html` | Verdict of the rules, the rule and ranking tables, what is not checked, findings with the origin of each limit, the audit of KiCad's rules, the file comparison, the limits used, and each ranking. |
| `dfm.json` | Schema `openpcbsi-dfm/2`, below. |
| `fab_files/` | The Gerber and drill files that were plotted and compared (without `--gerbers`). |
| `layer_<layer>.png` | Findings on each copper layer. |
| `reflow_peak_top.png` | Modelled peak temperature of the top layer. |
| `<board>_openpcbsi.kicad_pcb` | The board with a marker per warning and failure on User.9. |

### `dfm.json`

| Field | Content |
|-------|---------|
| `schema`, `tool`, `version`, `date`, `seconds` | Schema id, `openpcbsi`, its version, time and duration of the run. |
| `board`, `board_sha256`, `kicad`, `config`, `config_hash` | The inputs. |
| `fab` | `name`, `retrieved`, `sources`, `profile`, `profile_hash`, `notes`, and `limits`: a list of `{key, value, source: published / standard / assumed / config, used, note}`. |
| `limits_not_read_by_a_check_run` | Published and standard limits that no check read (a full run only). |
| `exit_code`, `status`, `counts` | Verdict of the rules (`fail`, `warn`, `pass`) and totals outside waivers. |
| `checks` | Per check: `check`, `requirement`, `tier` (`rule` or `ranking`), `status`, counts, `note`, `seconds`. |
| `findings` | Every finding in the model of `evaluation.json`, with `tier`, `basis` (the origin of its limit) and `field` (BGA) added. |
| `rules` | Per rule: `params`, `metrics`, its own `findings`; `rule_audit` has `rows`. |
| `rankings` | Per ranking: `params`, `summary`, `rows` (every row), and `findings` when `[dfm.judge]` applies. |
| `not_checked`, `bga_fields`, `cases_pass`, `solver_versions` | What the review does not cover; the BGA parts; whether the closed-form cases pass. |
| `observed`, `judge` | With `--observed`: the score of each ranking; the `[dfm.judge]` table in force. |
| `waivers_file`, `waivers_hash`, `waivers`, `images`, `marker_board` | As in `evaluation.json`. |
