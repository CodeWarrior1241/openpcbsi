# openpcbsi

Signal- and power-integrity checks for KiCad 10 boards: the electrical review that KiCad's DRC
does not do. openpcbsi reads a `.kicad_pcb` and reports return-path problems, impedance,
crosstalk, channel loss, supply drop, PDN impedance and board temperature, with every finding
placed on the board. It can also compare two layouts of the same design, and check a layout against
a fabricator's and an assembler's capabilities before the files are sent.

openpcbsi is research software. Its solvers are checked against closed-form reference cases, and
its signal-integrity model has been checked against a real, fabricated board: Antmicro's open
Signal Integrity Test Board, whose VNA measurements are public. On that board the model's stub
resonances land within 5 % of the measured ones (see
[Correlation with a measured board](#correlation-with-a-measured-board)). That is one board, one
frequency range and one kind of structure; the rest of the tool is not yet correlated with
measurements.

## What it checks

- **Return paths**: reference-plane continuity under each trace, crossings of plane gaps and
  splits, return vias at layer changes, stitching along fast routes and the board edge, wide pads
  over solid reference planes.
- **Routing**: differential-pair coupling and skew, spacing and parallel length to other nets,
  distance to plane and board edges, via stubs and their resonance (blind, buried and back-drilled
  vias included), dangling ends, breakout necks inside BGA and connector fields.
- **Stitching**: a board-wide map of the distance to the nearest via joining the ground planes, with
  its percentiles and the largest via-free areas; return paths on boards without planes.
- **Impedance**: a 2D field solver for every trace width and pair geometry on every layer, with
  dielectric and conductor loss, copper roughness, etch and solder mask, against the targets of
  each signal class, with manufacturing corners (etch, dielectric thickness, permittivity).
- **Crosstalk**: near- and far-end coupling for every victim and aggressor within range.
- **Signal integrity**: per-net SPICE and Touchstone models (lines, vias, pads), point-to-point,
  T and fly-by nets with one result per receiver, step response with a linear driver, channel loss
  and eye estimates, propagation delay and skew per matching group.
- **3D regions**: openEMS models of discontinuities (layer changes, AC-coupling capacitors, splits,
  connector launches, breakouts) cut from the board with transmission-line ports, mixed-mode
  S-parameters and TDR, each with its standard fix as a variant.
- **Supply drop**: DC drop from each regulator or current-sense shunt to the load pins, the
  ground return included, against the rail's tolerance window; where the drop is spent along the
  current path; track heating and current per via.
- **PDN impedance**: |Z(f)| at the loads against a target from ripple and step current inside each
  rail's band, an optional package and die model, the droop for a load step, plane resonances, and
  a decoupling-capacitor optimiser.
- **Thermal**: a 2.5D board temperature map from part dissipation and I²R heating.
- **EMC**: a ranked list of likely radiators, for review.
- **KiCad DRC**: KiCad's own design-rule check and schematic parity, run as part of the evaluation.
- **Manufacturing** (`dfm`, a separate command): against a fabricator profile, acute copper corners
  that hold etchant (acid traps), gaps below the same-net spacing, track width, solder-mask webs and
  copper inside or beside mask openings; for assembly, a reflow simulation of the board and, for
  each two-terminal chip, how far apart in time its two joints melt (tombstoning).

## Requirements

- KiCad 10, whose system `python3` provides the `pcbnew` module, and `kicad-cli` on the `PATH`.
- numpy, scipy and matplotlib for that same `python3`.
- Optional: [KiPI](https://github.com/drewsum/KiPI), a finite-element DC and PDN solver that
  openpcbsi can use in place of its own. `deps/setup_kipi.sh` installs the pinned version into
  `deps/` in its own environment; see [deps/README.md](deps/README.md).
- Optional: [openEMS](https://openems.de) for full-wave runs on small board regions. Point
  `OPENPCBSI_OPENEMS_PYTHON` at the Python of an environment that has openEMS and CSXCAD.
- Optional: [kicad-python](https://pypi.org/project/kicad-python/) (`kipy`) in its own environment,
  for `evaluate --live`. Point `OPENPCBSI_KIPY_PYTHON` at that environment's Python.

Run openpcbsi from the repository root, or install it with `pip install .` (into KiCad's Python, or a
venv made with `--system-site-packages` so that `pcbnew` is found) to get an `openpcbsi` command.
Caches and results go to `cache/` and `runs/` beside a checkout, or under `./openpcbsi/` for an
installed copy; `OPENPCBSI_CACHE`, `OPENPCBSI_RUNS` and `OPENPCBSI_DEPS` set them. The board is only
read, never modified.

openpcbsi reads KiCad 10.0 board files. A board saved by a newer KiCad is converted when the
`kicad-cli` on the `PATH` can downgrade it; otherwise openpcbsi stops with a message naming the file
version, and the board has to be saved as KiCad 10.0 first.

## Quick start

```sh
git clone https://github.com/CodeWarrior1241/openpcbsi.git
cd openpcbsi
python3 -m openpcbsi evaluate path/to/board.kicad_pcb
```

1. **Evaluate with no configuration.** openpcbsi detects the board's signal classes, BGA and
   connector fields, and power rails from the netlist and footprints. It runs every check and
   KiCad's DRC, and writes a scorecard to `runs/<board>_<hash>/evaluation.html`. The page lists
   everything detection had to assume.

2. **Draft a configuration and review it.**

   ```sh
   python3 -m openpcbsi init path/to/board.kicad_pcb      # writes boards/<board>.toml
   ```

   The draft holds what was detected, with the assumptions as comments. Fill in what the board
   cannot tell you:
   - edge rates of the signal classes;
   - rail currents (until you do, a rail is marked `current_assumed` and its DC drop is never
     judged worse than a warning);
   - tolerances and ripple targets per rail;
   - part dissipations for the thermal check.

   [boards/au7p-fmc.toml](boards/au7p-fmc.toml) is a complete, commented example for an 8-layer
   FPGA board.

3. **Evaluate against the reviewed configuration**, adding waivers for findings you have accepted:

   ```sh
   python3 -m openpcbsi evaluate board.kicad_pcb -c boards/board.toml --waivers waivers.toml
   ```

4. **Compare two layouts** of the same design:

   ```sh
   python3 -m openpcbsi compare a.kicad_pcb b.kicad_pcb -c boards/board.toml
   ```

`evaluate` exits with 0 when nothing fails or warns, 1 for warnings only, 2 when anything fails and
3 when the board could not be evaluated, so it can gate a CI job.

### During layout

openpcbsi can also be called while a board is being laid out (by a person, a script, or a layout
tool). It detects the layout stage: before routing is complete, connections that are not routed yet
are reported as `pending`, not as failures. Each finding has a stable ID, so successive runs can be
compared, and findings on geometric rules carry a first fix.

```sh
python3 -m openpcbsi evaluate board.kicad_pcb --families g1,rp -o runs/now        # rule checks: seconds
python3 -m openpcbsi evaluate board.kicad_pcb --since runs/before/evaluation.json  # new / resolved findings
python3 -m openpcbsi plan board.kicad_pcb --layer In2.Cu --zdiff 100 --gap 0.15     # width for a target
```

A tool can also call `openpcbsi.api.evaluate(snapshot, config, families, changed)` on a board
snapshot, without writing files. Solver results are cached by their inputs, so an unchanged geometry
is never solved twice.

The board open in KiCad can be evaluated as it stands, unsaved edits included, through KiCad's IPC
API:

```sh
python3 -m openpcbsi evaluate --live --families g1,rp
```

This needs KiCad running with its API server enabled (Preferences > Plugins > Enable IPC API server)
and a Python that has [kicad-python](https://pypi.org/project/kicad-python/) (`kipy`), named by
`OPENPCBSI_KIPY_PYTHON`. openpcbsi runs `kipy` in that Python, never in its own. A live board gives
the same findings as the same board saved, and acceptance test #32 checks this against a headless
KiCad.

### Sign-off

```sh
python3 -m openpcbsi signoff board.kicad_pcb -c boards/board.toml --waivers waivers.toml --reviewer NAME
python3 -m openpcbsi verify-signoff runs/<board>_<hash>/signoff.json
```

`signoff.json` ties the verdict to the SHA-256 of the board, its KiCad project files, the
configuration and the waivers; `verify-signoff` fails if any of them changed since.

### Manufacturing sign-off

```sh
python3 -m openpcbsi dfm board.kicad_pcb [-c boards/board.toml] [--fab jlcpcb] [--waivers waivers.toml]
```

`dfm` checks the layout against a fabricator profile and writes `runs/<board>_<hash>/dfm/dfm.html`,
with the same finding model, waiver format and exit codes as `evaluate`. The profile shipped,
[openpcbsi/fabs/jlcpcb.toml](openpcbsi/fabs/jlcpcb.toml), holds JLCPCB's published capabilities as
an exemplar; copy it for another fabricator and pass the file with `--fab`. Every limit in the
report is marked as published by the fabricator, assumed, or set in the board configuration.

| Check | What it finds |
|-------|---------------|
| `fab_track` | Tracks narrower than the fabricator's minimum for the layer count and copper weight. |
| `acid_trap` | Corners in copper of one net sharper than the angle below which liquid is held, with the depth of the pocket and its flush time. |
| `same_net_gap` | Slots and openings in copper of one net narrower than the fabricator's same-net spacing. |
| `mask_web` | Pads too close for a solder-mask web, webs below the minimum, openings that merge. |
| `mask_exposed` | Copper of another net inside a pad's mask opening, or too close to it. |
| `pad_balance` | Chips whose two pads differ in the copper leaving them, in a pour, or in a via. |
| `reflow` | Peak temperature, time above liquidus, ramp and lag of every part through a convection oven. |
| `tombstone` | The time between the two joints of a chip reaching liquidus. |

The etching and reflow models are checked against closed forms only, not against fabricated boards
or an assembler's defect data, and several limits are proposals to be calibrated; the scorecard says
which inputs were assumed. Clearance, annular rings, copper slivers and starved thermals stay with
KiCad's DRC, which `evaluate` runs. Models, limits and file formats: [docs/dfm.md](docs/dfm.md).

## Commands

| Command | What it does |
|---------|--------------|
| `evaluate BOARD [-c CFG] [--waivers FILE] [--families ...] [--checks ...] [--stage S] [--since EVAL] [--em none\|plan\|run] [--json OUT] [-o DIR] [--no-drc] [--no-markers]` | Every check plus KiCad DRC, a scorecard and an exit code. `--em plan` lists automatic 3D regions, `--em run` runs them (minutes each). |
| `signoff BOARD [evaluate options] [--reviewer NAME]` | Evaluate and write `signoff.json`, tied to the exact inputs. |
| `verify-signoff RECORD` | Check that a sign-off record still matches its board and files. |
| `dfm BOARD [-c CFG] [--fab NAME\|FILE] [--checks track,etch,mask,balance,reflow,tombstone] [--waivers FILE] [--json OUT] [-o DIR] [--no-markers]` | Manufacturing sign-off against a fabricator profile: etching, solder mask and reflow checks, a scorecard and an exit code. |
| `plan BOARD --layer L (--zdiff Z [--gap G] \| --z0 Z)` | Track width for a target impedance on a layer of the board's stackup. |
| `init BOARD [-o FILE] [--force]` | Draft a board configuration from the board. |
| `check BOARD [-c CFG] [--checks ...] [--rules ...] [-o DIR] [--no-markers]` | Run checks and write the report, without scorecard or DRC. |
| `compare A B -c CFG [--config-b CFG] [--checks ...]` | Run both boards with the same configuration and report what changed. |
| `rediff RUN_A RUN_B` | Compare two finished runs again without re-running them. |
| `whatif BOARD -c CFG SPEC` | Sweep a single line's geometry, driver and termination before routing it. |
| `region BOARD [-c CFG] SPEC` | Full-wave openEMS run on a region cut from the board (see below). |
| `validate [--3d]` | Run the solvers' reference cases; `--3d` adds the openEMS cases (minutes each). |
| `correlate [--fetch]` | Compare openpcbsi with a measured, published board (see [Correlation](#correlation-with-a-measured-board)). |
| `snapshot BOARD` | Re-read the board into the cache. |

All commands are run as `python3 -m openpcbsi <command>`. Without `-c`, `evaluate`, `check`,
`whatif` and `region` use the configuration detected from the board.

`--checks` takes a comma-separated list of `drc` (electrical rules), `impedance`, `crosstalk`,
`si`, `vias`, `timing`, `dc`, `pdn`, `thermal` and `emc`. `--families` selects by group instead: `g1`
and `rp` (rules and return paths), `sd` (supply delivery), `si`, `th` (thermal), `emc` and `ev`
(KiCad DRC).

A 3D region spec names the nets and a box; ports are placed where the nets cross the box edge:

```toml
[region]
name = "rx0_caps"
nets = ["RX0_P", "RX0_N", "RX0C_P", "RX0C_N"]
box = [118.0, 62.0, 126.0, 70.0]      # mm
fmax_ghz = 20
```

A region runs only after a pre-flight check of the model, on two grids, and is accepted only when
both runs settle and agree. The grids must agree on the shape of the TDR (each normalised to its own
port reading) within 2 % and on the line impedance within 3 %. The absolute level of a 3D run moves with
the mesh by up to about 2.5 % (on one pair, three grids gave 98.95, 101.09 and 99.9 ohm against 99.6 ohm
in 2D), so that level is checked against the 2D solver instead: the model's line is compared with 2D
for the same geometry, and the TDR is renormalised to it when they differ by more than 2 %.

Parts between region nets: a resistor is a lumped resistor; a capacitor too large to charge within the
run (any AC-coupling capacitor) is modelled as its metal body bridging the pads, since in band it is a
short (its reactance at the lowest reported frequency is recorded, and reported as a problem when it is
not small); a small capacitor or an inductor is a series R-L-C lumped element.

## Outputs

Each run writes a directory under `runs/` (or `-o DIR`):

| File | Content |
|------|---------|
| `evaluation.html`, `evaluation.json` | Scorecard per check family, every finding with a stable ID, value, limit and position, and waiver status. Format: [docs/evaluation.md](docs/evaluation.md). |
| `report.html`, `report.json` | Every metric and finding of every check. |
| `layer_<layer>.png` | Findings drawn on each copper layer. |
| `pdn_<rail>.png`, `thermal_top.png` | PDN impedance per rail; temperature map. |
| `nets/` | A SPICE subcircuit and a Touchstone file per extracted net. |
| `kicad_drc.json` | KiCad's DRC report. |
| `progress.json`, `report.partial.json` | What has run so far, written before each check. |
| `signoff.json` | The sign-off record (`signoff` only). |
| `em_jobs/` | The automatic 3D region specs (`--em plan` / `run`). |
| `<board>_openpcbsi.kicad_pcb` | A copy of the board with errors and warnings as markers on layer User.9. |

`compare` writes `compare.html` and `compare.json` to `runs/compare_<A>_vs_<B>/`. `dfm` writes
`dfm.html`, `dfm.json`, its own layer images, `reflow_peak_top.png` and a marker board to the `dfm/`
folder of the run ([docs/dfm.md](docs/dfm.md)).

## Configuration

A board configuration is a TOML file. It sets:
- signal classes (net patterns, edge rate, impedance target);
- exclusion areas;
- reference layers and the limits of each rule;
- power rails (net, voltage, tolerance, current, source, loads, ripple target);
- thermal parts and channels.

Every key has a default. To use KiPI for the DC and PDN solves:

```toml
[engines]
dc = "kipi"           # or "openpcbsi" (default)
pdn = "kipi"
cross_check = true    # also run the built-in solver and tabulate the difference
```

KiPI's PDN sweep is slow on large boards; leave `pdn` out of `--checks` for a quick run.

Without a configuration, openpcbsi uses the KiCad project itself: netclasses (pattern and explicit
assignments in the `.kicad_pro`), custom rules in the `.kicad_dru`, and footprint fields such as
`Current` (a load's current) and `Tj_max`. Limits are taken from, in this order: a `.kicad_dru` rule
covering the net, the board configuration, the interface standard the net's name maps to (PCIe, USB,
HDMI, MIPI, LVDS, Ethernet, SATA, multi-gigabit transceivers), and openpcbsi's defaults. Every
evaluation says which inputs were sourced and which were assumed.

## Correlation with a measured board

openpcbsi's models are compared with a board that was fabricated and measured by someone else, so
that the comparison does not depend on a measurement setup of our own.

### The reference data and where it comes from

| | |
|---|---|
| Board | [Antmicro Signal Integrity Test Board](https://github.com/antmicro/signal-integrity-test-board), rev 1.0.1: a 4-layer KiCad project (KiCad 6 format, which KiCad 10 reads), 0.185 mm microstrip on F.Cu over In1 through 0.12 mm of dielectric (εr 4.18), SMA edge connectors |
| Measurements | VNA S11 from 0.2 to 6 GHz (1000 points), published by Antmicro in the examples of their [gerber2ems](https://github.com/antmicro/gerber2ems) repository (`examples/stub_short/vna.csv`, `examples/stub_long/vna.csv`) |
| Licence | Apache-2.0, © Antmicro. openpcbsi does not copy these files into its repository |
| Pin | `deps/references.lock` names both repositories, the commits used (board `de0b592`, data `5001f1a`) and the SHA-256 of every file |

`python3 -m openpcbsi correlate --fetch` downloads the four files at those commits into the cache,
checks every hash, and runs the comparison. A file that does not match its hash stops the run.

### How the measured structures were matched to the board

Each gerber2ems example is a slice of the test board. Their part-position files and port positions
locate them on the KiCad board:

| Case | Net on the board | Structure | Measured from |
|---|---|---|---|
| `stub_short` | `Net-(J10-Pad1)` | 45.4 mm through-line J9 to J10 with an 8.35 mm stub ending at R4 | J9; J10 terminated in 50 Ω by the VNA's second port |
| `stub_long` | `Net-(J7-Pad1)` | 45.4 mm through-line J7 to J8 with a 20.9 mm stub ending at R3 | J7; J8 terminated in 50 Ω |

R3 and R4, 0 Ω resistors that would short the stubs to ground, carry a `DNP` field on the board, so
the measured stubs are open. That agrees with the measurement: near-total reflection where an open
stub is a quarter wavelength long. Antmicro's own simulation set-up for these slices (port 1
excited, port 2 a 50 Ω port) gives the termination.

### How openpcbsi is compared

openpcbsi models each case from the KiCad board alone, as it would any board. It extracts the
net's path from the measured connector, with the branch to R3 or R4 as a stub. The lines come from
its 2D field solver for the board's stackup (here Z0 52.3 Ω, εeff 3.16). The pads add their excess
capacitance against the reference copper actually under them. The far connector is a 50 Ω load.

The pass criterion was fixed before the first comparison. The quarter-wave reflection peak and the
half-wave null must fall within 5 % of the measured frequency, which corresponds to about 10 % in
εeff. These frequencies depend on the stubs' length and εeff, which the board model contains, and
barely on the SMA launches and cables, which it does not. The |S11| curve error is reported but not
judged. The measurement has a ripple with a period of about 260 MHz, from its cables, so both curves
are averaged over 300 MHz before the features are read.

### Results

| Case | Feature | Model | Measured | Error |
|---|---|---|---|---|
| `stub_short` | quarter-wave peak | 4.80 GHz | 5.03 GHz | −4.6 % |
| `stub_long` | quarter-wave peak | 1.99 GHz | 1.90 GHz | +4.6 % |
| `stub_long` | half-wave null | 4.22 GHz | 4.21 GHz | +0.3 % |

The smoothed |S11| curves differ by 1.8 dB rms up to 4 GHz. The sharp half-wave null is the most
reliable figure. The quarter-wave peaks are broad, and where their maximum is read depends on how
the cable ripple is averaged. With 50 MHz averaging instead of 300 MHz the `stub_long` peak is
within 0.9 %, but the `stub_short` peak, broad and near the top of the band, is 6.9 % off, which
fails.

`runs/correlation/antmicro_si.json` holds the full result and `antmicro_si.png` overlays the
curves. Acceptance test #8 runs this comparison.

### What the comparison changed

Running openpcbsi against this board found four faults, all fixed before the results above:
- **Do-not-populate fields.** The board marks unfitted parts with a `DNP` field (KiCad 6) rather
  than the attribute; openpcbsi read only the attribute and would have modelled the stubs as
  shorted.
- **T-junctions.** A stub branching part way along a track was not joined to it, so the first
  model had no stub at all.
- **Launch pads over voids.** The excess capacitance of the SMA pads was taken against solid In1,
  where the board has a void under them (about 2.8 pF instead of the real fraction), which made the
  whole line reflect.
- **Smoothing.** The first feature finder read the highest peak, not the first one.

### Limits

- One board, one fab, one measurement set, by a third party.
- 0.2–6 GHz, through SMA launches and cables that are not de-embedded.
- Outer-layer microstrip only: no stripline, vias or differential pairs.
- Signal integrity only. The DC drop, PDN, thermal and 3D results are checked against closed forms
  and other solvers, not against measurements.

## How it works

openpcbsi reads the board once through KiCad's `pcbnew` API into a cached snapshot. Boards saved by
a newer KiCad are first converted with `kicad-cli pcb downgrade`. Copper is rasterised per layer at
25 µm, and plane questions (reference copper, splits, antipads, edges) are answered on the rasters;
traces stay vector geometry.

Each check uses its own model:
- impedance and crosstalk: a 2D finite-volume field solver;
- DC drop: a resistive network of the copper;
- PDN impedance: a cavity and capacitor network;
- temperature: a 2.5D finite-volume model;
- reflow: the same model with heat capacity, stepped through the oven's zones, and a finer window
  around each chip for the lag of its two joints;
- etching: a ball rolled through the space between copper of one net, on a 6.25 µm raster drawn
  from the vector data, with a diffusion solve for each pocket it cannot reach.

Each finding carries a severity, a stable ID and a position, so waivers survive re-runs and two
layouts can be compared finding by finding.

## License

MIT, see [LICENSE](LICENSE). KiPI is not included; `deps/setup_kipi.sh` fetches it under its own
GPL-3.0 licence.
