# openpcbsi

Signal- and power-integrity checks for KiCad 10 boards: the electrical review that KiCad's DRC
does not do. openpcbsi reads a `.kicad_pcb` and reports return-path problems, impedance,
crosstalk, channel loss, supply drop, PDN impedance and board temperature, with every finding
placed on the board. It can also compare two layouts of the same design.

openpcbsi is research software. Its solvers are checked against closed-form reference cases, but
not yet against measured boards.

## What it checks

- **Return paths**: reference-plane continuity under each trace, crossings of plane gaps and
  splits, return vias at layer changes, stitching along fast routes and the board edge, wide pads
  over solid reference planes.
- **Routing**: differential-pair coupling and skew, spacing and parallel length to other nets,
  distance to plane and board edges, via stubs and their resonance, dangling ends.
- **Impedance**: a 2D field solver for every trace width and pair geometry on every layer, with
  dielectric and conductor loss, copper roughness, etch and solder mask, against the targets of
  each signal class.
- **Crosstalk**: near- and far-end coupling for every victim and aggressor within range.
- **Signal integrity**: per-net SPICE and Touchstone models (lines, vias, pads), step response
  with a linear driver, channel loss and eye estimates.
- **Supply drop**: DC drop from each regulator or current-sense shunt to the load pins, against
  the rail's tolerance window, with via currents.
- **PDN impedance**: |Z(f)| at the loads against a target from ripple and step current, plane
  resonances, and a decoupling-capacitor optimiser.
- **Thermal**: a 2.5D board temperature map from part dissipation and I²R heating.
- **EMC**: a ranked list of likely radiators, for review.
- **KiCad DRC**: KiCad's own design-rule check and schematic parity, run as part of the evaluation.

## Requirements

- KiCad 10, whose system `python3` provides the `pcbnew` module, and `kicad-cli` on the `PATH`.
- numpy, scipy and matplotlib for that same `python3`.
- Optional: [KiPI](https://github.com/drewsum/KiPI), a finite-element DC and PDN solver that
  openpcbsi can use in place of its own. `deps/setup_kipi.sh` installs the pinned version into
  `deps/` in its own environment; see [deps/README.md](deps/README.md).
- Optional: [openEMS](https://openems.de) for full-wave runs on small board regions. Point
  `OPENPCBSI_OPENEMS_PYTHON` at the Python of an environment that has openEMS and CSXCAD.

Run openpcbsi from the repository root. Caches and results go to `cache/` and `runs/` there. The
board is only read, never modified.

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

## Commands

| Command | What it does |
|---------|--------------|
| `evaluate BOARD [-c CFG] [--waivers FILE] [--families ...] [--checks ...] [--json OUT] [-o DIR] [--no-drc] [--no-markers]` | Every check plus KiCad DRC, a scorecard and an exit code. |
| `init BOARD [-o FILE] [--force]` | Draft a board configuration from the board. |
| `check BOARD [-c CFG] [--checks ...] [--rules ...] [-o DIR] [--no-markers]` | Run checks and write the report, without scorecard or DRC. |
| `compare A B -c CFG [--config-b CFG] [--checks ...]` | Run both boards with the same configuration and report what changed. |
| `rediff RUN_A RUN_B` | Compare two finished runs again without re-running them. |
| `whatif BOARD -c CFG SPEC` | Sweep a single line's geometry, driver and termination before routing it. |
| `region BOARD -c CFG SPEC` | Full-wave openEMS run on a region cut from the board. |
| `validate` | Run the solvers' reference cases. |
| `snapshot BOARD` | Re-read the board into the cache. |

All commands are run as `python3 -m openpcbsi <command>`. Without `-c`, `evaluate`, `check`,
`whatif` and `region` use the configuration detected from the board.

`--checks` takes a comma-separated list of `drc` (electrical rules), `impedance`, `crosstalk`,
`si`, `vias`, `dc`, `pdn`, `thermal` and `emc`. `--families` selects by group instead: `g1` and
`rp` (rules and return paths), `sd` (supply delivery), `si`, `th` (thermal), `emc` and `ev`
(KiCad DRC).

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
| `<board>_openpcbsi.kicad_pcb` | A copy of the board with errors and warnings as markers on layer User.9. |

`compare` writes `compare.html` and `compare.json` to `runs/compare_<A>_vs_<B>/`.

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

## How it works

openpcbsi reads the board once through KiCad's `pcbnew` API into a cached snapshot. Boards saved by
a newer KiCad are first converted with `kicad-cli pcb downgrade`. Copper is rasterised per layer at
25 µm, and plane questions (reference copper, splits, antipads, edges) are answered on the rasters;
traces stay vector geometry.

Each check uses its own model:
- impedance and crosstalk: a 2D finite-volume field solver;
- DC drop: a resistive network of the copper;
- PDN impedance: a cavity and capacitor network;
- temperature: a 2.5D finite-volume model.

Each finding carries a severity, a stable ID and a position, so waivers survive re-runs and two
layouts can be compared finding by finding.

## License

MIT, see [LICENSE](LICENSE). KiPI is not included; `deps/setup_kipi.sh` fetches it under its own
GPL-3.0 licence.
