# openlynx

Electrical checks for KiCad 10 boards of the kind a commercial SI/PI suite makes and KiCad's DRC
does not: return paths, impedance, crosstalk, channel loss, supply drop and noise, heat. It
implements `REQUIREMENTS.md`: G0 to G17, and the first step of the design evaluation in sections
5.10 to 5.16 (G18 to G21, RP1, SD1, EV1, EV5). The pilot is the AU7P FMC board.

The solvers pass their closed-form validation cases. None of it has been compared with a measured
board or with HyperLynx.

## Using it on a KiCad design

System `python3` of stock KiCad 10 (pcbnew, numpy, scipy, matplotlib) and `kicad-cli` on the PATH.
G8 also uses the openEMS venv in `Dev_Tools/OpenEMS/install/venv`; KiPI (G10/G11 engine) lives in
`deps/` (see `deps/README.md`). All caches and outputs stay in this folder (`cache/`, `runs/`). The
board is only read; it is never modified.

```
cd /media/fpgadev/Dev_Tools/Work/Open_Source_Hyperlynx
```

1. **First look, no config.** `evaluate` detects classes, fields and rails from the board itself
   (EV5) and runs every built check plus KiCad's DRC:

   ```
   python3 -m openlynx evaluate PATH/TO/BOARD.kicad_pcb
   ```

   The scorecard is printed and written to `runs/<board>_<sha8>/evaluation.html` and
   `evaluation.json`; the exit code is 0 (nothing fails or warns), 1 (warnings only), 2 (any fail)
   or 3 (could not evaluate). Everything detection assumed is listed at the top of the page.

2. **Draft a config and review it.**

   ```
   python3 -m openlynx init PATH/TO/BOARD.kicad_pcb        # writes boards/<board>.toml
   ```

   The draft has every detected class, field and rail, with the assumptions as comments. What the
   board cannot tell you has to be filled in by hand:
   - edge rates and Nyquist frequencies per class (detection uses 35 ps / 6.25 GHz for serial-lane
     names, 300 ps / 1 GHz for other pairs, 500 ps for clocks);
   - rail currents (detection assumes 0.5 A and marks the rail `current_assumed = true`: a DC
     result at an assumed current is never judged worse than `warn`);
   - tolerances, `ripple_pct` and `step_a` per rail (G11 runs only on rails that have them);
   - rails whose voltage is not in a net name or upstream net (listed under "Left out");
   - part dissipations for G13, G3 drivers and terminations, G7 channels.

   `boards/au7p-fmc.toml` is a reviewed config with every section commented.

3. **Evaluate against the reviewed config**, with waivers once findings have been reviewed:

   ```
   python3 -m openlynx evaluate BOARD.kicad_pcb -c boards/<board>.toml [--waivers waivers.toml]
   ```

   Findings carry stable IDs (G19), so a waiver keeps matching while its cause is unchanged, and
   one whose finding has gone is reported as stale. File formats: `docs/evaluation.md`.

4. **Compare two layouts** of the same design (G15): `compare A.kicad_pcb B.kicad_pcb -c CFG`.

### Commands

```
python3 -m openlynx evaluate BOARD [-c CFG] [--waivers FILE] [--families g1,rp,sd,si,th,emc,ev]
                                   [--checks ...] [--json OUT] [-o DIR] [--no-drc] [--no-markers]
python3 -m openlynx init     BOARD [-o boards/NAME.toml] [--force]          # EV5 draft config
python3 -m openlynx check    BOARD [-c CFG] [--checks drc,impedance] [--rules ref_plane,...] [-o DIR] [--no-markers]
python3 -m openlynx compare  A.kicad_pcb B.kicad_pcb -c CFG [--config-b CFG_B] [--checks ...]   # G15 A/B
python3 -m openlynx rediff   RUN_DIR_A RUN_DIR_B                            # A/B from two finished runs
python3 -m openlynx validate                                                # G17 cases
python3 -m openlynx whatif   BOARD -c CFG boards/whatif_se_fcu.toml         # G5
python3 -m openlynx region   BOARD -c CFG boards/regions/em1_fcu_line.toml  # G8 (openEMS)
python3 -m openlynx snapshot BOARD                                          # re-read the board
```

`check` runs the checks and writes the report without the scorecard, DRC or exit code. Without
`-c`, `check`, `whatif` and `region` use the detected config as `evaluate` does. `compare` needs a
config, so that both boards are judged by the same one.

Checks for `--checks`: `drc` (G1), `impedance` (G6), `dc` (G10), `pdn` (G11/G12), `crosstalk`
(G4), `thermal` (G13), `vias` (G9), `si` (G2/G3/G7), `emc` (G14). `--families` selects by the
families of REQUIREMENTS 5.11 to 5.15 instead: `g1` and `rp` run `drc`, `sd` runs `dc,pdn`, `si`
runs `impedance,crosstalk,si,vias`, `th` runs `thermal`, `emc` runs `pdn,emc`, `ev` runs KiCad's
DRC. `--rules` limits G1 to some rules (`ref_plane,return_via,...`).

A board with KiPI engines in its config (`[engines] dc/pdn = "kipi"`) spends most of its time in
KiPI's AC sweep (about 19 min for the pilot's 4 regulators); leave `pdn` out of `--checks` for a
quick run.

### Outputs

A run writes the following to `runs/<board>_<sha8>/` (or `-o DIR`):

- `evaluation.json`, `evaluation.html` (`evaluate` only): scorecard per family and check, every
  finding with ID, severity, value and limit, waivers and their status. Schema in `docs/evaluation.md`.
- `kicad_drc.json` (`evaluate` only): KiCad's own DRC report (EV1).
- `report.json`: every metric and finding, with positions.
- `report.html`: one page per run.
- `layer_<L>.png`: findings drawn on each copper layer.
- `pdn_<rail>.png`: Z(f) per rail.
- `thermal_top.png`: the temperature map.
- `nets/`: a SPICE subcircuit and a Touchstone file per net (characters other than letters, digits
  and `_.+-` in a net name, such as the `/` of hierarchical names, become `_`).
- `kipi/`: KiPI's setup, job and result, when KiPI ran.
- `<board>_openlynx.kicad_pcb`: a copy of the board with the errors and warnings as markers on
  User.9.

`compare` writes `runs/compare_<A>_vs_<B>/compare.html` and `compare.json`.

A board saved by a newer KiCad (Backplane) is converted first with
`kicad-cli pcb downgrade --target 10.0`; this is automatic when pcbnew cannot load the file, and the
report records the converted copy (`downgraded`).

## How it works

1. `snapshot.py` reads the board once through pcbnew and writes everything else needs to a gzip
   JSON keyed by the board's SHA-256: stackup (parsed from the file), outline, tracks, vias with
   flashed layers, zone fills, pads with polygons and pin functions, footprints, courtyards, net
   classes. It is the only module that imports pcbnew (apart from the marker board in `report.py`).
2. `model.py` builds numpy arrays and a 25 µm raster of net ids per copper layer (cached). Plane
   lookups (reference copper, antipads, splits, edges) are done on the rasters, distances by
   distance transforms. Track geometry stays vector.
3. Each check reads the board model and the per-board TOML config. Without a config file,
   `detect.py` builds one from the snapshot:
   - GND-like reference nets;
   - P/N pairs by name, grouped into serial lanes (by name), KiCad net classes, or other pairs,
     with Zdiff from a net class named like `90ohm`;
   - BGA and connector fields;
   - the supply tree, rooted at regulator outputs (pin functions where the footprint has them)
     and at power-entry connectors, through shunts, ferrites, inductors, fuses and pass FETs. A
     net fed through a series part inherits its supply's voltage.
4. `evaluate.py` turns the check results into G19 findings (`info` / `pass` / `warn` / `fail`,
   stable ID, value, limit, nets, references, position, evidence). It judges DC drop (SD1) and PDN
   impedance, runs KiCad's DRC (EV1), applies waivers (G20) and writes the scorecard (G18).

## Requirements and what was built

| ID | Status | Module | Notes |
|----|--------|--------|-------|
| G0 | done | `docs/G0_survey.md` | Survey by web search, 2026-11-14. Open tools have since appeared for G10, G11 and G13 (KiPI, padne, kicad-irdrop, KiCad_Thermal_Sim): see "Decisions open" below. |
| G1 | done | `drc.py`, `context.py` | Rules: reference-plane continuity (PC1), gap and split crossings (PC2), unreferenced segments, return vias (PC3, power-area references point at the nearest capacitor), via stubs with resonance (PC4), pads over planes (PC5, judged as RP1), pair coupling and skew (PC6), distance to other nets and parallel length (PC7), plane and board edges (PC8), stitching along corridors and the edge (PC9), decoupling via distance and distance to load (PC11), thermal vias (PC13), power areas overlapping across thin dielectric, dangling ends and T-branches. Exclusion areas from courtyards or rectangles; findings inside them are listed, not failed. |
| G2 | done (lumped vias, no package) | `si.py` | Point-to-point path per net or pair. Lines from G6 parameters, vias from G9, pad capacitance. Writes `.sp` (T lines + LC) and `.s2p`. |
| G3 | partial | `si.py` | Linear driver and termination, not IBIS: ngspice/KIBIS are not installed. Overshoot, ringback, monotonic edge, 50 % delay, suggested series R. |
| G4 | done | `crosstalk.py` | Every victim in scope against every aggressor within 0.6 mm. Kb/Kf from the 2D solver per layer, widths and spacing. NEXT saturating with length, FEXT linear. Ranked list; the top cases are listed for G8. |
| G5 | done (single line) | `si.py whatif` | TOML topology: layer, width, gap, length, driver, termination; cartesian sweeps. Multi-drop topologies are not yet supported. |
| G6 | done | `xsec.py`, `planner.py` | 2D finite volumes on a graded grid. Complex permittivity (Djordjevic-Sarkar), Wheeler conductor loss with Hammerstad roughness, etch trapezoid, conformal mask. Every width and pair geometry on every layer, against class targets and fabricator figures (`[xsec] fab`). |
| G7 | partial | `si.py channel` | Cascade of nets and series capacitors; masks: SDD11, ILD, IL at Nyquist; peak-distortion eye. Not yet: SCD21 (needs 4-port, G8), connector/package Touchstone in the cascade, IBIS-AMI, COM, a standards mask library. |
| G7a | not started | | Not needed for the pilot. |
| G8 | runs, not yet trustworthy | `region.py`, `em_runner.py` | Region cut from the snapshot, PEC copper sheets, via cylinders, lossy slabs, lumped ports at pads or points (nets under a port are clipped at the ports), two grid resolutions with a dB agreement test, cache by job hash, line impedance against G6. Pilot smoke test (`boards/regions/em1_fcu_line.toml`, 4 mm of 0.15 mm F.Cu): S21 -0.05 dB at 1 GHz, grids agree to 0.60 dB (limit 0.5: not accepted), extracted Z 80 ohm against 55 ohm from G6. The vertical lumped ports are not de-embedded, which dominates the extraction. Next: MSL/waveguide ports or port de-embedding (two line lengths), mixed mode, TDR. |
| G9 | done (closed form) | `vias.py` | Coaxial barrel model per via geometry, antipads measured on the plane rasters, nearest return via. Cached library `cache/via_library.json`. Not field-solved. |
| G10 | done: KiPI engine, openlynx cross-check | `kipi_bridge.py` + `deps/` (KiPI v1.0.1); `pdn.py dc_drop` | KiPI (FEM, `[engines] dc = "kipi"`) solves the whole power tree from the 12 V input: regulators, shunts and ferrites as series links, remote sense, and the ground net. The tree is derived from the rail config. Reported per rail: regulator-to-ball drop (shunts and ground included) against the window, PASS/FAIL, via currents, current density. openlynx's own raster solver runs beside it as a cross-check on a like-for-like figure (from the shunt or sense point, supply side): on the pilot the two agree within 0.64 mV on every rail, -0.4 to +13.7 % (board `d454fef3`; the larger percentages are on the 3-4 mV drops of the GTH rails), and within 1.2 mV, -5.8 to +8.7 % (board `e6b71042`). |
| G11 | done: KiPI engine, openlynx optimiser | `kipi_bridge.py`; `pdn.py pdn_impedance` | KiPI field-solves rail and ground per regulator and gives Z(f) at the load balls against V·ripple/ΔI up to `target_fmax_hz`. openlynx's own network (plane-pair cavity cells, track and via L, C-ESR-ESL capacitors) runs beside it for the cross-check and carries the greedy decoupling optimiser, which KiPI does not have. The report overlays both curves. KiPI AC is slow: about 19 min for the pilot's 4 regulators (tune `[kipi] ac_*`). |
| G12 | partial | `pdn.py` | Resonances (peaks of Z at the load, place of the largest response) and Z21 from the source within one rail. Coupling between different nets (e.g. +12V switcher to the GTH supplies) needs a multi-net cavity model, not built. |
| G13 | done | `thermal.py` | 2.5-D finite volumes, copper fraction per cell and layer, vias, convection both sides, I²R from G10. Tj = T_board + P θjb. Runs only when `[thermal.parts]` lists dissipations. |
| G14 | rule-based | `emc.py` | Ranked radiators from G12 resonances near the edge, PC9 edge stitching, pair skew/uncoupled length, PC2 crossings. Rule of thumb for review. |
| G15 | done | `compare.py` | Same config for both runs (a run is refused if solver sections differ); nets by name, then by pad set; better/worse per metric; overlaid Z(f); nets only in one board. |
| G16 | done | `report.py` | JSON, HTML, PNG per layer, marker board on User.9. |
| G17 | done for the solvers built | `validate.py`, `pdn.py`, `thermal.py` | 2D: Hammerstad-Jensen microstrip, Cohn stripline and coupled stripline (exact), homogeneous dielectric loss (exact), wide-strip conductor loss (exact limit), Pozar formula (informative). DC sheet resistance, capacitor SRF, cavity f10, thermal energy balance. The openEMS-against-G6 case (EM1-like) does not pass yet (see G8). A check whose solver has no passed case reports "unvalidated". |
| G18 | done for the built checks | `evaluate.py` | Scorecard per family and check: `fail` / `warn` / `pass`, `unvalidated` (solver without a passed G17 case), `assumed inputs` (DC at assumed currents), `not run` (inputs missing, e.g. no rail ripple targets), `not built`. A family that has unbuilt checks is marked `partial`. |
| G19 | done | `evaluate.py` | Severities `info` / `pass` / `warn` / `fail`. Stable ID: check plus a hash of nets, references and position rounded to 0.1 mm. Value, limit, unit, nets, references, layer, position, reason, evidence (layer image, plot, KiCad report). Findings in exclusion areas are counted separately. "Best achievable" is filled only where a check computes one (none yet: it comes with PL1). |
| G20 | done | `evaluate.py` | TOML waivers by finding ID, or by check with net / reference / area globs; reason and reviewer required. Waived findings stay in the report as waived; unused waivers are reported as stale. |
| G21 | done, without `--claims` and `--em` | `__main__.py`, `evaluate.py` | `evaluate BOARD [-c] [--waivers] [--families] [--json]`, exit codes 0 / 1 / 2 / 3, `[evaluate] warnings_exit` for CI. The claims file (EV3) and region jobs (FS9) are not built, so their options are not offered. |
| G22 | partial | `snapshot.py`, `evaluate.py` | KiCad 10 files read directly; newer files downgraded automatically and recorded in the report. The schematic is used for DRC parity when it sits next to the board; the project's `.kicad_pro` / `.kicad_dru` are used by KiCad's DRC, not yet by openlynx's own rules. |
| G23 | partial | `docs/evaluation.md` | `evaluation.json` is versioned (`openlynx-evaluation/1`) and documented. Marker output is the board copy only; plain JSON marker shapes for a viewer are not built. |
| G24 | not started | | Library API on a snapshot, kipy adapter, incremental runs. |
| G25 | done | `evaluate.py` | The evaluation records openlynx version, config hash, waiver hash, solver versions, KiPI pin and the KiCad DRC version. Caches are keyed by content hash. |
| PL1-PL10 | not started | | Placement family. Shown as `not built` in the scorecard. |
| SD1 | done | `evaluate.py` | Drop as a share of the window from either engine: `fail` above 100 %, `warn` above `[evaluate] dc_warn_pct` (50 %); engines disagreeing by more than 20 % is a `warn`. |
| SD2-SD9 | not started | | BGA supply entry, hot copper, current per via, shunt isolation, Kelvin sense, sense side, PDN target band, necks. |
| RP1 | done | `drc.py pad_plane` | On classes with Nyquist ≥ `fast_ghz` (3 GHz): solid reference under a wide pad is an error, partial a warning, cut-out info. Not yet: the `h_mm` depth limit and the void margin. |
| RP2-RP9 | as G1 | `drc.py` | The existing rules, reported under the RP family. RP2's per-class default limit, RP3's percentile map and RP5's cavity resonance are not built. |
| FS1-FS10 | not started | | FS1's 2D cases exist as G17; the gates themselves and the automatic regions do not. |
| EV1 | done | `evaluate.py` | `kicad-cli pcb drc --severity-all --refill-zones` (`--schematic-parity` when the schematic is next to the board). KiCad errors, unconnected items and parity issues are `fail`. KiCad warnings are `warn` by default (`[evaluate] drc_warning = "fail"` follows EV1 literally): on KiCad's own CM5 demo 172 of 176 violations are warnings, 112 of them silkscreen and text. |
| EV2-EV4 | not started | | DRC on the downgraded copy, claims file, A/B scorecard. |
| EV5 | done | `detect.py` | Fields, pairs, classes, shunts, regulators and rails from the netlist and footprints. Tested on the pilot (every rail source matches the reviewed config) and on KiCad's `cm5_minima` demo. |

## Limits to know

- With KiPI (default for the pilot) the ground net is solved; openlynx's own G10/G11 solvers still
  assume ideal ground. The PDN sees the board up to the balls (no package, no die).
- KiPI's input chain, regulators and ferrite DCRs come from the board config (`[kipi]`,
  `[kipi.input]`, `regulator =` per rail); the FET and fuse resistances there are estimates.
- The SI is linear-driver SI. Treat absolute overshoot and eye figures as estimates until IBIS
  (ngspice + KIBIS) is available; A/B differences are more trustworthy than absolute values.
- Copper roughness is a single RMS figure (config); the pilot's 0.4 µm is an assumption.
- Currents, dissipations, ripple and step currents in `boards/au7p-fmc.toml` are estimates, not
  Power Design Manager figures. The DC drop and PDN results scale with them.
- Detection (no config) is heuristic: net names and part values decide classes and the supply
  tree. Known gaps:
  - A rail named without its voltage (VCCO_65 from an adjustable regulator) is left out.
  - A module or connector that supplies a rail is not recognised as a source.
  - Serial-lane names include slower interfaces (Ethernet, MIPI), which get the 35 ps default.

  Its results are a first look; judge a design on a reviewed config.
- The `edge_distance` rule (PC8) is slow on large boards: 78 s of a 5 min evaluation of the CM5
  demo.

## Decisions open

- KiPI is wrapped for G10/G11 (research use, licence not a concern); see `deps/README.md` for the
  pin and the linkage. Still open: whether KiPI's electro-thermal mode should replace openlynx's G13,
  and whether padne is worth adding as a second independent DC check.
- Package location (REQUIREMENTS section 7): openlynx lives here as one package, run from this
  folder against any board path; board configs are in `boards/`. REQUIREMENTS still lists this as
  open.
