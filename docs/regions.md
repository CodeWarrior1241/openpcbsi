# 3D regions (openEMS)

A 3D region is a full-wave openEMS model of a small box cut from the board: a pair's layer change, a
connector launch, AC-coupling capacitors. `evaluate --em plan` lists the regions openpcbsi would model,
`evaluate --em run` solves them, and `region BOARD SPEC` runs one region you describe yourself. Each takes
minutes to hours. openEMS is found through `OPENPCBSI_OPENEMS_PYTHON` (see the README's Requirements).

## The region spec

A region spec names the nets and a box; ports are placed where the nets cross the box edge:

```toml
[region]
name = "rx0_caps"
nets = ["RX0_P", "RX0_N", "RX0C_P", "RX0C_N"]
box = [118.0, 62.0, 126.0, 70.0]      # mm
fmax_ghz = 20
```

Optional keys:

| Key | Default | Meaning |
|---|---|---|
| `rotate_deg`, `rotate_about` | 0, the box's centre | Turn the board before the box is applied, so that a 45-degree route runs along an axis and a side cuts it square. The planner sets `rotate_deg = -45` where the box's sides cut its nets at 45 degrees. |
| `include_mm` | 1.0 (at least 10 field heights) | How far from the region's own copper other nets' copper is kept in the model. Farther copper is left out, as the absorbing sides already stand for the board beyond them; 0 keeps everything. |
| `term_ohm`, `terminate` | 50, true | The load at a pad where a region net ends inside the box; `terminate = false` removes it. |
| `max_hours` | 6 (or `OPENPCBSI_EM_MAX_HOURS`) | Cap on the region's total run time. |

## Before solving: the pre-flight

The pre-flight refuses a model that would give no usable answer:

- a lane that does not run through the box to a port on another side or to a termination (an open end,
  or a break); routes are followed segment by segment, through vias, pads and series parts;
- a pair whose legs leave the box on different sides (a box corner cuts it), or a port whose track crosses
  the box edge more than 10 degrees off square;
- copper that ends in or runs along the absorbing boundary;
- a run predicted to stop on the time-step limit before the excitation has passed, or predicted to take
  longer than `max_hours` at this machine's median speed over its last finished grids (350 Mcells/s until
  it has some). The time step follows the smallest mesh cell (no two lines are closer than 5 um), and each
  run's step limit is 4 times the steps its excitation needs. The worst case (every run to its step limit,
  at the slowest measured speed) is logged.

## While it runs

- **Instability.** openEMS's field energy is watched. A grid whose energy grows a hundredfold after it has
  decayed (an FDTD instability, which only grows) is stopped at once and solved again, once, at 0.9 of its
  cell size. An instability belongs to one mesh: on the example FPGA board one region's grid was unstable at
  full size and stable at 0.9, and in a larger box with either absorber thickness.
- **Time cap.** A region still running at `max_hours` is stopped there; time spent waiting for a job slot
  does not count. The fine grid is not started when the coarse grid's measured time says it would pass
  `max_hours`: the cap holds for the region's grids together, so the message gives the time used, the time
  the next grid needs and the `max_hours` that would let it finish. A region stopped has no answer, and
  says why.
- **Parallel runs and fix variants.** `evaluate --em run` solves the regions and their fix variants as
  many at a time as there are openEMS job slots (`OPENPCBSI_JOB_SLOTS`, default one per 16 cores). A fix
  variant is planned only where it adds something: return vias only where the board has none within 1 mm
  of the layer change, placed 0.6-1.0 mm away and clear of every other net's via, track and pad (an added
  via has a pad only on its own net's planes and a clearance hole in every other plane); a cut-out only
  where the board has none.

## What is modelled

- **Terminations.** A region net that ends on a part's pad inside the box (a connector pin, a BGA ball)
  gets a 50 ohm lumped port from the pad to the nearest reference plane under it: the receiver's or
  cable's load (100 ohm differential for a pair) and a measured second port, instead of an open end that
  reflects everything. Its height is part of that load: a pad over a cut-out reaches a plane further down.
- **Copper thickness.** The region's own copper on inner layers has its real thickness, as on the outer
  layers (a zero-thickness sheet read an inner-layer pair 8.9 % high against 2D); other nets' inner copper
  stays a sheet.
- **Parts between region nets.** A resistor is a lumped resistor. A capacitor too large to charge within
  the run (any AC-coupling capacitor) is modelled as its metal body bridging the pads, since in band it is
  a short; its reactance at the lowest reported frequency is recorded, and reported as a problem when it
  is not small. A small capacitor or an inductor is a series R-L-C lumped element.

## Reading the results

- **Reflection floor.** The port leads are zero-thickness strips (openEMS's line port measures current on
  a sheet), so their impedance differs from the board's thick track: about 7 % on a microstrip, 13 to 16 %
  on pairs. That step at each box edge, and the lead's own mismatch to the 50 ohm-per-leg reference, put a
  floor under every S11 / Sdd11: a uniform line shows -20 to -23 dB. Each result states this floor
  (`reflection_floor_db`); a return loss near or below it is the model's, not the layout's.
- **TDR excursion.** `tdr_excursion` is read inside the box, past the near lead's step and short of the
  far one, against the model's calibrated line or else the 2D impedance of the same geometry. It is what
  the planner's FS9 findings judge, at the link's own edge when its class gives a rate. The time gate keeps
  about 1 % of the near and far leads' steps (their Gaussian tails), which bounds what it can resolve.
- **Absolute level.** The level of a 3D run moves with the mesh by up to about 2.5 % (on one pair, three
  grids gave 98.95, 101.09 and 99.9 ohm against 99.6 ohm in 2D), so it is checked against the 2D solver:
  the model's line is compared with 2D for the same geometry when the box holds a uniform run long enough.

## Convergence: two grids

A region runs on two grids, compared on what the region reports, never on the whole trace or on each
grid's own port reading (the port lead's level moves by up to about 3 % with the mesh; late time is
multiple reflections between the model's ports):

- the TDR excursion inside the box's window, both grids against the same reference (the 2D impedance of
  the same geometry), and the same in-box reflections time-gated and seen with the link's own edge (the
  class's `rise_ps`, or 0.15 / `nyquist_ghz`: 30 % of a bit, stated as an assumption);
- the transmission at the link's Nyquist (`nyquist_ghz`), in dB.

They agree when every excursion matches within the larger of 1.5 percentage points and 15 % of its size,
and the transmission within the larger of 0.1 dB and 15 %. When they do not, the region still stands if
both grids put it on the same side of the 10 % TDR limit (the spread is reported beside it); when they
straddle the limit, a third, finer grid decides, if `max_hours` allows. A region whose grids both stay
within 5 % and whose return loss is under the model's own floor is a pass: clean, below the model's
resolution.

**Runs cut on the step limit.** A run that stops on its step limit (a plane fragment that keeps ringing,
which the board's decoupling would damp) is accepted when every port's voltage over the last tenth of the
run is at least 50 dB below the largest port peak; the figure is recorded per grid. A model still ringing
there, but at least 40 dB down on both grids and within 6 dB of each other (the structure ringing down,
not one mesh), is judged on its TDR alone. The in-box window ends long before the cut, so the time record
it is read from is complete: records from the example FPGA board cut short to tails of -24 to -55 dB moved
the in-box TDR by 0.2 points at most. The spectrum of a cut ring is not complete, so no S-parameter result
is given and the grids are not compared on the transmission. The verdict says "usable for the TDR only".
