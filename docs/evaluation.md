# Design evaluation: files and schema

`python3 -m openpcbsi evaluate` writes `evaluation.json` and
`evaluation.html` into the run directory, beside `report.json`. This page describes the JSON, the
waiver file and the exit codes. The schema is versioned: a change that removes or renames a field
raises the version in `schema`.

## Exit codes

| Code | Meaning |
|------|---------|
| 0 | No `fail` and no `warn` outside exclusion areas and waivers. |
| 1 | Warnings only. Set `[evaluate] warnings_exit = false` in the board config for 0 instead. |
| 2 | At least one `fail`. |
| 3 | The board could not be evaluated: a board file newer than KiCad 10.0 (the message names both file versions), an unreadable board, or a crash in a check (the traceback is printed). |

Findings inside exclusion areas, waived findings, `pending` findings and `unvalidated` checks do not
change the code.

## `evaluation.json`, schema `openpcbsi-evaluation/1`

| Field | Content |
|-------|---------|
| `schema`, `tool`, `version`, `date` | Schema id, `openpcbsi`, openpcbsi version, ISO time of the run. |
| `board`, `board_sha256`, `downgraded` | Board path and hash; path of the KiCad 10 copy if the board was downgraded, else null. |
| `kicad_file_version`, `kicad` | The board file's format version and the KiCad that read it. |
| `stage`, `incomplete_nets`, `unrouted_pads` | Layout stage (`placement`, `routing`, `final`; detected, or `--stage`), the nets still in more than one copper island and the pads they leave unconnected. |
| `since` | With `--since`: `new`, `resolved` and `changed` findings against the earlier evaluation, by stable ID, and which of its waived findings are still present. |
| `em_jobs` | With `--em plan` / `run`: the automatic 3D regions (kind, nets, box, fix variant) and, with `run`, their results. |
| `config`, `config_hash` | Config file, or `"detected from the board"`; hash of the effective config. |
| `detected` | Without a config: `assumptions` and `left_out` (lists of sentences), `fields`, `reference_nets`. Else null. |
| `waivers_file`, `waivers_hash` | The waiver file and its SHA-256 (12 hex digits), or null. |
| `solver_versions`, `kipi` | Version of each solver (as in G17's validation cache); KiPI's pin (`deps/kipi.lock`) when KiPI ran. |
| `kicad_drc` | EV1: `ran`, `kicad` version, `parity` (schematic found), `command`, `violations`, `unconnected`, `parity_issues`. |
| `exit_code` | As above. |
| `families` | Per family: `family`, `title`, `status`, `partial` (some checks not built), `checks`, `fail`, `warn`. |
| `checks` | Per check: `check`, `family`, `requirement`, `status`, counts `fail` / `warn` / `pending` / `info` / `pass` / `in_area` / `waived`, `validation`, `why` (for `not run`, `not built`, `assumed inputs`). |
| `counts` | Totals outside areas and waivers (`fail`, `warn`, `pending`, `info`, `pass`), plus `in_area` and `waived`. |
| `waivers` | Each waiver as written, plus `status`: `applied`, `stale` or `invalid: ...`. |
| `findings` | Every finding (below), `info` included. |
| `report` | Path of the run's `report.json`. |

Status values, worst first: `fail`, `warn`, `pending` (only findings that wait for routing),
`unvalidated` (a solver check without a passed validation case), `not judged` (results computed but no
limit applies, e.g. impedances without a class target), `assumed inputs` (judged on assumed rail
currents), `nothing to check` (no item in scope), `not run` (inputs missing, e.g. no `ripple_pct` on
any rail), `not built`, `pass`. Only a check that judged something and found nothing wrong is a
`pass`. A family takes the worst status of the checks that ran.

Families: `G1` (general electrical rules), `RP` (return path and grounding), `SD` (supply
delivery: DC drop and PDN impedance), `SI` (impedance, crosstalk, net models, vias), `TH`
(thermal), `EMC`, `EV` (process: KiCad DRC), `PL` (placement) and `FS` (field-solver gates); the
last two are not built yet.

### A finding

| Field | Content |
|-------|---------|
| `id` | Stable ID `<check>-<10 hex>`: hash of the check, nets, references, layer, position rounded to 0.1 mm and the finding kind. A second finding with the same hash in one run gets `.1`, `.2`. |
| `severity` | `info`, `pass`, `pending`, `warn` or `fail`. G1 `error` maps to `fail`, `warning` to `warn`. Before the `final` stage, findings that only say "not connected yet" (KiCad's unconnected items, rails whose loads are not reached, dangling ends on incomplete nets) are `pending`. |
| `hint` | For geometric findings: a first fix (advice, not a constraint). |
| `check`, `family`, `requirement` | The G1 rule key or check name (`dc`, `pdn`, `kicad_drc`, ...), its family and the project's internal requirement ID (for example `RP1`, `SD1`). |
| `value`, `limit`, `unit` | Measured value and the limit it was judged against, where there is one. |
| `best` | Best achievable value, where a check computes one (none yet). |
| `nets`, `refs`, `layer`, `x`, `y` | Where it is (mm, board coordinates). |
| `reason` | One line. |
| `evidence` | The layer image, plot, `report.html` section or `kicad_drc.json` that shows it. |
| `area` | The exclusion area it lies in, or empty. Area findings are listed, not counted. |
| `waived` | Null, or `{reason, reviewer, waiver}`. |

## Waiver file

```toml
[[waiver]]
id = "split_crossing-3f9a0c1b2d"          # one finding, by its stable ID
reason = "L6 pair crosses the L5 split under the BGA breakout; reviewed in EM2"
reviewer = "board review 1"

[[waiver]]
check = "decap"                           # every finding of a check, narrowed by globs
refs = ["C19", "C2?"]
nets = ["3V3_AUX"]
area = "bga_U1"                           # optional: only findings in this exclusion area
reason = "3V3_AUX decoupling is placed for the monitor ICs, not the FPGA"
reviewer = "board review 1"
```

`reason` and `reviewer` are required; a waiver without them, or without `id` or `check`, is
reported as `invalid` and not applied. A waiver that matches nothing in the run is reported as
`stale`.

## Board config keys read by `evaluate`

```toml
[evaluate]
warnings_exit = true          # exit code 1 for warnings only
drc_warning = "warn"          # KiCad DRC warnings as "warn" or "fail" (EV1 literally: "fail")
dc_warn_pct = 50.0            # SD1: warn above this share of the window
engine_disagree_pct = 20.0    # SD1: KiPI and openpcbsi differing by more is a warn
drc_timeout_s = 900

trace_warn_c = 10.0           # SD3: track temperature rise (IPC-2221) warn / fail
trace_fail_c = 30.0
via_amps_per_0p2mm = 0.5      # SD4: current per via, scaled with the drill

[rules.pad_plane]
fast_ghz = 3.0                # RP1 is judged on classes with Nyquist at or above this

[thermal]
tj_max_c = 100.0              # junction limit when a part gives none (field Tj_max, or [thermal.parts])

[xsec.corners]                # manufacturing corners for every judged geometry
etch_um = 20.0                # width +/- 2 x etch (gap the other way)
thickness_pct = 10.0
er_pct = 5.0

[timing]
skew_ps = 10.0                # group skew limit (groups in [timing.groups])
pair_skew_ps = 1.0            # P/N skew limit (default: the interface standard's)

[em]
fast_ghz = 3.0                # --em plan / run: classes at or above this get automatic 3D regions
```

Limits are taken, in this order, from: a `.kicad_dru` rule that covers the net; the board config;
without a reviewed config, the interface standard the net's name maps to (`standards.py`); the
defaults above.

## Sign-off record

`python3 -m openpcbsi signoff BOARD [-c CFG] [--waivers FILE] [--reviewer NAME] [--sign-key KEY]` runs
the evaluation and writes `signoff.json` (schema `openpcbsi-signoff/2`) beside it: the SHA-256 of the
board, of each KiCad project file next to it (`.kicad_pro`, `.kicad_dru`, `.kicad_sch`), of the
config, the waivers and `evaluation.json`; the KiCad file version and openpcbsi version; every
check's status; the verdict (`release`, `release with warnings`, `blocked: failures`, `blocked: could
not evaluate`, or `blocked: layout stage ...` before the final stage); the blocking findings; and the
checks that were not judged (not built, not run, nothing to check, not judged, unvalidated, assumed
inputs). Paths in it are relative to its folder. With `--sign-key` (an SSH private key) the record is
signed: `ssh-keygen -Y sign`, namespace `openpcbsi-signoff`, into `signoff.json.sig`.

`python3 -m openpcbsi verify-signoff signoff.json [--allowed-signers FILE] [--rerun]` exits 1 when:

- any input's hash changed (a record never vouches for a board it did not judge);
- any field that follows from the hashed evaluation differs from it: the exit code, the counts, the
  checks, the families, the blockers or the verdict are recomputed from `evaluation.json`, never
  read from the record (an edited verdict or emptied blocker list does not verify);
- with `--allowed-signers` (an ssh allowed_signers file), the signature does not verify for the
  record's reviewer.

Only the signature authenticates the record as a whole and the reviewer's name: a record and its
evaluation edited together, consistently, verify by their hashes alone. `--rerun` evaluates the board
again with the recorded config and waivers and compares the verdict, which needs neither.
