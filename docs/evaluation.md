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
| 3 | The board could not be evaluated (unreadable board, a crash in a check). The traceback is printed. |

Findings inside exclusion areas, waived findings and `unvalidated` checks do not change the code.

## `evaluation.json`, schema `openpcbsi-evaluation/1`

| Field | Content |
|-------|---------|
| `schema`, `tool`, `version`, `date` | Schema id, `openpcbsi`, openpcbsi version, ISO time of the run. |
| `board`, `board_sha256`, `downgraded` | Board path and hash; path of the KiCad 10 copy if the board was downgraded, else null. |
| `config`, `config_hash` | Config file, or `"detected from the board"`; hash of the effective config. |
| `detected` | Without a config: `assumptions` and `left_out` (lists of sentences), `fields`, `reference_nets`. Else null. |
| `waivers_file`, `waivers_hash` | The waiver file and its SHA-256 (12 hex digits), or null. |
| `solver_versions`, `kipi` | Version of each solver (as in G17's validation cache); KiPI's pin (`deps/kipi.lock`) when KiPI ran. |
| `kicad_drc` | EV1: `ran`, `kicad` version, `parity` (schematic found), `command`, `violations`, `unconnected`, `parity_issues`. |
| `exit_code` | As above. |
| `families` | Per family: `family`, `title`, `status`, `partial` (some checks not built), `checks`, `fail`, `warn`. |
| `checks` | Per check: `check`, `family`, `requirement`, `status`, counts `fail` / `warn` / `info` / `pass` / `in_area` / `waived`, `validation`, `why` (for `not run`, `not built`, `assumed inputs`). |
| `counts` | Totals outside areas and waivers, plus `in_area` and `waived`. |
| `waivers` | Each waiver as written, plus `status`: `applied`, `stale` or `invalid: ...`. |
| `findings` | Every finding (below), `info` included. |
| `report` | Path of the run's `report.json`. |

Status values, worst first: `fail`, `warn`, `unvalidated` (a solver check without a passed G17
case), `assumed inputs` (judged on assumed rail currents), `not run` (inputs missing, e.g. no
`ripple_pct` on any rail), `not built`, `pass`. A family takes the worst status of the checks that
ran.

Families: `G1` (general electrical rules), `RP` (return path and grounding), `SD` (supply
delivery: DC drop and PDN impedance), `SI` (impedance, crosstalk, net models, vias), `TH`
(thermal), `EMC`, `EV` (process: KiCad DRC), `PL` (placement) and `FS` (field-solver gates); the
last two are not built yet.

### A finding

| Field | Content |
|-------|---------|
| `id` | Stable ID `<check>-<10 hex>`: hash of the check, nets, references, layer, position rounded to 0.1 mm and the finding kind. A second finding with the same hash in one run gets `.1`, `.2`. |
| `severity` | `info`, `pass`, `warn` or `fail`. G1 `error` maps to `fail`, `warning` to `warn`. |
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

[rules.pad_plane]
fast_ghz = 3.0                # RP1 is judged on classes with Nyquist at or above this
```
