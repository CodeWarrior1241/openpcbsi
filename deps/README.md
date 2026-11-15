# deps: external tools openpcbsi drives

Everything an external engine needs lives in this folder: its source, its Python environment,
pip's cache and its temp files. Nothing goes to the root disk.

## KiPI (G10 DC drop, G11/G12 PDN impedance)

| Path | What |
|------|------|
| `kipi.lock` | The pin: repo, tag (v1.0.1), commit, Python, venv. The single source of truth. |
| `setup_kipi.sh` | Recreates the clone and venv from the pin, then smoke-tests the entry points. `--check` only verifies. |
| `KiPI/` | Clone of drewsum/KiPI at the pinned commit. Do not edit; the setup script refuses a dirty clone. |
| `kipi-venv/` | KiPI's own environment (system Python 3.12, its `requirements-standalone.txt`), about 1.3 GB. |
| `kipi_driver.py` | Runs inside `kipi-venv`: loads the board and setup, calls KiPI, writes plain JSON. |
| `kipi-home/` | `KIPI_HOME` for runs from openpcbsi (KiPI would otherwise use `~/.local/share/kipi`). |
| `.pip-cache/`, `tmp/` | pip's cache and build directory for this folder. |

How openpcbsi uses it (`openpcbsi/kipi_bridge.py`):

1. It derives KiPI's power tree from the board config's `[rails.*]`. Each rail's `source` (the
   shunt) becomes a series link. The walk continues upstream through parts of 1 ohm or less and
   ferrites (DCR from `[kipi] ferrite_dcr_ohm`) to the rail's `regulator`. The `sense` net gives
   the VRM's sense pin. Loads are the `loads` parts' pins on the rail, with their ground pins.
   Probes and targets come from `ripple_pct` / `step_a`. Capacitor models are grouped by value
   and package.
2. It writes `<run>/kipi/<board>.kipi-setup.json` (KiPI reads JSON setups; the file opens in KiPI's
   GUI as well) and `<run>/kipi/kipi_job.json`.
3. It runs `kipi-venv/bin/python -I deps/kipi_driver.py <run>/kipi/kipi_job.json`, so openpcbsi never
   imports KiPI.
4. It maps `kipi_result.json` back into the report's G10 / G11 tables. openpcbsi's own solver runs
   beside it when `[engines] cross_check = true`, and the report shows the difference.

Which engine runs is set per board config:

```toml
[engines]
dc = "kipi"          # or "openpcbsi"
pdn = "kipi"         # or "openpcbsi"
cross_check = true   # also run the other engine and tabulate the difference
```

Any value other than `"kipi"` selects openpcbsi's own solver, so configs written before the rename
with `"openlynx"` behave the same.

KiPI's board snapshots are cached in `../cache/kipi/` by the board's SHA-256.

### Updating the pin

1. Edit `tag` and `commit` in `kipi.lock`.
2. Run `./setup_kipi.sh` to check out that commit and install its requirements.
3. Re-check `kipi_driver.py`. It calls internal API: `kipi.dc.pi.tree.solve_tree`,
   `kipi.ac.solver.sweep.solve_ac`, `kipi.sources.kicad_pcb.source.KiCadFileSource`,
   `kipi.core.model.snapshot_io`, `kipi.dc.pi.model.PISetup`, and the result dataclasses'
   field names.
4. Run `python3 -m openpcbsi check <pilot board> -c boards/au7p-fmc.toml --checks dc,pdn`
   and compare the cross-check table with the previous run.

### Viewing a setup in KiPI's GUI

```sh
cd /media/fpgadev/Dev_Tools/Work/Open_Source_Hyperlynx
KIPI_HOME=deps/kipi-home PYTHONPATH=deps/KiPI deps/kipi-venv/bin/python -m kipi BOARD.kicad_pcb
```

KiPI picks up a setup that sits next to the board as `<board>.kipi-setup.yml`, or as
`<board>.kipi-setup.json` (the older name, which it still reads). Copy `<run>/kipi/<board>.kipi-setup.json`
there, onto a copy of the board project, not into a board folder another tool owns (the pilot's
`kicad/` folder is read-only for openpcbsi).
