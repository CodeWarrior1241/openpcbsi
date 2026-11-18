# deps: external tools openpcbsi drives

Everything an external engine needs lives in this folder: its source, its Python environment,
pip's cache and its temp files. Nothing goes to the root disk. The setup scripts install into this
folder, or into `$OPENPCBSI_DEPS` when it is set (the folder openpcbsi then reads them from; the
scripts copy the lock files there).

Every Python distribution is pinned to an exact version and installed by SHA-256
(`pip install --require-hashes --only-binary :all:`), so every user gets the same files; no step
installs whatever PyPI serves that day, pip included. Both scripts have a `--check` mode that installs
nothing and fails loudly when the installed state differs from the locks.

## KiPI (G10 DC drop, G11/G12 PDN impedance)

| Path | What |
|------|------|
| `kipi.lock` | The pin: repo, tag (v1.0.1), commit, Python, venv, and the hashed requirement files below. |
| `requirements/kipi.in` | Every distribution in the venv, transitive ones included, at an exact version (taken from the working venv's `pip freeze`). They must satisfy KiPI's own ranges in `KiPI/requirements-standalone.txt`; the setup script checks that. |
| `requirements/kipi.txt` | Generated from `kipi.in`: the same pins with the SHA-256 of every wheel. What the setup script installs. |
| `requirements/kipi-pip.in`, `.txt` | pip 26.3 for the venv, by hash, installed over the venv's bundled pip (which comes with the Python install, not the network). |
| `setup_kipi.sh` | Recreates the clone and venv from the pins, then checks and smoke-tests them. `--check` only verifies: the commit, a clean clone, the venv's distributions exactly equal to the pins (none extra), every installed file equal to its RECORD hash, `pip check`, KiPI's ranges, the imports. |
| `KiPI/` | Clone of drewsum/KiPI at the pinned commit. Do not edit; the setup script refuses a dirty clone. |
| `kipi-venv/` | KiPI's own environment (system Python 3.12), about 1.3 GB. |
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
2. If KiPI's `requirements-standalone.txt` changed, or a dependency must move, edit the versions in
   `requirements/kipi.in` (see [Moving a Python pin](#moving-a-python-pin)).
3. Run `./setup_kipi.sh` to check out that commit and install the pinned environment.
4. Re-check `kipi_driver.py`. It calls internal API: `kipi.dc.pi.tree.solve_tree`,
   `kipi.ac.solver.sweep.solve_ac`, `kipi.sources.kicad_pcb.source.KiCadFileSource`,
   `kipi.core.model.snapshot_io`, `kipi.dc.pi.model.PISetup`, and the result dataclasses'
   field names.
5. Run `python3 -m openpcbsi check <pilot board> -c boards/au7p-fmc.toml --checks dc,pdn`
   and compare the cross-check table with the previous run.

### Viewing a setup in KiPI's GUI

```sh
cd <openpcbsi checkout>
KIPI_HOME=deps/kipi-home PYTHONPATH=deps/KiPI deps/kipi-venv/bin/python -m kipi BOARD.kicad_pcb
```

KiPI picks up a setup that sits next to the board as `<board>.kipi-setup.yml`, or as
`<board>.kipi-setup.json` (the older name, which it still reads). Copy `<run>/kipi/<board>.kipi-setup.json`
there, onto a copy of the board project, not into a board folder another tool owns (the pilot's
`kicad/` folder is read-only for openpcbsi).

## File parsers (reading the board)

| Path | What |
|------|------|
| `parsers.lock` | The pins: kicad-monkey 2026.9.7 (commit 4c029a39) and KiCadFiles 1.0.5 (tag v1.0.5), each with the URL and SHA-256 of its PyPI sdist, and the hashed requirement file for kicad-monkey's dependencies. |
| `requirements/parsers.in`, `.txt` | kicad-monkey's dependencies at exact versions; the `.txt` adds the wheels' SHA-256s. |
| `setup_parsers.sh` | Clones both at the pinned commits and checks each against its sdist: the sdist's SHA-256 must equal `sdist_sha256`, and the clone's tracked files under `clone_tree` must equal the sdist's files under `sdist_tree`, file for file and byte for byte (the git commit is what is imported; the sdist is the published release it must match). Installs the dependencies into `parsers-site/` by hash (`--no-deps`: the host's numpy is used), checks the site against the pins and RECORD hashes, smoke-tests the imports. `--check` does all the checks and installs nothing. The sdists are cached in `.pip-cache/sdist/` and re-hashed on every run. |
| `kicad_monkey/` | Clone of wavenumber-eng/kicad_monkey: the primary reader (openpcbsi/parse/monkey.py). |
| `KiCadFiles/` | Clone of Steffen-W/KiCadFiles: the independent cross-check (openpcbsi/parse/kicadfiles.py). |
| `parsers-site/` | kicad-monkey's dependencies (shapely, msgspec, ...). |

Only `openpcbsi/parse/` imports the parsers, so a change of their APIs stays there. Change a pin only
on purpose, then run acceptance tests 50 to 52: the corpus agreement with KiCad's own reading.
To move kicad-monkey or KiCadFiles, update `version`, `commit`, `sdist_url` and `sdist_sha256` together
(PyPI's JSON API, `https://pypi.org/pypi/<name>/<version>/json`, lists the sdist's URL and SHA-256);
`setup_parsers.sh` then refuses a commit whose Python source differs from that release.

## Moving a Python pin

The `.in` files under `requirements/` are the pins; the `.txt` files are generated from them and are
what the setup scripts install. After editing a version in a `.in` file:

```sh
python3 deps/pins.py hash           # rewrites requirements/*.txt with the hashes PyPI lists
python3 deps/pins.py hash --check   # (optional) regenerates in memory and reports any difference
deps/setup_kipi.sh                  # or deps/setup_parsers.sh: installs and verifies
```

`pins.py hash` reads each pin's files from PyPI's JSON API and keeps the SHA-256 of every wheel that
fits the Pythons named in the `.in` file's `# pythons:` line (3.12 and 3.13) on Linux (manylinux
x86_64, aarch64), macOS (x86_64, arm64) and Windows (x86-64). It refuses a pin without a Linux x86_64
wheel and lists the other gaps as `# no wheel:` comments at the top of the `.txt`. Sdists are left out
on purpose: building one would fetch build tools that are not hashed. KiPI's file must hold the whole
closure: `pip install --require-hashes` fails on any dependency that is missing from it, and
`setup_kipi.sh` fails on any installed distribution that is not pinned (remove it with
`kipi-venv/bin/pip uninstall`, or delete `kipi-venv/` and rerun). A dependency that only one platform
needs takes a marker, as the macOS-only `pyobjc` pins in `kipi.in` do (`; sys_platform == "darwin"`).

Covered and not covered: the Linux x86_64 / CPython 3.12 environments are built and tested from these
files. Windows x86-64 and macOS arm64 resolve from them with pip's `--platform` dry run (which cannot
evaluate the darwin markers, so the `pyobjc` pins are unchecked). There is no wheel of triangle
20250106 for macOS x86_64 or of triangle and pyamg for Linux aarch64, so KiPI's environment does not
install there.
