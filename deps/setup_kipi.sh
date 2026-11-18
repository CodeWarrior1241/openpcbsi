#!/usr/bin/env bash
# (Re)create openpcbsi's KiPI dependency, idempotently.
#   ./setup_kipi.sh            clone (if missing) at the pinned commit, build the venv by hash, smoke test
#   ./setup_kipi.sh --check    only verify: the commit, the venv against the pins, the smoke test
# Target: this folder, or $OPENPCBSI_DEPS when set (the folder openpcbsi runs KiPI from; kipi.lock and
# kipi_driver.py are copied there). The clone, the venv, pip's cache and temp files all stay in the
# target (KiPI's own ~/.local/share/kipi is not used).
# Pins: kipi.lock (the clone), requirements/kipi.txt and requirements/kipi-pip.txt (every
# distribution in the venv, with SHA-256s; regenerate with `python3 deps/pins.py hash` after editing
# requirements/kipi.in or kipi-pip.in).
set -euo pipefail
SRC=$(cd "$(dirname "$0")" && pwd)
DEPS=${OPENPCBSI_DEPS:-$SRC}
mkdir -p "$DEPS"; DEPS=$(cd "$DEPS" && pwd)
CHECK=0; [[ "${1:-}" == "--check" ]] && CHECK=1
cd "$DEPS"
get() { sed -n "s/^$1 *= *\"\([^\"]*\)\".*/\1/p" "$SRC/kipi.lock"; }
REPO=$(get repo); TAG=$(get tag); COMMIT=$(get commit); PY=$(get python); VENV=$(get venv)
PIPREQ="$SRC/$(get pip)"; REQ="$SRC/$(get requirements)"; UPSTREAM=$(get upstream)
export PIP_CACHE_DIR="$DEPS/.pip-cache" TMPDIR="$DEPS/tmp" KIPI_HOME="$DEPS/kipi-home"
mkdir -p "$PIP_CACHE_DIR" "$TMPDIR" "$KIPI_HOME"
if [[ "$DEPS" != "$SRC" && $CHECK == 0 ]]; then cp "$SRC/kipi.lock" "$SRC/kipi_driver.py" "$DEPS/"; fi
SITEPKG() { "$VENV/bin/python" -I -c "import sysconfig; print(sysconfig.get_paths()['purelib'])"; }

if [[ $CHECK == 0 ]]; then
  if [[ ! -d KiPI/.git ]]; then
    git clone --quiet --branch "$TAG" "$REPO" KiPI
  fi
  git -C KiPI fetch --quiet --tags origin "$TAG" 2>/dev/null || true
  git -C KiPI -c advice.detachedHead=false checkout --quiet "$COMMIT"
  if [[ ! -x "$VENV/bin/python" ]]; then
    "$PY" -m venv "$VENV"          # its bundled pip comes with the Python install, not the network
  fi
  if ! python3 -I "$SRC/pins.py" site "$(SITEPKG)" "$PIPREQ" "$REQ" >/dev/null 2>&1; then
    PIPI=("$VENV/bin/python" -m pip install --quiet --disable-pip-version-check --require-hashes --only-binary :all:)
    "${PIPI[@]}" --no-deps -r "$PIPREQ"
    "${PIPI[@]}" -r "$REQ"
  fi
fi

[[ -d KiPI/.git ]] || { echo "$DEPS/KiPI is missing: run setup_kipi.sh" >&2; exit 1; }
HEAD=$(git -C KiPI rev-parse HEAD)
if [[ "$HEAD" != "$COMMIT" ]]; then
  echo "KiPI is at $HEAD, kipi.lock pins $COMMIT" >&2; exit 1
fi
if [[ -n "$(git -C KiPI status --porcelain)" ]]; then
  echo "KiPI clone has local changes (git -C $DEPS/KiPI status)" >&2; exit 1
fi
[[ -x "$VENV/bin/python" ]] || { echo "$DEPS/$VENV is missing: run setup_kipi.sh" >&2; exit 1; }
# exactly the pinned distributions, files unmodified (an extra one: `$VENV/bin/pip uninstall` it)
python3 -I "$SRC/pins.py" site "$(SITEPKG)" "$PIPREQ" "$REQ"
OUT=$("$VENV/bin/python" -I -m pip check --disable-pip-version-check) || {
  echo "$OUT" >&2; echo "pip check fails in $DEPS/$VENV" >&2; exit 1; }
# the pins must satisfy KiPI's own ranges
"$VENV/bin/python" -I - "$UPSTREAM" <<'EOF'
import sys, importlib.metadata as md
from packaging.requirements import Requirement
bad = []
for line in open(sys.argv[1]):
    line = line.split("#")[0].strip()
    if not line:
        continue
    r = Requirement(line)
    if r.marker and not r.marker.evaluate():
        continue
    try:
        v = md.version(r.name)
    except md.PackageNotFoundError:
        bad.append(f"{r.name} missing")
        continue
    if not r.specifier.contains(v, prereleases=True):
        bad.append(f"{r.name} {v} outside KiPI's {r.specifier}")
if bad:
    sys.exit("pins do not satisfy " + sys.argv[1] + ": " + "; ".join(bad))
EOF
"$VENV/bin/python" -I -c "
import sys; sys.path.insert(0, 'KiPI')
from kipi.dc.pi.tree import solve_tree
from kipi.ac.solver.sweep import solve_ac
from kipi.sources.kicad_pcb.source import KiCadFileSource
print('KiPI', '$TAG', '$COMMIT'[:10], 'ok in', sys.executable)"
