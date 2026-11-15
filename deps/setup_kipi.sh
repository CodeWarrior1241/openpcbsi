#!/usr/bin/env bash
# (Re)create openlynx's KiPI dependency in this folder, idempotently.
#   ./setup_kipi.sh            clone (if missing) at the pinned commit, build the venv, smoke test
#   ./setup_kipi.sh --check    only verify the pin and the venv
# Everything stays in this folder: the clone, the venv, pip's cache and its temp files
# (nothing on the root disk; KiPI's own ~/.local/share/kipi is not used).
set -euo pipefail
cd "$(dirname "$0")"
get() { sed -n "s/^$1 *= *\"\([^\"]*\)\".*/\1/p" kipi.lock; }
REPO=$(get repo); TAG=$(get tag); COMMIT=$(get commit); PY=$(get python); VENV=$(get venv); REQ=$(get requirements)
export PIP_CACHE_DIR="$PWD/.pip-cache" TMPDIR="$PWD/tmp" KIPI_HOME="$PWD/kipi-home"
mkdir -p "$PIP_CACHE_DIR" "$TMPDIR" "$KIPI_HOME"

if [[ "${1:-}" != "--check" ]]; then
  if [[ ! -d KiPI/.git ]]; then
    git clone --quiet --branch "$TAG" "$REPO" KiPI
  fi
  git -C KiPI fetch --quiet --tags origin "$TAG" 2>/dev/null || true
  git -C KiPI -c advice.detachedHead=false checkout --quiet "$COMMIT"
  if [[ ! -x "$VENV/bin/python" ]]; then
    "$PY" -m venv "$VENV"
    "$VENV/bin/python" -m pip install --quiet --upgrade pip
  fi
  "$VENV/bin/python" -m pip install --quiet -r "$REQ"
fi

HEAD=$(git -C KiPI rev-parse HEAD)
if [[ "$HEAD" != "$COMMIT" ]]; then
  echo "KiPI is at $HEAD, kipi.lock pins $COMMIT" >&2; exit 1
fi
if [[ -n "$(git -C KiPI status --porcelain)" ]]; then
  echo "KiPI clone has local changes (git -C deps/KiPI status)" >&2; exit 1
fi
"$VENV/bin/python" -I -c "
import sys; sys.path.insert(0, 'KiPI')
from kipi.dc.pi.tree import solve_tree
from kipi.ac.solver.sweep import solve_ac
from kipi.sources.kicad_pcb.source import KiCadFileSource
print('KiPI', '$TAG', '$COMMIT'[:10], 'ok in', sys.executable)"
