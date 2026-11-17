#!/usr/bin/env bash
# (Re)create openpcbsi's file parsers in this folder, idempotently (plan phase 2).
#   ./setup_parsers.sh           clone (if missing) at the pinned commits, install the site, smoke test
#   ./setup_parsers.sh --check   only verify the pins and the site
# Everything stays in this folder: the clones, the site packages, pip's cache and temp files.
set -euo pipefail
cd "$(dirname "$0")"
export PIP_CACHE_DIR="$PWD/.pip-cache" TMPDIR="$PWD/tmp"
mkdir -p "$PIP_CACHE_DIR" "$TMPDIR"
val() { python3 -c "import tomllib,sys; d=tomllib.load(open('parsers.lock','rb')); print(d['$1']['$2'] if isinstance(d['$1']['$2'], str) else ' '.join(d['$1']['$2']))"; }
for p in kicad_monkey kicadfiles; do
  REPO=$(val $p repo); COMMIT=$(val $p commit); DIR=$(val $p dir)
  if [[ "${1:-}" != "--check" ]]; then
    [[ -d "$DIR/.git" ]] || git clone --quiet --filter=blob:none "$REPO" "$DIR"
    git -C "$DIR" fetch --quiet origin "$COMMIT" 2>/dev/null || true
    git -C "$DIR" -c advice.detachedHead=false checkout --quiet "$COMMIT"
  fi
  HEAD=$(git -C "$DIR" rev-parse HEAD)
  [[ "$HEAD" == "$COMMIT" ]] || { echo "$DIR is at $HEAD, parsers.lock pins $COMMIT" >&2; exit 1; }
  [[ -z "$(git -C "$DIR" status --porcelain)" ]] || { echo "$DIR has local changes" >&2; exit 1; }
done
SITE=$(val site dir)
if [[ "${1:-}" != "--check" ]]; then
  python3 -m pip install --quiet --no-deps --upgrade --target "$SITE" $(val site requirements)
fi
PYTHONPATH="$SITE:$(val kicad_monkey python_path):$(val kicadfiles python_path)" python3 -I -c "
import sys; sys.path[:0] = '$SITE:$(val kicad_monkey python_path):$(val kicadfiles python_path)'.split(':')
import kicad_monkey.kicad_pcb_parser, kicadfiles, shapely
print('kicad-monkey', '$(val kicad_monkey version)', 'kicadfiles', '$(val kicadfiles version)', 'shapely', shapely.__version__, 'ok')"
