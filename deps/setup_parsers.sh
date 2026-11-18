#!/usr/bin/env bash
# (Re)create openpcbsi's file parsers, idempotently (plan phase 2).
#   ./setup_parsers.sh           clone (if missing) at the pinned commits, check them against the PyPI
#                                sdists, install the site by hash, smoke test
#   ./setup_parsers.sh --check   only verify: the commits, the sdists, the installed site (no install)
# Target: this folder, or $OPENPCBSI_DEPS when set (the folder openpcbsi reads them from; the lock is
# copied there). The clones, the site packages, pip's cache and temp files all stay in the target.
# Pins: parsers.lock (the clones), requirements/parsers.txt (the site, with SHA-256s; regenerate with
# `python3 deps/pins.py hash` after editing requirements/parsers.in).
set -euo pipefail
SRC=$(cd "$(dirname "$0")" && pwd)
DEPS=${OPENPCBSI_DEPS:-$SRC}
mkdir -p "$DEPS"; DEPS=$(cd "$DEPS" && pwd)
CHECK=0; [[ "${1:-}" == "--check" ]] && CHECK=1
cd "$DEPS"
export PIP_CACHE_DIR="$DEPS/.pip-cache" TMPDIR="$DEPS/tmp"
mkdir -p "$PIP_CACHE_DIR" "$TMPDIR"
val() { python3 -I -c "import tomllib; print(tomllib.load(open('$SRC/parsers.lock','rb'))['$1']['$2'])"; }
if [[ "$DEPS" != "$SRC" && $CHECK == 0 ]]; then cp "$SRC/parsers.lock" "$DEPS/parsers.lock"; fi

for p in kicad_monkey kicadfiles; do
  REPO=$(val $p repo); COMMIT=$(val $p commit); DIR=$(val $p dir)
  if [[ $CHECK == 0 ]]; then
    [[ -d "$DIR/.git" ]] || git clone --quiet --filter=blob:none "$REPO" "$DIR"
    git -C "$DIR" fetch --quiet origin "$COMMIT" 2>/dev/null || true
    git -C "$DIR" -c advice.detachedHead=false checkout --quiet "$COMMIT"
  fi
  [[ -d "$DIR/.git" ]] || { echo "$DEPS/$DIR is missing: run setup_parsers.sh" >&2; exit 1; }
  HEAD=$(git -C "$DIR" rev-parse HEAD)
  [[ "$HEAD" == "$COMMIT" ]] || { echo "$DIR is at $HEAD, parsers.lock pins $COMMIT" >&2; exit 1; }
  [[ -z "$(git -C "$DIR" status --porcelain)" ]] || { echo "$DIR has local changes" >&2; exit 1; }
  # the source must be the published release: sdist SHA-256 and file-by-file equality with the clone
  python3 -I "$SRC/pins.py" sdist "$DEPS" "$p"
done

SITE=$(val site dir); REQ="$SRC/$(val site requirements)"
if [[ $CHECK == 0 ]] && ! python3 -I "$SRC/pins.py" site "$SITE" "$REQ" >/dev/null 2>&1; then
  # build a fresh site beside the old one and swap, so nothing unpinned survives
  rm -rf "$SITE.new"
  python3 -m pip install --quiet --disable-pip-version-check --no-deps --require-hashes \
    --only-binary :all: --target "$SITE.new" -r "$REQ"
  rm -rf "$SITE.old"; if [[ -d "$SITE" ]]; then mv "$SITE" "$SITE.old"; fi
  mv "$SITE.new" "$SITE"; rm -rf "$SITE.old"
fi
[[ -d "$SITE" ]] || { echo "$DEPS/$SITE is missing: run setup_parsers.sh" >&2; exit 1; }
python3 -I "$SRC/pins.py" site "$SITE" "$REQ"

python3 -I -c "
import sys; sys.path[:0] = '$SITE:$(val kicad_monkey python_path):$(val kicadfiles python_path)'.split(':')
import kicad_monkey.kicad_pcb_parser, kicadfiles, shapely, numpy
print('kicad-monkey', '$(val kicad_monkey version)', 'kicadfiles', '$(val kicadfiles version)',
      'shapely', shapely.__version__, 'numpy (host)', numpy.__version__, 'ok in', '$DEPS')"
