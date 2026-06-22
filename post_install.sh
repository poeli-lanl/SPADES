#!/bin/bash
set -euo pipefail

# Run this after installing Krona in the active conda environment.
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

log() {
  printf 'INFO: %s\n' "$*"
}

fail() {
  printf 'ERROR: %s\n' "$*" >&2
  exit 2
}

if [[ -z "${CONDA_PREFIX:-}" ]]; then
  fail "CONDA_PREFIX is not set. Activate the SPADES conda environment first."
fi

KRONA_TAXONOMY="$CONDA_PREFIX/opt/krona/taxonomy"
SPADES_TAXONOMY="$SCRIPT_DIR/data/taxonomy_db"
KRONA_DIR=$(dirname "$KRONA_TAXONOMY")

[[ -d "$SPADES_TAXONOMY" ]] || fail "SPADES taxonomy directory not found: $SPADES_TAXONOMY"
[[ -f "$SPADES_TAXONOMY/taxonomy.tab.gz" ]] || fail "Missing taxonomy file: $SPADES_TAXONOMY/taxonomy.tab.gz"

gzip -d $SPADES_TAXONOMY/taxonomy.tab.gz

if [[ -L "$KRONA_TAXONOMY" || -e "$KRONA_TAXONOMY" ]]; then
  rm -rf "$KRONA_TAXONOMY"
fi

ln -s "$SPADES_TAXONOMY" "$KRONA_TAXONOMY"
log "Linked Krona taxonomy: $KRONA_TAXONOMY -> $SPADES_TAXONOMY"
