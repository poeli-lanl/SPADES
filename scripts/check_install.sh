#!/bin/bash
set -uo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO_DIR=$(cd "$SCRIPT_DIR/.." && pwd)
STATUS=0

info() {
  printf 'INFO: %s\n' "$*"
}

ok() {
  printf 'OK: %s\n' "$*"
}

warn() {
  printf 'WARN: %s\n' "$*" >&2
}

err() {
  printf 'ERROR: %s\n' "$*" >&2
  STATUS=1
}

require_command() {
  local cmd="$1"
  if command -v "$cmd" >/dev/null 2>&1; then
    ok "command found: $cmd"
  else
    err "command not found: $cmd"
  fi
}

require_file() {
  local path="$1"
  if [[ -f "$path" ]]; then
    ok "file found: $path"
  else
    err "file missing: $path"
  fi
}

require_dir() {
  local path="$1"
  if [[ -d "$path" ]]; then
    ok "directory found: $path"
  else
    err "directory missing: $path"
  fi
}

info "Checking SPADES installation from $REPO_DIR"

for cmd in gottcha2 fastp fastplong samtools sylph ktImportTaxonomy minimap2 bcftools; do
  require_command "$cmd"
done

if command -v python >/dev/null 2>&1; then
  python - <<'PY'
import importlib
import sys

modules = ["pandas", "numpy", "pysam", "biom", "minify_html"]
missing = []
for module in modules:
    try:
        importlib.import_module(module)
    except Exception as exc:
        missing.append(f"{module} ({exc})")

if missing:
    print("ERROR: missing Python modules: " + ", ".join(missing), file=sys.stderr)
    sys.exit(1)

print("OK: required Python modules import")
PY
  if [[ $? -ne 0 ]]; then
    STATUS=1
  fi
else
  err "command not found: python"
fi

require_file "$REPO_DIR/data/pathogen.tsv"
require_dir "$REPO_DIR/data/taxonomy_db"
# require_file "$REPO_DIR/data/taxonomy_db/names.dmp"
# require_file "$REPO_DIR/data/taxonomy_db/nodes.dmp"

if [[ -n "${CONDA_PREFIX:-}" ]]; then
  KRONA_TAXONOMY="$CONDA_PREFIX/opt/krona/taxonomy"
  if [[ -e "$KRONA_TAXONOMY" ]]; then
    ok "Krona taxonomy path exists: $KRONA_TAXONOMY"
  else
    err "Krona taxonomy path missing: $KRONA_TAXONOMY. Run ./post_install.sh"
  fi
else
  warn "CONDA_PREFIX is not set; skipping Krona taxonomy path check"
fi

TEST_DB="$REPO_DIR/test/gottcha2_database_test/gottcha_db.species.fna"
if [[ -d "$REPO_DIR/test/gottcha2_database_test" ]]; then
  for ext in syldb zip stats tax.tsv; do
    if [[ -f "$TEST_DB.$ext" ]]; then
      ok "bundled test database asset found: $TEST_DB.$ext"
    else
      warn "bundled test database asset missing: $TEST_DB.$ext"
    fi
  done
else
  warn "test database directory not present; skipping bundled test database checks"
fi

if [[ $STATUS -eq 0 ]]; then
  info "SPADES installation check passed."
else
  err "SPADES installation check failed."
fi

exit "$STATUS"
