#!/bin/bash
CWD=$(dirname $(realpath "$0"))
RUN_COMMAND="$0 $*"
export PATH="$CWD/scripts:$PATH"
set -uo pipefail

usage() {
  cat <<'EOF'
Usage:
  run_SPADES.sh \
    -i <reads> \
    OR
    -1 <read1> -2 <read2> \
    -o <outdir> \
    -p <prefix> \
    -d <db_path> \
    -t <cpu> \
    [--js-external] \
    [--spades-data <spades_data_dir>] \
    [--min-depth <min_depth>] \
    [--ont] \
    [--ont-error-rate <float>] \
    [--clean]

Input:
  -i, --input         Single-end reads file (fastq/fq[.gz])
  -1, --read1         Paired-end read 1 (fastq/fq[.gz])
  -2, --read2         Paired-end read 2 (fastq/fq[.gz])

Required:
  -o, --outdir        Output directory
  -p, --prefix        Output prefix (base filename)
  -d, --db-path       GOTTCHA2 database path
  -t, --cpu           Threads/CPUs (positive integer)

Optional:
  --spades-data       Directory containing taxonomy_db/ and pathogen.tsv
  --min-depth         Minimum depth for variant calling (default: 10)
  --ont               Treat input as long reads; will split to 150bp and pass -np to gottcha2
  --ont-error-rate    Error rate for ONT reads (passed to gottcha2 -er), default: 0.03
  --js-external       Use CDN-hosted JavaScript/CSS assets in generated HTML reports
  --clean             Remove large intermediates after the run
  --version           Show script version and exit

EOF
}

# variables needed (declare + default)
VERSION="1.2.1"
INPUT=""
INPUT_QC=""
READ=""
OUTDIR=""
PREFIX=""
DB_PATH=""
CPU="1"
ONT="false"
SPADES_DATA="$CWD/data"
ONT_FLAG=""
ONT_ERROR_RATE="0.03"
CLEAN_FLAG=""
READ1=""
READ2=""
PAIRED="false"
INPUT_QC_R1=""
INPUT_QC_R2=""
JS_EXTERNAL="false"
MIN_DEPTH="10"
DB_LEVEL=""

# parse args
while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--input)        INPUT="${2:-}"; shift 2 ;;
    -1|--read1)        READ1="${2:-}"; shift 2 ;;
    -2|--read2)        READ2="${2:-}"; shift 2 ;;
    -o|--outdir)       OUTDIR="${2:-}"; shift 2 ;;
    -p|--prefix)       PREFIX="${2:-}"; shift 2 ;;
    -d|--db-path)      DB_PATH="${2:-}"; shift 2 ;;
    -t|--cpu)          CPU="${2:-}"; shift 2 ;;
    --spades-data)     SPADES_DATA="${2:-}"; shift 2 ;;
    --ont)             ONT="true"; shift ;;
    --ont-error-rate)  ONT_ERROR_RATE="${2:-}"; shift 2 ;;
    --js-external)     JS_EXTERNAL="true"; shift ;;
    --clean)           CLEAN_FLAG="true"; shift ;; 
    --min-depth)       MIN_DEPTH="${2:-}"; shift 2 ;;
    --version)         echo "$VERSION"; exit 0 ;;
    -h|--help)         usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage; exit 2 ;;
  esac
done

# validate input mode
if [[ -n "$READ1" || -n "$READ2" ]]; then
  [[ -z "$INPUT" ]] || {
    echo "ERROR: --input cannot be combined with --read1/--read2" >&2
    exit 2
  }
  [[ "$ONT" != "true" ]] || {
    echo "ERROR: --ont requires single-end input and cannot be combined with paired-end reads" >&2
    exit 2
  }
  [[ -n "$READ1" && -n "$READ2" ]] || {
    echo "ERROR: Both --read1 and --read2 are required for paired-end input" >&2
    exit 2
  }
  [[ -f "$READ1" && -f "$READ2" ]] || {
    echo "ERROR: Paired-end input files not found" >&2
    exit 2
  }
  PAIRED="true"
else
  [[ -n "$INPUT" ]] || {
    echo "ERROR: missing --input or --read1/--read2" >&2
    exit 2
  }
  [[ -f "$INPUT" ]] || {
    echo "ERROR: input not found: $INPUT" >&2
    exit 2
  }
fi

# validate required args
[[ -n "$OUTDIR" ]]          || { echo "ERROR: missing --outdir" >&2; usage; exit 2; }
[[ -n "$PREFIX" ]]          || { echo "ERROR: missing --prefix" >&2; usage; exit 2; }
[[ -n "$DB_PATH" ]]         || { echo "ERROR: missing --db-path" >&2; usage; exit 2; }
DB_PATH="${DB_PATH%.syldb}"
DB_PATH="${DB_PATH%.zip}"
DB_PATH="${DB_PATH%.stats}"
DB_PATH="${DB_PATH%.tax.tsv}"
[[ -f "$DB_PATH.syldb" ]]   || { echo "ERROR: GOTTCHA2 fast-profile database index not found: $DB_PATH.syldb" >&2; exit 2; }
[[ -f "$DB_PATH.zip" ]]     || { echo "ERROR: GOTTCHA2 signature archive not found: $DB_PATH.zip" >&2; exit 2; }
[[ -f "$DB_PATH.tax.tsv" ]] || { echo "ERROR: GOTTCHA2 taxonomy file not found: $DB_PATH.tax.tsv" >&2; exit 2; }
[[ -f "$DB_PATH.stats" ]]   || { echo "ERROR: GOTTCHA2 stats file not found: $DB_PATH.stats" >&2; exit 2; }
[[ "$CPU" =~ ^[0-9]+$ && "$CPU" -gt 0 ]]|| { echo "ERROR: --cpu must be a positive integer" >&2; exit 2; }
[[ -f "$SPADES_DATA/pathogen.tsv" ]]|| { echo "ERROR: missing $SPADES_DATA/pathogen.tsv" >&2; exit 2; }
[[ "$ONT_ERROR_RATE" =~ ^[0-9]+(\.[0-9]+)?$ ]] || { echo "ERROR: --ont-error-rate must be numeric (e.g., 0.03)" >&2; exit 2; }
[[ "$MIN_DEPTH" =~ ^[0-9]+$ && "$MIN_DEPTH" -gt 0 ]] || { echo "ERROR: --min-depth must be a positive integer" >&2; exit 2; }

detect_db_level() {
  local db_path="$1"
  local rank=""
  local base
  base=$(basename "$db_path")
  local IFS='.'
  local part
  for part in $base; do
    case "$part" in
      (superkingdom|phylum|class|order|family|genus|species|strain)
        rank="$part"
        break
        ;;
    esac
  done
  if [[ -z "$rank" ]]; then
    rank=$(awk -F '\t' 'NF && $1 !~ /^#/ && tolower($1) != "db_level" && tolower($1) != "level" {print $1; exit}' "$db_path.stats")
  fi
  [[ -n "$rank" ]] || {
    echo "ERROR: unable to detect GOTTCHA2 database level from $db_path or $db_path.stats" >&2
    return 2
  }
  printf '%s\n' "$rank"
}

DB_LEVEL=$(detect_db_level "$DB_PATH") || exit $?

log_start() {
  echo "=> [INFO] $*"
}

log_success() {
  echo "[INFO] $*"
}

log_error() {
  echo "[ERROR] $*"
}

create_placeholder_outputs() {
  local status=0
  local outputs=(
    "$OUTDIR/$PREFIX.krona.html"
    "$OUTDIR/$PREFIX.pathogen.tsv"
    "$OUTDIR/$PREFIX.pathogen.summary.txt"
    "$OUTDIR/$PREFIX.pathogen.full.tsv"
    "$OUTDIR/$PREFIX.pathogen.summary.tsv"
    "$OUTDIR/$PREFIX.pathogen.full.html"
    "$OUTDIR/$PREFIX.gottcha_${DB_LEVEL}.bam"
    "$OUTDIR/$PREFIX.gottcha_${DB_LEVEL}.coverage.tsv"
    "$OUTDIR/$PREFIX.coverage.html"
    "$OUTDIR/$PREFIX.info"
  )

  for outfile in "${outputs[@]}"; do
    touch "$outfile" || status=$?
  done

  if [[ $status -eq 0 ]]; then
    log_success "Placeholder outputs created for downstream steps."
  else
    log_error "Failed to create placeholder outputs (exit $status)."
  fi
  return $status
}

prepare_output_dir() {
  log_start "Preparing output directory..."
  mkdir -p "$OUTDIR"
  local status=$?
  if [[ $status -eq 0 ]]; then
    log_success "Output directory ready: '$OUTDIR'."
  else
    log_error "Unable to create output directory (exit $status)."
  fi
  return $status
}

fastp_qc() {
  log_start "Running fastp pre-processing..."
  INPUT_QC="$OUTDIR/${PREFIX}.qc.fastq.gz"
  if [[ "$PAIRED" == "true" ]]; then
    INPUT_QC_R1="$OUTDIR/${PREFIX}.R1.qc.fastq.gz"
    INPUT_QC_R2="$OUTDIR/${PREFIX}.R2.qc.fastq.gz"
    fastp -i "$READ1" -I "$READ2" -o "$INPUT_QC_R1" -O "$INPUT_QC_R2" -j "$INPUT_QC.json" -h "$INPUT_QC.html"
    local status=$?
    if [[ $status -eq 0 ]]; then
      log_success "Pre-processing completed: '$INPUT_QC_R1' and '$INPUT_QC_R2'."
    else
      log_error "Unable to pre-process input file (exit $status)."
    fi
    return $status
  else
    fastp -i "$INPUT" -o "$INPUT_QC" -j "$INPUT_QC.json" -h "$INPUT_QC.html"
    local status=$?
    if [[ $status -eq 0 ]]; then
      log_success "Pre-processing completed: '$INPUT_QC'."
    else
      log_error "Unable to pre-process input file (exit $status)."
    fi
    return $status
  fi
}

fastplong_qc() {
  log_start "Running fastplong pre-processing..."
  INPUT_QC="$OUTDIR/${PREFIX}.qc.fastq.gz"
  fastplong -i "$INPUT" -o "$INPUT_QC" -l 150 -j "$INPUT_QC.json" -h "$INPUT_QC.html"
  local status=$?
  if [[ $status -eq 0 ]]; then
    log_success "Pre-processing completed: '$INPUT_QC'."
  else
    log_error "Unable to pre-process input file (exit $status)."
  fi
  return $status
}

run_gottcha2() {
  log_start "Running GOTTCHA2 (fast-profile mode)..."
  if [[ "$PAIRED" == "true" ]]; then
      gottcha2 fast-profile \
              -i "$READ1" "$READ2" \
              -t "$CPU" \
              -o "$OUTDIR" \
              -p "$PREFIX" \
              -d "$DB_PATH" \
              -r "GENOMIC_CONTENT_EST" \
              -mf 0.9 \
              -mg 0 \
              --fast-min-kmer 5 \
              --mpa \
              --verbose \
              $ONT_FLAG
  else
      gottcha2 fast-profile \
              -i "$READ" \
              -t "$CPU" \
              -o "$OUTDIR" \
              -p "$PREFIX" \
              -d "$DB_PATH" \
              -r "GENOMIC_CONTENT_EST" \
              -mf 0.9 \
              -mg 0 \
              --fast-min-kmer 5 \
              --mpa \
              --verbose \
              $ONT_FLAG
  fi

  local status=$?
  if [[ $status -eq 0 ]]; then
    touch "$OUTDIR/$PREFIX.tsv"
    touch "$OUTDIR/$PREFIX.full.tsv"
    log_success "GOTTCHA2 run finished: results in '$OUTDIR'."
    if ! grep -q "superkingdom" "$OUTDIR/$PREFIX.tsv"; then
      log_start "No GOTTCHA2 results found. Creating placeholder outputs for downstream steps..."
      create_placeholder_outputs || exit $?
      log_success "Stopping pipeline."
      exit 0
    fi
  else
    log_error "GOTTCHA2 run failed (exit $status)."
  fi
  return $status
}

generate_krona_plot() {
  log_start "Generating Krona plot..."
  local status=0
  if grep -q "^species" "$OUTDIR/$PREFIX.tsv"; then
    grep "^species" "$OUTDIR/$PREFIX.tsv" \
      | ktImportTaxonomy -t 3 -m 9 -o "$OUTDIR/$PREFIX.krona.html" - \
      || status=$?
  else
    printf '<!DOCTYPE html><html><body><p>No species-level GOTTCHA2 rows available for Krona.</p></body></html>\n' > "$OUTDIR/$PREFIX.krona.html" \
      || status=$?
  fi
  if [[ $status -eq 0 ]]; then
    log_success "Krona plot generated: '$OUTDIR/$PREFIX.krona.html'."
  else
    log_error "Krona plot generation failed (exit $status)."
  fi
  return $status
}

annotate_pathogens_qualified() {
  log_start "Annotating pathogens for qualified taxa..."
  touch "$OUTDIR/$PREFIX.pathogen.tsv"
  pathogen.py \
    -i "$OUTDIR/$PREFIX.tsv" \
    -o "$OUTDIR/$PREFIX.pathogen.tsv" \
    -d "$SPADES_DATA/taxonomy_db" \
    -c "$DB_PATH.tax.tsv" \
    -p "$SPADES_DATA/pathogen.tsv" > "$OUTDIR/$PREFIX.pathogen.summary.txt"
  local status=$?
  if [[ $status -eq 0 ]]; then
    log_success "Pathogen annotation completed: '$OUTDIR/$PREFIX.pathogen.tsv'."
    log_success "Pathogen summary saved: '$OUTDIR/$PREFIX.pathogen.summary.txt'."
  else
    log_error "Pathogen annotation for qualified taxa failed (exit $status)."
  fi
  return $status
}

annotate_pathogens_full() {
  log_start "Annotating pathogens for all taxa..."
  touch "$OUTDIR/$PREFIX.pathogen.full.tsv"
  pathogen.py \
    -i "$OUTDIR/$PREFIX.full.tsv" \
    -o "$OUTDIR/$PREFIX.pathogen.full.tsv" \
    -d "$SPADES_DATA/taxonomy_db" \
    -c "$DB_PATH.tax.tsv" \
    -p "$SPADES_DATA/pathogen.tsv"
  local status=$?
  if [[ $status -eq 0 ]]; then
    log_success "Pathogen annotation completed: '$OUTDIR/$PREFIX.pathogen.full.tsv'."
  else
    log_error "Pathogen annotation for all taxa failed (exit $status)."
  fi
  return $status
}

extract_summary_levels() {
  log_start "Extracting genus and species level data for result viewer..."
  local summary="$OUTDIR/$PREFIX.pathogen.summary.tsv"
  local status=0

  if [[ ! -r "$OUTDIR/$PREFIX.pathogen.tsv" ]]; then
    log_error "Pathogen TSV not readable: '$OUTDIR/$PREFIX.pathogen.tsv'."
    return 1
  fi

  head -1 "$OUTDIR/$PREFIX.pathogen.tsv" > "$summary" || status=$?
  if [[ $status -eq 0 ]]; then
    grep "^species" "$OUTDIR/$PREFIX.pathogen.tsv" >> "$summary" || true
    grep "^genus" "$OUTDIR/$PREFIX.pathogen.tsv" >> "$summary" || true
  fi

  if [[ $status -eq 0 ]]; then
    log_success "Genus and species level data extracted: '$summary'."
  else
    log_error "Genus/species data extraction failed (exit $status)."
  fi
  return $status
}

generate_result_viewer() {
  log_start "Generating result viewer HTML..."

  local js_external_flag=""
  if [[ "$JS_EXTERNAL" == "true" ]]; then
    js_external_flag="--external"
  fi

  result_viewer.py \
    -i "$OUTDIR/$PREFIX.pathogen.full.tsv" \
    -o "$OUTDIR/$PREFIX.pathogen.full.html" \
    $js_external_flag
  local status=$?
  if [[ $status -eq 0 ]]; then
    log_success "Result viewer HTML generated: '$OUTDIR/$PREFIX.pathogen.full.html'."
  else
    log_error "Result viewer generation failed (exit $status)."
  fi
  return $status
}

generate_coverage_browser() {
  log_start "Generating coverage browser HTML..."
  local bam="$OUTDIR/$PREFIX.gottcha_${DB_LEVEL}.bam"
  local coverage="$OUTDIR/$PREFIX.gottcha_${DB_LEVEL}.coverage.tsv"
  local ref="$OUTDIR/$PREFIX.sylph_extracted.fa.gz"
  local vcf="$OUTDIR/$PREFIX.gottcha_${DB_LEVEL}.vcf.gz"
  local full_tsv="$OUTDIR/$PREFIX.full.tsv"
  local status=0

  local js_external_flag=""
  if [[ "$JS_EXTERNAL" == "true" ]]; then
    js_external_flag="--external"
  fi

  if [[ $status -eq 0 ]]; then
    samtools coverage "$bam" > "$coverage" || status=$?

    # convert FASTA to bgzip format
    gzip -dc "$ref" | bgzip -@ "$CPU" -c > "$OUTDIR/$PREFIX.sylph_extracted.fa.bgz"
    ref="$OUTDIR/$PREFIX.sylph_extracted.fa.bgz"
    samtools faidx "$ref"

    # Generate VCF
    # --ploidy 1 is usually appropriate for bacterial/viral haploid references.
    bcftools mpileup \
        -Ou \
        -f "$ref" \
        -q 20 \
        -Q 20 \
        -a FORMAT/DP,FORMAT/AD \
        --threads "$CPU" \
        "$bam" | \
    bcftools call \
        -mv \
        --ploidy 1 \
        -Oz \
        --threads "$CPU" \
        -o "$vcf" || status=$?
    
    # index the VCF
    bcftools index -t "$vcf" || status=$?
  fi
  if [[ $status -eq 0 ]]; then
    coverage_browser.py \
      -c "$coverage" \
      -f "$full_tsv" \
      --vcf "$vcf" \
      --min-depth "$MIN_DEPTH" \
      $js_external_flag \
      -o "$OUTDIR/$PREFIX.coverage.html" || status=$?
  fi

  if [[ $status -eq 0 ]]; then
    log_success "Coverage browser HTML generated: '$OUTDIR/$PREFIX.coverage.html'."
  else
    log_error "Coverage browser generation failed (exit $status)."
  fi
  return $status
}

after_run() {
  log_start "Writing pipeline metadata..."
  echo -e "SPADES_version=$VERSION" > "$OUTDIR/$PREFIX.info"
  fastp --version | sed 's/ /=/' >> "$OUTDIR/$PREFIX.info"
  fastplong --version | sed 's/ /=/' >> "$OUTDIR/$PREFIX.info"
  echo "GOTTCHA2_version=$(gottcha2 profile --version)" >> "$OUTDIR/$PREFIX.info"
  echo "GOTTCHA2_database=$DB_PATH" >> "$OUTDIR/$PREFIX.info"
  echo "ONT_mode=$ONT" >> "$OUTDIR/$PREFIX.info"
  echo "ONT_error_rate=$ONT_ERROR_RATE" >> "$OUTDIR/$PREFIX.info"
  samtools --version | awk '
    /^samtools / { print "samtools_version="$2 }
    /^Using htslib / { print "samtools_htslib_version="$3 }' >> "$OUTDIR/$PREFIX.info"
  gzip -dc "$OUTDIR/$PREFIX.gottcha_${DB_LEVEL}.vcf.gz" | grep "##bcftools" | sed 's/##//' | sed 's/; /\n/' >> "$OUTDIR/$PREFIX.info"

  local status=$?
  if [[ $status -eq 0 ]]; then
    log_success "Pipeline metadata recorded: '$OUTDIR/$PREFIX.info'."
  else
    log_error "Failed to record pipeline metadata (exit $status)."
  fi

  if [[ "$CLEAN_FLAG" == "true" ]]; then
    log_start "Cleaning up intermediate files..."
    rm -f \
      "$OUTDIR/$PREFIX.gottcha_${DB_LEVEL}.sam" \
      "$OUTDIR/$PREFIX.gottcha_${DB_LEVEL}.sam.temp" \
      "$OUTDIR/$PREFIX.split_reads.fasta.gz" \
      "$OUTDIR/$PREFIX.sylph_*" \
      "$INPUT_QC" \
      "$INPUT_QC_R1" \
      "$INPUT_QC_R2"

    log_success "Cleanup completed."
  fi

  return $status
}

run_pipeline() {
  log_start "Starting SPADES workflow..."

  ONT_FLAG=""

  if [[ "$PAIRED" == "true" ]]; then
    READ1="$READ1"
    READ2="$READ2"
  else
    READ="$INPUT"
  fi

  prepare_output_dir || exit $?

  if [[ "$ONT" == "true" ]]; then
    fastplong_qc || exit $?
    ONT_FLAG="-np -er $ONT_ERROR_RATE"
  else
    fastp_qc || exit $?
  fi

  if [[ "$PAIRED" == "true" ]]; then
    READ1="$INPUT_QC_R1"
    READ2="$INPUT_QC_R2"
  else
    READ="$INPUT_QC"
  fi

  run_gottcha2 || exit $?
  generate_krona_plot || exit $?
  annotate_pathogens_qualified || exit $?
  annotate_pathogens_full || exit $?
  extract_summary_levels || exit $?
  generate_result_viewer || exit $?
  generate_coverage_browser || exit $?
  after_run || exit $?

  log_success "Pipeline completed successfully. All results are in '$OUTDIR'."
}

run_pipeline