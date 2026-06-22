[![logo](https://gottcha2.readthedocs.io/en/latest/_images/gottcha_icon.png)](https://gottcha2.readthedocs.io/en/latest/_images/gottcha_icon.png)

# Genomic Origin Through Taxonomic CHAllenge (GOTTCHA2)

GOTTCHA2 is a gene-independent, signature-based metagenomic taxonomic profiler for sequencing reads. It is designed to reduce false discoveries while remaining practical to run on a workstation or laptop. Instead of relying on marker genes, it maps reads to precomputed unique signature fragments and estimates abundance from signature coverage and depth.

> GOTTCHA v1 databases are not compatible with GOTTCHA2.

---

## Table of contents

- [What's new](#whats-new)
- [Installation](#installation)
- [Dependencies](#dependencies)
- [Databases](#databases)
- [Quick start](#quick-start)
- [Command overview](#command-overview)
- [Profiling](#profiling)
- [Fast profile mode](#fast-profile-mode)
- [Read extraction](#read-extraction)
- [Output files](#output-files)
- [Thresholds and filtering](#thresholds-and-filtering)
- [Full report fields](#full-report-fields)
- [Troubleshooting](#troubleshooting)
- [License and citation](#license-and-citation)

---

## What's new

Recent GOTTCHA2 releases through v2.4.0 include several workflow changes that are worth knowing before you start:

- **Fast prefiltering mode**: `fast-profile` uses `sylph` to prefilter the reference set before read mapping while producing results comparable to the standard `profile` workflow. Depending on the sample and database, it often reduces runtime by about 5–10× and memory usage by roughly 2–10×.
- **Current CLI**: the supported entry points are `profile`, `fast-profile`, `extract`, `sam2bam`, `download`, and `version`.
- **Updated identity handling**: the reported `SNI_SCORE` is based on consensus identity rather than the legacy read-weighted identity metric.
- **BAM-based workflow**: runs use sorted and indexed BAM for downstream processing instead of keeping SAM as the main intermediate.
- **Legacy compatibility**: the older `gottcha2.py` workflow (SAM-based) is still available for compatibility, but it is frozen at v2.2.3.

---

## Installation

### Option A: Conda

```bash
conda install -c bioconda gottcha2
```

### Option B: Install from source

Install the external tools first, then install the Python package:

```bash
# required for profile
# add sylph as well if you plan to use fast-profile
git clone https://github.com/poeli/GOTTCHA2
cd GOTTCHA2
python -m pip install .

# development install
python -m pip install -e .
```

Confirm the installation:

```bash
gottcha2 version
gottcha2 profile --help
gottcha2 fast-profile --help
```

For containerized usage, see [DOCKER.md](../DOCKER.md).

---

## Dependencies

GOTTCHA2 requires Python 3.9+.

Runtime dependencies:

- `minimap2` for mapping
- `samtools` and `pysam` for BAM conversion and parsing
- `numpy` and `pandas`
- `requests`
- `tqdm`
- `biom-format` if you use `--format biom`
- `sylph` if you use `fast-profile`

A Conda environment file is provided as `environment.yml`.

---

## Databases

### Prebuilt databases

The default download target used by `gottcha2 download` is:

```text
https://ref-db.edgebioinformatics.org/gottcha2/latest/gottcha_db.species.tar
```

You can also download database bundles manually from the same host if you prefer.

### Database bundle contents

A standard profiling database should include these files with the same prefix:

- `gottcha_db.<level>.fna.mmi` for `profile`
- `gottcha_db.<level>.fna.tax.tsv` taxonomy mapping
- `gottcha_db.<level>.fna.stats` signature and genome statistics

Additional files used by fast mode:

- `gottcha_db.<level>.fna.syldb` `sylph` database for prefiltering
- `gottcha_db.<level>.fna.zip` archived signature sequences used to build the reduced reference

Pass either the shared database prefix or the database directory to `-d/--database`. GOTTCHA2 locates the required sidecar files from that path. For example:

```text
/path/to/db/gottcha_db.species.fna
```

or

```text
/path/to/db
```

### Download helper

[Not available yet] Use the built-in downloader to fetch the default database tarball into a new `database/` directory:

```bash
gottcha2 download
```

See available options with:

```bash
gottcha2 download --help
```

---

## Quick start

These examples use the most common workflows. Replace the example paths and filenames with your own sample and database locations.

### Example conventions

- `gottcha2 profile` maps reads to the selected signature database and produces taxonomic reports.
- `gottcha2 fast-profile` first narrow the reference set, then runs the profiling workflow on the reduced reference.
- `gottcha2 extract` pulls reads assigned to selected taxa from an existing GOTTCHA2 BAM file.
- `-d/--database` points to either a database prefix, such as `/path/to/db/gottcha_db.species.fna`, or to a directory that contains the matching database files.
- `-i/--input` supplies one or more read files. Use two files for paired-end Illumina reads and one file for single-end or Nanopore reads.
- `-b/--bam` reuses an existing sorted and indexed BAM instead of remapping reads.
- `-t/--threads` controls the number of CPU threads used by mapping and related processing steps.
- `-o/--outdir` chooses the output directory. GOTTCHA2 creates it if needed.
- `-p/--prefix` sets the output filename prefix. If you omit it, GOTTCHA2 derives a prefix from the input filename or BAM name.
- The backslash (`\`) at the end of a line lets long shell commands continue on the next line. You can also write each example as a single line.

### 1) Profile Illumina paired-end reads

Use this when your sample has forward and reverse FASTQ files.

```bash
gottcha2 profile \
  -d /path/to/db/gottcha_db.species.fna \
  -i sample_R1.fastq.gz sample_R2.fastq.gz \
  -t 8 \
  -o out \
  -p sample
```

What this command does:

- loads the species-level database specified with `-d`
- maps both paired-end read files supplied after `-i`
- uses 8 threads because of `-t 8`
- writes results into the `out/` directory
- names output files with the `sample` prefix, for example `sample.tsv`, `sample.full.tsv`, and `sample.gottcha_species.bam`

Use a sample-specific prefix whenever you process multiple samples into the same output directory.

### 2) Profile Illumina single-end reads

Use this when each sample has one FASTQ file.

```bash
gottcha2 profile \
  -d /path/to/db/gottcha_db.species.fna \
  -i sample.fastq.gz \
  -t 8 \
  -o out
```

Because `-p/--prefix` is omitted, GOTTCHA2 derives the output prefix from `sample.fastq.gz`. Add `-p sample_name` if you want a shorter or more explicit prefix.

### 3) Profile Oxford Nanopore reads

Nanopore mode expects exactly one input file. Add `--nanopore` so GOTTCHA2 uses long-read preprocessing and Nanopore-oriented default thresholds.

```bash
gottcha2 profile \
  -d /path/to/db/gottcha_db.species.fna \
  -i ont_reads.fastq.gz \
  --nanopore \
  -t 8 \
  -o out \
  -p ont_sample
```

The important difference from short-read mode is `--nanopore`. In this mode, GOTTCHA2 chunks long reads before mapping and uses more permissive default alignment and error-rate settings. See [Oxford Nanopore mode](#oxford-nanopore-mode) for details.

### 4) Re-run profiling from an existing BAM

Use this when you already have a sorted and indexed GOTTCHA2 BAM and want to re-aggregate results with different thresholds. This avoids the slower read-mapping step.

```bash
gottcha2 profile \
  -b sample.gottcha_species.bam \
  -d /path/to/db/gottcha_db.species.fna \
  -Mc 0.01 \
  -Mr 10 \
  -mi 0.95 \
  -t 8 \
  -o out \
  -p sample.refiltered
```

What the non-default options mean:

- `-b sample.gottcha_species.bam` reads alignments from an existing BAM instead of using `-i` input reads.
- `-Mc 0.01` requires at least 1% signature coverage for abundance calculation.
- `-Mr 10` requires at least 10 mapped reads.
- `-mi 0.95` keeps only matches with at least 95% alignment identity.
- `-p sample.refiltered` keeps the re-filtered output separate from the original run.

The BAM must be coordinate-sorted and indexed. Keep the database path consistent with the database used for the original mapping.

### 5) Run the faster prefiltering workflow

Use `fast-profile` when you want to reduce runtime and memory usage while producing results comparable to the standard `profile` workflow. It does this by preselecting likely reference sequences before read mapping.

```bash
gottcha2 fast-profile \
  -d /path/to/db/gottcha_db.species.fna \
  -i sample.fastq.gz \
  -t 8 \
  -o out \
  -p sample.fast
```

This command requires the standard database files plus fast-mode sidecars: `.syldb` and `.zip`. The outputs have the same general structure as `profile`, but mapping is performed against a reduced reference set selected by `sylph`.

### 6) Extract reads for a taxon from an existing BAM

Use `extract` after profiling when you want the reads assigned to one or more taxa. The example below extracts reads assigned to NCBI taxid `562`.

```bash
gottcha2 extract \
  -b sample.gottcha_species.bam \
  -e 562 \
  -o out \
  -p sample.ecoli
```

Key options:

- `-b` points to the GOTTCHA2 BAM created by `profile` or `fast-profile`.
- `-e` selects the taxon or taxa to extract. You can use taxids, names, or `@file` syntax.
- `-o` and `-p` control where the extracted FASTA or FASTQ output is written.

### Check the results

After any profiling run, start with these files:

```bash
ls out
less out/sample.gottcha_species.log
column -t -s $'\t' out/sample.tsv | less -S
column -t -s $'\t' out/sample.full.tsv | less -S
```

The summary report (`*.tsv`, `*.csv`, or `*.biom`) contains taxa that passed the selected filters. The full report (`*.full.tsv`) includes both passing and filtered taxa and records filtering reasons in the `NOTE` column.

---

## Command overview

GOTTCHA2 uses a subcommand-style CLI:

```text
gottcha2 <command> [options]
```

| Command | Use it when you need to... | Typical starting point |
| ------- | -------------------------- | ---------------------- |
| `profile` | Map reads or reuse a BAM and generate taxonomic profiles. | `gottcha2 profile -d DB -i reads.fastq.gz -o out` |
| `fast-profile` | Prefilter the database with `sylph`, then run profiling on a reduced reference set. | `gottcha2 fast-profile -d DB -i reads.fastq.gz -o out` |
| `extract` | Extract reads assigned to one or more taxa from an existing BAM. | `gottcha2 extract -b sample.bam -e 562` |
| `sam2bam` | Convert legacy GOTTCHA2 SAM output into sorted, indexed BAM. | `gottcha2 sam2bam -i sample.sam -o sample.bam` |
| `download` | Download the default database bundle, when supported by your build. | `gottcha2 download` |
| `version` | Print the installed GOTTCHA2 version. | `gottcha2 version` |

Use `--help` after any command to see command-specific options, defaults, and examples:

```bash
gottcha2 profile --help
gottcha2 extract --help
```

---

## Profiling

### Key concepts

GOTTCHA2 profiles metagenomic samples by mapping sequencing reads directly to taxon-specific signature fragments. GOTTCHA2 consolidates alignments across each genome's signature space to compute coverage and depth statistics, then derives an ANI-like metric called the signature nucleotide identity score (`SNI_SCORE`). Genome-level results are subsequently aggregated to higher taxonomic ranks.

### Oxford Nanopore mode

With `--nanopore`, GOTTCHA2 first converts the input read file into a temporary FASTA of fixed-length chunks to make long reads easier to map. In the current implementation, Nanopore mode:

- requires exactly one input file
- splits reads into non-overlapping 150 bp chunks
- drops the trailing remainder if it is shorter than 150 bp
- removes inconsistent chunk assignments after mapping

If you do not override them explicitly, Nanopore mode uses these defaults:

- `--matchIdentity 0.85`
- `--matchFraction 0.85`
- `--matchLength 100`
- `--errorRate 0.03`

### Accessions of interest

Use `--accList` to provide a text file containing one accession or signature ID per line. This is useful for plasmids, spike-ins, or other targets you want to track during profiling.

Use `--accListAction` to control how those reads are handled:

- `report_only` keeps all reads and reports the count in `AOI_READ_COUNT`
- `filter_out` removes reads matching listed accessions
- `filter_in` keeps only reads matching listed accessions

### Reporting level and database level

The database level is usually auto-detected from the database prefix or BAM name. For example, `gottcha_db.species.fna` implies `species`, and `sample.gottcha_species.bam` implies `species`.

If auto-detection is not possible, set it explicitly with `-l/--dbLevel`.

---

## Fast profile mode

`fast-profile` is a convenience wrapper for `profile --fast`. It adds a `sylph` prefiltering step before read mapping:

1. Query the `.syldb` database against the input sample.
2. Collect the subset of candidate signatures.
3. Extract those signatures from the `.zip` archive.
4. Map reads only against that reduced reference.

This mode is useful when you need faster execution with a smaller memory footprint. It still produces the standard GOTTCHA2 outputs, including the BAM and summary reports.

---

## Read extraction

GOTTCHA2 can extract reads for one or more taxa from an existing BAM file. Taxa can be provided as:

- comma-separated taxids, for example `-e "666,562"`
- comma-separated taxon names, for example `-e "Vibrio cholerae,Escherichia coli"`
- a file prefixed with `@`, for example `-e "@taxids.txt"`

The `extract` command is shorthand for running `profile` with `--extract` and `--extractOnly`.

### Example usages

Extract reads mapping to taxid `666`:

```bash
gottcha2 extract \
  -b sample.gottcha_species.bam \
  -e 666
```

Extract with explicit match thresholds:

```bash
gottcha2 extract \
  -b sample.gottcha_species.bam \
  -e 666 \
  -mi 0.9 \
  -mf 0.9
```

Extract multiple taxa:

```bash
gottcha2 extract -b sample.gottcha_species.bam -e "1234,5678"
gottcha2 extract -b sample.gottcha_species.bam -e "@taxids.txt"
```

Limit the number of reads per taxon and choose the output format with `:N:FORMAT`:

```bash
# up to 1000 reads per taxon, FASTQ output
gottcha2 extract -b sample.gottcha_species.bam -e "@taxids.txt:1000:fastq"
```

Extract up to 20 representative sequences per profiled reference:

```bash
gottcha2 extract -b sample.gottcha_species.bam -ef
```

### Extracted record format

Each extracted FASTA or FASTQ header encodes the matched reference, interval, taxon, and match statistics:

```text
>{READ_NAME}{MATE}|{REFERENCE}:{START}..{END} LEVEL={LEVEL} NAME={NAME} TAXID={TAXID} AOI={AOI} MG={MG} MI={MI} MF={MF}
```

Field definitions:

- `READ_NAME`: read identifier
- `MATE`: paired-end suffix (`.1`, `.2`, or empty)
- `REFERENCE`: matched reference sequence name
- `START..END`: mapped reference interval (1-based)
- `LEVEL`: extracted taxonomic rank
- `NAME`: extracted taxon name
- `TAXID`: extracted taxonomy ID
- `AOI`: accession-of-interest flag
- `MG`: alignment length
- `MI`: mapping identity
- `MF`: mapping fraction

Example:

```text
>read123.1|chrA|1|300|GCF10000:10..120 LEVEL=species NAME=Escherichia_coli TAXID=562 AOI=False MG=148 MI=98.65 MF=0.99
ACGT...
```

---

## Output files

By default, outputs go to `--outdir` and use a prefix derived from `--prefix`, the first input filename, or the BAM name.

Typical outputs:

- `*.tsv`, `*.csv`, or `*.biom` - summary report at the requested reporting level
- `*.full.tsv` - full report including filtered taxa and notes
- `*.lineage.tsv` - lineage table for qualified taxa
- `*.mpa.tsv` - MetaPhlAn-style output when `--mpa` is enabled
- `*.extract.fasta` or `*.extract.fastq` - extracted reads when `--extract` or `extract` is used
- `*.gottcha_<level>.bam` and `*.bai` - sorted BAM and index for reuse
- `*.gottcha_<level>.log` - run log including thresholds and processing steps

---

## Thresholds and filtering

Most taxonomic cutoffs default to `0` and are disabled unless you set them explicitly. Alignment thresholds are applied by default unless you lower them yourself.

Use `--noCutoff` to disable taxonomic profiling cutoffs. This is equivalent to:

```text
-Mc 0 -Mr 0 -Ml 0 -Mz 0 -ss 0,0,0
```

### Alignment thresholds

- `-mi, --matchIdentity <FLOAT>`
  Minimum alignment identity for a valid match. Default: `0.95` for short reads, `0.85` for Nanopore mode.

- `-mf, --matchFraction <FLOAT>`
  Minimum aligned fraction for a valid match. Default: `0.95` for short reads, `0.85` for Nanopore mode.

- `-mg, --matchLength <INT>`
  Minimum alignment length in bp. Default: `100`.

### Taxonomic profiling cutoffs

- `-er, --errorRate <FLOAT>`
  Estimated sequencing error rate. Default: `0.005` for short reads, `0.03` for Nanopore mode.

- `-ss, --sniScore <FLOAT>[,<FLOAT>,<FLOAT>]`
  SNI-score thresholds for `other,species,strain`. Default: `0.9,0.95,0.99`.

- `-Mc, --minCov <FLOAT>`
  Minimum signature coverage required for abundance calculation. Default: `0`.

- `-Mr, --minReads <INT>`
  Minimum number of mapped reads. Default: `0`.

- `-Ml, --minLen <INT>`
  Minimum covered signature length. Default: `0`.

- `-Mz, --maxZscore <FLOAT>`
  Maximum z-score for mapped-region depth distribution. Default: `0` (disabled).

Filtered taxa remain visible in `*.full.tsv`, with the reason recorded in `NOTE`.

---

## Full report fields

The full report (`<prefix>.full.tsv`) contains all computed metrics. The summary report contains the qualified rows shown at the requested reporting level.

| Field Name             | Description |
| ---------------------- | ----------- |
| LEVEL                  | Taxonomic rank (`superkingdom` through `strain`) |
| NAME                   | Taxon name |
| TAXID                  | NCBI taxonomy ID |
| READ_COUNT             | Reads mapped to this taxon |
| TOTAL_BP_MAPPED        | Total mapped bases across this taxon's signatures |
| SNI_SCORE              | Signature nucleotide identity used during filtering and aggregation |
| COVERED_SIG_LEN        | Total covered signature length |
| BEST_SIG_COV           | Highest signature coverage among rolled-up members |
| DEPTH                  | Depth of coverage (`TOTAL_BP_MAPPED / TOTAL_SIG_LEN`) |
| REL_ABUNDANCE_GC       | Relative abundance from genomic-content estimate |
| REL_ABUNDANCE          | Relative abundance from the field selected by `--relAbu` |
| PARENT_NAME            | Parent taxon name |
| PARENT_TAXID           | Parent taxonomy ID |
| AOI_READ_COUNT         | Reads matched to `--accList` entries |
| TOTAL_READ_LEN         | Total aligned read length |
| TOTAL_BP_MISMATCH      | Total mismatched bases |
| TOTAL_BP_INDEL         | Total inserted and deleted bases |
| READ_WT_SNI            | Read-weighted identity estimate |
| CONSENSUS_SEQ_SNI      | Consensus-sequence identity estimate |
| SNI_CI95_LH            | Low and high 95% confidence bounds for identity |
| SIG_COV                | Signature coverage (`COVERED_SIG_LEN / TOTAL_SIG_LEN`) |
| MAPPED_SIG_LEN         | Signature length with at least one mapped read |
| TOTAL_SIG_LEN          | Total signature length for the taxon |
| COVERED_SIG_DEPTH      | Depth across covered signature only |
| COVERED_MAPPED_SIG_COV | Covered fraction of mapped signature |
| ZSCORE                 | Depth-distribution z-score |
| GENOMIC_CONTENT_EST    | Genomic-content estimate |
| ABUNDANCE              | Raw abundance value from `--relAbu` |
| REL_ABUNDANCE_DEPTH    | Relative abundance computed from depth |
| SIG_LEVEL              | Signature rank used for mapping |
| GENOME_COUNT           | Number of rolled-up genomes |
| GENOME_SIZE            | Combined genome size used for GC normalization |
| NOTE                   | Filtering or rollup note |

---

## Troubleshooting

### BAM input must be sorted and indexed

If you provide `-b/--bam`, the BAM must already be coordinate-sorted and indexed.

For legacy GOTTCHA2 SAM output, convert it with:

```bash
gottcha2 sam2bam -i sample.sam -o sample.bam -t 8
```

### Database sidecar files are required

For `profile`, keep the database sidecar files next to the database prefix. At minimum, GOTTCHA2 expects:

```text
<db>.mmi
<db>.tax.tsv
<db>.stats
```

For `fast-profile`, it expects:

```text
<db>.syldb
<db>.zip
<db>.tax.tsv
<db>.stats
```

### `--nanopore` requires one input file

Nanopore mode only accepts a single FASTA or FASTQ input file. If you have multiple files, merge them first or process them separately.

### Python and external dependency checks happen at runtime

GOTTCHA2 checks for:

- Python 3.9+
- `minimap2`
- `samtools`
- `sylph` when `fast-profile` is used

If one of these tools is missing from `PATH`, the run stops before mapping begins. Install the missing tool or activate the environment that contains it, then rerun the command.

### No taxa are reported

If the summary report is empty, check the full report and the log before rerunning:

- `*.full.tsv` shows filtered taxa and the reason in `NOTE`.
- `*.gottcha_<level>.log` records the thresholds and database files used.
- Lowering `-Mc`, `-Mr`, `-mi`, or `-mf` can increase sensitivity, but may also increase false positives.

### Fast profile cannot find `.syldb` or `.zip`

`fast-profile` requires the standard database sidecar files plus the fast-mode `.syldb` and `.zip` files. If those files are missing, either download a fast-mode-compatible database bundle or use `profile` with the standard `.mmi` database.

### Identity and SNI changed from older releases

Modern GOTTCHA2 releases report `SNI_SCORE` from consensus identity. If you compare output against older `gottcha2.py` runs, expect differences in SNI-related columns and filtering behavior.

---

## License and citation

- License: [TBD]
- If you use GOTTCHA2 in publications, cite the GOTTCHA or GOTTCHA2 project, the database source, and the exact software version reported by `gottcha2 version`.
