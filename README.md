[![logo](docs/SPADES_logo_icon.png)](docs/SPADES_logo_icon.png)

# SPADES

The SPADES (Sequence-based Pathogen-Agnostic Diagnostics/Detection Solution) program is a full suite of  bioinformatics tools for detecting pathogens in complex metagenomic DNA sequencing datasets. The output of SPADES includes the detected organisms, highlighting the possible pathogens identified, with coverage-based validation. Under the hood, SPADES utilizes the GOTTCHA2 taxonomy profiling bioinformatics algorithms and databases, extracts pathogen hits, computes detailed signature/genome coverage summaries, and automatically generates an interactive report/plot for rapid interpretation and triage.

## Requirements

- Conda or mamba (recommended)
- Dependencies in `environment.yml` (python 3.12, minimap2, krona, samtools, etc.)
- A GOTTCHA2 database (`.syldb`, `.zip`, `.stats`, `.tax.tsv` for read input; `.stats` and `.tax.tsv` for existing-BAM input)
- Taxonomy and pathogen data (provided in `data/`)
- Minimal hardware: 4 CPU cores, 24 GB RAM, 100 GB disk space

## Setup

### Local Conda

```bash
mamba env create -f environment.yml
mamba activate spades
./post_install.sh
./scripts/check_install.sh
```

`post_install.sh` links the Krona taxonomy folder to `data/taxonomy_db`.
`scripts/check_install.sh` verifies command line tools, Python imports, bundled
data files, and Krona taxonomy wiring without running the sequencing workflow.

### Docker

```bash
docker build -t spades-g2:latest .
docker run --rm -it spades-g2:latest ./scripts/check_install.sh
```

See [Docker.md](Docker.md) for mounted input/database examples and multi-platform
build commands.

### Database Requirements

> **Warning:** If you downloaded the database before October 27, 2026 (`20261027`), we recommend downloading an updated version from the [SPADES database download page](https://ref-db.edgebioinformatics.org/spades-g2/).

The [**SPADES database**](https://ref-db.edgebioinformatics.org/spades-g2/) is available for download [here](https://ref-db.edgebioinformatics.org/spades-g2/). The script to start SPADES, `run_SPADES.sh`, uses `gottcha2 fast-profile`. The database path passed with
`-d/--db-path` accepts the database base path or a path ending in `.syldb`,
`.zip`, `.stats`, or `.tax.tsv`. For read input, these files must exist next
to the base path:

- `<db>.syldb`
- `<db>.zip`
- `<db>.stats`
- `<db>.tax.tsv`

Existing-BAM mode requires only `<db>.stats` and `<db>.tax.tsv`. The database
level is detected from the base filename or the `.stats` file.



## Usage

```bash
./run_SPADES.sh \
  -i reads.fastq.gz \
  -o outdir \
  -p sample \
  -d /path/to/gottcha_db.species.fna \
  -t 4 \
  [--spades-data data/] \
  [--ont] \
  [--ont-error-rate 0.03] \
  [--min-depth 10] \
  [--skip-qc] \
  [--auto-skip-qc-ratio 0.2] \
  [--clean]
```

For paired-end FASTQ input, replace `-i reads.fastq.gz` with
`-1 reads_R1.fastq.gz -2 reads_R2.fastq.gz`. Both mates are required;
paired-end input cannot be combined with `--ont`.

Re-profile an existing coordinate-sorted and indexed GOTTCHA2 BAM without
mapping the reads again:

```bash
./run_SPADES.sh \
  --bam sample.gottcha_species.bam \
  -o reprofiled \
  -p sample_reprofiled \
  -d /path/to/gottcha_db.species.fna \
  -t 4 \
  --ont
```

The BAM index must be named `<bam>.bai`. Existing-BAM mode runs `gottcha2
profile --bam`, regenerates the taxonomic/pathogen reports, and creates a
coverage browser without re-mapping reads. Read QC is bypassed in this mode.

For detailed instructions on using GOTTCHA2, please see its [wiki page](https://github.com/poeli/GOTTCHA2/wiki).

Options:

- `-i, --input` Single-end FASTA/FASTQ reads file (optionally gzipped)
- `-1, --read1` Paired-end read 1 (fastq/fq[.gz], etc.)
- `-2, --read2` Paired-end read 2 (fastq/fq[.gz], etc.)
- `-b, --bam` Existing coordinate-sorted GOTTCHA2 BAM; mutually exclusive with read input
- `-o, --outdir` Output directory
- `-p, --prefix` Output prefix (base filename)
- `-d, --db-path` GOTTCHA2 database base path or one of its companion files
- `-t, --cpu` Threads/CPUs (positive integer)
- `--spades-data` Directory containing `taxonomy_db/` and `pathogen.tsv` (default: `data/`)
- `--ont` Treat reads or BAM as ONT data; uses `fastplong` for FASTQ QC (unless skipped) and passes ONT profiling flags to GOTTCHA2
- `--ont-error-rate` Error rate for ONT reads passed to GOTTCHA2 (`-er`), default: `0.03`
- `--clean` Remove intermediate files after the run (see Outputs section for details)
- `--min-depth` Minimum depth for variant calling (default: 10)
- `--skip-qc` Skip fastp/fastplong quality control and preprocessing steps
- `--auto-skip-qc-ratio` Use original reads downstream when QC retains less than this fraction (decimal from `0` to `1`, default: `0.2`)
- `--js-external` Deprecated compatibility option; HTML resources are always embedded
- `--version` Show script version and exit
- `-h, --help` Show usage

Single-end FASTA input automatically bypasses FASTQ quality filtering. When QC
runs, its reports are retained even if `--auto-skip-qc-ratio` causes the pipeline
to use the original reads. Set the ratio to `0` to disable this automatic
fallback, or use `--skip-qc` to bypass QC entirely.

## Outputs

The pipeline keeps final tables in the main output directory, HTML reports in
`results_html/`, and supporting files in `intermediate/`.

### Final Outputs (in `<outdir>/`)

- `<prefix>.pathogen.tsv` Qualified GOTTCHA2 hits with pathogen annotations
- `<prefix>.pathogen.summary.txt` Pathogen detection summary text log
- `<prefix>.taxonomy.tsv` All qualified GOTTCHA2 taxa with pathogen annotations where available
- `<prefix>.taxonomy.unfiltered.tsv` All detected taxa, including those failing GOTTCHA2 filters, with pathogen annotations where available

### HTML Reports (in `<outdir>/results_html/`)

- `<prefix>.results_viewer.html` Interactive viewer of the unfiltered taxonomy table
- `<prefix>.results_krona.html` Interactive Krona taxonomic plot
- `<prefix>.coverage_browser.html` Interactive coverage browser with variant information
- `<prefix>.reads_qc.html` FastP/FastPlong QC report (when QC runs)

### Intermediate Files (in `<outdir>/intermediate/`)

- `<prefix>.info` Pipeline metadata (versions, database path, input mode, ONT settings)
- `<prefix>.pathogen.summary.tsv` Genus/species level pathogen data
- `<prefix>.tsv` GOTTCHA2 profiling result table (filtered)
- `<prefix>.full.tsv` GOTTCHA2 full profiling output
- `<prefix>.gottcha_<dbLevel>.bam` Alignment file (with index)
- `<prefix>.gottcha_<dbLevel>.sam` SAM alignment (if produced)
- `<prefix>.gottcha_<dbLevel>.vcf.gz` Variant calls (with index, when produced)
- `<prefix>.gottcha_<dbLevel>.coverage.tsv` Coverage summary table
- `<prefix>.qc.fastq.gz` (or `R1/R2`) Quality-controlled reads
- `<prefix>.qc.fastq.gz.json` FastP/FastPlong JSON report
- `<prefix>.mpa.tsv` MetaPhlAn-compatible output
- `<prefix>.lineage.tsv` Lineage information
- `<prefix>.sylph_*` Sylph extraction files
- `<prefix>.sylph_extracted.fa.bgz` Extracted reference FASTA used for alignment and coverage analysis (with indexes when produced)

Existing-BAM mode uses the supplied BAM and its index in place; it does not
copy them into the final output directory.

### Cleanup Behavior

- **Default** (no `--clean`): All files are preserved, with final tables in `<outdir>/`, HTML reports in `<outdir>/results_html/`, and supporting files in `<outdir>/intermediate/`.
- **With `--clean`**: The `intermediate/` directory and all its contents, including metadata and QC JSON, are deleted. Final tables and `results_html/` are retained.

If GOTTCHA2 produces no results, the script creates empty placeholder tables,
HTML reports, metadata, and (for read input) a BAM, then exits cleanly before
organizing the output directory. These placeholders remain directly in
`<outdir>/`; QC and other files already generated are also left in place.

### GOTTCHA2 Full Report

See [tutorial.md](GOTTCHA2-main/docs/tutorial.md) for details about GOTTCHA2 results.

## Example

```bash
cd test
./test_run.sh
```

This runs ONT, single-end, and paired-end examples using files under `test/`
and the bundled data in `data/`. Output can be found in the generated
directories under `test/`.

## Streaming Oxford Nanopore directories

Run `stream_spades.py` to continuously monitor a directory for newly generated, stable FASTA/FASTQ files (`.fa`, `.fasta`, `.fna`, `.fq`, or `.fastq`, optionally gzip-compressed) and automatically process them using the SPADES workflow.

Open `<outdir>/<prefix>.stream.html` in your web browser to monitor the workflow in real time. The report is a self-contained dashboard that automatically refreshes as new data and analysis results become available.

```bash
./stream_spades.py \
  --input-dir /data/ont/fastq_pass \
  --outdir /results/run_01_stream \
  --prefix run_01 \
  --db-path /db/gottcha_db.species.fna \
  --cpu 8
```

The monitor waits until a file's size and modification time have remained
unchanged for `--settle-seconds` (30 seconds by default). Subdirectories are
included by default; use `--no-recursive` to watch only the input directory,
`--once --settle-seconds 0` to process an existing batch and exit, or
`--max-files N` to stop after a fixed number of files.
When multiple stable files are present, they are processed strictly one at a
time in modification-time/name order. The next file does not start until the
current file's cumulative profile and state record are complete.

Use `--skip-qc` to bypass fastp/fastplong quality control and preprocessing for
each arriving file. In this mode, `stream_spades.py` counts records directly
from plain or gzip-compressed FASTA/FASTQ input, so the **Input reads** metric
remains available in the dashboard. Because no QC was performed, **Post-QC
reads** is shown as `—` instead of repeating the input count.

The monitor also accepts `--auto-skip-qc-ratio` (default: `0.2`) and forwards it
to `run_SPADES.sh`. If this fallback uses original reads, **Post-QC reads**
still reports the number retained by QC, rather than the number used downstream.

Per-batch working files are removed after successful processing by default.
Use `--keep-batches` when chunk and staged profile files are needed for
troubleshooting.

Streaming output is organized as follows:

```text
<outdir>/
  cumulative/<prefix>.gottcha_<level>.bam[.bai]
  cumulative/<prefix>.sylph_extracted.fa.bgz
  cumulative/<prefix>.pathogen.tsv
  cumulative/<prefix>.pathogen.summary.txt
  cumulative/<prefix>.taxonomy.tsv
  cumulative/<prefix>.taxonomy.unfiltered.tsv
  cumulative/results_html/<prefix>.results_viewer.html
  cumulative/results_html/<prefix>.results_krona.html
  cumulative/results_html/<prefix>.coverage_browser.html
  cumulative/<prefix>.gottcha_<level>.coverage.tsv
  cumulative/intermediate/<prefix>.info
  cumulative/...       # other files from the latest cumulative profile
  <prefix>.stream.html
  stream_state.json
  batches.tsv
  logs/
```

With `--keep-batches`, each batch also retains
`batches/batch_000001/chunk/` and
`batches/batch_000001/profile/`, including chunk QC reports. An alignment-free
run can publish the pipeline's root-level placeholders under `cumulative/`
instead of `cumulative/results_html/`.

`stream_state.json` makes restarts idempotent: an unchanged input file is not
processed twice. If a file at the same path changes size or modification time,
it is treated as a new timepoint. The cumulative BAM is replaced atomically
after the merged BAM and final report both validate.

`batches.tsv` records `raw_reads`, `filtered_reads`, and their cumulative
totals for every timepoint. When `--skip-qc` is active, `raw_reads` contains the
directly counted input records while `filtered_reads` and its cumulative total
are left empty. With `--keep-batches`, the manifest also records the retained
`qc_json` path when QC was run. `<prefix>.stream.html` is a self-contained,
auto-refreshing status dashboard. It shows the current and overall monitor
state, the latest cumulative screening finding, cumulative read totals, recent
timepoints, and one tile per identified human-pathogenic species. The footer
shows the version reported by `run_SPADES.sh --version`.
The latest-results section links directly to the cumulative
`cumulative/results_html/<prefix>.results_viewer.html` report and
`cumulative/results_html/<prefix>.coverage_browser.html` coverage report.
For compatibility with existing state files and dashboard readers, the manifest
and state keys `pathogen_full_tsv`, `pathogen_full_html`, and `coverage_html`
retain their names while pointing to the new filenames.
Current findings are visually distinguished from historical-only findings.
Clicking a species expands three line charts tracking supporting read
alignments, signature coverage (`SIG_COV`), and SNI score over time. The report
labels these as analytical screening evidence and keeps interpretation and
confirmation limitations visible for clinical users.

## Notes

- The script exports its own `scripts/` folder onto `PATH` and defaults `--spades-data` to `data/` next to the script.
- `stream_spades.py` always enables `--ont`; no separate `--ont` option is needed for streaming.
- `--clean` removes the `intermediate/` directory after a successful run, keeping final tables and `results_html/`. Without `--clean`, all intermediate files are preserved for debugging.
- `--min-depth` sets the minimum depth threshold for variant calling in the coverage browser (default: 10).
- Ensure the database path points to the base file (e.g. `gottcha_db.species.fna`) and that the corresponding `.syldb`, `.zip`, `.stats`, and `.tax.tsv` exist.

## Notice of Copyright Assertion (O5123)

This program is Open-Source under the BSD-3 License.
 
Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 
Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 
Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 
Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
