# SPADES

The SPADES (Sequence-based Pathogen-Agnostic Diagnostics/Detection Solution) program is a full suite of  bioinformatics tools for detecting pathogens in complex metagenomic DNA sequencing datasets. The output of SPADES includes the detected organisms, highlighting the possible pathogens identified, with coverage-based validation. Under the hood, SPADES utilizes the GOTTCHA2 taxonomy profiling bioinformatics algorithms and databases, extracts pathogen hits, computes detailed signature/genome coverage summaries, and automatically generates an interactive report/plot for rapid interpretation and triage.

## Requirements
- Conda or mamba (recommended)
- Dependencies in `environment.yml` (python 3.12, minimap2, krona, samtools, etc.)
- A GOTTCHA2 fast-profile database (base path with `.syldb`, `.zip`, `.stats`, `.tax.tsv`)
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

`run_SPADES.sh` uses `gottcha2 fast-profile`. The database path passed with
`-d/--db-path` should be the base path without extensions, and these files must
exist next to it:

- `<db>.syldb`
- `<db>.zip`
- `<db>.stats`
- `<db>.tax.tsv`

## Usage

```bash
./run_SPADES.sh \
  -i reads.fastq.gz \
  -o outdir \
  -p sample \
  -d /path/to/gottcha_db.species.fna \
  -t 4 \
  --js-external \
  [--spades-data data/] \
  [--ont] \
  [--ont-error-rate 0.03] \
  [--clean]
```

For detailed instructions on using GOTTCHA2, please see its [wiki page](https://github.com/poeli/GOTTCHA2/wiki).

Options:
- `-i, --input` Input reads file (fastq/fq[.gz], etc.)
- `-1, --read1` Paired-end read 1 (fastq/fq[.gz], etc.)
- `-2, --read2` Paired-end read 2 (fastq/fq[.gz], etc.)
- `-o, --outdir` Output directory
- `-p, --prefix` Output prefix (base filename)
- `-d, --db-path` GOTTCHA2 fast-profile database base path
- `-t, --cpu` Threads/CPUs (positive integer)
- `--spades-data` Directory containing `taxonomy_db/` and `pathogen.tsv` (default: `data/`)
- `--ont` Treat input as long reads; pre-processes with `fastplong`, splits to 150 bp, and passes relaxed error rate flags to GOTTCHA2
- `--ont-error-rate` Error rate for ONT reads passed to GOTTCHA2 (`-er`), default: `0.03`
- `--js-external` Use CDN-hosted JavaScript/CSS assets in generated HTML reports (required for local use)
- `--clean` Remove intermediate files after the run (see Outputs section for details)
- `--min-depth` Minimum depth for variant calling (default: 10)
- `--version` Show script version and exit
- `-h, --help` Show usage

## Outputs

The pipeline writes results to the output directory. Final outputs are kept in the main output directory, while intermediate files are organized in an `intermediate/` subdirectory for debugging purposes.

### Final Outputs (in `<outdir>/`)

- `<prefix>.info` Pipeline metadata (versions, database path, parameters)
- `<prefix>.krona.html` Interactive Krona taxonomic plot
- `<prefix>.qc.fastq.gz.html` FastP/FastPlong QC report
- `<prefix>.pathogen.tsv` Pathogen-annotated GOTTCHA2 hits (qualified taxa only)
- `<prefix>.pathogen.summary.txt` Pathogen detection summary text log
- `<prefix>.pathogen.summary.tsv` Genus/species level pathogen data for quick inspection
- `<prefix>.pathogen.full.tsv` Pathogen annotations for all detected taxa
- `<prefix>.pathogen.full.html` Interactive result viewer with filtering and visualization
- `<prefix>.coverage.html` Interactive coverage browser with variant information

### Intermediate Files (in `<outdir>/intermediate/`)

- `<prefix>.tsv` GOTTCHA2 profiling result table (filtered)
- `<prefix>.full.tsv` GOTTCHA2 full profiling output
- `<prefix>.gottcha_<dbLevel>.bam` Alignment file (with index)
- `<prefix>.gottcha_<dbLevel>.sam` SAM alignment (if produced)
- `<prefix>.gottcha_<dbLevel>.vcf.gz` Variant calls (with index)
- `<prefix>.gottcha_<dbLevel>.coverage.tsv` Coverage summary table
- `<prefix>.qc.fastq.gz` (or `R1/R2`) Quality-controlled reads
- `<prefix>.qc.fastq.gz.json` FastP/FastPlong JSON report
- `<prefix>.split_reads.fasta.gz` Split reads (ONT mode)
- `<prefix>.mpa.tsv` MetaPhlAn-compatible output
- `<prefix>.lineage.tsv` Lineage information
- `<prefix>.sylph_*` Sylph extraction files
- Reference FASTA files used for alignment

### Cleanup Behavior

- **Default** (no `--clean`): All files are preserved. Final outputs remain in `<outdir>/` for easy access, while intermediate files are moved to `<outdir>/intermediate/` for debugging and downstream analysis.
- **With `--clean`**: The `intermediate/` directory and all its contents are deleted after the pipeline completes, leaving only the final output files in `<outdir>/`.

If no qualifying GOTTCHA2 hits are found, the script creates placeholder outputs so downstream steps still receive expected files and then exits cleanly.

### GOTTCHA2 Full Report

See [tutorial.md](GOTTCHA2-main/docs/tutorial.md) for details about GOTTCHA2 results.

## Example

```bash
cd test
./test_run.sh
```

This runs a small ONT example using files under `test/` and the bundled data in
`data/`. Output can be found in the generated directories under `test/`.

## Notes

- The script exports its own `scripts/` folder onto `PATH` and defaults `--spades-data` to `data/` next to the script.
- Use `--js-external` when opening generated HTML reports directly from disk or from an environment that does not provide local `/publicdata` JavaScript/CSS assets.
- `--ont` uses `fastplong`, splits reads to 150 bp, and passes the error rate to GOTTCHA2 (`--ont-error-rate`, default `0.03`). Suggesting values: `0.05` for legacy R9.4.1 (HAC/SUP), `0.01` for R10.4.1 (Simplex/SUP), and `0.001` for R10.4.1 (Duplex/SUP).
- `--clean` removes the `intermediate/` directory after a successful run, keeping only final output files in the main output directory. Without `--clean`, all intermediate files are preserved in `intermediate/` for debugging.
- `--min-depth` sets the minimum depth threshold for variant calling in the coverage browser (default: 10).
- Ensure the database path points to the base file (e.g. `gottcha_db.species.fna`) and that the corresponding `.syldb`, `.zip`, `.stats`, and `.tax.tsv` exist.

## Notice of Copyright Assertion (O4958)

This program is Open-Source under the BSD-3 License.
 
Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 
Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 
Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 
Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.