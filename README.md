# SPADES

The SPADES (Sequence-based Pathogen-Agnostic Diagnostics/Detection Solution) program is a full suite of  bioinformatics tools for detecting pathogens in complex metagenomic DNA sequencing datasets. The output of SPADES includes the detected organisms, highlighting the possible pathogens identified, with coverage-based validation. Under the hood, SPADES utilizes the GOTTCHA2 taxonomy profiling bioinformatics algorithms and databases, extracts pathogen hits, computes detailed signature/genome coverage summaries, and automatically generates an interactive report/plot for rapid interpretation and triage.

## Requirements
- Conda or mamba (recommended)
- Dependencies in `environment.yml` (python 3.12, minimap2, krona, samtools, etc.)
- A GOTTCHA2 fast-profile database (base path with `.syldb`, `.zip`, `.stats`, `.tax.tsv`)
- Taxonomy and pathogen data (provided in `data/`)
- Minimal hardware: 4 CPU cores, 24 GB RAM, 100 GB disk space

## Setup

### Local Conda/Mamba

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
- `--js-external` Use CDN-hosted JavaScript/CSS assets in generated HTML reports (required for local used)
- `--clean` Remove large intermediates (SAM, QC FASTQ, split FASTA) after the run
- `--version` Show script version and exit
- `-h, --help` Show usage

## Outputs

The pipeline writes results to the output directory, including:
- `<prefix>.tsv` GOTTCHA2 profiling result table (filtered)
- `<prefix>.full.tsv` GOTTCHA2 full profiling output file
- `<prefix>.krona.html` Krona plot
- `<prefix>.pathogen.tsv` Pathogen-annotated hits
- `<prefix>.pathogen.summary.txt` Pathogen summary log
- `<prefix>.pathogen.summary.tsv` Genus/species subset for quick inspection
- `<prefix>.pathogen.full.tsv` Pathogen annotations for all taxa
- `<prefix>.pathogen.full.html` Interactive result viewer
- `<prefix>.coverage.html` Coverage browser
- `<prefix>.gottcha_<dbLevel>.sam/.bam` Alignment files (supplementary)
- `<prefix>.gottcha_<dbLevel>.coverage.tsv` Coverage summary (supplementary)
- `<prefix>.info` Pipeline metadata (GOTTCHA2 version)

If no qualifying GOTTCHA2 hits are found, the script creates placeholder outputs so downstream steps still receive expected files and then exits cleanly.

### GOTTCHA2 Full Report Fields

The GOTTCHA2 full report contains detailed information about the taxonomic profiling results. This document describes all 26 columns in the full report.

| Field Name             | Description                                                                                                                 |
| ---------------------- | --------------------------------------------------------------------------------------------------------------------------- |
| LEVEL                  | Taxonomic level (major levels from superkingdom to strain)                                                                  |
| NAME                   | Name of the taxon                                                                                                           |
| TAXID                  | NCBI Taxonomy ID                                                                                                            |
| READ_COUNT             | Number of reads mapped to this taxon                                                                                        |
| TOTAL_BP_MAPPED        | Total number of base pairs mapped to this taxon’s signatures                                                                |
| SNI_SCORE              | Signature Nucleotide Identity (SNI) score estimated from mapped signatures (default 0.95 for species, 0.99 for strains)     |
| COVERED_SIG_LEN        | Total length of mapped signature regions (formerly “linear coverage”)                                                       |
| BEST_SIG_COV           | Highest signature coverage observed (SIG_COV) at strain level                                                               |
| DEPTH                  | Depth of coverage (TOTAL_BP_MAPPED / TOTAL_SIG_LEN)                                                                         |
| REL_ABUNDANCE          | Relative abundance calculated using the specified abundance field (default: DEPTH)                                          |
| PARENT_NAME            | Name of the parent taxon                                                                                                    |
| PARENT_TAXID           | NCBI Taxonomy ID of the parent taxon                                                                                        |
| TOTAL_READ_LEN         | Total length of the mapped reads                                                                                            |
| READ_IDT               | Mapping identity (commonly (TOTAL_BP_MAPPED − TOTAL_BP_MISMATCH) / TOTAL_READ_LEN; may also account for indels if reported) |
| TOTAL_BP_MISMATCH      | Total number of base-pair mismatches in the alignments                                                                      |
| TOTAL_BP_INDEL         | Total number of base pairs involved in insertions/deletions (indels) in the alignments                                      |
| SNI_NAIVE              | “Naive” SNI estimate (unadjusted identity estimate, as reported by the tool)                                                |
| SNI_CI95_LH            | 95% confidence interval bounds for SNI (low/high, as reported)                                                              |
| SIG_COV                | Signature coverage (COVERED_SIG_LEN / TOTAL_SIG_LEN)                                                                        |
| MAPPED_SIG_LEN         | Total length of signature fragments that had at least one read mapped                                                       |
| TOTAL_SIG_LEN          | Total length of all signature sequences for this taxon                                                                      |
| COVERED_SIG_DEPTH      | Depth over covered regions only (TOTAL_BP_MAPPED / COVERED_SIG_LEN)                                                         |
| COVERED_MAPPED_SIG_COV | Ratio of covered signature length to mapped signature length (COVERED_SIG_LEN / MAPPED_SIG_LEN)                             |
| ZSCORE                 | Estimated Z-score for depth distribution (coverage uniformity indicator; lower is better)                                   |
| GENOMIC_CONTENT_EST    | Estimated genomic content scaled by genome size (commonly TOTAL_BP_MAPPED / TOTAL_SIG_LEN × GENOME_SIZE)                    |
| ABUNDANCE              | Raw abundance value (from the specified column used for REL_ABUNDANCE calculation)                                          |
| REL_ABUNDANCE_DEPTH    | Relative abundance calculated using depth of coverage                                                                       |
| REL_ABUNDANCE_GC       | Relative abundance calculated using genomic content estimate                                                                |
| SIG_LEVEL              | Taxonomic level of the signatures used for mapping                                                                          |
| GENOME_COUNT           | Number of reference genomes rolled up into this taxon                                                                       |
| GENOME_SIZE            | Combined size of the original genomes for this taxon                                                                        |
| NOTE                   | Additional information, including reasons for filtering if applicable                                                       |
| HUMAN_PATHOGEN         | Yes / No                                                                                                                    |
| PATHOGENIC_INFO        | Display pathogenic information of the taxa                                                                                  |

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
- `--clean` removes large intermediates (`*.sam`, QC FASTQ, split FASTA, and fast-profile Sylph extraction files) after a successful run.
- Ensure the database path points to the base file (e.g. `gottcha_db.species.fna`) and that the corresponding `.syldb`, `.zip`, `.stats`, and `.tax.tsv` exist.

## Notice of Copyright Assertion (O4958)

This program is Open-Source under the BSD-3 License.
 
Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 
Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 
Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 
Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.