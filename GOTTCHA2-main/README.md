[![logo](https://gottcha2.readthedocs.io/en/latest/_images/gottcha_icon.png)](https://gottcha2.readthedocs.io/en/latest/_images/gottcha_icon.png)

# Genomic Origin Through Taxonomic CHAllenge (GOTTCHA)

[![gottcha2](https://github.com/poeli/GOTTCHA2/actions/workflows/gottcha.yml/badge.svg?branch=master)](https://github.com/poeli/GOTTCHA2/actions/workflows/gottcha.yml)
[![bioconda](https://anaconda.org/bioconda/gottcha2/badges/version.svg)](https://anaconda.org/bioconda/gottcha2)


GOTTCHA is an application of a novel, gene-independent and signature-based metagenomic taxonomic profiling
method with significantly smaller false discovery rates (FDR) that is laptop deployable. Our algorithm was
tested and validated on twenty synthetic and mock datasets ranging in community composition and complexity,
was applied successfully to data generated from spiked environmental and clinical samples, and robustly
demonstrates superior performance compared with other available tools.

GOTTCHAv2 is currently under development in BETA stage. Pre-built databases for v1 are incompatible with v2.

-------------------------------------------------------------------
## DEPENDENCIES

GOTTCHA2 profiler is written in Python 3 and uses minimap2 to map reads to signature sequences. To run GOTTCHA2, install the following dependencies. A ready-to-use Conda environment file is provided at `environment.yml`.

- Python 3.9+
- minimap2 2.17+
- samtools
- numpy
- pandas
- requests
- biom-format
- pysam
- tqdm
- sylph

-------------------------------------------------------------------
## QUICK START

1. Install the package:

        via conda `conda install -c bioconda gottcha2`

        OR

        Download or git clone GOTTCHA2 from this repository and run:

        `python -m pip install .`

        (For development installs: `python -m pip install -e .`)

2. Download the latest version of the GOTTCHA2 database. (This step may take some time)

        https://ref-db.edgebioinformatics.org/gottcha2/RefSeq-r220/

3. Run GOTTCHA2:
        
        Regular mode:

        `gottcha2 profile -d /path/to/db/ -t 8 -i <FASTQ>`

        OR

        Fast-profiling mode:

        `gottcha2 fast-profile -d /path/to/db/ -t 8 -i <FASTQ>`

-------------------------------------------------------------------
## RESULT

GOTTCHA2 can output the profiling results in either CSV, TSV or BIOM format.
- summary (.tsv or .csv) - A summary of profiling results (10 columns) in taxonomic ranks breakdown
- full (.tsv or .csv) - A full profiling results including unfiltered profiling results and additional columns
- lineage (.lineage.tsv or .lineage.tsv) - output lineage and abundance of the profiled taxon per line
- extract (.extract[TAXID].fastq) - Extracted reads for a specific taxon.

-------------------------------------------------------------------
## DOCUMENTATION

Please refer to https://github.com/poeli/GOTTCHA2/wiki for more details.
