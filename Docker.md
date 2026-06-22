# Docker Usage

This repo includes a two-stage Dockerfile that builds the `spades` conda
environment, installs the vendored GOTTCHA2 package, and copies the runtime
files into a slim Debian image.

## Local Build

From the repo root:

```bash
docker build -t spades-g2:local .
```

Optional multi-platform publish examples:

```bash
docker buildx build --platform linux/arm64 -t spades-g2:[DATE]-arm64 .
docker buildx build --platform linux/amd64 -t spades-g2:[DATE]-amd64 .
docker buildx build --platform linux/amd64,linux/arm64 -t spades-g2:[DATE] -t spades-g2:latest .
```

## Check The Image

The image contains `data/`, `scripts/`, `run_SPADES.sh`, and `post_install.sh`.
Run the lightweight install checker before using a new image:

```bash
docker run --rm -it spades-g2:local ./scripts/check_install.sh
```

The checker does not run the sequencing workflow. It verifies required command
line tools, Python imports, SPADES data files, and Krona taxonomy wiring.

## Run

Mount your reads and database into the container so the pipeline can read inputs
and write outputs back to the host filesystem.

```bash
docker run --rm -it \
  -v "$PWD":/work \
  -v /path/to/gottcha_db:/db:ro \
  -w /work \
  spades-g2:local \
  /app/run_SPADES.sh \
    -i /work/reads.fastq.gz \
    -o /work/outdir \
    -p sample \
    -d /db/gottcha_db.species.fna \
    -t 8 \
    --spades-data /app/data \
    --js-external
```

Notes:

- `-d` expects the GOTTCHA2 fast-profile database base path without extensions.
- The database directory must contain:
  - `<db>.syldb`
  - `<db>.zip`
  - `<db>.stats`
  - `<db>.tax.tsv`
- The image already contains `data/taxonomy_db` and `data/pathogen.tsv`.
- Add `--ont` to `run_SPADES.sh` for long-read single-end inputs.
- Use `--js-external` unless the HTML reports will be served from an environment
  that provides the local `/publicdata` JavaScript/CSS assets.

## Run The Bundled Test

The runtime image does not bake in the full `test/` directory. Mount the repo if
you want to run the bundled test script:

```bash
docker run --rm -it \
  -v "$PWD":/work \
  -w /work/test \
  spades-g2:local \
  ./test_run.sh
```

