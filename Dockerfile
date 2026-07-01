FROM mambaorg/micromamba AS builder

ARG MAMBA_DOCKERFILE_ACTIVATE=1
ARG MAMBA_USER=mambauser
ENV MAMBA_ROOT_PREFIX=/opt/conda
WORKDIR /tmp/build

COPY --chown=${MAMBA_USER}:${MAMBA_USER} environment.yml ./environment.yml
COPY --chown=${MAMBA_USER}:${MAMBA_USER} GOTTCHA2-main/ ./GOTTCHA2-main/
RUN micromamba env create -y -n spades -f environment.yml \
    && micromamba clean -a -y

FROM debian:bookworm-slim

ENV CONDA_PREFIX=/opt/conda/envs/spades
ENV LANG=C.UTF-8
ENV LC_ALL=C.UTF-8
ENV PATH="/opt/conda/envs/spades/bin:/app/scripts:${PATH}"

COPY --from=builder /opt/conda /opt/conda

WORKDIR /app
COPY data/ ./data/
COPY scripts/ ./scripts/
COPY run_SPADES.sh post_install.sh ./
COPY README.md Docker.md ./

RUN ./post_install.sh

CMD ["bash"]
