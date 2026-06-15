# syntax=docker/dockerfile:1
#
# ultimaDSEcaller — production container.
#
# Stage 1 builds a fully static (musl) binary with all features.
# Stage 2 is a micromamba runtime that ALSO bundles the read aligners
# (minimap2, STAR, HISAT2) and samtools, so the image can run the complete
# FASTQ -> BAM -> differential-splicing-events pipeline on its own.
#
#   docker build -t ultimadsecaller:1.0.0 .
#   docker run --rm -v "$PWD":/data ultimadsecaller:1.0.0 \
#       run -s /data/samples.tsv -a /data/annotation.gtf -o /data/results \
#           --contrast group:treatment-control
#
# ---------------------------------------------------------------------------
FROM rust:1.85-bookworm AS builder

RUN apt-get update && apt-get install -y --no-install-recommends \
        musl-tools pkg-config \
    && rm -rf /var/lib/apt/lists/*
RUN rustup target add x86_64-unknown-linux-musl

WORKDIR /build
COPY Cargo.toml Cargo.lock ./
COPY src ./src
COPY benches ./benches
COPY tests ./tests

ENV CARGO_TARGET_X86_64_UNKNOWN_LINUX_MUSL_LINKER=musl-gcc
RUN cargo build --release --locked --all-features \
        --target x86_64-unknown-linux-musl --bin ultimaDSEcaller \
    && strip target/x86_64-unknown-linux-musl/release/ultimaDSEcaller

# ---------------------------------------------------------------------------
FROM mambaorg/micromamba:1.5-jammy AS runtime

LABEL org.opencontainers.image.title="ultimaDSEcaller" \
      org.opencontainers.image.description="Ultimate Differential Splicing Event caller for short- and long-read RNA-seq" \
      org.opencontainers.image.source="https://github.com/ebareke/ultimaDSEcaller" \
      org.opencontainers.image.url="https://ebareke.github.io/ultimaDSEcaller/" \
      org.opencontainers.image.licenses="MIT" \
      org.opencontainers.image.version="1.0.0"

# Bioinformatics toolchain for the FASTQ -> BAM front-end.
RUN micromamba install -y -n base -c conda-forge -c bioconda \
        samtools=1.19 \
        minimap2=2.26 \
        star=2.7.11b \
        hisat2=2.2.1 \
        htslib=1.19 \
    && micromamba clean --all --yes

# Static caller binary - works regardless of the conda env being active.
COPY --from=builder /build/target/x86_64-unknown-linux-musl/release/ultimaDSEcaller /usr/local/bin/ultimaDSEcaller

ENV PATH=/opt/conda/bin:$PATH
WORKDIR /data

# Quick self-check at build time.
RUN /usr/local/bin/ultimaDSEcaller --version

ENTRYPOINT ["/usr/local/bin/ultimaDSEcaller"]
CMD ["--help"]
