# syntax=docker/dockerfile:1
#
# ultimaDSEcaller — production container.
#
# Stage 1 builds the binary on Ubuntu 22.04 (jammy). Stage 2 — the micromamba
# runtime — is ALSO jammy-based, so the builder and runtime share an identical
# glibc. That avoids the dynamic-loader mismatch you hit when a musl binary
# (linked to /lib/ld-musl-*.so) lands on a glibc image, or when a binary built
# against a newer glibc runs on an older one. The C-dependency crates
# (bzip2 / xz, via noodles' CRAM support) vendor and statically link their own
# C code, so the only dynamic dependency that remains is the base glibc, which
# is present in the runtime.
#
# The runtime also bundles the read aligners (minimap2, STAR, HISAT2) and
# samtools, so the image runs the complete FASTQ -> BAM -> events pipeline.
#
#   docker build -t ultimadsecaller:1.0.0 .
#   docker run --rm -v "$PWD":/data ultimadsecaller:1.0.0 \
#       run -s /data/samples.tsv -a /data/annotation.gtf -o /data/results \
#           --contrast group:treatment-control
#
# ---------------------------------------------------------------------------
FROM ubuntu:22.04 AS builder

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates curl build-essential pkg-config cmake \
    && rm -rf /var/lib/apt/lists/*

# Pin the Rust toolchain to the project MSRV.
ENV RUSTUP_HOME=/usr/local/rustup CARGO_HOME=/usr/local/cargo \
    PATH=/usr/local/cargo/bin:$PATH
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs \
        | sh -s -- -y --no-modify-path --profile minimal --default-toolchain 1.85.0

WORKDIR /build
COPY Cargo.toml Cargo.lock ./
COPY src ./src
COPY benches ./benches
COPY tests ./tests

# Default (glibc) target, all features, locked deps.
RUN cargo build --release --locked --all-features --bin ultimaDSEcaller \
    && strip target/release/ultimaDSEcaller \
    && ./target/release/ultimaDSEcaller --version

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

# Caller binary (glibc; runs against the jammy base libc present in this image).
COPY --from=builder --chmod=755 /build/target/release/ultimaDSEcaller /usr/local/bin/ultimaDSEcaller

ENV PATH=/opt/conda/bin:$PATH
WORKDIR /data

# Quick self-check at build time.
RUN /usr/local/bin/ultimaDSEcaller --version

ENTRYPOINT ["/usr/local/bin/ultimaDSEcaller"]
CMD ["--help"]
