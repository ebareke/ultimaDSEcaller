# Container images

Both images bundle **ultimaDSEcaller** plus the read aligners (minimap2,
STAR, HISAT2) and samtools, so a single container runs the entire
**FASTQ → BAM → differential splicing events** workflow.

## Docker

```bash
# Build (from the repository root)
docker build -t ultimadsecaller:1.0.0 .

# Run the caller on BAMs you already have
docker run --rm -v "$PWD":/data ultimadsecaller:1.0.0 \
    run -s /data/samples.tsv -a /data/annotation.gtf -o /data/results \
        --contrast group:treatment-control

# Align FASTQ → BAM inside the same image (long-read example)
docker run --rm -v "$PWD":/data --entrypoint minimap2 ultimadsecaller:1.0.0 \
    -ax splice -uf /data/ref.fa /data/sample.fastq.gz \
  | docker run --rm -i -v "$PWD":/data --entrypoint samtools ultimadsecaller:1.0.0 \
        sort -o /data/sample.bam -
```

Published image (after a tagged release): `ghcr.io/ebareke/ultimadsecaller:1.0.0`.

## Apptainer / Singularity (HPC)

```bash
# Option A — pull the pre-built SIF published by the release workflow (ORAS)
apptainer pull ultimadsecaller.sif oras://ghcr.io/ebareke/ultimadsecaller-sif:1.0.0

# Option B — build the SIF from the published Docker image
apptainer build ultimadsecaller.sif containers/ultimadsecaller.def

# Option C — from a locally-built Docker image, no registry needed
docker build -t ultimadsecaller:1.0.0 .
apptainer build ultimadsecaller.sif docker-daemon://ultimadsecaller:1.0.0

# Run
apptainer run ultimadsecaller.sif run \
    -s samples.tsv -a annotation.gtf -o results \
    --contrast group:treatment-control
```

Each tagged release builds and pushes both artifacts automatically:

| Artifact | Location |
|---|---|
| Docker image | `ghcr.io/ebareke/ultimadsecaller:<version>` and `:latest` |
| Apptainer SIF | `oras://ghcr.io/ebareke/ultimadsecaller-sif:<version>` + GitHub Release asset |

## Image contents

| Tool            | Version | Purpose                                  |
|-----------------|---------|------------------------------------------|
| ultimaDSEcaller | 1.0.0   | Differential splicing event calling      |
| samtools        | 1.19    | BAM sort/index, CRAM handling            |
| minimap2        | 2.26    | Long-read (PacBio / ONT) spliced aligner |
| STAR            | 2.7.11b | Short-read (Illumina) spliced aligner    |
| HISAT2          | 2.2.1   | Alternative short-read spliced aligner   |

The container's caller binary is built on Ubuntu 22.04 (jammy), matching the
runtime's glibc. For a portable, statically-linked (musl) binary that runs on
any Linux x86-64 / aarch64 host, use the tarballs attached to each
[GitHub Release](https://github.com/ebareke/ultimaDSEcaller/releases).
