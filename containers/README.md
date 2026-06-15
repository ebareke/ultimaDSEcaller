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
# Option A — from the published image
apptainer build ultimadsecaller.sif containers/ultimadsecaller.def

# Option B — from a locally-built Docker image, no registry needed
docker build -t ultimadsecaller:1.0.0 .
apptainer build ultimadsecaller.sif docker-daemon://ultimadsecaller:1.0.0

# Run
apptainer run ultimadsecaller.sif run \
    -s samples.tsv -a annotation.gtf -o results \
    --contrast group:treatment-control
```

## Image contents

| Tool            | Version | Purpose                                  |
|-----------------|---------|------------------------------------------|
| ultimaDSEcaller | 1.0.0   | Differential splicing event calling      |
| samtools        | 1.19    | BAM sort/index, CRAM handling            |
| minimap2        | 2.26    | Long-read (PacBio / ONT) spliced aligner |
| STAR            | 2.7.11b | Short-read (Illumina) spliced aligner    |
| HISAT2          | 2.2.1   | Alternative short-read spliced aligner   |

The caller itself is a static musl binary, so it also works when copied out
of the image onto any Linux x86-64 host.
