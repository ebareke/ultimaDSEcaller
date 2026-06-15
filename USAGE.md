# `ultimaDSEcaller` — Tool Usage Guide

Complete reference for the v1.0.0 CLI. For an overview of capabilities,
see [`README.md`](README.md).

---

## Contents

1. [Starting from FASTQ](#starting-from-fastq)
2. [Subcommands](#subcommands)
3. [`run` — the full pipeline](#run--the-full-pipeline)
4. [`junctions` — junction extraction only](#junctions--junction-extraction-only)
5. [`dump-annotation` — splice-graph inspection](#dump-annotation--splice-graph-inspection)
6. [`report` — re-render HTML](#report--re-render-html)
7. [`pdf` — re-render PDF](#pdf--re-render-pdf)
8. [Sample sheet schema](#sample-sheet-schema)
9. [Design formula syntax](#design-formula-syntax)
10. [Consensus statistical engine](#consensus-statistical-engine)
11. [Environment variables](#environment-variables)
12. [Output directory layout](#output-directory-layout)
13. [Recipes](#recipes)
14. [Error codes](#error-codes)

---

## Starting from FASTQ

The `ultimaDSEcaller` binary itself consumes **BAM/CRAM** (the standard
interchange format). To go all the way from raw **FASTQ**, use either the
Nextflow pipeline or the container, both of which bundle the aligners:

```bash
# Nextflow (STAR for short-read, minimap2 for long-read; auto-routed per sample)
nextflow run nextflow/main.nf \
    --samplesheet samples.csv --gtf annotation.gtf --fasta genome.fa \
    --contrast 'group:treatment-control' -profile docker

# Or align by hand inside the container, then call:
docker run --rm -v "$PWD":/data --entrypoint minimap2 ultimadsecaller:1.0.0 \
    -ax splice -uf /data/genome.fa /data/sample.fastq.gz \
  | samtools sort -o sample.bam -
```

See [`nextflow/README.md`](nextflow/README.md) and
[`containers/README.md`](containers/README.md). The rest of this guide covers
the caller once you have BAMs.

---

## Subcommands

| Subcommand        | Purpose                                                            |
|-------------------|--------------------------------------------------------------------|
| `run`             | End-to-end pipeline: detect → quantify → test → report             |
| `junctions`       | Extract splice junctions from BAM/CRAM(s) to TSV (no event calls)  |
| `dump-annotation` | Parse GTF/GFF3 and dump per-gene splice graphs as JSON             |
| `report`          | Re-render the HTML report from an existing `results.json`          |
| `pdf`             | Render a PDF report from a `results.json` (requires `--features pdf`) |

Global flags: `--verbose`, `--threads N` (0 = all cores), `--help`, `--version`.

---

## `run` — the full pipeline

```text
ultimaDSEcaller run [OPTIONS]
```

### Required inputs

| Flag                  | Description                                                  |
|-----------------------|--------------------------------------------------------------|
| `-s, --sample-sheet`  | TSV/CSV with `sample`, `bam`, `group` columns                |
| `-a, --annotation`    | GTF or GFF3 (optionally bgzipped)                            |
| `-o, --out`           | Output directory (created if it doesn't exist)               |

### Optional inputs / config

| Flag                      | Description                                                  |
|---------------------------|--------------------------------------------------------------|
| `-c, --config <PATH>`     | YAML / TOML / JSON config file (CLI flags override it)       |
| `--reference <FASTA>`     | Required for CRAM; enables motif + protein analyses for BAM  |
| `--fusion-bedpe <PATH>`   | Fusion-caller BEDPE for FusionAssociated event detection     |

### Sequencing tech

| Flag             | Values                          | Effect                            |
|------------------|---------------------------------|-----------------------------------|
| `--tech`         | `short`, `pac-bio`, `ont`       | Sets defaults; long-read enables isoform reconstruction |
| `--min-mapq N`   | u8 (default 10 short, 1 long)   | Minimum mapping quality           |
| `--min-overhang N` | u32 (default 8 / 4)           | Anchor length per side of junction |
| `--junction-tolerance N` | u32 (default 0 / 5)      | Wobble bp for long-read collapsing  |
| `--multimap`     | `discard`, `primary`, `fractional` | NH>1 read handling             |

### Contrasts & formula

| Flag                  | Description                                                                 |
|-----------------------|-----------------------------------------------------------------------------|
| `--contrast <SPEC>`   | Single contrast: `variable:level1-level2` (e.g. `group:treatment-control`)  |
| `--contrasts <SPECS>` | Comma-separated multiple contrasts; each gets its own output subdirectory   |
| `--formula <FORMULA>` | Wilkinson-style design (e.g. `~ batch + group + batch:group`)               |
| `--random-effect <COL>` | Sample-sheet column used as the random intercept under `--test glmm`      |

### Statistical test

| Flag                 | Values                                       |
|----------------------|----------------------------------------------|
| `--test`             | `bb-lrt` (default), `glm`, `glmm`            |
| `--consensus`        | `stouffer` (default), `brown`, `weighted-fisher` — combines all available tests into one q-value |

### Filtering thresholds (high-confidence set)

| Flag                 | Default | Description                                  |
|----------------------|---------|----------------------------------------------|
| `--max-fdr`          | 0.05    | Maximum BH-adjusted p-value                  |
| `--min-delta-psi`    | 0.10    | Minimum \|ΔPSI\|                             |
| `--min-coverage`     | 10      | Minimum sample-mean coverage                 |

### Performance / IO

| Flag                | Effect                                                              |
|---------------------|---------------------------------------------------------------------|
| `--resume`          | Skip stages whose checkpoint files already exist                    |
| `--no-cache`        | Always re-parse the annotation (skip `.ultidse` cache)              |
| `--cache-path <P>`  | Override the default `<annotation>.ultidse` cache path              |
| `--threads N`       | rayon thread pool size (0 = all cores)                              |

### Output add-ons

| Flag        | Effect                                                                   |
|-------------|--------------------------------------------------------------------------|
| `--pdf`     | Also write `report.pdf` (requires the binary to be built with `--features pdf`) |

---

## `junctions` — junction extraction only

```text
ultimaDSEcaller junctions BAMS... -o junctions.tsv \
  [--reference FASTA] [--min-mapq N] [--min-overhang N]
```

Writes a TSV with columns `sample · chrom · donor_end · acceptor_start · count`.

---

## `dump-annotation` — splice-graph inspection

```text
ultimaDSEcaller dump-annotation -a annotation.gtf.gz [-o splice_graphs.json] [--gene ENSG00000123456]
```

Emits each gene's exons + transcripts + edge counts as JSON. Useful for
debugging parser behavior on idiosyncratic annotation files.

---

## `report` — re-render HTML

```text
ultimaDSEcaller report -i results/results.json -o results/report.html
```

Regenerates the HTML report from a previously-produced JSON summary.
Useful for re-rendering after a template change without re-running the
pipeline.

---

## `pdf` — re-render PDF

```text
ultimaDSEcaller pdf -i results/results.json -o results/report.pdf
```

PDF render from existing JSON. Requires `--features pdf` at build time.
The PDF includes:

* Cover page with run metadata and per-kind counts
* Top-20 events page per AS category
* Full high-confidence event set
* Consensus statistical results
* Protein-consequence breakdown

Plots are **not** embedded — for interactive plots, view `report.html`
in a browser and use the browser's "Print to PDF" feature if a
plot-rich PDF is needed.

---

## Sample sheet schema

TSV (or CSV — extension-detected). Required columns: `sample`, `bam`,
`group`. Any additional columns become covariates available to
`--formula` and `--random-effect`.

```tsv
sample	bam	group	batch	subject	rin
s01	/data/s01.bam	control	1	A	8.7
s02	/data/s02.bam	control	2	B	8.4
s03	/data/s03.bam	treatment	1	A	8.9
s04	/data/s04.bam	treatment	2	B	8.1
```

Numeric covariates (`rin`) are auto-detected and used as continuous
predictors; non-numeric covariates (`batch`, `subject`) are
treatment-coded with the first level as reference.

---

## Design formula syntax

Wilkinson notation, R/DESeq2-compatible subset:

| Formula                          | Means                                              |
|----------------------------------|----------------------------------------------------|
| `~ group`                        | Intercept + main effect of `group`                 |
| `~ group + batch`                | Additive main effects                              |
| `~ group + batch + group:batch`  | Main effects + their interaction                   |
| `~ group * batch`                | Sugar for the line above                           |
| `~ 0 + group`                    | No intercept (group means parameterization)        |
| `~ batch + group + rin`          | Mixed categorical + continuous covariates          |

Not yet supported: nested groupings (`a/b`), function transforms
(`I(x^2)`, `log(x)`), `factor()`, `^N` exponentiation.

---

## Consensus statistical engine

When `--consensus` is set, the pipeline runs **every** applicable test
per event (BB-LRT, GLM, GLMM, Fisher 2×2) and combines their p-values:

| Method               | Best for                                                       |
|----------------------|----------------------------------------------------------------|
| `stouffer` (default) | Weighted z-score mean; closed-form, robust to mild correlation |
| `brown`              | Fisher χ² with correlation correction; conservative            |
| `weighted-fisher`    | Classic Fisher with per-method weights                         |

Tunables (env vars):
- `ULTIMADSE_BROWN_RHO=0.3` — assumed inter-test correlation for Brown's method.

Consensus output:
- `combined_p` — pre-FDR consensus p-value
- `consensus_q` — BH-adjusted across all events
- `confidence` — composite (consensus, motif, coverage, reproducibility)
- `methods_combined` — how many tests contributed (1–4)

---

## Environment variables

| Variable                       | Effect                                                |
|--------------------------------|-------------------------------------------------------|
| `ULTIMADSE_NO_PROGRESS=1`   | Suppress progress bars (CI / non-TTY default)         |
| `ULTIMADSE_IR_READ_LEN=N`   | Override effective read length used in IR PSI         |
| `ULTIMADSE_BROWN_RHO=R`     | Inter-test correlation assumed by Brown's method      |
| `RUST_LOG=ultimadsecaller=debug`  | Fine-grained log filter (overrides `-v`)              |

---

## Output directory layout

```
out/
├── events.raw.tsv / .csv                  # full per-event table
├── events.high_confidence.tsv / .csv      # FDR + ΔPSI + confidence filtered
├── events.raw.parquet                     # if --features parquet
├── by_event_type/
│   ├── SE.raw.tsv      SE.high_confidence.tsv
│   ├── MXE.raw.tsv     MXE.high_confidence.tsv
│   ├── A5SS.raw.tsv    A5SS.high_confidence.tsv
│   ├── A3SS.raw.tsv    A3SS.high_confidence.tsv
│   ├── IR.raw.tsv      IR.high_confidence.tsv
│   ├── AFE.raw.tsv     AFE.high_confidence.tsv
│   └── ALE.raw.tsv     ALE.high_confidence.tsv
├── advanced_events.tsv                    # cryptic / MSE / MIR / recursive / ...
├── junctions.bin                          # checkpoint for --resume
├── results.json                           # machine-readable summary
├── report.html                            # interactive HTML report
├── report.pdf                             # if --pdf and --features pdf
└── contrast_<num>_vs_<denom>/             # for each --contrasts entry
    └── (same layout as above)
```

---

## Recipes

### Standard 2-group short-read run

```bash
ultimaDSEcaller run \
  --sample-sheet samples.tsv \
  --annotation gencode.v45.gtf.gz \
  --out results/ \
  --contrast group:treatment-control
```

### With reference FASTA (motif + protein consequence)

```bash
ultimaDSEcaller run \
  -s samples.tsv -a gencode.v45.gtf.gz \
  --reference GRCh38.fa --out results/ \
  --contrast group:treatment-control
```

### Multi-batch design via formula

```bash
ultimaDSEcaller run \
  -s samples.tsv -a gencode.v45.gtf.gz --out results/ \
  --contrast group:treatment-control \
  --test glm --formula '~ batch + group'
```

### Paired-subject GLMM

```bash
ultimaDSEcaller run \
  -s samples.tsv -a gencode.v45.gtf.gz --out results/ \
  --contrast group:post-pre \
  --test glmm --random-effect subject --formula '~ group'
```

### Multi-contrast cohort

```bash
ultimaDSEcaller run \
  -s samples.tsv -a gencode.v45.gtf.gz --out results/ \
  --contrasts 'tissue:liver-brain,tissue:lung-brain,tissue:kidney-brain' \
  --consensus stouffer
```

### Long-read isoform reconstruction

```bash
ultimaDSEcaller run \
  -s samples.tsv -a gencode.v45.gtf.gz --out results/ \
  --reference GRCh38.fa --tech ont \
  --contrast group:tumor-normal \
  --junction-tolerance 5
```

### Resume after BAM extraction

```bash
ultimaDSEcaller run ... --resume     # skips junction extraction if junctions.bin exists
```

### PDF + HTML in one run

```bash
# Built with: cargo build --release --features 'pdf parquet'
ultimaDSEcaller run ... --pdf
```

### Pipeline-only junction extraction

```bash
ultimaDSEcaller junctions data/*.bam -o cohort_junctions.tsv --min-mapq 30
```

---

## Error codes

All errors expose a stable machine-readable code. Tooling should pattern-match on the code, not the human message.

| Code   | Meaning                                |
|--------|----------------------------------------|
| `E0001` | I/O                                    |
| `E0010` | Annotation parse                       |
| `E0020` | BAM/CRAM parse                         |
| `E0030` | Reference FASTA                        |
| `E0040` | Configuration                          |
| `E0050` | Sample-sheet / design                  |
| `E0051` | Design formula parse                   |
| `E0060` | Statistical computation                |
| `E0070` | Cache                                  |
| `E0090` | Unsupported feature                    |
| `E0099` | Other                                  |
