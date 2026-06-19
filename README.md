# ultimaDSEcaller

**Ultimate Differential Splicing Event caller for short- and long-read RNA-seq.**

`ultimaDSEcaller` detects, quantifies, statistically tests, and reports
**differential alternative-splicing events** between groups of samples — from
Illumina short reads and PacBio / Oxford Nanopore long reads alike. It is a
single, fast, statically-linkable Rust binary, with a Nextflow pipeline and
container images that take you all the way from **FASTQ to publication-ready
results**.

[![CI](https://github.com/ebareke/ultimaDSEcaller/actions/workflows/ci.yml/badge.svg)](https://github.com/ebareke/ultimaDSEcaller/actions/workflows/ci.yml)
![Rust](https://img.shields.io/badge/rust-1.85%2B-orange)
![License](https://img.shields.io/badge/license-MIT-blue)
![Platforms](https://img.shields.io/badge/platform-linux%20%7C%20macOS-lightgrey)
[![Docs](https://img.shields.io/badge/docs-ebareke.github.io-1E6B4F)](https://ebareke.github.io/ultimaDSEcaller/)

Documentation: **<https://ebareke.github.io/ultimaDSEcaller/>**

---

## Features

- 🧬 **18 event classes** — the 7 canonical (SE, MXE, A5SS, A3SS, IR, AFE,
  ALE) plus 11 advanced: cryptic splice sites, multi-exon skipping,
  multi-intron retention, recursive & nested splicing, partial-exon
  inclusion, exonic-intronic hybrids, alternative promoter / polyadenylation,
  tandem UTRs, and fusion-associated junctions.
- 📏 **Short *and* long read** — junction-based PSI for Illumina; full-length
  isoform reconstruction and differential isoform usage for PacBio / ONT.
- 🧮 **Rigorous statistics** — beta-binomial likelihood-ratio test, logistic
  GLM / GLMM (random effects), Fisher's exact, and a **consensus engine**
  (Stouffer / Brown / weighted-Fisher) that fuses them into one calibrated
  q-value.
- 🧠 **Functional annotation** — GT-AG / GC-AG / AT-AC splice-motif scoring
  from a reference FASTA, plus protein-consequence and NMD prediction.
- 📊 **Self-contained HTML report** — volcano, MA, PCA, UMAP, heatmaps,
  per-position sashimi tracks, splice-graph views, isoform Sankey diagrams.
- 🗂️ **Many output formats** — TSV, CSV, JSON, Parquet, HTML, and PDF.
- ⚡ **Built for scale** — multithreaded (rayon), BAM-index aware,
  annotation caching, checkpoint/resume.
- 📦 **Reproducible** — Docker + Apptainer images and a Nextflow pipeline
  that run the whole FASTQ → events path.

## Install

### From source (Rust ≥ 1.85)

```bash
git clone https://github.com/ebareke/ultimaDSEcaller.git
cd ultimaDSEcaller
cargo build --release                 # binary at target/release/ultimaDSEcaller
cargo build --release --all-features  # adds Parquet + PDF output
```

### Container (no Rust needed; bundles aligners + samtools)

```bash
docker build -t ultimadsecaller:1.0.0 .
docker run --rm -v "$PWD":/data ultimadsecaller:1.0.0 --help
```

See [containers/README.md](containers/README.md) for Apptainer/Singularity.

## Quick start

### Already have BAMs

```bash
ultimaDSEcaller run \
    --sample-sheet samples.tsv \
    --annotation gencode.v45.gtf.gz \
    --reference GRCh38.fa \
    --out results/ \
    --contrast group:treatment-control \
    --tech short \
    --consensus stouffer
```

### Starting from FASTQ (Nextflow)

```bash
nextflow run nextflow/main.nf \
    --samplesheet samples.csv \
    --gtf annotation.gtf --fasta genome.fa \
    --contrast 'group:treatment-control' \
    -profile docker
```

### Try the bundled example (no aligner needed)

```bash
cargo build --release
bash example/run_example.sh
```

This generates a tiny synthetic dataset with a single exon-skipping event
(included in `control`, skipped in `treatment`) and calls it end-to-end —
proving the install in seconds. See [example/README.md](example/README.md).

## How it works

```
FASTQ ──▶ align ──▶ sorted BAM ──▶ ultimaDSEcaller ──▶ events + report
        STAR (short)            │   ├─ junction extraction
        minimap2 (long)         │   ├─ splice-graph event detection
                                │   ├─ PSI / ΔPSI quantification
                                │   ├─ BB-LRT / GLM / GLMM / consensus
                                │   ├─ motif + protein consequence
                                │   └─ TSV / CSV / JSON / Parquet / HTML / PDF
```

The Rust core consumes BAM/CRAM (the standard interchange format, like
rMATS / MAJIQ / SUPPA2); the Nextflow pipeline and container images provide
the alignment front-end so the end-to-end path is one command.

## Project layout

```
src/                Rust library + CLI (20 modules)
  annotation.rs       GTF/GFF3 parser + per-gene splice graphs
  junctions.rs        BAM/CRAM junction extraction (streaming + indexed)
  events.rs           Canonical event detection
  advanced.rs         Cryptic / MSE / MIR / recursive / nested / …
  quantify.rs         PSI / ΔPSI / coverage / confidence
  stats.rs glm.rs     BB-LRT, Fisher, BH-FDR, GLM, GLMM
  consensus.rs        Multi-test combination
  motif.rs protein.rs Splice motifs, protein consequence + NMD
  longread.rs         Isoform reconstruction + differential usage
  embedding.rs sashimi.rs output.rs report.rs pdf.rs   Reporting
tests/              integration + property tests
benches/            criterion benchmarks
example/            synthetic end-to-end example (BAM and FASTQ paths)
nextflow/           FASTQ → events pipeline (STAR / minimap2 + caller)
containers/         Apptainer definition + container docs
docs/               published documentation site
```

## Documentation

| File | Contents |
|---|---|
| [USAGE.md](USAGE.md) | Complete CLI reference, recipes, output schema |
| [example/README.md](example/README.md) | Runnable end-to-end example |
| [nextflow/README.md](nextflow/README.md) | FASTQ → events pipeline |
| [containers/README.md](containers/README.md) | Docker / Apptainer images |
| [CHANGELOG.md](CHANGELOG.md) | Release history |
| [ROADMAP.md](ROADMAP.md) | Planned work |
| [CONTRIBUTING.md](CONTRIBUTING.md) | How to contribute |
| [SECURITY.md](SECURITY.md) | Security model and reporting |
| [CITATION.md](CITATION.md) | How to cite ultimaDSEcaller |

## Authors

- **Eric B.** — <eb.bioinfo@pm.me>
- **Ethan B.** — <b.bioinfo@pm.me>
- **Conrad B.** — <b.bioinfo@pm.me>

## Citation

If you use ultimaDSEcaller in your research, please cite it — see
[CITATION.md](CITATION.md) (or the **Cite this repository** button on GitHub,
powered by [CITATION.cff](CITATION.cff)).

> Bareke, E., Bareke, E., & Bareke, C. (2026). *ultimaDSEcaller: Ultimate
> Differential Splicing Event caller for short- and long-read RNA-seq*
> (Version 1.0.0) [Computer software]. https://github.com/ebareke/ultimaDSEcaller

## License

[MIT](LICENSE) © 2026 Eric B., Ethan B., and Conrad B.
