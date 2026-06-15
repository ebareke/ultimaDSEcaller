# Changelog

All notable changes to ultimaDSEcaller are documented here. The format
follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the
project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

## [1.0.0] — 2026-06-15

First public release.

### Added

- **Event detection** for 18 alternative-splicing classes: the 7 canonical
  events (SE, MXE, A5SS, A3SS, IR, AFE, ALE) and 11 advanced events (cryptic
  splice sites, multi-exon skipping, multi-intron retention, recursive and
  nested splicing, partial-exon inclusion, exonic-intronic hybrids,
  alternative promoter and polyadenylation, tandem UTRs, fusion-associated
  junctions).
- **Short- and long-read support**: junction-based PSI/ΔPSI for Illumina;
  full-length isoform reconstruction with FLAIR-style annotation-anchored
  junction correction and differential isoform usage for PacBio / ONT.
- **Statistical engines**: beta-binomial likelihood-ratio test (Nelder-Mead
  MLE + χ²), logistic GLM (IRLS) and GLMM (PQL random intercept), Fisher's
  exact 2×2, Benjamini-Hochberg FDR, and a consensus engine
  (Stouffer / Brown / weighted-Fisher) combining all available tests.
- **Functional annotation** from a reference FASTA: GT-AG / GC-AG / AT-AC
  splice-motif classification, protein-consequence prediction
  (in-frame / frame-shift / PTC) and the 50-nt NMD rule.
- **Design support**: multi-group sample sheets, Wilkinson-style design
  formulas (`~ batch + group + batch:group`), and multiple contrasts in one
  run.
- **Reporting**: a self-contained interactive HTML report (volcano, MA, PCA,
  UMAP, PSI distributions, per-sample heatmap, coverage and chromosomal
  distributions, splice-junction graphs, per-position sashimi tracks, isoform
  Sankey diagrams, consensus and protein-consequence summaries) plus PDF.
- **Output formats**: TSV, CSV, JSON, Parquet (Snappy), HTML, PDF.
- **Performance & robustness**: rayon multithreading, BAM-index-aware region
  pileup, binary annotation cache, checkpoint/resume, stable machine-readable
  error codes, and progress reporting.
- **Reproducible deployment**: a multi-stage Docker image bundling the static
  caller with minimap2, STAR, HISAT2 and samtools; an Apptainer definition;
  and a Nextflow pipeline that runs the full FASTQ → events workflow.
- **Tests**: integration tests, property-based tests (proptest), and
  criterion benchmarks. A synthetic end-to-end example that runs without an
  aligner.

### Heritage

ultimaDSEcaller began as an internal alternative-splicing caller and was
consolidated, hardened, and rebranded for this release. Compared with that
prototype, the statistical engines, long-read support, functional
annotation, consensus combination, container/Nextflow deployment, and the
full report suite are new or substantially rewritten.
