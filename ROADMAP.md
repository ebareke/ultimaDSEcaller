# Roadmap

This roadmap is indicative, not a commitment. Items move as real-world
datasets and user feedback dictate.

## Near term

- **Benchmark suite** against rMATS, MAJIQ, SUPPA2, LeafCutter and FLAIR on
  simulated and validated datasets (recall@FDR, runtime, memory). This is
  experimental-design work needed for paper-level accuracy claims.
- **Per-event multi-test orchestration in the default path** — today the
  consensus engine runs every applicable test when `--consensus` is set; the
  goal is to make the combined q-value a first-class column in the headline
  tables.
- **Published container images** on GHCR and a Bioconda recipe so
  `conda install ultimadsecaller` and `docker pull ghcr.io/ebareke/...` work
  out of the box.

## Considered

- **Plot embedding in the PDF report** (currently tables-only; would require
  a headless-browser snapshot step for the Plotly figures).
- **GFF3 transcript→gene linkage** via `mRNA` records for annotations that
  omit `gene_id` on exon lines.
- **CRAM reference auto-resolution** from the BAM/CRAM header (`UR:`/`M5:`)
  rather than requiring `--reference`.
- **Dirichlet-multinomial** model for differential isoform usage (today: a
  χ² test on the K×2 table), better handling over-dispersion.
- **Second-order Laplace / adaptive Gauss-Hermite** for the GLMM (today:
  PQL), improving variance-component estimation in low-count regimes.
- **Full UMAP** (today: the spectral initializer) for large cohorts.
- **Streaming / out-of-core** junction and event tables for very large
  cohorts.
- **Splice-site strength scoring** (MaxEntScan-style) layered on top of the
  GT-AG / GC-AG / AT-AC motif classification.

## Explicitly out of scope

- Bundling or redistributing reference genomes / annotations.
- Acting as a general-purpose aligner — alignment is delegated to STAR and
  minimap2 in the pipeline and containers.
