//! # ultimaDSEcaller — alternative splicing analysis platform
//!
//! This crate provides a complete, production-grade Rust pipeline for
//! detecting, quantifying, and statistically testing alternative-splicing
//! events from short- and long-read RNA-seq data.
//!
//! ## Capabilities
//!
//! * All 7 canonical event types (SE, MXE, A5SS, A3SS, IR, AFE, ALE) plus
//!   11 advanced categories (cryptic, MSE, MIR, recursive, nested, partial
//!   exon, exonic-intronic hybrid, alt promoter, alt polyA, tandem UTR,
//!   fusion-associated).
//! * Multiple statistical tests (beta-binomial LRT, GLM, GLMM via PQL,
//!   Fisher 2×2) plus a Stouffer/Brown/weighted-Fisher consensus engine.
//! * Reference-FASTA splice-site motif annotation (GT-AG / GC-AG / AT-AC).
//! * Protein consequence prediction with NMD-rule scoring.
//! * Per-event sashimi tracks, sample heatmaps, PCA, UMAP, isoform Sankey
//!   diagrams — all in a single self-contained HTML report.
//! * Annotation binary cache, BAM-index acceleration, checkpoint/resume,
//!   and progress bars.
//!
//! ## Quick start
//!
//! See the crate's `README.md` for the full CLI guide. From Rust:
//!
//! ```no_run
//! use ultimadsecaller::{annotation, events, junctions, quantify, stats};
//! # fn main() -> anyhow::Result<()> {
//! let ann = annotation::parse(std::path::Path::new("genome.gtf"))?;
//! # let _ = ann;
//! # Ok(())
//! # }
//! ```
//!
//! ## Module map
//!
//! | Module | Purpose |
//! |--------|---------|
//! | [`annotation`] | GTF / GFF3 parser + per-gene splice graphs |
//! | [`cache`] | Binary annotation cache |
//! | [`junctions`] | BAM/CRAM junction extraction (streaming + indexed) |
//! | [`pileup`] | Region pileup for coverage and IR |
//! | [`events`] | Canonical event detection |
//! | [`advanced`] | Cryptic, MSE, MIR, recursive, nested, etc. |
//! | [`quantify`] | PSI / ΔPSI / coverage / confidence per event |
//! | [`stats`] | BB-LRT, BH-FDR, Fisher 2×2 |
//! | [`glm`] | GLM (IRLS) and GLMM (PQL) |
//! | [`consensus`] | Combine multiple tests into one q-value |
//! | [`formula`] | Wilkinson-style design formula parser |
//! | [`motif`] | Reference FASTA + splice-site motif classification |
//! | [`protein`] | Protein consequence + NMD prediction |
//! | [`embedding`] | PCA + UMAP-like sample embeddings |
//! | [`longread`] | Long-read isoform reconstruction + DIU |
//! | [`sashimi`] | Sashimi tracks + isoform Sankey data |
//! | [`output`] | TSV / CSV / JSON / Parquet writers + report payload |
//! | [`report`] | HTML report rendering |
//! | [`progress`] | Progress-bar facade |
//! | [`error`] | Crate-wide error type with stable error codes |
//!
//! ## Error codes
//!
//! All errors expose a stable machine-readable code via [`error::UltiError::code`].
//! See the [`error`] module for the full list.

pub mod advanced;
pub mod annotation;
pub mod cache;
pub mod cli;
pub mod config;
pub mod consensus;
pub mod embedding;
pub mod error;
pub mod events;
pub mod formula;
pub mod glm;
pub mod junctions;
pub mod longread;
pub mod motif;
pub mod output;
#[cfg(feature = "pdf")]
pub mod pdf;
pub mod pileup;
pub mod progress;
pub mod protein;
pub mod quantify;
pub mod report;
pub mod sashimi;
pub mod stats;

pub use error::{UltiError, UltiResult};

/// Genomic coordinate (1-based inclusive, GTF convention).
pub type Position = u64;

/// Strand of a feature.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, serde::Serialize, serde::Deserialize)]
pub enum Strand {
    Forward,
    Reverse,
    Unknown,
}

impl Strand {
    pub fn from_char(c: char) -> Self {
        match c {
            '+' => Strand::Forward,
            '-' => Strand::Reverse,
            _ => Strand::Unknown,
        }
    }
}

impl std::fmt::Display for Strand {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Strand::Forward => "+",
            Strand::Reverse => "-",
            Strand::Unknown => ".",
        })
    }
}
