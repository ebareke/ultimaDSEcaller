//! Command-line interface. The CLI is a thin wrapper that translates user
//! input into a [`crate::config::RunConfig`] and dispatches to the pipeline.

use std::path::PathBuf;

use clap::{Parser, Subcommand, ValueEnum};

/// `ultimaDSEcaller` — Alternative splicing detection, quantification and
/// differential analysis for short- and long-read RNA-seq.
#[derive(Parser, Debug)]
#[command(
    name = "ultimaDSEcaller",
    version,
    about,
    long_about = None,
    propagate_version = true,
)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Command,

    /// Increase log verbosity (-v, -vv, -vvv).
    #[arg(short, long, action = clap::ArgAction::Count, global = true)]
    pub verbose: u8,

    /// Number of worker threads for the rayon thread pool. 0 = all logical cores.
    #[arg(long, default_value_t = 0, global = true)]
    pub threads: usize,
}

#[derive(Subcommand, Debug)]
// clap parses this once at startup; variant size difference is irrelevant.
#[allow(clippy::large_enum_variant)]
pub enum Command {
    /// End-to-end pipeline: detect → quantify → test → report.
    Run(RunArgs),

    /// Parse GTF/GFF3 and dump the splice graph as JSON (debugging aid).
    DumpAnnotation(DumpAnnotationArgs),

    /// Extract junctions from one or more BAM/CRAM files without calling events.
    Junctions(JunctionsArgs),

    /// Render an HTML report from a previously-produced results JSON file.
    Report(ReportArgs),

    /// Render a PDF report from a previously-produced results JSON file.
    /// Requires the binary to have been built with `--features pdf`.
    Pdf(ReportArgs),
}

#[derive(Parser, Debug)]
pub struct RunArgs {
    /// Path to a YAML/TOML/JSON config file. CLI flags override config values.
    #[arg(short = 'c', long)]
    pub config: Option<PathBuf>,

    /// Sample sheet (TSV or CSV). Columns: sample, bam, group [, batch, ...].
    /// See README for the full schema.
    #[arg(short = 's', long)]
    pub sample_sheet: Option<PathBuf>,

    /// Annotation file (GTF or GFF3, optionally bgzipped). Format is detected by extension.
    #[arg(short = 'a', long)]
    pub annotation: Option<PathBuf>,

    /// Output directory. Created if it does not exist.
    #[arg(short = 'o', long)]
    pub out: Option<PathBuf>,

    /// Reference FASTA (required for CRAM input, optional for BAM).
    #[arg(long)]
    pub reference: Option<PathBuf>,

    /// Sequencing technology — affects junction tolerance, MAPQ defaults, and
    /// whether long-read isoform reconstruction is enabled.
    #[arg(long, value_enum)]
    pub tech: Option<Technology>,

    /// Contrast for differential analysis, e.g. "group:treatment-control".
    /// If omitted, all pairwise contrasts of the design column are emitted.
    #[arg(long)]
    pub contrast: Option<String>,

    /// Minimum MAPQ for a read to be counted toward junction support.
    #[arg(long)]
    pub min_mapq: Option<u8>,

    /// Minimum splice-junction overhang on each side (anchor length).
    #[arg(long)]
    pub min_overhang: Option<u32>,

    /// Strategy for multimapping reads.
    #[arg(long, value_enum)]
    pub multimap: Option<MultimapStrategy>,

    /// Junction wobble tolerance (bp). 0 disables collapsing; 3–5 typical for long reads.
    #[arg(long)]
    pub junction_tolerance: Option<u32>,

    /// Maximum FDR for the "high-confidence" output set.
    #[arg(long)]
    pub max_fdr: Option<f64>,

    /// Minimum |ΔPSI| for the "high-confidence" output set.
    #[arg(long)]
    pub min_delta_psi: Option<f64>,

    /// Minimum total inclusion+exclusion coverage per sample for inclusion in the test.
    #[arg(long)]
    pub min_coverage: Option<u32>,

    /// Statistical test to use. `bb-lrt` is the default. `glm` uses logistic
    /// GLM with a single binary contrast term. `glmm` adds a random
    /// intercept per level of a chosen grouping covariate (set via
    /// `--random-effect`).
    #[arg(long, value_enum)]
    pub test: Option<TestMethod>,

    /// Name of a sample-sheet column to use as the random-effect grouping
    /// factor when `--test glmm` is set (e.g. `subject`, `donor`, `batch`).
    #[arg(long)]
    pub random_effect: Option<String>,

    /// Optional fusion BEDPE for FusionAssociated event detection.
    #[arg(long)]
    pub fusion_bedpe: Option<PathBuf>,

    /// Wilkinson-style design formula (e.g. `~ batch + group + batch:group`).
    /// When set, the GLM/GLMM call uses this design instead of the default
    /// `[1, treatment]`. The right-hand side variables must exist in the
    /// sample sheet as `group`, `sample`, or covariate columns.
    #[arg(long)]
    pub formula: Option<String>,

    /// Multiple contrasts, each as `variable:level1-level2`, comma-separated.
    /// Each contrast produces its own per-contrast output subdirectory.
    /// Mutually exclusive with `--contrast`.
    #[arg(long, value_delimiter = ',')]
    pub contrasts: Vec<String>,

    /// Consensus combination method for multi-test events.
    #[arg(long, value_enum)]
    pub consensus: Option<ConsensusMethod>,

    /// Skip stages whose outputs already exist (junction matrix, results JSON, ...).
    #[arg(long)]
    pub resume: bool,

    /// Disable the annotation binary cache (always re-parse).
    #[arg(long)]
    pub no_cache: bool,

    /// Path for the annotation cache (default: `<annotation>.ultidse`).
    #[arg(long)]
    pub cache_path: Option<PathBuf>,

    /// Also emit a PDF report (`report.pdf`) alongside the HTML report.
    /// Requires the binary to have been built with `--features pdf`.
    #[arg(long)]
    pub pdf: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum TestMethod {
    /// Beta-binomial likelihood-ratio test (default).
    BbLrt,
    /// Logistic GLM via IRLS, single binary contrast.
    Glm,
    /// Logistic GLMM with random intercept (PQL).
    Glmm,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum ConsensusMethod {
    Stouffer,
    Brown,
    WeightedFisher,
}

impl From<ConsensusMethod> for crate::consensus::ConsensusMethod {
    fn from(c: ConsensusMethod) -> Self {
        match c {
            ConsensusMethod::Stouffer => crate::consensus::ConsensusMethod::Stouffer,
            ConsensusMethod::Brown => crate::consensus::ConsensusMethod::Brown,
            ConsensusMethod::WeightedFisher => crate::consensus::ConsensusMethod::WeightedFisher,
        }
    }
}

#[derive(Parser, Debug)]
pub struct DumpAnnotationArgs {
    #[arg(short = 'a', long)]
    pub annotation: PathBuf,

    /// Restrict to a single gene by gene_id.
    #[arg(long)]
    pub gene: Option<String>,

    /// Output path (default: stdout).
    #[arg(short = 'o', long)]
    pub out: Option<PathBuf>,
}

#[derive(Parser, Debug)]
pub struct JunctionsArgs {
    /// One or more BAM/CRAM files.
    #[arg(required = true)]
    pub bams: Vec<PathBuf>,

    /// Output TSV path.
    #[arg(short = 'o', long)]
    pub out: PathBuf,

    /// Reference FASTA (required for CRAM).
    #[arg(long)]
    pub reference: Option<PathBuf>,

    #[arg(long, default_value_t = 10)]
    pub min_mapq: u8,

    #[arg(long, default_value_t = 8)]
    pub min_overhang: u32,
}

#[derive(Parser, Debug)]
pub struct ReportArgs {
    /// Input results JSON (produced by `run`).
    #[arg(short = 'i', long)]
    pub input: PathBuf,

    /// Output HTML path.
    #[arg(short = 'o', long)]
    pub out: PathBuf,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum Technology {
    /// Illumina short reads (default behaviour).
    Short,
    /// Long reads, PacBio Iso-Seq.
    PacBio,
    /// Long reads, Oxford Nanopore.
    Ont,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "kebab-case")]
#[derive(Default)]
pub enum MultimapStrategy {
    /// Discard reads with NH > 1 (most conservative; mirrors rMATS default).
    Discard,
    /// Use a primary-alignment-only count (recommended for short reads).
    #[default]
    Primary,
    /// Weight each alignment by 1/NH (fractional counting).
    Fractional,
}
