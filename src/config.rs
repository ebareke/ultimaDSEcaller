//! Run configuration: merges defaults, optional config file (YAML / TOML /
//! JSON, detected by extension), and CLI flags into a single immutable
//! [`RunConfig`] passed to the rest of the pipeline.
//!
//! Also responsible for parsing the **sample sheet** — the tabular file that
//! describes which BAMs belong to which experimental groups and what
//! covariates apply. The schema is intentionally close to DESeq2's so users
//! can reuse existing files.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};

use crate::cli::{MultimapStrategy, RunArgs, Technology, TestMethod};
use crate::error::{UltiError, UltiResult};

/// One row of the sample sheet.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Sample {
    /// Unique sample identifier.
    pub id: String,
    /// Path to BAM or CRAM.
    pub bam: PathBuf,
    /// Primary experimental group (the column tested for differential usage).
    pub group: String,
    /// Any additional columns become covariates accessible by name. Reserved
    /// for future GLM / GLMM support.
    #[serde(default)]
    pub covariates: BTreeMap<String, String>,
}

/// User-facing thresholds for the "high-confidence" output set.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Filters {
    pub min_coverage: u32,
    pub min_supporting_reads: u32,
    pub max_fdr: f64,
    pub min_delta_psi: f64,
    pub min_confidence: f64,
}

impl Default for Filters {
    fn default() -> Self {
        Filters {
            min_coverage: 10,
            min_supporting_reads: 5,
            max_fdr: 0.05,
            min_delta_psi: 0.10,
            min_confidence: 0.5,
        }
    }
}

/// Read-extraction parameters.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReadParams {
    pub min_mapq: u8,
    pub min_overhang: u32,
    pub junction_tolerance: u32,
    pub multimap: MultimapStrategy,
}

impl ReadParams {
    pub fn for_tech(tech: Technology) -> Self {
        match tech {
            Technology::Short => ReadParams {
                min_mapq: 10,
                min_overhang: 8,
                junction_tolerance: 0,
                multimap: MultimapStrategy::Primary,
            },
            Technology::PacBio | Technology::Ont => ReadParams {
                min_mapq: 1,
                min_overhang: 4,
                junction_tolerance: 5,
                multimap: MultimapStrategy::Primary,
            },
        }
    }
}

/// Top-level immutable run configuration.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RunConfig {
    pub annotation: PathBuf,
    pub out: PathBuf,
    pub reference: Option<PathBuf>,
    pub tech: Technology,
    pub samples: Vec<Sample>,
    pub contrast: Option<Contrast>,
    pub filters: Filters,
    pub reads: ReadParams,
    pub threads: usize,
    pub test: TestMethod,
    pub random_effect: Option<String>,
    pub fusion_bedpe: Option<PathBuf>,
}

/// A differential-usage contrast: "group:level1-level2".
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Contrast {
    pub variable: String,
    pub numerator: String,
    pub denominator: String,
}

impl Contrast {
    pub fn parse(spec: &str) -> UltiResult<Self> {
        let (var, levels) = spec.split_once(':').ok_or_else(|| {
            UltiError::Config(format!(
                "contrast `{spec}` must look like `variable:level1-level2`"
            ))
        })?;
        let (num, denom) = levels.split_once('-').ok_or_else(|| {
            UltiError::Config(format!(
                "contrast `{spec}` must specify two levels separated by `-`"
            ))
        })?;
        Ok(Contrast {
            variable: var.trim().to_string(),
            numerator: num.trim().to_string(),
            denominator: denom.trim().to_string(),
        })
    }
}

/// Disk form of the config file. Every field is optional; missing fields fall
/// back to defaults or to CLI overrides.
#[derive(Debug, Default, Clone, Serialize, Deserialize)]
pub struct ConfigFile {
    pub annotation: Option<PathBuf>,
    pub sample_sheet: Option<PathBuf>,
    pub out: Option<PathBuf>,
    pub reference: Option<PathBuf>,
    pub tech: Option<Technology>,
    pub contrast: Option<String>,
    pub filters: Option<Filters>,
    pub reads: Option<ReadParams>,
    pub test: Option<TestMethod>,
    pub random_effect: Option<String>,
    pub fusion_bedpe: Option<PathBuf>,
}

impl ConfigFile {
    pub fn load(path: &Path) -> UltiResult<Self> {
        let raw = std::fs::read_to_string(path).map_err(|e| UltiError::io(path, e))?;
        let ext = path
            .extension()
            .and_then(|s| s.to_str())
            .unwrap_or("")
            .to_ascii_lowercase();
        match ext.as_str() {
            "yaml" | "yml" => {
                serde_yaml::from_str(&raw).map_err(|e| UltiError::Config(e.to_string()))
            }
            "toml" => toml::from_str(&raw).map_err(|e| UltiError::Config(e.to_string())),
            "json" => serde_json::from_str(&raw).map_err(|e| UltiError::Config(e.to_string())),
            other => Err(UltiError::Config(format!(
                "unknown config extension `{other}` — use .yaml, .toml, or .json"
            ))),
        }
    }
}

/// Sample-sheet parser. Auto-detects CSV vs TSV from the extension; falls
/// back to TSV.
pub fn read_sample_sheet(path: &Path) -> UltiResult<Vec<Sample>> {
    let delim = match path
        .extension()
        .and_then(|s| s.to_str())
        .unwrap_or("")
        .to_ascii_lowercase()
        .as_str()
    {
        "csv" => b',',
        _ => b'\t',
    };

    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(delim)
        .has_headers(true)
        .flexible(true)
        .from_path(path)
        .map_err(|e| UltiError::Design(format!("cannot open sample sheet {path:?}: {e}")))?;

    let headers = rdr
        .headers()
        .map_err(|e| UltiError::Design(format!("cannot read sample-sheet header: {e}")))?
        .clone();

    let required = ["sample", "bam", "group"];
    for req in required {
        if !headers.iter().any(|h| h.eq_ignore_ascii_case(req)) {
            return Err(UltiError::Design(format!(
                "sample sheet missing required column `{req}` (found: {:?})",
                headers.iter().collect::<Vec<_>>()
            )));
        }
    }

    let mut samples = Vec::new();
    for (row_idx, rec) in rdr.records().enumerate() {
        let rec = rec.map_err(|e| {
            UltiError::Design(format!("sample sheet row {row_idx} parse error: {e}"))
        })?;
        let mut sample_id = None;
        let mut bam = None;
        let mut group = None;
        let mut covariates = BTreeMap::new();
        for (h, v) in headers.iter().zip(rec.iter()) {
            match h.to_ascii_lowercase().as_str() {
                "sample" => sample_id = Some(v.to_string()),
                "bam" => bam = Some(PathBuf::from(v)),
                "group" => group = Some(v.to_string()),
                other => {
                    covariates.insert(other.to_string(), v.to_string());
                }
            }
        }
        samples.push(Sample {
            id: sample_id
                .ok_or_else(|| UltiError::Design(format!("row {row_idx}: missing `sample`")))?,
            bam: bam.ok_or_else(|| UltiError::Design(format!("row {row_idx}: missing `bam`")))?,
            group: group
                .ok_or_else(|| UltiError::Design(format!("row {row_idx}: missing `group`")))?,
            covariates,
        });
    }
    if samples.is_empty() {
        return Err(UltiError::Design("sample sheet is empty".into()));
    }
    Ok(samples)
}

/// Merge config-file values, CLI overrides, and defaults into a [`RunConfig`].
/// CLI flags take precedence over config-file values.
pub fn resolve(args: &RunArgs) -> UltiResult<RunConfig> {
    let cfg = match &args.config {
        Some(p) => ConfigFile::load(p)?,
        None => ConfigFile::default(),
    };

    let annotation = args
        .annotation
        .clone()
        .or(cfg.annotation)
        .ok_or_else(|| UltiError::Config("missing --annotation".into()))?;
    let out = args
        .out
        .clone()
        .or(cfg.out)
        .ok_or_else(|| UltiError::Config("missing --out".into()))?;
    let sample_sheet_path = args
        .sample_sheet
        .clone()
        .or(cfg.sample_sheet)
        .ok_or_else(|| UltiError::Config("missing --sample-sheet".into()))?;
    let samples = read_sample_sheet(&sample_sheet_path)?;

    let tech = args.tech.or(cfg.tech).unwrap_or(Technology::Short);

    let mut reads = cfg.reads.unwrap_or_else(|| ReadParams::for_tech(tech));
    if let Some(v) = args.min_mapq {
        reads.min_mapq = v;
    }
    if let Some(v) = args.min_overhang {
        reads.min_overhang = v;
    }
    if let Some(v) = args.junction_tolerance {
        reads.junction_tolerance = v;
    }
    if let Some(v) = args.multimap {
        reads.multimap = v;
    }

    let mut filters = cfg.filters.unwrap_or_default();
    if let Some(v) = args.max_fdr {
        filters.max_fdr = v;
    }
    if let Some(v) = args.min_delta_psi {
        filters.min_delta_psi = v;
    }
    if let Some(v) = args.min_coverage {
        filters.min_coverage = v;
    }

    let contrast = match args.contrast.clone().or(cfg.contrast) {
        Some(spec) => Some(Contrast::parse(&spec)?),
        None => None,
    };

    let test = args.test.or(cfg.test).unwrap_or(TestMethod::BbLrt);
    let random_effect = args.random_effect.clone().or(cfg.random_effect);

    Ok(RunConfig {
        annotation,
        out,
        reference: args.reference.clone().or(cfg.reference),
        tech,
        samples,
        contrast,
        filters,
        reads,
        threads: 0,
        test,
        random_effect,
        fusion_bedpe: args.fusion_bedpe.clone().or(cfg.fusion_bedpe),
    })
}
