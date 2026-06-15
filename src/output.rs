//! Output writers — TSV, CSV, JSON, and (behind the `parquet` feature) a
//! columnar Parquet table.

use std::collections::BTreeMap;
use std::fs::{self, File};
use std::io::BufWriter;
use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};

use crate::config::{Filters, RunConfig};
use crate::error::{UltiError, UltiResult};
use crate::events::EventKind;
use crate::quantify::EventQuant;
use crate::stats::{PValue, TestUsed};

/// Flat per-event result row — what gets written to the headline TSV / CSV.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ResultRow {
    pub event_id: String,
    pub gene_id: String,
    pub chrom: String,
    pub strand: String,
    pub event_type: String,
    pub coordinates: String,
    pub n_samples: usize,
    pub mean_psi_numerator: f64,
    pub mean_psi_denominator: f64,
    pub delta_psi: f64,
    pub mean_coverage_numerator: f64,
    pub mean_coverage_denominator: f64,
    pub inclusion_reads_total: f64,
    pub exclusion_reads_total: f64,
    pub complexity: u32,
    pub reproducibility: f64,
    pub confidence: f64,
    pub p_value: f64,
    pub adjusted_p_value: f64,
    pub test_used: String,
}

pub fn build_rows(quants: &[EventQuant], pvals: &[PValue]) -> Vec<ResultRow> {
    quants
        .iter()
        .zip(pvals.iter())
        .map(|(q, p)| {
            let cs = q.contrast_summary.as_ref();
            ResultRow {
                event_id: q.event.event_id.clone(),
                gene_id: q.event.gene_id.clone(),
                chrom: q.event.chrom.clone(),
                strand: q.event.strand.to_string(),
                event_type: q.event.kind.short().to_string(),
                coordinates: q
                    .event
                    .exons
                    .iter()
                    .map(|e| format!("{}-{}", e.start, e.end))
                    .collect::<Vec<_>>()
                    .join(","),
                n_samples: q.psi.len(),
                mean_psi_numerator: cs.map(|s| s.mean_psi_num).unwrap_or(f64::NAN),
                mean_psi_denominator: cs.map(|s| s.mean_psi_denom).unwrap_or(f64::NAN),
                delta_psi: cs.map(|s| s.delta_psi).unwrap_or(f64::NAN),
                mean_coverage_numerator: cs.map(|s| s.mean_coverage_num).unwrap_or(f64::NAN),
                mean_coverage_denominator: cs.map(|s| s.mean_coverage_denom).unwrap_or(f64::NAN),
                inclusion_reads_total: q.inclusion.iter().sum(),
                exclusion_reads_total: q.exclusion.iter().sum(),
                complexity: q.complexity,
                reproducibility: q.reproducibility,
                confidence: q.confidence,
                p_value: p.p_value,
                adjusted_p_value: p.adjusted_p_value,
                test_used: match p.test_used {
                    TestUsed::BetaBinomialLRT => "BB-LRT".into(),
                    TestUsed::FisherExact => "Fisher-2x2".into(),
                    TestUsed::Glm => "GLM".into(),
                    TestUsed::Glmm => "GLMM".into(),
                    TestUsed::Insufficient => "insufficient".into(),
                },
            }
        })
        .collect()
}

pub fn passes_filter(row: &ResultRow, f: &Filters) -> bool {
    if !row.adjusted_p_value.is_finite() || row.adjusted_p_value > f.max_fdr {
        return false;
    }
    if !row.delta_psi.is_finite() || row.delta_psi.abs() < f.min_delta_psi {
        return false;
    }
    if row.confidence < f.min_confidence {
        return false;
    }
    let mean_cov = 0.5
        * (row.mean_coverage_numerator + row.mean_coverage_denominator);
    if mean_cov < f.min_coverage as f64 {
        return false;
    }
    true
}

pub fn write_all(
    cfg: &RunConfig,
    rows: &[ResultRow],
) -> UltiResult<()> {
    write_all_with_payload(cfg, rows, |payload| payload)
}

/// Same as [`write_all`] but lets the caller augment the report payload
/// (PCA, UMAP, sashimi, isoform Sankey, advanced-event counts) before
/// serialization. The closure receives a freshly-built payload and
/// returns the version to write.
pub fn write_all_with_payload<F>(
    cfg: &RunConfig,
    rows: &[ResultRow],
    augment: F,
) -> UltiResult<()>
where
    F: FnOnce(ReportPayload) -> ReportPayload,
{
    let out = &cfg.out;
    fs::create_dir_all(out).map_err(|e| UltiError::io(out, e))?;

    // Headline tables — full set and high-confidence subset, in both TSV and CSV.
    write_table(&out.join("events.raw.tsv"), rows, b'\t')?;
    write_table(&out.join("events.raw.csv"), rows, b',')?;

    let hc: Vec<ResultRow> = rows
        .iter()
        .filter(|r| passes_filter(r, &cfg.filters))
        .cloned()
        .collect();
    write_table(&out.join("events.high_confidence.tsv"), &hc, b'\t')?;
    write_table(&out.join("events.high_confidence.csv"), &hc, b',')?;

    // Per-event-type breakouts. Each AS category gets its own TSV — both
    // raw and high-confidence. Easy for users to grep / load into R.
    let per_type_dir = out.join("by_event_type");
    fs::create_dir_all(&per_type_dir).map_err(|e| UltiError::io(&per_type_dir, e))?;
    for kind in [
        EventKind::SE,
        EventKind::MXE,
        EventKind::A5SS,
        EventKind::A3SS,
        EventKind::IR,
        EventKind::AFE,
        EventKind::ALE,
    ] {
        let subset_raw: Vec<_> = rows
            .iter()
            .filter(|r| r.event_type == kind.short())
            .cloned()
            .collect();
        let subset_hc: Vec<_> = hc
            .iter()
            .filter(|r| r.event_type == kind.short())
            .cloned()
            .collect();
        write_table(
            &per_type_dir.join(format!("{}.raw.tsv", kind.short())),
            &subset_raw,
            b'\t',
        )?;
        write_table(
            &per_type_dir.join(format!("{}.high_confidence.tsv", kind.short())),
            &subset_hc,
            b'\t',
        )?;
    }

    // Optional Parquet output.
    #[cfg(feature = "parquet")]
    {
        let path = out.join("events.raw.parquet");
        write_parquet(&path, rows)?;
    }

    // Machine-readable summary consumed by the HTML report.
    let summary = augment(build_report_payload(cfg, rows, &hc));
    let json_path = out.join("results.json");
    let f = File::create(&json_path).map_err(|e| UltiError::io(&json_path, e))?;
    // serde_json serializes non-finite floats (NaN/Inf — e.g. the mean ΔPSI of
    // an event kind with no events) as JSON `null`. That is valid for the
    // browser and for pandas, but the strict `f64` fields reject `null` on
    // read-back (the `report` / `pdf` subcommands). Sanitize numeric nulls to
    // 0.0 while preserving genuinely-optional fields, so the JSON round-trips
    // for every consumer.
    let mut value =
        serde_json::to_value(&summary).map_err(|e| UltiError::Other(anyhow::anyhow!(e)))?;
    sanitize_numeric_nulls(&mut value);
    let writer = BufWriter::new(f);
    serde_json::to_writer_pretty(writer, &value)
        .map_err(|e| UltiError::Other(anyhow::anyhow!(e)))?;

    Ok(())
}

/// Keys whose value is legitimately optional and may stay `null`.
const NULLABLE_KEYS: &[&str] = &[
    "contrast",
    "pca",
    "umap",
    "isoform_sankey",
    "heatmap",
    "explained",
];

/// Recursively replace numeric `null`s (from non-finite floats) with `0.0`,
/// leaving the known-nullable object fields untouched. Array elements that are
/// `null` (e.g. gaps in a heatmap PSI row) become `0.0`.
fn sanitize_numeric_nulls(value: &mut serde_json::Value) {
    use serde_json::Value;
    match value {
        Value::Object(map) => {
            for (k, v) in map.iter_mut() {
                if v.is_null() {
                    if !NULLABLE_KEYS.contains(&k.as_str()) {
                        *v = Value::from(0.0);
                    }
                } else {
                    sanitize_numeric_nulls(v);
                }
            }
        }
        Value::Array(arr) => {
            for v in arr.iter_mut() {
                if v.is_null() {
                    *v = Value::from(0.0);
                } else {
                    sanitize_numeric_nulls(v);
                }
            }
        }
        _ => {}
    }
}

#[cfg(feature = "parquet")]
fn write_parquet(path: &Path, rows: &[ResultRow]) -> UltiResult<()> {
    use parquet::basic::Compression;
    use parquet::data_type::ByteArray;
    use parquet::file::properties::WriterProperties;
    use parquet::file::writer::SerializedFileWriter;
    use parquet::schema::parser::parse_message_type;
    use std::sync::Arc;

    let schema_str = "
        message events {
            REQUIRED BYTE_ARRAY event_id (UTF8);
            REQUIRED BYTE_ARRAY gene_id (UTF8);
            REQUIRED BYTE_ARRAY chrom (UTF8);
            REQUIRED BYTE_ARRAY strand (UTF8);
            REQUIRED BYTE_ARRAY event_type (UTF8);
            REQUIRED BYTE_ARRAY coordinates (UTF8);
            REQUIRED INT64 n_samples;
            REQUIRED DOUBLE mean_psi_numerator;
            REQUIRED DOUBLE mean_psi_denominator;
            REQUIRED DOUBLE delta_psi;
            REQUIRED DOUBLE mean_coverage_numerator;
            REQUIRED DOUBLE mean_coverage_denominator;
            REQUIRED DOUBLE inclusion_reads_total;
            REQUIRED DOUBLE exclusion_reads_total;
            REQUIRED INT64 complexity;
            REQUIRED DOUBLE reproducibility;
            REQUIRED DOUBLE confidence;
            REQUIRED DOUBLE p_value;
            REQUIRED DOUBLE adjusted_p_value;
            REQUIRED BYTE_ARRAY test_used (UTF8);
        }
    ";
    let schema = Arc::new(parse_message_type(schema_str).map_err(|e| {
        UltiError::Other(anyhow::anyhow!("parquet schema parse: {e}"))
    })?);
    let props = Arc::new(
        WriterProperties::builder()
            // Snappy is the codec enabled in Cargo.toml (`snap` feature);
            // it is fast and widely supported by Parquet readers.
            .set_compression(Compression::SNAPPY)
            .build(),
    );
    let file = std::fs::File::create(path).map_err(|e| UltiError::io(path, e))?;
    let mut writer = SerializedFileWriter::new(file, schema, props)
        .map_err(|e| UltiError::Other(anyhow::anyhow!(e)))?;

    let mut row_group = writer
        .next_row_group()
        .map_err(|e| UltiError::Other(anyhow::anyhow!(e)))?;
    macro_rules! write_string_col {
        ($field:ident) => {{
            let values: Vec<ByteArray> = rows
                .iter()
                .map(|r| ByteArray::from(r.$field.as_bytes()))
                .collect();
            let mut col = row_group
                .next_column()
                .map_err(|e| UltiError::Other(anyhow::anyhow!(e)))?
                .ok_or_else(|| UltiError::Other(anyhow::anyhow!("missing column")))?;
            col.typed::<parquet::data_type::ByteArrayType>()
                .write_batch(&values, None, None)
                .map_err(|e| UltiError::Other(anyhow::anyhow!(e)))?;
            col.close().map_err(|e| UltiError::Other(anyhow::anyhow!(e)))?;
        }};
    }
    macro_rules! write_i64_col {
        ($expr:expr) => {{
            let values: Vec<i64> = rows.iter().map($expr).collect();
            let mut col = row_group
                .next_column()
                .map_err(|e| UltiError::Other(anyhow::anyhow!(e)))?
                .ok_or_else(|| UltiError::Other(anyhow::anyhow!("missing column")))?;
            col.typed::<parquet::data_type::Int64Type>()
                .write_batch(&values, None, None)
                .map_err(|e| UltiError::Other(anyhow::anyhow!(e)))?;
            col.close().map_err(|e| UltiError::Other(anyhow::anyhow!(e)))?;
        }};
    }
    macro_rules! write_f64_col {
        ($expr:expr) => {{
            let values: Vec<f64> = rows.iter().map($expr).collect();
            let mut col = row_group
                .next_column()
                .map_err(|e| UltiError::Other(anyhow::anyhow!(e)))?
                .ok_or_else(|| UltiError::Other(anyhow::anyhow!("missing column")))?;
            col.typed::<parquet::data_type::DoubleType>()
                .write_batch(&values, None, None)
                .map_err(|e| UltiError::Other(anyhow::anyhow!(e)))?;
            col.close().map_err(|e| UltiError::Other(anyhow::anyhow!(e)))?;
        }};
    }

    write_string_col!(event_id);
    write_string_col!(gene_id);
    write_string_col!(chrom);
    write_string_col!(strand);
    write_string_col!(event_type);
    write_string_col!(coordinates);
    write_i64_col!(|r| r.n_samples as i64);
    write_f64_col!(|r| r.mean_psi_numerator);
    write_f64_col!(|r| r.mean_psi_denominator);
    write_f64_col!(|r| r.delta_psi);
    write_f64_col!(|r| r.mean_coverage_numerator);
    write_f64_col!(|r| r.mean_coverage_denominator);
    write_f64_col!(|r| r.inclusion_reads_total);
    write_f64_col!(|r| r.exclusion_reads_total);
    write_i64_col!(|r| r.complexity as i64);
    write_f64_col!(|r| r.reproducibility);
    write_f64_col!(|r| r.confidence);
    write_f64_col!(|r| r.p_value);
    write_f64_col!(|r| r.adjusted_p_value);
    write_string_col!(test_used);

    row_group
        .close()
        .map_err(|e| UltiError::Other(anyhow::anyhow!(e)))?;
    writer
        .close()
        .map_err(|e| UltiError::Other(anyhow::anyhow!(e)))?;
    Ok(())
}

fn write_table(path: &Path, rows: &[ResultRow], delimiter: u8) -> UltiResult<()> {
    let mut wtr = csv::WriterBuilder::new()
        .delimiter(delimiter)
        .from_path(path)
        .map_err(|e| UltiError::Other(anyhow::anyhow!(e)))?;
    for r in rows {
        wtr.serialize(r)
            .map_err(|e| UltiError::Other(anyhow::anyhow!(e)))?;
    }
    wtr.flush().map_err(|e| UltiError::io(path, e))
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReportPayload {
    pub run: RunMeta,
    pub summary: SummaryStats,
    pub per_kind: BTreeMap<String, PerKindStats>,
    pub volcano: Vec<VolcanoPoint>,
    pub top_events: BTreeMap<String, Vec<ResultRow>>,
    pub high_confidence: Vec<ResultRow>,
    /// Extra data filled in by the pipeline driver after the headline tables
    /// are computed. Optional so the `report` subcommand on a partial JSON
    /// still works.
    #[serde(default)]
    pub ma_points: Vec<MaPoint>,
    #[serde(default)]
    pub pca: Option<crate::embedding::Embedding2D>,
    #[serde(default)]
    pub umap: Option<crate::embedding::Embedding2D>,
    #[serde(default)]
    pub sashimi: Vec<crate::sashimi::SashimiTrack>,
    #[serde(default)]
    pub isoform_sankey: Option<crate::sashimi::SankeyData>,
    #[serde(default)]
    pub advanced_event_counts: BTreeMap<String, usize>,
    #[serde(default)]
    pub diu: Vec<crate::longread::DiuRecord>,
    /// Optional heatmap of top events × samples (PSI values).
    #[serde(default)]
    pub heatmap: Option<HeatmapData>,
    /// Per-sample read-count summary for the coverage distribution plot.
    #[serde(default)]
    pub coverage_distribution: Vec<CoverageRow>,
    /// Counts of detected events per chromosome.
    #[serde(default)]
    pub chrom_distribution: Vec<ChromCount>,
    /// Compact splice-junction graphs for the top genes by significance.
    #[serde(default)]
    pub junction_graphs: Vec<JunctionGraph>,
    /// Consensus statistical results (optional — emitted when consensus is run).
    #[serde(default)]
    pub consensus: Vec<crate::consensus::ConsensusResult>,
    /// Per-event protein-consequence calls (optional — emitted when a
    /// reference FASTA and CDS catalog are available).
    #[serde(default)]
    pub protein_consequences: Vec<ProteinAnnotation>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HeatmapData {
    pub event_ids: Vec<String>,
    pub samples: Vec<String>,
    /// Row-major matrix: row = event, col = sample. NaN where uncovered.
    pub psi: Vec<Vec<f64>>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CoverageRow {
    pub sample: String,
    pub group: String,
    pub total_passing_reads: u64,
    pub low_mapq_reads: u64,
    pub mean_junction_support: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ChromCount {
    pub chrom: String,
    pub total: usize,
    pub significant: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct JunctionGraph {
    pub gene_id: String,
    pub chrom: String,
    pub nodes: Vec<JunctionGraphNode>,
    pub edges: Vec<JunctionGraphEdge>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct JunctionGraphNode {
    pub id: usize,
    pub start: u64,
    pub end: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct JunctionGraphEdge {
    pub source: usize,
    pub target: usize,
    pub transcripts: u32,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProteinAnnotation {
    pub event_id: String,
    pub gene_id: String,
    pub consequence: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MaPoint {
    pub event_id: String,
    pub event_type: String,
    pub gene_id: String,
    pub mean_log_total: f64,
    pub log_fold_change: f64,
    pub neg_log10_q: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RunMeta {
    pub n_samples: usize,
    pub n_groups: usize,
    pub contrast: Option<String>,
    pub annotation_path: PathBuf,
    pub tech: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SummaryStats {
    pub events_total: usize,
    pub events_high_confidence: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PerKindStats {
    pub total: usize,
    pub significant: usize,
    pub mean_delta_psi: f64,
    pub psi_numerator_histogram: Vec<(f64, usize)>,
    pub psi_denominator_histogram: Vec<(f64, usize)>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VolcanoPoint {
    pub event_id: String,
    pub event_type: String,
    pub gene_id: String,
    pub delta_psi: f64,
    pub neg_log10_q: f64,
}

pub fn build_report_payload(
    cfg: &RunConfig,
    rows: &[ResultRow],
    hc: &[ResultRow],
) -> ReportPayload {
    let mut per_kind: BTreeMap<String, PerKindStats> = BTreeMap::new();
    for k in ["SE", "MXE", "A5SS", "A3SS", "IR", "AFE", "ALE"] {
        let subset: Vec<&ResultRow> = rows.iter().filter(|r| r.event_type == k).collect();
        let sig = subset
            .iter()
            .filter(|r| r.adjusted_p_value.is_finite() && r.adjusted_p_value <= cfg.filters.max_fdr)
            .count();
        let mean_dpsi = mean_finite(subset.iter().map(|r| r.delta_psi));
        let hist_num = histogram(subset.iter().map(|r| r.mean_psi_numerator), 0.0, 1.0, 20);
        let hist_denom = histogram(subset.iter().map(|r| r.mean_psi_denominator), 0.0, 1.0, 20);
        per_kind.insert(
            k.to_string(),
            PerKindStats {
                total: subset.len(),
                significant: sig,
                mean_delta_psi: mean_dpsi,
                psi_numerator_histogram: hist_num,
                psi_denominator_histogram: hist_denom,
            },
        );
    }

    let volcano: Vec<VolcanoPoint> = rows
        .iter()
        .filter(|r| r.adjusted_p_value.is_finite())
        .map(|r| VolcanoPoint {
            event_id: r.event_id.clone(),
            event_type: r.event_type.clone(),
            gene_id: r.gene_id.clone(),
            delta_psi: r.delta_psi,
            neg_log10_q: -r.adjusted_p_value.max(1e-300).log10(),
        })
        .collect();

    let mut top_events: BTreeMap<String, Vec<ResultRow>> = BTreeMap::new();
    for k in ["SE", "MXE", "A5SS", "A3SS", "IR", "AFE", "ALE"] {
        let mut sorted: Vec<ResultRow> = rows
            .iter()
            .filter(|r| r.event_type == k && r.adjusted_p_value.is_finite())
            .cloned()
            .collect();
        sorted.sort_by(|a, b| {
            a.adjusted_p_value
                .partial_cmp(&b.adjusted_p_value)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        sorted.truncate(20);
        top_events.insert(k.to_string(), sorted);
    }

    let n_groups: usize = cfg
        .samples
        .iter()
        .map(|s| s.group.clone())
        .collect::<std::collections::HashSet<_>>()
        .len();
    let ma_points = rows
        .iter()
        .filter(|r| r.adjusted_p_value.is_finite())
        .map(|r| {
            let total =
                (r.mean_coverage_numerator + r.mean_coverage_denominator).max(1e-6);
            let p_num = r.mean_psi_numerator.max(1e-6).min(1.0 - 1e-6);
            let p_denom = r.mean_psi_denominator.max(1e-6).min(1.0 - 1e-6);
            let log_fc = (p_num / (1.0 - p_num)).ln() - (p_denom / (1.0 - p_denom)).ln();
            MaPoint {
                event_id: r.event_id.clone(),
                event_type: r.event_type.clone(),
                gene_id: r.gene_id.clone(),
                mean_log_total: total.ln(),
                log_fold_change: log_fc,
                neg_log10_q: -r.adjusted_p_value.max(1e-300).log10(),
            }
        })
        .collect();

    // Chromosomal distribution from headline rows.
    let mut by_chrom: BTreeMap<String, (usize, usize)> = BTreeMap::new();
    for r in rows {
        let entry = by_chrom.entry(r.chrom.clone()).or_insert((0, 0));
        entry.0 += 1;
        if r.adjusted_p_value.is_finite() && r.adjusted_p_value <= cfg.filters.max_fdr {
            entry.1 += 1;
        }
    }
    let chrom_distribution: Vec<ChromCount> = by_chrom
        .into_iter()
        .map(|(chrom, (total, sig))| ChromCount {
            chrom,
            total,
            significant: sig,
        })
        .collect();

    ReportPayload {
        run: RunMeta {
            n_samples: cfg.samples.len(),
            n_groups,
            contrast: cfg
                .contrast
                .as_ref()
                .map(|c| format!("{}:{}-{}", c.variable, c.numerator, c.denominator)),
            annotation_path: cfg.annotation.clone(),
            tech: format!("{:?}", cfg.tech),
        },
        summary: SummaryStats {
            events_total: rows.len(),
            events_high_confidence: hc.len(),
        },
        per_kind,
        volcano,
        top_events,
        high_confidence: hc.to_vec(),
        ma_points,
        pca: None,
        umap: None,
        sashimi: Vec::new(),
        isoform_sankey: None,
        advanced_event_counts: BTreeMap::new(),
        diu: Vec::new(),
        heatmap: None,
        coverage_distribution: Vec::new(),
        chrom_distribution,
        junction_graphs: Vec::new(),
        consensus: Vec::new(),
        protein_consequences: Vec::new(),
    }
}

/// Build a top-events × samples PSI heatmap. `top_k` events are chosen by
/// adjusted p-value across all kinds.
pub fn build_heatmap(
    quants: &[crate::quantify::EventQuant],
    rows: &[ResultRow],
    sample_ids: &[String],
    top_k: usize,
) -> Option<HeatmapData> {
    if sample_ids.len() < 2 {
        return None;
    }
    let mut ranked: Vec<&ResultRow> = rows
        .iter()
        .filter(|r| r.adjusted_p_value.is_finite())
        .collect();
    ranked.sort_by(|a, b| {
        a.adjusted_p_value
            .partial_cmp(&b.adjusted_p_value)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    let chosen: std::collections::HashSet<String> = ranked
        .iter()
        .take(top_k)
        .map(|r| r.event_id.clone())
        .collect();
    if chosen.is_empty() {
        return None;
    }
    let mut event_ids: Vec<String> = Vec::new();
    let mut psi_rows: Vec<Vec<f64>> = Vec::new();
    for q in quants {
        if chosen.contains(&q.event.event_id) {
            event_ids.push(q.event.event_id.clone());
            psi_rows.push(q.psi.clone());
        }
    }
    if event_ids.is_empty() {
        return None;
    }
    Some(HeatmapData {
        event_ids,
        samples: sample_ids.to_vec(),
        psi: psi_rows,
    })
}

/// Build a compact splice-junction graph for the top N genes by
/// significance count. Each node = exon, each edge = annotated intron.
pub fn build_junction_graphs(
    ann: &crate::annotation::Annotation,
    rows: &[ResultRow],
    top_n_genes: usize,
) -> Vec<JunctionGraph> {
    let mut counts: BTreeMap<String, usize> = BTreeMap::new();
    for r in rows {
        if r.adjusted_p_value.is_finite() && r.adjusted_p_value <= 0.05 {
            *counts.entry(r.gene_id.clone()).or_insert(0) += 1;
        }
    }
    let mut ranked: Vec<(String, usize)> = counts.into_iter().collect();
    ranked.sort_by(|a, b| b.1.cmp(&a.1));
    let mut out = Vec::new();
    for (gene_id, _) in ranked.into_iter().take(top_n_genes) {
        let Some(g) = ann.genes.get(&gene_id) else {
            continue;
        };
        let mut nodes: Vec<JunctionGraphNode> = Vec::new();
        let mut id_by_exon: std::collections::HashMap<crate::annotation::Exon, usize> =
            std::collections::HashMap::new();
        for (i, (ex, _)) in g.exon_index.iter().enumerate() {
            id_by_exon.insert(*ex, i);
            nodes.push(JunctionGraphNode {
                id: i,
                start: ex.start,
                end: ex.end,
            });
        }
        let mut edges: Vec<JunctionGraphEdge> = Vec::new();
        for edge_idx in g.graph.edge_indices() {
            let (src, tgt) = g.graph.edge_endpoints(edge_idx).unwrap();
            let src_ex = g.graph[src];
            let tgt_ex = g.graph[tgt];
            let w = &g.graph[edge_idx];
            if let (Some(&s), Some(&t)) = (id_by_exon.get(&src_ex), id_by_exon.get(&tgt_ex)) {
                edges.push(JunctionGraphEdge {
                    source: s,
                    target: t,
                    transcripts: w.transcripts.len() as u32,
                });
            }
        }
        out.push(JunctionGraph {
            gene_id: g.gene_id.clone(),
            chrom: g.chrom.clone(),
            nodes,
            edges,
        });
    }
    out
}

fn mean_finite<I: Iterator<Item = f64>>(it: I) -> f64 {
    let mut sum = 0.0;
    let mut n = 0u32;
    for x in it {
        if x.is_finite() {
            sum += x;
            n += 1;
        }
    }
    if n == 0 {
        f64::NAN
    } else {
        sum / n as f64
    }
}

fn histogram<I: Iterator<Item = f64>>(values: I, lo: f64, hi: f64, bins: usize) -> Vec<(f64, usize)> {
    let mut counts = vec![0usize; bins];
    let width = (hi - lo) / bins as f64;
    for v in values {
        if !v.is_finite() {
            continue;
        }
        let mut idx = ((v - lo) / width).floor() as isize;
        if idx < 0 {
            idx = 0;
        }
        if idx as usize >= bins {
            idx = bins as isize - 1;
        }
        counts[idx as usize] += 1;
    }
    counts
        .into_iter()
        .enumerate()
        .map(|(i, c)| (lo + (i as f64 + 0.5) * width, c))
        .collect()
}
