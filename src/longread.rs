//! Long-read isoform reconstruction.
//!
//! For each long read we extract its **junction-chain signature** — the
//! ordered list of `(donor_end, acceptor_start)` splice junctions implied
//! by the read's CIGAR. Reads with identical signatures (after optional
//! wobble collapsing — see [`crate::junctions::collapse_wobble`]) are
//! treated as evidence for the same isoform.
//!
//! Output: an [`IsoformCatalog`] with per-sample read support for each
//! isoform. Single-exon reads (no junctions) are accumulated into a
//! degenerate "unspliced" isoform per gene region.
//!
//! Limitations vs published tools (FLAIR / IsoQuant / TALON):
//! * No correction step against the annotation (e.g. trimming reads to
//!   nearest annotated TSS/PAS). Isoforms are taken at face value from
//!   the read signature.
//! * No quality / completeness scoring beyond raw read support.
//! * No collapsing across genes — a read is assigned by the chromosome
//!   range of its junctions, not by gene_id lookup against the GTF.
//!
//! These are well-defined extension points and the API is designed so
//! adding them is local to this module.

use std::collections::HashMap;

use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use crate::annotation::Annotation;
use crate::config::RunConfig;
use crate::error::{UltiError, UltiResult};
use crate::Position;

/// A reconstructed transcript isoform.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Isoform {
    pub id: String,
    pub gene_id: Option<String>,
    pub chrom: String,
    /// Ordered junction list as `(donor_end, acceptor_start)` tuples.
    pub junctions: Vec<(Position, Position)>,
    /// Per-sample read support, parallel to [`IsoformCatalog::samples`].
    pub support: Vec<f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct IsoformCatalog {
    pub samples: Vec<String>,
    pub isoforms: Vec<Isoform>,
}

impl IsoformCatalog {
    pub fn empty() -> Self {
        IsoformCatalog {
            samples: Vec::new(),
            isoforms: Vec::new(),
        }
    }
}

/// Reconstruct isoforms across the cohort. Returns an empty catalog when
/// the technology is short-read (long reads are required for isoform-level
/// reconstruction).
pub fn reconstruct(cfg: &RunConfig, ann: &Annotation) -> UltiResult<IsoformCatalog> {
    use crate::cli::Technology;
    if cfg.tech == Technology::Short {
        return Ok(IsoformCatalog::empty());
    }

    let tol = cfg.reads.junction_tolerance as u64;

    // **Annotation-anchored correction (FLAIR-style).**
    // Build a per-chromosome sorted list of all annotated junctions. Each
    // read's junctions are snapped to the nearest annotated one within
    // `tol`. Reads whose snapped chain matches an annotated isoform cluster
    // exactly, deterministically — no order-dependent wobble grouping.
    let annotated = AnnotatedJunctionCatalog::from(ann);

    // Per-sample: chain signature -> read count, post-snap.
    let per_sample: Vec<HashMap<Signature, u64>> = cfg
        .samples
        .par_iter()
        .map(|s| {
            read_chains(&s.bam, cfg.reads.min_mapq)
                .map(|raw| snap_chains(&raw, &annotated, tol))
                .unwrap_or_default()
        })
        .collect();

    // Post-snap, equal signatures cluster trivially. We still collapse any
    // residual wobble (reads whose closest annotated junction was just
    // outside `tol`) by the wobble-anchor pass.
    let mut anchors: Vec<Signature> = Vec::new();
    for chains in &per_sample {
        for sig in chains.keys() {
            if !anchors.iter().any(|a| signatures_close(a, sig, tol)) {
                anchors.push(sig.clone());
            }
        }
    }

    let n_samples = cfg.samples.len();
    let mut isoforms = Vec::with_capacity(anchors.len());
    for (i, sig) in anchors.iter().enumerate() {
        let mut support = vec![0.0_f64; n_samples];
        for (s_idx, chains) in per_sample.iter().enumerate() {
            for (other, &count) in chains {
                if signatures_close(sig, other, tol) {
                    support[s_idx] += count as f64;
                }
            }
        }
        let chrom = sig.chrom.clone();
        let gene_id = pick_gene(ann, &chrom, &sig.junctions);
        isoforms.push(Isoform {
            id: format!("ISO_{i:06}"),
            gene_id,
            chrom,
            junctions: sig.junctions.clone(),
            support,
        });
    }

    Ok(IsoformCatalog {
        samples: cfg.samples.iter().map(|s| s.id.clone()).collect(),
        isoforms,
    })
}

/// Annotated splice junctions, sorted per chromosome — fast nearest-neighbor
/// lookup for FLAIR-style snapping.
struct AnnotatedJunctionCatalog {
    /// chrom -> sorted vec of (donor_end, acceptor_start) tuples.
    by_chrom: HashMap<String, Vec<(Position, Position)>>,
}

impl AnnotatedJunctionCatalog {
    fn from(ann: &Annotation) -> Self {
        let mut by_chrom: HashMap<String, Vec<(Position, Position)>> = HashMap::new();
        for g in ann.genes.values() {
            let entry = by_chrom.entry(g.chrom.clone()).or_default();
            for j in g.junctions() {
                entry.push(j);
            }
        }
        for v in by_chrom.values_mut() {
            v.sort_unstable();
            v.dedup();
        }
        AnnotatedJunctionCatalog { by_chrom }
    }

    /// Snap a single junction to the nearest annotated one within `tol`.
    /// Returns `None` when no annotated junction is close enough.
    fn snap(
        &self,
        chrom: &str,
        donor: Position,
        acceptor: Position,
        tol: u64,
    ) -> Option<(Position, Position)> {
        let arr = self.by_chrom.get(chrom)?;
        // Binary search by donor — then linear scan a small window for the
        // best acceptor match. Annotated junctions adjacent in `donor` are
        // also adjacent in this sorted array.
        let idx = match arr.binary_search_by(|j| j.0.cmp(&donor)) {
            Ok(i) => i,
            Err(i) => i,
        };
        let lo = idx.saturating_sub(8);
        let hi = (idx + 8).min(arr.len());
        let mut best: Option<((Position, Position), u64)> = None;
        for &(jd, ja) in &arr[lo..hi] {
            let dd = (jd as i64 - donor as i64).unsigned_abs();
            let da = (ja as i64 - acceptor as i64).unsigned_abs();
            if dd <= tol && da <= tol {
                let score = dd + da;
                if best.as_ref().map(|(_, s)| score < *s).unwrap_or(true) {
                    best = Some(((jd, ja), score));
                }
            }
        }
        best.map(|(j, _)| j)
    }
}

/// Snap every junction in every read's chain to its nearest annotated
/// neighbor (within tolerance). Unsnappable junctions are kept verbatim.
fn snap_chains(
    raw: &HashMap<Signature, u64>,
    catalog: &AnnotatedJunctionCatalog,
    tol: u64,
) -> HashMap<Signature, u64> {
    let mut out: HashMap<Signature, u64> = HashMap::new();
    for (sig, &count) in raw {
        let snapped: Vec<(Position, Position)> = sig
            .junctions
            .iter()
            .map(|&(d, a)| catalog.snap(&sig.chrom, d, a, tol).unwrap_or((d, a)))
            .collect();
        let key = Signature {
            chrom: sig.chrom.clone(),
            junctions: snapped,
        };
        *out.entry(key).or_insert(0) += count;
    }
    out
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
struct Signature {
    chrom: String,
    junctions: Vec<(Position, Position)>,
}

fn signatures_close(a: &Signature, b: &Signature, tol: u64) -> bool {
    if a.chrom != b.chrom || a.junctions.len() != b.junctions.len() {
        return false;
    }
    for (ja, jb) in a.junctions.iter().zip(b.junctions.iter()) {
        let d_donor = (ja.0 as i64 - jb.0 as i64).unsigned_abs();
        let d_acc = (ja.1 as i64 - jb.1 as i64).unsigned_abs();
        if d_donor > tol || d_acc > tol {
            return false;
        }
    }
    true
}

/// Read BAM and extract one junction-chain per primary aligned read.
fn read_chains(path: &std::path::Path, min_mapq: u8) -> UltiResult<HashMap<Signature, u64>> {
    use noodles::bam;
    use noodles::sam::alignment::record::cigar::op::Kind;

    let mut reader = bam::io::reader::Builder
        .build_from_path(path)
        .map_err(|e| UltiError::Alignment {
            path: path.into(),
            message: format!("cannot open BAM: {e}"),
        })?;
    let header = reader.read_header().map_err(|e| UltiError::Alignment {
        path: path.into(),
        message: format!("cannot read header: {e}"),
    })?;

    let mut counts: HashMap<Signature, u64> = HashMap::new();
    for result in reader.records() {
        let Ok(record) = result else { continue };
        let flags = record.flags();
        if flags.is_unmapped() || flags.is_secondary() || flags.is_supplementary() {
            continue;
        }
        let mapq = record.mapping_quality().map(|q| q.get()).unwrap_or(0);
        if mapq < min_mapq {
            continue;
        }
        let ref_id = match record.reference_sequence_id() {
            Some(Ok(id)) => id,
            _ => continue,
        };
        let chrom = match header.reference_sequences().get_index(ref_id) {
            Some((name, _)) => name.to_string(),
            None => continue,
        };
        let mut pos = match record.alignment_start() {
            Some(Ok(p)) => p.get() as u64,
            _ => continue,
        };
        let mut junctions: Vec<(Position, Position)> = Vec::new();
        for op_result in record.cigar().iter() {
            let Ok(op) = op_result else {
                junctions.clear();
                break;
            };
            let len = op.len() as u64;
            match op.kind() {
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch | Kind::Deletion => {
                    pos += len;
                }
                Kind::Skip => {
                    let donor_end = pos.saturating_sub(1);
                    let acceptor_start = pos + len;
                    junctions.push((donor_end, acceptor_start));
                    pos += len;
                }
                _ => {}
            }
        }
        let sig = Signature {
            chrom,
            junctions,
        };
        *counts.entry(sig).or_insert(0) += 1;
    }
    Ok(counts)
}

/// Heuristic gene assignment: find a gene whose splice graph contains any
/// of the read's junctions. First match wins. Returns `None` for
/// completely intergenic / unannotated reads.
fn pick_gene(ann: &Annotation, chrom: &str, junctions: &[(Position, Position)]) -> Option<String> {
    if junctions.is_empty() {
        return None;
    }
    for gene_id in ann.gene_by_chrom.get(chrom)? {
        let g = ann.genes.get(gene_id)?;
        for j in g.junctions() {
            if junctions.iter().any(|x| x.0 == j.0 && x.1 == j.1) {
                return Some(gene_id.clone());
            }
        }
    }
    None
}

/// Mean read depth in the intron interior (a thin wrapper around
/// [`crate::pileup::pileup_regions`]).
pub fn intron_coverage(
    bam: &std::path::Path,
    chrom: &str,
    intron: (Position, Position),
) -> UltiResult<f64> {
    use crate::pileup::{pileup_regions, Region};
    let result = pileup_regions(
        bam,
        &[Region {
            chrom: chrom.to_string(),
            start: intron.0,
            end: intron.1,
        }],
        0,
    )?;
    Ok(result.mean_depth.first().copied().unwrap_or(0.0))
}

/// Differential isoform usage test. Compares per-sample isoform fractions
/// between the contrast's two groups using a per-gene chi-square test on
/// the rows-by-isoforms × cols-by-groups table.
///
/// Returns one record per gene that has at least 2 isoforms with non-zero
/// support across the cohort.
pub fn differential_usage(
    catalog: &IsoformCatalog,
    sample_groups: &[&str],
    numerator: &str,
    denominator: &str,
) -> Vec<DiuRecord> {
    // Bucket isoforms by gene.
    let mut by_gene: HashMap<String, Vec<&Isoform>> = HashMap::new();
    for iso in &catalog.isoforms {
        if let Some(g) = &iso.gene_id {
            by_gene.entry(g.clone()).or_default().push(iso);
        }
    }

    let mut out = Vec::new();
    for (gene, isos) in by_gene {
        if isos.len() < 2 {
            continue;
        }
        // Build per-group totals per isoform.
        let mut num_counts = vec![0.0_f64; isos.len()];
        let mut denom_counts = vec![0.0_f64; isos.len()];
        for (s_idx, g) in sample_groups.iter().enumerate() {
            for (iso_idx, iso) in isos.iter().enumerate() {
                let c = iso.support.get(s_idx).copied().unwrap_or(0.0);
                if *g == numerator {
                    num_counts[iso_idx] += c;
                } else if *g == denominator {
                    denom_counts[iso_idx] += c;
                }
            }
        }
        let total_num: f64 = num_counts.iter().sum();
        let total_denom: f64 = denom_counts.iter().sum();
        if total_num < 10.0 || total_denom < 10.0 {
            continue;
        }

        // Chi-square goodness of fit on the K × 2 table.
        let mut chi2 = 0.0_f64;
        let total = total_num + total_denom;
        for i in 0..isos.len() {
            let row_total = num_counts[i] + denom_counts[i];
            if row_total <= 0.0 {
                continue;
            }
            let e_num = row_total * total_num / total;
            let e_denom = row_total * total_denom / total;
            if e_num > 0.0 {
                chi2 += (num_counts[i] - e_num).powi(2) / e_num;
            }
            if e_denom > 0.0 {
                chi2 += (denom_counts[i] - e_denom).powi(2) / e_denom;
            }
        }
        let df = (isos.len() - 1).max(1) as f64;
        let p = chi_square_sf(chi2, df);

        // Per-isoform usage fractions.
        let mut iso_records = Vec::new();
        for i in 0..isos.len() {
            iso_records.push(IsoformUsage {
                isoform_id: isos[i].id.clone(),
                num_count: num_counts[i],
                denom_count: denom_counts[i],
                num_fraction: if total_num > 0.0 { num_counts[i] / total_num } else { 0.0 },
                denom_fraction: if total_denom > 0.0 {
                    denom_counts[i] / total_denom
                } else {
                    0.0
                },
            });
        }
        out.push(DiuRecord {
            gene_id: gene,
            n_isoforms: isos.len(),
            chi2,
            df,
            p_value: p,
            isoforms: iso_records,
        });
    }
    out
}

/// Differential-isoform-usage record per gene.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DiuRecord {
    pub gene_id: String,
    pub n_isoforms: usize,
    pub chi2: f64,
    pub df: f64,
    pub p_value: f64,
    pub isoforms: Vec<IsoformUsage>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct IsoformUsage {
    pub isoform_id: String,
    pub num_count: f64,
    pub denom_count: f64,
    pub num_fraction: f64,
    pub denom_fraction: f64,
}

fn chi_square_sf(stat: f64, df: f64) -> f64 {
    use statrs::distribution::{ChiSquared, ContinuousCDF};
    if !stat.is_finite() || stat <= 0.0 || df <= 0.0 {
        return 1.0;
    }
    match ChiSquared::new(df) {
        Ok(d) => (1.0 - d.cdf(stat)).clamp(0.0, 1.0),
        Err(_) => 1.0,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn signatures_collapse_with_wobble() {
        let a = Signature {
            chrom: "chr1".into(),
            junctions: vec![(100, 200), (300, 400)],
        };
        let b = Signature {
            chrom: "chr1".into(),
            junctions: vec![(102, 199), (301, 402)],
        };
        let c = Signature {
            chrom: "chr1".into(),
            junctions: vec![(150, 250)],
        };
        assert!(signatures_close(&a, &b, 5));
        assert!(!signatures_close(&a, &b, 1));
        assert!(!signatures_close(&a, &c, 100));
    }

    #[test]
    fn snap_chain_to_annotated() {
        // Annotated junctions: (100, 200), (300, 400).
        let mut by_chrom = std::collections::HashMap::new();
        by_chrom.insert("chr1".to_string(), vec![(100, 200), (300, 400)]);
        let catalog = AnnotatedJunctionCatalog { by_chrom };

        // Read with slightly off coords — should snap.
        let mut raw = std::collections::HashMap::new();
        raw.insert(
            Signature {
                chrom: "chr1".into(),
                junctions: vec![(102, 199), (298, 401)],
            },
            5_u64,
        );
        // Read whose junction is too far from any annotated — should not snap.
        raw.insert(
            Signature {
                chrom: "chr1".into(),
                junctions: vec![(150, 250)],
            },
            3_u64,
        );

        let snapped = snap_chains(&raw, &catalog, 5);
        // Snapped read collapses to the canonical (100,200),(300,400).
        let canonical = Signature {
            chrom: "chr1".into(),
            junctions: vec![(100, 200), (300, 400)],
        };
        assert_eq!(snapped.get(&canonical), Some(&5));
        // Non-snappable read is preserved verbatim.
        let unsnappable = Signature {
            chrom: "chr1".into(),
            junctions: vec![(150, 250)],
        };
        assert_eq!(snapped.get(&unsnappable), Some(&3));
    }

    #[test]
    fn diu_detects_clear_switch() {
        let catalog = IsoformCatalog {
            samples: vec!["s1".into(), "s2".into(), "s3".into(), "s4".into()],
            isoforms: vec![
                Isoform {
                    id: "iso_a".into(),
                    gene_id: Some("G1".into()),
                    chrom: "chr1".into(),
                    junctions: vec![(100, 200)],
                    support: vec![100.0, 110.0, 5.0, 8.0],
                },
                Isoform {
                    id: "iso_b".into(),
                    gene_id: Some("G1".into()),
                    chrom: "chr1".into(),
                    junctions: vec![(100, 300)],
                    support: vec![10.0, 8.0, 95.0, 100.0],
                },
            ],
        };
        let groups: Vec<&str> = vec!["ctrl", "ctrl", "treat", "treat"];
        let recs = differential_usage(&catalog, &groups, "treat", "ctrl");
        assert_eq!(recs.len(), 1);
        assert!(recs[0].p_value < 1e-10, "p = {}", recs[0].p_value);
    }
}
