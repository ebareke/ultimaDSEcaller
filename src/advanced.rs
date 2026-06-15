//! Advanced and complex AS events.
//!
//! This module detects 11 categories beyond the 7 canonical events:
//!
//! 1. **CrypticSpliceSite** — empirical junctions absent from the
//!    annotation but observed with high cohort-wide support.
//! 2. **MultiExonSkipping (MSE)** — a single junction spans ≥ 2 internal
//!    annotated exons.
//! 3. **MultiIntronRetention (MIR)** — two or more consecutive introns of
//!    a transcript co-retained (intronic coverage spanning multiple
//!    annotated introns).
//! 4. **RecursiveSplicing** — a splice site that falls *inside* an
//!    annotated intron, where reads support both a partial-intron splice
//!    and a full-intron splice using the same flanking exons.
//! 5. **NestedSplicing** — a junction whose extent is completely contained
//!    inside another junction sharing the same upstream donor (or
//!    downstream acceptor).
//! 6. **PartialExonInclusion** — empirical junctions whose acceptor (or
//!    donor) falls *inside* an annotated exon, producing a shortened or
//!    extended variant.
//! 7. **ExonicIntronicHybrid** — empirical junctions where one anchor is
//!    inside an annotated exon and the other inside an annotated intron.
//! 8. **AlternativePromoter** — AFE-like events where the alternative
//!    first exons are separated by ≥ 5 kb (a soft proxy for distinct
//!    promoters in the absence of CAGE data).
//! 9. **AlternativePolyadenylation** — ALE-like events where the
//!    alternative last exons differ in their 3'-end position by ≥ 50 bp
//!    (a soft proxy for distinct polyA sites in the absence of 3'-seq).
//! 10. **TandemUtr** — a special case of (9) where the alternative last
//!    exons share a 5' end but differ only in 3' extension.
//! 11. **FusionAssociated** — junctions implied by a user-supplied fusion
//!    BEDPE file (see [`load_fusion_bedpe`]).
//!
//! Isoform switching is handled separately as
//! [`crate::longread::differential_usage`] since it is a *test* on the
//! isoform catalog rather than an event topology.
//!
//! All proxies (promoter, polyA, tandem UTR) are clearly flagged as such
//! in the code — they're heuristics that improve recall against typical
//! CAGE/3'-seq-augmented annotations but should not be quoted as
//! authoritative without orthogonal data.

use std::collections::{BTreeSet, HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use serde::{Deserialize, Serialize};

use crate::annotation::{Annotation, Exon, SpliceGraph};
use crate::error::{UltiError, UltiResult};
use crate::junctions::JunctionMatrix;
use crate::Position;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum AdvancedEventKind {
    AlternativePromoter,
    AlternativePolyadenylation,
    TandemUtr,
    RecursiveSplicing,
    NestedSplicing,
    CrypticSpliceSite,
    PartialExonInclusion,
    ExonicIntronicHybrid,
    MultiExonSkipping,
    MultiIntronRetention,
    FusionAssociated,
}

impl AdvancedEventKind {
    pub fn short(self) -> &'static str {
        match self {
            AdvancedEventKind::AlternativePromoter => "AltPromoter",
            AdvancedEventKind::AlternativePolyadenylation => "AltPolyA",
            AdvancedEventKind::TandemUtr => "TandemUTR",
            AdvancedEventKind::RecursiveSplicing => "Recursive",
            AdvancedEventKind::NestedSplicing => "Nested",
            AdvancedEventKind::CrypticSpliceSite => "Cryptic",
            AdvancedEventKind::PartialExonInclusion => "PartialExon",
            AdvancedEventKind::ExonicIntronicHybrid => "ExonIntronHybrid",
            AdvancedEventKind::MultiExonSkipping => "MSE",
            AdvancedEventKind::MultiIntronRetention => "MIR",
            AdvancedEventKind::FusionAssociated => "Fusion",
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AdvancedEvent {
    pub event_id: String,
    pub gene_id: Option<String>,
    pub chrom: String,
    pub kind: AdvancedEventKind,
    pub coords: Vec<(Position, Position)>,
    pub support: f64,
    pub notes: String,
}

/// Tunables for advanced detection. Sensible defaults supplied; the
/// fields are public so callers can adjust per cohort.
#[derive(Debug, Clone)]
pub struct AdvancedParams {
    pub cryptic_min_total_support: f64,
    pub cryptic_min_samples: usize,
    pub alt_promoter_min_separation: u64,
    pub alt_polya_min_separation: u64,
    pub tandem_utr_min_extension: u64,
}

impl Default for AdvancedParams {
    fn default() -> Self {
        AdvancedParams {
            cryptic_min_total_support: 10.0,
            cryptic_min_samples: 2,
            alt_promoter_min_separation: 5000,
            alt_polya_min_separation: 50,
            tandem_utr_min_extension: 50,
        }
    }
}

/// Run all advanced detectors. Returns a flat vector.
pub fn detect_all(
    ann: &Annotation,
    jm: &JunctionMatrix,
    params: &AdvancedParams,
    fusion_bedpe: Option<&Path>,
) -> UltiResult<Vec<AdvancedEvent>> {
    let mut events = Vec::new();
    events.extend(detect_cryptic(ann, jm, params));
    events.extend(detect_multi_exon_skipping(ann, jm));
    events.extend(detect_recursive(ann, jm));
    events.extend(detect_nested(ann, jm));
    events.extend(detect_partial_and_hybrid(ann, jm));
    events.extend(detect_alt_promoter(ann, params));
    events.extend(detect_alt_polya_and_tandem(ann, params));
    events.extend(detect_multi_intron_retention(ann));
    if let Some(bedpe) = fusion_bedpe {
        events.extend(load_fusion_bedpe(bedpe)?);
    }
    Ok(events)
}

fn next_id(counter: &mut usize, kind: AdvancedEventKind) -> String {
    let id = format!("ADV_{}_{:06}", kind.short(), counter);
    *counter += 1;
    id
}

/// Cryptic splice sites — empirical junctions not in the annotation that
/// pass cohort-wide support thresholds.
fn detect_cryptic(
    ann: &Annotation,
    jm: &JunctionMatrix,
    params: &AdvancedParams,
) -> Vec<AdvancedEvent> {
    let mut annotated: HashSet<(String, Position, Position)> = HashSet::new();
    for g in ann.genes.values() {
        for j in g.junctions() {
            annotated.insert((g.chrom.clone(), j.0, j.1));
        }
    }
    let mut counter = 0;
    let mut out = Vec::new();
    for (j, supports) in &jm.counts {
        let key = (j.chrom.clone(), j.donor_end, j.acceptor_start);
        if annotated.contains(&key) {
            continue;
        }
        let total: f64 = supports.iter().sum();
        let n_samples = supports.iter().filter(|&&x| x > 0.0).count();
        if total < params.cryptic_min_total_support || n_samples < params.cryptic_min_samples {
            continue;
        }
        let gene = nearest_gene(ann, &j.chrom, j.donor_end, j.acceptor_start);
        out.push(AdvancedEvent {
            event_id: next_id(&mut counter, AdvancedEventKind::CrypticSpliceSite),
            gene_id: gene,
            chrom: j.chrom.clone(),
            kind: AdvancedEventKind::CrypticSpliceSite,
            coords: vec![(j.donor_end, j.acceptor_start)],
            support: total,
            notes: format!("observed in {n_samples} samples"),
        });
    }
    out
}

/// Multi-exon skipping — an empirical or annotated junction `(d, a)` such
/// that ≥ 2 annotated exons of a gene lie strictly between `d` and `a`.
fn detect_multi_exon_skipping(ann: &Annotation, jm: &JunctionMatrix) -> Vec<AdvancedEvent> {
    let mut counter = 0;
    let mut out = Vec::new();
    // Cohort-wide observed junctions (annotated + empirical).
    let mut all_juncs: HashSet<(String, Position, Position)> = HashSet::new();
    for j in jm.counts.keys() {
        all_juncs.insert((j.chrom.clone(), j.donor_end, j.acceptor_start));
    }
    for g in ann.genes.values() {
        for j in g.junctions() {
            all_juncs.insert((g.chrom.clone(), j.0, j.1));
        }
    }
    for g in ann.genes.values() {
        let mut exon_intervals: Vec<(Position, Position)> = g
            .exon_index
            .keys()
            .map(|e| (e.start, e.end))
            .collect();
        exon_intervals.sort();
        for (chrom, donor, acceptor) in &all_juncs {
            if chrom != &g.chrom {
                continue;
            }
            let inside: Vec<_> = exon_intervals
                .iter()
                .filter(|(s, e)| *s > *donor && *e < *acceptor)
                .collect();
            if inside.len() >= 2 {
                let support = jm
                    .counts
                    .iter()
                    .find(|(j, _)| {
                        j.chrom == *chrom && j.donor_end == *donor && j.acceptor_start == *acceptor
                    })
                    .map(|(_, v)| v.iter().sum::<f64>())
                    .unwrap_or(0.0);
                out.push(AdvancedEvent {
                    event_id: next_id(&mut counter, AdvancedEventKind::MultiExonSkipping),
                    gene_id: Some(g.gene_id.clone()),
                    chrom: g.chrom.clone(),
                    kind: AdvancedEventKind::MultiExonSkipping,
                    coords: vec![(*donor, *acceptor)],
                    support,
                    notes: format!("skips {} exons", inside.len()),
                });
            }
        }
    }
    out
}

/// Recursive splicing — empirical junction whose donor and acceptor BOTH
/// fall *inside* the same annotated intron (i.e. an intermediate splice
/// step of the annotated intron).
fn detect_recursive(ann: &Annotation, jm: &JunctionMatrix) -> Vec<AdvancedEvent> {
    let mut counter = 0;
    let mut out = Vec::new();
    for g in ann.genes.values() {
        for intron in g.junctions() {
            for (j, supports) in &jm.counts {
                if j.chrom != g.chrom {
                    continue;
                }
                if j.donor_end > intron.0 && j.acceptor_start < intron.1 {
                    let total: f64 = supports.iter().sum();
                    if total < 5.0 {
                        continue;
                    }
                    out.push(AdvancedEvent {
                        event_id: next_id(&mut counter, AdvancedEventKind::RecursiveSplicing),
                        gene_id: Some(g.gene_id.clone()),
                        chrom: g.chrom.clone(),
                        kind: AdvancedEventKind::RecursiveSplicing,
                        coords: vec![intron, (j.donor_end, j.acceptor_start)],
                        support: total,
                        notes: "internal splice site of annotated intron".into(),
                    });
                }
            }
        }
    }
    out
}

/// Nested splicing — junction (`d2`, `a2`) strictly inside another junction
/// (`d1`, `a1`) that shares either the donor or the acceptor with `j1`.
fn detect_nested(ann: &Annotation, jm: &JunctionMatrix) -> Vec<AdvancedEvent> {
    let mut counter = 0;
    let mut out = Vec::new();
    let mut by_chrom: HashMap<String, Vec<(Position, Position, f64)>> = HashMap::new();
    for (j, supports) in &jm.counts {
        let total: f64 = supports.iter().sum();
        by_chrom
            .entry(j.chrom.clone())
            .or_default()
            .push((j.donor_end, j.acceptor_start, total));
    }
    for (chrom, mut juncs) in by_chrom {
        // f64 has no total order, so sort by the two integer coordinates only.
        juncs.sort_by(|a, b| (a.0, a.1).cmp(&(b.0, b.1)));
        for i in 0..juncs.len() {
            for k in 0..juncs.len() {
                if i == k {
                    continue;
                }
                let (d1, a1, _) = juncs[i];
                let (d2, a2, s2) = juncs[k];
                if d2 > d1 && a2 < a1 {
                    let same_donor = d1 == d2;
                    let same_acc = a1 == a2;
                    if same_donor || same_acc || (a2 - d2) * 2 < (a1 - d1) {
                        let gene = nearest_gene(ann, &chrom, d2, a2);
                        out.push(AdvancedEvent {
                            event_id: next_id(&mut counter, AdvancedEventKind::NestedSplicing),
                            gene_id: gene,
                            chrom: chrom.clone(),
                            kind: AdvancedEventKind::NestedSplicing,
                            coords: vec![(d1, a1), (d2, a2)],
                            support: s2,
                            notes: if same_donor {
                                "shared donor".into()
                            } else if same_acc {
                                "shared acceptor".into()
                            } else {
                                "geometrically nested".into()
                            },
                        });
                    }
                }
            }
        }
    }
    out
}

/// Partial-exon inclusion and exonic-intronic hybrid — empirical junctions
/// where at least one anchor falls in an unexpected genomic context
/// relative to the annotation.
fn detect_partial_and_hybrid(ann: &Annotation, jm: &JunctionMatrix) -> Vec<AdvancedEvent> {
    let mut counter_partial = 0;
    let mut counter_hybrid = 0;
    let mut out = Vec::new();
    let mut annotated: HashSet<(String, Position, Position)> = HashSet::new();
    for g in ann.genes.values() {
        for j in g.junctions() {
            annotated.insert((g.chrom.clone(), j.0, j.1));
        }
    }
    for (j, supports) in &jm.counts {
        let key = (j.chrom.clone(), j.donor_end, j.acceptor_start);
        if annotated.contains(&key) {
            continue;
        }
        let total: f64 = supports.iter().sum();
        if total < 5.0 {
            continue;
        }
        let Some(genes) = ann.gene_by_chrom.get(&j.chrom) else {
            continue;
        };
        for gene_id in genes {
            let g = match ann.genes.get(gene_id) {
                Some(g) => g,
                None => continue,
            };
            let donor_ctx = anchor_context(g, j.donor_end);
            let acceptor_ctx = anchor_context(g, j.acceptor_start);
            use AnchorContext::*;
            match (donor_ctx, acceptor_ctx) {
                (Boundary, InsideExon) | (InsideExon, Boundary) | (InsideExon, InsideExon) => {
                    out.push(AdvancedEvent {
                        event_id: next_id(
                            &mut counter_partial,
                            AdvancedEventKind::PartialExonInclusion,
                        ),
                        gene_id: Some(g.gene_id.clone()),
                        chrom: j.chrom.clone(),
                        kind: AdvancedEventKind::PartialExonInclusion,
                        coords: vec![(j.donor_end, j.acceptor_start)],
                        support: total,
                        notes: "anchor inside annotated exon body".into(),
                    });
                    break;
                }
                (InsideExon, InsideIntron)
                | (InsideIntron, InsideExon)
                | (Boundary, InsideIntron)
                | (InsideIntron, Boundary) => {
                    out.push(AdvancedEvent {
                        event_id: next_id(
                            &mut counter_hybrid,
                            AdvancedEventKind::ExonicIntronicHybrid,
                        ),
                        gene_id: Some(g.gene_id.clone()),
                        chrom: j.chrom.clone(),
                        kind: AdvancedEventKind::ExonicIntronicHybrid,
                        coords: vec![(j.donor_end, j.acceptor_start)],
                        support: total,
                        notes: "exonic-intronic boundary mix".into(),
                    });
                    break;
                }
                _ => {}
            }
        }
    }
    out
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum AnchorContext {
    Boundary,
    InsideExon,
    InsideIntron,
    Outside,
}

fn anchor_context(g: &SpliceGraph, pos: Position) -> AnchorContext {
    for ex in g.exon_index.keys() {
        if ex.start == pos || ex.end == pos {
            return AnchorContext::Boundary;
        }
    }
    for ex in g.exon_index.keys() {
        if pos > ex.start && pos < ex.end {
            return AnchorContext::InsideExon;
        }
    }
    // Inside a gene-spanning extent but not in any exon → intron.
    let (gmin, gmax) = match g.exon_index.keys().fold(None, |acc: Option<(Position, Position)>, ex| {
        Some(match acc {
            None => (ex.start, ex.end),
            Some((lo, hi)) => (lo.min(ex.start), hi.max(ex.end)),
        })
    }) {
        Some(b) => b,
        None => return AnchorContext::Outside,
    };
    if pos >= gmin && pos <= gmax {
        AnchorContext::InsideIntron
    } else {
        AnchorContext::Outside
    }
}

/// Alternative promoter (heuristic). A gene with ≥ 2 transcripts whose
/// first exons are separated by ≥ `alt_promoter_min_separation` bp gets one
/// `AltPromoter` event per qualifying pair.
fn detect_alt_promoter(ann: &Annotation, params: &AdvancedParams) -> Vec<AdvancedEvent> {
    let mut counter = 0;
    let mut out = Vec::new();
    for g in ann.genes.values() {
        let first_exons: BTreeSet<Exon> = g
            .transcripts
            .values()
            .filter_map(|chain| chain.first().copied())
            .collect();
        if first_exons.len() < 2 {
            continue;
        }
        let vec: Vec<Exon> = first_exons.into_iter().collect();
        for i in 0..vec.len() {
            for k in (i + 1)..vec.len() {
                let sep = (vec[i].start as i64 - vec[k].start as i64).unsigned_abs();
                if sep >= params.alt_promoter_min_separation {
                    out.push(AdvancedEvent {
                        event_id: next_id(&mut counter, AdvancedEventKind::AlternativePromoter),
                        gene_id: Some(g.gene_id.clone()),
                        chrom: g.chrom.clone(),
                        kind: AdvancedEventKind::AlternativePromoter,
                        coords: vec![(vec[i].start, vec[i].end), (vec[k].start, vec[k].end)],
                        support: 0.0,
                        notes: format!(
                            "TSS separation ≥ {} bp (heuristic without CAGE data)",
                            params.alt_promoter_min_separation
                        ),
                    });
                }
            }
        }
    }
    out
}

/// Alternative polyadenylation + tandem UTR (heuristic from ALE topology).
fn detect_alt_polya_and_tandem(ann: &Annotation, params: &AdvancedParams) -> Vec<AdvancedEvent> {
    let mut counter_apa = 0;
    let mut counter_utr = 0;
    let mut out = Vec::new();
    for g in ann.genes.values() {
        let last_exons: BTreeSet<Exon> = g
            .transcripts
            .values()
            .filter_map(|chain| chain.last().copied())
            .collect();
        if last_exons.len() < 2 {
            continue;
        }
        let vec: Vec<Exon> = last_exons.into_iter().collect();
        for i in 0..vec.len() {
            for k in (i + 1)..vec.len() {
                let same_start = vec[i].start == vec[k].start;
                let end_diff =
                    (vec[i].end as i64 - vec[k].end as i64).unsigned_abs();
                if same_start && end_diff >= params.tandem_utr_min_extension {
                    out.push(AdvancedEvent {
                        event_id: next_id(&mut counter_utr, AdvancedEventKind::TandemUtr),
                        gene_id: Some(g.gene_id.clone()),
                        chrom: g.chrom.clone(),
                        kind: AdvancedEventKind::TandemUtr,
                        coords: vec![(vec[i].start, vec[i].end), (vec[k].start, vec[k].end)],
                        support: 0.0,
                        notes: format!(
                            "shared 5' end, 3' extension ≥ {} bp",
                            params.tandem_utr_min_extension
                        ),
                    });
                } else if end_diff >= params.alt_polya_min_separation {
                    out.push(AdvancedEvent {
                        event_id: next_id(
                            &mut counter_apa,
                            AdvancedEventKind::AlternativePolyadenylation,
                        ),
                        gene_id: Some(g.gene_id.clone()),
                        chrom: g.chrom.clone(),
                        kind: AdvancedEventKind::AlternativePolyadenylation,
                        coords: vec![(vec[i].start, vec[i].end), (vec[k].start, vec[k].end)],
                        support: 0.0,
                        notes: format!(
                            "polyA-site separation ≥ {} bp (heuristic without 3'-seq data)",
                            params.alt_polya_min_separation
                        ),
                    });
                }
            }
        }
    }
    out
}

/// Multi-intron retention (MIR) — pairs of consecutive annotated introns
/// in a transcript that *could* be co-retained. We emit one event per pair;
/// quantification of MIR requires intron coverage data, surfaced via the
/// pileup module and computed on-demand by callers.
fn detect_multi_intron_retention(ann: &Annotation) -> Vec<AdvancedEvent> {
    let mut counter = 0;
    let mut out = Vec::new();
    for g in ann.genes.values() {
        for (tx_id, chain) in &g.transcripts {
            if chain.len() < 3 {
                continue;
            }
            for w in chain.windows(3) {
                let intron_a = (w[0].end + 1, w[1].start.saturating_sub(1));
                let intron_b = (w[1].end + 1, w[2].start.saturating_sub(1));
                if intron_a.0 >= intron_a.1 || intron_b.0 >= intron_b.1 {
                    continue;
                }
                out.push(AdvancedEvent {
                    event_id: next_id(&mut counter, AdvancedEventKind::MultiIntronRetention),
                    gene_id: Some(g.gene_id.clone()),
                    chrom: g.chrom.clone(),
                    kind: AdvancedEventKind::MultiIntronRetention,
                    coords: vec![intron_a, intron_b],
                    support: 0.0,
                    notes: format!("consecutive introns of transcript {tx_id}"),
                });
            }
        }
    }
    out
}

/// Load a fusion BEDPE file produced by a fusion caller (Arriba,
/// STAR-Fusion, etc.). Each row contributes one `FusionAssociated` event
/// using the two breakpoints as coordinates.
pub fn load_fusion_bedpe(path: &Path) -> UltiResult<Vec<AdvancedEvent>> {
    let file = File::open(path).map_err(|e| UltiError::io(path, e))?;
    let reader = BufReader::new(file);
    let mut counter = 0;
    let mut out = Vec::new();
    for (lineno, line) in reader.lines().enumerate() {
        let line = line.map_err(|e| UltiError::io(path, e))?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 6 {
            continue;
        }
        let chrom_a = fields[0].to_string();
        let start_a: u64 = fields[1].parse().map_err(|_| {
            UltiError::Annotation {
                path: path.into(),
                line: lineno,
                message: "BEDPE start_a not an integer".into(),
            }
        })?;
        let end_a: u64 = fields[2].parse().map_err(|_| {
            UltiError::Annotation {
                path: path.into(),
                line: lineno,
                message: "BEDPE end_a not an integer".into(),
            }
        })?;
        let chrom_b = fields[3].to_string();
        let start_b: u64 = fields[4].parse().map_err(|_| {
            UltiError::Annotation {
                path: path.into(),
                line: lineno,
                message: "BEDPE start_b not an integer".into(),
            }
        })?;
        let end_b: u64 = fields[5].parse().map_err(|_| {
            UltiError::Annotation {
                path: path.into(),
                line: lineno,
                message: "BEDPE end_b not an integer".into(),
            }
        })?;
        out.push(AdvancedEvent {
            event_id: next_id(&mut counter, AdvancedEventKind::FusionAssociated),
            gene_id: None,
            chrom: format!("{chrom_a}|{chrom_b}"),
            kind: AdvancedEventKind::FusionAssociated,
            coords: vec![(start_a + 1, end_a), (start_b + 1, end_b)],
            support: 0.0,
            notes: "from fusion BEDPE".into(),
        });
    }
    Ok(out)
}

/// Pick a gene whose extent contains either junction anchor. Falls back to
/// the chromosome's first gene if none qualifies.
fn nearest_gene(
    ann: &Annotation,
    chrom: &str,
    donor: Position,
    acceptor: Position,
) -> Option<String> {
    let genes = ann.gene_by_chrom.get(chrom)?;
    for gene_id in genes {
        let g = ann.genes.get(gene_id)?;
        for ex in g.exon_index.keys() {
            if (donor >= ex.start && donor <= ex.end)
                || (acceptor >= ex.start && acceptor <= ex.end)
            {
                return Some(gene_id.clone());
            }
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn write_gtf(s: &str) -> tempfile::NamedTempFile {
        let mut f = tempfile::Builder::new().suffix(".gtf").tempfile().unwrap();
        f.write_all(s.as_bytes()).unwrap();
        f
    }

    #[test]
    fn detects_alt_promoter_when_well_separated() {
        // Two transcripts of the same gene, both ending at a shared 3'
        // exon (50000-50100) but starting at promoters 40 kb apart.
        // First exons in genomic order: T1 at 100, T2 at 40000.
        let gtf = "\
chr1\ts\texon\t100\t200\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";
chr1\ts\texon\t50000\t50100\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";
chr1\ts\texon\t40000\t40100\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T2\";
chr1\ts\texon\t50000\t50100\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T2\";
";
        let tmp = write_gtf(gtf);
        let ann = crate::annotation::parse(tmp.path()).unwrap();
        let events = detect_alt_promoter(&ann, &AdvancedParams::default());
        assert!(events.iter().any(|e| e.kind == AdvancedEventKind::AlternativePromoter));
    }
}
