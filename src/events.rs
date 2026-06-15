//! Event detection — enumerates the 7 canonical alternative-splicing events
//! from each gene's splice graph.
//!
//! For every event we record:
//! * Stable coordinates (used to produce a deterministic `event_id`)
//! * The set of *inclusion* junctions — reads supporting the alternative form
//! * The set of *exclusion* junctions — reads supporting the constitutive form
//!
//! PSI for an event then becomes simply:
//!
//! ```text
//! psi = inclusion_reads / (inclusion_reads + exclusion_reads)
//! ```
//!
//! Intron retention (IR) is the one structural exception: the "retained"
//! form contributes no junction, only intronic coverage. The IR coords are
//! emitted so the quantification stage can compute coverage from a pileup
//! (currently TODO — see [`crate::quantify`]).

use std::collections::HashSet;

use petgraph::visit::EdgeRef;
use serde::{Deserialize, Serialize};

use crate::annotation::{Annotation, Exon, SpliceGraph};
use crate::junctions::JunctionMatrix;
use crate::{Position, Strand};

/// Canonical AS event categories.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum EventKind {
    SE,
    MXE,
    A5SS,
    A3SS,
    IR,
    AFE,
    ALE,
}

impl EventKind {
    pub fn short(&self) -> &'static str {
        match self {
            EventKind::SE => "SE",
            EventKind::MXE => "MXE",
            EventKind::A5SS => "A5SS",
            EventKind::A3SS => "A3SS",
            EventKind::IR => "IR",
            EventKind::AFE => "AFE",
            EventKind::ALE => "ALE",
        }
    }
}

/// One enumerated AS event. The `event_id` is deterministic from coordinates
/// so the same event has the same ID across runs and across samples.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ASEvent {
    pub event_id: String,
    pub gene_id: String,
    pub chrom: String,
    pub strand: Strand,
    pub kind: EventKind,
    /// Exons involved, in canonical order (see per-kind comments below).
    pub exons: Vec<Exon>,
    /// Junctions whose presence indicates the inclusion / alternative form.
    pub inclusion_junctions: Vec<(Position, Position)>,
    /// Junctions whose presence indicates the exclusion / constitutive form.
    pub exclusion_junctions: Vec<(Position, Position)>,
    /// For IR events: the genomic interval of the intron whose retention is
    /// being measured. `None` for all other kinds.
    pub retained_intron: Option<(Position, Position)>,
}

impl ASEvent {
    fn new(
        gene_id: &str,
        chrom: &str,
        strand: Strand,
        kind: EventKind,
        exons: Vec<Exon>,
        inclusion_junctions: Vec<(Position, Position)>,
        exclusion_junctions: Vec<(Position, Position)>,
        retained_intron: Option<(Position, Position)>,
    ) -> Self {
        let id = format!(
            "{gene_id}|{kind}|{chrom}:{coords}",
            kind = kind.short(),
            coords = exons
                .iter()
                .map(|e| format!("{}-{}", e.start, e.end))
                .collect::<Vec<_>>()
                .join("_"),
        );
        ASEvent {
            event_id: id,
            gene_id: gene_id.to_string(),
            chrom: chrom.to_string(),
            strand,
            kind,
            exons,
            inclusion_junctions,
            exclusion_junctions,
            retained_intron,
        }
    }
}

/// Enumerate all events across the annotation, optionally augmenting with
/// novel junctions empirically observed in the cohort.
pub fn detect_all(ann: &Annotation, junctions: &JunctionMatrix) -> Vec<ASEvent> {
    let empirical: HashSet<(String, Position, Position)> = junctions
        .counts
        .keys()
        .map(|j| (j.chrom.clone(), j.donor_end, j.acceptor_start))
        .collect();

    let mut all = Vec::new();
    for gene in ann.genes.values() {
        all.extend(detect_gene(gene, &empirical));
    }
    // De-duplicate by event_id — different traversals can rediscover the
    // same SE/MXE/etc.
    let mut seen = HashSet::new();
    all.retain(|e| seen.insert(e.event_id.clone()));
    all
}

fn detect_gene(g: &SpliceGraph, empirical: &HashSet<(String, Position, Position)>) -> Vec<ASEvent> {
    let mut events = Vec::new();
    detect_se(g, empirical, &mut events);
    detect_mxe(g, empirical, &mut events);
    detect_a5ss(g, &mut events);
    detect_a3ss(g, &mut events);
    detect_ir(g, &mut events);
    detect_afe_ale(g, &mut events);
    events
}

/// SE: upstream U → middle M → downstream D, with a direct U → D junction
/// (the "skip" form).
fn detect_se(
    g: &SpliceGraph,
    empirical: &HashSet<(String, Position, Position)>,
    out: &mut Vec<ASEvent>,
) {
    for u_idx in g.graph.node_indices() {
        let u = g.graph[u_idx];
        for um_edge in g.graph.edges(u_idx) {
            let m_idx = um_edge.target();
            let m = g.graph[m_idx];
            for md_edge in g.graph.edges(m_idx) {
                let d_idx = md_edge.target();
                let d = g.graph[d_idx];
                let u_d_annotated = g
                    .graph
                    .edges_connecting(u_idx, d_idx)
                    .next()
                    .is_some();
                let u_d_empirical =
                    empirical.contains(&(g.chrom.clone(), u.end, d.start));
                if u_d_annotated || u_d_empirical {
                    let inc = vec![(u.end, m.start), (m.end, d.start)];
                    let exc = vec![(u.end, d.start)];
                    out.push(ASEvent::new(
                        &g.gene_id,
                        &g.chrom,
                        g.strand,
                        EventKind::SE,
                        vec![u, m, d],
                        inc,
                        exc,
                        None,
                    ));
                }
            }
        }
    }
}

/// MXE: upstream U → either A1 → downstream D or A2 → D, with no path
/// A1 ↔ A2 inside the gene.
fn detect_mxe(
    g: &SpliceGraph,
    _empirical: &HashSet<(String, Position, Position)>,
    out: &mut Vec<ASEvent>,
) {
    for u_idx in g.graph.node_indices() {
        let u = g.graph[u_idx];
        let mids: Vec<_> = g
            .graph
            .edges(u_idx)
            .map(|e| (e.target(), g.graph[e.target()]))
            .collect();
        for i in 0..mids.len() {
            for j in (i + 1)..mids.len() {
                let (a1_idx, a1) = mids[i];
                let (a2_idx, a2) = mids[j];
                // Skip overlap or contiguity — MXE requires non-overlapping alts.
                if a1.end >= a2.start && a2.end >= a1.start {
                    continue;
                }
                // Both must lead to a common downstream exon D, and there
                // must be no edge between A1 and A2 in either direction.
                if g.graph.contains_edge(a1_idx, a2_idx)
                    || g.graph.contains_edge(a2_idx, a1_idx)
                {
                    continue;
                }
                let a1_targets: HashSet<_> = g.graph.neighbors(a1_idx).collect();
                let a2_targets: HashSet<_> = g.graph.neighbors(a2_idx).collect();
                for d_idx in a1_targets.intersection(&a2_targets) {
                    let d = g.graph[*d_idx];
                    // Order alternates so A1 is the 5' one in genomic coords.
                    let (alt1, alt2) = if a1.start < a2.start { (a1, a2) } else { (a2, a1) };
                    let inc = vec![(u.end, alt1.start), (alt1.end, d.start)];
                    let exc = vec![(u.end, alt2.start), (alt2.end, d.start)];
                    out.push(ASEvent::new(
                        &g.gene_id,
                        &g.chrom,
                        g.strand,
                        EventKind::MXE,
                        vec![u, alt1, alt2, d],
                        inc,
                        exc,
                        None,
                    ));
                }
            }
        }
    }
}

/// A5SS: two exons sharing the *acceptor* (start) but with different *donors*
/// (ends) on the + strand. Conceptually: shorter and longer donor variants of
/// the same upstream exon, both joining to the same downstream exon.
///
/// We detect this as: two distinct edges from two distinct exon nodes that
/// share `start` (the donor exon body) to the same downstream exon, where
/// the two donor exons differ only in `end`.
fn detect_a5ss(g: &SpliceGraph, out: &mut Vec<ASEvent>) {
    detect_alt_ss(g, out, /* alt_on_end = */ true);
}

/// A3SS: mirror — two acceptor exons sharing `end` but with different `start`,
/// both reachable from the same upstream donor.
fn detect_a3ss(g: &SpliceGraph, out: &mut Vec<ASEvent>) {
    detect_alt_ss(g, out, /* alt_on_end = */ false);
}

fn detect_alt_ss(g: &SpliceGraph, out: &mut Vec<ASEvent>, alt_on_end: bool) {
    // Group exons by the "anchor" side. For A5SS the anchor is `start`
    // (acceptor end stays fixed at the body's start), so we group donor
    // exons by start and find pairs with different ends.
    use std::collections::HashMap;
    let mut by_anchor: HashMap<Position, Vec<Exon>> = HashMap::new();
    for ex in g.exon_index.keys() {
        let anchor = if alt_on_end { ex.start } else { ex.end };
        by_anchor.entry(anchor).or_default().push(*ex);
    }

    for (_anchor, mut exons) in by_anchor {
        if exons.len() < 2 {
            continue;
        }
        exons.sort();
        for i in 0..exons.len() {
            for j in (i + 1)..exons.len() {
                let e1 = exons[i];
                let e2 = exons[j];
                // Identify the "long" and "short" variants.
                let (short, long) = if alt_on_end {
                    if e1.end < e2.end { (e1, e2) } else { (e2, e1) }
                } else {
                    if e1.start > e2.start { (e1, e2) } else { (e2, e1) }
                };
                let s_idx = match g.exon_index.get(&short) {
                    Some(x) => *x,
                    None => continue,
                };
                let l_idx = match g.exon_index.get(&long) {
                    Some(x) => *x,
                    None => continue,
                };

                // For A5SS: shared downstream exon reached from both via outgoing edge.
                // For A3SS: shared upstream exon reached from both via incoming edge.
                let s_neighbors: HashSet<_> = if alt_on_end {
                    g.graph.neighbors(s_idx).collect()
                } else {
                    g.graph
                        .neighbors_directed(s_idx, petgraph::Direction::Incoming)
                        .collect()
                };
                let l_neighbors: HashSet<_> = if alt_on_end {
                    g.graph.neighbors(l_idx).collect()
                } else {
                    g.graph
                        .neighbors_directed(l_idx, petgraph::Direction::Incoming)
                        .collect()
                };
                for partner_idx in s_neighbors.intersection(&l_neighbors) {
                    let partner = g.graph[*partner_idx];
                    let (long_donor_end, short_donor_end, partner_start_or_end) = if alt_on_end {
                        // A5SS, + strand convention: junction (donor.end → partner.start)
                        (long.end, short.end, partner.start)
                    } else {
                        // A3SS: junction (partner.end → acceptor.start)
                        (long.start, short.start, partner.end)
                    };
                    let (inc, exc) = if alt_on_end {
                        // Inclusion = short isoform retains more of exon; exclusion = long form
                        (
                            vec![(short_donor_end, partner_start_or_end)],
                            vec![(long_donor_end, partner_start_or_end)],
                        )
                    } else {
                        (
                            vec![(partner_start_or_end, short_donor_end)],
                            vec![(partner_start_or_end, long_donor_end)],
                        )
                    };
                    out.push(ASEvent::new(
                        &g.gene_id,
                        &g.chrom,
                        g.strand,
                        if alt_on_end { EventKind::A5SS } else { EventKind::A3SS },
                        vec![short, long, partner],
                        inc,
                        exc,
                        None,
                    ));
                }
            }
        }
    }
}

/// IR: an annotated intron between U and D where a third exon completely
/// covers that intron in some other transcript (i.e. the intron is "exonic"
/// in that other isoform). Inclusion is detected by intronic coverage; the
/// exclusion junction is the annotated intron.
fn detect_ir(g: &SpliceGraph, out: &mut Vec<ASEvent>) {
    let mut seen = HashSet::new();
    for edge in g.graph.edge_references() {
        let intron = edge.weight();
        let donor = g.graph[edge.source()];
        let acceptor = g.graph[edge.target()];
        let intron_lo = intron.donor_end + 1;
        let intron_hi = intron.acceptor_start.saturating_sub(1);
        if intron_lo >= intron_hi {
            continue;
        }
        // Look for any exon whose extent strictly covers the intron region.
        for ex in g.exon_index.keys() {
            if ex.start <= intron_lo && ex.end >= intron_hi {
                if ex.start == donor.start && ex.end == acceptor.end {
                    continue; // same as covering the whole donor+acceptor; ignore degenerate
                }
                if seen.insert((intron_lo, intron_hi)) {
                    let exc = vec![(intron.donor_end, intron.acceptor_start)];
                    let inc = vec![]; // measured via coverage in quantify stage
                    out.push(ASEvent::new(
                        &g.gene_id,
                        &g.chrom,
                        g.strand,
                        EventKind::IR,
                        vec![donor, *ex, acceptor],
                        inc,
                        exc,
                        Some((intron_lo, intron_hi)),
                    ));
                }
            }
        }
    }
}

/// AFE / ALE: identify exons that have no incoming edges (AFE candidates) or
/// no outgoing edges (ALE candidates), then pair those that share a common
/// neighbor on the other side.
fn detect_afe_ale(g: &SpliceGraph, out: &mut Vec<ASEvent>) {
    use petgraph::Direction;
    let mut firsts: Vec<_> = Vec::new();
    let mut lasts: Vec<_> = Vec::new();
    for idx in g.graph.node_indices() {
        if g.graph.neighbors_directed(idx, Direction::Incoming).count() == 0
            && g.graph.neighbors_directed(idx, Direction::Outgoing).count() > 0
        {
            firsts.push(idx);
        }
        if g.graph.neighbors_directed(idx, Direction::Outgoing).count() == 0
            && g.graph.neighbors_directed(idx, Direction::Incoming).count() > 0
        {
            lasts.push(idx);
        }
    }

    for i in 0..firsts.len() {
        for j in (i + 1)..firsts.len() {
            let f1 = g.graph[firsts[i]];
            let f2 = g.graph[firsts[j]];
            if f1.end >= f2.start && f2.end >= f1.start {
                continue; // overlapping = same exon really
            }
            let n1: HashSet<_> = g.graph.neighbors(firsts[i]).collect();
            let n2: HashSet<_> = g.graph.neighbors(firsts[j]).collect();
            for d_idx in n1.intersection(&n2) {
                let d = g.graph[*d_idx];
                let (a, b) = if f1.start < f2.start { (f1, f2) } else { (f2, f1) };
                let inc = vec![(a.end, d.start)];
                let exc = vec![(b.end, d.start)];
                out.push(ASEvent::new(
                    &g.gene_id,
                    &g.chrom,
                    g.strand,
                    EventKind::AFE,
                    vec![a, b, d],
                    inc,
                    exc,
                    None,
                ));
            }
        }
    }

    for i in 0..lasts.len() {
        for j in (i + 1)..lasts.len() {
            let l1 = g.graph[lasts[i]];
            let l2 = g.graph[lasts[j]];
            if l1.end >= l2.start && l2.end >= l1.start {
                continue;
            }
            let p1: HashSet<_> = g
                .graph
                .neighbors_directed(lasts[i], petgraph::Direction::Incoming)
                .collect();
            let p2: HashSet<_> = g
                .graph
                .neighbors_directed(lasts[j], petgraph::Direction::Incoming)
                .collect();
            for u_idx in p1.intersection(&p2) {
                let u = g.graph[*u_idx];
                let (a, b) = if l1.start < l2.start { (l1, l2) } else { (l2, l1) };
                let inc = vec![(u.end, a.start)];
                let exc = vec![(u.end, b.start)];
                out.push(ASEvent::new(
                    &g.gene_id,
                    &g.chrom,
                    g.strand,
                    EventKind::ALE,
                    vec![u, a, b],
                    inc,
                    exc,
                    None,
                ));
            }
        }
    }
}
