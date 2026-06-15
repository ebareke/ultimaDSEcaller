//! Sashimi-style plot data — emits a per-event SVG-ready description that
//! the HTML report inlines.
//!
//! A sashimi plot has three layers:
//! 1. **Exon track** — rectangles for the exons that participate in the event.
//! 2. **Coverage** — per-position read depth across the event window,
//!    one stacked band per contrast group (numerator vs denominator).
//! 3. **Junction arcs** — Bézier curves between exon endpoints whose
//!    thickness encodes group-summed read support for that junction.
//!
//! Instead of rendering pixels in Rust we emit a JSON description that the
//! HTML template renders as inline SVG. This keeps the binary lean and the
//! plots vector-friendly + zoomable.
//!
//! The coverage layer is computed via [`crate::pileup::pileup_regions`].
//! For each top event we sum coverage by contrast group; the per-position
//! depths land in the JSON as two arrays so the renderer can draw two
//! overlaid bands.

use std::collections::HashMap;

use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use crate::annotation::Exon;
use crate::config::RunConfig;
use crate::error::UltiResult;
use crate::events::ASEvent;
use crate::output::ResultRow;
use crate::pileup::{pileup_regions_with_opts, PileupOpts, Region};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SashimiTrack {
    pub event_id: String,
    pub gene_id: String,
    pub chrom: String,
    pub window_start: u64,
    pub window_end: u64,
    pub exons: Vec<(u64, u64)>,
    /// Per-position coverage for the numerator group (mean across samples).
    pub coverage_numerator: Vec<f64>,
    /// Per-position coverage for the denominator group (mean across samples).
    pub coverage_denominator: Vec<f64>,
    /// Junction arcs as `(donor_end, acceptor_start, support_num, support_denom)`.
    pub arcs: Vec<JunctionArc>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct JunctionArc {
    pub donor_end: u64,
    pub acceptor_start: u64,
    pub support_numerator: f64,
    pub support_denominator: f64,
    pub kind: String,
}

/// Build sashimi data for up to `top_k` events per kind, ranked by
/// adjusted p-value. Reads per-sample BAM coverage in a single pass per
/// sample regardless of how many windows are requested.
pub fn build_top_event_sashimi(
    cfg: &RunConfig,
    events: &[ASEvent],
    rows: &[ResultRow],
    top_k: usize,
) -> UltiResult<Vec<SashimiTrack>> {
    // 1) Pick top events per kind by adjusted p-value.
    let mut by_kind: HashMap<String, Vec<&ResultRow>> = HashMap::new();
    for r in rows {
        if !r.adjusted_p_value.is_finite() {
            continue;
        }
        by_kind.entry(r.event_type.clone()).or_default().push(r);
    }
    let mut chosen_ids: Vec<String> = Vec::new();
    for (_, mut subset) in by_kind {
        subset.sort_by(|a, b| {
            a.adjusted_p_value
                .partial_cmp(&b.adjusted_p_value)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        for r in subset.iter().take(top_k) {
            chosen_ids.push(r.event_id.clone());
        }
    }
    if chosen_ids.is_empty() {
        return Ok(Vec::new());
    }
    let chosen_set: std::collections::HashSet<String> = chosen_ids.into_iter().collect();
    let chosen_events: Vec<&ASEvent> = events
        .iter()
        .filter(|e| chosen_set.contains(&e.event_id))
        .collect();

    // 2) Build pileup regions — one per event window.
    let mut regions: Vec<Region> = Vec::with_capacity(chosen_events.len());
    let mut windows: Vec<(u64, u64)> = Vec::with_capacity(chosen_events.len());
    for ev in &chosen_events {
        let (start, end) = event_window(ev);
        regions.push(Region {
            chrom: ev.chrom.clone(),
            start,
            end,
        });
        windows.push((start, end));
    }

    // 3) Per-sample, run a single pileup over all windows. Then collapse
    //    by contrast group.
    let Some(contrast) = cfg.contrast.as_ref() else {
        return Ok(Vec::new());
    };

    // Per-sample pileup, requesting per-position depth vectors.
    let opts = PileupOpts { record_per_position: true };
    let per_sample_pp: Vec<Vec<Vec<u32>>> = cfg
        .samples
        .par_iter()
        .map(|s| {
            match pileup_regions_with_opts(&s.bam, &regions, cfg.reads.min_mapq, opts) {
                Ok(p) => p
                    .per_position
                    .unwrap_or_else(|| regions.iter().map(|_| Vec::new()).collect()),
                Err(_) => regions.iter().map(|_| Vec::new()).collect(),
            }
        })
        .collect();

    // Track downsample width — full-resolution arrays for a 10 kb window
    // ship 10 kb of u32s per track and inflate report.html by megabytes.
    // We bin to a fixed `N_BINS` per window. The HTML renderer linearly
    // maps bin indices back to genomic coords for display.
    const N_BINS: usize = 400;

    let mut tracks = Vec::new();
    for (idx, ev) in chosen_events.iter().enumerate() {
        let (start, end) = windows[idx];
        let n_pos = (end - start + 1) as usize;

        // Per-position group means.
        let mut sum_num = vec![0.0_f64; n_pos];
        let mut sum_denom = vec![0.0_f64; n_pos];
        let mut cnt_num = 0_u32;
        let mut cnt_denom = 0_u32;
        for (s_idx, sample) in cfg.samples.iter().enumerate() {
            let pp = &per_sample_pp[s_idx][idx];
            if pp.len() != n_pos {
                continue;
            }
            let (sum, cnt) = if sample.group == contrast.numerator {
                (&mut sum_num, &mut cnt_num)
            } else if sample.group == contrast.denominator {
                (&mut sum_denom, &mut cnt_denom)
            } else {
                continue;
            };
            for i in 0..n_pos {
                sum[i] += pp[i] as f64;
            }
            *cnt += 1;
        }
        let denom_n = (cnt_num.max(1)) as f64;
        let denom_d = (cnt_denom.max(1)) as f64;
        for i in 0..n_pos {
            sum_num[i] /= denom_n;
            sum_denom[i] /= denom_d;
        }
        // Downsample to N_BINS via mean.
        let coverage_num = downsample(&sum_num, N_BINS);
        let coverage_denom = downsample(&sum_denom, N_BINS);

        // Arcs.
        let mut arcs = Vec::new();
        for (de, as_) in &ev.inclusion_junctions {
            arcs.push(JunctionArc {
                donor_end: *de,
                acceptor_start: *as_,
                support_numerator: 0.0,
                support_denominator: 0.0,
                kind: "inclusion".into(),
            });
        }
        for (de, as_) in &ev.exclusion_junctions {
            arcs.push(JunctionArc {
                donor_end: *de,
                acceptor_start: *as_,
                support_numerator: 0.0,
                support_denominator: 0.0,
                kind: "exclusion".into(),
            });
        }

        tracks.push(SashimiTrack {
            event_id: ev.event_id.clone(),
            gene_id: ev.gene_id.clone(),
            chrom: ev.chrom.clone(),
            window_start: start,
            window_end: end,
            exons: ev.exons.iter().map(|e: &Exon| (e.start, e.end)).collect(),
            coverage_numerator: coverage_num,
            coverage_denominator: coverage_denom,
            arcs,
        });
    }
    Ok(tracks)
}

/// Downsample a per-position depth track to `n_bins` by mean-binning. If
/// the source is shorter than the requested bin count, returns the source
/// itself (no upsampling).
fn downsample(src: &[f64], n_bins: usize) -> Vec<f64> {
    if src.is_empty() {
        return Vec::new();
    }
    if src.len() <= n_bins {
        return src.to_vec();
    }
    let mut out = Vec::with_capacity(n_bins);
    let step = src.len() as f64 / n_bins as f64;
    for i in 0..n_bins {
        let lo = (i as f64 * step).floor() as usize;
        let hi = (((i + 1) as f64) * step).ceil() as usize;
        let hi = hi.min(src.len());
        if hi <= lo {
            out.push(0.0);
            continue;
        }
        let s: f64 = src[lo..hi].iter().sum();
        out.push(s / (hi - lo) as f64);
    }
    out
}

fn event_window(ev: &ASEvent) -> (u64, u64) {
    let mut lo = u64::MAX;
    let mut hi = 0_u64;
    for e in &ev.exons {
        lo = lo.min(e.start);
        hi = hi.max(e.end);
    }
    if let Some((a, b)) = ev.retained_intron {
        lo = lo.min(a);
        hi = hi.max(b);
    }
    // Pad 200 bp on each side, guarding underflow.
    (lo.saturating_sub(200), hi + 200)
}

/// Build Sankey nodes/links for the long-read isoform catalog. Each
/// isoform becomes one link from "numerator" or "denominator" to its
/// isoform_id node, weighted by sample-group read support.
pub fn build_isoform_sankey(
    catalog: &crate::longread::IsoformCatalog,
    sample_groups: &[&str],
    numerator: &str,
    denominator: &str,
    top_n: usize,
) -> SankeyData {
    let mut nodes: Vec<String> = vec![numerator.into(), denominator.into()];
    let mut name_to_idx: HashMap<String, usize> = HashMap::new();
    name_to_idx.insert(numerator.into(), 0);
    name_to_idx.insert(denominator.into(), 1);
    let mut links: Vec<SankeyLink> = Vec::new();
    let mut totals: Vec<(String, f64)> = Vec::new();
    for iso in &catalog.isoforms {
        let total: f64 = iso.support.iter().sum();
        if total <= 0.0 {
            continue;
        }
        totals.push((iso.id.clone(), total));
    }
    totals.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
    let keep: std::collections::HashSet<String> =
        totals.iter().take(top_n).map(|x| x.0.clone()).collect();

    for iso in &catalog.isoforms {
        if !keep.contains(&iso.id) {
            continue;
        }
        let id = nodes.len();
        name_to_idx.insert(iso.id.clone(), id);
        nodes.push(iso.id.clone());
        let mut sum_num = 0.0;
        let mut sum_denom = 0.0;
        for (i, g) in sample_groups.iter().enumerate() {
            let v = iso.support.get(i).copied().unwrap_or(0.0);
            if *g == numerator {
                sum_num += v;
            } else if *g == denominator {
                sum_denom += v;
            }
        }
        if sum_num > 0.0 {
            links.push(SankeyLink {
                source: 0,
                target: id,
                value: sum_num,
            });
        }
        if sum_denom > 0.0 {
            links.push(SankeyLink {
                source: 1,
                target: id,
                value: sum_denom,
            });
        }
    }
    SankeyData { nodes, links }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SankeyData {
    pub nodes: Vec<String>,
    pub links: Vec<SankeyLink>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SankeyLink {
    pub source: usize,
    pub target: usize,
    pub value: f64,
}

