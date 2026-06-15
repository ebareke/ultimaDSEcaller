//! Junction → event quantification. For each detected event, look up the
//! event's inclusion / exclusion junctions in the cohort-wide junction
//! matrix, build per-sample inclusion/exclusion counts, and compute
//! summary statistics (PSI, ΔPSI, coverage, confidence, complexity,
//! reproducibility).
//!
//! Group aggregation happens here too — the statistical test downstream
//! only needs to know "k successes out of n trials in group A vs group B".

use std::collections::HashMap;

use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use crate::cli::Technology;
use crate::config::{Contrast, RunConfig};
use crate::events::{ASEvent, EventKind};
use crate::junctions::JunctionMatrix;
use crate::pileup::{pileup_regions, Region};

/// Per-event quantification for the whole cohort.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EventQuant {
    pub event: ASEvent,
    /// Per-sample PSI in the same order as `JunctionMatrix::samples`.
    /// `NaN` where total coverage is zero.
    pub psi: Vec<f64>,
    /// Per-sample inclusion + exclusion counts (the "n" in beta-binomial).
    pub inclusion: Vec<f64>,
    pub exclusion: Vec<f64>,
    /// Group-level summaries for the active contrast (None if no contrast set).
    pub contrast_summary: Option<ContrastSummary>,
    /// Number of distinct junction structures involved (≥ 2 for any AS event).
    pub complexity: u32,
    /// 1 − mean(CV of PSI within group). Higher = more reproducible.
    pub reproducibility: f64,
    /// Composite confidence in (0, 1) — see [`compute_confidence`].
    pub confidence: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ContrastSummary {
    pub numerator_group: String,
    pub denominator_group: String,
    pub mean_psi_num: f64,
    pub mean_psi_denom: f64,
    pub delta_psi: f64,
    pub mean_coverage_num: f64,
    pub mean_coverage_denom: f64,
    /// Per-sample inclusion / exclusion split by group; passed to the stats engine.
    pub inclusion_num: Vec<f64>,
    pub exclusion_num: Vec<f64>,
    pub inclusion_denom: Vec<f64>,
    pub exclusion_denom: Vec<f64>,
}

/// Quantify every detected event across the cohort.
pub fn quantify(
    cfg: &RunConfig,
    events: &[ASEvent],
    jm: &JunctionMatrix,
) -> Vec<EventQuant> {
    // Pre-index the junction matrix for O(1) lookup.
    let chrom_lookup: HashMap<(String, u64, u64), &Vec<f64>> = jm
        .counts
        .iter()
        .map(|(j, v)| ((j.chrom.clone(), j.donor_end, j.acceptor_start), v))
        .collect();

    let sample_to_group: Vec<&str> = cfg
        .samples
        .iter()
        .map(|s| s.group.as_str())
        .collect();
    let n = cfg.samples.len();

    // For every IR event with a retained-intron interval, compute per-sample
    // mean depth in the intron via pileup. We pay one BAM pass per sample but
    // batch all IR intervals into that single pass.
    let ir_coverage = compute_ir_coverage(cfg, events);

    events
        .par_iter()
        .map(|ev| {
            // Sum the supports of all inclusion / exclusion junctions.
            let mut inclusion = vec![0.0_f64; n];
            let mut exclusion = vec![0.0_f64; n];
            for (de, as_) in &ev.inclusion_junctions {
                if let Some(v) = chrom_lookup.get(&(ev.chrom.clone(), *de, *as_)) {
                    for i in 0..n {
                        inclusion[i] += v[i];
                    }
                }
            }
            for (de, as_) in &ev.exclusion_junctions {
                if let Some(v) = chrom_lookup.get(&(ev.chrom.clone(), *de, *as_)) {
                    for i in 0..n {
                        exclusion[i] += v[i];
                    }
                }
            }

            // IR events: replace junction-based inclusion (always zero) with
            // a normalized intronic-coverage estimate. The normalization
            // converts mean read depth into an "effective number of
            // inclusion reads" comparable to junction read counts.
            //
            //   inclusion ≈ mean_depth · intron_length / effective_read_length
            //
            // The effective read length defaults to 100 (short) or 1000
            // (long) per `cfg.tech`. Users tuning for unusual library types
            // can override via the IR_READ_LEN env var (handled in
            // [`effective_read_length`]).
            if ev.kind == EventKind::IR {
                if let Some(intron) = ev.retained_intron {
                    let len = (intron.1 - intron.0 + 1) as f64;
                    let read_len = effective_read_length(cfg.tech);
                    let key = (ev.chrom.clone(), intron.0, intron.1);
                    if let Some(depths) = ir_coverage.get(&key) {
                        for i in 0..n {
                            inclusion[i] = depths[i] * len / read_len;
                        }
                    }
                }
            }

            let psi: Vec<f64> = inclusion
                .iter()
                .zip(exclusion.iter())
                .map(|(i, e)| {
                    let t = i + e;
                    if t > 0.0 {
                        i / t
                    } else {
                        f64::NAN
                    }
                })
                .collect();

            let contrast_summary = cfg.contrast.as_ref().map(|c| {
                summarize_contrast(c, &sample_to_group, &inclusion, &exclusion, &psi)
            });

            let complexity =
                (ev.inclusion_junctions.len() + ev.exclusion_junctions.len()) as u32;
            let reproducibility = within_group_reproducibility(&sample_to_group, &psi);
            let confidence = compute_confidence(
                &inclusion,
                &exclusion,
                contrast_summary.as_ref().map(|s| s.delta_psi).unwrap_or(0.0),
                reproducibility,
            );

            EventQuant {
                event: ev.clone(),
                psi,
                inclusion,
                exclusion,
                contrast_summary,
                complexity,
                reproducibility,
                confidence,
            }
        })
        .collect()
}

fn summarize_contrast(
    c: &Contrast,
    sample_groups: &[&str],
    inclusion: &[f64],
    exclusion: &[f64],
    psi: &[f64],
) -> ContrastSummary {
    let (mut inc_num, mut inc_denom) = (Vec::new(), Vec::new());
    let (mut exc_num, mut exc_denom) = (Vec::new(), Vec::new());
    let (mut psi_num, mut psi_denom) = (Vec::new(), Vec::new());
    let (mut cov_num, mut cov_denom) = (Vec::new(), Vec::new());

    for (i, g) in sample_groups.iter().enumerate() {
        let cov = inclusion[i] + exclusion[i];
        if *g == c.numerator {
            inc_num.push(inclusion[i]);
            exc_num.push(exclusion[i]);
            psi_num.push(psi[i]);
            cov_num.push(cov);
        } else if *g == c.denominator {
            inc_denom.push(inclusion[i]);
            exc_denom.push(exclusion[i]);
            psi_denom.push(psi[i]);
            cov_denom.push(cov);
        }
    }

    let mean_psi_num = nanmean(&psi_num);
    let mean_psi_denom = nanmean(&psi_denom);

    ContrastSummary {
        numerator_group: c.numerator.clone(),
        denominator_group: c.denominator.clone(),
        mean_psi_num,
        mean_psi_denom,
        delta_psi: mean_psi_num - mean_psi_denom,
        mean_coverage_num: nanmean(&cov_num),
        mean_coverage_denom: nanmean(&cov_denom),
        inclusion_num: inc_num,
        exclusion_num: exc_num,
        inclusion_denom: inc_denom,
        exclusion_denom: exc_denom,
    }
}

fn nanmean(xs: &[f64]) -> f64 {
    let mut sum = 0.0;
    let mut n = 0u32;
    for x in xs {
        if x.is_finite() {
            sum += *x;
            n += 1;
        }
    }
    if n == 0 {
        f64::NAN
    } else {
        sum / n as f64
    }
}

fn within_group_reproducibility(sample_groups: &[&str], psi: &[f64]) -> f64 {
    let mut by_group: HashMap<&str, Vec<f64>> = HashMap::new();
    for (g, p) in sample_groups.iter().zip(psi.iter()) {
        if p.is_finite() {
            by_group.entry(*g).or_default().push(*p);
        }
    }
    if by_group.is_empty() {
        return 0.0;
    }
    let mut cvs = Vec::new();
    for vals in by_group.values() {
        if vals.len() < 2 {
            continue;
        }
        let m = nanmean(vals);
        if !m.is_finite() || m.abs() < 1e-9 {
            continue;
        }
        let var = vals.iter().map(|x| (x - m).powi(2)).sum::<f64>() / (vals.len() - 1) as f64;
        let sd = var.sqrt();
        cvs.push(sd / m.abs());
    }
    if cvs.is_empty() {
        return 0.5; // not enough info — neutral
    }
    let mean_cv = cvs.iter().copied().sum::<f64>() / cvs.len() as f64;
    (1.0 - mean_cv.min(1.0)).max(0.0)
}

/// Composite per-event confidence in (0, 1). Components:
///
/// * **coverage_term**: 1 − exp(−mean_coverage / 20). Reaches ~0.6 at 20×,
///   ~0.86 at 40×, asymptotes to 1.0.
/// * **effect_term**: tanh(2 · |ΔPSI|). Penalizes events whose PSI shift is
///   indistinguishable from noise even before the formal test.
/// * **reproducibility**: 0..1, already on the right scale.
///
/// The three are combined as a weighted geometric mean to penalize any
/// single weak component.
///
/// > **Open design choice (flagged for user input):** the *exact* functional
/// > form and the relative weighting are reasonable defaults, not a derived
/// > optimum. Different cohorts may want different mixes (e.g. low-coverage
/// > pilot studies may want coverage to dominate; high-replicate cohorts may
/// > want reproducibility to dominate). The function is `pub` so external
/// > callers can override.
pub fn compute_confidence(
    inclusion: &[f64],
    exclusion: &[f64],
    delta_psi: f64,
    reproducibility: f64,
) -> f64 {
    let mean_cov = inclusion
        .iter()
        .zip(exclusion.iter())
        .map(|(i, e)| i + e)
        .sum::<f64>()
        / inclusion.len().max(1) as f64;
    let coverage_term = 1.0 - (-mean_cov / 20.0).exp();
    let effect_term = (2.0 * delta_psi.abs()).tanh();
    let repro = reproducibility.clamp(0.0, 1.0);

    // Geometric mean with weights — penalizes any weak component.
    let (wc, we, wr) = (1.0, 1.0, 1.0);
    let total_w = wc + we + wr;
    (coverage_term.max(1e-6).powf(wc / total_w))
        * (effect_term.max(1e-6).powf(we / total_w))
        * (repro.max(1e-6).powf(wr / total_w))
}

/// Effective read length used to convert intronic depth into an "effective
/// inclusion read count" for IR PSI. Reasonable defaults per technology.
/// Override via the `ULTIMADSE_IR_READ_LEN` env var if your library is
/// atypical.
pub fn effective_read_length(tech: Technology) -> f64 {
    if let Ok(s) = std::env::var("ULTIMADSE_IR_READ_LEN") {
        if let Ok(v) = s.parse::<f64>() {
            if v > 0.0 {
                return v;
            }
        }
    }
    match tech {
        Technology::Short => 100.0,
        Technology::PacBio | Technology::Ont => 1000.0,
    }
}

/// Run pileup over every IR event's retained intron, for every sample.
/// Returns a map `(chrom, intron_start, intron_end) -> per-sample mean depth`.
fn compute_ir_coverage(
    cfg: &RunConfig,
    events: &[ASEvent],
) -> HashMap<(String, u64, u64), Vec<f64>> {
    let regions: Vec<Region> = events
        .iter()
        .filter_map(|e| {
            if e.kind == EventKind::IR {
                e.retained_intron.map(|(s, t)| Region {
                    chrom: e.chrom.clone(),
                    start: s,
                    end: t,
                })
            } else {
                None
            }
        })
        .collect();

    if regions.is_empty() {
        return HashMap::new();
    }

    // Deduplicate regions while keeping the order stable.
    let mut seen: HashMap<(String, u64, u64), usize> = HashMap::new();
    let mut unique_regions: Vec<Region> = Vec::new();
    for r in &regions {
        let k = (r.chrom.clone(), r.start, r.end);
        if !seen.contains_key(&k) {
            seen.insert(k, unique_regions.len());
            unique_regions.push(r.clone());
        }
    }

    let n_samples = cfg.samples.len();
    // Per-sample pileup in parallel.
    let per_sample: Vec<Vec<f64>> = cfg
        .samples
        .par_iter()
        .map(|s| {
            match pileup_regions(&s.bam, &unique_regions, cfg.reads.min_mapq) {
                Ok(p) => p.mean_depth,
                Err(_) => vec![0.0; unique_regions.len()],
            }
        })
        .collect();

    // Transpose into the (region -> per-sample) shape.
    let mut out: HashMap<(String, u64, u64), Vec<f64>> = HashMap::new();
    for (region_idx, r) in unique_regions.iter().enumerate() {
        let mut v = vec![0.0_f64; n_samples];
        for s in 0..n_samples {
            v[s] = per_sample[s][region_idx];
        }
        out.insert((r.chrom.clone(), r.start, r.end), v);
    }
    out
}

/// For tests / future work — drop events that don't pass the configured
/// per-sample coverage minimum in *both* contrast groups.
pub fn meets_coverage(quant: &EventQuant, min_cov_per_sample: u32) -> bool {
    let Some(cs) = quant.contrast_summary.as_ref() else {
        return true;
    };
    let min_cov_f = min_cov_per_sample as f64;
    cs.inclusion_num
        .iter()
        .zip(cs.exclusion_num.iter())
        .all(|(i, e)| i + e >= min_cov_f)
        && cs
            .inclusion_denom
            .iter()
            .zip(cs.exclusion_denom.iter())
            .all(|(i, e)| i + e >= min_cov_f)
}

