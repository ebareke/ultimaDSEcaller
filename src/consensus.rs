//! Consensus statistical engine (spec §8.3).
//!
//! Combines per-event p-values from multiple, partially-correlated tests
//! (BB-LRT, GLM, GLMM, Fisher, DIU) into a single consensus q-value.
//! Three combination methods are supported:
//!
//! * **Stouffer's Z** (default) — weighted mean of z-scores.
//! * **Brown's method** — Fisher with a correlation correction.
//! * **Weighted Fisher** — classic χ² = −2 Σ wᵢ ln pᵢ.
//!
//! Final BH-FDR is applied across consensus p-values.

use serde::{Deserialize, Serialize};
use statrs::distribution::{ChiSquared, ContinuousCDF, Normal};

use crate::motif::SpliceMotif;
use crate::stats::bh_fdr;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
pub enum ConsensusMethod {
    #[default]
    Stouffer,
    Brown,
    WeightedFisher,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Weights {
    pub bb_lrt: f64,
    pub glm: f64,
    pub glmm: f64,
    pub fisher: f64,
    pub diu: f64,
}

impl Default for Weights {
    fn default() -> Self {
        Weights {
            bb_lrt: 1.0,
            glm: 1.0,
            glmm: 1.0,
            fisher: 0.5,
            diu: 0.5,
        }
    }
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct EventEvidence {
    pub p_bb_lrt: Option<f64>,
    pub p_glm: Option<f64>,
    pub p_glmm: Option<f64>,
    pub p_fisher: Option<f64>,
    pub p_diu: Option<f64>,
    pub motif: Option<SpliceMotif>,
    pub mean_coverage: f64,
    pub replicate_reproducibility: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ConsensusResult {
    pub combined_p: f64,
    pub consensus_q: f64,
    pub confidence: f64,
    pub methods_combined: u32,
    pub method: ConsensusMethod,
}

pub fn combine(
    evidence: &[EventEvidence],
    weights: &Weights,
    method: ConsensusMethod,
) -> Vec<ConsensusResult> {
    let combined: Vec<(f64, u32)> = evidence
        .iter()
        .map(|e| match method {
            ConsensusMethod::Stouffer => stouffer(e, weights),
            ConsensusMethod::Brown => brown(e, weights),
            ConsensusMethod::WeightedFisher => weighted_fisher(e, weights),
        })
        .collect();
    let ps: Vec<f64> = combined.iter().map(|(p, _)| *p).collect();
    let qs = bh_fdr(&ps);
    evidence
        .iter()
        .zip(combined.iter())
        .zip(qs.iter())
        .map(|((ev, &(p, k)), &q)| ConsensusResult {
            combined_p: p,
            consensus_q: q,
            confidence: confidence_score(ev, q),
            methods_combined: k,
            method,
        })
        .collect()
}

fn collect_inputs(ev: &EventEvidence, w: &Weights) -> Vec<(f64, f64)> {
    let mut out: Vec<(f64, f64)> = Vec::new();
    for (p, weight) in [
        (ev.p_bb_lrt, w.bb_lrt),
        (ev.p_glm, w.glm),
        (ev.p_glmm, w.glmm),
        (ev.p_fisher, w.fisher),
        (ev.p_diu, w.diu),
    ] {
        if let Some(p) = p {
            if p.is_finite() && weight > 0.0 {
                out.push((p.clamp(1e-300, 1.0), weight));
            }
        }
    }
    out
}

fn stouffer(ev: &EventEvidence, w: &Weights) -> (f64, u32) {
    let inputs = collect_inputs(ev, w);
    if inputs.is_empty() {
        return (f64::NAN, 0);
    }
    let normal = Normal::standard();
    let mut z_sum = 0.0_f64;
    let mut w2_sum = 0.0_f64;
    for (p, wi) in &inputs {
        let z = normal.inverse_cdf((1.0 - p / 2.0).clamp(1e-12, 1.0 - 1e-12));
        z_sum += wi * z;
        w2_sum += wi * wi;
    }
    let z_combined = z_sum / w2_sum.sqrt().max(1e-12);
    let p_combined = 2.0 * (1.0 - normal.cdf(z_combined.abs()));
    (p_combined.clamp(0.0, 1.0), inputs.len() as u32)
}

fn weighted_fisher(ev: &EventEvidence, w: &Weights) -> (f64, u32) {
    let inputs = collect_inputs(ev, w);
    if inputs.is_empty() {
        return (f64::NAN, 0);
    }
    let k = inputs.len() as f64;
    let w_sum: f64 = inputs.iter().map(|(_, w)| *w).sum();
    let scale = k / w_sum;
    let chi2: f64 = inputs
        .iter()
        .map(|(p, wi)| -2.0 * (wi * scale) * p.ln())
        .sum();
    let df = 2.0 * k;
    let dist = match ChiSquared::new(df) {
        Ok(d) => d,
        Err(_) => return (1.0, inputs.len() as u32),
    };
    let p = (1.0 - dist.cdf(chi2)).clamp(0.0, 1.0);
    (p, inputs.len() as u32)
}

fn brown(ev: &EventEvidence, w: &Weights) -> (f64, u32) {
    let inputs = collect_inputs(ev, w);
    if inputs.is_empty() {
        return (f64::NAN, 0);
    }
    let k = inputs.len() as f64;
    let rho = std::env::var("ULTIMADSE_BROWN_RHO")
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(0.3)
        .clamp(0.0, 0.99);
    let chi2: f64 = inputs.iter().map(|(p, _)| -2.0 * p.ln()).sum();
    let var = 4.0 * k + 2.0 * k * (k - 1.0) * (3.25 + 0.75 * rho) * rho;
    let c = var / (4.0 * k);
    let f = (2.0 * k).powi(2) / var;
    let scaled = chi2 / c;
    let dist = match ChiSquared::new(f) {
        Ok(d) => d,
        Err(_) => return (1.0, inputs.len() as u32),
    };
    let p = (1.0 - dist.cdf(scaled)).clamp(0.0, 1.0);
    (p, inputs.len() as u32)
}

fn confidence_score(ev: &EventEvidence, q: f64) -> f64 {
    let sig_term = if q.is_finite() {
        1.0 - q.clamp(0.0, 1.0)
    } else {
        0.0
    };
    let motif_term = match ev.motif {
        Some(SpliceMotif::GtAg | SpliceMotif::GtAgReverse) => 1.0,
        Some(SpliceMotif::GcAg | SpliceMotif::GcAgReverse) => 0.85,
        Some(SpliceMotif::AtAc | SpliceMotif::AtAcReverse) => 0.8,
        Some(SpliceMotif::NonCanonical) => 0.4,
        _ => 0.7,
    };
    let cov_term = 1.0 - (-ev.mean_coverage / 20.0).exp();
    let repro = ev.replicate_reproducibility.clamp(0.0, 1.0);
    let pieces = [sig_term, motif_term, cov_term, repro];
    let n = pieces.len() as f64;
    pieces
        .iter()
        .map(|x| x.max(1e-6))
        .product::<f64>()
        .powf(1.0 / n)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn stouffer_combines_strong_evidence() {
        let ev = EventEvidence {
            p_bb_lrt: Some(1e-6),
            p_glm: Some(1e-6),
            p_glmm: Some(1e-6),
            p_fisher: Some(1e-5),
            motif: Some(SpliceMotif::GtAg),
            mean_coverage: 40.0,
            replicate_reproducibility: 0.9,
            ..Default::default()
        };
        let res = combine(&[ev], &Weights::default(), ConsensusMethod::Stouffer);
        assert!(res[0].combined_p < 1e-8, "got {}", res[0].combined_p);
        assert_eq!(res[0].methods_combined, 4);
    }

    #[test]
    fn stouffer_combines_weak_evidence() {
        let ev = EventEvidence {
            p_bb_lrt: Some(0.5),
            p_glm: Some(0.4),
            p_glmm: Some(0.6),
            p_fisher: Some(0.5),
            ..Default::default()
        };
        let res = combine(&[ev], &Weights::default(), ConsensusMethod::Stouffer);
        assert!(res[0].combined_p > 0.05, "got {}", res[0].combined_p);
    }

    #[test]
    fn weighted_fisher_classical_recovery() {
        let ev = EventEvidence {
            p_bb_lrt: Some(0.001),
            p_glm: Some(0.01),
            ..Default::default()
        };
        let res = combine(&[ev], &Weights::default(), ConsensusMethod::WeightedFisher);
        assert!(res[0].combined_p < 1e-3, "got {}", res[0].combined_p);
    }

    #[test]
    fn brown_at_least_as_conservative_as_fisher() {
        let ev = EventEvidence {
            p_bb_lrt: Some(1e-4),
            p_glm: Some(1e-4),
            p_glmm: Some(1e-4),
            ..Default::default()
        };
        let fisher = combine(
            std::slice::from_ref(&ev),
            &Weights::default(),
            ConsensusMethod::WeightedFisher,
        );
        let brown = combine(&[ev], &Weights::default(), ConsensusMethod::Brown);
        assert!(brown[0].combined_p >= fisher[0].combined_p - 1e-9);
    }
}
