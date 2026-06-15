//! Statistical engine.
//!
//! ## With replicates: beta-binomial likelihood-ratio test
//!
//! For each event we have, per group, a vector of (inclusion, exclusion)
//! counts across replicates. Model each sample's inclusion count as
//! Beta-Binomial(n_i, α, β) where n_i = inclusion + exclusion. The
//! Beta-Binomial captures over-dispersion that a plain Binomial cannot —
//! crucial because PSI estimates have biological variability beyond the
//! purely combinatorial.
//!
//! Test statistic:
//!
//! ```text
//!   D = 2 · (LL_alt − LL_null)        ~  χ²(df = 2)
//! ```
//!
//! where LL_null fits one (α, β) to both groups pooled and LL_alt fits
//! independent (α, β) per group.
//!
//! ## Without replicates: Fisher's exact 2 × 2
//!
//! When a group has only one sample, the beta-binomial dispersion is
//! unidentifiable. We fall back to a two-sided Fisher's exact on the
//! 2 × 2 table of (inclusion_A, exclusion_A, inclusion_B, exclusion_B),
//! which is the textbook exact test for differential splice-junction usage.
//!
//! ## Multiple-testing correction: Benjamini-Hochberg FDR.
//!
//! ## Optimization
//!
//! MLE is found by Nelder-Mead in (log α, log β) space — keeps the
//! parameters positive without needing bounded optimizers and converges
//! reliably for the relatively flat BB log-likelihood surface.

use rayon::prelude::*;
use statrs::distribution::{ChiSquared, ContinuousCDF, Discrete, Hypergeometric};
use statrs::function::gamma::ln_gamma;

use crate::cli::TestMethod;
use crate::config::RunConfig;
use crate::glm::{fit_glm, fit_glmm, wald_test};
use crate::quantify::EventQuant;

#[derive(Debug, Clone, Copy)]
pub struct PValue {
    pub p_value: f64,
    pub adjusted_p_value: f64,
    pub effect_size: f64,
    pub test_used: TestUsed,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TestUsed {
    BetaBinomialLRT,
    FisherExact,
    Glm,
    Glmm,
    Insufficient,
}

impl Default for TestUsed {
    fn default() -> Self {
        TestUsed::Insufficient
    }
}

/// Run the statistical test on every quantified event in parallel, then
/// adjust p-values across the whole set.
pub fn test_all(quants: &[EventQuant]) -> Vec<PValue> {
    test_all_with_method(quants, None)
}

/// Variant that lets the caller select the test method. When `cfg` is
/// supplied with `cfg.test != BbLrt`, GLM or GLMM is run instead of the
/// beta-binomial LRT. The Fisher fallback for under-replicated rows still
/// applies.
pub fn test_all_with_method(quants: &[EventQuant], cfg: Option<&RunConfig>) -> Vec<PValue> {
    let method = cfg.map(|c| c.test).unwrap_or(TestMethod::BbLrt);
    let random_groups: Option<Vec<usize>> = match (method, cfg) {
        (TestMethod::Glmm, Some(c)) => Some(build_random_groups(c)),
        _ => None,
    };

    let raw: Vec<(f64, f64, TestUsed)> = quants
        .par_iter()
        .map(|q| {
            let Some(cs) = q.contrast_summary.as_ref() else {
                return (f64::NAN, 0.0, TestUsed::Insufficient);
            };
            let (p, used) = match method {
                TestMethod::BbLrt => test_one(
                    &cs.inclusion_num,
                    &cs.exclusion_num,
                    &cs.inclusion_denom,
                    &cs.exclusion_denom,
                ),
                TestMethod::Glm => test_one_glm(
                    &cs.inclusion_num,
                    &cs.exclusion_num,
                    &cs.inclusion_denom,
                    &cs.exclusion_denom,
                ),
                TestMethod::Glmm => test_one_glmm(
                    &cs.inclusion_num,
                    &cs.exclusion_num,
                    &cs.inclusion_denom,
                    &cs.exclusion_denom,
                    random_groups.as_deref(),
                ),
            };
            (p, cs.delta_psi, used)
        })
        .collect();

    let p_values: Vec<f64> = raw.iter().map(|t| t.0).collect();
    let adjusted = bh_fdr(&p_values);

    raw.into_iter()
        .zip(adjusted.into_iter())
        .map(|((p, eff, used), q)| PValue {
            p_value: p,
            adjusted_p_value: q,
            effect_size: eff,
            test_used: used,
        })
        .collect()
}

/// Multi-test result: every applicable test method's p-value for one event.
/// Fields are `None` when the method couldn't be applied (e.g. GLMM without
/// a random-effect column, or any test where the contrast is absent).
#[derive(Debug, Clone, Default)]
pub struct MultiTestResult {
    pub p_bb_lrt: Option<f64>,
    pub p_glm: Option<f64>,
    pub p_glmm: Option<f64>,
    pub p_fisher: Option<f64>,
    pub primary: TestUsed,
    pub primary_p: f64,
}

/// Run **all** applicable test methods on every event simultaneously.
/// This is what the consensus engine consumes when you want to combine
/// independent statistical signals into one q-value.
///
/// Methods that can't apply to the cohort (e.g. GLMM with no random
/// effect column, BB-LRT with `n_a < 2`) are skipped, not faked.
pub fn test_all_multi(quants: &[EventQuant], cfg: Option<&RunConfig>) -> Vec<MultiTestResult> {
    let random_groups: Option<Vec<usize>> = cfg
        .filter(|c| c.random_effect.is_some())
        .map(build_random_groups);
    let primary_method = cfg.map(|c| c.test).unwrap_or(TestMethod::BbLrt);

    quants
        .par_iter()
        .map(|q| {
            let mut out = MultiTestResult::default();
            let Some(cs) = q.contrast_summary.as_ref() else {
                out.primary = TestUsed::Insufficient;
                out.primary_p = f64::NAN;
                return out;
            };
            let (inc_a, exc_a, inc_b, exc_b) = (
                &cs.inclusion_num,
                &cs.exclusion_num,
                &cs.inclusion_denom,
                &cs.exclusion_denom,
            );
            // BB-LRT (only when ≥ 2 replicates per group).
            if inc_a.len() >= 2 && inc_b.len() >= 2 {
                let (p, used) = test_one(inc_a, exc_a, inc_b, exc_b);
                if used == TestUsed::BetaBinomialLRT && p.is_finite() {
                    out.p_bb_lrt = Some(p);
                }
            }
            // GLM (always applicable when there's at least one sample/group).
            {
                let (p, used) = test_one_glm(inc_a, exc_a, inc_b, exc_b);
                if used == TestUsed::Glm && p.is_finite() {
                    out.p_glm = Some(p);
                }
            }
            // GLMM (only when --random-effect supplies grouping).
            if let Some(g) = random_groups.as_deref() {
                let (p, used) = test_one_glmm(inc_a, exc_a, inc_b, exc_b, Some(g));
                if used == TestUsed::Glmm && p.is_finite() {
                    out.p_glmm = Some(p);
                }
            }
            // Fisher 2x2 (always defined — collapses replicates).
            let p_f = fisher_exact_2x2(
                total_inc(inc_a),
                total_inc(exc_a),
                total_inc(inc_b),
                total_inc(exc_b),
            );
            if p_f.is_finite() {
                out.p_fisher = Some(p_f);
            }
            // Primary = the method the user selected via --test.
            let (p_primary, used_primary) = match primary_method {
                TestMethod::BbLrt => test_one(inc_a, exc_a, inc_b, exc_b),
                TestMethod::Glm => test_one_glm(inc_a, exc_a, inc_b, exc_b),
                TestMethod::Glmm => {
                    test_one_glmm(inc_a, exc_a, inc_b, exc_b, random_groups.as_deref())
                }
            };
            out.primary = used_primary;
            out.primary_p = p_primary;
            out
        })
        .collect()
}

pub fn test_one(
    inc_a: &[f64],
    exc_a: &[f64],
    inc_b: &[f64],
    exc_b: &[f64],
) -> (f64, TestUsed) {
    let n_a = inc_a.len();
    let n_b = inc_b.len();
    if n_a == 0 || n_b == 0 {
        return (f64::NAN, TestUsed::Insufficient);
    }

    // Coverage check: at least one informative sample per group.
    let total_a: f64 = inc_a.iter().zip(exc_a.iter()).map(|(i, e)| i + e).sum();
    let total_b: f64 = inc_b.iter().zip(exc_b.iter()).map(|(i, e)| i + e).sum();
    if total_a < 1.0 || total_b < 1.0 {
        return (1.0, TestUsed::Insufficient);
    }

    if n_a < 2 || n_b < 2 {
        let p = fisher_exact_2x2(
            total_inc(inc_a),
            total_inc(exc_a),
            total_inc(inc_b),
            total_inc(exc_b),
        );
        return (p, TestUsed::FisherExact);
    }

    // Beta-binomial LRT.
    let (k_a, n_obs_a) = inc_n(inc_a, exc_a);
    let (k_b, n_obs_b) = inc_n(inc_b, exc_b);
    let mut k_all = k_a.clone();
    k_all.extend_from_slice(&k_b);
    let mut n_all = n_obs_a.clone();
    n_all.extend_from_slice(&n_obs_b);

    let (a_null, b_null, ll_null) = fit_beta_binomial(&k_all, &n_all);
    let (_, _, ll_a) = fit_beta_binomial(&k_a, &n_obs_a);
    let (_, _, ll_b) = fit_beta_binomial(&k_b, &n_obs_b);
    let ll_alt = ll_a + ll_b;

    let d = 2.0 * (ll_alt - ll_null);
    if !d.is_finite() || d <= 0.0 {
        let _ = (a_null, b_null);
        return (1.0, TestUsed::BetaBinomialLRT);
    }
    let chisq = ChiSquared::new(2.0).expect("df=2 chi-squared");
    let p = 1.0 - chisq.cdf(d);
    (p.clamp(0.0, 1.0), TestUsed::BetaBinomialLRT)
}

/// GLM-based per-event test: logistic regression of (k, n-k) on a single
/// binary indicator for `numerator` vs `denominator` group membership.
pub fn test_one_glm(
    inc_a: &[f64],
    exc_a: &[f64],
    inc_b: &[f64],
    exc_b: &[f64],
) -> (f64, TestUsed) {
    let n_a = inc_a.len();
    let n_b = inc_b.len();
    if n_a == 0 || n_b == 0 {
        return (f64::NAN, TestUsed::Insufficient);
    }
    let (y, n, x) = build_glm_inputs(inc_a, exc_a, inc_b, exc_b);
    match fit_glm(&y, &n, &x) {
        Ok(fit) => {
            let w = wald_test(&fit.beta, &fit.vcov, &[0.0, 1.0]);
            (w.p_value.clamp(0.0, 1.0), TestUsed::Glm)
        }
        Err(_) => (1.0, TestUsed::Insufficient),
    }
}

/// GLMM-based per-event test. `random_groups` is the per-sample
/// random-intercept group index, ordered as `(samples_in_group_A then
/// samples_in_group_B)`. If `None`, falls back to a plain GLM.
pub fn test_one_glmm(
    inc_a: &[f64],
    exc_a: &[f64],
    inc_b: &[f64],
    exc_b: &[f64],
    random_groups: Option<&[usize]>,
) -> (f64, TestUsed) {
    let n_a = inc_a.len();
    let n_b = inc_b.len();
    if n_a == 0 || n_b == 0 {
        return (f64::NAN, TestUsed::Insufficient);
    }
    let (y, n, x) = build_glm_inputs(inc_a, exc_a, inc_b, exc_b);
    let groups: Vec<usize> = match random_groups {
        Some(g) if g.len() == y.len() => g.to_vec(),
        _ => return test_one_glm(inc_a, exc_a, inc_b, exc_b),
    };
    match fit_glmm(&y, &n, &x, &groups) {
        Ok(fit) => {
            let w = wald_test(&fit.beta, &fit.vcov, &[0.0, 1.0]);
            (w.p_value.clamp(0.0, 1.0), TestUsed::Glmm)
        }
        Err(_) => (1.0, TestUsed::Insufficient),
    }
}

fn build_glm_inputs(
    inc_a: &[f64],
    exc_a: &[f64],
    inc_b: &[f64],
    exc_b: &[f64],
) -> (Vec<f64>, Vec<f64>, Vec<Vec<f64>>) {
    let mut y = Vec::with_capacity(inc_a.len() + inc_b.len());
    let mut n = Vec::with_capacity(inc_a.len() + inc_b.len());
    let mut x = Vec::with_capacity(inc_a.len() + inc_b.len());
    for i in 0..inc_a.len() {
        y.push(inc_a[i].max(0.0));
        n.push((inc_a[i] + exc_a[i]).max(0.0));
        x.push(vec![1.0, 0.0]); // intercept + treatment=0
    }
    for i in 0..inc_b.len() {
        y.push(inc_b[i].max(0.0));
        n.push((inc_b[i] + exc_b[i]).max(0.0));
        x.push(vec![1.0, 1.0]); // intercept + treatment=1
    }
    (y, n, x)
}

/// Build the per-sample random-effect group vector, restricted to the
/// samples participating in the active contrast, in the same order
/// (`numerator` samples first, `denominator` samples second) that
/// `build_glm_inputs` produces.
pub fn build_random_groups(cfg: &RunConfig) -> Vec<usize> {
    let re_col = match &cfg.random_effect {
        Some(s) => s,
        None => return Vec::new(),
    };
    let contrast = match &cfg.contrast {
        Some(c) => c,
        None => return Vec::new(),
    };
    let mut level_ids: std::collections::HashMap<String, usize> = std::collections::HashMap::new();
    let mut next = 0_usize;
    let mut a_groups = Vec::new();
    let mut b_groups = Vec::new();
    for s in &cfg.samples {
        let level = if re_col.eq_ignore_ascii_case("sample") {
            s.id.clone()
        } else if re_col.eq_ignore_ascii_case("group") {
            s.group.clone()
        } else {
            match s.covariates.get(re_col.as_str()) {
                Some(v) => v.clone(),
                None => continue,
            }
        };
        let id = *level_ids.entry(level).or_insert_with(|| {
            let v = next;
            next += 1;
            v
        });
        if s.group == contrast.numerator {
            a_groups.push(id);
        } else if s.group == contrast.denominator {
            b_groups.push(id);
        }
    }
    a_groups.extend(b_groups);
    a_groups
}

fn total_inc(xs: &[f64]) -> u64 {
    xs.iter().map(|x| x.max(0.0).round() as u64).sum()
}

fn inc_n(inc: &[f64], exc: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let k: Vec<f64> = inc.iter().map(|x| x.max(0.0)).collect();
    let n: Vec<f64> = inc
        .iter()
        .zip(exc.iter())
        .map(|(i, e)| (i + e).max(0.0))
        .collect();
    (k, n)
}

/// Beta-binomial log-likelihood for observations (k_i, n_i) given (α, β).
/// The binomial-coefficient term is dropped because it's constant in (α, β).
fn bb_log_lik(alpha: f64, beta: f64, k: &[f64], n: &[f64]) -> f64 {
    if alpha <= 0.0 || beta <= 0.0 {
        return f64::NEG_INFINITY;
    }
    let mut ll = 0.0;
    for (&k_i, &n_i) in k.iter().zip(n.iter()) {
        if n_i <= 0.0 {
            continue;
        }
        ll += ln_gamma(k_i + alpha) + ln_gamma(n_i - k_i + beta)
            - ln_gamma(n_i + alpha + beta)
            + ln_gamma(alpha + beta)
            - ln_gamma(alpha)
            - ln_gamma(beta);
    }
    ll
}

/// Fit (α, β) by Nelder-Mead on the log-likelihood surface. Returns
/// `(α̂, β̂, ll̂)`.
fn fit_beta_binomial(k: &[f64], n: &[f64]) -> (f64, f64, f64) {
    // Method-of-moments initial guess gives Nelder-Mead a strong start.
    let p_hat = k.iter().sum::<f64>() / n.iter().sum::<f64>().max(1.0);
    let p_clamped = p_hat.clamp(1e-3, 1.0 - 1e-3);
    let init = (p_clamped * 10.0, (1.0 - p_clamped) * 10.0);

    let f = |x: (f64, f64)| -> f64 {
        // Optimize on log-scale to keep α, β > 0.
        let (la, lb) = x;
        -bb_log_lik(la.exp(), lb.exp(), k, n)
    };

    let x0 = (init.0.ln(), init.1.ln());
    let (la, lb) = nelder_mead_2d(f, x0, 1e-6, 200);
    let (a, b) = (la.exp(), lb.exp());
    let ll = bb_log_lik(a, b, k, n);
    (a, b, ll)
}

/// Plain 2-D Nelder-Mead. Returns the best simplex vertex after at most
/// `max_iter` iterations or when the simplex contracts below `tol`.
fn nelder_mead_2d<F: Fn((f64, f64)) -> f64>(
    f: F,
    x0: (f64, f64),
    tol: f64,
    max_iter: usize,
) -> (f64, f64) {
    let mut simplex = [x0, (x0.0 + 0.5, x0.1), (x0.0, x0.1 + 0.5)];
    let mut vals = [f(simplex[0]), f(simplex[1]), f(simplex[2])];

    for _ in 0..max_iter {
        // Order by value ascending so simplex[0] is best.
        let mut idx = [0, 1, 2];
        idx.sort_by(|&i, &j| vals[i].partial_cmp(&vals[j]).unwrap_or(std::cmp::Ordering::Equal));
        let (b, m, w) = (idx[0], idx[1], idx[2]);
        if (simplex[w].0 - simplex[b].0).hypot(simplex[w].1 - simplex[b].1) < tol {
            return simplex[b];
        }
        // Centroid of best two.
        let c = ((simplex[b].0 + simplex[m].0) / 2.0, (simplex[b].1 + simplex[m].1) / 2.0);
        // Reflect worst through centroid.
        let r = (2.0 * c.0 - simplex[w].0, 2.0 * c.1 - simplex[w].1);
        let fr = f(r);
        if fr < vals[b] {
            // Expand.
            let e = (3.0 * c.0 - 2.0 * simplex[w].0, 3.0 * c.1 - 2.0 * simplex[w].1);
            let fe = f(e);
            if fe < fr {
                simplex[w] = e;
                vals[w] = fe;
            } else {
                simplex[w] = r;
                vals[w] = fr;
            }
        } else if fr < vals[m] {
            simplex[w] = r;
            vals[w] = fr;
        } else {
            // Contract.
            let cc = ((c.0 + simplex[w].0) / 2.0, (c.1 + simplex[w].1) / 2.0);
            let fc = f(cc);
            if fc < vals[w] {
                simplex[w] = cc;
                vals[w] = fc;
            } else {
                // Shrink toward best.
                for &i in &[m, w] {
                    simplex[i] = (
                        (simplex[b].0 + simplex[i].0) / 2.0,
                        (simplex[b].1 + simplex[i].1) / 2.0,
                    );
                    vals[i] = f(simplex[i]);
                }
            }
        }
    }
    let mut best = 0;
    for i in 1..3 {
        if vals[i] < vals[best] {
            best = i;
        }
    }
    simplex[best]
}

/// Two-sided Fisher's exact test on a 2 × 2 contingency table:
///
/// |          | inclusion | exclusion |
/// |----------|-----------|-----------|
/// | group A  |    a      |    b      |
/// | group B  |    c      |    d      |
///
/// Uses the hypergeometric distribution: under independence, the count `a`
/// is Hypergeometric(N = a+b+c+d, K = a+c, n = a+b). The p-value is the
/// sum of probabilities of all outcomes at least as extreme as observed.
pub fn fisher_exact_2x2(a: u64, b: u64, c: u64, d: u64) -> f64 {
    let n = a + b + c + d;
    if n == 0 {
        return 1.0;
    }
    let col1 = a + c;
    let row1 = a + b;
    let dist = match Hypergeometric::new(n, col1, row1) {
        Ok(d) => d,
        Err(_) => return 1.0,
    };
    let observed = dist.pmf(a);
    let lo = row1.saturating_sub(n - col1);
    let hi = row1.min(col1);
    // Two-sided p (minimum-likelihood method, as in R's `fisher.test`): sum
    // the probability of every table at least as extreme as observed. A
    // *relative* tolerance keeps the comparison numerically stable so the
    // result is invariant under row/column swaps and transposition — an
    // absolute epsilon (e.g. 1e-15) is far too tight for pmf values of order
    // 1e-2 and makes a borderline table flip in or out asymmetrically.
    let threshold = observed * (1.0 + 1e-7);
    let mut p = 0.0;
    for k in lo..=hi {
        let pmf = dist.pmf(k);
        if pmf <= threshold {
            p += pmf;
        }
    }
    p.clamp(0.0, 1.0)
}

/// Benjamini-Hochberg FDR adjustment. NaN p-values are preserved as NaN in
/// the output.
pub fn bh_fdr(p_values: &[f64]) -> Vec<f64> {
    let n = p_values.len();
    if n == 0 {
        return Vec::new();
    }
    // Pair (p, original index), keeping NaNs aside.
    let mut indexed: Vec<(usize, f64)> = p_values
        .iter()
        .enumerate()
        .filter(|(_, p)| p.is_finite())
        .map(|(i, p)| (i, *p))
        .collect();
    indexed.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));

    // BH adjusted p-value at rank i (1-based, ascending) is
    //   min over j >= i of (p_j * m / j), capped at 1.
    let m_total = indexed.len();
    let m = m_total as f64;
    let mut adjusted = vec![f64::NAN; n];
    let mut current_min = 1.0_f64;
    for i in (0..m_total).rev() {
        let (orig_idx, p) = indexed[i];
        let rank = (i + 1) as f64;
        let q = (p * m / rank).min(1.0);
        if q < current_min {
            current_min = q;
        }
        adjusted[orig_idx] = current_min;
    }
    adjusted
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bh_known_values() {
        // R: p.adjust(c(0.001, 0.008, 0.039, 0.041, 0.042, 0.060, 0.074, 0.205), method="BH")
        // -> 0.00800 0.03200 0.06720 0.06720 0.06720 0.08000 0.08457 0.20500
        // (verified against R 4.4 base::p.adjust)
        let p = [0.001, 0.008, 0.039, 0.041, 0.042, 0.060, 0.074, 0.205];
        let adj = bh_fdr(&p);
        let expected = [0.00800, 0.03200, 0.06720, 0.06720, 0.06720, 0.08000, 0.08457, 0.20500];
        for (a, e) in adj.iter().zip(expected.iter()) {
            assert!((a - e).abs() < 1e-3, "got {a}, expected {e}");
        }
    }

    #[test]
    fn fisher_no_difference() {
        let p = fisher_exact_2x2(50, 50, 50, 50);
        assert!(p > 0.9, "should be ~1.0, got {p}");
    }

    #[test]
    fn fisher_extreme_difference() {
        let p = fisher_exact_2x2(100, 0, 0, 100);
        assert!(p < 1e-10, "should be tiny, got {p}");
    }

    #[test]
    fn bb_lrt_no_difference() {
        // Same generative parameters in both groups → large p-value.
        let inc_a = vec![50.0, 55.0, 48.0];
        let exc_a = vec![50.0, 45.0, 52.0];
        let inc_b = vec![49.0, 51.0, 53.0];
        let exc_b = vec![51.0, 49.0, 47.0];
        let (p, used) = test_one(&inc_a, &exc_a, &inc_b, &exc_b);
        assert_eq!(used, TestUsed::BetaBinomialLRT);
        assert!(p > 0.05, "no difference should give p > 0.05, got {p}");
    }
}
