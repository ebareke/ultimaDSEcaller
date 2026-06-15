//! Generalized linear models for splice-event differential testing.
//!
//! Two estimators are provided:
//!
//! * **GLM** — binomial family with logit link, fit by IRLS (iteratively
//!   reweighted least squares). This is the textbook logistic regression
//!   for an event's per-sample (inclusion, n) counts against an arbitrary
//!   design matrix. Suitable when samples are independent.
//!
//! * **GLMM** — the same model with one Gaussian random-intercept per
//!   grouping factor, fit by PQL (penalized quasi-likelihood). This is
//!   appropriate when samples are clustered (subjects, donors, time
//!   points) and you want to share strength across observations within a
//!   group rather than pretending they're independent.
//!
//! Both produce a fitted coefficient vector, its variance-covariance matrix,
//! and a log-likelihood. The Wald test on a contrast vector is a one-liner
//! given the vcov ([`wald_test`]).
//!
//! PQL is a *first-order* approximation that's known to underestimate
//! variance components in low-count regimes. It's what lme4 falls back to
//! when full Laplace fails, and it's good enough for the
//! splice-junction-count signal-to-noise typical of bulk RNA-seq. A future
//! upgrade to second-order Laplace or adaptive Gauss-Hermite quadrature
//! would land entirely within this file.

use nalgebra::{DMatrix, DVector};
use serde::{Deserialize, Serialize};
use statrs::distribution::{ContinuousCDF, Normal};

use crate::error::{UltiError, UltiResult};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GlmFit {
    pub beta: Vec<f64>,
    pub vcov: Vec<Vec<f64>>,
    pub loglik: f64,
    pub iterations: usize,
    pub converged: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GlmmFit {
    pub beta: Vec<f64>,
    pub vcov: Vec<Vec<f64>>,
    pub variance_component: f64,
    pub random_effects: Vec<f64>,
    pub loglik: f64,
    pub iterations: usize,
    pub converged: bool,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct WaldResult {
    pub estimate: f64,
    pub se: f64,
    pub z: f64,
    pub p_value: f64,
}

/// Fit a binomial GLM with logit link.
///
/// * `y` — successes per observation (must be ≤ n[i])
/// * `n` — trials per observation
/// * `x` — n×p design matrix (caller is responsible for including a column
///   of ones if an intercept is wanted)
///
/// Returns an error if the design matrix is rank-deficient.
pub fn fit_glm(y: &[f64], n: &[f64], x: &[Vec<f64>]) -> UltiResult<GlmFit> {
    let n_obs = y.len();
    if n_obs == 0 || n.len() != n_obs || x.len() != n_obs {
        return Err(UltiError::Stats(
            "fit_glm: mismatched input lengths".into(),
        ));
    }
    let p = x[0].len();
    if p == 0 {
        return Err(UltiError::Stats("fit_glm: empty design matrix".into()));
    }

    let y_v = DVector::from_vec(y.to_vec());
    let n_v = DVector::from_vec(n.to_vec());
    let x_m = matrix_from_rows(x, p)?;

    // Starting values via Haldane-adjusted empirical logit.
    let mut beta = DVector::<f64>::zeros(p);
    {
        // Use intercept estimate = logit(mean p_hat); other betas = 0.
        let mut k_sum = 0.0;
        let mut n_sum = 0.0;
        for i in 0..n_obs {
            k_sum += y_v[i];
            n_sum += n_v[i];
        }
        if n_sum > 0.0 {
            let p_hat = ((k_sum + 0.5) / (n_sum + 1.0)).clamp(1e-6, 1.0 - 1e-6);
            beta[0] = (p_hat / (1.0 - p_hat)).ln();
        }
    }

    let mut iter = 0;
    let max_iter = 50;
    let tol = 1e-7;
    let mut converged = false;
    let mut vcov: DMatrix<f64> = DMatrix::zeros(p, p);

    while iter < max_iter {
        iter += 1;
        let eta: DVector<f64> = &x_m * &beta;
        let p_hat: DVector<f64> = eta.map(|e| 1.0 / (1.0 + (-e).exp()));
        // Working weights w_i = n_i * p_i * (1 - p_i)
        let w: DVector<f64> = DVector::from_iterator(
            n_obs,
            (0..n_obs).map(|i| (n_v[i] * p_hat[i] * (1.0 - p_hat[i])).max(1e-12)),
        );
        // Working response z_i = eta_i + (y_i - n_i p_i) / w_i
        let z: DVector<f64> = DVector::from_iterator(
            n_obs,
            (0..n_obs).map(|i| eta[i] + (y_v[i] - n_v[i] * p_hat[i]) / w[i]),
        );

        // β_new = (Xᵀ W X)^(-1) Xᵀ W z, solved as a weighted normal equation.
        let xtw_x = weighted_xtwx(&x_m, &w);
        let xtw_z = weighted_xtwz(&x_m, &w, &z);
        let lu = xtw_x.clone().lu();
        let beta_new = match lu.solve(&xtw_z) {
            Some(v) => v,
            None => {
                return Err(UltiError::Stats(
                    "fit_glm: design matrix appears singular".into(),
                ))
            }
        };
        // Store vcov = (XᵀWX)^(-1)
        match xtw_x.try_inverse() {
            Some(inv) => vcov = inv,
            None => {
                return Err(UltiError::Stats(
                    "fit_glm: cannot invert weighted normal matrix".into(),
                ))
            }
        }

        let diff = (&beta_new - &beta).norm();
        beta = beta_new;
        if diff < tol {
            converged = true;
            break;
        }
    }

    // Log-likelihood at MLE (binomial, ignoring the constant n choose k term).
    let mut loglik = 0.0;
    let eta: DVector<f64> = &x_m * &beta;
    for i in 0..n_obs {
        let p = 1.0 / (1.0 + (-eta[i]).exp());
        if p > 0.0 && p < 1.0 {
            loglik += y_v[i] * p.ln() + (n_v[i] - y_v[i]) * (1.0 - p).ln();
        }
    }

    Ok(GlmFit {
        beta: beta.iter().copied().collect(),
        vcov: matrix_to_vec(&vcov),
        loglik,
        iterations: iter,
        converged,
    })
}

/// Fit a binomial GLMM with one random intercept per level of `group_id`.
///
/// `group_id[i]` is the integer level (0..G) of observation i. The model
/// is `logit(p_i) = X_i β + u_{g(i)}`, `u_g ~ N(0, σ²)`.
///
/// Algorithm: PQL — alternates between
/// 1. Fix u, fit β by IRLS (a single weighted least squares step).
/// 2. Fix β, update u_g by ridge regression with penalty 1/σ².
/// 3. Re-estimate σ² from u via method of moments.
pub fn fit_glmm(
    y: &[f64],
    n: &[f64],
    x: &[Vec<f64>],
    group_id: &[usize],
) -> UltiResult<GlmmFit> {
    let n_obs = y.len();
    if n_obs == 0 || n.len() != n_obs || x.len() != n_obs || group_id.len() != n_obs {
        return Err(UltiError::Stats(
            "fit_glmm: mismatched input lengths".into(),
        ));
    }
    let p = x[0].len();
    let n_groups = group_id.iter().copied().max().unwrap_or(0) + 1;

    let y_v = DVector::from_vec(y.to_vec());
    let n_v = DVector::from_vec(n.to_vec());
    let x_m = matrix_from_rows(x, p)?;

    let mut beta = DVector::<f64>::zeros(p);
    let mut u = DVector::<f64>::zeros(n_groups);
    let mut sigma2 = 1.0_f64;

    let mut iter = 0;
    let max_iter = 30;
    let tol = 1e-5;
    let mut converged = false;
    let mut vcov: DMatrix<f64> = DMatrix::zeros(p, p);

    while iter < max_iter {
        iter += 1;
        // Linear predictor with current β and u.
        let mut eta = DVector::<f64>::zeros(n_obs);
        for i in 0..n_obs {
            let mut s = 0.0;
            for j in 0..p {
                s += x_m[(i, j)] * beta[j];
            }
            eta[i] = s + u[group_id[i]];
        }
        let p_hat = eta.map(|e| 1.0 / (1.0 + (-e).exp()));
        let w = DVector::from_iterator(
            n_obs,
            (0..n_obs).map(|i| (n_v[i] * p_hat[i] * (1.0 - p_hat[i])).max(1e-12)),
        );
        let z = DVector::from_iterator(
            n_obs,
            (0..n_obs).map(|i| eta[i] + (y_v[i] - n_v[i] * p_hat[i]) / w[i]),
        );

        // Step 1: β update — weighted normal equation against (z - u_offset).
        let mut z_minus_u = z.clone();
        for i in 0..n_obs {
            z_minus_u[i] -= u[group_id[i]];
        }
        let xtw_x = weighted_xtwx(&x_m, &w);
        let xtw_z = weighted_xtwz(&x_m, &w, &z_minus_u);
        let beta_new = match xtw_x.clone().lu().solve(&xtw_z) {
            Some(v) => v,
            None => {
                return Err(UltiError::Stats(
                    "fit_glmm: design matrix appears singular".into(),
                ))
            }
        };
        vcov = match xtw_x.try_inverse() {
            Some(inv) => inv,
            None => DMatrix::zeros(p, p),
        };

        // Step 2: u update — for each group, ridge-regress residual on intercept.
        let mut u_new = DVector::<f64>::zeros(n_groups);
        let mut group_w_sum = vec![0.0_f64; n_groups];
        let mut group_wr_sum = vec![0.0_f64; n_groups];
        for i in 0..n_obs {
            let mut xbeta = 0.0;
            for j in 0..p {
                xbeta += x_m[(i, j)] * beta_new[j];
            }
            let r = z[i] - xbeta;
            group_w_sum[group_id[i]] += w[i];
            group_wr_sum[group_id[i]] += w[i] * r;
        }
        let lambda = 1.0 / sigma2.max(1e-8);
        for g in 0..n_groups {
            u_new[g] = group_wr_sum[g] / (group_w_sum[g] + lambda);
        }

        // Step 3: σ² update by MoM.
        let sigma2_new = if n_groups >= 2 {
            let mean = u_new.iter().sum::<f64>() / n_groups as f64;
            let var = u_new.iter().map(|x| (x - mean).powi(2)).sum::<f64>()
                / (n_groups - 1) as f64;
            var.max(1e-6)
        } else {
            sigma2
        };

        let diff = (&beta_new - &beta).norm() + (&u_new - &u).norm()
            + (sigma2_new - sigma2).abs();
        beta = beta_new;
        u = u_new;
        sigma2 = sigma2_new;
        if diff < tol {
            converged = true;
            break;
        }
    }

    // Approximate log-likelihood — quasi-likelihood (no integration constant).
    let mut loglik = 0.0;
    for i in 0..n_obs {
        let mut s = u[group_id[i]];
        for j in 0..p {
            s += x_m[(i, j)] * beta[j];
        }
        let p_i = 1.0 / (1.0 + (-s).exp());
        if p_i > 0.0 && p_i < 1.0 {
            loglik += y_v[i] * p_i.ln() + (n_v[i] - y_v[i]) * (1.0 - p_i).ln();
        }
    }
    // Penalty term for random effects.
    let half_log_sigma = 0.5 * sigma2.ln();
    for ug in u.iter() {
        loglik -= 0.5 * ug * ug / sigma2.max(1e-8) + half_log_sigma;
    }

    Ok(GlmmFit {
        beta: beta.iter().copied().collect(),
        vcov: matrix_to_vec(&vcov),
        variance_component: sigma2,
        random_effects: u.iter().copied().collect(),
        loglik,
        iterations: iter,
        converged,
    })
}

/// Wald test for a linear contrast cᵀβ. Returns the estimate, SE, z, and
/// two-sided p-value.
pub fn wald_test(beta: &[f64], vcov: &[Vec<f64>], contrast: &[f64]) -> WaldResult {
    let p = beta.len();
    if vcov.len() != p || contrast.len() != p {
        return WaldResult {
            estimate: f64::NAN,
            se: f64::NAN,
            z: f64::NAN,
            p_value: f64::NAN,
        };
    }
    let estimate: f64 = beta.iter().zip(contrast.iter()).map(|(b, c)| b * c).sum();
    let mut var = 0.0_f64;
    for i in 0..p {
        for j in 0..p {
            var += contrast[i] * vcov[i][j] * contrast[j];
        }
    }
    let se = var.max(0.0).sqrt();
    let z = if se > 0.0 { estimate / se } else { f64::NAN };
    let p_value = if z.is_finite() {
        let normal = Normal::standard();
        2.0 * (1.0 - normal.cdf(z.abs()))
    } else {
        f64::NAN
    };
    WaldResult {
        estimate,
        se,
        z,
        p_value,
    }
}

fn matrix_from_rows(rows: &[Vec<f64>], p: usize) -> UltiResult<DMatrix<f64>> {
    let n = rows.len();
    let mut m = DMatrix::<f64>::zeros(n, p);
    for (i, r) in rows.iter().enumerate() {
        if r.len() != p {
            return Err(UltiError::Stats(format!(
                "design matrix row {i} has {} cols, expected {p}",
                r.len()
            )));
        }
        for j in 0..p {
            m[(i, j)] = r[j];
        }
    }
    Ok(m)
}

fn weighted_xtwx(x: &DMatrix<f64>, w: &DVector<f64>) -> DMatrix<f64> {
    let (n, p) = x.shape();
    let mut out = DMatrix::<f64>::zeros(p, p);
    for i in 0..n {
        let wi = w[i];
        for a in 0..p {
            let xa = x[(i, a)];
            for b in 0..p {
                out[(a, b)] += wi * xa * x[(i, b)];
            }
        }
    }
    out
}

fn weighted_xtwz(x: &DMatrix<f64>, w: &DVector<f64>, z: &DVector<f64>) -> DVector<f64> {
    let (n, p) = x.shape();
    let mut out = DVector::<f64>::zeros(p);
    for i in 0..n {
        let wi = w[i];
        for a in 0..p {
            out[a] += wi * x[(i, a)] * z[i];
        }
    }
    out
}

fn matrix_to_vec(m: &DMatrix<f64>) -> Vec<Vec<f64>> {
    let (r, c) = m.shape();
    let mut out = vec![vec![0.0; c]; r];
    for i in 0..r {
        for j in 0..c {
            out[i][j] = m[(i, j)];
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn glm_recovers_strong_signal() {
        // Two groups, very different success rates.
        // Group 0 (treatment indicator = 0): p ≈ 0.2
        // Group 1 (treatment indicator = 1): p ≈ 0.8
        let y = vec![20.0, 22.0, 18.0, 80.0, 78.0, 82.0];
        let n = vec![100.0; 6];
        let x = vec![
            vec![1.0, 0.0],
            vec![1.0, 0.0],
            vec![1.0, 0.0],
            vec![1.0, 1.0],
            vec![1.0, 1.0],
            vec![1.0, 1.0],
        ];
        let fit = fit_glm(&y, &n, &x).unwrap();
        assert!(fit.converged);
        // Treatment coefficient should be ~ logit(0.8) - logit(0.2) ≈ 2.77.
        assert!(
            (fit.beta[1] - 2.77).abs() < 0.3,
            "got coef {} expected ≈ 2.77",
            fit.beta[1]
        );
        let w = wald_test(&fit.beta, &fit.vcov, &[0.0, 1.0]);
        assert!(w.p_value < 1e-10, "p = {}", w.p_value);
    }

    #[test]
    fn glm_null_signal_returns_high_p() {
        let y = vec![50.0; 6];
        let n = vec![100.0; 6];
        let x = vec![
            vec![1.0, 0.0],
            vec![1.0, 0.0],
            vec![1.0, 0.0],
            vec![1.0, 1.0],
            vec![1.0, 1.0],
            vec![1.0, 1.0],
        ];
        let fit = fit_glm(&y, &n, &x).unwrap();
        let w = wald_test(&fit.beta, &fit.vcov, &[0.0, 1.0]);
        assert!(w.p_value > 0.05);
    }

    #[test]
    fn glmm_handles_grouped_data() {
        // 4 subjects, 2 samples each. Treatment effect is in column 1.
        let y = vec![20.0, 25.0, 80.0, 75.0, 22.0, 18.0, 78.0, 82.0];
        let n = vec![100.0; 8];
        let x = vec![
            vec![1.0, 0.0],
            vec![1.0, 0.0],
            vec![1.0, 1.0],
            vec![1.0, 1.0],
            vec![1.0, 0.0],
            vec![1.0, 0.0],
            vec![1.0, 1.0],
            vec![1.0, 1.0],
        ];
        let group = vec![0, 0, 0, 0, 1, 1, 1, 1];
        let fit = fit_glmm(&y, &n, &x, &group).unwrap();
        let w = wald_test(&fit.beta, &fit.vcov, &[0.0, 1.0]);
        assert!(w.p_value < 1e-5, "p = {}", w.p_value);
    }
}
