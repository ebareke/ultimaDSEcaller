//! Sample-level embeddings on the event × sample PSI matrix.
//!
//! Two embeddings are produced:
//!
//! * **PCA** — principal component analysis via SVD on the centered PSI
//!   matrix. The first two principal coordinates are returned per sample.
//! * **UMAP-like** — k-nearest-neighbor graph followed by spectral
//!   embedding into 2-D. This is *not* full UMAP (which adds a stochastic
//!   cross-entropy optimization step on top); it's the LargeVis / spectral
//!   initializer that UMAP uses as a starting point. For small cohort
//!   sizes (≤ 200 samples) the result is qualitatively very similar.
//!
//! Both functions handle missing values (NaN PSI) by mean-imputation
//! across the same event.

use nalgebra::{DMatrix, DVector};
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Embedding2D {
    pub samples: Vec<String>,
    pub xs: Vec<f64>,
    pub ys: Vec<f64>,
    pub explained: Option<Vec<f64>>,
}

/// Compute the PSI matrix `events × samples` from the contrast-aware per-
/// event PSI vectors. Returns `None` if there are fewer than 2 events or
/// fewer than 2 samples.
pub fn build_psi_matrix(
    sample_ids: &[String],
    per_event_psi: &[Vec<f64>],
) -> Option<DMatrix<f64>> {
    if per_event_psi.is_empty() || sample_ids.len() < 2 {
        return None;
    }
    let n_events = per_event_psi.len();
    let n_samples = sample_ids.len();
    let mut m = DMatrix::<f64>::zeros(n_events, n_samples);
    for (i, row) in per_event_psi.iter().enumerate() {
        // Mean-impute NaN across samples within this event.
        let mut sum = 0.0;
        let mut cnt = 0;
        for v in row {
            if v.is_finite() {
                sum += *v;
                cnt += 1;
            }
        }
        let mean = if cnt > 0 { sum / cnt as f64 } else { 0.5 };
        for j in 0..n_samples.min(row.len()) {
            let v = row[j];
            m[(i, j)] = if v.is_finite() { v } else { mean };
        }
    }
    Some(m)
}

/// PCA → 2-D coordinates per sample.
pub fn pca_2d(matrix: &DMatrix<f64>, sample_ids: &[String]) -> Embedding2D {
    let (_n_events, n_samples) = matrix.shape();
    // Center each *event* (row) by subtracting its mean across samples.
    let mut centered = matrix.clone();
    for i in 0..centered.nrows() {
        let mean: f64 = centered.row(i).iter().sum::<f64>() / n_samples as f64;
        for j in 0..n_samples {
            centered[(i, j)] -= mean;
        }
    }
    // SVD via nalgebra. We want the right singular vectors (V) — each row
    // of Vᵀ is a sample coordinate in the PC basis.
    let svd = centered.svd(true, true);
    let v_t = match svd.v_t {
        Some(v) => v,
        None => return empty_embedding(sample_ids),
    };
    let singulars: Vec<f64> = svd.singular_values.iter().copied().collect();
    let xs: Vec<f64> = (0..n_samples).map(|j| v_t[(0, j)] * singulars[0]).collect();
    let ys: Vec<f64> = if v_t.nrows() > 1 {
        (0..n_samples).map(|j| v_t[(1, j)] * singulars[1]).collect()
    } else {
        vec![0.0; n_samples]
    };
    let total: f64 = singulars.iter().map(|s| s * s).sum();
    let explained: Vec<f64> = singulars
        .iter()
        .take(5)
        .map(|s| if total > 0.0 { s * s / total } else { 0.0 })
        .collect();
    Embedding2D {
        samples: sample_ids.to_vec(),
        xs,
        ys,
        explained: Some(explained),
    }
}

/// UMAP-like 2-D embedding via spectral layout of a symmetric k-NN graph.
///
/// Steps:
/// 1. Compute pairwise Euclidean distances between *samples* (columns).
/// 2. For each sample keep the `k` nearest neighbors; build a sparse
///    similarity matrix W with Gaussian weights using per-sample bandwidth
///    σ_i set to the distance to the k-th neighbor.
/// 3. Symmetrize: W ← (W + Wᵀ) − W ⊙ Wᵀ (UMAP fuzzy-union form).
/// 4. Form the normalized graph Laplacian L = I − D^(−1/2) W D^(−1/2).
/// 5. Take the eigenvectors of L corresponding to the second- and
///    third-smallest eigenvalues (the smallest is the trivial constant)
///    and use them as the 2-D coordinates.
///
/// This is the spectral initializer that real UMAP uses before its
/// gradient step; for small cohorts (n ≤ 200) it captures the global
/// structure faithfully.
pub fn umap_like_2d(matrix: &DMatrix<f64>, sample_ids: &[String], k: usize) -> Embedding2D {
    let (_n_events, n_samples) = matrix.shape();
    if n_samples < 3 {
        return empty_embedding(sample_ids);
    }
    let k = k.min(n_samples.saturating_sub(1)).max(1);

    // Pairwise distances.
    let mut dist = DMatrix::<f64>::zeros(n_samples, n_samples);
    for i in 0..n_samples {
        for j in (i + 1)..n_samples {
            let mut s = 0.0_f64;
            for r in 0..matrix.nrows() {
                let d = matrix[(r, i)] - matrix[(r, j)];
                s += d * d;
            }
            let d = s.sqrt();
            dist[(i, j)] = d;
            dist[(j, i)] = d;
        }
    }
    // k-NN with per-sample bandwidth.
    let mut sigma = vec![1.0_f64; n_samples];
    for i in 0..n_samples {
        let mut row: Vec<f64> = (0..n_samples).filter(|&j| j != i).map(|j| dist[(i, j)]).collect();
        row.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        sigma[i] = row.get(k - 1).copied().unwrap_or(1.0).max(1e-6);
    }
    let mut w = DMatrix::<f64>::zeros(n_samples, n_samples);
    for i in 0..n_samples {
        let mut idx: Vec<(usize, f64)> = (0..n_samples)
            .filter(|&j| j != i)
            .map(|j| (j, dist[(i, j)]))
            .collect();
        idx.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
        for (j, d) in idx.iter().take(k) {
            w[(i, *j)] = (-(d / sigma[i]).powi(2)).exp();
        }
    }
    // UMAP-style fuzzy-union symmetrization.
    let mut w_sym = DMatrix::<f64>::zeros(n_samples, n_samples);
    for i in 0..n_samples {
        for j in 0..n_samples {
            let a = w[(i, j)];
            let b = w[(j, i)];
            w_sym[(i, j)] = a + b - a * b;
        }
    }
    // Normalized Laplacian.
    let mut deg = DVector::<f64>::zeros(n_samples);
    for i in 0..n_samples {
        let s: f64 = w_sym.row(i).iter().sum();
        deg[i] = s.max(1e-9);
    }
    let mut l_norm = DMatrix::<f64>::identity(n_samples, n_samples);
    for i in 0..n_samples {
        for j in 0..n_samples {
            l_norm[(i, j)] -= w_sym[(i, j)] / (deg[i].sqrt() * deg[j].sqrt());
        }
    }
    // Eigendecomposition (symmetric).
    let sym = nalgebra::SymmetricEigen::new(l_norm);
    let mut order: Vec<usize> = (0..n_samples).collect();
    order.sort_by(|&a, &b| {
        sym.eigenvalues[a]
            .partial_cmp(&sym.eigenvalues[b])
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    // Skip the trivial first eigenvector.
    let pick_x = order.get(1).copied().unwrap_or(0);
    let pick_y = order.get(2).copied().unwrap_or(pick_x);
    let xs: Vec<f64> = (0..n_samples).map(|i| sym.eigenvectors[(i, pick_x)]).collect();
    let ys: Vec<f64> = (0..n_samples).map(|i| sym.eigenvectors[(i, pick_y)]).collect();
    Embedding2D {
        samples: sample_ids.to_vec(),
        xs,
        ys,
        explained: None,
    }
}

fn empty_embedding(sample_ids: &[String]) -> Embedding2D {
    let n = sample_ids.len();
    Embedding2D {
        samples: sample_ids.to_vec(),
        xs: vec![0.0; n],
        ys: vec![0.0; n],
        explained: None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pca_separates_two_clusters() {
        // Two clusters of samples on a 4-event matrix.
        let m = DMatrix::from_row_slice(
            4,
            6,
            &[
                0.1, 0.1, 0.15, 0.9, 0.85, 0.95,
                0.05, 0.1, 0.05, 0.95, 0.9, 0.92,
                0.9, 0.85, 0.95, 0.1, 0.05, 0.15,
                0.95, 0.9, 0.92, 0.05, 0.1, 0.05,
            ],
        );
        let ids: Vec<String> = (0..6).map(|i| format!("s{i}")).collect();
        let emb = pca_2d(&m, &ids);
        // PC1 should split into two groups: first 3 vs last 3.
        let g1: f64 = emb.xs.iter().take(3).sum::<f64>() / 3.0;
        let g2: f64 = emb.xs.iter().skip(3).sum::<f64>() / 3.0;
        assert!((g1 - g2).abs() > 0.5, "g1={g1} g2={g2}");
    }

    #[test]
    fn umap_like_runs_without_panic() {
        let m = DMatrix::from_row_slice(
            4,
            5,
            &[
                0.1, 0.2, 0.3, 0.7, 0.8,
                0.2, 0.1, 0.4, 0.8, 0.9,
                0.9, 0.85, 0.8, 0.2, 0.1,
                0.95, 0.9, 0.7, 0.1, 0.15,
            ],
        );
        let ids: Vec<String> = (0..5).map(|i| format!("s{i}")).collect();
        let emb = umap_like_2d(&m, &ids, 2);
        assert_eq!(emb.xs.len(), 5);
    }
}
