//! Benchmark suite — run with `cargo bench`. Each benchmark group
//! targets a hot path:
//!
//! * `bb_lrt`  — beta-binomial likelihood-ratio test
//! * `bh_fdr`  — Benjamini-Hochberg FDR
//! * `pca`     — PCA on the PSI matrix

use criterion::{Criterion, black_box, criterion_group, criterion_main};

fn bench_bb_lrt(c: &mut Criterion) {
    let inc_a = vec![50.0_f64, 55.0, 48.0, 52.0, 51.0];
    let exc_a = vec![50.0, 45.0, 52.0, 48.0, 49.0];
    let inc_b = vec![80.0, 78.0, 82.0, 81.0, 79.0];
    let exc_b = vec![20.0, 22.0, 18.0, 19.0, 21.0];
    c.bench_function("bb_lrt", |b| {
        b.iter(|| {
            ultimadsecaller::stats::test_one(
                black_box(&inc_a),
                black_box(&exc_a),
                black_box(&inc_b),
                black_box(&exc_b),
            )
        })
    });
}

fn bench_bh_fdr(c: &mut Criterion) {
    let p: Vec<f64> = (0..10_000).map(|i| (i as f64 + 1.0) / 20_001.0).collect();
    c.bench_function("bh_fdr_10k", |b| {
        b.iter(|| ultimadsecaller::stats::bh_fdr(black_box(&p)))
    });
}

fn bench_pca(c: &mut Criterion) {
    use nalgebra::DMatrix;
    let n_events = 500;
    let n_samples = 20;
    let m = DMatrix::from_fn(n_events, n_samples, |i, j| {
        let bias = if j < n_samples / 2 { 0.3 } else { 0.7 };
        ((i * 7 + j * 3) % 100) as f64 / 1000.0 + bias
    });
    let ids: Vec<String> = (0..n_samples).map(|i| format!("s{i}")).collect();
    c.bench_function("pca_500x20", |b| {
        b.iter(|| ultimadsecaller::embedding::pca_2d(black_box(&m), black_box(&ids)))
    });
}

criterion_group!(benches, bench_bb_lrt, bench_bh_fdr, bench_pca);
criterion_main!(benches);
