//! Property tests via `proptest`. Generative testing complements the
//! example-based tests by exploring a much wider input space and finding
//! counter-examples we wouldn't have thought to write.
//!
//! ## What we assert
//!
//! 1. **BH-FDR is monotone non-decreasing in original p-value rank**:
//!    sorting p-values ascending and applying BH must give a sequence
//!    where adjusted values never decrease.
//!
//! 2. **BH-FDR adjusted p-values are bounded** in `[0, 1]`.
//!
//! 3. **Fisher 2×2** is symmetric under row swap: swapping (a,b) and
//!    (c,d) must give the same p-value.
//!
//! 4. **Splice motif classification** is total — every valid 2-byte donor
//!    × 2-byte acceptor pair classifies to *some* variant (never panics).
//!
//! Hooks for miri runs (`cargo miri test --lib stats::tests`) are in
//! `.github/workflows/ci.yml`; valgrind is enabled by running the binary
//! under `valgrind ./target/release/ultimaDSEcaller …` and is documented in
//! the README's "Memory safety" subsection.

use proptest::prelude::*;

use ultimadsecaller::motif::SpliceMotif;
use ultimadsecaller::stats::{bh_fdr, fisher_exact_2x2};

proptest! {
    #[test]
    fn bh_is_monotone_non_decreasing(
        ps in prop::collection::vec(0.0f64..1.0, 1..200),
    ) {
        let adj = bh_fdr(&ps);
        // Sort indices by original p-value ascending.
        let mut idx: Vec<usize> = (0..ps.len()).collect();
        idx.sort_by(|&a, &b| ps[a].partial_cmp(&ps[b]).unwrap_or(std::cmp::Ordering::Equal));
        let mut prev = -f64::INFINITY;
        for i in idx {
            let q = adj[i];
            prop_assert!(q.is_finite(), "expected finite, got {q}");
            prop_assert!(q >= prev - 1e-12, "BH not monotone: {prev} -> {q}");
            prop_assert!((0.0..=1.0).contains(&q));
            prev = q;
        }
    }

    #[test]
    fn fisher_symmetric_under_row_swap(
        a in 0u64..200,
        b in 0u64..200,
        c in 0u64..200,
        d in 0u64..200,
    ) {
        let p1 = fisher_exact_2x2(a, b, c, d);
        let p2 = fisher_exact_2x2(c, d, a, b);
        prop_assert!((p1 - p2).abs() < 1e-9, "p1={p1} p2={p2}");
    }

    #[test]
    fn motif_classification_total(
        d in proptest::array::uniform2(0u8..=255u8),
        a in proptest::array::uniform2(0u8..=255u8),
    ) {
        // Should never panic for any 2-byte input.
        let _ = SpliceMotif::from_dinucleotides(&d, &a);
    }
}

#[test]
fn miri_smoke() {
    // Tiny example that runs cleanly under `cargo miri test` — used to catch
    // any unsafe-block regressions in the codebase.
    let p = [0.01, 0.05, 0.1, 0.5];
    let adj = bh_fdr(&p);
    assert_eq!(adj.len(), 4);
}
