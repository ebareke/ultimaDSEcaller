# Contributing

Issues and pull requests are welcome.

## Getting started

```bash
git clone https://github.com/ebareke/ultimaDSEcaller.git
cd ultimaDSEcaller
cargo build --release
cargo test --all-features
bash example/run_example.sh     # end-to-end sanity check (needs samtools, python3)
```

## Ground rules

- **Correctness first.** Statistical and genomic logic must be backed by a
  test. Bugs in this kind of tool are silent and propagate into other
  people's biology — new behaviour comes with a test that would have caught
  the bug.
- **No fake results.** A code path either computes the real quantity or
  returns a typed error / `NaN` — never a plausible-looking placeholder.
- **The core stays input-format-honest.** The binary consumes BAM/CRAM;
  alignment belongs in the Nextflow pipeline and containers, not in the Rust
  core.

## Where contributions help most

- **New event types and edge cases** in `src/events.rs` / `src/advanced.rs` —
  each detector should come with a unit test on a small synthetic GTF.
- **Statistical methods** in `src/stats.rs`, `src/glm.rs`, `src/consensus.rs`
  — validate against R (`fisher.test`, `p.adjust`, `glm`) and pin the
  expected values in a test.
- **Aligner / pipeline modules** in `nextflow/modules/` — keep them container-
  pinned and reproducible.
- **Real-data parser robustness** — Ensembl, Gencode, RefSeq and FlyBase
  GTF/GFF3 all have quirks; failing inputs make great regression tests.

## Pull-request checklist

- [ ] `cargo build --release --all-features` succeeds.
- [ ] `cargo test --all-features` passes; new behaviour has tests.
- [ ] `cargo fmt --all` and `cargo clippy --all-targets --all-features`
      are clean.
- [ ] `bash example/run_example.sh` still calls the expected event.
- [ ] User-facing changes are reflected in `USAGE.md` and `CHANGELOG.md`.

## Commit style

Short imperative subject lines ("Add tandem-UTR detector", "Fix Fisher
two-sided symmetry"). Reference issues where relevant.

## Maintainers

- Eric B. — <eb.bioinfo@pm.me>
- Ethan B. — <eb.bioinfo@pm.me>
- Conrad B. — <eb.bioinfo@pm.me>
