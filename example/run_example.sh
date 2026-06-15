#!/usr/bin/env bash
#
# End-to-end example for ultimaDSEcaller — no aligner required.
#
# 1. Generates a tiny synthetic dataset (reference, GTF, 4 samples) whose
#    `control` samples include a cassette exon and whose `treatment`
#    samples skip it.
# 2. Converts the synthetic SAM files to sorted, indexed BAMs (samtools).
# 3. Runs the caller and prints the called exon-skipping event.
#
# For the full *FASTQ → BAM → events* path (real alignment with
# minimap2 / STAR), use the Nextflow pipeline in ../nextflow/ or the
# container images, which bundle the aligners.
#
# Requirements: samtools, python3, and the ultimaDSEcaller binary
# (built with `cargo build --release` or available on PATH).

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SYN="$HERE/synthetic"
OUT="$HERE/results"

# Locate the binary: prefer a release build, fall back to PATH.
BIN="${ULTIMADSE_BIN:-}"
if [[ -z "$BIN" ]]; then
  if [[ -x "$HERE/../target/release/ultimaDSEcaller" ]]; then
    BIN="$HERE/../target/release/ultimaDSEcaller"
  elif command -v ultimaDSEcaller >/dev/null 2>&1; then
    BIN="ultimaDSEcaller"
  else
    echo "ERROR: ultimaDSEcaller binary not found. Build it with 'cargo build --release'." >&2
    exit 1
  fi
fi

echo ">> [1/3] Generating synthetic dataset"
python3 "$HERE/scripts/make_synthetic.py"

echo ">> [2/3] Converting SAM -> sorted, indexed BAM"
for sam in "$SYN"/*.sam; do
  base="$(basename "$sam" .sam)"
  samtools sort -O bam -o "$SYN/$base.sorted.bam" "$sam"
  samtools index "$SYN/$base.sorted.bam"
done

echo ">> [3/3] Running ultimaDSEcaller"
# The sample sheet lists BAMs by basename; run from the synthetic dir so the
# relative paths resolve.
cd "$SYN"
"$BIN" run \
  --sample-sheet samples.tsv \
  --annotation annotation.gtf \
  --reference reference.fa \
  --out "$OUT" \
  --contrast group:treatment-control \
  --tech short \
  --consensus stouffer \
  --min-coverage 5 \
  -v

echo
echo ">> Done. Key outputs:"
echo "   $OUT/events.raw.tsv"
echo "   $OUT/events.high_confidence.tsv"
echo "   $OUT/report.html"
echo
echo ">> High-confidence events:"
column -t -s$'\t' "$OUT/events.high_confidence.tsv" 2>/dev/null || cat "$OUT/events.high_confidence.tsv"
