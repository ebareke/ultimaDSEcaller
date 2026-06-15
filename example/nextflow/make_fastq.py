#!/usr/bin/env python3
"""
Emit genuine spliced-cDNA FASTQ for the ultimaDSEcaller Nextflow `test`
profile.

Reads the synthetic reference (../synthetic/reference.fa) and builds the two
transcript sequences of the toy gene:

    T1 (inclusion) = exon1(100-200) + exon2(400-500) + exon3(700-800)
    T2 (skipping)  = exon1(100-200)            + exon3(700-800)

It then writes per-sample single-end FASTQ where `control` samples are
inclusion-heavy and `treatment` samples are skipping-heavy. Because these are
real cDNA sequences, `minimap2 -ax splice` against the genome recovers the
exon-skipping junctions — making the Nextflow test a true FASTQ → events run.
"""

import gzip
import os
import random

random.seed(11)

HERE = os.path.dirname(os.path.abspath(__file__))
SYN = os.path.normpath(os.path.join(HERE, "..", "synthetic"))
FQDIR = os.path.join(HERE, "fastq")
os.makedirs(FQDIR, exist_ok=True)

# 1-based inclusive exon coordinates (must match make_synthetic.py).
EXON1 = (100, 200)
EXON2 = (400, 500)
EXON3 = (700, 800)


def read_reference(path):
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].strip().split()[0]
            else:
                seq.append(line.strip())
    return name, "".join(seq)


def sub(seq, a, b):
    return seq[a - 1 : b]  # 1-based inclusive


def main():
    ref_path = os.path.join(SYN, "reference.fa")
    if not os.path.exists(ref_path):
        raise SystemExit(
            f"{ref_path} not found — run ../scripts/make_synthetic.py first."
        )
    _, genome = read_reference(ref_path)

    e1 = sub(genome, *EXON1)
    e2 = sub(genome, *EXON2)
    e3 = sub(genome, *EXON3)
    t1 = e1 + e2 + e3  # inclusion
    t2 = e1 + e3       # skipping

    samples = [
        ("control_rep1", 40, 5, "control"),
        ("control_rep2", 38, 6, "control"),
        ("treatment_rep1", 6, 40, "treatment"),
        ("treatment_rep2", 5, 44, "treatment"),
    ]

    sheet = ["sample,group,platform,fastq_1,fastq_2"]
    qual = "I"

    for name, n_inc, n_skip, group in samples:
        recs = []
        for i in range(n_inc):
            recs.append((f"{name}_T1_{i}", t1))
        for i in range(n_skip):
            recs.append((f"{name}_T2_{i}", t2))
        random.shuffle(recs)
        fq_path = os.path.join(FQDIR, f"{name}.fastq.gz")
        with gzip.open(fq_path, "wt") as fh:
            for rid, seq in recs:
                fh.write(f"@{rid}\n{seq}\n+\n{qual * len(seq)}\n")
        sheet.append(f"{name},{group},ont,fastq/{name}.fastq.gz,")

    with open(os.path.join(HERE, "samplesheet.csv"), "w") as fh:
        fh.write("\n".join(sheet) + "\n")

    print(f"Wrote FASTQ for {len(samples)} samples to {FQDIR}")
    print(f"Wrote samplesheet to {os.path.join(HERE, 'samplesheet.csv')}")


if __name__ == "__main__":
    main()
