#!/usr/bin/env python3
"""
Generate a tiny, self-contained synthetic dataset for ultimaDSEcaller.

It produces, under ``example/synthetic/``:

* ``reference.fa``  — one 2 kb contig ``chrT`` with GT-AG canonical splice
  motifs placed at the intron boundaries of a single test gene.
* ``annotation.gtf`` — one gene ``GENE1`` with two transcripts that differ
  by a single cassette exon (a textbook exon-skipping event).
* ``samples.tsv``   — a 4-sample sheet: 2 ``control`` + 2 ``treatment``.
* ``*.sam``         — per-sample alignments. ``control`` samples are
  inclusion-heavy, ``treatment`` samples are skipping-heavy, so the tool
  must call a significant exon-skipping (SE) event with a large ΔPSI.

The companion ``run_example.sh`` converts the SAM files to sorted, indexed
BAMs with samtools and then runs the caller — no aligner required, so the
example runs anywhere samtools is installed.

Gene model (genomic, 1-based inclusive):
    exon1 100-200   cassette 400-500   exon3 700-800
    inclusion junctions: (200->400) and (500->700)
    skipping  junction:  (200->700)
"""

import os
import random

random.seed(7)

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.normpath(os.path.join(HERE, "..", "synthetic"))
os.makedirs(OUT, exist_ok=True)

CONTIG = "chrT"
LENGTH = 2000

# ----------------------------------------------------------------------
# 1. Reference FASTA with canonical GT-AG motifs at every intron boundary.
# ----------------------------------------------------------------------
seq = list("".join(random.choice("ACGT") for _ in range(LENGTH)))


def place(pos1, bases):
    """Write `bases` starting at 1-based genomic position `pos1`."""
    for i, b in enumerate(bases):
        seq[pos1 - 1 + i] = b


# Intron 1: 201..399  (donor GT at 201, acceptor AG at 398..399)
place(201, "GT")
place(398, "AG")
# Intron 2: 501..699  (donor GT at 501, acceptor AG at 698..699)
place(501, "GT")
place(698, "AG")
# Skip intron: 201..699 shares the same donor (201 GT) and acceptor (698 AG).

with open(os.path.join(OUT, "reference.fa"), "w") as fh:
    fh.write(f">{CONTIG}\n")
    s = "".join(seq)
    for i in range(0, len(s), 60):
        fh.write(s[i : i + 60] + "\n")

# ----------------------------------------------------------------------
# 2. Annotation GTF — one gene, two transcripts.
# ----------------------------------------------------------------------
gtf_lines = []


def gtf(feature, start, end, attrs):
    gtf_lines.append(
        f"{CONTIG}\tsynthetic\t{feature}\t{start}\t{end}\t.\t+\t.\t{attrs}"
    )


# Transcript T1 — inclusion (exon1 + cassette + exon3)
for (s, e) in [(100, 200), (400, 500), (700, 800)]:
    gtf("exon", s, e, 'gene_id "GENE1"; transcript_id "GENE1.T1"; gene_name "TestGene";')
    gtf("CDS", s, e, 'gene_id "GENE1"; transcript_id "GENE1.T1"; gene_name "TestGene";')
# Transcript T2 — skipping (exon1 + exon3)
for (s, e) in [(100, 200), (700, 800)]:
    gtf("exon", s, e, 'gene_id "GENE1"; transcript_id "GENE1.T2"; gene_name "TestGene";')
    gtf("CDS", s, e, 'gene_id "GENE1"; transcript_id "GENE1.T2"; gene_name "TestGene";')

with open(os.path.join(OUT, "annotation.gtf"), "w") as fh:
    fh.write("\n".join(gtf_lines) + "\n")

# ----------------------------------------------------------------------
# 3. Synthetic alignments.
# ----------------------------------------------------------------------
READLEN_ANCHOR = 50  # bases of exon on each side of a junction


def sam_header():
    return [f"@HD\tVN:1.6\tSO:coordinate", f"@SQ\tSN:{CONTIG}\tLN:{LENGTH}"]


def spliced_read(qname, donor_end, intron_len, left_anchor=50, right_anchor=51):
    """Build a single spliced SAM record spanning one junction."""
    pos = donor_end - left_anchor + 1
    cigar = f"{left_anchor}M{intron_len}N{right_anchor}M"
    seqlen = left_anchor + right_anchor
    seqstr = "A" * seqlen
    qual = "I" * seqlen
    return f"{qname}\t0\t{CONTIG}\t{pos}\t60\t{cigar}\t*\t0\t0\t{seqstr}\t{qual}\tNH:i:1"


def write_sample(name, n_inclusion, n_skipping):
    lines = sam_header()
    rid = 0
    # Inclusion evidence: reads across junction1 (200->400) and junction2 (500->700)
    for _ in range(n_inclusion):
        lines.append(spliced_read(f"{name}_inc1_{rid}", 200, 199))  # 201..399
        rid += 1
        lines.append(spliced_read(f"{name}_inc2_{rid}", 500, 199))  # 501..699
        rid += 1
    # Skipping evidence: reads across the skip junction (200->700)
    for _ in range(n_skipping):
        lines.append(spliced_read(f"{name}_skip_{rid}", 200, 499))  # 201..699
        rid += 1
    path = os.path.join(OUT, f"{name}.sam")
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    return path


# control: inclusion-heavy ; treatment: skipping-heavy
samples = [
    ("control_rep1", 40, 4, "control"),
    ("control_rep2", 36, 6, "control"),
    ("treatment_rep1", 5, 38, "treatment"),
    ("treatment_rep2", 7, 42, "treatment"),
]

sheet = ["sample\tbam\tgroup"]
for name, n_inc, n_skip, group in samples:
    write_sample(name, n_inc, n_skip)
    sheet.append(f"{name}\t{name}.sorted.bam\t{group}")

with open(os.path.join(OUT, "samples.tsv"), "w") as fh:
    fh.write("\n".join(sheet) + "\n")

print(f"Wrote synthetic dataset to {OUT}")
print("  reference.fa, annotation.gtf, samples.tsv, and 4 *.sam files")
