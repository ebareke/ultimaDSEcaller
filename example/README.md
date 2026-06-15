# Example — end-to-end in seconds

Two self-contained demonstrations of ultimaDSEcaller on a tiny synthetic
gene with a single cassette-exon (SE) event that is **included in `control`**
and **skipped in `treatment`**.

## 1. BAM → events (no aligner required)

Needs `samtools` and `python3` plus a built binary
(`cargo build --release`).

```bash
bash example/run_example.sh
```

It will:

1. generate `synthetic/reference.fa`, `annotation.gtf`, `samples.tsv` and
   four per-sample alignments (`scripts/make_synthetic.py`);
2. convert them to sorted, indexed BAMs with samtools;
3. run `ultimaDSEcaller run … --contrast group:treatment-control`;
4. print the high-confidence events.

Expected result: one significant **SE** event with a large negative ΔPSI
(treatment skips the exon), e.g.

```
event_type  delta_psi  adjusted_p_value  test_used
SE          -0.71      2.1e-04           BB-LRT
```

## 2. FASTQ → events (Nextflow)

Needs Nextflow and a container engine (Docker / Apptainer); the aligners
are bundled in the image.

```bash
python3 example/scripts/make_synthetic.py    # reference + GTF
python3 example/nextflow/make_fastq.py       # spliced-cDNA FASTQ + samplesheet
nextflow run nextflow/main.nf -profile test,docker
```

The `make_fastq.py` generator emits genuine spliced cDNA reads from the two
transcripts, so `minimap2 -ax splice` recovers the junctions and the
pipeline reproduces the same SE call starting from raw reads.

## Files

```
scripts/make_synthetic.py    reference + GTF + per-sample SAM (BAM path)
run_example.sh               BAM → events driver
nextflow/make_fastq.py       spliced-cDNA FASTQ + samplesheet (FASTQ path)
```

Generated artifacts (BAMs, FASTQ, results) are git-ignored — the scripts
reproduce them deterministically.
