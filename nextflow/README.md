# ultimaDSEcaller Nextflow pipeline

A complete **FASTQ → BAM → differential splicing events** workflow.

* **Short-read (Illumina)** samples are aligned with **STAR** (2-pass).
* **Long-read (PacBio / ONT)** samples are aligned with **minimap2 -ax splice**.
* Aligned BAMs are handed to **ultimaDSEcaller** for detection,
  quantification, statistical testing, and reporting.

All processes run inside the project container, so the only host
requirements are **Nextflow ≥ 23.04** and a container engine
(Docker / Singularity / Apptainer).

## Samplesheet

CSV with a header row:

```csv
sample,group,platform,fastq_1,fastq_2
ctrl1,control,illumina,reads/ctrl1_R1.fastq.gz,reads/ctrl1_R2.fastq.gz
ctrl2,control,illumina,reads/ctrl2_R1.fastq.gz,reads/ctrl2_R2.fastq.gz
trt1,treatment,ont,reads/trt1.fastq.gz,
trt2,treatment,pacbio,reads/trt2.fastq.gz,
```

* `platform` ∈ `illumina` | `pacbio` | `ont`.
* `fastq_2` is empty for single-end and long-read samples.
* You can mix short- and long-read samples in one run; each is routed to
  the appropriate aligner and merged before calling.

## Run

```bash
nextflow run nextflow/main.nf \
    --samplesheet samples.csv \
    --gtf  annotation.gtf \
    --fasta genome.fa \
    --contrast 'group:treatment-control' \
    --outdir results \
    -profile docker
```

Swap `-profile docker` for `-profile singularity` / `-profile apptainer`
on HPC.

## Built-in test

A tiny synthetic dataset is wired into the `test` profile. First generate
its inputs, then run:

```bash
python3 example/scripts/make_synthetic.py     # reference + GTF
python3 example/nextflow/make_fastq.py        # spliced-cDNA FASTQ + samplesheet
nextflow run nextflow/main.nf -profile test,docker
```

The test data contains a single cassette-exon (SE) event that is included
in `control` and skipped in `treatment`; the pipeline should report it as
a significant event with a large negative ΔPSI.

## Parameters

| Parameter        | Default                                  | Description                                   |
|------------------|------------------------------------------|-----------------------------------------------|
| `--samplesheet`  | —                                        | CSV described above (required)                |
| `--gtf`          | —                                        | GTF / GFF3 annotation (required)              |
| `--fasta`        | —                                        | Genome FASTA (required)                       |
| `--contrast`     | `group:treatment-control`                | Differential contrast                         |
| `--outdir`       | `results`                                | Output directory                              |
| `--tech_extra`   | `''`                                     | Extra flags forwarded to `ultimaDSEcaller run` |
| `--container`    | `ghcr.io/ebareke/ultimadsecaller:1.0.0`  | Container image used by every process         |

## Outputs

```
results/
├── bam/                       aligned, sorted, indexed BAMs
├── ultimadse_results/         full ultimaDSEcaller output (TSV/CSV/JSON/HTML)
├── samples.tsv                generated caller sample sheet
└── pipeline_info/             timeline, report, trace
```
