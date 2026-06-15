#!/usr/bin/env nextflow
/*
 * ultimaDSEcaller — Nextflow pipeline
 * FASTQ -> (align) -> sorted BAM -> differential splicing events
 *
 * Short-read (Illumina) samples are aligned with STAR; long-read
 * (PacBio / ONT) samples with minimap2 -ax splice. Aligned BAMs are then
 * handed to ultimaDSEcaller for event detection, quantification, testing
 * and reporting.
 *
 * Usage:
 *   nextflow run main.nf \
 *     --samplesheet samples.csv \
 *     --gtf annotation.gtf \
 *     --fasta genome.fa \
 *     --contrast 'group:treatment-control' \
 *     --outdir results \
 *     -profile docker
 *
 * Samplesheet (CSV) columns:
 *   sample,group,platform,fastq_1,fastq_2
 * where platform is one of: illumina | pacbio | ont
 * and fastq_2 is empty for single-end / long-read data.
 */

nextflow.enable.dsl = 2

include { STAR_INDEX     } from './modules/star_index.nf'
include { STAR_ALIGN     } from './modules/star_align.nf'
include { MINIMAP2_ALIGN } from './modules/minimap2_align.nf'
include { ULTIMADSE      } from './modules/ultimadse.nf'

workflow {

    if (!params.samplesheet) { error "Please provide --samplesheet" }
    if (!params.gtf)         { error "Please provide --gtf" }
    if (!params.fasta)       { error "Please provide --fasta" }

    gtf   = file(params.gtf,   checkIfExists: true)
    fasta = file(params.fasta, checkIfExists: true)

    // Parse the samplesheet into a channel of per-sample maps + reads.
    reads = Channel
        .fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def meta = [ id: row.sample, group: row.group, platform: row.platform.toLowerCase() ]
            def r1 = file(row.fastq_1, checkIfExists: true)
            def r2 = (row.fastq_2 && row.fastq_2.trim()) ? file(row.fastq_2, checkIfExists: true) : null
            tuple(meta, r1, r2)
        }

    // Split by platform.
    short_reads = reads.filter { meta, r1, r2 -> meta.platform == 'illumina' }
    long_reads  = reads.filter { meta, r1, r2 -> meta.platform in ['pacbio', 'ont'] }

    // --- Short-read path: STAR ---
    star_index = STAR_INDEX(fasta, gtf)
    star_bams  = STAR_ALIGN(short_reads, star_index.collect(), gtf)

    // --- Long-read path: minimap2 ---
    mm2_bams = MINIMAP2_ALIGN(long_reads, fasta)

    // Merge aligned BAMs from both paths.
    all_bams = star_bams.bam.mix(mm2_bams.bam)

    // Build the ultimaDSEcaller sample sheet (sample<TAB>bam<TAB>group) and
    // collect every BAM + index for the single caller invocation.
    sheet_rows = all_bams.map { meta, bam, bai -> "${meta.id}\t${bam.name}\t${meta.group}" }
    bam_files  = all_bams.map { meta, bam, bai -> [bam, bai] }.flatten().collect()

    ULTIMADSE(
        sheet_rows.collect(),
        bam_files,
        gtf,
        fasta,
        params.contrast
    )
}

workflow.onComplete {
    log.info ( workflow.success
        ? "\n[ultimaDSEcaller] Done. Results in: ${params.outdir}\n"
        : "\n[ultimaDSEcaller] Pipeline failed. See .nextflow.log\n" )
}
