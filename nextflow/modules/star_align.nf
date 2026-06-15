process STAR_ALIGN {
    label 'align'
    tag "${meta.id}"
    publishDir "${params.outdir}/bam", mode: 'copy', pattern: "*.bam*"

    input:
    tuple val(meta), path(r1), path(r2)
    path index
    path gtf

    output:
    tuple val(meta), path("${meta.id}.bam"), path("${meta.id}.bam.bai"), emit: bam

    script:
    def reads = r2 ? "${r1} ${r2}" : "${r1}"
    def gzipped = r1.name.endsWith('.gz') ? '--readFilesCommand zcat' : ''
    """
    STAR --runMode alignReads \\
        --genomeDir ${index} \\
        --readFilesIn ${reads} \\
        ${gzipped} \\
        --runThreadN ${task.cpus} \\
        --outSAMtype BAM SortedByCoordinate \\
        --outSAMattributes NH HI AS nM \\
        --twopassMode Basic \\
        --outFileNamePrefix ${meta.id}.

    mv ${meta.id}.Aligned.sortedByCoord.out.bam ${meta.id}.bam
    samtools index -@ ${task.cpus} ${meta.id}.bam
    """
}
