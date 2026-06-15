process MINIMAP2_ALIGN {
    label 'align'
    tag "${meta.id}"
    publishDir "${params.outdir}/bam", mode: 'copy', pattern: "*.bam*"

    input:
    tuple val(meta), path(r1), path(r2)
    path fasta

    output:
    tuple val(meta), path("${meta.id}.bam"), path("${meta.id}.bam.bai"), emit: bam

    script:
    // PacBio Iso-Seq and ONT cDNA/dRNA both use the spliced preset; the
    // strand flag differs (PacBio is typically forward-stranded).
    def preset = (meta.platform == 'pacbio') ? '-ax splice:hq -uf' : '-ax splice'
    """
    minimap2 ${preset} -t ${task.cpus} ${fasta} ${r1} \\
      | samtools sort -@ ${task.cpus} -o ${meta.id}.bam -
    samtools index -@ ${task.cpus} ${meta.id}.bam
    """
}
