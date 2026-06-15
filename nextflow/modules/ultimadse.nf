process ULTIMADSE {
    tag "ultimaDSEcaller"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    val  sheet_rows
    path bams          // every BAM + .bai, staged into the work dir
    path gtf
    path fasta
    val  contrast

    output:
    path "ultimadse_results/**", emit: results
    path "samples.tsv"

    script:
    // Pick a tech flag: if any long-read platform was present the caller is
    // still fine in short mode for junction extraction, but we expose the
    // knob through params.tech_extra for explicit control.
    def header = "sample\\tbam\\tgroup"
    def rows   = sheet_rows.join('\\n')
    """
    printf '${header}\\n${rows}\\n' > samples.tsv

    ultimaDSEcaller run \\
        --sample-sheet samples.tsv \\
        --annotation ${gtf} \\
        --reference ${fasta} \\
        --out ultimadse_results \\
        --contrast '${contrast}' \\
        --consensus stouffer \\
        --threads ${task.cpus} \\
        ${params.tech_extra}
    """
}
