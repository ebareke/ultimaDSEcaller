process STAR_INDEX {
    label 'index'
    tag "STAR index"

    input:
    path fasta
    path gtf

    output:
    path "star_index"

    script:
    def overhang = params.sjdb_overhang ?: 100
    """
    mkdir -p star_index
    # Small-genome safe sjdbOverhang + genomeSAindexNbases auto-scaling.
    GENOME_LEN=\$(grep -v '^>' ${fasta} | tr -d '\\n' | wc -c)
    SA_INDEX=\$(awk -v L=\$GENOME_LEN 'BEGIN{n=int((log(L)/log(2))/2 - 1); if(n>14)n=14; if(n<4)n=4; print n}')
    STAR --runMode genomeGenerate \\
        --genomeDir star_index \\
        --genomeFastaFiles ${fasta} \\
        --sjdbGTFfile ${gtf} \\
        --sjdbOverhang ${overhang} \\
        --genomeSAindexNbases \$SA_INDEX \\
        --runThreadN ${task.cpus}
    """
}
