process BAM_COMBINE {
    label 'process_low'

    container "hukai916/miniconda3_bio:0.3"

    input:
    path bam

    output:
    tuple val([id:"id_token"]), path("combined/combined.bam"),     emit: bam

    script:

    """
    mkdir combined # in case some raw bam files are named combined
    samtools merge combined/combined.bam *.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS

    """
}
