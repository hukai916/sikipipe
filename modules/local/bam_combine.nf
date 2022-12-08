process BAM_COMBINE {
    label 'process_low'

    container "hukai916/miniconda3_bio:0.3"

    input:
    path bam

    output:
    tuple val("token"), path("combined/combined.bam"),     emit: bam

    script:

    """
    mkdir combined # in case some raw bam files are named combined
    touch combined/test.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS

    """
}
