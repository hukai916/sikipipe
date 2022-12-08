process BAM_COMBINE {
    tag "$meta.id"
    label 'process_low'

    container "hukai916/miniconda3_bio:0.3"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("combined/combined.bam"),     emit: bam

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
