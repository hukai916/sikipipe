process CLASSIFY_INDEL {
    tag "$meta.id"
    label 'process_low'

    container "hukai916/miniconda3_bio:0.3"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("3a_classify_indel/no_indel/*.fastq.gz"),    emit: reads_no_indel
    tuple val(meta), path("3a_classify_indel/indel_5p/*.fastq.gz"),    emit: reads_indel_5p
    tuple val(meta), path("3a_classify_indel/indel_3p/*.fastq.gz"),    emit: reads_indel_3p
    tuple val(meta), path("3a_classify_indel/indel_5p_3p/*.fastq.gz"), emit: reads_indel_5p_3p
    path "3a_classify_indel/stat/*.csv",             emit: stat
    path  "versions.yml",                            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p 3a_classify_indel/no_indel 3a_classify_indel/indel_5p 3a_classify_indel/indel_3p 3a_classify_indel/indel_5p_3p 3a_classify_indel/stat

    classify_indel.py ${prefix}_1.fastq.gz ${prefix}_2.fastq.gz 3a_classify_indel/no_indel 3a_classify_indel/indel_5p 3a_classify_indel/indel_3p 3a_classify_indel/indel_5p_3p 3a_classify_indel/stat ${prefix} $args

    gzip 3a_classify_indel/no_indel/*
    gzip 3a_classify_indel/indel_5p/*
    gzip 3a_classify_indel/indel_3p/*
    gzip 3a_classify_indel/indel_5p_3p/*

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS

    """
}
