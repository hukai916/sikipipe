process CLASSIFY_READTHROUGH {
    tag "$meta.id"
    label 'process_low'

    container "hukai916/miniconda3_bio:0.3"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("4a_classify_readthrough/readthrough/*.fastq.gz"),        emit: reads_through
    tuple val(meta), path("4a_classify_readthrough/non_readthrough/*.fastq.gz"),    emit: reads_nonethrough
    path "4a_classify_readthrough/stat/*.csv",                    emit: stat
    path  "versions.yml",                                         emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p 4a_classify_readthrough/readthrough 4a_classify_readthrough/non_readthrough 4a_classify_readthrough/stat

    classify_readthrough.py ${prefix}_1.fastq.gz ${prefix}_2.fastq.gz 4a_classify_readthrough/readthrough 4a_classify_readthrough/non_readthrough 4a_classify_readthrough/stat ${prefix} $args

    gzip 4a_classify_readthrough/readthrough/*
    gzip 4a_classify_readthrough/non_readthrough/*

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS

    """
}
