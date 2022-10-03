process BBMERGE {
    tag "$meta.id"
    label 'process_low'

    container "hukai916/bbmap:0.1"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("4b_bbmerge/merged/*.fastq.gz"),        emit: reads_merged
    tuple val(meta), path("4b_bbmerge/non_merged/*.fastq.gz"),    emit: reads_non_merged
    path "4b_bbmerge/stat/*.csv",               emit: stat
    path  "versions.yml",                       emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p 4b_bbmerge/merged 4b_bbmerge/non_merged 4b_bbmerge/stat

    bbmerge.sh in1=${prefix}_1.fastq.gz in2=${prefix}_2.fastq.gz out=4b_bbmerge/merged/${prefix}.fastq.gz outu1=4b_bbmerge/non_merged/${prefix}_1.fastq.gz outu2=4b_bbmerge/non_merged/${prefix}_2.fastq.gz

    c_merge=\$(echo \$(zcat 4b_bbmerge/merged/${prefix}.fastq.gz | wc -l)/4|bc)
    c_non_merge=\$(echo \$(zcat 4b_bbmerge/non_merged/${prefix}_1.fastq.gz | wc -l)/4|bc)

    echo "${prefix},\$c_merge,\$c_non_merge" > 4b_bbmerge/stat/${prefix}.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS

    """
}
