process READ_LENGTH_DIST {
    tag "$meta.id"
    label 'process_low'

    container "hukai916/miniconda3_bio:0.3"

    input:
    tuple val(meta), path(reads)
    val outdir

    output:
    tuple val(meta), path("*/length/*.csv"),       emit: length
    tuple val(meta), path("*/input/*.fastq.gz"),   emit: reads
    path  "versions.yml",                          emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p ${outdir}/count ${outdir}/input
    cp -P $reads ${outdir}/input/

    read_length_dist.py ${prefix}.fastq.gz ${prefix} ${outdir} $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS

    """
}
