process UMI_GROUP_STAT {
    tag "$meta.id"
    label 'process_low'

    container "hukai916/miniconda3_bio:0.3"

    input:
    tuple val(meta), path(csv)
    tuple val(meta), path(reads)
    val outdir

    output:
    tuple val(meta), path("*/*.csv"),            emit: stat
    tuple val(meta), path("*/input/*.fastq.gz"), emit: reads
    path  "versions.yml",                        emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    # examine repeat length distribution for UMI group 3, 10, and 100, each select 10 UMI to plot

    mkdir -p ${outdir} ${outdir}/input/
    cp -P $reads ${outdir}/input/

    umi_group_stat.py $csv $prefix ${outdir} $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS

    """
}
