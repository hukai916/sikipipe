process READ_LENGTH_DIST_WITHIN_UMI_GROUP {
    tag "$meta.id"
    label 'process_low'

    container "hukai916/miniconda3_bio:0.3"

    input:
    tuple val(meta), path(csv)
    val outdir

    output:
    path "*/UMI_3/stat/*.csv",        emit: umi_3_stat
    path "*/UMI_3/plot/*",            emit: umi_3_plot
    path "*/UMI_10/stat/*.csv",       emit: umi_10_stat
    path "*/UMI_10/plot/*",           emit: umi_10_plot
    path "*/UMI_30/stat/*.csv",       emit: umi_30_stat
    path "*/UMI_30/plot/*",           emit: umi_30_plot
    path  "versions.yml",             emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    # examine repeat length distribution for UMI group 3, 10, and 30, each select 10 UMI to plot

    mkdir -p ${outdir}/UMI_3/stat ${outdir}/UMI_3/plot ${outdir}/UMI_10/stat ${outdir}/UMI_10/plot ${outdir}/UMI_30/stat ${outdir}/UMI_30/plot

    read_length_dist_within_umi_group.py $csv $prefix 3 ${outdir}/UMI_3/stat ${outdir}/UMI_3/plot $args
    read_length_dist_within_umi_group.py $csv $prefix 10 ${outdir}/UMI_10/stat ${outdir}/UMI_10/plot $args
    read_length_dist_within_umi_group.py $csv $prefix 30 ${outdir}/UMI_30/stat ${outdir}/UMI_30/plot $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS

    """
}
