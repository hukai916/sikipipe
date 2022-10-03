process CAT_STAT {
    label 'process_low'

    container "hukai916/miniconda3_bio:0.3"

    input:
    path stat
    val outdir
    val header

    output:
    path "**/all_sample.csv", emit: stat
    path  "versions.yml",         emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    mkdir -p $outdir

    (echo -e "$header" && cat *.csv | sort -n) > $outdir/all_sample.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS

    """
}
