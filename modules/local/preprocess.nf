process PREPROCESS {
    label 'process_low'

    container "hukai916/miniconda3_bio:0.3"

    input:
    tuple val(meta), path(reads)
    val outdir

    output:
    tuple val(meta), path("*/reads_normal/*.fastq.gz"),     emit: reads_normal
    tuple val(meta), path("*/reads_abnormal/*.fastq.gz"),   emit: reads_abnormal
    path "*/stat/*.csv",                                    emit: stat
    path "versions.yml",                                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def univ_5p = task.ext.univ_5p ?: ''
    def univ_3p = task.ext.univ_3p ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p $outdir/reads_normal $outdir/reads_abnormal $outdir/stat
    touch $outdir/stat/${prefix}.csv

    prep_data.py $reads $prefix ${outdir}/reads_normal/${prefix}.fastq.gz ${outdir}/reads_abnormal/${prefix}.fastq.gz $outdir/stat/${prefix}.csv $univ_5p $univ_3p 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS

    """
}
