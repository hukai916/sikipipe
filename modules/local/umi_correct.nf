process UMI_CORRECT {
    tag "$meta.id"
    label 'process_low'

    container "hukai916/miniconda3_bio:0.3"

    input:
    tuple val(meta), path(csv)
    tuple val(meta), path(reads)
    val outdir

    output:
    tuple val(meta), path("*/cutoff_1/*.fastq.gz"),   emit: reads_cutoff1
    tuple val(meta), path("*/cutoff_5/*.fastq.gz"),   emit: reads_cutoff5
    tuple val(meta), path("*/cutoff_30/*.fastq.gz"),  emit: reads_cutoff30
    path "*/stat/*.csv",                              emit: stat
    path  "versions.yml",                             emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p ${outdir}/cutoff_1 ${outdir}/cutoff_5 ${outdir}/cutoff_30 ${outdir}/stat

    umi_correct.py $csv $prefix $reads ${outdir}/cutoff_1/${prefix}.fastq 1 $args
    umi_correct.py $csv $prefix $reads ${outdir}/cutoff_5/${prefix}.fastq 5 $args
    umi_correct.py $csv $prefix $reads ${outdir}/cutoff_30/${prefix}.fastq 30 $args

    gzip ${outdir}/cutoff_1/${prefix}.fastq
    gzip ${outdir}/cutoff_5/${prefix}.fastq
    gzip ${outdir}/cutoff_30/${prefix}.fastq

    count1=\$(zcat < ${outdir}/cutoff_1/${prefix}.fastq.gz | awk 'NR%4==0{print 1}' | wc -l)
    count2=\$(zcat < ${outdir}/cutoff_5/${prefix}.fastq.gz | awk 'NR%4==0{print 1}' | wc -l)
    count3=\$(zcat < ${outdir}/cutoff_30/${prefix}.fastq.gz | awk 'NR%4==0{print 1}' | wc -l)

    echo "${prefix}",\$count1,\$count2,\$count3 > ${outdir}/stat/${prefix}.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed -e "s/python //g" )
    END_VERSIONS

    """
}
