process UMI_PATTERN {
    tag "$meta.id"
    label 'process_low'

    container "hukai916/miniconda3_bio:0.3"

    input:
    tuple val(meta), path(reads)
    val outdir

    output:
    path "*/stat/*.csv",    emit: stat
    path "*/plot/*",        emit: plot
    path  "versions.yml",   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    if (meta.single_end) {
      """
      mkdir -p ${outdir}/stat ${outdir}/plot

      zcat ${prefix}.fastq.gz | awk 'NR%4==1 {n = split(\$1, array, "_"); print array[n]}' | sort | uniq -c | awk '{print \$1}' | sort | uniq -c | awk '{print \$2","\$1}' | sort -n > tem.txt

      (echo -e "${prefix},count" && cat tem.txt | sort -n) > ${outdir}/stat/${prefix}_UMI_distribution.csv
      rm tem.txt

      plot_umi_pattern.py ${outdir}/stat/${prefix}_UMI_distribution.csv ${prefix} ${outdir}/plot/${prefix}_UMI_distribution.png $args

      cat <<-END_VERSIONS > versions.yml
      "${task.process}":
          python: \$( python --version | sed -e "s/python //g" )
      END_VERSIONS

      """

    } else {
      """
      mkdir -p ${outdir}/stat ${outdir}/plot

      zcat ${prefix}_1.fastq.gz | awk 'NR%4==1 {n = split(\$1, array, "_"); print array[n]}' | sort | uniq -c | awk '{print \$1}' | sort | uniq -c | awk '{print \$2","\$1}' | sort -n > tem.txt

      (echo -e "${prefix},count" && cat tem.txt | sort -n) > ${outdir}/stat/${prefix}_UMI_distribution.csv
      rm tem.txt

      plot_umi_pattern.py ${outdir}/stat/${prefix}_UMI_distribution.csv ${prefix} ${outdir}/plot/${prefix}_UMI_distribution.png $args

      cat <<-END_VERSIONS > versions.yml
      "${task.process}":
          python: \$( python --version | sed -e "s/python //g" )
      END_VERSIONS

      """
    }
}
