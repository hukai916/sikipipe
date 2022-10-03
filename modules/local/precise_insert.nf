process PRECISE_INSERT {
    tag "$meta.id"
    label 'process_low'

    container "hukai916/miniconda3_bio:0.3"

    input:
    tuple val(meta), path(reads)
    tuple val(meta), path(bam)
    val outdir

    output:
    path "*/grep_precise_insert_fastq/*.fastq.gz",     emit: reads_grep_precise
    path "*/grep_precise_insert_cigar/*.csv",          emit: cigar_grep_precise
    path "*/stat/*.csv",                               emit: stat
    path  "versions.yml",                              emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    if (meta.single_end) {
      """
      mkdir -p ${outdir}/grep_precise_insert_fastq ${outdir}/grep_precise_insert_cigar ${outdir}/stat

      # grep precise reads:
      zcat < ${reads} | grep -A 2 -B 1 "$args" | grep -v "^--\$" > ${outdir}/grep_precise_insert_fastq/${prefix}.fastq || true # otherwise the job exit with 1 if empty which breaks the Nextflow pipe

      # grep precise read CIGAR from BAM:
      touch ${outdir}/grep_precise_insert_cigar/${prefix}.precise_insert_cigar.csv
      samtools view $bam | awk 'NR == FNR { cigar_dict["@"\$1] = \$6; next } { if (FNR % 4 == 1) { if (\$1 in cigar_dict) { print \$1","cigar_dict[\$1] } else { print \$1",""NA" } } }' - ${outdir}/grep_precise_insert_fastq/${prefix}.fastq > ${outdir}/grep_precise_insert_cigar/${prefix}.precise_insert_cigar.csv || true
      gzip ${outdir}/grep_precise_insert_fastq/${prefix}.fastq

      # stat: how many precise insertion get it to the BAM file:
      total=\$(cat ${outdir}/grep_precise_insert_cigar/${prefix}.precise_insert_cigar.csv | wc -l)
      cigar=\$(cat ${outdir}/grep_precise_insert_cigar/${prefix}.precise_insert_cigar.csv | grep -v ",NA" | wc -l) || true
      echo "${prefix},\$total,\$cigar" > ${outdir}/stat/${prefix}.precise_insert_cigar_ratio.csv

      cat <<-END_VERSIONS > versions.yml
      "${task.process}":
          python: \$( python --version | sed -e "s/python //g" )
      END_VERSIONS

      """

    }
}
