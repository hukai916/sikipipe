process PREP_REF {
    label 'process_low'

    container "hukai916/miniconda3_bio:0.3"

    input:
    path ref
    val ref_rc
    val outdir

    output:
    path "*/ref/ref.fasta", emit: ref
    path  "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    if (ref_rc) {
      """
      mkdir -p $outdir/ref
      touch $outdir/ref/ref.fasta
      prep_ref.py $ref $outdir/ref/ref.fasta

      cat <<-END_VERSIONS > versions.yml
      "${task.process}":
          python: \$( python --version | sed -e "s/python //g" )
      END_VERSIONS

      """
    } else {

      """
      mkdir -p $outdir/ref
      cp -P $ref $outdir/ref/ref.fasta

      cat <<-END_VERSIONS > versions.yml
      "${task.process}":
          python: \$( python --version | sed -e "s/python //g" )
      END_VERSIONS

      """
    }


}
