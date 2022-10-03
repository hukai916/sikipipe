process BWA_INDEX {
    tag "$fasta"
    label 'process_single'

    container "hukai916/bwa_xenial:0.7.17"

    input:
    path fasta
    val outdir

    output:
    path "*/bwa",        emit: index
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mkdir -p ${outdir}/bwa
    bwa \\
        index \\
        $args \\
        -p ${outdir}/bwa/${fasta.baseName} \\
        $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
    END_VERSIONS
    """
}
