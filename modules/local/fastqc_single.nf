process FASTQC_SINGLE {
    tag "$meta.id"
    label 'process_medium'

    // conda (params.enable_conda ? "bioconda::fastqc=0.11.9" : null)
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/fastqc:0.11.9--0' :
    //     'quay.io/biocontainers/fastqc:0.11.9--0' }"
    container "hukai916/fastqc_0.11.9:0.1"

    input:
    tuple val(meta), path(reads)
    val outdir

    output:
    tuple val(meta), path("*/*.html"), emit: html
    tuple val(meta), path("*/*.zip") , emit: zip
    path  "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // Add soft-links to original FastQs for consistent naming in pipeline
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (meta) {
        """
        [ ! -f  ${prefix}.fastq.gz ] && ln -s $reads ${prefix}.fastq.gz
        mkdir $outdir
        fastqc $args --threads $task.cpus ${prefix}.fastq.gz
        mv *.html *.zip $outdir/

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastqc: \$( fastqc --version | sed -e "s/FastQC v//g" )
        END_VERSIONS
        """
    }
}
