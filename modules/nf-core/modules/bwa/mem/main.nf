process BWA_MEM {
    // tag "$meta.id"
    label 'process_medium'

    container "hukai916/bwa_xenial:0.7.17"

    input:
    tuple val(meta), path(reads)
    path index
    val sort_bam
    path insert_fasta
    val outdir

    output:
    tuple val(meta), path("*/*.bam"),                   emit: bam
    tuple val(meta), path("*/mapped/insert/*.bam"),     emit: insert
    tuple val(meta), path("*/mapped/non_insert/*.bam"), emit: non_insert
    tuple val(meta), path("*/unmapped/*"),              emit: unmapped
    tuple val(meta), path("*/stat/*.csv"),              emit: stat
    path  "versions.yml",                               emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def samtools_command = sort_bam ? 'sort' : 'view'
    """
    mkdir $outdir
    mkdir ${outdir}/mapped
    mkdir ${outdir}/mapped/insert
    mkdir ${outdir}/mapped/non_insert
    mkdir ${outdir}/unmapped
    mkdir ${outdir}/stat

    INDEX=`find -L ./ -name "*.amb" | sed 's/.amb//'`

    bwa mem \\
        $args \\
        -t $task.cpus \\
        \$INDEX \\
        $reads \\
        | samtools $samtools_command $args2 --threads $task.cpus -o ${outdir}/${prefix}.bam -

    # seperate unmapped reads:
    samtools view -f 4 ${outdir}/${prefix}.bam -o ${outdir}/unmapped/${prefix}.bam
    samtools view -F 4 ${outdir}/${prefix}.bam -o ${outdir}/mapped/${prefix}.bam
    count_unmapped=\$(samtools view -c ${outdir}/unmapped/${prefix}.bam)
    count_mapped=\$(samtools view -c ${outdir}/mapped/${prefix}.bam)
    echo "${prefix},\$count_mapped,\$count_unmapped" > ${outdir}/stat/${prefix}.stat.csv

    bedtools bamtofastq -i ${outdir}/unmapped/${prefix}.bam -fq ${outdir}/unmapped/${prefix}.fastq
    bedtools bamtofastq -i ${outdir}/mapped/${prefix}.bam -fq ${outdir}/mapped/${prefix}.fastq

    # seperate mappable bam into insert and non_insert bam
    separate_insert.py ${outdir}/mapped/${prefix}.bam $insert_fasta ${outdir}/mapped/insert/${prefix}.bam ${outdir}/mapped/non_insert/${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
