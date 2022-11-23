process BLAT {
    tag "$meta.id"
    label 'process_low'

    container "hukai916/miniconda3_blat:0.2.1"

    input:
    tuple val(meta), path(reads)
    path ref
    val outdir

    output:
    path "*/blat_bam/*",                                           emit: blat_bam
    tuple val(meta), path("*/blat_bam_tuned/*.tuned.bam"),         emit: blat_bam_tuned
    path "*/fasta/*.fasta",                                        emit: fasta
    tuple val(meta), path("*/input/*.fastq.gz"),                   emit: reads
    tuple val(meta), path("*/mapped/*"),                           emit: mapped
    tuple val(meta), path("*/unmapped/*"),                         emit: unmapped
    path  "versions.yml",                                          emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    if (meta.single_end) {
      """
      mkdir -p ${outdir}/blat_bam ${outdir}/blat_bam_tuned ${outdir}/fasta ${outdir}/input
      cp -P $reads ${outdir}/input/

      # convert fastq to fasta:
      zcat < $reads | awk '{ if (NR%4==1) {sub(/^@/, ">", \$0); print \$0} else if (NR%4==2) {print \$0} }' > ${outdir}/fasta/${prefix}.fasta

      # blat to psl, convert to bam, add header:
      blat $ref ${outdir}/fasta/${prefix}.fasta ${outdir}/blat_bam/${prefix}.psl
      psl2sam.pl ${outdir}/blat_bam/${prefix}.psl > ${outdir}/blat_bam/${prefix}.sam
      ref_name=\$(awk 'NR==1 {print substr(\$0,2,length(\$0))}' $ref)
      # remove all trailing spaces in the reference name because samtools is very sensitive about it, any trailing spaces will invalidate the reference check, thereby, mapping the third column a * (some ref.fasta contains a special ending character)
      ref_name2=\$(echo \$ref_name | awk '{ gsub(/[[:space:]]+\$/,""); print }')
      echo -e "@HD\\tVN:1.0\\tGO:query" > tem.txt
      echo -e "@SQ\\tSN:\${ref_name2}\\tLN:X" >> tem.txt
      cat tem.txt ${outdir}/blat_bam/${prefix}.sam > ${outdir}/blat_bam/${prefix}.header.sam

      # sort:
      samtools sort ${outdir}/blat_bam/${prefix}.header.sam -o ${outdir}/blat_bam/${prefix}.bam
      samtools index ${outdir}/blat_bam/${prefix}.bam

      ## hard clip to soft clip:
      samtools view -h ${outdir}/blat_bam/${prefix}.bam | awk -F'\t' -v OFS='\t' 'NR>=4 {gsub(/H/, "S", \$6)} 1' | samtools view -o ${outdir}/blat_bam_tuned/${prefix}.soft.bam

      # tweak BAM to add seq and qual:

      ## version2: biopython
      tweak_blat_bam.py ${outdir}/blat_bam_tuned/${prefix}.soft.bam $reads ${outdir}/blat_bam_tuned/${prefix}.tuned.sam
      samtools view ${outdir}/blat_bam_tuned/${prefix}.tuned.sam -o ${outdir}/blat_bam_tuned/${prefix}.tuned.bam

      ## version1: can not reverse complement reads if flag contains 16, this leads to plotting issue for CrisprVariants because the seq portion in BAM should be reverse complemented if 16.
      ## add seq and qual:
      #zcat < $reads | awk -v ORS="" '{if (NR%4 == 0) {print \$0"\\n"} else {print \$0"\t"} }' > ${outdir}/blat_bam_tuned/${prefix}.fastq.tsv
      #samtools view -h ${outdir}/blat_bam_tuned/${prefix}.soft.bam | awk 'BEGIN{FS=OFS="\t"} NR==FNR { seq_dict[substr(\$1, 2)] = \$2; qual_dict[substr(\$1, 2)] = \$4; next } { if (FNR < 4) {print \$0} else { if (\$1 in seq_dict) { \$10 = seq_dict[\$1]; \$11 = qual_dict[\$1]; print \$0 } } } ' ${outdir}/blat_bam_tuned/${prefix}.fastq.tsv - | samtools view -o ${outdir}/blat_bam_tuned/${prefix}.tuned.bam

      ## seperate unmapped reads:
      mkdir ${outdir}/mapped
      mkdir ${outdir}/unmapped

      samtools view -f 4 ${outdir}/blat_bam_tuned/${prefix}.tuned.bam -o ${outdir}/unmapped/${prefix}.tuned.bam
      samtools view -F 4 ${outdir}/blat_bam_tuned/${prefix}.tuned.bam -o ${outdir}/mapped/${prefix}.tuned.bam

      samtools index ${outdir}/mapped/${prefix}.tuned.bam
      samtools index ${outdir}/unmapped/${prefix}.tuned.bam

      bedtools bamtofastq -i ${outdir}/unmapped/${prefix}.tuned.bam -fq ${outdir}/unmapped/${prefix}.tuned.fastq
      bedtools bamtofastq -i ${outdir}/mapped/${prefix}.tuned.bam -fq ${outdir}/mapped/${prefix}.tuned.fastq

      cat <<-END_VERSIONS > versions.yml
      "${task.process}":
          python: \$( python --version | sed -e "s/python //g" )
      END_VERSIONS

      """

    }
}
