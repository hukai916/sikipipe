process CRISPRVARIANTS {
    tag "$meta.id"
    label 'process_low'

    container "hukai916/crisprvariants:0.1"

    input:
    tuple val(meta), path(bam)
    path ref
    path guide_bed
    val outdir

    output:
    path "*/vc/*.rds",        emit: vc
    path "*/plot/*.png",      emit: plot

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def extend_by = task.ext.extend_by ?: ''
    def chimera_to_target = task.ext.chimera_to_target ?: ''
    def top_n = task.ext.top_n ?: ''

    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env Rscript

    library(CrispRVariants)
    library(rtracklayer)
    library(GenomicFeatures)
    library(BSgenome.Drerio.UCSC.danRer11)
    library(TxDb.Drerio.UCSC.danRer11.refGene)

    dir.create("${outdir}/plot", recursive = TRUE)
    dir.create("${outdir}/vc", recursive = TRUE)
    file.create("${outdir}/vc/test.rds")
    file.create("${outdir}/plot/test.png")


      # Step1: obtain gdl
    gd_fname <- "$guide_bed"
    gd <- rtracklayer::import(gd_fname)
    extend_by = $extend_by
    gdl <- GenomicRanges::resize(gd, width(gd) + extend_by, fix = "center")
    zero_coordinate <- gd\$thick@start - start(gd) + extend_by / 2

      # Step2: obtain reference
    reference <- system(sprintf("samtools faidx $ref %s:%s-%s", seqnames(gdl)[1], start(gdl)[1], end(gdl)[1]), intern = TRUE)
    reference <- paste(reference[2:length(reference)], collapse = '')
    reference <- Biostrings::reverseComplement(Biostrings::DNAString(reference))

      # Step3: create crispr_set
    bam_fname <- "$bam"

    tryCatch(
      {
        crispr_set <- readsToTarget(bam_fname,
                                    target = gdl,
                                    reference = reference,
                                    chimera.to.target = $chimera_to_target, # allow larger gaps
                                    target.loc = zero_coordinate)
        vc <- variantCounts(crispr_set)

          # Step4: plot variants
        top.n <- $top_n
        png(file="${outdir}/plot/${prefix}.png", res = 120, width = 2500, height = 1000)
        p <- plotVariants(crispr_set,
                          gene.text.size = 8,
                          row.ht.ratio = c(1,8),
                          col.wdth.ratio = c(4,2),
                          plotAlignments.args = list(line.weight = 0.5, ins.size = 2, legend.symbol.size = 4, top.n = top.n),
                          plotFreqHeatmap.args = list(plot.text.size = 3,
                                                      x.size = 5, top.n = top.n,
                                                      legend.text.size = 8,
                                                      legend.key.height = grid::unit(0.5, "lines")))

        saveRDS(vc, "${outdir}/vc/${prefix}.rds")
      },
      error = function(cond) {
        message("CrisprVRiants failed!")
        return(NA)
      }
    )


    """
}
