/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowSikipipe.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.fasta ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CAT_FASTQ                   } from '../modules/nf-core/modules/cat/fastq/main'
include { PREPROCESS                  } from '../modules/local/preprocess'
include { CAT_STAT as CAT_STAT_PREPROCESS } from '../modules/local/cat_stat'
include { CUTADAPT                    } from '../modules/nf-core/modules/cutadapt/main'
include { UMI_EXTRACT                 } from '../modules/local/umi_extract'
include { FASTQC                      } from '../modules/nf-core/modules/fastqc/main'
include { UMI_PATTERN                 } from '../modules/local/umi_pattern'
include { READ_LENGTH_DIST            } from '../modules/local/read_length_dist'
include { READ_LENGTH_DIST_WITHIN_UMI_GROUP } from '../modules/local/read_length_dist_within_umi_group'
include { UMI_CORRECT                 } from '../modules/local/umi_correct'
include { CAT_STAT as CAT_STAT_UMI_CORRECT } from '../modules/local/cat_stat'
include { CAT_STAT as CAT_STAT_PRECISE_INSERT_UMI_1 } from '../modules/local/cat_stat'
include { CAT_STAT as CAT_STAT_PRECISE_INSERT_UMI_5 } from '../modules/local/cat_stat'

include { PREP_REF                    } from '../modules/local/prep_ref'
include { PREP_REF as PREP_REF_2      } from '../modules/local/prep_ref'
include { BWA_INDEX                   } from '../modules/nf-core/modules/bwa/index/main'
include { BWA_MEM as BWA_UMI_1        } from '../modules/nf-core/modules/bwa/mem/main'
include { BWA_MEM as BWA_UMI_5        } from '../modules/nf-core/modules/bwa/mem/main'
include { BLAT as BLAT_UMI_1          } from '../modules/local/blat'
include { BLAT as BLAT_UMI_5          } from '../modules/local/blat'
include { PRECISE_INSERT as PRECISE_INSERT_UMI_1 } from '../modules/local/precise_insert'
include { PRECISE_INSERT as PRECISE_INSERT_UMI_5 } from '../modules/local/precise_insert'
include { CRISPRVARIANTS as CRISPRVARIANTS_UMI_1 } from '../modules/local/crisprvariants'
include { CRISPRVARIANTS as CRISPRVARIANTS_UMI_5 } from '../modules/local/crisprvariants'
include { CRISPRVARIANTS as CRISPRVARIANTS_BWA_UMI_1 } from '../modules/local/crisprvariants'
include { CRISPRVARIANTS as CRISPRVARIANTS_BWA_UMI_1_INSERT } from '../modules/local/crisprvariants'
include { CRISPRVARIANTS as CRISPRVARIANTS_BWA_UMI_1_NON_INSERT } from '../modules/local/crisprvariants'
include { CRISPRVARIANTS as CRISPRVARIANTS_BWA_UMI_5 } from '../modules/local/crisprvariants'
include { CRISPRVARIANTS as CRISPRVARIANTS_BWA_UMI_5_INSERT } from '../modules/local/crisprvariants'
include { CRISPRVARIANTS as CRISPRVARIANTS_BWA_UMI_5_NON_INSERT } from '../modules/local/crisprvariants'

include { INSERT_STAT as INSERT_STAT_BLAT_1      } from '../modules/local/insert_stat'
include { INSERT_STAT as INSERT_STAT_BLAT_5      } from '../modules/local/insert_stat'
include { INSERT_STAT as INSERT_STAT_BWA_1       } from '../modules/local/insert_stat'
include { INSERT_STAT as INSERT_STAT_BWA_5       } from '../modules/local/insert_stat'

// include { MAP_LOCUS                   } from '../modules/local/map_locus'
// include { CLASSIFY_INDEL              } from '../modules/local/classify_indel'
// include { CLASSIFY_READTHROUGH        } from '../modules/local/classify_readthrough'
// include { BBMERGE                     } from '../modules/local/bbmerge'
// include { FASTQC_SINGLE               } from '../modules/local/fastqc_single'

// include { REPEAT_DIST_WITHIN_UMI_GROUP} from '../modules/local/repeat_dist_within_umi_group'
include { UMI_GROUP_STAT              } from '../modules/local/umi_group_stat'

include { MULTIQC                     } from '../modules/nf-core/modules/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow SIKIPIPE {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    // INPUT_CHECK.out.reads.view()

    // MODULE: Cat Fastq: to merge different lanes if multiple
    CAT_FASTQ (
      INPUT_CHECK.out.reads
      )
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions)

    // MODULE: preprocess: to rc and get rid of abnormal reads
    PREPROCESS (
      INPUT_CHECK.out.reads,
      "01a_preprocessed"
      )
    ch_versions = ch_versions.mix(PREPROCESS.out.versions)

    CAT_STAT_PREPROCESS (
      PREPROCESS.out.stat.collect(),
      "01a_preprocessed/stat",
      "sample_name,normal_read_count,abnormal_read_count" // header to be added
      )
    ch_versions = ch_versions.mix(CAT_STAT_PREPROCESS.out.versions)

    // MODULE: cutadapt: trim off 5p universal primer
    CUTADAPT (
      PREPROCESS.out.reads_normal,
      "01b_cutadapt_5p"
      )
    ch_versions = ch_versions.mix(CUTADAPT.out.versions)

    // MODULE: FastQC
    FASTQC (
      CUTADAPT.out.reads,
      "01c_fastqc"
      )
    ch_versions = ch_versions.mix(FASTQC.out.versions)

    // MODULE: UMI extract: extact UMI and append to read name
    UMI_EXTRACT (
      CUTADAPT.out.reads,
      "02a_umi_extract"
      )
    ch_versions = ch_versions.mix(UMI_EXTRACT.out.versions)

    // MODULE: UMI pattern:
    UMI_PATTERN (
      UMI_EXTRACT.out.reads,
      "02b_umi_pattern"
      )
    ch_versions = ch_versions.mix(UMI_PATTERN.out.versions)

    // MODULE: read length distribution: for each read, figure out the read length
    READ_LENGTH_DIST (
      UMI_EXTRACT.out.reads,
      "02c_read_length_dist"
      )
    ch_versions = ch_versions.mix(READ_LENGTH_DIST.out.versions)

    // MODULE: read length distribution within umi group
    READ_LENGTH_DIST_WITHIN_UMI_GROUP (
      READ_LENGTH_DIST.out.length,
      "02d_read_length_dist_within_umi_group"
      )
    ch_versions = ch_versions.mix(READ_LENGTH_DIST_WITHIN_UMI_GROUP.out.versions)

    // MODULE: UMI group stat: UMI read_count mean mode: 5c
    UMI_GROUP_STAT (
      READ_LENGTH_DIST.out.length,
      READ_LENGTH_DIST.out.reads,
      "02e_umi_group_stat"
      )
    ch_versions = ch_versions.mix(UMI_GROUP_STAT.out.versions)

    // MODULE: umi correct with mode
    UMI_CORRECT (
      UMI_GROUP_STAT.out.stat,
      UMI_GROUP_STAT.out.reads,
      "02f_umi_correct"
      )
    ch_versions = ch_versions.mix(UMI_CORRECT.out.versions)

    CAT_STAT_UMI_CORRECT (
      UMI_CORRECT.out.stat.collect(),
      "02f_umi_correct/stat",
      "sample_name,cutoff_1,cutoff_5,cutoff_30" // header to be added
      )
    ch_versions = ch_versions.mix(CAT_STAT_PREPROCESS.out.versions)

    // MODULE: rc reference sequence
    PREP_REF (
      params.ref,
      params.ref_need_rc,
      "03a_prep_ref"
    )
    ch_versions = ch_versions.mix(PREP_REF.out.versions)

    PREP_REF_2 (
      params.ref2,
      params.ref_need_rc,
      "03a_prep_ref_2"
    )
    ch_versions = ch_versions.mix(PREP_REF.out.versions)


    // MODULE: bwa index
    BWA_INDEX (
      PREP_REF.out.ref,
      "03b_bwa_index"
    )
    ch_versions = ch_versions.mix(BWA_INDEX.out.versions)

    // BWA_INDEX.out.index.view()
    // UMI_CORRECT.out.reads_cutoff1.view()

    // MODULE: bwa_mem
    BWA_UMI_1 (
      UMI_CORRECT.out.reads_cutoff1,
      BWA_INDEX.out.index,
      "sort",
      params.insert_fasta,
      params.insert_frac_size,
      "03c_bwa_cutoff1"
    )
    ch_versions = ch_versions.mix(BWA_UMI_1.out.versions)

    // MODULE: bwa_mem
    BWA_UMI_5 (
      UMI_CORRECT.out.reads_cutoff5,
      BWA_INDEX.out.index,
      "sort",
      params.insert_fasta,
      params.insert_frac_size,
      "03c_bwa_cutoff5"
    )
    ch_versions = ch_versions.mix(BWA_UMI_5.out.versions)

    // MODULE: BLAT
    BLAT_UMI_1 (
      UMI_CORRECT.out.reads_cutoff1,
      PREP_REF.out.ref,
      "04a_blat_cutoff1"
    )
    ch_versions = ch_versions.mix(BLAT_UMI_1.out.versions)

    // MODULE: BLAT
    BLAT_UMI_5 (
      UMI_CORRECT.out.reads_cutoff1,
      PREP_REF.out.ref,
      "04a_blat_cutoff5"
    )
    ch_versions = ch_versions.mix(BLAT_UMI_5.out.versions)

    // MODULE: how many precise insertions
    PRECISE_INSERT_UMI_1 (
      BLAT_UMI_1.out.reads,
      BLAT_UMI_1.out.blat_bam_tuned,
      "04b_precise_insert_blat_cutoff1"
    )

    CAT_STAT_PRECISE_INSERT_UMI_1 (
      PRECISE_INSERT_UMI_1.out.stat.collect(),
      "04b_precise_insert_blat_cutoff1/stat",
      "sample_name,grep,cigar" // header to be added
      )
    ch_versions = ch_versions.mix(CAT_STAT_PRECISE_INSERT_UMI_1.out.versions)

    // MODULE: how many preciesinserction
    PRECISE_INSERT_UMI_5 (
      BLAT_UMI_5.out.reads,
      BLAT_UMI_5.out.blat_bam_tuned,
      "04b_precise_insert_blat_cutoff5"
    )

    CAT_STAT_PRECISE_INSERT_UMI_5 (
      PRECISE_INSERT_UMI_5.out.stat.collect(),
      "04b_precise_insert_blat_cutoff5/stat",
      "sample_name,grep,cigar" // header to be added
      )
    ch_versions = ch_versions.mix(CAT_STAT_PRECISE_INSERT_UMI_5.out.versions)

    // MODULE: crisprvariants
    CRISPRVARIANTS_UMI_1 (
      BLAT_UMI_1.out.blat_bam_tuned,
      PREP_REF.out.ref,
      params.guide_bed,
      "05a_crisprvariants_blat_umi_1"
    )

    CRISPRVARIANTS_UMI_5 (
      BLAT_UMI_5.out.blat_bam_tuned,
      PREP_REF.out.ref,
      params.guide_bed,
      "05a_crisprvariants_blat_umi_5"
    )

    CRISPRVARIANTS_BWA_UMI_1 (
      BWA_UMI_1.out.bam,
      PREP_REF.out.ref,
      params.guide_bed,
      "05b_crisprvariants_bwa_umi_1/all_reads"
    )

    // Mappable reads and with insert
    CRISPRVARIANTS_BWA_UMI_1_INSERT (
      BWA_UMI_1.out.bam_insert,
      PREP_REF.out.ref,
      params.guide_bed,
      "05b_crisprvariants_bwa_umi_1/insert_reads"
    )

    // Mappable reads and without insert
    CRISPRVARIANTS_BWA_UMI_1_NON_INSERT (
      BWA_UMI_1.out.bam_non_insert,
      PREP_REF_2.out.ref,
      params.guide_bed,
      "05b_crisprvariants_bwa_umi_1/non_insert_reads"
    )



    CRISPRVARIANTS_BWA_UMI_5 (
      BWA_UMI_5.out.bam,
      PREP_REF.out.ref,
      params.guide_bed,
      "05b_crisprvariants_bwa_umi_5/all_reads"
    )

    // Mappable reads and with insert
    CRISPRVARIANTS_BWA_UMI_5_INSERT (
      BWA_UMI_5.out.bam_insert,
      PREP_REF.out.ref,
      params.guide_bed,
      "05b_crisprvariants_bwa_umi_5/insert_reads"
    )

    // Mappable reads and without insert
    CRISPRVARIANTS_BWA_UMI_5_NON_INSERT (
      BWA_UMI_5.out.bam_non_insert,
      PREP_REF_2.out.ref,
      params.guide_bed,
      "05b_crisprvariants_bwa_umi_5/non_insert_reads"
    )


    // MODULE: large insert stat
    INSERT_STAT_BLAT_1 (
      BLAT_UMI_1.out.blat_bam_tuned,
      "06a_large_insert_stat_blat_umi_1"
    )

    INSERT_STAT_BLAT_5 (
      BLAT_UMI_5.out.blat_bam_tuned,
      "06a_large_insert_stat_blat_umi_5"
    )

    INSERT_STAT_BWA_1 (
      BWA_UMI_1.out.bam,
      "06a_large_insert_stat_bwa_umi_1"
    )

    INSERT_STAT_BWA_5 (
      BWA_UMI_5.out.bam,
      "06a_large_insert_stat_bwa_umi_5"
    )

    // MODULE: read length distribution within umi group
    // READ_LENGTH_DIST_WITHIN_UMI_GROUP (
    //   REPEAT_DIST_DISTANCE.out.count_r1,
    //   "5b_r1_repeat_dist_within_umi_group"
    //   )
    // ch_versions = ch_versions.mix(REPEAT_DIST_WITHIN_UMI_GROUP.out.versions)
    //

    // //
    // // MODULE: map locus
    // //
    // MAP_LOCUS (
    //   CUTADAPT.out.reads
    //   )
    // ch_versions = ch_versions.mix(MAP_LOCUS.out.versions)
    //
    // //
    // // MODULE: combine MAP_LOCUS.out.stat into one file
    // //
    // CAT_STAT (
    //   MAP_LOCUS.out.stat.collect(),
    //   "1a_map_locus/stat",
    //   "sample_name,locus,percent,misprimed,percent,problem,percent" // header to be added
    //   )
    // ch_versions = ch_versions.mix(CAT_STAT.out.versions)
    //
    // //
    // // MODULE: classify INDEL
    // //
    // CLASSIFY_INDEL (
    //   MAP_LOCUS.out.reads_locus
    //   )
    // ch_versions = ch_versions.mix(CLASSIFY_INDEL.out.versions)
    //
    // //
    // // MODULE: combine CLASSIFY_INDEL.out.stat into one file
    // //
    // CAT_STAT2 (
    //   CLASSIFY_INDEL.out.stat.collect(),
    //   "3a_classify_indel/stat",
    //   "sample_name,no_indel,no_indel_percent,indel_5p,indel_5p_percent,indel_3p,indel_3p_percent,indel_5p_3p,indel_5p_3p_percent" // header to be added
    //   )
    // ch_versions = ch_versions.mix(CAT_STAT2.out.versions)
    //
    // //
    // // MODULE: classify_readthrough
    // //
    // CLASSIFY_READTHROUGH (
    //   CLASSIFY_INDEL.out.reads_no_indel
    //   )
    // ch_versions = ch_versions.mix(CLASSIFY_READTHROUGH.out.versions)
    //
    // //
    // // MODULE: combine CLASSIFY_READTHROUGH.out.stat into one file
    // //
    // CAT_STAT3 (
    //   CLASSIFY_READTHROUGH.out.stat.collect(),
    //   "4a_classify_readthrough/stat",
    //   "sample_name,count_readthrough,count_readthrough_percent,count_non_readthrough,p_count_non_readthrough_percent" // header to be added
    //   )
    // ch_versions = ch_versions.mix(CAT_STAT3.out.versions)
    //
    // //
    // // MODULE: BBmerge
    // //
    // BBMERGE (
    //   CLASSIFY_READTHROUGH.out.reads_through
    //   )
    // ch_versions = ch_versions.mix(BBMERGE.out.versions)
    //
    // //
    // // MODULE: combine CLASSIFY_READTHROUGH.out.stat into one file
    // //
    // CAT_STAT4 (
    //   BBMERGE.out.stat.collect(),
    //   "4b_bbmerge/stat",
    //   "sample_name,count_merge,count_non_merge" // header to be added
    //   )
    // ch_versions = ch_versions.mix(CAT_STAT4.out.versions)
    //
    // //
    // // MODULE: FastQC
    // //
    // FASTQC1 (
    //   CLASSIFY_READTHROUGH.out.reads_through,
    //   "4c_fastqc_r1_r2"
    //   )
    // ch_versions = ch_versions.mix(FASTQC1.out.versions)
    //
    // //
    // // MODULE: FastQC
    // //
    // FASTQC_SINGLE (
    //   BBMERGE.out.reads_merged,
    //   "4c_fastqc_merged"
    //   )
    // ch_versions = ch_versions.mix(FASTQC_SINGLE.out.versions)
    //
    // //
    // // MODULE: repeat distribution distance with R1/R2 reads
    // //
    // REPEAT_DIST_DISTANCE (
    //   CLASSIFY_READTHROUGH.out.reads_through,
    //   "4d_repeat_distribution_distance"
    //   )
    // ch_versions = ch_versions.mix(REPEAT_DIST_DISTANCE.out.versions)
    //
    // //
    // // MODULE: repeat distribution distance with merged reads
    // //
    // REPEAT_DIST_DISTANCE_MERGED (
    //   BBMERGE.out.reads_merged,
    //   "4d_repeat_distribution_distance_merged"
    //   )
    // ch_versions = ch_versions.mix(REPEAT_DIST_DISTANCE_MERGED.out.versions)
    //
    // //
    // // MODULE: UMI pattern: 5a
    // //
    // UMI_PATTERN2 (
    //   CLASSIFY_READTHROUGH.out.reads_through,
    //   "5a_umi_pattern"
    //   )
    // ch_versions = ch_versions.mix(UMI_PATTERN2.out.versions)
    //
    //

    // //
    // // MODULE: repeat distribution within umi group
    // //
    // REPEAT_DIST_WITHIN_UMI_GROUP2 (
    //   REPEAT_DIST_DISTANCE.out.count_r2,
    //   "5b_r2_repeat_dist_within_umi_group"
    //   )
    // ch_versions = ch_versions.mix(REPEAT_DIST_WITHIN_UMI_GROUP2.out.versions)
    //
    //
    // //
    // // MODULE: UMI group stat: UMI read_count mean mode: 5c
    // //
    // UMI_GROUP_STAT (
    //   REPEAT_DIST_DISTANCE.out.count_r1,
    //   "5c_r1_umi_group_stat"
    //   )
    // ch_versions = ch_versions.mix(UMI_GROUP_STAT.out.versions)
    //
    // //
    // // MODULE: repeat dist UMI corrected: 5d
    // //
    // REPEAT_DIST_UMI_CORRECT (
    //   UMI_GROUP_STAT.out.stat,
    //   "5d_r1_repeat_dist_umi_correct"
    //   )
    // ch_versions = ch_versions.mix(REPEAT_DIST_UMI_CORRECT.out.versions)
    //
    //
    // //
    // // MODULE: combine CLASSIFY_READTHROUGH.out.stat into one file
    // //
    // // REPEAT_DIST_DISTANCE.out.frac_r1.collect().view()
    // CAT_STAT5 (
    //   REPEAT_DIST_DISTANCE.out.frac_r1.collect(),
    //   "4d_repeat_distribution_distance/frac_r1",
    //   "blow,above" // header to be added
    //   )
    // ch_versions = ch_versions.mix(CAT_STAT5.out.versions)
    //
    // //
    // // MODULE: combine CLASSIFY_READTHROUGH.out.stat into one file
    // //
    // CAT_STAT6 (
    //   REPEAT_DIST_DISTANCE.out.frac_r2.collect(),
    //   "4d_repeat_distribution_distance/frac_r2",
    //   "blow,above" // header to be added
    //   )
    // ch_versions = ch_versions.mix(CAT_STAT6.out.versions)
    //
    // //
    // // MODULE: combine REPEAT_DIST_UMI_CORRECT.out.frac_x into one file
    // //
    // // REPEAT_DIST_UMI_CORRECT.out.frac_1.collect().view()
    // CAT_STAT7 ( REPEAT_DIST_UMI_CORRECT.out.frac_1.collect(), "5d_r1_repeat_dist_umi_correct/frac_1", "blow,above" )
    // ch_versions = ch_versions.mix(CAT_STAT7.out.versions)
    // CAT_STAT8 ( REPEAT_DIST_UMI_CORRECT.out.frac_3.collect(), "5d_r1_repeat_dist_umi_correct/frac_3", "blow,above" )
    // ch_versions = ch_versions.mix(CAT_STAT8.out.versions)
    // CAT_STAT9 ( REPEAT_DIST_UMI_CORRECT.out.frac_10.collect(), "5d_r1_repeat_dist_umi_correct/frac_10", "blow,above" )
    // ch_versions = ch_versions.mix(CAT_STAT9.out.versions)
    // CAT_STAT10 ( REPEAT_DIST_UMI_CORRECT.out.frac_30.collect(), "5d_r1_repeat_dist_umi_correct/frac_30", "blow,above" )
    // ch_versions = ch_versions.mix(CAT_STAT10.out.versions)
    // CAT_STAT11 ( REPEAT_DIST_UMI_CORRECT.out.frac_100.collect(), "5d_r1_repeat_dist_umi_correct/frac_100", "blow,above" )
    // ch_versions = ch_versions.mix(CAT_STAT11.out.versions)

    //
    // MODULE: repeat distribution R1 distance
    //
    // REPEAT_DIST_DISTANCE_MERGED (
    //   BBMERGE.out.reads_merged,
    //   "4c_merge_fastqc"
    //   )
    // ch_versions = ch_versions.mix(FASTQC_SINGLE.out.versions)

    //
    // MODULE: repeat distribution R1 distance
    //
    // REPEAT_DIST_ALIGNMENT (
    //   BBMERGE.out.reads_merged,
    //   "4c_merge_fastqc"
    //   )
    // ch_versions = ch_versions.mix(FASTQC_SINGLE.out.versions)
    //
    // //
    // // MODULE: repeat distribution R1 distance
    // //
    // REPEAT_DIST_ALIGNMENT_MERGED (
    //   BBMERGE.out.reads_merged,
    //   "4c_merge_fastqc"
    //   )
    // ch_versions = ch_versions.mix(FASTQC_SINGLE.out.versions)


    // MAP_LOCUS.out.stat.collect()
    //
    // MODULE: umi distribution statistics
    //
    // UMI_STAT (
    //   CUTADAPT.out.reads
    //   )
    // ch_versions = ch_versions.mix(UMI_STAT.out.versions)


    //
    // MODULE: Run FastQC
    //
    // FASTQC (
    //     INPUT_CHECK.out.reads
    // )
    // ch_versions = ch_versions.mix(FASTQC.out.versions.first())
    //
    // CUSTOM_DUMPSOFTWAREVERSIONS (
    //     ch_versions.unique().collectFile(name: 'collated_versions.yml')
    // )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowSikipipe.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    // ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    // ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    // MULTIQC (
    //     ch_multiqc_files.collect()
    // )
    // multiqc_report = MULTIQC.out.report.toList()
    // ch_versions    = ch_versions.mix(MULTIQC.out.versions)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
