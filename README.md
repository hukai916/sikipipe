## Introduction

**sikipipe** is a Nextflow pipeline designed to analyze CRISPR/Cas variants with PacBio sequencing.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

## Pipeline summary

Currently, sikipipe has the following modules:
1. **Preprocess:**  
  a. 00a_merge_lane: merge laned if multiple lanes  
  b. 01a_preprocessed: preprocess: reverse complement if needed and filter out abnormal length  
  c. 01b_cutadap_5p: trim off 5p universal primer   
  d. 01c_fastqc: FastQC on 1c results
2. **UMI-extract:**   
  a. 02a_umi_extract: extract UMI and append to read name   
  b. 02b_umi_pattern: plot how many reads per UMI
  c. 02c_read_length_dist: read length distribution   
  d. 02d_read_length_dist_within_umi_group: read length distribution within UMI group   
  e. 02e_umi_group_stat: umi group stat: mean and mode   
  f. 02f_umi_correct: UMI-correction(collapse) with mode using different cutoffs
3. **Map with BWA:**    
  a. 03a_prep_ref: prepare reference, reverse complement if specified   
  b. 03b_bwa_index: build BWA index   
  c. 03c_bwa_cutoff1: map 02f result with bwa, umi-corrected with cutoff_1  
  d. 03c_bwa_cutoff5: map 02f result with bwa, umi-corrected with cutoff_5    
4.  **Map with BLAT:** by looking at the CIGAR, many reads only mapped a small fraction   
  a. 04a_blat_cutoff1: map 02f result with BLAT, umi-corrected with cutoff_1    
  b. 04a_blat_cutoff5: map 02f result with BLAT, umi-corrected with cutoff_5    
  c. 04b_precise_insert_blat_cutoff1: how many precise insert by using grep method    
  d. 04b_precise_insert_blat_cutoff5: how many precise insert by using grep method
5. **CrispRVariants:** need to supply reverse complemented reference if desired
  a. 05a_crisprvariants_blat_umi_1: run CrisprVRiants and plot variants: BLAT, UMI_cutoff1    
  b. 05a_crisprvariants_blat_umi_5: run CrisprVRiants and plot variants: BLAT, UMI_cutoff5    
  c: 05b_crisprvariants_bwa_umi_1: run CrisprVRiants and plot variants: BWA, UMI_cutoff1    
  d: 05b_crisprvariants_bwa_umi_5: run CrisprVRiants and plot variants: BWA, UMI_cutoff5    
6. **Large insert stat**: since CrispRVariants drops the naems of the reads, parse BAM to figure out large inserts
  a. 06a_large_insert_stat_blat_umi_1: stat/all_sample.csv    
  Header: sample_name, CIGAR, insert_size, pos_ref, insert_seq, read_seq    
  b. 06a_large_insert_stat_blat_umi_5   
  c. 06a_large_insert_stat_bwa_umi_1    
  d. 06a_large_insert_stat_bwa_umi_5    

## TODO:
1. Analyze other samples in addition to the 6 samples.
2. Debug the BLAT output to bam issue.
