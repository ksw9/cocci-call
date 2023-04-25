#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
INSERT PIPELINE DESCRIPTION
*/

// ----------------Workflow---------------- //

include { TrimFastQ } from './modules/trimgalore.nf'
include { Kraken } from './modules/kraken.nf'
include { MapReads_BWA as MapReads_Immitis } from './modules/map_reads_bwa.nf'
include { MapReads_BWA as MapReads_Posadasii } from './modules/map_reads_bwa.nf'
include { CheckMappedFiles } from './modules/check_mapping.nf'
include { GATK } from './workflows/gatk_calling.nf'
include { LOFREQ } from './workflows/lofreq_calling.nf'
include { SummarizeRun } from './modules/make_run_summary.nf'

workflow {

  // LOAD GENOME RESOURCES ---------------- //

  // Channels for immitis genome reference fasta
  immitis_reference_fasta = Channel.fromPath("${params.resources_dir}/${params.immitis_fasta_path}")
  
  // Channel for immitis genome reference fasta index
  immitis_reference_fasta_index = Channel.fromPath("${params.resources_dir}/${params.immitis_fasta_index_path}")

  // Channels for immitis BWA index
  Channel.fromPath("${params.resources_dir}/${params.immitis_bwa_index_path}/*{amb,ann,bwt,pac,sa}")
    .collect()
    .set{immitis_bwa_index}

  // Channel for immitis GATK dictionary (absolute path from params won't do since it has to be present in the dir where GATK is launched)
  immitis_gatk_dictionary = Channel.fromPath("${params.resources_dir}/${params.immitis_gatk_dictionary_path}")

  // Channels for posadasii genome reference fasta
  posadasii_reference_fasta = Channel.fromPath("${params.resources_dir}/${params.posadasii_fasta_path}")

  // Channel for posadasii genome reference fasta index
  posadasii_reference_fasta_index = Channel.fromPath("${params.resources_dir}/${params.posadasii_fasta_index_path}")

  // Channels for posadasii BWA index
  Channel.fromPath("${params.resources_dir}/${params.posadasii_bwa_index_path}/*{amb,ann,bwt,pac,sa}")
    .collect()
    .set{posadasii_bwa_index}

  // Channel for posadasii GATK dictionary (absolute path from params won't do since it has to be present in the dir where GATK is launched)
  posadasii_gatk_dictionary = Channel.fromPath("${params.resources_dir}/${params.posadasii_gatk_dictionary_path}")

  // Channel for Kraken2 database
  Channel.fromPath("${params.resources_dir}/${params.kraken_database_path}/*{kmer_distrib,k2d,txt,map}")
    .collect()
    .set{kraken_database}

  // CREATING RAW-READS CHANNEL ----------- //

  Channel
    .fromPath("${params.resources_dir}/${params.reads_list}")
    .splitCsv(header: true, sep: '\t')
    .map{row -> tuple(row.sample, row.fastq_1, row.fastq_2, row.batch)}
    .set{raw_reads}

  // TRIMGALORE --------------------------- //

  TrimFastQ(raw_reads)

  // KRAKEN ------------------------------- //

  // Running Kraken2
  Kraken(kraken_database, TrimFastQ.out.trimmed_fastq_files)

  // MAPPING READS TO C. IMMITIS ---------- //

  // Mapping and removing duplicates
  MapReads_Immitis("immitis", immitis_reference_fasta, immitis_bwa_index, Kraken.out.kraken_filtered_files)

  // MAPPING READS TO C. IMMITIS ---------- //

  // Mapping and removing duplicates
  MapReads_Posadasii("posadasii", posadasii_reference_fasta, posadasii_bwa_index, Kraken.out.kraken_filtered_files)

  // ASSIGN SPECIES (IF NOT DONE BY K2) --- //

  // Merge MapReads_Immitis and MapReads_Posadasii bam_files channels
  MapReads_Immitis.out.bam_files
    .join(MapReads_Posadasii.out.bam_files, by: [0,1], remainder: false)
    .set{ immitis_and_posadasii_bam_files }

  // Compare mappings and assign species if not done already by Kraken2
  CheckMappedFiles(immitis_reference_fasta, immitis_reference_fasta_index, immitis_gatk_dictionary, posadasii_reference_fasta, posadasii_reference_fasta_index, posadasii_gatk_dictionary, immitis_and_posadasii_bam_files)

  // VARIANT CALLING ---------------------- //

  // GATK variant calling, consensus fasta generation, and cvs file annotation
  GATK(CheckMappedFiles.out.bam_files, CheckMappedFiles.out.reference_fasta, CheckMappedFiles.out.reference_fasta_index, CheckMappedFiles.out.gatk_dictionary)
  
  // Running LoFreq variant calling and cvs file annotation, if desired

  if (params.run_lofreq == true) {

    LOFREQ(CheckMappedFiles.out.bam_files, CheckMappedFiles.out.reference_fasta, CheckMappedFiles.out.reference_fasta_index)

  }

  // MAKING SUMMARY REPORT ---------------- //

  // Creating channel for reads_list file (needed to parse trimming_reports)
  Channel
    .fromPath("${params.resources_dir}/${params.reads_list}")
    .set{ reads_list_file }

  SummarizeRun(reads_list_file, TrimFastQ.out.trimming_reports.flatten().collect(), Kraken.out.kraken_reports.collect(), mapping_reports.collect(), coverage_stats.collect(), dup_metrics.collect())

}