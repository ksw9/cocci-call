#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
INSERT PIPELINE DESCRIPTION
*/

// ----------------Workflow---------------- //

include { TrimFastQ } from './modules/trimgalore.nf'
include { Kraken } from './modules/kraken.nf'

workflow {

  // CREATING RAW-READS CHANNEL ----------- //

  Channel
  .fromPath("${params.resources_dir}/${params.reads_list}")
  .splitCsv(header: true, sep: '\t')
  .map{row -> tuple(row.sample, row.fastq_1, row.fastq_2, row.batch)}
  .set{raw_reads}

  // TRIMGALORE --------------------------- //

  TrimFastQ(raw_reads)

  // KRAKEN ------------------------------- //

  // Channel for Kraken2 database
  Channel.fromPath("${params.resources_dir}/${params.kraken_database_path}/*{kmer_distrib,k2d,txt,map}")
  .collect()
  .set{kraken_database}

  // Running Kraken2
  Kraken(kraken_database, TrimFastQ.out.trimmed_fastq_files)

}