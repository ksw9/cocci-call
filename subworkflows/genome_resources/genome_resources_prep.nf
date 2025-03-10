
// ----------------Workflow---------------- //

include { DownloadRefs } from '../../modules/download_refs/download_ref.nf'
include { RunRepeatMasker } from '../../modules/mask_refs/repeatmasker.nf'
include { RunNucmerMasker } from '../../modules/mask_refs/nucmer.nf'
include { MergeMasking } from '../../modules/mask_refs/merge_masking.nf'
include { FastaIndex } from '../../modules/ref_indexing/fasta_indexing.nf'
include { BwaIndex } from '../../modules/mapping/bwa_indexing.nf'
include { MakeGatkDict } from '../../modules/variant_calling/make_gatk_dict.nf'
include { SnpeffInputPrep } from '../../modules/snpeff_prep/snpeff_input_prep.nf'

workflow GENOMERESOURCES {

  take:
  assembly_identifier
  species
	
  main:
  // Channel for scripts directory
  scripts_dir = Channel.fromPath("${projectDir}/scripts")

  // DOWNLOAD REFERENCE GENOMES ----------- //

  DownloadRefs(assembly_identifier, species)

  // MASKING ------------------------------ //

  // RepeatMasker
  if (params.repeatmasker_mask) {

    // Channel for repeats database
    if ("${species}" == "immitis") {

      masking_file = Channel.fromPath("${projectDir}/masking_files/immitis_repeats.fa")

    }
    else {

      masking_file = Channel.fromPath("${projectDir}/masking_files/posadasii_repeats.fa")

    }
    
    RunRepeatMasker(DownloadRefs.out.fasta, masking_file)

  }

  // Nucmer
  if (params.nucmer_mask) {
    
    RunNucmerMasker(scripts_dir, DownloadRefs.out.fasta)

  }

  // Set masked fasta channel
  if (params.repeatmasker_mask & params.nucmer_mask) {

    MergeMasking(scripts_dir, RunRepeatMasker.out.masked_fasta, RunNucmerMasker.out.masked_fasta)
    masked_fasta = MergeMasking.out.masked_fasta

  }
  else if (params.repeatmasker_mask) {

    masked_fasta = RunRepeatMasker.out.masked_fasta

  }
  else if (params.nucmer_mask) {

    masked_fasta = RunNucmerMasker.out.masked_fasta

  }
  else {

    masked_fasta = DownloadRefs.out.fasta

  }

  // FASTA INDEXING

  // Index fasta
  FastaIndex(masked_fasta)

  // Merge channels
  masked_fasta
  .join(FastaIndex.out.masked_fasta_index, by: 0, remainder: false)
  .set{ masked_indexed_fasta }

  // BWA INDEXING ------------------------- //
  
  BwaIndex(masked_indexed_fasta)

  // GATK DICTIONARY GENERATION ----------- //

  MakeGatkDict(masked_indexed_fasta)

  // SNPEFF INPUT PREP

  // Merge channels
  masked_fasta
  .join(DownloadRefs.out.gff, by: 0, remainder: false)
  .set{ masked_indexed_fasta }

  SnpeffInputPrep(masked_fasta_and_gff)

  emit:
  snpeff_input = SnpeffInputPrep.out.snpeff_input

}
