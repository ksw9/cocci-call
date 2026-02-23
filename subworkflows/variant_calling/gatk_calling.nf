
// ----------------Workflow---------------- //

include { VariantsGATK } from '../../modules/variant_calling/variants_gatk.nf'
include { IndexVCF as IndexRawVCF } from '../../modules/variant_calling/index_vcf.nf'
include { IndexVCF as IndexFilteredVCF } from '../../modules/variant_calling/index_vcf.nf'
include { FilterVCF } from '../../modules/variant_calling/filter_vcf.nf'
include { ConvertVCF } from '../../modules/variant_calling/vcf2fasta.nf'
include { AnnotateVCF } from '../../modules/variant_calling/annotate_vcf.nf'

workflow GATK {

  take:
  scripts_dir
  reference_fasta
  reference_fasta_index
  gatk_dictionary
  snpeff_dir
  snpeff_datapath
  bam_files
	
  main:
  // GATK VARIANT CALLER ------------------ //

  // Join channels by sample_id and batch
  reference_fasta
    .join(reference_fasta_index, by: [0,1], remainder: false)
    .join(gatk_dictionary, by: [0,1], remainder: false)
    .join(bam_files, by: [0,1], remainder: false)
    .set{ gatk_input }

  // Variant calling
  VariantsGATK(gatk_input)

  // Index vcf
  IndexRawVCF(VariantsGATK.out.gatk_vcf_unfiltered)

  // Filtering
  FilterVCF(scripts_dir, VariantsGATK.out.gatk_vcf_unfiltered)

  // Index vcf
  IndexFilteredVCF(FilterVCF.out.filtered_vcf)

  // CONVERTING VCF TO FASTA -------------- //

  // Define input
  if (params.use_filtered_vcf_for_fasta) {

    reference_fasta
      .join(FilterVCF.out.filtered_vcf, by: [0,1], remainder: false)
      .join(IndexFilteredVCF.out.vcf_index, by: [0,1], remainder: false)
      .set{ vcf_to_fasta_input }

  }
  else {

    reference_fasta
      .join(VariantsGATK.out.gatk_vcf_unfiltered, by: [0,1], remainder: false)
      .join(IndexRawVCF.out.vcf_index, by: [0,1], remainder: false)
      .set{ vcf_to_fasta_input }

  }

  ConvertVCF("gatk", vcf_to_fasta_input)

  // ANNOTATE GATK VCF -------------------- //

  // Join reference fasta and filtered vcf channels by sample_id and batch
  reference_fasta
    .join(FilterVCF.out.filtered_vcf, by: [0,1], remainder: false)
    .join(IndexFilteredVCF.out.vcf_index, by: [0,1], remainder: false)
    .set{ snpeff_input }

  // Annotate
  AnnotateVCF("gatk", snpeff_dir, snpeff_datapath, snpeff_input)

}
