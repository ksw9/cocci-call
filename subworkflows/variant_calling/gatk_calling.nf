
// ----------------Workflow---------------- //

include { VariantsGATK } from '../../modules/variant_calling/variants_gatk.nf'
include { FilterVCF } from '../../modules/variant_calling/filter_vcf.nf'
include { ConvertVCF } from '../../modules/variant_calling/vcf2fasta.nf'
include { AnnotateVCF } from '../../modules/variant_calling/annotate_vcf.nf'

workflow GATK {

  take:
  scripts_dir
  bam_files
  reference_fasta
  reference_fasta_index
  gatk_dictionary
	
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

  // Filtering
  FilterVCF(scripts_dir, VariantsGATK.out.gatk_vcf_unfiltered)

  // Join reference fasta and filtered vcf channels by sample_id and batch
  reference_fasta
    .join(FilterVCF.out.filtered_vcf, by: [0,1], remainder: false)
    .set{ gatk_filtered_vcf }

  // CONVERTING VCF TO FASTA -------------- //

  ConvertVCF("gatk", gatk_filtered_vcf)

  // ANNOTATE GATK VCF -------------------- //

  // Channels for snpEff resources
  Channel.fromPath("${params.resources_dir}/${params.snpeff_dir}")
    .set{ snpeff_dir }

  Channel.fromPath("${params.resources_dir}/${params.snpeff_datapath}")
    .set{ snpeff_datapath }

  // Annotation
  AnnotateVCF("gatk", snpeff_dir, snpeff_datapath, gatk_filtered_vcf)

}
