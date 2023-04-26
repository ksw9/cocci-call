
// ----------------Workflow---------------- //

include { VariantsGATK } from '../modules/variants_gatk.nf'
include { ConvertVCF } from '../modules/vcf2fasta.nf'
include { AnnotateVCF } from '../modules/annotate_vcf.nf'

workflow GATK {

  take:
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

  // CONVERTING VCF TO FASTA -------------- //

  // Join channels by sample_id and batch
  reference_fasta
    .join(VariantsGATK.out.gatk_vcf_filt, by: [0,1], remainder: false)
    .set{ gatk_vcf_filt }

  ConvertVCF("gatk", gatk_vcf_filt)

  // ANNOTATE GATK VCF -------------------- //

  // Channels for snpEff resources
  Channel.fromPath("${params.resources_dir}/${params.snpeff_dir}")
    .set{ snpeff_dir }

  Channel.fromPath("${params.resources_dir}/${params.snpeff_datapath}")
    .set{ snpeff_datapath }

  // Join channels by sample_id and batch
  reference_fasta
    .join(VariantsGATK.out.gatk_vcf_filt, by: [0,1], remainder: false)
    .set{ gatk_vcf_filt }

  // Annotation
  AnnotateVCF("gatk", snpeff_dir, snpeff_datapath, gatk_vcf_filt)

}
