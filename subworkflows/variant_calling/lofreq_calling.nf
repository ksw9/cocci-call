
// ----------------Workflow---------------- //

include { VariantsLoFreq } from '../../modules/variant_calling/variants_lofreq.nf'
include { IndexVCF as IndexFilteredVCF } from '../../modules/variant_calling/index_vcf.nf'
include { AnnotateVCF } from '../../modules/variant_calling/annotate_vcf.nf'

workflow LOFREQ {

  take:
  reference_fasta
  reference_fasta_index
  snpeff_dir
  snpeff_datapath
  bam_files
	
  main:
  // LOFREQ VARIANT CALLER ---------------- //

  // Join channels by sample_id and batch
  reference_fasta
    .join(reference_fasta_index, by: [0,1], remainder: false)
    .join(bam_files, by: [0,1], remainder: false)
    .set{ lofreq_input }

  // Variant calling
  VariantsLoFreq(lofreq_input)

  // Index vcf
  IndexFilteredVCF(VariantsLoFreq.out.lofreq_vcf_filtered)

  // ANNOTATE GATK VCF -------------------- //

  // Join channels by sample_id and batch
  reference_fasta
    .join(VariantsLoFreq.out.lofreq_vcf_filtered, by: [0,1], remainder: false)
    .join(IndexFilteredVCF.out.vcf_index, by: [0,1], remainder: false)
    .set{ snpeff_input }

  // Annotation
  AnnotateVCF("LoFreq", snpeff_dir, snpeff_datapath, snpeff_input)

}