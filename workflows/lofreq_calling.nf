
// ----------------Workflow---------------- //

include { VariantsLoFreq } from '../modules/variants_lofreq.nf'
include { AnnotateVCF } from '../modules/annotate_vcf.nf'

workflow LOFREQ {

  take:
  bam_files
  reference_fasta
  reference_fasta_index
	
  main:
  // LOFREQ VARIANT CALLER ---------------- //

  // Join channels by sample_id and batch
  reference_fasta
    .join(reference_fasta_index, by: [0,1], remainder: false)
    .join(bam_files, by: [0,1], remainder: false)
    .set{ lofreq_input }

  // Variant calling
  VariantsLoFreq(lofreq_input)

  // ANNOTATE GATK VCF -------------------- //

  // Join channels by sample_id and batch
  reference_fasta
    .join(VariantsLoFreq.out.lofreq_vcf_filt, by: [0,1], remainder: false)
    .set{ lofreq_vcf_filt }

  // Annotation
  AnnotateVCF("LoFreq", lofreq_vcf_filt)

}