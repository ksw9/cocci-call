
// ----------------Workflow---------------- //

include { VariantsLoFreq } from '../../modules/variant_calling/variants_lofreq.nf'
include { AnnotateVCF } from '../../modules/variant_calling/annotate_vcf.nf'

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

  Channel.fromPath("${params.resources_dir}/${params.snpeff_dir}")
    .set{ snpeff_dir }

  Channel.fromPath("${params.resources_dir}/${params.snpeff_datapath}")
    .set{ snpeff_datapath }

  // Join channels by sample_id and batch
  reference_fasta
    .join(VariantsLoFreq.out.lofreq_vcf_filtered, by: [0,1], remainder: false)
    .set{ lofreq_vcf_filtered }

  // Annotation
  AnnotateVCF("LoFreq", snpeff_dir, snpeff_datapath, lofreq_vcf_filtered)

}