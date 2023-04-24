process ConvertVCF {

  // Convert single sample VCF to fasta

  label 'slurm'

  publishDir "${projectDir}/results/${batch}/${sample_id}/fasta", mode: "copy", pattern: "*.fa"

  input:
  each variant_caller
  tuple val(sample_id), val(batch), path(reference), path(vcf)

  output:
  path "${sample_id}_${variant_caller}.fa", emit: unmasked_fasta

  """
  # Index vcf
  tabix -p vcf ${vcf}

  # N.B. The vcf files come from individual samples, so no need to specify --sample in bcftools consensus (also, LoFreq does not store sample name info in the vcf).

  # Output - Consensus without ppe masking, but with quality filters applied (exclude indels)
  bcftools consensus --include "(TYPE!='indel' && FORMAT/DP >= ${params.depth_threshold}) && (QUAL >= ${params.qual_threshold}|| RGQ  >= ${params.qual_threshold})" --fasta-ref ${reference} --missing 'N' --absent 'N' ${vcf}  > ${sample_id}_${variant_caller}.fa
  """

}
