process ConvertVCF {

  // Convert single sample VCF to fasta

  label 'variantcalling'

  publishDir "${projectDir}/results/${batch}/${sample_id}/fasta", mode: "copy", pattern: "*.fa"

  input:
  each variant_caller
  tuple val(sample_id), val(batch), path(reference), path(vcf)

  output:
  path "${sample_id}_${variant_caller}.fa", emit: unmasked_fasta

  """
  # Index vcf
  tabix -p vcf ${vcf}

  # Output - Consensus (exclude indels)
  # N.B. The vcf files come from individual samples, so no need to specify --sample in bcftools consensus (also, LoFreq does not store sample name info in the vcf).
  if [ ${params.fasta_calls_only} == false ]
  then

    bcftools consensus --include "TYPE!='indel'" --fasta-ref ${reference} ${vcf}  > ${sample_id}_${variant_caller}.fa

  else

    bcftools consensus --include "TYPE!='indel'" --fasta-ref ${reference} --missing 'N' --absent 'N' ${vcf}  > ${sample_id}_${variant_caller}.fa

  fi
  """

}
