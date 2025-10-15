process VariantsGATK {

  // Variant calling with GATK
  
  label 'variantcalling'

  publishDir "${projectDir}/results/${batch}/${sample_id}/vars", mode: "copy", pattern: "*_{gatk.g,gatk_unfilt,gatk_unfilt_norm}.vcf.gz"

  input:
  tuple val(sample_id), val(batch), path(reference), path(reference_index), path(dictionary), path(bam)

  output:
  //tuple val(sample_id), val(batch), path("${sample_id}_gatk.g.vcf.gz"), emit: gatk_gvcf
  tuple val(sample_id), val(batch), path("${sample_id}_gatk_unfilt_norm.vcf.gz"), emit: gatk_vcf_unfiltered

  """
  # Indexing bam
  samtools index ${bam}

  if [ ${params.variants_only} == false ]
  then 
  
    # Call variants with GATK, output GVCF
    # ERC: Reference model emitted with condensed non-variant blocks, i.e. the GVCF format
    gatk --java-options "-Xmx4g" HaplotypeCaller \
    -R ${reference} \
    -ploidy ${params.ploidy} \
    -I ${bam} \
    -ERC BP_RESOLUTION \
    --output-mode EMIT_ALL_CONFIDENT_SITES \
    -O ${sample_id}_gatk.g.vcf.gz

    # GVCF to VCF. Min base quality score is 10 by default. Including non-variant sites in order to differentiate between consensus call and no-call sites.
    gatk --java-options '-Xmx100g' GenotypeGVCFs \
    -R ${reference} \
    --variant ${sample_id}_gatk.g.vcf.gz \
    -ploidy ${params.ploidy} \
    --include-non-variant-sites true \
    --output ${sample_id}_gatk_unfilt.vcf.gz

    # vcf normalization
    bcftools norm \
    -f ${reference} \
    -m -both \
    -O z \
    -o ${sample_id}_gatk_unfilt_norm.vcf.gz \
    ${sample_id}_gatk_unfilt.vcf.gz
  
  else

    # Call variants with GATK, output GVCF
    # ERC: Reference model emitted with condensed non-variant blocks, i.e. the GVCF format
    gatk --java-options "-Xmx4g" HaplotypeCaller \
    -R ${reference} \
    -ploidy ${params.ploidy} \
    -I ${bam} \
    -ERC GVCF \
    --output-mode EMIT_VARIANTS_ONLY \
    -O ${sample_id}_gatk.g.vcf.gz

    # GVCF to VCF. Min base quality score is 10 by default. Including non-variant sites in order to differentiate between consensus call and no-call sites.
    gatk --java-options '-Xmx100g' GenotypeGVCFs \
    -R ${reference} \
    --variant ${sample_id}_gatk.g.vcf.gz \
    -ploidy ${params.ploidy} \
    --include-non-variant-sites false \
    --output ${sample_id}_gatk_unfilt.vcf.gz

    # vcf normalization
    bcftools norm \
    -f ${reference} \
    -m -both \
    -O z \
    -o ${sample_id}_gatk_unfilt_norm.vcf.gz \
    ${sample_id}_gatk_unfilt.vcf.gz
  
  fi
  """

}
