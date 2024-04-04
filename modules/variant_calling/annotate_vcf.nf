process AnnotateVCF {

  // Annotate VCF for variant examination
  
  label 'variantcalling'

  publishDir "${projectDir}/results/${batch}/${sample_id}/vars", mode: "copy", pattern: "${sample_id}_${variant_caller}_ann.vcf.gz"

  input:
  each variant_caller
  each path(snpeff_dir)
  each path(snpeff_datapath)
  tuple val(sample_id), val(batch), path(reference), path(vcf)

  output:
  path "${sample_id}_${variant_caller}_ann.vcf.gz"

  """
  # Index vcf
  tabix ${vcf}

  if [[ "${reference}" == "immitis.fasta" ]]
  then

    java -jar -Xmx8g ${snpeff_dir}/snpEff.jar eff ${params.snpeff_immitis_db} ${vcf} -dataDir ${snpeff_datapath} -noStats -no-downstream -no-upstream -canon | bgzip > ${sample_id}_${variant_caller}_ann.vcf.gz

  else

    java -jar -Xmx8g ${snpeff_dir}/snpEff.jar eff ${params.snpeff_posadasii_db} ${vcf} -dataDir ${snpeff_datapath} -noStats -no-downstream -no-upstream -canon | bgzip > ${sample_id}_${variant_caller}_ann.vcf.gz
  
  fi
  """

}
