process AnnotateVCF {

  // Annotate VCF for variant examination
  
  label 'slurm'

  publishDir "${projectDir}/results/${batch}/${sample_id}/vars", mode: "copy", pattern: "*_${variant_caller}_ann.vcf.gz"

  input:
  each variant_caller
  tuple val(sample_id), val(batch), path(reference), path(vcf)

  output:
  path "${sample_id}_${variant_caller}_ann.vcf.gz"

  """
  # Index vcf
  tabix ${vcf}

  if [[ "${reference}" == "immitis.fasta" ]]
  then

    java -jar -Xmx8g ${params.resources_dir}/${params.snpeff} eff ${params.snpeff_immitis_db} ${sample_id}_renamed.vcf.gz -c ${params.resources_dir}/${params.snpeff_config} -noStats -no-downstream -no-upstream -canon | bgzip > ${sample_id}_${variant_caller}_ann.vcf.gz

  else

    # Run snpEff
    java -jar -Xmx8g ${params.resources_dir}/${params.snpeff} eff ${params.snpeff_posadasii_db} ${sample_id}_renamed.vcf.gz -c ${params.resources_dir}/${params.snpeff_config} -noStats -no-downstream -no-upstream -canon | bgzip > ${sample_id}_${variant_caller}_ann.vcf.gz
  
  fi
  """

}
