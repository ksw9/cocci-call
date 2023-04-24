process VariantsLoFreq {

  // Variant calling with LoFreq
  
  label 'slurm'

  publishDir "${projectDir}/results/${batch}/${sample_id}/vars", mode: "copy", pattern: "*{_lofreq_unfilt,_lofreq_filt}.vcf.gz"

  input:
  tuple val(sample_id), val(batch), path(reference), path(reference_index), path(bam)

  output:
  tuple val(sample_id), val(batch), path("${sample_id}_lofreq_unfilt.vcf.gz"), emit: lofreq_vcf_unfilt
  tuple val(sample_id), val(batch), path("${sample_id}_lofreq_filt.vcf.gz"), emit: lofreq_vcf_filt

  """
  # Indexing bam
  samtools index ${bam}
  
  # Call variants with LoFreq, no filter
  lofreq call-parallel --call-indels --pp-threads \$SLURM_CPUS_ON_NODE --no-default-filter -f ${reference} -o ${sample_id}_lofreq_unfilt.vcf ${bam}
  
  # Bgzipping
  bgzip ${sample_id}_lofreq_unfilt.vcf
  
  # Call variants with LoFreq, default filter
  lofreq call-parallel --call-indels --pp-threads \$SLURM_CPUS_ON_NODE -f ${reference} -o ${sample_id}_lofreq_filt.vcf ${bam}
  
  # Bgzipping
  bgzip ${sample_id}_lofreq_filt.vcf
  """

}
