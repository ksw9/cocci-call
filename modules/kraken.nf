process Kraken {

  // Filter reads taxonomically with Kraken
  
  label 'slurm'

  //publishDir "${projectDir}/results/${batch}/${sample_id}/kraken", mode: "copy", pattern: "*_kr_{1,2}.fq.gz"
  publishDir "${projectDir}/results/${batch}/${sample_id}/kraken", mode: "copy", pattern: "*_kraken.report"

  input:
  path(kraken_db)
  tuple val(sample_id), path(read1), path(read2), val(batch)

  output:
  path "*_kraken.report", emit: kraken_reports
  tuple val("${sample_id}"), path("${sample_id}_kr_1.fq.gz"), path("${sample_id}_kr_2.fq.gz"), val("${batch}"), emit: kraken_filtered_files

  """
  # run kraken to taxonomically classify paired-end reads and write output file.
  kraken2 --db . --paired --gzip-compressed --threads \$SLURM_CPUS_ON_NODE --report ${sample_id}_kraken.report --report-minimizer-data --use-names ${read1} ${read2} --output ${sample_id}.out

  # Remove Illumina suffixes from read names (Kraken reads list does not include suffixes) 
  zcat ${read1} | sed 's|/1\$||' | bgzip > ${sample_id}_plain_1.fq.gz
  zcat ${read2} | sed 's|/2\$||' | bgzip > ${sample_id}_plain_2.fq.gz

  # Select for reads directly assigned to Coccidioides genus (G) (taxid 5500), reads assigned directly to Coccidioides immitis complex (taxid 5501), and reads assigned directly to Coccidioides posadasii complex (taxid 199306).
  grep -E 'Coccidioides (taxid 5500)|Coccidioides immitis|Coccidioides posadasii' ${sample_id}.out | awk '{print \$2}' > ${sample_id}_reads.list

  # Use bbmap to select reads corresponding to taxa of interest.
  filterbyname.sh int=false in1=${sample_id}_plain_1.fq.gz in2=${sample_id}_plain_2.fq.gz out1=${sample_id}_kr_1.fq.gz out2=${sample_id}_kr_2.fq.gz names=${sample_id}_reads.list include=true overwrite=true
  
  # Summarize Kraken statistics. 
  #${projectDir}/scripts/kraken_stats.sh ${sample_id}_kraken.report
  """

}
