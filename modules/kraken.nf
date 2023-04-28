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
  tuple val("${sample_id}"), path("${sample_id}_*_kr_1.fq.gz"), path("${sample_id}_*_kr_2.fq.gz"), val("${batch}"), emit: kraken_filtered_files

  """
  # run kraken to taxonomically classify paired-end reads and write output file.
  kraken2 --db . --paired --gzip-compressed --threads \$SLURM_CPUS_ON_NODE --report ${sample_id}_kraken.report --report-minimizer-data --use-names ${read1} ${read2} --output ${sample_id}.out

  # Remove Illumina suffixes from read names (Kraken reads list does not include suffixes) 
  zcat ${read1} | sed 's|/1\$||' | bgzip > ${sample_id}_plain_1.fq.gz
  zcat ${read2} | sed 's|/2\$||' | bgzip > ${sample_id}_plain_2.fq.gz

  # Select for reads directly assigned to Coccidioides genus (G) (taxid 5500), reads assigned directly to Coccidioides immitis complex (taxid 5501), and reads assigned directly to Coccidioides posadasii complex (taxid 199306).
  #grep -E 'Coccidioides (taxid 5500)|Coccidioides immitis|Coccidioides posadasii' ${sample_id}.out | awk '{print \$2}' > ${sample_id}_reads.list
  grep -E 'Coccidioides (taxid 5500)|Coccidioides immitis|Coccidioides posadasii' ${sample_id}.out | awk '{print \$2" "}' > ${sample_id}_reads.list

  # Use bbmap to select reads corresponding to taxa of interest.
  #filterbyname.sh int=false in1=${sample_id}_plain_1.fq.gz in2=${sample_id}_plain_2.fq.gz out1=${sample_id}_kr_1.fq.gz out2=${sample_id}_kr_2.fq.gz names=${sample_id}_reads.list include=true overwrite=true

  # bbmap tends to glitch, so the following code is meant to replace it
  bgzip -d ${sample_id}_plain_1.fq.gz
  awk '{ if(\$0 ~ /@/) { print \$0" " } else { print \$0 } }' ${sample_id}_plain_1.fq > ${sample_id}_plain_1_mod.fq
  awk '{ print }' ${sample_id}_reads.list | grep -f - -A 3 ${sample_id}_plain_1_mod.fq | bgzip > ${sample_id}_kr_1.fq.gz
  
  bgzip -d ${sample_id}_plain_2.fq.gz
  awk '{ if(\$0 ~ /@/) { print \$0" " } else { print \$0 } }' ${sample_id}_plain_2.fq > ${sample_id}_plain_2_mod.fq
  awk '{ print }' ${sample_id}_reads.list | grep -f - -A 3 ${sample_id}_plain_2_mod.fq | bgzip > ${sample_id}_kr_2.fq.gz
  
  # Summarize Kraken statistics. 
  #${projectDir}/scripts/kraken_stats.sh ${sample_id}_kraken.report

  # Assign species based on number of distinct minimizers
  # If the log ratio of the immitis vs posadasii is above/below a certain absolute threshold, the species is assigned, otherwise it's left as "coccidioides"
  # N.B. Log is calculated using awk, which uses natural logs, so has to be transformed to base 2
  immitis_minimizers=\$(grep "Coccidioides immitis" ${sample_id}_kraken.report | head -1 | awk '{ score = log(\$5) / log(2) ; print score }')
  posadasii_minimizers=\$(grep "Coccidioides posadasii" ${sample_id}_kraken.report | head -1 | awk '{ score = log(\$5) / log(2) ; print score }')
  assignment_score=\$(echo "\$immitis_minimizers \$posadasii_minimizers" | awk '{ print \$1-\$2 }')
  species=\$(echo "\$assignment_score ${params.minimizers_log_ratio_thr}" | awk '{ if(\$1 > \$2) { print "immitis" } else if(\$1 < -\$2) { print "posadasii" } else { print "coccidioides" } }')

  # Modify kr_2.fq.gz files to incorporate species name
  mv ${sample_id}_kr_1.fq.gz ${sample_id}_\${species}_kr_1.fq.gz
  mv ${sample_id}_kr_2.fq.gz ${sample_id}_\${species}_kr_2.fq.gz
  """

}
