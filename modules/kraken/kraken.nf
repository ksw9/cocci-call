process Kraken {

  // Filter reads taxonomically with Kraken
  
  label 'kraken2'

  //publishDir "${projectDir}/results/${batch}/${sample_id}/kraken", mode: "copy", pattern: "*_kr_{1,2}.fq.gz"
  publishDir "${projectDir}/results/${batch}/${sample_id}/kraken", mode: "copy", pattern: "*_kraken.report"

  input:
  path(kraken_db)
  tuple val(sample_id), val(batch), path(read1), path(read2)

  output:
  path "*_kraken.report", emit: kraken_reports
  tuple val(sample_id), val(batch), path("${sample_id}_*_kr_1.fq.gz"), path("${sample_id}_*_kr_2.fq.gz"), emit: kraken_filtered_files

  """
  # run kraken to taxonomically classify paired-end reads and write output file.
  kraken2 --db . --paired --gzip-compressed --threads \$SLURM_CPUS_ON_NODE --report ${sample_id}_kraken.report --report-minimizer-data --use-names ${read1} ${read2} --output ${sample_id}.out

  if [ ${params.kraken_filtering} == true ]
  then

    # Get list of reads that belong to Coccidioides
    grep -E 'Coccidioides (taxid 5500)|Coccidioides immitis|Coccidioides posadasii' ${sample_id}.out | awk '{print \$2}' > ${sample_id}_reads.list
    
    # Use seqtk to select reads corresponding to the Coccidioides genus
    seqtk subseq ${read1} ${sample_id}_reads.list | bgzip > ${sample_id}_kr_1.fq.gz
    seqtk subseq ${read2} ${sample_id}_reads.list | bgzip > ${sample_id}_kr_2.fq.gz

  else

    cp ${read1} ${sample_id}_kr_1.fq.gz
    cp ${read2} ${sample_id}_kr_2.fq.gz

  fi
    
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
