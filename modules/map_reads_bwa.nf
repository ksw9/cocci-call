process MapReads_BWA {

  // Map reads and remove duplicates in same step
  
  label 'slurm'

  publishDir "${projectDir}/results/${batch}/${sample_id}/bams", mode: "copy", pattern: "*_merged_mrkdup_${target_species}.bam"
  publishDir "${projectDir}/results/${batch}/${sample_id}/stats", mode: "copy", pattern: "*_mapping_${target_species}.log"
  publishDir "${projectDir}/results/${batch}/${sample_id}/stats", mode: "copy", pattern: "*_coverage_stats_${target_species}.txt"
  publishDir "${projectDir}/results/${batch}/${sample_id}/stats", mode: "copy", pattern: "*_marked_dup_metrics_${target_species}.txt"

  input:
  each target_species
  each path(reference_fasta)
  path(bwa_index)
  tuple val(sample_id), path(read1), path(read2), val(batch)

  output:
  tuple val(sample_id), val(batch), path("${sample_id}_merged_mrkdup_*.bam"), path("${sample_id}_mapping_*.log"), emit: bam_files
  path "${sample_id}_mapping_*.log", emit: mapping_reports
  path "${sample_id}_coverage_stats_*.txt", emit: coverage_stats
  path "${sample_id}_marked_dup_metrics_*.txt", emit: dup_metrics

  """
  # Check species assignment from Kraken2
  if [[ "${read1}" == *"immitis"* ]]
  then

    species="immitis"

  elif [[ "${read1}" == *"posadasii"* ]]
  then

    species="posadasii"

  else

    species="coccidioides"

  fi

  # Determine whether to run alignment or not
  if [[ "\${species}" == "${target_species}" ]]
  then

    run_toggle=true

  elif [[ "\${species}" == "coccidioides" ]]
  then

    run_toggle=true

  else

    run_toggle=false

  fi

  # Run alignment only if species is target_species or coccidioides (i.e. unclear)
  if \$run_toggle
  then

    # Get machine id_lane from SRA read identifier (old Illumina fastq format)
    read_name=\$(zcat ${read1} | head -n 1)
    read_name=\${read_name/@}
    read_name=\$(echo \${read_name} | cut -d' ' -f1)
    flowcell="\$(echo \${read_name} | cut -d: -f1-2)"
    barcode="\$(echo \${read_name} | cut -d: -f3)"
    lane="\$(echo \${read_name} | cut -d: -f4)"
    if [ "\$flowcell" == "\$lane" ]; then 
      ID=\${read_name}
      PU=\${read_name}
    else 
      ID="\${flowcell}'.'\${lane}"
      PU="\${flowcell}'.'\${barcode}'.'\${lane}"
    fi
  
    # Mapping with BWA
    bwa mem \
    -t \$SLURM_CPUS_ON_NODE \
    ${reference_fasta} \
    ${read1} ${read2} > temp.sam

    # Sort and convert to bam
    gatk SortSam \
    -I temp.sam \
    -O temp.bam \
    -SORT_ORDER coordinate

    # Extracting mapping stats
    sambamba flagstat -t \$SLURM_CPUS_ON_NODE temp.bam > ${sample_id}_mapping_${target_species}.log

    # Collect coverage stats with Picard
    picard CollectWgsMetrics \
    R=${reference_fasta} \
    I=temp.bam \
    O=${sample_id}_coverage_stats_${target_species}.txt

    # Add/replace read groups for post-processing with GATK
    picard AddOrReplaceReadGroups \
    INPUT=temp.bam \
    OUTPUT=temp_rg.bam \
    RGID=\${ID} \
    RGLB=${params.library} \
    RGPU=\${PU} \
    RGPL=${params.seq_platform} \
    RGSM=${sample_id}

    # Mark duplicates & produce library complexity metrics
    gatk MarkDuplicates \
    -I temp_rg.bam \
    -O ${sample_id}_merged_mrkdup_${target_species}.bam  \
    -M ${sample_id}_marked_dup_metrics_${target_species}.txt  
      
    # Mark duplicates with Spark (also sorts BAM) & produce library complexity metrics
    # Need to use Java 7 or Java 11 (https://gatk.broadinstitute.org/hc/en-us/community/posts/4417665825307-java-lang-reflect-InaccessibleObjectException-Unable-to-make-field-transient-java-lang-Object-java-util-ArrayList-elementData-accessible-module-java-base-does-not-opens-java-util-to-unnamed-module-3bf44630)
    #gatk MarkDuplicatesSpark \
    #-I temp_rg.bam \
    #-O ${sample_id}_merged_mrkdup_${target_species}.bam \
    #-M ${sample_id}_marked_dup_metrics.txt  
  
    # Base (Quality Score) Recalibration: not done because no 'known variants.'

  else

    # Creating mock outputs
    touch ${sample_id}_merged_mrkdup_mock.bam
    touch ${sample_id}_mapping_mock.log
    touch ${sample_id}_coverage_stats_mock.txt
    touch ${sample_id}_marked_dup_metrics_mock.txt

  fi
  """

}