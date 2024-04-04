process CheckMappedFiles {

  // Check BWA mapping and assign species if not already doen so in Kraken2 step
  
  label 'kraken2'

  input:
  each path(immitis_reference_fasta)
  each path(immitis_reference_fasta_index)
  each path(immitis_gatk_dictionary)
  each path(posadasii_reference_fasta)
  each path(posadasii_reference_fasta_index)
  each path(posadasii_gatk_dictionary)
  tuple val(sample_id), val(batch), path(immitis_bam), path(immitis_mapping_report), path(posadasii_bam), path(posadasii_mapping_report)

  output:
  tuple val(sample_id), val(batch), path("${sample_id}_*.bam"), emit: bam_files
  tuple val(sample_id), val(batch), path("*.fasta"), emit: reference_fasta
  tuple val(sample_id), val(batch), path("*.fasta.fai"), emit: reference_fasta_index
  tuple val(sample_id), val(batch), path("*.dict"), emit: gatk_dictionary

  """
  # Check species assignment from Kraken2
  if [[ "${posadasii_bam}" == *"mock"* ]]
  then

    # Kraken2 assignment to C. immitis
    cp ${immitis_bam} ${sample_id}_immitis.bam
    cp ${immitis_reference_fasta} immitis.fasta
    cp ${immitis_reference_fasta_index} immitis.fasta.fai
    cp ${immitis_gatk_dictionary} immitis.dict

  elif [[ "${immitis_bam}" == *"mock"* ]]
  then

    # Kraken2 assignment to C. posadasii
    cp ${posadasii_bam} ${sample_id}_posadasii.bam
    cp ${posadasii_reference_fasta} posadasii.fasta
    cp ${posadasii_reference_fasta_index} posadasii.fasta.fai
    cp ${posadasii_gatk_dictionary} posadasii.dict

  else

    # Kraken2 assignment failed. Assigning to best % mapped species
    immitis_mapping=\$(awk '(NR == 5) { print \$0 }' ${immitis_mapping_report} | awk '{ print \$5 }' | sed 's/(//g' | sed 's/%.*//g')
    posadasii_mapping=\$(awk '(NR == 5) { print \$0 }' ${posadasii_mapping_report} | awk '{ print \$5 }' | sed 's/(//g' | sed 's/%.*//g')
    species=\$(echo "\$immitis_mapping \$posadasii_mapping" | awk '{ if(\$1 >= \$2) { print "immitis" } else { print "posadasii" } }')

    # Determine whether to run alignment or not
    if [[ "\${species}" == "immitis" ]]
    then

      cp ${immitis_bam} ${sample_id}_immitis.bam
      cp ${immitis_reference_fasta} immitis.fasta
      cp ${immitis_reference_fasta_index} immitis.fasta.fai
      cp ${immitis_gatk_dictionary} immitis.dict

    else

      cp ${posadasii_bam} ${sample_id}_posadasii.bam
      cp ${posadasii_reference_fasta} posadasii.fasta
      cp ${posadasii_reference_fasta_index} posadasii.fasta.fai
      cp ${posadasii_gatk_dictionary} posadasii.dict

    fi

  fi
  """

}