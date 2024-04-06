process DownloadNCBI {
	
  // Download NCBI fasta files for bacteria, fungi, viral, human complete genomes
  label 'kraken2'

  errorStrategy 'ignore'

  input:
  tuple val(group), val(taxid), val(target_file), val(ftp_link)

  output:
  path "*.fna", optional: true, emit: fasta

  """
  # Test that the link is good
  wget_test=\$(wget -S --spider "${ftp_link}" 2>&1)

  # Download if good link
  if [[ \${wget_test} == *"Remote file exists."* ]]
  then

    wget -w 1 --tries 20 --retry-connrefused --retry-on-host-error "${ftp_link}"

    # Uncompress file
    gzip -d *.gz

    # Modify fasta header
    mv ${target_file} temp.fna
    awk -v taxid=${taxid} '{ if(\$0 ~ /^>.*/) { print \$1"|kraken:taxid|"taxid } else { print \$0 } }' temp.fna > ${target_file}
    rm temp.fna

    # Add prefix "nomasking_" if Human fasta
    if [[ "${group}" == "vertebrate_mammalian/Homo_sapiens/" ]]
    then

      mv ${target_file} nomasking_${target_file}

    fi

  fi
  """

}
