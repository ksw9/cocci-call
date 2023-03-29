#!/bin/bash

# This script builds a Kraken2 database with bacteria, fungi, viral, human complete genomes and EuPathDB46 (eukaryotic pathogen genomes)

# Make folder for database
kraken_db=${1}
mkdir -p ${kraken_db}
cd ${kraken_db}

# Define threads to use for building the database
threads=${2}

# Check if manual download is desired for Kraken2 standard resources
download_type=${3:-"kraken_implementation"}

# Download kraken2 taxonomy
kraken2-build --db . --download-taxonomy

# Download kraken2 standard resources
if [[ ${download_type} == "kraken_implementation" ]]
then

	# Download via Kraken2 --download-library (N.B. using ftp since rsync seems to have issues with some databases)
	for db in bacteria fungi viral human
	do

		kraken2-build --use-ftp --db . --download-library ${db}

	done

else

	# Download genome assemblies, then change headers for compatibility with Kraken2
	for group in bacteria fungi viral "vertebrate_mammalian/Homo_sapiens/"
	do

		# Download assemblies list
		wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/${group}/assembly_summary.txt

		# Create arrays of taxids and ftp links and other info
		taxid_arr=($(awk -F'\t' -v header_row=2 -v taxid_column=6 '(NR > header_row) { print $taxid_column }' assembly_summary.txt))
		assembly_accession=($(awk -F'\t' -v header_row=2 -v accession_column=1 '(NR > header_row) { print $accession_column }' assembly_summary.txt))
		assembly_level=($(awk -F'\t' -v header_row=2 -v level_column=12 '(NR > header_row) { print $level_column }' assembly_summary.txt | sed 's/ /_/g')) # N.B. Spaces are removed or they will break up lines into multiple variables
		asm_name=($(awk -F'\t' -v header_row=2 -v name_column=16 '(NR > header_row) { print $name_column }' assembly_summary.txt | sed 's/ /_/g')) # Sometimes it's the species name, so there's a space which needs to be removed
		ftp_arr=($(awk -F'\t' -v header_row=2 -v ftp_column=20 '(NR > header_row) { print $ftp_column }' assembly_summary.txt))

		# Download assemblies and change headers to be Kraken2 compatible
		for index in "${!taxid_arr[@]}"
		do
    	
    		# N.B. only complete genomes are downloaded (too many otherwise)
    		if [[ ${assembly_level[index]} == "Complete_Genome" ]]
    		then

    			# Download assembly
    			target_file="${assembly_accession[index]}_${asm_name[index]}_genomic.fna"
	    		ftp_link="${ftp_arr[index]}/${target_file}.gz"
    			wget ${ftp_link}

    			# Uncompress file
    			gzip -d ${target_file}

	    		# Modify fasta header
    			mv ${target_file} temp.fna
    			awk -v taxid=${taxid_arr[index]} '{ if($0 ~ /^>.*/) { print $1"|kraken:taxid|"taxid } else { print $0 } }' temp.fna > ${target_file}
    			rm temp.fna

	    		# Add sequence to library
    			if [[ ${group} == "vertebrate_mammalian/Homo_sapiens/" ]]
				then

					kraken2-build --db . --add-to-library ${target_file} --no-masking

				else
			
					kraken2-build --db . --add-to-library ${target_file}

				fi

				# Remove assembly fasta file (useless now that it was added to K2 library)
				rm ${target_file}

			fi

		done

		# Remove assembly_summary.txt
		rm assembly_summary.txt

	done

fi

# Download EuPathDB46 seqid2taxid.map
wget ftp://ftp.ccb.jhu.edu/pub/data/EuPathDB46/seqid2taxid.map

# Download EuPathDB46 (eukaryotic pathogen genomes) sequences to library folder
for db in AmoebaDB46 CryptoDB46 FungiDB46 GiardiaDB46 MicrosporidiaDB46 PiroplasmaDB46 PlasmoDB46 ToxoDB46 TrichDB46 TriTrypDB46
do

	# Download sequences
	wget ftp://ftp.ccb.jhu.edu/pub/data/EuPathDB46/${db}.tgz
	tar -xzvf ${db}.tgz
	mv $(basename -s '46' ${db}) ${db} # Databases are named with '46' in gzipped file, but once uncompressed, folder is missing '46'
	rm ${db}.tgz

	# Make headers compatible with Kraken2, then add sequence to library
	for fasta in ${db}/*fna
	do

		# Get the first fasta sequence header so as to find the taxid
		sample_header=$(grep -E '^>.*' ${fasta} | head -1)

		# Find taxid for species
		if [[ "${sample_header}" == *"kraken:taxid|"* ]]
		then

			taxid=$(echo ${sample_header} | sed 's/>//g' | grep -f - seqid2taxid.map | head -1 | awk '{ print $2 }')

		else

			taxid=$(echo ${sample_header} | sed 's/>//g' | sed 's/ //g' | sed 's/|//g' | grep -f - seqid2taxid.map | head -1 | awk '{ print $2 }')

		fi

		# Modify sequence headers to be compatible with Kraken2
		awk -v taxid=${taxid} '{ if($0 ~ /^>.*/ && $0 !~ />kraken.*/) { print $1"|kraken:taxid|"taxid } else { print $0 } }' ${fasta} > corrected_fasta.fna

		# Replace fasta
		mv corrected_fasta.fna ${fasta}

		# Add sequence to library
		kraken2-build --db . --add-to-library ${fasta}

	done

	# Removed db folder
	rm -r ${db}

done

# Remove EuPathDB46 seqid2taxid.map (was only needed to convert sequence names to be Kraken2 compatible)
rm seqid2taxid.map

# Build database
kraken2-build --build --threads ${threads} --db .

# Clean up
kraken2-build --clean --db .
