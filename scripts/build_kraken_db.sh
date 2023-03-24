#!/bin/bash

# This script builds a database with standard (archaea, bacteria, viral, plasmid, human, UniVec_Core) and EuPathDB46 

# Make folder for database
kraken_db=${1}
mkdir -p ${kraken_db}
cd ${kraken_db}

# Define threads to use for building the database
threads=${2}

# Download kraken2 taxonomy
kraken2-build --db . --download-taxonomy

# Download EuPathDB46 seqid2taxid.map
wget ftp://ftp.ccb.jhu.edu/pub/data/EuPathDB46/seqid2taxid.map

# Making library folder
mkdir library

# Download EuPathDB46 sequences to library folder
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

		# Modify header
		python ../convert_fasta_headers.py --fasta ${fasta} --seqid2taxid seqid2taxid.map

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

# Download kraken2 standard resources
for db in archaea bacteria fungi viral plasmid human UniVec_Core
do

	if [[ ${db} == "human" ]]
	then

		kraken2-build --db . --download-library ${db} --no-masking

	else

		kraken2-build --db . --download-library ${db}

	fi

done

# Build database
kraken2-build --build --threads ${threads} --db .
