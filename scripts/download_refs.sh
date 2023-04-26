#!/bin/bash

#################################################################################################
#### Download references & index reference genomes for Coccidioides variant calling pipeline ####
#################################################################################################

# Define local container tool: docker (default), podman, or singularity.
container=${1:-"docker"}

# Define locations
res_dir=resources/
immitis_ref_dir=immitis_refs/
immitis_bwa_index_dir=immitis_bwa_index/
immitis_gatk_dictionary_dir=immitis_gatk_dictionary/
posadasii_ref_dir=posadasii_refs/
posadasii_bwa_index_dir=posadasii_bwa_index/
posadasii_gatk_dictionary_dir=posadasii_gatk_dictionary/
kraken2_db_dir=kraken_db/
input_dir=input/

# Define main variables for Coccidioides immitis and posadasii
snpeff_url=https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
immitis_ncbi_id="GCF_000149335.2"
immitis_ref_name="${immitis_ncbi_id}.fasta"
immitis_gatk_dictionary_name="${immitis_ncbi_id}.dict"
posadasii_ref_version="GCA_018416015.2"
posadasii_ref_name="${posadasii_ref_version}.fasta"
posadasii_gatk_dictionary_name="${posadasii_ref_version}.dict"

# Define run command and options
if [ "$container" = "docker" ]
then

	run_command="run"
	bind_option="-v $(pwd):/home"
	other_options="--rm"
	image="ksw9/mtb-call"

elif [ "$container" = "podman" ]
then

	run_command="run"
	bind_option="-v $(pwd):/home"
	other_options="--rm"
	image="ksw9/mtb-call"

else

	run_command="exec"
	bind_option="--bind $(pwd)"
	other_options="--cleanenv"
	image="docker://ksw9/mtb-call"
	
fi

# Make necessary directories
mkdir -p ${res_dir}
cd ${res_dir}
mkdir -p ${immitis_ref_dir}
mkdir -p ${immitis_bwa_index_dir}
mkdir -p ${immitis_gatk_dictionary_dir}
mkdir -p ${posadasii_ref_dir}
mkdir -p ${posadasii_bwa_index_dir}
mkdir -p ${posadasii_gatk_dictionary_dir}
#mkdir -p ${kraken2_db_dir} # Will be made by the script invoked at point 3
mkdir -p ${input_dir}

# Requires: entrez-direct (conda), bwa, GATK, samtools, kraken2, all present in the container.
# Run each step within the image. Need to mount local directory so that the resources are downloaded locally, not just in the container.

## 1. Reference genome ##

# Download reference fasta for C. immitis
${container} ${run_command} ${bind_option} ${other_options} ${image} /bin/bash -c "esearch -db nucleotide -query ${immitis_ncbi_id} | efetch -format fasta > ${immitis_ref_dir}${immitis_ref_name}"

# Download reference fasta and GFF for C. posadasii
${container} ${run_command} ${bind_option} ${other_options} ${image} curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/${posadasii_ref_version}/download?include_annotation_type=GENOME_FASTA,GENOME_GFF,RNA_FASTA,CDS_FASTA,PROT_FASTA,SEQUENCE_REPORT&filename=${posadasii_ref_version}.zip" -H "Accept: application/zip"
unzip ${posadasii_ref_version}.zip
mv ncbi_dataset/data/${posadasii_ref_version}/${posadasii_ref_version}*.fna ${posadasii_ref_dir}${posadasii_ref_name}
mv ncbi_dataset/data/${posadasii_ref_version}/*.gff ${posadasii_ref_dir}genes.gff
rm -r README.md ncbi_dataset ${posadasii_ref_version}.zip

# Add fasta sequence to C. posadasii GFF
echo "" >> ${posadasii_ref_dir}genes.gff
echo "##FASTA" >> ${posadasii_ref_dir}genes.gff
cat ${posadasii_ref_dir}${posadasii_ref_name} >> ${posadasii_ref_dir}genes.gff

# bwa index references
${container} ${run_command} ${bind_option} ${other_options} ${image} bwa index ${immitis_ref_dir}${immitis_ref_name}
mv ${immitis_ref_dir}*.{amb,ann,bwt,pac,sa} ${immitis_bwa_index_dir}

${container} ${run_command} ${bind_option} ${other_options} ${image} bwa index ${posadasii_ref_dir}${posadasii_ref_name}
mv ${posadasii_ref_dir}*.{amb,ann,bwt,pac,sa} ${posadasii_bwa_index_dir}

# Create fasta index files
${container} ${run_command} ${bind_option} ${other_options} ${image} samtools faidx ${immitis_ref_dir}${immitis_ref_name}

${container} ${run_command} ${bind_option} ${other_options} ${image} samtools faidx ${posadasii_ref_dir}${posadasii_ref_name}

# create GATK reference dictionaries
${container} ${run_command} ${bind_option} ${other_options} ${image} gatk CreateSequenceDictionary -R ${immitis_ref_dir}${immitis_ref_name}
mv ${immitis_ref_dir}${immitis_gatk_dictionary_name} ${immitis_gatk_dictionary_dir}

${container} ${run_command} ${bind_option} ${other_options} ${image} gatk CreateSequenceDictionary -R ${posadasii_ref_dir}${posadasii_ref_name}
mv ${posadasii_ref_dir}${posadasii_gatk_dictionary_name} ${posadasii_gatk_dictionary_dir}

## 2. Annotation information ##

# Download SnpEff for gene annotation.
${container} ${run_command} ${bind_option} ${other_options} ${image} wget ${snpeff_url}

# Unzip file
unzip snpEff_latest_core.zip
rm snpEff_latest_core.zip

# Download the updated C. immitis annotations
${container} ${run_command} ${bind_option} ${other_options} ${image} java -jar snpEff/snpEff.jar download Coccidioides_immitis_rs_gca_000149335

# Building annotations for C. posadasii
# Step 1. Configure a new genome
echo "" >> snpEff/snpEff.config
echo "# Coccidioides posadasii genome, version ${posadasii_ref_version}" >> snpEff/snpEff.config
echo "${posadasii_ref_version}.genome : Coccidioides_posadasii_${posadasii_ref_version}" >> snpEff/snpEff.config

# Step 2. Building a database from GFF
mkdir -p snpEff/data/${posadasii_ref_version}
cp ${posadasii_ref_dir}genes.gff snpEff/data/${posadasii_ref_version}
${container} ${run_command} ${bind_option} ${other_options} ${image} java -jar snpEff/snpEff.jar build -gff3 -v -noCheckCds -noCheckProtein ${posadasii_ref_version}

## 3. Kraken2 Database ##

# Build Kraken2 database
${container} ${run_command} ${bind_option} ${other_options} ${image} /bin/bash ../build_kraken_db.sh ${kraken2_db_dir} ${SLURM_CPUS_ON_NODE}
