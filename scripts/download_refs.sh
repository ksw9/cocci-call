#!/bin/bash

#################################################################################################
#### Download references & index reference genomes for Coccidioides variant calling pipeline ####
#################################################################################################

# Define local container tool: docker (default), podman, or singularity.
container=${1:-"docker"}

# Define locations
kraken2_db_dir=kraken_db/

# Make necessary directories
mkdir -p ${kraken2_db_dir}

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

# Build Kraken2 database
${container} ${run_command} ${bind_option} ${other_options} ${image} /bin/bash ./build_kraken_db.sh ${kraken2_db_dir} ${SLURM_CPUS_ON_NODE}
