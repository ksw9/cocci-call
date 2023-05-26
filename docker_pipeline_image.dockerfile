FROM continuumio/miniconda3:4.12.0

### UPDATING CONDA ------------------------- ###

RUN conda update -y conda

### INSTALLING PIPELINE PACKAGES ----------- ###

# Adding bioconda to the list of channels
RUN conda config --add channels bioconda

# Adding conda-forge to the list of channels
RUN conda config --add channels conda-forge

# Installing mamba
RUN conda install -y mamba

# Installing packages
RUN mamba install -y \
    bbmap=39.01 \
    bcftools=1.12 \
    bwa=0.7.17 \
    cutadapt=4.4 \
    entrez-direct=16.2 \
    fastqc=0.12.1 \
    gatk4=4.3.0.0 \
    kraken2=2.1.2 \
    lofreq=2.1.5 \
    mummer=3.23 \
    mykrobe=0.12.1 \
    picard=2.27.5 \
    sambamba=0.8.1 \
    seqtk=1.3 \
    tabix=1.11 \
    trim-galore=0.6.10 && \
    conda clean -afty

RUN mamba install --force-reinstall -y java-jdk # Needed to fix "java: symbol lookup error: java: undefined symbol: JLI_StringDup" error

### FIX KRAKEN2-BUILD ---------------------- ###

# Fixes the "rsync_from_ncbi.pl: unexpected FTP path (new server?)" error
# Thanks to Bent Petersen, PhD (https://www.bpetersen.dk/post/kraken2-rsync_from_ncbi-pl-unexpected-ftp-path-new-server-for)
#RUN mv /opt/conda/libexec/rsync_from_ncbi.pl /opt/conda/libexec/rsync_from_ncbi.pl.bak && \
#    sed '46 s/ftp/https/' /opt/conda/libexec/rsync_from_ncbi.pl.bak > /opt/conda/libexec/rsync_from_ncbi.pl && \
#    chmod 775 /opt/conda/libexec/rsync_from_ncbi.pl

### SETTING WORKING ENVIRONMENT ------------ ###

# Set workdir to /home/
WORKDIR /home/

# Launch bash automatically
CMD ["/bin/bash"]
