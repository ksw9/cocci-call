# *Coccidioides* species and variant identification pipeline

Pipeline for *Coccidioides* *immitis* and *posadasii* species and variant identification from short-read data for epidemiology and phylogenetics. Briefly, this pipeline takes raw short-read data, pre-processes reads (read adapter trimming and taxonomic filtering), maps them to a provided reference genome, calls variants, and outputs consensus FASTA sequences for downstream applications.

The pipeline is designed to be run on High Performance Computing (HPC) systems using the Slurm Workload Manager. Tasks are parallelized thanks to the Nextflow workflow system and the pipeline uses Docker containers to enhance reproducibility and portability.

/// ---------------------------------------- ///

## Overview:

### RESOURCESPREP workflow

Run this workflow by selecting the 'download_refs' profile (see 'RUN COMMAND' section below).

* Reference FASTA and GFF download

* Masking (optional) with RepeatMasker, or Nucmer, or both

* FASTA indexing

* BWA index generation

* GATK dictionary generation

* snpEff download and databases generation

* Kraken2 database generation (bacteria, fungi, viral, Homo sapiens, and EuPathDB46 genomes)

### VARIANTCALLING workflow

Run this workflow by selecting the 'variant_calling' profile (see 'RUN COMMAND' section below).

* Reads trimming and QC with TrimGalore

* Kraken2 filtering and species determination (*immitis* vs *posadasii*) based on number of unique minimizers

* Mapping with BWA

* Marking duplicates with GATK

* Variant calling with GATK and (optionally) LoFreq

* Annotation of vcf files with snpEff

* Converting vcf to consensus FASTA

* Summarizing of run info

Note that species determination is determined by a user-defined threshold for the log2 ratio of unique minimizers for *immitis* and *posadasii*. If species could not be determined at the Kraken2 step (*e.g.* threshold set too high), then the species is assigned at the BWA step, as the one for which mapping percentage is higher.

/// ---------------------------------------- ///

## SETUP:

1. Clone Github repo.

```
git clone https://github.com/ksw9/cocci-call.git
```

2. Make sure that Nextflow and your container platform of choice (Docker, Podman, or Singularity) are installed. If using Lmod, load the necessary modules: *e.g.*

```
module load docker java nextflow 
```

3. Unzip the masking files in masking_files folder.

```
gzip -d masking_files/*.fa.gz 
```

4. Run RESOURCESPREP workflow to populate your resources directory with required references and databases. Note that Kraken2 database generation can take days, so make sure to submit a Slurm job with appropriate quality of service specifications.

```
cd cocci-call
nextflow run main.nf -profile singularity,download_refs --repeatmasker_mask true --nucmer_mask true --resources_dir "$(pwd)/resources"
```

5. Modify the config file (nextflow.config):

* Update resources_dir (full path to directory resources)

* Verify resources paths (relative to resources_dir)

* Set clusterOptions parameters according to your HPC settings

/// ---------------------------------------- ///

## INPUT:

The variant calling workflow reads samples' info from a tab-delimited file, placed in "/path/to/resources_dir/input"
The file has the following columns:

* **sample**: unique sample identifier

* **fastq_1**: full path to fastq mate 1

* **fastq_2**: full path to fastq mate 2

* **batch**: batch name

/// ---------------------------------------- ///

## OUTPUTS:

### RESOURCESPREP workflow

A folder named "resources" is created by default in the project directory. The directory has the following structure:

<pre>
<b>resources</b>
│
├── <b>immitis_bwa_index</b>
│
├── <b>immitis_gatk_dictionary</b>
│
├── <b>immitis_refs</b>
│
├── <b>input</b>
│
├── <b>kraken_db</b>
│
├── <b>posadasii_bwa_index</b>
│
├── <b>posadasii_gatk_dictionary</b>
│
├── <b>posadasii_refs</b>
│
└── <b>snpEff</b>
</pre>

### VARIANTCALLING workflow

All outputs are stored in the results directory, within the project directory. Directory structure mirrors the input reads file, with directories organized by batch, then sample.

<pre>
<b>results</b>
│
└── <b>batch_0</b>
    │
    ├── <b>sample</b>
    |   │
    |   ├── <b>bams</b>
    |   │
    |   ├── <b>fasta</b>
    |   │
    |   ├── <b>kraken</b>
    |   │
    |   ├── <b>stats</b>
    |   │
    |   ├── <b>trim</b>
    |   │
    |   └── <b>vars</b>
    |
    └── <b>batch_n</b>
</pre>

/// ---------------------------------------- ///

## RUN COMMAND:

```
nextflow run main.nf -profile [PROFILES] [OPTIONS]
```

### PROFILES:

* **standard**: runs the pipeline using Docker

* **docker**: runs the pipeline using Docker

* **podman**: runs the pipeline using Podman

* **singularity**: runs the pipeline using Singularity

* **download_refs**: runs the RESOURCESPREP workflow, which populates the "resources" directory with the necessary files for the analysis workflow

* **variant_calling**: runs the VARIANTCALLING workflow

### OPTIONS:

* **repeatmasker_mask**: set to **true** if you want to mask the genome with RepeatMasker, **false** otherwise

* **nucmer_mask**: set to **true** if you want to mask the genome with RepeatMasker, **false** otherwise

* **minimizers_log_ratio_thr**: threshold for the log ratio of Kraken 2 *C. immitis* and *C. posadasii* minimizers

* **run_lofreq**: set to **true** if you want to run LoFreq, **false** otherwise. Note that GATK is run anyway

* **seq_platform**: name of sequencing platform to be added to reads by picard AddOrReplaceReadGroups prior to GATK processing

* **library**: name of library to be added to reads by picard AddOrReplaceReadGroups prior to GATK processing

* **vcf_filter**: filter for vcf files, based on any header column. Parameters need to be specified for variants to be kept, *e.g.* 'QUAL > 20' will result in any variant with quality below 20 to be discarded. Multiple parameters can be combined by "&&" (*i.e.* *and*) or "||" (*i.e.* *or*) operators, *e.g.* 'QD > 2.0 && FS < 60.0 && MQ > 40.0 && DP > 10 && GQ > 50 && QUAL > 20 && QUAL != "." && RGQ > 20 && (TYPE == "SNP" || TYPE == "indel" || TYPE != "insert")'

* **ploidy**: Ploidy to be used by GATK

* **nextseq**: set to **true** if data comes from an Illumina NextSeq system, **false** otherwise

* **nextseq_qual_threshold**: read quality threshold for Illumina NextSeq data

* **variants_only**: set to **true** if GATK HaplotypeCaller should only output variant sites, **false** to output also invariant ones

/// ---------------------------------------- ///

## DEPENDENCIES:

* Nextflow 23.10.1+

* Container platform, one of
    * Docker 20.10.21+
    * Podman
    * Singularity 3.8.5+

* Slurm

/// ---------------------------------------- ///

## DAG

### RESOURCESPREP workflow

<p align="center">
  <img width="1585" height="273" src="https://github.com/ksw9/cocci-call/blob/main/images/pipeline_dag_resourcesprep.png">
</p>

### VARIANTCALLING workflow

<p align="center">
  <img width="731" height="494" src="https://github.com/ksw9/cocci-call/blob/main/images/pipeline_dag_variantcalling.png">
</p>
