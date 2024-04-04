process DownloadRefs {

  // Download assembly fasta and gff from NCBI

  label 'download_refs'

  publishDir "${params.resources_dir}/${species}_refs", mode: "copy", pattern: "*.fasta"

  input:
  val assembly_identifier
  val species

  output:
  tuple val("${species}"), path("${assembly_identifier}.fasta"), emit: fasta
  tuple val("${species}"), path("${assembly_identifier}.gff"), emit: gff

  """
  # Download reference fasta and GFF for C. posadasii
  curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/${assembly_identifier}/download?include_annotation_type=GENOME_FASTA,GENOME_GFF&filename=${assembly_identifier}.zip" -H "Accept: application/zip"
  
  unzip ${assembly_identifier}.zip
  
  mv ncbi_dataset/data/${assembly_identifier}/*.fna ${assembly_identifier}.fasta
  mv ncbi_dataset/data/${assembly_identifier}/*.gff ${assembly_identifier}.gff
  """

}
