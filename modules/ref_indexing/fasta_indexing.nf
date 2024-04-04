process FastaIndex {

  // Index a fasta file

  label 'mapping'

  publishDir "${params.resources_dir}/${species}_refs", mode: "copy", pattern: "*.fai"

  input:
  tuple val(species), path(fasta)

  output:
  tuple val("${species}"), path("${fasta}.fai"), emit: masked_fasta_index

  """
  samtools faidx ${fasta}
  """

}
