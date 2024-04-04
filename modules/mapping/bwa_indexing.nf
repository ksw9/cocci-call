process BwaIndex {

  // Generate a BWA index

  label 'mapping'

  publishDir "${params.resources_dir}/${species}_bwa_index", mode: "copy", pattern: "*.{amb,ann,bwt,pac,sa}"

  input:
  tuple val(species), path(fasta), path(fasta_index)

  output:
  tuple val("${species}"), path("*.amb"), path("*.ann"), path("*.bwt"), path("*.pac"), path("*.sa"), emit: bwa_index

  """
  bwa index ${fasta}
  """

}
