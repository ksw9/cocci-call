process RunNucmerMasker {

  // Mask a fasta genome using Nucmer

  label 'nucmer'

  publishDir "${params.resources_dir}/${species}_refs", mode: "copy", pattern: "*_nucmer.fasta"

  input:
  path scripts_dir
  tuple val(species), path(fasta)

  output:
  tuple val("${species}"), path("${species}_nucmer.fasta"), emit: masked_fasta

  """
  nucmer --maxmatch --nosimplify --prefix=${species} ${fasta} ${fasta}

  show-coords -r ${species}.delta > ${species}.coords
  
  python ${scripts_dir}/coords_masking.py --fasta_file ${fasta} --coords_file ${species}.coords

  mv masked.fasta ${species}_nucmer.fasta
  """

}
