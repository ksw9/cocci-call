process RunRepeatMasker {

  // Mask a fasta genome unisng RepeatMasker
  label 'repeatmasker'

  publishDir "${params.resources_dir}/${species}_refs", mode: "copy", pattern: "*_repeatmasker.fasta"

  input:
  tuple val(species), path(fasta)
  path masking_file

  output:
  tuple val("${species}"), path("${species}_repeatmasker.fasta"), emit: masked_fasta

  """
  RepeatMasker -lib ${masking_file} --norna ${fasta}

  mv ${fasta}.masked ${species}_repeatmasker.fasta
  """

}
