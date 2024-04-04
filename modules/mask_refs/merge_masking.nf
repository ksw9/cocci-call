process MergeMasking {

  // Convert single sample VCF to fasta

  label 'repeatmasker'

  publishDir "${params.resources_dir}/${species_repetmasker}_refs", mode: "copy", pattern: "*_fullymasked.fasta"

  input:
  path scripts_dir
  tuple val(species_repetmasker), path(fasta_repetmasker)
  tuple val(species_nucmer), path(fasta_nucmer)

  output:
  tuple val("${species_repetmasker}"), path("${species_repetmasker}_fullymasked.fasta"), emit: masked_fasta

  """
  python ${scripts_dir}/merge_masking.py --fasta_files "${fasta_repetmasker};${fasta_nucmer}"

  mv merged.fasta ${species_repetmasker}_fullymasked.fasta
  """

}
