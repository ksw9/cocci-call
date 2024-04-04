process SnpeffInputPrep {

  // Prep input gff for SnpEff database creation

  label 'download_refs'

  input:
  tuple val(species), path(fasta)
  tuple val(species_copy), path(gff)

  output:
  tuple val("${species}"), path("${species}_genes.gff"), emit: snpeff_input

  """
  # Add fasta sequence to gff
  cp ${gff} ${species}_genes.gff
  echo "" >> ${species}_genes.gff
  echo "##FASTA" >> ${species}_genes.gff
  cat ${fasta} >> ${species}_genes.gff
  """

}