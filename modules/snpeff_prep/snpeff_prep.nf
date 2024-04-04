process SnpeffPrep {

  // Prep databases for SnpEff

  label 'download_refs'

  publishDir "${params.resources_dir}", mode: "copy", pattern: "snpEff"

  input:
  tuple val(immitis), path(gff_immitis)
  tuple val(posadasii), path(gff_posadasii)

  output:
  path "snpEff", emit: snpeff

  """
  # Download SnpEff for gene annotation.
  wget ${params.snpeff_url}

  # Unzip file
  unzip snpEff_latest_core.zip

  ### Building annotations for C. immitis
  # Step 1. Configure a new genome
  echo "" >> snpEff/snpEff.config
  echo "# Coccidioides immitis genome, version ${params.immitis_assembly_identifier}" >> snpEff/snpEff.config
  echo "${params.immitis_assembly_identifier}.genome : Coccidioides_immitis_${params.immitis_assembly_identifier}" >> snpEff/snpEff.config

  # Step 2. Building a database from GFF
  mkdir -p snpEff/data/${params.immitis_assembly_identifier}
  cp ${gff_immitis} snpEff/data/${params.immitis_assembly_identifier}/genes.gff
  java -jar snpEff/snpEff.jar build -gff3 -v -noCheckCds -noCheckProtein ${params.immitis_assembly_identifier}

  ### Building annotations for C. posadasii
  # Step 1. Configure a new genome
  echo "" >> snpEff/snpEff.config
  echo "# Coccidioides posadasii genome, version ${params.posadasii_assembly_identifier}" >> snpEff/snpEff.config
  echo "${params.posadasii_assembly_identifier}.genome : Coccidioides_posadasii_${params.posadasii_assembly_identifier}" >> snpEff/snpEff.config

  # Step 2. Building a database from GFF
  mkdir -p snpEff/data/${params.posadasii_assembly_identifier}
  cp ${gff_posadasii} snpEff/data/${params.posadasii_assembly_identifier}/genes.gff
  java -jar snpEff/snpEff.jar build -gff3 -v -noCheckCds -noCheckProtein ${params.posadasii_assembly_identifier}
  """

}
