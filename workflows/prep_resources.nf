
// ----------------Workflow---------------- //

include { GENOMERESOURCES as GENOMERESOURCES_IMMITIS } from '../subworkflows/genome_resources/genome_resources_prep.nf'
include { GENOMERESOURCES as GENOMERESOURCES_POSADASII } from '../subworkflows/genome_resources/genome_resources_prep.nf'
include { SnpeffPrep } from '../modules/snpeff_prep/snpeff_prep.nf'
include { GenerateKraken2DB } from '../modules/kraken/generate_kraken2_db.nf'

workflow RESOURCESPREP {

  // DOWNLOAD AND PREP GENOME REFS -------- //

  GENOMERESOURCES_IMMITIS(params.immitis_assembly_identifier, "immitis")
  GENOMERESOURCES_POSADASII(params.posadasii_assembly_identifier, "posadasii")

  // SNPEFF SETUP ------------------------- //

  SnpeffPrep(GENOMERESOURCES_IMMITIS.out.snpeff_input, GENOMERESOURCES_POSADASII.out.snpeff_input)

  // GENERATE KRAKEN2 DATABASE ------------ //

  GenerateKraken2DB()

}