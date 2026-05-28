library(tidyverse)
library(study.wrangler)

# Only set variable metadata if the column exists (some cols may be dropped as empty)
set_var_meta_if_present <- function(entity, col, ...) {
  if (col %in% names(entity@data)) {
    set_variable_metadata(entity, col, ...)
  } else {
    entity
  }
}

wrangle_phi_base <- function(filename) {
  genePhenotype = entity_from_file(filename)

  genePhenotype <- genePhenotype %>%
    set_entity_metadata(
      name = "genePhenotypeData",
      display_name = "Gene Phenotype Data",
      stable_id = "GENE_PHENOTYPE_DATA_ENTITY",
      display_name_plural = "Gene Phenotype Data"
    )

  genePhenotype <- genePhenotype %>%
    set_variable_display_names_from_provider_labels()

  genePhenotype <- genePhenotype %>%
    modify_data(mutate(ID = row_number())) %>%
    sync_variable_metadata() %>%
    redetect_column_as_id('ID') %>%
    redetect_columns_as_variables('phenotype_id')

  genePhenotype <- genePhenotype %>%
    set_variable_metadata('phenotype_id', display_order = 1, display_name = "phenotype_id", hidden=list('variableTree')) %>%
    set_variable_metadata('gene', stable_id = "VEUPATHDB_GENE_ID", display_order = 2, display_name = "Gene", definition = "Gene", hidden=list('variableTree')) %>%
    set_var_meta_if_present('PHI.base.entry', display_order = 3, display_name = "PHI-base Entry", definition = "PHI.base.entry") %>%
    set_var_meta_if_present('Essential.gene', display_order = 4, display_name = "Essential Gene", definition = "Essential.gene") %>%
    set_var_meta_if_present('Multiple.mutations', display_order = 5, display_name = "Multiple Mutations", definition = "Multiple.mutations") %>%
    set_var_meta_if_present('Pathogen.species', display_order = 6, display_name = "Pathogen Species", definition = "Pathogen.species") %>%
    set_var_meta_if_present('Pathogen.strain', display_order = 7, display_name = "Pathogen Strain", definition = "Pathogen.strain") %>%
    set_var_meta_if_present('Host.species', display_order = 8, display_name = "Host Species", definition = "Host.species") %>%
    set_var_meta_if_present('Host.strain', display_order = 9, display_name = "Host Strain", definition = "Host.strain") %>%
    set_var_meta_if_present('Tissue', display_order = 10, display_name = "Tissue", definition = "Tissue") %>%
    set_var_meta_if_present('Mutant.Phenotype', display_order = 11, display_name = "Mutant Phenotype", definition = "Mutant.Phenotype") %>%
    set_var_meta_if_present('Disease', display_order = 12, display_name = "Disease", definition = "Disease") %>%
    set_var_meta_if_present('Mating.defect', display_order = 13, display_name = "Mating Defect", definition = "Mating.defect") %>%
    set_var_meta_if_present('Prepenetration.defect', display_order = 14, display_name = "Prepenetration Defect", definition = "Prepenetration.defect") %>%
    set_var_meta_if_present('Penetration.defect', display_order = 15, display_name = "Penetration Defect", definition = "Penetration.defect") %>%
    set_var_meta_if_present('Postpenetration.defect', display_order = 16, display_name = "Postpenetration Defect", definition = "Postpenetration.defect") %>%
    set_var_meta_if_present('Disease.manifestation', display_order = 17, display_name = "Disease Manifestation", definition = "Disease.manifestation") %>%
    set_var_meta_if_present('Vegetative.spores', display_order = 18, display_name = "Vegetative Spores", definition = "Vegetative.spores") %>%
    set_var_meta_if_present('Sexual.spores', display_order = 19, display_name = "Sexual Spores", definition = "Sexual.spores") %>%
    set_var_meta_if_present('Invitro.growth', display_order = 20, display_name = "Invitro Growth", definition = "Invitro.growth") %>%
    set_var_meta_if_present('Spore.germination', display_order = 21, display_name = "Spore Germination", definition = "Spore.germination") %>%
    set_var_meta_if_present('Gene.inducer', display_order = 22, display_name = "Gene Inducer", definition = "Gene.inducer") %>%
    set_var_meta_if_present('Experimental.technique', display_order = 23, display_name = "Experimental Technique", definition = "Experimental.technique") %>%
    set_var_meta_if_present('PMID', display_order = 24, display_name = "PMID", definition = "PMID") %>%
    set_var_meta_if_present('Comments', display_order = 25, display_name = "Comments", definition = "Comments")

  study = study(name="TEMP_STUDY_NAME", genePhenotype)

  return(study)
}
