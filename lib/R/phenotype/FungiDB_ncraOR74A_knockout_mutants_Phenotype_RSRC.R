library(tidyverse)
library(study.wrangler)

wrangle <- function() {
  rm(list = ls())

  # read in file
  genePhenotype = entity_from_file("ko_mutants_final_09_12_18.txt")

  # Set meta data for entity
 genePhenotype <- genePhenotype %>% set_entity_metadata(name = "genePhenotypeData", display_name = "Gene Phenotype Data", stable_id="genePhenotypeData", display_name_plural="Gene Phenotype Data")

  #default column/variable labels
 genePhenotype <- genePhenotype %>%  set_variable_display_names_from_provider_labels()

  #  deal with the primary Key (gene variable). boilerplate
  genePhenotype <- genePhenotype %>%
    modify_data(mutate(ID = row_number())) %>%
    sync_variable_metadata() %>%
    redetect_column_as_id('ID')


  genePhenotype <- genePhenotype %>%
    set_variable_metadata('gene', display_name = "Gene", provider_label=list("gene"), display_order=1, hidden=list('variableTree')) %>%
    set_variable_metadata('pubmed', display_name = "Pubmed", display_order=2, definition = "Pubmed", data_type = "string", data_shape = "categorical") %>%
    set_variable_metadata('gene_classification', display_name = "Gene Classification", display_order=3,definition = "Gene Classification", hidden=list('variableTree')) %>%
    set_variable_metadata('fgsc', display_name = "FGSC ID", display_order=4,definition = "FGSC ID", hidden=list('variableTree')) %>%
    set_variable_metadata('gene_name', display_name = "Gene Name", display_order=5, definition = "Gene Name") %>%
    set_variable_metadata('mating_type', display_name = "Mating Type", display_order=6,definition = "Mating Type") %>%
    set_variable_metadata('basal_hyphae_growth_rate', display_name = "Basal Hyphae Growth Rate", display_order=7,definition = "Basal Hyphae Growth Rate") %>%
    set_variable_metadata('aerial_hyphae_height', display_name = "Aerial Hyphae Height", display_order=8, definition = "Aerial Hyphae Height") %>%
    set_variable_metadata('conidia_number', display_name = "Conidia Number", display_order=9, definition = "Conidia Number") %>%
    set_variable_metadata('conidia_morphology', display_name = "Conidia Morphology", display_order=10, definition = "Conidia Morphology") %>%
    set_variable_metadata('protoperithecia_number', display_name = "Protoperithecia Number", display_order=11, definition = "Protoperithecia Number") %>%
    set_variable_metadata('protoperithecial_morphology', display_name = "Protoperithecia Morphology", display_order=12, definition = "Protoperithecia Morphology") %>%
    set_variable_metadata('perithecia_number', display_name = "Perithecia Number", display_order=13, definition = "Perithecia Number") %>%
    set_variable_metadata('perithecia_morphology', display_name = "Perithecia Morphology", display_order=14, definition = "Perithecia Morphology") %>%
    set_variable_metadata('ascospore_number', display_name = "Ascospore Number", display_order=15, definition = "Ascospore Number") %>%
    set_variable_metadata('ascospore_morphology', display_name = "Ascospore Morphology", display_order=16, definition = "Ascospore Morphology") %>%
    modify_data(mutate(pubmed = as.character(pubmed)))
  
  study = study(name="TEMP_STUDY_NAME", genePhenotype)

  return(study)

}
