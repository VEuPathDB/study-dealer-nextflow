library(tidyverse)
library(study.wrangler)


wrangle <- function() {
  rm(list = ls())

  # read in file
  genePhenotype = entity_from_file("phenotype_results.txt")

  # Skip Sexual_development column
  genePhenotype <- genePhenotype %>% modify_data(select(-Sexual_Development))

  # Set meta data for entity
 genePhenotype <- genePhenotype %>% set_entity_metadata(name = "genePhenotypeImageData", display_name = "Phenotype Image Collection", stable_id="genePhenotypeImageCollection", display_name_plural="Gene Phenotype Image Collection")

   # make gene column a variable
 genePhenotype <- genePhenotype %>% redetect_columns_as_variables('gene')
 
  #default column/variable labels
 genePhenotype <- genePhenotype %>%  set_variable_display_names_from_provider_labels()

  #  deal with the primary Key (gene variable). boilerplate
 genePhenotype <- genePhenotype %>%
   modify_data(mutate(ID = row_number())) %>%
   sync_variable_metadata() %>%
   redetect_column_as_id('ID')  %>%
   redetect_columns_as_variables(c('Image'))

  genePhenotype <- genePhenotype %>%
    set_variable_metadata('gene', display_name = "Gene", provider_label=list("gene"), display_order=1, hidden=list('variableTree')) %>%
    set_variable_metadata('Image', display_order=2, display_name = "Image", definition = "Image") %>%
    set_variable_metadata('location', display_order=3, display_name = "location", definition = "location") %>%
    set_variable_metadata('media', display_order=4, display_name = "media", definition = "media") %>%
    set_variable_metadata('morphology', display_order=5, display_name = "morphology", definition = "morphology") %>%
    set_variable_metadata('physiology', display_order=6, display_name = "physiology", definition = "physiology") %>%
    set_variable_metadata('temperature', display_order=7, display_name = "temperature", definition = "temperature") %>%
    set_variable_metadata('timepoint', display_order=8, display_name = "timepoint", definition = "timepoint")

  study = study("ncraOR74A_phenotype_GeneImage_NAFeaturePhenotypeImage", genePhenotype)

  return(study)

}
