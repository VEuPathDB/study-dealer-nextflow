library(tidyverse)
library(study.wrangler)


wrangle <- function() {
  rm(list = ls())

  # read in file
  genePhenotype = entity_from_file("phenotype_results.txt")

  # Remove ImageName and Orthologs columns
  genePhenotype <- genePhenotype %>% modify_data(select(-Imagename, -Orthologs))

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
    set_variable_metadata('Temperature', display_order=3, display_name = "Temperature", definition = "Temperature") %>%
    set_variable_metadata('Timepoint', display_order=4, display_name = "Timepoint", definition = "Timepoint") %>%
    set_variable_metadata('Media', display_order=5, display_name = "Media", definition = "Media")

  study = study("afumA1163_TF_KO_Image_NAFeaturePhenotypeImage", genePhenotype)

  return(study)

}
