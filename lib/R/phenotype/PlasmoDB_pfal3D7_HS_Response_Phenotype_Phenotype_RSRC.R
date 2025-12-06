library(tidyverse)
library(study.wrangler)


wrangle <- function() {
  rm(list = ls())

  # read in file
  genePhenotype = entity_from_file("fc_growth_fc_hs_pfs.txt")

  # Set meta data for entity
 genePhenotype <- genePhenotype %>% set_entity_metadata(name = "genePhenotypeData", display_name = "Gene Phenotype Data", stable_id="genePhenotypeData", display_name_plural="Gene Phenotype Data")

   # make gene column a variable
 genePhenotype <- genePhenotype %>% redetect_columns_as_variables('gene')
 
  #default column/variable labels
 genePhenotype <- genePhenotype %>%  set_variable_display_names_from_provider_labels()

  #  deal with the primary Key (gene variable). boilerplate
 genePhenotype <- genePhenotype %>%
   modify_data(mutate(ID = row_number())) %>%
   sync_variable_metadata() %>%
   redetect_column_as_id('ID')

  genePhenotype <- genePhenotype %>%
    set_variable_metadata('gene', display_name = "Gene", provider_label=list("gene"), display_order=1, hidden=list('variableTree'))

  study = study("pfal3D7_phenotype_HS_Response_Phenotype", genePhenotype)

  return(study)

}
