library(tidyverse)
library(study.wrangler)


wrangle <- function() {
  rm(list = ls())
  
  # read in file
  genePhenotype = entity_from_file("phenotype_with_names.txt")
  
  # Set meta data for entity
  genePhenotype <- genePhenotype %>% set_entity_metadata(name = "genePhenotypeData", display_name = "Gene Phenotype Data", stable_id="GENE_PHENOTYPE_DATA_ENTITY", display_name_plural="Gene Phenotype Data")
  
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
    set_variable_metadata('gene', stable_id = "VEUPATHDB_GENE_ID", display_name = "Gene", provider_label=list("gene"), display_order=1, hidden=list('variableTree')) %>%
    set_variable_metadata('pubmedId', display_order=2, display_name = "Pubmed ID", definition = "Pubmed ID", data_type = "string", data_shape = "categorical") %>%
    set_variable_metadata('modType', display_order=3, display_name = "Mod Type", definition = "Mod Type") %>%
    set_variable_metadata('isSuccessful', display_order=4, display_name = "Is Successful", definition = "Is Successful") %>%
    set_variable_metadata('qualityTerm', display_order=5, display_name = "Quality Term ID", definition = "Quality Term ID", hidden=list('variableTree')) %>%
    set_variable_metadata('entityTerm', display_order=6, display_name = "Entity Term ID", definition = "Entity Term ID", hidden=list('variableTree')) %>%
    set_variable_metadata('timing', display_order=7, display_name = "Timing", definition = "Timing")  %>%
    set_variable_metadata('lifeCycleTerm', display_order=8, display_name = "Life Cycle Term ID", definition = "Life Cycle Term ID", hidden=list('variableTree'))  %>%
    set_variable_metadata('evidenceTerm', display_order=9, display_name = "Evidence Term ID", definition = "Evidence Term ID", hidden=list('variableTree'))  %>%
    set_variable_metadata('note', display_order=10, display_name = "Note", definition = "Note") %>%
    set_variable_metadata('qualityTerm_name', display_order=5, display_name = "Quality Term", definition = "Quality Term") %>%
    set_variable_metadata('entityTerm_name', display_order=6, display_name = "Entity Term", definition = "Entity Term") %>%
    set_variable_metadata('lifeCycleTerm_name', display_order=8, display_name = "Life Cycle Term", definition = "Life Cycle Term")  %>%
    set_variable_metadata('evidenceTerm_name', display_order=9, display_name = "Evidence Term", definition = "Evidence Term")  %>%
    modify_data(mutate(pubmedId = as.character(pubmedId)))
  
  
  study = study(name="TEMP_STUDY_NAME", genePhenotype)

  return(study)
  
}
