library(tidyverse)
library(study.wrangler)

wrangle <- function() {
  rm(list = ls())
  
  # read in file
  genePhenotype = entity_from_file("phenotype_growth_rates.txt")
  
  # Set meta data for entity
  genePhenotype <- genePhenotype %>% set_entity_metadata(name = "genePhenotypeData", display_name = "Gene Phenotype Data", stable_id="genePhenotypeData", display_name_plural="Gene Phenotype Data")
  
  # default column/variable labels
  genePhenotype <- genePhenotype %>%  set_variable_display_names_from_provider_labels()
  
  # make gene column a variable
  genePhenotype <- genePhenotype %>% redetect_columns_as_variables(c("gene", "construct"))
  
  #  deal with the primary Key (gene variable). boilerplate
  genePhenotype <- genePhenotype %>%
    modify_data(mutate(ID = row_number())) %>%
    sync_variable_metadata() %>%
    redetect_column_as_id('ID')
  
  
  # update columns as needed
  genePhenotype <- genePhenotype %>%
    set_variable_metadata('gene', display_name = "Gene", provider_label=list("gene"), display_order=1, hidden=list('variableTree'))  %>%
    set_variable_metadata('phenotype', display_order=2, display_name = "Phenotype", definition = "Phenotype") %>%
    set_variable_metadata('relative_growth_rate', display_order=3, display_name = "Relative Growth Rate", definition = "Relative Growth Rate") %>%
    set_variable_metadata('confidence', display_order=4, display_name = "Confidence", definition = "Confidence") %>%
    set_variable_metadata('expected_variance', display_order=5, display_name = "Expected Variance", definition = "Expected Variance") %>%
    set_variable_metadata('RGR_CI_low', display_order=6, display_name = "RGR CI Low", definition = "RGR CI Low") %>%
    set_variable_metadata('RGR_CI_high', display_order=7, display_name = "RGR CI High", definition = "RGR CI High") %>%
    set_variable_metadata('times_analyzed', display_order=8, display_name = "Times Analyzed", definition = "Times Analyzed") %>%
    set_variable_metadata('construct', display_order=9, display_name = "Targeting Vector", definition = "Targeting Vector") %>%
    set_variable_metadata('notes', display_order=10, display_name = "Notes", definition = "Notes")
  
  study = study(name="TEMP_STUDY_NAME", genePhenotype)

  return(study)
}

