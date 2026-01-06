library(tidyverse)
library(study.wrangler)


wrangle <- function() {
  rm(list = ls())

  ## read in file
  genePhenotype = entity_from_file("phenotype_results.txt")

  ## Set meta data for entity
  genePhenotype <- genePhenotype %>%
    set_entity_metadata(name = "genePhenotypeData", display_name = "Gene Phenotype Data", stable_id="genePhenotypeData", display_name_plural="Gene Phenotype Data")

  ## make gene column a variable
  genePhenotype <- genePhenotype %>%
    redetect_columns_as_variables('gene')
  
  ##default column/variable labels
  genePhenotype <- genePhenotype %>%
    set_variable_display_names_from_provider_labels()

  ##  deal with the primary Key (gene variable). boilerplate
  genePhenotype <- genePhenotype %>%
    modify_data(mutate(ID = row_number())) %>%
    sync_variable_metadata() %>%
    redetect_column_as_id('ID')

  genePhenotype <- genePhenotype %>%
    set_variable_metadata('gene', display_name = "Gene", provider_label=list("gene"), display_order=1, hidden=list('variableTree')) %>%
    set_variable_metadata('Feature_Type', display_order=2, display_name = "Feature Type", definition = "Feature Type", hidden=list('variableTree')) %>%
    set_variable_metadata('Gene_Name', display_order=3, display_name = "Gene Name", definition = "Gene Name", hidden=list('variableTree')) %>%
    set_variable_metadata('CGDID', display_order=4, display_name = "CGDID", definition = "CGDID", hidden=list('variableTree')) %>%
    set_variable_metadata('Reference', display_order=5, display_name = "Reference", definition = "Reference") %>%
    set_variable_metadata('Experiment_Type', display_order=6, display_name = "Experiment Type", definition = "Experiment Type") %>%
    set_variable_metadata('Mutant_Type', display_order=7, display_name = "Mutant Type", definition = "Mutant Type") %>%
    set_variable_metadata('Allele', display_order=8, display_name = "Allele", definition = "Allele") %>%
    set_variable_metadata('Strain_background', display_order=9, display_name = "Strain background", definition = "Strain background") %>%
    set_variable_metadata('Phenotype', display_order=10, display_name = "Phenotype", definition = "Phenotype") %>%
    set_variable_metadata('Chemical', display_order=11, display_name = "Chemical", definition = "Chemical") %>%
    set_variable_metadata('Condition', display_order=12, display_name = "Condition", definition = "Condition") %>%
    set_variable_metadata('Details', display_order=13, display_name = "Details", definition = "Details") %>%
    set_variable_metadata('Reporter', display_order=14, display_name = "Reporter", definition = "Reporter") %>%
    set_variable_metadata('Virulence_Model', display_order=16, display_name = "Virulence Model", definition = "Virulence Model")
  
  study = study("calbSC5314_phenotype_CGD_pheno", genePhenotype)

  return(study)
}
