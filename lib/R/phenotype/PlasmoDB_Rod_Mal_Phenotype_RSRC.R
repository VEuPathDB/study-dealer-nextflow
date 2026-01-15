library(tidyverse)
library(study.wrangler)

wrangle <- function() {
  rm(list = ls())

  # read in file
  genePhenotype = entity_from_file("rmgm_modifications.txt")

  # Set meta data for entity
  genePhenotype <- genePhenotype %>% set_entity_metadata(name = "genePhenotypeData", display_name = "Gene Phenotype Data", stable_id="GENE_PHENOTYPE_DATA_ENTITY", display_name_plural="Gene Phenotype Data")

  # default column/variable labels
  genePhenotype <- genePhenotype %>%  set_variable_display_names_from_provider_labels()

  # make gene column a variable
  genePhenotype <- genePhenotype %>%
    redetect_columns_as_variables('gene') %>%
    set_variable_metadata('gene', stable_id = "VEUPATHDB_GENE_ID", display_name = "Gene", provider_label=list("gene"), display_order=1, hidden=list('variableTree'))

  #  deal with the primary Key (ID is rmgmid from the file)
  genePhenotype <- genePhenotype %>%
    modify_data(mutate(ID = row_number())) %>%
    sync_variable_metadata() %>%
    redetect_column_as_id('ID')

  # update columns as needed
  genePhenotype <- genePhenotype %>%
    set_variable_metadata('rmgmid', display_order=2, display_name = "RMGM ID", definition = "Rodent Malaria genetically modified parasite database ID", data_type = "string", data_shape = "categorical", hidden=list('variableTree')) %>%
    redetect_columns_as_variables(c('rmgmid')) %>%
    set_variable_metadata('success_of_genetic_modification', display_order=3, display_name = "Success of Genetic Modification", definition = "Whether the genetic modification was successful (yes/no)") %>%
    set_variable_metadata('pubmedId', display_order=4, display_name = "Pubmed ID", definition = "Pubmed ID", data_type = "string", data_shape = "categorical") %>%
    set_variable_metadata('mod_types', display_order=5, display_name = "Modification Type", definition = "Type of genetic modification (e.g., disrupted, mutated)") %>%
    set_variable_metadata('species', display_order=6, display_name = "Species", definition = "Plasmodium species studied") %>%
    set_variable_metadata('phenotype_found_in', display_order=7, display_name = "Phenotype Found In", definition = "Life stage(s) where phenotype was observed") %>%
    set_variable_metadata('asexualPhenotype', display_order=8, display_name = "Asexual Phenotype", definition = "Phenotype observed in asexual blood stage") %>%
    set_variable_metadata('gametocytePhenotype', display_order=9, display_name = "Gametocyte Phenotype", definition = "Phenotype observed in gametocyte stage") %>%
    set_variable_metadata('ookinetePhenotype', display_order=10, display_name = "Ookinete Phenotype", definition = "Phenotype observed in ookinete stage") %>%
    set_variable_metadata('oocystPhenotype', display_order=11, display_name = "Oocyst Phenotype", definition = "Phenotype observed in oocyst stage") %>%
    set_variable_metadata('sporozoitePhenotype', display_order=12, display_name = "Sporozoite Phenotype", definition = "Phenotype observed in sporozoite stage") %>%
    set_variable_metadata('liverstagePhenotype', display_order=13, display_name = "Liver Stage Phenotype", definition = "Phenotype observed in liver stage") %>%
    set_variables_multivalued('gene' = ';', 'phenotype_found_in' = ';', 'mod_types' = ';', 'species' = ';')

  study = study(name="TEMP_STUDY_NAME", genePhenotype)

  return(study)
}
