library(tidyverse)
library(study.wrangler)

wrangle <- function() {
  rm(list = ls())

  # read in file
  genePhenotype = entity_from_file("MFS_MIS.txt") ##, preprocess_fn=removeColumns)

  # Set meta data for entity
  genePhenotype <- genePhenotype %>% set_entity_metadata(name = "genePhenotypeData", display_name = "Gene Phenotype Data", stable_id="genePhenotypeData", display_name_plural="Gene Phenotype Data")

  #default column/variable labels
  genePhenotype <- genePhenotype %>%  set_variable_display_names_from_provider_labels()

   # make gene column a variable
   genePhenotype <- genePhenotype %>%
    redetect_columns_as_variables('gene') %>%
    set_variable_metadata('gene', display_name = "Gene", provider_label=list("gene"), display_order=1, hidden=list('variableTree'))
 
  #  deal with the primary Key (gene variable). boilerplate
  genePhenotype <- genePhenotype %>%
    modify_data(mutate(ID = row_number())) %>%
    sync_variable_metadata() %>%
    redetect_column_as_id('ID')

  ## update columns as needed
  genePhenotype <- genePhenotype %>%
    set_variable_metadata('mutant.fitness.score', display_order=2, display_name = "Mutant Fitness Score", definition="Mutant fitness score") %>%
    set_variable_metadata('mutagenesis.index.score', display_order=3, display_name = "Mutagenesis Index Score", definition="Mutagenesis index score")

  study = study("pfal3D7_phenotype_pB_mutagenesis_MIS_MFS", genePhenotype)

  return(study)

}
