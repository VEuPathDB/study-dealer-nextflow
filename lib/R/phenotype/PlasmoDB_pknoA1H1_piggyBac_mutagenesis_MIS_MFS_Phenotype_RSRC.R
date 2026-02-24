library(tidyverse)
library(study.wrangler)


wrangle <- function() {
  rm(list = ls())

  # read in file
  genePhenotype = entity_from_file("MIS_MFS_MISplus_mutability.txt")

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

  ## create rank variables and update metadata
  genePhenotype <- genePhenotype %>%
    modify_data(mutate(
      MIS.plus.rank = rank(MIS.plus),
      MIS.rank = rank(MIS),
      MFS.rank = rank(MFS)
    )) %>%
    sync_variable_metadata() %>%

    set_variable_metadata('gene', stable_id = "VEUPATHDB_GENE_ID", display_name = "Gene", provider_label=list("gene"), display_order=1, hidden=list('variableTree')) %>%
    set_variable_metadata('mutability', display_order=2, display_name = "Mutability in CDS", definition = "Mutability in CDS") %>%
    set_variable_metadata('MIS.plus.rank', display_order=3, display_name = "MIS Plus Rank", hidden=list('variableTree'), definition="Rank of MIS Plus") %>%
    set_variable_metadata('MIS.plus', display_order=4, display_name = "MIS Plus", definition = "MIS Plus") %>%
    set_variable_metadata('MIS.rank', display_order=5, display_name = "MIS Rank", hidden=list('variableTree'), definition="Rank of MIS") %>%
    set_variable_metadata('MIS', display_order=6, display_name = "MIS", definition = "MIS") %>%
    set_variable_metadata('MFS.rank', display_order=7, display_name = "MFS Rank", hidden=list('variableTree'), definition="Rank of MFS") %>%
    set_variable_metadata('MFS', display_order=8, display_name = "MFS", definition = "MFS")
  
  study = study(name="TEMP_STUDY_NAME", genePhenotype)

  return(study)

}
