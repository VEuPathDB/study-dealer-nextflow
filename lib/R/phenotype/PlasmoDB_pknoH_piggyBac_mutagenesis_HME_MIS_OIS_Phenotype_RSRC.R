wrangle <- function() {
  rm(list = ls())

  # read in file
  genePhenotype = entity_from_file("HMS_OIS.txt")

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
    set_variable_metadata('gene', display_name = "Gene", provider_label=list("gene"), display_order=1, hidden=list('variableTree')) %>%
    set_variable_metadata('Hybrid.model.score', display_order=2, display_name = "Hybrid Model Score", definition = "Hybrid Model Score") %>%
    set_variable_metadata('Occupancy.index.score', display_order=3, display_name = "Occupancy Index Score", definition = "Occupancy Index Score")
  
  crisprStudy = study("pknoH_phenotype_piggyBac_mutagenesis_HME_MIS_OIS_RSRC", genePhenotype)

  return(crisprStudy)

}
