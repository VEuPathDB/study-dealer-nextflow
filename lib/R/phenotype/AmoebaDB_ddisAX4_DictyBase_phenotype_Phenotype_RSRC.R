wrangle <- function() {
  rm(list = ls())
  
  # read in file
  genePhenotype = entity_from_file("phenotype_results.txt")
  
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
    set_variable_metadata('PAN_NAME', display_order=2, display_name = "PAN NAME", definition = "PAN NAME") %>%
    set_variable_metadata('PAN_URI', display_order=3, display_name = "PAN URI", definition = "PAN URI") %>%
    set_variable_metadata('Genename', display_order=4, display_name = "Gene Name", definition = "Gene Name") %>%
    set_variable_metadata('Strain', display_order=5, display_name = "Strain", definition = "Strain") %>%
    set_variable_metadata('Mutantphenotype', display_order=6, display_name = "Mutant Phenotype", definition = "QMutant Phenotype")

  crisprStudy = study("ddisAX4_DictyBase_phenotype_NAFeaturePhenotypeGeneric_RSRC", genePhenotype)
  
  return(crisprStudy)
  
}
