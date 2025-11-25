wrangle <- function() {
  rm(list = ls())
  
  # read in file
  genePhenotype = entity_from_file("concatenated_phenotype.txt")
  
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
    set_variable_metadata('name', display_order=2, display_name = "Name", definition = "Name") %>%
    set_variable_metadata('pubmedId', display_order=3, display_name = "Pubmed ID", definition = "Pubmed ID") %>%
    set_variable_metadata('modType', display_order=4, display_name = "Mod Type", definition = "Mod Type") %>%
    set_variable_metadata('organism', display_order=5, display_name = "Organism", definition = "Organism") %>%
    set_variable_metadata('qualityTerm', display_order=6, display_name = "Quality Term", definition = "Quality Term") %>%
    set_variable_metadata('entityTerm', display_order=7, display_name = "Entity Term", definition = "Entity Term") %>%
    set_variable_metadata('note', display_order=8, display_name = "Note", definition = "Note") %>%
    set_variable_metadata('experimentType', display_order=9, display_name = "Experiment Type", definition = "Experiment Type") %>%
    set_variable_metadata('allele', display_order=10, display_name = "Allele", definition = "Allele") %>%
    set_variable_metadata('chebiAnnotationExtension', display_order=11, display_name = "chebi Annotation Extension", definition = "chebi Annotation Extension") 

  crisprStudy = study("afumAf293_phenotype_VEuPathDB_curated_phenotype_RSRC", genePhenotype)
  
  return(crisprStudy)
  
}
