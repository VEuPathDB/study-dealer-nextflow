library(tidyverse)
library(study.wrangler)


wrangle <- function() {
  rm(list = ls())
  
  # read in file
  genePhenotype = entity_from_file("concatenated_phenotype_with_names.txt")
  
  # Set meta data for entity
  genePhenotype <- genePhenotype %>% set_entity_metadata(name = "genePhenotypeData", display_name = "Gene Phenotype Data", stable_id="GENE_PHENOTYPE_DATA_ENTITY", display_name_plural="Gene Phenotype Data")
  
  # make gene column a variable
  genePhenotype <- genePhenotype %>% redetect_columns_as_variables('geneSourceId')
  
  #default column/variable labels
  genePhenotype <- genePhenotype %>%  set_variable_display_names_from_provider_labels()
  
  #  deal with the primary Key (gene variable). boilerplate
  genePhenotype <- genePhenotype %>%
    modify_data(mutate(ID = row_number())) %>%
    sync_variable_metadata() %>%
    redetect_column_as_id('ID')
  
  genePhenotype <- genePhenotype %>%
    set_variable_metadata('geneSourceId', stable_id = "VEUPATHDB_GENE_ID", display_name = "Gene", provider_label=list("geneSourceId"), display_order=1, hidden=list('variableTree')) %>%
    set_variable_metadata('name', display_order=2, display_name = "Name", definition = "Name", hidden=list('variableTree')) %>%
    set_variable_metadata('pubmedId', display_order=2, display_name = "Pubmed ID", definition = "Pubmed ID", data_type = "string", data_shape = "categorical") %>%
    set_variable_metadata('modType', display_order=4, display_name = "Modification Type", definition = "Modification Type") %>%
    set_variable_metadata('organism', display_order=5, display_name = "Organism", definition = "Organism") %>%
    set_variable_metadata('qualityTerm', display_order=5, display_name = "Quality Term ID", definition = "Quality Term ID", hidden=list('variableTree')) %>%
    set_variable_metadata('entityTerm', display_order=6, display_name = "Entity Term ID", definition = "Entity Term ID", hidden=list('variableTree')) %>%
    set_variable_metadata('qualityTerm_name', display_order=5, display_name = "Phenotype Quality", definition = "A PATO term describing how the phenotype deviates (viable, lethal, increased amount, delayed, abnormal). This is the modifier/attribute of the phenotype.") %>%
    set_variable_metadata('entityTerm_name', display_order=6, display_name = "Affected Process/Structure", definition = "A GO/APO/CHEBI/FAO term naming what biological thing is affected (cell growth, conidium formation, virulence, RNA accumulation).") %>%
    set_variable_metadata('note', display_order=8, display_name = "Description", definition = "Description") %>%
    set_variable_metadata('experimentType', display_order=9, display_name = "Experiment Type", definition = "Experiment Type") %>%
    set_variable_metadata('allele', display_order=10, display_name = "Allele", definition = "Allele") %>%
    set_variable_metadata('AnnotationExtension', display_order=11, display_name = "Phenotype Target", definition = "A ChEBI compound, gene ID, or other qualifier that specifies the context of the observation (e.g. the specific chemical being accumulated, or the gene whose
  RNA changes).") %>%
    modify_data(mutate(pubmedId = as.character(pubmedId)))

  study = study(name="TEMP_STUDY_NAME", genePhenotype)

  return(study)
  
}
