library(tidyverse)
library(study.wrangler)


removeColumns = function(data) {
  return(select(data, -c(standard_error, gene_fdr, sg_fdr)))
}

wrangle <- function() {
  rm(list = ls())

  ## read in file. remove unwanted columns
  genePhenotype = entity_from_file("crispr_phenotype.txt", preprocess_fn=removeColumns)

  ## Set meta data for entity
  genePhenotype <- genePhenotype %>% set_entity_metadata(name = "genePhenotypeData", display_name = "Gene Phenotype Data", stable_id="genePhenotypeData", display_name_plural="Gene Phenotype Data")

  ## default column/variable labels
  genePhenotype <- genePhenotype %>%  set_variable_display_names_from_provider_labels()

  ## make gene column a variable
  genePhenotype <- genePhenotype %>%
    redetect_columns_as_variables('gene') %>%
    set_variable_metadata('gene', display_name = "Gene", provider_label=list("gene"), display_order=1, hidden=list('variableTree'))

  ##  deal with the primary Key (gene variable). boilerplate
  genePhenotype <- genePhenotype %>%
    modify_data(mutate(ID = row_number())) %>%
    sync_variable_metadata() %>%
    redetect_column_as_id('ID')


  genePhenotype <- genePhenotype %>%
    set_variable_metadata('mean_phenotype', display_order=2, display_name = "Mean Phenotype Score", definition="In CRISPR screens, a \"mean phenotype score\" refers to the average effect of a specific gene perturbation (e.g., knockout or gene silencing) on a measured phenotype, often calculated across multiple guide RNAs targeting that gene. A lower or more negative score indicates a stronger growth defect.") %>%
    set_variable_metadata('rank', display_order=3, hidden=list('variableTree'),  display_name = "Phenotype Rank", definition="Order/Rank or the Mean Phenotype Score")


  crisprStudy = study("tgonGT1_crisprPhenotype_CrisprScreen_RSRC", genePhenotype)

  return(crisprStudy)

}
