library(tidyverse)
library(study.wrangler)


removeColumns = function(data) {
  return(select(data, -c("PE p-value", "Lung p-value", "Liver p-value", "Spleen p-value")))
}

wrangle <- function() {
  rm(list = ls())

  ## read in file
  genePhenotype = entity_from_file("ToxoDB_MouseScreenTable-Susanne-working.txt", preprocess_fn=removeColumns)

  ## Set meta data for entity
  genePhenotype <- genePhenotype %>% set_entity_metadata(name = "genePhenotypeData", display_name = "Gene Phenotype Data", stable_id="GENE_PHENOTYPE_DATA_ENTITY", display_name_plural="Gene Phenotype Data")

  ## default column/variable labels
  genePhenotype <- genePhenotype %>%  set_variable_display_names_from_provider_labels()


  ## make gene column a variable
  genePhenotype <- genePhenotype %>%
    redetect_columns_as_variables('gene') %>%
    set_variable_metadata('gene', stable_id = "VEUPATHDB_GENE_ID", display_name = "Gene", provider_label=list("gene"), display_order=1, hidden=list('variableTree'))


  ##  deal with the primary Key (gene variable). boilerplate
  genePhenotype <- genePhenotype %>%
    modify_data(mutate(ID = row_number())) %>%
    sync_variable_metadata() %>%
    redetect_column_as_id('ID')


  ## update columns as needed
  genePhenotype <- genePhenotype %>%

  ##modify_data(mutate(across(c(PE.p.value, Liver.p.value, Lung.p.value, Spleen.p.value), ~ log(.x, base = 10)))) %>%
  modify_data(mutate(
      PE.Fitness.rank = rank(PE.Fitness),
      Liver.Fitness.rank = rank(Liver.Fitness),
      Lung.Fitness.rank = rank(Lung.Fitness),
      Spleen.Fitness.rank = rank(Spleen.Fitness),
      PE.Composite.Score.rank = rank(PE.Composite.Score),
      Liver.Composite.Score.rank = rank(Liver.Composite.Score),
      Lung.Composite.Score.rank = rank(Lung.Composite.Score),
      Spleen.Composite.Score.rank = rank(Spleen.Composite.Score)
    )) %>%

  sync_variable_metadata() %>%


  set_variable_metadata('Liver.Fitness.rank', display_name = "Liver Fitness Rank",display_order=9, hidden=list('variableTree')) %>%
  set_variable_metadata('Liver.Composite.Score.rank', display_name = "Liver Composite Score Rank",display_order=10, hidden=list('variableTree')) %>%
  set_variable_metadata('Liver.Fitness', display_name = "Liver Fitness",display_order=11, definition="This is a measure of how much the gRNA in the liver dropped out relative to the start of the screen - the average fold change in gRNA abundance for gRNAs targeting a given gene. Fitness scores quantify the effect of gene knockout on parasite growth, with lower scores indicating a greater fitness defect.") %>%
  set_variable_metadata('Liver.Composite.Score', display_name = "Liver Composite Score",display_order=12, definition="The composite score is the product of multiplying the (Differential Fitness Score)( –log10(p-value of fitness score)).") %>%
  set_variable_metadata('Liver.vs..Cell.Culture.Differential.Fitness', display_name = "Liver vs. Cell Culture Differential Fitness",display_order=13, definition="Differential fitness is a metric that compares the gene fitness score in a tissue vs. cell culture.") %>%


  set_variable_metadata('Lung.Fitness.rank', display_name = "Lung Fitness Rank",display_order=14, hidden=list('variableTree')) %>%
  set_variable_metadata('Lung.Composite.Score.rank', display_name = "Lung Composite Score Rank",display_order=15, hidden=list('variableTree')) %>%
  set_variable_metadata('Lung.Fitness', display_name = "Lung Fitness",display_order=16, definition="This is a measure of how much the gRNA in the lung dropped out relative to the start of the screen - the average fold change in gRNA abundance for gRNAs targeting a given gene. Fitness scores quantify the effect of gene knockout on parasite growth, with lower scores indicating a greater fitness defect.") %>%
  set_variable_metadata('Lung.Composite.Score', display_name = "Lung Composite Score",display_order=17, definition="The composite score is the product of multiplying the (Differential Fitness Score)( –log10(p-value of fitness score)).") %>%
  set_variable_metadata('Lung.vs..Cell.Culture.Differential.Fitness', display_name = "Lung vs. Cell Culture Differential Fitness",display_order=18, definition="Differential fitness is a metric that compares the gene fitness score in a tissue vs. cell culture.") %>%

  set_variable_metadata('PE.Fitness.rank', display_name = "Peritoneum Fitness Rank",display_order=19, hidden=list('variableTree')) %>%
  set_variable_metadata('PE.Composite.Score.rank', display_name = "Peritoneum Composite Score Rank",display_order=20, hidden=list('variableTree')) %>%
  set_variable_metadata('PE.Fitness', display_name = "Peritoneum Fitness",display_order=21, definition="This is a measure of how much the gRNA in the peritoneum dropped out relative to the start of the screen - the average fold change in gRNA abundance for gRNAs targeting a given gene. Fitness scores quantify the effect of gene knockout on parasite growth, with lower scores indicating a greater fitness defect.") %>%
  set_variable_metadata('PE.Composite.Score', display_name = "Peritoneum Composite Score",display_order=22, definition="The composite score is the product of multiplying the (Differential Fitness Score)( –log10(p-value of fitness score)).") %>%
  set_variable_metadata('PE.vs..Cell.Culture.Differential.Fitness', display_name = "Peritoneum vs. Cell Culture Differential Fitness",display_order=23, definition="Differential fitness is a metric that compares the gene fitness score in a tissue vs. cell culture.") %>%


  set_variable_metadata('Spleen.Fitness.rank', display_name = "Spleen Fitness Rank",display_order=24, hidden=list('variableTree')) %>%
  set_variable_metadata('Spleen.Composite.Score.rank', display_name = "Spleen Composite Score Rank",display_order=25, hidden=list('variableTree')) %>%
  set_variable_metadata('Spleen.Fitness', display_name = "Spleen Fitness",display_order=26, definition="This is a measure of how much the gRNA in the spleen dropped out relative to the start of the screen - the average fold change in gRNA abundance for gRNAs targeting a given gene. Fitness scores quantify the effect of gene knockout on parasite growth, with lower scores indicating a greater fitness defect.") %>%
  set_variable_metadata('Spleen.Composite.Score', display_name = "Spleen Composite Score",display_order=27, definition="The composite score is the product of multiplying the (Differential Fitness Score)( –log10(p-value of fitness score)).") %>%
  set_variable_metadata('Spleen.vs..Cell.Culture.Differential.Fitness', display_name = "Spleen vs. Cell Culture Differential Fitness",display_order=28, definition="Differential fitness is a metric that compares the gene fitness score in a tissue vs. cell culture.") %>%

  set_variable_metadata('P4', display_name = "Cell Culture P4", display_order=2, definition="Fitness scores of genes at cell culture passage #4- the average fold change in gRNA abundance for gRNAs targeting a given gene. Fitness scores quantify the effect of gene knockout on parasite growth, with lower scores indicating a greater fitness defect.") %>%
  set_variable_metadata('P5', display_name = "Cell Culture P5", display_order=3, definition="Fitness scores of genes at cell culture passage #5- the average fold change in gRNA abundance for gRNAs targeting a given gene. Fitness scores quantify the effect of gene knockout on parasite growth, with lower scores indicating a greater fitness defect.") %>%
  set_variable_metadata('P6', display_name = "Cell Culture P6", display_order=4, definition="Fitness scores of genes at cell culture passage #6- the average fold change in gRNA abundance for gRNAs targeting a given gene. Fitness scores quantify the effect of gene knockout on parasite growth, with lower scores indicating a greater fitness defect.") %>%
  set_variable_metadata('P7', display_name = "Cell Culture P7", display_order=5, definition="Fitness scores of genes at cell culture passage #7- the average fold change in gRNA abundance for gRNAs targeting a given gene. Fitness scores quantify the effect of gene knockout on parasite growth, with lower scores indicating a greater fitness defect.") %>%
  set_variable_metadata('P8', display_name = "Cell Culture P8", display_order=6, definition="Fitness scores of genes at cell culture passage #8- the average fold change in gRNA abundance for gRNAs targeting a given gene. Fitness scores quantify the effect of gene knockout on parasite growth, with lower scores indicating a greater fitness defect.") %>%
  set_variable_metadata('P9', display_name = "Cell Culture P9", display_order=7, definition="Fitness scores of genes at cell culture passage #9- the average fold change in gRNA abundance for gRNAs targeting a given gene. Fitness scores quantify the effect of gene knockout on parasite growth, with lower scores indicating a greater fitness defect.") %>%
  set_variable_metadata('P10', display_name = "Cell Culture P10", display_order=8, definition="Fitness scores of genes at cell culture passage #10- the average fold change in gRNA abundance for gRNAs targeting a given gene. Fitness scores quantify the effect of gene knockout on parasite growth, with lower scores indicating a greater fitness defect.") %>%

  set_variable_metadata('gene', display_name = "Gene")

  crisprStudy = study(name="TEMP_STUDY_NAME", genePhenotype)

  return(crisprStudy)
}
