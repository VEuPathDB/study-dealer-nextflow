library(tidyverse)
library(study.wrangler)

wrangle <- function() {
  rm(list = ls())

  # read in both files
  mutants_ddb_g <- read.table("all-mutants-ddb_g.txt", header = TRUE, sep = "\t", quote = "", comment.char = "", stringsAsFactors = FALSE)
  mutants <- read.table("all-mutants.txt", header = TRUE, sep = "\t", quote = "", comment.char = "", stringsAsFactors = FALSE)

  # remove specified columns
  mutants_ddb_g <- mutants_ddb_g %>%
    select(-`Associated.gene.s.`, -Strain_Descriptor, -Phenotypes)

  mutants <- mutants %>%
    select(-`Associated.gene.s.`)

  # rename Systematic.Name to Systematic_Name in mutants file to match mutants_ddb_g
  mutants <- mutants %>%
    rename(Systematic_Name = Systematic.Name)

  # convert to tibbles for dplyr operations
  mutants_ddb_g <- as_tibble(mutants_ddb_g)
  mutants <- as_tibble(mutants)

  # join on Systematic Name (many-to-many relationship expected)
  joined_data <- mutants_ddb_g %>%
    inner_join(mutants, by = "Systematic_Name", suffix = c("_ddb_g", "_mutants"), relationship = "many-to-many")

  # rename Phenotypes_ddb_g to gene
  joined_data <- joined_data %>%
    rename(gene = DDB_G_ID)

  # create entity from joined data
  genePhenotype = entity_from_tibble(joined_data)

  # Set meta data for entity
  genePhenotype <- genePhenotype %>% set_entity_metadata(name = "genePhenotypeData", display_name = "Gene Phenotype Data", stable_id="genePhenotypeData", display_name_plural="Gene Phenotype Data")

  # make gene column a variable (using DDB_G_ID as the gene identifier)
  #genePhenotype <- genePhenotype %>% redetect_columns_as_variables('DDB_G_ID')

  #default column/variable labels
  genePhenotype <- genePhenotype %>%  set_variable_display_names_from_provider_labels()

  #  deal with the primary Key (gene variable). boilerplate
  genePhenotype <- genePhenotype %>%
    modify_data(mutate(ID = row_number())) %>%
    sync_variable_metadata() %>%
    redetect_column_as_id('ID')

  genePhenotype <- genePhenotype %>%
    set_variable_metadata('gene', display_name = "Gene", provider_label=list("DDB_G_ID"), display_order=1, hidden=list('variableTree')) %>%
    set_variable_metadata('Systematic_Name', display_order=2, display_name = "Systematic Name", definition = "Systematic Name") %>%
    set_variable_metadata('Strain.Descriptor', display_order=3, display_name = "Strain Descriptor", definition = "Strain Descriptor") %>%
    set_variables_multivalued('gene' = ' | ', 'Phenotypes' = ' | ')

  study = study("ddisAX4_DictyBase_phenotype_NAFeaturePhenotypeGeneric", genePhenotype)
  
  return(study)
  
}
