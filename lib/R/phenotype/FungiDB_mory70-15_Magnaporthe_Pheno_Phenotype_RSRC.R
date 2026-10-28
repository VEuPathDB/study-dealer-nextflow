library(tidyverse)
library(study.wrangler)


wrangle <- function() {
  rm(list = ls())
  
  # read in file
  genePhenotype = entity_from_file("phenotype_with_names.txt")
  
  # Set meta data for entity
  genePhenotype <- genePhenotype %>% set_entity_metadata(name = "genePhenotypeData", display_name = "Gene Phenotype Data", stable_id="GENE_PHENOTYPE_DATA_ENTITY", display_name_plural="Gene Phenotype Data")
  
  # make gene column a variable; PubMed IDs are categorical strings, not numbers
  genePhenotype <- genePhenotype %>%
    modify_data(mutate(pubmedId = as.character(pubmedId))) %>%
    # MGG_04163's note is 1578 chars (a full abstract); EDA backend limit is 1000,
    # so substitute a precis of it
    modify_data(mutate(note = if_else(
      gene == "MGG_04163" & str_detect(note, "^Deletion of DES1"),
      "Deletion of DES1 had no apparent effect on developmental morphogenesis, appressorium formation or penetration, and the mutant showed no cell wall integrity defects. The Deltades1 mutant was hypersensitive to oxidative stress, with strongly reduced extracellular peroxidase and laccase activity. Rice cells inoculated with the mutant showed strong defense responses (ROS accumulation, autofluorescence, PR gene induction) and infectious hyphal extension was significantly reduced, decreasing pathogenicity. Suppressing ROS with the NADPH oxidase inhibitor DPI restored normal infectious growth, suggesting DES1 is a pathogenicity gene that compromises ROS-mediated plant defense. Chi et al., 2009 [Precis from 1,578 character note; Claude Sonnet 5.5]",
      note))) %>%
    redetect_columns_as_variables(c('gene', 'pubmedId'))
  
  #default column/variable labels
  genePhenotype <- genePhenotype %>%  set_variable_display_names_from_provider_labels()
  
  #  deal with the primary Key (gene variable). boilerplate
  genePhenotype <- genePhenotype %>%
    modify_data(mutate(ID = row_number())) %>%
    sync_variable_metadata() %>%
    redetect_column_as_id('ID')
  
  genePhenotype <- genePhenotype %>%
    set_variable_metadata('gene', stable_id = "VEUPATHDB_GENE_ID", display_name = "Gene", provider_label=list("gene"), display_order=1, hidden=list('variableTree')) %>%
    set_variable_metadata('pubmedId', display_order=2, display_name = "Pubmed ID", definition = "Pubmed ID") %>%
    set_variable_metadata('organism', display_order=3, display_name = "Organism", definition = "Organism") %>%
    set_variable_metadata('qualityTerm', display_order=5, display_name = "Quality Term ID", definition = "Quality Term ID", hidden=list('variableTree')) %>%
    set_variable_metadata('entityTerm', display_order=6, display_name = "Entity Term ID", definition = "Entity Term ID", hidden=list('variableTree')) %>%
    set_variable_metadata('qualityTerm_name', display_order=5, display_name = "Quality Term", definition = "Quality Term") %>%
    set_variable_metadata('entityTerm_name', display_order=6, display_name = "Entity Term", definition = "Entity Term") %>%
    set_variable_metadata('note', display_order=6, display_name = "Note", definition = "Note")

  # Sort PubMed IDs numerically (default is lexical). Any non-integer value
  # goes at the end in lexical order, so unexpected levels fall back quietly.
  sortPmids <- function(pmids) {
    pmids[order(suppressWarnings(as.integer(pmids)), pmids)]
  }

  genePhenotype <- genePhenotype %>%
    set_variable_vocabulary_order("pubmedId", order = sortPmids)
  
  study = study(name="TEMP_STUDY_NAME", genePhenotype)

  return(study)
  
}
