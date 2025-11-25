wrangle <- function() {
  rm(list = ls())

  # read in file
  genePhenotype = entity_from_file("concatenated_phenotype_results.txt")

  # Set meta data for entity
  genePhenotype <- genePhenotype %>% 
    set_entity_metadata(
      name = "genePhenotypeData", 
      display_name = "Gene Phenotype Data", 
      stable_id = "genePhenotypeData", 
      display_name_plural = "Gene Phenotype Data"
    )

  # default column/variable labels
  genePhenotype <- genePhenotype %>% 
    set_variable_display_names_from_provider_labels()

  # deal with the primary key (gene variable). boilerplate
  genePhenotype <- genePhenotype %>%
    modify_data(mutate(ID = row_number())) %>%
    sync_variable_metadata() %>%
    redetect_column_as_id('ID')

  # Set variable metadata with display_order, display_name, and definition
  genePhenotype <- genePhenotype %>%
    set_variable_metadata('PAN_NAME', display_order = 1, display_name = "PAN_NAME", definition = "PAN_NAME") %>%
    set_variable_metadata('GENE_SOURCE_ID', display_order = 2, display_name = "Gene", definition = "GENE_SOURCE_ID") %>%
    set_variable_metadata('PHI.base.entry', display_order = 3, display_name = "PHI-base Entry", definition = "PHI.base.entry") %>%
    set_variable_metadata('Essential.gene', display_order = 4, display_name = "Essential Gene", definition = "Essential.gene") %>%
    set_variable_metadata('Multiple.mutations', display_order = 5, display_name = "Multiple Mutations", definition = "Multiple.mutations") %>%
    set_variable_metadata('Pathogen.species', display_order = 6, display_name = "Pathogen Species", definition = "Pathogen.species") %>%
    set_variable_metadata('Pathogen.strain', display_order = 7, display_name = "Pathogen Strain", definition = "Pathogen.strain") %>%
    set_variable_metadata('Host.species', display_order = 8, display_name = "Host Species", definition = "Host.species") %>%
    set_variable_metadata('Host.strain', display_order = 9, display_name = "Host Strain", definition = "Host.strain") %>%
    set_variable_metadata('Tissue', display_order = 10, display_name = "Tissue", definition = "Tissue") %>%
    set_variable_metadata('Mutant.Phenotype', display_order = 11, display_name = "Mutant Phenotype", definition = "Mutant.Phenotype") %>%
    set_variable_metadata('Disease', display_order = 12, display_name = "Disease", definition = "Disease") %>%
    set_variable_metadata('Mating.defect', display_order = 13, display_name = "Mating Defect", definition = "Mating.defect") %>%
    set_variable_metadata('Prepenetration.defect', display_order = 14, display_name = "Prepenetration Defect", definition = "Prepenetration.defect") %>%
    set_variable_metadata('Penetration.defect', display_order = 15, display_name = "Penetration Defect", definition = "Penetration.defect") %>%
    set_variable_metadata('Postpenetration.defect', display_order = 16, display_name = "Postpenetration Defect", definition = "Postpenetration.defect") %>%
    set_variable_metadata('Disease.manifestation', display_order = 17, display_name = "Disease Manifestation", definition = "Disease.manifestation") %>%
    set_variable_metadata('Vegetative.spores', display_order = 18, display_name = "Vegetative Spores", definition = "Vegetative.spores") %>%
    set_variable_metadata('Sexual.spores', display_order = 19, display_name = "Sexual Spores", definition = "Sexual.spores") %>%
    set_variable_metadata('Invitro.growth', display_order = 20, display_name = "Invitro Growth", definition = "Invitro.growth") %>%
    set_variable_metadata('Spore.germination', display_order = 21, display_name = "Spore Germination", definition = "Spore.germination") %>%
    set_variable_metadata('Gene.inducer', display_order = 22, display_name = "Gene Inducer", definition = "Gene.inducer") %>%
    set_variable_metadata('Experimental.technique', display_order = 23, display_name = "Experimental Technique", definition = "Experimental.technique") %>%
    set_variable_metadata('PMID', display_order = 24, display_name = "PMID", definition = "PMID") %>%
    set_variable_metadata('Comments', display_order = 25, display_name = "Comments", definition = "Comments")

  # Create study entity
  crisprStudy = study("fgraPH-1_PHI-base_curated_phenotype_NAFeaturePhenotypeGeneric_RSRC", genePhenotype)

  return(crisprStudy)
}

