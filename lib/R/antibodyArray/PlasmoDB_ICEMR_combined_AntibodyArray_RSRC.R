library(tidyverse)
library(study.wrangler)

wrangle <- function() {
  sampleEntityFile = "samples.txt"

  # Multiple profile files for this combined ICEMR dataset
  profileFiles <- c(
    "profiles_Amazonia_Brazil.txt",
    "profiles_Amazonia_Peru.txt",
    "profiles_Amazonia.txt",
    "profiles_India.txt",
    "profiles_Malawi.txt",
    "profiles_South_Africa_Zambia_Dec.txt",
    "profiles_South_Africa_Zambia_Jun.txt",
    "profiles_Southeast_Asia.txt",
    "profiles_Southern_Africa.txt",
    "profiles_South_Pacific_PNG.txt",
    "profiles_South_Pacific.txt",
    "profiles_Uganda_East_Africa.txt"
  )

  my_r_lib <- Sys.getenv("MY_R_LIB")

  core_protein_array_env <- new.env()
  coreProteinMicroarrayScript <- paste0(my_r_lib, "/antibodyArray/coreProteinMicroarray.R");

  source(coreProteinMicroarrayScript, local = core_protein_array_env)

  sample_entity = core_protein_array_env$createProteinArraySampleEntity(sampleEntityFile);

  # Remove completely empty columns (all NA)
  sample_entity <- sample_entity %>%
    modify_data({
      empty_cols <- (.) %>%
        select(-SampleName) %>%
        select(where(~all(is.na(.)))) %>%
        names()

      if (length(empty_cols) > 0) {
        message("Removing ", length(empty_cols), " empty columns: ", paste(empty_cols, collapse=", "))
        (.) %>% select(-all_of(empty_cols))
      } else {
        (.)
      }
    })

  # Apply display names and definitions from ontology mapping file
  sample_entity <- applyOntologyMapping(sample_entity, "header_ontology_mapping_deduplicated.txt")

  sample_entity <- sample_entity %>%
    redetect_column_as_id("SampleName") %>%
    set_variable_metadata('SampleName', entity_name = 'Sample') %>%
    set_variable_metadata('Temperature', unit = 'degrees C') %>%
    set_variable_metadata('P_vivax_parasite_density_by_microscope', unit = 'microliter') %>%
    set_variable_metadata('P_falciparum_parasite_density_by_microscope', unit = 'microliter') %>%
    set_variable_metadata('parasite_density_by_pcr', unit = 'microliter') %>%
    set_variable_metadata('P_vivax_gametocyte_density_by_microscope', unit = 'microliter') %>%
    set_variable_metadata('P_falciparum_gametocyte_density_by_microscope', unit = 'microliter') %>%
    set_variable_metadata('parasite_density_by_microscope', unit = 'microliter') %>%
    set_variable_metadata('p_vivax_parasite_density_by_pcr', unit = 'microliter') %>%
    set_variable_metadata('gametocyte_density_by_microscope', unit = 'microliter') %>%
    set_variable_metadata('p_falciparum_parasite_density_by_pcr', unit = 'microliter') %>%
    set_variable_metadata('p_falciparum_parasite_density_by_rt_pcr', unit = 'microliter') %>%
    set_variable_metadata('hemoglobin_level', unit = 'g/dL')

                                        # Inspect the sample entity
  message("\nSample entity summary:")
  inspect(sample_entity)

  # Process each profile file individually and combine
  array_entity <- core_protein_array_env$createCombinedProteinArrayAssayEntity(profileFiles, my_r_lib)

  # ===== CREATE STUDY =====
  message("\n=== Creating Study ===")
  study <- study_from_entities(
    entities = list(sample_entity, array_entity),
    name = "TEMP_STUDY_NAME"
  )
  return(study)
}


#' Apply ontology mapping to set display names and definitions
#'
#' @param entity The entity to apply mappings to
#' @param mapping_file Path to the header_ontology_mapping file
#' @return Entity with updated display names and definitions
applyOntologyMapping <- function(entity, mapping_file) {
  message("\n=== Applying Ontology Mapping ===")

  # Read the mapping file
  mappings <- read_tsv(mapping_file, col_types = cols(.default = "c"))

  # Get current variables in the entity
  var_metadata <- entity %>% get_variable_metadata()

  # Apply mappings for each variable that has a mapping
  for (i in 1:nrow(mappings)) {
    header_name <- mappings$headerName[i]
    display_name <- mappings$DisplayName[i]
    definition <- mappings$Definition[i]

    # Check if this variable exists in the entity
    if (header_name %in% var_metadata$variable) {
      message("  Mapping: ", header_name, " -> ", display_name)

      # Set display name and definition if they are not empty
      if (!is.na(display_name) && display_name != "") {
        entity <- entity %>%
          set_variable_metadata(header_name, display_name = display_name)
      }

      if (!is.na(definition) && definition != "") {
        entity <- entity %>%
          set_variable_metadata(header_name, definition = definition)
      }
    }
  }

  entity <- entity %>% sync_variable_metadata()

  return(entity)
}
