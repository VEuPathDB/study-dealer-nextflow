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
  array_entity <- createCombinedProteinArrayAssayEntity(profileFiles, my_r_lib)

  # ===== VALIDATE ENTITIES =====
  message("\n=== Validating Entities ===")
  sample_validation <- validate(sample_entity, profile = "eda")
  if (!is.null(sample_validation) && length(sample_validation) > 0) {
    warning("Sample entity validation issues:")
    print(sample_validation)
  } else {
    message("Sample entity validation passed!")
  }

  array_validation <- validate(array_entity, profile = "eda")
  if (!is.null(array_validation) && length(array_validation) > 0) {
    warning("Antibody Array entity validation issues:")
    print(array_validation)
  } else {
    message("Antibody Array entity validation passed!")
  }

  # ===== CREATE STUDY =====
  message("\n=== Creating Study ===")
  study <- study_from_entities(
    entities = list(sample_entity, array_entity),
    name = "TEMP_STUDY_NAME"
  )
  # Validate study
  message("\n=== Validating Study ===")
  study_validation <- validate(study, profile = "eda")
  if (!is.null(study_validation) && length(study_validation) > 0) {
    warning("Study validation issues:")
    print(study_validation)
  } else {
    message("Study validation passed!")
  }

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


#' Create combined protein array assay entity from multiple profile files
#'
#' Reads each profile file individually, adds a dataset variable based on the
#' file name, and combines them into a single entity.
#'
#' @param profileFiles Vector of profile file names
#' @param my_r_lib Path to the R library directory
#' @return Combined array entity
createCombinedProteinArrayAssayEntity <- function(profileFiles, my_r_lib) {
  wrangler_utils_env <- new.env()
  wranglerUtilsScript <- paste0(my_r_lib, "/utils.R")
  source(wranglerUtilsScript, local = wrangler_utils_env)

  message("\n=== Creating Combined Antibody Microarray Data Entity ===")

  # Process each profile file and collect transposed data
  all_profiles <- list()

  for (file in profileFiles) {
    message("Processing: ", file)

    # Extract dataset name from filename (e.g., "profiles_India.txt" -> "India")
    dataset_name <- gsub("^profiles_|\\.txt$", "", file)

    # Load and preprocess profiles data
    profiles_wide <- entity_from_file(
      file_path = file,
      preprocess_fn = wrangler_utils_env$averageDuplicateGenes,
      skip_type_convert = TRUE
    )

    profiles_data <- profiles_wide@data

    # Transpose: convert from genes x samples to samples x genes
    profiles_transposed <- profiles_data %>%
      column_to_rownames("gene") %>%
      t() %>%
      as_tibble(rownames = "SampleName") %>%
      mutate(dataset = dataset_name)

    message("  - ", nrow(profiles_transposed), " samples with dataset = '", dataset_name, "'")

    all_profiles[[file]] <- profiles_transposed
  }

  # Combine all transposed profiles
  combined_profiles <- bind_rows(all_profiles)

  message("\nCombined format: ", nrow(combined_profiles), " samples with ",
          ncol(combined_profiles) - 2, " gene variables (plus SampleName and dataset)")

  # Create entity from the combined data
  array_entity <- entity_from_tibble(
    combined_profiles,
    name = "AntibodyArray",
    display_name = "Antibody Microarray Assay",
    display_name_plural = "Antibody Microarray Assay",
    description = "P. falciparum protein microarray antibody response data",
    stable_id = prefixed_alphanumeric_id(prefix = "ENT_", length = 8, seed_string = "AntibodyArray"),
    skip_type_convert = TRUE
  )

  # Set variable metadata
  array_entity <- array_entity %>%
    set_variable_metadata('SampleName',
                          display_name = "Sample Name",
                          definition = "Reference to parent Sample entity",
                          display_order = 1,
                          provider_label = list("Sample.Name")) %>%
    set_variable_metadata('dataset',
                          display_name = "Dataset",
                          definition = "Source profile dataset") %>%
    modify_data(mutate(ID = row_number())) %>%
    sync_variable_metadata() %>%
    redetect_column_as_id('ID') %>%
    set_variable_display_names_from_provider_labels() %>%
    set_parents("Sample", "SampleName")

  # Get gene variables (exclude SampleName, dataset, and ID)
  gene_vars <- array_entity %>%
    get_variable_metadata() %>%
    filter(!variable %in% c("SampleName", "dataset", "ID")) %>%
    pull(variable)

  array_entity <- array_entity %>%
    create_variable_category(
      category_name = 'normalized_intensity',
      children = gene_vars,
      display_name = 'Gene Antibody Array Intensities',
      definition = 'Gene Antibody Array Intensities'
    ) %>%
    create_variable_collection(
      'normalized_intensity',
      member = 'gene',
      member_plural = 'genes',
      is_proportion = FALSE,
      is_compositional = TRUE
    )

  return(array_entity)
}
