
createProteinArraySampleEntity <- function(file) {
  # ===== CREATE SAMPLE ENTITY =====
  message("\n=== Creating Sample Entity ===")

  sample_entity <- entity_from_file(
    file_path = file
  )

  # Set entity metadata
  sample_entity <- sample_entity %>%
    set_entity_metadata(
      name = "Sample",
      display_name = "Sample",
      display_name_plural = "Samples",
      description = "Blood serum samples from subjects",
      stable_id = prefixed_alphanumeric_id(prefix = "ENT_", length = 8, seed_string = "Sample")
    )

  # Set all variables to use provider labels as defaults
  sample_entity <- sample_entity %>%
    set_variable_display_names_from_provider_labels()

  return(sample_entity)
}




createProteinArrayAssayEntity <- function(file, strip_x_prefix = TRUE) {
  wrangler_utils_env <- new.env()
  wranglerUtilsScript <- paste0(my_r_lib, "/utils.R");

  source(wranglerUtilsScript, local = wrangler_utils_env)


  message("\n=== Creating Antibody Microarray Data Entity (Tall Format) ===")
  # Load and preprocess profiles data
  # skip_type_convert since preprocessing already handles type conversion
  message("Loading profiles data...")
  profiles_wide <- entity_from_file(
    file_path = file,
    preprocess_fn = wrangler_utils_env$averageDuplicateGenes,
    skip_type_convert = TRUE
  )

  profiles_data <- profiles_wide@data

  # Convert to tall format: each row is a sample-gene combination
  message("Converting to tall format...")
  profiles_tall <- profiles_data %>%
    pivot_longer(
      cols = -gene,
      names_to = "SampleName",
      values_to = "Score"
    ) %>%
    select(SampleName, gene, Score)

  # Optionally remove X prefix added by R for numeric names
  if (strip_x_prefix) {
    profiles_tall <- profiles_tall %>%
      mutate(SampleName = str_remove(SampleName, "^X(?=\\d)"))
  }

  message("Tall format: ", nrow(profiles_tall), " rows (sample-gene combinations)")
  message("  Samples: ", n_distinct(profiles_tall$SampleName))
  message("  Genes: ", n_distinct(profiles_tall$gene))

  # Create entity from the tall data
  array_entity <- entity_from_tibble(
    profiles_tall,
    name = "AntibodyArray",
    display_name = "Antibody Microarray Assay",
    display_name_plural = "Antibody Microarray Assay",
    description = "P. falciparum protein microarray antibody response data",
    stable_id = "GENE_ANTIBODY_ARRAY_DATA",
    skip_type_convert = TRUE
  )

  # Deal with the primary key (ID variable). Boilerplate
  array_entity <- array_entity %>%
    modify_data(mutate(ID = row_number())) %>%
    sync_variable_metadata() %>%
    redetect_column_as_id('ID')

  # Set variable metadata
  array_entity <- array_entity %>%
    set_variable_metadata('SampleName',
                          display_name = "Sample Name",
                          definition = "Reference to parent Sample entity",
                          display_order = 1,
                          provider_label = list("Sample.Name")) %>%
    set_variable_metadata('gene',
                          display_name = "Gene",
                          stable_id = "VEUPATHDB_GENE_ID",
                          definition = "Gene identifier from the microarray",
                          display_order = 2) %>%
    set_variable_metadata('Score',
                          display_name = "Normalized Intensity",
                          stable_id = "NORMALIZED_INTENSITY",
                          definition = "Normalized antibody array intensity score",
                          display_order = 3,
                          data_shape = "continuous") %>%
    set_parents(names = c('Sample'), id_columns = c('SampleName'))


  return(array_entity)
}



#' Create combined protein array assay entity from multiple profile files
#'
#' Reads each profile file individually, adds a dataset variable based on the
#' file name (e.g., "profiles_India.txt" -> "India"), and combines them into
#' a single entity.
#'
#' @param profileFiles Vector of profile file names
#' @param my_r_lib Path to the R library directory
#' @return Combined array entity
createCombinedProteinArrayAssayEntity <- function(profileFiles, my_r_lib, strip_x_prefix = TRUE) {
  wrangler_utils_env <- new.env()
  wranglerUtilsScript <- paste0(my_r_lib, "/utils.R")
  source(wranglerUtilsScript, local = wrangler_utils_env)

  message("\n=== Creating Combined Antibody Microarray Data Entity (Tall Format) ===")

  # Process each profile file and collect tall format data
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

    # Convert to tall format: each row is a sample-gene combination
    profiles_tall <- profiles_data %>%
      pivot_longer(
        cols = -gene,
        names_to = "SampleName",
        values_to = "Score"
      ) %>%
      mutate(
        # Add dataset identifier
        dataset = dataset_name
      ) %>%
      select(SampleName, gene, Score, dataset)

    # Optionally remove X prefix added by R for numeric names
    if (strip_x_prefix) {
      profiles_tall <- profiles_tall %>%
        mutate(SampleName = str_remove(SampleName, "^X(?=\\d)"))
    }

    message("  - ", nrow(profiles_tall), " rows (sample-gene combinations) with dataset = '", dataset_name, "'")

    all_profiles[[file]] <- profiles_tall
  }

  # Combine all tall format profiles
  combined_profiles <- bind_rows(all_profiles)

  message("\nCombined tall format: ", nrow(combined_profiles), " rows (sample-gene combinations)")
  message("  Samples: ", n_distinct(combined_profiles$SampleName))
  message("  Genes: ", n_distinct(combined_profiles$gene))
  message("  Datasets: ", n_distinct(combined_profiles$dataset))

  # Create entity from the combined tall data
  array_entity <- entity_from_tibble(
    combined_profiles,
    name = "AntibodyArray",
    display_name = "Antibody Microarray Assay",
    display_name_plural = "Antibody Microarray Assay",
    description = "P. falciparum protein microarray antibody response data",
    stable_id = "GENE_ANTIBODY_ARRAY_DATA",
    skip_type_convert = TRUE
  )

  # Deal with the primary key (ID variable). Boilerplate
  array_entity <- array_entity %>%
    modify_data(mutate(ID = row_number())) %>%
    sync_variable_metadata() %>%
    redetect_column_as_id('ID')

  # Set variable metadata
  array_entity <- array_entity %>%
    set_variable_metadata('SampleName',
                          display_name = "Sample Name",
                          definition = "Reference to parent Sample entity",
                          display_order = 1,
                          provider_label = list("Sample.Name")) %>%
    set_variable_metadata('gene',
                          stable_id = "VEUPATHDB_GENE_ID",
                          display_name = "Gene",
                          definition = "Gene identifier from the microarray",
                          display_order = 2) %>%
    set_variable_metadata('Score',
                          display_name = "Normalized Intensity",
                          stable_id = "NORMALIZED_INTENSITY",
                          definition = "Normalized antibody array intensity score",
                          display_order = 3,
                          data_shape = "continuous") %>%
    set_variable_metadata('dataset',
                          display_name = "Dataset",
                          definition = "Source profile dataset",
                          display_order = 4) %>%
    set_parents(names = c('Sample'), id_columns = c('SampleName'))


  return(array_entity)
}
