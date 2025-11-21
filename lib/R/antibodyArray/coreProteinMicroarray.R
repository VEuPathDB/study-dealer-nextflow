
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




createProteinArrayAssayEntity <- function(file) {
  wrangler_utils_env <- new.env()
  wranglerUtilsScript <- paste0(my_r_lib, "/utils.R");

  source(wranglerUtilsScript, local = wrangler_utils_env)

  
  message("\n=== Creating Antibody Microarray Data Entity ===")
  # Load and preprocess profiles data
  # skip_type_convert since preprocessing already handles type conversion
  message("Loading profiles data...")
  profiles_wide <- entity_from_file(
    file_path = file,
    preprocess_fn = wrangler_utils_env$averageDuplicateGenes,
    skip_type_convert = TRUE
  )

  # Transpose: convert from genes x samples to samples x genes
  message("Transposing profiles from genes x samples to samples x genes...")
  profiles_data <- profiles_wide@data

  # Transpose the data
  # First set gene column as rownames, transpose, then convert back to tibble
  profiles_transposed <- profiles_data %>%
    column_to_rownames("gene") %>%
    t() %>%
    as_tibble(rownames = "SampleName")

  message("Transposed format: ", nrow(profiles_transposed), " samples with ",
          ncol(profiles_transposed) - 1, " gene variables")

  # Create entity from the transposed data
  array_entity <- entity_from_tibble(
    profiles_transposed,
    name = "AntibodyArray",
    display_name = "Antibody Microarray Assay",
    display_name_plural = "Antibody Microarray Assay",
    description = "P. falciparum protein microarray antibody response data",
    stable_id = prefixed_alphanumeric_id(prefix = "ENT_", length = 8, seed_string = "AntibodyArray"),
    skip_type_convert = TRUE
  )

  # Set variable metadata for SampleName
  array_entity <- array_entity %>%
    set_variable_metadata('SampleName',
                          display_name = "Sample Name",
                          definition = "Reference to parent Sample entity",
                          display_order = 1,
                          provider_label = list("Sample.Name"))  %>%
    modify_data(mutate(ID = row_number())) %>%
    sync_variable_metadata() %>%
    redetect_column_as_id('ID') %>%
    set_variable_display_names_from_provider_labels()  %>%
    set_parents("Sample", "SampleName")

  
    array_entity <- array_entity %>%
      create_variable_category(
        category_name = 'normalized_intensity',
        children = array_entity %>% get_variable_metadata() %>% pull(variable),
        display_name = 'Gene Antibody Array Intensities',
        definition = 'Gene Antibody Array Intensities'
      ) %>%
      create_variable_collection(
        'normalized_intensity',
        member = 'gene',
        member_plural = 'genes',
        is_proportion =  FALSE,
        is_compositional = TRUE
      )

  return(array_entity)
}
