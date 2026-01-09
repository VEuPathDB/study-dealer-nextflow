
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
    mutate(
      # Remove X prefix added by R for numeric names
      SampleName = str_remove(SampleName, "^X(?=\\d)")
    ) %>%
    select(SampleName, gene, Score)

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
    stable_id = prefixed_alphanumeric_id(prefix = "ENT_", length = 8, seed_string = "AntibodyArray"),
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
                          definition = "Gene identifier from the microarray",
                          display_order = 2) %>%
    set_variable_metadata('Score',
                          display_name = "Normalized Intensity",
                          definition = "Normalized antibody array intensity score",
                          display_order = 3,
                          data_shape = "continuous") %>%
    set_parents(names = c('Sample'), id_columns = c('SampleName'))

  # Create variable collection for the Score variable
  array_entity <- array_entity %>%
    create_variable_category(
      category_name = 'normalized_intensity',
      children = 'Score',
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
