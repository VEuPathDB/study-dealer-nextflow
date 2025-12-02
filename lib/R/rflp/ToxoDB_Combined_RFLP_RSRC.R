library(tidyverse)
library(study.wrangler)

wrangle <- function() {
  rm(list = ls())

  message("\n=== Processing Combined RFLP Data ===")

  # Get all .txt files from the directory
  data_dir <- "."
  rflp_files <- list.files(path = data_dir, pattern = "\\.txt$", full.names = TRUE)

  # Exclude the HTML reference file if present
  rflp_files <- rflp_files[!grepl("\\.html$", rflp_files)]

  message("Found ", length(rflp_files), " data files to combine")

  # Read and combine all files
  all_isolates <- list()

  for (file in rflp_files) {
    message("Processing: ", basename(file))

    # Read the file
    data <- read_tsv(file, col_types = cols(.default = "c"), show_col_types = FALSE)

    # Standardize the isolate column name
    # Different files use different names: "Isolates", "Isolate ID", "Isolate designation", "Isolate"
    isolate_col <- names(data)[1]
    if (isolate_col %in% c("Isolates", "Isolate ID", "Isolate designation", "Isolate")) {
      data <- data %>% rename(Isolate_ID = !!isolate_col)
    }

    # Add source file column for mapping
    data <- data %>% mutate(source_file_name = basename(file))

    message("  - ", nrow(data), " isolates from ", basename(file))

    all_isolates[[file]] <- data
  }

  # Combine all data - bind_rows will include all columns from all files
  # Missing columns in individual files will be filled with NA
  combined_data <- bind_rows(all_isolates)

  # Map file names to dataset names
  combined_data <- combined_data %>%
    mutate(Dataset = case_when(
      source_file_name == "TgDogBA22.txt" ~ "Canine RFLP genotypes",
      source_file_name == "Menezes_chicken.txt" ~ "Chicken RFLP genotypes",
      source_file_name == "Scherrer_Eurasian_lynx_beaver.txt" ~ "Eurasian lynx and beaver RFLP genotype",
      source_file_name == "ChunleiSu.txt" ~ "RFLP genotypes",
      source_file_name == "Bird.txt" ~ "RFLP genotypes of Brazilian isolates from birds for human consumption",
      source_file_name == "Pig_Goat.txt" ~ "RFLP genotypes of Brazilian isolates from pig and goat for human consumption",
      source_file_name == "Tgon_Cat_Brazil.txt" ~ "Molecular markers in cats RFLP phenotype",
      TRUE ~ source_file_name  # fallback to filename if not matched
    )) %>%
    select(-source_file_name)  # Remove the temporary column

  message("\nCombined total: ", nrow(combined_data), " isolates from ", length(rflp_files), " files")
  message("Total columns: ", ncol(combined_data))

  # Show dataset distribution
  message("\n=== Dataset Distribution ===")
  dataset_counts <- combined_data %>% count(Dataset)
  print(dataset_counts)

  # Create entity from combined data
  rflp_entity <- entity_from_tibble(
    combined_data,
    name = "RFLPIsolate",
    display_name = "RFLP Isolate",
    display_name_plural = "RFLP Isolates",
    description = "Restriction Fragment Length Polymorphism (RFLP) genotyping data for Toxoplasma gondii isolates",
    stable_id = prefixed_alphanumeric_id(prefix = "ENT_", length = 8, seed_string = "RFLPIsolate")
  )

  # Set variable display names from provider labels
  rflp_entity <- rflp_entity %>%
    set_variable_display_names_from_provider_labels()

  # Create unique ID and set as primary key
  rflp_entity <- rflp_entity %>%
    modify_data(mutate(ID = row_number())) %>%
    sync_variable_metadata() %>%
    redetect_column_as_id('ID')

  # Set specific variable metadata for key columns
  rflp_entity <- rflp_entity %>%
    set_variable_metadata('Dataset',
                          display_name = "Dataset",
                          definition = "Dataset source for the RFLP isolate",
                          display_order = 1) %>%
    set_variable_metadata('Isolate_ID',
                          display_name = "Isolate ID",
                          definition = "Unique identifier for the RFLP isolate",
                          display_order = 2) %>%
    set_variable_metadata('Host_Common_Name',
                          display_name = "Host Common Name",
                          definition = "Common name of the host organism",
                          display_order = 3) %>%
    set_variable_metadata('Host_Species',
                          display_name = "Host Species",
                          definition = "Scientific name of the host species",
                          display_order = 4) %>%
    set_variable_metadata('Country',
                          display_name = "Country",
                          definition = "Country where the isolate was collected",
                          display_order = 5) %>%
    set_variable_metadata('Continent',
                          display_name = "Continent",
                          definition = "Continent where the isolate was collected",
                          display_order = 6) %>%
    set_variable_metadata('City_Region',
                          display_name = "City/Region",
                          definition = "City or region where the isolate was collected",
                          display_order = 7) %>%
    set_variable_metadata('ToxoDB_Genotype',
                          display_name = "PCR-RFLP genotype #",
                          definition = "PCR-RFLP genotype #",
                          display_order = 8) %>%
    set_variable_metadata('Organc',
                          display_name = "Organ",
                          display_order = 9) %>%
    create_variable_category(
      "Genetic_Marker",
      display_name = "Genetic Marker",
      children = c("SAG1","X5..3..SAG2","alt..SAG2","SAG3","CS3","BTUB","GRA6","c22.8","c29.2","L358","PK1","Apico")
    )


  # Set PMID metadata if the column exists
  if ("PMID" %in% names(rflp_entity@data)) {
    rflp_entity <- rflp_entity %>%
      set_variable_metadata('PMID',
                            display_name = "PubMed ID",
                            definition = "PubMed identifier for the publication",
                            data_type = "string")
    message("PMID variable set to data_type: string")
  }

  # Validate entity
  message("\n=== Validating Entity ===")
  entity_validation <- validate(rflp_entity, profile = "eda")
  if (!is.null(entity_validation) && length(entity_validation) > 0) {
    warning("Entity validation issues:")
    print(entity_validation)
  } else {
    message("Entity validation passed!")
  }

  # Create study
  message("\n=== Creating Study ===")
  rflp_study <- study(
    name = "ToxoDB_Combined_RFLP_RSRC",
    rflp_entity
  )

  # Validate study
  message("\n=== Validating Study ===")
  study_validation <- validate(rflp_study, profile = "eda")
  if (!is.null(study_validation) && length(study_validation) > 0) {
    warning("Study validation issues:")
    print(study_validation)
  } else {
    message("Study validation passed!")
  }

  message("\n=== Wrangle Complete ===")

  return(rflp_study)
}
