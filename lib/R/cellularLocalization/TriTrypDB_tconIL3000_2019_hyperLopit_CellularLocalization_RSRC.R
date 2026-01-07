library(tidyverse)
library(study.wrangler)

wrangle <- function() {
  rm(list = ls())

  # Read in the two hyperLOPIT files for bloodstream and procyclic forms
  bsf_file = "20220329_TcBSF_n3_TriTrypDB_tables.csv"
  pcf_file = "20220329_TcPCF_n3_TriTrypDB_tables.csv"

  # Read both files - first column contains gene IDs
  bsf_data <- read_csv(bsf_file, col_types = cols(.default = "c"))
  pcf_data <- read_csv(pcf_file, col_types = cols(.default = "c"))

  # Add suffixes to distinguish BSF vs PCF values before merging
  # Remove organism column as it's not needed
  bsf_data <- bsf_data %>%
    select(-organism) %>%
    rename_with(~paste0(., ".BSF"), -gene)

  pcf_data <- pcf_data %>%
    select(-organism) %>%
    rename_with(~paste0(., ".PCF"), -gene)

  # Merge the two files by gene
  merged_data <- bsf_data %>%
    full_join(pcf_data, by = "gene")

  # Add Dataset variable
  merged_data <- merged_data %>%
    mutate(dataset = "hyperLOPIT Cellular Localization")

  # Create entity from the merged data
  lopitEntity <- entity_from_tibble(
    merged_data,
    name = "hyperLopitData",
    display_name = "hyperLOPIT Data",
    display_name_plural = "hyperLOPIT Data",
    stable_id = "hyperLopitData"
  )

  # Set default column/variable labels from provider labels
  lopitEntity <- lopitEntity %>%
    set_variable_display_names_from_provider_labels()

  # Make gene column a variable
  lopitEntity <- lopitEntity %>%
    redetect_columns_as_variables('gene') %>%
    set_variable_metadata('gene',
                         display_name = "Gene",
                         provider_label = list("gene"),
                         display_order = 1,
                         hidden = list('variableTree'))

  # Set dataset variable metadata
  lopitEntity <- lopitEntity %>%
    set_variable_metadata('dataset',
                         display_name = "Dataset",
                         display_order = 2) %>%
    set_variable_metadata('tagm.map.localisation.prediction.BSF', display_name = "Bloodstream Form Prediction") %>%
    set_variable_metadata('tagm.map.localisation.probability.BSF', display_name = "Bloodstream Form Probability") %>%
    set_variable_metadata('tagm.map.localisation.prediction.PCF', display_name = "Procyclic Form Prediction") %>%
    set_variable_metadata('tagm.map.localisation.probability.PCF', display_name = "Procyclic Form Probability")

  # Deal with the primary Key (ID variable). Boilerplate
  lopitEntity <- lopitEntity %>%
    modify_data(mutate(ID = row_number())) %>%
    sync_variable_metadata() %>%
    redetect_column_as_id('ID')

  # Get all variable names for categorization
  var_metadata <- lopitEntity %>% get_variable_metadata()

  # Get BSF and PCF variables (excluding gene, dataset, and ID)
  bsf_vars <- var_metadata %>%
    filter(str_ends(variable, "\\.BSF")) %>%
    pull(variable)

  pcf_vars <- var_metadata %>%
    filter(str_ends(variable, "\\.PCF")) %>%
    pull(variable)

  # Create variable categories for BSF and PCF life stages
  lopitEntity <- lopitEntity %>%
    create_variable_category(
      category_name = "bsf_localization",
      display_name = "Bloodstream Form Localization",
      definition = "TAGM-MAP cellular localization predictions and probabilities for T. congolense bloodstream form",
      children = bsf_vars
    ) %>%
    create_variable_category(
      category_name = "pcf_localization",
      display_name = "Procyclic Form Localization",
      definition = "TAGM-MAP cellular localization predictions and probabilities for T. congolense procyclic form",
      children = pcf_vars
    )

  # Create the study
  lopitStudy <- study(name="TEMP_STUDY_NAME", lopitEntity)

  return(lopitStudy)
}
