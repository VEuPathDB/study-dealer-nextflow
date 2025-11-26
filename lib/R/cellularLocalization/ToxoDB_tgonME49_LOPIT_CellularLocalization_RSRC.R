library(tidyverse)
library(study.wrangler)

wrangle <- function() {
  rm(list = ls())

  # Read in the two LOPIT files
  map_file = "TAGM_MAP_Joint_Probability.tab"
  mcmc_file = "TAGM_MCMC_Joint_Probability.tab"

  # Read both files - first column contains gene IDs
  map_data <- read_tsv(map_file, col_types = cols(.default = "c"))
  mcmc_data <- read_tsv(mcmc_file, col_types = cols(.default = "c"))

  # Rename the first column to "gene" if it's not already named
  if (names(map_data)[1] != "gene") {
    names(map_data)[1] <- "gene"
  }
  if (names(mcmc_data)[1] != "gene") {
    names(mcmc_data)[1] <- "gene"
  }

  # Merge the two files by gene, adding suffixes to distinguish MAP vs MCMC values
  merged_data <- map_data %>%
    left_join(mcmc_data, by = "gene", suffix = c(".MAP", ".MCMC"))

  # Add Dataset variable
  merged_data <- merged_data %>%
    mutate(dataset = "LOPIT Cellular Localization")

  # Create entity from the merged data
  lopitEntity <- entity_from_tibble(
    merged_data,
    name = "lopitData",
    display_name = "LOPIT Data",
    display_name_plural = "LOPIT Data",
    stable_id = "lopitData"
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
                         display_order = 2)

  # Deal with the primary Key (ID variable). Boilerplate
  lopitEntity <- lopitEntity %>%
    modify_data(mutate(ID = row_number())) %>%
    sync_variable_metadata() %>%
    redetect_column_as_id('ID')

  # Get all variable names for categorization
  var_metadata <- lopitEntity %>% get_variable_metadata()

  # Get MAP and MCMC variables (excluding gene, dataset, and ID)
  map_vars <- var_metadata %>%
    filter(str_ends(variable, "\\.MAP")) %>%
    pull(variable)

  mcmc_vars <- var_metadata %>%
    filter(str_ends(variable, "\\.MCMC")) %>%
    pull(variable)

  # Create variable categories for MAP and MCMC methods
  lopitEntity <- lopitEntity %>%
    create_variable_category(
      category_name = "map_probabilities",
      display_name = "TAGM MAP Joint Probabilities",
      definition = "Joint probability estimates for cellular localization using the TAGM MAP (Maximum A Posteriori) method",
      children = map_vars
    ) %>%
    create_variable_category(
      category_name = "mcmc_probabilities",
      display_name = "TAGM MCMC Joint Probabilities",
      definition = "Joint probability estimates for cellular localization using the TAGM MCMC (Markov Chain Monte Carlo) method",
      children = mcmc_vars
    )

  # Create the study
  lopitStudy <- study("tgonME49_cellularLocalization_LOPIT_CellularLocalization_RSRC", lopitEntity)

  return(lopitStudy)
}
