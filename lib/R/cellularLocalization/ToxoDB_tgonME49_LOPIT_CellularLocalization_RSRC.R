library(tidyverse)
library(study.wrangler)

wrangle <- function() {
  rm(list = ls())

  # Read in the two LOPIT files
  map_file = "TAGM_MAP_Joint_Probability.tab"
  mcmc_file = "TAGM_MCMC_Joint_Probability.tab"

  # Read both files - first column contains gene IDs
  # Using read.delim for more reliable parsing of tab-delimited files
  map_data <- read.delim(map_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  mcmc_data <- read.delim(mcmc_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

  # Rename the first column to "gene" if it's not already named
  if (names(map_data)[1] != "gene") {
    names(map_data)[1] <- "gene"
  }
  if (names(mcmc_data)[1] != "gene") {
    names(mcmc_data)[1] <- "gene"
  }

  # Convert to tibble for tidyverse operations
  map_data <- as_tibble(map_data)
  mcmc_data <- as_tibble(mcmc_data)

  # Convert MAP data to tall format: find max probability compartment for each gene
  map_tall <- map_data %>%
    pivot_longer(cols = -gene, names_to = "compartment", values_to = "probability") %>%
    mutate(probability = as.numeric(probability)) %>%
    group_by(gene) %>%
    slice_max(probability, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(gene,
           tagm.map.localisation.prediction = compartment,
           tagm.map.localisation.probability = probability)

  # Convert MCMC data to tall format: find max probability compartment for each gene
  mcmc_tall <- mcmc_data %>%
    pivot_longer(cols = -gene, names_to = "compartment", values_to = "probability") %>%
    mutate(probability = as.numeric(probability)) %>%
    group_by(gene) %>%
    slice_max(probability, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(gene,
           tagm.mcmc.localisation.prediction = compartment,
           tagm.mcmc.localisation.probability = probability)

  # Add suffixes to distinguish MAP vs MCMC values before merging
  map_tall <- map_tall %>%
    rename_with(~paste0(., ".MAP"), -gene)

  mcmc_tall <- mcmc_tall %>%
    rename_with(~paste0(., ".MCMC"), -gene)

  # Merge the two datasets by gene
  merged_data <- map_tall %>%
    full_join(mcmc_tall, by = "gene")

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

  # Set dataset variable metadata and individual variable metadata
  lopitEntity <- lopitEntity %>%
    set_variable_metadata('dataset',
                         display_name = "Dataset",
                         display_order = 2) %>%
    set_variable_metadata('tagm.map.localisation.prediction.MAP', display_name = "MAP Prediction") %>%
    set_variable_metadata('tagm.map.localisation.probability.MAP', display_name = "MAP Probability") %>%
    set_variable_metadata('tagm.mcmc.localisation.prediction.MCMC', display_name = "MCMC Prediction") %>%
    set_variable_metadata('tagm.mcmc.localisation.probability.MCMC', display_name = "MCMC Probability")

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
      category_name = "map_localization",
      display_name = "TAGM MAP Localization",
      definition = "Cellular localization predictions and probabilities using the TAGM MAP (Maximum A Posteriori) method",
      children = map_vars
    ) %>%
    create_variable_category(
      category_name = "mcmc_localization",
      display_name = "TAGM MCMC Localization",
      definition = "Cellular localization predictions and probabilities using the TAGM MCMC (Markov Chain Monte Carlo) method",
      children = mcmc_vars
    )

  # Create the study
  lopitStudy <- study("tgonME49_cellularLocalization_LOPIT_CellularLocalization_RSRC", lopitEntity)

  return(lopitStudy)
}
