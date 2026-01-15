library(tidyverse)
library(study.wrangler)

wrangle <- function() {
  sampleEntityFile = "samples.txt"

  # Three profile files for this dataset
  profileFiles <- c(
    "profiles_bioinf.txt",
    "profiles_random.txt",
    "profiles_vaccine.txt"
  )

  my_r_lib <- Sys.getenv("MY_R_LIB")

  core_protein_array_env <- new.env()
  coreProteinMicroarrayScript <- paste0(my_r_lib, "/antibodyArray/coreProteinMicroarray.R");

  source(coreProteinMicroarrayScript, local = core_protein_array_env)

  sample_entity = core_protein_array_env$createProteinArraySampleEntity(sampleEntityFile);

  sample_entity <- sample_entity %>%
    redetect_column_as_id("SampleName") %>%
    set_variable_metadata('SampleName', entity_name = 'Sample') %>%
    set_variable_metadata('sex', display_name = "Sex") %>%
    set_variable_metadata('age.years', display_name = "Age (years)") %>%
    set_variable_metadata('PCR.result', display_name = "PCR result") %>%
    set_variable_metadata('fever', display_name = "Febrile") %>%
    set_variable_metadata('parasite.organism', display_name = "Parasite organism") %>%
    set_variable_metadata('parasite_density_by_microscope_microliter', display_name = "Parasite Density by Microscope (microliter)") %>%
    create_variable_category(
      "participant.information",
      display_name = "Participant Information",
      children = c("age.years", "sex", "fever")
    )   %>%
    create_variable_category(
      "parasite.information",
      display_name = "Parasite Information",
      children = c("parasite.organism")
    )   %>%
  create_variable_category(
      "laboratory.findings",
      display_name = "Laboratory Findings",
      children = c("parasite_density_by_microscope_microliter", "PCR.result")
    )   %>%
  sync_variable_metadata()

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
