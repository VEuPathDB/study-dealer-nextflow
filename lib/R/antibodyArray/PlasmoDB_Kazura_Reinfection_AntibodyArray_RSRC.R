library(tidyverse)
library(study.wrangler)

wrangle <- function() {
  sampleEntityFile = "samples.txt"
  antibodyArrayEntityFile = "profiles.txt"

  my_r_lib <- Sys.getenv("MY_R_LIB")

  core_protein_array_env <- new.env()
  coreProteinMicroarrayScript <- paste0(my_r_lib, "/antibodyArray/coreProteinMicroarray.R");

  source(coreProteinMicroarrayScript, local = core_protein_array_env)

  sample_entity = core_protein_array_env$createProteinArraySampleEntity(sampleEntityFile);

  
  
  # Add the Dataset Variable (here there is only one but adding to be consistent with other ab array studies)
  sample_entity <- sample_entity %>%
    modify_data(mutate(dataset = "Treatment-time to reinfection cohort from Kisumu area, Kenya collected in 2003")) %>%
    sync_variable_metadata() %>%
    set_variable_metadata('dataset', display_name = "Dataset")


  
  sample_entity <- sample_entity %>%
    redetect_column_as_id("SampleName") %>%
    set_variable_metadata('SampleName', entity_name = 'Sample') %>%
    set_variable_metadata('specimen',
                          display_name = "Sample Type") %>%
    set_variable_metadata('parasite.organism',
                          display_name = "Parasite Organism") %>%
    set_variable_metadata('parasite.strain',
                          display_name = "Parasite Strain") %>%
    set_variable_metadata('age_years',
                          display_name = "Age (years)") %>%
    set_variable_metadata('country',
                          display_name = "Country") %>%
    set_variable_metadata('time_to_first_malaria_dx',
                          display_name = "Time to First Malaria Diagnosis") %>%
    set_variable_metadata('time_to_reinfection',
                          display_name = "Time to Reinfection") %>%
    set_variable_metadata('end_of_observation_period',
                          display_name = "End of Observation Period") %>%
    set_variable_metadata('Specimen.collection.date',
                          display_name = "Collection Date") %>%
    set_variable_metadata('microscopy_result',
                          display_name = "Asexual Parasites Present (Microscopy)") %>%
    set_variable_metadata('blood_smear_summary',
                          display_name = "Blood Smear Summary") %>%
    set_variable_metadata('parasite_density_by_microscope_microliter',
                          display_name = "Parasite Density by Microscope (microliter)") %>%
    create_variable_category(
      "participant.information",
      display_name = "Participant Information",
      children = c("age_years", "country", "time_to_reinfection", "time_to_first_malaria_dx", "end_of_observation_period")
    )   %>%
    create_variable_category(
      "parasite.information",
      display_name = "Parasite Information",
      children = c("parasite.organism", "parasite.strain")
    )   %>%
  create_variable_category(
      "laboratory.findings",
      display_name = "Laboratory Findings",
      children = c("microscopy_result", "blood_smear_summary", "parasite_density_by_microscope_microliter")
    )   %>%
    create_variable_category(
      "sample.collection",
      display_name = "Sample Collection",
      children = c("specimen", "Specimen.collection.date")
    )   %>%
    sync_variable_metadata()

  # Inspect the sample entity
  message("\nSample entity summary:")
  inspect(sample_entity)

  array_entity <- core_protein_array_env$createProteinArrayAssayEntity(antibodyArrayEntityFile);


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
