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
    modify_data(mutate(dataset = "Protein targets of serum antibodies in response to infection")) %>%
    sync_variable_metadata() %>%
    set_variable_metadata('dataset', display_name = "Dataset")

  sample_entity <- sample_entity %>%
    set_variable_metadata('SampleName', entity_name = 'Sample') %>%
    set_variable_metadata('Subject.ID',
                          display_name = "Host ID") %>%
    set_variable_metadata('age.years',
                          display_name = "Age (years)") %>%
    set_variable_metadata('organism',
                          display_name = "Host Organism") %>%
    set_variable_metadata('sex',
                          display_name = "Sex") %>%
    set_variable_metadata('health_status',
                          display_name = "Health Status") %>%
    set_variable_metadata('time_to_onset.days',
                          display_name = "Time to Onset") %>%
    set_variable_metadata('Collection.date',
                          display_name = "Collection Date") %>%
    set_variable_metadata('specimen',
                          display_name = "Sample Type") %>%
    set_variable_metadata('iPCR.result',
                          display_name = "iPCR Result") %>%
    set_variable_metadata('technical.replicate',
                          display_name = "Technical Replicate") %>%
    set_variable_metadata('parasite.organism',
                          display_name = "Parasite Organism") %>%
    set_variable_metadata('parasite.strain',
                          display_name = "Parasite Strain") %>%
    redetect_column_as_id("SampleName") %>%
    create_variable_category(
      "participant.information",
      display_name = "Participant Information",
      children = c("organism", "time_to_onset.days", "sex", "age.years", "Subject.ID", "health_status")
    )   %>%
    create_variable_category(
      "parasite.information",
      display_name = "Parisite Information",
      children = c("parasite.organism", "parasite.strain")
    )   %>%
    create_variable_category(
      "laboratory.findings",
      display_name = "Laboratory Findings",
      children = c("iPCR.result")
    )   %>%
    create_variable_category(
      "sample.collection",
      display_name = "Sample Collection",
      children = c("Collection.date", "specimen", "technical.replicate")
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
    name = "PlasmoDB_Crompton_Mali_AntibodyArray_RSRC"
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

  ## # Print summary
  ## message("\n=== Curation Summary ===")
  ## message("Study name: ", study@name)
  ## message("Number of entities: 2")
  ## message("Sample entity: ", nrow(sample_entity@data), " observations, ",
  ##         ncol(sample_entity@data), " variables")
  ## message("Antibody Array entity: ", nrow(array_entity@data), " observations (gene-sample pairs), ",
  ##         ncol(array_entity@data), " variables")

  return(study)
}
