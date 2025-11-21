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
    modify_data(mutate(dataset = "Antibody array of recent exposure to Plasmodium infection")) %>%
    sync_variable_metadata() %>%
    set_variable_metadata('dataset', display_name = "Dataset")



  sample_entity <- sample_entity %>%
    redetect_column_as_id("SampleName") %>%
    set_variable_metadata('SampleName', entity_name = 'Sample') %>%
    # Setting variable display names
    set_variable_metadata('PNAS_paper', display_name = "Study group") %>%
    set_variable_metadata('study_type', display_name = "Study cohort") %>%
    set_variable_metadata('malCount_Before', display_name = "Cumulative malaria episode count") %>%

    set_variable_metadata('city', display_name = "City, village, or region") %>%
    set_variable_metadata('country', display_name = "Country") %>%
    set_variable_metadata('sex', display_name = "Sex") %>%
    set_variable_metadata('age.at.collection', display_name = "Age") %>%
    set_variable_metadata('enrlDate', display_name = "Enrollment date") %>%
    set_variable_metadata('wthdrDate', display_name = "Withdrawal date") %>%
    set_variable_metadata('lastDate', display_name = "Last date observed") %>%
    set_variable_metadata('plasmaDate', display_name = "Sample collection date") %>%
    set_variable_metadata('mal30', display_name = "Malaria diagnosed in 30 days before sample collection") %>%
    set_variable_metadata('mal90', display_name = "Malaria diagnosed in 90 days before sample collection") %>%
    set_variable_metadata('mal180', display_name = "Malaria diagnosed in 180 days before sample collection") %>%
    set_variable_metadata('mal365', display_name = "Malaria diagnosed in 365 days before sample collection") %>%
    set_variable_metadata('par30', display_name = "P. falciparum detected in 30 days before sample collection") %>%
    set_variable_metadata('par90', display_name = "P. falciparum detected in 90 days before sample collection") %>%
    set_variable_metadata('par180', display_name = "P. falciparum detected in 180 days before sample collection") %>%
    set_variable_metadata('par365', display_name = "P. falciparum detected in 365 days before sample collection") %>%
    set_variable_metadata('parasite.organism', display_name = "Parasite organism") %>%
    set_variable_metadata('subject', display_name = "Participant ID") %>%
    set_variable_metadata('hhid', display_name = "Household ID") %>%
    set_variable_metadata('nightlyMosq', display_name = "Nightly female Anophelese count") %>%
    set_variable_metadata('cumMosq', display_name = "Cumulative female Anopheles count") %>%
    set_variable_metadata('incMal', display_name = "Malaria incidence in the previous 365 days") %>%
    set_variable_metadata('mal_at_plasma', display_name = "Malaria diagnosed at sample collection") %>%
    set_variable_metadata('par_at_plasma', display_name = "P. falciparum detected at sample collection") %>%
    set_variable_metadata('mal_7dBefore', display_name = "Malaria diagnosed in 7 days before sample collection") %>%
    set_variable_metadata('mal_7dAfter', display_name = "Malaria diagnosed in 7 days after sample collection") %>%
    set_variable_metadata('timesince_mal', display_name = "Days since last malaria diagnosis") %>%
    set_variable_metadata('timesince_par', display_name = "Days since last P. falciparum detection") %>%
    set_variable_metadata('sincemal', display_name = "Days since last malaria diagnosis, adjusted") %>%
    set_variable_metadata('sincepar', display_name = "Days since last P. falciparum detection, adjusted") %>%
    set_variable_metadata('sincepar_log', display_name = "Analysis metric: Days since last P. falciparum detection, adjusted & log transformed") %>%
    set_variable_metadata('incMal_log', display_name = "Analysis metric: Malaria incidence in the previous 365 days, log transformed") %>%
    set_variable_metadata('pardens', display_name = "P. falciparum density at sample collection") %>%
    create_variable_category(
      "participant.information",
      display_name = "Participant Information",
      children = c(
        "age_years",
        "country",
        "time_to_reinfection",
        "time_to_first_malaria_dx",
        "end_of_observation_period",
        "study_type",
        "city",
        "sex",
        "age.at.collection",
        "subject"
      )
    )  %>%
  create_variable_category(
    "data.set",
    display_name = "Data Set",
    children = c(
      "PNAS_paper"
    )
  ) %>%
  create_variable_category(
    "entomology",
    display_name = "Entomology",
    children = c(
      "cumMosq",  # Cumulative female Anopheles count
      "nightlyMosq" # Nightly female Anophelese count
    )
  ) %>%
  create_variable_category(
    "geographic.location",
    display_name = "Geographic Location",
    children = c(
      "city", # City, village, or region
      "country" # Country
    )
  ) %>%
  create_variable_category(
    "organism.under.investigation",
    display_name = "Organism Under Investigation",
    children = c(
      "parasite.organismy" # Parasite organism
    )
  ) %>%
  create_variable_category(
    "sample.collection",
    display_name = "Sample Collection",
    children = c(
      "plasmaDate" # Sample collection date
    )
  ) %>%
  create_variable_category(
    "sample.source",
    display_name = "Sample Source",
    children = c(
      "age.at.collection", # Age
      "sincepar_log", # Analysis metric: Days since last P. falciparum detection, adjusted & log transformed
      "incMal_log", # Analysis metric: Malaria incidence in the previous 365 days, log transformed
      "malCount_Before", # Cumulative malaria episode count
      "timesince_mal", # Days since last malaria diagnosis
      "sincemal", # Days since last P. falciparum detection
      "timesince_par", # Days since last P. falciparum detection, adjusted
      "sincepar", # Days since last P. falciparum detection, adjusted
      "hhid", # Household ID
      "mal_at_plasma", # Malaria diagnosed at sample collection
      "mal_7dAfter", # Malaria diagnosed in 7 days after sample collection
      "mal_7dBefore", # Malaria diagnosed in 7 days before sample collection
      "mal30", # Malaria diagnosed in 30 days before sample collection
      "mal90", # Malaria diagnosed in 90 days before sample collection
      "mal180", # Malaria diagnosed in 180 days before sample collection
      "mal365", # Malaria diagnosed in 365 days before sample collection
      "incMal", # Malaria incidence in the previous 365 days
      "pardens", # P. falciparum density at sample collection
      "par_at_plasma", # P. falciparum detected at sample collection
      "par30", # P. falciparum detected in 30 days before sample collection
      "par90", # P. falciparum detected in 90 days before sample collection
      "par180", # P. falciparum detected in 180 days before sample collection
      "par365" # P. falciparum detected in 365 days before sample collection
    )
  ) %>%
  create_variable_category(
    "participant.study.details",
    display_name = "Participant Study Details",
    children = c(
      "enrlDate", # Enrollment date
      "wthdrDate", # Withdrawal date
      "lastDate" # Last date observed
    )
  ) %>%
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


  return(study)
}
