wrangle <- function() {
  rm(list = ls())

  # read in file
  genePhenotype = entity_from_file("merged_phenotypes_scores.txt")

  # Set meta data for entity
 genePhenotype <- genePhenotype %>% set_entity_metadata(name = "genePhenotypeData", display_name = "Gene Phenotype Data", stable_id="genePhenotypeData", display_name_plural="Gene Phenotype Data")

   # make gene column a variable
 genePhenotype <- genePhenotype %>% redetect_columns_as_variables('gene')
 
  #default column/variable labels
 genePhenotype <- genePhenotype %>%  set_variable_display_names_from_provider_labels()

  #  deal with the primary Key (gene variable). boilerplate
 genePhenotype <- genePhenotype %>%
   modify_data(mutate(ID = row_number())) %>%
   sync_variable_metadata() %>%
   redetect_column_as_id('ID')

genePhenotype <- genePhenotype %>%
  set_variable_metadata('gene', display_order=1, display_name="Gene", definition="Gene") %>%
  set_variable_metadata('OD600_biofilm_no_dox', display_order=2, display_name="1: OD600 of biofilm growth RPMI-MOPS 37C WITHOUT DOX [pmid:31235750]", definition="1: OD600 of biofilm growth RPMI-MOPS 37C WITHOUT DOX [pmid:31235750]") %>%
  set_variable_metadata('OD600_biofilm_with_dox', display_order=3, display_name="2: OD600 of biofilm growth RPMI-MOPS 37C WITH DOX [pmid:31235750]", definition="2: OD600 of biofilm growth RPMI-MOPS 37C WITH DOX [pmid:31235750]") %>%
  set_variable_metadata('OD600_diff_biofilm', display_order=4, display_name="3: OD600 differential of biofilm growth ([RPMI-MOPS + DOX] - [RPMI-MOPS]) [pmid:31235750]", definition="3: OD600 differential of biofilm growth ([RPMI-MOPS + DOX] - [RPMI-MOPS]) [pmid:31235750]") %>%
  set_variable_metadata('freq_opaque_white_switch', display_order=5, display_name="4: Frequency of opaque-white switching in 5% CO2 [pmid:22621896]", definition="4: Frequency of opaque-white switching in 5% CO2 [pmid:22621896]") %>%
  set_variable_metadata('doubling_time_SC_37C', display_order=6, display_name="5: Doubling time (min) grown in SC media at 37C [pmid:20543849]", definition="5: Doubling time (min) grown in SC media at 37C [pmid:20543849]") %>%
  set_variable_metadata('abundance_BALBc', display_order=7, display_name="6: Abundance of mutant recovered infection of BALBc mice [pmid:20543849]", definition="6: Abundance of mutant recovered infection of BALBc mice [pmid:20543849]") %>%
  set_variable_metadata('colonies_gal_antimycin', display_order=8, display_name="7: Number of colonies growth with Galactose and Antimycin A [pmid:27614020]", definition="7: Number of colonies growth with Galactose and Antimycin A [pmid:27614020]") %>%
  set_variable_metadata('colonies_glc_antimycin', display_order=9, display_name="8: Number of colonies growth with Glucose and Antimycin A [pmid:27614020]", definition="8: Number of colonies growth with Glucose and Antimycin A [pmid:27614020]") %>%
  set_variable_metadata('LDH_Caco2_1', display_order=10, display_name="9: LDH release during infecttion of Confluent Caco-2 cells (CO2 [5%]) [pmid:29871918]", definition="9: LDH release during infecttion of Confluent Caco-2 cells (CO2 [5%]) [pmid:29871918]") %>%
  set_variable_metadata('LDH_Caco2_2', display_order=11, display_name="10: LDH release during infection of Confluent Caco-2 cells (CO2 [5%])  [pmid:29871918]", definition="10: LDH release during infection of Confluent Caco-2 cells (CO2 [5%])  [pmid:29871918]") %>%
  set_variable_metadata('LDH_Caco2_3', display_order=12, display_name="11: LDH release during infection of Confluent Caco-2 cells (CO2 [5%])  [pmid:29871918]", definition="11: LDH release during infection of Confluent Caco-2 cells (CO2 [5%])  [pmid:29871918]") %>%
  set_variable_metadata('LDH_Caco2_4', display_order=13, display_name="12: LDH release during infection of Confluent Caco-2 cells (CO2 [5%])  [pmid:29871918]", definition="12: LDH release during infection of Confluent Caco-2 cells (CO2 [5%])  [pmid:29871918]") %>%
  set_variable_metadata('colony_size_spider', display_order=14, display_name="13: Colony size growth in Spider Media at 30C [pmid:28211844]", definition="13: Colony size growth in Spider Media at 30C [pmid:28211844]") %>%
  set_variable_metadata('morphology_spider', display_order=15, display_name="14: Morphology in Spider Media at 30C [pmid:28211844]", definition="14: Morphology in Spider Media at 30C [pmid:28211844]") %>%
  set_variable_metadata('invasion_ring_spider', display_order=16, display_name="15: Invasion ring size in Spider Media at 30C [pmid:28211844]", definition="15: Invasion ring size in Spider Media at 30C [pmid:28211844]") %>%
  set_variable_metadata('growth_YPD_1M_NaCl', display_order=17, display_name="16: Growth on solid YPD media with 1M NaCl [pmid:27385340]", definition="16: Growth on solid YPD media with 1M NaCl [pmid:27385340]") %>%
  set_variable_metadata('growth_YPD_300M_Menadione', display_order=18, display_name="17: Growth on solid YPD media with 300_M Menadione [pmid:27385340]", definition="17: Growth on solid YPD media with 300_M Menadione [pmid:27385340]") %>%
  set_variable_metadata('survival_GlcNAc', display_order=19, display_name="18: Survival in GlcNAc [2%] at 30C [pmid:26350972]", definition="18: Survival in GlcNAc [2%] at 30C [pmid:26350972]") %>%
  set_variable_metadata('fitness_3d_allKO_1', display_order=20, display_name="19: Fitness after 3d of infection of mice along with all other KO mutants  [pmid:30870623]", definition="19: Fitness after 3d of infection of mice along with all other KO mutants  [pmid:30870623]") %>%
  set_variable_metadata('fitness_3d_allKO_2', display_order=21, display_name="20: Fitness after 3d of infection of mice along with all other KO mutants  [pmid:30870623]", definition="20: Fitness after 3d of infection of mice along with all other KO mutants  [pmid:30870623]") %>%
  set_variable_metadata('fitness_10d_allKO_1', display_order=22, display_name="21: Fitness after 10d of infection of mice along with all other KO mutants [pmid:30870623]", definition="21: Fitness after 10d of infection of mice along with all other KO mutants [pmid:30870623]") %>%
  set_variable_metadata('fitness_10d_allKO_2', display_order=23, display_name="22: Fitness after 10d of infection of mice along with all other KO mutants  [pmid:30870623]", definition="22: Fitness after 10d of infection of mice along with all other KO mutants  [pmid:30870623]") %>%
  set_variable_metadata('fitness_10d_allKO_3', display_order=24, display_name="23: Fitness after 10d of infection of mice along with all other KO mutants  [pmid:30870623]", definition="23: Fitness after 10d of infection of mice along with all other KO mutants  [pmid:30870623]") %>%
  set_variable_metadata('fitness_5d_except_efg1_1', display_order=25, display_name="24: Fitness after 5d of infection of mice along with all other KO mutants except efg1 KO  [pmid:30870623]", definition="24: Fitness after 5d of infection of mice along with all other KO mutants except efg1 KO  [pmid:30870623]") %>%
  set_variable_metadata('fitness_5d_except_efg1_2', display_order=26, display_name="25: Fitness after 5d of infection of mice along with all other KO mutants except efg1 KO  [pmid:30870623]", definition="25: Fitness after 5d of infection of mice along with all other KO mutants except efg1 KO  [pmid:30870623]") %>%
  set_variable_metadata('fitness_5d_except_efg1_3', display_order=27, display_name="26: Fitness after 5d of infection of mice along with all other KO mutants except efg1 KO  [pmid:30870623]", definition="26: Fitness after 5d of infection of mice along with all other KO mutants except efg1 KO  [pmid:30870623]") %>%
  set_variable_metadata('fitness_10d_except_efg1_1', display_order=28, display_name="27: Fitness after 10d of infection of mice along with all other KO mutants except efg1 KO [pmid:30870623]", definition="27: Fitness after 10d of infection of mice along with all other KO mutants except efg1 KO [pmid:30870623]") %>%
  set_variable_metadata('fitness_10d_except_efg1_2', display_order=29, display_name="28: Fitness after 10d of infection of mice along with all other KO mutants except efg1 KO [pmid:30870623]", definition="28: Fitness after 10d of infection of mice along with all other KO mutants except efg1 KO [pmid:30870623]") %>%
  set_variable_metadata('fitness_10d_except_efg1_3', display_order=30, display_name="29: Fitness after 10d of infection of mice along with all other KO mutants except efg1 KO [pmid:30870623]", definition="29: Fitness after 10d of infection of mice along with all other KO mutants except efg1 KO [pmid:30870623]") %>%
  set_variable_metadata('fitness_3d_except_efg1_rob1_brg1_1', display_order=31, display_name="30: Fitness after 3d of infection of mice along with all other KO mutants except efg1 rob1 and brg1 KOs [pmid:30870623]", definition="30: Fitness after 3d of infection of mice along with all other KO mutants except efg1 rob1 and brg1 KOs [pmid:30870623]") %>%
  set_variable_metadata('fitness_3d_except_efg1_rob1_brg1_2', display_order=32, display_name="31: Fitness after 3d of infection of mice along with all other KO mutants except efg1 rob1 and brg1 KOs [pmid:30870623]", definition="31: Fitness after 3d of infection of mice along with all other KO mutants except efg1 rob1 and brg1 KOs [pmid:30870623]") %>%
  set_variable_metadata('fitness_10d_except_efg1_rob1_brg1_1', display_order=33, display_name="32: Fitness after 10d of infection of mice along with all other KO mutants except the efg1 rob1 and brg1 KOs [pmid:30870623]", definition="32: Fitness after 10d of infection of mice along with all other KO mutants except the efg1 rob1 and brg1 KOs [pmid:30870623]") %>%
  set_variable_metadata('fitness_10d_except_efg1_rob1_brg1_2', display_order=34, display_name="33: Fitness after 10d of infection of mice along with all other KO mutants except the efg1 rob1 and brg1 KOs [pmid:30870623]", definition="33: Fitness after 10d of infection of mice along with all other KO mutants except the efg1 rob1 and brg1 KOs [pmid:30870623]") %>%
  set_variable_metadata('growth_butyric_acid_ypd_pH55_37C', display_order=35, display_name="34: Growth in Butyric Acid in YPD (pH 5.5) at 37C [pmid:26297702]", definition="34: Growth in Butyric Acid in YPD (pH 5.5) at 37C [pmid:26297702]") %>%
  set_variable_metadata('growth_propionic_acid_ypd_pH55_37C', display_order=36, display_name="35: Growth in Propionic Acid in YPD (pH 5.5) at 37C [pmid:26297702]", definition="35: Growth in Propionic Acid in YPD (pH 5.5) at 37C [pmid:26297702]") %>%
  set_variable_metadata('growth_acetic_acid_ypd_pH55_37C', display_order=37, display_name="36: Growth in Acetic Acid in YPD (pH 5.5) at 37C [pmid:26297702]", definition="36: Growth in Acetic Acid in YPD (pH 5.5) at 37C [pmid:26297702]") %>%
  set_variable_metadata('growth_lactic_acid_ypd_pH55_37C', display_order=38, display_name="37: Growth in Lactic Acid in YPD (pH 5.5) at 37C [pmid:26297702]", definition="37: Growth in Lactic Acid in YPD (pH 5.5) at 37C [pmid:26297702]") %>%
  set_variable_metadata('growth_nacl_ypd_pH55_37C', display_order=39, display_name="38: Growth in NaCl in YPD (pH 5.5) at 37C [pmid:26297702]", definition="38: Growth in NaCl in YPD (pH 5.5) at 37C [pmid:26297702]") %>%
  set_variable_metadata('growth_ypd_pH55_37C', display_order=40, display_name="39: Growth in YPD (pH 5.5) at 37C [pmid:26297702]", definition="39: Growth in YPD (pH 5.5) at 37C [pmid:26297702]") %>%
  set_variable_metadata('diff_growth_butyric_acid_ypd_pH55_37C', display_order=41, display_name="40: Differential growth significance in Butyric Acid in YPD (pH 5.5) at 37C [pmid:26297702]", definition="40: Differential growth significance in Butyric Acid in YPD (pH 5.5) at 37C [pmid:26297702]") %>%
  set_variable_metadata('diff_growth_propionic_acid_ypd_pH55_37C', display_order=42, display_name="41: Differential growth significance in Propionic Acid in YPD (pH 5.5) at 37C [pmid:26297702]", definition="41: Differential growth significance in Propionic Acid in YPD (pH 5.5) at 37C [pmid:26297702]") %>%
  set_variable_metadata('diff_growth_acetic_acid_ypd_pH55_37C', display_order=43, display_name="42: Differential growth significance in Acetic Acid in YPD (pH 5.5) at 37C [pmid:26297702]", definition="42: Differential growth significance in Acetic Acid in YPD (pH 5.5) at 37C [pmid:26297702]") %>%
  set_variable_metadata('diff_growth_lactic_acid_ypd_pH55_37C', display_order=44, display_name="43: Differential growth significance in Lactic Acid in YPD (pH 5.5) at 37C [pmid:26297702]", definition="43: Differential growth significance in Lactic Acid in YPD (pH 5.5) at 37C [pmid:26297702]") %>%
  set_variable_metadata('diff_growth_nacl_ypd_pH55_37C', display_order=45, display_name="44: Differential growth significance in NaCl in YPD (pH 5.5) at 37C [pmid:26297702]", definition="44: Differential growth significance in NaCl in YPD (pH 5.5) at 37C [pmid:26297702]") %>%
  set_variable_metadata('diff_growth_ypd_pH55_37C', display_order=46, display_name="45: Differential growth significance in YPD (pH 5.5) at 37C [pmid:26297702]", definition="45: Differential growth significance in YPD (pH 5.5) at 37C [pmid:26297702]") %>%
  set_variable_metadata('cell_size_YPD_30C_1', display_order=47, display_name="46: Cell size following growth in YPD at 30C [pmid:30921326]", definition="46: Cell size following growth in YPD at 30C [pmid:30921326]") %>%
  set_variable_metadata('cell_size_YPD_30C_2', display_order=48, display_name="47: Cell size following growth in YPD at 30C  [pmid:30921326]", definition="47: Cell size following growth in YPD at 30C  [pmid:30921326]") %>%
  set_variable_metadata('cell_size_YPD_30C_3', display_order=49, display_name="48: Cell size following growth in YPD at 30C [pmid:30921326]", definition="48: Cell size following growth in YPD at 30C [pmid:30921326]") %>%
  set_variable_metadata('colony_size_undecanoic_1', display_order=50, display_name="49: Colony size on YPD at 30C with undecanoic acid  [pmid:31243082]", definition="49: Colony size on YPD at 30C with undecanoic acid  [pmid:31243082]") %>%
  set_variable_metadata('colony_size_var_undecanoic_1', display_order=51, display_name="50: Colony size variation on YPD at 30C with undecanoic acid  [pmid:31243082]", definition="50: Colony size variation on YPD at 30C with undecanoic acid  [pmid:31243082]") %>%
  set_variable_metadata('colony_size_undecanoic_2', display_order=52, display_name="51: Colony size on YPD at 30C with undecanoic acid [pmid:31243082]", definition="51: Colony size on YPD at 30C with undecanoic acid [pmid:31243082]") %>%
  set_variable_metadata('colony_size_var_undecanoic_2', display_order=53, display_name="52: Colony size variation on YPD at 30C with undecanoic acid [pmid:31243082]", definition="52: Colony size variation on YPD at 30C with undecanoic acid [pmid:31243082]") %>%
  set_variable_metadata('colony_fitness_undecanoic_1', display_order=54, display_name="53: Colony fitness on YPD at 30C with undecanoic acid [pmid:31243082]", definition="53: Colony fitness on YPD at 30C with undecanoic acid [pmid:31243082]") %>%
  set_variable_metadata('colony_fitness_var_undecanoic_1', display_order=55, display_name="54: Colony fitness variation on YPD at 30C with undecanoic acid [pmid:31243082]", definition="54: Colony fitness variation on YPD at 30C with undecanoic acid [pmid:31243082]") %>%
  set_variable_metadata('colony_fitness_sig_undecanoic', display_order=56, display_name="55: Colony fitness significance on YPD at 30C with undecanoic acid [pmid:31243082]", definition="55: Colony fitness significance on YPD at 30C with undecanoic acid [pmid:31243082]") %>%
  set_variable_metadata('adherence_no_dox', display_order=57, display_name="56: Adherence to RPMI + Bovine Serum-primed 6-well plates without DOX [pmid:27870871]", definition="56: Adherence to RPMI + Bovine Serum-primed 6-well plates without DOX [pmid:27870871]") %>%
  set_variable_metadata('adherence_with_dox', display_order=58, display_name="57: Adherence to RPMI + Bovine Serum-primed 6-well plates with 0.5ug_mL DOX [pmid:27870871]", definition="57: Adherence to RPMI + Bovine Serum-primed 6-well plates with 0.5ug_mL DOX [pmid:27870871]")


 crisprStudy = study("calbSC5314_phenotype_Shapiro_2021_Calb_RSRC", genePhenotype)

  return(crisprStudy)

}

