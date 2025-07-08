# ==============================================================================
#
# Format output files from CoreMS to comply with ESS-DIVE CSV reporting 
# format (Velliquette et al. 2021) and consistency with other data packages
#
# ==============================================================================
#
# Author: Brieanne Forbes (brieanne.forbes@pnnl.gov)
# 8 July 2025
#
# ==============================================================================

require(pacman)
p_load(tidyverse) #install and/or library necessary packages

# remove anything in the environment
rm(list=ls(all=T))

# ================================= set dir ================================

current_path <- rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path))

# ============================= update column names ============================

 './Processed_clean_CoreMS_XML_Int_4_reps_1p5ppm_cal_Data.csv'%>%
  read_csv() %>%
  rename_with(~ str_remove(.x, "_[^_]*$")) %>% # remove IAT from sample names
  rename(Calibrated_Mass = 1) %>%
  write_csv('./Processed_clean_CoreMS_XML_Int_4_reps_1p5ppm_cal_Data.csv')

'./FTICR_combined_unique_formulas_Data.csv'%>%
  read_csv() %>%
  rename_with(~ str_remove(.x, "_[^_]*$")) %>% # remove IAT from sample names
  rename(Calibrated_Mass = 1) %>%
  write_csv('./FTICR_combined_unique_formulas_Data.csv')


'./Processed_clean_CoreMS_XML_Int_1p5ppm_cal_pts_Mol.csv' %>%
  read_csv() %>%
  select(-`Molecular.Formula`) %>% #remove column as it is the same as MolForm column 
  rename(Calibrated_Mass = 1,
         Is_Isotopologue = `Is.Isotopologue`,
         Heteroatom_Class = `Heteroatom.Class`,
         Calculated_Mass = `Calculated.m.z`,
         Error_ppm = `m.z.Error..ppm.`) %>%
  write_csv('./Processed_clean_CoreMS_XML_Int_1p5ppm_cal_pts_Mol.csv')

