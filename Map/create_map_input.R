# ==============================================================================
#
# Create input file to map
# 
# Status: In progress
#
# ==============================================================================
#
# Author: Brieanne Forbes (brieanne.forbes@pnnl.gov)
# 17 March 2025

# ==============================================================================

# remove all files
rm(list=ls(all=TRUE))

#set working directory to parent github folder
current_path <- rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path))
setwd("./../")

library(tidyverse)

# ============================ find and read files =============================

metadata <- read_csv('./EC_Data_Package/EC_Field_Metadata.csv') %>%
  select(Parent_ID, Sample_Latitude, Sample_Longitude) %>% # select columns with coordinates 
  rename(Latitude = Sample_Latitude, # rename for simplicity
         Longitude = Sample_Longitude)

# get list of sites used for analysis in this paper 
sites_used <- read_csv('./Data/Medians_of Median_molecular_properties_per_site_and_treatment_unique_formulas.csv') %>%
  select(site) %>%
  distinct() %>%
  pull()


input_file <- metadata %>%
  mutate(Type = case_when(Parent_ID %in% c('EC_011','EC_012', 'EC_023','EC_052','EC_053','EC_057') ~ 'Removed from analysis',
                          Parent_ID %in% sites_used ~ 'Used in analysis',
                          TRUE ~ 'Not used in analysis'))

write_csv(input_file, './Map/Map_Input_File.csv')

# This input file will not be included in the data package. Rerun this script in order to 
# reproduce the map. Use the "Repair Data Source" function to select the file created by this script. 
