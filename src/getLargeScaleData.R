rm(list=ls(all=T))
library(tidyverse)
library(dplyr) # For reorganization
library(stringr) # For string manipulation
# ==== Defining paths and working directories ======
github = 'C:/Users/gara009/OneDrive - PNNL/Documents/GitHub/Geospatial_variables/'
# ====== Read in data ======
comids =  read.csv("C:/Users/gara009/OneDrive - PNNL/Documents/GitHub/Geospatial_variables/Example_Code/v4_RCSFA_Geospatial_Data_Package/v4_RCSFA_Geospatial_Site_Information.csv")%>%
dplyr::select(site = Site_ID, comid = COMID)

sample_data = read_csv(paste0('EC_Data_Package/EC_Field_Metadata.csv'),comment = '#', na = c('N/A', -9999)) %>%
  dplyr::select(Parent_ID,site = Site_ID)

geospatial = read.csv(paste0(github,'v4_RCSFA_Extracted_Geospatial_Data_2025-01-31.csv'))%>%
  mutate(PctFst = pctmxfst2019ws + pctdecid2019ws + pctconif2019ws,
         PctAg = pctcrop2019ws + pcthay2019ws) %>%
  dplyr::select(site,slope,elevation = elevws,AridityWs,Pct_shrub= pctshrb2019ws,precipitation = PrecipSite,PctFst,PctAg) %>%
  distinct(site, .keep_all = TRUE)

# et_npp = read.csv(paste0(github,'Example_Code/v4_RCSFA_Extracted_NPP_ET_Correct_2026-01-02.csv')) %>%
#   group_by(site) %>%
#   summarise(
#     NPP_mean_kgC_m2_yr = mean(NPP_mean_kgC_m2_yr, na.rm = TRUE),
#     ET_mean_kg_m2_yr   = mean(ET_mean_kg_m2_yr, na.rm = TRUE),
#     .groups = "drop"
#   )

et_npp = read.csv(paste0(github,'Example_Code/watershed_pet_et_npp_samplingwindow_mean_2026-01-04.csv')) %>%
  dplyr::select(comid,
    NPP_mean_kgC_m2_yr = NPP_mean,
    PET_mean_mm_yr   = PET_mean,
    ET_mean_mm_yr   = ET_mean
  ) %>%
  left_join(comids, by='comid') %>%
  dplyr::select(-comid)

# ==== Merge data  ===
# data_merged = sample_data %>%
#   left_join(geospatial, by='site') %>%
#   left_join(et_npp %>% dplyr::select(site,NPP_mean_kgC_m2_yr, ET_mean_kg_m2_yr), by='site')
data_merged = sample_data %>%
  left_join(geospatial, by='site') %>%
  left_join(et_npp , by='site')
# ==== Save data ====
write.csv(data_merged, 'data/EC_Field_Geospatial_NPP_PET_VGC.csv', row.names = FALSE)
