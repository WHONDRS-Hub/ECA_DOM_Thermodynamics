# ==== Loading libraries =========
rm(list=ls(all=T))
library(stringr); library(devtools);  library("plyr")
library("readr");  library(readxl);library(crayon); library(vegan)
# Load in necessary libraries first
library(reshape2)
library(ggpubr) # For to combine plots
library(dplyr) # For reorganization
library(stringr) # For string manipulation
# ==== Defining paths and working directories ======
github = 'C:/Users/gara009/OneDrive - PNNL/Documents/GitHub/ECA_DOM_Thermodynamics/'
data_path = paste0(github,'Data/')
figure_path = paste0(github,'Figures/')
# ====== Read in data ======
# Processed ICR Data
data = read.csv(list.files(path = data_path, pattern = "*Median_Summary_Unique_Peaks.csv", full.names = T))
names(data)[1] = 'Sample_Name'

# ===== Calculate a median of these metrics per site and treatment ====

df.stats <- data %>%
  select(Sample_Name, contains("median"), number.of.peaks, peaks.with.formula) %>%
  select(-AI.median,-DBE_O.median,-KenMass.median,-KenDef.median,-HtoC.median,-"OtoC.median",-"NtoC.median",-"PtoC.median",-"NtoP.median",-"GFE_0.median",-GFEperC_0.median,-lambdaO2_0.median,-"mass.median",-"number.of.peaks")%>%
  mutate(site = str_extract(Sample_Name, "^[^_]+_[^_]+")) %>%
  mutate(Treatment = if_else(str_detect(Sample_Name, "W"), "Wet",if_else(str_detect(Sample_Name, "D"), "Dry", NA_character_)))


df.medians = df.stats %>%
  group_by(site,Treatment) %>%
  summarise(
    Median_AI_mod = median(AI_Mod.median, na.rm = TRUE),
    Median_NOSC = median(NOSC.median, na.rm = TRUE),
    Median_DBE = median(DBE.median, na.rm = TRUE),
    Median_delGcoxPerCmol = median(GFEperC_7.median, na.rm = TRUE),
    Median_delGcoxPerCompmol = median(GFE_7.median, na.rm = TRUE),
    Median_Lambda = median(lambdaO2_7.median, na.rm = TRUE),
    Median_Formulas = median(peaks.with.formula,na.rm = TRUE)
  )

write.csv(df.medians, 'Data/Medians_of Median_molecular_properties_per_site_and_treatment_unique_formulas.csv', row.names = F)
