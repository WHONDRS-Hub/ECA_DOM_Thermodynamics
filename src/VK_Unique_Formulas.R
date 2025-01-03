# ==== Loading libraries =========
rm(list=ls(all=T))

library(stringr); library(devtools);  library("plyr")
library("readr");  library(readxl);library(crayon); library(vegan)
# Load in necessary libraries first
library(reshape2)
library(ggpubr) # For to combine plots
library(dplyr) # For reorganization
library(stringr) # For string manipulation
library(tidyr)

# ==== Defining paths and working directories ======
github = 'C:/Users/gara009/OneDrive - PNNL/Documents/GitHub/ECA_DOM_Thermodynamics/'
# ====== Read in data ======
# Processed ICR Data
data = read.csv(list.files(path = github, pattern = "*unique_formulas_Data.csv", full.names = T),row.names = 1)
mol = read.csv(list.files(path = github, pattern = "*cal_pts_Mol.csv"), row.names = 1)
# Fixing colnames 
colnames(data) = gsub('SIR.','SIR-',colnames(data))

sample_data = read_csv(paste0(github,'EC_Data_Package/Sample_Data/EC_Sediment_Sample_Data_Summary.csv'),comment = '#', na = c('N/A', -9999)) %>%
  slice(-(1:11))%>%
  mutate_at(vars(-Sample_Name,-Field_Name,-IGSN,-Material), as.numeric) %>%
  dplyr::select(-Field_Name)

effect_size = read_csv(paste0(github,'EC_Data_Package/Sample_Data/EC_Sediment_Effect_Size.csv'),comment = '#', na = c('N/A', -9999)) %>%
  slice(-(1:11))%>%
  dplyr::select('Sample_Name',"Effect_Size_Respiration_Rate_mg_DO_per_L_per_H","Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H","Effect_Size_Initial_Gravimetric_Moisture_g_per_g","Effect_Size_Final_Gravimetric_Moisture_g_per_g","Effect_Size_Extractable_NPOC_mg_per_kg","Effect_Size_Extractable_TN_mg_per_L")%>%
  mutate_at(vars(-Sample_Name), as.numeric)


effect_size$site = effect_size$Sample_Name
effect_size$site = gsub('_all','',effect_size$site)

factors = data.frame(Sample_Name = sample_data$Sample_Name, site = sample_data$Sample_Name, Treatment = sample_data$Sample_Name)
factors$site = str_extract(factors$site, "EC_0[0-9]{2}|EC_([A-Za-z0-9]+)")
factors$Treatment = str_extract(factors$Treatment, "W|D|Blk")

factors$Treatment = gsub('W','Wet', factors$Treatment)
factors$Treatment = gsub('D','Dry', factors$Treatment)

sample_data = merge(factors,sample_data, by = 'Sample_Name')
# ======= Calculate Relative abundance of VK classes =======
mol2 = mol %>% dplyr::select(bs1_class)

merged_data = merge(data, mol2, by = 0)

# Calculate the relative abundance of Van Krevelen classes per sample
relative_abundance <- merged_data %>%
  gather(key = "sample", value = "intensity", -Row.names, -bs1_class) %>%
  group_by(sample, bs1_class) %>%
  summarise(total_intensity = sum(intensity, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(sample) %>%
  mutate(relative_abundance = total_intensity / sum(total_intensity) * 100) %>%
  ungroup()

relative_abundance$Sample_Name = str_extract(relative_abundance$sample,    "EC_\\d{3}_[^_]*-(W|D)")
relative_abundance$Sample_Name = gsub('_SIR','',relative_abundance$Sample_Name)

relative_abundance = merge(relative_abundance,factors, by = 'Sample_Name', all = T)


# Create the ggplot
plots <- relative_abundance %>%
  filter(bs1_class == "Amino Sugar") %>%
  ggplot(aes(x = as.factor(site), y = relative_abundance, fill = Treatment)) +
  geom_boxplot() +
  facet_wrap(~ bs1_class, scales = "free_y") +
  labs(
    x = " ",
    y = "Percent Relative Abundance",
    color = "Treatment"
  ) +
  theme_bw() + theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Print the plot
print(plots)
