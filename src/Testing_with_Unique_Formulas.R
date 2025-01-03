rm(list=ls(all=T))

# ===== Load necessary libraries ======
library(dplyr)
library(tidyr)
library(stringr) # For string manipulation
library(ggplot2)
library(readr)
# ==== Working directories =====
input_path = 'C:/Users/gara009/OneDrive - PNNL/Documents - Core Richland and Sequim Lab-Field Team/Data Generation and Files/ECA/FTICR/03_ProcessedData/CoreMS/EC_Data_Processed_FTICR/Processed_with_XML/'
github = 'C:/Users/gara009/OneDrive - PNNL/Documents/GitHub/ECA_DOM_Thermodynamics/'

# ====== Read in and clean up data ======
# Processed ICR Data
data = read.csv(list.files(path = github, pattern = "*comparison_results", full.names = T))

# Read in summary file from the data package
sample_data = read_csv(paste0(github,'EC_Data_Package/Sample_Data/EC_Sediment_Sample_Data_Summary.csv'),comment = '#', na = c('N/A', -9999)) %>%
  slice(-(1:11))%>%
  mutate_at(vars(-Sample_Name,-Field_Name,-IGSN,-Material), as.numeric)

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

# ====== Looking at ratios ======
summ = data %>%
  group_by(site, status) %>%
  summarise(Number_of_Formulas = n()) %>%
  filter(status != 'shared') %>%
  spread(key = status, value = Number_of_Formulas) %>%
  mutate(Ratio_Wet_to_Dry = `unique_to_wet` / `unique_to_dry`,
         Category = ifelse(Ratio_Wet_to_Dry < 1, "Lower_than_1", "Higher_than_1")) %>%
  select(site, Ratio_Wet_to_Dry, Category)


ggplot(summ, aes(x = Ratio_Wet_to_Dry)) +
  geom_histogram(color = "black", fill = "skyblue", bins = 30) +
  geom_vline(aes(xintercept = 1), linetype = "dashed", color = "red") +
  labs(title = "Ratio of unique formulas per treatment", 
       x = "Ratio of Wet to Dry", 
       y = "Count") +
  theme_bw()

summary = summ %>%
  group_by(Category) %>%
  summarise(Total = n())

# ==== Looking at some correlations ====
df = merge(summ,effect_size, by = 'site') %>%
  mutate(Cubic_Root_Effect_Size_Respiration = abs(Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H)^(1/3),
         Moment = ifelse(Cubic_Root_Effect_Size_Respiration < 5, "Neutral_Moment", "Cold_Moment"))


ggplot(df, aes(x = Category, y = Cubic_Root_Effect_Size_Respiration, fill = Moment)) +
  geom_boxplot(color = "black") +
  labs(title = "Effect Size of Respiration Rate by Category", 
       x = "Category", 
       y = "Cubic Root Effect Size Respiration (mg DO per kg per H)")  +
  theme_bw() + theme(legend.position = 'top')


ggplot(df, aes(x = Ratio_Wet_to_Dry, y = Cubic_Root_Effect_Size_Respiration, color = Moment)) +
  geom_point(size = 3, alpha = 0.7) +
  labs(title = "Effect Size of Respiration Rate by Category", 
       x = "Category", 
       y = "Cubic Root Effect Size Respiration (mg DO per kg per H)") +
  theme_bw() +
  theme(legend.position = 'top')

# ====== Other plots =====
df2 <- merge(summ, sample_data, by = 'site', all = TRUE) %>%
  mutate(
    Median_Respiration_Rate_mg_DO_per_kg_per_H = ifelse(        Median_Respiration_Rate_mg_DO_per_kg_per_H == 0,
      -min(abs(Median_Respiration_Rate_mg_DO_per_kg_per_H[
        Median_Respiration_Rate_mg_DO_per_kg_per_H != 0 & 
          !is.na(Median_Respiration_Rate_mg_DO_per_kg_per_H)
      ])) / 2, 
      Median_Respiration_Rate_mg_DO_per_kg_per_H
    ),                                       
    Cubic_Root_Respiration = abs(Median_Respiration_Rate_mg_DO_per_kg_per_H)^(1/3)
  )%>%
  filter(!is.na(Cubic_Root_Respiration))%>% filter(!is.na(Category))


ggplot(df2, aes(x = Category, y = Cubic_Root_Respiration)) +
  geom_boxplot(color = "black") +
  labs(title = "Median Respiration Rates per site by Category", 
       x = "Category", 
       y = "Cubic Root Respiration (mg DO per kg per H)")  +
  theme_bw() + theme(legend.position = 'top')


ggplot(df2, aes(x = Ratio_Wet_to_Dry, y = Cubic_Root_Respiration, color = Treatment)) +
  geom_point(size = 3, alpha = 0.7) +
  labs(title = " ", 
       x = "Unique Formulas W/D", 
       y = "Cubic Root Respiration (mg DO per kg per H)") +
  theme_bw() +
  theme(legend.position = 'top')


num_sites <- length(unique(df2$site))


shapes <- c(0:25, 32:47)  # These are ggplot2's shape codes for a variety of shapes

# If you have more than 38 sites, it will cycle through the shapes
unique_sites <- unique(df2$site)
if(length(unique_sites) > length(shapes)) {
  warning("More sites than available shapes. Shapes will be reused.")
  shapes <- rep(shapes, length.out = length(unique_sites))
}

# Create a named vector for the shapes based on site names
site_shapes <- setNames(shapes, unique_sites)

# Plot with ggplot2
ggplot(df2, aes(x = Median_62948_Final_Gravimetric_Moisture_g_per_g, y = Cubic_Root_Respiration, color = Treatment, shape = site)) +
  geom_point(size = 3, alpha = 0.7) +
  labs(x = "Final Gravimetric Moisture g per g", 
       y = "Cubic Root Respiration (mg DO per kg per H)") +
  scale_shape_manual(values = site_shapes) +
  theme_bw() +
  theme(legend.position = 'top')


ggplot(df2, aes(x = Median_62948_Final_Gravimetric_Moisture_g_per_g, y = Ratio_Wet_to_Dry, color = Treatment, shape = site)) +
  geom_point(size = 3, alpha = 0.7) +
  labs(x = "Final Gravimetric Moisture g per g", 
       y = "Ratio Unique formulas wet to dry") +
  scale_shape_manual(values = site_shapes) +
  theme_bw() +
  theme(legend.position = 'top')


ggplot(df2, aes(y = Median_62948_Final_Gravimetric_Moisture_g_per_g, x = Category)) +
  geom_boxplot() +
  labs(y = "Final Gravimetric Moisture g per g", 
       x = " ") +
  theme_bw() +
  theme(legend.position = 'top')

# ====== Testing some things directly with the ICR data ====
