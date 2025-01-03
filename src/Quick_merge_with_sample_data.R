rm(list=ls(all=T))

library(dplyr) # For reorganization
library(stringr) # For string manipulation
# ==== Defining paths and working directories ======
github = 'C:/Users/gara009/OneDrive - PNNL/Documents/GitHub/ECA_DOM_Thermodynamics/'
# ====== Read in data ======
sample_data = read_csv(paste0(github,'EC_Data_Package/Sample_Data/EC_Sediment_Sample_Data_Summary.csv'),comment = '#', na = c('N/A', -9999)) %>%
  slice(-(1:11))%>%
  mutate_at(vars(-Sample_Name,-Field_Name,-IGSN,-Material), as.numeric)

grainsize_data = read_csv(paste0(github,'v4_CM_SSS_Data_Package/Sample_Data/v3_CM_SSS_Sediment_Grain_Size.csv'),comment = '#', na = c('N/A', -9999)) %>%
  slice(-(1:11))%>%
  mutate_at(vars(-Sample_Name,-Field_Name,-IGSN,-Material,-Methods_Deviation), as.numeric)%>%
  dplyr::select(-Field_Name,-IGSN,-Material,-Methods_Deviation)

grainsize_data$site = grainsize_data$Sample_Name
grainsize_data$site = gsub('CM','EC',grainsize_data$site)
grainsize_data$site = gsub('_GRN','',grainsize_data$site)

grainsize_data = grainsize_data %>%
              dplyr::select(-Sample_Name)


effect_size = read_csv(paste0(github,'EC_Data_Package/Sample_Data/EC_Sediment_Effect_Size.csv'),comment = '#', na = c('N/A', -9999)) %>%
  slice(-(1:11))%>%
  dplyr::select('Sample_Name',"Effect_Size_Respiration_Rate_mg_DO_per_L_per_H","Effect_Size_Respiration_Rate_mg_DO_per_kg_per_H","Effect_Size_Initial_Gravimetric_Moisture_g_per_g","Effect_Size_Final_Gravimetric_Moisture_g_per_g","Effect_Size_Extractable_NPOC_mg_per_kg","Effect_Size_Extractable_TN_mg_per_L")%>%
  mutate_at(vars(-Sample_Name), as.numeric)


effect_size$site = effect_size$Sample_Name
effect_size$site = gsub('_all','',effect_size$site)

effect_size = merge(effect_size,grainsize_data, by = 'site')

factors = data.frame(Sample_Name = sample_data$Sample_Name, site = sample_data$Sample_Name, Treatment = sample_data$Sample_Name)
factors$site = str_extract(factors$site, "EC_0[0-9]{2}|EC_([A-Za-z0-9]+)")
factors$Treatment = str_extract(factors$Treatment, "W|D|Blk")

factors$Treatment = gsub('W','Wet', factors$Treatment)
factors$Treatment = gsub('D','Dry', factors$Treatment)

sample_data = merge(factors,sample_data, by = 'Sample_Name')
sample_data = merge(sample_data,grainsize_data, by = 'site')

Total_MF = read.csv('Median_Site_Total_MF_order.csv')

unique_MF = read.csv('Median_Site_order.csv')

# ==== Merge data and correlate ===

moisture_ratio = sample_data %>%
  group_by(site,Treatment) %>%
  summarise(
    median_moisture = median(Median_62948_Final_Gravimetric_Moisture_g_per_g, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  spread(key = Treatment, value = median_moisture) %>%  # Spread the data for Wet and Dry
  mutate(
    moisture_ratio = Wet / Dry  # Calculate the ratio of Wet to Dry median moisture
  )

total_mf_merge = merge(Total_MF,moisture_ratio, by = 'site')

unique_mf_merge = merge(unique_MF,moisture_ratio, by = 'site')

# ===== merge with grain size ====

total_mf_merge = merge(Total_MF,grainsize_data, by = 'site')

unique_mf_merge = merge(unique_MF,grainsize_data, by = 'site')

plot((total_mf_merge$ratio)~total_mf_merge$Percent_Tot_Sand)

plot((unique_mf_merge$ratio)~unique_mf_merge$Percent_Tot_Sand)


# ===== Loading medians mol prop =====
mol.prop = read.csv('Medians_of_Weighted_averages_for_molecular_properties_per_site_and_treatment.csv')

mol.merge = merge(mol.prop,grainsize_data, by = 'site')

plot((mol.merge$Median_Weighted_Avg_Lambda),mol.merge$Percent_Tot_Sand)

ggplot(mol.merge, aes(x = Median_Weighted_Avg_Lambda, y = Percent_Tot_Sand, color = Treatment)) +
  geom_point(size = 3, alpha = 0.8) +  # Scatter plot with customized point size and transparency
  labs(
    x = "Median Weighted Average Lambda",
    y = "Percent Total Sand",
    color = "Treatment"
  ) +
  scale_color_manual(values = c("Wet" = "lightblue", "Dry" = "darkorange")) +  # Custom colors for Treatment
  theme_bw() +  # Minimal theme
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.position = "top"
  )
