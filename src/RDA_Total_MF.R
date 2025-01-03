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
library(vegan)
# ==== Defining paths and working directories ======
github = 'C:/Users/gara009/OneDrive - PNNL/Documents/GitHub/ECA_DOM_Thermodynamics/'
# ====== Read in data ======
data = read.csv('Medians_of_Weighted_averages_for_molecular_properties_per_site_and_treatment_Total_Formulas.csv')

df = read.csv('Median_Site_Total_MF_order.csv') %>%
  dplyr::select(site,Dry,Wet) %>%
  pivot_longer(
    cols = c("Wet","Dry"),  # Choose columns to pivot 
    names_to = "Treatment",            
    values_to = "Median_Number_of_Molecular_Formulas"                  # New column for the values of the properties
  )

data = merge(data,df, by = c('site','Treatment'))
row.names(data) = paste0(data$site,'_',data$Treatment)

sample_data = read_csv(paste0(github,'EC_Data_Package/Sample_Data/EC_Sediment_Sample_Data_Summary.csv'),comment = '#', na = c('N/A', -9999)) %>%
  slice(-(1:11)) %>%  # Remove the first 11 rows
  mutate_at(vars(-Sample_Name, -Field_Name, -IGSN, -Material), as.numeric) %>%  # Convert specified columns to numeric
  mutate(
    Treatment = case_when(
      str_detect(Sample_Name, "-W") ~ "Wet",    # If Sample_Name contains "-W", assign "Wet"
      str_detect(Sample_Name, "-D") ~ "Dry",    # If Sample_Name contains "-D", assign "Dry"
      TRUE ~ NA_character_                       # Otherwise, assign NA
    )
  )

sample_data$site = gsub('-W|-D','',sample_data$Sample_Name)


grainsize_data = read_csv(paste0(github,'v4_CM_SSS_Data_Package/Sample_Data/v3_CM_SSS_Sediment_Grain_Size.csv'),comment = '#', na = c('N/A', -9999)) %>%
  slice(-(1:11))%>%
  mutate_at(vars(-Sample_Name,-Field_Name,-IGSN,-Material,-Methods_Deviation), as.numeric)%>%
  dplyr::select(-Field_Name,-IGSN,-Material,-Methods_Deviation)

grainsize_data$site = grainsize_data$Sample_Name
grainsize_data$site = gsub('CM','EC',grainsize_data$site)
grainsize_data$site = gsub('_GRN','',grainsize_data$site)

grainsize_data = grainsize_data %>%
  dplyr::select(site,Percent_Fine_Sand)
# ==== Set up data =====
response_data = data %>%
  dplyr::filter(site != "EC_023") %>%  # Remove rows where site is EC_023 because it has NA in some of the explanatory data
  dplyr::select(-site,-Treatment)

explanatory_data = data %>%
  dplyr::select(site,Treatment)

explanatory_data = merge(explanatory_data,sample_data, by = c('site','Treatment'))

explanatory_data = explanatory_data %>%
  dplyr::select(-Field_Name,-Sample_Name,-IGSN, -Material,-"Median_Respiration_Rate_mg_DO_per_L_per_H", -"Median_62948_Initial_Gravimetric_Moisture_g_per_g",- "Median_Fe_mg_per_L",-"Median_ATP_nanomoles_per_L",-"Median_Extractable_NPOC_mg_per_L",-"Median_Extractable_TN_mg_per_L",-"Median_Missing_Reps")

explanatory_data = merge(explanatory_data,grainsize_data, by = 'site')

explanatory_data$Treatment <- as.factor(explanatory_data$Treatment)
explanatory_data$site <- as.factor(explanatory_data$site)

# ===== Check for co-linearity =====
# To help with model overfitting for the number of samples we have, assess co-linearity and pick between 5-8 explanatory variables to put into the model. 
# Check for multicollinearity
vif_values <- vif(lm(Median_SpC_microsiemens_per_cm ~ Median_pH + Median_Temperature_degC + 
                       Median_Respiration_Rate_mg_DO_per_kg_per_H + Median_62948_Final_Gravimetric_Moisture_g_per_g + 
                       Median_Fe_mg_per_kg + Median_ATP_picomoles_per_g + Median_Extractable_NPOC_mg_per_kg + 
                       Median_Extractable_TN_mg_per_kg + Median_01395_C_percent_per_mg + 
                       Median_01397_N_percent_per_mg + Percent_Fine_Sand, data = explanatory_data))
print(vif_values)


# Need to keep variables with vif< 5

library(corrplot)
corr_matrix <- cor(select(explanatory_data, Median_SpC_microsiemens_per_cm, Median_pH, 
                          Median_Temperature_degC, Median_Respiration_Rate_mg_DO_per_kg_per_H, 
                          Median_62948_Final_Gravimetric_Moisture_g_per_g, Median_Fe_mg_per_kg, 
                          Median_ATP_picomoles_per_g, Median_Extractable_NPOC_mg_per_kg, 
                          Median_Extractable_TN_mg_per_kg, Median_01395_C_percent_per_mg, 
                          Median_01397_N_percent_per_mg, Percent_Fine_Sand), use = "complete.obs")
corrplot(corr_matrix, method = "color", addCoef.col = "black", number.cex = 0.7)

# Identify and remove high VIF variables and also other variables that might not have theoretical relevance to avoid model over fitting 
predictors_reduced <- explanatory_data %>%
  select(-site, -Treatment,
         -Median_Respiration_Rate_mg_DO_per_kg_per_H,
         -Median_Temperature_degC,
         -Median_pH,
         -Median_SpC_microsiemens_per_cm)

lm_reduced <- lm(Median_62948_Final_Gravimetric_Moisture_g_per_g  ~ ., data = predictors_reduced)
vif_reduced <- vif(lm_reduced)
print(vif_reduced)

final_predictors <- explanatory_data %>%
  select(site,Treatment,colnames(predictors_reduced))%>%
  drop_na() # site 23 has NA in some of geochem data

final_predictors_scaled <- final_predictors %>%
  mutate_at(vars(-Treatment, -site), ~ scale(.))


# ==== Run RDA =====
# Partial RDA controlling for site

rda_partial <- rda(response_data ~ Treatment + Median_62948_Final_Gravimetric_Moisture_g_per_g + Median_Fe_mg_per_kg + Median_ATP_picomoles_per_g + Median_Extractable_NPOC_mg_per_kg + Median_Extractable_TN_mg_per_kg + Median_01395_C_percent_per_mg + Median_01397_N_percent_per_mg + Percent_Fine_Sand + Condition(site), data = final_predictors_scaled)

summary(rda_partial)
anova(rda_partial, permutations = 999)
anova(rda_partial, by = "axis", permutations = 999)  # Axis significance
anova(rda_partial, by = "terms", permutations = 999)  # Variable significance

# Plot RDA
plot(rda_partial, scaling = 2, main = "Partial RDA Biplot")


# Calculate the proportion of variance explained
variance_explained <- summary(rda_partial)$cont$importance


# Check for influential samples
cooks.distance <- cooks.distance(rda_partial)
plot(cooks.distance, type = "h", main = "Cooks Distance for Partial RDA")
abline(h = 4/(length(cooks.distance)), col = "red")  # Common threshold


site_scores <- as.data.frame(scores(rda_partial, display = "sites", scaling = 2))
site_scores$Treatment <- final_predictors_scaled$Treatment  # Add treatment info

biplot_scores <- as.data.frame(scores(rda_partial, display = "bp", scaling = 2))
biplot_scores$Variable <- rownames(biplot_scores)  # Add variable names

ggplot() +
  # Plot site scores (samples)
  geom_point(data = site_scores, aes(x = RDA1, y = RDA2, color = Treatment), size = 3) +
  # Add arrows for biplot variables
  geom_segment(data = biplot_scores, aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  # Add labels for biplot variables
  geom_text(data = biplot_scores, aes(x = RDA1, y = RDA2, label = Variable),
            hjust = 1.2, vjust = 1.2, color = "black") +
  # Customize the theme
  theme_bw() +
  labs(title = "Partial RDA Biplot Colored by Treatment",
       x = "RDA1", y = "RDA2") +
  scale_color_manual(values = c("Dry" = "darkorange", "Wet" = "lightblue"))  # Adjust colors

# ===== Explore interactions ======
rda_with_interaction <- rda(response_data ~ Median_62948_Final_Gravimetric_Moisture_g_per_g * Treatment + 
                              Median_Fe_mg_per_kg + 
                              Median_ATP_picomoles_per_g + 
                              Median_Extractable_NPOC_mg_per_kg + 
                              Median_Extractable_TN_mg_per_kg + 
                              Median_01395_C_percent_per_mg + 
                              Median_01397_N_percent_per_mg + 
                              Percent_Fine_Sand + 
                              Condition(site), 
                            data = final_predictors_scaled)

summary(rda_with_interaction)
anova(rda_with_interaction, permutations = 999)

# ===== Separate the treatments =====
# Dry Treatments
dry_data <- final_predictors_scaled %>%
  filter(Treatment == "Dry")
response_data_dry =  data %>%
  filter(Treatment == "Dry") %>%
  dplyr::filter(site != "EC_023") %>%  # Remove rows where site is EC_023 because it has NA in some of the explanatory data
  dplyr::select(-site,-Treatment)

rda_dry <- rda(response_data_dry ~ Median_62948_Final_Gravimetric_Moisture_g_per_g + 
                 Median_Fe_mg_per_kg + 
                 Median_Extractable_NPOC_mg_per_kg + 
                 Percent_Fine_Sand + 
                 Condition(site), 
               data = dry_data)

plot(rda_dry, scaling = 2, main = "Partial RDA - Dry Treatments")

# Wet Treatments
wet_data <- final_predictors_scaled %>%
  filter(Treatment == "Wet")

rda_wet <- rda(response_data ~ Median_62948_Final_Gravimetric_Moisture_g_per_g + 
                 Median_Fe_mg_per_kg + 
                 Median_ATP_picomoles_per_g + 
                 Median_Extractable_NPOC_mg_per_kg + 
                 Median_Extractable_TN_mg_per_kg + 
                 Median_01395_C_percent_per_mg + 
                 Median_01397_N_percent_per_mg + 
                 Percent_Fine_Sand + 
                 Condition(site), 
               data = wet_data)

plot(rda_wet, scaling = 2, main = "Partial RDA - Wet Treatments")
