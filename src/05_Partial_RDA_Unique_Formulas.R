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
data_path = paste0(github,'Data/')
figure_path = paste0(github,'Figures/')
# ====== Read in data ======
data = read.csv(paste0(data_path,'Medians_of_Weighted_averages_for_molecular_properties_per_site_and_treatment_unique_formulas.csv'))
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
# Check for multicolinearity
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

# Identify and remove high VIF variables
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
variance_explained <- summary(rda_partial)$cont$proportion
print(variance_explained)

# Check for influential samples
cooks.distance <- cooks.distance(rda_partial)
plot(cooks.distance, type = "h", main = "Cooks Distance for Partial RDA")
abline(h = 4/(length(cooks.distance)), col = "red")  # Common threshold

# ==== Plots ===
site_scores <- as.data.frame(scores(rda_partial, display = "sites", scaling = 2))
sites_df <- as.data.frame(site_scores)
sites_df$Treatment <- final_predictors_scaled$Treatment  # Add treatment info

biplot_scores <- as.data.frame(scores(rda_partial, display = "bp", scaling = 2))
signif_vars <- c("Median_62948_Final_Gravimetric_Moisture_g_per_g",
                 "Median_Extractable_NPOC_mg_per_kg")  #

# Filter species_scores for significant variables
species_scores_signif <- biplot_scores[rownames(biplot_scores) %in% signif_vars, ]

# Create a data frame for significant variables
species_df_signif <- as.data.frame(species_scores_signif)
species_df_signif$Variable <- rownames(species_scores_signif)

species_df_signif <- species_df_signif %>%
  mutate(Variable = case_when(
    Variable == "Median_Extractable_NPOC_mg_per_kg" ~ "Extractable NPOC",
    Variable == "Median_62948_Final_Gravimetric_Moisture_g_per_g" ~ "Gravimetric Moisture",
    TRUE ~ Variable  # Keep other variables as they are
  ))

library(vegan)
library(ggplot2)
library(ggrepel)

# Plot Partial RDA using ggplot2
partial_rda_plot <- ggplot() +
  # Plot the sites (samples)
  geom_point(data = sites_df, 
             aes(x = RDA1, y = RDA2, color = Treatment), 
             size = 3, alpha = 0.8) +
  # Plot arrows for significant variables
  geom_segment(data = species_df_signif, 
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(length = unit(0.3, "cm")), 
               color = "black") +
  # Label significant variables
  geom_text_repel(data = species_df_signif, 
                  aes(x = RDA1, y = RDA2, label = Variable),
                  size = 4, 
                  box.padding = 0.3, 
                  point.padding = 0.5,
                  segment.color = "grey50") +
  # Customize colors for treatments
  scale_color_manual(values = c("Dry" = "darkorange", "Wet" = "lightblue")) +
  # Add labels and theme adjustments
  labs(title = " ",
       x = "RDA1 (33.33% Variance Explained)",
       y = "RDA2 (6.95% Variance Explained)",
       color = "Treatment") +
  theme_bw() +
  theme(legend.position = "top",
        plot.title = element_text(hjust = 0.5, face = "bold"))

# Save as PDF
ggsave(paste0(figure_path,"Figure3_Partial-RDA_plots.pdf"), plot = partial_rda_plot, width = 8, height = 8)

# Save as PNG
ggsave(paste0(figure_path,"Figure3_Partial-RDA_plots.png"), plot = partial_rda_plot, width = 8, height = 8)
