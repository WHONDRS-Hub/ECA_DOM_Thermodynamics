# ==== Loading libraries =========
rm(list=ls(all=T))

library(stringr); library(devtools);  library("plyr")
library("readr"); library(tidyverse); library(readxl);library(crayon); library(vegan)
# Load in necessary libraries first
library(reshape2)
library(ggpubr) # For to combine plots
library(dplyr) # For reorganization
library(stringr) # For string manipulation
# ==== Defining paths and working directories ======
setwd('C:/Users/gara009/OneDrive - PNNL/Documents - Core Richland and Sequim Lab-Field Team/Data Generation and Files/ECA/FTICR/03_ProcessedData/EC_Data_Processed_FTICR/')

github = 'C:/Users/gara009/OneDrive - PNNL/Documents/GitHub/ECA_DOM_Thermodynamics/'
# ====== Read in data ======
# Processed ICR Data
data = read.csv(list.files(pattern = "*four_reps_Intensity_Data.csv"), row.names = 1)
mol = read.csv(list.files(pattern = "*clean_Intensity_Mol.csv"), row.names = 1)
# Fixing colnames 
colnames(data) = gsub('SIR.','SIR-',colnames(data))

effect.size = read.csv(paste0(github,'ECA_Effect_Size_ReadyForBoye_2023-11-08.csv'))

effect.size$Effect_Size = gsub('-9999',NA,effect.size$Effect_Size)
effect.size$Effect_Size = as.numeric(effect.size$Effect_Size)
effect.size$Location = str_extract(effect.size$Sample_Name, "EC_0[0-9]{2}|EC_([A-Za-z0-9]+)")
effect.size$Moments = NA
# Assign values based on the condition to the new column
effect.size$Moments[effect.size$Effect_Size < 70] <- "Neutral"
effect.size$Moments[effect.size$Effect_Size >= 70] <- "Cold"

effect = effect.size %>% dplyr::select(Location,Moments)

# Load rates
rates = read.csv(paste0(github,'rates_clean.csv'))
rates$Samples = rates$Sample_Name
rates$Samples = gsub('INC', 'SIR', rates$Samples)

# ========= Plots ======
# Removing peaks that were not assigned a molecular formula
data = data[!is.na(mol$MolForm),]
mol = mol[!is.na(mol$MolForm),]

# Creating factors sheet
factors = data.frame(Samples = colnames(data), Location = colnames(data), Treatment = colnames(data))
factors$Location = str_extract(factors$Location, "EC_0[0-9]{2}|EC_([A-Za-z0-9]+)")
factors$Treatment = str_extract(factors$Treatment, "W|D|Blk")

# Select the mol variables of interest
mol2 = mol %>% dplyr::select(AI_Mod,DBE_1,NOSC,delGcoxPerCmol,lamO2)

# Calculating weighed avg metris and thermodynamics per sample
df.merge = merge(data,mol2, by = 'row.names', all = TRUE)

df.stats = as.data.frame(matrix(NA, nrow = ncol(data), ncol = 6))
colnames(df.stats) = c('Samples','AI_mod','NOSC','DBE','Gibbs','Lambda')

for (i in 2:(ncol(data)+1)){
  df.stats$Samples[i-1] = colnames(df.merge[i])
  df.stats$Gibbs[i-1] = weighted.mean(df.merge$delGcoxPerCmol[which(df.merge[, i] > 0)],df.merge[,i][which(df.merge[, i] > 0)])
  df.stats$Lambda[i-1] = weighted.mean(df.merge$lamO2[which(df.merge[, i] > 0)],df.merge[,i][which(df.merge[, i] > 0)])
  df.stats$AI_mod[i-1] = weighted.mean(df.merge$AI_Mod[which(df.merge[, i] > 0)],df.merge[,i][which(df.merge[, i] > 0)])
  df.stats$NOSC[i-1] = weighted.mean(df.merge$NOSC[which(df.merge[, i] > 0)],df.merge[,i][which(df.merge[, i] > 0)])
  df.stats$DBE[i-1] = weighted.mean(df.merge$DBE_1[which(df.merge[, i] > 0)],df.merge[,i][which(df.merge[, i] > 0)])
}

df <- merge(merge(df.stats, factors, by = "Samples", all = TRUE),rates, by = "Samples", all = TRUE)



df = na.omit(df[!is.na(df$Gibbs), ])
df = na.omit(df[!is.na(df$Respiration_Rate_mg_DO_per_L_per_H), ])

ggplot(df, aes(x = Location, y = Gibbs, fill = Treatment)) +
  geom_boxplot() +
  labs(x = "Location", y = "Gibbs") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

df$log_rate = log10(df$Respiration_Rate_mg_DO_per_L_per_H+1)

# Calculate median per Location and Treatment
median_data <- df %>%
  group_by(Location, Treatment) %>%
  summarise_at(median)

# Calculate ratio between W and D
median_data <- median_data %>%
  mutate(ratio_W_to_D = W / D)

# Plot histogram of the ratios
ggplot(median_data, aes(x = ratio_W_to_D)) +
  geom_histogram(binwidth = 0.1, fill = "skyblue", color = "black") +
  labs(x = "Ratio of W to D", y = "Frequency") +
  theme_minimal()

cols_of_interest <- 2:6  # Columns of interest
plots <- list()  # Initialize a list to store plots

for (col in cols_of_interest) {
  plot <- ggplot(df, aes_string(x = colnames(df)[col], y = "log_rate", color = "Treatment")) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) +  # Add linear regression line
    facet_wrap(~ Treatment) +  # Separate plots by Treatment
    labs(x = colnames(df)[col], y = "Log Rate") +  # Label axes
    theme_minimal()
  
  plots[[as.character(col)]] <- plot  # Store each plot in the list
}

# Composite plot with all individual plots
composite_plot <- cowplot::plot_grid(plotlist = plots, align = "hv")

df2 = subset(df, df$log_rate<= 1)
for (col in cols_of_interest) {
  plot <- ggplot(df2, aes_string(x = colnames(df)[col], y = "log_rate", color = "Treatment")) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) +  # Add linear regression line
    facet_wrap(~ Treatment) +  # Separate plots by Treatment
    labs(x = colnames(df2)[col], y = "Log Rate") +  # Label axes
    theme_minimal()
  
  plots[[as.character(col)]] <- plot  # Store each plot in the list
}

# Composite plot with all individual plots
composite_plot <- cowplot::plot_grid(plotlist = plots, align = "hv")


for (col in cols_of_interest) {
  plot <- ggplot(df2, aes_string(x = colnames(df)[col], y = 'Respiration_Rate_mg_DO_per_L_per_H', color = "Treatment")) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) +  # Add linear regression line
    facet_wrap(~ Treatment) +  # Separate plots by Treatment
    labs(x = colnames(df2)[col], y = "Rate") +  # Label axes
    theme_minimal()
  
  plots[[as.character(col)]] <- plot  # Store each plot in the list
}

# Composite plot with all individual plots
composite_plot <- cowplot::plot_grid(plotlist = plots, align = "hv")
