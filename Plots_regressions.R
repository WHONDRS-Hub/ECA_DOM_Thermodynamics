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
data = read.csv(list.files(pattern = "*four_reps_Data.csv"), row.names = 1)
mol = read.csv(list.files(pattern = "*clean_Mol.csv"), row.names = 1)
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

# Renaming bs1_classes
mol$bs1_class[grep(";", mol$bs1_class)] = "Multi-class"

# Adding mass
mol$Mass = as.numeric(as.character(row.names(mol)))

# Creating factors sheet
factors = data.frame(Samples = colnames(data), Location = colnames(data), Treatment = colnames(data))
factors$Location = str_extract(factors$Location, "EC_0[0-9]{2}|EC_([A-Za-z0-9]+)")
factors$Treatment = str_extract(factors$Treatment, "W|D|Blk")

factors = merge(factors,effect, by = 'Location')

# Calculating mean and median thermodynamics per sample
df.merge = merge(data,mol, by = 'row.names', all = TRUE)

df.stats = as.data.frame(matrix(NA, nrow = ncol(data), ncol = 5))
colnames(df.stats) = c('Samples','Gibbs_mean', 'Gibbs_median','Lambda_mean','Lambda_median')


for (i in 2:(ncol(data)+1)){
  df.stats$Samples[i-1] = colnames(df.merge[i])
  df.stats$Gibbs_mean[i-1] = mean(na.omit(df.merge$delGcoxPerCmol[which(df.merge[, i] > 0)]))
  df.stats$Gibbs_median[i-1] = median(na.omit(df.merge$delGcoxPerCmol[which(df.merge[, i] > 0)]))
  df.stats$Lambda_mean[i-1] = mean(na.omit(df.merge$lamO2[which(df.merge[, i] > 0)]))
  df.stats$Lambda_median[i-1] = median(na.omit(df.merge$lamO2[which(df.merge[, i] > 0)]))
}

df <- merge(merge(df.stats, factors, by = "Samples", all = TRUE),rates, by = "Samples", all = TRUE)

df = na.omit(df[!is.na(df$Gibbs_mean), ])
df = na.omit(df[!is.na(df$Respiration_Rate_mg_DO_per_L_per_H), ])

df$log_rate = log10(df$Respiration_Rate_mg_DO_per_L_per_H+1)

cols_of_interest <- 2:5  # Columns of interest
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
