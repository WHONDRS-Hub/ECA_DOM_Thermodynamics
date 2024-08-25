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
setwd('C:/Users/gara009/OneDrive - PNNL/Documents - Core Richland and Sequim Lab-Field Team/Data Generation and Files/ECA/FTICR/03_ProcessedData/Formularity/EC_Data_Processed_FTICR/')

github = 'C:/Users/gara009/OneDrive - PNNL/Documents/GitHub/ECA_DOM_Thermodynamics/'
# ====== Read in data ======
# Processed ICR Data
data = read.csv(list.files(pattern = "*three_reps_four_cal_points_Data.csv"), row.names = 1)
mol = read.csv(list.files(pattern = "*clean_four_cal_points_Mol.csv"), row.names = 1)
# Fixing colnames 
colnames(data) = gsub('SIR.','SIR-',colnames(data))

# ========= Data processing  ======
# Removing peaks that were not assigned a molecular formula
data = data[!is.na(mol$MolForm),]
mol = mol[!is.na(mol$MolForm),]

# Creating factors sheet
factors = data.frame(Samples = colnames(data), Location = colnames(data), Treatment = colnames(data))
factors$Location = str_extract(factors$Location, "EC_0[0-9]{2}|EC_([A-Za-z0-9]+)")
factors$Treatment = str_extract(factors$Treatment, "W|D|Blk")

# Select the mol variables of interest
mol2 = mol %>% dplyr::select(AI_Mod,DBE_1,NOSC,delGcoxPerCmol,delGd,lamO2)

df.merge = merge(data,mol2, by = 'row.names', all = TRUE)
# if lambda is higher than 0.3 or negative turn into NA
df.merge <- df.merge %>%
  mutate(lamO2 = ifelse(lamO2 < 0 | lamO2 > 0.3, NA, lamO2))

# ==== Calculating median metrics and thermodynamics per sample ====

df.stats = as.data.frame(matrix(NA, nrow = ncol(data), ncol = 7))
colnames(df.stats) = c('Samples','Median_AI_mod','Median_NOSC','Median_DBE','Median_Gibbs','Median_Gibbs_Compound','Median_Lambda')

# Calculate medians
for (i in 2:(ncol(data)+1)){
  df.stats$Samples[i-1] = colnames(df.merge[i])
  df.stats$Median_Gibbs[i-1] = median(df.merge$delGcoxPerCmol[which(df.merge[, i] > 0)], na.rm = T)
  df.stats$Median_Gibbs_Compound[i-1] = median(df.merge$delGd[which(df.merge[, i] > 0)], na.rm = T)
  df.stats$Median_Lambda[i-1] = median(df.merge$lamO2[which(df.merge[, i] > 0)],na.rm = T)
  df.stats$Median_AI_mod[i-1] = median(df.merge$AI_Mod[which(df.merge[, i] > 0)],na.rm = T)
  df.stats$Median_NOSC[i-1] = median(df.merge$NOSC[which(df.merge[, i] > 0)],na.rm = T)
  df.stats$Median_DBE[i-1] = median(df.merge$DBE_1[which(df.merge[, i] > 0)],na.rm = T)
}

df <- merge(df.stats, factors, by = "Samples", all = TRUE)

# ==== Export medians per sample ===
write.csv(df, paste0(github,'Median_metrics_Formularity_3_reps_4_cal_pts.csv'),row.names = F)

# ===== Medians per treatment ====
location = unique(factors$Location)

df.stats.wet = as.data.frame(matrix(NA, nrow = (length(location)), ncol = 9))
colnames(df.stats.wet) = c('Samples','Treatment','Number_of_reps','Median_AI_mod','Median_NOSC','Median_DBE','Median_Gibbs','Median_Gibbs_Compound','Median_Lambda')

for (i in 1:length(location)){
  wet.reps = grep(pattern = paste0(location[i],'_SIR-W'), colnames(df.merge))
  df.mol = df.merge[,557:562]

a <- list()
  for (j in 1:length(wet.reps)){
    a[[j]] <- df.mol[which(df.merge[, wet.reps[j]] > 0),]
    
  }

combined_df <- do.call(rbind, a)
    df.stats.wet$Samples[i] = paste0(location[i],'_SIR-W')
    df.stats.wet$Treatment[i] = 'Wet'
    df.stats.wet$Number_of_reps[i] = length(wet.reps)
    df.stats.wet$Median_Gibbs[i] = median(combined_df$delGcoxPerCmol, na.rm = T)
    df.stats.wet$Median_Gibbs_Compound[i] = median(combined_df$delGd, na.rm = T)
    df.stats.wet$Median_Lambda[i] = median(combined_df$lamO2,na.rm = T)
    df.stats.wet$Median_AI_mod[i] = median(combined_df$AI_Mod,na.rm = T)
    df.stats.wet$Median_NOSC[i] = median(combined_df$NOSC,na.rm = T)
    df.stats.wet$Median_DBE[i] = median(combined_df$DBE_1,na.rm = T)

}


df.stats.dry = as.data.frame(matrix(NA, nrow = (length(location)), ncol = 9))
colnames(df.stats.dry) = c('Samples','Treatment','Number_of_reps','Median_AI_mod','Median_NOSC','Median_DBE','Median_Gibbs','Median_Gibbs_Compound','Median_Lambda')

for (i in 1:length(location)){
  dry.reps = grep(pattern = paste0(location[i],'_SIR-D'), colnames(df.merge))
  df.mol = df.merge[,557:562]
  
  a <- list()
  for (j in 1:length(dry.reps)){
    a[[j]] <- df.mol[which(df.merge[, dry.reps[j]] > 0),]
    
  }
  
  combined_df <- do.call(rbind, a)
  df.stats.dry$Samples[i] = paste0(location[i],'_SIR-D')
  df.stats.dry$Treatment[i] = 'Dry'
  df.stats.dry$Number_of_reps[i] = length(dry.reps)
  df.stats.dry$Median_Gibbs[i] = median(combined_df$delGcoxPerCmol, na.rm = T)
  df.stats.dry$Median_Gibbs_Compound[i] = median(combined_df$delGd, na.rm = T)
  df.stats.dry$Median_Lambda[i] = median(combined_df$lamO2,na.rm = T)
  df.stats.dry$Median_AI_mod[i] = median(combined_df$AI_Mod,na.rm = T)
  df.stats.dry$Median_NOSC[i] = median(combined_df$NOSC,na.rm = T)
  df.stats.dry$Median_DBE[i] = median(combined_df$DBE_1,na.rm = T)
  
}

df_medians = rbind(df.stats.dry,df.stats.wet)

# ==== Export medians per sample ===
write.csv(df_medians, paste0(github,'Median_per_site_Formularity_3_reps_4_cal_pts.csv'),row.names = F)
