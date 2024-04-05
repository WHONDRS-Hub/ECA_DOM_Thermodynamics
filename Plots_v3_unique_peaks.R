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


# ====== Read in data ======
# Processed ICR Data
data = read.csv(list.files(pattern = "*four_reps_Data.csv"), row.names = 1)
mol = read.csv(list.files(pattern = "*clean_Mol.csv"), row.names = 1)
# Fixing colnames 
colnames(data) = gsub('SIR.','SIR-',colnames(data))

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

# If a peak is present in at least 50% of the wet samples and not present in at least 50% of the dry samples is considered unique
data.wet = data[which(grepl('W',colnames(data)),)]
data.dry = data[which(grepl('D',colnames(data)),)]

unique.wet = as.data.frame(row.names(data.wet)[(rowSums(data.wet)/ncol(data.wet)) > 0.5])
names(unique.wet)[1] = 'Wet'
unique.dry = as.data.frame(row.names(data.dry)[(rowSums(data.dry)/ncol(data.dry)) > 0.5])
names(unique.dry)[1] = 'Dry'

unique.wet$Keep = NA
for (i in 1:nrow(unique.wet)){
  if(unique.wet$Wet[i] %in% unique.dry$Dry == TRUE){
    unique.wet$Keep[i] = FALSE
  }else{
    unique.wet$Keep[i] = TRUE
  }
}

unique.wet_2 = subset(unique.wet, unique.wet$Keep == TRUE)
row.names(unique.wet_2) = unique.wet_2$Wet

unique.dry$Keep = NA
for (i in 1:nrow(unique.dry)){
  if(unique.dry$Dry[i] %in% unique.wet$Wet == TRUE){
    unique.dry$Keep[i] = FALSE
  }else{
    unique.dry$Keep[i] = TRUE
  }
}

unique.dry_2 = subset(unique.dry, unique.dry$Keep == TRUE)
row.names(unique.dry_2) = unique.dry_2$Dry

mol.wet = merge(mol, unique.wet_2, by = 'row.names')
mol.wet = mol.wet %>% dplyr::select('DBE_1','NOSC','AI_Mod')
mol.w = melt(mol.wet)
mol.w$Treatment = 'W'
mol.dry = merge(mol, unique.dry_2, by = 'row.names')
mol.dry = mol.dry %>% dplyr::select('DBE_1','NOSC','AI_Mod')
mol.d = melt(mol.dry)
mol.d$Treatment = 'D'

melt.mol.by.type = rbind(mol.d,mol.w)
# Plotting metrics
ggplot(melt.mol.by.type, aes(x = value, group = Treatment))+
  geom_density(aes(fill = Treatment), alpha = 0.5)+
  facet_wrap(variable~., scales = "free", ncol = 1)+
  theme_bw() + theme(text = element_text(size=12, color="black"),
                     axis.text = element_text(color = "black"),
                     axis.ticks = element_line(color = "black"),
                     panel.background = element_blank(),
                     panel.grid = element_blank())

