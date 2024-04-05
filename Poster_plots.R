# ==== Loading libraries =========
rm(list=ls(all=T))

library(stringr); library(devtools);  library("plyr")
library("readr"); library(tidyverse); library(readxl);library(crayon); library(vegan)

# ==== Defining paths and working directories ======
setwd('C:/Users/gara009/OneDrive - PNNL/Documents - Core Richland and Sequim Lab-Field Team/Data Generation and Files/ECA/FTICR/03_ProcessedData/EC_Data_Processed_FTICR/')


# ====== Read in data ======
# Processed ICR Data
data = read.csv(list.files(pattern = "*four_reps_Data.csv"), row.names = 1)
mol = read.csv(list.files(pattern = "*Mol.csv"), row.names = 1)
# Fixing colnames 
colnames(data) = gsub('SIR.','SIR-',colnames(data))

# ========= Plots ======

# Selecting only molecular formulae
data = data[!is.na(mol$MolForm),]
mol = mol[!is.na(mol$MolForm),]

# Calculating StoC ratios
mol$StoC_ratio = with(mol, S/C)

# Fixing compound class definitions
mol$bs1_class[grep(";", mol$bs1_class)] = "Multi-class"


# ################################################################ #
#### Creating factors ####
# ################################################################ #
# Creating factors sheet
factors = data.frame(Samples = colnames(data), Location = colnames(data), Treatment = colnames(data))
factors$Location = str_extract(factors$Location, "EC_0[0-9]{2}|EC_([A-Za-z0-9]+)")
factors$Treatment = str_extract(factors$Treatment, "W|D|Blk")

# ####################################################### #
#### Calculating average derived statistics by sample #####
# ####################################################### #

char = data.frame(NOSC = NA, DBE = NA, AI_Mod = NA, NtoC = NA, PtoC = NA, StoC = NA,Peaks = NA, Location = factors$Location,Treatment = factors$Treatment, Samples = factors$Samples,
                  row.names = colnames(data), stringsAsFactors = F)

for(i in 1:ncol(data)){
  temp = mol[which(data[,i] > 0),]
  
  char$NOSC[i] = mean(temp$NOSC, na.rm = T)
  char$DBE[i] = mean(temp$DBE_1, na.rm = T)
  char$AI_Mod[i] = mean(temp$AI_Mod, na.rm = T)
  char$NtoC[i] = mean(temp$NtoC_ratio, na.rm = T)
  char$PtoC[i] = mean(temp$PtoC_ratio, na.rm = T)
  char$StoC[i] = mean(temp$StoC_ratio, na.rm = T)
  char$Peaks[i] = nrow(temp)
  
  rm(temp)
} # Looping through samples and averaging characteristics

rm(i)


################################################# #
#### Calculating elemental composition by sample ####
# ################################################# #

el.comp = matrix(data = 0, nrow = ncol(data), ncol = length(unique(mol$El_comp)), 
                 dimnames = list(colnames(data), unique(mol$El_comp)))

for(i in 1:nrow(el.comp)){
  temp = mol[which(data[,i] > 0),] # Mol data for a given sample
  
  for(j in 1:ncol(el.comp)){
    el.comp[i,j] = length(which(temp$El_comp %in% colnames(el.comp)[j]))
  }
} # Counting the number of times a given elemental composition appears in a dataset

rm(temp, i, j)

el.comp = as.data.frame(t(apply(el.comp, 1, function(x) (x/sum(x))*100)))

# #################################################### #
#### Calculating compound classifications by sample ####
# #################################################### #

comp.class = matrix(data = 0, nrow = ncol(data), ncol = length(unique(mol$bs1_class)), 
                    dimnames = list(colnames(data), unique(mol$bs1_class)))

for(i in 1:nrow(comp.class)){
  temp = mol[which(data[,i] > 0),] # Mol data for a given sample
  
  for(j in 1:ncol(comp.class)){
    comp.class[i,j] = length(which(temp$bs1_class %in% colnames(comp.class)[j]))
  }
} # Counting the number of times a given elemental composition appears in a dataset

rm(temp, i, j)

comp.class = as.data.frame(t(apply(comp.class, 1, function(x) (x/sum(x))*100)))
