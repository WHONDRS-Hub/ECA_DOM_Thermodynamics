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
 
# === Density plots of OM properties ====
temp.data = data[which(grepl('W',colnames(data))),]
temp.mol = mol[,c("AI_Mod", "DBE_1", "NOSC")] 

df.merge = merge(temp.data,temp.mol, by = 'row.names', all = TRUE)

AI_wet = NA
NOSC_wet = NA
DBE_wet = NA
for (i in 1:178){
   temp.AI = df.merge$AI_Mod[which(df.merge[, i] > 0)]
   temp.NOSC = df.merge$NOSC[which(df.merge[, i] > 0)] 
   temp.DBE = df.merge$DBE_1[which(df.merge[, i] > 0)] 
  AI_wet = c(AI_wet,temp.AI)
  NOSC_wet = c(NOSC_wet,temp.NOSC) 
  DBE_wet = c(DBE_wet,temp.DBE) 
  
}
df.wet = as.data.frame(cbind(AI = AI_wet,NOSC = NOSC_wet,DBE = DBE_wet))
df.wet = df.wet[-1, ]
df.mol = melt(as.matrix(df.wet)) # Melting as a matrix to get the Var1/Var2 melt format
df.mol$Treatment = "W" # Adding on a qualifier to know what sample type the data is coming from

melt.mol.by.type = df.mol # Creating the object to eventually go into ggplot

temp.data = data[which(grepl('D',colnames(data))),]
temp.mol = mol[,c("AI_Mod", "DBE_1", "NOSC")] 

df.merge = merge(temp.data,temp.mol, by = 'row.names', all = TRUE)

AI_dry = NA
NOSC_dry = NA
DBE_dry = NA
for (i in 1:178){
  temp.AI = df.merge$AI_Mod[which(df.merge[, i] > 0)]
  temp.NOSC = df.merge$NOSC[which(df.merge[, i] > 0)] 
  temp.DBE = df.merge$DBE_1[which(df.merge[, i] > 0)] 
  AI_dry = c(AI_dry,temp.AI)
  NOSC_dry = c(NOSC_dry,temp.NOSC) 
  DBE_dry = c(DBE_dry,temp.DBE) 
  
}
df.dry = as.data.frame(cbind(AI = AI_dry,NOSC = NOSC_dry,DBE = DBE_dry))
df.dry = df.dry[-1, ]
df.mol = melt(as.matrix(df.dry))
df.mol$Treatment = "D" # Adding on a qualifier to know what sample type the data is coming from

melt.mol.by.type = rbind(melt.mol.by.type, df.mol)

rm(temp.mol)

# Adding in molecular formula-by-sample count into our melt.mol.by.type object
melt.mol.by.type = rbind(melt.mol.by.type, data.frame(Var1 = colnames(data), 
                                                      Var2 = "Molecular Formula Count",
                                                      value = colSums(na.omit(data)), 
                                                      Treatment = factors$Treatment))


# === Repeat fo only thermodynamic properties ======
temp.data = data[which(grepl('W',colnames(data))),]
temp.mol = mol[,c("delGcoxPerCmol", "lamO2", "delGd")] 
temp.mol = temp.mol %>% filter(lamO2 < 0.3) %>%
  mutate(lamO2 = ifelse(lamO2 < 0, 0, lamO2)) %>%
  filter(!is.na(lamO2))
df.merge = merge(temp.data,temp.mol, by = 'row.names', all = TRUE)

GFE_wet = NA
lambda_wet = NA
dGd_wet = NA
for (i in 1:178){
  temp.GFE = df.merge$delGcoxPerCmol[which(df.merge[, i] > 0)]
  temp.lambda = df.merge$lamO2[which(df.merge[, i] > 0)] 
  temp.dGd = df.merge$delGd[which(df.merge[, i] > 0)] 
  GFE_wet = c(GFE_wet,temp.GFE)
  lambda_wet = c(lambda_wet,temp.lambda) 
  dGd_wet = c(dGd_wet,temp.dGd) 
  
}
df.wet = as.data.frame(cbind(Gibbs_per_C_mol = GFE_wet,Gibbs_per_compound_mol = dGd_wet,Lambda = lambda_wet))
df.wet = df.wet[-1, ]
df.mol = melt(as.matrix(df.wet)) # Melting as a matrix to get the Var1/Var2 melt format
df.mol$Treatment = "W" # Adding on a qualifier to know what sample type the data is coming from

melt.mol.by.type_thermo = df.mol # Creating the object to eventually go into ggplot

temp.data = data[which(grepl('D',colnames(data))),]

df.merge = merge(temp.data,temp.mol, by = 'row.names', all = TRUE)

GFE_dry = NA
lambda_dry = NA
dGd_dry = NA
for (i in 1:178){
  temp.GFE = df.merge$delGcoxPerCmol[which(df.merge[, i] > 0)]
  temp.lambda = df.merge$lamO2[which(df.merge[, i] > 0)] 
  temp.dGd = df.merge$delGd[which(df.merge[, i] > 0)] 
  GFE_dry = c(GFE_dry,temp.GFE)
  lambda_dry = c(lambda_dry,temp.lambda) 
  dGd_dry = c(dGd_dry,temp.dGd) 
  
}
df.dry = as.data.frame(cbind(Gibbs_per_C_mol = GFE_dry,Gibbs_per_compound_mol = dGd_dry,Lambda = lambda_dry))
df.dry = df.dry[-1, ]
df.mol = melt(as.matrix(df.dry)) # Melting as a matrix to get the Var1/Var2 melt format
df.mol$Treatment = "D" # Adding on a qualifier to know what sample type the data is coming from

melt.mol.by.type_thermo = rbind(melt.mol.by.type_thermo, df.mol)

rm(temp.mol)

melt.mol.by.type_thermo$Var2 = gsub('Gibbs_per_compound_mol','01_Gibbs_per_compound_mol',melt.mol.by.type_thermo$Var2)

# Plotting metrics
ggplot(melt.mol.by.type, aes(x = value, group = Treatment))+
  geom_histogram(aes(fill = Treatment), alpha = 1, color = 'black')+
  scale_fill_manual(values = c("darkorange","lightblue")) +  # Set colors for W and D
  facet_wrap(Var2~., scales = "free", ncol = 1)+
  theme_bw() + theme(text = element_text(size=18, color="black"),
                     axis.text = element_text(color = "black"),
                     axis.ticks = element_line(color = "black"),
                     panel.background = element_blank(),
                     panel.grid = element_blank())+
  guides(fill = FALSE)+
  labs(x = " ", y = " ")
  
ggplot(melt.mol.by.type_thermo, aes(x = value, group = Treatment))+
  geom_histogram(aes(fill = Treatment), alpha = 1, color = 'black')+
  scale_fill_manual(values = c("darkorange","lightblue")) +  # Set colors for W and D
  facet_wrap(Var2~., scales = "free", ncol = 1)+
  theme_bw() + theme(text = element_text(size=18, color="black"),
                     axis.text = element_text(color = "black"),
                     axis.ticks = element_line(color = "black"),
                     panel.background = element_blank(),
                     panel.grid = element_blank())+
  guides(fill = FALSE)+
  labs(x = " ", y = " ")
