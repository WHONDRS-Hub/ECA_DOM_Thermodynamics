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
# OM data
om = read.csv('Median_Thermodynamics_DOM.csv')
# all the other data 
data = read.csv(paste0(github,'ECA_Means_ESS_PI.csv'))
names(data)[1] = 'Location'
names(data)[2] = 'Treatment'
# ==== Merging by Site and Treatment ====
df = merge(om,data, by = c('Location','Treatment'))

# ========= Correlations ===========
# Calculate ranks
df$rank_Mean_OutliersRemoved_Respiration_Rate_mg_DO_per_kg_per_H <- rank(df$Mean_OutliersRemoved_Respiration_Rate_mg_DO_per_kg_per_H)
df$rank_Gibbs_per_C <- rank(df$Gibbs_per_C)
df$rank_Gibbs_per_compound <- rank(df$Gibbs_per_compound)
df$rank_Lambda <- rank(df$Lambda)
df$rank_ATP <- rank(df$Mean_ATP_picomol_per_g)

# Calculate regression and add to data frame
model <- lm(rank_Mean_OutliersRemoved_Respiration_Rate_mg_DO_per_kg_per_H ~ rank_Gibbs_per_C, data = df)
df$predicted <- predict(model, df)

p <- ggplot(df, aes(x = rank_Mean_OutliersRemoved_Respiration_Rate_mg_DO_per_kg_per_H, y = rank_Gibbs_per_C)) +
  geom_point() +
  stat_cor(method = "spearman",cor.coef.name = c( "rho"),
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")), 
           label.x = 17,label.y = 36,color='black',size=4)+ 
  labs(title = " ",
       x = "Rank of Mean_OutliersRemoved_Respiration_Rate_mg_DO_per_kg_per_H",
       y = "Rank of Gibbs_per_C ") +
  facet_wrap(~ Treatment)+ theme_bw()

p1 <- ggplot(df, aes(x = rank_Mean_OutliersRemoved_Respiration_Rate_mg_DO_per_kg_per_H, y = rank_Gibbs_per_compound)) +
  geom_point() +
  stat_cor(method = "spearman",cor.coef.name = c( "rho"),
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")), 
           label.x = 17,label.y = 36,color='black',size=4)+ 
  labs(title = " ",
       x = "Rank of Mean_OutliersRemoved_Respiration_Rate_mg_DO_per_kg_per_H",
       y = "Rank of Gibbs_per_Comp ") +
  facet_wrap(~ Treatment)+ theme_bw()

p2 <- ggplot(df, aes(x = rank_Mean_OutliersRemoved_Respiration_Rate_mg_DO_per_kg_per_H, y = rank_Lambda)) +
  geom_point() +
  stat_cor(method = "spearman",cor.coef.name = c( "rho"),
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")), 
           label.x = 17,label.y = 36,color='black',size=4)+ 
  labs(title = " ",
       x = "Rank of Mean_OutliersRemoved_Respiration_Rate_mg_DO_per_kg_per_H",
       y = "Rank of Lambda ") +
  facet_wrap(~ Treatment)+ theme_bw()


p3 <- ggplot(df, aes(x = rank_Mean_OutliersRemoved_Respiration_Rate_mg_DO_per_kg_per_H, y = rank_ATP)) +
  geom_point() +
  stat_cor(method = "spearman",cor.coef.name = c( "rho"),
           aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")), 
           label.x = 17,label.y = 36,color='black',size=4)+ 
  labs(title = " ",
       x = "Rank of Mean_OutliersRemoved_Respiration_Rate_mg_DO_per_kg_per_H",
       y = "Rank of ATP ") +
  facet_wrap(~ Treatment)+ theme_bw()


