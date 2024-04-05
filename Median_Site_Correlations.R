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

# ========= 