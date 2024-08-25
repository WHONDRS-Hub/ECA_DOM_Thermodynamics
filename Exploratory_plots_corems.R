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
input = 'C:/Users/gara009/OneDrive - PNNL/Documents - Core Richland and Sequim Lab-Field Team/Data Generation and Files/ECA/FTICR/03_ProcessedData/CoreMS/EC_Data_Processed_FTICR/'

github = 'C:/Users/gara009/OneDrive - PNNL/Documents/GitHub/ECA_DOM_Thermodynamics/'
out_plots = paste0(github,'CoreMS/Plots/')
out_data = paste0(github,'CoreMS/Data/')

# ====== Read in data ======
# Processed ICR Data
data = read.csv(list.files(pattern = "*Processed_Data.csv", path = input, full.names = T), row.names = 1)
mol = read.csv(list.files(pattern = "*Processed_Mol.csv", path = input, full.names = T), row.names = 1)

# Load the mapping file
mapping =  read.csv(list.files(pattern = "*Mapping.csv", path = input, full.names = T),)

# Fixing and Mapping sample names
colnames(data) = gsub('SIR.','SIR-',colnames(data))
colnames(data) = gsub('_[0-9]+\\.corems','',colnames(data))
colnames(data) = gsub('.corems','',colnames(data))


for (i in 141:ncol(data)) {
  # Find the corresponding Randomized_ID for the current column name
  names(data)[i] <- mapping$Sample_ID[which(mapping$Randomized_ID == names(data)[i])]
}
for (i in 1:140){
  names(data)[i] = ifelse(grepl("EC_[0-9]{3}_SIR",   names(data)[i]),   names(data)[i], sub("EC_([0-9]{1,2})_SIR", "EC_0\\1_SIR",   names(data)[i]))
}
# ===== Data cleaning ==============
# Factor creation
factors = data.frame(Samples = colnames(data), 
                     Sample_Type = case_when(grepl("Blk", colnames(data)) ~ "Process Blank",
                                             grepl("W", colnames(data)) ~ "Wet",
                                             grepl("D", colnames(data)) ~ "Dry",TRUE~"Samples"))


# Converting data to presence/absence
int.data = data
data[data > 0] = 1


#Generate clean data

#
### Moving onto the permissive dataset
# Peaks in >50% of the blanks
peaks.to.remove = row.names(data)[(rowSums(data[,which(factors$Sample_Type %in% "Process Blank")])/
                                     ncol(data[,which(factors$Sample_Type %in% "Process Blank")])) > 0.5]

# Removing peaks
# data_clean = data[-which(row.names(data) %in% peaks.to.remove),]
# mol_clean = mol[-which(row.names(mol) %in% peaks.to.remove),]
data_clean = data
mol_clean = mol

# Confirming order
if(!identical(row.names(data_clean), row.names(mol_clean))){
  stop("Your conservative dataset has row name issues.")
}

# Remove process blanks from the data
data_clean = data_clean[,-grep("Blk", colnames(data_clean))]

# ===== Some quick checks ======
number_of_peaks <- sapply(data, function(column) sum(column != 0))
total_peaks_data <- data.frame(Samples = names(number_of_peaks), Total_Peaks = number_of_peaks)

number_of_peaks_clean <- sapply(data_clean, function(column) sum(column != 0))
total_peaks_data_clean <- data.frame(Samples = names(number_of_peaks_clean), Total_Peaks_Clean = number_of_peaks_clean)

df = merge(total_peaks_data,total_peaks_data_clean, by = 'Samples', all = T)
df$Difference = df$Total_Peaks- df$Total_Peaks_Clean

df2 = merge(df,factors, by = 'Samples', all = T)
df2$Location = str_extract(df2$Samples, "EC_[A-Z0-9]+")

# ==== PCA =====
factors$Location = str_extract(factors$Samples, "EC_[A-Z0-9]+")
factors$Treatment = factors$Samples
factors$Treatment = str_extract(factors$Treatment, "W|D|Blk")

pca = prcomp(x = t(data)) # Calculating PCA
pca = as.data.frame(scores(pca)) # Converting to PCA scores in order to plot using ggplot
pca = cbind(factors, pca)


pca %>%
  ggplot(aes(x = PC1, y = PC2))+
  geom_point(aes(color =  Treatment),size = 2) +
  # geom_point() +
  #geom_label(aes(label = row.names(pca)))+
  scale_shape_manual(values = c(0,1,2,3,4,5,6,7,8,9,10))+ theme_bw()

# xlim(-25, 35) + ylim(-25, 35)+
# hori_x_theme 
#+ theme(legend.position = "top")+
#scale_color_manual(values=as.vector(alphabet(57)))
ggsave(paste0("PCA_all_sites","_",Sys.Date(),".pdf"))
