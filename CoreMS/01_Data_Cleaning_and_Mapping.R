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

### Moving onto the permissive dataset
# Peaks in >50% of the blanks
peaks.to.remove = row.names(data)[(rowSums(data[,which(factors$Sample_Type %in% "Process Blank")])/
                                     ncol(data[,which(factors$Sample_Type %in% "Process Blank")])) > 0.25]

# Removing peaks
data_clean_int = int.data[-which(row.names(int.data) %in% peaks.to.remove),]
data_clean = data[-which(row.names(data) %in% peaks.to.remove),]
mol_clean = mol[-which(row.names(mol) %in% peaks.to.remove),]

# Confirming order
if(!identical(row.names(data_clean), row.names(mol_clean))){
  stop("Your conservative dataset has row name issues.")
}


# ===== Some quick checks ======
number_of_peaks <- sapply(data, function(column) sum(column != 0))
total_peaks_data <- data.frame(Samples = names(number_of_peaks), Total_Peaks = number_of_peaks)

number_of_peaks_clean <- sapply(data_clean, function(column) sum(column != 0))
total_peaks_data_clean <- data.frame(Samples = names(number_of_peaks_clean), Total_Peaks_Clean = number_of_peaks_clean)

df = merge(total_peaks_data,total_peaks_data_clean, by = 'Samples', all = T)
df$Difference = df$Total_Peaks- df$Total_Peaks_Clean

df2 = merge(df,factors, by = 'Samples', all = T)
names(df2)[2] = 'Total_Peaks_CoreMS_Raw'
names(df2)[3] = 'Total_Peaks_Clean_CoreMS_Raw'
names(df2)[4] = 'Difference_CoreMS'

write.csv(df2, paste0(github,'CoreMS_Raw_Peaks.csv'), row.names = F)

# ===== Load Total Peaks formularity =====
peaks = read.csv(paste0(github,'Total_peaks_Formularity.csv'))
test = merge(peaks,df2, by = c('Samples','Sample_Type'))

library(viridis)
library(ggpmisc)
library(segmented)

lm_model <- lm(Total_Peaks_Clean ~ Total_Peaks_Clean_CoreMS, data = test)

# Get the p-value
p_value <- summary(lm_model)$coefficients[2, 4]
breaks <- seq(min(test$final), max(test$final), by = 30)

ggplot(test, aes(y = Total_Peaks_Clean, x = Total_Peaks_Clean_CoreMS, color = final)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +  
  labs(y = "Total Peaks Clean Formularity", x = "Total Peaks Clean CoreMS", color = "Final Calibration points Formularity") +
  scale_color_viridis(option = "D", direction = -1, breaks = breaks) +  # Continuous color scale with specified breaks and viridis color palette
  theme_bw() +
  geom_smooth(method = "lm", se = FALSE, aes(group = 1), formula = y ~ x) +  # Adding linear regression line
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~~")), 
               label.x = "left", label.y = 60, size = 4)+
  annotate("text", x = 2000, y = 8000, 
           label = paste("p-value:", signif(p_value, digits = 3)), vjust = 1, hjust = 1, color = "black")


# ====== Remove samples below a Calibration threshold ======

test2 = filter(test, test$final >15)

lm_model <- lm(Total_Peaks_Clean ~ Total_Peaks_Clean_CoreMS, data = test2)

# Get the p-value
p_value <- summary(lm_model)$coefficients[2, 4]
ggplot(test2, aes(y = Total_Peaks_Clean, x = Total_Peaks_Clean_CoreMS, color = final)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +  
  labs(y = "Total Peaks Clean Formularity", x = "Total Peaks Clean CoreMS", color = "Final Calibration points Formularity") +
  scale_color_viridis(option = "D", direction = -1, breaks = breaks) +  # Continuous color scale with specified breaks and viridis color palette
  theme_bw() +
  geom_smooth(method = "lm", se = FALSE, aes(group = 1), formula = y ~ x) +  # Adding linear regression line
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~~")), 
               label.x = "left", label.y = 60, size = 4)+
  annotate("text", x = 2000, y = 8000, 
           label = paste("p-value:", signif(p_value, digits = 3)), vjust = 1, hjust = 1, color = "black")


ggplot(test2, aes(y = Total_Peaks_Clean, x = Total_Peaks_Clean_CoreMS, color = as.factor(final))) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +  
  labs(y = "Total Peaks Clean Formularity", x = "Total Peaks Clean CoreMS", color = "Final Calibration points Formularity") +
  theme_bw() +
  geom_vline(xintercept = 1000, linetype = "dotted", color = "red")


# ==== Export Data =====

# Remove process blanks from the data
data_clean = data_clean[,-grep("Blk", colnames(data_clean))]
data_clean_int = data_clean_int[,-grep("Blk", colnames(data_clean_int))]

write.csv(data_clean_int,paste0(input,'Processed_CoreMS_Clean_Data.csv'))

write.csv(mol_clean,paste0(input,'Processed_CoreMS_Clean_Mol.csv'))
