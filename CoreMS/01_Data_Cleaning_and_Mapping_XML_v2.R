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
input = 'C:/Users/gara009/OneDrive - PNNL/Documents - Core Richland and Sequim Lab-Field Team/Data Generation and Files/ECA/FTICR/03_ProcessedData/CoreMS/EC_Data_Processed_FTICR/Processed_with_XML/'

github = 'C:/Users/gara009/OneDrive - PNNL/Documents/GitHub/ECA_DOM_Thermodynamics/'
out_plots = paste0(github,'CoreMS/Plots/')
out_data = paste0(github,'CoreMS/Data/')

# ====== Read in data ======
# Processed ICR Data
data = read.csv(list.files(pattern = "*Processed_Data.csv", path = input, full.names = T), row.names = 1)
mol = read.csv(list.files(pattern = "*Processed_Mol.csv", path = input, full.names = T), row.names = 1)

# Fixing and Mapping sample names
colnames(data) = gsub('SIR.','SIR-',colnames(data))
colnames(data) = gsub('-Blk.','_Blk-',colnames(data))
colnames(data) = gsub('_[0-9]+\\.corems','',colnames(data))
colnames(data) = gsub('.corems','',colnames(data))

# ====== Removing poorly calibrated samples ====
# filter out poorly calibrated samples
cal = read.csv(list.files(pattern = "*Calibration_Results.csv", path = input, full.names = T),check.names = F) %>%
  filter(!is.na(`Cal. Points`) & `Cal. Points` >= 1) %>%
  filter(!`Cal. RMS Error (ppm)` > 1.5) %>%
  dplyr::select(Samples = 'Sample',Cal_pts = 'Cal. Points') %>%
  mutate(Calibration_Status = 'Good')

data_clean = data[,which(colnames(data) %in% cal$Samples)] # Removing poorly calibrated samples but keep the original data present

# Factor creation
factors = data.frame(Samples = colnames(data), 
                     Sample_Type = case_when(grepl("Blk", colnames(data)) ~ "Process Blank",
                                             grepl("W", colnames(data)) ~ "Wet",
                                             grepl("D", colnames(data)) ~ "Dry",TRUE~"Samples"))

factors$Location = str_extract(factors$Samples, "EC_[A-Z0-9]+")


factors = merge(factors,cal, by = 'Samples', all = T)

factors$Calibration_Status[is.na(factors$Calibration_Status)] <- 'Bad'

# Export the factors sheet 
write.csv(factors,'Factor_sheet.csv', row.names = F)
# ===== Data cleaning ==============

#Generate clean data

# Remove peaks in >50% of the blanks
blk_cols <- grep("Blk", colnames(data_clean), ignore.case = TRUE)

peaks.to.remove = row.names(data_clean)[(rowSums(data_clean[,blk_cols])/ ncol(data_clean[,blk_cols])) > 0.50]

# Removing peaks
data_clean = data[-which(row.names(data) %in% peaks.to.remove),]
mol_clean = mol[-which(row.names(mol) %in% peaks.to.remove),]

# Confirming order
if(!identical(row.names(data_clean), row.names(mol_clean))){
  stop("Your conservative dataset has row name issues.")
}


# ==== Export Data =====

# Remove process blanks from the data

write.csv(data_clean,paste0(input,'Processed_CoreMS_XML_No_poorly-calibrated_Clean_Data.csv'))

write.csv(mol_clean,paste0(input,'Processed_CoreMS_XML_No_poorly-calibrated_Clean_Mol.csv'))

# ===== Exporting data per reps ======
# Deciding how many samples to keep based on a minimum number of resps that passed QAQC 
df2 = factors %>% dplyr::select(Samples, Location, Calibration_Status)
df2$Site = str_extract(df2$Samples, "EC_[A-Z0-9]+_SIR-(D|W)")
df2$Result = df2$Calibration_Status
df2$Result = gsub('Good','Pass',df2$Result)
df2$Result = gsub('Bad','Fail',df2$Result)
df2 = df2[!is.na(df2$Site), ]
df2 = df2[!is.na(df2$Result), ]
df2 = df2 %>% dplyr::select(-'Calibration_Status')

# Keeping at least 4 cal points
# For min of 4 reps
result_4_reps <- df2 %>%
  group_by(Location, Site) %>%
  dplyr::summarise(pass_count = sum(Result == "Pass"))%>%
  group_by(Location) %>%
  mutate(Keep = ifelse(any(pass_count >= 4), all(pass_count >= 4), FALSE))

# For min of 3 reps
result_3_reps <- df2 %>%
  group_by(Location, Site) %>%
  dplyr::summarise(pass_count = sum(Result == "Pass"))%>%
  group_by(Location) %>%
  mutate(Keep = ifelse(any(pass_count >= 3), all(pass_count >= 3), FALSE))

# Pull sample names to extract
samples_to_extract_4_reps <- result_4_reps %>%
  filter(Keep == TRUE) %>%
  pull(Location) # Extracting Location values where Keep is TRUE

filtered_df_4_reps <- df2 %>%
  filter(grepl(paste(samples_to_extract_4_reps, collapse = "|"), Location))

# Pull sample names to extract
samples_to_extract_3_reps <- result_3_reps %>%
  filter(Keep == TRUE) %>%
  pull(Location) # Extracting Location values where Keep is TRUE

filtered_df_3_reps <- df2 %>%
  filter(grepl(paste(samples_to_extract_3_reps, collapse = "|"), Location))

# Now filter the data_clean to only include the samples that met QAQC  threshold and also the minimum number of reps

data_4_reps <- data_clean[, intersect(names(data_clean), filtered_df_4_reps$Samples)]

data_3_reps <- data_clean[, intersect(names(data_clean), filtered_df_3_reps$Samples)]

write.csv(data_4_reps,'Processed_clean_CoreMS_XML_Int_4_reps_1p5ppm_cal_Data.csv')

write.csv(data_3_reps,'Processed_clean_CoreMS_XML_Int_3_reps_1p5ppm_cal_Data.csv')

write.csv(mol_clean,'Processed_clean_CoreMS_XML_Int_1p5ppm_cal_pts_Mol.csv')
# Rows are the masses
###########################################################################
# ====== Exploratory checks =====
# ===== Some quick checks ======
number_of_peaks <- sapply(data, function(column) sum(column != 0))
total_peaks_data <- data.frame(Samples = names(number_of_peaks), Total_Peaks = number_of_peaks)

number_of_peaks_clean <- sapply(data_clean, function(column) sum(column != 0))
total_peaks_data_clean <- data.frame(Samples = names(number_of_peaks_clean), Total_Peaks_Clean = number_of_peaks_clean)

df = merge(total_peaks_data,total_peaks_data_clean, by = 'Samples', all = T)
df$Difference = df$Total_Peaks- df$Total_Peaks_Clean

df2 = merge(df,factors, by = 'Samples', all = T)
names(df2)[2] = 'Total_Peaks_CoreMS_XML'
names(df2)[3] = 'Total_Peaks_Clean_CoreMS_XML'
names(df2)[4] = 'Difference_CoreMS_XML'
#write.csv(df2,'CoreMS_XML_Peaks.csv',row.names = F)
# ===== Load Total Peaks formularity =====
peaks = read.csv(paste0(github,'Formularity_Peaks.csv'))
core_peaks = read.csv(paste0(github,'CoreMS_Raw_Peaks.csv'))
test = merge(peaks,df2, by = c('Samples','Sample_Type'))
test2 = merge(core_peaks, df2, by = c('Samples','Sample_Type'))
library(viridis)
library(ggpmisc)
library(segmented)

lm_model <- lm(Total_Peaks_with_MF_Form ~ Total_Peaks_CoreMS_XML, data = test)

# Get the p-value
#p_value <- summary(lm_model)$coefficients[2, 4]
breaks <- seq(min(test$final.y), max(test$final.y), by = 30)

ggplot(test, aes(y = Total_Peaks_with_MF_Form, x = Total_Peaks_CoreMS_XML)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +  
  labs(y = "Total Peaks with MF Formularity", x = "Total Peaks with MF CoreMS XML") +
  theme_bw() +
  geom_smooth(method = "lm", se = FALSE, aes(group = 1), formula = y ~ x) +  # Adding linear regression line
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~~")), 
               label.x = "left", label.y = 60, size = 4)+
  annotate("text", x = 2000, y = 8000, 
           label = paste("p-value< 0.001"), vjust = 1, hjust = 1, color = "black")


ggplot(test, aes(y = Total_Peaks_with_MF_Clean_Form, x = Total_Peaks_Clean_CoreMS_XML)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +  
  labs(y = "Total Clean Peaks with MF Formularity", x = "Total Clean Peaks with MF CoreMS XML") +
  theme_bw() +
  geom_smooth(method = "lm", se = FALSE, aes(group = 1), formula = y ~ x) +  # Adding linear regression line
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~~")), 
               label.x = "left", label.y = 60, size = 4)+
  annotate("text", x = 2000, y = 8000, 
           label = paste("p-value< 0.001"), vjust = 1, hjust = 1, color = "black")


hist(test$final.x)

hist(test$final.x,
     breaks = 30,                        # Set number of bins
     col = "#4F81BD",                    # Use a nice color
     border = "white",                   # Set border color
     main = "Calibration Points in Formularity",  # Title of the plot
     xlab = "Calibration Points",        # Label for x-axis
     ylab = "Frequency",                 # Label for y-axis
     cex.main = 1.2,                     # Adjust title size
     cex.lab = 1.1,                      # Adjust axis label size
     cex.axis = 0.9)                     # Adjust 


hist(test$final.y,
     breaks = 30,                        # Set number of bins
     col = "#4F81BD",                    # Use a nice color
     border = "white",                   # Set border color
     main = "Calibration Points in CoreMS with XMLs",  # Title of the plot
     xlab = "Calibration Points",        # Label for x-axis
     ylab = "Frequency",                 # Label for y-axis
     cex.main = 1.2,                     # Adjust title size
     cex.lab = 1.1,                      # Adjust axis label size
     cex.axis = 0.9)                     # Adjust 


lm_model <- lm(Total_Peaks_CoreMS_Raw ~ Total_Peaks_CoreMS_XML, data = test2)

# Get the p-value
p_value <- summary(lm_model)$coefficients[2, 4]

ggplot(test2, aes(y = Total_Peaks_CoreMS_Raw, x = Total_Peaks_CoreMS_XML)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +  
  labs(y = "Total Peaks with MF CoreMS Raw", x = "Total Peaks with MF CoreMS XML") +
  theme_bw() +
  geom_smooth(method = "lm", se = FALSE, aes(group = 1), formula = y ~ x) +  # Adding linear regression line
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~~")), 
               label.x = "left", label.y = 60, size = 4)+
  annotate("text", x = 2000, y = 8000, 
           label = paste("p-value < 0.001 "), vjust = 1, hjust = 1, color = "black")


ggplot(test, aes(y = Total_Peaks_with_MF_Clean_Form, x = Total_Peaks_Clean_CoreMS_XML)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +  
  labs(y = "Total Clean Peaks with MF Formularity", x = "Total Clean Peaks with MF CoreMS XML") +
  theme_bw() +
  geom_smooth(method = "lm", se = FALSE, aes(group = 1), formula = y ~ x) +  # Adding linear regression line
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~~")), 
               label.x = "left", label.y = 60, size = 4)+
  annotate("text", x = 2000, y = 8000, 
           label = paste("p-value< 0.001"), vjust = 1, hjust = 1, color = "black")



