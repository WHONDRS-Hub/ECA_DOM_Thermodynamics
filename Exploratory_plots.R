# ==== Loading libraries =========
rm(list=ls(all=T))

library(stringr); library(devtools);  library("plyr")
library("readr"); library(tidyverse); library(readxl);library(crayon); library(vegan)

# ==== Defining paths and working directories ======
setwd('C:/Users/gara009/OneDrive - PNNL/Documents - Core Richland and Sequim Lab-Field Team/Data Generation and Files/ECA/FTICR/03_ProcessedData/EC_Data_Processed_FTICR/')


# ====== Read in data ======
# Processed ICR Data
data = read.csv(list.files(pattern = "*EC_Data.csv"), row.names = 1)
mol = read.csv(list.files(pattern = "EC_Mol.csv"), row.names = 1)
# Fixing colnames 
colnames(data) = gsub('SIR.','SIR-',colnames(data))
colnames(data) = gsub('SIR-Blk.','SIR_Blk-',colnames(data))

# All Calibration file
calib = read.csv(list.files(pattern = "*_All_Calibrations.csv"))
# Load Poorly calibrated file
poor.cal = read.csv(list.files(pattern = "*_Poorly_Calibrated_Samples.csv"))

# ======= Cleaning up data ======
# Check the blanks and remove peaks from the data that are in at least 50% of the blanks 

file.name = "EC_Clean" # 

# Fixing column names if they begin with numbers
if(length(grep("^X", colnames(data))) > 0){
  colnames(data) = gsub("^X", "", colnames(data))
} # R begins integer column names with X's - this fixes that

# Factor creation
factors = data.frame(Samples = colnames(data), 
                     Sample_Type = case_when(grepl("Blk", colnames(data)) ~ "Process Blank",
                                             grepl("W", colnames(data)) ~ "Wet",
                                             grepl("D", colnames(data)) ~ "Dry",TRUE~"Samples"),
                     Cal_Quality = "Good")

factors$Cal_Quality[which(factors$Samples %in% poor.cal$samples)] = "Poor"



# Converting data to presence/absence
int.data = data
data[data > 0] = 1


#Generate clean data

#
### Moving onto the permissive dataset
# Peaks in >50% of the blanks
peaks.to.remove = row.names(data)[(rowSums(data[,which(factors$Sample_Type %in% "Process Blank")])/
                                     ncol(data[,which(factors$Sample_Type %in% "Process Blank")])) > 0.25]

# Removing peaks
data_clean = data[-which(row.names(data) %in% peaks.to.remove),]
mol_clean = mol[-which(row.names(mol) %in% peaks.to.remove),]

# Confirming order
if(!identical(row.names(data_clean), row.names(mol_clean))){
  stop("Your conservative dataset has row name issues.")
}

# Remove process blanks from the data
data_clean = data_clean[,-grep("Blk", colnames(data_clean))]

# ===== Some quick checks ======
cal = calib %>% dplyr::select(Samples = samples,'final')
factors = merge(factors,cal, by = 'Samples', all = T)

number_of_peaks <- sapply(data, function(column) sum(column != 0))
total_peaks_data <- data.frame(Samples = names(number_of_peaks), Total_Peaks = number_of_peaks)

number_of_peaks_clean <- sapply(data_clean, function(column) sum(column != 0))
total_peaks_data_clean <- data.frame(Samples = names(number_of_peaks_clean), Total_Peaks_Clean = number_of_peaks_clean)

df = merge(total_peaks_data,total_peaks_data_clean, by = 'Samples', all = T)
df$Difference = df$Total_Peaks- df$Total_Peaks_Clean

df2 = merge(df,factors, by = 'Samples', all = T)
# ====== Calibration Threshold ======
df2$Site = str_extract(df2$Samples, "EC_[A-Z0-9]+_SIR-(D|W)")
df2$Result = NA
df2$Result = ifelse(df2$final < 15, "Fail", "Pass")
df2$Total = 1
df2 = df2[!is.na(df2$Difference), ]
df2$Location = str_extract(df2$Site, "EC_[A-Z0-9]+")

# Convert Result to a factor for better plotting
df2$Result = factor(df2$Result)

# Create a stacked bar plot
ggplot(df2, aes(x = Site, y = Total, fill = Result)) +
  geom_bar(stat = "identity") +
  labs(x = "Site", y = "Total", fill = "Result") +
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1))

# Deciding how many samples to keep based on a minimum number of resps that passed QAQC 

# For min of 4 reps
result_4_reps <- df2 %>%
  mutate(Location = str_extract(Site, "EC_[A-Z0-9]+")) %>%
group_by(Location, Site) %>%
  summarise(pass_count = sum(Result == "Pass"))%>%
  group_by(Location) %>%
  mutate(Keep = ifelse(any(pass_count >= 4), all(pass_count >= 4), FALSE))


# For min of 3 reps
result_3_reps <- df2 %>%
  mutate(Location = str_extract(Site, "EC_[A-Z0-9]+")) %>%
  group_by(Location, Site) %>%
  summarise(pass_count = sum(Result == "Pass"))%>%
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

write.csv(data_4_reps,'Processed_EC_four_reps_Data.csv')

write.csv(data_3_reps,'Processed_EC_three_reps_Data.csv')

write.csv(mol_clean,'Processed_EC_clean_Mol.csv')
