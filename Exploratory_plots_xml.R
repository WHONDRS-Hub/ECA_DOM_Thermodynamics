# ==== Loading libraries =========
rm(list=ls(all=T))

library(stringr); library(devtools);  library("plyr")
library("readr"); library(tidyverse); library(readxl);library(crayon); library(vegan)

# ==== Defining paths and working directories ======
setwd('C:/Users/gara009/OneDrive - PNNL/Documents - Core Richland and Sequim Lab-Field Team/Data Generation and Files/ECA/FTICR/03_ProcessedData/CoreMS/EC_Data_Processed_FTICR/Processed_with_XML/')


# ====== Read in data ======
# Processed ICR Data
data = read.csv(list.files(pattern = "*Data.csv"), row.names = 1)
mol = read.csv(list.files(pattern = "Mol.csv"), row.names = 1)
# Fixing colnames 
colnames(data) = gsub('SIR.','SIR-',colnames(data))
colnames(data) = gsub('SIR-Blk.','SIR_Blk-',colnames(data))
colnames(data) = gsub('.corems','',colnames(data))

# All Calibration file
calib = read.csv(list.files(pattern = "*Calibration_Results.csv"))
# Load Poorly calibrated file

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
cal = calib %>% dplyr::select(Samples = Sample,final = 'Cal..Points')
factors = merge(factors,cal, by = 'Samples', all = T)

number_of_peaks <- sapply(data, function(column) sum(column != 0))
total_peaks_data <- data.frame(Samples = names(number_of_peaks), Total_Peaks = number_of_peaks)

number_of_peaks_clean <- sapply(data_clean, function(column) sum(column != 0))
total_peaks_data_clean <- data.frame(Samples = names(number_of_peaks_clean), Total_Peaks_Clean = number_of_peaks_clean)

df = merge(total_peaks_data,total_peaks_data_clean, by = 'Samples', all = T)
df$Difference = df$Total_Peaks- df$Total_Peaks_Clean

df2 = merge(df,factors, by = 'Samples', all = T)
df2$Location = str_extract(df2$Samples, "EC_[A-Z0-9]+")

# ====== Calibration Threshold ======
library(ggpubr); library(viridis)

lm_model <- lm(Total_Peaks ~ final, data = df2)

num_categories <- length(unique(df2$Location))
palette <- magma(num_categories)

ggplot(df2, aes(color = factor(Location), y = Total_Peaks, x = final)) +
  geom_point() +
  scale_color_manual(values = palette) +
  geom_smooth(method = "lm", se = FALSE, color = "black", formula = y ~ x) +
  annotate("text", x = 25, y = 2000, 
           label = paste("R² =", round(summary(lm_model)$r.squared, 2), "\n", "p < 0.001"), hjust = 1, vjust = 1) +
  labs(x = "Calibration Points", y = "Total peaks") +
  theme_bw()

lm_model <- lm(final ~ Total_Peaks, data = df2)
palette <- viridis(num_categories)

ggplot(df2, aes(color = factor(Location), x = Total_Peaks, y = final)) +
  geom_point() +
  scale_color_manual(values = palette) +
  geom_smooth(method = "lm", se = FALSE, color = "black", formula = y ~ x) +
  annotate("text", y = 75, x = 500, 
           label = paste("R² =", round(summary(lm_model)$r.squared, 2), "\n", "p < 0.001"), hjust = 1, vjust = 1) +
  labs(y = "Calibration Points", x = "Total peaks") +
  theme_bw()+
  geom_hline(yintercept = c(4, 15), linetype = "dashed", color = c("red",'blue'))


df2$Site = str_extract(df2$Samples, "EC_[A-Z0-9]+_SIR-(D|W)")
df2$Result = NA
df2$Result = ifelse(df2$final < 15, "Fail", "Pass")
df2$Total = 1
df2 = df2[!is.na(df2$Difference), ]

# Convert Result to a factor for better plotting
df2$Result = factor(df2$Result)



# Create a stacked bar plot
ggplot(df2, aes(x = Site, y = Total, fill = Result)) +
  geom_bar(stat = "identity") +
  labs(x = "Site", y = "Total", fill = "Result") +
  theme_bw()+
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1))


# Create a stacked bar plot
ggplot(df2, aes(x = Location, y = Total, fill = Result)) +
  geom_bar(stat = "identity") +
  labs(x = "Site", y = "Total", fill = "Result") +
  theme_bw()+
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1))

# 4 cal points
df2$Result = NA
df2$Result = ifelse(df2$final < 4, "Fail", "Pass")
df2$Total = 1
df2 = df2[!is.na(df2$Difference), ]

# Convert Result to a factor for better plotting
df2$Result = factor(df2$Result)



# Create a stacked bar plot
ggplot(df2, aes(x = Site, y = Total, fill = Result)) +
  geom_bar(stat = "identity") +
  labs(x = "Site", y = "Total", fill = "Result") +
  theme_bw()+
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1))


# Create a stacked bar plot
ggplot(df2, aes(x = Location, y = Total, fill = Result)) +
  geom_bar(stat = "identity") +
  labs(x = "Site", y = "Total", fill = "Result") +
  theme_bw()+
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1))


# ====== Multivariate Stats =======
factors$Location = str_extract(factors$Samples, "EC_[A-Z0-9]+")
factors$Treatment = factors$Samples
factors$Treatment = str_extract(factors$Treatment, "W|D|Blk")

pca = prcomp(x = t(data)) # Calculating PCA
pca = as.data.frame(scores(pca)) # Converting to PCA scores in order to plot using ggplot
pca = cbind(factors, pca)

pca <- pca %>%
  mutate(final_category = ifelse(final > 15, "greater_than_15", "less_or_equal_15"))%>%
  mutate(other_category = ifelse(final > 4, "greater_than_4", "less_or_equal_4"))

# Define colors for the new categories
color_palette <- c("greater_than_4" ="blue", "less_or_equal_4" = "red")

pca %>%
  ggplot(aes(x = PC1, y = PC2))+
  geom_point(aes(color = other_category, shape = Treatment),size = 2) +
  # geom_point() +
 # geom_label(aes(label = row.names(pca)))+
  scale_shape_manual(values = c(0,1,2,3,4,5,6,7,8,9,10))+
  scale_color_manual(values = color_palette) +
  theme_bw()


### Beta-diversity
# Creating distance matrix
dist = vegdist(x = t(data), method = "jaccard") # Using Jaccard for historical reasons (ICR data is often analyzed using it)
library(reshape2)
# Plotting a Jaccard heatmap
dist.melt = melt(as.matrix(dist))


ggplot(data = dist.melt, aes(x = Var1, y = Var2, fill = value))+
  geom_tile() + scale_fill_gradient2(low = "gray100", mid = "gray80", high = "darkred", midpoint = 0.4)+
  xlab(NULL) + ylab(NULL)

ggsave(paste0("Jacad_Heat_map_All_sites","_",Sys.Date(),".pdf"))

# Plotting Jaccard NMDS
nms = metaMDS(dist, trymax = 1000) # Determining NMDS
nms = as.data.frame(scores(nms)) # Conveting to scores
nms = cbind(factors, nms)

library(pals)
nms %>%
  #filter(!grepl("YDE_Raw", Location)) %>% filter(!grepl("rep1", rep))%>%
  ggplot(aes(x = NMDS1, y = NMDS2))+
  #geom_point() +
  geom_point(aes(color = Treatment),size = 2) +
  theme_bw()
# scale_color_manual(values=as.vector(alphabet(26)))
ggsave(paste0("NMDS_all_sites","_",Sys.Date(),".pdf"))

# ===== Exporting data ====
# Deciding how many samples to keep based on a minimum number of resps that passed QAQC 
# Within Core MS, anything that has calibration points is technically good so exporting data with min of 4 points
df2 = df2 %>%
  mutate(Result = ifelse(is.na(final), "Fail", "Pass"))
# For min of 4 reps
result_4_reps <- df2 %>%
group_by(Location, Site) %>%
  summarise(pass_count = sum(Result == "Pass"))%>%
  group_by(Location) %>%
  mutate(Keep = ifelse(any(pass_count >= 4), all(pass_count >= 4), FALSE))


# For min of 3 reps
result_3_reps <- df2 %>%
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

write.csv(data_4_reps,'Processed_EC_four_reps_four_cal_pts_Data.csv')

write.csv(data_3_reps,'Processed_EC_three_reps_four_cal_pts_Data.csv')

write.csv(mol_clean,'Processed_EC_clean_Mol.csv')
