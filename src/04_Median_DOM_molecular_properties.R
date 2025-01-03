# ==== Loading libraries =========
rm(list=ls(all=T))
library(stringr); library(devtools);  library("plyr")
library("readr");  library(readxl);library(crayon); library(vegan)
# Load in necessary libraries first
library(reshape2)
library(ggpubr) # For to combine plots
library(dplyr) # For reorganization
library(stringr) # For string manipulation
# ==== Defining paths and working directories ======
github = 'C:/Users/gara009/OneDrive - PNNL/Documents/GitHub/ECA_DOM_Thermodynamics/'
data_path = paste0(github,'Data/')
figure_path = paste0(github,'Figures/')
# ====== Read in data ======
# Processed ICR Data
data = read.csv(list.files(path = data_path, pattern = "*unique_formulas_Data.csv", full.names = T),row.names = 1)
mol = read.csv(list.files(path = data_path, pattern = "*cal_pts_Mol.csv", full.names = T), row.names = 1)
# Fixing colnames 
colnames(data) = gsub('SIR.','SIR-',colnames(data))

# ========= Calculate weighted Average ======
# Select the mol variables of interest
mol2 = mol %>% dplyr::select(AI_Mod,DBE_1,NOSC,delGcoxPerCmol,delGd,lamO2)

# Calculating weighed avg metris and thermodynamics per sample
df.merge = merge(data,mol2, by = 'row.names')

df.stats = as.data.frame(matrix(NA, nrow = ncol(data), ncol = 7))
colnames(df.stats) = c('Sample_Name','Weighted_Avg_AI_mod','Weighted_Avg_NOSC','Weighted_Avg_DBE','Weighted_Avg_delGcoxPerCmol','Weighted_Avg_delGcoxPerCompmol','Weighted_Avg_Lambda')

for (i in 2:(ncol(data)+1)){
  df.stats$Sample_Name[i-1] = colnames(df.merge[i])
  df.stats$Weighted_Avg_delGcoxPerCmol[i-1] = weighted.mean(df.merge$delGcoxPerCmol[which(df.merge[, i] > 0)],df.merge[,i][which(df.merge[, i] > 0)])
  df.stats$Weighted_Avg_delGcoxPerCompmol[i-1] = weighted.mean(df.merge$delGd[which(df.merge[, i] > 0)],df.merge[,i][which(df.merge[, i] > 0)])
  df.stats$Weighted_Avg_Lambda[i-1] = weighted.mean(df.merge$lamO2[which(df.merge[, i] > 0)],df.merge[,i][which(df.merge[, i] > 0)])
  df.stats$Weighted_Avg_AI_mod[i-1] = weighted.mean(df.merge$AI_Mod[which(df.merge[, i] > 0)],df.merge[,i][which(df.merge[, i] > 0)])
  df.stats$Weighted_Avg_NOSC[i-1] = weighted.mean(df.merge$NOSC[which(df.merge[, i] > 0)],df.merge[,i][which(df.merge[, i] > 0)])
  df.stats$Weighted_Avg_DBE[i-1] = weighted.mean(df.merge$DBE_1[which(df.merge[, i] > 0)],df.merge[,i][which(df.merge[, i] > 0)])
}

# Adding factors manually to df.stats to not complicate merging
df.stats$site = str_extract(df.stats$Sample_Name, "EC_0[0-9]{2}|EC_([A-Za-z0-9]+)")
df.stats$Treatment = str_extract(df.stats$Sample_Name, "W|D|Blk")

df.stats$Treatment = gsub('W','Wet', df.stats$Treatment)
df.stats$Treatment = gsub('D','Dry', df.stats$Treatment)

write.csv(df.stats,'Data/Weighted_averages_for_molecular_properties_per_sample_unique_formulas.csv', row.names = F)

# ===== Stats ====
wilcox_test <- function(df, variable, group_var) {
  results <- data.frame()
  unique_sites <- unique(df[[group_var]])
  for (site in unique_sites) {
    site_data <- df[df[[group_var]] == site, ]
    wet_data <- site_data[site_data$Treatment == 'Wet', variable]
    dry_data <- site_data[site_data$Treatment == 'Dry', variable]
    
    if (length(wet_data) > 0 & length(dry_data) > 0) {
      if (length(wet_data) > 1 & length(dry_data) > 1) {
        test_result <- wilcox.test(wet_data, dry_data, exact = FALSE)
        results <- rbind(results, data.frame(Site = site, Variable = variable, 
                                             W = test_result$statistic, 
                                             p.value = test_result$p.value))
      } else {
        results <- rbind(results, data.frame(Site = site, Variable = variable, 
                                             W = NA, 
                                             p.value = NA))
      }
    } else {
      results <- rbind(results, data.frame(Site = site, Variable = variable, 
                                           W = NA, 
                                           p.value = NA))
    }
  }
  return(results)
}

variables_to_test <- c('Weighted_Avg_AI_mod', 'Weighted_Avg_NOSC', 'Weighted_Avg_DBE', 'Weighted_Avg_delGcoxPerCmol','Weighted_Avg_delGcoxPerCompmol', 'Weighted_Avg_Lambda')
wilcox_results <- data.frame()

# Perform Wilcoxon tests for each variable and combine the results
for (variable in variables_to_test) {
  wilcox_results <- rbind(wilcox_results, wilcox_test(df.stats, variable, 'site'))
}

#Filter for significant p-values (e.g., p < 0.05)
significant_results <- wilcox_results %>% filter(p.value < 0.05)

# Count how many sites have significant differences across multiple mol values
significant_count <- significant_results %>%
  group_by(Site) %>%
  summarise(Num_Significant_Variables = n())

df_filtered <- df.stats %>% filter(site %in% significant_count$Site)



# ===== Calculate a median of these metrics per site and treatment ====
df.medians = df.stats %>%
  group_by(site,Treatment) %>%
  summarise(
    Median_Weighted_Avg_AI_mod = median(Weighted_Avg_AI_mod, na.rm = TRUE),
    Median_Weighted_Avg_NOSC = median(Weighted_Avg_NOSC, na.rm = TRUE),
    Median_Weighted_Avg_DBE = median(Weighted_Avg_DBE, na.rm = TRUE),
    Median_Weighted_Avg_delGcoxPerCmol = median(Weighted_Avg_delGcoxPerCmol, na.rm = TRUE),
    Median_Weighted_Avg_delGcoxPerCompmol = median(Weighted_Avg_delGcoxPerCmol, na.rm = TRUE),
    Median_Weighted_Avg_Lambda = median(Weighted_Avg_Lambda, na.rm = TRUE)
  )

write.csv(df.medians, 'Data/Medians_of_Weighted_averages_for_molecular_properties_per_site_and_treatment_unique_formulas.csv', row.names = F)
