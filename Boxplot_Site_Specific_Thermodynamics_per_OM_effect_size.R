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
# Processed ICR Data
data = read.csv(list.files(pattern = "*four_reps_Intensity_Data.csv"), row.names = 1)
mol = read.csv(list.files(pattern = "*clean_Intensity_Mol.csv"), row.names = 1)
# Fixing colnames 
colnames(data) = gsub('SIR.','SIR-',colnames(data))

# ========= Data set-up ======
# Creating factors sheet
factors = data.frame(Samples = colnames(data), Location = colnames(data), Treatment = colnames(data))
factors$Location = str_extract(factors$Location, "EC_0[0-9]{2}|EC_([A-Za-z0-9]+)")
factors$Treatment = str_extract(factors$Treatment, "W|D|Blk")

# Need to remove lambdas with no biological meaning (i.e. lambda > 0.3). Also if lambda is negative it will be set to zero. This also removes NA in lambda. Alls per communication with Hyun

# Select the mol variables of interest
mol2 = mol %>% dplyr::select(delGcoxPerCmol,lamO2,delGd)%>% filter(lamO2 < 0.3) %>%
  mutate(lamO2 = ifelse(lamO2 < 0, 0, lamO2)) %>%
  filter(!is.na(lamO2))

# ==== Calculating avg metrics and thermodynamics per sample ====
df.merge = merge(data,mol2, by = 'row.names', all = TRUE)

df.stats = as.data.frame(matrix(NA, nrow = ncol(data), ncol = 4))
colnames(df.stats) = c('Samples','Gibbs_per_C','Gibbs_per_compound','Lambda')

for (i in 2:(ncol(data)+1)){
  df.stats$Samples[i-1] = colnames(df.merge[i])
  df.stats$Gibbs_per_C[i-1] = median(na.omit(df.merge$delGcoxPerCmol[which(df.merge[, i] > 0)]))
  df.stats$Lambda[i-1] = median(na.omit(df.merge$lamO2[which(df.merge[, i] > 0)]))
  df.stats$Gibbs_per_compound[i-1] = median(na.omit(df.merge$delGd[which(df.merge[, i] > 0)]))
}

df <- merge(df.stats, factors, by = "Samples", all = TRUE)

df = na.omit(df[!is.na(df$Gibbs_per_C), ])

df_melt = melt(df)

my_colors <- c("W" = "lightblue", "D" = "darkorange")

# ==== Order the sites increasing effect size of drying over gibbs per Comp mol ====
# Calculate median per Location and Treatment and then calculate a ratio of W/D as teh effect size
df2 = df[,-1]
medians <- df2 %>%
  group_by(Location, Treatment) %>%
  summarise(across(where(is.numeric), median, na.rm = TRUE))

write.csv(medians,'Median_Thermodynamics_DOM.csv', row.names = F)
# Calculate the ratio of medians for Treatment 'W' divided by Treatment 'D' within each Location
ratios <- medians %>%
  group_by(Location) %>%
  summarise(across(where(is.numeric), function(x) {
    w_median <- x[Treatment == "W"]
    d_median <- x[Treatment == "D"]
    var_name <- cur_column()
    
    if (!is.null(w_median) && !is.null(d_median) && d_median != 0) {
      data.frame(
        Ratio = w_median / d_median
      )
    } else {
      data.frame(
        Ratio = NA
      )
    }
  })) %>%
  tidyr::unnest(cols = everything(), names_sep = "_")

ratios_melted = melt(ratios)
write.csv(ratios,'Ratios_thermodynamics.csv',row.names = F)
write.csv(ratios_melted,'Ratios_thermodynamics_melted.csv',row.names = F)

# Assign position for increasing order of Gibbs per comp
# order the data
ratios2 = ratios %>% dplyr::select('Location','Gibbs_per_compound_Ratio')
ratio_order = ratios2[order(ratios$Gibbs_per_compound_Ratio), ]
ratio_order$position <- sprintf("%02d", 1:nrow(ratio_order))

ratio_order$Name = paste0(ratio_order$position,'_',ratio_order$Location)

order = ratio_order %>% dplyr::select(Location,Name)

df_melt2 = merge(df_melt,order, by = 'Location')
df_melt2$variable = gsub('Gibbs_per_compound','01_Gibbs_per_compound',df_melt2$variable)

# === Boxplot W/D per site ====
# Create the boxplot with free scales. All the variables in one
ggplot(df_melt2, aes(x = Name, y = value, fill = Treatment)) +
  geom_boxplot() +
  facet_grid(variable ~ ., scales = "free", switch = "y") +
  scale_fill_manual(values = my_colors) +
  labs(x = " ", y = " ", fill = "Treatment") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 16)) + theme(axis.text.y = element_text(size = 16)) +  guides(fill = FALSE)
ggsave("Box_plots_effect_size.pdf",width = 800, height = 600, units = "px", dpi = 72)

# ===== Plot effect siize ratios ====
# Define custom breaks for the x-axis
# Define custom breaks for the x-axis
custom_breaks <- seq(0.8, 1.5, by = 0.05)

#################################
p1 = ggplot(subset(ratios_melted, variable == var[2]), aes(x = value)) +
  geom_histogram(binwidth = 0.012, fill = "grey", color = "black") +
  geom_text(aes(label = paste(var[2])), x = 1, y = 6, size = 10, color = 'black') +
  #ggtitle(paste(var[2])) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red", size = 1) +
  scale_x_continuous(breaks = seq(0.8, 1.5, by = 0.02)) +
  labs(x = " ", y = " ") +
  theme_bw()+
  theme(axis.text.x = element_text(size = 16)) + theme(axis.text.y = element_text(size = 16)) 

p2 = ggplot(subset(ratios_melted, variable == var[1]), aes(x = value)) +
  geom_histogram(binwidth = 0.0012, fill = "grey", color = "black") +
  geom_text(aes(label = paste(var[1])), x = 1, y = 6, size = 10, color = 'black') +
 # ggtitle(paste(var[1])) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red", size = 1) +
  scale_x_continuous(breaks = seq(0.8, 1.5, by = 0.002)) +
  labs(x = " ", y = " ") +
  theme_bw()+
  theme(axis.text.x = element_text(size = 16)) + theme(axis.text.y = element_text(size = 16)) 

p3 = ggplot(subset(ratios_melted, variable == var[3]), aes(x = value)) +
  geom_histogram(binwidth = 0.012, fill = "grey", color = "black") +
  geom_text(aes(label = paste(var[3])), x = 1, y = 6, size = 10, color = 'black') +
 # ggtitle(paste(var[3])) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red", size = 1) +
  scale_x_continuous(breaks = seq(0.8, 1.5, by = 0.02)) +
  labs(x = "Ratio of Medians (W / D)", y = " ") +
  theme_bw()+
  theme(axis.text.x = element_text(size = 16)) + theme(axis.text.y = element_text(size = 16)) 

#ggsave('lambda_hist.pdf',width = 20, height = 15)
####
# Create separate kernel density plots for each variable
plots_list <- lapply(unique(ratios_melted$variable), function(var) {
  ggplot(subset(ratios_melted, variable == var), aes(x = value)) +
    geom_density(fill = "grey", color = "black") +
    ggtitle(paste(var)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "red", size = 1) +
    scale_x_continuous(breaks = custom_breaks) +
    labs(x = "Ratio of Medians (W / D)", y = "Density") +
    theme_bw() 
})

library(gridExtra)

final_plot <- grid.arrange(p1, p2, p3, ncol = 1,nrow = 3)

