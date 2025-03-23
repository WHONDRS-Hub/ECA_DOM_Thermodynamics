# ==== Loading libraries =========
rm(list=ls(all=T))

library(stringr); library(devtools);  library("plyr")
library("readr"); library(tidyr); library(readxl);library(crayon); library(vegan)
library(reshape2)
library(ggpubr) # For to combine plots
library(dplyr) # For reorganization
library(stringr) # For string manipulation
# ==== Defining paths and working directories ======
github = 'C:/Users/gara009/OneDrive - PNNL/Documents/GitHub/ECA_DOM_Thermodynamics/'
data_path = paste0(github,'Data/')
figure_path = paste0(github,'Figures/')
# ====== Read in and clean up data ======
# Processed ICR Data
data = read.csv('Data/FTICR_combined_unique_formulas_Data.csv', row.names = 1, check.names = F)
data$Mass = rownames(data) # doing it like this to keep as much of the decimals as possible

mol = read.csv('Data/Processed_clean_CoreMS_XML_Int_1p5ppm_cal_pts_Mol.csv', row.names = 1)
mol$Mass = rownames(mol)


# Fixing colnames
colnames(data) = gsub('SIR.','SIR-',colnames(data))

# store site sample count
site.count = table(gsub("_ICR.*", "", colnames(data)))


# ========= Plots ======
# Adding mass
mol$Mass = as.numeric(as.character(row.names(mol)))

# Creating factors sheet
factors = data.frame(Samples = colnames(data), Location = colnames(data), Treatment = colnames(data))
factors$Location = str_extract(factors$Location, "EC_0[0-9]{2}|EC_([A-Za-z0-9]+)")
factors$Treatment = str_extract(factors$Treatment, "W|D|Blk")

sample_info <- data.frame(
  sample = colnames(data),
  site = str_extract(colnames(data), "EC_[A-Z0-9]+"),
  treatment = case_when(grepl("W", colnames(data)) ~ "Wet",
                        grepl("D", colnames(data)) ~ "Dry",TRUE~"Samples"))


# === Function prepare the data for Wet and Dry treatments plots ===
process_treatment_data <- function(data, mol, treatment, property_cols) {
  
  # Subset data based on Wet or Dry treatment
  treatment_data <- data %>% select(matches(treatment))
  treatment_mol <- mol %>% select(all_of(property_cols))
  
  # Merge treatment data with molecular properties
  df_merge <- merge(treatment_data,treatment_mol, by = 'row.names', all = FALSE)
  
  # Clean data by removing rows with all zeros (or NAs in some cases)
  df_cleaned <- df_merge %>%
    pivot_longer(cols = starts_with("EC"), names_to = "Sample", values_to = "Intensity") %>%
    filter(Intensity > 0) %>%
    pivot_longer(cols = c(property_cols), names_to = "Property", values_to = "Value") %>%
    mutate(Treatment = treatment)
  
  return(df_cleaned)
}

# === Process both Wet and Dry treatments molecular and thermodynamic properties ===
om_properties <- c("delGcoxPerCmol", "lamO2")

# Process Wet and Dry treatment data for OM properties
wet_om_data <- process_treatment_data(data, mol, "W", om_properties)
dry_om_data <- process_treatment_data(data, mol, "D", om_properties)

# Combine Wet and Dry data
om_combined_data <- bind_rows(wet_om_data, dry_om_data)

final_combined_data <- om_combined_data # for code compatibility

#fixing names
final_combined_data$Property = gsub('DBE_1','DBE', final_combined_data$Property)
final_combined_data$Property = gsub('delGd','delGcoxPerCompmol', final_combined_data$Property)


final_combined_data$Treatment = gsub('W','Wet',final_combined_data$Treatment)
final_combined_data$Treatment = gsub('D','Dry',final_combined_data$Treatment)

metabolite_count_data <- data %>%
  dplyr::summarise(dplyr::across(everything(), ~sum(. != 0))) %>%  # Summing non-zero values
  pivot_longer(cols = everything(), names_to = "Sample", values_to = "MetaboliteCount")

# Merge metabolite count with treatment information (Wet vs Dry)
# Merge metabolite count with treatment information (Wet vs Dry)
metabolite_count_data <- metabolite_count_data %>%
  left_join(sample_info, by = c("Sample" = "sample"))

# === Combine Metabolite Count with other data for final plot ===
final_combined_data_with_metabolite <- bind_rows(final_combined_data, 
                                                 data.frame(Property = "Metabolite Count",
                                                            Value = metabolite_count_data$MetaboliteCount,
                                                            Sample = metabolite_count_data$Sample,
                                                            Treatment = metabolite_count_data$treatment))

# === Plot Histograms for each Property ===
property_levels <- unique(final_combined_data$Property)
data_A <- final_combined_data %>% filter(Property == property_levels[1])
data_B <- final_combined_data %>% filter(Property == property_levels[2])

plot_A <- ggplot(data_A, aes(x = Value, fill = Treatment)) + 
  geom_histogram(alpha = 1, color = 'black') +
  scale_fill_manual(values = c("darkorange","lightblue")) +
  theme_bw() +
  labs(x = expression(paste(Delta, G[cox], ~(kJ~Cmol^{-1}))), 
       y = "Frequency", 
       title = "A") +
  theme(legend.position = "none")+ theme(aspect.ratio = 1)

plot_B <- ggplot(data_B, aes(x = Value, fill = Treatment)) + 
  geom_histogram(alpha = 1, color = 'black') +
  scale_fill_manual(values = c("darkorange","lightblue")) +
  theme_bw() +
  labs(x = expression(paste(lambda)), 
       y = " ", 
       title = "B") +
  theme(legend.position = "none")+ theme(aspect.ratio = 1)

library(gridExtra)
library(grid)

legend_plot <- ggplot(final_combined_data, aes(x = Value, fill = Treatment)) + 
  geom_histogram() +  # No need to customize much
  scale_fill_manual(values = c("darkorange","lightblue")) +
  theme_bw() +
  theme(legend.position = "bottom")  # Make sure legend is shown

# Extract the legend grob
tmp <- ggplot_gtable(ggplot_build(legend_plot))
leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
legend_grob <- tmp$grobs[[leg]]

# Arrange the plots with the legend at the bottom
histogram_plot = grid.arrange(
  arrangeGrob(plot_A, plot_B, ncol = 2),
  legend_grob,
  heights = c(10, 1)
)

# ===== Save just the histograms ====
# Save the plot as PNG
# ggsave(filename = file.path(figure_path, "Figure1_histogram_plot_unique_peaks.png"),
#        plot = histogram_plot,
#        width = 10, height = 6, dpi = 300, device = "png")
# 
# # Save the plot as PDF
# ggsave(filename = file.path(figure_path, "Figure1_histogram_plot_unique_peaks.pdf"),
#        plot = histogram_plot,
#        width = 10, height = 6, device = "pdf")


 # ==== Read in other data and make additional plots ====
data = read.csv(paste0(data_path,'ECA_Median_Summary_Unique_Peaks.csv')) %>% dplyr::select(Sample_Name = X, Median_delGcoxPerCmol = GFEperC_7.median,Median_Lambda = lambdaO2_7.median)

# Fixing sample names

data$Sample_Name = substr(data$Sample_Name, 1, 13)
data$Sample_Name = gsub('SIR','INC',data$Sample_Name)

sample_data = read_csv(paste0(github,'EC_Data_Package/Sample_Data/EC_Sediment_Gravimetric_Moisture.csv'),comment = '#', na = c('N/A', -9999)) %>%
  slice(-(1:11)) %>%  # Remove the first 11 rows
  dplyr::select(Sample_Name,'62948_Final_Gravimetric_Moisture_g_per_g') %>%
  mutate(`62948_Final_Gravimetric_Moisture_g_per_g` = as.numeric(`62948_Final_Gravimetric_Moisture_g_per_g`))

respiration_data =  read_csv(paste0(github,'EC_Data_Package/Sample_Data/EC_Sediment_SpC_pH_Temp_Respiration.csv'),comment = '#', na = c('N/A', -9999)) %>%
  slice(-(1:11)) %>%  # Remove the first 11 rows
  dplyr::select(Sample_Name,'Respiration_Rate_mg_DO_per_kg_per_H') %>%
  mutate(`Respiration_Rate_mg_DO_per_kg_per_H` = as.numeric(`Respiration_Rate_mg_DO_per_kg_per_H`))

data = merge(data,sample_data, by = 'Sample_Name')
data = merge(data,respiration_data, by = 'Sample_Name')

# ==== Set up data 
# Removing stes with issues
data <- data %>%
  dplyr::filter(!grepl("EC_023|EC_011|EC_012|EC_052|EC_053|EC_057", Sample_Name))

# ===== Add treatment column ====
data <- data %>%
  mutate(Treatment = case_when(
    str_detect(Sample_Name, "W") ~ "Inundated",
    str_detect(Sample_Name, "D") ~ "Dry",
    TRUE ~ "Unknown"
  ))

# ==== Make additional plots ===
# Density plot
moisture_density <- ggplot(data, aes(x = `62948_Final_Gravimetric_Moisture_g_per_g`, fill = Treatment)) +
  geom_density(alpha = 0.6, color = "black") +
  scale_fill_manual(values = c("Dry" = "darkorange", "Inundated" = "lightblue")) +
  theme_bw() +
  labs(
    x = "Gravimetric Moisture (g/g)",
    y = "Density", title = "C") +
  theme(legend.position = 'bottom',aspect.ratio = 1)

# Scatter plots

gibbs_scatter <- ggplot(data, aes(x = `62948_Final_Gravimetric_Moisture_g_per_g`, y = Median_delGcoxPerCmol, color = Treatment)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("Dry" = "darkorange", "Inundated" = "lightblue")) +
  theme_bw() +
  labs(
    x = "Gravimetric Moisture (g/g)",
    y = expression(paste(Delta, G[cox], " (kJ Cmol"^{-1}, ")")), title = "D") +
  theme(legend.position = 'none',aspect.ratio = 1)

lambda_scatter <- ggplot(data, aes(x = `62948_Final_Gravimetric_Moisture_g_per_g`, y = Median_Lambda, color = Treatment)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("Dry" = "darkorange", "Inundated" = "lightblue")) +
  theme_bw() +
  labs(
    x = "Gravimetric Moisture (g/g)",
    y = expression(lambda),title = "E")+
  theme(legend.position= 'none', aspect.ratio = 1)

legend_plot <- ggplot(data, aes(x = `62948_Final_Gravimetric_Moisture_g_per_g`, fill = Treatment)) +
  geom_density() +
  scale_fill_manual(values = c("Dry" = "darkorange", "Inundated" = "lightblue")) +
  theme_bw() +
  theme(
    legend.position = "right",
    legend.direction = "vertical",
    legend.box = "vertical",
  )

# Extract legend
tmp <- ggplot_gtable(ggplot_build(legend_plot))
leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
legend_grob <- tmp$grobs[[leg]]

# ===== Combine all plots ====
# old code
# library(cowplot)  # For plot_grid with labels
# 
# combined_plots <- plot_grid(
#   plot_A,
#   plot_B,
#   moisture_density,
#   gibbs_scatter,
#   lambda_scatter,
#   legend_grob,
#   nrow = 2,
#   labels = NULL
# )
# 
# # Add legend at bottom
# final_plot <- grid.arrange(
#   combined_plots,
#   nrow = 2,
#  heights = c(10, 1)
# )


library(patchwork)
final_plot =  plot_A / plot_B / moisture_density
# ==== Save the plot =====

ggsave(filename = file.path(figure_path, "Figure1_histogram_plot_unique_peaks_and_others.png"),
       plot = final_plot,
       width = 6, height = 10, dpi = 300, device = "png")

# Save the plot as PDF
ggsave(filename = file.path(figure_path, "Figure1_histogram_plot_unique_peaks_and_others.pdf"),
       plot = final_plot,
       width = 6, height = 10, device = "pdf")
