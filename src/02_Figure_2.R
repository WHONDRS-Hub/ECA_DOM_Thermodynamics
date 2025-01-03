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
mol = read.csv(list.files(path = data_path, pattern = "*cal_pts_Mol.csv", full.names = T), row.names = 1)
data = read.csv(list.files(path = data_path, pattern = "*Int_4_reps_1p5ppm_cal_Data", full.names = T), row.names = 1)

# Fixing colnames
colnames(data) = gsub('SIR.','SIR-',colnames(data))

# clean up missing peaks
data = data[-which(rowSums(data) == 0),]

# removing singletons (formulas found only in one site)
singletons = apply(data, 1, function(x) length(which(x > 0))) # identify
data = data[-which(singletons == 1),]

# store site sample count
site.count = table(gsub("_ICR.*", "", colnames(data)))

# clean up
rm(singletons, site.count)

# === Calculate total number of formulas =====

number_of_MF <- sapply(data, function(column) sum(column != 0))
total_MF_data <- data.frame(Samples = names(number_of_MF), Total_Number_of_MF = number_of_MF)

# Extract site identifiers from column names
sample_info <- data.frame(
  Samples = colnames(data),
  site = str_extract(colnames(data), "EC_[A-Z0-9]+"),
  treatment = case_when(grepl("W", colnames(data)) ~ "Wet",
                        grepl("D", colnames(data)) ~ "Dry",TRUE~"Samples"))


df = merge(total_MF_data,sample_info, by = 'Samples')

# ==== Box plot of all formula per treatment =====
wilcox_results <- df %>%
  summarise(
    p_value = wilcox.test(Total_Number_of_MF ~ treatment)$p.value
  ) %>%
  mutate(significant = ifelse(p_value < 0.05, TRUE, FALSE))
p_value <- signif(wilcox_results$p_value, digits = 2)

boxplot = ggplot(df, aes(x = treatment, y = Total_Number_of_MF, fill = treatment)) +
  # Add boxplot layer
  geom_boxplot(color = "black", outlier.shape = NA) +
  # Add points for individual data with transparency
  geom_jitter(position = position_jitter(width = 0.2), alpha = 0.8, size = 2, color = "darkgrey") +
  # Add labels and titles
  labs(x = " ", 
       y = "Total Molecular Formula Count",
       fill = "Treatment") +
  # Set custom fill colors for treatments
  scale_fill_manual(values = c("Wet" = "lightblue", "Dry" = "darkorange")) +
  # Apply theme
  theme_bw() +
  theme(legend.position = "top") +
  # Add p-value annotation in the top-left corner
  annotate("text", x = 0.5, y = max(df$Total_Number_of_MF, na.rm = TRUE), 
           label = paste("p =", p_value), 
           hjust = 0, vjust = 1, size = 4, fontface = "italic")


# ===== Calculate median values within sites and treatments and error as IQR ===
error_data <- df %>%
  group_by(site, treatment) %>%
  summarise(
    median = median(Total_Number_of_MF),
    lower_quartile = quantile(Total_Number_of_MF, 0.25),  # 25th percentile (Q1)
    upper_quartile = quantile(Total_Number_of_MF, 0.75),  # 75th percentile (Q3)
    .groups = "drop"
  )


# Ordering plot based on the treatment that is higher
median_data <- df  %>%
  group_by(site, treatment) %>%
  summarise(median = median(Total_Number_of_MF), .groups = "drop")

#Calculate mean Wet and Dry values per site
median_values <- median_data %>%
  spread(key = treatment, value = median) %>%
  mutate(order_group = if_else(Wet > Dry, 1, 2)) %>%  # Group sites by Wet > Dry or Dry >= Wet
  arrange(order_group, desc(Wet), desc(Dry))  # Order by group, then by descending Wet or Dry mean


# Reorder the site factor based on the new order
site_order <- median_values %>%
  pull(site)  

# Reorder the site factor based on the new order
error_data <- error_data %>%
  mutate(site = factor(site, levels = site_order))  # Apply the new order to the site factor

#Create the bar plot with the reordered sites
barplot = ggplot(error_data, aes(x = site, y = median, fill = treatment)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), color = "black") +  # Bar plot
  geom_errorbar(
    aes(ymin = lower_quartile, ymax = upper_quartile),  # Error bars based on IQR (Q1 and Q3)
    position = position_dodge(width = 0.9),  # Match dodge width of bars
    width = 0.2  # Error bar width
  ) +
  labs(
    x = " ",
    y = "Median Total Molecular Formula Count\nwithin a site and treatment",
    fill = "Treatment"
  ) +
  scale_fill_manual(values = c("Wet" = "lightblue", "Dry" = "darkorange")) +  # Custom colors
  theme_bw() +  # Minimal theme
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")  # Rotate x-axis labels for readability


median_values = median_values %>%
  group_by(site)%>%
  mutate(effect_size = Wet - Dry,
         ratio = Wet/Dry)


write.csv(median_values,'Archive/Figure2_Median_Site_Total_MF_order.csv',row.names = F)

# ===== Merge plots and export ====
library(patchwork)

# Combine the two plots and add labels
combined_plot <- (boxplot + labs(tag = "A")) | (barplot + labs(tag = "B"))

combined_plot = combined_plot + plot_layout(widths = c(1, 2))

# Save as PDF
ggsave(paste0(figure_path,"Figure2_combined_plots.pdf"), plot = combined_plot, width = 12, height = 6)

# Save as PNG
ggsave(paste0(figure_path,"Figure2_combined_plots.png"), plot = combined_plot, width = 12, height = 6, dpi = 300)
