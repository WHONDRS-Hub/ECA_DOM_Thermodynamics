rm(list=ls(all=T))

# ===== Load necessary libraries ======
library(dplyr)
library(tidyr)
library(stringr) # For string manipulation

# ==== Working directories =====
input_path = 'C:/Users/gara009/OneDrive - PNNL/Documents - Core Richland and Sequim Lab-Field Team/Data Generation and Files/ECA/FTICR/03_ProcessedData/CoreMS/EC_Data_Processed_FTICR/Processed_with_XML/'
github = 'C:/Users/gara009/OneDrive - PNNL/Documents/GitHub/ECA_DOM_Thermodynamics/'

# ====== Read in and clean up data ======
# Processed ICR Data
data = read.csv(list.files(path = github, pattern = "*Int_4_reps_1p5ppm_cal_Data", full.names = T), row.names = 1)

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

number_of_MF <- sapply(data, function(column) sum(column != 0))
total_MF_data <- data.frame(Samples = names(number_of_MF), Total_Number_of_MF = number_of_MF)

# Extract site identifiers from column names
sample_info <- data.frame(
  Samples = colnames(data),
  site = str_extract(colnames(data), "EC_[A-Z0-9]+"),
  treatment = case_when(grepl("W", colnames(data)) ~ "Wet",
                        grepl("D", colnames(data)) ~ "Dry",TRUE~"Samples"))


df = merge(total_MF_data,sample_info, by = 'Samples')


error_data <- df %>%
  group_by(site, treatment) %>%
  summarise(
    median = median(Total_Number_of_MF),
    lower_quartile = quantile(Total_Number_of_MF, 0.25),  # 25th percentile (Q1)
    upper_quartile = quantile(Total_Number_of_MF, 0.75),  # 75th percentile (Q3)
    .groups = "drop"
  )

# Step 2: Create the bar plot with error bars using IQR
ggplot(error_data, aes(x = site, y = median, fill = treatment)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), color = "black") +  # Bar plot
  geom_errorbar(
    aes(ymin = lower_quartile, ymax = upper_quartile),  # Error bars based on IQR (Q1 and Q3)
    position = position_dodge(width = 0.9),  # Match dodge width of bars
    width = 0.2  # Error bar width
  ) +
  labs(
    x = "Site",
    y = "Median Total Molecular Formula Count\nwithin a site and treatment",
    fill = "Treatment"
  ) +
  scale_fill_manual(values = c("Wet" = "lightblue", "Dry" = "darkorange")) +  # Custom colors
  theme_bw() +  # Minimal theme
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")  # Rotate x-axis labels for readability

# make a bar plot with error info as well as ordering based on the treatment that is higher
median_data <- df  %>%
  group_by(site, treatment) %>%
  summarise(median = median(Total_Number_of_MF), .groups = "drop")

avg_data <- df  %>%
  group_by(site, treatment) %>%
  summarise(mean = mean(Total_Number_of_MF), .groups = "drop")
# Step 2: Calculate mean Wet and Dry values per site
median_values <- median_data %>%
  spread(key = treatment, value = median) %>%
  mutate(order_group = if_else(Wet > Dry, 1, 2)) %>%  # Group sites by Wet > Dry or Dry >= Wet
  arrange(order_group, desc(Wet), desc(Dry))  # Order by group, then by descending Wet or Dry mean


# Step 3: Reorder the site factor based on the new order
site_order <- median_values %>%
  pull(site)  

# Step 3: Reorder the site factor based on the new order
error_data <- error_data %>%
  mutate(site = factor(site, levels = site_order))  # Apply the new order to the site factor

# Step 4: Create the bar plot with the reordered sites
ggplot(error_data, aes(x = site, y = median, fill = treatment)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), color = "black") +  # Bar plot
  geom_errorbar(
    aes(ymin = lower_quartile, ymax = upper_quartile),  # Error bars based on IQR (Q1 and Q3)
    position = position_dodge(width = 0.9),  # Match dodge width of bars
    width = 0.2  # Error bar width
  ) +
  labs(
    x = "Site",
    y = "Median Total Molecular Formula Count\nwithin a site and treatment",
    fill = "Treatment"
  ) +
  scale_fill_manual(values = c("Wet" = "lightblue", "Dry" = "darkorange")) +  # Custom colors
  theme_bw() +  # Minimal theme
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")  # Rotate x-axis labels for readability


median_values = median_values %>%
  group_by(site)%>%
  mutate(effect_size = Wet - Dry,
         ratio = Wet/Dry)


write.csv(median_values,'Median_Site_Total_MF_order.csv',row.names = F)
