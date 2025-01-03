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

# Extract site identifiers from column names
sample_info <- data.frame(
  sample = colnames(data),
  site = str_extract(colnames(data), "EC_[A-Z0-9]+"),
  treatment = case_when(grepl("W", colnames(data)) ~ "Wet",
                        grepl("D", colnames(data)) ~ "Dry",TRUE~"Samples"))


# Split data by site
split_data <- split(sample_info, sample_info$site)

# Function to compare treatments within a site
compare_treatments <- function(site_data, full_data) {
  site_name <- unique(site_data$site)
  message("Processing site: ", site_name)
  
  # Get columns for wet and dry treatments
  wet_cols <- site_data$sample[site_data$treatment == "Wet"]
  dry_cols <- site_data$sample[site_data$treatment == "Dry"]
  
  # Subset the data
  wet_data <- full_data[, wet_cols, drop = FALSE]
  dry_data <- full_data[, dry_cols, drop = FALSE]
  
  # Summarize presence/absence
  wet_presence <- rowSums(wet_data > 0) > 0
  dry_presence <- rowSums(dry_data > 0) > 0
  
  # Identify differences
  unique_to_wet <- rownames(full_data)[wet_presence & !dry_presence]
  unique_to_dry <- rownames(full_data)[dry_presence & !wet_presence]
  shared <- rownames(full_data)[wet_presence & dry_presence]
  
  # Return a summary
  list(
    site = site_name,
    unique_to_wet = unique_to_wet,
    unique_to_dry = unique_to_dry,
    shared = shared
  )
}

# Apply the function across all sites
results <- lapply(split_data, compare_treatments, full_data = data)

# Combine results into a summary dataframe
summary_results <- lapply(results, function(res) {
  data.frame(
    site = res$site,
    formula = c(res$unique_to_wet, res$unique_to_dry, res$shared),
    status = c(
      rep("unique_to_wet", length(res$unique_to_wet)),
      rep("unique_to_dry", length(res$unique_to_dry)),
      rep("shared", length(res$shared))
    )
  )
}) %>%
  bind_rows()
write.csv(summary_results, "Molecular_Formula_site_treatment_comparison_results.csv", row.names = FALSE)

# ===== Plot to visualize results =====

# Load necessary libraries
library(ggplot2)

# Summarize results into counts for plotting
plot_data <- summary_results %>%
  group_by(site, status) %>%
  summarise(count = n(), .groups = "drop")

# Create the stacked bar plot
ggplot(plot_data, aes(x = site, y = count, fill = status)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    title = "Compare Molecular Formulas Between Wet and Dry Treatments",
    x = "Site",
    y = "Count of Molecular Formulas",
    fill = "Status"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), # Rotate x-axis labels
    text = element_text(size = 12)
  )

# ===== Subset the ICR data to only keep those unique formulas ====
library(dplyr)

# Function to subset the data for a specific site and treatment
subset_unique_formulas <- function(site_data, full_data, treatment_type) {
  site_name <- unique(site_data$site)
  message("Subsetting for site: ", site_name, " | Treatment: ", treatment_type)
  
  # Separate wet and dry samples
  wet_cols <- site_data$sample[site_data$treatment == "Wet"]
  dry_cols <- site_data$sample[site_data$treatment == "Dry"]
  
  # Subset the data for both treatments
  wet_data <- full_data[, wet_cols, drop = FALSE]
  dry_data <- full_data[, dry_cols, drop = FALSE]
  
  # Determine unique formulas
  wet_presence <- rowSums(wet_data > 0) > 0
  dry_presence <- rowSums(dry_data > 0) > 0
  
  if (treatment_type == "Wet") {
    unique_formulas <- wet_presence & !dry_presence
  } else if (treatment_type == "Dry") {
    unique_formulas <- dry_presence & !wet_presence
  } else {
    stop("Invalid treatment type. Use 'Wet' or 'Dry'.")
  }
  
  # Subset the original data for the unique formulas
  unique_formula_data <- full_data[unique_formulas, , drop = FALSE]
  
  # Further filter to keep only columns matching the current treatment
  treatment_cols <- site_data$sample[site_data$treatment == treatment_type]
  unique_formula_data <- unique_formula_data[, treatment_cols, drop = FALSE]
  
  return(unique_formula_data)
}

# Initialize an empty dataframe to combine all subsets (using list for efficient merging)
combined_data_list <- list()

# Loop through each site and treatment to subset the data
for (site in unique(sample_info$site)) {
  site_data <- sample_info %>% filter(site == !!site)
  
  for (treatment_type in c("Wet", "Dry")) {
    subset_data <- subset_unique_formulas(site_data, data, treatment_type)
    
    # Add row names as a column for later merging
    subset_data <- subset_data %>%
      tibble::rownames_to_column("mass")
    
    # Store the subset data in the list
    combined_data_list[[paste(site, treatment_type, sep = "_")]] <- subset_data
  }
}

# Step 2: Combine all subsets into one dataframe
# We will use `dplyr::full_join` to merge them by the "mass" column (molecular formulas)
combined_data <- combined_data_list[[1]]  # Start with the first dataset

for (i in 2:length(combined_data_list)) {
  combined_data <- full_join(combined_data, combined_data_list[[i]], by = "mass", suffix = c("", paste("_", names(combined_data_list)[i], sep = "")))
}

# Step 3: Ensure all formulas are present as rownames, replace non-existing formulas with 0
# Replace mass column back into rownames and remove it from the data
rownames(combined_data) <- combined_data$mass
combined_data$mass <- NULL  # Remove the mass column after setting it as rownames

# Step 4: Replace missing values (NA) with 0 (indicating no presence of that MF in that sample)
combined_data[is.na(combined_data)] <- 0

# Save the final dataframe
write.csv(combined_data, "combined_unique_formulas_Data.csv", row.names = TRUE, quote = FALSE)

# ===== Some checks =====
count_non_zero_rows_data1 <- combined_data %>%
  summarise(across(everything(), ~sum(. != 0)))  %>% pivot_longer(cols = everything(), names_to = "Column", values_to = "NonZeroCount") %>%
  left_join(sample_info, by = c("Column" = "sample"))

average_non_zero_counts_data1 <- count_non_zero_rows_data1 %>%
  group_by(site, treatment) %>%          # Group by site and treatment
  summarise(AverageNonZeroCount = mean(NonZeroCount), .groups = "drop")  


median_non_zero_counts_data1 <- count_non_zero_rows_data1 %>%
  group_by(site, treatment) %>%          # Group by site and treatment
  summarise(MedianNonZeroCount = median(NonZeroCount), .groups = "drop")  

library(ggplot2)

# Create a barplot to visualize the average non-zero counts per site and treatment
ggplot(average_non_zero_counts_data1, aes(x = site, y = AverageNonZeroCount, fill = treatment)) +
  geom_bar(stat = "identity", position = "dodge", color = 'black') +    # Create the bar plot with dodged bars
  labs(x = "Site", y = "Average Unique Molecular Formula Count\nwithin a site and treatment", fill = "Treatment") +  # Labels for the plot
  scale_fill_manual(values = c("Wet" = "lightblue", "Dry" = "darkorange")) +   # Custom colors for Wet and Dry
  theme_bw() +      # Use a minimal theme
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")

# Make a bar plot with error bars
error_data <- count_non_zero_rows_data1 %>%
  group_by(site, treatment) %>%
  summarise(
    mean = mean(NonZeroCount),
    sd = sd(NonZeroCount),
    se = sd / sqrt(n()),  # Standard error of the mean
    .groups = "drop"
  )

# Step 2: Create the bar plot with error bars
ggplot(error_data, aes(x = site, y = mean, fill = treatment)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), color = "black") +  # Bar plot
  geom_errorbar(
    aes(ymin = mean - se, ymax = mean + se),  # Error bars based on SE
    position = position_dodge(width = 0.9),  # Match dodge width of bars
    width = 0.2  # Error bar width
  ) +
  labs(
    x = "Site",
    y = "Average Unique Molecular Formula Count\nwithin a site and treatment",
    fill = "Treatment"
  ) +
  scale_fill_manual(values = c("Wet" = "lightblue", "Dry" = "darkorange")) +  # Custom colors
  theme_bw() +  # Minimal theme
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")  # Rotate x-axis labels for readability

# make a bar plot with error info as well as ordering based on the treatment that is higher
# Step 1: Calculate mean values for each site-treatment combination
error_data <- count_non_zero_rows_data1 %>%
  group_by(site, treatment) %>%
  summarise(
    mean = mean(NonZeroCount),
    sd = sd(NonZeroCount),
    se = sd / sqrt(n()),  # Standard error of the mean
    .groups = "drop"
  )


mean_data <- count_non_zero_rows_data1  %>%
  group_by(site, treatment) %>%
  summarise(mean = mean(NonZeroCount), .groups = "drop")

# Step 2: Calculate mean Wet and Dry values per site
mean_values <- mean_data %>%
  spread(key = treatment, value = mean) %>%
  mutate(order_group = if_else(Wet > Dry, 1, 2)) %>%  # Group sites by Wet > Dry or Dry >= Wet
  arrange(order_group, desc(Wet), desc(Dry))  # Order by group, then by descending Wet or Dry mean

write.csv(mean_values,'Site_order.csv',row.names = F)
# Step 3: Reorder the site factor based on the new order
site_order <- mean_values %>%
  pull(site)  

# Step 3: Reorder the site factor based on the new order
error_data <- error_data %>%
  mutate(site = factor(site, levels = site_order))  # Apply the new order to the site factor

# Step 4: Create the bar plot with the reordered sites
ggplot(error_data, aes(x = site, y = mean, fill = treatment)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), color = "black") +  # Bar plot
  geom_errorbar(
    aes(ymin = mean - se, ymax = mean + se),  # Error bars based on SE
    position = position_dodge(width = 0.9),  # Match dodge width of bars
    width = 0.2  # Error bar width
  ) +
  labs(
    x = "Site",
    y = "Average Unique Molecular Formula Count\nwithin a site and treatment",
    fill = "Treatment"
  ) +
  scale_fill_manual(values = c("Wet" = "lightblue", "Dry" = "darkorange")) +  # Custom colors
  theme_bw() +  # Minimal theme
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")  # Rotate x-axis labels for readability

# Create a boxplot to visualize the distribution of non-zero counts per site and treatment
# Step 1: Perform pairwise Wilcoxon tests per site
wilcox_results <- count_non_zero_rows_data1 %>%
  group_by(site) %>%
  summarise(
    p_value = wilcox.test(NonZeroCount ~ treatment)$p.value
  ) %>%
  mutate(significant = ifelse(p_value < 0.05, TRUE, FALSE))

# Step 2: Merge significance results with original data for plotting
count_non_zero_rows_data1 <- count_non_zero_rows_data1 %>%
  left_join(wilcox_results, by = "site")

# Step 3: Add stars to the plot
ggplot(count_non_zero_rows_data1, aes(x = site, y = NonZeroCount, fill = treatment)) +
  geom_boxplot(color = "black", outlier.shape = 16, outlier.colour = "darkgreen") +
  labs(x = "Site", 
       y = "Unique MF Count",
       fill = "Treatment") +
  scale_fill_manual(values = c("Wet" = "lightblue", "Dry" = "darkorange")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top") +
  # Add stars for significant comparisons
  geom_text(data = wilcox_results %>% filter(significant),
            aes(x = site, y = 0, label = "*"), 
            inherit.aes = FALSE, size = 6, vjust = 0.5, color = "red")


# ===== Repeat with median ====
# Create a barplot to visualize the median non-zero counts per site and treatment
ggplot(median_non_zero_counts_data1, aes(x = site, y = MedianNonZeroCount, fill = treatment)) +
  geom_bar(stat = "identity", position = "dodge", color = 'black') +    # Create the bar plot with dodged bars
  labs(x = "Site", y = "Median Unique Molecular Formula Count\nwithin a site and treatment", fill = "Treatment") +  # Labels for the plot
  scale_fill_manual(values = c("Wet" = "lightblue", "Dry" = "darkorange")) +   # Custom colors for Wet and Dry
  theme_bw() +      # Use a minimal theme
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")

# Make a bar plot with error bars
# Since we are using medians now we will have IQR for the error bars which span from the lower Q to the upper Q
error_data <- count_non_zero_rows_data1 %>%
  group_by(site, treatment) %>%
  summarise(
    median = median(NonZeroCount),
    lower_quartile = quantile(NonZeroCount, 0.25),  # 25th percentile (Q1)
    upper_quartile = quantile(NonZeroCount, 0.75),  # 75th percentile (Q3)
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
    y = "Median Unique Molecular Formula Count\nwithin a site and treatment",
    fill = "Treatment"
  ) +
  scale_fill_manual(values = c("Wet" = "lightblue", "Dry" = "darkorange")) +  # Custom colors
  theme_bw() +  # Minimal theme
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")  # Rotate x-axis labels for readability

# make a bar plot with error info as well as ordering based on the treatment that is higher
median_data <- count_non_zero_rows_data1  %>%
  group_by(site, treatment) %>%
  summarise(median = median(NonZeroCount), .groups = "drop")

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
    y = "Median Unique Molecular Formula Count\nwithin a site and treatment",
    fill = "Treatment"
  ) +
  scale_fill_manual(values = c("Wet" = "lightblue", "Dry" = "darkorange")) +  # Custom colors
  theme_bw() +  # Minimal theme
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")  # Rotate x-axis labels for readability

# ==== effect size of unique peaks ====

median_values = median_values %>%
  group_by(site)%>%
  mutate(effect_size = Wet - Dry,
         ratio = Wet/Dry)
  

write.csv(median_values,'Median_Site_order.csv',row.names = F)
# Create a boxplot to visualize the distribution of non-zero counts per site and treatment
# Step 1: Perform pairwise Wilcoxon tests per site
wilcox_results <- count_non_zero_rows_data1 %>%
  group_by(site) %>%
  summarise(
    p_value = wilcox.test(NonZeroCount ~ treatment)$p.value
  ) %>%
  mutate(significant = ifelse(p_value < 0.05, TRUE, FALSE))

# Step 2: Merge significance results with original data for plotting
count_non_zero_rows_data1 <- count_non_zero_rows_data1 %>%
  left_join(wilcox_results, by = "site")

# Step 3: Add stars to the plot
ggplot(count_non_zero_rows_data1, aes(x = site, y = NonZeroCount, fill = treatment)) +
  geom_boxplot(color = "black", outlier.shape = 16, outlier.colour = "darkgreen") +
  labs(x = "Site", 
       y = "Unique MF Count",
       fill = "Treatment") +
  scale_fill_manual(values = c("Wet" = "lightblue", "Dry" = "darkorange")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top") +
  # Add stars for significant comparisons
  geom_text(data = wilcox_results %>% filter(significant),
            aes(x = site, y = 0, label = "*"), 
            inherit.aes = FALSE, size = 6, vjust = 0.5, color = "red")

# Create a boxplot to visualize the distribution of non-zero counts per treatment
wilcox_test <- wilcox.test(NonZeroCount ~ treatment, data = count_non_zero_rows_data1)

# Extract the p-value
p_value <- wilcox_test$p.value

# Add p-value to the plot
ggplot(count_non_zero_rows_data1, aes(x = treatment, y = NonZeroCount, fill = treatment)) +
  geom_boxplot(color = "black", outlier.shape = 16, outlier.colour = "darkgreen") +
  labs(x = " ", 
       y = "Unique MF Count",
       fill = "Treatment") +
  scale_fill_manual(values = c("Wet" = "lightblue", "Dry" = "darkorange")) +
  theme_bw() +
  theme(legend.position = "top") +
  annotate("text", x = 1.5, y = max(count_non_zero_rows_data1$NonZeroCount, na.rm = TRUE), 
           label = paste("p =", signif(p_value, digits = 3)), 
           hjust = 0, vjust = -0.5, size = 4, fontface = "italic")  # Add p-value at top-left # Place the legend at the top
# Boxplot for avg counts
wilcox_test <- wilcox.test(AverageNonZeroCount ~ treatment, data = average_non_zero_counts_data1)

# Extract the p-value
p_value <- wilcox_test$p.value

# Add p-value to the plot
ggplot(average_non_zero_counts_data1, aes(x = treatment, y = AverageNonZeroCount, fill = treatment)) +
  geom_boxplot(color = "black", outlier.shape = 16, outlier.colour = "darkgreen") +
  labs(x = " ", 
       y = "Avg Unique MF Count",
       fill = "Treatment") +
  scale_fill_manual(values = c("Wet" = "lightblue", "Dry" = "darkorange")) +
  theme_bw() +
  theme(legend.position = "top") +
  annotate("text", x = 1.5, y = max(average_non_zero_counts_data1$AverageNonZeroCount, na.rm = TRUE), 
           label = paste("p =", signif(p_value, digits = 3)), 
           hjust = 0, vjust = -0.5, size = 4, fontface = "italic")  # Add p-value at top-left

# Original data 
count_non_zero_rows_data <- data %>%
  summarise(across(everything(), ~sum(. != 0)))  %>% pivot_longer(cols = everything(), names_to = "Column", values_to = "NonZeroCount") %>%
  left_join(sample_info, by = c("Column" = "sample"))

average_non_zero_counts_data <- count_non_zero_rows_data%>%
  group_by(site, treatment) %>%          # Group by site and treatment
  summarise(AverageNonZeroCount = mean(NonZeroCount), .groups = "drop")  

library(ggplot2)

# Create a barplot to visualize the average non-zero counts per site and treatment
ggplot(average_non_zero_counts_data, aes(x = site, y = AverageNonZeroCount, fill = treatment)) +
  geom_bar(stat = "identity", position = "dodge", color = 'black') +    # Create the bar plot with dodged bars
  labs(x = "Site", y = "Average Total Molecular Formula Count\nwithin a site and treatment", fill = "Treatment") +  # Labels for the plot
  scale_fill_manual(values = c("Wet" = "lightblue", "Dry" = "darkorange")) +   # Custom colors for Wet and Dry
  theme_bw() +      # Use a minimal theme
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")


# Create a boxplot to visualize the distribution of non-zero counts per site and treatment
ggplot(count_non_zero_rows_data, aes(x = site, y = NonZeroCount, fill = treatment)) +
  geom_boxplot(color = "black", outlier.shape = 16, outlier.colour = "darkgreen") +  # Boxplot with black border and red outliers
  labs(x = "Site", 
       y = "Total MF Count",  # Y-axis label
       fill = "Treatment") +  # Fill legend label
  scale_fill_manual(values = c("Wet" = "lightblue", "Dry" = "darkorange")) +   # Custom colors for Wet and Dry
  theme_bw() +      # Use a minimal theme
  theme(axis.text.x = element_text(angle = 45, hjust = 1),   # Rotate x-axis labels for readability
        legend.position = "top")  # Place the legend at the top

