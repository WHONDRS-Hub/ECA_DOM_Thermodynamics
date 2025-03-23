rm(list=ls(all=T))

# ===== Load necessary libraries ======
library(dplyr)
library(tidyr)
library(stringr) # For string manipulation

# ==== Working directories =====
github = 'C:/Users/gara009/OneDrive - PNNL/Documents/GitHub/ECA_DOM_Thermodynamics/'
data_path = paste0(github,'Data/')
figure_path = paste0(github,'Figures/')
# ====== Read in and clean up data ======
# Processed ICR Data
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

# Extract site identifiers from column names
sample_info <- data.frame(
  sample = colnames(data),
  site = str_extract(colnames(data), "EC_[A-Z0-9]+"),
  treatment = case_when(grepl("W", colnames(data)) ~ "Wet",
                        grepl("D", colnames(data)) ~ "Dry",TRUE~"Samples"))


# Split data by site
split_data <- split(sample_info, sample_info$site)

# ===== Function to compare treatments within a site =====
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

# ===== Calculate unique and shared formulas ====
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
write.csv(summary_results, "Data/Unique_Molecular_Formulas_per_site_per_treatment.csv", row.names = FALSE)

# ===== Plot to visualize results =====

# Load necessary libraries
library(ggplot2)

# Summarize results into counts for plotting
plot_data <- summary_results %>%
  group_by(site, status) %>%
  summarise(count = n(), .groups = "drop")

# Create the stacked bar plot
barplot = ggplot(plot_data, aes(x = site, y = count, fill = status)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    title = " ",
    x = " ",
    y = "Count of Molecular Formulas",
    fill = "Molecular Formula Status"
  ) +
  theme_bw() +
  theme(legend.position = "top") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), # Rotate x-axis labels
    text = element_text(size = 12)
  )

# Save as PDF
ggsave(paste0(figure_path,"FigureS1_unique_formula_plots.pdf"), plot = barplot, width = 12, height = 6)

# Save as PNG
ggsave(paste0(figure_path,"FigureS1_unique_formula_plots.png"), plot = barplot, width = 12, height = 6, dpi = 300)
# ===== Subset the ICR data to only keep those unique formulas for each treatment ====
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

# Combine all subsets into one dataframe. Merge them by the "mass" column (molecular formulas)
combined_data <- combined_data_list[[1]]  # Start with the first dataset

for (i in 2:length(combined_data_list)) {
  combined_data <- full_join(combined_data, combined_data_list[[i]], by = "mass", suffix = c("", paste("_", names(combined_data_list)[i], sep = "")))
}

# Ensure all formulas are present as rownames, replace non-existing formulas with 0

rownames(combined_data) <- combined_data$mass
combined_data$mass <- NULL  # Remove the mass column after setting it as rownames

#Replace missing values (NA) with 0 (indicating no presence of that MF in that sample)
combined_data[is.na(combined_data)] <- 0

# Save the final dataframe
write.csv(combined_data, "Data/FTICR_combined_unique_formulas_Data.csv", row.names = TRUE, quote = FALSE)

