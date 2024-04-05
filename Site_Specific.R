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
# ====== Add in respiration, effect size and moisture====

# Read in data
effect.size = read.csv(paste0(github,'ECA_Effect_Size_ReadyForBoye_2023-11-08.csv'))

effect.size$Effect_Size = gsub('-9999',NA,effect.size$Effect_Size)
effect.size$Effect_Size = as.numeric(effect.size$Effect_Size)
effect.size$Location = str_extract(effect.size$Sample_Name, "EC_0[0-9]{2}|EC_([A-Za-z0-9]+)")
effect.size$Moments = NA
# Assign values based on the condition to the new column
effect.size$Moments[effect.size$Effect_Size < 70] <- "Neutral"
effect.size$Moments[effect.size$Effect_Size >= 70] <- "Cold"

effect = effect.size %>% dplyr::select(Location,Moments)

# Load rates
rates = read.csv(paste0(github,'ECA_Sediment_Incubations_mg_kg_rates_laan208_on_2023-12-01.csv'))
rates$Samples = rates$Sample_Name
rates$Samples = gsub('INC', 'SIR', rates$Samples)

moisture = read.csv(paste0(github,'ECA_Drying_Masses_Summary_merged_by_laan208_on_2023-12-01.csv'))
moisture$Samples = moisture$Sample_Name
moisture$Samples = gsub('INC', 'SIR', moisture$Samples)

# Merging and keeping important variables
m1 = merge(rates,moisture, by = 'Samples', all = TRUE)
m1 = m1 %>% dplyr::select('Samples','Respiration_Rate_mg_DO_per_kg_per_H','Moisture')

m1$Respiration_Rate_mg_DO_per_kg_per_H = gsub('-9999',NA,m1$Respiration_Rate_mg_DO_per_kg_per_H)

# Removing NA from respiration
m1 = na.omit(m1[!is.na(m1$Respiration_Rate_mg_DO_per_kg_per_H), ])

m1$log_rate = log10(abs(as.numeric(m1$Respiration_Rate_mg_DO_per_kg_per_H))+1)


# ========= Data set-up ======
# Creating factors sheet
factors = data.frame(Samples = colnames(data), Location = colnames(data), Treatment = colnames(data))
factors$Location = str_extract(factors$Location, "EC_0[0-9]{2}|EC_([A-Za-z0-9]+)")
factors$Treatment = str_extract(factors$Treatment, "W|D|Blk")

# Need to remove lambdas with no biological meaning (i.e. lambda > 0.3). Also if lambda is negative it will be set to zero. This also removes NA in lambda. Alls per communication with Hyun

# Select the mol variables of interest
# mol2 = mol %>% dplyr::select(AI_Mod,DBE_1,NOSC,delGcoxPerCmol,lamO2,delGd)%>% filter(lamO2 < 0.3) %>%
#   mutate(lamO2 = ifelse(lamO2 < 0, 0, lamO2)) %>%
#   filter(!is.na(lamO2))

mol2 = mol %>% dplyr::select(delGcoxPerCmol,lamO2,delGd)%>% filter(lamO2 < 0.3) %>%
  mutate(lamO2 = ifelse(lamO2 < 0, 0, lamO2)) %>%
  filter(!is.na(lamO2))

# ==== Calculating avg metrics and thermodynamics per sample ====
df.merge = merge(data,mol2, by = 'row.names', all = TRUE)

# df.stats = as.data.frame(matrix(NA, nrow = ncol(data), ncol = 7))
# colnames(df.stats) = c('Samples','AI_mod','NOSC','DBE','Gibbs_per_C','Gibbs_per_compound','Lambda')

df.stats = as.data.frame(matrix(NA, nrow = ncol(data), ncol = 4))
colnames(df.stats) = c('Samples','Gibbs_per_C','Gibbs_per_compound','Lambda')
for (i in 2:(ncol(data)+1)){
  df.stats$Samples[i-1] = colnames(df.merge[i])
  df.stats$Gibbs_per_C[i-1] = median(na.omit(df.merge$delGcoxPerCmol[which(df.merge[, i] > 0)]))
  df.stats$Lambda[i-1] = median(na.omit(df.merge$lamO2[which(df.merge[, i] > 0)]))
  # df.stats$AI_mod[i-1] = median(na.omit(df.merge$AI_Mod[which(df.merge[, i] > 0)]))
  # df.stats$NOSC[i-1] = median(na.omit(df.merge$NOSC[which(df.merge[, i] > 0)]))
  # df.stats$DBE[i-1] = median(na.omit(df.merge$DBE_1[which(df.merge[, i] > 0)]))
  df.stats$Gibbs_per_compound[i-1] = median(na.omit(df.merge$delGd[which(df.merge[, i] > 0)]))
}

df <- merge(df.stats, factors, by = "Samples", all = TRUE)

df = na.omit(df[!is.na(df$Gibbs_per_C), ])
df.all = merge(m1,df, by = 'Samples')

df_melt = melt(df.all)

my_colors <- c("W" = "lightblue", "D" = "darkorange")
# === Boxplot W/D per site ====
# Create the boxplot with free scales. All the variables in one
ggplot(df_melt, aes(x = Location, y = value, fill = Treatment)) +
  geom_boxplot() +
  facet_grid(variable ~ ., scales = "free", switch = "y") +
  scale_fill_manual(values = my_colors) +
  labs(x = "Location", y = "Value", fill = "Treatment") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Separating Variables into chemistry and thermodynamics 
# Subset the data into two groups based on the variable
df_melt_subset1 <- df_melt[df_melt$variable %in% c("AI_mod", "DBE", "NOSC"), ]
df_melt_subset2 <- df_melt[df_melt$variable %in% c("Gibbs_per_C", "Gibbs_per_compound", "Lambda"), ]
df_melt_subset3 <- df_melt[df_melt$variable %in% c("Moisture", "log_rate"), ]

# Create plot for variables 1 to 3
ggplot(df_melt_subset1, aes(x = Location, y = value, fill = Treatment)) +
  geom_boxplot() +
  facet_grid(variable ~., scales = "free", switch = "y") +
  scale_fill_manual(values = my_colors) +
  labs(x = "Location", y = "Value", fill = "Treatment") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('W_vs_D_Properties_per_site.pdf',width = 20, height = 15)
plot3 <- ggplot(df_melt_subset3, aes(x = Location, y = value, fill = Treatment)) +
  geom_boxplot() +
  facet_grid(variable ~., scales = "free", switch = "y") +
  scale_fill_manual(values = my_colors) +
  labs(x = "Location", y = "Value", fill = "Treatment") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Create plot for variables 4 to 6
ggplot(df_melt_subset2, aes(x = Location, y = value, fill = Treatment)) +
  geom_boxplot() +
  facet_grid(variable ~., scales = "free", switch = "y") +
  scale_fill_manual(values = my_colors) +
  labs(x = "Location", y = "Value", fill = "Treatment") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('W_vs_D_Thermodynamics_per_site.pdf',width = 20, height = 15)


# Stats - Wilcox two sided
library(dplyr)
# Filter the dataframe for Treatment W and D separately
df_W <- df_melt_subset2 %>% filter(Treatment == "W")
df_D <- df_melt_subset2 %>% filter(Treatment == "D")

# Get unique locations and variables
unique_locations <- unique(df_melt_subset2$Location)
unique_variables <- c("Gibbs_per_C", "Gibbs_per_compound", "Lambda")

# Initialize an empty list to store test results
wilcox_results <- list()

# Loop through unique locations and variables to perform Wilcoxon tests
for (loc in unique_locations) {
  for (var in unique_variables) {
    # Subset data for the current location and variable
    data_loc_var_W <- df_W %>% filter(Location == loc, variable == var)
    data_loc_var_D <- df_D %>% filter(Location == loc, variable == var)
    
    # Perform Wilcoxon rank sum test for the current location and variable
    test_result <- wilcox.test(data_loc_var_W$value, data_loc_var_D$value, alternative = "two.sided", exact = F)
    
    # Store the results in a list
    result <- tibble(
      Location = loc,
      variable = var,
      p_value = test_result$p.value
    )
    
    # Append the result to the list
    wilcox_results[[length(wilcox_results) + 1]] <- result
  }
}

# Combine results into a single dataframe
wilcox_results_df <- do.call(rbind, wilcox_results)
write.csv(wilcox_results_df,'Wilcox_results_W_vs_D_thermodynamics.csv',row.names = F)
wilcox_results_df$p_value = as.numeric(wilcox_results_df$p_value)
library(ggsignif)
# Filter out non-numeric or missing p-values
wilcox_results_df <- wilcox_results_df %>%
  filter(!is.na(p_value) & is.numeric(p_value))

# saving plot 2 to add significance by hand


# === Median per location and treatment and ratio ===
# Calculate median per Location and Treatment and then calculate a ratio
df2 = df[,-1]
medians <- df2 %>%
  group_by(Location, Treatment) %>%
  summarise(across(where(is.numeric), median, na.rm = TRUE))

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

# Define custom breaks for the x-axis
# Define custom breaks for the x-axis
custom_breaks <- seq(0.8, 1.5, by = 0.05)

# Define custom binwidths for each variable
custom_binwidths <- c(AI_mod = 0.02, NOSC = 0.04, DBE = 0.02,
                      Gibbs_per_C = 0.005, Gibbs_per_compound = 0.02,
                      Lambda = 0.02)  # Assign the binwidths as needed

# Create separate plots with variable binwidths for each variable
plots_list <- lapply(unique(ratios_melted$variable), function(var) {
  ggplot(subset(ratios_melted, variable == var), aes(x = value)) +
    geom_histogram(binwidth = custom_binwidths[var], fill = "grey", color = "black") +
    ggtitle(paste(var)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "red", size = 1) +
    scale_x_continuous(breaks = custom_breaks) +
    labs(x = "Ratio of Medians (W / D)", y = "Frequency") +
    theme_bw()
  #+
   # theme(strip.text = element_text(size = 8)) +
   # coord_cartesian(expand = FALSE)
})

# Arrange plots in a grid
do.call(gridExtra::grid.arrange, c(plots_list, ncol = 3, nrow = 2))

custom_breaks <- seq(0.8, 1.5, by = 0.05)


#################################
var = unique(ratios_melted$variable)[3]
ggplot(subset(ratios_melted, variable == var), aes(x = value)) +
  geom_histogram(binwidth = 0.012, fill = "grey", color = "black") +
  ggtitle(paste(var)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red", size = 1) +
  scale_x_continuous(breaks = seq(0.8, 1.5, by = 0.02)) +
  labs(x = "Ratio of Medians (W / D)", y = "Frequency") +
  theme_bw()
ggsave('lambda_hist.pdf',width = 20, height = 15)
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

# Arrange plots in a grid
do.call(gridExtra::grid.arrange, c(plots_list, ncol = 3))


library(ggplot2)
library(gridExtra)
columns_to_plot <- names(df.all)[5:10]

# Create a PDF file to save the plots
pdf("scatterplots_by_location_2.pdf")

# Loop through each unique location and create a page with plots arranged in a grid
for (loc in unique(df.all$Location)) {
  df_location <- subset(df.all, Location == loc)
  
  # Create a list to store plots for columns 5:10
  plots_list <- lapply(columns_to_plot, function(col) {
    ggplot(df_location, aes_string(x = col, y = "log_rate", color = "Moisture")) +
      geom_point() +
      labs(x = col, y = "log_rate", color = "Moisture") +
      ggtitle(paste("Location:", loc, "-", col)) +
      theme_bw()
  })
  
  # Arrange plots in a grid for columns 5:10
  plot_grid = do.call(grid.arrange, c(plots_list, ncol = 2))
  
  # Print the grid to the PDF
  print(plot_grid)
}

# Close the PDF device
dev.off()

# Pick sites which ratios are higher than 1 for lambda
Location_higher_one = ratios$Location[which(ratios$Lambda_Ratio>1)]

df_higher = df.all[which(df.all$Location %in% Location_higher_one),]

df_higher$Log_lambda = log10(df_higher$Lambda)

columns_to_plot <- names(df_higher)[5:10]
# Create a vector of colors for the color scale with a specific breakpoint


# Create a PDF file to save the plots
pdf("scatterplots_columns_5_to_10_higher_one_2.pdf")

# Create a list to store plots for columns 5:10
plots_list <- lapply(columns_to_plot, function(col) {
  ggplot(df_higher, aes_string(x = col, y = "log_rate", color = "Moisture")) +
    geom_point() +
    labs(x = col, y = "log_rate", color = "Moisture") +
    ggtitle(paste(col)) +
    scale_color_gradient2(midpoint = 40, low = "orange", mid = "grey", high = "lightblue") +
    theme_bw()
})

# Arrange plots in a grid for columns 5:10
plot_grid = do.call(grid.arrange, c(plots_list, ncol = 2))

# Print the grid to the PDF
print(plot_grid)

# Close the PDF device
dev.off()

Location_lower_one = ratios$Location[which(ratios$Lambda_Ratio<1)]

df_lower = df.all[which(df.all$Location %in% Location_lower_one),]

df_lower$Log_lambda = log10(df_lower$Lambda)

columns_to_plot <- names(df_lower)[5:10]
# Create a vector of colors for the color scale with a specific breakpoint


# Create a PDF file to save the plots
pdf("scatterplots_columns_5_to_10_lower_one_2.pdf")

# Create a list to store plots for columns 5:10
plots_list <- lapply(columns_to_plot, function(col) {
  ggplot(df_lower, aes_string(x = col, y = "log_rate", color = "Moisture")) +
    geom_point() +
    labs(x = col, y = "log_rate", color = "Moisture") +
    ggtitle(paste(col)) +
    scale_color_gradient2(midpoint = 100, low = "orange", mid = "grey", high = "lightblue") +
    theme_bw()
})

# Arrange plots in a grid for columns 5:10
plot_grid = do.call(grid.arrange, c(plots_list, ncol = 2))

# Print the grid to the PDF
print(plot_grid)

# Close the PDF device
dev.off()

# Looking at locations again only at thermodynamics
# 
# library(dplyr)
# library(ggplot2)
# library(ggpubr)
# 
# 
# # Create a PDF file to save the plots
# pdf("scatterplots_by_location_thermo_2.pdf")
# 
# # Loop through each unique location and create a page with plots arranged in a grid
# plot_list <- list()
# 
# for (loc in unique(df.all$Location)) {
#   df_location <- subset(df.all, Location == loc)
#   moisture_threshold <- median(df_location$Moisture)  
#   df_high_moisture <- df_location %>%
#     filter(Moisture > moisture_threshold)
#   
#   df_low_moisture <- df_location %>%
#     filter(Moisture <= moisture_threshold)
#   
#   for (i in 8:10){
#    p1 = ggplot(df_location, aes_string(x = names(df_location)[i], y = "log_rate", color = "Moisture")) +
#       geom_point() +
#       labs(x = names(df_location)[i], y = "log_rate", color = "Moisture") +
#       ggtitle(paste(loc, "-", names(df_location)[i])) +
#       geom_smooth(method = "lm", se = FALSE) +  # Add linear regression line
#      geom_smooth(data = df_high_moisture, aes_string(x = col, y = "log_rate"), method = "lm", se = FALSE) +  # Add linear regression line for high moisture
#      geom_smooth(data = df_low_moisture, aes_string(x = col, y = "log_rate"), method = "lm", se = FALSE) +   
#      stat_poly_eq(data = df_high_moisture, aes_string(label = paste(..eq.label.., sep = "~~~~")), formula = y ~ x, parse = TRUE) +  # Add R-squared and p-value for high moisture
#      stat_poly_eq(data = df_low_moisture, aes_string(label = paste(..eq.label.., sep = "~~~~")), formula = y ~ x, parse = TRUE) +  # Add R-squared and p-value for low moisture
#      theme_bw() +
#      theme_bw() +
#      theme(plot.margin = margin(1, 1, 1.5, 1.5, "cm"))
#       stat_cor(label.x = median(df_location[,i]), label.y = 0.2, aes(label = paste(..rr.label.., ..p.label.., sep = " - "))) +  # Add R-squared and p-value as text
#       theme_bw()
#    plot_list[[i-7]] <- p1
#   }
# }
#   
#   # Arrange plots in a grid for columns 5:10
#   plot_grid = do.call(grid.arrange, c(plots_list, ncol = 3))
#   
#   # Print the grid to the PDF
#   print(plot_grid)
# 
# 
# # Close the PDF device
# dev.off()
# 
# library(dplyr)
# library(ggpmisc)
# library(gridExtra)
# library(ggpubr)
columns_to_plot <- names(df.all)[8:10]

# Create a PDF file to save the plots
pdf("scatterplots_by_location_thermo_high_low_moisture_with_stats_2.pdf")

# Loop through each unique location and create a page with plots arranged in a grid
for (loc in unique(df.all$Location)) {
  df_location <- subset(df.all, Location == loc)
  
  # Define a threshold for moisture (you can set your threshold here)
  moisture_threshold <- median(df_location$Moisture)  # For example, using median
  
  # Separate data into high and low moisture groups
  df_high_moisture <- df_location %>%
    filter(Moisture > moisture_threshold)
  
  df_low_moisture <- df_location %>%
    filter(Moisture <= moisture_threshold)
  
  # Create a list to store plots for columns 8 to 10
  plots_list <- lapply(columns_to_plot, function(col) {
    p <- ggplot(df_location, aes_string(x = col, y = "log_rate", color = "Moisture")) +
      geom_point(data = df_high_moisture, aes_string(x = col, y = "log_rate"), alpha = 1) +
      geom_point(data = df_low_moisture, aes_string(x = col, y = "log_rate"), alpha = 1) +
      labs(x = col, y = "log_rate", color = "Moisture") +
      ggtitle(paste(loc, "-", col)) +
      geom_smooth(data = df_high_moisture, aes_string(x = col, y = "log_rate"), method = "lm", se = FALSE, formula = y ~ x) +  # Add linear regression line for high moisture
      geom_smooth(data = df_low_moisture, aes_string(x = col, y = "log_rate"), method = "lm", se = FALSE, formula = y ~ x) +   # Add linear regression line for low moisture
     # stat_poly_eq(data = df_high_moisture, aes(label = paste(..eq.label.., sep = "~~~~")), formula = y ~ x, parse = TRUE, geom = "text", size = 3, color = "black") +  # Add R-squared and p-value for high moisture
    #  stat_poly_eq(data = df_low_moisture, aes(label = paste(..eq.label.., sep = "~~~~")), formula = y ~ x, parse = TRUE, geom = "text", size = 3, color = "black") +  # Add R-squared and p-value for low moisture
      scale_color_gradient2(midpoint = moisture_threshold, low = "orange", mid = "grey", high = "lightblue")+
      theme_bw() 
    #+ theme(aspect.ratio=1)
      #theme(plot.margin = margin(1, 1, 1.5, 1.5, "cm"))
    
    return(p)
  })
  
  # Arrange plots in a grid for columns 8 to 10
  plot_grid = do.call(grid.arrange, c(plots_list, ncol = 1, nrow = 3))
  
  # Print the grid to the PDF
  print(plot_grid)
}

# Close the PDF device
dev.off()

df.all$Treatment = str_extract(df.all$Samples, "W|D")

columns_to_plot <- names(df.all)[8:10]

# Create a PDF file to save the plots
pdf("scatterplots_by_location_thermo_high_low_treatment_with_stats.pdf")

# Loop through each unique location and create a page with plots arranged in a grid
for (loc in unique(df.all$Location)) {
  df_location <- subset(df.all, Location == loc)
  
  # Separate data into different Treatment groups
  df_treatment_W <- df_location %>%
    filter(Treatment == "W")
  
  df_treatment_D <- df_location %>%
    filter(Treatment == "D")
  
  # Create a list to store plots for columns 8 to 10
  plots_list <- lapply(columns_to_plot, function(col) {
    p <- ggplot(df_location, aes_string(x = col, y = "log_rate", color = "Treatment")) +
      geom_point(data = df_treatment_W, aes_string(x = col, y = "log_rate"), alpha = 1, color = 'lightblue') +
      geom_point(data = df_treatment_D, aes_string(x = col, y = "log_rate"), alpha = 1, color = 'orange') +
      labs(x = col, y = "log_rate", color = "Treatment") +
      ggtitle(paste(loc, "-", col)) +
      geom_smooth(data = df_treatment_W, aes_string(x = col, y = "log_rate"), method = "lm", se = FALSE, formula = y ~ x, color = 'lightblue') +  # Add linear regression line for Treatment W
      geom_smooth(data = df_treatment_D, aes_string(x = col, y = "log_rate"), method = "lm", se = FALSE, formula = y ~ x, color = 'orange') +   # Add linear regression line for Treatment D
      #stat_poly_eq(data = df_treatment_W, aes(label = paste(..eq.label.., sep = "~~~~")), formula = y ~ x, parse = TRUE, geom = "text", size = 3, color = "black") +  # Add R-squared and p-value for Treatment W
      #stat_poly_eq(data = df_treatment_D, aes(label = paste(..eq.label.., sep = "~~~~")), formula = y ~ x, parse = TRUE, geom = "text", size = 3, color = "black") +  # Add R-squared and p-value for Treatment D
      theme_bw() 
    
    return(p)
  })
  
  # Arrange plots in a grid for columns 8 to 10
  plot_grid = do.call(grid.arrange, c(plots_list, ncol = 1, nrow = 3))
  
  # Print the grid to the PDF
  print(plot_grid)
}

# Close the PDF device
dev.off()

# === Pulling regression coefficients for the plots above =====
# By treatment

# Initialize an empty dataframe to store regression results
regression_results <- data.frame()

# Loop through each unique combination of Location, Treatment, and predictor columns
for (loc in unique(df.all$Location)) {
  for (treatment in unique(df.all$Treatment)) {
    for (col in names(df.all)[8:10]) {
      df_subset <- df.all %>%
        filter(Location == loc, Treatment == treatment)
      
      # Fit linear regression model for the current column
      lm_model <- lm(paste("log_rate ~", col), data = df_subset)
      
      # Extract coefficients, R-squared, and p-value
      coefficients <- coef(lm_model)
      rsquared <- summary(lm_model)$r.squared
      p_value <- summary(lm_model)$coefficients[, 4][2]  # P-value for the predictor variable
      
      # Store results in the dataframe
      result <- data.frame(
        Location = loc,
        Treatment = treatment,
        Column = col,
        Slope = coefficients[2],  # Extracting the slope for the predictor variable
        R_squared = rsquared,
        P_value = p_value
      )
      
      # Append to the main dataframe
      regression_results <- bind_rows(regression_results, result)
    }
  }
}

# Export the dataframe to a CSV file
write.csv(regression_results, "regression_results_per_column.csv", row.names = FALSE)

# Initialize an empty dataframe to store correlation results
correlation_results <- data.frame()

# Loop through each unique combination of Location, Treatment, and predictor columns
for (loc in unique(df.all$Location)) {
  for (treatment in unique(df.all$Treatment)) {
    for (col in names(df.all)[8:10]) {
      df_subset <- df.all %>%
        filter(Location == loc, Treatment == treatment)
      
      # Calculate Pearson correlation and its p-value for the current column
      correlation_test <- cor.test(df_subset$log_rate, df_subset[[col]], method = "pearson")
      
      # Extract correlation and p-value
      correlation <- correlation_test$estimate
      p_value <- correlation_test$p.value
      
      # Store results in the dataframe
      result <- data.frame(
        Location = loc,
        Treatment = treatment,
        Column = col,
        Correlation = correlation,
        P_value = p_value
      )
      
      # Append to the main dataframe
      correlation_results <- bind_rows(correlation_results, result)
    }
  }
}

# Export the dataframe to a CSV file
write.csv(correlation_results, "correlation_results_per_column_with_pvalues.csv", row.names = FALSE)
library(corrplot)

# === Picking now only a couple of sites for the regression plots ====

columns_to_plot <- names(df.all)[8:10]

# Create a PDF file to save the plots
pdf("scatterplots_by_selected_location_thermo_high_low_treatment.pdf")

# Loop through each unique location and create a page with plots arranged in a grid
selected_loc = c('EC_095','EC_064','EC_094','EC_037','EC_009')
for (loc in selected_loc) {
  df_location <- subset(df.all, Location == loc)
  
  # Separate data into different Treatment groups
  df_treatment_W <- df_location %>%
    filter(Treatment == "W")
  
  df_treatment_D <- df_location %>%
    filter(Treatment == "D")
  
  # Create a list to store plots for columns 8 to 10
  plots_list <- lapply(columns_to_plot, function(col) {
    p <- ggplot(df_location, aes_string(x = col, y = "log_rate", color = "Treatment")) +
      geom_point(data = df_treatment_W, aes_string(x = col, y = "log_rate"), alpha = 1, color = 'lightblue', size = 2) +
      geom_point(data = df_treatment_D, aes_string(x = col, y = "log_rate"), alpha = 1, color = 'orange', size = 2) +
      labs(x = col, y = "log_rate", color = "Treatment") +
      ggtitle(paste(loc)) +
      geom_smooth(data = df_treatment_W, aes_string(x = col, y = "log_rate"), method = "lm", se = FALSE, formula = y ~ x, color = 'lightblue') +  # Add linear regression line for Treatment W
      geom_smooth(data = df_treatment_D, aes_string(x = col, y = "log_rate"), method = "lm", se = FALSE, formula = y ~ x, color = 'orange') +   # Add linear regression line for Treatment D
      #stat_poly_eq(data = df_treatment_W, aes(label = paste(..eq.label.., sep = "~~~~")), formula = y ~ x, parse = TRUE, geom = "text", size = 3, color = "black") +  # Add R-squared and p-value for Treatment W
      #stat_poly_eq(data = df_treatment_D, aes(label = paste(..eq.label.., sep = "~~~~")), formula = y ~ x, parse = TRUE, geom = "text", size = 3, color = "black") +  # Add R-squared and p-value for Treatment D
      theme_bw() +
      theme(aspect.ratio=1)
      
    
    return(p)
  })
  
  # Arrange plots in a grid for columns 8 to 10
  plot_grid = do.call(grid.arrange, c(plots_list, ncol = 3, nrow = 2))
  
  # Print the grid to the PDF
  print(plot_grid)
}

# Close the PDF device
dev.off()
