# ==== Loading libraries =========
rm(list=ls(all=T))

# organizational
library(tidyverse)
# library(ggpubr)
# library(egg)
# library(ggthemes)
# library(ggridges)

# technical
library(readxl)
library(Rfast)

# analytical
library(rstatix)
library(vegan)
library(canprot)
#install.packages("WGCNA")
#install.packages("BiocManager")
#BiocManager::install("GO.db")
#BiocManager::install("impute")
#BiocManager::install("preprocessCore")

library(GO.db)
library('WGCNA')


# ==== Defining paths and working directories ======
github = 'C:/Users/gara009/OneDrive - PNNL/Documents/GitHub/ECA_DOM_Thermodynamics/'
input_path = 'C:/Users/gara009/OneDrive - PNNL/Documents - Core Richland and Sequim Lab-Field Team/Data Generation and Files/ECA/FTICR/03_ProcessedData/CoreMS/EC_Data_Processed_FTICR/Processed_with_XML/'
out_plots = paste0(github,'CoreMS/Plots/')
out_data = paste0(github,'CoreMS/Data/')
#setwd(github)
# ====== Read in and clean up data ======
# Processed ICR Data
data = read.csv(list.files(path = github, pattern = "*Int_4_reps_1p5ppm_cal_Data", full.names = T), row.names = 1)
#%>% rename_with(function(x) gsub(".corems", "", x))
mol = read.csv(list.files(path = github, pattern = "*_cal_pts_Mol", full.names = T), row.names = 1)

# clean up missing peaks
mol = mol[-which(rowSums(data) == 0),]
data = data[-which(rowSums(data) == 0),]

# removing singletons (formulas found only in one site)
singletons = apply(data, 1, function(x) length(which(x > 0))) # identify
data = data[-which(singletons == 1),]
mol = mol[-which(singletons == 1),]

# store site sample count
site.count = table(gsub("_ICR.*", "", colnames(data)))

# clean up
rm(singletons, site.count)

# ====== WGCNA =====
# Transpose the data for WGCNA
dat = data

# dat$Mass = row.names(data)
# # Note that this is only an error because we have presence absence. Intensities wouldn't have an issue
# dat <- as.data.frame(lapply(dat, function(x) {
#   if (is.integer(x)) {
#     as.numeric(x)  # Convert integers to numeric
#   } else {
#     x  # Leave other columns unchanged
#   }
# }))
# row.names(dat) = dat$Mass
# dat = dat %>% dplyr::select(-'Mass')

transposed_data <- t(dat)

# ==== Normalize data with log 2 ====
constant <- 1
normalized_data <- log2(transposed_data + constant)


# WGCNA options
thresh <- 0.15
sign <- "signed"
export <- TRUE

allowWGCNAThreads()

# ==== Calculate soft thresholding ======
# Thresholding
sft <- pickSoftThreshold(normalized_data, powerVector = c(1:20, seq(22, 30, by = 2)),
                         verbose = 5, RsquaredCut = 0.8,
                         networkType = sign)

# Plot results 
par(mfrow = c(1,2));
cex1 = 0.9;

plot(sft$fitIndices[, 1],
     sft$fitIndices[, 2],
     #-sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     main = paste("Scale independence")
)

abline(h = 0.90, col = "red")
plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     main = paste("Mean connectivity")
)

#


p <- sft$powerEstimate
# Note higher threshold value helps get more meaningful relationships

#p <- 10  # Override if needed

# Automatic module detection
block <- blockwiseModules(normalized_data, power = p, networkType = sign,
                          TOMType = sign, minModuleSize = 20,
                          numericLabels = TRUE, corType = "pearson",
                          maxBlockSize = 26000)
# Increase minModuleSize to 100, do this before doing any stats, do this. The goal is to try to get the number of modules to 100-75 something higher than that might be hard to capture dynamics. 
# Try ax blocksize 24k. Chops data to make things more manageable but some peaks that might belong to a module get trapped to another block. You never get a complete module. Module merging step correlate all modules to themselves. If the maxblock is larger than your data it ensures that everything stays within ta block
# ==== Visualize the modules ====

# Convert labels to colors for plotting
mergedColors = labels2colors(block$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(
  block$dendrograms[[1]],
  mergedColors[block$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 )

# ==== Relate the modules to the treatments ====
module_df <- data.frame(
  metabolite = names(block$colors),
  colors = labels2colors(block$colors),
  Modules = paste0("ME", block$colors)
)

write_delim(module_df,
            file = "metabolite_modules_4_reps_1p5ppm_cal_intensities_11-24.txt",
            delim = "\t")

write.csv(module_df,
            file = "metabolite_modules_4_reps_1p5_cal_intensities_11-24.csv",row.names = F)

# Count features in each module and visualize
module_summary <- module_df %>%
  count(colors, name = "Feature_Count") %>%
  arrange(desc(Feature_Count))

# Print the summary
print(module_summary)

# Visualize with a bar plot
module_summary %>%
  ggplot(aes(x = reorder(colors, -Feature_Count), y = Feature_Count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(x = "Module", y = "Number of Features", title = "Feature Count per Module") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# ==== Plot treatment relationships =====
# Get Module Eigengenes per cluster
MEs0 <- moduleEigengenes(normalized_data, mergedColors)$eigengenes

write.csv(MEs0,'Module_Eigengenes_11-24.csv')
# Reorder modules so similar modules are next to each other
MEs0 <- orderMEs(MEs0)
module_order = names(MEs0) %>% gsub("ME","", .)

# Add treatment names
MEs0$treatment = row.names(MEs0)
MEs0$treatment = if_else(grepl("W", MEs0$treatment), "Wet", "Dry")
# Edits from talking with Bob
library(rstatix)
# Functions written to be friendly with tidyverse. corr test wilcox etc but tidy

# tidy & plot data
mME = MEs0 %>%
  pivot_longer(-treatment) %>%
  mutate(
    name = gsub("ME", "", name),
    name = factor(name, levels = module_order)
  )
write.csv(mME,
          file = "Module_color_per_treatment_11-24.csv",row.names = F)

stats = mME %>% 
  group_by(name) %>%
  wilcox_test(value ~ treatment) %>%
  ungroup() %>%
  mutate(padj = p.adjust(p,method = 'fdr')) # False discovery rate, permissible but doesn't get rid of 
# Pick anything that is significant from a p-value perspective 
# significantly different clusters

top_clusters = stats %>%
              filter(p < 0.05)

# It is debatable if we need a p-value correction
# PCA with the modules that I have minus the grey (No scale) and see if there is separation between wet dry and look at the biplot to see if there is something driving the split. After fishing the main modules on the PCA find which principal component demonstrates the difference between wt and dry it might not be PC1 or PC2 it might be PC10. Find that principal component when wet and dry separate and the correlate the modules to that principal component. Then look at correlation strength between the module and the PC that showed a split and then 
# Do a simple t test for each PC between wet and dry and then pick the one with the split and then do the correlation mentioned above. 
# The most variation in the data is not related to wet and dry, the hope is that one of the other PCs captures that

# Run the modules through OPLSA to pull differences after removing the gray module. Maybe gray is some of the noise or the complexity of the data to be better suited to find the differences. OPLSA is supervised ML. This should give an ordination saying which modules are the most important and then I can pull the MFs from the modules that are relevant 

# WGCNA is partially a network analysis so we could run separate networks one for dry and one for wet and then we can look at how those look and see what they are driving each other. We can't do with this is compare a module across wet and dry but we can still say something about it. Look for MFs and see where they wind up in the treatment and it will indicate how the biogeochem is changing for example how is NOSC changing on the treatments. If we want to see that NOSC is important to OM composition, we can regenerate networks by removing specific MFs that match specific NOSC criteria and then we remove all formulas outside the bin and then regenerate the network then if we see more coordination at higher nosc and we see no coordination then those NOSC were important to structuring the whole OM

#WGCNA needs at least 15 samples per treatment, 12 total seems 


mME %>% ggplot(., aes(x=treatment, y=name, fill=value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1,1)) +
  theme(axis.text.x = element_text(angle=90)) +
  labs(title = "Module-trait Relationships", y = "Modules", fill="corr")

mME %>% filter(name == 'indianred') %>%
  ggplot(., aes(x=treatment, y=value, fill=treatment)) +
  geom_boxplot() +
  theme_bw() # to check where the modules are preferred

# I can also do a ratio, average or max of the values. Avg indian red in wet / avg indian red in dry treatments if above 1 it is a wet if below 1 it is more indicative of dry. Log 2 + 1

# there is a relatioonship between the module eigen value and where it is peaking in the treatment
# Filter for top modules
threshold <- 0.8  # Adjust this threshold as needed
top_modules <- mME %>%
  filter(abs(value) > threshold)

# Save the most correlated modules
write.csv(most_correlated_modules, file = "Top_Modules_per_treatment_11-17.csv", row.names = FALSE)

top_modules %>%
  ggplot(aes(x = reorder(name, abs(value)), y = abs(value), fill = value > 0)) +  # Use absolute value for the height
  geom_bar(stat = "identity", show.legend = FALSE) +
  scale_fill_manual(values = c("red", "blue")) +
  theme_minimal() +
  labs(x = "Module", y = "Absolute Correlation with Treatment", title = "Top 10 Most Correlated Modules") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
# ===== Modules of interest to cytoscape ====
# pick out a few modules of interest here
modules_of_interest = top_modules$name

genes_of_interest = module_df %>%
  subset(colors %in% modules_of_interest)

expr_of_interest = dat[genes_of_interest$metabolite,]

# Only recalculate TOM for modules of interest (faster, altho there's some online discussion if this will be slightly off)
TOM = TOMsimilarity(adjacency(t(expr_of_interest), power = p, type = "signed"),
                    TOMType = sign) # calculates TOM; the correlations used in

# exporting network for cytoscape

dimnames(TOM) <- list(colnames(t(expr_of_interest)), colnames(t(expr_of_interest)))
# Export network for Cytoscape
cyt <- exportNetworkToCytoscape(TOM, edgeFile = paste0("ECA_Top_Modules_Edges_", thresh, "-thresh_", sign, ".txt"),
 nodeFile = paste0("ECA_Top_Modules_Nodes_", thresh, "-thresh_", sign, ".txt"),threshold = thresh)
                            
# Create a data frame with module colors
colors <- genes_of_interest

# Generate a larger set of colors to match the number of unique modules
num_modules = length(unique(colors$Modules))
library(viridis)

# Generate 13 distinct colors using the 'viridis' color palette
alt_colors <- data.frame(Modules = unique(colors$Modules),
                         AltColors = viridis(num_modules))

colors = merge(colors,alt_colors, by = 'Modules')


# Write module list without duplicates
write.table(colors, paste0("Cytoscape_Top_Modules_", sign, ".txt"), 
            sep = "\t", row.names = FALSE, quote = FALSE)


                            
# ===== Export to Cytoscape =====
# Export TOM if needed
  TOM <- TOMsimilarity(adjacency(transposed_data, power = p, type = sign), TOMType = sign)
  dimnames(TOM) <- list(colnames(transposed_data), colnames(transposed_data))
  
  # Export network for Cytoscape
  cyt <- exportNetworkToCytoscape(TOM,
                                  edgeFile = paste0("ECA_Edges_", thresh, "-thresh_", sign, ".txt"),
                                  nodeFile = paste0("ECA_Nodes_", thresh, "-thresh_", sign, ".txt"),
                                  threshold = thresh)
  
  # Create a module list
  colors <- data.frame(user_genome = names(block$colors),
                       Modules = paste0("ME", block$colors),
                       colors = labels2colors(block$colors))
  
  # Add alternative colors
  # Generate 131 distinct colors using the rainbow function
  unique_colors <- c("gray", rainbow(130))
  
  # Create the data frame
  mod.cols <- data.frame(
    Modules = paste0("ME", unique(block$colors)[order(unique(block$colors))]), 
    AltColors = unique_colors
  )
  
  colors <- colors %>% left_join(mod.cols, by = "Modules")
  
  # Write module list without duplicates
  write.table(colors, paste0("Cytoscape_All_Modules_", sign, ".txt"), 
              sep = "\t", row.names = FALSE, quote = FALSE)

# ===== Stop running here =====
# Pulling out modules
MEs <- block$MEs
AEs <- moduleEigengenes(input, block$colors)$averageExpr # average value of the intensities

# Identify strongest module membership.If a peak is significantly related to a module
ME.strength <- NULL

for (curr.mod in colnames(MEs)) {
  temp.meta <- colors %>% filter(Modules %in% curr.mod)
  
  temp.data <- data.frame(input[which(rownames(input) %in% temp.meta$user_genome),], 
                          row.names = 1)
  
  temp.corr <- cor(t(temp.data), MEs[, which(colnames(MEs) %in% curr.mod)])
  
  temp.corr <- data.frame(user_genome = row.names(temp.corr), r_value = temp.corr)
  
  ME.strength <- rbind(ME.strength, temp.corr)
}

# Add module correlations to the metadata
tax <- colors %>% left_join(ME.strength)

# Clean up
rm(block, sft, p, curr.mod, ME.strength, colors, cyt, TOM, input, export, sign, thresh, mod.cols)
