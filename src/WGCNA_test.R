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
install.packages("WGCNA")
install.packages("BiocManager")
BiocManager::install("GO.db")
library(GO.db)
library('WGCNA')


# ==== Defining paths and working directories ======
github = 'C:/Users/vaneg/OneDrive/Documents/GitHub/ECA_DOM_Thermodynamics/'
out_plots = paste0(github,'CoreMS/Plots/')
out_data = paste0(github,'CoreMS/Data/')
setwd(github)
# ====== Read in and clean up data ======
# Processed ICR Data
data = read.csv(list.files(path = 'data/', pattern = "Processed_EC_three_reps", full.names = T), row.names = 1)
#%>% rename_with(function(x) gsub(".corems", "", x))
mol = read.csv(list.files(path = 'data/', pattern = "Processed_EC_clean", full.names = T), row.names = 1)

# mol = read.csv("FTICR_Data/GROWsurf-Processed_Mol.csv",
#                check.names = F)
# 
# # filter out poorly calibrated samples
# cal = read.csv("FTICR_Data/GROWsurf_CoreMS_Calibration_Results.csv",
#                check.names = F) %>%
#   filter(!`Cal. Points` < 3) %>%
#   filter(!`Cal. RMS Error (ppm)` > 1.5)
# data = data[,which(colnames(data) %in% cal$Sample)]
# 
# rm(cal)

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
rm(singletons, data, site.count)

# ====== WGCNA =====
# Transpose the data for WGCNA
dat = data
dat$Mass = row.names(data)
dat <- as.data.frame(lapply(dat, function(x) {
  if (is.integer(x)) {
    as.numeric(x)  # Convert integers to numeric
  } else {
    x  # Leave other columns unchanged
  }
}))
row.names(dat) = dat$Mass
dat = dat %>% dplyr::select(-'Mass')
input <- t(dat)


#Q for Bob
#Data normalization? also should I use intensities?


# WGCNA options
thresh <- 0.15
sign <- "signed"
export <- TRUE

allowWGCNAThreads()

# Thresholding
sft <- pickSoftThreshold(input, powerVector = c(1:20, seq(22, 30, by = 2)),
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
block <- blockwiseModules(input, power = p, networkType = sign,
                          TOMType = sign, minModuleSize = 20,
                          numericLabels = TRUE, corType = "pearson",
                          maxBlockSize = 12000)

# Export TOM if needed
if (export) {
  TOM <- TOMsimilarity(adjacency(input, power = p, type = sign), TOMType = sign)
  dimnames(TOM) <- list(colnames(input), colnames(input))
  
  # Export network for Cytoscape
  cyt <- exportNetworkToCytoscape(TOM,
                                  edgeFile = paste0("GROWsurf Information/GROWsurf_Edges_", thresh, "-thresh_", sign, ".txt"),
                                  nodeFile = paste0("GROWsurf Information/GROWsurf_Nodes_", thresh, "-thresh_", sign, ".txt"),
                                  threshold = thresh)
  
  # Create a module list
  colors <- data.frame(user_genome = names(block$colors),
                       Modules = paste0("ME", block$colors),
                       colors = labels2colors(block$colors))
  
  # Add alternative colors
  mod.cols <- data.frame(Modules = paste0("ME", unique(block$colors)[order(unique(block$colors))]), 
                         AltColors = c("gray", RColorBrewer::brewer.pal("Reds", n = 5)[-1], 
                                       RColorBrewer::brewer.pal("Blues", n = 5)[-1], 
                                       RColorBrewer::brewer.pal("Greens", n = 5)[-1],
                                       RColorBrewer::brewer.pal("Purples", n = 5)[-1]))
  
  colors <- colors %>% left_join(mod.cols, by = "Modules")
  
  # Write module list without duplicates
  write.table(colors, paste0("GROWsurf Information/GROWsurf_Modules_", sign, ".txt"), 
              sep = "\t", row.names = FALSE, quote = FALSE)
}

# Pulling out modules
MEs <- block$MEs
AEs <- moduleEigengenes(input, block$colors)$averageExpr

# Identify strongest module membership
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
