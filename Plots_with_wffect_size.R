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
data = read.csv(list.files(pattern = "*four_reps_Data.csv"), row.names = 1)
mol = read.csv(list.files(pattern = "*clean_Mol.csv"), row.names = 1)
# Fixing colnames 
colnames(data) = gsub('SIR.','SIR-',colnames(data))

effect.size = read.csv(paste0(github,'ECA_Effect_Size_ReadyForBoye_2023-11-08.csv'))

effect.size$Effect_Size = gsub('-9999',NA,effect.size$Effect_Size)
effect.size$Effect_Size = as.numeric(effect.size$Effect_Size)
effect.size$Location = str_extract(effect.size$Sample_Name, "EC_0[0-9]{2}|EC_([A-Za-z0-9]+)")
effect.size$Moments = NA
# Assign values based on the condition to the new column
effect.size$Moments[effect.size$Effect_Size < 70] <- "Neutral"
effect.size$Moments[effect.size$Effect_Size >= 70] <- "Cold"

effect = effect.size %>% dplyr::select(Location,Moments)
# ========= Plots ======
# Removing peaks that were not assigned a molecular formula
data = data[!is.na(mol$MolForm),]
mol = mol[!is.na(mol$MolForm),]

# Renaming bs1_classes
mol$bs1_class[grep(";", mol$bs1_class)] = "Multi-class"

# Adding mass
mol$Mass = as.numeric(as.character(row.names(mol)))

# Creating factors sheet
factors = data.frame(Samples = colnames(data), Location = colnames(data), Treatment = colnames(data))
factors$Location = str_extract(factors$Location, "EC_0[0-9]{2}|EC_([A-Za-z0-9]+)")
factors$Treatment = str_extract(factors$Treatment, "W|D|Blk")

factors = merge(factors,effect, by = 'Location')

# Principal component analysis
pca = prcomp(x = t(data))

# Everything below this line should not change very much whether we use PCA or NMDS
ordination.scores = scores(pca) # Works with both PCA and NMDS, change the object accordingly
ordination.scores = as.data.frame(ordination.scores) # ggplot doesn't like matrices - needs to be converted to a data frame
ordination.scores$Moments= factors$Moments # Adding in sample type to our ordination scores object

# We have everything necessary for ggplot - we want to plot PC1 and PC2
ggplot(data = ordination.scores, aes(x = PC1, y = PC2, color = Moments))+
  xlab(paste0("PC1 (", summary(pca)$importance[2,1]*100, "%)"))+
  ylab(paste0("PC2 (", summary(pca)$importance[2,2]*100, "%)"))+
  geom_point(size = 2) + theme_bw()+
  theme(text = element_text(color = "black", size = 12),
        axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(color = "black", size = 14),
        legend.text = element_text(color = "black", size = 14),
        panel.border = element_rect(color = "black"),
        axis.ticks = element_line(colour = "black"), 
        panel.grid = element_blank(),
        panel.background = element_blank())

### Beta-diversity
# Creating distance matrix
dist = vegdist(x = t(data), method = "jaccard") # Using Jaccard for historical reasons (ICR data is often analyzed using it)

# Plotting a Jaccard heatmap
dist.melt = melt(as.matrix(dist))


ggplot(data = dist.melt, aes(x = Var1, y = Var2, fill = value))+
  geom_tile() + scale_fill_gradient2(low = "gray100", mid = "gray80", high = "darkred", midpoint = 0.4)+
  xlab(NULL) + ylab(NULL)+
  vert_x_theme_2
ggsave(paste0("Jacad_Heat_map_All_sites","_",Sys.Date(),".pdf"))

# Plotting Jaccard NMDS
nms = metaMDS(dist, trymax = 1000) # Determining NMDS
nms = as.data.frame(scores(nms)) # Conveting to scores
nms = cbind(factors, nms)

library(pals)
nms %>%
  #filter(!grepl("YDE_Raw", Location)) %>% filter(!grepl("rep1", rep))%>%
  ggplot(aes(x = NMDS1, y = NMDS2))+
  #geom_point() +
  geom_point(aes(color = Moments),size = 4) 

