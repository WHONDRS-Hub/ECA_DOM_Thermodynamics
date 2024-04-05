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


# ====== Read in data ======
# Processed ICR Data
data = read.csv(list.files(pattern = "*four_reps_Data.csv"), row.names = 1)
mol = read.csv(list.files(pattern = "*clean_Mol.csv"), row.names = 1)
# Fixing colnames 
colnames(data) = gsub('SIR.','SIR-',colnames(data))

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
 
# === Density plots of OM properties ====
temp.data = data[which(grepl('W',colnames(data))),]
temp.mol = mol[,c("AI_Mod", "DBE_1", "NOSC")] 

df.merge = merge(temp.data,temp.mol, by = 'row.names', all = TRUE)

AI_wet = NA
NOSC_wet = NA
DBE_wet = NA
for (i in 1:178){
   temp.AI = df.merge$AI_Mod[which(df.merge[, i] > 0)]
   temp.NOSC = df.merge$NOSC[which(df.merge[, i] > 0)] 
   temp.DBE = df.merge$DBE_1[which(df.merge[, i] > 0)] 
  AI_wet = c(AI_wet,temp.AI)
  NOSC_wet = c(NOSC_wet,temp.NOSC) 
  DBE_wet = c(DBE_wet,temp.DBE) 
  
}
df.wet = as.data.frame(cbind(AI = AI_wet,NOSC = NOSC_wet,DBE = DBE_wet))
df.wet = df.wet[-1, ]
df.mol = melt(as.matrix(df.wet)) # Melting as a matrix to get the Var1/Var2 melt format
df.mol$Treatment = "W" # Adding on a qualifier to know what sample type the data is coming from

melt.mol.by.type = df.mol # Creating the object to eventually go into ggplot

temp.data = data[which(grepl('D',colnames(data))),]
temp.mol = mol[,c("AI_Mod", "DBE_1", "NOSC")] 

df.merge = merge(temp.data,temp.mol, by = 'row.names', all = TRUE)

AI_dry = NA
NOSC_dry = NA
DBE_dry = NA
for (i in 1:178){
  temp.AI = df.merge$AI_Mod[which(df.merge[, i] > 0)]
  temp.NOSC = df.merge$NOSC[which(df.merge[, i] > 0)] 
  temp.DBE = df.merge$DBE_1[which(df.merge[, i] > 0)] 
  AI_dry = c(AI_dry,temp.AI)
  NOSC_dry = c(NOSC_dry,temp.NOSC) 
  DBE_dry = c(DBE_dry,temp.DBE) 
  
}
df.dry = as.data.frame(cbind(AI = AI_dry,NOSC = NOSC_dry,DBE = DBE_dry))
df.dry = df.dry[-1, ]
df.mol = melt(as.matrix(df.dry))
df.mol$Treatment = "D" # Adding on a qualifier to know what sample type the data is coming from

melt.mol.by.type = rbind(melt.mol.by.type, df.mol)

rm(temp.mol)

# Adding in molecular formula-by-sample count into our melt.mol.by.type object
melt.mol.by.type = rbind(melt.mol.by.type, data.frame(Var1 = colnames(data), 
                                                      Var2 = "Molecular Formula Count",
                                                      value = colSums(na.omit(data)), 
                                                      Treatment = factors$Treatment))

# Statistics
metric.stats = NULL
uniq.met = unique(as.character(melt.mol.by.type$Var2))

for(i in 1:length(uniq.met)){
  w = which(melt.mol.by.type$Var2 %in% uniq.met[i])
  wil.test = wilcox.test(value~Treatment, data = melt.mol.by.type[w,], alternative = "two.sided")
  wil.test = data.frame(Comparison = uniq.met[i], W = wil.test$statistic, P.value = wil.test$p.value,
                        Test = wil.test$alternative)
  metric.stats = rbind(metric.stats, wil.test)
}

metric.stats$P.value = p.adjust(metric.stats$P.value, method = "fdr")

rm(wil.test, w, i)

# Plotting metrics
ggplot(melt.mol.by.type, aes(x = value, group = Treatment))+
  geom_density(aes(fill = Treatment), alpha = 1, color = 'black')+
  scale_fill_manual(values = c("darkorange","lightblue")) +  # Set colors for W and D
  facet_wrap(Var2~., scales = "free", ncol = 1)+
  theme_bw() + theme(text = element_text(size=18, color="black"),
                     axis.text = element_text(color = "black"),
                     axis.ticks = element_line(color = "black"),
                     panel.background = element_blank(),
                     panel.grid = element_blank())
ggsave('W_vs_D_properties_per_overall_density.pdf',width = 20, height = 15)
# === Repeat fo only thermodynamic properties ======
temp.data = data[which(grepl('W',colnames(data))),]
temp.mol = mol[,c("delGcoxPerCmol", "lamO2", "delGd")] 
temp.mol = temp.mol %>% filter(lamO2 < 0.3) %>%
  mutate(lamO2 = ifelse(lamO2 < 0, 0, lamO2)) %>%
  filter(!is.na(lamO2))
df.merge = merge(temp.data,temp.mol, by = 'row.names', all = TRUE)

GFE_wet = NA
lambda_wet = NA
dGd_wet = NA
for (i in 1:178){
  temp.GFE = df.merge$delGcoxPerCmol[which(df.merge[, i] > 0)]
  temp.lambda = df.merge$lamO2[which(df.merge[, i] > 0)] 
  temp.dGd = df.merge$delGd[which(df.merge[, i] > 0)] 
  GFE_wet = c(GFE_wet,temp.GFE)
  lambda_wet = c(lambda_wet,temp.lambda) 
  dGd_wet = c(dGd_wet,temp.dGd) 
  
}
df.wet = as.data.frame(cbind(Gibbs_per_C_mol = GFE_wet,Gibbs_per_compound_mol = dGd_wet,Lambda = lambda_wet))
df.wet = df.wet[-1, ]
df.mol = melt(as.matrix(df.wet)) # Melting as a matrix to get the Var1/Var2 melt format
df.mol$Treatment = "W" # Adding on a qualifier to know what sample type the data is coming from

melt.mol.by.type = df.mol # Creating the object to eventually go into ggplot

temp.data = data[which(grepl('D',colnames(data))),]

df.merge = merge(temp.data,temp.mol, by = 'row.names', all = TRUE)

GFE_dry = NA
lambda_dry = NA
dGd_dry = NA
for (i in 1:178){
  temp.GFE = df.merge$delGcoxPerCmol[which(df.merge[, i] > 0)]
  temp.lambda = df.merge$lamO2[which(df.merge[, i] > 0)] 
  temp.dGd = df.merge$delGd[which(df.merge[, i] > 0)] 
  GFE_dry = c(GFE_dry,temp.GFE)
  lambda_dry = c(lambda_dry,temp.lambda) 
  dGd_dry = c(dGd_dry,temp.dGd) 
  
}
df.dry = as.data.frame(cbind(Gibbs_per_C_mol = GFE_dry,Gibbs_per_compound_mol = dGd_dry,Lambda = lambda_dry))
df.dry = df.dry[-1, ]
df.mol = melt(as.matrix(df.dry)) # Melting as a matrix to get the Var1/Var2 melt format
df.mol$Treatment = "D" # Adding on a qualifier to know what sample type the data is coming from

melt.mol.by.type = rbind(melt.mol.by.type, df.mol)

rm(temp.mol)

# Statistics
metric.stats = NULL
uniq.met = unique(as.character(melt.mol.by.type$Var2))

for(i in 1:length(uniq.met)){
  w = which(melt.mol.by.type$Var2 %in% uniq.met[i])
  wil.test = wilcox.test(value~Treatment, data = melt.mol.by.type[w,], alternative = "two.sided")
  wil.test = data.frame(Comparison = uniq.met[i], W = wil.test$statistic, P.value = wil.test$p.value,
                        Test = wil.test$alternative)
  metric.stats = rbind(metric.stats, wil.test)
}

metric.stats$P.value = p.adjust(metric.stats$P.value, method = "fdr")

rm(wil.test, w, i)

# Plotting metrics
ggplot(melt.mol.by.type, aes(x = value, group = Treatment))+
  geom_histogram(aes(fill = Treatment), alpha = 1, color = 'black')+
  scale_fill_manual(values = c("darkorange","lightblue")) +  # Set colors for W and D
  facet_wrap(Var2~., scales = "free", ncol = 1)+
  theme_bw() + theme(text = element_text(size=18, color="black"),
                     axis.text = element_text(color = "black"),
                     axis.ticks = element_line(color = "black"),
                     panel.background = element_blank(),
                     panel.grid = element_blank())

test = subset(melt.mol.by.type, melt.mol.by.type$Treatment == 'D')
ggplot(test, aes(x = value, group = Treatment))+
  geom_histogram(aes(fill = Treatment), alpha = 1, color = 'black')+
  scale_fill_manual(values = c("darkorange","lightblue")) +  # Set colors for W and D
  facet_wrap(Var2~., scales = "free", ncol = 1)+
  theme_bw() + theme(text = element_text(size=18, color="black"),
                     axis.text = element_text(color = "black"),
                     axis.ticks = element_line(color = "black"),
                     panel.background = element_blank(),
                     panel.grid = element_blank())
ggsave('W_vs_D_Thermodynamics_per_overall.pdf',width = 20, height = 15)

# === Repeating only counting molecular formula once =====
# Unique sample types (i.e., water and sediment)
uniq.sample.type = unique(factors$Treatment)

# Shorthand way of doing those 6 commands above
data.by.type = as.data.frame(matrix(nrow = nrow(data), ncol = length(uniq.sample.type), 
                                    dimnames = list(row.names(data), uniq.sample.type)))

# Using a for-loop to run through our data and identify peaks present in each sample type
for(i in 1:length(uniq.sample.type)){
  # Find temporary datasets for sample types
  temp = data[,which(factors$Treatment %in% uniq.sample.type[i])]
  
  # Add in sums into the parent data frame
  data.by.type[,i] = rowSums(temp)
}

rm(temp, i)

# Resetting data back to presence/absence
data.by.type[data.by.type > 0] = 1

# Partitioning the molecular information based upon sample
temp.mol = mol[which(data.by.type$W > 0),] # Finding molecular information for wet
temp.mol = mol[,c("delGcoxPerCmol", "lamO2", "delGd")] 
temp.mol = temp.mol %>% filter(lamO2 < 0.3) %>%
  mutate(lamO2 = ifelse(lamO2 < 0, 0, lamO2)) %>%
  filter(!is.na(lamO2))
temp.mol = melt(as.matrix(temp.mol)) # Melting as a matrix to get the Var1/Var2 melt format
temp.mol$Treatment = "W" # Adding on a qualifier to know what sample type the data is coming from

melt.mol.by.type = temp.mol # Creating the object to eventually go into ggplot

temp.mol = mol[which(data.by.type$D > 0),]
temp.mol = mol[,c("delGcoxPerCmol", "lamO2", "delGd")] 
temp.mol = temp.mol %>% filter(lamO2 < 0.3) %>%
  mutate(lamO2 = ifelse(lamO2 < 0, 0, lamO2)) %>%
  filter(!is.na(lamO2))
temp.mol = melt(as.matrix(temp.mol)) # Melting as a matrix to get the Var1/Var2 melt format
temp.mol$Treatment = "D" # Adding on a qualifier to know what sample type the data is coming from

melt.mol.by.type = rbind(melt.mol.by.type, temp.mol)
# ################################################ #
#### Comparing el. comp. and bs1 between groups ####
# ################################################ #

### Elemental Composition by sample
el.comp = matrix(data = 0, nrow = ncol(data), ncol = length(unique(mol$El_comp)), 
                 dimnames = list(colnames(data), unique(mol$El_comp)))

for(i in 1:nrow(el.comp)){
  temp = mol[which(data[,i] > 0),] # Mol data for a given sample
  
  for(j in 1:ncol(el.comp)){
    el.comp[i,j] = length(which(temp$El_comp %in% colnames(el.comp)[j]))
  }
} # Counting the number of times a given elemental composition appears in a dataset

el.comp = as.data.frame(t(apply(el.comp, 1, function(x) (x/sum(x))*100))) # Relative abundance
el.comp = cbind(factors, el.comp)
el.comp = melt(el.comp, id.vars = colnames(factors))

### Compound Class by sample
comp.class = matrix(data = 0, nrow = ncol(data), ncol = length(unique(mol$bs1_class)), 
                    dimnames = list(colnames(data), unique(mol$bs1_class)))

for(i in 1:nrow(comp.class)){
  temp = mol[which(data[,i] > 0),] # Mol data for a given sample
  
  for(j in 1:ncol(comp.class)){
    comp.class[i,j] = length(which(temp$bs1_class %in% colnames(comp.class)[j]))
  }
} # Counting the number of times a given elemental composition appears in a dataset

rm(temp, i, j)

comp.class = as.data.frame(t(apply(comp.class, 1, function(x) (x/sum(x))*100)))
comp.class = cbind(factors, comp.class)
comp.class = melt(comp.class, id.vars = colnames(factors))
comp.class$variable = gsub("Hydrocarbon", "HC", comp.class$variable) # Shortening names
comp.class$variable = gsub("Carbohydrate", "Carb.", comp.class$variable)

# Performing stats on el. comp
el.stats = NULL

for(curr.el in unique(el.comp$variable)){
  temp = el.comp[which(el.comp$variable %in% curr.el),]
  max.val = max(temp$value)
  temp = wilcox.test(value~Treatment, data = temp, alternative = "two.sided")
  temp = data.frame(Comparison = curr.el, W = temp$statistic, p.value = temp$p.value, max.val = max.val)
  el.stats = rbind(el.stats, temp)
}

el.stats$p.value = p.adjust(el.stats$p.value, method = "fdr")
el.stats$Symbol = symnum(el.stats$p.value, corr = FALSE, na = FALSE, 
                         cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", ".", " "))

rm(curr.el, temp)

# Performing stats on comp. class
comp.stats = NULL

for(curr.comp in unique(comp.class$variable)){
  temp = comp.class[which(comp.class$variable %in% curr.comp),]
  max.val = max(temp$value)
  temp = wilcox.test(value~Treatment, data = temp, alternative = "two.sided")
  temp = data.frame(Comparison = curr.comp, W = temp$statistic, p.value = temp$p.value, max.val = max.val)
  comp.stats = rbind(comp.stats, temp)
}

comp.stats$p.value = p.adjust(comp.stats$p.value, method = "fdr")
comp.stats$Symbol = symnum(comp.stats$p.value, corr = FALSE, na = FALSE, 
                           cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", ".", " "))

rm(curr.comp, temp)

# Making boxplots for compound class and elem. comp.
el.plot = ggplot(data = el.comp)+
  geom_boxplot(aes(x = variable, y = value, color = Treatment))+ 
  theme_bw() + xlab(NULL) + ylab("Relative Abundance (%)") + labs(color = "Treatment:")+
  geom_text(data = el.stats, aes(x = Comparison, y = max.val+2, label = as.character(Symbol)), size = 7)+
  theme(text = element_text(color = "black", size = 12),
        axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(color = "black", size = 14),
        legend.text = element_text(color = "black", size = 14),
        legend.title = element_text(color = "black", size = 14),
        panel.border = element_rect(color = "black"),
        axis.ticks = element_line(colour = "black"), 
        panel.grid = element_blank(),
        panel.background = element_blank())

comp.plot = ggplot(data = comp.class)+
  geom_boxplot(aes(x = variable, y = value, color = Treatment))+ 
  theme_bw() + xlab(NULL) + ylab("Relative Abundance (%)") + labs(color = "Treatment:")+
  geom_text(data = comp.stats, aes(x = Comparison, y = max.val+2, label = as.character(Symbol)), size = 7)+
  theme(text = element_text(color = "black", size = 12),
        axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(color = "black", size = 14),
        legend.text = element_text(color = "black", size = 14),
        legend.title = element_text(color = "black", size = 14),
        panel.border = element_rect(color = "black"),
        axis.ticks = element_line(colour = "black"), 
        panel.grid = element_blank(),
        panel.background = element_blank())

ggarrange(el.plot, comp.plot, ncol = 1, common.legend = T)


# ########################### #
#### Multivariate Analyses ####
# ########################### #

# Principal component analysis
pca = prcomp(x = t(data))

# Everything below this line should not change very much whether we use PCA or NMDS
ordination.scores = scores(pca) # Works with both PCA and NMDS, change the object accordingly
ordination.scores = as.data.frame(ordination.scores) # ggplot doesn't like matrices - needs to be converted to a data frame
ordination.scores$Treatment = factors$Treatment # Adding in sample type to our ordination scores object

# We have everything necessary for ggplot - we want to plot PC1 and PC2
ggplot(data = ordination.scores, aes(x = PC1, y = PC2, color = Treatment))+
  xlab(paste0("PC1 (", summary(pca)$importance[2,1]*100, "%)"))+
  ylab(paste0("PC2 (", summary(pca)$importance[2,2]*100, "%)"))+
  geom_point(size = 2) + 
  scale_color_manual(values = c("orange","lightblue")) +  # Colors for Treatment levels
  
  theme_bw()+
  theme(text = element_text(color = "black", size = 12),
        axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(color = "black", size = 14),
        legend.text = element_text(color = "black", size = 14),
        panel.border = element_rect(color = "black"),
        axis.ticks = element_line(colour = "black"), 
        panel.grid = element_blank(),
        panel.background = element_blank())
ggsave('W_vs_D_PCA.pdf',width = 20, height = 15)
# PERMANOVA
dist = vegdist(t(data), method = "euclidean", binary = T)
perm = adonis(dist~factors$Treatment, permutations = 999)

# Beta-dispersion analysis
beta.disp = betadisper(dist, group = factors$Treatment)
beta.disp = data.frame(Type = as.character(beta.disp$group), Distance = as.numeric(beta.disp$distances),
                       stringsAsFactors = F)
beta.stats = wilcox.test(Distance~Type, data = beta.disp)

# Plotting beta-diserpsion
ggplot(data = beta.disp, aes(x = Type, y = Distance))+
  geom_boxplot(aes(color = Type))+
  xlab(NULL)+ylab("Distance to Centroid")+
  theme_bw()+theme(legend.position = "none",
                   axis.text = element_text(color = "black", size = 11),
                   axis.title = element_text(color = "black", size = 13),
                   axis.ticks = element_line(color = "black"),
                   panel.border = element_rect(color = "black"),
                   panel.background = element_blank(),
                   panel.grid = element_blank())
