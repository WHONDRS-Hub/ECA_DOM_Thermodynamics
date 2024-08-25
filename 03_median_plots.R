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
github_path = 'C:/Users/gara009/OneDrive - PNNL/Documents/GitHub/ECA_DOM_Thermodynamics/'
path = 'C:/Users/gara009/OneDrive - PNNL/Documents/GitHub/ECA_Multireactor_Incubations/Data/Cleaned Data/'

# ====== Read in data ======
# Processed ICR Data
icr.data = read.csv(list.files(pattern = "*Median_per_site_", path = github_path, full.names = T))
icr.data$Samples = gsub('_SIR-D|_SIR-W','_INC',icr.data$Samples)

icr.data.2 = read.csv(list.files(pattern = "*Median_metrics_", path = github_path, full.names = T)) 


other.data = read.csv(list.files(pattern = "2024-05-29_Medians_ECA", path = path, full.names = T), row.names = 1)
names(other.data)[1] = 'Samples'
names(other.data)[2] = 'Treatment'

effect.data = read.csv(list.files(pattern = "2024-05-29_Effect_Median_ECA", path = path, full.names = T), row.names = 1)
names(effect.data)[1] = 'Samples'
names(effect.data)[2] = 'Treatment'
# ==== Merge data ===
df = merge(icr.data, other.data, by = c('Samples','Treatment'), all = T)

df.effect = merge(icr.data, effect.data, by = c('Samples','Treatment'), all = T)

cube_root <- function(x) sign(x) * (abs(x))^(1/3)
cube_data = other.data %>% 
  mutate(across(where(is.numeric), cube_root)) %>% 
  rename_with(.cols = c(median_SpC:median_mean_ssa), .fn = ~ paste0("cube_", .x))

df_cube = merge(icr.data, cube_data, by = c('Samples','Treatment'), all = T)
# ===== Plots =====
# === Scatter plots ===
ggplot(df_cube, aes(y = cube_median_Respiration_Rate_mg_DO_per_L_per_H, x = Median_Gibbs, color = as.factor(Treatment))) +
  geom_point() +
  #geom_smooth(method = "lm", se = FALSE, color = "black") +
  theme_bw() +
  labs(
    title = " ",
    y = "cube root Median Respiration Rate (mg DO/L/hr)",
    x = "Median Gibbs (kJ/C mol)",
    color = 'Treatment'
  )

ggplot(df_cube, aes(y = cube_median_Respiration_Rate_mg_DO_per_L_per_H, x = Median_Lambda, color = as.factor(Treatment))) +
  geom_point() +
  #geom_smooth(method = "lm", se = FALSE, color = "black") +
  theme_bw() +
  labs(
    title = " ",
    y = "cube root Median Respiration Rate (mg DO/L/hr)",
    x = "Median Lambda",
    color = 'Treatment'
  )

ggplot(df_cube, aes(y = cube_median_Respiration_Rate_mg_DO_per_L_per_H, x = Median_Gibbs_Compound, color = as.factor(Treatment))) +
  geom_point() +
  #geom_smooth(method = "lm", se = FALSE, color = "black") +
  theme_bw() +
  labs(
    title = " ",
    y = "cube root Median Respiration Rate (mg DO/L/hr)",
    x = "Median Gibbs per compound",
    color = 'Treatment'
  )


# ==== boxplots with medians per sample  ====

ggplot(icr.data.2, aes(y = Median_Gibbs, x = Treatment, fill = as.factor(Treatment))) +
  geom_boxplot() +
  theme_bw() +
  labs(
    title = " ",
    x = "Median Gibbs (kJ/C mol)",
    fill = 'Treatment'
  )

ggplot(icr.data.2, aes(y = Median_Lambda, x = Treatment, fill = as.factor(Treatment))) +
  geom_boxplot() +
  theme_bw() +
  labs(
    title = " ",
    y = "Median Lambda",
    fill = 'Treatment'
  )

ggplot(icr.data.2, aes(y = Median_Gibbs_Compound, x = Treatment, fill = as.factor(Treatment))) +
  geom_boxplot() +
  theme_bw() +
  labs(
    title = " ",
    y = "Median Gibbs per compound",
    fill = 'Treatment'
  )

# ==== Boxplots with medians per treatment and site  ====

ggplot(df, aes(y = Median_Gibbs, x = Treatment, fill = as.factor(Treatment))) +
  geom_boxplot() +
  theme_bw() +
  labs(
    title = " ",
    x = "Median Gibbs (kJ/C mol)",
    fill = 'Treatment'
  )

ggplot(df, aes(y = Median_Lambda, x = Treatment, fill = as.factor(Treatment))) +
  geom_boxplot() +
  theme_bw() +
  labs(
    title = " ",
    y = "Median Lambda",
    fill = 'Treatment'
  )

ggplot(df, aes(y = Median_Gibbs_Compound, x = Treatment, fill = as.factor(Treatment))) +
  geom_boxplot() +
  theme_bw() +
  labs(
    title = " ",
    y = "Median Gibbs per compound",
    fill = 'Treatment'
  )

# ==== PCA with medians of all reps =====
factors = icr.data.2 %>%
  dplyr::select(c(Samples,Location,Treatment))

data.pca = icr.data.2 %>%
  column_to_rownames(var = 'Samples') %>%
  dplyr::select(c(-Location,-Treatment))


# Check for missing values
if (sum(is.na(data.pca)) > 0) {
  data.pca <- na.omit(data.pca)  # Remove rows with missing values
}

data.pca <- prcomp(data.pca, scale = TRUE,
                          center = TRUE, retx = T)

# 
# Extract the proportion of variance explained
explained_variance <- summary(data.pca)$importance[2, ]

# Create axis labels with the percentage of explained variance
x_label <- paste0("PC1 (", round(explained_variance[1] * 100, 1), "%)")
y_label <- paste0("PC2 (", round(explained_variance[2] * 100, 1), "%)")

# Extract loadings (rotation) for plotting arrows
loadings <- as.data.frame(data.pca$rotation)

# Scale loadings for better visualization in the plot
loading_scale_factor <- 5
loadings_scaled <- loadings * loading_scale_factor

# === Making the PCA plot ====
pca = as.data.frame(scores(data.pca)) # Converting to PCA scores in order to plot using ggplot
pca = cbind(factors, pca)

pca %>%
  ggplot(aes(x = PC1, y = PC2))+
  geom_point(aes(color = Treatment),size = 2) +
  theme_bw() +
  labs(
    x = x_label,
    y = y_label)+
  # Add arrows for the loadings
  geom_segment(data = loadings_scaled, aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.3, "cm")), color = "red") +
  # Add labels for the arrows
  geom_text(data = loadings_scaled, aes(x = PC1, y = PC2, label = rownames(loadings_scaled)), 
            vjust = -0.5, color = "red")


# ===== Look at a heatmap for the PCA loadings =====

loadings = data.pca$rotation

print(loadings)

loadings_df <- as.data.frame(loadings) %>% 
  select(c(PC1, PC2, PC3, PC4, PC5)) 

# Generate heatmap

row_names <- rownames(loadings_df)
loadings_df$Variable <- row_names

# Melt the dataframe for plotting
loadings_melted <- reshape2::melt(loadings_df, id.vars = "Variable")

# Plotting the heatmap
ggplot(loadings_melted, aes(x = variable, y = Variable, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(x = "Principal Component", y = "Variable", fill = "Loadings") +
  ggtitle("PCA Loadings Heatmap")+
  geom_text(aes(label = round(value, 2)), color = "black", size = 3)

# ==== PCA with medians per treatment =====
factors = icr.data %>%
  dplyr::select(c(Samples,Treatment))

data.pca = icr.data %>%
  dplyr::select(c(-Samples,-Treatment,-Number_of_reps))


# Check for missing values
if (sum(is.na(data.pca)) > 0) {
  data.pca <- na.omit(data.pca)  # Remove rows with missing values
}

data.pca <- prcomp(data.pca, scale = TRUE,
                   center = TRUE, retx = T)

# 
# Extract the proportion of variance explained
explained_variance <- summary(data.pca)$importance[2, ]

# Create axis labels with the percentage of explained variance
x_label <- paste0("PC1 (", round(explained_variance[1] * 100, 1), "%)")
y_label <- paste0("PC2 (", round(explained_variance[2] * 100, 1), "%)")

# Extract loadings (rotation) for plotting arrows
loadings <- as.data.frame(data.pca$rotation)

# Scale loadings for better visualization in the plot
loading_scale_factor <- 5
loadings_scaled <- loadings * loading_scale_factor

# === Making the PCA plot ====
pca = as.data.frame(scores(data.pca)) # Converting to PCA scores in order to plot using ggplot
pca = cbind(factors, pca)

pca %>%
  ggplot(aes(x = PC1, y = PC2))+
  geom_point(aes(color = Treatment),size = 2) +
  theme_bw() +
  labs(
    x = x_label,
    y = y_label)+
  # Add arrows for the loadings
  geom_segment(data = loadings_scaled, aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.3, "cm")), color = "red") +
  # Add labels for the arrows
  geom_text(data = loadings_scaled, aes(x = PC1, y = PC2, label = rownames(loadings_scaled)), 
            vjust = -0.5, color = "red")

# ==== Make a classification for low-mid-high respiration ===

# Define the breaks and labels
breaks <- c(-Inf, 5, 50, Inf)
labels <- c("low", "mid", "high")

# Use the cut function to create the categorical variable
df$respiration_rate_category <- cut(df$median_Respiration_Rate_mg_DO_per_L_per_H, 
                                    breaks = breaks, labels = labels, right = F)

ggplot(df, aes(y = Median_Lambda, x = respiration_rate_category)) +
  geom_boxplot() +
  theme_bw()

ggplot(df, aes(y = median_NPOC_mg_C_per_L, x = respiration_rate_category)) +
  geom_boxplot() +
  theme_bw()

ggplot(df, aes(y = median_NPOC_mg_C_per_L, x = respiration_rate_category, fill = Treatment)) +
  geom_boxplot() +
  theme_bw()

ggplot(df, aes(y = median_NPOC_mg_C_per_L, x = Treatment, fill = Treatment)) +
  geom_boxplot() +
  theme_bw()

test2 <- df %>%
  filter(Samples %in% test$Samples)

ggplot(test2, aes(y = median_NPOC_mg_C_per_L, x = Samples, fill = Treatment)) +
  geom_boxplot() +
  theme_bw()

test3 <- df %>%
  filter(!(Samples %in% test$Samples))
#%>%
 # filter(respiration_rate_category != 'high')

ggplot(test3, aes(y = median_NPOC_mg_C_per_L, x = Samples, fill = Treatment)) +
  geom_boxplot() +
  theme_bw().

ggplot(df, aes(y = median_Final_Gravimetric_Moisture, x = Treatment, fill = Treatment)) +
  geom_boxplot() +
  theme_bw()
# === Filter out high resp ====

df.test = filter(df, df$respiration_rate_category!= 'high')

ggplot(df.test, aes(x = (Median_Gibbs), y = log10(median_NPOC_mg_C_per_L), color = as.factor(Treatment))) +
  geom_point() +
  #geom_smooth(method = "lm", se = FALSE, color = "black") +
  theme_bw() 
