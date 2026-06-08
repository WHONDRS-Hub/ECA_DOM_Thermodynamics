rm(list=ls(all = TRUE))

#libraries
library(tidyverse)
# read resppiration data
resp_data = read_csv(paste0('EC_Data_Package/Sample_Data/EC_Sediment_SpC_pH_Temp_Respiration.csv'), 
                       comment = '#', na = c('N/A', -9999)) %>%
  slice(-(1:11)) %>%
  mutate(
    Theoretical = dplyr::case_when(
      !is.na(Methods_Deviation) & stringr::str_detect(Methods_Deviation, "RATE_005") ~ "yes",
      TRUE ~ "no"
    ),
    Theoretical = factor(Theoretical, levels = c("no", "yes"))
  ) %>% 
  #mutate(across(-Sample_Name:-Material, as.numeric)) %>%
  mutate(Treatment = case_when(
    str_detect(Sample_Name, "-W") ~ "Wet",
    str_detect(Sample_Name, "-D") ~ "Dry",
    TRUE ~ NA_character_)) %>%
  mutate(site = gsub('-W|-D', '', Sample_Name)) %>%
  dplyr::select(Sample_Name, Treatment,  Theoretical,Respiration_Rate_mg_DO_per_kg_per_H, DO_Concentration_At_Incubation_Time_Zero)

# read DOsta at P and temp duirng calibration 
do_sat = readxl::read_excel(paste0('EC_Data_Package/ec_fast_rate_calculations.xlsx')) %>%
  dplyr::select(Sample_Name = Sample_ID, DO_sat = DO_sat_mg_L)
#read moisture
moi_data = read_csv(paste0('EC_Data_Package/Sample_Data/EC_Sediment_Gravimetric_Moisture.csv'), 
                    comment = '#', na = c('N/A', -9999)) %>%
  slice(-(1:11)) %>%
  mutate(across(-Sample_Name:-Material, as.numeric)) %>%
  dplyr::select(Sample_Name, Final_Water_Mass_g,Dry_Sediment_Mass_g) %>%
  #merge with do sat via sample name
  left_join(do_sat, by = 'Sample_Name') %>%
  mutate(Displacement = Dry_Sediment_Mass_g/2.65,
         Water_added = 50 - Displacement - Final_Water_Mass_g,
         Mixed_DO_Calculated = ((0*Final_Water_Mass_g)+(DO_sat*Water_added))/(Water_added + Final_Water_Mass_g))

# Join moi data with resp data
full_data = resp_data %>%
  left_join(moi_data, by = 'Sample_Name') %>%
  dplyr::filter(
    !stringr::str_detect(
      Sample_Name,
      "^(EC_023|EC_012|EC_011|EC_053|EC_057|EC_052)"
    )
  )
  #Remove EC_011, EC_012, EC_023, EC_052, EC_053, and EC_057 for having too much water added (11 and 12), no mg/kg calculation (23), and being duplicated NEON sites (52, 53, and 57)

full_data = full_data %>%
  #remove rows in NA in respiration rates
  dplyr::filter(!is.na(Respiration_Rate_mg_DO_per_kg_per_H))


library(dplyr)
library(ggplot2)

# Fit linear model
mod <- lm(Mixed_DO_Calculated ~ DO_sat, data = full_data)
mod_sum <- summary(mod)

# Extract regression stats
r2 <- mod_sum$r.squared
p_val <- coef(mod_sum)[2, 4]

# Make label text
stats_label <- paste0(
  "R² = ", round(r2, 3),
  "\np = ", signif(p_val, 3)
)

# Plot
ggplot(full_data, aes(x = DO_sat, y = Mixed_DO_Calculated, color =  Theoretical)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  annotate(
    "text",
    x = Inf, y = Inf,
    label = stats_label,
    hjust = 1.1, vjust = 1.5,
    size = 5
  ) +
  labs(
    x = "DO_sat",
    y = "Mixed_DO_Calculated",
    color = "Theorethical",
    title = "DO_sat vs Mixed_DO_Calculated"
  ) +
  theme_bw()


# Plot
ggplot(full_data, aes(x = as.numeric(DO_Concentration_At_Incubation_Time_Zero), y = Mixed_DO_Calculated, color =  Theoretical)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  annotate(
    "text",
    x = Inf, y = Inf,
    label = stats_label,
    hjust = 1.1, vjust = 1.5,
    size = 5
  ) +
  labs(
    x = "Measured DO",
    y = "Mixed_DO_Calculated",
    color = "Theorethical",
    title = "DO measured time zero vs Mixed_DO_Calculated"
  ) +
  theme_bw()

full_data = full_data %>%
  mutate(DO_diff = Mixed_DO_Calculated - DO_sat)

#Plot histogram of do_diff
ggplot(full_data, aes(x = DO_diff, fill = Theoretical)) +
  geom_histogram(position = "dodge", bins = 20, alpha = 0.8) +
  labs(
    x = "Difference (Mixed_DO_Calculated - DO_sat)",
    y = "Count",
    fill = "Theoretical",
    title = "Histogram of Difference between Mixed_DO_Calculated and DO_sat"
  ) +
  theme_bw()


full_data = full_data %>%
  mutate(DO_diff_measured = Mixed_DO_Calculated - as.numeric(DO_Concentration_At_Incubation_Time_Zero))

#Plot histogram of do_diff
ggplot(full_data, aes(x = DO_diff_measured, fill = Theoretical)) +
  geom_histogram(position = "dodge", bins = 20, alpha = 0.8) +
  labs(
    x = "Difference (Mixed_DO_Calculated - DOmeasured)",
    y = "Count",
    fill = "Theoretical",
  ) +
  theme_bw()


ggplot(full_data, aes(x = as.numeric(DO_Concentration_At_Incubation_Time_Zero), y = DO_diff_measured, color =  Theoretical)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(
    x = "DO measured",
    y = "Difference (Mixed_DO_Calculated - DO measured)",
    color = "Theoretical",
  ) +
  theme_bw()  



library(dplyr)
library(ggplot2)

# make x variable once
full_data2 <- full_data %>%
  mutate(Measured_DO = as.numeric(DO_Concentration_At_Incubation_Time_Zero))

# calculate separate regression stats for each Theoretical group
stats_df <- full_data2 %>%
  group_by(Theoretical) %>%
  do({
    model <- lm(Mixed_DO_Calculated ~ Measured_DO, data = .)
    data.frame(
      r2 = summary(model)$r.squared,
      p = summary(model)$coefficients[2, 4]
    )
  }) %>%
  ungroup() %>%
  mutate(
    label = paste0(
      "Theoretical = ", Theoretical,
      "\nR² = ", round(r2, 3),
      "\np = ", signif(p, 3)
    ),
    x = c(12.8, 12.8),   # adjust if needed
    y = c(13.0, 11.8)    # separate label positions
  )

p1 = ggplot(
  full_data2,
  aes(x = Measured_DO, y = Mixed_DO_Calculated, color = Theoretical)
) +
  geom_point(size = 3, alpha = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  geom_smooth(method = "lm", se = FALSE) +
  geom_text(
    data = stats_df,
    aes(x = x, y = y, label = label, color = Theoretical),
    hjust = 1,
    size = 5,
    show.legend = FALSE
  ) +
  labs(
    x = "Measured DO",
    y = "Mixed_DO_Calculated",
    color = "Theoretical",
    title = "DO measured time zero vs Mixed_DO_Calculated"
  ) +
  scale_x_continuous(limits = c(0, 13)) +
  scale_y_continuous(limits = c(5.5, 13.5)) +
  theme_bw()


# subset to theoretical only
theory_only <- full_data2 %>%
  filter(Theoretical == "yes")

# fit model
mod_theory <- lm(Mixed_DO_Calculated ~ Measured_DO, data = theory_only)

# extract R2
r2_theory <- summary(mod_theory)$r.squared

# make label
r2_label <- paste0("R² = ", round(r2_theory, 3))

# plot
p_theory <- ggplot(
  theory_only,
  aes(
    x = Measured_DO,
    y = Mixed_DO_Calculated
  )
) +
  
  # points
  geom_point(
    color = "#1f9aa0",
    size = 3,
    alpha = 0.9
  ) +
  
  # regression line
  geom_smooth(
    method = "lm",
    se = FALSE,
    color = "black",
    linewidth = 1
  ) +
  
  # 1:1 line
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "red",
    alpha = 0.7
  ) +
    annotate(
    "text",
    x = 8.6,
    y = 0.6,
    label = r2_label,
    size = 5,
    fontface = "bold",
    hjust = 1
  ) +
  
  labs(
    x = expression("Measured O"[2]*" (mg L"^-1*")"),
    y = expression("Mixed O"[2]*" Calculated (mg L"^-1*")")
  ) +
  
  coord_cartesian(
    xlim = c(-0.5, 9),
    ylim = c(-0.5, 9)
  ) +
  
  scale_x_continuous(
    breaks = seq(0, 10, by = 2.5)
  ) +
  
  scale_y_continuous(
    breaks = seq(0, 10, by = 2.5)
  ) +
  
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

p_theory

ggsave(
  "Figures/FigureS4.png",
  p_theory,
  width = 10,
  height = 10,
  dpi = 300
)

ggsave(
  "Figures/FigureS4.pdf",
  p_theory,
  width = 10,
  height = 10,
  device = cairo_pdf
)

