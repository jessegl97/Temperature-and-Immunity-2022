#Temperature and Immunity 2022 Analysis
  #Updated Oct 2025
  #JGL
  #Variability Models

rm(list=ls())
####read in + format data####
setwd('/Users/jesse/Documents/GitHub/Temperature-and-Immunity-2022/meanModels')

#Automatically install & load if missing
#install.packages("pacman")
if (!require("pacman")) install.packages("pacman")
library(pacman)

# Core tidyverse
p_load(tidyverse, patchwork, dplyr)

# Modeling
p_load(lme4, glmmTMB, DHARMa, effects, ordinal, AICcmodavg, emmeans, parameters, car, optimx, multcompView, ordinal, performance, mgcv)

#Writing
p_load(writexl, broom, broom.mixed)

# Visualization
p_load(ggplot2, gridExtra, gtsummary)


####Format theme####
#set colors for temperature
temp_colors <- c("#277DA1", "#F94144")
#set color scheme MG treatment
inf_colors <- c("black", "#BC6C25")

# full factorial
#treat_colors <- c("#277DA1", "#90BE6D", "#F94144", "#F9C74F")
treat_colors <- c("#90BE6D", "#277DA1", "#F9C74F", "#F94144")
#set theme
theme_set(
  theme_bw() +
    theme(
      axis.title.y = element_text(color = "black", size = 15),
      axis.text.y = element_text(color = "black", size = 15),
      axis.title.x = element_text(color = "black", size = 15),
      axis.text.x = element_text(color = "black", size = 15),
      legend.background = element_rect(linewidth = 0.25, linetype = "solid"),
      legend.position = "top",
      legend.direction = "horizontal",
      legend.title = element_text(size=15),
      legend.text = element_text(size=15)
    )
)

#Import Variability Functions:
source("Variability_Functions.R")

n_boot <- 1000

####Antibodies####
source("dataCleaning_antibody.R")

unique(ti.pi$dpi)
ti.pi$dpi <- as.factor(ti.pi$dpi)

#Set reference categories "Warm" and "Sham"
ti.pi$temp <- relevel(as.factor(ti.pi$temp), ref = "Warm")
ti.pi$treatment <- relevel(as.factor(ti.pi$treatment), ref = "Sham")

#omit seropositve birds from quarantine
unique(Qinf$band_number)

ti.pi <- ti.pi %>%
  filter(!band_number %in% c(2667, 2750))

ti.pi %>%
  filter(dpi ==28)%>%
  dplyr::select(treatment, temp, inf_temp)%>%
  tbl_summary(
  )%>%
  modify_header(
    label ~ "**All Birds**"
  )

#missing elisa data from 10 samples:
#DPI 9: 2684, 2893, 2632, 2842, 2847
#DPI 28: 2684, 2893, 2846, 2869, 2926

#2926 only on dpi 9; missing 28
ti.pi %>%
  filter(elisa_od>0.08) %>%
  dplyr::select(dpi, elisa_od, band_number, treatment)

ti.pi.full <- ti.pi

ti.pi <- ti.pi %>%
  filter(band_number != 2926)

# Remove rows with NA values in elisa_od
ti.pi.clean <- ti.pi %>%
  filter(!is.na(elisa_od))

ggplot(ti.pi.clean, aes(x=dpi, y=elisa_od, color=temp))+
  geom_point()

ti.pi.clean.count <- ti.pi.clean %>%
  group_by(dpi, temp, treatment)%>%
  summarise(n = length(elisa_od))

# Variability analysis for elisa_od
elisa_var <- ti.pi.clean %>%
  group_by(dpi, groups) %>%
  dplyr::reframe(
    temp = temp,
    band_number = band_number,
    elisa_od = elisa_od,
    treatment = treatment,
    groups = groups,
    
    # Calculations for elisa_od
    mean_elisa = mean(elisa_od, na.rm = TRUE),
    median_elisa = median(elisa_od, na.rm = TRUE),
    bird_cv_elisa = calculate_cv(elisa_od),
    bird_pv_elisa = calculate_pv(elisa_od),
    bird_v2_elisa = calculate_v2(elisa_od),
    
    # Bootstrapping for elisa_od
    cv_bootstrap_elisa = list(replicate(n_boot, calculate_cv(sample(elisa_od, replace = TRUE)))),
    pv_bootstrap_elisa = list(replicate(n_boot, calculate_pv(sample(elisa_od, replace = TRUE)))),
    v2_bootstrap_elisa = list(replicate(n_boot, calculate_v2(sample(elisa_od, replace = TRUE))))
  ) %>%
  mutate(
    # Confidence intervals for elisa_od
    cv_lower_ci_elisa = map_dbl(cv_bootstrap_elisa, ~ quantile(.x, 0.025, na.rm = TRUE)),
    cv_upper_ci_elisa = map_dbl(cv_bootstrap_elisa, ~ quantile(.x, 0.975, na.rm = TRUE)),
    pv_lower_ci_elisa = map_dbl(pv_bootstrap_elisa, ~ quantile(.x, 0.025, na.rm = TRUE)),
    pv_upper_ci_elisa = map_dbl(pv_bootstrap_elisa, ~ quantile(.x, 0.975, na.rm = TRUE)),
    v2_lower_ci_elisa = map_dbl(v2_bootstrap_elisa, ~ quantile(.x, 0.025, na.rm = TRUE)),
    v2_upper_ci_elisa = map_dbl(v2_bootstrap_elisa, ~ quantile(.x, 0.975, na.rm = TRUE))
  ) %>%
  dplyr::select(-cv_bootstrap_elisa, -pv_bootstrap_elisa, -v2_bootstrap_elisa) %>%
  ungroup()

# Summary statistics for elisa_od
summary_elisa_tibble <- elisa_var %>%
  group_by(dpi, groups) %>%
  dplyr::reframe(
    mean_elisa = mean(mean_elisa, na.rm = TRUE),
    median_elisa = mean(median_elisa, na.rm = TRUE),
    mean_bird_cv_elisa = mean(bird_cv_elisa, na.rm = TRUE),
    mean_bird_pv_elisa = mean(bird_pv_elisa, na.rm = TRUE),
    mean_bird_v2_elisa = mean(bird_v2_elisa, na.rm = TRUE),
    mean_cv_lower_ci_elisa = mean(cv_lower_ci_elisa, na.rm = TRUE),
    mean_cv_upper_ci_elisa = mean(cv_upper_ci_elisa, na.rm = TRUE),
    mean_pv_lower_ci_elisa = mean(pv_lower_ci_elisa, na.rm = TRUE),
    mean_pv_upper_ci_elisa = mean(pv_upper_ci_elisa, na.rm = TRUE),
    mean_v2_lower_ci_elisa = mean(v2_lower_ci_elisa, na.rm = TRUE),
    mean_v2_upper_ci_elisa = mean(v2_upper_ci_elisa, na.rm = TRUE),
    n_individuals = n_distinct(band_number) #calculate number of individuals in each group
  )

ggplot(summary_elisa_tibble, aes(x = groups, color = groups)) +
  geom_point(aes(y = mean_bird_pv_elisa, shape = "PV"), size = 3) +
  geom_errorbar(aes(y = mean_bird_pv_elisa, ymin = mean_pv_lower_ci_elisa, ymax = mean_pv_upper_ci_elisa), width = 0.1) +
  geom_point(aes(y = mean_bird_cv_elisa, shape = "CV"), size = 2) +
  
  
  labs(
    y = "Variability",
    x = "Temperature",
    shape = "Metric Type",
    color = "Temperature",
    title = "Variability in Antibody Levels"
  ) +
  scale_color_manual(values = c(treat_colors))+
  scale_shape_manual(name = "Variability Metric",
                     values = c("PV" = 16, "CV" = 1)) +
  facet_wrap(~dpi)+
  theme_bw()+
  theme(
    axis.text.x = element_text(angle = 45, hjust=1)
  )

ggplot(ti.pi.clean, aes(x = groups, y=elisa_od, color = groups)) +
  geom_boxplot(width=0.5, outlier.shape=3)+
  geom_jitter(size = 2, width=0.25, alpha=0.5) +

  labs(
    y = "ELISA OD",
    x = "Treatment Group",
    color = "Treatment Group",
    title = "Mean Antibody Levels"
  ) +
  scale_color_manual(values = c(treat_colors))+
  facet_wrap(~dpi)+
  theme_bw()+
  theme(
    axis.text.x = element_text(angle = 45, hjust=1)
  )

ti.pi.clean$groups <- as.factor(ti.pi.clean$groups)
#Brown-Forsythe
bf_ab_treat <- leveneTest(elisa_od ~ treatment, data = ti.pi.clean %>% filter(dpi == 28), center = median)
bf_ab_temp <- leveneTest(elisa_od ~ temp, data = ti.pi.clean %>% filter(dpi == 28), center = median)
bf_ab_groups <- leveneTest(elisa_od ~ groups, data = ti.pi.clean %>% filter(dpi == 28), center = median)

leveneTest(elisa_od ~ interaction(groups, dpi), data = ti.pi.clean, center = median)

