#Temperature and Immunity 2022 Analysis
  #Updated Oct 2025
  #JGL
  #Variability Models

rm(list=ls())
####read in + format data####
setwd('/Users/jesse/Documents/GitHub/Temperature-and-Immunity-2022/Analysis/')

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
p_load(ggplot2, gridExtra, gtsummary, cowplot)


####Format theme####
#set colors for temperature
temp_colors <- c(
  "Cold" = "#277DA1", 
  "Warm" = "#F94144")
#set color scheme MG treatment
inf_colors <- c("black", "#BC6C25")

# full factorial
#treat_colors <- c("#277DA1", "#90BE6D", "#F94144", "#F9C74F")
treat_colors <- c(
  "Cold Inoculated" = "#277DA1",  # blue
  "Cold Control"  = "#90BE6D",  # green
  "Warm Inoculated" = "#F94144",  # red
  "Warm Control"  = "#F9C74F"   # yellow
)
#set theme
theme_set(
  theme_bw() +
    theme(
      axis.title.y = element_text(color = "black", size = 15),
      axis.text.y = element_text(color = "black", size = 15),
      axis.title.x = element_text(color = "black", size = 15),
      axis.text.x = element_text(color = "black", size = 15),
      legend.background = element_rect(linewidth = 0.25, linetype = "solid"),
      legend.position = "right",
      legend.direction = "vertical",
      legend.title = element_text(size=15),
      legend.text = element_text(size=15)
    )
)


#Import Variability Functions:
source("r_scripts/Variability_Functions.R")

n_boot <- 1000

####Antibodies####
source("r_scripts/dataCleaning_antibody.R")
ti %>%
  dplyr::select(dpi, treatment, temp, groups)%>%
  tbl_summary(
    by="dpi"
  )%>%
  modify_header(
    label ~ "**All Birds**"
  )

unique(ti.rm$band_number)

#Antibody analysis sample sizes; see dataCleaning_antibody.R for removal breakdown
ti.rm %>%
  dplyr::select(dpi, treatment, temp, groups)%>%
  tbl_summary(
    by="dpi"
  )%>%
  modify_header(
    label ~ "**Birds for Antibodies**"
  )

ti.ab <- ti.rm

ti.ab <- ti.ab %>% 
  mutate(groups = factor(groups,
                levels = c("Warm Control", "Warm Inoculated",
                           "Cold Control", "Cold Inoculated")))

unique(ti.ab$dpi)
ti.ab$dpi <- as.factor(ti.ab$dpi)

#Set reference categories "Warm" and "Sham"
ti.ab$temp <- relevel(as.factor(ti.ab$temp), ref = "Warm")
ti.ab$treatment <- relevel(as.factor(ti.ab$treatment), ref = "Control")
ti.ab$dpi.f <- as.factor(ti.ab$dpi)

# Variability analysis for elisa_od
elisa_var <- ti.ab %>%
  group_by(dpi, groups, sex) %>%
  dplyr::reframe(
    temp = temp,
    sex = sex,
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
  group_by(dpi, groups, sex, temp) %>%
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

dodge <- position_dodge(width=0.75)

var.ab <- ggplot(summary_elisa_tibble, aes(x = temp, color = groups)) +
  geom_errorbar(aes(y = mean_bird_pv_elisa, ymin = mean_pv_lower_ci_elisa, ymax = mean_pv_upper_ci_elisa, group = dpi), width = 0, color="black", position=dodge) +
  geom_point(aes(y = mean_bird_pv_elisa, shape = "PV", group = dpi), size = 3, shape=1, color="black", stroke=1.5, position=dodge) +
  geom_point(aes(y = mean_bird_pv_elisa, shape = "PV", group = dpi, alpha=dpi), size = 3, shape=16, position=dodge) +

  #geom_point(aes(y = mean_bird_cv_elisa, shape = "CV", group=dpi), size = 2, position = dodge) +
  
  
  labs(
    y = "Variability (PV)",
    x = "Temperature",
    shape = "Metric Type",
    color = "Temperature",
    #title = "Variability in Antibody Levels"
  ) +
  scale_color_manual(values = c(treat_colors))+
  scale_alpha_manual(values = c(0.5, 0.75, 1))+
  theme(
    axis.text.x = element_text(angle = 45, hjust=1)
  )+
  facet_wrap(~sex)

var.ab

sex_colors <- c("#A76F98", "#578E3F")

ggplot(summary_elisa_tibble, aes(x = sex, color = sex, shape=dpi)) +
  geom_errorbar(aes(y = mean_bird_pv_elisa, ymin = mean_pv_lower_ci_elisa, ymax = mean_pv_upper_ci_elisa, group = dpi), width = 0, color="black", position=dodge) +
  geom_point(aes(y = mean_bird_pv_elisa, group = dpi), size = 4, color="black", stroke=1, position=dodge) +
  geom_point(aes(y = mean_bird_pv_elisa, group = dpi), size = 3, position=dodge) +
  
  #geom_point(aes(y = mean_bird_cv_elisa, shape = "CV", group=dpi), size = 2, position = dodge) +
  
  
  labs(
    y = "Variability (PV)",
    x = "Treatment Group",
    shape = "Metric Type",
    color = "Temperature",
    #title = "Variability in Antibody Levels"
  ) +
  scale_color_manual(values = c(sex_colors))+
  scale_shape_manual(values = c(1, 17, 16))+
  # theme(
  #   axis.text.x = element_text(angle = 45, hjust=1)
  # )+
  facet_wrap(~temp)


mean.ab <- ggplot(ti.ab, aes(x = sex, y=elisa_od, color = groups)) +
  geom_boxplot(width=0.5, outlier.shape=3)+
  geom_jitter(size = 2, width=0.25, alpha=0.5) +

  labs(
    y = "ELISA OD",
    x = "Sex",
    color = "Treatment Group",
    #title = "Mean Antibody Levels"
  ) +
  scale_color_manual(values = c(treat_colors))+
  facet_wrap(~sex)+
  theme_bw()+
  theme(
    axis.text.x = element_text(angle = 45, hjust=1)
  )


var.ab + theme(legend.position = "none")+
  mean.ab

ti.ab$groups <- as.factor(ti.ab$groups)
#Brown-Forsythe; are antibody levels equally variable across treatments?
bf_ab_treat <- leveneTest(elisa_od ~ treatment, data = ti.ab %>% filter(dpi == 28), center = median)
bf_ab_temp <- leveneTest(elisa_od ~ temp, data = ti.ab %>% filter(dpi == 28), center = median)
bf_ab_groups <- leveneTest(elisa_od ~ groups, data = ti.ab %>% filter(dpi == 28), center = median)


# coerce to data frames with a term column
df_treat  <- bf_ab_treat  %>% as.data.frame() %>% rownames_to_column("term")
df_temp   <- bf_ab_temp   %>% as.data.frame() %>% rownames_to_column("term")
df_groups <- bf_ab_groups %>% as.data.frame() %>% rownames_to_column("term")

# (optional) add context columns
df_treat$comparison  <- "treatment"
df_temp$comparison   <- "temp"
df_groups$comparison <- "groups"
df_treat$dpi <- df_temp$dpi <- df_groups$dpi <- "28"

# write to Excel with three sheets
# write_xlsx(
#   list(
#     by_treatment = df_treat,
#     by_temp      = df_temp,
#     by_groups    = df_groups
#   ),
#   path = "/Users/jesse/Documents/GitHub/Temperature-and-Immunity-2022/Analysis/Results/Variability/bf_elisa_dpi28.xlsx"
# )

####Phagocytosis####
source("r_scripts/dataCleaning_phago.R")
lvl <- c("Warm Control", "Warm Inoculated", "Cold Control", "Cold Inoculated")

ti.p %>%
  dplyr::select(dpi, treatment, temp, groups)%>%
  tbl_summary(
  )%>%
  modify_header(
    label ~ "**All Birds**"
  )


# Remove rows with NA values in phago_score
ti.p.clean <- ti.p %>%
  filter(!is.na(phago_score)) %>%
  mutate(
    groups = paste(temp, ifelse(treatment == "Inoculated", "Inoculated", "Control")),
    groups = str_trim(groups),
    groups = factor(groups, levels = lvl)
  )

# Count entries for dpi == 28 by temp and treatment
ti.p.clean.count <- ti.p.clean %>%
  group_by(temp, treatment) %>%
  summarize(n = length(phago_score))%>%
  ungroup()

# Variability analysis for phago_score
phago_var <- ti.p.clean %>%
  group_by(groups, temp, treatment) %>%
  dplyr::reframe(
    groups = groups,
    temp = temp,
    band_number = band_number,
    phago_score = phago_score,
    treatment = treatment,
    
    # Calculations for phago_score
    mean_phago = mean(phago_score, na.rm = TRUE),
    median_phago = median(phago_score, na.rm = TRUE),
    bird_cv_phago = calculate_cv(phago_score),
    bird_pv_phago = calculate_pv(phago_score),
    bird_v2_phago = calculate_v2(phago_score),
    
    # Bootstrapping for phago_score
    cv_bootstrap_phago = list(replicate(n_boot, calculate_cv(sample(phago_score, replace = TRUE)))),
    pv_bootstrap_phago = list(replicate(n_boot, calculate_pv(sample(phago_score, replace = TRUE)))),
    v2_bootstrap_phago = list(replicate(n_boot, calculate_v2(sample(phago_score, replace = TRUE))))
  ) %>%
  mutate(
    # Confidence intervals for phago_score
    cv_lower_ci_phago = map_dbl(cv_bootstrap_phago, ~ quantile(.x, 0.025, na.rm = TRUE)),
    cv_upper_ci_phago = map_dbl(cv_bootstrap_phago, ~ quantile(.x, 0.975, na.rm = TRUE)),
    pv_lower_ci_phago = map_dbl(pv_bootstrap_phago, ~ quantile(.x, 0.025, na.rm = TRUE)),
    pv_upper_ci_phago = map_dbl(pv_bootstrap_phago, ~ quantile(.x, 0.975, na.rm = TRUE)),
    v2_lower_ci_phago = map_dbl(v2_bootstrap_phago, ~ quantile(.x, 0.025, na.rm = TRUE)),
    v2_upper_ci_phago = map_dbl(v2_bootstrap_phago, ~ quantile(.x, 0.975, na.rm = TRUE))
  ) %>%
  dplyr::select(-cv_bootstrap_phago, -pv_bootstrap_phago, -v2_bootstrap_phago) %>%
  ungroup()

# Summary statistics for phago_score
summary_phago_tibble <- phago_var %>%
  group_by(groups, temp, treatment) %>%
  dplyr::reframe(
    mean_phago = mean(mean_phago, na.rm = TRUE),
    median_phago = mean(median_phago, na.rm = TRUE),
    mean_bird_cv_phago = mean(bird_cv_phago, na.rm = TRUE),
    mean_bird_pv_phago = mean(bird_pv_phago, na.rm = TRUE),
    mean_bird_v2_phago = mean(bird_v2_phago, na.rm = TRUE),
    mean_cv_lower_ci_phago = mean(cv_lower_ci_phago, na.rm = TRUE),
    mean_cv_upper_ci_phago = mean(cv_upper_ci_phago, na.rm = TRUE),
    mean_pv_lower_ci_phago = mean(pv_lower_ci_phago, na.rm = TRUE),
    mean_pv_upper_ci_phago = mean(pv_upper_ci_phago, na.rm = TRUE),
    mean_v2_lower_ci_phago = mean(v2_lower_ci_phago, na.rm = TRUE),
    mean_v2_upper_ci_phago = mean(v2_upper_ci_phago, na.rm = TRUE)
  )

# Plot variability for phago_score
ggplot(summary_phago_tibble, aes(x = groups, color = groups)) +
  geom_errorbar(aes(y = mean_bird_pv_phago, ymin = mean_pv_lower_ci_phago, ymax = mean_pv_upper_ci_phago), width = 0, color="black") +
  geom_point(aes(y = mean_bird_pv_phago), size = 3) +
  geom_point(aes(y = mean_bird_pv_phago), size = 3, shape=1, stroke=1.5, color="black") +
  
  # geom_point(aes(y = mean_bird_cv_phago, shape = "CV"), size = 2) +
  # geom_point(aes(y = mean_bird_v2_phago, shape = "V2"), size = 2) +
  
  labs(
    y = "Variability in Phagocytosis Score (PV)",
    x = "Treatment Group",
    # shape = "Metric Type",
    color = "Temperature",
    #title = "Variability in Phagocytosis"
  ) +
  scale_color_manual(values = treat_colors) +
  theme(
    axis.text.x = element_text(angle=45, hjust=1)
  )
  # scale_shape_manual(name = "Variability Metric",
  #                    values = c("PV" = 16, "CV" = 2, "V2" = 0)) +
  #coord_flip() +


#Brown-Forsythe; are phagocytosis levels equally variable across treatments?
bf_phago_treat <- leveneTest(phago_score ~ treatment, data = ti.p.clean, center = median)
bf_phago_temp <- leveneTest(phago_score ~ temp, data = ti.p.clean, center = median)
bf_phago_groups <- leveneTest(phago_score ~ groups, data = ti.p.clean, center = median)

#just inoculated
bf_phago_inoc <- leveneTest(phago_score ~ temp, data = ti.p.clean %>% filter(treatment == "Inoculated"), center = median)

# coerce to data frames with a term column
df_treat  <- bf_phago_treat  %>% as.data.frame() %>% rownames_to_column("term")
df_temp   <- bf_phago_temp   %>% as.data.frame() %>% rownames_to_column("term")
df_groups <- bf_phago_groups %>% as.data.frame() %>% rownames_to_column("term")

# (optional) add context columns
df_treat$comparison  <- "treatment"
df_temp$comparison   <- "temp"
df_groups$comparison <- "groups"
df_treat$dpi <- df_temp$dpi <- df_groups$dpi <- "28"

# write to Excel with three sheets
# write_xlsx(
#   list(
#     by_treatment = df_treat,
#     by_temp      = df_temp,
#     by_groups    = df_groups
#   ),
#   path = "/Users/jesse/Documents/GitHub/Temperature-and-Immunity-2022/Results/Variability/bf_phagocytosis.xlsx"
# )


####Eye Score####
source("dataCleaning_eyeScore.R")

ti.mg %>%
  dplyr::select(dpi, treatment, temp, groups)%>%
  tbl_summary(
    by="dpi"
  )%>%
  modify_header(
    label ~ "**All Birds**"
  )

ti.mg$band_number <- as.factor(ti.mg$band_number)

ti.mg <- ti.mg %>%
  mutate(
    groups = paste(temp, ifelse(treatment == "Inoculated", "Inoculated", "Control")),
    groups = str_trim(groups),
    groups = factor(groups, levels = lvl)
  )


#add small constant
ti.mg$total_eye_score1 <- ti.mg$total_eye_score+0.001
#All Birds
#No eye scores dpi -12, 3, 35 so get rid 
eye_score_var <- ti.mg %>%
  group_by(dpi, temp, treatment, groups) %>%
  filter(!dpi %in% c(-12, 3, 35) & treatment=="Inoculated")%>%
  dplyr::reframe(
    dpi = dpi,
    temp = temp,
    treatment = treatment,
    band_number = band_number,
    total_eye_score = total_eye_score,
    
    # Calculations for total_eye_score
    mean_eye_score = mean(total_eye_score, na.rm = TRUE),
    median_eye_score = median(total_eye_score, na.rm = TRUE),
    bird_cv_eye_score = calculate_cv(total_eye_score1),
    bird_pv_eye_score = calculate_pv(total_eye_score1),
    bird_v2_eye_score = calculate_v2(total_eye_score1),
    
    # Bootstrapping for total_eye_score
    cv_bootstrap_eye_score = list(replicate(n_boot, calculate_cv(sample(total_eye_score1, replace = TRUE)))),
    pv_bootstrap_eye_score = list(replicate(n_boot, calculate_pv(sample(total_eye_score1, replace = TRUE)))),
    v2_bootstrap_eye_score = list(replicate(n_boot, calculate_v2(sample(total_eye_score1, replace = TRUE))))
  ) %>%
  mutate(
    # Confidence intervals for total_eye_score
    cv_lower_ci_eye_score = map_dbl(cv_bootstrap_eye_score, ~ quantile(.x, 0.025, na.rm = TRUE)),
    cv_upper_ci_eye_score = map_dbl(cv_bootstrap_eye_score, ~ quantile(.x, 0.975, na.rm = TRUE)),
    pv_lower_ci_eye_score = map_dbl(pv_bootstrap_eye_score, ~ quantile(.x, 0.025, na.rm = TRUE)),
    pv_upper_ci_eye_score = map_dbl(pv_bootstrap_eye_score, ~ quantile(.x, 0.975, na.rm = TRUE)),
    v2_lower_ci_eye_score = map_dbl(v2_bootstrap_eye_score, ~ quantile(.x, 0.025, na.rm = TRUE)),
    v2_upper_ci_eye_score = map_dbl(v2_bootstrap_eye_score, ~ quantile(.x, 0.975, na.rm = TRUE))
  ) %>%
  dplyr::select(-cv_bootstrap_eye_score, -pv_bootstrap_eye_score, -v2_bootstrap_eye_score) %>%
  ungroup()

# Summary statistics for total_eye_score
summary_eye_score_tibble <- eye_score_var %>%
  group_by(dpi, temp, treatment, groups) %>%
  dplyr::reframe(
    mean_eye_score = mean(mean_eye_score, na.rm = TRUE),
    median_eye_score = mean(median_eye_score, na.rm = TRUE),
    mean_bird_cv_eye_score = mean(bird_cv_eye_score, na.rm = TRUE),
    mean_bird_pv_eye_score = mean(bird_pv_eye_score, na.rm = TRUE),
    mean_bird_v2_eye_score = mean(bird_v2_eye_score, na.rm = TRUE),
    mean_cv_lower_ci_eye_score = mean(cv_lower_ci_eye_score, na.rm = TRUE),
    mean_cv_upper_ci_eye_score = mean(cv_upper_ci_eye_score, na.rm = TRUE),
    mean_pv_lower_ci_eye_score = mean(pv_lower_ci_eye_score, na.rm = TRUE),
    mean_pv_upper_ci_eye_score = mean(pv_upper_ci_eye_score, na.rm = TRUE),
    mean_v2_lower_ci_eye_score = mean(v2_lower_ci_eye_score, na.rm = TRUE),
    mean_v2_upper_ci_eye_score = mean(v2_upper_ci_eye_score, na.rm = TRUE)
  )

# Plot variability for total_eye_score
ggplot(summary_eye_score_tibble, aes(x = dpi, color = temp)) +
  geom_point(aes(y = mean_bird_pv_eye_score), size = 2, position = position_dodge(width=0.75)) +
  geom_errorbar(aes(y = mean_bird_pv_eye_score, ymin = mean_pv_lower_ci_eye_score, ymax = mean_pv_upper_ci_eye_score), 
                width = 0.1, position = position_dodge(width=0.75)) +
  # geom_point(aes(y = mean_bird_cv_eye_score, shape = "CV"), size = 2) +
  # geom_point(aes(y = mean_bird_v2_eye_score, shape = "V2"), size = 2) +
  
  labs(
    y = "Variability in Eye Score (PV)",
    x = "Days Post-Inoculation",
    shape = "Metric Type",
    color = "Temperature",
    #title = "Variability in Total Eye Score by DPI and Temperature"
  ) +
  scale_y_continuous(limits=c(0,1))+
  scale_color_manual(values = temp_colors) 
  # scale_shape_manual(name = "Variability Metric",
  #                    values = c("PV" = 16, "CV" = 2, "V2" = 0)) +
  #facet_wrap(~groups)

ti.mg.m <- ti.mg %>%
  group_by(band_number, ever_diseased, temp, treatment, groups) %>%
  reframe(
    max_tes = max(total_eye_score1))

ggplot(ti.mg.m, aes(x = temp, y = max_tes, color = temp, fill=temp)) +
  geom_jitter(width=0.1, height=0, size = 2.6) +
  geom_boxplot(alpha=0.3, color="black", width=0.5)+
  stat_summary(aes(group=temp, y=max_tes, x=temp), fun="mean", geom = "point", shape = 1, size=3)+
  scale_color_manual(values = temp_colors) +
  scale_fill_manual(values = temp_colors) +
  labs(x = "Group", y = "Max Eyescore", color="Temperature", fill="Temperature")

#title = "Variability of Peak Eye Score Across Birds") +
#theme(axis.text.x = element_text(angle = 45, hjust = 1))

max_tes_var <- ti.mg.m %>%
  group_by(temp, treatment, groups) %>%             # add ever_diseased here if desired
  reframe(
    # raw values (for bootstrap sampling)
    vals = list(max_tes),
    
    # point estimates
    mean_max_tes = mean(max_tes, na.rm = TRUE),
    cv_max_tes   = calculate_cv(max_tes),
    pv_max_tes   = calculate_pv(max_tes),
    v2_max_tes   = calculate_v2(max_tes),
    
    # bootstrap replicates (n_boot defined elsewhere)
    cv_boot_max  = list(replicate(n_boot, calculate_cv(sample(max_tes, replace = TRUE)))),
    pv_boot_max  = list(replicate(n_boot, calculate_pv(sample(max_tes, replace = TRUE)))),
    v2_boot_max  = list(replicate(n_boot, calculate_v2(sample(max_tes, replace = TRUE))))
  ) %>%
  mutate(
    cv_max_lower = map_dbl(cv_boot_max, ~ quantile(.x, 0.025, na.rm = TRUE)),
    cv_max_upper = map_dbl(cv_boot_max, ~ quantile(.x, 0.975, na.rm = TRUE)),
    pv_max_lower = map_dbl(pv_boot_max, ~ quantile(.x, 0.025, na.rm = TRUE)),
    pv_max_upper = map_dbl(pv_boot_max, ~ quantile(.x, 0.975, na.rm = TRUE)),
    v2_max_lower = map_dbl(v2_boot_max, ~ quantile(.x, 0.025, na.rm = TRUE)),
    v2_max_upper = map_dbl(v2_boot_max, ~ quantile(.x, 0.975, na.rm = TRUE))
  ) %>%
  #select(-vals, -cv_boot_max, -pv_boot_max, -v2_boot_max) %>%
  ungroup()

ggplot(max_tes_var, aes(x = temp, y = pv_max_tes, color = temp)) +
  geom_errorbar(aes(ymin = pv_max_lower, ymax = pv_max_upper), width = 0.0) +
  geom_point(size = 2.6) +
  scale_color_manual(values = temp_colors) +
  scale_y_continuous(limits = c(0,1))+
  labs(x = "Group", y = "Variability in Max Eyescore (PV)", color="Temperature")+
       #title = "Variability of Peak Eye Score Across Birds") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Brown-Forsythe; are max eye scores equally variable across temperatures?
bf_max_ex_temp <- leveneTest(max_tes ~ temp, data = ti.mg.m, center = median)


####Fever **Needs Work :(**####
source("dataCleaning_fever.R")

ti.f %>%
  dplyr::select(dpi, treatment, temp, groups)%>%
  tbl_summary(
    by="dpi"
  )%>%
  modify_header(
    label ~ "**All Birds**"
  )

ti.f.clean <- ti.f %>%
  filter(!is.na(fever_score))

# Count entries by temp and treatment
ti.f.clean.count <- ti.f.clean %>%
  group_by(dpi, temp, treatment) %>%
  summarize(n = length(fever_score)) %>%
  ungroup()

lvl <- c("Warm Control", "Warm Inoculated”, “Cold Control", "Cold Inoculated")

ti.f.clean <- ti.f.clean %>%
  mutate(
    groups = paste(temp, ifelse(treatment == "Inoculated", "Inoculated", "Control")),
    groups = str_trim(groups)
  )

# Variability in fever_score over time
fever_var <- ti.f.clean %>%
  group_by(dpi, temp, treatment) %>%
  dplyr::reframe(
    temp = temp,
    band_number = band_number,
    fever_score = fever_score,
    treatment = treatment,
    groups = groups,
    
    # Calculations for fever_score
    mean_fever = mean(fever_score, na.rm = TRUE),
    median_fever = median(fever_score, na.rm = TRUE),
    bird_cv_fever = calculate_cv(fever_score),
    bird_pv_fever = calculate_pv(fever_score),
    bird_v2_fever = calculate_v2(fever_score),
    
    # Bootstrapping for fever_score
    cv_bootstrap_fever = list(replicate(n_boot, calculate_cv(sample(fever_score, replace = TRUE)))),
    pv_bootstrap_fever = list(replicate(n_boot, calculate_pv(sample(fever_score, replace = TRUE)))),
    v2_bootstrap_fever = list(replicate(n_boot, calculate_v2(sample(fever_score, replace = TRUE))))
  ) %>%
  mutate(
    # Confidence intervals for fever_score
    cv_lower_ci_fever = map_dbl(cv_bootstrap_fever, ~ quantile(.x, 0.025, na.rm = TRUE)),
    cv_upper_ci_fever = map_dbl(cv_bootstrap_fever, ~ quantile(.x, 0.975, na.rm = TRUE)),
    pv_lower_ci_fever = map_dbl(pv_bootstrap_fever, ~ quantile(.x, 0.025, na.rm = TRUE)),
    pv_upper_ci_fever = map_dbl(pv_bootstrap_fever, ~ quantile(.x, 0.975, na.rm = TRUE)),
    v2_lower_ci_fever = map_dbl(v2_bootstrap_fever, ~ quantile(.x, 0.025, na.rm = TRUE)),
    v2_upper_ci_fever = map_dbl(v2_bootstrap_fever, ~ quantile(.x, 0.975, na.rm = TRUE))
  ) %>%
  dplyr::select(-cv_bootstrap_fever, -pv_bootstrap_fever, -v2_bootstrap_fever) %>%
  ungroup()

# Summary statistics for fever_score
summary_fever_tibble <- fever_var %>%
  group_by(dpi, groups, temp, treatment) %>%
  dplyr::reframe(
    mean_fever = mean(mean_fever, na.rm = TRUE),
    median_fever = mean(median_fever, na.rm = TRUE),
    mean_bird_cv_fever = mean(bird_cv_fever, na.rm = TRUE),
    mean_bird_pv_fever = mean(bird_pv_fever, na.rm = TRUE),
    mean_bird_v2_fever = mean(bird_v2_fever, na.rm = TRUE),
    mean_cv_lower_ci_fever = mean(cv_lower_ci_fever, na.rm = TRUE),
    mean_cv_upper_ci_fever = mean(cv_upper_ci_fever, na.rm = TRUE),
    mean_pv_lower_ci_fever = mean(pv_lower_ci_fever, na.rm = TRUE),
    mean_pv_upper_ci_fever = mean(pv_upper_ci_fever, na.rm = TRUE),
    mean_v2_lower_ci_fever = mean(v2_lower_ci_fever, na.rm = TRUE),
    mean_v2_upper_ci_fever = mean(v2_upper_ci_fever, na.rm = TRUE)
  )

# Plot variability for fever_score over time
ggplot(summary_fever_tibble %>% filter(dpi != 3), aes(x = dpi, color = groups)) +
  geom_point(aes(y = mean_bird_pv_fever, shape = "PV"), size = 2, position=position_dodge(width=2)) +
  geom_errorbar(aes(y = mean_bird_pv_fever, ymin = mean_pv_lower_ci_fever, ymax = mean_pv_upper_ci_fever), width = 0.5, position=position_dodge(width=2)) +
  geom_point(aes(y = mean_bird_cv_fever, shape = "CV"), size = 2, position=position_dodge(width=2)) +
  #geom_point(aes(y = mean_bird_v2_fever, shape = "V2"), size = 2) +
  
  labs(
    y = "Variability",
    x = "Temperature",
    shape = "Metric Type",
    color = "Temperature",
    title = "Variability in Fever Score"
  ) +
  scale_color_manual(values = treat_colors) +
  scale_shape_manual(name = "Variability Metric",
                     values = c("PV" = 16, "CV" = 2, "V2" = 0)) +
  #coord_flip() +
  theme_bw()+
  facet_wrap(~fct_rev(treatment), ncol=2)+
  theme(axis.text.x = element_text(angle = 90))


#Variability in fever change (baseline to max)
ti.f <- ti.f %>%
  group_by(band_number)%>%
  mutate(fever_peak = max(fever_score),
         fever_high = max(fever_high = max(fever_score[dpi != 0])))

ti.fc <- ti.f %>%
  dplyr::select(band_number, treatment, temp, dpi, fever_score, fever_peak, mass, sex) %>%
  group_by(band_number, treatment, temp, sex) %>%
  summarise(
    baseline   = fever_score[dpi == 0],
    peak       = max(fever_score, na.rm = TRUE),
    high       = max(fever_high = max(fever_score[dpi != 0])),
    magnitude  = peak - baseline,
    mag_high   = high - baseline,
    .groups    = "drop"
  ) %>%
  mutate(group = paste(temp, ifelse(treatment == "Inoculated", "Inoculated", "Control")))

ti.fc.clean <- ti.fc %>%
  filter(!is.na(baseline | peak))

# Count entries by temp and treatment
ti.fc.clean_count <- ti.fc.clean %>%
  group_by(temp, treatment, group) %>%
  summarize(n = length(baseline)) %>%
  ungroup()

lvl <- c("Warm Control", "Warm Inoculated”, “Cold Control", "Cold Inoculated")

# ti.fc.clean <- ti.fc.clean %>%
#   mutate(
#     groups = paste(temp, ifelse(treatment == "Inoculated", "Inoculated", "Control")),
#     groups = str_trim(groups)
#   )

# Variability in mag_high 
  #Doesn't work right becuase there are negative numbers
fever_c_var <- ti.fc.clean %>%
  group_by(group, temp, treatment) %>%
  dplyr::reframe(
    temp = temp,
    band_number = band_number,
    baseline = baseline,
    peak = peak,
    high = high,
    magnitude = magnitude,
    mag_high = mag_high,
    treatment = treatment,
    group = group,
    fever_score = fever_score,
    
    # Calculations for mag_high
    mean_fever = mean(mag_high, na.rm = TRUE),
    median_fever = median(mag_high, na.rm = TRUE),
    bird_cv_fever = calculate_cv(mag_high),
    bird_pv_fever = calculate_pv(mag_high),
    bird_v2_fever = calculate_v2(mag_high),
    
    # Bootstrapping for mag_high
    cv_bootstrap_fever = list(replicate(n_boot, calculate_cv(sample(mag_high, replace = TRUE)))),
    pv_bootstrap_fever = list(replicate(n_boot, calculate_pv(sample(mag_high, replace = TRUE)))),
    v2_bootstrap_fever = list(replicate(n_boot, calculate_v2(sample(mag_high, replace = TRUE))))
  ) %>%
  mutate(
    # Confidence intervals for mag_high
    cv_lower_ci_fever = map_dbl(cv_bootstrap_fever, ~ quantile(.x, 0.025, na.rm = TRUE)),
    cv_upper_ci_fever = map_dbl(cv_bootstrap_fever, ~ quantile(.x, 0.975, na.rm = TRUE)),
    pv_lower_ci_fever = map_dbl(pv_bootstrap_fever, ~ quantile(.x, 0.025, na.rm = TRUE)),
    pv_upper_ci_fever = map_dbl(pv_bootstrap_fever, ~ quantile(.x, 0.975, na.rm = TRUE)),
    v2_lower_ci_fever = map_dbl(v2_bootstrap_fever, ~ quantile(.x, 0.025, na.rm = TRUE)),
    v2_upper_ci_fever = map_dbl(v2_bootstrap_fever, ~ quantile(.x, 0.975, na.rm = TRUE))
  ) %>%
  dplyr::select(-cv_bootstrap_fever, -pv_bootstrap_fever, -v2_bootstrap_fever) %>%
  ungroup()

fever_var_summary <- ti.fc.clean %>%
  group_by(group, temp, treatment) %>%
  reframe(
    across(
      c(baseline, peak, high, magnitude, mag_high),
      list(
        mean = ~mean(.x, na.rm = TRUE),
        cv   = ~calculate_cv(.x),
        pv   = ~calculate_pv(.x),
        v2   = ~calculate_v2(.x)
      ),
      .names = "{.col}_{.fn}"
    ),
    
    # Bootstrapping for each variable + metric
    # baseline_cv_boot = list(replicate(n_boot, calculate_cv(sample(baseline, replace = TRUE)))),
    # baseline_pv_boot = list(replicate(n_boot, calculate_pv(sample(baseline, replace = TRUE)))),
    # baseline_v2_boot = list(replicate(n_boot, calculate_v2(sample(baseline, replace = TRUE)))),
    # 
    # peak_cv_boot = list(replicate(n_boot, calculate_cv(sample(peak, replace = TRUE)))),
    # peak_pv_boot = list(replicate(n_boot, calculate_pv(sample(peak, replace = TRUE)))),
    # peak_v2_boot = list(replicate(n_boot, calculate_v2(sample(peak, replace = TRUE)))),
    # 
    # high_cv_boot = list(replicate(n_boot, calculate_cv(sample(high, replace = TRUE)))),
    # high_pv_boot = list(replicate(n_boot, calculate_pv(sample(high, replace = TRUE)))),
    # high_v2_boot = list(replicate(n_boot, calculate_v2(sample(high, replace = TRUE)))),
    # 
    # magnitude_cv_boot = list(replicate(n_boot, calculate_cv(sample(magnitude, replace = TRUE)))),
    # magnitude_pv_boot = list(replicate(n_boot, calculate_pv(sample(magnitude, replace = TRUE)))),
    # magnitude_v2_boot = list(replicate(n_boot, calculate_v2(sample(magnitude, replace = TRUE)))),
    
    mag_high_cv_boot = list(replicate(n_boot, calculate_cv(sample(mag_high, replace = TRUE)))),
    mag_high_pv_boot = list(replicate(n_boot, calculate_pv(sample(mag_high, replace = TRUE)))),
    mag_high_v2_boot = list(replicate(n_boot, calculate_v2(sample(mag_high, replace = TRUE))))
  ) %>%
  mutate(
    # Confidence intervals for each metric
    across(
      ends_with("_cv_boot"),
      list(
        lower = ~quantile(unlist(.x), 0.025, na.rm = TRUE),
        upper = ~quantile(unlist(.x), 0.975, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    ),
    across(
      ends_with("_pv_boot"),
      list(
        lower = ~quantile(unlist(.x), 0.025, na.rm = TRUE),
        upper = ~quantile(unlist(.x), 0.975, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    ),
    across(
      ends_with("_v2_boot"),
      list(
        lower = ~quantile(unlist(.x), 0.025, na.rm = TRUE),
        upper = ~quantile(unlist(.x), 0.975, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    )
  ) %>%
  dplyr::select(-ends_with("_boot")) %>%
  ungroup()

ggplot(fever_var_summary, aes(x = group, color = group)) +
  geom_point(aes(y = mag_high_pv, shape = "PV"), size = 2, position=position_dodge(width=2)) +
  geom_errorbar(aes(y = mag_high_pv, ymin = mag_high_pv_boot_lower, ymax = mag_high_pv_boot_upper), width = 0.5, position=position_dodge(width=2)) +
  #geom_point(aes(y = mean_bird_cv_fever, shape = "CV"), size = 2, position=position_dodge(width=2)) +
  #geom_point(aes(y = mean_bird_v2_fever, shape = "V2"), size = 2) +
  
  labs(
    y = "Variability",
    x = "Temperature",
    shape = "Metric Type",
    color = "Temperature",
    title = "Variability in Magnitude of Fever Change"
  ) +
  scale_color_manual(values = treat_colors) +
  scale_shape_manual(name = "Variability Metric",
                     values = c("PV" = 16, "CV" = 2, "V2" = 0)) +
  #coord_flip() +
  theme_bw()+
  facet_wrap(~fct_rev(treatment), ncol=2)+
  theme(axis.text.x = element_text(angle = 90))


# fever_c_var <- ti.fc.clean %>%
#   group_by(group, temp, treatment) %>%
#   dplyr::reframe(
#     temp = temp,
#     band_number = band_number,
#     baseline = baseline,
#     peak = peak,
#     high = high,
#     magnitude = magnitude,
#     mag_high = mag_high,
#     treatment = treatment,
#     group = group,
#     fever_score = fever_score,
#     
#     # Calculations for mag_high
#     mean_fever = mean(mag_high, na.rm = TRUE),
#     median_fever = median(mag_high, na.rm = TRUE),
#     bird_cv_fever = calculate_cv(mag_high),
#     bird_pv_fever = calculate_pv(mag_high),
#     bird_v2_fever = calculate_v2(mag_high),
#     
#     # Bootstrapping for mag_high
#     cv_bootstrap_fever = list(replicate(n_boot, calculate_cv(sample(mag_high, replace = TRUE)))),
#     pv_bootstrap_fever = list(replicate(n_boot, calculate_pv(sample(mag_high, replace = TRUE)))),
#     v2_bootstrap_fever = list(replicate(n_boot, calculate_v2(sample(mag_high, replace = TRUE))))
#   ) %>%
#   mutate(
#     # Confidence intervals for mag_high
#     cv_lower_ci_fever = map_dbl(cv_bootstrap_fever, ~ quantile(.x, 0.025, na.rm = TRUE)),
#     cv_upper_ci_fever = map_dbl(cv_bootstrap_fever, ~ quantile(.x, 0.975, na.rm = TRUE)),
#     pv_lower_ci_fever = map_dbl(pv_bootstrap_fever, ~ quantile(.x, 0.025, na.rm = TRUE)),
#     pv_upper_ci_fever = map_dbl(pv_bootstrap_fever, ~ quantile(.x, 0.975, na.rm = TRUE)),
#     v2_lower_ci_fever = map_dbl(v2_bootstrap_fever, ~ quantile(.x, 0.025, na.rm = TRUE)),
#     v2_upper_ci_fever = map_dbl(v2_bootstrap_fever, ~ quantile(.x, 0.975, na.rm = TRUE))
#   ) %>%
#   dplyr::select(-cv_bootstrap_fever, -pv_bootstrap_fever, -v2_bootstrap_fever) %>%
#   ungroup()

ti.f <- ti.f %>%
  dplyr::mutate(sex = dplyr::recode(sex,
                                    "M" = "Male",
                                    "F" = "Female"))
ti.f <- ti.f %>%
  dplyr::mutate(groups = dplyr::recode(groups,
                                       "Cold Control" = "Cold Control",
                                       "Cold Infected" = "Cold Inoculated",
                                       "Warm Control" = "Warm Control",
                                       "Warm Infected" = "Warm Inoculated"))

fever_var <- ti.f %>%
  group_by(dpi, temp, treatment, groups, sex) %>%
  filter(dpi %in% c(-0, 3, 7, 14, 18, 24, 28, 35))%>%
  dplyr::reframe(
    dpi = dpi,
    temp = temp,
    treatment = treatment,
    sex = sex,
    band_number = band_number,
    fever_score = fever_score,
    
    # Calculations for total_eye_score
    mean_fever = mean(fever_score, na.rm = TRUE),
    median_fever = median(fever_score, na.rm = TRUE),
    bird_cv_fever = calculate_cv(fever_score),
    bird_pv_fever = calculate_pv(fever_score),
    #bird_v2_fever = calculate_v2(fever_score),
    
    # Bootstrapping for total_eye_score
    cv_bootstrap_fever = list(replicate(n_boot, calculate_cv(sample(fever_score, replace = TRUE)))),
    pv_bootstrap_fever = list(replicate(n_boot, calculate_pv(sample(fever_score, replace = TRUE)))),
    v2_bootstrap_fever = list(replicate(n_boot, calculate_v2(sample(fever_score, replace = TRUE))))
  ) %>%
  mutate(
    # Confidence intervals for total_eye_score
    cv_lower_ci_fever = map_dbl(cv_bootstrap_fever, ~ quantile(.x, 0.025, na.rm = TRUE)),
    cv_upper_ci_fever = map_dbl(cv_bootstrap_fever, ~ quantile(.x, 0.975, na.rm = TRUE)),
    pv_lower_ci_fever = map_dbl(pv_bootstrap_fever, ~ quantile(.x, 0.025, na.rm = TRUE)),
    pv_upper_ci_fever = map_dbl(pv_bootstrap_fever, ~ quantile(.x, 0.975, na.rm = TRUE)),
    v2_lower_ci_fever = map_dbl(v2_bootstrap_fever, ~ quantile(.x, 0.025, na.rm = TRUE)),
    v2_upper_ci_fever = map_dbl(v2_bootstrap_fever, ~ quantile(.x, 0.975, na.rm = TRUE))
  ) %>%
  dplyr::select(-cv_bootstrap_fever, -pv_bootstrap_fever, -v2_bootstrap_fever) %>%
  ungroup()

# Summary statistics for total_eye_score
summary_fever_tibble <- fever_var %>%
  group_by(dpi, temp, treatment, groups, sex) %>%
  dplyr::reframe(
    mean_fever = mean(mean_fever, na.rm = TRUE),
    median_fever = mean(median_fever, na.rm = TRUE),
    mean_bird_cv_fever = mean(bird_cv_fever, na.rm = TRUE),
    mean_bird_pv_fever = mean(bird_pv_fever, na.rm = TRUE),
    #mean_bird_v2_fever = mean(bird_v2_fever, na.rm = TRUE),
    mean_cv_lower_ci_fever = mean(cv_lower_ci_fever, na.rm = TRUE),
    mean_cv_upper_ci_fever = mean(cv_upper_ci_fever, na.rm = TRUE),
    mean_pv_lower_ci_fever = mean(pv_lower_ci_fever, na.rm = TRUE),
    mean_pv_upper_ci_fever = mean(pv_upper_ci_fever, na.rm = TRUE),
    #mean_v2_lower_ci_fever = mean(v2_lower_ci_fever, na.rm = TRUE),
    #mean_v2_upper_ci_fever = mean(v2_upper_ci_fever, na.rm = TRUE)
  )

# Plot variability for total_eye_score
dodge = position_dodge(width=2)

ggplot(summary_fever_tibble, aes(x = dpi, color = groups)) +
  geom_point(aes(y = mean_bird_pv_fever, color=groups), size = 2, position = dodge) +
  geom_errorbar(aes(y = mean_bird_pv_fever, ymin = mean_pv_lower_ci_fever, ymax = mean_pv_upper_ci_fever, color=groups), 
                width = 0.1, position = dodge) +
   #geom_point(aes(y = mean_bird_cv_fever, shape = "CV"), size = 2, shape=1) +
  # geom_point(aes(y = mean_bird_v2_eye_score, shape = "V2"), size = 2) +
  
  labs(
    y = "Variability in Fever Score (PV)",
    x = "Days Post-Inoculation",
    shape = "Metric Type",
    color = "Temperature",
    #title = "Variability in Total Eye Score by DPI and Temperature"
  ) +
  #scale_y_continuous(limits=c(0,1))+
  scale_color_manual(values = treat_colors) +
  facet_grid(~fct_rev(temp)~sex)+
  theme(strip.text = element_text(size=12))
# scale_shape_manual(name = "Variability Metric",
#                    values = c("PV" = 16, "CV" = 2, "V2" = 0)) +
#facet_wrap(~groups)

leveneTest(fever_score ~ treatment, data=ti.f, center=median)
leveneTest(fever_score ~ temp, data=ti.f, center=median)
leveneTest(fever_score ~ sex, data=ti.f, center=median)

#
leveneTest(fever_score ~ sex, data=ti.f %>% filter(treatment == "Inoculated"), center=median)
leveneTest(fever_score ~ sex, data=ti.f %>% filter(treatment == "Sham"), center=median)

ggplot(ti.f, aes(x= sex, y=fever_score, color= treatment))+
  geom_jitter(width=0.2)+
  #scale_color_manual(values=treat_colors)+
  facet_grid(~treatment)+
  theme(axis.text.x = element_text(angle=45, hjust=1))+
  labs(x="Treatment Groups", y="Fever Score (All Days Combined B-H Test)", color="Treatment Groups")

leveneTest(fever_score ~ groups, data=ti.f, center=median)


ggplot(ti.f, aes(x= groups, y=fever_score, color= groups))+
  geom_jitter(width=0.2)+
  scale_color_manual(values=treat_colors)+
  facet_grid(~sex~treatment)+
  theme(axis.text.x = element_text(angle=45, hjust=1))+
  labs(x="Treatment Groups", y="Fever Score (All Days Combined B-H Test)", color="Treatment Groups")



####Mass####
source("r_scripts/dataCleaning_TI22.R")
ti %>%
  dplyr::select(dpi, treatment, temp, groups)%>%
  tbl_summary(
    by="dpi"
  )%>%
  modify_header(
    label ~ "**All Birds**"
  )

ti.m <- ti %>%
  filter(dpi %in% c(-12, 3, 14, 21, 28, 35),
         !is.na(mass))

#Antibody analysis sample sizes; see dataCleaning_antibody.R for removal breakdown
ti.m %>%
  dplyr::select(dpi, treatment, temp, groups)%>%
  tbl_summary(
    by="dpi"
  )%>%
  modify_header(
    label ~ "**Birds for Mass**"
  )

ti.m <- ti.m %>% 
  mutate(groups = factor(groups,
                         levels = c("Warm Control", "Warm Inoculated",
                                    "Cold Control", "Cold Inoculated")))

unique(ti.m$dpi)
ti.m$dpi <- as.factor(ti.m$dpi)

#Set reference categories "Warm" and "Sham"
ti.m$temp <- relevel(as.factor(ti.m$temp), ref = "Warm")
ti.m$treatment <- relevel(as.factor(ti.m$treatment), ref = "Control")
ti.m$dpi.f <- as.factor(ti.m$dpi)

# Variability analysis for mass
mass_var <- ti.m %>%
  group_by(dpi, groups, ever_diseased) %>%
  dplyr::reframe(
    temp = temp,
    band_number = band_number,
    mass = mass,
    treatment = treatment,
    groups = groups,
    ever_diseased = ever_diseased,
    
    # Calculations for mass
    mean_mass = mean(mass, na.rm = TRUE),
    median_mass = median(mass, na.rm = TRUE),
    bird_cv_mass = calculate_cv(mass),
    bird_pv_mass = calculate_pv(mass),
    bird_v2_mass = calculate_v2(mass),
    
    # Bootstrapping for mass
    cv_bootstrap_mass = list(replicate(n_boot, calculate_cv(sample(mass, replace = TRUE)))),
    pv_bootstrap_mass = list(replicate(n_boot, calculate_pv(sample(mass, replace = TRUE)))),
    v2_bootstrap_mass = list(replicate(n_boot, calculate_v2(sample(mass, replace = TRUE))))
  ) %>%
  mutate(
    # Confidence intervals for mass
    cv_lower_ci_mass = map_dbl(cv_bootstrap_mass, ~ quantile(.x, 0.025, na.rm = TRUE)),
    cv_upper_ci_mass = map_dbl(cv_bootstrap_mass, ~ quantile(.x, 0.975, na.rm = TRUE)),
    pv_lower_ci_mass = map_dbl(pv_bootstrap_mass, ~ quantile(.x, 0.025, na.rm = TRUE)),
    pv_upper_ci_mass = map_dbl(pv_bootstrap_mass, ~ quantile(.x, 0.975, na.rm = TRUE)),
    v2_lower_ci_mass = map_dbl(v2_bootstrap_mass, ~ quantile(.x, 0.025, na.rm = TRUE)),
    v2_upper_ci_mass = map_dbl(v2_bootstrap_mass, ~ quantile(.x, 0.975, na.rm = TRUE))
  ) %>%
  dplyr::select(-cv_bootstrap_mass, -pv_bootstrap_mass, -v2_bootstrap_mass) %>%
  ungroup()

# Summary statistics for mass
summary_mass_tibble <- mass_var %>%
  group_by(dpi, groups, temp, ever_diseased) %>%
  dplyr::reframe(
    mean_mass = mean(mean_mass, na.rm = TRUE),
    median_mass = mean(median_mass, na.rm = TRUE),
    mean_bird_cv_mass = mean(bird_cv_mass, na.rm = TRUE),
    mean_bird_pv_mass = mean(bird_pv_mass, na.rm = TRUE),
    mean_bird_v2_mass = mean(bird_v2_mass, na.rm = TRUE),
    mean_cv_lower_ci_mass = mean(cv_lower_ci_mass, na.rm = TRUE),
    mean_cv_upper_ci_mass = mean(cv_upper_ci_mass, na.rm = TRUE),
    mean_pv_lower_ci_mass = mean(pv_lower_ci_mass, na.rm = TRUE),
    mean_pv_upper_ci_mass = mean(pv_upper_ci_mass, na.rm = TRUE),
    mean_v2_lower_ci_mass = mean(v2_lower_ci_mass, na.rm = TRUE),
    mean_v2_upper_ci_mass = mean(v2_upper_ci_mass, na.rm = TRUE),
    n_individuals = n_distinct(band_number) #calculate number of individuals in each group
  )

dodge <- position_dodge(width=0.75)

var.mass <- ggplot(summary_mass_tibble, aes(x = dpi, color = groups)) +
  geom_errorbar(aes(x = dpi, y = mean_bird_pv_mass, ymin = mean_pv_lower_ci_mass, ymax = mean_pv_upper_ci_mass, group = groups),
                width = 0, color="black", position=dodge) +
  geom_point(aes(y = mean_bird_pv_mass, shape = "PV", group = groups),
             size = 3, shape=1, color="black", stroke=1.5, position=dodge) +
  geom_point(aes(y = mean_bird_pv_mass, shape = "PV", group = groups),
             size = 3, shape=16, position=dodge) +
  geom_line(aes(x = dpi, y = mean_bird_pv_mass, group = groups),
            position=dodge)+
  
  #geom_point(aes(y = mean_bird_cv_mass, shape = "CV", group=dpi), size = 2, position = dodge) +
  
  labs(
    y = "Variability (PV)",
    x = "Days Post Inoculation",
    shape = "Metric Type",
    color = "Treatment Group"
    #title = "Variability in Mass"
  ) +
  scale_color_manual(values = c(treat_colors))+
  facet_grid(~temp~ever_diseased)+
  theme(
    strip.text = element_text(size=12)
  )
var.mass


bf_mass_treat <- leveneTest(mass ~ treatment, data = ti.m, center = median)
bf_mass_temp <- leveneTest(mass ~ temp, data = ti.m, center = median)
bf_mass_groups <- leveneTest(mass ~ groups, data = ti.m, center = median)
bf_mass_dis <- leveneTest(mass ~ as.factor(ever_diseased), data = ti.m, center = median)

bf_mass_dis_cold <- leveneTest(mass ~ as.factor(ever_diseased), data = ti.m %>% filter(temp == "Cold"), center = median)
bf_mass_dis_warm <- leveneTest(mass ~ as.factor(ever_diseased), data = ti.m %>% filter(temp == "Warm"), center = median)

bf_mass_treat_cold <- leveneTest(mass ~ treatment, data = ti.m %>% filter(temp == "Cold" & dpi == 21), center = median)
bf_mass_treat_warm <- leveneTest(mass ~ treatment, data = ti.m %>% filter(temp == "Warm" & dpi == 21), center = median)

