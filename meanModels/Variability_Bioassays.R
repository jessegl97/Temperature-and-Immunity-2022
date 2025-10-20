####Variability in Bioassays###
rm(list=ls())

####read in + format data####
setwd('/Users/jesse/Documents/GitHub/Temperature-and-Immunity-2022/meanModels/')
library(ggplot2)
library(tidyverse)
library(lme4)
library(glmmTMB)
library(effects)
library(DHARMa)
library(gtsummary)

#set theme
theme_set(theme_minimal(base_size=12))
temp_cols <- c("#669BBC", "#C1121F")

#Antibodies
source("dataCleaning_antibody.R")

ti.pi %>%
  filter(dpi ==28)%>%
  dplyr::select(treatment, temp, inf_temp)%>%
  tbl_summary(
  )%>%
  modify_header(
    label ~ "**All Birds**"
  )

elisa_counts <- ti.pi %>%
  group_by(dpi, temp, band_number) %>%
  summarise(counts = sum(!is.na(elisa_od)), .groups="drop")

#missing elisa data from 10 samples:
  #DPI 9: 2684, 2893, 2632, 2842, 2847
  #DPI 28: 2684, 2893, 2846, 2869, 2926

ti.pi$dpi <- as.factor(ti.pi$dpi)
ti.pi$band_number <- as.factor(ti.pi$band_number)

#2926 has high antibody levels but no plasma from dpi 28 which sucks
ti.pi %>%
  filter(elisa_od>0.08) %>%
  dplyr::select(dpi, elisa_od, band_number, treatment)

ggplot(ti.pi, aes(x=dpi, y=elisa_od, color=temp))+
  #geom_jitter(height=0, width=0.15, alpha=0.5)+
  geom_point(alpha=0.5)+
  geom_line(aes(group=band_number), alpha=0.5)+
  scale_color_manual(values = c("Cold" = "#669BBC", "Warm" = "#C1121F"))

#I guess omit 2926 because it doesn't have ab levels dpi 28?
#for analysis with 2926
ti.pi.full <- ti.pi

ti.pi <- ti.pi %>%
  filter(band_number != 2926)
#mean model
lm1<-glm(elisa_od~treatment+temp+dpi + treatment:dpi,
         data=ti.pi, 
         family=Gamma(link = "inverse"))
summary(lm1)
simulateResiduals(lm1, plot=T)

lm1.f <- glm(elisa_od~treatment+temp+dpi + treatment:dpi,
             data=ti.pi.full, 
             family=Gamma(link = "inverse"))

summary(lm1.f)

#no interaction
lm1.5<- glm(elisa_od ~ treatment+temp + dpi, data=ti.pi, family=Gamma(link="inverse"))
summary(lm1.5)
plot(allEffects(lm1.5))
simulateResiduals(lm1.5, plot=T)

lm1.5.f<- glm(elisa_od ~ treatment+temp + dpi, data=ti.pi.full, family=Gamma(link="inverse"))
summary(lm1.5.f)

#antibody levels were elevated on dpi 28 in infected birds, but there was no effect of temperature (p=0.5)
#including 2926 does not change this result.

mod <- lm1.5
dat.new=expand.grid(elisa_od=unique(ti.pi$elisa_od),
                    temp=unique(ti.pi$temp),
                    dpi = unique(ti.pi$dpi),
                    treatment= unique(ti.pi$treatment))

dat.new$yhat=predict(mod, type="response", newdata = dat.new)
#prediction intervals
preds = predict(mod, type = "link", newdata = dat.new, se.fit =T)
#bind se's and fitted points
dat.new = cbind(dat.new, preds)
#inverse link function
ilink <- family(mod)$linkinv
#back transform CIs
dat.new <- transform(dat.new,
                     Fitted = ilink(fit),
                     Upper = ilink(fit + (2*se.fit)),
                     Lower = ilink(fit - (2*se.fit)))

ggplot(ti.pi, aes(x=dpi, y=elisa_od, color=temp))+
  geom_hline(yintercept = 0.061, alpha=2, lty="dashed")+
  geom_point(alpha=0.5)+
  geom_line(aes(group=band_number), alpha=2)+
  geom_point(data=dat.new, aes(x=dpi, y=yhat, color=temp), size=2, shape=15)+
  geom_line(data=dat.new, aes(x=dpi, y=yhat, group = temp, color=temp), size=1)+
  geom_errorbar(data=dat.new, aes(x=dpi, y=yhat, ymin=Lower, ymax=Upper, color=temp), width = 0.1)+
  scale_color_manual(values = c("Cold" = "#669BBC", "Warm" = "#C1121F"))+
  labs(x="Days Post Inoculation", y="Antibody Levels", color="Temperature")+
  facet_wrap(~treatment)

dat.new_summary <- dat.new %>%
  group_by(temp, dpi, treatment) %>%
  summarize(
    yhat = mean(yhat, na.rm = TRUE),  # Mean of predicted values
    Upper = mean(Upper, na.rm = TRUE),  # Mean of Upper CI
    Lower = mean(Lower, na.rm = TRUE))%>%  # Mean of Lower CI
    ungroup()

#write.csv(dat.new_summary, "/Users/jesse/Documents/GitHub/Temperature-and-Immunity-2022/Variability/anibodies_model_preds.csv", row.names = FALSE)
#Within-group variability - just look dpi 28

#Are antibody responses more variable in cold temperatures?
#install.packages("remotes")
#remotes::install_github("T-Engel/CValternatives")
library(CValternatives)
# Define functions for calculating CV, PV, and V2
calculate_cv <- function(data) {
  return(sd(data) / mean(data))
}

calculate_pv <- function(data) {
  return(PV(data))
}

calculate_v2 <- function(data) {
  mean_x <- mean(data)
  sd_x <- sd(data)
  v2 <- sd_x^2 / (sd_x^2 + mean_x^2)
  return(v2)
}

n_boot <- 1000

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
  group_by(dpi, temp, treatment) %>%
  dplyr::reframe(
    temp = temp,
    band_number = band_number,
    elisa_od = elisa_od,
    treatment = treatment,
    
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
  group_by(dpi, temp, treatment) %>%
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
    mean_v2_upper_ci_elisa = mean(v2_upper_ci_elisa, na.rm = TRUE)
  )

ggplot(summary_elisa_tibble, aes(x = interaction(temp, dpi, treatment), color = temp)) +
  geom_point(aes(y = mean_bird_pv_elisa, shape = "PV"), size = 2) +
  geom_errorbar(aes(y = mean_bird_pv_elisa, ymin = mean_pv_lower_ci_elisa, ymax = mean_pv_upper_ci_elisa), width = 0) +
  geom_point(aes(y = mean_bird_cv_elisa, shape = "CV"), size = 2) +
  geom_point(aes(y = mean_bird_v2_elisa, shape = "V2"), size = 2) +
  
  
  labs(
    y = "Variability",
    x = "Temperature",
    shape = "Metric Type",
    color = "Temperature",
    title = "Variability in Antibody Levels DPI 28"
  ) +
  scale_color_manual(values = c("#669BBC", "#C1121F"))+
  scale_shape_manual(name = "Variability Metric",
                     values = c("PV" = 16, "CV" = 2, "V2" = 0)) +
  theme_minimal()+
  theme(
    axis.text.x = element_text(angle = 45, hjust=1)
  )

#write.csv(summary_elisa_tibble, "/Users/jesse/Documents/GitHub/Temperature-and-Immunity-2022/Variability/Elisa_Variability.csv")

####Phagocytosis####

source("dataCleaning_phago.R")

ti.p %>%
  dplyr::select(dpi, treatment, temp, inf_temp)%>%
  tbl_summary(
  )%>%
  modify_header(
    label ~ "**All Birds**"
  )

phago_counts <- ti.p %>%
  group_by(dpi, temp, band_number) %>%
  summarise(counts = sum(!is.na(phago_score)), .groups="drop")

## Binomial model and the
#    weights represent the total number of trials, since the response is a
#    pre-calculated proportion. The interaction between temp and treatment was
#    not significant.

glm1 <- glmer(phago_score~temp+treatment + (1|band_number), 
              weights=wbc_total+phago_total, 
              data=ti.p, family="binomial")

simulateResiduals(glm1, plot = T)

summary(glm1)

mod <- glm1
dat.new=expand.grid(phago_score=unique(ti.p$phago_score),
                    temp=unique(ti.p$temp),
                    treatment= unique(ti.p$treatment),
                    band_number = unique(ti.p$band_number))

dat.new$yhat=predict(mod, type="response", newdata = dat.new, re.form = NA)
# #prediction intervals
preds = predict(mod, type = "link", newdata = dat.new, se.fit =T, re.form = NA)
# #bind se's and fitted points
 dat.new = cbind(dat.new, preds)
# #inverse link function
 ilink <- family(mod)$linkinv
# #back transform CIs
dat.new <- transform(dat.new,
                     Fitted = ilink(fit),
                     Upper = ilink(fit + (2*se.fit)),
                     Lower = ilink(fit - (2*se.fit)))

ggplot(ti.p, aes(x=interaction(temp,treatment), y=phago_score, color=temp))+
  geom_jitter(height=0, width=0.2, alpha=0.5)+
  scale_color_manual(values = temp_cols)+
  geom_point(data=dat.new, aes(x=interaction(temp,treatment), y=yhat, color=temp), size=2)+
  geom_errorbar(data=dat.new, aes(x=interaction(temp, treatment), ymin=Lower, ymax=Upper, y=yhat), width=0, alpha=0.75)

dat.new_summary <- dat.new %>%
  group_by(temp, treatment) %>%
  summarize(
    yhat = mean(yhat, na.rm = TRUE),  # Mean of predicted values
    Upper = mean(Upper, na.rm = TRUE),  # Mean of Upper CI
    Lower = mean(Lower, na.rm = TRUE))%>%  # Mean of Lower CI
  ungroup()

#write.csv(dat.new_summary, "/Users/jesse/Documents/GitHub/Temperature-and-Immunity-2022/Variability/phago_model_preds.csv", row.names = FALSE)

#Variability
# Remove rows with NA values in phago_score
ti.p.clean <- ti.p %>%
  filter(!is.na(phago_score))

# Count entries for dpi == 28 by temp and treatment
ti.p.clean.count <- ti.p.clean %>%
  group_by(temp, treatment) %>%
  summarize(n = length(phago_score))%>%
  ungroup()

# Variability analysis for phago_score
phago_var <- ti.p.clean %>%
  group_by(temp, treatment) %>%
  dplyr::reframe(
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
  group_by(temp, treatment) %>%
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
ggplot(summary_phago_tibble, aes(x = interaction(temp, treatment), color = temp)) +
  geom_point(aes(y = mean_bird_pv_phago, shape = "PV"), size = 2) +
  geom_errorbar(aes(y = mean_bird_pv_phago, ymin = mean_pv_lower_ci_phago, ymax = mean_pv_upper_ci_phago), width = 0) +
  geom_point(aes(y = mean_bird_cv_phago, shape = "CV"), size = 2) +
  geom_point(aes(y = mean_bird_v2_phago, shape = "V2"), size = 2) +
  
  labs(
    y = "Variability",
    x = "Temperature",
    shape = "Metric Type",
    color = "Temperature",
    title = "Variability in Phagocytosis"
  ) +
  scale_color_manual(values = c("#669BBC", "#C1121F")) +
  scale_shape_manual(name = "Variability Metric",
                     values = c("PV" = 16, "CV" = 2, "V2" = 0)) +
  #coord_flip() +
  theme_minimal()

# Write the summarized tibble to a CSV file
#write.csv(summary_phago_tibble, "/Users/jesse/Documents/GitHub/Temperature-and-Immunity-2022/Variability/Phagocytosis_Variability.csv")

####Fever####
source("dataCleaning_fever.R")

ti.f %>%
  dplyr::select(dpi, treatment, temp, groups)%>%
  tbl_summary(
    by="dpi"
  )%>%
  modify_header(
    label ~ "**All Birds**"
  )

fever_counts <- ti.f %>%
  group_by(dpi, groups) %>%
  summarise(counts = sum(!is.na(fever_score)), .groups="drop")

#model
lm2 <- glmmTMB(fever_change ~ temp+treatment*as.factor(dpi) + 
                   (1|band_number),
                 data=ti.f)
summary(lm2)
simulateResiduals(lm2, plot=T)

library(emmeans)
emresults<- emmeans(lm2, ~treatment*dpi)
summary(emresults)


ti.fn <- ti.f %>%
  mutate(
    norm_fever = case_when(
      temp == "Cold" ~ fever_score + 2,  # Add 2 degrees for Cold birds
      temp == "Warm" ~ fever_score - 2,  # Subtract 2 degrees for Warm birds
    )
  )

ggplot(ti.fn, aes(x=dpi, y=norm_fever, color=temp))+
  geom_jitter(height=0, width=0.2, alpha=0.1)+
  stat_summary(aes(group = groups, x=dpi, y=norm_fever, color=temp, shape=groups), fun="mean", geom="point", size=2)+
  stat_summary(aes(group = groups, x=dpi, y=norm_fever, color=temp, shape=groups), fun="mean", geom="line")+
  scale_color_manual(values=c(temp_cols))+
  facet_wrap(~treatment)+
  theme_bw()

lm2.5 <- glmmTMB(norm_fever ~ treatment+temp+as.factor(dpi)+ (1|band_number), data=ti.fn)
summary(lm2.5)

lm2.75 <- glmmTMB(norm_fever ~ temp + treatment + as.factor(dpi) + 
                 treatment:as.factor(dpi) + temp:as.factor(dpi) + ar1(as.factor(dpi) + 0|band_number),
               data=ti.fn)
summary(lm2.75)
#Treatment and DPI, but not temperature affect fever change.

# mod <- lm2
# dat.new=expand.grid(fever_change=unique(ti.f$fever_change),
#                     temp=unique(ti.f$temp),
#                     dpi = unique(ti.f$dpi),
#                     treatment= unique(ti.f$treatment),
#                     band_number = unique(ti.f$band_number))
# 
# dat.new$yhat=predict(mod, type="response", newdata = dat.new)
# #prediction intervals
# preds = predict(mod, type = "link", newdata = dat.new, se.fit =T)
# #bind se's and fitted points
# dat.new = cbind(dat.new, preds)
# #inverse link function
# ilink <- family(mod)$linkinv
# #back transform CIs
# dat.new <- transform(dat.new,
#                      Fitted = ilink(fit),
#                      Upper = ilink(fit + (2*se.fit)),
#                      Lower = ilink(fit - (2*se.fit)))
# 
ggplot(ti.f, aes(x=dpi, y=fever_change, color=temp))+
  geom_jitter(height=0, width=0.2)+
  geom_line(aes(group=band_number), alpha=0.2)+
  facet_grid(~temp~treatment)+
  scale_color_manual(values=c(temp_cols))
  geom_point(data=dat.new, x=dpi, y=yhat, color=groups)

  
  
  
#variability in fever_score (raw temperature) by temperature
  # Remove rows with NA values in fever_score
  ti.f.clean <- ti.f %>%
    filter(!is.na(fever_score))
  
  # Count entries by temp and treatment
  ti.f.clean.count <- ti.f.clean %>%
    group_by(dpi, temp, treatment) %>%
    summarize(n = length(fever_score)) %>%
    ungroup()
  
  # Variability analysis for fever_score
  fever_var <- ti.f.clean %>%
    group_by(dpi, temp, treatment) %>%
    dplyr::reframe(
      temp = temp,
      band_number = band_number,
      fever_score = fever_score,
      treatment = treatment,
      
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
    group_by(dpi, temp, treatment) %>%
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
    )%>%
    # Add a new column 'norm_fever' based on the temperature adjustment
    mutate(
      norm_fever = case_when(
        temp == "Cold" ~ mean_fever + 2,  # Add 2 degrees for Cold birds
        temp == "Warm" ~ mean_fever - 2,  # Subtract 2 degrees for Warm birds
        TRUE ~ mean_fever                # Default (no change)
      )
    )
  
  # Plot variability for fever_score
  ggplot(summary_fever_tibble, aes(x = interaction(temp, dpi, treatment), color = temp)) +
    geom_point(aes(y = mean_bird_pv_fever, shape = "PV"), size = 2) +
    geom_errorbar(aes(y = mean_bird_pv_fever, ymin = mean_pv_lower_ci_fever, ymax = mean_pv_upper_ci_fever), width = 0) +
    geom_point(aes(y = mean_bird_cv_fever, shape = "CV"), size = 2) +
    #geom_point(aes(y = mean_bird_v2_fever, shape = "V2"), size = 2) +
    
    labs(
      y = "Variability",
      x = "Temperature",
      shape = "Metric Type",
      color = "Temperature",
      title = "Variability in Fever Score"
    ) +
    scale_color_manual(values = c("#669BBC", "#C1121F")) +
    scale_shape_manual(name = "Variability Metric",
                       values = c("PV" = 16, "CV" = 2, "V2" = 0)) +
    #coord_flip() +
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 90))

  # Plot mean for fever_score
  ggplot(summary_fever_tibble, aes(x = interaction(temp, dpi, treatment), color = temp)) +
    geom_point(aes(y = norm_fever), size = 2) +
    labs(
      y = "Mean",
      x = "Temperature",
      color = "Temperature",
      title = "Mean Fever Score"
    ) +
    scale_color_manual(values = c("#669BBC", "#C1121F")) +
    #coord_flip() +
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 90))
  
  
  # Write the summarized tibble to a CSV file
  #write.csv(summary_fever_tibble, "/Users/jesse/Documents/GitHub/Temperature-and-Immunity-2022/Variability/Fever_score_Variability.csv")
  
  
ti.fm <- ti.f %>%
  group_by(dpi, groups, temp, treatment)%>%
  reframe(max_fc = max(fever_change),
            min_fc = abs(min(fever_change)),
          max_fever_score = max(fever_score),
            diseased = diseased,
            ever_diseased = ever_diseased,
            infected = infected,
            ever_infected = ever_infected)

ggplot(ti.fm, aes(x=groups, color=temp))+
  geom_point(aes(y=max_fc))


ti.fm <- ti.f %>%
    group_by(band_number) %>%
    mutate(shifted_fever_change = fever_change + min(fever_change, na.rm = TRUE)) %>%
    reframe(
      max_fc = max(fever_change),
      max_shifted_fc = max(shifted_fever_change),
      max_fever_score = max(fever_score),
      diseased = diseased,
      ever_diseased=ever_diseased,
      temp=temp,
      treatment=treatment,
      groups=groups,
      dpi=dpi,
      band_number = band_number)
  
ggplot(ti.fm, aes(x=groups, color=temp))+
  geom_point(aes(y=max_fever_score))


#Does disease predict fever change?
lm3 <- glmmTMB(max_shifted_fc ~ diseased + temp, data=ti.fm)
summary(lm3)
simulateResiduals(lm3, plot=T)
#Birds with active pathology had higher fever scores

#That's kinda problematic though above and below 0 matters.
#
lm3.5 <- glmmTMB(fever_score ~ treatment + temp + dpi + (1|band_number), data=ti.f)
summary(lm3.5) 
simulateResiduals(lm3.5, plot=T)

#so what about peak fever score?
lm3.75 <- glmmTMB(norm_fever ~ treatment + temp, data=ti.fn %>% filter(dpi == 14))
summary(lm3.75)
simulateResiduals(lm3.75, plot=T)


mod <- lm3.75
dat.new=expand.grid(norm_fever=unique(ti.fn$norm_fever),
                    temp=unique(ti.fn$temp),
                    treatment= unique(ti.fn$treatment))

dat.new$yhat=predict(mod, type="response", newdata = dat.new, re.form = NA)
# #prediction intervals
preds = predict(mod, type = "link", newdata = dat.new, se.fit =T, re.form = NA)
# #bind se's and fitted points
dat.new = cbind(dat.new, preds)
# #inverse link function
ilink <- family(mod)$linkinv
# #back transform CIs
dat.new <- transform(dat.new,
                     Fitted = ilink(fit),
                     Upper = ilink(fit + (2*se.fit)),
                     Lower = ilink(fit - (2*se.fit)))

ggplot(ti.fn %>% filter(dpi ==14), aes(x=interaction(temp,treatment), y=norm_fever, color=temp))+
  geom_jitter(height=0, width=0.2, alpha=0.5)+
  scale_color_manual(values = temp_cols)+
  geom_point(data=dat.new, aes(x=interaction(temp,treatment), y=yhat, color=temp), size=2)+
  geom_errorbar(data=dat.new, aes(x=interaction(temp, treatment), ymin=Lower, ymax=Upper, y=yhat), width=0, alpha=0.75)

dat.new_summary <- dat.new %>%
  group_by(temp, treatment) %>%
  summarize(
    yhat = mean(yhat, na.rm = TRUE),  # Mean of predicted values
    Upper = mean(Upper, na.rm = TRUE),  # Mean of Upper CI
    Lower = mean(Lower, na.rm = TRUE))%>%  # Mean of Lower CI
  ungroup()

#write.csv(dat.new_summary, "/Users/jesse/Documents/GitHub/Temperature-and-Immunity-2022/Variability/max_fever_score_model_preds.csv", row.names = FALSE)


####eye score####
source("dataCleaning_eyeScore.R")

ti.mg %>%
  dplyr::select(dpi, treatment, temp, groups)%>%
  tbl_summary(
    by="dpi"
  )%>%
  modify_header(
    label ~ "**All Birds**"
  )

eyescore_counts <- ti.mg %>%
  group_by(dpi, groups) %>%
  summarise(counts = sum(!is.na(total_eye_score)), .groups="drop")

ti.mg$band_number <- as.factor(ti.mg$band_number)
#model
lm1 <- glmmTMB(tes ~ temp*dpi + (1|band_number),
               data=ti.mg, 
               ziformula = ~ temp,
               family = ziGamma(link = "log"))
#zero-inflated model allowing for temperature to dictate the probability of zero.
#for the non-zero birds, we used a gamma model with a log link function, 
#including a two way interaction between temp and dpi and all lower order effects

ggplot(ti.mg, aes(x=dpi, y=total_eye_score, color=temp))+
  geom_point(alpha=0.75)+
  geom_line(aes(group = band_number), alpha=0.5)+
  stat_summary(data=ti.mg %>% filter(ever_infected == 1), aes(x=dpi, y=total_eye_score, group = ever_infected, color=temp),
               geom="line", fun="mean", size=1.5)+
  scale_color_manual(values=c(temp_cols))

#variability
# Variability analysis for total_eye_score
#add small constant
ti.mg$total_eye_score1 <- ti.mg$total_eye_score+0.001
#All Birds
  #No eye scores dpi -12, 3, 35 so get rid 
nboot=1000

eye_score_var <- ti.mg %>%
  group_by(dpi, temp, treatment) %>%
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
  group_by(dpi, temp, treatment) %>%
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
ggplot(summary_eye_score_tibble, aes(x = interaction(temp, dpi), color = temp)) +
  geom_point(aes(y = mean_bird_pv_eye_score, shape = "PV"), size = 2) +
  geom_errorbar(aes(y = mean_bird_pv_eye_score, ymin = mean_pv_lower_ci_eye_score, ymax = mean_pv_upper_ci_eye_score), width = 0) +
  geom_point(aes(y = mean_bird_cv_eye_score, shape = "CV"), size = 2) +
  geom_point(aes(y = mean_bird_v2_eye_score, shape = "V2"), size = 2) +
  
  labs(
    y = "Variability",
    x = "DPI, Temperature, and Treatment",
    shape = "Metric Type",
    color = "Temperature",
    title = "Variability in Total Eye Score by DPI and Temperature"
  ) +
  scale_color_manual(values = c("#669BBC", "#C1121F")) +
  scale_shape_manual(name = "Variability Metric",
                     values = c("PV" = 16, "CV" = 2, "V2" = 0)) +
  #coord_flip() +
  theme_minimal()


#Write the summarized tibble to a CSV file
#write.csv(summary_eye_score_tibble, "/Users/jesse/Documents/GitHub/Temperature-and-Immunity-2022/Variability/Eye_Score_Inoculated_Variability.csv")

#Compare inoculated, infected, and diseased
esc<- read.csv("~/Documents/GitHub/Temperature-and-Immunity-2022/Variability/Eye Score/Eye_Score_Compare_Variability.csv")
esc
#Compare variability in inoculated, infected, and diseased on day 14
ggplot(esc %>% filter(dpi == 14), aes(x=status, y=mean_bird_pv_eye_score, color=temp))+
  geom_point(aes(shape="PV"), position=position_dodge(width=0.2))+
  geom_point(aes(x=status, y=mean_bird_cv_eye_score, color=temp, shape="CV"),
             position=position_dodge(width=0.1))+
  geom_point(aes(x=status, y=mean_bird_v2_eye_score, color=temp, shape="V2"),
             position=position_dodge(width=0.35))+
  geom_errorbar(aes(x=status, y=mean_bird_pv_eye_score,
                    ymin=mean_pv_lower_ci_eye_score, ymax=mean_pv_upper_ci_eye_score),
                position=position_dodge(width=0.2), width=0.01)+
  geom_errorbar(aes(x=status, y=mean_bird_cv_eye_score,
                    ymin=mean_cv_lower_ci_eye_score, ymax=mean_cv_upper_ci_eye_score),
                position=position_dodge(width=0.1), width=0.01, lty="dashed")+
  geom_errorbar(aes(x=status, y=mean_bird_v2_eye_score,
                    ymin=mean_v2_lower_ci_eye_score, ymax=mean_v2_upper_ci_eye_score),
                position=position_dodge(width=0.35), width=0.01, lty="dotted")+
  labs(x="Birds Included", y="Variability", color="Temperature", title="Comparison of Variability in Peak Eye Score")+
  scale_color_manual(values = c("#669BBC", "#C1121F")) +
  scale_shape_manual(
    name="Variability Metric",
    values=c("CV" = 1, "PV"=16, "V2" = 2)
  )

#compare means in diseased infected inoculated
lm4 <- glm(mean_eye_score ~ temp + status, data = esc%>% filter(dpi == 14))
summary(lm4)
simulateResiduals(lm4, plot=T)

lm4.5 <- glm(total_eye_score ~ temp, data=ti.mg %>% filter(dpi == 14 & treatment=="Inoculated"))
summary(lm4.5)
simulateResiduals(lm4.5, plot=T)

mod <- lm4.5
dat.new=expand.grid(total_eye_score=unique(ti.mg$total_eye_score),
                    temp=unique(ti.mg$temp))

dat.new$yhat=predict(mod, type="response", newdata = dat.new, re.form = NA)
# #prediction intervals
preds = predict(mod, type = "link", newdata = dat.new, se.fit =T, re.form = NA)
# #bind se's and fitted points
dat.new = cbind(dat.new, preds)
# #inverse link function
ilink <- family(mod)$linkinv
# #back transform CIs
dat.new <- transform(dat.new,
                     Fitted = ilink(fit),
                     Upper = ilink(fit + (2*se.fit)),
                     Lower = ilink(fit - (2*se.fit)))

peak_tes_fig <- ggplot(ti.mg %>% filter(dpi == 14), aes(x=temp, y=total_eye_score, color=temp))+
  geom_jitter(width=0.1, alpha=0.5)+
  geom_point(data=dat.new, aes(x=temp, y=yhat, color=temp))+
  geom_errorbar(data=dat.new, aes(x=temp, y=yhat, ymin=Lower, ymax=Upper, color=temp), width=0.0, alpha=0.75)+
  scale_y_continuous(limits = c(0,6))+
  scale_color_manual(values = c("#669BBC","#C1121F"))+
  labs(x="Temperature", y="Peak Eye Score", color="Temperature")

#mean eye score at peak infection
peak_tes_fig<- ggplot(esc %>% filter(dpi == 14), aes(x=status, y=mean_eye_score, color=temp))+
  geom_point(alpha=0.5)+
  geom_errorbar(data=dat.new, aes(x=status, y=yhat, ymin=Lower, ymax=Upper, color=temp), width=0.05, alpha=0.75)+
  scale_y_continuous(limits = c(0,4))

#calculate max tes
ti.mg.m <- ti.mg %>%
  group_by(band_number) %>%
  reframe(
    max_tes = max(total_eye_score1),
    diseased = diseased,
    ever_diseased=ever_diseased,
    temp=temp,
    treatment=treatment,
    groups=groups,
    dpi=dpi,
    band_number = band_number)

#what about total eye score?
lm5 <- glmmTMB(max_tes ~ temp, data = ti.mg.m%>% filter(treatment == "Inoculated" & dpi == 3))
summary(lm5)
simulateResiduals(lm5, plot=T)



mod <- lm5
dat.new=expand.grid(max_tes=unique(ti.mg.m$max_tes),
                    temp=unique(ti.mg.m$temp))

dat.new$yhat=predict(mod, type="response", newdata = dat.new, re.form = NA)
# #prediction intervals
preds = predict(mod, type = "link", newdata = dat.new, se.fit =T, re.form = NA)
# #bind se's and fitted points
dat.new = cbind(dat.new, preds)
# #inverse link function
ilink <- family(mod)$linkinv
# #back transform CIs
dat.new <- transform(dat.new,
                     Fitted = ilink(fit),
                     Upper = ilink(fit + (2*se.fit)),
                     Lower = ilink(fit - (2*se.fit)))

ggplot(ti.mg.m %>% filter(treatment=="Infected" & dpi ==3), aes(x=temp, y=max_tes, color=temp))+
  geom_jitter(height=0, width=0.2, alpha=0.5)+
  scale_color_manual(values = c("#C1121F","#669BBC"))+
  geom_point(data=dat.new, aes(x=temp, y=yhat, color=temp), size=2)+
  geom_errorbar(data=dat.new, aes(x=temp, ymin=Lower, ymax=Upper, y=yhat), width=0, alpha=0.75)

dat.new_summary <- dat.new %>%
  group_by(temp) %>%
  summarize(
    yhat = mean(yhat, na.rm = TRUE),  # Mean of predicted values
    Upper = mean(Upper, na.rm = TRUE),  # Mean of Upper CI
    Lower = mean(Lower, na.rm = TRUE))%>%  # Mean of Lower CI
  ungroup()

#write.csv(dat.new_summary, "/Users/jesse/Documents/GitHub/Temperature-and-Immunity-2022/Variability/max_tes_model_preds.csv", row.names = FALSE)


glmmTMB(tes ~ temp*dpi + (1|band_number),
        data=ti.mg, 
        ziformula = ~ temp,
        family = ziGamma(link = "log"))

####Pathogen Load####
ti.q <- ti %>%
  filter(!is.na(quantity))
unique(ti.q$dpi)

ti.q %>%
  dplyr::select(dpi, treatment, temp, groups)%>%
  tbl_summary(
    by="dpi"
  )%>%
  modify_header(
    label ~ "**Quantity**"
  )

quant_counts <- ti.q %>%
  group_by(dpi, groups) %>%
  summarise(counts = sum(!is.na(quantity)), .groups="drop")

ti.q$quantity1 <- ti.q$quantity + 1
hist(ti.q$quantity1)

#dpi 9 max path load
ggplot(ti.q, aes(x=dpi, y=quantity1, color=temp))+
  geom_point()+
  stat_summary(aes(x=dpi, y=quantity1, group = treatment), geom="line", fun="mean")

#only dpi 9 for peak path load for all inoculated birds
glm1<- lm(quantity1 ~ temp, data=ti.q %>% filter(dpi == 9 & treatment == "Inoculated"))
summary(glm1)

mod <- glm1
dat.new=expand.grid(quantity1=unique(ti.q$quantity1),
                    temp=unique(ti.q$temp))

dat.new$yhat=predict(mod, type="response", newdata = dat.new, re.form = NA)
# #prediction intervals
preds = predict(mod, type = "response", newdata = dat.new, se.fit =T, re.form = NA)
# #bind se's and fitted points
dat.new = cbind(dat.new, preds)
# #inverse link function
ilink <- family(mod)$linkinv
# #back transform CIs
dat.new <- transform(dat.new,
                     Fitted = ilink(fit),
                     Upper = ilink(fit + (2*se.fit)),
                     Lower = ilink(fit - (2*se.fit)))

peak_path_fig <- ggplot(ti.q %>% filter(treatment=="Inoculated" & dpi ==9), aes(x=temp, y=quantity1, color=temp))+
  geom_jitter(height=0, width=0.1, alpha=0.5)+
  scale_color_manual(values = c("#669BBC","#C1121F"))+
  geom_point(data=dat.new, aes(x=temp, y=yhat, color=temp), size=2)+
  geom_errorbar(data=dat.new, aes(x=temp, ymin=Lower, ymax=Upper, y=yhat), width=0, alpha=0.75)+
  labs(x="Temperature", y="Peak Pathogen Load", color="Temperature")


mean_figure <- wrap_plots(
  peak_tes_fig +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.title = element_text(face = "bold"),
      axis.text = element_text(face = "bold")
    ),
  peak_path_fig +
    theme_bw() +
    theme(
      axis.title = element_text(face = "bold"),
      axis.text = element_text(face = "bold")
    ),
  nrow = 1
)
print(mean_figure)
#model - needs work
glm2 <- glmmTMB(quantity1 ~ temp + treatment + dpi + (1|band_number), data=ti.q, family=Gamma(link="log"))
summary(glm2)
simulateResiduals(glm2, plot=T)

ti.q.m <- ti.q %>%
  group_by(band_number) %>%
  reframe(
    max_quant = max(quantity1),
    diseased = diseased,
    ever_diseased=ever_diseased,
    temp=temp,
    treatment=treatment,
    groups=groups,
    dpi=dpi,
    band_number = band_number)

#what about max path load?
lm6 <- glm(log10(max_quant) ~ temp, data = ti.q.m%>% filter(treatment == "Inoculated" & dpi == 3))
summary(lm6)
simulateResiduals(lm6, plot=T)

lm7<- glm(max_quant ~ temp, data = ti.q.m%>% filter(treatment == "Inoculated" & dpi == 3))
summary(lm7)
simulateResiduals(lm7, plot=T)

mod <- lm7
dat.new=expand.grid(max_quant=unique(ti.q.m$max_quant),
                    temp=unique(ti.q.m$temp))

dat.new$yhat=predict(mod, type="response", newdata = dat.new, re.form = NA)
# #prediction intervals
preds = predict(mod, type = "link", newdata = dat.new, se.fit =T, re.form = NA)
# #bind se's and fitted points
dat.new = cbind(dat.new, preds)
# #inverse link function
ilink <- family(mod)$linkinv
# #back transform CIs
dat.new <- transform(dat.new,
                     Fitted = ilink(fit),
                     Upper = ilink(fit + (2*se.fit)),
                     Lower = ilink(fit - (2*se.fit)))

ggplot(ti.q.m %>% filter(treatment=="Inoculated" & dpi ==3), aes(x=temp, y=max_quant, color=temp))+
  geom_jitter(height=0, width=0.2, alpha=0.5)+
  scale_color_manual(values = c("#C1121F","#669BBC"))+
  geom_point(data=dat.new, aes(x=temp, y=yhat, color=temp), size=2)+
  geom_errorbar(data=dat.new, aes(x=temp, ymin=Lower, ymax=Upper, y=yhat), width=0, alpha=0.75)

dat.new_summary <- dat.new %>%
  group_by(temp) %>%
  summarize(
    yhat = mean(yhat, na.rm = TRUE),  # Mean of predicted values
    Upper = mean(Upper, na.rm = TRUE),  # Mean of Upper CI
    Lower = mean(Lower, na.rm = TRUE))%>%  # Mean of Lower CI
  ungroup()

#write.csv(dat.new_summary, "/Users/jesse/Documents/GitHub/Temperature-and-Immunity-2022/Variability/max_quant_model_preds.csv", row.names = FALSE)

#variability
# Variability analysis for quantity1
ti.q.clean <- ti.q %>%
  filter(!is.na(quantity1) & treatment == "Inoculated")
ti.q.clean$log10_quantity <- log10(ti.q.clean$quantity1)+0.01

quantity_var <- ti.q.clean %>%
  group_by(dpi, temp) %>%
  dplyr::reframe(
    temp = temp,
    treatment = treatment,
    band_number = band_number,
    quantity1 = quantity1,
    log10_quantity=log10_quantity,
    
    # Calculations for quantity1
    mean_quantity = mean(log10_quantity, na.rm = TRUE),
    median_quantity = median(log10_quantity, na.rm = TRUE),
    bird_cv_quantity = calculate_cv(log10_quantity),
    bird_pv_quantity = calculate_pv(log10_quantity),
    bird_v2_quantity = calculate_v2(log10_quantity),
    
    # Bootstrapping for quantity1
    cv_bootstrap_quantity = list(replicate(n_boot, calculate_cv(sample(log10_quantity, replace = TRUE)))),
    pv_bootstrap_quantity = list(replicate(n_boot, calculate_pv(sample(log10_quantity, replace = TRUE)))),
    v2_bootstrap_quantity = list(replicate(n_boot, calculate_v2(sample(log10_quantity, replace = TRUE))))
  ) %>%
  mutate(
    # Confidence intervals for quantity1
    cv_lower_ci_quantity = map_dbl(cv_bootstrap_quantity, ~ quantile(.x, 0.025, na.rm = TRUE)),
    cv_upper_ci_quantity = map_dbl(cv_bootstrap_quantity, ~ quantile(.x, 0.975, na.rm = TRUE)),
    pv_lower_ci_quantity = map_dbl(pv_bootstrap_quantity, ~ quantile(.x, 0.025, na.rm = TRUE)),
    pv_upper_ci_quantity = map_dbl(pv_bootstrap_quantity, ~ quantile(.x, 0.975, na.rm = TRUE)),
    v2_lower_ci_quantity = map_dbl(v2_bootstrap_quantity, ~ quantile(.x, 0.025, na.rm = TRUE)),
    v2_upper_ci_quantity = map_dbl(v2_bootstrap_quantity, ~ quantile(.x, 0.975, na.rm = TRUE))
  ) %>%
  dplyr::select(-cv_bootstrap_quantity, -pv_bootstrap_quantity, -v2_bootstrap_quantity) %>%
  ungroup()

# Summary statistics for quantity_var
summary_quantity_tibble <- quantity_var %>%
  group_by(dpi, temp) %>%
  dplyr::reframe(
    mean_quant = mean(mean_quantity, na.rm = TRUE),
    median_eye_score = mean(median_quantity, na.rm = TRUE),
    mean_bird_cv_quantity = mean(bird_cv_quantity, na.rm = TRUE),
    mean_bird_pv_quantity = mean(bird_pv_quantity, na.rm = TRUE),
    mean_bird_v2_quantity = mean(bird_v2_quantity, na.rm = TRUE),
    mean_cv_lower_ci_quantity = mean(cv_lower_ci_quantity, na.rm = TRUE),
    mean_cv_upper_ci_quantity = mean(cv_upper_ci_quantity, na.rm = TRUE),
    mean_pv_lower_ci_quantity = mean(pv_lower_ci_quantity, na.rm = TRUE),
    mean_pv_upper_ci_quantity = mean(pv_upper_ci_quantity, na.rm = TRUE),
    mean_v2_lower_ci_quantity = mean(v2_lower_ci_quantity, na.rm = TRUE),
    mean_v2_upper_ci_quantity = mean(v2_upper_ci_quantity, na.rm = TRUE)
  )

ggplot(summary_quantity_tibble, aes(x = interaction(temp, dpi), color = temp)) +
  geom_point(aes(y = mean_bird_pv_quantity, shape = "PV"), size = 2) +
  geom_errorbar(aes(y = mean_bird_pv_quantity, ymin = mean_pv_lower_ci_quantity, ymax = mean_pv_upper_ci_quantity), width = 0) +
  geom_point(aes(y = mean_bird_cv_quantity, shape = "CV"), size = 2) +
  geom_point(aes(y = mean_bird_v2_quantity, shape = "V2"), size = 2) +
  
  labs(
    y = "Variability",
    x = "Temperature",
    shape = "Metric Type",
    color = "Temperature",
    title = "Variability in Log10(Quantity) by Temperature"
  ) +
  scale_color_manual(values = c("#669BBC", "#C1121F")) +
  scale_shape_manual(name = "Variability Metric",
                     values = c("PV" = 16, "CV" = 2, "V2" = 0)) +
 # coord_flip() +
  theme_minimal()+
  theme(
    axis.text.x = element_text(angle = 45, hjust=1)
  )

# Write the summarized tibble to a CSV file
#write.csv(summary_quantity_tibble, "/Users/jesse/Documents/GitHub/Temperature-and-Immunity-2022/Variability/Quantity_Variability.csv")


#Max Quantity
# Compute max_quantity for each band_number
ti.q.max <- ti.q %>%
  group_by(band_number) %>%
  reframe(max_quantity = max(quantity1, na.rm = TRUE),
          max_log10_quantity = max(log10_quantity))%>%
  ungroup()

# Merge max_quantity back into the dataset
ti.q <- merge(ti.q, ti.q.max, by = "band_number")

# Variability analysis for max_quantity
max_quantity_var <- ti.q %>%
  filter(treatment == "Inoculated")%>%
  group_by(temp) %>%
  dplyr::reframe(
    temp = temp,
    treatment = treatment,
    band_number = band_number,
    max_quantity = max_quantity,
    max_log10_quantity = max_log10_quantity,
    
    # Calculations for max_quantity
    mean_max_quantity = mean(max_log10_quantity, na.rm = TRUE),
    median_max_quantity = median(max_log10_quantity, na.rm = TRUE),
    bird_cv_max_quantity = calculate_cv(max_log10_quantity),
    bird_pv_max_quantity = calculate_pv(max_log10_quantity),
    bird_v2_max_quantity = calculate_v2(max_log10_quantity),
    
    # Bootstrapping for max_quantity
    cv_bootstrap_max_quantity = list(replicate(n_boot, calculate_cv(sample(max_log10_quantity, replace = TRUE)))),
    pv_bootstrap_max_quantity = list(replicate(n_boot, calculate_pv(sample(max_log10_quantity, replace = TRUE)))),
    v2_bootstrap_max_quantity = list(replicate(n_boot, calculate_v2(sample(max_log10_quantity, replace = TRUE))))
  ) %>%
  mutate(
    # Confidence intervals for max_quantity
    cv_lower_ci_max_quantity = map_dbl(cv_bootstrap_max_quantity, ~ quantile(.x, 0.025, na.rm = TRUE)),
    cv_upper_ci_max_quantity = map_dbl(cv_bootstrap_max_quantity, ~ quantile(.x, 0.975, na.rm = TRUE)),
    pv_lower_ci_max_quantity = map_dbl(pv_bootstrap_max_quantity, ~ quantile(.x, 0.025, na.rm = TRUE)),
    pv_upper_ci_max_quantity = map_dbl(pv_bootstrap_max_quantity, ~ quantile(.x, 0.975, na.rm = TRUE)),
    v2_lower_ci_max_quantity = map_dbl(v2_bootstrap_max_quantity, ~ quantile(.x, 0.025, na.rm = TRUE)),
    v2_upper_ci_max_quantity = map_dbl(v2_bootstrap_max_quantity, ~ quantile(.x, 0.975, na.rm = TRUE))
  ) %>%
  dplyr::select(-cv_bootstrap_max_quantity, -pv_bootstrap_max_quantity, -v2_bootstrap_max_quantity) %>%
  ungroup()

# Summary statistics for max_quantity
summary_max_quantity_tibble <- max_quantity_var %>%
  group_by(temp) %>%
  dplyr::reframe(
    mean_max_quantity = mean(mean_max_quantity, na.rm = TRUE),
    median_max_quantity = mean(median_max_quantity, na.rm = TRUE),
    mean_bird_cv_max_quantity = mean(bird_cv_max_quantity, na.rm = TRUE),
    mean_bird_pv_max_quantity = mean(bird_pv_max_quantity, na.rm = TRUE),
    mean_bird_v2_max_quantity = mean(bird_v2_max_quantity, na.rm = TRUE),
    mean_cv_lower_ci_max_quantity = mean(cv_lower_ci_max_quantity, na.rm = TRUE),
    mean_cv_upper_ci_max_quantity = mean(cv_upper_ci_max_quantity, na.rm = TRUE),
    mean_pv_lower_ci_max_quantity = mean(pv_lower_ci_max_quantity, na.rm = TRUE),
    mean_pv_upper_ci_max_quantity = mean(pv_upper_ci_max_quantity, na.rm = TRUE),
    mean_v2_lower_ci_max_quantity = mean(v2_lower_ci_max_quantity, na.rm = TRUE),
    mean_v2_upper_ci_max_quantity = mean(v2_upper_ci_max_quantity, na.rm = TRUE)
  )

# Plot variability for max_quantity
ggplot(summary_max_quantity_tibble, aes(x = temp, color = temp)) +
  geom_point(aes(y = mean_bird_pv_max_quantity, shape = "PV"), size = 2) +
  geom_errorbar(aes(y = mean_bird_pv_max_quantity, ymin = mean_pv_lower_ci_max_quantity, ymax = mean_pv_upper_ci_max_quantity), width = 0) +
  geom_point(aes(y = mean_bird_cv_max_quantity, shape = "CV"), size = 2) +
  geom_point(aes(y = mean_bird_v2_max_quantity, shape = "V2"), size = 2) +
  
  labs(
    y = "Variability",
    x = "Temperature",
    shape = "Metric Type",
    color = "Temperature",
    title = "Variability in Max Quantity by Temperature"
  ) +
  scale_color_manual(values = c("#669BBC", "#C1121F")) +
  scale_shape_manual(name = "Variability Metric",
                     values = c("PV" = 16, "CV" = 2, "V2" = 0)) +
  #coord_flip() +
  theme_minimal()

# Write the summarized tibble to a CSV file
#write.csv(summary_max_quantity_tibble, "/Users/jesse/Documents/GitHub/Temperature-and-Immunity-2022/Variability/Max_log_Quantity_Variability.csv")

####Synthesis####
setwd("/Users/jesse/Documents/GitHub/Temperature-and-Immunity-2022/Variability/")
var <- read.csv("Variability_Synthesis.csv", row.names = NULL)
colnames(var)

unique(var$test)

var_fig <- var %>% filter(test %in% c("eye_score_dpi14", "quantity1_dpi9", "eye_score_dpi8", "eye_score_dpi9"))

peak_path_var <- ggplot(data=var_fig %>% filter(treatment=="Inoculated" & test=="quantity1_dpi9"), aes(x=temp, y=quantity1, color=temp))+
  scale_color_manual(values = c("#669BBC","#C1121F"))+
  # geom_point(aes(x=temp, y=cv, color=temp), size=2, shape = 2, position=position_dodge(width=1))+
  # geom_errorbar(aes(x=temp, ymin=cv_lower_ci, ymax=cv_upper_ci, y=cv), width=0.1, alpha=0.75, position=position_dodge(width=1), lty="dashed")+
  geom_point(aes(x=temp, y=pv, color=temp), size=2)+
  geom_errorbar(aes(x=temp, ymin=pv_lower_ci, ymax=pv_upper_ci, y=pv), width=0, alpha=0.75)+
  #scale_y_continuous(limits=c(0,3))+
  scale_y_continuous(limits=c(0,1))+
   labs(x="Temperature", y="Peak Pathogen Load PV", color="Temperature")

peak_tes_var <- ggplot(var_fig %>% filter(treatment=="Inoculated" & test == "eye_score_dpi14"), aes(x=temp, y=pv, color=temp))+
  scale_color_manual(values = c("#669BBC","#C1121F"))+
  # geom_point(aes(x=temp, y=cv, color=temp), size=2, shape = 2, position=position_dodge(width=1))+
  # geom_errorbar(aes(x=temp, ymin=cv_lower_ci, ymax=cv_upper_ci, y=cv), width=0.1, alpha=0.75, position=position_dodge(width=1), lty="dashed")+
  geom_point(aes(x=temp, y=pv, color=temp), size=2)+
  geom_errorbar(aes(x=temp, ymin=pv_lower_ci, ymax=pv_upper_ci, y=pv), width=0, alpha=0.75)+
  #scale_y_continuous(limits=c(0,1.25))+
  scale_y_continuous(limits=c(0,1))+
  labs(x="Temperature", y="Peak Eye Score PV", color="Temperature")
peak_tes_var


var_figure <- wrap_plots(
  peak_tes_var +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.title = element_text(face = "bold"),
      axis.text = element_text(face = "bold")
    ),
  peak_path_var +
    theme_bw() +
    theme(
      axis.title = element_text(face = "bold"),
      axis.text = element_text(face = "bold")
    ),
  nrow = 1
)

print(var_figure)
peak_tes_var_tg <- ggplot(var_fig %>% filter(treatment=="Inoculated" & test == "eye_score_dpi8"), aes(x=temp, y=pv, color=temp))+
  scale_color_manual(values = c("#669BBC","#C1121F"))+
  # geom_point(aes(x=temp, y=cv, color=temp), size=2, shape = 2, position=position_dodge(width=1))+
  # geom_errorbar(aes(x=temp, ymin=cv_lower_ci, ymax=cv_upper_ci, y=cv), width=0.1, alpha=0.75, position=position_dodge(width=1), lty="dashed")+
  geom_point(aes(x=temp, y=pv, color=temp), size=2)+
  geom_errorbar(aes(x=temp, ymin=pv_lower_ci, ymax=pv_upper_ci, y=pv), width=0.1, alpha=0.75)+
  #scale_y_continuous(limits=c(0,1.25))+
  scale_y_continuous(limits=c(0,1))+
  labs(x="Temperature", y="Peak Eye Score PV", color="Temperature", title = "TG23 D8")
peak_tes_var_tg

peak_tes_var9 <- ggplot(var_fig %>% filter(treatment=="Inoculated" & test == "eye_score_dpi9"), aes(x=temp, y=pv, color=temp))+
  scale_color_manual(values = c("#669BBC","#C1121F"))+
  # geom_point(aes(x=temp, y=cv, color=temp), size=2, shape = 2, position=position_dodge(width=1))+
  # geom_errorbar(aes(x=temp, ymin=cv_lower_ci, ymax=cv_upper_ci, y=cv), width=0.1, alpha=0.75, position=position_dodge(width=1), lty="dashed")+
  geom_point(aes(x=temp, y=pv, color=temp), size=2)+
  geom_errorbar(aes(x=temp, ymin=pv_lower_ci, ymax=pv_upper_ci, y=pv), width=0.1, alpha=0.75)+
  #scale_y_continuous(limits=c(0,1.25))+
  scale_y_continuous(limits=c(0,1))+
  labs(x="Temperature", y="Peak Eye Score PV", color="Temperature", title = "TI22 D9")

var_figure_tg <- wrap_plots(
                            #peak_tes_var +  theme(legend.position = "none"),
                            peak_tes_var9 +  theme(legend.position = "none"),
                            peak_tes_var_tg + theme(legend.position = "none"),
                         peak_path_var,
                         nrow=1)

var_bio <- var %>%
  filter(test %in% c("elisa_dpi28", "phagocytosis", "eye_score_dpi14", "quantity1_dpi9", "peak_fever_score"))

ggplot(var_bio %>% filter(treatment == "Inoculated"), aes(x=test, color=temp))+
  geom_point(aes(y=pv), size=1, position=position_dodge(width=0.5))+
  geom_errorbar(aes(y=pv, ymin=pv_lower_ci, ymax=pv_upper_ci), width=0, position=position_dodge(width=0.5))+
  scale_color_manual(values=c(temp_cols))+
  scale_y_continuous(limits = c(0,1))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face="bold"))+
  labs(x="Assay", y="PV", color="Temperature")



mean_bio <- var %>%
  filter(test %in% c("elisa_dpi28", "phagocytosis", "max_eye_score", "max_quantity1", "max_log10_quantity1", "max_fever_score"))
var_bio_mean <- mean_bio %>%
  filter(treatment == "Inoculated")%>%
  group_by(temp) %>%
  summarize(mean_pv = mean(pv),
            mean_cv = mean(cv),
            mean_v2 = mean(v2))

ggplot(mean_bio %>% filter(treatment == "Inoculated"), aes(x=test, color=temp))+
  geom_point(aes(y=pv, shape = "PV"), size=2, position=position_dodge(width=0.6))+
  #geom_point(aes(y=v2, shape="V2"), size=2, position=position_dodge(width=1))+
  geom_point(aes(y=cv, shape="CV"), size =2, position=position_dodge(width=0.2))+
  #geom_errorbar(aes(y=v2, ymin=v2_lower_ci, ymax=v2_upper_ci), width=0, position=position_dodge(width=1), alpha=0.5, lty="dashed")+
  geom_errorbar(aes(y=cv, ymin=cv_lower_ci, ymax=cv_upper_ci), width=0, position=position_dodge(width=0.2), alpha=0.5, lty="dashed")+
  geom_errorbar(aes(y=pv, ymin=pv_lower_ci, ymax=pv_upper_ci), width=0, position=position_dodge(width=0.6))+
  #scale_y_continuous(limits = c(0,1))+
  scale_color_manual(values=c(temp_cols))+
  # geom_hline(data=var_bio_mean %>% filter(temp=="Cold"), aes(yintercept=mean_pv), color="#669BBC", lty="dashed", alpha=0.5)+
  # geom_hline(data=var_bio_mean %>% filter(temp=="Warm"), aes(yintercept=mean_pv), color="#C1121F", lty="dashed", alpha=0.5)+
  coord_flip()+
  scale_shape_manual(name="Variability Metric",
                     values=c("PV" = 16, "V2" = 1, "CV" = 0))+
  labs(x="Test", y="Variability", color="Temperature")


ab<-ggplot(var_bio %>% filter(test == "elisa_dpi28"), aes(x="Antibody Levels DPI 28", y=yhat, color= temp, shape=treatment))+
  geom_point(size=2.5, position=position_dodge(2))+
  geom_errorbar(aes(x="Antibody Levels DPI 28", y=yhat, ymin=Lower, ymax=Upper), position=position_dodge(2), width=0)+
  scale_color_manual(values=c(temp_cols))+
  geom_point(data=ti.pi %>% filter(dpi ==28), aes(x="Antibody Levels DPI 28", y=elisa_od, color=temp), position=position_dodge(2), alpha=0.5, size=1)+
  scale_y_continuous(limits=c(0.04, 0.08))+
  #labs(x="", y="Antibody Levels")+
  labs(x="", y="")+
  coord_flip()
ab

tes<-ggplot(var_bio %>% filter(test == "max_eye_score"), aes(x="Max Eye Score", y=yhat, color= temp, shape=treatment))+
  geom_point(size=2.5, position=position_dodge(2))+
  geom_errorbar(aes(x="Max Eye Score", y=yhat, ymin=Lower, ymax=Upper), position=position_dodge(2), width=0)+
  geom_point(data=ti.mg %>% filter(dpi ==28), aes(x="Max Eye Score", y=max_eye_score, color=temp), position=position_dodge(2), alpha=0.5, size=1)+
  #scale_y_continuous(limits=c(0, 2))+
  scale_color_manual(values=c(temp_cols))+
  #labs(x="", y="Max Eye Score")+
  labs(x="", y="")#+
  coord_flip()
tes
quant<-ggplot(var_bio %>% filter(test == "max_quantity1" & treatment == "Infected"), aes(x="Max Pathogen Load", y=yhat, color= temp, shape=treatment))+
 geom_point(size=2.5, position=position_dodge(2))+
  geom_errorbar(aes(x="Max Pathogen Load", y=yhat, ymin=Lower, ymax=Upper), position=position_dodge(2), width=0)+
  geom_point(data=ti.q %>% filter(treatment == "Infected" & dpi == 28), aes(x="Max Pathogen Load", y=max_quantity, color=temp), position=position_dodge(2), alpha=0.5, size=1)+
  #scale_y_continuous(limits=c(0, 120000))+
  scale_color_manual(values=c(temp_cols))+
  #labs(x="", y="Max Pathogen Load")+
  labs(x="", y="")#+
  coord_flip()
quant
quant10<-ggplot(var_bio %>% filter(test == "max_log10_quantity1"), aes(x="Max Log10(Pathogen Load)", y=yhat, color= temp, shape=treatment))+
  geom_point(size=2.5, position=position_dodge(2))+
  geom_errorbar(aes(x="Max Log10(Pathogen Load)", y=yhat, ymin=Lower, ymax=Upper), position=position_dodge(2), width=0)+
  geom_point(data=ti.q %>% filter(treatment == "Infected" & dpi == 28), aes(x="Max Log10(Pathogen Load)", y=max_log10_quantity, color=temp), position=position_dodge(2), alpha=0.5, size=1)+
  #scale_y_continuous(limits=c(0, 5))+
  scale_color_manual(values=c(temp_cols))+
  labs(x="", y="",color = "Temperature", shape="Treatment")#+
  coord_flip()
quant10
phago<-ggplot(var_bio %>% filter(test == "phagocytosis"), aes(x="Phagocytosis", y=yhat, color= temp, shape=treatment))+
  geom_point(size=2.5, position=position_dodge(2))+
  geom_errorbar(aes(x="Phagocytosis", y=yhat, ymin=Lower, ymax=Upper), position=position_dodge(2), width=0)+
  geom_point(data=ti.p, aes(x="Phagocytosis", y=phago_score, color=temp), position=position_dodge(2), alpha=0.5, size=1)+
  #scale_y_continuous(limits=c(0, 0.2))+
  scale_color_manual(values=c(temp_cols))+
  labs(x="", y="", color="Temperature", shape="Treatment")#+
  coord_flip()
phago
fever<-ggplot(var_bio %>% filter(test == "max_fever_score"), aes(x="Max Fever", y=yhat, color= temp, shape=treatment))+
  geom_point(size=2.5, position=position_dodge(2))+
  geom_errorbar(aes(x="Max Fever", y=yhat, ymin=Lower, ymax=Upper), position=position_dodge(2), width=0)+
  geom_point(data=ti.fm %>% filter(dpi==28), aes(x="Max Fever", y=max_fever_score, color=temp), position=position_dodge(2), alpha=0.5, size=1)+
  #scale_y_continuous(limits=c(30, 40))+
  scale_color_manual(values=c(temp_cols))+
  labs(x="", y="", shape = "Treatment", color = "Temperature")+
  theme_minimal()#+
  coord_flip()
fever

# Combine plots theme minimal and no legend
combined_plot_vert <- wrap_plots(
  tes + coord_flip() + theme(legend.position = "none"),
  quant + coord_flip() +  theme(legend.position = "none"),
  quant10 + coord_flip() +  theme(legend.position = "none"),
  ab + coord_flip() +  theme(legend.position = "none"),
  phago + coord_flip() +  theme(legend.position = "none"),
  fever + coord_flip() + theme(legend.position = "none"), # Keep the legend for `fever`
  #nrow = 2
  ncol=1
) 

# Display the combined plot
print(combined_plot_vert)

# Combine plots theme minimal and horizontal no legend
combined_plot_ho <- wrap_plots(
  tes + theme_bw() + theme(legend.position = "none"),
  quant +  theme_bw() +  theme(legend.position = "none"),
  quant10 + theme_bw() + theme(legend.position = "none"),
  ab + theme_bw() + theme(legend.position = "none"),
  phago + theme_bw() + theme(legend.position = "none"),
  fever + theme_bw() + theme(legend.position = "none"), # Keep the legend for `fever`
  nrow = 3
  #ncol=1
) 

# Display the combined plot
print(combined_plot_ho)

ab<-ggplot(var_bio %>% filter(test == "elisa_dpi28"), aes(x=interaction(temp,treatment), y=yhat, color= temp, shape=treatment))+
  geom_point(size=2.5, position=position_dodge(2))+
  geom_errorbar(aes(x=interaction(temp,treatment), y=yhat, ymin=Lower, ymax=Upper), position=position_dodge(2), width=0)+
  scale_color_manual(values=c(temp_cols))+
  geom_point(data=ti.pi %>% filter(dpi ==28), aes(x=interaction(temp,treatment), y=elisa_od, color=temp), position=position_dodge(2), alpha=0.5, size=1)+
  scale_y_continuous(limits=c(0.04, 0.08))+
  #labs(x="", y="Antibody Levels")+
  labs(x="", y="", title = "Antibody Levels")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))#+
  coord_flip()
ab

tes<-ggplot(var_bio %>% filter(test == "max_eye_score"), aes(x=temp, y=yhat, color= temp, shape=treatment))+
  geom_point(size=2.5, position=position_dodge(2))+
  geom_errorbar(aes(x=temp, y=yhat, ymin=Lower, ymax=Upper), position=position_dodge(2), width=0)+
  geom_point(data=ti.mg %>% filter(dpi ==28), aes(x=temp, y=max_eye_score, color=temp), position=position_dodge(2), alpha=0.5, size=1)+
  #scale_y_continuous(limits=c(0, 2))+
  scale_color_manual(values=c(temp_cols))+
  #labs(x="", y="Max Eye Score")+
  labs(x="", y="", title="Max Eye Score")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))#+
coord_flip()
tes
quant<-ggplot(var_bio %>% filter(test == "max_quantity1" & treatment == "Infected"), aes(x=temp, y=yhat, color= temp, shape=treatment))+
  geom_point(size=2.5, position=position_dodge(2))+
  geom_errorbar(aes(x=temp, y=yhat, ymin=Lower, ymax=Upper), position=position_dodge(2), width=0)+
  geom_point(data=ti.q %>% filter(treatment == "Infected" & dpi == 28), aes(x=temp, y=max_quantity, color=temp), position=position_dodge(2), alpha=0.5, size=1)+
  #scale_y_continuous(limits=c(0, 120000))+
  scale_color_manual(values=c(temp_cols))+
  labs(x="", y="", title="Max Pathogen Load")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))#+
coord_flip()
quant

quant10<-ggplot(var_bio %>% filter(test == "max_log10_quantity1"), aes(x=temp, y=yhat, color= temp, shape=treatment))+
  geom_point(size=2.5, position=position_dodge(2))+
  geom_errorbar(aes(x=temp, y=yhat, ymin=Lower, ymax=Upper), position=position_dodge(2), width=0)+
  geom_point(data=ti.q %>% filter(treatment == "Infected" & dpi == 28), aes(x=temp, y=max_log10_quantity, color=temp), position=position_dodge(2), alpha=0.5, size=1)+
  #scale_y_continuous(limits=c(0, 5))+
  scale_color_manual(values=c(temp_cols))+
  labs(x="", y="",color = "Temperature", shape="Treatment", title = "Max Log10(Pathogen Load)")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))#+
coord_flip()
quant10
phago<-ggplot(var_bio %>% filter(test == "phagocytosis"), aes(x=interaction(temp,treatment), y=yhat, color= temp, shape=treatment))+
  geom_point(size=2.5, position=position_dodge(2))+
  geom_errorbar(aes(x=interaction(temp,treatment), y=yhat, ymin=Lower, ymax=Upper), position=position_dodge(2), width=0)+
  geom_point(data=ti.p, aes(x=interaction(temp,treatment), y=phago_score, color=temp), position=position_dodge(2), alpha=0.5, size=1)+
  #scale_y_continuous(limits=c(0, 0.2))+
  scale_color_manual(values=c(temp_cols))+
  labs(x="", y="", color="Temperature", shape="Treatment", title = "Phagocytosis")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))#+
coord_flip()
phago
fever<-ggplot(var_bio %>% filter(test == "max_fever_score"), aes(x=interaction(temp,treatment), y=yhat, color= temp, shape=treatment))+
  geom_point(size=2.5, position=position_dodge(2))+
  geom_errorbar(aes(x=interaction(temp,treatment), y=yhat, ymin=Lower, ymax=Upper), position=position_dodge(2), width=0)+
  geom_point(data=ti.fm %>% filter(dpi==28), aes(x=interaction(temp,treatment), y=max_fever_score, color=temp), position=position_dodge(2), alpha=0.5, size=1)+
  #scale_y_continuous(limits=c(30, 40))+
  scale_color_manual(values=c(temp_cols))+
  labs(x="", y="", shape = "Treatment", color = "Temperature", title= "Max Fever")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))#+
coord_flip()
fever

# Combine plots theme minimal and horizontal no legend
combined_plot_ho_lab <- wrap_plots(
  tes + theme_bw() + theme(legend.position = "none"),
  ab + theme_bw() + theme(legend.position = "none"),
  quant +  theme_bw() +  theme(legend.position = "none"),
  quant10 + theme_bw() + theme(legend.position = "none"),
  phago + theme_bw() + theme(legend.position = "none"),
  fever + theme_bw() + theme(legend.position = "none"), # Keep the legend for `fever`
  nrow = 3
  #ncol=1
) 

# Display the combined plot
print(combined_plot_ho_lab)

combined_plot_small <- wrap_plots(
  tes + theme_bw() + theme(legend.position = "none"),
  quant +  theme_bw() +  theme(legend.position = "none"),
  #quant10 + theme_bw() + theme(legend.position = "none"),
  nrow = 1
  #ncol=1
) 
print(combined_plot_small)

####Can I use a LMM to test whether temperature affects variability in each bioassay?
mod <- lmer(cv ~ temp + treatment + (1|test), data=var_bio)
summary(mod)
simulateResiduals(mod, plot=T)
plot(allEffects(mod))

ggplot(var_bio, aes(x = interaction(temp, treatment), y = pv, color = temp, shape=test)) +
  geom_point(size = 2) +
  labs(title = "Variability (PV) by Temperature and Bioassay", 
       x = "Temperature", 
       y = "Proportional Variability") +
  theme_minimal()
