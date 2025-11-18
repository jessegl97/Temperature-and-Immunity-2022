#Temperature and Immunity 2022 Analysis
  #Updated Oct 2025
  #JGL
  #Mean Models

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
p_load(ggplot2, gridExtra, gtsummary)


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

####Antibody Analysis####
#Do differences in temperatures affect antibody responses to MG infection?

source("r_scripts/dataCleaning_antibody.R")

ti %>%
  dplyr::select(dpi, treatment, temp, groups)%>%
  tbl_summary(
    by="dpi"
  )%>%
  modify_header(
    label ~ "**All Birds**"
  )

unique(ti.ab$band_number)

#Antibody analysis sample sizes; see dataCleaning_antibody.R for removal breakdown
ti.ab %>%
  dplyr::select(dpi, treatment, temp, groups)%>%
  tbl_summary(
    by="dpi"
  )%>%
  modify_header(
    label ~ "**All Birds**"
  )

unique(ti.ab$dpi)
ti.ab$dpi <- as.factor(ti.ab$dpi)

#Set reference categories "Warm" and "Sham"
ti.ab$temp <- relevel(as.factor(ti.ab$temp), ref = "Warm")
ti.ab$treatment <- relevel(as.factor(ti.ab$treatment), ref = "Control")
ti.ab$dpi.f <- as.factor(ti.ab$dpi)

ti.ab.mod <- ti.ab %>%
  filter(dpi == 28)

hist(ti.ab.mod$elisa_od)

#three way interaction
a <- glm(elisa_od ~ treatment * temp,
             data=ti.ab.mod,
             family=Gamma(link="log"))

simulateResiduals(a, plot=T)

#treat x temp
b <- glm(elisa_od ~ treatment + temp,
             data=ti.ab.mod,
             family=Gamma(link="log"))

c <- glm(elisa_od ~ treatment * temp * sex,
         data=ti.ab.mod,
         family=Gamma(link="log"))

d <- glm(elisa_od ~ treatment * temp + sex,
         data=ti.ab.mod,
         family=Gamma(link="log"))

e <- glm(elisa_od ~ treatment,
         data=ti.ab.mod,
         family=Gamma(link="log"))

f <- glm(elisa_od ~ temp,
       data=ti.ab.mod,
       family=Gamma(link="log"))
#Null
g <- glm(elisa_od ~ 1,
             data=ti.ab,
             family=Gamma(link="log"))

aictab(cand.set=list(a, b, c,  d, e, f, g), 
       modnames=c("a",  "b", "c", "d",  "e", "f", "g"))



#Model Selection
#three way interaction
a <- glmmTMB(elisa_od ~ treatment * temp * dpi.f + (1|band_number),
             data=ti.ab,
             family=Gamma(link="log"))

#treat x temp
b <- glmmTMB(elisa_od ~ treatment * temp + dpi.f + (1|band_number),
             data=ti.ab,
             family=Gamma(link="log"))

#treat + temp
c <- glmmTMB(elisa_od ~ treatment + temp + dpi.f + (1|band_number),
             data=ti.ab,
             family=Gamma(link="log"))

#treat * dpi
d <- glmmTMB(elisa_od ~ treatment * dpi.f + temp + (1|band_number),
             data=ti.ab,
             family=Gamma(link="log"))

#temp * dpi
e <- glmmTMB(elisa_od ~ temp * dpi.f + treatment + (1|band_number),
             data=ti.ab,
             family=Gamma(link="log"))

e.5 <- glmmTMB(elisa_od ~ temp * dpi.f + (1|band_number),
               data=ti.ab,
               family=Gamma(link="log"))

f <- glmmTMB(elisa_od~treatment+temp+dpi.f + treatment:dpi.f + (1|band_number),
             data=ti.ab, 
             family=Gamma(link = "log"))

fs <- glmmTMB(elisa_od~treatment+temp+dpi.f + treatment:dpi.f + sex + (1|band_number),
              data=ti.ab, 
              family=Gamma(link = "log"))

#Null
g <- glmmTMB(elisa_od ~ 1 + (1|band_number),
             data=ti.ab,
             family=Gamma(link="log"))

aictab(cand.set=list(a, b, c,  d, e, e.5, f, fs, g), 
       modnames=c("a",  "b", "c", "d",  "e", "e.5", "f", "fs", "g"))

BIC(f, d)
summary(d)
summary(f)


#Final Model: Are antibodies predicted by treatment, temp, dpi, or the interaction between temp and dpi while controlling for individual band number
lm1<-glmmTMB(elisa_od~treatment+temp+dpi + treatment:dpi + (1|band_number),
                 data=ti.ab, 
                 family=Gamma(link = "log"))

simulateResiduals(lm1, plot=T)

lm1s<-glmmTMB(elisa_od~treatment+temp+dpi + treatment:dpi + sex + (1|band_number),
             data=ti.ab, 
             family=Gamma(link = "log"))

#sex not significant
anova(lm1, lm1s, type="Chisq")
summary(lm1s)


Anova(lm1, type = "II")

emm <- emmeans(lm1, ~ treatment | temp*dpi, type = "response")

pairs(emm, by=c("temp","dpi"), adjust = "tukey")

emm_df <- as.data.frame(emm)

dodge = position_dodge(0.75)
ggplot(emm_df, aes(x = groups, y = response, color = groups, shape=dpi, groups=dpi)) +
  geom_jitter(
    data = ti.ab,
    aes(y = elisa_od, color = groups, shape = dpi),
    alpha = 0.75, size = 3,
    position = position_jitterdodge(dodge.width = 0.75, jitter.width=0.5)
  ) +
  geom_point(size = 3, color="black", stroke = 1,
             position = dodge) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.0, color="black",
                position = dodge) +
  scale_color_manual(values = treat_colors, breaks = lvl, name = "Treatment") +
  scale_shape_manual(values = c(0, 1, 16))+
  #facet_wrap(~ dpi, nrow = 2) +
  labs(x = "Treatment Group", y = "ELISA OD", shape = "Days Post Inoculation") +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggplot(emm_df, aes(x = dpi, y = response, color = groups, groups=groups)) +
  geom_jitter(
    data = ti.ab,
    aes(y = elisa_od, color = groups, shape = dpi),
    alpha = 0.75, size = 3,
    position = position_jitterdodge(dodge.width = 0.75, jitter.width=0.2)
  ) +
  geom_point(size = 3, color="black",
             position = dodge) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.0, color="black",
                position = dodge) +
  scale_color_manual(values = treat_colors, breaks = lvl, name = "Treatment") +
  scale_shape_manual(values = c(0, 1, 16))+
  #facet_wrap(~ dpi, nrow = 2) +
  labs(x = "Treatment Group", y = "ELISA OD", shape = "Days Post Inoculation") +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggplot(ti.ab, aes(x=dpi, y=elisa_od, color=groups))+
  geom_point(aes(shape=as.factor(ever_diseased)))+
  geom_line(aes(group = as.factor(band_number)))+
  scale_color_manual(values=treat_colors)+
  facet_grid(~ever_infected~ever_diseased)

#set contrasts so that results do not depend on reference-level
#only for type II Anova: 
  #Intercept corresponds to overall mean across all factor levels
  #each coefficient represents how far a given level is from the overall mean
  #the model fit is the same to default setting, but effects are represented differently
    #Have to change for type II Anova so that there is no reference-category. 
    #Every level contributes equally
  #contr.sum = changes ordered factors into dummy variables; symmetric coding = no reference
  #contr.poly = encodes polynomial trends for ordered factors
#options(contrasts = c("contr.sum", "contr.poly"))

#Anova type II asks whether there is a main effect overall, averaged across the other variables
lm1_anova_II <- car::Anova(lm1, type = "II")

# Analysis of Deviance Table (Type II Wald chisquare tests)
# 
# Response: elisa_od
#                   Chisq Df Pr(>Chisq)    
#   treatment     12.1200  1  0.0004988 ***
#   temp           1.7131  1  0.1905797    
#   dpi            7.4981  1  0.0061766 ** 
#   treatment:dpi  4.2154  1  0.0400585 *   

#Type III asks whether there is a maen effect that is consistent across the other variables, even after accounting for the interaction
#change contrasts back to baseline - back to reference categories
  #contr.treatment = treatment contrasts where the first (reference) level for each factor is baseline
  #contr.poly stays the same


emm <- emmeans(lm1, ~ treatment | temp*dpi, type = "response")

pairs(emm, by=c("temp","dpi"), adjust = "tukey")

emm_df <- as.data.frame(emm)

lvl <- c("Warm Control", "Warm Inoculated", "Cold Control", "Cold Inoculated")

emm_df <- as.data.frame(emm) %>%
  mutate(
    groups = paste(temp, ifelse(treatment == "Inoculated", "Inoculated", "Control")),
    groups = factor(groups, levels = lvl)
  )

ti.ab <- ti.ab %>%
  mutate(
    groups = paste(temp, ifelse(treatment == "Inoculated", "Inoculated", "Control")),
    groups = factor(groups, levels = lvl),
    dpi    = factor(dpi)
  )

dodge = position_dodge(0.75)
ggplot(emm_df, aes(x = groups, y = response, color = groups, shape=dpi, groups=dpi)) +
  geom_jitter(
    data = ti.ab,
    aes(y = elisa_od, color = groups, shape = dpi),
    alpha = 0.75, size = 2,
    position = position_jitterdodge(dodge.width = 0.75, jitter.width=0.2)
  ) +
  geom_point(size = 2.5, color="black",
             position = dodge) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.0, color="black",
                position = dodge) +
  scale_color_manual(values = treat_colors, breaks = lvl, name = "Treatment") +
  scale_shape_manual(values = c(0, 1, 16))+
  #facet_wrap(~ dpi, nrow = 2) +
  labs(x = "Treatment Group", y = "ELISA OD", shape = "Days Post Inoculation") +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


#options(contrasts = c("contr.treatment", "contr.poly"))
lm1_anova_III <- car::Anova(lm1, type = "III")

summary(lm1)


plot(allEffects(lm1))

# Table S1
# write_xlsx(
#   list(
#     Fixed_Effects = fixed_effects,
#     Dispersion_Model = dispersion_effects,
#     Random_Effects = random_effects,
#     Pairwise_Comparisons = pairwise_df
#   ),
#   "/Users/jesse/Documents/GitHub/EEID_1A_Mechanistic_Link/EEID_1A_Antibody_Analysis_files/Results/Supplementary/TableS1.xlsx"
# )

 #Model Predictions
mod <- lm1
# dat.new=expand.grid(#band_number = unique(ti.ab$band_number),
#                     temp = unique(ti.ab$temp),
#                     treatment= unique(ti.ab$treatment),
#                     dpi = unique(ti.ab$dpi))
dat.new <- dplyr::distinct(ti.ab[, c("temp","treatment","dpi")])

dat.new$yhat=predict(mod, type="response", newdata = dat.new, re.form = NA)
#prediction intervals
preds = predict(mod, type = "link", newdata = dat.new, se.fit =T, re.form= NA)

pr <- predict(mod, newdata = dat.new, type = "link", se.fit = TRUE, re.form = NA)
#bind se's and fitted points
dat.new = cbind(dat.new, preds)
#inverse link function
ilink <- family(mod)$linkinv
#back transform CIs
dat.new <- transform(dat.new,
                     Fitted = ilink(fit),
                     Upper = ilink(fit + (2*se.fit)),
                     Lower = ilink(fit - (2*se.fit)))

ggplot(ti.ab, aes(x=dpi, y=elisa_od, color=treatment))+
  geom_jitter(shape = 16, alpha=0.35, size = 2, position=position_jitterdodge(jitter.height=0, jitter.width=0.1, dodge.width=0.5))+
  geom_point(data=dat.new, aes(x=dpi, y=yhat, color=treatment), shape = 16, size=2.5, position=position_dodge(width=0.5))+
  geom_errorbar(data=dat.new, aes(x=dpi, y=yhat, ymax=Upper, ymin=Lower, color=treatment),
                position=position_dodge(width=0.5), width=0.1)+
  labs(x="Days Post-Inoculation", y= "ELISA OD", color="Treatment")+
  scale_color_manual(values=c(inf_colors))+
  facet_wrap(~temp)+
  theme(
    strip.text = element_text(size = 12, face="bold")
  )


# tidy or convert results
coef_fixed   <- broom.mixed::tidy(lm1, effects = "fixed", conf.int = TRUE)
coef_random  <- broom.mixed::tidy(lm1, effects = "ran_pars", conf.int = TRUE)
fit_glance   <- broom.mixed::glance(lm1) 
anova_II_df  <- broom::tidy(lm1_anova_II)
anova_III_df <- broom::tidy(lm1_anova_III)
emm_df       <- as.data.frame(emm)
pairs_df     <- as.data.frame(pairs(emm, by = c("temp", "dpi"), adjust = "tukey"))


# combine into one named list
out_list <- list(
  "Model_Fixed_Coeffs"     = coef_fixed,
  "Model_Random_Vars"      = coef_random,
  "Model_Fit_Summary"      = fit_glance,
  "ANOVA_Type_II"          = anova_II_df,
  "ANOVA_Type_III"         = anova_III_df,
  "Estimated_Marginal_Means" = emm_df,
  "Pairwise_Comparisons"   = pairs_df
)

# write to one Excel file
#write_xlsx(out_list, "/Users/jesse/Documents/GitHub/Temperature-and-Immunity-2022/Results/Antibody_Results.xlsx")

####Fever Score####
source("dataCleaning_fever.R")
#fever_change = score - baseline
#fever_diff = change from previous score

#New thought: Look only at baseline versus peak fever score. This simplifies models a ton and still looks at magnitude of response
ti.f$dpi <- as.factor(ti.f$dpi)

#peak score
ti.f <- ti.f %>%
  group_by(band_number)%>%
  mutate(fever_peak = max(fever_score),
         fever_high = max(fever_high = max(fever_score[dpi != 0])))

ti.f$fever_peak  

#Lood at differences between baseline (DPI 0) and peak fever score regardless of when it occurs.
ti.fc <- ti.f %>%
  dplyr::select(band_number, treatment, temp, dpi, fever_score, fever_peak, mass, sex, ever_infected) %>%
  group_by(band_number, treatment, temp, sex, ever_infected) %>%
  summarise(
    baseline   = fever_score[dpi == 0],
    peak       = max(fever_score, na.rm = TRUE),
    high       = max(fever_high = max(fever_score[dpi != 0])),
    magnitude  = peak - baseline,
    mag_high   = high - baseline,
    .groups    = "drop"
  ) %>%
  mutate(groups = paste(temp, ifelse(treatment == "Inoculated", "Inoculated", "Control")))

ggplot(ti.fc, aes(x=groups, y = mag_high, fill=groups))+
  geom_boxplot(outlier.shape = 4)+
  scale_fill_manual(values=c(treat_colors))+
  geom_jitter(width=0.1, height=0, alpha=1, shape=1)+
  labs(x="Treatment Group", y="Magnitude of Fever Change (Peak - Baseline)", fill = "Treatment Group", alpha="Treatment Group")+
  theme(
    axis.text.x = element_text(angle=45, hjust=1),
    legend.position = "right",
    legend.direction = "vertical"
  )+
  facet_wrap(~sex)



ti.long <- pivot_longer(ti.fc,
                        cols = c(baseline, high),
                        names_to = "fever_type",
                        values_to = "fever_value")

#Ocular temperature baseline vs peak
ggplot(ti.long, aes(x=groups, y=fever_value, alpha=fever_type, fill=groups), color="black")+
  geom_boxplot(position=position_dodge(width= 1), aes(alpha=fever_type))+
  geom_jitter(position = position_jitterdodge(dodge.width = 1, jitter.width=0.2), shape=1)+
  scale_color_manual(values=c(treat_colors))+
  scale_fill_manual(values=c(treat_colors))+
  scale_alpha_manual(values = c(0.5, 1))+
  labs(x="Treatment", y="Ocular Temperature (Degrees C)", fill = "Treatment", alpha="Fever Type", color="Treatment")+
  theme(
    axis.text.x = element_text(angle=25, hjust=1, size=12),
    legend.position = "right",
    legend.direction = "vertical"
  )+
  facet_wrap(~sex)


ti.long$temp <- relevel(as.factor(ti.long$temp), ref = "Warm")
ti.long$treatment <- relevel(as.factor(ti.long$treatment), ref = "Sham")
ti.long$fever_type <- relevel(as.factor(ti.long$fever_type), ref = "baseline")

ti.long <- ti.long %>%
  mutate(fever_type = dplyr::recode(fever_type,
                             "baseline" = "Baseline",
                             "high" = "Peak"))

ti.long <- ti.long %>%
  mutate(sex = dplyr::recode(sex,
                                  "M" = "Male",
                                  "F" = "Female"))


ggplot(ti.long, aes(x=fever_value, fill=temp))+
  geom_histogram(position="dodge", binwidth=0.5, alpha=0.75)+
  facet_wrap(~fever_type)+
  scale_fill_manual(values=temp_colors)

#model selection
fc1 <- glmmTMB(fever_value ~ fever_type * temp * treatment + (1|band_number), data=ti.long)
fc1s <- glmmTMB(fever_value ~ fever_type * temp * treatment + sex + (1|band_number), data=ti.long)

fc2 <- glmmTMB(fever_value ~ fever_type + temp * treatment + (1|band_number), data=ti.long)
fc2s <- glmmTMB(fever_value ~ fever_type + temp * treatment + sex + (1|band_number), data=ti.long) 

fc3 <-glmmTMB(fever_value ~ fever_type + temp:treatment + (1|band_number), data=ti.long)
fc3s <-glmmTMB(fever_value ~ fever_type + temp:treatment + sex + (1|band_number), data=ti.long)

fc4 <-glmmTMB(fever_value ~ fever_type * temp + treatment + (1|band_number), data=ti.long)
fc4s <-glmmTMB(fever_value ~ fever_type * temp + treatment + sex + (1|band_number), data=ti.long)

fc5 <-glmmTMB(fever_value ~ fever_type * treatment + temp + (1|band_number), data=ti.long)
fc5s <-glmmTMB(fever_value ~ fever_type * treatment + temp + sex + (1|band_number), data=ti.long)
fc5si <-glmmTMB(fever_value ~ fever_type * treatment + temp * sex + (1|band_number), data=ti.long)
fc5si2 <-glmmTMB(fever_value ~ fever_type * treatment * sex + temp + (1|band_number), data=ti.long)
fc5si2i <-glmmTMB(fever_value ~ fever_type * treatment * sex * temp + (1|band_number), data=ti.long)

fc6 <-glmmTMB(fever_value ~ 1 + (1|band_number), data=ti.long)

aictab(cand.set=list(fc1, fc1s, fc2, fc2s, fc3, fc3s, fc4, fc4s, fc5, fc5s, fc5si, fc5si2, fc5si2i, fc6), 
       modnames=c("fc1",  "fc1s", "fc2",  "fc2s", "fc3", "fc3s", "fc4", "fc4s", "fc5", "fc5s",  "fc5si", "fc5si2", "fc5si2i", "fc6"))

#Does the interaction between fever type and treatment or temperature predict fever_value?
  #Is fever_value affected by the interaction between fever type and treatment, the interaction between temperature and sex,
  #or each main effect while accounting for repeated measures
lm2<- glmmTMB(fever_value ~ fever_type * treatment + temp * sex + (1|band_number), data=ti.long)

drop1(lm2, test = "Chisq")

simulateResiduals(lm2, plot=T)
summary(lm2)
plot(allEffects(lm2))

#Anova type II asks whether there is a main effect overall, averaged across the other variables
lm2_anova_II <- car::Anova(lm2, type = "II")

lm2_anova_III <- car::Anova(lm2, type = "III")

#emmeans
#estimate fever_type means within each combination of temperature x sex x treatment
emm <-emmeans(lm2, ~fever_type | temp*sex*treatment)
tests <- contrast(emm, method = "pairwise", adjust = "tukey")


emm_df <- as.data.frame(emm)  # has emmean, SE, lower.CL, upper.CL

lvl_f <- c("Cold Control", "Cold Inoculated", "Warm Control", "Warm Inoculated")

emm_df <- as.data.frame(emm) %>%
  mutate(
    groups = paste(temp, ifelse(treatment == "Inoculated", "Inoculated", "Control")),
    groups = factor(groups, levels = lvl_f)
  )


#Plot Emmeans
ggplot(ti.long, aes(x = groups,
                    y = fever_value,
                    color = groups,
                    shape = fever_type)) +
  
  geom_jitter(
    position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.9),
    alpha = 0.75, size = 2.5, stroke=1
  ) +
  # model means — DODGE BY fever_type
  geom_point(
    data = emm_df,
    aes(x = groups, y = emmean, group = fever_type, shape = fever_type, color=groups),
    position = position_dodge(width = 0.9),
    size = 2.5, 
    stroke = 1,
    color = "black",
    inherit.aes = FALSE
  ) +
  geom_errorbar(
    data = emm_df,
    aes(x = groups, ymin = lower.CL, ymax = upper.CL,
        group = fever_type, color=groups),
    width = 0.0,
    position = position_dodge(width = 0.9),
    color = "black",
    inherit.aes = FALSE
  ) +

  scale_shape_manual(values = c(1, 16)) +
  scale_color_manual(values = c(treat_colors))+
  labs(x = "Treatment Groups", y = "Ocular Temperature (C)",
       color = "Treatment Group", shape = "Fever Phase") +
  theme(
    axis.text.x = element_text(size=12, angle=45, hjust=1),
    strip.text = element_text(size=12)
  )+
  facet_wrap(~sex)

# tidy or convert results
lm2_coef_fixed   <- broom.mixed::tidy(lm2, effects = "fixed", conf.int = TRUE)
lm2_coef_random  <- broom.mixed::tidy(lm2, effects = "ran_pars", conf.int = TRUE)
lm2_fit_glance   <- broom.mixed::glance(lm2) 
lm2_anova_II_df  <- broom::tidy(lm2_anova_II)
lm2_anova_III_df <- broom::tidy(lm2_anova_III)
lm2_emm_df       <- as.data.frame(emm)
lm2_pairs_df     <- as.data.frame(pairs(emm, by = c("temp", "treatment"), adjust = "tukey"))


# combine into one named list
lm2_out_list <- list(
  "Model_Fixed_Coeffs"     = lm2_coef_fixed,
  "Model_Random_Vars"      = lm2_coef_random,
  "Model_Fit_Summary"      = lm2_fit_glance,
  "ANOVA_Type_II"          = lm2_anova_II_df,
  "ANOVA_Type_III"         = lm2_anova_III_df,
  "Estimated_Marginal_Means" = lm2_emm_df,
  "Pairwise_Comparisons"   = lm2_pairs_df
)

# write to one Excel file
#write_xlsx(lm2_out_list, "/Users/jesse/Documents/GitHub/Temperature-and-Immunity-2022/Results/Fever_Change_Results.xlsx")



bp_df <- ti.f %>%
  dplyr::select(band_number, treatment, temp, dpi, fever_score, sex) %>%
  group_by(band_number, treatment, temp, sex) %>%
  summarise(
    baseline   = fever_score[dpi == 0],
    peak       = max(fever_score, na.rm = TRUE),
    high       = max(fever_high = max(fever_score[dpi != 0])),
    magnitude  = peak - baseline,
    mag_high   = high - baseline,
    .groups  = "drop"
  ) %>%
  mutate(treat.temp = interaction(treatment, temp, drop = TRUE)) %>%
  pivot_longer(
    c(baseline, high), #Change peak to high; high = max after baseline > allows negative numbers; peak = max including baseline > anything negative will be 0
    names_to  = "phase",
    values_to = "fever_value"
  ) %>%
  drop_na(fever_value) %>%
  mutate(
    phase      = factor(phase, levels = c("baseline","high"), labels = c("Baseline","High")),
    treat.temp = factor(treat.temp,
                        levels = c("Sham.Cold", "Inoculated.Cold", "Sham.Warm", "Inoculated.Warm"))  # set your preferred order if needed
  )


#model selection
fm1 <-glmmTMB(fever_value ~ phase * treatment * temp * sex + (1|band_number),
              data = bp_df)

fm1.2 <-glmmTMB(fever_value ~ phase * treatment * temp + sex + (1|band_number),
              data = bp_df)
fm1.3 <-glmmTMB(fever_value ~ phase * treatment * temp + (1|band_number),
              data = bp_df)


fm2 <-glmmTMB(fever_value ~ phase * treatment + temp + (1|band_number),
              data = bp_df)
fm2.1 <-glmmTMB(fever_value ~ phase * treatment + temp * sex + (1|band_number),
              data = bp_df)

fm3 <-glmmTMB(fever_value ~ phase * treatment + (1|band_number),
              data = bp_df)
fm3.1 <-glmmTMB(fever_value ~ phase * treatment + sex + (1|band_number),
              data = bp_df)


fm4 <-glmmTMB(fever_value ~ 1 + (1|band_number),
              data = bp_df)


aictab(cand.set=list( fm1, fm1.2, fm1.3, fm2, fm2.1, fm3, fm3.1, fm4), 
       modnames=c("fm1", "fm1.2", "fm1.3", "fm2", "fm2.1", "fm3", "fm3.1", "fm4"))

simulateResiduals(fm2.1, plot=T)

m_phase <- glmmTMB(fever_value ~ phase * treatment + temp + sex + (1|band_number),
                   data = bp_df)
summary(m_phase)
simulateResiduals(m_phase, plot=T)
car::Anova(m_phase, type = "II")
emmeans(m_phase, pairwise ~ treatment | phase)
emmeans(m_phase, pairwise ~ temp | phase)
emmeans(m_phase, pairwise ~ sex | phase)

emm <- emmeans(m_phase, ~ phase | treatment + temp, type = "response")
emm_df <- as.data.frame(emm) %>%
  mutate(treat.temp = interaction(treatment, temp, drop = TRUE))

emm_df <- emm_df %>%
  mutate(
    phase = factor(phase, levels = c("Baseline","High")),
    treat.temp = factor(treat.temp, levels = levels(bp_df$treat.temp)),
    x_num = as.numeric(treat.temp),
    x_pos = x_num + ifelse(phase == "Baseline", -0.18, 0.18)
  )

seg_emm <- emm_df %>%
  dplyr::select(treat.temp, phase, x_pos, emmean) %>%
  pivot_wider(names_from = phase, values_from = c(x_pos, emmean)) %>%
  filter(!is.na(emmean_Baseline), !is.na(emmean_High))

nudge <- 0.18
pos_df <- bp_df %>%
  mutate(x_num = as.numeric(treat.temp),
         x_pos = x_num + ifelse(phase == "Baseline", -nudge, +nudge))

ggplot() +
  # raw paired lines per bird
  geom_segment(
    data = pos_df %>% 
      dplyr::select(band_number, treat.temp, phase, fever_value, x_pos) %>%
      pivot_wider(names_from = phase, values_from = c(fever_value, x_pos)) %>%
      filter(!is.na(fever_value_Baseline), !is.na(fever_value_High)),
    aes(x = x_pos_Baseline, xend = x_pos_High,
        y = fever_value_Baseline, yend = fever_value_High,
        color = treat.temp),
    alpha = 0.15
  ) +
  # raw points
  geom_point(
    data = pos_df,
    aes(x = x_pos, y = fever_value, color = treat.temp, shape = phase),
    size = 3, alpha = 1
  ) +
  # EMMEANS segment (Baseline -> Peak)
  geom_segment(
    data = seg_emm,
    aes(x = x_pos_Baseline, xend = x_pos_High,
        y = emmean_Baseline, yend = emmean_High),
    color="black",
    linewidth = 0.75
  ) +
  # EMMEANS points
  geom_point(
    data = emm_df,
    aes(x = x_pos, y = emmean, shape = phase),
    color="black",
    size = 3,
  ) +
  # EMMEANS CIs
  geom_errorbar(
    data = emm_df,
    aes(x = x_pos, ymin = lower.CL, ymax = upper.CL),
    color="black",
    width = 0.05, linewidth = 0.75
  ) +
  # axis & scales
  scale_x_continuous(
    breaks = sort(unique(pos_df$x_num)),
    labels = levels(pos_df$treat.temp)
  ) +
  scale_shape_manual(values = c(Baseline = 1, High = 16)) +
  scale_color_manual(values=c(treat_colors))+
  labs(x = "Treatment × Temp", y = "Ocular Temperature (C)",
       color = "Group", shape = "Phase") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


##Fever by dpi##
ggplot(ti.f, aes(x=dpi, y=fever_score, color=groups))+
  geom_point()

ti.f$dpi.f <- as.factor(ti.f$dpi)
ti.f$temp <- relevel(as.factor(ti.f$temp), ref = "Warm")
ti.f$treatment <- relevel(as.factor(ti.f$treatment), ref = "Sham")

#did fever differ between treatments at baseline?
t.test(fever_score ~ treatment, data=ti.f %>% filter(dpi.f == 0))

#Model Selection
#because sampled_first is a function of dpi (it is fixed for each dpi), it is collinear with dpi.f, so do not include.
lmf.a <- glmmTMB(fever_score ~ temp * treatment + sex + dpi.f + (1|band_number), data=ti.f)
lmf.b <- glmmTMB(fever_score ~ temp * treatment +  dpi.f +(1|band_number), data=ti.f)


lmf.c <- glmmTMB(fever_score ~ temp + treatment + sex + dpi.f +(1|band_number), data=ti.f)
lmf.cd <- glmmTMB(fever_score ~  sex + treatment  + temp * dpi.f +(1|band_number), data=ti.f)
lmf.ce <- glmmTMB(fever_score ~  sex + temp + treatment  * dpi.f +(1|band_number), data=ti.f)
lmf.cec <- glmmTMB(fever_score ~  sex + temp + treatment  * as.numeric(dpi) +(1|band_number), data=ti.f)
lmf.cei <- glmmTMB(fever_score ~  treatment + temp * sex + dpi.f +(1|band_number), data=ti.f)
lmf.cf <- glmmTMB(fever_score ~  sex + temp * treatment  * dpi.f +(1|band_number), data=ti.f)

lmf.d <- glmmTMB(fever_score ~ temp + treatment +  dpi.f +(1|band_number), data=ti.f)
lmf.dc <- glmmTMB(fever_score ~ temp + treatment +  as.numeric(dpi) +(1|band_number), data=ti.f)
lmf.dd <- glmmTMB(fever_score ~ treatment + temp *  dpi.f +(1|band_number), data=ti.f)

lmf.e <- glmmTMB(fever_score ~ temp + sex + dpi.f +(1|band_number), data=ti.f)
lmf.f <- glmmTMB(fever_score ~ treatment + sex + dpi.f +(1|band_number), data=ti.f)

lmf.g <- glmmTMB(fever_score ~ temp * treatment + dpi.f +(1|band_number), data=ti.f)
lmf.h <- glmmTMB(fever_score ~ temp * treatment + sex + dpi.f +(1|band_number), data=ti.f)
lmf.i <- glmmTMB(fever_score ~ 1 + (1|band_number), data=ti.f)

aictab(cand.set=list(lmf.a, lmf.b, lmf.c, lmf.cd, lmf.ce, lmf.cec, lmf.cei, lmf.cf, lmf.d, lmf.dc, lmf.dd, lmf.e, lmf.f, lmf.g, lmf.h, lmf.i), 
       modnames=c("lmf.a",  "lmf.b", "lmf.c", "lmf.cd", "lmf.ce", "lmf.cec", "lmf.cei", "lmf.cf", "lmf.d", "lmf.dc", "lmf.dd", "lmf.e", "lmf.f", "lmf.g", "lmf.h", "lmf.i"))

#second round
a <- glmmTMB(fever_score ~ temp + treatment + sex + dpi.f + (1|band_number), data=ti.f)
b <- glmmTMB(fever_score ~  treatment  + temp : dpi.f + sex : dpi.f + (1|band_number), data=ti.f)
c <- glmmTMB(fever_score ~  treatment : dpi.f + temp + sex : dpi.f + (1|band_number), data=ti.f)
c.5 <- glmmTMB(fever_score ~  treatment : dpi.f + temp : dpi.f + sex : dpi.f + (1|band_number), data=ti.f)
c.75 <- glmmTMB(fever_score ~ (treatment + temp + sex) * as.numeric(dpi.f) + (1|band_number), data=ti.f)
c.8 <- glmmTMB(fever_score ~  sex + temp : dpi.f + treatment  : dpi.f + (1|band_number), data=ti.f)
c.9 <- glmmTMB(fever_score ~ temp * dpi.f + treatment * dpi.f + sex + (1|band_number),
        data = ti.f)
c.10 <- glmmTMB(fever_score ~ temp * dpi.f + treatment * dpi.f + sex * dpi.f + (1|band_number),
                data = ti.f)
d <- glmmTMB(fever_score ~  treatment + temp + sex * dpi.f +(1|band_number), data=ti.f)
e <- glmmTMB(fever_score ~  temp * sex * dpi.f + (1|band_number), data=ti.f)
f <- glmmTMB(fever_score ~  treatment * sex * dpi.f +(1|band_number), data=ti.f)
null <- glmmTMB(fever_score ~ 1 + (1|band_number), data=ti.f)

aictab(cand.set=list(a, b, c, c.5, c.75, c.8, c.9, c.10, d, e, f, null), 
       modnames=c("a", "b", "c", "c.5", "c.75", "c.8", "c.9", "c.10", "d", "e", "f", "null"))


drop1(c.9)

anova(c.9, c.5, test = "Chisq")


simulateResiduals(c.9, plot=T)
hist(resid(c.9))
summary(c.9)
summary(c.8)


#glm_fever <- glmmTMB(fever_score ~  sex + temp + treatment  * dpi.f +(1|band_number), data=ti.f)
glm_fever <-   glmmTMB(fever_score ~ temp * dpi.f + treatment * dpi.f + sex + (1|band_number),
                       data = ti.f)

simulateResiduals(glm_fever, plot=T)

summary(glm_fever)
car::Anova(glm_fever, type = "III")
car::Anova(glm_fever, type = "II")

# Get emmeans (predicted marginal means)
emm_df <-  emmeans(glm_fever, ~ sex * treatment * dpi.f * temp, type = "response") %>%
  as.data.frame() %>%
  left_join(
    ti.f %>%
      distinct(temp, treatment, dpi.f, sex, groups),
    by = c("temp", "treatment", "dpi.f", "sex")
  ) %>%
  dplyr::mutate(groups = dplyr::recode(groups,
                                       "Warm Control" = "Warm Control",
                                       "Cold Control" = "Cold Control",
                                       "Cold Infected" = "Cold Inoculated",
                                       "Warm Infected" = "Warm Inoculated"))

emm_df <- emm_df %>%
  mutate(sex = dplyr::recode(sex,
                             "M" = "Male",
                             "F" = "Female"))

ti.f.g <- ti.f %>%
  dplyr::mutate(groups = dplyr::recode(groups,
                                       "Cold Infected" = "Cold Inoculated",
                                       "Warm Infected" = "Warm Inoculated"))

ti.f.g <- ti.f.g %>%
  dplyr::mutate(sex = dplyr::recode(sex,
                                       "M" = "Male",
                                       "F" = "Female"))

dodge <- position_dodge(width = .5)

#Model predictions by temp
ggplot() +
  geom_point(data = ti.f.g,
             aes(x = dpi.f, y = fever_score, color = groups),
             position =dodge,
             alpha = 0.5, size = 2) +
  # geom_line(data=ti.m.g,
  #           aes(x=dpi.f, y=mass, color=groups, group = as.factor(band_number)),
  #           position = dodge, alpha=0.1)+
  
  geom_errorbar(data = emm_df,
                aes(x = dpi.f, ymin = lower.CL, ymax = upper.CL, groups = treatment),
                color="black", position = dodge, width = 0.1) +
  geom_line(data = emm_df,
            aes(x = as.numeric(dpi.f), y = emmean, linetype = treatment),
            color="black", position = dodge, size =0.25) +
  
  geom_point(data = emm_df,
             aes(x = dpi.f, y = emmean, color = groups),
            position = dodge, size = 2.8) +
  geom_point(data = emm_df,
             aes(x = dpi.f, y = emmean, group=groups),
             position = dodge, size = 2.8,
             shape=1, color="black") +
  
  scale_color_manual(values = treat_colors) +
  scale_shape_manual(values = c(1, 1, 1, 1))+
  labs(x = "Days Post Inoculation", y = "Ocular Temperature (C)", color = "Treatment Group", shape="Treatment Group", linetype="Inoculation Type") +
  scale_linetype_manual(
    name   = "Inoculation Type",
    values = c("dashed", "solid"),
    labels = c("Control", "Inoculated")
  )+
  facet_grid(~temp~sex)+ #, scales = "free_y"
  theme(strip.text = element_text(size=12))



######Magnitude of fever change; probably don't use this one####
# ti.fc$temp <- relevel(as.factor(ti.fc$temp), ref = "Warm")
# ti.fc$treatment <- relevel(as.factor(ti.fc$treatment), ref = "Sham")
# 
# ti.fc$magnitude <- ti.fc$magnitude+0.001
# 
# #model selection
# #Magnitude does not allow for negative numbers as it subtracts baseline from the peak. If the peak was at baseline, then magnitude is 0
# fm1 <- glmmTMB(magnitude ~ temp*treatment, data=ti.fc, family=Gamma(link="log"))
# 
# fm2 <- glmmTMB(magnitude ~ temp+treatment, data=ti.fc, family=Gamma(link="log"))
# #inverse does not converge
# fm2.i <- glmmTMB(magnitude ~ temp+treatment, data=ti.fc, family=Gamma(link="inverse"))
# fm2.1 <- glmmTMB(magnitude ~ temp+treatment+sex, data=ti.fc, family=Gamma(link="log"))
# fm2.2 <- glmmTMB(magnitude ~ temp+treatment*sex, data=ti.fc, family=Gamma(link="log"))
# fm2.3 <- glmmTMB(magnitude ~ temp*sex+treatment, data=ti.fc, family=Gamma(link="log"))
# fm2.4 <- glmmTMB(magnitude ~ temp*treatment*sex, data=ti.fc, family=Gamma(link="log"))
# 
# fm3 <- glmmTMB(magnitude ~ treatment, data=ti.fc, family=Gamma(link="log"))
# fm4 <- glmmTMB(magnitude ~ temp, data=ti.fc, family=Gamma(link="log"))
# fm5 <- glmmTMB(magnitude ~ 1, data=ti.fc, family=Gamma(link="log"))
# 
# aictab(cand.set=list( fm1, fm2, fm2.i, fm2.1, fm2.2, fm2.3, fm2.4, fm3, fm4, fm5), 
#        modnames=c("fm1", "fm2", "fm2.i", "fm2.1", "fm2.2", "fm2.3", "fm2.4", "fm3", "fm4", "fm5"))
# 
# #Sex does not affect
# summary(fm2.1)
# anova(fm2, fm2.1, type="Chisq")
# 
# hist(ti.fc$magnitude)
# hist(resid(fm2))
# 
# #Does the magnitude of change in eye temperature between baseline and peak differ between treatments?
# lm3 <- glmmTMB(magnitude ~ temp+treatment, data=ti.fc, family=Gamma(link="log"))
# simulateResiduals(lm3, plot=T)
# summary(lm3)
# plot(allEffects(lm3))
# 
# #Anova type II asks whether there is a main effect overall, averaged across the other variables
# car::Anova(lm3, type = "II")
# 
# car::Anova(lm3, type = "III")
# 
# #emmeans
# emm <-emmeans(lm3, ~temp+treatment)
# tests <- contrast(emm, method = "pairwise", adjust = "tukey")
# 
# 
# emm_df <- as.data.frame(emm)  # has emmean, SE, lower.CL, upper.CL
#  
# emm_df <- as.data.frame(emm) %>%
#   mutate(
#     groups = paste(temp, ifelse(treatment == "Inoculated", "Inoculated", "Control")),
#     groups = factor(groups, levels = lvl_f)
#   )
# 
# 
# ggplot(ti.fc, aes(x=groups, y=magnitude, color = groups))+
#   geom_jitter(width = 0.1, height=0, shape =16, size=2)+
#   geom_point(data=emm_df, aes(y=emmean, x=groups))+
#   geom_errorbar(data=emm_df, aes(y=emmean, ymin=asymp.LCL, ymax=asymp.UCL, x=groups), width=0.1)+
#   labs(x="Treatment x Temp", y= "Magnitude of Fever Change (Peak - Baseline)", color="Treatment")+
#   scale_color_manual(values=c(treat_colors))+
#   theme(
#     axis.text.x = element_text(angle = 45, hjust=1),
#     legend.position = "right",
#     legend.direction = "vertical"
#   )
# 
# 
# mod <- lm3
# dat.new=expand.grid(magnitude=unique(ti.fc$magnitude),
#                     temp = unique(ti.fc$temp),
#                     treatment= unique(ti.fc$treatment))
# 
# dat.new$yhat=predict(mod, type="response", newdata = dat.new)
# #prediction intervals
# preds = predict(mod, type = "response", newdata = dat.new, se.fit =T)
# #bind se's and fitted points
# dat.new = cbind(dat.new, preds)
# #inverse link function
# ilink <- family(mod)$linkinv
# #back transform CIs
# dat.new <- transform(dat.new,
#                      Fitted = ilink(fit),
#                      Upper = ilink(fit + (2*se.fit)),
#                      Lower = ilink(fit - (2*se.fit)))

#######Try mag_high which is max ocular temperature excluding baseline to allow for negative numbers####
hist(ti.fc$mag_high)

fh1 <- lm(mag_high ~ temp*treatment, data=ti.fc)
fh2 <- lm(mag_high ~ temp+treatment, data=ti.fc)
fh2.1 <- lm(mag_high ~ temp+treatment+sex, data=ti.fc)
fh2.2 <- lm(mag_high ~ temp+treatment*sex, data=ti.fc)

fh3 <- lm(mag_high ~ treatment, data=ti.fc)
fh4 <- lm(mag_high ~ temp, data=ti.fc)
fh5 <- lm(mag_high ~ 1, data=ti.fc)

aictab(cand.set=list(fh1, fh2, fh2.1, fh2.2, fh3, fh4, fh5), 
       modnames=c("fh1", "fh2", "fh2.1", "fh2.2", "fh3", "fh4", "fh5"))

#no effect of sex
summary(fh2.1)
#Does the magnitude of change in eye temperature between baseline and peak differ between treatments?
lm3 <- lm(mag_high ~ temp+treatment, data=ti.fc)
simulateResiduals(lm3, plot=T)
summary(lm3)
plot(allEffects(lm3))

#Anova type II asks whether there is a main effect overall, averaged across the other variables
lm3_anova_II <- car::Anova(lm3, type = "II")

lm3_anova_III <- car::Anova(lm3, type = "III")

#emmeans
emm <-emmeans(lm3, ~temp+treatment)
tests <- contrast(emm, method = "pairwise", adjust = "tukey")


emm_df <- as.data.frame(emm)  # has emmean, SE, lower.CL, upper.CL

emm_df <- as.data.frame(emm) %>%
  mutate(
    groups = paste(temp, ifelse(treatment == "Inoculated", "Inoculated", "Control")),
    groups = factor(groups, levels = lvl)
  )


ggplot(ti.fc, aes(x=groups, y=mag_high, color = groups))+
  geom_jitter(width = 0.1, height=0, shape =16, size=2.5, alpha=0.75)+
  geom_point(data=emm_df, aes(y=emmean, x=groups), color="black", size=3)+
  geom_errorbar(data=emm_df, aes(y=emmean, ymin=lower.CL, ymax=upper.CL, x=groups), width=0., color="black", size=.75)+
  labs(x="Treatment Group", y= "Magnitude of Fever Change (Peak - Baseline [C])", color="Treatment Group")+
  scale_color_manual(values=c(treat_colors))+
  theme(
    axis.text.x = element_text(size=13, angle=45, hjust=1),
    legend.position = "right",
    legend.direction = "vertical",
  )

# tidy or convert results
lm3_coef_fixed   <- broom.mixed::tidy(lm3, effects = "fixed", conf.int = TRUE)
lm3_coef_random  <- broom.mixed::tidy(lm3, effects = "ran_pars", conf.int = TRUE)
lm3_fit_glance   <- broom.mixed::glance(lm3) 
lm3_anova_II_df  <- broom::tidy(lm3_anova_II)
lm3_anova_III_df <- broom::tidy(lm3_anova_III)
lm3_emm_df       <- as.data.frame(emm)
lm3_pairs_temp_df     <- as.data.frame(pairs(emm, by = c("temp"), adjust = "tukey"))
lm3_pairs_treatment_df     <- as.data.frame(pairs(emm, by = c("treatment"), adjust = "tukey"))

# combine into one named list
lm3_out_list <- list(
  "Model_Fixed_Coeffs"     = lm3_coef_fixed,
  "Model_Random_Vars"      = lm3_coef_random,
  "Model_Fit_Summary"      = lm3_fit_glance,
  "ANOVA_Type_II"          = lm3_anova_II_df,
  "ANOVA_Type_III"         = lm3_anova_III_df,
  "Estimated_Marginal_Means" = lm3_emm_df,
  "Pairwise_Comparisons_Temp"   = lm3_pairs_temp_df,
  "Pairwise_Comparisons_Treatment"   = lm3_pairs_treatment_df
)

# write to one Excel file
#write_xlsx(lm3_out_list, "/Users/jesse/Documents/GitHub/Temperature-and-Immunity-2022/Results/Fever_Magnitude_Results.xlsx")


#Model Predictions > lm only so different 
mod <- lm3
# build prediction grid from the same factor levels
dat.new <- ti.fc %>%
  distinct(temp, treatment, groups)

# predictions (confidence interval for the mean)
pr <- predict(lm3, newdata = dat.new, se.fit = TRUE)
dat.new <- dat.new %>%
  mutate(
    yhat = pr$se.fit*0 + pr$fit,     # just to keep the naming obvious
    Lower  = pr$fit - 1.96 * pr$se.fit,
    Upper  = pr$fit + 1.96 * pr$se.fit
  )

#Model Predictions
ggplot(ti.fc, aes(x=groups, y=mag_high, color = groups))+
  geom_jitter(width = 0.1, height=0, shape =16, size=2)+
  geom_point(data=dat.new, aes(y=yhat, x=groups), color="black")+
  geom_errorbar(data=dat.new, aes(y=yhat, ymin=Lower, ymax=Upper, x=groups), width=0.1, color="black")+
  labs(x="Treatment x Temp", y= "Magnitude of Fever Change (Peak - Baseline)", color="Treatment")+
  scale_color_manual(values=c(treat_colors))+
  theme(
    axis.text.x = element_text(angle = 45, hjust=1),
    legend.position = "right",
    legend.direction = "vertical"
  )

####Phagocytosis####
source("dataCleaning_phago.R")

##   Binomial model and the
#    weights represent the total number of trials, since the response is a
#    pre-calculated proportion. The interaction between temp and treatment was
#    not significant.

ti.p$temp <- as.factor(ti.p$temp)
ti.p$treatment <- as.factor(ti.p$treatment)

ti.p <- ti.p %>%
  mutate(
    temp      = relevel(temp, ref = "Warm"),        # or whichever baseline you want
    treatment = relevel(treatment, ref = "Sham") # change as needed
  )

glm4 <- glmmTMB(phago_score~temp+treatment + (1|band_number), 
                weights=wbc_total+phago_total, 
                data=ti.p, family="binomial")

simulateResiduals(glm4, plot = T)

summary(glm4)

plot(allEffects(glm4))

#Anova type II asks whether there is a main effect overall, averaged across the other variables
glm4_anova_II <- car::Anova(glm4, type = "II")

glm4_anova_III <- car::Anova(glm4, type = "III")


#emmeans
emm <- emmeans(glm4, ~ temp * treatment, type = "response", re.form = NA)
emm_df <- as.data.frame(emm) %>%
  mutate(
    groups = paste(temp, ifelse(treatment == "Inoculated", "Inoculated", "Control")),
    groups = factor(groups, levels = lvl)
  )

ti.p <- ti.p %>%
  mutate(
    groups = paste(temp, ifelse(treatment == "Inoculated", "Inoculated", "Control")),
    groups = factor(groups, levels = lvl)
  )

ggplot(emm_df, aes(x = groups, y = prob, color=groups))+
  geom_jitter(data = ti.p, aes(x = groups, y = phago_score, color = groups),
              width = 0.2, alpha = 0.75, size=2.5)+
  geom_point(size=3, color="black")+
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL, y=prob), width=0., color="black", size=0.75)+
  labs(x="Treatment", y="Phagocytosis Score", color="Treatment")+
  scale_color_manual(values=treat_colors)+
  theme(
    axis.text.x = element_text(size=13, angle=45, hjust=1)
  )
  #facet_wrap(~temp)

ggplot(ti.p, aes(x=groups, y=phago_score, color=groups))+
  geom_jitter()+
  facet_wrap(~infected)+
  scale_color_manual(values=treat_colors)+
  theme(
    axis.text.x = element_text(size=13, angle=45, hjust=1)
  )

# tidy or convert results
glm4_coef_fixed   <- broom.mixed::tidy(glm4, effects = "fixed", conf.int = TRUE)
glm4_coef_random  <- broom.mixed::tidy(glm4, effects = "ran_pars", conf.int = TRUE)
glm4_fit_glance   <- broom.mixed::glance(glm4) 
glm4_anova_II_df  <- broom::tidy(glm4_anova_II)
glm4_anova_III_df <- broom::tidy(glm4_anova_III)
glm4_emm_df       <- as.data.frame(emm)
glm4_pairs_df     <- as.data.frame(pairs(emm, by = c("temp"), adjust = "tukey"))

# combine into one named list
glm4_out_list <- list(
  "Model_Fixed_Coeffs"     = glm4_coef_fixed,
  "Model_Random_Vars"      = glm4_coef_random,
  "Model_Fit_Summary"      = glm4_fit_glance,
  "ANOVA_Type_II"          = glm4_anova_II_df,
  "ANOVA_Type_III"         = glm4_anova_III_df,
  "Estimated_Marginal_Means" = glm4_emm_df,
  "Pairwise_Comparisons"   = glm4_pairs_df
)

# write to one Excel file
#write_xlsx(glm4_out_list, "/Users/jesse/Documents/GitHub/Temperature-and-Immunity-2022/Results/Phagocytosis_Results.xlsx")

#####Eye Score####
source("dataCleaning_TI22.R")

range(ti.cont$quantity) #highest control quantity = 95.38
ti$quant_cutoff = 100
ti$seropos_cutoff = 0.061
ti$sympt_cutoff = 0.1

#Sample sizes
ti %>%
  filter(dpi ==0)%>%
  dplyr::select(treatment, temp, inf_temp)%>%
  tbl_summary(
  )%>%
  modify_header(
    label ~ "**All Birds**"
  )

ti %>%
  dplyr::select(dpi, treatment, temp, inf_temp)%>%
  tbl_summary(
    by=dpi
  )%>%
  modify_header(
    label ~ "**All Days**"
  )


#data frame with only inoculated birds on days with eye scores
ti.inoc <- ti%>%
  filter(treatment == "Inoculated" & dpi %in% c(0, 3,7,9,14,18,21,24,28,35))

hist(ti.inoc$total_eye_score)
hist(ti$total_eye_score)

ti.inoc$sympt_cutoff = 0.5
ti.inoc$quant_cutoff = 22.5

#diseased = eye score
#infected = eye score or pathology >= cutoff
ti.inoc <- ti.inoc %>%
  mutate(diseased= ifelse(total_eye_score>=sympt_cutoff, 1, 0),
    sick = ifelse(total_eye_score>=sympt_cutoff | quantity >= quant_cutoff, 1, 0),
    infected = ifelse(quantity >= quant_cutoff, 1, 0))

#new column with whether the bird is ever diseased
ti.inoc <- ti.inoc %>%
  group_by(band_number) %>%
  mutate(ever_diseased = ifelse(any(coalesce(diseased, 0) == 1), 1, 0),
         ever_infected = ifelse(any(coalesce(infected, 0) == 1), 1, 0),
         ever_sick = ifelse(any(coalesce(sick, 0) == 1), 1,0))%>%
  ungroup()
  

#Sample sizes
  #ever_diseased: cold n=6, warm n=8
  #ever_infected: cold n=12, warm n=13
    #Asymptomatic: cold n=6, warm n=5
ti.inoc %>%
  filter(dpi ==28)%>%
  dplyr::select(temp, ever_diseased, ever_infected, ever_sick)%>%
  tbl_summary(
    by=temp
  )%>%
  modify_header(
    label ~ "**Inoculated Birds**"
  )

#Source Eye Score Formatting
source("dataCleaning_eyeScore.R")
#Ordinal Mixed Effects Model
#Sham birds did not get infected therefore did not have eyescores. Remove from analysis so that model converges.
  #No eye score on DPI 3 so also remove. Keep DPI 0 as reference category.
ti.mod <- ti %>%
  dplyr::select(date, dpi, band_number, bird_ID, treatment, temp, total_eye_score)%>%
  filter(treatment != "Sham" & dpi %in% c(0, 7,9,14,18,21,24,28))%>%
  mutate(temp = relevel(factor(temp), ref = "Warm"),
         dpi_f  = factor(dpi, levels = c(0,7,9,14,18,21,24,28)))

#We need to define thresholds based off of the observed data. The technical bounds are [0,6], but we only have [0,4]
levels_obs <- sort(unique(ti.mod$total_eye_score))

#turn observed numeric scores into an ordered factor with factors
ti.mod$eye_score_ord <- ordered(
  ti.mod$total_eye_score,
  levels=levels_obs
)

with(ti.mod, xtabs(~ eye_score_ord + temp + dpi_f))

ti.mod <- ti.mod %>%
  mutate(eye5 = fct_collapse(eye_score_ord,
                             "0"="0",
                             "1"=c("0.5","1"),
                             "2"=c("1.5","2"),
                             "3"=c("2.5","3"),
                             "4"=c("3.5","4")))


ti.mod <- ti.mod %>%
  mutate(
    dpi_f = relevel(dpi_f, ref = "0"),
    temp  = relevel(temp,  ref = "Warm")
  )

#Model Comparisons
clm.a <- clmm(eye_score_ord ~ temp * dpi_f + (1|band_number),
              data=ti.mod,
              link = "logit")

clm.b <- clmm(eye5 ~ temp * dpi_f + (1|band_number),
              data=ti.mod,
              link = "logit")

clm.c <- clmm(eye_score_ord ~ temp + dpi_f + (1|band_number),
              data=ti.mod,
              link = "logit")

clm.d <- clmm(eye5 ~ temp + dpi_f + (1|band_number),
              data=ti.mod,
              link = "logit")

summary(clm.a)
summary(clm.b)
summary(clm.c)
summary(clm.d)

aictab(cand.set=list(clm.a, clm.c),
      modnames= c("clm.a", "clm.c"))

aictab(cand.set=list(clm.b, clm.d),
       modnames= c("clm.b", "clm.d"))

#Does collapsing eye score into 5 categories hurt model?

#can't directly compare, but lower AIC = good
AIC(clm.c, clm.d)

#Lower logLikelihood = good
logLik(clm.c)
logLik(clm.d)

#R2 same if not slightly better = good
r2_nakagawa(clm.c)
r2_nakagawa(clm.d)

#predictions look the same
emmip(clm.c, temp ~ dpi_f)
emmip(clm.d, temp ~ dpi_f)

#use eye5 in final model
clm1 <- clmm(eye5 ~ temp + dpi_f + (1|band_number),
              data=ti.mod,
              link = "logit")

summary(clm1)

ti.mod <- ti.mod %>%
  mutate(
    eye5_num = as.numeric(as.character(eye5)),               # 0..4
    dpi_num  = as.numeric(as.character(dpi_f))               # e.g., 0,7,9,14,...
  )

# 2) Model predictions: expected category (mean class) + 95% CIs
emm_mean <- emmeans(clm1, ~ temp * dpi_f,
                    type = "response", mode = "mean.class") %>%
  as.data.frame() %>%
  mutate(
    yhat = estimate - 1,
    ymin = asymp.LCL - 1,
    ymax = asymp.UCL - 1,
    dpi_num = as.numeric(as.character(dpi_f))
  )

ggplot() +
  # raw observations
  geom_jitter(
    data = ti.mod,
    aes(x = dpi_num, y = eye5_num, color = temp),
    width = 0.6, height = 0.08, alpha = 0.35, size = 1.7
  ) +
  # model mean class with CI
  geom_errorbar(
    data = emm_mean,
    aes(x = dpi_num, ymin = ymin, ymax = ymax, color = temp),
    width = 0.8, size = 0.5
  ) +
  geom_line(
    data = emm_mean,
    aes(x = dpi_num, y = yhat, color = temp, group = temp),
    size = 0.7
  ) +
  geom_point(
    data = emm_mean,
    aes(x = dpi_num, y = yhat, color = temp),
    size = 2.5
  ) +
  scale_y_continuous(breaks = 0:4, limits = c(-1, 4.1)) +
  labs(x = "Days Post Inoculation", y = "Total Eye Score",
       color = "Temperature") +
  theme_minimal(base_size = 12)


ggplot(ti.inoc, aes(x=dpi, y=total_eye_score, color=temp))+
  geom_jitter(width=0.5, height=0.1, alpha=0.25)+
  #geom_line(aes(group = band_number), alpha=0.25)+
  stat_summary(aes(x=dpi, y=total_eye_score, color = temp), geom = "point", size=3)+
  stat_summary(aes(x=dpi, y=total_eye_score, color = temp), geom = "errorbar", width=0.25)+
  stat_summary(aes(x=dpi, y=total_eye_score, color = temp), geom = "line", width=0.25, alpha=0.5)+
  scale_color_manual(values = temp_colors)+
  labs(x="Days Post Inoculation", y="Total Eye Score", color="Temperature")




#GAM
mod.e <- ti.inoc %>%
  dplyr::filter(treatment == "Inoculated" & dpi > 0)%>%
  dplyr::mutate(total_eye_score1 = total_eye_score+0.001, #add small constant for gamma model
                dpi = as.numeric(dpi),
                temp = factor(temp),
                band_number = factor(band_number)
  )
                
#Use GAM to ask how temperature affects eye score
# Fit a GAM with smooth function of dpi, separate smooths for each temperature
#7 dpi, so need k to be less than 7

##Gamma Model - Old
# gam1 <- gam(
#   total_eye_score1 ~ temp + s(dpi, by = temp, k = 5, bs = "cr") +
#     s(band_number, bs="re"),
#   data   = mod.e,
#   family = Gamma(link = "log"),
#   method = "REML",
#   select = TRUE #Select = TRUE tells the model to add extra penalty to each smooth term; if smoothing doesn't help, the model drops it
# )
# 
# 
# summary(gam1)
# gam.check(gam1)   # look at k-index & residuals
# plot(gam1, pages = 1, shade = TRUE)
# 
# #Test whether temperature smooth improves fit - null model
# gam_null <- gam(
#   total_eye_score1 ~ temp + s(dpi, k = 5, bs = "cr") + s(band_number, bs="re"),
#   data   = mod.e,
#   family = Gamma(link="log"),
#   method = "REML",
#   select = TRUE
# )
# 
# #Does including temperature  smooth improve model fit?
# anova(gam_null, gam1, test = "Chisq")
# 
# newdat <- expand.grid(
#   dpi  = seq(min(mod.e$dpi, na.rm=TRUE), max(mod.e$dpi, na.rm=TRUE), length.out = 200),
#   temp = levels(mod.e$temp)
# )
# 
# # dummy level (any level from the fit is fine)
# newdat$band_number <- levels(mod.e$band_number)[1]
# 
# # exclude the RE smooth when predicting
# pr <- predict(gam1, newdat, type = "link", se.fit = TRUE, 
#               exclude = "s(band_number)")
# 
# pr   <- predict(gam1, newdat, type = "link", se.fit = TRUE)
# newdat$fit <- family(gam1)$linkinv(pr$fit)
# newdat$lwr <- family(gam1)$linkinv(pr$fit - 1.96*pr$se.fit)
# newdat$upr <- family(gam1)$linkinv(pr$fit + 1.96*pr$se.fit)
# 
# ggplot(newdat, aes(dpi, fit, color = temp, fill = temp)) +
#   geom_jitter(data=mod.e, aes(x=dpi, y=total_eye_score, color= temp), height=0.05, width=0)+
#   geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.15, linewidth = 0) +
#   geom_line(linewidth = 1.1) +
#   scale_fill_manual(values=c(temp_colors))+
#   scale_color_manual(values=c(temp_colors))+
#   labs(x = "Days Post Infection", y = "Total Eye Score", fill= "Temperature", color = "Temperature")+
#   theme_bw()

#tweedie GAM to handle zero inflation instead of adding small constant
gam1 <- gam(
  total_eye_score ~ temp + s(dpi, by = temp, k = 5, bs = "cr") +
    s(band_number, bs = "re"),
  data = mod.e,
  family = tw(link = "log"),
  method = "REML"
)

summary(gam1)
gam.check(gam1)   # look at k-index & residuals
plot(gam1, pages = 1, shade = TRUE)

#Test whether temperature smooth improves fit - null model
gam_null <- gam(
  total_eye_score ~ temp + s(dpi, k = 5, bs = "cr") + s(band_number, bs="re"),
  data   = mod.e,
  family = tw(link = "log"),
  method = "REML",
  select = TRUE
)

#Does including temperature  smooth improve model fit?
anova(gam_null, gam1, test = "Chisq")

newdat <- expand.grid(
  dpi  = seq(min(mod.e$dpi, na.rm=TRUE), max(mod.e$dpi, na.rm=TRUE), length.out = 200),
  temp = levels(mod.e$temp)
)

# dummy level (any level from the fit is fine)
newdat$band_number <- levels(mod.e$band_number)[1]

# exclude the RE smooth when predicting
pr <- predict(gam1, newdat, type = "link", se.fit = TRUE, 
              exclude = "s(band_number)")

pr   <- predict(gam1, newdat, type = "link", se.fit = TRUE)
newdat$fit <- family(gam1)$linkinv(pr$fit)
newdat$lwr <- family(gam1)$linkinv(pr$fit - 1.96*pr$se.fit)
newdat$upr <- family(gam1)$linkinv(pr$fit + 1.96*pr$se.fit)

ggplot(newdat, aes(dpi, fit, color = temp, fill = temp)) +
  geom_jitter(data=mod.e, aes(x=dpi, y=total_eye_score, color= temp), height=0.05, width=0.5, size=2, alpha=0.5)+
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.15, linewidth = 0) +
  geom_line(linewidth = 1.1) +
  scale_fill_manual(values=c(temp_colors))+
  scale_color_manual(values=c(temp_colors))+
  labs(x = "Days Post Infection", y = "Total Eye Score", fill= "Temperature", color = "Temperature")+
  theme_bw()

#Save statistics
s <- summary(gam1)

# Parametric (fixed) terms
param_df <- as_tibble(s$p.table, rownames = "term") %>%
  rename(
    estimate  = `Estimate`,
    std_error = `Std. Error`
  ) %>%
  mutate(
    statistic = `t value`,
    p_value   = `Pr(>|t|)`,
    model     = "Tweedie_GAM"
  ) %>%
  dplyr::select(model, term, estimate, std_error, statistic, p_value)

# Smooth terms 
smooth_df <- as_tibble(s$s.table, rownames = "smooth") %>%
  mutate(
    edf    = edf,
    ref_df = `Ref.df`,
    stat   = `F`,
    p_value = `p-value`,
    model  = "Tweedie_GAM"
  ) %>%
  dplyr::select(model, smooth, edf, ref_df, stat, p_value)

# Model-level fit metrics
fit_df <- tibble(
  model         = "Tweedie_GAM",
  family        = gam1$family$family,
  link          = gam1$family$link,
  shape_p       = if (!is.null(gam1$family$variance)) s$family$p else NA_real_,  # Tweedie p
  AIC           = AIC(gam1),
  dev_expl_frac = s$dev.expl,                      # 0–1
  dev_expl_pct  = round(100 * s$dev.expl, 1),
  r2_adj        = s$r.sq %||% NA_real_,
  scale_est     = s$scale,
  n             = s$n
)

cmp_tbl <- anova(gam_null, gam1, test = "Chisq") %>%
  as_tibble() %>%
  mutate(comparison = "Tweedie: null vs temp-specific")

#write to excel
# write_xlsx(
#   list(
#     tweedie_param   = param_df,
#     tweedie_smooth  = smooth_df,
#     tweedie_fit     = fit_df,
#     tweedie_compare = cmp_tbl
#   ),
#   path = "/Users/jesse/Documents/GitHub/Temperature-and-Immunity-2022/Results/Eye_Score_tweedie_gam_Results.xlsx"
# )

ti.inoc$dpi.f <- as.factor(ti.inoc$dpi)
##Zero-Inflated Gamma Model##
lmes <- glmmTMB(total_eye_score ~ temp*dpi + (1|band_number),
               data=ti.inoc, 
               ziformula = ~ temp,
               family = ziGamma(link = "log"))

#zero-inflated model allowing for temperature to dictate the probability of zero.
#for the non-zero birds, we used a gamma model with a log link function, including a two way interaction between temp and dpi and all lower order effects

simulateResiduals(lmes, plot=T)
summary(lmes)
car::Anova(lmes, type = "III")

# Pairwise comparisons (Tukey-adjusted)
dpi_vals <- sort(unique(ti.inoc$dpi))

# 1) EMMs for temp at each dpi (population-level: drop REs)
emm <- emmeans(
  lmes,
  specs = ~ temp | dpi,              # simple effects of temp by dpi
  at    = list(dpi = dpi_vals),      # evaluate at each observed dpi
  type  = "response",
  re.form = NA
)

library(multcomp)
# Extract the emmeans table (not the contrasts)
means <- cld(emmeans_results$emmeans,
             adjust = "tukey",    # adjust letters for multiple comparisons
             Letters = letters,   # use a,b,c... labeling
             alpha = 0.05)        # significance threshold

# Convert to data frame
means_df <- as.data.frame(means)

head(means_df)

ggplot(means_df, aes(x = dpi, y = response, color = fct_rev(temp)))+
  geom_jitter(data = ti.inoc, aes(x = dpi, y = total_eye_score), 
              width = 0.1, alpha = 0.5)+
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),     
                width = 0.1)+
  stat_summary(data= ti.inoc, aes(x=dpi, y=total_eye_score, group= temp), geom="line", fun = "mean")+
  #geom_line(aes(group = temp), linewidth = 1.2)+
  geom_point(size = 3)+
  theme_bw()#+
facet_wrap(~temp)


#Max Eye Score
ti.mod.m <- ti.mod %>%
  filter(dpi > 0)%>%
  group_by(band_number, temp, treatment) %>%
  reframe(
    max_tes = max(total_eye_score))

ggplot(ti.mod.m, aes(x=max_tes, fill=temp))+
  geom_histogram(aes(groups=temp), position = "dodge")

wilcox.test(max_tes ~ temp, data = ti.mod.m)

# data:  max_tes by temp
# W = 108.5, p-value = 0.6225

ggplot(ti.mod.m, aes(x = temp, y = max_tes, color = temp, fill=temp)) +
  geom_jitter(width=0.1, height=0, size = 2.5, alpha=0.75) +
  #geom_boxplot(alpha=0.3, color="black", width=0.05)+
  stat_summary(geom="point", fun="mean", shape = 16, size=3, color="black")+
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0., color = "black") +
  # geom_point(data=dat.new, aes(x=temp, y=yhat), size=3)+
  # geom_errorbar(data=dat.new, aes(x=temp, y=yhat, ymin=Lower, ymax = Upper), width=0.05)+
  scale_color_manual(values = temp_colors) +
  scale_fill_manual(values = temp_colors) +
  labs(x = "Group", y = "Max Eyescore", color="Temperature", fill="Temperature")


####Pathogen Load####
source("dataCleaning_pathogenLoad.R")

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

ti.q$quantity1 <- ti.q$quantity + 1
hist(ti.q$quantity1)


dat.q <- ti.q %>%
  filter(treatment == "Inoculated" & dpi > 0) %>%
  mutate(
    dpi  = as.numeric(dpi),
    temp = factor(temp),
    band_number = factor(band_number)
  )


#Use GAM to ask how temperature affects pathogen load
# Fit a GAM with smooth function of dpi, separate smooths for each temperature
  #7 dpi, so need k to be less than 7
gam2 <- gam(
  quantity1 ~ temp + s(dpi, by = temp, k = 5, bs = "cr") +
    s(band_number, bs="re"),
  data   = dat.q,
  family = Gamma(link = "log"),
  method = "REML",
  select = TRUE #Select = TRUE tells the model to add extra penalty to each smooth term; if smoothing doesn't help, the model drops it
)


summary(gam2)
gam.check(gam2)   # look at k-index & residuals
plot(gam2, pages = 1, shade = TRUE)

newdat <- expand.grid(
  dpi  = seq(min(dat.q$dpi, na.rm=TRUE), max(dat.q$dpi, na.rm=TRUE), length.out = 200),
  temp = levels(dat.q$temp)
)

# dummy level (any level from the fit is fine)
newdat$band_number <- levels(dat.q$band_number)[1]

# exclude the RE smooth when predicting
pr <- predict(gam2, newdat, type = "link", se.fit = TRUE, 
              exclude = "s(band_number)")

pr   <- predict(gam2, newdat, type = "link", se.fit = TRUE)
newdat$fit <- family(gam2)$linkinv(pr$fit)
newdat$lwr <- family(gam2)$linkinv(pr$fit - 1.96*pr$se.fit)
newdat$upr <- family(gam2)$linkinv(pr$fit + 1.96*pr$se.fit)

ggplot(newdat, aes(dpi, fit, color = temp, fill = temp)) +
  geom_jitter(data=dat.q, aes(x=dpi, y=quantity1, color= temp),height=0.05, width=0.5, size=2, alpha=0.5)+
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.15, linewidth = 0) +
  geom_line(linewidth = 1.1) +
  scale_fill_manual(values=c(temp_colors))+
  scale_color_manual(values=c(temp_colors))+
  labs(x = "Days Post Infection", y = "Pathogen Load", fill= "Temperature", color = "Temperature")+
  scale_y_log10()
#Test whether temperature smooth improves fit - null model
gam_null <- gam(
  quantity1 ~ temp + s(dpi, k = 5, bs = "cr") + s(band_number, bs="re"),
  data   = dat.q,
  family = Gamma(link="log"),
  method = "REML",
  select = TRUE
)

#Does including temperature  smooth improve model fit?
anova(gam_null, gam2, test = "Chisq")


s <- summary(gam2)

# Parametric (fixed) terms
param_df <- as_tibble(s$p.table, rownames = "term") %>%
  rename(
    estimate  = `Estimate`,
    std_error = `Std. Error`
  ) %>%
  mutate(
    statistic = `t value`,
    p_value   = `Pr(>|t|)`,
    model     = "Gamma_GAM"
  ) %>%
  dplyr::select(model, term, estimate, std_error, statistic, p_value)

# Smooth terms 
smooth_df <- as_tibble(s$s.table, rownames = "smooth") %>%
  mutate(
    edf    = edf,
    ref_df = `Ref.df`,
    stat   = `F`,
    p_value = `p-value`,
    model  = "Gamma_GAM"
  ) %>%
  dplyr::select(model, smooth, edf, ref_df, stat, p_value)

# Model-level fit metrics
fit_df <- tibble(
  model         = "Gamma_GAM",
  family        = gam2$family$family,
  link          = gam2$family$link,
  shape_p       = if (!is.null(gam2$family$variance)) s$family$p else NA_real_,  # Tweedie p
  AIC           = AIC(gam2),
  dev_expl_frac = s$dev.expl,                      # 0–1
  dev_expl_pct  = round(100 * s$dev.expl, 1),
  r2_adj        = s$r.sq %||% NA_real_,
  scale_est     = s$scale,
  n             = s$n
)

cmp_tbl <- anova(gam_null, gam2, test = "Chisq") %>%
  as_tibble() %>%
  mutate(comparison = "Gamma: null vs temp-specific")

#write to excel
# write_xlsx(
#   list(
#     gamma_param   = param_df,
#     gamma_smooth  = smooth_df,
#     gamma_fit     = fit_df,
#     gamma_compare = cmp_tbl
#   ),
#   path = "/Users/jesse/Documents/GitHub/Temperature-and-Immunity-2022/Results/Path_Load_gamma_gam_Results.xlsx"
# )


####Tolerance####
ti.t.max <- ti %>%
  dplyr::select(band_number, dpi, tes, treatment, temp, mass, sex, quantity, log10_quantity, groups, ever_diseased, ever_infected) %>%
  filter(treatment == "Inoculated")%>%
  group_by(band_number) %>%
  summarise(
    max_tes = if (all(is.na(tes))) NA_real_ else max(tes, na.rm = TRUE),
    max_quantity1 = if (all(is.na(quantity))) NA_real_ else max(quantity+1, na.rm = TRUE),
    max_mass = if (all(is.na(mass))) NA_real_ else max(mass, na.rm = TRUE),
    tolerance = -(max_tes / log10(max_quantity1)),
    treatment = first(treatment),
    temp = first(temp),
    mass = first(mass),
    sex = first(sex),
    ever_diseased = first(ever_diseased),
    ever_infected = first(ever_infected),
    groups = first(groups)
  ) %>%
  ungroup()

ggplot(ti.t.max, aes(x=log10(max_quantity1), y=max_tes, color=temp))+
  geom_point()+
  facet_wrap(~sex)

#tissue-specific tolerance
ggplot(ti.t.max, aes(x=temp, y=tolerance, color=temp))+
  geom_boxplot(width=0.5, aes(fill=temp), alpha=0.1, color="black")+
  geom_jitter(size=3, width=0.25, alpha=0.75)+
  scale_color_manual(values=temp_colors)+
  scale_fill_manual(values=temp_colors)+
  labs(x="Temperature", y="Tolerance: \n -(Max Pathology / Max Pathogen Load)", color="Temperature", fill="Temperature")+
  facet_wrap(~sex, labeller = labeller(sex = c('F'="Female", 'M' = "Male")))+
  theme(strip.text = element_text(size = 14))

ti.wx <- ti.t.max %>%
  filter(ever_infected == 1)

#tissue-specific tolerance - only infected
ggplot(ti.wx, aes(x=temp, y=tolerance, color=temp))+
  geom_boxplot(width=0.25, aes(fill=temp), alpha=0.1, color="black")+
  geom_jitter(size=3, width=0.25, alpha=0.75)+
  scale_color_manual(values=temp_colors)+
  scale_fill_manual(values=temp_colors)+
  labs(x="Temperature", y="Tolerance: \n -(Max Pathology / Max Pathogen Load)", color="Temperature", fill="Temperature")+
  facet_wrap(~sex)+
  theme(strip.text = element_text(size = 14))

#Test: Does temperature affect tolerance in birds that got infected?
wilcox.test(tolerance ~ temp, data=ti.wx)
wilcox.test(tolerance ~ sex, data=ti.t.max)

lmt <- lm(tolerance ~ temp + sex, data= ti.t.max)
simulateResiduals(lmt, plot = T)
summary(lmt)

# data:  tolerance by temp
# W = 84, p-value = 0.7545

lmg <- lm(tolerance ~ temp + sex, data=ti.wx)

simulateResiduals(lmg, plot=T)
summary(lmg)

#Brown-Forsythe; is tolerance equally variable across temperatures?
bf_tol <- leveneTest(tolerance ~ temp, data = ti.wx, center = median)

# Levene's Test for Homogeneity of Variance (center = median)
#       Df F value Pr(>F)
# group  1  0.1664 0.6871
#       23              


####Mass####
#Sample sizes
ti.t.max %>%
  dplyr::select(treatment, temp, groups)%>%
  tbl_summary(
  )%>%
  modify_header(
    label ~ "**All Birds**"
  )

ti.m <- ti %>%
  filter(dpi %in% c(-28, -12, 3, 14, 21, 28, 35))

ggplot(data=ti.m, aes(x=dpi, y=mass, color=temp))+
  geom_point()+
  geom_line(aes(group = as.factor(band_number)), size=0.1)


ggplot(ti.m, aes(x=dpi, y=mass, color=temp))+
  geom_jitter()+
  facet_wrap(~sampled_first)

hist(ti.m$mass)

ti.m$dpi.f <- as.factor(ti.m$dpi)
ti.m$temp <- relevel(as.factor(ti.m$temp), ref = "Warm")
ti.m$treatment <- relevel(as.factor(ti.m$treatment), ref = "Sham")

#did mass differ at baseline?
t.test(mass ~ temp, data=ti.m %>% filter(dpi.f == -28))

t.test(mass ~ temp, data=ti.m %>% filter(dpi.f == -12))

#Model Selection
#because sampled_first is a function of dpi (it is fixed for each dpi), it is collinear with dpi.f, so do not include.
lm.a <- glmmTMB(mass ~ temp * treatment + sex + dpi.f + (1|band_number), data=ti.m)
lm.b <- glmmTMB(mass ~ temp * treatment +  dpi.f +(1|band_number), data=ti.m)


lm.c <- glmmTMB(mass ~ temp + treatment + sex + dpi.f +(1|band_number), data=ti.m)
lm.d <- glmmTMB(mass ~ temp + treatment +  dpi.f +(1|band_number), data=ti.m)

lm.dc <- glmmTMB(mass ~ temp + treatment +  dpi +(1|band_number), data=ti.m)
lm.dd <- glmmTMB(mass ~ treatment + temp *  dpi.f +(1|band_number), data=ti.m)

lm.e <- glmmTMB(mass ~ temp + sex + dpi.f +(1|band_number), data=ti.m)
lm.f <- glmmTMB(mass ~ treatment + sex + dpi.f +(1|band_number), data=ti.m)

lm.g <- glmmTMB(mass ~ temp * treatment + dpi.f +(1|band_number), data=ti.m)
lm.h <- glmmTMB(mass ~ temp * treatment + sex + dpi.f +(1|band_number), data=ti.m)
lm.i <- glmmTMB(mass ~ 1 + (1|band_number), data=ti.m)

aictab(cand.set=list(lm.a, lm.b, lm.c,  lm.d, lm.dc, lm.dd, lm.e, lm.f, lm.g, lm.h, lm.i), 
       modnames=c("lm.a",  "lm.b", "lm.c", "lm.d", "lm.dc", "lm.dd", "lm.e", "lm.f", "lm.g", "lm.h", "lm.i"))

lm.dd <- glmmTMB(mass ~ treatment + temp *  dpi.f + (1|band_number), data=ti.m)
lm.de <- glmmTMB(mass ~ treatment * dpi.f + temp *  dpi.f +(1|band_number), data=ti.m)
lm.df <- glmmTMB(mass ~ treatment * dpi.f + temp *  dpi.f + sex + (1|band_number), data=ti.m)
lm.dg <- glmmTMB(mass ~ treatment * dpi.f + temp *  dpi.f + sex *dpi.f + (1|band_number), data=ti.m)
lm.dh <- glmmTMB(mass ~ treatment + temp *  dpi.f + sex + (1|band_number), data=ti.m)

aictab(cand.set=list(lm.dd, lm.de, lm.df, lm.dg, lm.dh), 
       modnames=c("lm.dd", "lm.de", "lm.df", "lm.dg", "lm.dh"))

#sex not significant
summary(lm.dh)

#lm.dd best supported by AICc
simulateResiduals(lm.dd, plot=T)
summary(lm.dd)
plot(allEffects(lm.dd))
hist(resid(lm.dd))

simulateResiduals(lm.e, plot=T)
summary(lm.e)

lm7 <- glmmTMB(mass ~ treatment + temp * dpi.f +(1|band_number), data=ti.m)

DHARMa::simulateResiduals(lm7, plot=T)

summary(lm7)
plot(allEffects(lm7))

hist(resid(lm7))

car::Anova(lm7, type = "III")
car::Anova(lm7, type = "II")

#What if we compare only baseline after acclimization to temperatures (dpi -12 instead of -28)
lm7.no.bl <- glmmTMB(mass ~ treatment + temp * dpi.f +(1|band_number), data=ti.m %>% filter(dpi != -28))

summary(lm7.no.bl)

emm_df <- emmeans(lm7, ~ temp * treatment * dpi.f, re.form = NA) %>%
  as.data.frame()

emm_df <- emmeans(lm7, ~ temp * treatment * dpi.f, re.form = NA) %>%
  as.data.frame() %>%
  left_join(
    ti.m %>%
      distinct(temp, treatment, dpi.f, groups),
    by = c("temp", "treatment", "dpi.f")
  ) %>%
  dplyr::mutate(groups = dplyr::recode(groups,
                         "Cold Infected" = "Cold Inoculated",
                         "Warm Infected" = "Warm Inoculated"))


ti.m.g <- ti.m %>%
  dplyr::mutate(groups = dplyr::recode(groups,
                                       "Cold Infected" = "Cold Inoculated",
                                       "Warm Infected" = "Warm Inoculated"))

dodge <- position_dodge(width = .5)

#Model predictions by temp
ggplot() +
  geom_point(data = ti.m.g,
              aes(x = dpi.f, y = mass, color = groups),
              position =dodge,
              alpha = 0.25, size = 2) +
  # geom_line(data=ti.m.g,
  #           aes(x=dpi.f, y=mass, color=groups, group = as.factor(band_number)),
  #           position = dodge, alpha=0.1)+
  
  geom_errorbar(data = emm_df,
                aes(x = dpi.f, ymin = lower.CL, ymax = upper.CL, groups = treatment),
                color="black", position = dodge, width = 0.1) +
  geom_line(data = emm_df,
            aes(x = as.numeric(dpi.f), y = emmean, linetype = treatment),
            color="black", position = dodge, size =0.25) +
  
  geom_point(data = emm_df,
             aes(x = dpi.f, y = emmean, color = groups),
             shape=16, position = dodge, size = 2.8) +
  geom_point(data = emm_df,
             aes(x = dpi.f, y = emmean, group=groups),
             position = dodge, size = 2.8,
             shape=1, color="black") +
  
  scale_color_manual(values = treat_colors) +
  scale_shape_manual(values = c(1, 1, 1, 1))+
  labs(x = "Days Post Inoculation", y = "Mass (g)", color = "Treatment Group", shape="Treatment Group", linetype="Inoculation Type") +
  scale_linetype_manual(
    name   = "Inoculation Type",
    values = c("dashed", "solid"),
    labels = c("Control", "Inoculated")
  )+
    facet_wrap(~temp, ncol=2)+
  theme(strip.text = element_text(size=12))


ggplot() +
  geom_point(data = ti.m.g,
             aes(x = dpi.f, y = mass, color = groups),
             position =dodge,
             alpha = 0.25, size = 2) +
  # geom_line(data=ti.m.g,
  #           aes(x=dpi.f, y=mass, color=groups, group = as.factor(band_number)),
  #           position = dodge, alpha=0.1)+
  
  geom_errorbar(data = emm_df,
                aes(x = dpi.f, ymin = lower.CL, ymax = upper.CL, groups = groups),
                color="black", position = dodge, width = 0.1) +
  geom_line(data = emm_df,
            aes(x = as.numeric(dpi.f), y = emmean, linetype = groups),
            color="black", position = dodge, size =0.25) +
  
  geom_point(data = emm_df,
             aes(x = dpi.f, y = emmean, color = groups),
             shape=16, position = dodge, size = 2.8) +
  geom_point(data = emm_df,
             aes(x = dpi.f, y = emmean, group=groups),
             position = dodge, size = 2.8,
             shape=1, color="black") +
  
  scale_color_manual(values = treat_colors) +
  scale_shape_manual(values = c(1, 1, 1, 1))+
  labs(x = "Days Post Inoculation", y = "Mass (g)", color = "Treatment Group", shape="Treatment Group", linetype="Treatment Group") +
  scale_linetype_manual(values= c("dashed", "solid", "dashed", "solid"))+
  # scale_linetype_manual(
  #   name   = "Inoculation Type",
  #   values = c("dashed", "solid"),
  #   labels = c("Control", "Inoculated")
  # )+
  facet_wrap(~treatment, ncol=2)

#Predictions
newdata <- ti.m %>%
  distinct(temp, treatment, dpi.f)

## 2) Predict (population-level: re.form = NA) and add 95% CIs
pp <- predict(lm7, newdata = newdata, se.fit = TRUE, re.form = NA, type = "response")
pred_df <- newdata %>%
  mutate(
    emmean   = as.numeric(pp$fit),
    se.fit   = as.numeric(pp$se.fit),
    lower.CL = emmean - 1.96 * se.fit,
    upper.CL = emmean + 1.96 * se.fit
  )

## 3) (Optional) attach your 'groups' label just for plotting aesthetics
pred_df <- pred_df %>%
  left_join(
    ti.m %>% distinct(temp, treatment, dpi.f, groups),
    by = c("temp", "treatment", "dpi.f")
  )

## 4) Plot: raw data + predicted means/CI from predict()
dodge <- position_dodge(width = 0.7)

ggplot() +
  # raw points
  geom_jitter(
    data = ti.m,
    aes(x = dpi.f, y = mass, color = groups),
    position = position_jitterdodge(dodge.width = 0.7, jitter.width = 0.15, jitter.height = 0),
    alpha = 0.6, size = 1.8
  ) +
  # predicted CIs
  geom_errorbar(
    data = pred_df,
    aes(x = dpi.f, ymin = lower.CL, ymax = upper.CL, group = groups),
    position = dodge, width = 0.12, color="black"
  ) +
  # predicted mean lines (optional)
  geom_line(
    data = pred_df,
    aes(x = dpi.f, y = emmean, linetype = treatment, group = interaction(temp, treatment, groups)),
    position = dodge, linewidth = 0.4, color = "black"
  ) +
  # predicted mean points
  geom_point(
    data = pred_df,
    aes(x = dpi.f, y = emmean, color = groups),
    position = dodge, size = 2.6
  ) +
  # black outline for emphasis
  geom_point(
    data = pred_df,
    aes(x = dpi.f, y = emmean, group = groups),
    position = dodge, size = 2.6, shape = 1, color = "black"
  ) +
  scale_color_manual(values = treat_colors) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  labs(x = "Days Post Inoculation (factor)", y = "Mass (g)",
       color = "Treatment Group", linetype = "Inoculation Type",
       title = "Observed mass with model predictions") +
      facet_wrap(~ temp, ncol = 2) 


##Mass vs fever
ti.mf <- ti.f %>%
  filter(dpi %in% c(3, 14, 28, 35))%>%
  dplyr::select(dpi, dpi.f, band_number, mass, fever_score, total_eye_score, diseased, infected, ever_diseased, ever_infected, groups, temp, treatment, sex, quantity)


ggplot(ti.mf, aes(y=fever_score, x=mass, color=groups))+
  geom_point()+
  scale_color_manual(values = treat_colors)+
  facet_grid(~dpi~temp)

hist(ti.mf$fever_score)

#model selection
a <- glmmTMB(fever_score ~ mass * sex * temp * dpi.f + (1|band_number), data=ti.mf%>% filter(treatment=="Inoculated"))
b <- glmmTMB(fever_score ~ mass + sex * temp * dpi.f + (1|band_number), data=ti.mf%>% filter(treatment=="Inoculated"))
c <- glmmTMB(fever_score ~ mass * sex + temp * dpi.f + (1|band_number), data=ti.mf%>% filter(treatment=="Inoculated"))
d <- glmmTMB(fever_score ~ mass * sex * temp + dpi.f + (1|band_number), data=ti.mf%>% filter(treatment=="Inoculated"))
e <- glmmTMB(fever_score ~ mass + sex + temp * dpi.f + (1|band_number), data=ti.mf%>% filter(treatment=="Inoculated"))
f <- glmmTMB(fever_score ~ mass + sex + temp + dpi.f + (1|band_number), data=ti.mf%>% filter(treatment=="Inoculated"))
g <- glmmTMB(fever_score ~ mass + temp + sex * dpi.f + (1|band_number), data=ti.mf%>% filter(treatment=="Inoculated"))
h <- glmmTMB(fever_score ~ mass + temp * dpi.f + (1|band_number), data=ti.mf%>% filter(treatment=="Inoculated"))

null <- glmmTMB(fever_score ~ 1 + (1|band_number), data=ti.mf%>% filter(treatment=="Inoculated"))

aictab(cand.set=list(a, b, c,  d, e, f, g, h, null), 
       modnames=c("a",  "b", "c", "d",  "e", "f", "g", "h", "null"))

drop1(f)
simulateResiduals(f, plot=T)

#Stepwise Hypotheses: 
#1: Males and females differ in their relationships between body temperature and mass.
  #H1:Males will have higher mass and body temperature.

ti.mf$sex <- as.factor(ti.mf$sex)
glmfm.i <- glmmTMB(fever_score ~ mass + sex + temp + dpi.f + (1|band_number), data=ti.mf %>% filter(treatment=="Inoculated"))
glmfm.c <- glmmTMB(fever_score ~ mass + sex + temp + dpi.f + (1|band_number), data=ti.mf %>% filter(treatment=="Sham"))
glmfm.f <- glmmTMB(fever_score ~ mass + sex + temp + treatment * dpi.f + (1|band_number), data=ti.mf)

#Mass, sex, and temperature affect fever score in inoculated birds
simulateResiduals(glmfm.i, plot=T)
summary(glmfm.i)
Anova(glmfm.i, type="III")

#Only temperature affects fever score in control birds
simulateResiduals(glmfm.c, plot=T)
summary(glmfm.c)
Anova(glmfm.c, type="III")

#Sex and temperature affect fever score when looking at all birds.
simulateResiduals(glmfm.f, plot=T)
summary(glmfm.f)
Anova(glmfm.f, type="III")

#So mass only appears to play a role in predicting fever score during infection.

ggplot(ti.mf %>% filter(dpi == 14), aes(x=fever_score, y=diseased, color=sex))+
  geom_point()+
  facet_wrap(~temp)

#What about infection versus disease (asymptomatic vs symptomatic infections)
  #These seem like they should be slightly different as the fever score is a measure of temperature at the sight of infection
  #Symptomatic infections likely will have higher fever scores


#First, cutoffs
ti.mf$quant_cutoff = 50
ti.mf$sympt_cutoff = 0.1

#diseased = eye score
#infected = eye score or pathology >= cutoff
ti.mf <- ti.mf %>%
  mutate(diseased= ifelse(total_eye_score>=sympt_cutoff, 1, 0),
         sick = ifelse(total_eye_score>=sympt_cutoff | quantity >= quant_cutoff, 1, 0),
         infected = ifelse(quantity >= quant_cutoff, 1, 0))

#new column with whether the bird is ever diseased
ti.mf <- ti.mf %>%
  group_by(band_number) %>%
  mutate(ever_diseased = ifelse(any(coalesce(diseased, 0) == 1), 1, 0),
         ever_infected = ifelse(any(coalesce(infected, 0) == 1), 1, 0),
         ever_sick = ifelse(any(coalesce(sick, 0) == 1), 1,0))%>%
  ungroup()


#So first, let's look at infections by sex
lmsi <- glm(as.factor(ever_infected) ~ sex, data=ti.mf %>% filter(treatment == "Inoculated"), family=binomial())

simulateResiduals(lmsi, plot=T)
summary(lmsi)

ggplot(ti.mf %>% filter(treatment == "Inoculated" & dpi == 14), aes(x=sex, y=fever_score, color=as.factor(ever_infected)))+
  geom_jitter(width=0.1)+
  facet_wrap(~ever_diseased)

#And now disease
lmsd <- glm(ever_diseased ~ sex, data=ti.mg %>% filter(treatment == "Inoculated"), family=binomial())

simulateResiduals(lmsd, plot=T)
summary(lmsd)
#huge effect of sex
plot(allEffects(lmsd))

#So there is no difference in infection by pathogen load, but males were much more likely to develop disease

glm.d <- glmmTMB(ever_diseased ~ sex * temp, data=ti.mf %>% filter(dpi == 14 & treatment == "Inoculated"), family=binomial())
glm.i <- glmmTMB(ever_infected ~ sex * temp, data=ti.mf %>% filter(dpi == 14 & treatment == "Inoculated"), family=binomial())

simulateResiduals(glm.d, plot=T)
summary(glm.d)
plot(allEffects(glm.d))
hist(resid(glm.d))

simulateResiduals(glm.i, plot=T)
summary(glm.i)
plot(allEffects(glm.i))
hist(resid(glm.i))

ti.mf %>%
  filter(dpi == 14)%>%
  dplyr::select(temp, treatment, groups, sex, ever_infected, ever_diseased)%>%
  tbl_summary(
    by=groups
  )%>%
  modify_header(
    label ~ "**All Birds**"
  )

ti.mf %>%
  filter(dpi == 14) %>%
  dplyr::select(sex, groups, ever_infected, ever_diseased) %>%
  tbl_strata(
    strata = sex,
    .tbl_fun = ~ tbl_summary(
      data = .x,
      by = groups,
      statistic = everything() ~ "{n} ({p}%)"
    )
  ) %>%
  modify_header(label ~ "**All Birds**")

library(scales)  # for percent_format

# Make a long dataset with both outcomes
ti_long <- ti.mf %>%
  dplyr::select(sex, temp, treatment, ever_infected) %>%
  pivot_longer(
    cols = c(ever_infected),
             #,ever_diseased),
    names_to = "outcome",
    values_to = "status"
  ) %>%
  dplyr::mutate(
    status = ifelse(status == 1, "Yes", "No"),
    outcome = dplyr::recode(outcome,
                     "ever_infected" = "Ever infected (qPCR+)")
                     #"ever_diseased" = "Ever diseased (eye score)")
  )

# Summarize to get proportions
summary_df <- ti_long %>%
  count(sex, temp, treatment, outcome, status) %>%
  group_by(sex, temp, treatment, outcome) %>%
  mutate(
    total = sum(n),
    prop = n / total,
    label = paste0(n, "/", total)        # text label
  ) %>%
  ungroup()

# Plot: stacked bars, proportions by sex
ggplot(summary_df %>% filter(treatment == "Inoculated"), aes(x = sex, y = prop, fill = status)) +
  geom_col() +
  geom_text(data = summary_df %>%
              filter(treatment == "Inoculated", status == "Yes"),
            aes(label = label, y = 1.5), size = 4)+
  facet_grid(~ temp ~ outcome) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_fill_manual(values = c("gray50", "black"))+
  labs(x = "Sex",
       y = "Proportion of Birds",
       fill = "Status") +
  theme(strip.text = element_text(size=12))


ti_sub <- ti.mf %>%
  ungroup() %>%
  mutate(infected_no_path = case_when(
    ever_infected == 1 & ever_diseased == 0 ~ "Yes",
    TRUE ~ "No"
  ))

summary_no_path <- ti_sub %>%
  count(sex, temp, infected_no_path) %>%
  group_by(sex, temp) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

ggplot(summary_no_path, aes(x = sex, y = prop, fill = infected_no_path)) +
  geom_col() +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_fill_manual(values = c("gray50", "black"))+
  labs(x = "Sex",
       y = "Proportion of birds",
       fill = "Infected but no pathology")+
  facet_wrap(~temp)
