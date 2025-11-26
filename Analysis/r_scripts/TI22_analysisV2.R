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
    label ~ "**Birds for Antibodies**"
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
lm1<-glmmTMB(elisa_od~temp+dpi + treatment*dpi + (1|band_number),
                 data=ti.ab, 
                 family=Gamma(link = "log"))

simulateResiduals(lm1, plot=T)
summary(lm1)

lm1s<-glmmTMB(elisa_od~temp+dpi + treatment*dpi + sex + (1|band_number),
             data=ti.ab, 
             family=Gamma(link = "log"))

#sex not significant
anova(lm1, lm1s, type="Chisq")
summary(lm1s)

Anova(lm1, type = "II")

emm <- emmeans(lm1, ~ treatment | temp*dpi, type = "response")

pairs(emm, by=c("temp","dpi"), adjust = "tukey")

emm_df <- as.data.frame(emm)

emm_df <- emm_df %>%
  mutate(
    groups = case_when(
      temp == "Warm" & treatment == "Control"    ~ "Warm Control",
      temp == "Warm" & treatment == "Inoculated" ~ "Warm Inoculated",
      temp == "Cold" & treatment == "Control"    ~ "Cold Control",
      temp == "Cold" & treatment == "Inoculated" ~ "Cold Inoculated"
    ),
    groups = factor(groups,
                    levels = c("Warm Control", "Warm Inoculated",
                               "Cold Control", "Cold Inoculated"))
  )

dodge = position_dodge(0.8)
ggplot(emm_df, aes(x = groups, y = response, color = groups, shape=dpi)) +
  geom_jitter(
    data = ti.ab,
    aes(y = elisa_od, color = groups, shape = dpi),
    size = 2.5, alpha=1,
    position = position_jitterdodge(dodge.width = 0.8, jitter.width=0.5)
  ) +
  geom_point(size = 3, color="black", stroke = 1,
             position = dodge) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.0, color="black",
                position = dodge) +
  scale_color_manual(values = treat_colors, name = "Treatment") +
  scale_shape_manual(values = c(1, 17, 16))+
  #facet_wrap(~ ever_diseased, nrow = 1) +
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
    position = position_jitterdodge(dodge.width = 0.75, jitter.width=0.15)
  ) +
  geom_point(size = 3, color="black",
             position = dodge) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.0, color="black",
                position = dodge) +
  scale_color_manual(values = treat_colors, name = "Treatment") +
  scale_shape_manual(values = c(0, 1, 16))+
  #facet_wrap(~ dpi, nrow = 2) +
  labs(x = "Treatment Group", y = "ELISA OD", shape = "Days Post Inoculation") +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggplot(ti.ab, aes(x=dpi, y=elisa_od, color=groups))+
  geom_point(aes(shape=sex), size=2)+
  geom_line(aes(group = as.factor(band_number)), linewidth=0.1)+
  scale_color_manual(values=treat_colors)+
  facet_grid(
    ever_infected ~ ever_diseased,
    labeller = labeller(
      ever_infected = c(`0` = "Never Infected", `1` = "Ever Infected"),
      ever_diseased = c(`0` = "Never Diseased", `1` = "Ever Diseased")
    )
  )+
  theme(strip.text = element_text(size=12))+
  labs(x="Days Post Inoculation", y="ELISA OD", color="Treatment Groups", shape="Sex")

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

dat.new <- dat.new %>%
  mutate(
    groups = case_when(
      temp == "Warm" & treatment == "Control"    ~ "Warm Control",
      temp == "Warm" & treatment == "Inoculated" ~ "Warm Inoculated",
      temp == "Cold" & treatment == "Control"    ~ "Cold Control",
      temp == "Cold" & treatment == "Inoculated" ~ "Cold Inoculated"
    ),
    groups = factor(groups,
                    levels = c("Warm Control", "Warm Inoculated",
                               "Cold Control", "Cold Inoculated"))
  )


dodge = position_dodge(0.8)
ggplot(dat.new, aes(x = groups, y = yhat, color = groups, shape=dpi)) +
  geom_jitter(
    data = ti.ab,
    aes(y = elisa_od, color = groups, shape = dpi),
    size = 3, alpha=0.75,
    position = position_jitterdodge(dodge.width = 0.8, jitter.width=0.6)
  ) +
  geom_point(size = 3, color="black", stroke = 1,
             position = dodge) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.0, color="black",
                position = dodge) +
  scale_color_manual(values = treat_colors, name = "Treatment") +
  scale_shape_manual(values = c(1, 17, 16))+
  #facet_wrap(~ ever_diseased, nrow = 1) +
  labs(x = "Treatment Group", y = "ELISA OD", shape = "Days Post Inoculation") +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
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
#write_xlsx(out_list, "/Users/jesse/Documents/GitHub/Temperature-and-Immunity-2022/Analysis/Results/Antibody_Results.xlsx")

####Sex Effects on Disease####
source("r_scripts/dataCleaning_TI22.R")

ti %>%
  filter(dpi == 0)%>%
  dplyr::select(groups, ever_infected, ever_diseased)%>%
  tbl_summary(
    by=groups
  )%>%
  modify_header(
    label ~ "**All Birds**"
  )

ti.inoc <- ti %>%
  filter(dpi == 0 & treatment == "Inoculated")

ti.inoc$ever_diseased <- as.factor(ti.inoc$ever_diseased)
ti.inoc$ever_infected <- as.factor(ti.inoc$ever_infected)


da <- glm(ever_diseased ~ sex + temp, data=ti.inoc, family=binomial())
db <- glm(ever_diseased ~ sex * temp, data=ti.inoc, family=binomial())
dc <- glm(ever_diseased ~ sex, data=ti.inoc, family=binomial())
dd <- glm(ever_diseased ~ temp, data=ti.inoc, family=binomial())
dn <- glm(ever_diseased ~ 1, data=ti.inoc, family=binomial())

aictab(cand.set=list(da, db, dc, dd, dn), 
       modnames=c("da", "db", "dc", "dd", "null"))

lm.dis <- glm(ever_diseased ~ sex + temp, data=ti.inoc, family=binomial())
simulateResiduals(lm.dis, plot=T)
summary(lm.dis)
car::Anova(lm.dis, type = "II")

#glm is not better than null model: so try Fisher's Exact Test for each temperature
# ftd.m <- data.frame(
#   "Cold" = c(0, 4, 3),
#   "Warm" = c(0, 1, 6),
#   row.names = c("Uninfected", "Asymptomatic", "Diseased"),
#   stringsAsFactors = FALSE)
# colnames(ftd.m)
# 
# mosaicplot(ftd.m,
#            main = "Mosaic",
#            color = TRUE)
# 
# chisq.test(ftd.m)$expected
# 
# #Is there a difference in the breakdown between statuses in males between temperatures?
# ftm <- fisher.test(ftd.m)
# ftm
# 
# ftd.f <- data.frame(
#   "Cold" = c(2, 2, 3),
#   "Warm" = c(1, 4, 2),
#   row.names = c("Uninfected", "Asymptomatic", "Diseased"),
#   stringsAsFactors = FALSE)
# colnames(ftd.f)
# 
# mosaicplot(ftd.f,
#            main = "Mosaic",
#            color = TRUE)
# 
# chisq.test(ftd.f)$expected
# 
# #Is there a difference in the breakdown between statuses in females between temperatures?
# ftf <- fisher.test(ftd.f)
# ftf
# 
# #Within temperatures
# ftd.c <- data.frame(
#   "Female" = c(2, 2, 3),
#   "Male" = c(0, 4, 3),
#   row.names = c("Uninfected", "Asymptomatic", "Diseased"),
#   stringsAsFactors = FALSE)
# colnames(ftd.c)
# 
# mosaicplot(ftd.c,
#            main = "Mosaic",
#            color = TRUE)
# 
# chisq.test(ftd.c)$expected
# 
# #Is there a difference in the breakdown between statuses between sexes in the cold room?
# ftc <- fisher.test(ftd.c)
# ftc
# 
# #Within temperatures - warm
# ftd.w <- data.frame(
#   "Female" = c(1, 4, 2),
#   "Male" = c(0, 6, 1),
#   row.names = c("Uninfected", "Asymptomatic", "Diseased"),
#   stringsAsFactors = FALSE)
# colnames(ftd.w)
# 
# mosaicplot(ftd.w,
#            main = "Mosaic",
#            color = TRUE)
# 
# chisq.test(ftd.w)$expected
# 
# #Is there a difference in the breakdown between statuses between sexes in the cold room?
# ftw <- fisher.test(ftd.w)
# ftw

a<-glm(ever_infected ~ sex + temp, data=ti.inoc, family=binomial())
null <- glm(ever_infected ~ 1, data=ti.inoc, family=binomial())

aictab(cand.set=list(a, null), 
       modnames=c("a", "null"))

lm.inf <- glm(ever_infected ~ sex + temp, data=ti.inoc, family=binomial())
simulateResiduals(lm.inf, plot=T)
summary(lm.inf)
car::Anova(lm.inf, type = "III")

ti_inc <- ti %>%
  filter(treatment == "Inoculated",
         dpi == 0) %>%
  mutate(
    status = dplyr::case_when(
      ever_infected == 1 & ever_diseased == 1 ~ "Diseased",
      ever_infected == 1 & ever_diseased == 0 ~ "Asymptomatic",
      ever_infected == 0 & ever_diseased == 1 ~ "Diseased Only",
      ever_infected == 0 & ever_diseased == 0 ~ "Uninfected"
    ),
    status = factor(
      status,
      levels = c(
        "Uninfected",
        "Diseased Only",
        "Asymptomatic",
        "Diseased"
      )
    )
  )

sum_status <- ti_inc %>%
  group_by(sex, status) %>%
  summarise(n = n(), .groups = "drop")

ggplot(sum_status,
       aes(x = sex, y = n, fill = status)) +
  geom_col(position = "stack", color="black") +
  #facet_wrap(~ temp) +
  scale_fill_manual(values=c("grey70", "grey40","black"))+
  geom_text(aes(label = n),
            position = position_stack(vjust=0.5),
            color="grey90", size=5, fontface="bold")+
  labs(
    x = "Sex",
    y = "Number of individuals",
    fill = "Status",
    title = "Disease Status by Sex (Inoculated Birds)"
  )+
  theme(strip.text = element_text(size=12))

sum_status_t <- ti_inc %>%
  group_by(temp, sex, status) %>%
  summarise(n = n(), .groups = "drop")

ggplot(sum_status_t,
       aes(x = sex, y = n, fill = status)) +
  geom_col(position = "stack", color="black", size=0.75) +
  facet_wrap(~ temp) +
  scale_fill_manual(values=c("grey70", "grey40","black"))+
  geom_text(aes(label = n),
            position = position_stack(vjust=0.5),
            color="grey90", size=5, fontface="bold")+
  labs(
    x = "Sex",
    y = "Number of individuals",
    fill = "Status",
    title = "Disease Status by Sex and Temperature (Inoculated Birds)"
  )+
  theme(strip.text = element_text(size=12))

#do temperature or sex predict whether a bird develops an infection when inoculated?
res <- glm(ever_infected ~ temp * sex, data=ti.inoc, family=binomial())
simulateResiduals(res, plot=T)
summary(res)
car::Anova(res, type="III")

#do temperature or sex predict whether a bird develops pathology when inoculated?
dis <- glm(ever_diseased ~ temp * sex, data=ti.inoc, family=binomial())
simulateResiduals(dis, plot=T)
summary(dis)
car::Anova(dis, type="III")

#do temperature or sex predict whether a bird develops an eye score when infected?
dis.i <- glm(ever_diseased ~ temp * sex, data=ti.inoc %>% filter(ever_infected == 1), family=binomial())
simulateResiduals(dis.i, plot=T)
summary(dis.i)

####Fever Score####
source("r_scripts/dataCleaning_fever.R")
#fever_change = score - baseline
#fever_diff = change from previous score

ti.f$dpi <- as.factor(ti.f$dpi)

ggplot(ti.f, aes(x=dpi, y=fever_score, color=groups))+
  geom_line(aes(group=as.factor(band_number)))+
  scale_color_manual(values=treat_colors)+
  facet_wrap(~temp~ever_diseased)

#peak score
ti.f <- ti.f %>%
  group_by(band_number)%>%
  mutate(fever_peak = max(fever_score),
         fever_high = max(fever_high = max(fever_score[dpi != 0])))

ti.f$dpi.f <- as.factor(ti.f$dpi)


#####Fever 1) Do kinetics of fever response differ between temperature, treatment, or sex?####

##Fever by dpi##
ggplot(ti.f, aes(x=dpi, y=fever_score, color=groups))+
  geom_point()

ti.f$dpi.f <- as.factor(ti.f$dpi)
ti.f$temp <- relevel(as.factor(ti.f$temp), ref = "Warm")
ti.f$treatment <- relevel(as.factor(ti.f$treatment), ref = "Control")

#did fever differ between treatments at baseline?
t.test(fever_score ~ treatment, data=ti.f %>% filter(dpi.f == 0))

#did temperature differ between sexes between the rooms at baseline?
#model selection: Does sex predict body temperature when controlling for room temperature?
ti.f.cont <- ti.f %>%
  filter(treatment == "Control")
c1 <- glmmTMB(fever_score ~ temp + sex + dpi.f + (1|band_number), data=ti.f.cont)
c2 <- glmmTMB(fever_score ~ temp + sex * dpi.f + (1|band_number), data=ti.f.cont)
c3 <- glmmTMB(fever_score ~ temp * sex + dpi.f + (1|band_number), data=ti.f.cont)
c4 <- glmmTMB(fever_score ~ temp * sex * dpi.f + (1|band_number), data=ti.f.cont)
c5 <- glmmTMB(fever_score ~ temp + dpi.f + (1|band_number), data=ti.f.cont)
c6 <- glmmTMB(fever_score ~ sex + dpi.f + (1|band_number), data=ti.f.cont)
null_fc <- glmmTMB(fever_score ~ 1 + (1|band_number), data=ti.f.cont)

aictab(list(c1, c2, c3, c4, c5, c6, null_fc),
       modnames = c("c1","c2","c3","c4","c5","c6","null"))

#So it isn't actually a sex effect
summary(c1)

#It's sex and infection
lm.ts <- glmmTMB(fever_score ~ temp + treatment + sex + dpi.f + (1|band_number), data=ti.f)
simulateResiduals(lm.ts, plot=T)
summary(lm.ts)
plot(allEffects(lm.ts))
car::Anova(lm.ts, type="III")

#Model Selection
#because sampled_first is a function of dpi (it is fixed for each dpi), it is collinear with dpi.f, so do not include.
#We treat dpi as a factor (dpi.f) because this is a linear model therefore we must look only at linear differences
#between dpis (no smoothing)

#I want to know whether treatment or temperature affect fever_score or shape of fever response
#while accounting for sex
m0 <- glmmTMB(fever_score ~ temp + treatment + sex + dpi.f + (1|band_number), data = ti.f)
m1 <- glmmTMB(fever_score ~ temp * treatment + sex + dpi.f + (1|band_number), data = ti.f)
m2 <- glmmTMB(fever_score ~ temp * dpi.f + treatment + sex + (1|band_number), data = ti.f)
m3 <- glmmTMB(fever_score ~ treatment * dpi.f + temp + sex + (1|band_number), data = ti.f)
m4 <- glmmTMB(fever_score ~ temp * treatment * dpi.f + sex + (1|band_number), data = ti.f)
m5 <- glmmTMB(fever_score ~ temp * treatment * dpi.f * sex + (1|band_number), data = ti.f)
m6 <- glmmTMB(fever_score ~ (treatment + temp) * dpi.f + sex + (1|band_number), data = ti.f)
m6s <- glmmTMB(fever_score ~ (treatment + temp) * dpi.f * sex + (1|band_number), data = ti.f)
m7 <- glmmTMB(fever_score ~ (treatment + temp + sex) * dpi.f + (1|band_number), data = ti.f)
m_null <- glmmTMB(fever_score ~ 1 + (1|band_number), data = ti.f)

aictab(list(m0, m1, m2, m3, m4, m5, m6, m6s, m7, m_null),
       modnames = c("m0","m1","m2","m3","m4","m5","m6","m6s", "m7", "null"))

simulateResiduals(m6, plot=T)
hist(resid(m6))
summary(m6)
car::Anova(m6, type="III")

#Interaction between temp and dpi does not have one significant dpi, but overall interaction is significant
m_full  <- glmmTMB(fever_score ~ (treatment + temp) * dpi.f + sex + (1|band_number), data=ti.f)
m_drop  <- glmmTMB(fever_score ~ treatment * dpi.f + temp + dpi.f + sex + (1|band_number), data=ti.f)

#m_full has lower AIC = do not drop the interaction
AIC(m_full, m_drop)

glm_fever <- glmmTMB(fever_score ~ (treatment + temp) * dpi.f + sex + (1|band_number), data=ti.f)

simulateResiduals(glm_fever, plot=T)

#Residuals deviate a little - how is residual distribution?
hist(resid(glm_fever))
#Normally distributed more or less

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
  
  geom_errorbar(data = emm_df,
                aes(x = dpi.f, ymin = lower.CL, ymax = upper.CL, groups = groups),
                color="black", position = dodge, width = 0.1) +
  
  # geom_ribbon(data = emm_df,
  #               aes(x = as.numeric(dpi.f), ymin = lower.CL, ymax = upper.CL, groups = groups, fill=groups),
  #               position = dodge, alpha=0.2) +
  # 
  geom_line(data = emm_df,
            aes(x = as.numeric(dpi.f), y = emmean, linetype = treatment, groups = groups),
            color="black", position = dodge, size =0.25) +
  
  geom_point(data = emm_df,
             aes(x = dpi.f, y = emmean, color = groups),
             position = dodge, size = 2.8) +
  geom_point(data = emm_df,
             aes(x = dpi.f, y = emmean, group=groups),
             position = dodge, size = 2.8,
             shape=1, color="black") +
  # stat_summary(data = ti.f.g, aes(x=dpi.f, y= fever_score, group = interaction(treatment, temp), color=groups), geom="line", fun = "mean")+
  
  scale_color_manual(values = treat_colors) +
  #scale_fill_manual(values=treat_colors)+
  scale_shape_manual(values = c(1, 1, 1, 1))+
  labs(x = "Days Post Inoculation", y = "Ocular Temperature (C)", color = "Treatment Group", shape="Treatment Group", linetype="Inoculation Type") +
  scale_linetype_manual(
    name   = "Inoculation Type",
    values = c("dashed", "solid"),
    labels = c("Control", "Inoculated")
  )+
  facet_grid(~sex)+ #, scales = "free_y"
  theme(strip.text = element_text(size=12))


#Results:

#                       Chisq Df Pr(>Chisq)    
# (Intercept)     13305.2657  1  < 2.2e-16 ***
# treatment           0.0085  1  0.9263886    
# temp              140.9492  1  < 2.2e-16 ***
# dpi.f              27.5258  6  0.0001154 ***
# sex                 7.0116  1  0.0080985 ** 
# treatment:dpi.f    41.0611  6  2.817e-07 ***
# temp:dpi.f         20.1471  6  0.0026072 ** 

#Cold temp = lower eye temperature
# temps differed across dpi
# Temps differed between sexes (M higher)
# Different kinetics between treatments (inoculated increased ocular temperature)
# Different kinetics between temps


#Model predictions by temp comparing sex directly in inoculated individuals
#We know there is a trend for males to be more likely to develop pathology when infected in warm temperatures
#Therefore, we would expect higher ocular temps at high temperatures, so we need a sex interaction
#But subset to only inocualted individuals

ti.f.i <- ti.f %>% 
  filter(treatment == "Inoculated")

#Model selection
m0i <- glmmTMB(fever_score ~ temp + sex + dpi.f + (1|band_number), data = ti.f.i)
m1i <- glmmTMB(fever_score ~ temp * sex + dpi.f + (1|band_number), data = ti.f.i)
m2i <- glmmTMB(fever_score ~ sex * dpi.f + temp + (1|band_number), data = ti.f.i)
m3i <- glmmTMB(fever_score ~ temp * dpi.f + sex + (1|band_number), data = ti.f.i)
m4i <- glmmTMB(fever_score ~ dpi.f * temp * sex + (1|band_number), data = ti.f.i)
m_nulli <- glmmTMB(fever_score ~ 1 + (1|band_number), data = ti.f.i)

aictab(list(m0i, m1i, m2i, m3i, m4i, m_nulli),
       modnames = c("m0","m1","m2","m3","m4", "null"))

simulateResiduals(m1i, plot=T)

glm_fever.i <- glmmTMB(fever_score ~ temp * sex + dpi.f + (1|band_number), data = ti.f.i)
simulateResiduals(glm_fever.i, plot=T)
summary(glm_fever.i)
car::Anova(glm_fever.i, type = "III")

# Get emmeans (predicted marginal means)
emm_df <-  emmeans(glm_fever.i, ~ sex * dpi.f * temp, type = "response") %>%
  as.data.frame() %>%
  left_join(
    ti.f.i %>%
      distinct(temp, dpi.f, sex, groups),
    by = c("temp", "dpi.f", "sex")
  ) 

emm_df <- emm_df %>%
  mutate(sex = dplyr::recode(sex,
                             "M" = "Male",
                             "F" = "Female"))

ti.f.i <- ti.f.i %>%
  dplyr::mutate(sex = dplyr::recode(sex,
                                    "M" = "Male",
                                    "F" = "Female"))

ggplot() +
  geom_point(data = ti.f.i,
             aes(x = dpi.f, y = fever_score, color = sex, shape = sex),
             position =dodge,
             alpha = 0.5, size = 2) +
  # geom_line(data=ti.m.g,
  #           aes(x=dpi.f, y=mass, color=groups, group = as.factor(band_number)),
  #           position = dodge, alpha=0.1)+
  
  geom_errorbar(data = emm_df,
                aes(x = dpi.f, ymin = lower.CL, ymax = upper.CL, groups = sex),
                color="black", position = dodge, width = 0.1) +
  geom_line(data = emm_df,
            aes(x = as.numeric(dpi.f), y = emmean, linetype = sex),
            color="black", position = dodge, size =0.25) +
  
  geom_point(data = emm_df,
             aes(x = dpi.f, y = emmean, color = sex, groups = sex, shape = sex),
             position = dodge, size = 2.8) +
  # geom_point(data = emm_df,
  #            aes(x = dpi.f, y = emmean, group=groups, shape ),
  #            position = dodge, size = 2.8,
  #            color="black") +
  
  scale_color_manual(values = c("black", "gray50")) +
  scale_shape_manual(values = c(16, 17))+
  labs(x = "Days Post Inoculation", y = "Ocular Temperature (C)", color = "Treatment Group", shape="Treatment Group", linetype="Inoculation Type") +
  # scale_linetype_manual(
  #   name   = "Inoculation Type",
  #   values = c("dashed", "solid"),
  #   labels = c("Control", "Inoculated")
  # )+
  facet_grid(~temp)+ #, scales = "free_y"
  theme(strip.text = element_text(size=12))


#Males have higher ocular temperature in warm rooms - is this because more develop pathology?
#Males more likely to be diseased (almost), if we account for this, does effect go away?

#Model selection
m1id <- glmmTMB(fever_score ~ temp * sex + dpi.f + ever_diseased + (1|band_number), data = ti.f.i)
m2id <- glmmTMB(fever_score ~ sex * ever_diseased * dpi.f + temp + (1|band_number), data = ti.f.i)
m3id <- glmmTMB(fever_score ~ sex * ever_diseased + dpi.f + temp + (1|band_number), data = ti.f.i)
m4id <- glmmTMB(fever_score ~ sex + ever_diseased * dpi.f + temp + (1|band_number), data = ti.f.i)
m4ids <- glmmTMB(fever_score ~ sex * ever_diseased * dpi.f + temp + (1|band_number), data = ti.f.i)
m5id <- glmmTMB(fever_score ~ sex + ever_diseased * dpi.f * temp + (1|band_number), data = ti.f.i)
m6id <- glmmTMB(fever_score ~ dpi.f * temp * sex * ever_diseased + (1|band_number), data = ti.f.i)
m7id <- glmmTMB(fever_score ~ ever_diseased * dpi.f + temp * sex + (1|band_number), data = ti.f.i)
m8id <- glmmTMB(fever_score ~ (sex + ever_diseased) * dpi.f + temp + (1|band_number), data = ti.f.i)

m_nulli <- glmmTMB(fever_score ~ 1 + (1|band_number), data = ti.f.i)

aictab(list(m1id, m2id, m3id, m4id, m4ids, m5id, m6id, m7id, m8id, m_nulli),
       modnames = c("m1","m2","m3","m4", "m4ids", "m5", "m6", "m7", "m8", "null"))

simulateResiduals(m4id, plot=T)
summary(m4id)
car::Anova(m4id, type="III")

glm_fever.id <- glmmTMB(fever_score ~ sex + ever_diseased * dpi.f + temp + (1|band_number), data = ti.f.i)
simulateResiduals(glm_fever.id, plot=T)
summary(glm_fever.id)
car::Anova(glm_fever.id, type = "III")

#compare to treatment
glm_fd <- glmmTMB(fever_score ~ sex + ever_diseased * dpi.f + temp + (1|band_number), data = ti.f)
glm_ft <- glmmTMB(fever_score ~ sex + treatment * dpi.f + temp + (1|band_number), data = ti.f)
glm_fi <- glmmTMB(fever_score ~ sex + ever_infected * dpi.f + temp + (1|band_number), data = ti.f)


aictab(list(glm_fd, glm_ft, glm_fi),
       modnames = c("glm_fd","glm_ft","glm_fi"))

# Get emmeans (predicted marginal means)
emm_df <-  emmeans(glm_fever.id, ~ sex * dpi.f * temp * ever_diseased, type = "response") %>%
  as.data.frame() 

emm_df <- emm_df %>%
  mutate(sex = dplyr::recode(sex,
                             "M" = "Male",
                             "F" = "Female"))



ti.f.g <- ti.f.i %>%
  dplyr::mutate(sex = dplyr::recode(sex,
                                    "M" = "Male",
                                    "F" = "Female"))

dodge <- position_dodge(width = 0.5)

#Model predictions by temp
ggplot() +
  geom_point(data = ti.f.g,
             aes(x = dpi.f, y = fever_score, color = temp, shape=sex, fill=temp),
             position =dodge,
             alpha = 0.75, size = 2) +
  # geom_line(data=ti.m.g,
  #           aes(x=dpi.f, y=mass, color=groups, group = as.factor(band_number)),
  #           position = dodge, alpha=0.1)+
  
  geom_errorbar(data = emm_df,
                aes(x = dpi.f, ymin = lower.CL, ymax = upper.CL, group = interaction(sex, temp)),
                color="black", position = dodge, width = 0.1) +
  geom_line(data = emm_df,
            aes(x = as.numeric(dpi.f), y = emmean, linetype = sex, group = interaction(sex, temp)),
            position = dodge) +
  
  geom_point(data = emm_df,
             aes(x = dpi.f, y = emmean, fill = temp, shape=sex),
             position = dodge, size = 2.8, stroke = 1) +
  # stat_summary(data=ti.f.g, aes(x=dpi.f, y=fever_score, color = sex, groups= interaction(ever_diseased, sex, temp)),
  #              geom="point", fun = mean)+
  
  scale_color_manual(values = temp_colors) +
  scale_fill_manual(values = temp_colors)+
  scale_shape_manual(values = c(21, 24))+
  labs(x = "Days Post Inoculation", y = "Ocular Temperature (C)", color = "Temperature", shape="Sex", linetype="Sex") +
  guides(fill=FALSE)+
  scale_linetype_manual(
    values = c("dashed", "solid")
  )+
  facet_grid(~ever_diseased, 
             labeller=labeller(
               ever_diseased = c('0' = "Never Diseased", '1' = "Ever Diseased")))+
  theme(strip.text = element_text(size=12))



####Fever 2) Baseline versus peak fever score####
#Look only at baseline versus peak fever score. This simplifies models a ton and still looks at magnitude of response
#Lood at differences between baseline (DPI 0) and peak fever score regardless of when it occurs.
  #one observation per bird
ti.fc <- ti.f %>%
  dplyr::select(band_number, treatment, temp, groups, dpi, fever_score, fever_peak, mass, sex, ever_infected, ever_diseased) %>%
  group_by(band_number, treatment, temp, groups, sex, ever_infected, ever_diseased) %>%
  summarise(
    baseline   = fever_score[dpi == 0],
    peak       = max(fever_score, na.rm = TRUE),
    high       = max(fever_high = max(fever_score[dpi != 0])),
    magnitude  = peak - baseline,
    mag_high   = high - baseline,
    .groups    = "drop"
  ) 

#What if we omit DPI 18?
ti.fc.omit <- ti.f %>%
  filter(dpi != 18)%>%
  dplyr::select(band_number, treatment, temp, groups, dpi, fever_score, fever_peak, mass, sex, ever_infected, ever_diseased) %>%
  group_by(band_number, treatment, temp, groups, sex, ever_infected, ever_diseased) %>%
  summarise(
    baseline   = fever_score[dpi == 0],
    peak       = max(fever_score, na.rm = TRUE),
    high       = max(fever_high = max(fever_score[dpi != 0])),
    magnitude  = peak - baseline,
    mag_high   = high - baseline,
    .groups    = "drop"
  ) 

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
  facet_wrap(~ever_diseased)


ti.long <- pivot_longer(ti.fc, #To check robustness of analysis, change to ti.fc.omit and run analysis to exclude DPI 18
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
ti.long$treatment <- relevel(as.factor(ti.long$treatment), ref = "Control")
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
fc6 <- glmmTMB(fever_value ~ fever_type * treatment + fever_type * sex + temp * sex + (1|band_number), data=ti.long)

fc7 <-glmmTMB(fever_value ~ 1 + (1|band_number), data=ti.long)

aictab(cand.set=list(fc1, fc1s, fc2, fc2s, fc3, fc3s, fc4, fc4s, fc5, fc5s, fc5si, fc5si2, fc5si2i, fc6, fc7), 
       modnames=c("fc1",  "fc1s", "fc2",  "fc2s", "fc3", "fc3s", "fc4", "fc4s", "fc5", "fc5s",  "fc5si", "fc5si2", "fc5si2i", "fc6", "fc7"))

#Does the interaction between fever type and treatment or temperature predict fever_value?
  #Is fever_value affected by the interaction between fever type and treatment, the interaction between temperature and sex,
  #or each main effect while accounting for repeated measures
lm2<- glmmTMB(fever_value ~ fever_type * treatment + temp * sex + (1|band_number), data=ti.long)

simulateResiduals(lm2, plot=T)
summary(lm2)
plot(allEffects(lm2))

####Omitting DPI 18 Results:

#                                     Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                        36.120018   0.344555  104.83  < 2e-16 ***
# fever_typePeak                     -0.143750   0.365794   -0.39 0.694334    
# treatmentInoculated                 0.002291   0.350067    0.01 0.994779    
# tempCold                           -3.030851   0.363427   -8.34  < 2e-16 ***
# sexMale                             1.255208   0.355191    3.53 0.000409 ***
# fever_typePeak:treatmentInoculated  1.358929   0.478936    2.84 0.004548 ** 
# tempCold:sexMale                   -1.291523   0.503814   -2.56 0.010363 *


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

#sex side by side
ti.long <- ti.long %>%
  mutate(
    sex_phase = interaction(sex, fever_type),
    sex_phase = recode_factor(
      sex_phase,
      "Female.Baseline" = "Female Baseline",
      "Male.Baseline" = "Male Baseline",
      "Female.Peak"     = "Female Peak",
      "Male.Peak"     = "Male Peak"
    )
  )

emm_df <- emm_df %>%
  mutate(
    sex_phase = interaction(sex, fever_type),
    sex_phase = recode_factor(
      sex_phase,
      "Female.Baseline" = "Female Baseline",
      "Male.Baseline" = "Male Baseline",
      "Female.Peak"     = "Female Peak",
      "Male.Peak"     = "Male Peak"
    )
  )

#Plot Emmeans by sex
ggplot(ti.long, aes(x = groups,
                    y = fever_value,
                    color = groups,
                    shape = sex_phase)) +
  
  geom_jitter(
    position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.75),
    alpha = 0.75, size = 2.5, stroke=1
  ) +
  # model means — DODGE BY fever_type
  geom_point(
    data = emm_df,
    aes(x = groups, y = emmean, group = sex_phase, shape = sex_phase, color=groups),
    position = position_dodge(width = 0.75),
    size =3, 
    stroke = 1,
    color = "black",
    inherit.aes = FALSE
  ) +
  geom_errorbar(
    data = emm_df,
    aes(x = groups, ymin = lower.CL, ymax = upper.CL,
        group = sex_phase, color=groups),
    width = 0.0,
    position = position_dodge(width = 0.75),
    color = "black",
    inherit.aes = FALSE
  ) +
  
  scale_shape_manual(values = c(1,  2,16, 17)) +
  scale_color_manual(values = c(treat_colors))+
  labs(x = "Treatment Groups", y = "Ocular Temperature (C)",
       color = "Treatment Group", shape = "Fever Phase") +
  theme(
    axis.text.x = element_text(size=12, angle=45, hjust=1),
    strip.text = element_text(size=12)
  )

lm2_anova_III

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


#We expect that males have more pathology when infected
  #So subset to only inoculated birds to look at effect of sex and temp and disease on fever
  #Subset to 

ti.long.i <- ti.long %>% 
  filter(treatment=="Inoculated")

#What about disease?
d1 <- glmmTMB(fever_value ~ (fever_type + temp * sex) * ever_diseased + (1|band_number), data=ti.long.i)
d2 <- glmmTMB(fever_value ~ (fever_type + temp + sex) * ever_diseased + (1|band_number), data=ti.long.i)
d3 <- glmmTMB(fever_value ~ fever_type * temp * sex * ever_diseased + (1|band_number), data=ti.long.i)
d4 <- glmmTMB(fever_value ~ fever_type + temp + sex + ever_diseased + (1|band_number), data=ti.long.i)
d5 <- glmmTMB(fever_value ~ fever_type * ever_diseased + sex + temp + (1|band_number), data=ti.long.i)
null <- glmmTMB(fever_value ~ 1 + (1|band_number), data=ti.long.i)

aictab(cand.set=list(d1, d2, d3, d4, d5, null), 
       modnames=c("d1", "d2", "d3", "d4", "d5", "null"))

#
lm2d <- glmmTMB(fever_value ~ fever_type * ever_diseased + sex + temp + (1|band_number), data=ti.long.i)
simulateResiduals(lm2d, plot=T)
summary(lm2d)
car::Anova(lm2d, type = "III")

lm2dis <- glmmTMB(fever_value ~ fever_type * ever_diseased + sex + temp + (1|band_number), data=ti.long)
lm2t <- glmmTMB(fever_value ~ fever_type * treatment + sex + temp + (1|band_number), data=ti.long)
lm2i <- glmmTMB(fever_value ~ fever_type * ever_infected + sex + temp + (1|band_number), data=ti.long)

#Ever diseased is a better predictor of fever_value than treatment (change lm2d data to ti.long)
aictab(cand.set=list(lm2dis, lm2t, lm2i), 
       modnames=c("lm2dis", "lm2t", "lm2i"))

#sex significant if you include controls
simulateResiduals(lm2dis, plot=T)
summary(lm2dis)
car::Anova(lm2dis, type="III")

#emmeans including ever_disaesed
#estimate fever_type means within each combination of temperature x sex x ever_diseased
emm <-emmeans(lm2d, ~fever_type | temp*sex*ever_diseased)
tests <- contrast(emm, method = "pairwise", adjust = "tukey")


emm_df <- as.data.frame(emm)  # has emmean, SE, lower.CL, upper.CL

ti.long <- ti.long %>%
  mutate(
    sex_phase = interaction(sex, fever_type),
    sex_phase = recode_factor(
      sex_phase,
      "Female.Baseline" = "Female Baseline",
      "Male.Baseline" = "Male Baseline",
      "Female.Peak"     = "Female Peak",
      "Male.Peak"     = "Male Peak"
    )
  )

emm_df <- emm_df %>%
  mutate(
    sex_phase = interaction(sex, fever_type),
    sex_phase = recode_factor(
      sex_phase,
      "Female.Baseline" = "Female Baseline",
      "Male.Baseline" = "Male Baseline",
      "Female.Peak"     = "Female Peak",
      "Male.Peak"     = "Male Peak"
    )
  )

#Plot Emmeans sex comparison
ggplot(ti.long.i, aes(x = temp,
                    y = fever_value,
                    color = temp,
                    shape = sex_phase)) +
  
  geom_jitter(
    position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.9),
    alpha = 0.75, size = 2.5, stroke=1
  ) +
  # model means — DODGE BY fever_type
  geom_point(
    data = emm_df,
    aes(x = temp, y = emmean, group = sex_phase, shape = sex_phase, color=groups),
    position = position_dodge(width = 0.9),
    size = 3, 
    stroke = 1,
    color = "black",
    inherit.aes = FALSE
  ) +
  geom_errorbar(
    data = emm_df,
    aes(x = temp, ymin = lower.CL, ymax = upper.CL,
        group = sex_phase, color=groups),
    width = 0.0,
    position = position_dodge(width = 0.9),
    color = "black",
    inherit.aes = FALSE
  ) +
  
  scale_shape_manual(values = c(1, 2, 16, 17)) +
  scale_color_manual(values = c(temp_colors))+
  labs(x = "Temperature Groups", y = "Ocular Temperature (C)",
       color = "Temperature Groups", shape = "Fever Phase") +
  theme(
    axis.text.x = element_text(size=12, angle=45, hjust=1),
    strip.text = element_text(size=12)
  )+
  facet_grid(~ever_diseased, 
             labeller=labeller(
               ever_diseased = c('0' = "Never Diseased", '1' = "Ever Diseased")))




####Fever 3) Magnitude of Fever Change#####
hist(ti.fc$mag_high)

#Question: Does the magnitude of change differ based on treatment, temperature, or sex
fh1 <- lm(mag_high ~ temp*treatment, data=ti.fc)
fh2 <- lm(mag_high ~ temp+treatment, data=ti.fc)
fh2.1 <- lm(mag_high ~ temp+treatment+sex, data=ti.fc)
fh2.2 <- lm(mag_high ~ temp+treatment*sex, data=ti.fc)
null_fh <- lm(mag_high ~ 1, data=ti.fc)


aictab(cand.set=list(fh1, fh2, fh2.1, fh2.2, null_fh), 
       modnames=c("fh1", "fh2", "fh2.1", "fh2.2", "null"))

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
    groups = paste(temp, ifelse(treatment == "Inoculated", "Inoculated", "Control"))
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

#Is this all just being driven by who gets disease?
ti.fci <- ti.fc %>%
  filter(treatment == "Inoculated")
ti.fci$ever_diseased <- as.factor(ti.fci$ever_diseased)

ti.fci <- ti.fci %>%
  dplyr::mutate(ever_diseased = dplyr::recode(ever_diseased,
                              "0" = "Never Diseased",
                              "1" = "Ever Diseased"))

fi1 <- lm(mag_high ~ ever_diseased, data=ti.fci)
fi1.1 <- lm(mag_high ~ ever_diseased + sex, data=ti.fci)
fi1.2 <- lm(mag_high ~ ever_diseased * sex, data=ti.fci)
fi2 <- lm(mag_high ~ temp*ever_diseased, data=ti.fci)
fi3 <- lm(mag_high ~ temp*ever_diseased+sex, data=ti.fci)
fi4 <- lm(mag_high ~ temp*ever_diseased*sex, data=ti.fci)
fi5 <- lm(mag_high ~ temp+ever_diseased, data=ti.fci)
fi6 <- lm(mag_high ~ temp+ever_diseased+sex, data=ti.fci)
null_fi <- lm(mag_high ~ 1, data=ti.fci)

aictab(cand.set=list(fi1, fi1.1, fi1.2, fi2, fi3, fi4, fi5, fi6, null_fi), 
       modnames=c("fi1", "fi1.1", "fi1.2", "fi2", "fi3", "fi4", "fi5", "fi6", "null"))

summary(fi1.1)

simulateResiduals(fi1, plot=T)
hist(resid(fi1))
summary(fi1)
car::Anova(fi1, type="II")

#Compare to treatment
fc1 <- lm(mag_high ~ treatment, data=ti.fc)
fc2 <- lm(mag_high ~ ever_diseased, data = ti.fc)
fc3 <- lm(mag_high ~ ever_infected, data = ti.fc)

#ever_diseased best model
aictab(cand.set=list(fc1, fc2, fc3), 
       modnames=c("fc1", "fc2", "fc3"))

#emmeans
emm <-emmeans(fi1, ~ever_diseased)
tests <- contrast(emm, method = "pairwise", adjust = "tukey")


emm_df <- as.data.frame(emm)  # has emmean, SE, lower.CL, upper.CL

emm_df$ever_diseased <- as.factor(emm_df$ever_diseased)

#Ever diseased predicts fever change magnitude
ggplot(ti.fci, aes(x=ever_diseased, y=mag_high, color = temp))+
  geom_jitter(size=2.5, alpha=0.75, aes(shape=sex, group= sex),
              position=position_jitterdodge(dodge.width=0.25, jitter.width = 0.1, jitter.height = 0))+
  geom_point(data=emm_df, aes(y=emmean, x=ever_diseased), color="black", size=3)+
  geom_errorbar(data=emm_df, aes(y=emmean, ymin=lower.CL, ymax=upper.CL, x=ever_diseased), width=0., color="black", size=.75)+
  labs(x="Disease Category", y= "Magnitude of Fever Change (Peak - Baseline [C])", color="Temperature", shape = "Sex")+
  scale_color_manual(values=c(temp_colors))+
  theme(
    axis.text.x = element_text(size=13),
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

#####Eye Score####
source("r_scripts/dataCleaning_TI22.R")

range(ti.cont$quantity) #highest control quantity = 95.38
ti$quant_cutoff = 50
ti$seropos_cutoff = 0.061
ti$sympt_cutoff = 0.1

#Sample sizes
ti %>%
  filter(dpi ==0)%>%
  dplyr::select(treatment, temp, groups)%>%
  tbl_summary(
  )%>%
  modify_header(
    label ~ "**All Birds**"
  )

ti %>%
  dplyr::select(dpi, treatment, temp, groups)%>%
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
source("r_scripts/dataCleaning_eyeScore.R")


####Eye Score: Ordinal Mixed Effects Model####
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




#####Eye score: GAM####
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
####Eye Score: Zero-Inflated Gamma Model####
ti.inoc <- ti.inoc %>%
  mutate(
    temp      = relevel(temp, ref = "Warm"),
    sex = relevel(sex, ref = "F") 
  )

unique(ti.inoc$dpi)

lmes <- glmmTMB(total_eye_score ~ temp*dpi + (1|band_number),
                data=ti.inoc, 
                ziformula = ~ temp,
                family = ziGamma(link = "log"))


simulateResiduals(lmes, plot=T)
summary(lmes)


lm1 <- glmmTMB(total_eye_score ~ temp*dpi + (1|band_number),
               data=ti.inoc, 
               ziformula = ~ temp,
               family = ziGamma(link = "log"))

AIC(lmes, lm1)

library(multcomp)
means <- cld(emmeans(lm1, pairwise ~ dpi*temp, adjust = "tukey", type = "response"))

# Obtain estimated marginal means (emmeans) with pairwise comparisons
emmeans_results <- emmeans(lm1, pairwise ~ dpi * temp, adjust = "tukey", type = "response")

# Extract the 'emmeans' part and apply cld()
means <- cld(emmeans_results$emmeans)

# Convert to data frame if needed
means_df <- as.data.frame(means)

ggplot(means_df, aes(x = dpi, y = response, color = fct_rev(temp)))+
  geom_jitter(data = ti.mg, aes(x = dpi, y = tes), 
              width = 0.1, alpha = 0.5)+
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),     
                width = 0.1)+
  stat_summary(data= ti.mg, aes(x=dpi, y=tes, group= temp), geom="line", fun = "mean")+
  #geom_line(aes(group = temp), linewidth = 1.2)+
  geom_point(size = 3)+
  theme_bw()#+
facet_wrap(~temp)


#zero-inflated model allowing for temperature to dictate the probability of zero.
#for the non-zero birds, we used a gamma model with a log link function, including a two way interaction between temp and dpi and all lower order effects

#model selection
l1 <- glmmTMB(total_eye_score ~ temp*dpi + (1|band_number),
                 data=ti.inoc, 
                 ziformula = ~ temp,
                 family = ziGamma(link = "log"))

l2 <- glmmTMB(total_eye_score ~ temp*dpi + (1|band_number),
              data=ti.inoc, 
              ziformula = ~ sex,
              family = ziGamma(link = "log"))

l3 <- glmmTMB(total_eye_score ~ temp*dpi + sex + (1|band_number),
              data=ti.inoc, 
              ziformula = ~ temp,
              family = ziGamma(link = "log"))

l4 <- glmmTMB(total_eye_score ~ temp*dpi * sex + (1|band_number),
              data=ti.inoc, 
              ziformula = ~ temp,
              family = ziGamma(link = "log"))

aictab(cand.set=list(l1, l2, l3, l4), 
       modnames=c("l1", "l2", "l3", "l4"))

simulateResiduals(l2, plot=T)
summary(l2)
car::Anova(l2, type = "III")


# 1. No zero inflation (baseline)
m0 <- glmmTMB(
  total_eye_score ~ temp * poly(dpi, 2) + sex + (1|band_number),
  family = Gamma(link = "log"),
  data = ti.inoc
)

# 2. Simple zi, constant over birds
m1 <- glmmTMB(
  total_eye_score ~ temp * poly(dpi, 2) + sex + (1|band_number),
  ziformula = ~ 1,
  family = ziGamma(link = "log"),
  data = ti.inoc
)

# 3. Additive zi
m2 <- glmmTMB(
  total_eye_score ~ temp * poly(dpi, 2) + sex + (1|band_number),
  ziformula = ~ temp + poly(dpi, 2) + sex,
  family = ziGamma(link = "log"),
  data = ti.inoc
)

# 4. Interactive zi
m3 <- glmmTMB(
  total_eye_score ~ temp * poly(dpi, 2) + sex + (1|band_number),
  ziformula = ~ temp * poly(dpi, 2) + sex,
  family = ziGamma(link = "log"),
  data = ti.inoc
)

AIC(m1, m2, m3)


l1_poly <- glmmTMB(total_eye_score ~ temp * poly(dpi, 2) + (1|band_number),
              data=ti.inoc, 
              ziformula = ~ temp,
              family = ziGamma(link = "log"))

l2_poly <- glmmTMB(
  total_eye_score ~ sex + temp * poly(dpi, 2) + (1 | band_number),
  data      = ti.inoc,
  ziformula = ~ sex,
  family    = ziGamma(link = "log")
)

l3_poly <- glmmTMB(total_eye_score ~ temp* poly(dpi, 2) + sex + (1|band_number),
              data=ti.inoc, 
              ziformula = ~ temp,
              family = ziGamma(link = "log"))

l4_poly <- glmmTMB(total_eye_score ~ temp*poly(dpi, 2) * sex + (1|band_number),
              data=ti.inoc, 
              ziformula = ~ temp,
              family = ziGamma(link = "log"))

l5_poly <- glmmTMB(total_eye_score ~ sex + temp * poly(dpi, 2) + (1|band_number),
                   data=ti.inoc, 
                   ziformula = ~ temp * poly(dpi, 2) + sex, #zero inflation: 
                   family = ziGamma(link = "log"))

l6_poly <- glmmTMB(total_eye_score ~ sex * temp * poly(dpi, 2) + (1|band_number),
                   data=ti.inoc, 
                   ziformula = ~ temp * poly(dpi, 2) + sex, #zero inflation: 
                   family = ziGamma(link = "log"))

l7_poly <- glmmTMB(total_eye_score ~ sex + temp * poly(dpi, 2) + (1|band_number),
                   data=ti.inoc, 
                   ziformula = ~ temp + poly(dpi, 2) + sex, #zero inflation: 
                   family = ziGamma(link = "log"))

l8_poly <- glmmTMB(total_eye_score ~ sex * temp * poly(dpi, 2) + (1|band_number),
                   data=ti.inoc, 
                   ziformula = ~ temp + poly(dpi, 2) + sex, #zero inflation: 
                   family = ziGamma(link = "log"))

l9_poly <- glmmTMB(total_eye_score ~  temp * poly(dpi, 2) + (1|band_number),
                    data=ti.inoc, 
                    ziformula = ~ temp + poly(dpi, 2) + sex, #zero inflation: 
                    family = ziGamma(link = "log"))

l10_poly <- glmmTMB(total_eye_score ~  temp + poly(dpi, 2) + (1|band_number),
                   data=ti.inoc, 
                   ziformula = ~ temp + poly(dpi, 2) + sex, #zero inflation: 
                   family = ziGamma(link = "log"))

aictab(cand.set=list(l1_poly, l2_poly, l3_poly, l4_poly, l5_poly, l6_poly, l7_poly, l8_poly, l9_poly, l10_poly, l2), 
       modnames=c("l1_poly", "l2_poly", "l3_poly", "l4_poly", "l5_poly", "l6_poly", "l7_poly", "l8_poly", "l9_poly", "l10_poly", "l2"))

lp1 <- glmmTMB(total_eye_score ~ sex + temp * poly(dpi, 2) + (1|band_number),
                   data=ti.inoc, 
                   ziformula = ~ temp + poly(dpi, 2) + sex, #zero inflation: 
                   family = ziGamma(link = "log"))

lp2 <- glmmTMB(total_eye_score ~ sex + temp * poly(dpi, 2) + (1|band_number),
               data=ti.inoc, 
               ziformula = ~ temp + poly(dpi, 2), #zero inflation: 
               family = ziGamma(link = "log"))

lp3 <- glmmTMB(total_eye_score ~  temp * poly(dpi, 2) + (1|band_number),
               data=ti.inoc, 
               ziformula = ~ temp + poly(dpi, 2) , #zero inflation: 
               family = ziGamma(link = "log"))
lp4 <- glmmTMB(total_eye_score ~  temp * poly(dpi, 2) + (1|band_number),
               data=ti.inoc, 
               ziformula = ~ temp + poly(dpi, 2) + sex, #zero inflation: 
               family = ziGamma(link = "log"))

AIC(lp1, lp2, lp3, lp4) #keep sex

simulateResiduals(l9_poly, plot=T)
simulateResiduals(l7_poly, plot=T)

summary(l9_poly)
summary(l7_poly)
car::Anova(l9_poly, type="III")
car::Anova(l7_poly, type="III")

emm_poly <- emmeans(
  l7_poly,
  ~ temp | dpi*sex,
  at = list(dpi = sort(unique(ti.inoc$dpi))),
  component = "response", #averaging in zero component
  #component = "cond", #only infected individuals
  type = "response"
)

emm_poly_df <- as.data.frame(emm_poly)

#For emmeans, we can look at conditional model (severity among birds that have pathology), or overall mean

sex_colors <- c("#A76F98", "#578E3F")

ggplot() +
  # raw data
  geom_jitter(
    data = ti.inoc,
    aes(x = dpi, y = total_eye_score, color = sex, shape=sex),
    width = 0., height = 0.,
    alpha = 0.75, size = 2
  ) +
  geom_line(data = ti.inoc,
            aes(x=dpi, y=total_eye_score, color=sex, linetype = sex, groups=as.factor(band_number)),
            alpha=0.5)+
  # model means: component = response
  geom_line(
    data = emm_poly_df,
    aes(x = dpi, y = emmean, color = sex, linetype = sex),
    size = 0.9
  ) +
  # 95% CI around means: component = response
  geom_ribbon(
    data = emm_poly_df,
    aes(x = dpi,
        ymin = asymp.LCL, ymax = asymp.UCL, y = emmean,
        fill = sex),
    alpha = 0.25
  ) +
  # model means: component = cond
  # geom_line(
  #   data = emm_poly_df,
  #   aes(x = dpi, y = response, color = sex, linetype = sex),
  #   size = 0.9
  # ) +
  # # 95% CI around means: component = cond
  # geom_ribbon(
  #   data = emm_poly_df,
  #   aes(x = dpi,
  #       ymin = asymp.LCL, ymax = asymp.UCL, y = response,
  #       fill = sex),
  #   alpha = 0.25
  #) +
  labs(
    x = "Days Post Inoculation",
    y = "Ocular Pathology",
    color = "Sex",
    fill  = "Sex",
    shape = "Sex",
    linetype="Sex"
  )+
  scale_linetype_manual(values = c("dashed", "solid"))+
   scale_color_manual(values = sex_colors)+
   scale_fill_manual(values = sex_colors)+
  facet_grid(~temp)

emm_poly <- emmeans(
  l9_poly,
  ~ temp | dpi,
  at = list(dpi = sort(unique(ti.inoc$dpi))),
  #component = "response", #averaging in zero component
  component = "cond", #only infected individuals
  type = "response"
)

emm_poly_df <- as.data.frame(emm_poly)


#color temp
ggplot() +
  # raw data
  geom_jitter(
    data = ti.inoc,
    aes(x = dpi, y = total_eye_score, color = temp, shape = sex, group = sex),
    position=position_jitterdodge(jitter.height=0.1, jitter.width=0.1, dodge.width = 1),
    alpha = 0.75, size = 2
  ) +
  # geom_line(data = ti.inoc,
  #           aes(x=dpi, y=total_eye_score, color=temp, groups=as.factor(band_number)),
  #           alpha=0.1)+
  # model means: component = response
  geom_line(
    data = emm_poly_df,
    aes(x = dpi, y = response, color=temp),
    size = 0.9
  ) +
  # 95% CI around means: component = response
  geom_ribbon(
    data = emm_poly_df,
    aes(x = dpi,
        ymin = asymp.LCL, ymax = asymp.UCL, y = response,
        fill = temp),
        linetype="dotted",
        color="black",
        size=0.25,
    alpha = 0.15
  ) +
  labs(
    x = "Days Post Inoculation",
    y = "Ocular Pathology",
    color = "Temperature",
    fill  = "Temperature",
    shape = "Sex",
    linetype="Temperature"
  )+
  scale_shape_manual(
    values = c(F = 16, M = 17),
    labels = c(F = "Female", M = "Male")
  )+
  scale_color_manual(values = temp_colors)+
  scale_fill_manual(values = temp_colors)#+
  facet_wrap(~sex, labeller = labeller(sex = c('F'="Female", 'M' = "Male")))
  

#not polynomial - doesn't allow for changing shape 
emm_df <- emmeans(
  lm1,
  ~ temp | dpi,
  at = list(dpi = sort(unique(ti.inoc$dpi))),
  #component = "response",
  component = "cond", #
  type = "response"
)

emm_df <- as.data.frame(emm_df)

ggplot() +
  # raw data
  geom_jitter(
    data = ti.inoc,
    aes(x = dpi, y = total_eye_score, color = temp),
    width = 0., height = 0.,
    alpha = 0.75, size = 2
  ) +
  # model means
  geom_line(
    data = emm_df,
    aes(x = dpi, y = response, color = temp),
    size = 0.9
  ) +
  geom_line(data = ti.inoc,
            aes(x=dpi, y=total_eye_score, color=temp, groups=as.factor(band_number)),
            alpha=0.5)+
  # 95% CI around means
  geom_ribbon(
    data = emm_df,
    aes(x = dpi,
        ymin = asymp.LCL, ymax = asymp.UCL, y = response,
        fill = temp),
    alpha = 0.25
  ) +
  labs(
    x = "Days Post Inoculation",
    y = "Ocular Pathology",
    color = "Temperature",
    fill  = "Temperature",
    shape = "Temperature"
  )+
  scale_linetype_manual(values = c("dashed", "solid"))+
  scale_color_manual(values = temp_colors)+
  scale_fill_manual(values = temp_colors)+
  facet_wrap(~temp)

#####Max Eye Score####


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
