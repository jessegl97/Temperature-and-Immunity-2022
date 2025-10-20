################ Jesse G-L Mean Models #####################
# 
# This file contains the "final versions of Jesse's mean models
#   for the temperature and immunity project. Each modeling section
#   contains details about the models fit and why they were chosen.
#   There are also graphs to use as a starting place for publication
#   ready figures if needed/ wanted.
#
############################################################
#
# Set-Up
#
# I used the beginning sections of each code file to pull together
#   data set up files that can be sourced into this one for cleaner
#   reading. See the sourced files to make changes to data cleaning
#   procedure.

library(emmeans)
library(multcomp)
library(ggplot2)
library(tidyverse)
library(lme4)
library(glmmTMB)
library(effects)
#library(AICcmodavg)
library(DHARMa)
library(optimx)

# Antibody Assay

source("dataCleaning_antibody.R")

unique(ti.pi$dpi)

ti.pi$dpi <- as.factor(ti.pi$dpi)

# For this model, I started with the three-way interaction of treatment,
#    temp, and dpi and then removed higher level interactions sequentially
#    when they were not significant. The final result only had one two-way
#    interaction in it. I chose the inverse link over the log link based
#    on BIC metrics. There is some evidence of heteroscedasticity in the
#    DHARMa residuals, but the model we have is the best one I could figure.
#    I think it has to do with how small the values are and the fact that 
#    we don't have a ton of data. I may think about this more later, but
#    I am comfortable that the results are reasonable given the data.

ti.pi$dpif <- as.factor(ti.pi$dpi)
lm1<-glm(elisa_od~treatment+temp+dpi + treatment:dpi,
             data=ti.pi, 
             family=Gamma(link = "inverse"))

simulationOutput <-simulateResiduals(lm1, plot = T)

car::Anova(lm1, type = "III")
summary(lm1)

means <- as.data.frame(cld(emmeans(lm1, pairwise ~ treatment:dpi, type = "response",
                                   adjust = "tukey")))


ggplot(means, aes(x = dpi, y = response))+
  geom_crossbar(aes(ymin = lower.CL, ymax = upper.CL, fill = treatment),
                alpha = 0.5, width = 0.75)+
  geom_jitter(data = ti.pi, aes(x = dpi, y = elisa_od, color = treatment),
              width = 0.1)+
  scale_color_manual(values = c("#F4A460", "#9370DB"))+
  scale_fill_manual(values = c("#F4A460", "#9370DB"))+
  labs(color = "Treatment", fill = "Treatment", x="Days Post Infection", y="Antibody Levels")+
  facet_wrap(~temp)+
  theme_bw()

# Fever Score

source("dataCleaning_fever.R")

## For this model, I used a linear mixed model to model and I tested all 
#     possible interactions and removed non-significant ones in a stepwise
#     selection process.

ti.f$dpi <- as.factor(ti.f$dpi)

lm2 <- glmmTMB(fever_change ~ temp*treatment*dpi + 
                 ar1(dpi +0 |band_number),
                data=ti.f)

lm2og <- glmmTMB(fever_change ~ temp*treatment*dpi + 
                 (1|band_number),
               data=ti.f)

BIC(lm2, lm2og) #AR1 covariance structure better; BIC less biased small sample sizes, AICc also works

hist(residuals(lm2og))
plot(residuals(lm2), predict(lm2))
car::Anova(lm2, type = "III")

### stepwise selection:
car::Anova(lm2, type = "III")

lm2 <- glmmTMB(fever_change ~ temp + treatment + dpi + 
              treatment:dpi + temp:dpi + ar1(dpi + 0|band_number),
            data=ti.f)
#p = 1.00 = low effect size 
car::Anova(lm2, type = "III")
summary(lm2)
#lm2 shows difference b/t infected and control, then use lm2a to look for differences b/t temps

##### Keep two-way interactions
#look only at infected birds
lm2a <- glmmTMB(fever_change ~ temp * dpi + ar1(dpi + 0|band_number), data=ti.f %>%filter(treatment=="Infected"))
summary(lm2a)
means <- as.data.frame(cld(emmeans(lm2a, pairwise ~ dpi:temp, adjust = "tukey")))

ggplot(means, aes(x = dpi, y = emmean, color = temp))+
  geom_line(aes(group = temp))+
  scale_color_manual(values = c("#F4A460", "#9370DB"))+
  scale_fill_manual(values = c("#F4A460", "#9370DB"))+
  geom_point()
#no difference in kinetics

ti.f$pred <- predict(lm2)

ti.fMeans <- ti.f %>% group_by(dpi, treatment, temp) %>%
  summarize(meanFevCh = mean(pred))

ggplot(ti.fMeans, aes(x = dpi, y = meanFevCh))+
  geom_line(aes(group = treatment, color = treatment),
            linewidth = 1.2)+
  geom_point(aes(color = treatment), size = 3)+
  geom_jitter(data = ti.f, aes(x = dpi, y = fever_change, color = treatment), 
              alpha = 0.2,
              width = 0.1)+
  scale_color_manual(values = c("#F4A460", "#9370DB"))+
  scale_fill_manual(values = c("#F4A460", "#9370DB"))+
  labs(y="Mean Fever Change", x= "Days Post Infection", color="Treatment")+
  theme_bw()+
  facet_wrap(~temp, nrow=2)

# Phagocytosis Assay

source("dataCleaning_phago.R")

## Jesse already had this model in great shape! It is a binomial model and the
#    weights represent the total number of trials, since the response is a
#    pre-calculated proportion. The interaction between temp and treatment was
#    not significant.

glm1 <- glmer(phago_score~temp+treatment + (1|band_number), 
              weights=wbc_total+phago_total, 
              data=ti, family="binomial")

simulateResiduals(glm1, plot = T)

summary(glm1)

means <- as.data.frame(cld(emmeans(glm1, pairwise ~ treatment,
                                   type = "response", adjust = "tukey")))

ggplot(means, aes(x = treatment, y = prob))+
  geom_crossbar(aes(ymin = asymp.LCL, ymax = asymp.UCL))+
  geom_jitter(data = ti, aes(x = treatment, y = phago_score, color = fct_rev(temp)),
              width = 0.1, alpha = 0.5, size=2)+
  labs(x="Treatment", y="Phagocytosis Score", color="Temperature")+
  theme_bw()


# Eye Score

source("dataCleaning_eyeScore.R")

ti.mg.mod$band_number <- as.factor(ti.mg.mod$band_number)

lm1 <- glmmTMB(tes ~ temp*dpi + (1|band_number),
             data=ti.mg.mod, 
             ziformula = ~ temp,
             family = ziGamma(link = "log"))
#zero-inflated model allowing for temperature to dictate the probability of zero.
#for the non-zero birds, we used a gamma model with a log link function, including a two way interaction between temp and dpi and all lower order effects

lm1ar<- glmmTMB(tes ~ temp*dpi + ar1(dpi + 0|band_number),
                data=ti.mg.mod, 
                ziformula = ~1,
                family = ziGamma(link = "log"))

BIC(lm1, lm1ar) #Difference < 2 effectively the same

simulateResiduals(lm1, plot = T)

car::Anova(lm1, type = "III")
#zero-inflated model: (logistic regression model) there was a tendency towards warm birds being more likely to develop eye scores than cold birds (p = 0.078, etc.)
#conditional model: (gamma model) the interaction effect was not significant, but was informative to the model (AIC better w/ it included). DPI was significant, etc.
summary(lm1)

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

#model predictions
newdata <- expand.grid(
  temp = unique(ti.mg.mod$temp),
  dpi = unique(ti.mg.mod$dpi),
  band_number = unique(ti.mg.mod$band_number)
)

# Make predictions using the fitted model
newdata$pred_tes <- predict(lm1, newdata = newdata, type = "response", re.form=NA)
# Zoom In to Show DPI differences combine temps w/ error

means2 <- cld(emmeans(lm1, pairwise ~ dpi, adjust = "tukey", type = "response"))

ggplot(means2, aes(x = dpi, y = response))+
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),     
                width = 0.1)+
  geom_line(aes(group = 1), linewidth = 1.2)+
  geom_point(size = 3)+
  theme_bw()


