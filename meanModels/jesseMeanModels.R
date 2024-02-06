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
library(AICcmodavg)
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
                alpha = 0.5)+
  geom_jitter(data = ti.pi, aes(x = dpi, y = elisa_od, color = treatment),
              width = 0.1)+
  facet_wrap(~temp)+
  theme_bw()

# Fever Score

source("dataCleaning_fever.R")

## For this model, I used a linear mixed model to model and I tested all 
#     possible interactions and removed non-significant ones in a stepwise
#     selection process.

ti.f$dpi <- as.factor(ti.f$dpi)

lm2 <- lmer(fever_change ~ temp*treatment*dpi + (1|band_number),data=ti.f)

hist(residuals(lm2))
plot(residuals(lm2), predict(lm2))
car::Anova(lm2, type = "III")

### stepwise selection:
car::Anova(lm2, type = "III")

lm2 <- lmer(fever_change ~ temp + treatment + dpi + 
              treatment:dpi + temp:dpi + 
              temp:treatment + (1|band_number),data=ti.f)

car::Anova(lm2, type = "III")

##### Keep two-way interactions

means <- as.data.frame(cld(emmeans(lm2, pairwise ~ dpi:treatment, adjust = "tukey")))

ggplot(means, aes(x = dpi, y = emmean, color = treatment))+
  geom_line(aes(group = treatment))+
  geom_point()

means2 <- as.data.frame(cld(emmeans(lm2, pairwise ~ dpi:temp, adjust = "tukey")))

ggplot(means2, aes(x = dpi, y = emmean, color = temp))+
  geom_line(aes(group = temp))+
  geom_point()

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
  theme_bw()+
  facet_wrap(~temp)

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
  geom_jitter(data = ti, aes(x = treatment, y = phago_score, color = temp),
              width = 0.1, alpha = 0.5)+
  theme_bw()


# Eye Score

source("dataCleaning_eyeScore.R")

ti.mg.mod$band_number <- as.factor(ti.mg.mod$band_number)

lm1 <- glmer(tes ~ temp*dpi + (1|band_number),
             data=ti.mg.mod, 
             family = poisson,
             control = glmerControl(optimizer ='optimx', optCtrl=list(method='L-BFGS-B')))

AF1 <- allFit(lm1, verbose=F)
AF1_lliks <- sort(sapply(AF1,logLik))

# Model parameters do not depend on optimizer - all within thousandths of a decimal

bind_rows(AF1_lliks) %>%
  remove_rownames(.) %>%
  mutate(model = c("NB Mixed Effects Model")) %>%
  select(model, everything()) %>%
  group_by(model) %>%
  gather(., Optimizer, llik, 2:ncol(.)) %>%
  ggplot(.,aes(Optimizer, llik)) + geom_point() +
  facet_wrap(~model) + coord_flip() +
  ylab("Log-Likelihood") +
  labs(title = "The Log-Likelihoods of Seven Different Optimizers in Our Model")

simulateResiduals(lm1, plot = T)

car::Anova(lm1, type = "III")

summary(lm1)

means <- cld(emmeans(lm1, pairwise ~ dpi*temp, adjust = "tukey", type = "response"))

ggplot(means, aes(x = dpi, y = rate, color = temp))+
  geom_jitter(data = ti.mg, aes(x = dpi, y = tes), 
              width = 0.1, alpha = 0.5)+
  geom_line(aes(group = temp), linewidth = 1.2)+
  geom_point(size = 3)+
  theme_bw()

# Zoom In to Show DPI differences

means2 <- cld(emmeans(lm1, pairwise ~ dpi, adjust = "tukey", type = "response"))

ggplot(means2, aes(x = dpi, y = rate))+
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), 
                width = 0.1)+
  geom_line(aes(group = 1), linewidth = 1.2)+
  geom_point(size = 3)+
  theme_bw()


