################ Jesse G-L Mean Models #####################
# 
# This file contains the "final" versions of Jesse's variability models
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
ti.pi$band_number <- as.factor(ti.pi$band_number)


#### Within-Bird Variability

# There are only two data points per bird so PV or CV may not be 
#   the best option. I am going to try those and try using
#   just a simple difference.

# Make the day 9 and 28 OD values their own columns

# Simplify data set so pivot_wider works
ti.pi <- ti.pi %>% select(dpi, band_number, treatment, temp,
                          sex, elisa_od)

wibird <- ti.pi %>% pivot_wider(names_from = dpi,
                                names_glue = "{.value}_{dpi}",
                                values_from = elisa_od)

wibird <- wibird %>% mutate(elisa_diff = abs(elisa_od_9 - elisa_od_28))

# Calcuate a PV for each bird

library(CValternatives)

pvs <- c()

for (i in 1:nrow(wibird)){
  
  vec <- c(wibird$elisa_od_9[i], wibird$elisa_od_28[i])
  
  pvs[i] <- PV(vec)
  
}

wibird$pv <- pvs

#Graph to see how they differ - not too different when absolute value of diff used

ggplot(wibird, aes(x = treatment, y = pv))+
  geom_point()+
  facet_wrap(~temp)

ggplot(wibird, aes(x = treatment, y = elisa_diff))+
  geom_point()+
  facet_wrap(~temp)

# Model PVs

hist(wibird$elisa_diff)
hist(wibird$pv) #looks pretty skewed - gamma regression probably ideal

# No random effect needed because only one PV per bird.

# Add teeny constant to the values to make the gamma work

wibird$pv_plus <- wibird$pv + 0.0001

mod <- glm(pv_plus ~ temp*treatment, data = wibird,
               family = Gamma())

summary(mod)

simulateResiduals(mod, plot = T) #Assumptions check out

car::Anova(mod, type = "III")

#### Between-Bird Variability

# Make filtered data sets

coldcon <- ti.pi %>% filter(temp == "Cold" & treatment == "Control")
warmcon <- ti.pi %>% filter(temp == "Warm" & treatment == "Control")
coldinf <- ti.pi %>% filter(temp == "Cold" & treatment == "Infected")
warminf <- ti.pi %>% filter(temp == "Warm" & treatment == "Infected")

# Fit models

modcc <- glmmTMB(elisa_od ~ 1 + (1|band_number), data = coldcon,
             family = Gamma(link = "log"))

modwc <- glmmTMB(elisa_od ~ 1 + (1|band_number), data = warmcon,
                 family = Gamma(link = "log"))

modci <- glmmTMB(elisa_od ~ 1 + (1|band_number), data = coldinf,
                 family = Gamma(link = "log"))

modwi <- glmmTMB(elisa_od ~ 1 + (1|band_number), data = warminf,
                 family = Gamma(link = "log"))

# Get variance components

vcCC<- ranef(modcc)
vcWC<- ranef(modwc)
vcCI<- ranef(modci)
vcWI<- ranef(modwi)

# Combine into data set

varcomps <- c(vcCC$cond$band_number$`(Intercept)`,
                  vcCI$cond$band_number$`(Intercept)`,
                  vcWC$cond$band_number$`(Intercept)`,
                  vcWI$cond$band_number$`(Intercept)`)

varcompsABS <- abs(varcomps)

btwbird <- data.frame(temp = c(rep("Cold",20), rep("Warm", 24)),
                      treatment = c(rep("Control", 6), rep("Infected", 14),
                                    rep("Control", 10), rep("Infected", 14)),
                      varcomp = varcompsABS)

head(btwbird)

hist(varcomps) # should be normalish because that's where they were drawn from

hist(varcompsABS) #now all positive - use gamma?

modbtw <- lm(varcomp ~ temp*treatment, data = btwbird)

hist(residuals(modbtw))
plot(residuals(modbtw), predict(modbtw))

summary(modbtw)

means <- emmeans(modbtw, pairwise ~ temp:treatment, adjust = "tukey",
        type ="response")

cld(means)


ggplot(btwbird, aes(x = temp, y = varcomp))+
  geom_jitter()+
  facet_wrap(~treatment)


ggplot(ti.pi, aes(x = temp, y = elisa_od))+
  geom_jitter(width = 0.1)+
  facet_wrap(~treatment)

# Weird stuff happening - keep getting mean VC estimates of 0

# New plan! Get error for each birds obs as it differs from group mean.

# Calculate means - need to add separate ones for each dpi if interested in that

medelisa  <- median(ti.pi$elisa_od, na.rm = T)
meanelisa  <- mean(ti.pi$elisa_od, na.rm = T)


ggplot(data = ti.pi, aes(x = dpi, y = elisa_od))+
  geom_jitter(width = 0.1, aes(color = temp))+
  geom_hline(yintercept = medelisa)+
  geom_hline(yintercept = meanelisa, color = "green")+
  facet_wrap(~treatment)

ti.pi$res <- abs(ti.pi$elisa_od - medelisa)


# Model - there is no random effect because it didn't help the model at all

mod <- glmmTMB(res ~ temp*treatment*dpi,
               data = ti.pi,
               family = Gamma())
summary(mod)

simulateResiduals(mod, plot = T)

car::Anova(mod, type = "III")

# Stepwise

modstep <- stepAIC(mod)

car::Anova(modstep, type = "III")

emmeans(modstep, pairwise ~ treatment, adjust = "tukey",
        type = "response")

emmeans(modstep, pairwise ~ temp:dpi, adjust = "tukey",
        type = "response")

#These results will take some nuance to understand

#### Eye Score - I think I need a categorical measure of variability

#Variability in categorical variables can be calculated as 
# 1-sum((proportions in eachgroup)^2)

# Antibody Assay

source("dataCleaning_eyeScore.R")

unique(ti.mg$dpi)

ti.mg$dpi <- as.factor(ti.mg$dpi)
ti.mg$band_number <- as.factor(ti.mg$band_number)


#### Within-Bird Variability

# Make a little table that shows counts for each bird in each category

wibird <- ti.mg %>% group_by(band_number, treatment, temp, total_eye_score) %>%
  summarize(n = n())

wibird <- wibird %>% na.omit(total_eye_score)

wibird <- wibird %>% group_by(band_number) %>% mutate(prop = n/sum(n)) %>%
  ungroup()

# Make even tinier table with the variability measure
#    - closer to 1 = more variable :)

wibird$propsq <- wibird$prop^2

wibirdVar <- wibird %>% group_by(band_number, treatment, temp) %>%
  summarize(catVar = 1-sum(propsq))

# Model variability

hist(wibirdVar$catVar) #weird bimodal distribution. Gonna use quasibinomial

# No random effect needed because only one PV per bird.

mod <- glm(catVar ~ temp, data = wibirdVar, family = quasibinomial())

summary(mod)

car::Anova(mod, type = "III")

#### Between-Bird Variability

# New plan! Get error for each birds obs as it differs from group mean.

# Add difference column
meanes <- mean(ti.mg$total_eye_score, na.rm = T)
medianes <- median(ti.mg$total_eye_score, na.rm = T)

ti.mg$res <- abs(ti.mg$total_eye_score-meanes)

# Model

hist(ti.mg$res)

range(ti.mg$res, na.rm =T)

mod <- glmmTMB(res ~ temp+dpi + (1|band_number),
               data = ti.mg, 
               family = Gamma())

# Tried AR1 covariance structure and it was too complicated to calculate
# Also tried zi model and nothing was significant and fit was no better

summary(mod)

simulateResiduals(mod, plot = T) #Pretty weird residuals
# I think this is because the outcome is actually categorical.

car::Anova(mod, type = "III")

cld(emmeans(mod, pairwise ~ dpi, adjust = "tukey", type = "response"))

# All in all, this analysis is weird and really since this is a poisson variable
#   the mean tells us about the variance too. I think this is giving us sensible
#   results, so I think it's fine to use as long as we lean on the conservative
#   side for stating results.

###### Fever
source("dataCleaning_fever.R")

#### Within-Bird Variability

# Make the different day eye score values their own columns

# Simplify data set so pivot_wider works
ti.f <- ti.f %>% select(dpi, band_number, treatment, temp,
                          sex, fever_score)

unique(ti.mg$dpi)

# Add tiny constant to values so I can use PVs

wibird <- ti.f %>% pivot_wider(names_from = dpi,
                                names_glue = "{.value}_{dpi}",
                                values_from = fever_score)

# Calcuate a PV for each bird

library(CValternatives)

pvs <- c()

for (i in 1:nrow(wibird)){
  
  vec <- c(wibird$`fever_score_0`[i], wibird$fever_score_3[i],
           wibird$`fever_score_7`[i], wibird$fever_score_14[i],
           wibird$`fever_score_18`[i], wibird$fever_score_24[i],
           wibird$`fever_score_28`[i], wibird$fever_score_35[i])
  
  pvs[i] <- PV(vec)
  
}

wibird$pv <- pvs

#Graph to see how they look

ggplot(wibird, aes(x = treatment, y = pv))+
  geom_jitter()+
  facet_wrap(~temp)


# Model PVs

hist(wibird$pv)

# No random effect needed because only one PV per bird.

mod <- glmmTMB(pv ~ temp+treatment, data = wibird, family = Gamma())

summary(mod)

simulateResiduals(mod, plot = T) #Assumptions check out

car::Anova(mod, type = "III")

#### Between-Bird Variability

# New plan! Get error for each birds obs as it differs from group mean.

# Add difference column
meanfev <- mean(ti.f$fever_score, na.rm = T)

ti.f$res <- abs(ti.f$fever_score-meanfev)

# Model

hist(ti.f$res)

ti.f$dpi <- as.factor(ti.f$dpi)

mod <- glmmTMB(res ~ temp*treatment*dpi + (1|band_number),
               data = ti.f, family = gaussian())

summary(mod)

hist(residuals(mod))
plot(residuals(mod), predict(mod)) #Assumptions look okay

car::Anova(mod, type = "III") # Gotta keep the fat model

cld(emmeans(mod, pairwise ~ temp:treatment:dpi, adjust = "tukey", type = "response"))
#Godspeed....
# JK I'll help you sort through these lol



############### Some notes ########################
#
# I think we definitely have the data to talk about within bird
# variability. For between bird variability, we *technically*
# don't have replication, so it is a little weird to use tests.
# That being said, I don't see anything wrong with reporting these
# and just being conservative about reporting the results.
#