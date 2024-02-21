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
meancc <- mean(coldcon$elisa_od, na.rm =T)
meanci <- mean(coldinf$elisa_od, na.rm =T)
meanwc <- mean(warmcon$elisa_od, na.rm =T)
meanwi <- mean(warminf$elisa_od, na.rm =T)

coldcon$res <- abs(coldcon$elisa_od - meancc)
coldinf$res <- abs(coldinf$elisa_od - meanci)
warmcon$res <- abs(warmcon$elisa_od - meanwc)
warminf$res <- abs(warminf$elisa_od - meanwi)

dat <- rbind(coldcon, coldinf, warmcon, warminf)


# Model

mod <- glm(res ~ temp*treatment*dpi, data = dat,
           family = Gamma())

hist(residuals(mod))
plot(residuals(mod), predict(mod))

car::Anova(mod, type = "III")

# Stepwise

modstep <- stepAIC(mod)

car::Anova(modstep, type = "III")

emmeans(modstep, pairwise ~ treatment, adjust = "tukey",
        type = "response")

ggplot(dat, aes(x = dpi, y = res, color = temp))+
  geom_jitter(width = 0.1)+
  facet_wrap(~treatment)


ggplot(dat, aes(x = elisa_od, y = res))+
  geom_point()



#### Eye Score

# Antibody Assay

source("dataCleaning_eyeScore.R")

unique(ti.mg$dpi)

ti.mg$dpi <- as.factor(ti.mg$dpi)
ti.mg$band_number <- as.factor(ti.mg$band_number)


#### Within-Bird Variability

# Make the different day eye score values their own columns

# Simplify data set so pivot_wider works
ti.mg <- ti.mg %>% select(dpi, band_number, treatment, temp,
                          sex, total_eye_score)

ti.mg <- ti.mg %>% na.omit(total_eye_score)

unique(ti.mg$dpi)

# Add tiny constant to values so I can use PVs
ti.mg$total_eye_score <- ti.mg$total_eye_score +0.01

wibird <- ti.mg %>% pivot_wider(names_from = dpi,
                                names_glue = "{.value}_{dpi}",
                                values_from = total_eye_score)





# Calcuate a PV for each bird

library(CValternatives)

pvs <- c()

for (i in 1:nrow(wibird)){
  
  vec <- c(wibird$`total_eye_score_-12`[i], wibird$total_eye_score_3[i],
           wibird$`total_eye_score_7`[i], wibird$total_eye_score_9[i],
           wibird$`total_eye_score_14`[i], wibird$total_eye_score_18[i],
           wibird$`total_eye_score_21`[i], wibird$total_eye_score_24[i],
           wibird$`total_eye_score_28`[i], wibird$total_eye_score_35[i])
  
  pvs[i] <- PV(vec)
  
}

wibird$pv <- pvs

#Graph to see how they differ

ggplot(wibird, aes(x = treatment, y = pv))+
  geom_jitter()+
  facet_wrap(~temp)

# Model PVs

hist(wibird$pv) #weird bimodal distribution... ZI necessary

# No random effect needed because only one PV per bird.

mod <- glmmTMB(pv ~ temp, data = wibird, family = ziGamma(),
               ziformula = ~temp)

summary(mod)

simulateResiduals(mod, plot = T) #Assumptions check out

car::Anova(mod, type = "III")

#### Between-Bird Variability

# New plan! Get error for each birds obs as it differs from group mean.

# Get rid of my tiny constant from before
ti.mg$total_eye_score <- ti.mg$total_eye_score - 0.01

# Make separate data sets
coldneg12 <- ti.mg %>% filter(temp == "Cold"& dpi == "-12")
cold3 <- ti.mg %>% filter(temp == "Cold"& dpi == "3")
cold7 <- ti.mg %>% filter(temp == "Cold"& dpi == "7")
cold9 <- ti.mg %>% filter(temp == "Cold"& dpi == "9")
cold14 <- ti.mg %>% filter(temp == "Cold"& dpi == "14")
cold18 <- ti.mg %>% filter(temp == "Cold"& dpi == "18")
cold21 <- ti.mg %>% filter(temp == "Cold"& dpi == "21")
cold24 <- ti.mg %>% filter(temp == "Cold"& dpi == "24")
cold28 <- ti.mg %>% filter(temp == "Cold"& dpi == "28")
cold35 <- ti.mg %>% filter(temp == "Cold"& dpi == "35")

warmneg12 <- ti.mg %>% filter(temp == "warm"& dpi == "-12")
warm3 <- ti.mg %>% filter(temp == "warm"& dpi == "3")
warm7 <- ti.mg %>% filter(temp == "warm"& dpi == "7")
warm9 <- ti.mg %>% filter(temp == "warm"& dpi == "9")
warm14 <- ti.mg %>% filter(temp == "warm"& dpi == "14")
warm18 <- ti.mg %>% filter(temp == "warm"& dpi == "18")
warm21 <- ti.mg %>% filter(temp == "warm"& dpi == "21")
warm24 <- ti.mg %>% filter(temp == "warm"& dpi == "24")
warm28 <- ti.mg %>% filter(temp == "warm"& dpi == "28")
warm35 <- ti.mg %>% filter(temp == "warm"& dpi == "35")

# Calculate means
meancneg12 <- mean(coldneg12$total_eye_score, na.rm =T)
meanc3 <- mean(cold3$total_eye_score, na.rm =T)
meanc7 <- mean(cold7$total_eye_score, na.rm =T)
meanc9 <- mean(cold9$total_eye_score, na.rm =T)
meanc14 <- mean(cold14$total_eye_score, na.rm =T)
meanc18 <- mean(cold18$total_eye_score, na.rm =T)
meanc21 <- mean(cold21$total_eye_score, na.rm =T)
meanc24 <- mean(cold2412$total_eye_score, na.rm =T)
meanc28 <- mean(cold28$total_eye_score, na.rm =T)
meanc35 <- mean(cold35$total_eye_score, na.rm =T)

meanwneg12 <- mean(warmneg12$total_eye_score, na.rm =T)
meanw3 <- mean(warm3$total_eye_score, na.rm =T)
meanw7 <- mean(warm7$total_eye_score, na.rm =T)
meanw9 <- mean(warm9$total_eye_score, na.rm =T)
meanw14 <- mean(warm14$total_eye_score, na.rm =T)
meanw18 <- mean(warm18$total_eye_score, na.rm =T)
meanw21 <- mean(warm21$total_eye_score, na.rm =T)
meanw24 <- mean(warm2412$total_eye_score, na.rm =T)
meanw28 <- mean(warm28$total_eye_score, na.rm =T)
meanw35 <- mean(warm35$total_eye_score, na.rm =T)

coldneg12$res <- abs(coldneg12$total_eye_score - meanc)
cold3$res <- abs(cold3$total_eye_score - meanc)
cold7$res <- abs(cold7$total_eye_score - meanc)
cold9$res <- abs(cold9$total_eye_score - meanc)
cold14$res <- abs(cold14$total_eye_score - meanc)
cold18$res <- abs(cold18$total_eye_score - meanc)
cold21$res <- abs(cold21$total_eye_score - meanc)
cold24$res <- abs(cold24$total_eye_score - meanc)
cold28$res <- abs(cold28$total_eye_score - meanc)
cold35$res <- abs(cold35$total_eye_score - meanc)

warmneg12$res <- abs(warmneg12$total_eye_score - meanc)
warm3$res <- abs(warm3$total_eye_score - meanc)
warm7$res <- abs(warm7$total_eye_score - meanc)
warm9$res <- abs(warm9$total_eye_score - meanc)
warm14$res <- abs(warm14$total_eye_score - meanc)
warm18$res <- abs(warm18$total_eye_score - meanc)
warm21$res <- abs(warm21$total_eye_score - meanc)
warm24$res <- abs(warm24$total_eye_score - meanc)
warm28$res <- abs(warm28$total_eye_score - meanc)
warm35$res <- abs(warm35$total_eye_score - meanc)


dat <- rbind(coldneg12, cold3, cold7, cold9, cold14,
             cold18, cold21, cold24, cold28, cold35,
             warmneg12, warm3, warm7, warm9, warm14,
             warm18, warm21, warm24, warm28, warm35)


# Model

hist(dat$res)

mod <- glm(res ~ temp+dpi, data = dat, family = Gamma())

summary(mod)

hist(residuals(mod))
plot(residuals(mod), predict(mod))

car::Anova(mod, type = "III")

