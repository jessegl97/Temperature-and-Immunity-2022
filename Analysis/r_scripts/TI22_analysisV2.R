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
sex_colors <- c("#A76F98", "#578E3F")
sex_shapes <- c(
  "Female" = 16,
  "Male" = 17
)

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


####Sex Effects on Disease####
source("r_scripts/dataCleaning_TI22.R")

ti %>%
  filter(dpi == 0)%>%
  dplyr::select(groups, ever_infected, ever_diseased, sex)%>%
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

#Do temperature or sex affect probability of birds developing disease when inoculated?
da <- glm(ever_diseased ~ sex + temp, data=ti.inoc, family=binomial())
db <- glm(ever_diseased ~ sex * temp, data=ti.inoc, family=binomial())
dc <- glm(ever_diseased ~ sex, data=ti.inoc, family=binomial())
dd <- glm(ever_diseased ~ temp, data=ti.inoc, family=binomial())
dn <- glm(ever_diseased ~ 1, data=ti.inoc, family=binomial())

#Null model is the best so none of these really make sense
aictab(cand.set=list(da, db, dc, dd, dn), 
       modnames=c("da", "db", "dc", "dd", "null"))

# lm.dis <- glm(ever_diseased ~ sex + temp, data=ti.inoc, family=binomial())
# simulateResiduals(lm.dis, plot=T)
# summary(lm.dis)
# car::Anova(lm.dis, type = "II")

#glm is not better than null model: so try Fisher's Exact Test for each temperature
ftd.m <- data.frame(
  "Cold" = c(0, 4, 3),
  "Warm" = c(0, 1, 6),
  row.names = c("Uninfected", "Asymptomatic", "Diseased"),
  stringsAsFactors = FALSE)
colnames(ftd.m)

mosaicplot(ftd.m,
           main = "Mosaic",
           color = TRUE)

chisq.test(ftd.m)$expected

#Is there a difference in the breakdown between statuses in males between temperatures?
ftm <- fisher.test(ftd.m)
ftm

ftd.f <- data.frame(
  "Cold" = c(2, 2, 3),
  "Warm" = c(1, 4, 2),
  row.names = c("Uninfected", "Asymptomatic", "Diseased"),
  stringsAsFactors = FALSE)
colnames(ftd.f)

mosaicplot(ftd.f,
           main = "Mosaic",
           color = TRUE)

chisq.test(ftd.f)$expected

#Is there a difference in the breakdown between statuses in females between temperatures?
ftf <- fisher.test(ftd.f)
ftf

#Within temperatures
ftd.c <- data.frame(
  "Female" = c(2, 2, 3),
  "Male" = c(0, 4, 3),
  row.names = c("Uninfected", "Asymptomatic", "Diseased"),
  stringsAsFactors = FALSE)
colnames(ftd.c)

mosaicplot(ftd.c,
           main = "Mosaic",
           color = TRUE)

chisq.test(ftd.c)$expected

#Is there a difference in the breakdown between statuses between sexes in the cold room?
ftc <- fisher.test(ftd.c)
ftc

#Within temperatures - warm
ftd.w <- data.frame(
  "Female" = c(1, 4, 2),
  "Male" = c(0, 6, 1),
  row.names = c("Uninfected", "Asymptomatic", "Diseased"),
  stringsAsFactors = FALSE)
colnames(ftd.w)

mosaicplot(ftd.w,
           main = "Mosaic",
           color = TRUE)

chisq.test(ftd.w)$expected

#Is there a difference in the breakdown between statuses between sexes in the cold room?
ftw <- fisher.test(ftd.w)
ftw

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
  filter(treatment == "Inoculated" & dpi %in% c(0, 3, 7,9,14,18,21,24,28,35))

ti.inoc %>%
  dplyr::select(dpi, treatment, temp, groups)%>%
  tbl_summary(
    by=dpi
  )%>%
  modify_header(
    label ~ "**All Days Inoculated**"
  )

#Sample sizes for inoculated
  #Warm Inoculated n = 14 per DPI
  #Cold Inoculated n = 14 per DPI


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


ti.inoc$dpi.f <- as.factor(ti.inoc$dpi)
ti.inoc$temp <- as.factor(ti.inoc$temp)
ti.inoc$sex <- as.factor(ti.inoc$sex)
####Eye Score: Zero-Inflated Gamma Model####
ti.inoc <- ti.inoc %>%
  mutate(
    temp  = relevel(temp, ref = "Warm"),
    sex = relevel(sex, ref = "Female") 
  )

#Remove DPI 0, 3, and 35 because they lack any variability in eye score
ti.inoc.mod <- ti.inoc %>%
  filter(!dpi.f %in% c(0, 3, 35))

unique(ti.inoc.mod$dpi)

#zero-inflated model allowing for temperature or sex or their interaction to dictate the probability of zero.
#for the non-zero birds, we used a gamma model with a log link function, including a two way interaction between temp and dpi and all lower order effects

#model selection
l1 <- glmmTMB(total_eye_score ~ temp*dpi.f + (1|band_number),
              data=ti.inoc.mod, 
              ziformula = ~ temp,
              family = ziGamma(link = "log"))

l2 <- glmmTMB(total_eye_score ~ temp*dpi.f + (1|band_number),
              data=ti.inoc.mod, 
              ziformula = ~ sex,
              family = ziGamma(link = "log"))

l3 <- glmmTMB(total_eye_score ~ temp*dpi.f + sex + (1|band_number),
              data=ti.inoc.mod, 
              ziformula = ~ temp,
              family = ziGamma(link = "log"))

l4 <- glmmTMB(total_eye_score ~ temp*dpi.f * sex + (1|band_number),
              data=ti.inoc.mod, 
              ziformula = ~ temp,
              family = ziGamma(link = "log"))

l5 <- glmmTMB(total_eye_score ~ dpi.f + (1|band_number),
              data=ti.inoc.mod, 
              ziformula = ~ sex * temp,
              family = ziGamma(link = "log"))

l6 <- glmmTMB(total_eye_score ~ dpi.f + temp + (1|band_number),
              data=ti.inoc.mod, 
              ziformula = ~ sex * temp,
              family = ziGamma(link = "log"))

l7 <- glmmTMB(total_eye_score ~ dpi.f + (1|band_number),
              data=ti.inoc.mod, 
              ziformula = ~ temp,
              family = ziGamma(link = "log"))

l8 <- glmmTMB(total_eye_score ~ dpi.f + (1|band_number),
              data=ti.inoc.mod, 
              ziformula = ~ sex,
              family = ziGamma(link = "log"))

aictab(cand.set=list(l1, l2, l3, l4, l5, l6, l7, l8), 
       modnames=c("l1", "l2", "l3", "l4", "l5", "l6", "l7", "l8"))

simulateResiduals(l5, plot=T)
summary(l5)
car::Anova(l5, type = "III")

#No effect of temperature
simulateResiduals(l6, plot=T)
summary(l6)
car::Anova(l6, type = "III")

glm.zi.es <-glmmTMB(total_eye_score ~ dpi.f + (1|band_number),
        data=ti.inoc.mod, 
        ziformula = ~ sex * temp,
        family = ziGamma(link = "log"))

simulateResiduals(glm.zi.es, plot=T)

summary(glm.zi.es)
car::Anova(glm.zi.es, type="III")

emm <- emmeans(
  glm.zi.es,
  ~ temp | dpi.f,
  at = list(dpi = sort(unique(ti.inoc$dpi))),
  component = "response", #averaging in zero component
  #component = "cond", #only infected individuals
  type = "response"
)

emm_df <- as.data.frame(emm)

#merge back dpi 0, 3, and 35 for emmean visualization
emm_df <- emm_df %>%
  mutate(
    dpi.f = factor(dpi.f, levels = levels(ti.inoc$dpi.f)),
    temp  = factor(temp,  levels = levels(ti.inoc$temp))
  )

full_grid <- expand_grid(
  temp  = levels(ti.inoc$temp),
  dpi.f = levels(ti.inoc$dpi.f)
)

emm_df_full <- full_grid %>%
  left_join(emm_df, by = c("temp", "dpi.f"))

dodge = position_dodge(width = 0.2)

g.tes <- ggplot(ti.inoc, aes(x = dpi.f, y = total_eye_score, color = temp)) +
  # raw data points
  # geom_point(
  #   alpha  = 0.5,
  #   size   = 2, 
  #   position = dodge,
  #   aes(shape = sex, group = temp)
  # ) +
  #Jittered
  geom_jitter(
    alpha  = 0.5,
    size   = 2, 
    position=position_jitterdodge(jitter.width = 0.1, jitter.height = 0.1, dodge.width = 0.2),
    aes(shape = sex, group = temp, color=temp)
  ) +
  # individual bird trajectories
  # geom_line(
  #   aes(group = interaction(band_number, temp)),
  #   alpha = 0.25,
  #   position=dodge
  # ) +
  # 95% CI around means
  geom_errorbar(
    data = emm_df_full,
    aes(
      x    = dpi.f,
      ymin = asymp.LCL,
      ymax = asymp.UCL,
      color = temp,
      group = fct_rev(temp)
    ),
    position = dodge,
    alpha        = 1,
    inherit.aes  = FALSE,
    na.rm        = TRUE,
    width=0.1
  ) +
  #model means (component = "response")
  geom_point(
    data = emm_df_full,
    aes(x = dpi.f, y = emmean, color = temp, group = fct_rev(temp)),
    size = 3,
    position=dodge,
    na.rm = TRUE
  ) +
  geom_point(
    data = emm_df_full,
    aes(x = dpi.f, y = emmean, group = fct_rev(temp)),
    size = 3.2,
    color="black",
    shape=1,
    position=dodge,
    na.rm = TRUE
  ) +
  geom_line(
    data = emm_df_full,
    aes(x = dpi.f, y = emmean, color = temp, group = fct_rev(temp)),
    linewidth = 0.5,
    position=dodge,
    na.rm = TRUE
  ) +
  labs(
    x     = "Days Post Inoculation",
    y     = "Total Eye Score",
    color = "Temperature",
    fill  = "Temperature",
    shape = "Sex",
    title = "Gamma Model"
  ) +
  scale_color_manual(values = temp_colors) +
  scale_fill_manual(values  = temp_colors)+
  scale_shape_manual(values = sex_shapes)


# Zero-inflation component: probability of extra zeros
emm_zi <- emmeans(
  glm.zi.es,
  ~ sex * temp,
  component = "zi",    
  type      = "response"
)

emm_zi_df <- as.data.frame(emm_zi)

# emm_zi_df will have columns like:
# sex, temp, prob, SE, asymp.LCL, asymp.UCL, etc.

g.zi<- ggplot(emm_zi_df,
       aes(x = temp, y = response, color = sex, group = sex, shape = sex)) +
  geom_point(size = 3) +
  geom_line() +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),
                width = 0.03) +
  labs(
    x = "Temperature",
    y = "Probability of No Pathology [Pr(structural zero)]",
    color = "Sex",
    shape = "Sex",
    title = "Zero-inflation Probability"
  ) +
  scale_color_manual(values= sex_colors)+
  scale_shape_manual(values = sex_shapes)

(g.tes + g.zi) + plot_layout(widths = c(3,1))+
  plot_annotation(tag_levels = 'I', tag_prefix = "", tag_suffix = ")")

####Pathogen Load####
source("r_scripts/dataCleaning_pathogenLoad.R")
#make quantity at dpi 0 = 0 as baseline for models
ti$quantity[ti$dpi == 0] <- 0

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

#Add small constant for transformation
ti.q$quantity1 <- ti.q$quantity + 1
hist(ti.q$quantity1)


dat.q <- ti.q %>%
  filter(treatment == "Inoculated" & dpi > -1) %>%
  mutate(
    dpi  = as.numeric(dpi),
    temp = factor(temp),
    band_number = factor(band_number),
    ever_infected = factor(ever_infected),
    eve_diseased = factor(ever_diseased)
  )


#Is gamma distributed glm better?
dat.q <- ti.q %>%
  filter(treatment == "Inoculated" & dpi > -1) %>%
  mutate(
    dpi  = as.numeric(dpi),
    temp = factor(temp),
    band_number = factor(band_number),
    ever_infected = factor(ever_infected),
    eve_diseased = factor(ever_diseased)
  )

#Quantity at dpi 0 was 0 



dat.q$sex <- as.factor(dat.q$sex)
dat.q <- dat.q %>%
  mutate(sex = dplyr::recode(sex,
                             "M" = "Male",
                             "F" = "Female"))

dat.q <- dat.q %>%
  mutate(
    temp  = relevel(temp, ref = "Warm"),
    sex = relevel(sex, ref = "Female") 
  )

dat.q$dpi.f <- as.factor(dat.q$dpi)
dat.q$log10_quantity1 <- log10(dat.q$quantity+1)
hist(dat.q$log10_quantity1)


#Zero-inflated model
#model selection
dat.q.mod <- dat.q %>%
  filter(dpi.f != 0)

q1 <- glmmTMB(log10_quantity1 ~ temp*dpi.f + (1|band_number),
              data=dat.q.mod, 
              ziformula = ~ temp,
              family = ziGamma(link = "log"))

q2 <- glmmTMB(log10_quantity1 ~ temp*dpi.f + (1|band_number),
              data=dat.q.mod, 
              ziformula = ~ sex,
              family = ziGamma(link = "log"))

q3 <- glmmTMB(log10_quantity1 ~ temp*dpi.f + sex + (1|band_number),
              data=dat.q.mod, 
              ziformula = ~ temp,
              family = ziGamma(link = "log"))

q4 <- glmmTMB(log10_quantity1 ~ temp*dpi.f * sex + (1|band_number),
              data=dat.q.mod, 
              ziformula = ~ temp,
              family = ziGamma(link = "log"))

q5 <- glmmTMB(log10_quantity1 ~ dpi.f + (1|band_number),
              data=dat.q.mod, 
              ziformula = ~ sex * temp,
              family = ziGamma(link = "log"))

q6 <- glmmTMB(log10_quantity1 ~ dpi.f + temp + (1|band_number),
              data=dat.q.mod, 
              ziformula = ~ sex * temp,
              family = ziGamma(link = "log"))

q7 <- glmmTMB(log10_quantity1 ~ dpi.f + (1|band_number),
              data=dat.q.mod, 
              ziformula = ~ temp,
              family = ziGamma(link = "log"))

q8 <- glmmTMB(log10_quantity1 ~ dpi.f + (1|band_number),
              data=dat.q.mod, 
              ziformula = ~ sex,
              family = ziGamma(link = "log"))

q9 <- glmmTMB(log10_quantity1 ~ dpi.f * sex + (1|band_number),
              data=dat.q.mod, 
              ziformula = ~ sex * temp,
              family = ziGamma(link = "log"))

aictab(cand.set=list(q1, q2, q3, q4, q5, q6, q7, q8, q9), 
       modnames=c("1", "2", "3", "4", "5", "6", "7", "8", "9"))

BIC(q5, q6)

glm.zi.pl <- glmmTMB(log10_quantity1 ~ dpi.f + (1|band_number),
                     data=dat.q.mod, 
                     ziformula = ~ sex * temp,
                     family = ziGamma(link = "log"))

simulateResiduals(glm.zi.pl, plot=T)
hist(resid(glm.zi.pl))
summary(glm.zi.pl)
car::Anova(glm.zi.pl, type = "III")

emm <- emmeans(
  glm.zi.pl,
  ~ temp | dpi.f,
  at = list(dpi = sort(unique(dat.q.mod$dpi))),
  component = "response", #averaging in zero component
  #component = "cond", #only infected individuals
  type = "response"
)

emm_df <- as.data.frame(emm)

#merge back dpi 0, 3, and 35 for emmean visualization
emm_df <- emm_df %>%
  mutate(
    dpi.f = factor(dpi.f, levels = levels(dat.q.mod$dpi.f)),
    temp  = factor(temp,  levels = levels(dat.q.mod$temp))
  )

ti.q <- ti.q %>%
  filter(dpi != "-12")

ti.q$dpi.f <- as.factor(ti.q$dpi)
ti.q$temp <- as.factor(ti.q$temp)

full_grid <- expand_grid(
  temp  = levels(ti.q$temp),
  dpi.f = levels(ti.q$dpi.f)
)

emm_df_full <- full_grid %>%
  left_join(emm_df, by = c("temp", "dpi.f"))

dodge = position_dodge(width = 0.2)

g.pl <- ggplot(dat.q, aes(x = dpi.f, y = log10_quantity1), color = temp) +
   #Jittered
  geom_jitter(
    alpha  = 0.5,
    size   = 2, 
    position=position_jitterdodge(jitter.width = 0.1, jitter.height = 0.1, dodge.width = 0.2),
    aes(shape = sex, group = temp, color=temp)
  ) +
  # 95% CI around means
  geom_errorbar(
    data = emm_df_full,
    aes(
      x    = dpi.f,
      ymin = asymp.LCL,
      ymax = asymp.UCL,
      color = temp,
      group = temp
    ),
    position = dodge,
    alpha        = 1,
    inherit.aes  = FALSE,
    na.rm        = TRUE,
    width=0.1
  ) +
  #model means (component = "response")
  geom_point(
    data = emm_df_full,
    aes(x = dpi.f, y = emmean, color = temp, group = temp),
    size = 3,
    position=dodge,
    na.rm = TRUE
  ) +
  geom_point(
    data = emm_df_full,
    aes(x = dpi.f, y = emmean, group = temp),
    size = 3.2,
    color="black",
    shape=1,
    position=dodge,
    na.rm = TRUE
  ) +
  geom_line(
    data = emm_df_full,
    aes(x = dpi.f, y = emmean, color = temp, group =temp),
    linewidth = 0.5,
    position=dodge,
    na.rm = TRUE
  ) +
  labs(
    x     = "Days Post Inoculation",
    y     = "Log10(Pathogen Load +1)",
    color = "Temperature",
    fill  = "Temperature",
    shape = "Sex",
    title = "Gamma Model"
  ) +
  scale_color_manual(values = temp_colors) +
  scale_fill_manual(values  = temp_colors)+
  scale_shape_manual(values = sex_shapes)
g.pl

# Zero-inflation component: probability of extra zeros
emm_zi <- emmeans(
  glm.zi.pl,
  ~ sex * temp,
  component = "zi",    
  type      = "response"
)

emm_zi_df <- as.data.frame(emm_zi)

# emm_zi_df will have columns like:
# sex, temp, prob, SE, asymp.LCL, asymp.UCL, etc.

g.zi.pl<- ggplot(emm_zi_df,
              aes(x = temp, y = response, color = sex, group = sex, shape = sex)) +
  geom_point(size = 3) +
  geom_line() +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),
                width = 0.03) +
  labs(
    x = "Temperature",
    y = "Probability of No Pathogen Load [Pr(structural zero)]",
    color = "Sex",
    shape = "Sex",
    title = "Zero-inflation Probability"
  ) +
  scale_color_manual(values= sex_colors)+
  scale_shape_manual(values = sex_shapes)

(g.pl + g.zi.pl) + plot_layout(widths = c(3,1))+
  plot_annotation(tag_levels = 'I', tag_prefix = "", tag_suffix = ")")


####Tolerance####
ti.t.max <- ti %>%
  dplyr::select(band_number, dpi, total_eye_score, treatment, temp, mass, sex, quantity, log10_quantity, groups, ever_diseased, ever_infected) %>%
  filter(treatment == "Inoculated")%>%
  group_by(band_number) %>%
  summarise(
    max_tes = if (all(is.na(total_eye_score))) NA_real_ else max(total_eye_score, na.rm = TRUE),
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

ti.wx <- ti.wx %>%
  mutate(sex = dplyr::recode(sex,
                             "F" = "Female",
                             "M" = "Male"))

ti.wx <- ti.wx %>%
  mutate(temp = fct_relevel(temp,
                            "Warm", "Cold"))

#identify seropositive birds from quaran

#tissue-specific tolerance - only infected
ggplot(ti.wx, aes(x=sex, y=tolerance, color=temp))+
  geom_boxplot(width=0.25, aes(fill=temp), alpha=0.1, color="black", position=position_dodge(width=0.5))+
  geom_jitter(size=3, alpha=0.75, 
              position=position_jitterdodge(dodge.width=0.5, jitter.width = 0.25, jitter.height = 0))+
  scale_color_manual(values=temp_colors)+
  scale_fill_manual(values=temp_colors)+
  labs(x="Sex", y="Tolerance: \n -(Max Pathology / Max Pathogen Load)", color="Temperature", fill="Temperature")+
  theme(strip.text = element_text(size = 14))

#Test: Does temperature affect tolerance in birds that got infected?
wilcox.test(tolerance ~ temp, data=ti.wx)
wilcox.test(tolerance ~ sex, data=ti.wx %>% filter(temp == "Warm"))
wilcox.test(tolerance ~ sex, data=ti.wx %>% filter(temp == "Cold"))

wilcox.test(tolerance ~ temp, data=ti.t.max %>% filter(sex == "M"))
wilcox.test(tolerance ~ temp, data=ti.t.max %>% filter(sex == "F"))

lmt <- lm(tolerance ~ temp * sex, data= ti.t.max)
simulateResiduals(lmt, plot = T)
summary(lmt)

# data:  tolerance by temp
# W = 84, p-value = 0.7545

lmg <- lm(tolerance ~ temp * sex, data=ti.wx)

simulateResiduals(lmg, plot=T)
summary(lmg)

#Brown-Forsythe; is tolerance equally variable across temperatures?
bf_tol <- leveneTest(tolerance ~ temp, data = ti.wx, center = median)

# Levene's Test for Homogeneity of Variance (center = median)
#       Df F value Pr(>F)
# group  1  0.1664 0.6871
#       23              


####Fever Score####
source("r_scripts/dataCleaning_fever.R")
#fever_change = score - baseline
#fever_diff = change from previous score

ti.f$dpi <- as.factor(ti.f$dpi)

unique(ti.f$dpi)

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
##Control Only##
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

#So it isn't a sex effect without stimulation
summary(c1)
summary(c5)
simulateResiduals(c5, plot=T)
plot(allEffects(c5))

#What were the ranges of fever scores between rooms?
frange <- ti.f %>%
  group_by(temp, treatment) %>%
  summarise(tmax = max(fever_score),
            tmin = min(fever_score))
frange

#Model Selection
#because sampled_first is a function of dpi (it is fixed for each dpi), it is collinear with dpi.f, so do not include.
#We treat dpi as a factor (dpi.f) because this is a linear model therefore we must look only at linear differences
#between dpis (no smoothing)

#DPI 18 maybe is an issue - uncomment and re-run to see exclusion analysis
# ti.f <- ti.f %>%
#   filter(dpi.f != 18)

##All Birds##
#I want to know whether treatment or temperature affect fever_score or shape of fever response
#while accounting for sex
m0 <- glmmTMB(fever_score ~ temp + treatment + sex + dpi.f + (1|band_number), data = ti.f)
m0.5 <- glmmTMB(fever_score ~ treatment + temp * sex + dpi.f + (1|band_number), data = ti.f)
m1 <- glmmTMB(fever_score ~ temp * treatment + sex + dpi.f + (1|band_number), data = ti.f)
m2 <- glmmTMB(fever_score ~ temp * dpi.f + treatment + sex + (1|band_number), data = ti.f)
m3 <- glmmTMB(fever_score ~ treatment * dpi.f + temp + sex + (1|band_number), data = ti.f)
m4 <- glmmTMB(fever_score ~ temp * treatment * dpi.f + sex + (1|band_number), data = ti.f)
m5 <- glmmTMB(fever_score ~ temp * treatment * dpi.f * sex + (1|band_number), data = ti.f)
m6 <- glmmTMB(fever_score ~ (treatment + temp) * dpi.f + sex + (1|band_number), data = ti.f)
m6s <- glmmTMB(fever_score ~ (treatment + temp) * dpi.f * sex + (1|band_number), data = ti.f)
m7 <- glmmTMB(fever_score ~ (treatment + temp + sex) * dpi.f + (1|band_number), data = ti.f)
m_null <- glmmTMB(fever_score ~ 1 + (1|band_number), data = ti.f)

aictab(list(m0, m0.5, m1, m2, m3, m4, m5, m6, m6s, m7, m_null),
       modnames = c("m0","m0.5","m1","m2","m3","m4","m5","m6","m6s", "m7", "null"))

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
plot(allEffects(glm_fever))

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

#Is there an effect of which room was sampled first?
ti.f$sampled_first <- as.factor(ti.f$sampled_first)
glm_fever.cont.sf <- glmmTMB(fever_score ~ temp + sex + sampled_first + (1|band_number), 
                          data=ti.f %>% filter(treatment == "Control"))

simulateResiduals(glm_fever.cont.sf, plot=T)
summary(glm_fever.cont.sf)
plot(allEffects(glm_fever.cont.sf))

glm_fever.cont <- glmmTMB(fever_score ~ temp + sex + dpi.f + (1|band_number), 
                          data=ti.f %>% filter(treatment == "Control"))

#Sampling order is nowhere near as important as dpi. Since they are collinear, they cannot both
  #be included in the model. Therefore, use dpi.
AIC(glm_fever.cont.sf, glm_fever.cont)

simulateResiduals(glm_fever.cont, plot=T)
summary(glm_fever.cont)
car::Anova(glm_fever.cont)
plot(allEffects(glm_fever.cont))

# Get emmeans
emm_df <-  emmeans(glm_fever.cont, ~  sampled_first * temp, type = "response") %>%
  as.data.frame() %>%
  left_join(
    ti.f %>%
      distinct(temp, sampled_first, groups),
    by = c("temp", "sampled_first")
  )

# emm_df <- emm_df %>%
#   mutate(sex = dplyr::recode(sex,
#                              "M" = "Male",
#                              "F" = "Female"))

ggplot() +
  geom_point(data = ti.f.g %>% filter(treatment == "Control"),
             aes(x = temp, y = fever_score, color = sampled_first),
             position =dodge,
             alpha = 0.5, size = 2) +
  
  geom_errorbar(data = emm_df,
                aes(x = temp, ymin = lower.CL, ymax = upper.CL, groups = sampled_first),
                color="black", position = dodge, width = 0.1) +
  
  # geom_ribbon(data = emm_df,
  #               aes(x = as.numeric(dpi.f), ymin = lower.CL, ymax = upper.CL, groups = groups, fill=groups),
  #               position = dodge, alpha=0.2) +
  # 
  geom_line(data = emm_df,
            aes(x = temp, y = emmean, linetype = temp, groups = sampled_first),
            color="black", position = dodge, size =0.25) +
  
  geom_point(data = emm_df,
             aes(x = temp, y = emmean, color = sampled_first),
             position = dodge, size = 2.8) +
  geom_point(data = emm_df,
             aes(x = temp, y = emmean, group=sampled_first),
             position = dodge, size = 2.8,
             shape=1, color="black") +
  # stat_summary(data = ti.f.g, aes(x=dpi.f, y= fever_score, group = interaction(treatment, temp), color=groups), geom="line", fun = "mean")+
  
  #scale_color_manual(values = sex_colors) +
  #scale_fill_manual(values=treat_colors)+
  #scale_shape_manual(values = c(1, 1, 1, 1))+
  labs(x = "Room", y = "Ocular Temperature (C)", color = "Sampled First", shape="Sampled First") +
  # scale_linetype_manual(
  #   name   = "Inoculation Type",
  #   values = c("dashed", "solid"),
  #   labels = c("Control", "Inoculated")
  # )+
  #facet_grid(~t)+ #, scales = "free_y"
  theme(strip.text = element_text(size=12))

# Get emmeans
emm_df <-  emmeans(glm_fever.cont, ~ sex * dpi.f * temp, type = "response") %>%
  as.data.frame() %>%
  left_join(
    ti.f %>%
      distinct(temp, dpi.f, sex, groups),
    by = c("temp", "dpi.f", "sex")
  )

emm_df <- emm_df %>%
  mutate(sex = dplyr::recode(sex,
                             "M" = "Male",
                             "F" = "Female"))

ggplot() +
  geom_point(data = ti.f.g %>% filter(treatment == "Control"),
             aes(x = dpi.f, y = fever_score, color = sex),
             position =dodge,
             alpha = 0.5, size = 2) +
  
  geom_errorbar(data = emm_df,
                aes(x = dpi.f, ymin = lower.CL, ymax = upper.CL, groups = sex),
                color="black", position = dodge, width = 0.1) +
  
  # geom_ribbon(data = emm_df,
  #               aes(x = as.numeric(dpi.f), ymin = lower.CL, ymax = upper.CL, groups = groups, fill=groups),
  #               position = dodge, alpha=0.2) +
  # 
  geom_line(data = emm_df,
            aes(x = as.numeric(dpi.f), y = emmean, linetype = temp, groups = sex),
            color="black", position = dodge, size =0.25) +
  
  geom_point(data = emm_df,
             aes(x = dpi.f, y = emmean, color = sex),
             position = dodge, size = 2.8) +
  geom_point(data = emm_df,
             aes(x = dpi.f, y = emmean, group=sex),
             position = dodge, size = 2.8,
             shape=1, color="black") +
  # stat_summary(data = ti.f.g, aes(x=dpi.f, y= fever_score, group = interaction(treatment, temp), color=groups), geom="line", fun = "mean")+
  
  scale_color_manual(values = sex_colors) +
  #scale_fill_manual(values=treat_colors)+
  #scale_shape_manual(values = c(1, 1, 1, 1))+
  labs(x = "Days Post Inoculation", y = "Ocular Temperature (C)", color = "Treatment Group", shape="Treatment Group", linetype="Inoculation Type") +
  # scale_linetype_manual(
  #   name   = "Inoculation Type",
  #   values = c("dashed", "solid"),
  #   labels = c("Control", "Inoculated")
  # )+
  facet_grid(~temp)+ #, scales = "free_y"
  theme(strip.text = element_text(size=12))


#Model predictions by temp comparing sex directly in inoculated individuals
#We know there is a trend for males to be more likely to develop pathology when infected in warm temperatures
#Therefore, we would expect higher ocular temps at high temperatures, so we need a sex interaction
#But subset to only inoculated individuals

##Inoculated Only##
ti.f.i <- ti.f %>% 
  filter(treatment == "Inoculated")

#Model selection
m0i <- glmmTMB(fever_score ~ temp + sex + dpi.f + (1|band_number), data = ti.f.i)
m0.5i <- glmmTMB(fever_score ~ temp + sex * dpi.f + (1|band_number), data = ti.f.i)
m1i <- glmmTMB(fever_score ~ temp * sex + dpi.f + (1|band_number), data = ti.f.i)
m2i <- glmmTMB(fever_score ~ sex * dpi.f + temp + (1|band_number), data = ti.f.i)
m3i <- glmmTMB(fever_score ~ temp * dpi.f + sex + (1|band_number), data = ti.f.i)
m4i <- glmmTMB(fever_score ~ dpi.f * temp * sex + (1|band_number), data = ti.f.i)
m_nulli <- glmmTMB(fever_score ~ 1 + (1|band_number), data = ti.f.i)

aictab(list(m0i, m0.5i, m1i, m2i, m3i, m4i, m_nulli),
       modnames = c("m0", "m0.5i", "m1","m2","m3","m4", "null"))

simulateResiduals(m1i, plot=T)

glm_fever.i <- glmmTMB(fever_score ~ temp * sex + dpi.f + (1|band_number), data = ti.f.i)
simulateResiduals(glm_fever.i, plot=T)
summary(glm_fever.i)
car::Anova(glm_fever.i, type = "III")
plot(allEffects(glm_fever.i))

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
  geom_point(data = ti.f.i %>% filter(ever_diseased == 1),
             aes(x = dpi.f, y = fever_score, color = sex, shape = sex),
             position =dodge,
             alpha = 0.5, size = 2) +
  geom_line(data=ti.f.i %>% filter(ever_diseased == 1),
            aes(x=dpi.f, y=fever_score, color=sex, group = as.factor(band_number)),
            position = dodge, alpha=0.1)+
    
  geom_errorbar(data = emm_df,
                aes(x = dpi.f, ymin = lower.CL, ymax = upper.CL, group = sex),
                color="black", position = dodge, width = 0.1) +
  geom_line(data = emm_df,
            aes(x = as.numeric(dpi.f), y = emmean, linetype = sex),
            color="black", position = dodge, size =0.25) +
  
  geom_point(data = emm_df,
             aes(x = dpi.f, y = emmean, color = sex, groups = sex, shape = sex),
             position = dodge, size = 2.8) +
  scale_color_manual(values = sex_colors) +
  scale_shape_manual(values = c(16, 17))+
  labs(x = "Days Post Inoculation", y = "Ocular Temperature (C)", color = "Sex", shape="Sex", linetype="Sex") +
  facet_grid(~temp)+ #, scales = "free_y"
  theme(strip.text = element_text(size=12))


#Males have higher ocular temperature in warm rooms - is this because more develop pathology?
#Males more likely to be diseased (almost), if we account for this, does effect go away?
ti.f.i$ever_diseased <- as.factor(ti.f.i$ever_diseased)
ti.f.i$ever_infected <- as.factor(ti.f.i$ever_infected)

#Model selection
m1id <- glmmTMB(fever_score ~ temp * sex + dpi.f + ever_diseased + (1|band_number), data = ti.f.i)
m2id <- glmmTMB(fever_score ~ sex * ever_diseased * dpi.f + temp + (1|band_number), data = ti.f.i)
m3id <- glmmTMB(fever_score ~ sex * ever_diseased + dpi.f + temp + (1|band_number), data = ti.f.i)
m4id <- glmmTMB(fever_score ~ sex + ever_diseased * dpi.f + temp + (1|band_number), data = ti.f.i)
m5id <- glmmTMB(fever_score ~ sex + ever_diseased * dpi.f * temp + (1|band_number), data = ti.f.i)
m6id <- glmmTMB(fever_score ~ dpi.f * temp * sex * ever_diseased + (1|band_number), data = ti.f.i)
m7id <- glmmTMB(fever_score ~ ever_diseased * dpi.f + temp * sex + (1|band_number), data = ti.f.i)
m8id <- glmmTMB(fever_score ~ (sex + ever_diseased) * dpi.f + temp + (1|band_number), data = ti.f.i)
m_nulli <- glmmTMB(fever_score ~ 1 + (1|band_number), data = ti.f.i)

aictab(list(m1id, m2id, m3id, m4id, m5id, m6id, m7id, m8id, m_nulli),
       modnames = c("m1","m2","m3","m4", "m5", "m6", "m7", "m8", "null"))

simulateResiduals(m4id, plot=T)
summary(m4id)
car::Anova(m4id, type="III")
summary(m7id)
car::Anova(m7id, type ="III")

#compare to treatment
glm_fd <- glmmTMB(fever_score ~ sex + ever_diseased * dpi.f + temp + (1|band_number), data = ti.f)
glm_ft <- glmmTMB(fever_score ~ sex + treatment * dpi.f + temp + (1|band_number), data = ti.f)
glm_fi <- glmmTMB(fever_score ~ sex + ever_infected * dpi.f + temp + (1|band_number), data = ti.f)


aictab(list(glm_fd, glm_ft, glm_fi),
       modnames = c("glm_fd","glm_ft","glm_fi"))

glm_fdi <- glmmTMB(fever_score ~ sex + ever_diseased * dpi.f + temp + (1|band_number), data = ti.f.i)
glm_fii <- glmmTMB(fever_score ~ sex + ever_infected * dpi.f + temp + (1|band_number), data = ti.f.i)

aictab(list(glm_fdi, glm_fii),
       modnames = c("glm_fdi","glm_fii"))

#Final model
glm_fever.id <- glmmTMB(fever_score ~ sex + ever_diseased * dpi.f + temp + (1|band_number),
                        data = ti.f.i)

simulateResiduals(glm_fever.id, plot=T)
summary(glm_fever.id)
car::Anova(glm_fever.id, type = "III")

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
             alpha = 0.5, size = 2) +
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
      "Female.Peak"     = "Female Peak",
      "Male.Baseline" = "Male Baseline",
      "Male.Peak"     = "Male Peak"
    )
  )

emm_df <- emm_df %>%
  mutate(
    sex_phase = interaction(sex, fever_type),
    sex_phase = recode_factor(
      sex_phase,
      "Female.Baseline" = "Female Baseline",
      "Female.Peak"     = "Female Peak",
      "Male.Baseline" = "Male Baseline",
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
  
  scale_shape_manual(values = c(1,  16, 2, 17)) +
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

ti.long$ever_diseased <- as.factor(ti.long$ever_diseased)
ti.long$ever_infected <- as.factor(ti.long$ever_infected)

ti.long$sex <- as.factor(ti.long$sex)

#What about disease?
d1 <- glmmTMB(fever_value ~ fever_type * ever_diseased + temp * sex + (1|band_number), data=ti.long)
d2 <- glmmTMB(fever_value ~ (fever_type + temp + sex) * ever_diseased + (1|band_number), data=ti.long)
d3 <- glmmTMB(fever_value ~ fever_type * temp * sex * ever_diseased + (1|band_number), data=ti.long)
d4 <- glmmTMB(fever_value ~ fever_type + temp + sex + ever_diseased + (1|band_number), data=ti.long)
d5 <- glmmTMB(fever_value ~ fever_type * ever_diseased + sex + temp + (1|band_number), data=ti.long)
null <- glmmTMB(fever_value ~ 1 + (1|band_number), data=ti.long)

aictab(cand.set=list(d1, d2, d3, d4, d5, null), 
       modnames=c("d1", "d2", "d3", "d4", "d5", "null"))

simulateResiduals(d1, plot=T)
summary(d1)
plot(allEffects(d1))
car::Anova(d1, type="III")

lm2d <- glmmTMB(fever_value ~ fever_type * ever_diseased + sex * temp + (1|band_number), data=ti.long)
simulateResiduals(lm2d, plot=T)
summary(lm2d)
car::Anova(lm2d, type = "III")

lm2dis <- glmmTMB(fever_value ~ fever_type * ever_diseased + sex * temp + (1|band_number), data=ti.long)
lm2t <- glmmTMB(fever_value ~ fever_type * treatment + sex * temp + (1|band_number), data=ti.long)
lm2i <- glmmTMB(fever_value ~ fever_type * ever_infected + sex * temp + (1|band_number), data=ti.long)

#Ever diseased is a better predictor of fever_value than treatment (change lm2d data to ti.long)
aictab(cand.set=list(lm2dis, lm2t, lm2i), 
       modnames=c("lm2dis", "lm2t", "lm2i"))

#sex significant if you include controls
simulateResiduals(lm2dis, plot=T)
summary(lm2dis)
car::Anova(lm2dis, type="III")

#emmeans including ever_disaesed
#estimate fever_type means within each combination of temperature x sex x ever_diseased
emm <-emmeans(lm2dis, ~fever_type | temp*sex*ever_diseased)
tests <- contrast(emm, method = "pairwise", adjust = "tukey")


emm_df <- as.data.frame(emm)  # has emmean, SE, lower.CL, upper.CL

ti.long <- ti.long %>%
  mutate(
    sex_phase = interaction(sex, fever_type),
    sex_phase = recode_factor(
      sex_phase,
      "Female.Baseline" = "Female Baseline",
      "Female.Peak"     = "Female Peak",
      "Male.Baseline" = "Male Baseline",
      "Male.Peak"     = "Male Peak"
    )
  )

emm_df <- emm_df %>%
  mutate(
    sex_phase = interaction(sex, fever_type),
    sex_phase = recode_factor(
      sex_phase,
      "Female.Baseline" = "Female Baseline",
      "Female.Peak"     = "Female Peak",
      "Male.Baseline" = "Male Baseline",
      "Male.Peak"     = "Male Peak"
    )
  )

#Plot Emmeans sex comparison
ggplot(ti.long, aes(x = temp,
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
  
  scale_shape_manual(values = c(1,  16,2, 17)) +
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
ti.fci <- ti.fc #%>%
  filter(treatment == "Inoculated")
ti.fci$ever_diseased <- as.factor(ti.fci$ever_diseased)

ti.fci <- ti.fci %>%
  dplyr::mutate(ever_diseased = dplyr::recode(ever_diseased,
                                              "0" = "Never Diseased",
                                              "1" = "Ever Diseased"),
                sex = dplyr::recode(sex,
                "F" = "Female",
                "M" = "Male"))

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

summary(fi1.2)

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
              position=position_jitterdodge(dodge.width=0.4, jitter.width = 0.2, jitter.height = 0))+
  geom_point(data=emm_df, aes(y=emmean, x=ever_diseased), color="black", size=3)+
  geom_errorbar(data=emm_df, aes(y=emmean, ymin=lower.CL, ymax=upper.CL, x=ever_diseased), width=0., color="black", size=.75)+
  labs(x="Disease Category", y= "Magnitude of Fever Change\n(Peak - Baseline [C])", color="Temperature", shape = "Sex")+
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

####Mass####
#Sample sizes


ti.m <- ti %>%
  filter(dpi %in% c(-28, -12, 3, 14, 21, 28, 35),
         !is.na(mass))

ti.m %>%
  dplyr::select(treatment, temp, groups, dpi)%>%
  tbl_summary(
    by=dpi
  )%>%
  modify_header(
    label ~ "**All Birds**"
  )

ggplot(data=ti.m, aes(x=dpi, y=mass, color=temp))+
  geom_point()+
  geom_line(aes(group = as.factor(band_number)), size=0.1)


ggplot(ti.m, aes(x=dpi, y=mass, color=temp))+
  geom_jitter()+
  facet_wrap(~sampled_first)

hist(ti.m$mass)

ti.m$dpi.f <- as.factor(ti.m$dpi)
ti.m$temp <- relevel(as.factor(ti.m$temp), ref = "Warm")
ti.m$treatment <- relevel(as.factor(ti.m$treatment), ref = "Control")

#did mass differ at baseline?
t.test(mass ~ temp, data=ti.m %>% filter(dpi.f == -28))

t.test(mass ~ temp, data=ti.m %>% filter(dpi.f == -12))

#Model Selection
#because sampled_first is a function of dpi (it is fixed for each dpi), it is collinear with dpi.f, so do not include.
lm.a <- glmmTMB(mass ~ temp * treatment + sex + dpi.f + (1|band_number), data=ti.m)
lm.b <- glmmTMB(mass ~ temp * treatment +  dpi.f +(1|band_number), data=ti.m)


lm.c <- glmmTMB(mass ~ temp + treatment + sex + dpi.f +(1|band_number), data=ti.m)
lm.d <- glmmTMB(mass ~ temp + treatment +  dpi.f +(1|band_number), data=ti.m)

lm.e <- glmmTMB(mass ~ temp + treatment *  dpi.f + sex +(1|band_number), data=ti.m)
lm.f <- glmmTMB(mass ~ treatment + temp *  dpi.f +(1|band_number), data=ti.m)

lm.g<- glmmTMB(mass ~ temp + sex + dpi.f +(1|band_number), data=ti.m)
lm.h <- glmmTMB(mass ~ treatment + sex + dpi.f +(1|band_number), data=ti.m)

lm.i <- glmmTMB(mass ~ temp * treatment + dpi.f +(1|band_number), data=ti.m)
lm.j <- glmmTMB(mass ~ temp * treatment + sex + dpi.f +(1|band_number), data=ti.m)
lm.k <- glmmTMB(mass ~ 1 + (1|band_number), data=ti.m)

aictab(cand.set=list(lm.a, lm.b, lm.c,  lm.d, lm.e, lm.f, lm.g, lm.h, lm.i, lm.j, lm.k), 
       modnames=c("lm.a",  "lm.b", "lm.c", "lm.d", "lm.e", "lm.f", "lm.g", "lm.h", "lm.i", "lm.j", "lm.k"))

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
car::Anova(lm.dd, type = "III")

lm7 <- glmmTMB(mass ~ treatment + temp * dpi.f +(1|band_number), data=ti.m)

DHARMa::simulateResiduals(lm7, plot=T)

summary(lm7)
plot(allEffects(lm7))

hist(resid(lm7))

car::Anova(lm7, type = "III")
car::Anova(lm7, type = "II")

#What if we compare only baseline after acclimization to temperatures (dpi -12 instead of -28)
lm7.no.bl <- glmmTMB(mass ~ treatment + temp * dpi.f +(1|band_number), data=ti.m %>% filter(dpi != -28))

lm7.dis <- glmmTMB(mass ~ ever_diseased + temp * dpi.f +(1|band_number), data=ti.m %>% filter(dpi != -28))

summary(lm7.dis)

summary(lm7.no.bl)
simulateResiduals(lm7.no.bl, plot=T)
car::Anova(lm7.no.bl, type="III")

emm_df <- emmeans(lm7, ~ temp * treatment * dpi.f, re.form = NA) %>%
  as.data.frame()

emm_df <- emmeans(lm7.no.bl, ~ temp * treatment * dpi.f, re.form = NA) %>%
  as.data.frame() %>%
  left_join(
    ti.m %>% filter(dpi != -28)%>%
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
  geom_point(data = ti.m.g %>% filter(dpi != -28),
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
lm1<-glmmTMB(elisa_od~treatment * dpi.f + temp + (1|band_number),
             data=ti.ab, 
             family=Gamma(link = "log"))

simulateResiduals(lm1, plot=T)
summary(lm1)

lm1s<-glmmTMB(elisa_od~treatment * dpi.f + temp + sex + (1|band_number),
              data=ti.ab, 
              family=Gamma(link = "log"))

simulateResiduals(lm1s, plot=T)

#sex not significant
AIC(lm1, lm1s)
summary(lm1s)

Anova(lm1, type = "III")

emm <- emmeans(lm1, ~ treatment | temp*dpi.f, type = "response")

pairs(emm, by=c("temp","dpi.f"), adjust = "tukey")

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
ggplot(emm_df, aes(x = groups, y = response, color = groups, shape=dpi.f)) +
  geom_jitter(
    data = ti.ab,
    aes(y = elisa_od, color = groups, shape = dpi.f),
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


####Phagocytosis####
source("/Users/jesse/Documents/GitHub/Temperature-and-Immunity-2022/Analysis/r_scripts/dataCleaning_phago.R")

##   Binomial model and the
#    weights represent the total number of trials, since the response is a
#    pre-calculated proportion. The interaction between temp and treatment was
#    not significant.

ti.p$temp <- as.factor(ti.p$temp)
ti.p$treatment <- as.factor(ti.p$treatment)

ti.p <- ti.p %>%
  mutate(
    temp      = relevel(temp, ref = "Warm"),        # or whichever baseline you want
    treatment = relevel(treatment, ref = "Control") # change as needed
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

ggplot(emm_df, aes(x = groups, y = prob, color=groups))+
  geom_jitter(data = ti.p, aes(x = groups, y = phago_score, color = groups),
              width = 0.2, alpha = 0.75, size=2.5)+
  geom_point(size=3, color="black")+
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL, y=prob), width=0., color="black", size=0.75)+
  labs(x="Treatment", y="Phagocytosis Score", color="Treatment")+
  scale_color_manual(values=treat_colors)+
  theme(
    axis.text.x = element_text(size=13, angle=45, hjust=1)
  )+
facet_wrap(~dis_9)

ggplot(ti.p, aes(x=groups, y=phago_score, color=groups))+
  geom_jitter()+
  facet_wrap(~dis_9)+
  scale_color_manual(values=treat_colors)+
  theme(
    axis.text.x = element_text(size=13, angle=45, hjust=1)
  )

glm4.d <- glmmTMB(phago_score~temp+ dis_9 + (1|band_number), 
                weights=wbc_total+phago_total, 
                data=ti.p %>% filter(treatment == "Inoculated"), family="binomial")

glm4.i <- glmmTMB(phago_score~temp+ inf_9 + (1|band_number), 
                  weights=wbc_total+phago_total, 
                  data=ti.p%>% filter(treatment == "Inoculated"), family="binomial")

AIC(glm4, glm4.d, glm4.i)

simulateResiduals(glm4.i, plot=T)
summary(glm4.i)
car::Anova(glm4.d, type="III")

emm <- emmeans(glm4.i, ~ temp * inf_9, type = "response", re.form = NA)
emm_df <- as.data.frame(emm)  

i<- ggplot(emm_df, aes(x = temp, y = prob, color=temp))+
  geom_jitter(data = ti.p%>% filter(treatment == "Inoculated"), aes(x = temp, y = phago_score, color = groups, shape= sex),
              width = 0.2, alpha = 0.75, size=2.5)+
  geom_point(size=3, color="black")+
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL, y=prob), width=0., color="black", size=0.75)+
  labs(x="Temperature", y="Phagocytosis Score", color="Temperature")+
  scale_color_manual(values=treat_colors)+
  labs(title = "Infected DPI 9")+
  theme(
    axis.text.x = element_text(size=13, angle=45, hjust=1)
  )+
  facet_wrap(~inf_9, labeller = labeller(inf_9 = c('0'="Uninfected", '1' = "Infected")))

emm <- emmeans(glm4.d, ~ temp * dis_9, type = "response", re.form = NA)
emm_df <- as.data.frame(emm)  

d<- ggplot(emm_df, aes(x = temp, y = prob, color=temp))+
  geom_jitter(data = ti.p%>% filter(treatment == "Inoculated"), aes(x = temp, y = phago_score, color = groups, shape = sex),
              width = 0.2, alpha = 0.75, size=2.5)+
  geom_point(size=3, color="black")+
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL, y=prob), width=0., color="black", size=0.75)+
  labs(x="Temperature", y="Phagocytosis Score", color="Temperature")+
  scale_color_manual(values=treat_colors)+
  theme(
    axis.text.x = element_text(size=13, angle=45, hjust=1)
  )+
  labs(title = "Diseased DPI 9")+
  facet_wrap(~dis_9, labeller = labeller(dis_9 = c('0'="No Pathology", '1' = "Pathology")))


emm <- emmeans(glm4, ~ temp * treatment, type = "response", re.form = NA)
emm_df <- as.data.frame(emm) 

t <-ggplot(emm_df, aes(x = temp, y = prob, color=temp))+
  geom_jitter(data = ti.p, aes(x = temp, y = phago_score, color = groups, shape = sex),
              width = 0.2, alpha = 0.75, size=2.5)+
  geom_point(size=3, color="black")+
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL, y=prob), width=0., color="black", size=0.75)+
  labs(x="Temperature", y="Phagocytosis Score", color="Temperature")+
  scale_color_manual(values=treat_colors)+
  labs(title = "Treatment")+
  theme(
    axis.text.x = element_text(size=13, angle=45, hjust=1)
  )+
  facet_wrap(~treatment)

i + theme(legend.position = "none")+
d + theme(legend.position = "none")+
t
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
