#T+I 2022 MG IgY Antibody Analysis
rm(list=ls())
# setwd("/Users/jesse/Documents/Virginia Tech/Research/Temp + Immunity 2022/RAW DATA")
library(ggplot2)
library(tidyverse)
library(glmmTMB)
library(lme4)
library(DHARMa)

ti <- read.csv("ti_merged_data.csv")

#combine dpi 9 and 10
ti$dpi <- ifelse(ti$dpi %in% c(9, 10), "9", as.integer(ti$dpi))
ti$dpi <- as.integer(ti$dpi)
unique(ti$dpi)
ti$quant
#add infected and seropositivity thresholds
ti.cont <- ti%>%
  filter(treatment == "Control")%>%
  drop_na(quantity)

ti.cont
range(ti.cont$quantity) #highest control quantity = 95.38
ti$threshold_cutoff = 96
ti.cont$elisa_od
ti$seropos_cutoff = 0.061
ti$sympt_cutoff = 0.1


#Add column for total eye score by adding l and r eye score
ti <- ti %>%
  group_by(bird_ID) %>%
  mutate(total_eye_score=(l_eye_score+r_eye_score))

ti <- ti%>%
  mutate(diseased= ifelse(total_eye_score>sympt_cutoff, 1, 0))

#new column with whether the bird is ever diseased
ti <- ti %>%
  group_by(band_number) %>%
  mutate(ever_diseased = ifelse(any(coalesce(diseased, 0) == 1), 1, 0)) %>%
  ungroup()


ti <- ti %>%
  mutate(infected = ifelse(quantity>threshold_cutoff, 1, 0))  #generate infection data; if path load > 50 copies = infected

#new column with whether the bird is ever infected
ti <- ti %>%
  group_by(band_number) %>%
  mutate(ever_infected = ifelse(any(coalesce(infected, 0) == 1), 1, 0)) %>%
  ungroup()

#identify seropositive birds from quarantine in dataset
Qinf <- ti %>%
  filter(dpi==-28)%>%
  filter(elisa_od >=0.061)%>%
  select(band_number, elisa_od, groups, dpi)
Qinf


#remove birds that were seropositive dpi -28 from antibody analysis
ti <- ti %>%
  filter(band_number != 2667 & band_number != 2750)

#remove dpi that didn't have ELISA data run
ti <- ti %>%
  filter(dpi == -28 | dpi == 9 | dpi == 28)

####Analysis####
### Overall Hypotheses: Pathogen-induced antibody response will differ between hot and cold rooms
#Ho: Antibody levels were the same between hot and cold infected groups
#Ha: Antibody levels differed between hot and cold infected groups

#Across all days post-infection (dpi 9 and 28)
ti.pi <- ti %>%
  filter(dpi > 0)
unique(ti.pi$dpi)

ti.pi$dpi <- as.factor(ti.pi$dpi)

lm1<-glmmTMB(elisa_od~treatment+temp+dpi + treatment:dpi + 
               (1|band_number), data=ti.pi, 
           family=Gamma(link = "inverse"))


summary(lm1)
car::Anova(lm1, type = "III")

simulationOutput <-simulateResiduals(lm1, plot = T)

library(emmeans)
library(multcomp)
means <- as.data.frame(cld(emmeans(lm1, pairwise ~ treatment:dpi, type = "response",
            adjust = "tukey")))

ggplot(means, aes(x = dpi, y = response))+
  geom_crossbar(aes(ymin = asymp.LCL, ymax = asymp.UCL, fill = treatment),
                alpha = 0.5)+
  geom_jitter(data = ti.pi, aes(x = dpi, y = elisa_od, color = treatment),
              width = 0.1)+
  facet_wrap(~temp)+
  theme_bw()


hist(residuals(lm1))

plot(residuals(lm1), predict(lm1, type = "response"))

hist(predict(lm1, type = "response"))


ggplot(data = ti.pi, aes(x = elisa_od))+
  geom_histogram(binwidth = 0.001, fill = "white", color = "black")+
  facet_wrap(~treatment)

ggplot(data = ti.pi, aes(x = treatment, y = elisa_od))+
  geom_boxplot()+
  geom_jitter(width = 0.1)


ggplot(data = ti.pi, aes(x = dpi, y = elisa_od, color = treatment))+
  geom_boxplot()+
  geom_jitter(width = 0.1)+
  facet_wrap(~temp)


hist((ti.pi$elisa_od), breaks = 100, col = ti.pi$treatment)
range(ti.pi$elisa_od, na.rm = T)
lm1<-lmer(elisa_od~treatment+temp + (1|band_number), data=ti.pi)

simulateResiduals(lm1, plot = T)

plot(effects::allEffects(lm1))

lm1.5<-glm(elisa_od~treatment+temp, data=ti.pi%>% filter(dpi==28), family=Gamma())
summary(lm1.5)

hist(resid(lm1))
hist(resid(lm1))
qqnorm(resid(lm1))
qqline(resid(lm1))
plot(residuals(lm1), predict(lm1))

#model selection
a1<-glmmTMB(elisa_od~treatment*temp + (1|band_number), data=ti.pi, family=Gamma(log))
a2<-glmmTMB(elisa_od~treatment+temp + (1|band_number), data=ti.pi, family=Gamma(log))
a3<-glmmTMB(elisa_od~temp + (1|band_number), data=ti.pi, family=Gamma(log))
a4<-glmmTMB(elisa_od~treatment + (1|band_number), data=ti.pi, family=Gamma(log))
a5<-glmmTMB(elisa_od~1 + (1|band_number), data=ti.pi, family=Gamma(log))
a6<-glmmTMB(elisa_od~treatment+temp + sex + (1|band_number), data=ti.pi, family=Gamma(log))
#a6<-glmmTMB(elisa_od~treatment + temp + (1|dpi) + (1|band_number), data=ti.pi, family=Gamma())
#a7<-glmmTMB(elisa_od~treatment + (1|dpi) + (1|band_number), data=ti.pi, family=Gamma())
#a8<-glmmTMB(elisa_od~temp + (1|dpi) + (1|band_number), data=ti.pi, family=Gamma())
#a9<-glmmTMB(elisa_od~treatment * temp + (1|dpi) + (1|band_number), data=ti.pi, family=Gamma())

#temperature is correlated with band number
#band_number only - more conservative and don't really need both

a12 <-glmmTMB(elisa_od~treatment+temp + (1|band_number), data=ti.pi, family=Gamma())
a13 <-glmmTMB(elisa_od~treatment+temp + (1|dpi), data=ti.pi, family=Gamma())
summary(a12)
summary(a13)
summary(a11)
resid <- simulateResiduals(a11)
plot(resid)
summary(a11)
hist(ti.pi$elisa_od)
aictab(cand.set=list(a1, a2, a3, a4, a5, a6),
       modnames=c("a1", "a2", "a3", "a4", "a5", "a6"))
summary(a2)
plot(allEffects(a6))



#model 7 is the best, but model 6 is second best - use to keep temp

#day 28
ti28<- ti%>%
  filter(dpi == 28)

lm2<- glm(elisa_od~treatment+temp, data=ti28, family=Gamma())
summary(lm2)
plot(allEffects(lm2))

#model selection
b1<-glm(elisa_od~treatment*temp, data=ti28, family=Gamma())
b2<-glm(elisa_od~treatment+temp, data=ti28, family=Gamma())
b3<-glm(elisa_od~temp, data=ti28, family=Gamma())
b4<-glm(elisa_od~treatment, data=ti28, family=Gamma())
b5<-glm(elisa_od~1, data=ti28, family=Gamma())


aictab(cand.set=list(b1, b2, b3, b4, b5), modnames=c("b1", "b2", "b3", "b4", "b5"))



#b4 is the best model

summary(b4)
#MG treatment had a significant impact on antibody levels (Estimate = 18.67, SE = 0.59, p=0.0113).

summary(b2)
plot(allEffects(b2))
#MG treatment (Estimate = 19.15, SE=0.79), but not temperature (Estimate = -0.966, SE = 1.03) had a significant impact on antibody levels.

####Model Predictions####
dat.new=expand.grid(temp=unique(ti$temp),
                    treatment=unique(ti$treatment))#new grid to put predictions into
dat.new$yhat = predict(b2, type="response", newdata=dat.new, re.form=NA) #predicted values based off glm.phago model
head(dat.new)

dat.new$inf_temp <-  paste(dat.new$treatment, dat.new$temp, sep = "_")

#plot predicted values over raw data
ab.pred <- ggplot(data=ti28, aes(x=treatment, y=elisa_od, shape=temp))+
  geom_jitter(size=2, width=0.1)+
  geom_line(data=dat.new, aes(x=treatment, y=yhat, linetype=temp, group=temp), size=1)+ #model predictions
  stat_summary(aes(group=treatment), fun=mean, alpha=1, size=2, shape="-")+
  geom_hline(aes(yintercept=0.061), linetype="dashed", alpha=0.75)+
  labs(x="Infection Treatment", y="Antibody Levels PID 28", linetype="Temperature", shape="Temperature")+
  scale_shape_manual(values=c(15, 0))+
  scale_linetype_manual(values = c("solid", "dashed"))+
  scale_x_discrete(labels = c("Control", "Infected")) +
  #scale_colour_manual(values=c("cadetblue4", "blue", "hotpink2","red"))
  theme_minimal()

ab.pred

#are there differences in elisa_od between all treatment groups?
g.ab <- ggplot(data=ti %>% filter(dpi == 28), aes(x=groups, y=elisa_od, shape=groups))+
  geom_jitter(width = .1)+
  stat_summary(fun.y=mean,
               fun.min = function(y) mean(y)-sd(y),
               fun.max = function(y) mean(y)+sd(y),
               geom= "errorbar", size=0.5, width=0.1, alpha=0.5)+ #error bars
  stat_summary(aes(group=groups), fun=mean, alpha=1, size=2, shape="-")+
  geom_hline(aes(yintercept=0.061), linetype="dashed", alpha=0.75)+
  labs(x="Treatment", y="Antibody Levels PID 28", shape="Treatment")+
  scale_shape_manual(values = c(1, 16, 0, 15))+
  theme_minimal()

g.ab

#split violin plot attempt
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin,
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])

                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                         1))
                               quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           })

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ...,
                              draw_quantiles = NULL, trim = FALSE, scale = "area", na.rm = FALSE,
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin,
        position = position, show.legend = show.legend, inherit.aes = inherit.aes,
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

#split violin plot
ggplot(ti %>% filter(dpi != -28 & treatment=="Infected"), aes(x=0, y=elisa_od, fill = fct_rev(temp)))+
  geom_split_violin(alpha=0.75)+
  geom_point(
    aes(group = fct_rev(temp)),  # Group by temp
    position = position_dodge(width = 0.025),  # Adjust width as needed
    size = 5,
    color="black",
    alpha=0.5,
    shape = "-") +
  stat_summary(aes(group = fct_rev(temp)),  # Group by temp
               fun = mean, color="black", geom = "point", size = 10, shape = "-", alpha=1,
               position = position_dodge(width = 0.04))+  # Adjust width as needed
  scale_fill_manual(values = c("Cold" = "white", "Warm" = "gray30")) +  # Set color for each level
  theme_minimal()+
  labs(title = "ELISA", x = "Infection", y = "ELISA OD", fill = "Temperature", color= "Temperature")+
  geom_hline(yintercept=0.061, alpha =0.75)+
  facet_wrap(~dpi)


#are there differe
g.ab.mg <- ggplot(data=ti28, aes(x=groups, y=elisa_od, shape=groups, color=temp))+
  geom_jitter(size=2, width=0.1)+
  stat_summary(aes(group=groups), fun=mean, alpha=1, size=2, shape="-")+
  stat_summary(fun.y=mean,
               fun.min = function(y) mean(y)-sd(y),
               fun.max = function(y) mean(y)+sd(y),
               geom= "errorbar", size=0.1, width=0.1, alpha=0.75)+ #error bars
  geom_hline(aes(yintercept=0.061), linetype="dashed", alpha=0.75)+
  labs(x="Infection Treatment", y="Antibody Levels PID 28", linetype="Temperature", shape="Temperature")+
  #scale_shape_manual(values=c(1, 16))+
  scale_color_manual(values=c("blue","red"))+
  theme_minimal()

g.ab.mg

ti%>%
  filter(elisa_od > 0.061 & ever.inf==0)


#Change in ab over all dpi
 g.ab.change <- ggplot(data=ti, aes(x=as.factor(dpi), y=elisa_od), color = groups, groups = band_number)+
   geom_point(aes(color=groups),size=0.4)+
  geom_line(aes(group = band_number, color = groups), size=0.5, linetype="solid", alpha=0.25)+
  stat_summary(aes(group=groups, color=groups), fun=mean, geom="line", linetype="solid", alpha=1, size=1)+
  stat_summary(aes(color=groups),fun.y=mean,
               #fun.min = function(x) mean(x)-sd(x),
               #fun.max = function(x) mean(x)+sd(x),
               #fun.min = function(x) mean(x)-(sd(x)/mean(x)), #CV
               #fun.max = function(x) mean(x)+(sd(x)/mean(x)), #CV
               geom= "pointrange", size=1.5, shape="+")+
  geom_hline(aes(yintercept=0.061), linetype='dashed')+
  ylab("MG Ab OD")+
  xlab("Days Post Infection")+
  labs(y = "Antibody Levels", x="Days Post Infection", shape= "Treatment Group", color="Treatment Group")+
  scale_color_manual(values=c("lightskyblue", "blue", "lightpink2", "red"))+
   theme_gray()+
  facet_wrap(~groups~ever_infected)
g.ab.change

#####Variance#####

#COME BACK TO THIS ONE

#subset antibodies by groups to look at their variance
#subset by groups (number doesn't matter) with 1|band_number as random effect
#variance of random effect should be highest in highest heterogeneity groups
library(lme4)
ti <- ti %>%
  drop_na(elisa_od)
#across all  days, which groups are most heterogeneous?
glm.all <- glmmTMB(elisa_od~1 +(1|band_number), data=ti, family=Gamma())
glm.tnmg <- glmmTMB(elisa_od~1 +(1|band_number), data=subset(ti, groups == "Warm Infected"), family=Gamma())
glm.tnc <- glmmTMB(elisa_od~1 +(1|band_number), data=subset(ti, groups == "Warm Control"), family=Gamma())
glm.stmg <- glmmTMB(elisa_od~1 +(1|band_number), data=subset(ti, groups == "Cold Infected"), family=Gamma())
glm.stc <- glmmTMB(elisa_od~1 +(1|band_number), data=subset(ti, groups == "Cold Control"), family=Gamma())

#summary(glm.all.nest)
summary(glm.all); summary(glm.tnmg); summary(glm.tnc); summary(glm.stmg); summary(glm.stc)
#options - use SD to get a CI around each four variance estimates OR
#add up the log lik from the 4models and calculate AIC

#glm.all: V = 1.092; AIC = -976.1; SD = 1.045; LogLik = 491.1
#glm.tnmg: V = 1.813; AIC = -288.6; SD = 1.346; LogLik = 147.3
#glm.tnc:  V = 5.099e-09; AIC = -215.9; SD = 7.141e-05; LogLik = 111.0
#glm.stmg: V = 0.6598; AIC = -309.6; SD = 0.8123; LogLik = 157.8
#glm.stc: V = 0.147; AIC = -173.2; SD = 0.3834; LogLik = 89.6

AIC(glm.all, glm.tnmg, glm.tnc, glm.stmg, glm.stc)

#quick graph showing differences in variance by group
var <- c(1.092, 1.813, 5.099e-09, 0.6598, 0.147) #variance from table above
ci <- c(1.045, 1.346, 7.141e-05, 0.8123, 0.3834) #Standard Deviation from table above
g <- c("All Groups", "Warm Infected", "Warm Control", "Cold Infected", "Cold Control") #Groups
temp <- c("all", "Warm", "Warm", "Cold", "Cold")

#new df with variance and error for error bars (min/max)
variability <- data.frame(g,var, ci, temp)
variability$min <- variability$var - variability$ci #new column for lowerbound error
variability$max <- variability$var + variability$ci #new column for upperbound error

g.var <- ggplot(variability %>% filter(g != "All Groups"), aes(x=g, y=var, shape=g))+
  geom_point(aes(shape=g), size=3)+
  #scale_shape_manual(values = c( 0, 15, 1, 16))+
  scale_shape_manual(values = c( 0, 1, 2, 3))+
  #scale_color_manual(values=c("blue", "red"))+
  #geom_errorbar(aes(ymin=var, ymax=max), width=0.5)+ #add errorbars representing +/- 1 SD
  labs(x="Treatment Groups", y="Variance", shape="Treatment Groups", color="Temperature")+
  ylim(c(0,4))+
  #scale_fill_manual(values=c("black", "gray75", "white", "gray35", "gray20"))+
  theme_minimal()
g.var


####Histograms####

treat_names_n <- c(
  'Cold Infected' = "Infected Cold (n= x)",
  'Warm Control' = "Control Warm (n= x)",
  'Warm Infected' = "Infected Warm (n= x)",
  'Cold Control' = "Control Cold (n= x)"
)
#new df with labels for graphs - treatment only
treat_names <- c(
  'Cold Infected' = "Infected Cold",
  'Warm Control' = "Control Warm",
  'Warm Infected' = "Infected Warm",
  'Cold Control' = "Control Cold"
)

#Add variable indicating the elisa_od mean for each treatment group on each dpi
ti.ab <- ti%>%
  group_by(groups, dpi)%>%
  drop_na(elisa_od)%>%
  mutate(groupmean.dpi = mean(elisa_od))

ti.ab$groupmean.dpi

ti.ab <- ti.ab %>%
  group_by(groups)%>%
  drop_na(elisa_od)%>%
  mutate(groupmean = mean(elisa_od))
ti.ab$groupmean

#histograms of treatment groups across all sample days
ggplot(ti.ab, aes(x=elisa_od, fill=groups))+
  geom_histogram(binwidth = 0.002, position="identity", color="black", alpha=1)+
  #scale_fill_manual(values=c("lightskyblue", "blue", "lightpink2", "red"))+
  scale_fill_manual(values=c("gray75", "white", "gray35", "black"))+
  geom_vline(xintercept=0.061, linetype="dashed")+
  labs(fill="Treatment Group", x= "ELISA OD", y= "Count")+
  geom_vline(data=ti.ab,aes(xintercept=groupmean, color="red"), linetype="dotted", alpha=1, size=0.75, show.legend=FALSE)+
  facet_grid(c("groups"))

#density plot of treatment groups across all sample days
ggplot(ti.ab, aes(x=elisa_od, fill=groups))+
  geom_density()+
  #scale_fill_manual(values=c("lightskyblue", "blue", "lightpink2", "red"))+
  scale_fill_manual(values=c("gray75", "white", "gray35", "black"))+
  geom_vline(xintercept=0.061, linetype="dashed")+
  geom_vline(data=ti.ab,aes(xintercept=groupmean, color="red"), linetype="dotted", alpha=1, size=0.75, show.legend=FALSE)+
  facet_grid(c("groups"))



  #Histogram showing distribution of elisa_od by dpi and treatment groups
  ggplot(ti.ab, aes(x=elisa_od, fill=groups))+
    #geom_density(color="black", alpha=0.5)+
    geom_histogram(binwidth = 0.002, position="identity", alpha=2, color="black", size=0.35)+
    geom_vline(xintercept=0.061, linetype="dashed", alpha=0.5)+
    geom_vline(data = ti.ab, aes(xintercept = groupmean.dpi, color = "red"), linetype = "dotted",alpha = 1, show.legend=FALSE) +
    labs(y= "Count", x= "ELISA OD", fill="Treatment Groups", color= "Average OD")+
    scale_fill_manual(values=c("gray75", "white", "gray35", "black"))+
    scale_color_manual(values = c("red"), guide = guide_legend(title= NULL))+
    facet_grid(c("dpi","groups"))+
    theme_bw()

  #Density plot showing distribution of elisa_od by dpi and treatment groups
  ggplot(ti.ab, aes(x=elisa_od, fill=groups))+
    geom_density()+
    #geom_histogram(binwidth = 0.002, position="identity", alpha=2, color="black", size=0.35)+
    geom_vline(xintercept=0.061, linetype="dashed", alpha=0.5)+
    geom_vline(data = ti.ab, aes(xintercept = groupmean.dpi, color = "red"), linetype = "dotted",alpha = 1, show.legend=FALSE) +
    labs(y= "Count", x= "ELISA OD", fill="Treatment Groups", color= "Average OD")+
    scale_fill_manual(values=c("gray75", "white", "gray35", "black"))+
    scale_color_manual(values = c("red"), guide = guide_legend(title= NULL))+
    facet_grid(c("dpi","groups"))+
    theme_bw()

####Calculate CV####
  cv.all <- ti.ab %>%
    filter(dpi !=-28)%>%
    group_by(groups, temp)%>%
    summarise(bird_cv = sd(elisa_od)/mean(elisa_od),
              bird_sd = sd(elisa_od))

  #generate error bar values
  cv.all$max <- cv.all$bird_cv + cv.all$bird_sd
  cv.all$min <- cv.all$bird_cv - cv.all$bird_sd

  #Graph of CV calculated from dpi 9 to 28 with error bars +/- 1 SD
  cv.all.primary <- ggplot(cv.all, aes(x=groups, y=bird_cv, shape=groups))+
    geom_point(aes(shape=groups), size=3)+
    scale_shape_manual(values = c( 0, 15, 1, 16))+
    scale_color_manual(values=c("blue", "red"))+
    labs(x="Treatment Groups", y="CV", shape="Treatment Groups", color="Temperature")+
    ylim(c(0, 0.3))+
    theme_minimal()

  cv.all.primary

ggplot(ti.ab %>% filter(dpi == 9 | dpi == 28), aes(x=groups, y=elisa_od, shape=groups))+
  geom_jitter(width = 0.25)+
  facet_wrap(~ever_diseased)

#CV by day
cv.all.day <- ti.ab %>%
  group_by(dpi, groups, treatment, temp) %>%
  summarise(
    bird_cv = sd(elisa_od) / mean(elisa_od),
    bird_sd = sd(elisa_od),
    max = bird_cv + bird_sd,
    min = bird_cv - bird_sd
  ) %>%
  ungroup()  # Remove grouping for further operations

  cv.primary <- ggplot(cv.all.day, aes(x=dpi, y=bird_cv, shape=groups, color=temp))+
    geom_point(size=3)+
    geom_line(aes(linetype=treatment))+
    scale_shape_manual(values = c( 0, 15, 1, 16))+
    scale_color_manual(values=c("blue", "red"))+
    scale_linetype_manual(values = c("dashed", "solid"))+
    labs(x="Days Post Infection", y="CV", shape="Infection Treatment",
         linetype="Infection Treatment", color="Temperature")+
    #scale_fill_manual(values=c("gray75", "white", "gray35", "black"))+
    #                  labels = c("High", "Low", "Sham"))+
    theme_bw()

  cv.primary

  ggplot(ti.ab, aes(x=groups, y=elisa_od, shape=groups))+
    geom_jitter(width = 0.25)+
  facet_wrap(~dpi)

unique_counts <- ti %>%
  #filter(!is.na(elisa_od)) %>%
  group_by(groups) %>%
  group_by(dpi) %>%
  summarise(unique_band_num = n_distinct(band_number))

unique_counts


unique_group_counts <- ti %>%
  filter(!is.na(elisa_od)) %>%
  count(dpi, band_number) %>%
  split(.$n)

unique_group_counts

#calculate difference from baseline to get difference in Ab OD
ti.t <- ti %>%
  drop_na(band_number, elisa_od, dpi, treatment)

ti.t$band_number_c = as.character(ti.t$band_number)

unique(ti.t$dpi)
str(ti.t$dpi)
unique(ti.t$band_number_c)
ti.t$timepoint[ti.t$dpi==(-28)]=1
ti.t$timepoint[ti.t$dpi ==(9) | ti.t$dpi==(10)]=2
ti.t$timepoint[ti.t$dpi==(28)]=3
ti.t$timepoint

ti.t = ti.t %>%
  drop_na(timepoint) %>%
  dplyr::group_by(band_number_c)%>%
  dplyr::arrange(band_number, timepoint)%>%
  mutate(ab.time.point.neg28 = dplyr::lag(elisa_od))%>%
  mutate(ab.diff = dplyr::lag(elisa_od)-elisa_od)


View(
  ti.t %>%
    dplyr::select(c(band_number,dpi,timepoint,elisa_od, ab.time.point.neg28, ab.diff))
)

ti.t$ab_diff <- -1*ti.t$ab.diff

#change in antibodies from baseline to days 9+10, and from days 9+10 to 28
g.ab.change <- ggplot(data=ti.t %>%
                        filter(elisa_cv<=10) %>%
                        filter(timepoint>1)%>%
                 filter(dpi %in% c(-28, 9,10, 28)),
               aes(x=groups, y=ab_diff, color = groups))+
  stat_summary(aes(group=groups, color=groups), fun=mean, geom="line", alpha=1, size=2)+
  stat_summary(aes(group=groups), fun.y=mean,
               fun.min = function(x) mean(x)-sd(x),
               fun.max = function(x) mean(x)+sd(x),
               geom= "pointrange", color="black", size=0.5)+
  geom_jitter(size=1)+
  geom_point(size=2, shape=2, stroke=0)+
  ylab("Delta MG Antibodies OD")+
  xlab("Treatment Group")+
  labs(color="Treatment Group")+
  scale_colour_manual(values=c("violetred",  "red", "cadetblue","blue"))+
  facet_wrap(~timepoint)

g.ab.change

#calculate difference from baseline to get overall difference in Ab OD
ti.tr <- ti %>%
  drop_na(band_number, elisa_od, dpi, treatment)

ti.tr$band_number_c = as.character(ti.tr$band_number)

unique(ti.tr$dpi)
str(ti.tr$dpi)
unique(ti.tr$band_number_c)
ti.tr$timepoint[ti.tr$dpi==(-28)]=1
ti.tr$timepoint[ti.tr$dpi==(28)]=2
ti.tr$timepoint

ti.tr = ti.tr%>%
  drop_na(timepoint) %>%
  dplyr::group_by(band_number_c)%>%
  dplyr::arrange(band_number, timepoint)%>%
  mutate(ab.time.point.neg28 = dplyr::lag(elisa_od))%>%
  mutate(ab.diff = dplyr::lag(elisa_od)-elisa_od)


View(
  ti.tr %>%
    dplyr::select(c(band_number,dpi,timepoint,elisa_od, ab.time.point.neg28, ab.diff))
)

ti.tr$ab_diff <- -1*ti.tr$ab.diff

#total difference in ab OD from baseline to pid28
g.ab.c.total <- ggplot(data=ti.tr %>%
                        filter(elisa_cv<=10) %>%
                        filter(dpi %in% c(-28, 28)),
                      aes(x=groups, y=ab_diff, shape = groups))+
  stat_summary(aes(group=groups), fun=mean, geom="point", shape = "-", alpha=1, size=10)+
  #stat_summary(aes(group=groups), fun.y=mean,
  #             fun.min = function(x) mean(x)-sd(x),
  #             fun.max = function(x) mean(x)+sd(x),
  #             geom= "pointrange", color="black", size=0.5)+
  geom_jitter(size=3, width=0.25)+
  ylab("Delta MG Antibodies OD")+
  xlab("Treatment Group")+
  labs(color="Treatment Group", shape="Treatment Group", x="Treatment Group", y="Change in Antibody Levels from Baseline")+
  scale_shape_manual(values = c(0, 15, 1, 16))

g.ab.c.total

library(effects)
library(DHARMa)
#models
lm1 <- lm(elisa_od ~ room*dpi+treatment*dpi, data=ti)
summary(lm1)
plot(allEffects(lm1))
res<-simulateResiduals(lm1)
plot(res)

hist(ti$elisa_od)
hist(log10(ti$elisa_od))
library(glmmTMB)
lm2<-glmmTMB(elisa_od ~ treatment*dpi + room*dpi + (1|band_number), family=gaussian, data=ti %>% filter(elisa_cv<=10))
summary(lm2)
plot(allEffects(lm2))

lm2.5<-glmmTMB(elisa_od ~ treatment*room*dpi + (1|band_number), family=gaussian, data=ti %>% filter(elisa_cv<=10))
summary(lm2.5)
plot(allEffects(lm2.5))

phago.pred
lm3 <- lm(ab_diff~treatment+room, data=ti.t %>% filter(elisa_cv<=10))
summary(lm3)
plot(allEffects(lm3))
hist(ti.t$ab_diff)
lm3.5 <- lm(ab_diff~treatment+room + (1|band_number), data=ti.t %>% filter(elisa_cv<=10))
summary(lm3)
plot(allEffects(lm3))

lm4 <- lm(elisa_od ~ dpi*treatment+dpi*room, data=ti %>% filter(elisa_cv<=10))
summary(lm4)
plot(allEffects(lm4))

lm4 <- lm(elisa_od ~ as.factor(dpi)*treatment+as.factor(dpi)*room, data=ti %>% filter(elisa_cv<=10))
summary(lm4)
plot(allEffects(lm4))

hist(ti.t$ab_diff)
lm5 <-glmmTMB(ab_diff~room*timepoint + treatment*timepoint +(1|band_number), data=ti.t %>% filter(elisa_cv<=10), family=gaussian)
summary(lm5)
plot(allEffects(lm5))

dat.new=expand.grid(room=unique(ti.t$room),
                    treatment=unique(ti.t$treatment),
                    timepoint = unique(ti.t$timepoint),
                    band_number=unique(ti.t$band_number),
                    groups=unique(ti.t$groups))
dat.new$yhat = predict(lm5, type="response", newdata=dat.new, re.form=NA)
head(dat.new)

#Plot predicted values
ab.pred <- ggplot(data=ti.t, aes(x=treatment, y=ab_diff, color=groups))+
  geom_jitter(size=2, width=0.1)+
  geom_line(data=dat.new, aes(x=treatment, y=yhat, col=groups, group=groups))+
  ylab("Antibody Change")+
  xlab("Treatment")+
  labs(x="Infection Treatment", y="Antibody Change", color="Temperature")+
  scale_colour_manual(values=c("violetred",  "red", "cadetblue","blue"), labels = c("Warm Sham", "Warm MG", "Cold Sham", "Cold MG"))
ab.pred

