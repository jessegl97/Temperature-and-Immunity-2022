#T+I 2022 Phagocytosis Analysis
rm(list=ls())
####read in + format data####
setwd('/Users/jesse/Documents/Virginia Tech/Research/Temp + Immunity 2022/RAW DATA/')
library(ggplot2)
library(tidyverse)
library(lme4)

ti <- read.csv("ti_merged_data.csv")


#select only phagocytosis assay and make new wbc total, phago total, and phago score columns
ti <- ti %>%
  filter(Well.1.WBC != "NA") %>% #select only phagocytosis assay
  mutate(wbc_total=Well.1.WBC+Well.2.WBC+Well.3.WBC+Well.4.WBC)%>% #sum of wbc quadruplicates
  mutate(phago_total=Well.1.Phago+Well.2.Phago+Well.3.Phago+Well.4.Phago) %>% #sum of phago quadruplicates
  mutate(phago_score = phago_total/(wbc_total+phago_total)) #phago score = phagocytic / total

####Analysis####
### Overall Hypotheses: Both baseline and pathogen-induced immune response will differ
                       #between hot and cold rooms

par(mfrow=c(2,2))
pr <- function(m) printCoefmat(coef(summary(m)), digits=3,signif.stars=FALSE) #function for model interpretation

ti <- ti%>%
  filter(dpi != 2) %>% #remove dpi 2 because only a small subset were run day 2
  drop_na(phago_score)

#PHAGOCYTOSIS MODEL#
#Binomial model testing whether temperature or MG treatment affected phagocytosis 
#Each cell is either phagocytic or not
#Therefore binomial model

#null hypothesis: H0: phagocytosis score is not affected by temperature or MG infection
#alternative hypothesis: Ha1: phagocytosis score is higher in MG infected birds
#Ha2: phagocytosis score is lower in cold birds
glm.phago <- glmer(phago_score~temp+treatment + (1|band_number), weights=wbc_total+phago_total, data=ti, family="binomial")
summary(glm.phago)

library(glmmTMB)
glm1 <- glmmTMB(phago_score ~ temp + treatment + (1|band_number), weights=wbc_total+phago_total, data=ti, family="binomial")
summary(glm1)
#Model Selection AICc

p1 <- glmer(phago_score~temp+treatment + (1|band_number), weights=wbc_total+phago_total, data=ti, family="binomial")
p2<- glmer(phago_score~temp*treatment + (1|band_number), weights=wbc_total+phago_total, data=ti, family="binomial")
p3<- glmer(phago_score~temp + (1|band_number), weights=wbc_total+phago_total, data=ti, family="binomial")
p4<- glmer(phago_score~treatment + (1|band_number), weights=wbc_total+phago_total, data=ti, family="binomial")
p5<- glmer(phago_score~1 + (1|band_number), weights=wbc_total+phago_total, data=ti, family="binomial")
p6 <- glmer(phago_score~groups + (1|band_number), weights=wbc_total+phago_total, data=ti, family="binomial")
p7 <- glmer(phago_score~temp+treatment + sex + (1|band_number), weights=wbc_total+phago_total, data=ti, family="binomial")
aictab(cand.set=list(p1, p2, p3, p4, p5, p6, p7), modnames=c("p1", "p2", "p3", "p4", "p5", "p6", "p7"))

#Model selection based on AICc:
#  
#  K   AICc Delta_AICc AICcWt Cum.Wt      LL
#p4 3 482.64       0.00   0.35   0.35 -238.05
#p1 4 483.08       0.44   0.28   0.62 -237.07
#p5 2 484.72       2.08   0.12   0.74 -240.22
#p3 3 485.19       2.56   0.10   0.84 -239.32
#p6 5 485.58       2.94   0.08   0.92 -237.07
#p2 5 485.58       2.94   0.08   1.00 -237.07

#p4 best model, but use p1 to keep both temp and treatment in model
summary(p1)
summary(p7)
plot(resid(glm.phago))

summary(glm.phago)

#Fixed effects:
#               Estimate Std. Error z value Pr(>|z|)    
#(Intercept)    -2.2724     0.1941 -11.709   <2e-16 ***
#tempTN          0.3294     0.2349   1.403   0.1608    
#treatmentSham  -0.5186     0.2384  -2.175   0.0296 *  

#infection with MG, but not cold temperatures increase phagocytosis

library(DHARMa)
p1.resid<- simulateResiduals(glm.phago)
plot(p1.resid) #residuals normal; deviation ns; good QQ

plot(allEffects(glm.phago))
#trend towards decreased phagocytosis in ST group, but not significant

#interactive model temp*treatment
summary(p2)
#Fixed effects:
#  Estimate Std. Error z value Pr(>|z|)    
#(Intercept)          -2.273206   0.218143 -10.421   <2e-16 ***
#  tempTN                0.331082   0.306384   1.081    0.280    
#treatmentSham        -0.516567   0.339014  -1.524    0.128    
#tempTN:treatmentSham -0.004014   0.476925  -0.008    0.993  

#there is no significant affect of temperature and treatment on phago score
p2.resid <- simulateResiduals(p2)
plot(p2.resid) #residuals normal; deviation ns; good QQ
plot(allEffects(p2)) 
#it certainly looks like infection increases phagocytosis, and that the increase is slightly abrogated by cold temperatures

#install.packages("betareg")
library(betareg)

#beta regression
beta_mod <- betareg(phago_score ~ temp + treatment, data=ti)
summary(beta_mod)
plot(allEffects(beta_mod))
plot(resid(beta_mod))

emm_results <- emmeans(beta_mod, ~ temp + treatment, scale = "response")
pairs(emm_results)

####Model Predictions####
dat.new=expand.grid(temp=unique(ti$temp),
                   treatment=unique(ti$treatment))#new grid to put predictions into
dat.new$yhat = predict(beta_mod, type="response", newdata=dat.new, re.form=NA) #predicted values based off glm.phago model

dat.new$inf_temp <-  paste(dat.new$treatment, dat.new$temp, sep = "_")

#plot predicted values over raw data
phago.pred <- ggplot(data=ti, aes(x=treatment, y=phago_score, shape=temp))+
  geom_jitter(size=2, width=0.1)+
  geom_line(data=dat.new, aes(x=treatment, y=yhat, linetype=temp, group=temp), size=1)+ #model predictions
  stat_summary(aes(group=treatment), fun=mean, alpha=1, size=3, shape="-")+
  stat_summary(data=ti, aes(group=treatment), fun.y=mean,
               fun.min = function(y) mean(y)-sd(y),
               fun.max = function(y) mean(y)+sd(y),
               geom= "errorbar", size=0.5, width=0.1, alpha=0.75)+ #error bars
  labs(x="Infection Treatment", y="Phagocytosis Score", linetype="Temperature", shape="Temperature")+
  scale_shape_manual(values=c(0, 15),  labels= c("Cold", "Warm"))+
  scale_linetype_manual(values = c("dotted", "solid"), labels= c("Cold", "Warm"))+
  #scale_colour_manual(values=c("cadetblue4", "blue", "hotpink2","red"))
  theme_minimal()

phago.pred

####Graphs####

#theme_set(
#  theme_bw() +
#theme(
#  text = element_text(size = 50),
#  axis.title.y = element_text(size = 100, face = "bold"),
#  axis.text.y = element_text(size = 75,  face = "bold"),
#  axis.title.x = element_text(size = 100, face = "bold"),
#  axis.text.x = element_text(size = 75,  face = "bold"),
#  strip.text = element_text(size = 100,  face = "bold")
#)
#)

#are there differences in phago_score between all treatment groups?
g.phago <- ggplot(data=ti, aes(x=groups, y=phago_score, shape=groups))+
  geom_jitter(width = .1, size=2, stroke=1)+
  stat_summary(fun.y=mean,
               fun.min = function(y) mean(y)-sd(y),
               fun.max = function(y) mean(y)+sd(y),
               geom= "errorbar", size=0.2, width=0.1, alpha=0.75)+ #error bars
  stat_summary(aes(group=groups), fun=mean, alpha=1, size=3, shape="-")+
  labs(x="Treatment", y="Phagocytosis Score", shape="Treatment Groups")+
  scale_shape_manual(values=c(1, 2, 3, 4))
  #scale_color_manual(values=c("cadetblue4", "blue", "hotpink2","red"))+ #in case I want color later

g.phago

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
ggplot(ti, aes(x=treatment, y=phago_score, fill = fct_rev(temp)))+ 
  geom_split_violin()+
  geom_point(
    aes(group = fct_rev(temp)),  # Group by temp
    position = position_dodge(width = 0.04),  # Adjust width as needed
    size = 5,
    alpha=0.4,
    shape = "-") +
  stat_summary(aes(group = fct_rev(temp)),  # Group by temp
              fun = mean, geom = "point", size = 7, shape = "-",
               position = position_dodge(width = 0.06))+  # Adjust width as needed
  scale_fill_manual(values = c("Warm" = "gray45", "Cold" = "white")) +  # Set color for each level
  theme_minimal()+
  labs(title = "Phagocytosis", x = "Infection", y = "Phagocytosis Score", fill = "Temperature", color= "Temperature")

ti%>%
  group_by(groups)%>%
  dplyr::summarise(sd(phago_score)) #double checking stat_summary error bars; look good

inf_names <- c(
  'MG' = "Infected",
  'Sham' = "Control"
)

#are there differences in phago_score between infected and uninfected groups regardless of temperature?
g.phago.mg <- ggplot(data=ti %>% filter(treatment == "Infected"), aes(x=temp, y=phago_score, color=temp))+
  geom_jitter(width = .1, size=3, stroke = 2)+
  stat_summary(fun.y=mean,
               fun.min = function(y) mean(y)-sd(y),
               fun.max = function(y) mean(y)+sd(y),
               geom= "errorbar", size=0.5, width=0.1, alpha=0.75)+ #error bars
  stat_summary(aes(group=temp), fun=mean, alpha=1, size=3, shape="-")+
  labs(x="\nTemperature", y="Phagocytosis Score\n", shape="Temperature", color="Temperature")+
  scale_shape_manual(values=c(16, 16))+
  scale_color_manual(values=c("blue", "red"))+
  facet_grid(~temp, scales="free_x")
g.phago.mg

####Heterogeneity####

#subset phagocytosis assay by groups to look at their variance
#variance of random effect should be highest in highest heterogeneity groups
glm.phago2 <- glmer(phago_score~1 + (1|band_number), weights=wbc_total+phago_total, data=ti, family="binomial")
glm.phago2a <- glmer(phago_score~1 + (1|band_number), weights=wbc_total+phago_total, data=subset(ti, inf_temp=="MG_ST"), family="binomial")
glm.phago2b <- glmer(phago_score~1 + (1|band_number), weights=wbc_total+phago_total, data=subset(ti, inf_temp=="Sham_ST"), family="binomial")
glm.phago2c <- glmer(phago_score~1 + (1|band_number), weights=wbc_total+phago_total, data=subset(ti, inf_temp=="MG_TN"), family="binomial")
glm.phago2d <- glmer(phago_score~1 + (1|band_number), weights=wbc_total+phago_total, data=subset(ti, inf_temp=="Sham_TN"), family="binomial")
#options - use SD to get a CI around each four variance estimates OR
#add up the log lik from the 4models and calculate AIC

summary(glm.phago2a); summary(glm.phago2b); summary(glm.phago2c); summary(glm.phago2d)

#MG_ST: Variance = 0.9908; AIC = 147.9; SD = 0.9954; logLik = -72.0
#Sham_ST: Variance = 0.7769; AIC = 97.3; SD = 0.8814; logLik = -46.7
#MG_TN: Variance = 0.3096; AIC = 142.3; SD = 0.5564; logLik = -69.2
#Sham_TN: Variance = 0.5205; AIC = 98.3; SD = 0.7215; logLik = -47.2

glm.phago3 <- glm(phago_score~1 , weights=wbc_total+phago_total, data=ti, family="binomial")
summary(glm.phago2)
AIC(glm.phago2a,glm.phago2b, glm.phago2c, glm.phago2d, glm.phago3, glm.phago2)

#Cold temperatures increase heterogeneity in immune response, both baseline and pathogen induced compared to warm temperatures

#quick graph showing differences in variance by group
var <- c(0.9908, 0.7769, 0.3096, 0.5205)
sd <- c(0.9954, 0.8814, 0.5564, 0.7215)
g <- c("Infected Cold", "Control Cold", "Infected Warm", "Control Warm")
temp <- c("Cold", "Cold", "Warm", "Warm")


variability <- data.frame(g,var, sd, temp)
variability$min <- variability$var - variability$ci
variability$max <- variability$var + variability$ci

variability$g <- factor(variability$g, levels = c("Control Cold", "Infected Cold", "Control Warm", "Infected Warm"))

#infected only
g.var.i <- ggplot(variability %>% filter(g == "Infected Cold" | g == "Infected Warm"), aes(x=temp, y=var, color=temp))+
  geom_point(size=5)+
  scale_shape_manual(values = c(16, 16))+
  ylim(c(0,1.25))+
  scale_color_manual(values=c("blue", "red"))+
  labs(x="\nTemperature", y="Variance\n", shape="Temperature", color="Temperature")+
  facet_grid(~temp, scales="free_x")

g.var.i

g.var <- ggplot(variability, aes(x=g, y=var, shape=g))+
  geom_point(size=2)+
  ylim(c(0,1.25))+
  scale_shape_manual(values=c(1, 2, 3, 4))+
  labs(x="\nTemperature", y="Variance\n", shape="Temperature", color="Temperature")

g.var
####Heirarchical Models####
glm1 <- glmmTMB(phago_score ~ (1|band_number), weights=wbc_total+phago_total, data=ti, family=binomial)

glm1

ti$resid <- resid(glm1)
hist(ti$resid)
lm1 <- glm(resid ~ temp*treatment, data=ti)

anova(lm1)

int <- aov(resid ~ temp*treatment, data=ti)
summary(int)

ggplot(ti, aes(x=groups, y=resid))+
  geom_point()+
  theme_bw()

####Histograms####

#new df with labels for graphs - treatment + n
#count nubmer of observations per group
ti %>%
  group_by(groups)%>%
  count(length(band_number))

treat_names_n <- c(
  'Cold Infected' = "Infected Cold (n=14)",
  'Warm Control' = "Control Warm (n=10)",
  'Warm Infected' = "Infected Warm (n=14)",
  'Cold Control' = "Control Cold (n=10)"
)
#new df with labels for graphs - treatment only
treat_names <- c(
  'Cold Infected' = "Infected Cold",
  'Warm Control' = "Control Warm",
  'Warm Infected' = "Infected Warm",
  'Cold Control' = "Control Cold"
)

#Add variable indicating the phago_score mean for each treatment group
p.ti <- ti%>%
  group_by(groups)%>%
  mutate(groupmean = mean(phago_score))%>%
  ungroup()
p.ti$groupmean
m.ti <- p.ti %>%
  mutate(meanall = mean(phago_score)) %>%
  ungroup()
m.ti$meanall

#histogram showing variability in phagocytosis scores across treatment groups
ggplot(p.ti, aes(x=phago_score, fill=groups))+
  geom_histogram(binwidth = 0.01, position="identity", alpha=1, color="black")+
  geom_vline(data = p.ti, aes(xintercept = groupmean, color = "red"), linetype = "dotted",alpha = 1, show.legend=FALSE) + #mean of each group to show differences
  geom_vline(data=m.ti, aes(xintercept = meanall, color = "black"), linetype="dashed", alpha=1, show.legend=FALSE)+ #mean of all groups combined to compare
  scale_fill_manual(values=c( "gray75", "white", "gray35","black"))+
  scale_color_manual(values = c("black", "red"), guide = guide_legend(title= NULL))+
  facet_wrap(~groups, nrow=4, labeller = as_labeller(treat_names_n), scales="fixed")+
  labs(x="Phagocytosis Score", y="Count", fill="Treatment")+
  theme_minimal()
table(p.ti$dpi)
#density plot showing variability in phagocytosis scores across treatment groups
ggplot(p.ti, aes(x=phago_score, fill=groups))+
  geom_density()+
  geom_vline(data = p.ti, aes(xintercept = groupmean, color = "red"), linetype = "dotted",alpha = 1, show.legend=FALSE) + #mean of each group to show differences
  geom_vline(data=m.ti, aes(xintercept = meanall, color = "black"), linetype="dashed", alpha=1, show.legend=FALSE)+ #mean of all groups combined to compare
  scale_fill_manual(values=c( "gray75", "white", "gray35","black"), labels=treat_names)+
  scale_color_manual(values = c("black", "red"), guide = guide_legend(title= NULL))+
  facet_wrap(~groups, nrow=4, labeller = as_labeller(treat_names_n), scales="fixed")+
  labs(x="Phagocytosis Score", y="Count", fill="Treatment")+
  theme_minimal()

####CVs####
#calculate CV Kate's way
p.ti$N = p.ti$wbc_total + p.ti$phago_total #add column with total number of cells per band_number

m.cv.all <- p.ti %>% 
  group_by(groups)%>%
  summarise(mean.phago = mean(phago_score, na.rm=T), #mean
            all.Ns = sum(N, na.rm = T), #total number of cells
            sd.binom = sqrt(all.Ns*mean.phago*(1-mean.phago)), #standard deviation = square root(N * mean * (1-mean))
            bird_cv = (sd.binom/mean.phago)) #CV = standard deviation / mean * 100

m.cv.all
#generate error bar values
m.cv.all$max <- m.cv.all$bird_cv + m.cv.all$sd.binom
m.cv.all$min <- m.cv.all$bird_cv - m.cv.all$sd.binom
m.cv.all$groups <- factor(m.cv.all$groups, levels = c("Cold Control", "Cold Infected", "Warm Control", "Warm Infected"))
#Graph of CVs
cv.all.primary <- ggplot(m.cv.all, aes(x=groups, y=bird_cv, shape=groups))+
  labs(x="Treatment Groups", y="CV", shape="Treatment Groups", color="Temperature")+
  geom_point(aes(shape=groups), size=3)+
  #scale_shape_manual(values = c(0, 15, 1, 16))+
  #scale_color_manual(values=c("blue", "red"))+
  ylim(c(0, 300))

cv.all.primary

#Jesse's CV Attempt
cvs <- p.ti %>%
  group_by(groups) %>%
  summarise(n = sum(N),
            mean_p = sum(N * phago_score),
            p = mean(phago_score),
            var_p = sum(N * phago_score * (1 - phago_score))) %>%
  mutate(cv = (sqrt(var_p / n) / p)*100)
cvs

cvs$max <- cvs$cv + cvs$var_p
cvs$min <- cvs$cv - cvs$var_p

g.cvs <- ggplot(cvs, aes(x=fct_rev(groups), y=cv, shape=groups))+
  #geom_errorbar(aes(ymax=max, ymin=min), color="black", size=0.5, width=0.5)+
  labs(x="Treatment Groups", y="CV Jesse", shape="Treatment Groups")+
  geom_point(aes(shape=groups), size=3)+
  scale_shape_manual(values = c(1, 16, 0, 15))+
  ylim(c(0, 400))+
  theme_minimal()

g.cvs

library(patchwork)
both.cvs <- ggplot() + cv.all.primary + g.cvs
both.cvs <- cv.all.primary + g.cvs + plot_layout(ncol=2)
print(both.cvs)


variability <- ggplot() + cv.all.primary + g.cvs + g.var
variability <- cv.all.primary + g.cvs + g.var + plot_layout(ncol=3)
variability

#CV of Phagocytosis Scores (not binomial)
cv.all <- p.ti %>% 
  group_by(groups)%>%
  summarise(bird_cv = sd(phago_score)/mean(phago_score),
            bird_sd = sd(phago_score))
cv.all
#generate error bar values
cv.all$max <- cv.all$bird_cv + cv.all$bird_sd
cv.all$min <- cv.all$bird_cv - cv.all$bird_sd

#Graph of CV calculated from dpi 9 to 28 with error bars +/- 1 SD
cv.all.score <- ggplot(cv.all, aes(x=groups, y=bird_cv, shape=groups))+
  geom_point(aes(shape=groups))+
  scale_shape_manual(values = c(1, 16, 0, 15))+
  #geom_col(aes(fill=groups), color="black", size=0.5)+
  #geom_errorbar(aes(ymax=max, ymin=min), color="black", size=0.5, width=0.5)+
  labs(x="Treatment Groups", y="CV", fill="Treatment")+
  ylim(c(0, 1))+
  #scale_fill_manual(values=c("gray75", "white", "gray35", "black"))+
  #                  labels = c("High", "Low", "Sham"))+
  theme_minimal()

cv.all.score

