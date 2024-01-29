#T+I 23 Eye Score Analysis
rm(list=ls())
####read in + format data####
setwd('/Users/jesse/Documents/Virginia Tech/Research/Temp + Immunity 2022/RAW DATA/')
library(ggplot2)
library(tidyverse)
library(lme4)

ti <- read.csv("ti_merged_data.csv")

theme_set(theme_bw())

#Add column for total eye score by adding l and r eye score
ti <- ti %>%
  group_by(bird_ID) %>%
  mutate(total_eye_score=(l_eye_score+r_eye_score))

#make total eye score an integer by multiplying by two for models
ti$tes <- ti$total_eye_score*2

#data frame with only inoculated birds
ti.mg <- ti%>%
  filter(treatment == "Infected")

hist(ti.mg$total_eye_score)
hist(ti$total_eye_score)
hist(ti.mg$tes)

ti.mg$sympt_cutoff = 0.1

ti.mg <- ti.mg %>%
  mutate(diseased= ifelse(total_eye_score>sympt_cutoff, 1, 0))


#new column with whether the bird is ever diseased
ti.mg <- ti.mg %>%
  group_by(band_number) %>%
  mutate(ever_diseased = ifelse(any(coalesce(diseased, 0) == 1), 1, 0)) %>%
  ungroup()



####Analysis####
### Overall Hypotheses: Both baseline and pathogen-induced immune response will differ
#between hot and cold rooms

par(mfrow=c(2,2))
pr <- function(m) printCoefmat(coef(summary(m)), digits=3,signif.stars=FALSE) #function for model interpretation

#H0: Ambient temperature did not affect eye score
#Ha: Ambient temperature did affect eye score

t1 <- ti.mg%>% 
  group_by(treatment, temp, dpi)%>%
  summarize(count = n())
t1

ti.mg <- ti.mg%>%
  dplyr::select(dpi, temp, treatment, band_number, bird_ID, sex, tes, total_eye_score, groups)

summary(ti.mg)
ti.mg<- na.omit(ti.mg)
####Eye Score Models####
lm1 <- glmmTMB(tes ~ temp + (1|band_number), data=ti.mg, family=nbinom2)
summary(lm1)
plot(allEffects(lm1))

#zero inflated model
#library(pscl)
#lm2 <- zeroinfl(2*total_eye_score ~ temp + (1|band_number), data=ti.mg, dist="negbin")
summary(lm2)
head(ti.mg)
summary(ti.mg)
ti.mg$total_eye_score_f <- factor(ti.mg$total_eye_score)
#tmblm2 <- glmmTMB(total_eye_score ~ temp + (1|band_number),
#                  ziformula = ~temp + (1|band_number),
#                  family = nbinom2, data = ti.mg)
summary(tmblm2)
plot(allEffects(tmblm2))
ti.mg$tes
#sum eye score for each bird over all dpi
sum_score <- ti.mg %>%
  group_by(band_number, temp, treatment, sex)%>%
  summarise(score = sum(total_eye_score, na.rm =TRUE))
hist(sum_score$score)
summary(ti.mg)
summary(sum_score)
table(sum_score$score)
sum_score$ts <- sum_score$score * 2
#does temperature predict total sum of eye score?
library(MASS)
hist(sum_score$ts)
table(sum_score$ts)
lm3 <- glm.nb(ts ~ temp , data=sum_score)
summary(lm3)
plot(allEffects(lm3))

#Temperature does not affect total eye score
#dotplot of total eye score in each temp
ggplot(sum_score, aes(x=temp, y=score, shape = temp)) +
  stat_summary(aes(group=temp), fun=mean, geom="point", alpha=1, size=5, shape="-")+
  stat_summary(aes(group=temp, shape = temp), fun.y=mean,
               fun.min = function(x) mean(x)-sd(x),
               fun.max = function(x) mean(x)+sd(x),
               #fun.min = function(x) mean(x)-(sd(x)/mean(x)), #CV
               #fun.max = function(x) mean(x)+(sd(x)/mean(x)), #CV
               geom= "errorbar", size=0.25, width = 0.25)+
  geom_jitter(width=0.15)+
  labs(x="Temperature", y="Sum of All Eye Scores", shape="Temperature")+
  scale_shape_manual(values=c(0, 15),  labels= c("Cold", "Warm"))+
  scale_linetype_manual(values = c("dotted", "solid"), labels= c("Cold", "Warm"))+
  scale_x_discrete(labels = c("Cold", "Warm")) +
  theme_minimal()



#calculate mean of each temp
sum_score<- sum_score %>%
  group_by(temp)%>%
  mutate(mean=mean(score))

####Heterogeneity####

#subset eye score by temperatures to look at their variance
#variance of random effect should be highest in highest heterogeneity groups
glm.eye.all <- glmmTMB(tes ~ 1 + (1|band_number), data=ti.mg, family=nbinom2)
glm.eye.cold <- glmmTMB(tes ~ 1 + (1|band_number), data=subset(ti.mg, temp=="Cold"), family=nbinom2)
glm.eye.warm <- glmmTMB(tes ~ 1 + (1|band_number), data=subset(ti.mg, temp=="Warm"), family=nbinom2)
#options - use SD to get a CI around each four variance estimates OR
#add up the log lik from the 4models and calculate AIC

#quantify w stats
glm <- glmmTMB(tes ~ dpi + (1|band_number), data=ti.mg, family=nbinom2)

glm
summary(glm)
glm1 <-glmmTMB(tes~ (1+dpi|band_number), data=ti.mg, family=nbinom2)

ti.mg$resid <- resid(glm)
hist(ti.mg$resid)
lm1 <- lm(resid ~ temp, data=ti.mg)

anova(lm1)

ggplot(ti.mg, aes(x=groups, y=resid, color=fct_rev(groups)))+
  geom_point()+
  labs(x="Groups", y="|Residuals|", color= "Groups", title = "Residuals of Total Eye Score")


summary(glm.eye.all); summary(glm.eye.cold); summary(glm.eye.warm)

#all: Variance = 7.074; AIC = 523.6; SD = 2.66; logLik = -258.8
#cold: Variance = 9.917; AIC = 234.6; SD = 3.149; logLik = -114.3
#warm: Variance = 5.486; AIC = 293.2; SD = 2.342; logLik = -143.6

AIC(glm.eye.all, glm.eye.cold, glm.eye.warm)

#Cold temperatures increase heterogeneity in immune response, both baseline and pathogen induced compared to warm temperatures

#quick graph showing differences in variance by group
var <- c(7.074, 9.917, 5.486)
ci <- c(2.66, 3.149, 2.342)
g <- c("All", "Cold", "Warm")

variability <- data.frame(g,var, ci)
variability$min <- variability$var - variability$ci
variability$max <- variability$var + variability$ci

g.var <- ggplot(variability %>% filter(g != "All"), aes(x=g, y=var, shape = g, color=g))+
  geom_point(size= 5, alpha = 1)+
  #geom_col(color="black", size=0.5, alpha=0.9)+
  labs(x="\nTemperature", y="Variance\n", shape="Temperature", color="Temperature")+
  #geom_errorbar(aes(ymax=max, ymin = min), size=0.5, width=0.25)+
  scale_shape_manual(values=c(15, 16))+
  scale_color_manual(values = c("blue", "red"))+
  ylim(c(0,12))+
  theme_minimal()

g.var

#for poster
g.var <- ggplot(variability %>% filter(g != "All"), aes(x=g, y=var, shape = g, color=g))+
  geom_point(size=40, alpha = 1)+
  #geom_col(color="black", size=0.5, alpha=0.9)+
  labs(x="\nTemperature", y="Variance\n", shape="Temperature", color="Temperature")+
  #geom_errorbar(aes(ymax=max, ymin = min), size=0.5, width=0.25)+
  scale_shape_manual(values=c(16, 16))+
  scale_color_manual(values = c("blue", "red"))+
  ylim(c(0,12))+
  facet_grid(~g, scales="free_x")

g.var

unique(ti.mg$dpi)
#calculate CVs
cv.all <- ti.mg %>%
  group_by(temp)%>%
  summarise(bird_cv = sd(total_eye_score)/mean(total_eye_score),
            bird_sd = sd(total_eye_score))

#generate error bar values
cv.all$max <- cv.all$bird_cv + cv.all$bird_sd
cv.all$min <- cv.all$bird_cv - cv.all$bird_sd
cv.all$min
cv.all$max
cv.all$k <- (1/(cv.all$bird_cv^2))

#Graph of CVs calculated all days with error bars +/- 1 SD
cv.all.g <- ggplot(cv.all, aes(x=temp, y=bird_cv, color=temp))+
  #geom_col(aes(fill=temp), color="black", size=0.5)+
  geom_point(size= 5, alpha = 1)+
  #geom_errorbar(aes(ymax=max, ymin=min), color="black", size=0.5, width=0.25)+
  #scale_fill_manual(values=c("white", "gray35"), labels = c("Cold", "Warm"))+
  #scale_shape_manual(values=c(1, 16))+
  scale_color_manual(values=c("blue", "red"))+
  ylim(c(0, 3.5))+
  labs(x="Treatment Groups", y="CV", color="Treatment")+
  theme_minimal()

cv.all.g

k.all.g <- ggplot(cv.all, aes(x=temp, y=k, color=temp))+
  #geom_col(aes(fill=temp), color="black", size=0.5)+
  geom_point(size= 5, alpha = 1)+
  #geom_errorbar(aes(ymax=max, ymin=min), color="black", size=0.5, width=0.25)+
  #scale_fill_manual(values=c("white", "gray35"), labels = c("Cold", "Warm"))+
  #scale_shape_manual(values=c(1, 16))+
  scale_color_manual(values=c("blue", "red"))+
  ylim(c(0, 1))+
  labs(x="Treatment Groups", y="CV", color="Treatment")+
  theme_minimal()

k.all.g

#Calculate CV of individual birds and the average of the bird_cvs by temp
cv.ind <- ti.mg %>%
  dplyr::select(band_number, temp, treatment, total_eye_score)%>%
  group_by(band_number, temp, treatment)%>%
  mutate(
    bird_cv = ifelse(mean(total_eye_score) != 0, sd(total_eye_score) / mean(total_eye_score), 0),
    bird_sd = ifelse(mean(total_eye_score) != 0, sd(total_eye_score), 0)
  ) %>%
  ungroup() %>%
  group_by(temp)%>%
  mutate(mean_cv = mean(bird_cv),
         count = n())

table(cv.ind$temp, cv.ind$bird_cv)

ggplot(cv.ind%>% filter(treatment == "Infected"), aes(x=temp, y= bird_cv, color=as.factor(band_number)))+
  geom_jitter(height=0)+
  geom_point(aes(x=temp, y=mean_cv), shape="-", color="black", size=5)+
  theme_bw()
  #scale_color_manual(values= c("blue", "red"))

ggplot(cv.ind%>% filter(treatment == "Infected"), aes(x=temp, y= bird_cv, color=fct_rev(temp)))+
  geom_point()+
  geom_point(aes(x=temp, y=mean_cv), shape="-", color="black", size=5)+
  theme_bw()
  #scale_color_manual(values= c("blue", "red"))

#calculate the means of the CVs
cv.ind <- cv.ind %>%
  filter(treatment == "Infected")%>%
  group_by(temp, treatment, band_number)%>%
  summarise(mean_cv = mean(bird_cv))

cv.ind

ggplot(cv.ind, aes(x=temp, y= mean_cv, color=temp))+
  geom_point()+
  theme_bw()

#CV of means
#calculate means
cv.all.m <- ti.mg %>%
  group_by(band_number, temp)%>%
  summarise(bird_mean = mean(total_eye_score))
#calculate CVs
cv.all.m <- cv.all.m %>%
  group_by(temp)%>%
  summarise(bird_cv = sd(bird_mean)/mean(bird_mean),
            bird_sd = sd(bird_mean))

#generate error bar values
cv.all.m$max <- cv.all.m$bird_cv + cv.all$bird_sd
cv.all.m$min <- cv.all.m$bird_cv - cv.all$bird_sd

#Graph of CVs calculated all days with error bars +/- 1 SD
cv.all.m.g <- ggplot(cv.all.m, aes(x=temp, y=bird_cv, shape=temp))+
  #geom_col(aes(fill=temp), color="black", size=0.5)+
  geom_point(color="black", size= 5, alpha = 1)+
  geom_errorbar(aes(ymax=max, ymin=min), color="black", size=0.5, width=0.25)+
  #scale_fill_manual(values=c("white", "gray35"), labels = c("Cold", "Warm"))+
  scale_shape_manual(values=c(1, 16))+
  ylim(c(0, 3.5))+
  labs(x="Treatment Groups", y="CV of Means", shape="Treatment")+
  theme_minimal()

cv.all.m.g

#mean tes combined
ggplot(cv.all.m, aes(x=temp, y= bird_mean, shape=temp))+
  geom_jitter(width = 0.25, height=0.15)+
  scale_shape_manual(values=c(0, 15))+
  labs(x="Temperature", y="Mean of Total Eye Score", shape = "Temperature")


#all tes combined
ggplot(ti.mg, aes(x=temp, y= total_eye_score, shape=temp))+
  geom_jitter(width = 0.25, height=0)+
  scale_shape_manual(values=c(0, 15))+
  labs(x="Temperature", y="Total Eye Score", shape = "Temperature")+
  theme_bw()

#CV by day
m.cv <- ti.mg %>% 
  group_by(dpi, temp)%>%
  summarise(bird_cv = sd(total_eye_score)/mean(total_eye_score),
            bird_sd = sd(total_eye_score))
m.cv
unique(m.cv$dpi)
#generate error bar values
m.cv$max <- m.cv$bird_cv + m.cv$bird_sd
m.cv$min <- m.cv$bird_cv - m.cv$bird_sd

cv.primary <- ggplot(m.cv %>% filter(dpi != -12 & dpi != 3 & dpi !=35), aes(x=dpi, y=bird_cv, color=temp))+
  geom_point(aes(color=temp), size= 3, alpha = 1)+
  geom_line(aes(color=temp))+
  #geom_errorbar(aes(ymax=max, ymin=min), color="black", size=0.25, width=0.5)+
  #geom_jitter(ti.mg, aes(x=temp, y=total_eye_score))+
  labs(x="Days Post Infection", y="CV", linetype="Temperature", color="Temperature", shape = "Temperature")+
  #scale_shape_manual(values=c(1, 16))+
  scale_color_manual(values=c("blue","red"))+
  #scale_linetype_manual(values = c("dashed", "solid"))+
  ylim(c(0, 3.5))+
  scale_x_continuous(breaks=unique(m.cv$dpi), labels=paste0(unique(m.cv$dpi)))+
  theme_minimal()

cv.primary

#all tes combined
ggplot(ti.mg %>% filter(dpi >3 & dpi !=35), aes(x=temp, y= total_eye_score, shape=temp))+
  geom_jitter(width = 0.25, height=0.15)+
  scale_shape_manual(values=c(0, 15))+
  labs(x="Temperature", y="Total Eye Score", shape = "Temperature")+
  facet_wrap(~dpi, nrow=1)


unique_counts <- ti %>%
  group_by(groups) %>%
  group_by(dpi) %>%
  summarise(unique_band_num = n_distinct(band_number))

unique_counts
treat_names_n <- c(
  'Cold' = "Cold (n= 14)",
  'Warm' = "Warm (n= 14)"
)



#Add variable indicating the elisa_od mean for each treatment group on each dpi
ti.mg <- ti.mg%>%
  group_by(temp, dpi)%>%
  drop_na(total_eye_score)%>%
  mutate(groupmean.dpi = mean(total_eye_score))

ti.mg$groupmean.dpi

ti.mg <- ti.mg %>%
  group_by(temp)%>%
  drop_na(total_eye_score)%>%
  mutate(groupmean = mean(total_eye_score))
ti.mg$groupmean

#histograms of treatment groups across all sample days
ggplot(ti.mg%>% filter(dpi > 3 & dpi !=35), aes(x=total_eye_score, fill=temp))+
  geom_histogram(binwidth = 0.5, position="identity", color="black", alpha=1)+
  #scale_fill_manual(values=c("lightskyblue", "blue", "lightpink2", "red"))+
  scale_fill_manual(values=c("white", "gray35"))+
  labs(y="Count", x= "Total Eye Score", fill = "Temperature")+
  geom_vline(data=ti.mg,aes(xintercept=groupmean, color="red"), linetype="dotted", alpha=1, size=0.75, show.legend=FALSE)+
  facet_grid(c("temp"))+
  theme_bw()

#histograms of treatment groups across all sample days
ggplot(ti.mg %>% filter(dpi > 3 & dpi !=35), aes(x=total_eye_score, fill=temp))+
  geom_histogram(binwidth = 0.5, position="identity", color="black", alpha=1)+
  #scale_fill_manual(values=c("lightskyblue", "blue", "lightpink2", "red"))+
  scale_fill_manual(values=c("white", "gray35"), labels = c("Cold", "Warm"))+
  geom_vline(data=ti.mg %>% filter(dpi > 3 & dpi !=35),aes(xintercept=groupmean.dpi, color="red"),
             linetype="dotted", alpha=1, size=0.75, show.legend=FALSE)+
  facet_grid(c("temp", "dpi"))+
  theme_bw()


animated_plot <- ggplot(ti.mg, aes(x = total_eye_score, fill = fct_rev(temp))) +
  geom_histogram(position = "identity", alpha = 0.75, binwidth = 0.5, color="black") +
  transition_states(dpi, transition_length = 6, state_length = 12) +
  geom_vline(xintercept=0, alpha=0.5)+
  labs(title= "Day Post Inoculation: {closest_state}", fill = "Temperature", x= "Total Eye Score", y= "Count")+
  scale_fill_manual(values = c("Warm" = "gray15", "Cold" = "white"))+
  facet_wrap(~fct_rev(temp)~treatment)+
  enter_fade() +
  exit_fade()

#animate(animated_plot, renderer = gifski_renderer())
#anim_save("/Users/jesse/Documents/Virginia Tech/Research/Temp + Immunity 2022/Animated Plots/eye_score_hist.gif",
#          animation = last_animation(), renderer = gifski_renderer(), format = "gif")

#histograms of treatment groups across all sample days
ggplot(ti.mg, aes(x=total_eye_score, fill=temp))+
  geom_density()+
  #geom_histogram(binwidth = 0.5, position="identity", color="black", alpha=1)+
  #scale_fill_manual(values=c("lightskyblue", "blue", "lightpink2", "red"))+
  scale_fill_manual(values=c("white", "gray35"), labels = c("Cold", "Warm"))+
  geom_vline(data=ti.mg,aes(xintercept=groupmean, color="red"), linetype="dotted", alpha=1, size=0.75, show.legend=FALSE)+
  facet_grid(c("temp"))

#histograms of treatment groups across all sample days
ggplot(ti.mg %>% filter(dpi != -12 & dpi != 3), aes(x=total_eye_score, fill=temp))+
  #geom_density()+
  geom_histogram(binwidth = 0.75, position="identity", color="black", alpha=1)+
  #scale_fill_manual(values=c("lightskyblue", "blue", "lightpink2", "red"))+
  scale_fill_manual(values=c("white", "gray35"), labels = c("Cold", "Warm"))+
  geom_vline(data=ti.mg%>% filter(dpi != -12 & dpi != 3),aes(xintercept=groupmean.dpi, color="red"), linetype="dotted", alpha=1, size=0.75, show.legend=FALSE)+
  facet_grid(c("temp", "dpi"))+
  labs(y="Count", x= "Total Eye Score", fill = "Temperature")

#split violin plot
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
ggplot(ti.mg, aes(x=treatment, y=total_eye_score, fill = fct_rev(temp)))+ 
  geom_split_violin()+
  geom_point(
    aes(group = fct_rev(temp)),  # Group by temp
    position = position_dodge(width = 0.04),  # Adjust width as needed
    size = 5,
    alpha=0.4,
    shape = "-")+
  geom_count(
    aes(group = fct_rev(temp), color=temp),
    position = position_dodge(width = 0.04),
    shape = "-",
    show.legend = TRUE) +
  stat_summary(aes(group = fct_rev(temp)),  # Group by temp
               fun = mean, geom = "point", size = 10, shape = "-",
               position = position_dodge(width = 0.06))+  # Adjust width as needed
  scale_fill_manual(values = c("Cold" = "gray85", "Warm" = "gray35"))+
  ylim(-0, 5)+
  theme_minimal()+
  labs(title = "Total Eye Score", x = "Infection",
       y = "Total Eye Score", fill = "Temperature", color= "Temperature")

ti %>%
  dplyr::select(quantity, bird_ID, band_number, dpi, treatment, temp)%>%
  filter(treatment=="Control" & quantity > 50)



g.e.change <- ggplot(data=ti.mg, aes(x=as.factor(dpi), y=total_eye_score), color = temp, groups = band_number)+
  geom_jitter(aes(color=temp),size=0.4, width=0.1)+
  geom_smooth(aes(group = band_number, color = temp), method="loess", se=FALSE, size=0.2, span=0.5, linetype="solid", alpha=0.25)+
  stat_summary(aes(group=temp, color=temp), fun=mean, geom="line", linetype="solid", alpha=1, size=1)+
  stat_summary(aes(color=temp),fun.y=mean,
               #fun.min = function(x) mean(x)-sd(x),
               #fun.max = function(x) mean(x)+sd(x),
               #fun.min = function(x) mean(x)-(sd(x)/mean(x)), #CV
               #fun.max = function(x) mean(x)+(sd(x)/mean(x)), #CV
               geom= "pointrange", size=1.5, shape="+")+
  labs(y = "Total Eye Score", x="Days Post Infection", shape= "Treatment Group", color="Treatment Group")+
  scale_color_manual(values=c("blue","red"))+
  facet_wrap(~temp)+
  theme_bw()
g.e.change

#for poster
g.change <- ggplot(data=ti.mg %>% filter(dpi > 0), aes(x=dpi, y=total_eye_score), color = temp, groups = band_number)+
  geom_point(aes(color=temp),size=3)+
  geom_smooth(aes(group = band_number, color = temp), method="loess", se=FALSE, size=1, span=0.4, linetype="solid", alpha=0.25)+
  stat_summary(aes(group=temp, color=temp), fun=mean, geom="line", linetype="solid", alpha=1, size=2)+
  stat_summary(aes(group=temp), fun=mean, geom="line", linetype="solid", color = "black", alpha=1, size=0.5)+
  labs(y = "Conjunictivitis Severity Score\n", x="\nDays Post Inoculation", shape= "Treatment Group", color="Treatment Group")+
  scale_color_manual(values=c("blue","red"))+
  ylim(c(0,5))+
  facet_wrap(~temp)
g.change

#Simpson's Index
#treat each eye score as a different species
#Diversity index value increases when the number of types increases and the evenness increases
table(ti.mg$total_eye_score, ti.mg$dpi)
table(ti.mg$temp, ti.mg$total_eye_score, ti.mg$dpi)

install.packages("abdiv")
library(abdiv)

#cold first
cd7 <- c(11,2,1,0,0,0,0,0,0)
cd9 <- c(10,0,1,1,0,1,0,0,1)
cd14 <-c(8,1,0,2,2,0,0,0,1)
cd18 <- c(8,3,0,1,0,1,0,1,0)
cd21 <- c(10,3,0,0,1,0,0,0,0)
cd24 <- c(9,4,0,0,0,1,0,0,0)
cd28 <- c(11,1,1,1,0,0,0,0,0)
cd35 <- c(13,1,0,0,0,0,0,0,0)

wd7 <- c(8,5,1,0,0,0,0,0,0)
wd9 <- c(6,1,2,2,1,2,0,0,0)
wd14 <- c(7,2,0,3,0,0,2,0,0)
wd18 <- c(8,0,0,3,1,0,2,0,0)
wd21 <- c(8,3,0,0,1,1,1,0,0)
wd24 <- c(8,3,1,1,0,1,0,0,0)
wd28 <- c(10,3,0,0,1,0,0,0,0)
wd35 <- c(14,0,0,0,0,0,0,0,0)

#simpson
c7<- simpson(cd7)
c9 <-simpson(cd9)
c14 <-simpson(cd14)
c18 <-simpson(cd18)
c21<- simpson(cd21)
c28<- simpson(cd28)
c35<- simpson(cd35)

w7<- simpson(wd7)
w9 <-simpson(wd9)
w14 <-simpson(wd14)
w18 <-simpson(wd18)
w21<- simpson(wd21)
w28<- simpson(wd28)
w35<- simpson(wd35)

#simpsons for all observations over all days
table(ti.mg$total_eye_score, ti.mg$temp)
c.all <- c(108,15,3,5,3,3,0,1,2)
w.all <- c(97,17,4,9,4,4,5,0,0)
ct<- simpson(c.all)
wt <- simpson(w.all)

#create data frame
cold <- c(c7, c9, c14, c18, c21, c28, c35)
warm <- c(w7, w9, w14, w18, w21, w28, w35)
dpi <- c(7, 9, 14, 18, 21, 28, 35)
cold.all <- c(ct, ct, ct, ct, ct, ct, ct)
warm.all <- c(wt, wt, wt, wt, wt, wt, wt)
all
# Create the 'simp' data frame
simp <- data.frame(cold, warm, dpi, cold.all, warm.all)

simp$cold <- 1/simp$cold
simp$warm <- 1/simp$warm
simp$cold.all <- 1/simp$cold.all
simp$warm.all <- 1/simp$warm.all
#show simpson's diversity index
ggplot(simp %>% filter(dpi != 35), aes(x=dpi))+
  geom_point(aes(x=dpi, y=cold), color="blue", size=3)+
  geom_line(aes(x=dpi, y=cold), color="blue")+
  geom_point(aes(x=dpi, y=warm), color="red", size=3)+
  geom_line(aes(x=dpi, y=warm), color="red")+
  geom_path(aes(x=dpi, y=cold.all), color="blue")+
  geom_path(aes(x=dpi, y=warm.all), color="red")+
  labs(x="Days Post Infection", y= "Simpson's Diversity Index: Eye Score", color="Temperature")+
  theme_minimal()

#higher values indicate lower diversity/heterogeneity <- take the inverse to show higher diversity
