##Haptoglobin##
rm(list=ls())

library(ggplot2)
library(tidyverse)
install.packages(glmmTMB)
library(glmmTMB)

setwd('/Users/jesse/Documents/Virginia Tech/Research/Temp + Immunity 2022/RAW DATA/')

ti <- read.csv("ti_merged_data.csv")

ti$hapto_ext
ti$hapto <- as.numeric(ti$hapto_ext)
ti$hapto_int <- as.numeric(ti$hapto_int)
ti$red <- as.numeric(ti$red)
ti$hapto_int <- as.numeric(ti$hapto_int)


ti.h<-ti %>%
  filter(hapto_plate == 2 | hapto_plate == 3) %>%
  drop_na(hapto_ext, treatment, hapto_plate)

ti.h <- ti.h %>%
  select(band_number, dpi, bird_ID, hapto_int, hapto_abs, hapto.oor, red, hapto_plate, hapto_cv, hapto_ext, inf_temp,
         groups, temp, treatment, hapto_plate, sex)
ti.h%>%
  group_by(dpi, groups)%>%
  summarize(count = n())

#Group counts
#dpi groups        count  
#<int> <chr>         <int>
#1   -12 Cold Control      9
#2   -12 Cold Infected    12
#3   -12 Warm Control      7
#4   -12 Warm Infected    11
#5     2 Cold Control      8
#6     2 Cold Infected    11
#7     2 Warm Control      8
#8     2 Warm Infected    12

ti.h$hapto_ext_log <- log10(ti.h$hapto_ext)
hist(ti.h$hapto_ext_log)

#compare plates
ggplot(data=ti.h, aes(x=log10(hapto_int), y=log10(hapto_ext), color=as.factor(hapto_plate)))+
  geom_point()

#compare plates raw values
ggplot(data=ti.h, aes(x=hapto_int, y=hapto_ext, color=as.factor(hapto_plate)))+
  geom_point()

#compare plates int to od
ggplot(data=ti.h, aes(x=hapto_int, y=hapto_abs, color=as.factor(hapto_plate)))+
  geom_point()

#compare plates ext to od
ggplot(data=ti.h, aes(x=hapto_ext, y=hapto_abs, color=as.factor(hapto_plate)))+
  geom_point()

ti.h2 <- ti.h %>%
  filter(dpi == 2)

hist(ti.h2$hapto_ext)
#does the assay plate predict haptoglobin concentration?
lm1 <- glm(hapto_ext~hapto_plate + (1|band_number), data=ti.h, family=Gamma())
summary(lm1)
plot(allEffects(lm1))
hist(resid(lm1))

#plate does not affect hapto_ext

#On day 2 post-infection, do treatment + temp predict haptoglobin levels?
lm2 <- glm(hapto_ext~treatment + temp + red , data=ti.h2, family=Gamma())
summary(lm2)
plot(allEffects(lm2))
hist(resid(lm2))

lm2.5 <- glm(hapto_ext~treatment + temp + red, data=ti.h2, family=Gamma())
summary(lm2.5)
hist(resid(lm2.5))

lm2.75 <- glm(hapto_ext~treatment * temp + red, data=ti.h2, family=Gamma())
summary(lm2.75)

#the residuals of lm2 are more normally distributed than lm2.5 so I think I should use lm2?
#model selection
a1 <- glm(hapto_ext~treatment + temp , data=ti.h2, family=Gamma())
a2 <- glm(hapto_ext~treatment, data=ti.h2, family=Gamma())
a3 <- glm(hapto_ext~temp, data=ti.h2, family=Gamma())
a4 <- glm(hapto_ext~treatment * temp, data=ti.h2, family=Gamma())
a5 <- glm(hapto_ext~treatment + temp + red, data=ti.h2, family=Gamma())
a6 <- glm(hapto_ext~treatment + red, data=ti.h2, family=Gamma())
a7 <- glm(hapto_ext~temp + red, data=ti.h2, family=Gamma())
a8 <- glm(hapto_ext~treatment * temp + red, data=ti.h2, family=Gamma())
a9 <- glm(hapto_ext~treatment + temp + red + sex , data=ti.h2, family=Gamma())
a10 <- glm(hapto_ext~1, data=ti.h2, family = Gamma())

#AICc 
aictab(cand.set=list(a1, a2, a3, a4, a5, a6, a7, a8, a9, a10),
       modnames=c("a1","a2", "a3", "a4", "a5", "a6", "a7", "a8", "a9", "a10"))

#Model selection based on AICc:

#  K   AICc Delta_AICc AICcWt Cum.Wt    LL
#a6  4 -15.24       0.00   0.47   0.47 12.21
#a8  6 -13.25       2.00   0.17   0.65 13.94
#a9  6 -13.16       2.08   0.17   0.81 13.89
#a5  5 -12.96       2.28   0.15   0.96 12.39
#a7  4 -10.19       5.06   0.04   1.00  9.68
#a2  3  10.15      25.40   0.00   1.00 -1.73
#a1  4  10.32      25.56   0.00   1.00 -0.57
#a4  5  12.46      27.71   0.00   1.00 -0.32
#a10 2  18.62      33.87   0.00   1.00 -7.15
#a3  3  19.03      34.27   0.00   1.00 -6.17

#a6 is best model but keep in temp so use a8
summary(a6)
plot(allEffects(a6))

summary(a8)
plot(allEffects(a8))

summary(a9)
plot(allEffects(a9))

#are there differences in hapto_ext between all treatment groups?
g.hp <- ggplot(data=ti.h2, aes(x=groups, y=hapto_ext, shape=groups))+
  geom_jitter(width = .1)+
  stat_summary(fun.y=mean,
               fun.min = function(y) mean(y)-sd(y),
               fun.max = function(y) mean(y)+sd(y),
               geom= "errorbar", size=0.5, width=0.1, alpha=0.5)+ #error bars
  stat_summary(aes(group=groups), fun=mean, alpha=1, size=2, shape="-")+
  geom_hline(aes(yintercept=0), linetype="solid", alpha=0.75)+
  labs(x="Treatment", y="Haptoglobin Concentration Day 2", shape="Treatment")+
  scale_shape_manual(values=c(1, 2, 3, 4))+
  theme_minimal()
g.hp

####Model Predictions####
dat.new=expand.grid(temp=unique(ti.h$temp),
                    treatment=unique(ti.h$treatment),
                    red=unique(ti.h$red),
                    band_number=unique(ti.h$band_number))#new grid to put predictions into
dat.new$yhat = predict(a5, type="response", newdata=dat.new, re.form=NA) #predicted values based off a5 model
head(dat.new)
dat.new$inf_temp <-  paste(dat.new$treatment, dat.new$temp, sep = "_")

#plot predicted values over raw data
hapto.pred <- ggplot(data=ti.h2, aes(x=inf_temp, y=hapto_ext, shape=temp))+
  geom_jitter(size=2, width=0.1)+
  #geom_line(data=dat.new, aes(x=inf_temp, y=yhat, linetype=treatment, group=treatment), size=1)+ #model predictions
  labs(x="Infection Treatment", y="Haptoglobin Concentration", linetype="Temperature", shape="Temperature")+
  scale_shape_manual(values=c(0, 15),  labels= c("Cold", "Warm"))+
  scale_linetype_manual(values = c("dotted", "solid"), labels= c("Cold", "Warm"))+
  #scale_colour_manual(values=c("cadetblue4", "blue", "hotpink2","red"))
  theme_minimal()
hapto.pred

#are there differences in hapto_ext between all treatment groups?
g.mg.hp <- ggplot(data=ti.h2, aes(x=treatment, y=hapto_ext, shape=treatment))+
  geom_jitter(width = .1)+
  stat_summary(fun.y=mean,
               fun.min = function(y) mean(y)-sd(y),
               fun.max = function(y) mean(y)+sd(y),
               geom= "errorbar", size=0.5, width=0.1, alpha=0.5)+ #error bars
  stat_summary(aes(group=treatment), fun=mean, alpha=1, size=2, shape="-")+
  geom_hline(aes(yintercept=0), linetype="solid", alpha=0.75)+
  labs(x="Treatment", y="Haptoglobin Concentration Day 2", shape="Treatment")+
  scale_shape_manual(values=c(16, 1), labels= c("Infected", "Control"))+
  scale_x_discrete(labels = c("Infected", "Control")) +
  theme_minimal()
g.mg.hp

#model selection
b1 <- glm(hapto_ext~treatment + temp + (1|band_number) , data=ti.h, family=Gamma())
b2 <- glm(hapto_ext~treatment + (1|band_number), data=ti.h, family=Gamma())
b3 <- glm(hapto_ext~temp + (1|band_number), data=ti.h, family=Gamma())
b4 <- glm(hapto_ext~treatment * temp + (1|band_number), data=ti.h, family=Gamma())
b5 <- glm(hapto_ext~treatment + temp + red + (1|band_number), data=ti.h, family=Gamma())
b6 <- glm(hapto_ext~treatment + red + (1|band_number), data=ti.h, family=Gamma())
b7 <- glm(hapto_ext~temp + red + (1|band_number), data=ti.h, family=Gamma())
b8 <- glm(hapto_ext~treatment * temp + red + (1|band_number), data=ti.h, family=Gamma())
b9 <- glm(hapto_ext~treatment + temp + red + sex + (1|band_number), data=ti.h, family=Gamma())
b10 <- glm(hapto_ext~1 + (1|band_number), data=ti.h, family = Gamma())

#AICc 
aictab(cand.set=list(b1, b2, b3, b4, b5, b6, b7, b8, b9, b10),
       modnames=c("b1","b2", "b3", "b4", "b5", "b6", "b7", "b8", "b9", "b10"))

#Model selection based on AICc:
  
#  K   AICc Delta_AICc AICcWt Cum.Wt    LL
#b7  4 -14.58       0.00   0.30   0.30 11.57
#b6  4 -14.57       0.02   0.30   0.60 11.56
#b9  6 -13.56       1.02   0.18   0.78 13.37
#b5  5 -12.61       1.97   0.11   0.89 11.72
#b8  6 -12.55       2.03   0.11   1.00 12.87
#b1  4  17.97      32.55   0.00   1.00 -4.71
#b3  3  19.65      34.23   0.00   1.00 -6.66
#b4  5  19.94      34.53   0.00   1.00 -4.55
#b2  3  21.56      36.14   0.00   1.00 -7.62
#b10 2  23.00      37.58   0.00   1.00 -9.42

summary(b7)
summary(b6) #neither treatment nor temperature predict haptoglobin concentration across days
summary(b9)

#haptoglobin concentrations day -12 and 2 post-infeciton
ggplot(data=ti.h, aes(x=treatment, y=hapto_ext, shape =treatment))+
  geom_jitter(width=0.1)+
  stat_summary(aes(group=treatment), fun=mean, alpha=1, size=2, shape="-")+
  stat_summary(fun.y=mean,
               fun.min = function(y) mean(y)-sd(y),
               fun.max = function(y) mean(y)+sd(y),
               geom= "errorbar", size=0.25, width=0.1, alpha=0.75)+ #error bars
  labs(y="log10(1+[Hp])", x="Treatment Group", shape = "Treatment Group")+
  scale_shape_manual(values=c(16, 1), labels= c("Infected", "Control"))+
  scale_x_discrete(labels = c("Infected", "Control")) +
  theme_minimal()+
  facet_wrap(~dpi)
       
ti.hod <- ti %>%
  filter(dpi == -12 | dpi == 2)

#model selection
bo1 <- glm(hapto_abs~treatment + temp + (1|band_number) , data=ti.hod, family=Gamma())
bo2 <- glm(hapto_abs~treatment + (1|band_number), data=ti.hod, family=Gamma())
bo3 <- glm(hapto_abs~temp + (1|band_number), data=ti.hod, family=Gamma())
bo4 <- glm(hapto_abs~treatment * temp + (1|band_number), data=ti.hod, family=Gamma())
bo5 <- glm(hapto_abs~treatment + temp + red + (1|band_number), data=ti.hod, family=Gamma())
bo6 <- glm(hapto_abs~treatment + red + (1|band_number), data=ti.hod, family=Gamma())
bo7 <- glm(hapto_abs~temp + red + (1|band_number), data=ti.hod, family=Gamma())
bo8 <- glm(hapto_abs~treatment * temp + red + (1|band_number), data=ti.hod, family=Gamma())
bo9 <- glm(hapto_abs~treatment + temp + red + sex + (1|band_number), data=ti.hod, family=Gamma())
bo10 <- glm(hapto_abs~1 + (1|band_number), data=ti.hod, family = Gamma())

#AICc 
aictab(cand.set=list(bo1, bo2, bo3, bo4, bo5, bo6, bo7, bo8, bo9, bo10),
       modnames=c("bo1","bo2", "bo3", "bo4", "bo5", "bo6", "bo7", "bo8", "bo9", "bo10"))

summary(bo7)
summary(bo5)
#haptoglobin absorbance day -12 and 2 post-infeciton
ggplot(data=ti.hod, aes(x=treatment, y=hapto_abs, shape =treatment))+
  geom_jitter(width=0.1)+
  stat_summary(aes(group=treatment), fun=mean, alpha=1, size=2, shape="-")+
  stat_summary(fun.y=mean,
               fun.min = function(y) mean(y)-sd(y),
               fun.max = function(y) mean(y)+sd(y),
               geom= "errorbar", size=0.25, width=0.1, alpha=0.75)+ #error bars
  labs(y="Haptoglobin RAW OD", x="Treatment Group", shape = "Treatment Group")+
  scale_shape_manual(values=c(16, 1), labels= c("Infected", "Control"))+
  scale_x_discrete(labels = c("Infected", "Control")) +
  theme_minimal()+
  facet_wrap(~dpi)



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
ggplot(ti.h, aes(x=treatment, y=hapto_ext_log, fill = fct_rev(temp)))+ 
  geom_split_violin()+
  geom_point(
    aes(group = fct_rev(temp), color=temp),  # Group by temp
    position = position_dodge(width = 0.04),  # Adjust width as needed
    size = 5,
    alpha=1,
    shape = "-") +
  stat_summary(aes(group = fct_rev(temp), color=temp),  # Group by temp
               fun = mean, geom = "point", size = 10, shape = "-",
               position = position_dodge(width = 0.15))+  # Adjust width as needed
  scale_color_manual(values = c("ST" = "black", "TN" = "black")) +  # Set color for each level
  scale_fill_manual(values = c("ST" = "white", "TN" = "gray"))+
  theme_minimal()+
  labs(title = "Haptoglobin Concentration", x = "Infection", y = "log10([Hp])", fill = "Temperature", color= "Temperature")+
  facet_wrap(~dpi)

#Change in haptoglobin over all dpi
g.hp.change <- ggplot(data=ti.h, aes(x=as.factor(dpi), y=hapto_ext), color = groups, groups = band_number)+
  geom_point(aes(shape=groups, color=groups),size=1)+
  geom_line(aes(group = band_number, color = groups), size=0.5, linetype="solid", alpha=0.25)+
  stat_summary(aes(group=groups, color=groups), fun=mean, geom="line", linetype="solid", alpha=1, size=1)+
  stat_summary(aes(color=groups),fun.y=mean,
               fun.min = function(x) mean(x)-sd(x),
               fun.max = function(x) mean(x)+sd(x),
               #fun.min = function(x) mean(x)-(sd(x)/mean(x)), #CV
               #fun.max = function(x) mean(x)+(sd(x)/mean(x)), #CV
               geom= "pointrange", size=1.5, shape="+")+
  geom_hline(aes(yintercept=0), linetype='solid')+
  labs(y = "Haptoglobin Concentration", x="Days Post Infection", shape= "Treatment Group", color="Treatment Group")+
  scale_color_manual(values=c("lightskyblue", "blue", "lightpink2", "red"))+
  scale_shape_manual(values=c(1, 2, 3, 4))+
  theme_minimal()+
facet_wrap(~groups)
g.hp.change
