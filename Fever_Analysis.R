#T+I 22 Fever Analysis
rm(list=ls())


####read in + format data####
setwd('/Users/jesse/Documents/Virginia Tech/Research/Temp + Immunity 2022')
library(ggplot2)
library(tidyverse)
library(lme4)
#install.packages("glmmTMB")
#install.packages('TMB', type = 'source')
library(glmmTMB)
library(effects)

ti <- read.csv("RAW DATA/ti_merged_data.csv")

#poster theme
#theme_set(
#  theme_bw() +
#    theme(
#      text = element_text(size = 50),
#      axis.title.y = element_text(size = 100, face = "bold"),
#      axis.text.y = element_text(size = 75,  face = "bold"),
#      axis.title.x = element_text(size = 100, face = "bold"),
#      axis.text.x = element_text(size = 75,  face = "bold"),
#      strip.text = element_text(size = 100,  face = "bold")
#    )
#)

theme_set(theme_bw())
#Add column for total eye score by adding l and r eye score
ti <- ti %>%
  group_by(bird_ID) %>%
  mutate(total_eye_score=(l_eye_score+r_eye_score))

#make total eye score an integer by multiplying by two for models
ti$tes <- ti$total_eye_score*2

ti$threshold_cutoff = 100
ti$sympt_cutoff = 0.1

ti <- ti %>%
  mutate(infected = ifelse(quantity>threshold_cutoff, 1, 0))  #generate infection data; if path load > 50 copies = infected

#new column with whether the bird is ever infected
ti <- ti %>%
  group_by(band_number) %>%
  mutate(ever_infected = ifelse(any(coalesce(infected, 0) == 1), 1, 0)) %>%
  ungroup()
ti$ever_infected
#new column with whether the bird is diseased on a particular day
ti <- ti %>%
  mutate(diseased= ifelse(total_eye_score>sympt_cutoff, 1, 0))
ti$diseased

#new column with whether the bird is ever diseased
ti <- ti %>%
  group_by(band_number) %>%
  mutate(ever_diseased = ifelse(any(coalesce(diseased, 0) == 1), 1, 0)) %>%
  ungroup()

#table w/ number of infected birds 
t1<-ti%>%  
  group_by(groups, dpi, diseased, ever_diseased)%>%
  summarize(count = n())
t1
####Format Fever Data####
#Add column for fever score by finding average score
ti <- ti %>%
  group_by(bird_ID) %>%
  mutate(fever_score=((((l_min + l_max)/2) + (r_min + r_max)/2)/2))
hist(ti$fever_score)
table(ti$dpi)
 ti.tr <- ti %>%
  drop_na (band_number, fever_score, dpi, treatment, temp, groups, sex)

 #fever change from baseline
ti.f <- ti.tr %>%
  group_by(band_number) %>%
  mutate(fever_change = fever_score - first(fever_score))

table(ti.f$dpi)

summary(ti.f$fever_change)
summary(ti.f$fever_score)

ggplot(ti.f, aes(x=fever_score, y=l_max + r_max, color=groups))+
  geom_point()+
  facet_wrap(~dpi)


#fever change from previous score
ti.f <- ti.f %>%
  group_by(band_number) %>%
  mutate(fever_diff = fever_score - lag(fever_score))

print(ti.f)

#only days w/fever
ti.f <- ti.f%>%
  filter(dpi %in% c(0, 3, 7, 14, 18, 24, 28, 35))

p <- ggplot(ti.f, aes(x=fever_score, fill = fct_rev(temp)))+
         geom_histogram(position = "identity", alpha=0.75)+
  facet_wrap(~dpi)
p

ggplot(ti.f, aes(x=fever_change, fill = fct_rev(temp)))+
  geom_histogram(position = "identity", alpha=0.75, binwidth = 1, color="black")+
  geom_vline(aes(xintercept=0))+
  facet_wrap(~ever_diseased~temp~dpi, ncol=7)

####Animated Histograms####
#install.packages("gganimate")
library(gganimate)
#install.packages("gifski")
library(gifski)
#install.packages("magick")
library(magick)

temp_names <- c(
  'Cold' = "Cold",
  'Warm' = "Warm"
)

dis_names <- c(
  "1" = "Eye Score",
  "0" = "No Eye Score"
)
treat_names <- c(
  "Infected" = "Infected",
  "Control" = "Control"
)

ggplot(ti.f, aes(x = fever_change, fill = fct_rev(temp))) +
  geom_histogram(position = "identity", alpha = 0.75, binwidth = 0.75, color="black") +
  #transition_states(dpi, transition_length = 2, state_length = 1) +
  labs(title= "Day Post Inoculation: {closest_state}", fill = "Temperature")+
  scale_fill_manual(values = c("TN" = "red", "ST" = "blue"), labels= c("Infected", "Control"))+
  facet_wrap(~fct_rev(temp)~treatment~ever_diseased, labeller=as_labeller(c(temp_names, dis_names, treat_names)))

animated_plot <- ggplot(ti.f, aes(x = fever_change, fill = fct_rev(temp))) +
  geom_histogram(position = "identity", alpha = 0.75, binwidth = 0.75, color="black") +
  transition_states(dpi, transition_length = 2, state_length = 2) +
  geom_vline(xintercept=0, alpha=0.5)+
  labs(title= "Day Post Inoculation: {closest_state}", fill = "Temperature")+
  scale_fill_manual(values = c("Warm" = "red", "Cold" = "blue"))+
  facet_wrap(~fct_rev(temp)~treatment~ever_diseased, labeller=as_labeller(c(temp_names, dis_names, treat_names)))+
  enter_fade() +
  exit_fade()

#animate(animated_plot, renderer = gifski_renderer())
#anim_save("/Users/jesse/Documents/Virginia Tech/Research/Temp + Immunity 2022/Animated Plots/fever_score_hist.gif",
#          animation = last_animation(), renderer = gifski_renderer(), format = "gif")

####fever change vs fever difference####
#which should I use? Fever change or fever difference?
check<-ti.f %>%
  dplyr::select(fever_score, fever_change, bird_ID, band_number, fever_diff)

#compare a subset of each
ggplot(ti.f %>% filter(band_number>2900), aes(x=dpi, y=fever_diff, groups=as.factor(band_number), color=as.factor(band_number))) +
  geom_point()+
  geom_line()+
  geom_line(aes(x=dpi, y=fever_change, groups=as.factor(band_number), color=as.factor(band_number)), linetype="dashed")+
  geom_point(aes(x=dpi, y=fever_change, groups=as.factor(band_number), color=as.factor(band_number)), shape=1)+
  facet_wrap(~temp)

#comapre to fever score
ggplot(ti.f %>% filter(band_number>2900), aes(x=dpi, y=fever_score, color=as.factor(band_number), groups=as.factor(band_number)))+
  geom_line()+
  geom_point()+
  facet_wrap(~temp)

#what's going on on dpi 24 in the warm room?
ggplot(ti.f, aes(x=dpi, y=fever_diff, color=treatment, groups=as.factor(band_number))) +
  geom_point()+
  geom_line()+
  scale_x_continuous(breaks=unique(ti$dpi), labels=paste0(unique(ti$dpi)))+
  facet_wrap(~temp)

#remove dpi 24 from analysis
ti.f <- ti.f %>%
  filter(dpi != 24)

t1 <- ti.f%>% 
  group_by(groups, dpi)%>%
  summarize(count = n(),
            fever = mean(fever_score))

t1
unique_counts <- ti.f %>%  
  #filter(!is.na(elisa_od)) %>%  
  group_by(groups) %>%
  group_by(dpi) %>%
  summarise(unique_band_num = n_distinct(band_number))

unique_counts


# Get all possible band_number values
all_band_numbers <- ti.f %>% pull(band_number) %>% unique()

# Create a data frame with all combinations of dpi and band_number
all_combinations <- expand.grid(dpi = unique(ti.f$dpi), band_number = all_band_numbers)

# Merge with the actual data to identify missing band_numbers
missing_band_numbers <- all_combinations %>%
  anti_join(ti.f, by = c("dpi", "band_number"))

# Display the missing band_numbers
print(missing_band_numbers)

#missing 2684 after it died on dpi 28
#2648 dpi 35 - pictures not taken
#2632 dpi 28 - pictures not taken

unique_group_counts <- ti.f %>%
  filter(!is.na(fever_score)) %>%
  count(dpi, band_number) %>%
  split(.$n)
#view(unique_group_counts)

####Models####
##Linear Model All Days##
hist(ti.f$fever_change)
#looks normally distributed - use gaussian distribution
lm1 <- glmmTMB(fever_change ~temp+treatment + (1|band_number),data=ti.f)
summary(lm1)
hist(resid(lm1))
plot(allEffects(lm1))

Anova(lm1, "III")
table(ti.f$fever_change)
lm2 <- glmmTMB(fever_change ~temp*treatment +  (1|band_number),data=ti.f)
hist(resid(lm2))
summary(lm2)

#Conditional model:
#                       Estimate Std. Error z value Pr(>|z|)    
#(Intercept)            1.0198     0.2250   4.533 5.83e-06 ***
#tempTN                -0.2611     0.2197  -1.188 0.234638    
#treatmentSham         -0.8395     0.2423  -3.465 0.000531 ***
#tempTN:treatmentSham  -0.3144     0.3423  -0.919 0.358313    

plot(allEffects(lm2))
Anova(lm2, "III")


range(ti.f$fever_change)
shapiro.test(resid(lm2))
#Treatment, but not temperature or their interaction had a significant effect on fever score across all dpi

#Model Selection <- this isn't working my models are too complex
e1 <- glmmTMB(fever_change ~temp+treatment + (1|band_number),data=ti.f)
e2 <- glmmTMB(fever_change ~temp+treatment + sex + (1|band_number),data=ti.f)
e3 <- glmmTMB(fever_change ~temp*treatment + (1|band_number),data=ti.f) 
e4 <- glmmTMB(fever_change ~temp*treatment + sex + (1|band_number),data=ti.f) 
e5 <- glmmTMB(fever_change ~temp + (1|band_number),data=ti.f)
e6 <- glmmTMB(fever_change ~treatment + (1|band_number),data=ti.f) 
e7 <- glmmTMB(fever_change ~1 + (1|band_number),data=ti.f)

aictab(cand.set=list(e1, e2,e3,e4,e6,e7),
       modnames=c("e1", "e2","e3","e4", "e6","e7"))

summary(e1)
summary(e2)
summary(e3)

f1<- glmmTMB(fever_score ~temp+treatment + (1|band_number),data=ti.f) 
summary(f1)

temp_names <- c(
  'Cold' = "Cold",
  'Warm' = "Warm"
)

ever_inf_names <- c(
  '0' = "Never Infected",
  '1' = "Infected"
)
table(ti.f$groups)
group_names <- c(
  'Cold Control' = "Cold Control",
  'Cold Infected' = "Cold Infected",
  'Warm Control' = "Warm Control",
  'Warm Infected' = "Warm Infected"
)

ti.f%>%
  filter(fever_change > 5) %>%
 dplyr::select(dpi, band_number, fever_change, fever_score)



#graph showing fever change over the course of infection for poster
ggplot(data=ti.f %>% filter(treatment == "Infected"), aes(x=dpi, y=fever_change, color = temp))+
  geom_jitter(aes(shape=temp), size=1, width = 0.15, stroke=2)+
  #geom_line(aes(group=band_number, linetype=temp), method="loess", se=FALSE, size=1) +
  geom_smooth(aes(group=band_number, linetype=temp), method="loess", se=FALSE, span = 0.71, size=0.5)+ 
  stat_summary(aes(group=temp, linetype=temp), fun=mean, geom="line", alpha=1, size=2)+
  stat_summary(aes(group=temp, linetype=temp), fun=mean, geom="line", color="black", alpha=1, size=.5)+
  #stat_summary(aes(group=groups), fun=mean,
  #             fun.min = function(x) mean(x)-sd(x),
  #             fun.max = function(x) mean(x)+sd(x),
  #             geom= "errorbar", size=0.25, width=0.5)+
  scale_color_manual(values=c("blue", "red"))+
  #scale_x_continuous(breaks=unique(ti$dpi), labels=paste0(unique(ti$dpi)), limits = c(0,37))+
  scale_shape_manual(values=c(16, 16))+
  scale_linetype_manual(values=c("solid", "solid"))+
  geom_hline(yintercept = 0, alpha=0.25)+
  labs(x="\nDays Post Infection", y="Fever Change From Baseline\n", 
       linetype= "Temperature", shape= "Temperature", color="Temperature")+
  facet_grid(~temp, labeller = as_labeller(c(temp_names)), scales="free_x")

ggplot(data=ti.f, aes(x=dpi, y=fever_change, color = groups))+
  geom_jitter(aes(shape=treatment), size=1, width = 0.15)+
  #geom_path(aes(x = dpi, y = fever_change, group = groups), alpha = 0.3, size = 0.5) +  
  geom_line(aes(x=dpi, y=fever_change, linetype=treatment, group=(band_number)), alpha=0.5, size=0.5)+
  #geom_smooth(aes(group=band_number, linetype=treatment), method="loess", se=FALSE, size=0.2, span=0.65) +
  stat_summary(aes(group=groups, linetype=treatment), fun=mean, geom="line", alpha=1, size=1)+
  stat_summary(aes(group=groups, linetype=treatment), fun=mean, geom="line", color="black", alpha=1, size=.5)+
  stat_summary(aes(group=groups), fun=mean,
               fun.min = function(x) mean(x)-sd(x),
               fun.max = function(x) mean(x)+sd(x),
               geom= "errorbar", size=0.5, width=0.25)+
  scale_color_manual(values=c("cyan3", "blue", "violet", "red"))+
  scale_x_continuous(breaks=unique(ti$dpi), labels=paste0(unique(ti$dpi)), limits = c(0,37))+
  scale_shape_manual(values=c(16, 1))+
  scale_linetype_manual(values=c("dashed", "solid"))+
  labs(x="Days Post Infection", y="Fever Score", 
       linetype= "Treatment", shape= "Treatment", color="Temperature")+
  facet_wrap(~groups, labeller = as_labeller(c(group_names)), scales="fixed")+
  #facet_wrap(~groups~ever_diseased, labeller = as_labeller(c(group_names, dis_names)), scales="fixed")+
  theme_bw()

lmfs<- glmmTMB(fever_change~temp + (1|band_number), data=ti.f)

summary(lmfs)
plot(allEffects(lmfs))
#graph showing fever change over the course of infection
ggplot(data=ti.f %>% filter(dpi >0) , aes(x=groups, y=fever_change, shape = groups, color=groups))+
  geom_jitter(aes(shape=groups), size=1, width = 0.15)+
  stat_summary(aes(group=groups), fun=mean, shape="-", size=1)+
  stat_summary(aes(group=groups), fun=mean,
               fun.min = function(x) mean(x)-sd(x),
               fun.max = function(x) mean(x)+sd(x),
               geom= "errorbar", size=0.25,
               width=0.25)+
  scale_shape_manual(values=c(1, 2, 3, 4))+
  scale_color_manual(values=c("lightskyblue", "blue", "lightpink2", "red"))+
  #scale_x_continuous(breaks=unique(ti$dpi), labels=paste0(unique(ti$dpi)), limits = c(0,37))+
  geom_hline(yintercept = 0, alpha=0.25)+
  labs(x="Days Post Infection", y="Fever Change From Baseline)", shape= "Treatment")+
  #facet_wrap(~diseased, labeller = as_labeller(dis_names), scales="fixed")+
  theme_bw()+
  facet_wrap(~dpi)

#install.packages("transformr")
library(transformr)

f<- ggplot(data = ti.f %>% filter(ever_diseased == 1), aes(x = dpi, y = fever_change, color = as.factor(band_number))) +
  geom_jitter(
    aes(shape = treatment, alpha=cumsum(!is.na(fever_change))),
    size = 1, alpha = 0.5) +
  geom_line(aes(x = dpi, y = fever_change, group=as.factor(band_number)), alpha = 0.75, size = 0.5) +
  stat_summary(aes(group=groups, linetype=treatment), fun=mean, geom="line", alpha=1, size=0.25)+
  stat_summary(aes(group=groups, linetype=treatment), fun=mean, geom="line", color="black", alpha=1, size=0.25)+
  #scale_color_manual(values = c("lightskyblue", "blue", "lightpink2", "red")) +
  scale_x_continuous(breaks = unique(ti.f$dpi), labels = paste0(unique(ti.f$dpi)), limits = c(0, 37)) +
  scale_shape_manual(values = c(16, 1)) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  geom_hline(yintercept = 0, alpha = 0.25) +
  scale_alpha_continuous(range = c(0.2, 1)) +
  labs(
    x = "Days Post Infection",
    y = "Fever Change (Day(x) Fever - Baseline Fever)",
    linetype = "Treatment",
    shape = "Treatment",
    color = "Temperature"
  ) +
  facet_wrap(~temp~ever_diseased, labeller = as_labeller(c(temp_names, dis_names)), scales = "fixed") +
  theme_bw() 

e<- ggplot(data = ti.f %>% filter(ever_diseased == 1), aes(x = dpi, y = total_eye_score, color = as.factor(band_number))) +
  geom_jitter(
    aes(shape = treatment, alpha=cumsum(!is.na(total_eye_score))),
    size = 1, alpha = 0.5) +
  geom_line(aes(x = dpi, y = total_eye_score, group=as.factor(band_number)), alpha = 0.75, size = 0.5) +
  stat_summary(aes(group=groups, linetype=treatment), fun=mean, geom="line", alpha=1, size=0.25)+
  stat_summary(aes(group=groups, linetype=treatment), fun=mean, geom="line", color="black", alpha=1, size=0.25)+
  #scale_color_manual(values = c("lightskyblue", "blue", "lightpink2", "red")) +
  scale_x_continuous(breaks = unique(ti.f$dpi), labels = paste0(unique(ti.f$dpi)), limits = c(0, 37)) +
  scale_shape_manual(values = c(16, 1)) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  geom_hline(yintercept = 0, alpha = 0.25) +
  scale_alpha_continuous(range = c(0.2, 1)) +
  labs(
    x = "Days Post Infection",
    y = "Total Eye Score)",
    linetype = "Treatment",
    shape = "Treatment",
    color = "Temperature"
  ) +
  facet_wrap(~temp~ever_diseased, labeller = as_labeller(c(temp_names, dis_names)), scales = "fixed") +
  theme_bw() 
e
library(patchwork)  

f + e

ggplot(data = ti.f, aes(x = dpi, y = fever_change, color = groups)) +
  geom_jitter(
    aes(shape = treatment, alpha = 1),
    size = 1, alpha = 0.5) +
  #geom_path(aes(group = groups), alpha = 0.75, size = 0.5) +
  stat_summary(aes(group = groups, linetype = treatment), fun = mean, geom = "line", alpha = 1, size = 1) +
  scale_color_manual(values=c("cyan3", "blue", "violet", "red"))+
  scale_x_continuous(breaks = unique(ti.f$dpi), labels = paste0(unique(ti.f$dpi)), limits = c(0, 37)) +
  scale_shape_manual(values = c(16, 1)) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  geom_hline(yintercept = 0, alpha = 0.25) +
  scale_alpha_continuous(range = c(0.2, 1)) +
  labs(x = "Days Post Infection", y = "Fever Change (Day(x) Fever - Baseline Fever)",
       linetype = "Treatment", shape = "Treatment", color = "Temperature") +
  facet_wrap(~treatment)+
  #facet_wrap(~ever_diseased, labeller = as_labeller(dis_names), scales = "fixed") +
  theme_bw() 
#install.packages("rust")
#library(rust)
animated_plot2 <- ggplot(data = ti.f, aes(x = dpi, y = fever_change, color = groups)) +
  stat_summary(data=ti.f,
    aes(group = groups, linetype = treatment),
    fun = mean,
    geom = "line",
    alpha = 1,
    size = 1) +
  geom_jitter(
    aes(shape = treatment, alpha = 1), size = 1, alpha = 0.5) +
  scale_color_manual(values = c("lightskyblue", "blue", "lightpink2", "red")) +
  scale_x_continuous(breaks = unique(ti.f$dpi), labels = paste0(unique(ti.f$dpi)), limits = c(0, 37)) +
  scale_shape_manual(values = c(16, 1)) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  geom_hline(yintercept = 0, alpha = 0.25) +
  labs(x = "Days Post Infection", y = "Fever Change (Day(x) Fever - Baseline Fever)",
       linetype = "Treatment", shape = "Treatment", color = "Temperature") +
  facet_wrap(~ever_diseased, labeller = as_labeller(dis_names), scales = "fixed") +
  theme_bw() +
  transition_states(dpi, transition_length = 2, state_length = 2) +
  shadow_mark(past = TRUE)

#animate(animated_plot2, renderer = gifski_renderer())


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
ggplot(ti.f, aes(x=0, y=fever_change, fill = temp))+ 
  geom_split_violin()+
  geom_point(
    aes(group = groups),  # Group by temp
    position = position_dodge(width = 0.04),  # Adjust width as needed
    size = 5,
    alpha=0.4,
    shape = "-") +
  stat_summary(aes(group = groups),  # Group by temp
               fun = mean, geom = "point", size = 7, shape = "-",
               position = position_dodge(width = 0.06))+  # Adjust width as needed
  #scale_color_manual(values = c("ST" = "blue", "TN" = "red")) +  # Set color for each level
  theme_minimal()+
  labs(title = "Fever Change from Baseline", x = "Probability", y = "Fever Change", fill = "Temperature", color= "Temperature")+
  facet_wrap(~treatment~dpi, ncol=7)+
  theme_bw()

anim_viol<-ggplot(ti.f, aes(x=0, y=fever_change, fill = groups))+ 
  geom_split_violin()+
  geom_point(
    aes(group = groups),  # Group by temp
    position = position_dodge(width = 0.04),  # Adjust width as needed
    size = 5,
    alpha=0.4,
    shape = "-") +
  stat_summary(aes(group = groups),  # Group by temp
               fun = mean, geom = "point", size = 7, shape = "-",
               position = position_dodge(width = 0.06))+  # Adjust width as needed
  #scale_color_manual(values = c("ST" = "blue", "TN" = "red")) +  # Set color for each level
  theme_minimal()+
  labs(title = "Fever: DPI {closest_state}", x = "Probability", y = "Fever Change", fill = "Temperature", color= "Temperature")+
  facet_wrap(~treatment~ever_diseased, labeller = as_labeller(c(treat_names,dis_names)))+
  transition_states(dpi, transition_length = 2, state_length = 2) +
  exit_fade()

#animate(anim_viol, renderer = gifski_renderer())
#anim_save("/Users/jesse/Documents/Virginia Tech/Research/Temp + Immunity 2022/Animated Plots/fever_score_viol.gif",
#          animation = last_animation(), renderer = gifski_renderer(), format = "gif")

####Variance####
#across all  days, which groups are most heterogeneous?
glm.all <- glmmTMB(fever_change~1 +(1|band_number), data=ti.f)
glm.tnmg <- glmmTMB(fever_change~1 +(1|band_number), data=subset(ti.f, groups == "Warm Infected"))
glm.tnc <- glmmTMB(fever_change~1 +(1|band_number), data=subset(ti.f, groups == "Warm Control"))
glm.stmg <- glmmTMB(fever_change~1 +(1|band_number), data=subset(ti.f, groups == "Cold Infected"))
glm.stc <- glmmTMB(fever_change~1 +(1|band_number), data=subset(ti.f, groups == "Cold Control"))

summary(glm.all); summary(glm.tnmg); summary(glm.tnc); summary(glm.stmg); summary(glm.stc)

#without dpi 24
#glm.all: V = 0.7357; AIC = 986.6; SD = 0.8577; LogLik = -490.3
#glm.tnmg: V = 0.5979; AIC = 347.2; SD = 0.7732; LogLik = -169.6
#glm.tnc:  V = 0.3482; AIC = 199.2; SD = 0.5901; LogLik = -95.6
#glm.stmg: V = 1.035; AIC = 392.7; SD = 1.017; LogLik = -192.3
#glm.stc: V = 0.3932; AIC = 181.1; SD = 0.6022; LogLik = -86.5 

#including dpi 24
#glm.all: V = 0.7096; AIC = 1312.7; SD = 0.8424; LogLik = -652.3
#glm.tnmg: V = 0.7322; AIC = 396.2; SD = 0.8557; LogLik = -194.1
#glm.tnc:  V = 0.1622; AIC = 230.1; SD = 0.4028; LogLik = -111.0
#glm.stmg: V = 1.2286; AIC = 438.2; SD = 1.1084; LogLik = -215.1
#glm.stc: V = 2.571e-02; AIC = 206.0; SD = 0.1603470; LogLik = -99.0 

AIC(glm.all, glm.tnmg, glm.tnc, glm.stmg, glm.stc)

#quick graph showing differences in variance by group including dpi 24
#var <- c(0.7096, 0.7322, 0.1622, 1.2286, 2.571e-02) #variance from table above
#ci <- c(0.8424, 0.8557, 0.4028, 1.1084, 0.1603470) #Standard Deviation from table above
#g <- c("All Groups", "Warm Infected", "Warm Control", "Cold Infected", "Cold Control") #Groups

#differences in variance by group excluding dpi 24
var <- c(0.7357, 0.5979, 0.3482, 1.035, 0.3932) #variance from table above
ci <- c(0.8577, 0.7732, 0.5901, 1.017, 0.6022) #Standard Deviation from table above
g <- c("All Groups", "Warm Infected", "Warm Control", "Cold Infected", "Cold Control") #Groups
temp <- c("all", "Warm", "Warm", "Cold", "Cold")

#new df with variance and error for error bars (min/max)
variability <- data.frame(g,var, ci, temp)
variability$min <- variability$var - variability$ci #new column for lowerbound error
variability$max <- variability$var + variability$ci #new column for upperbound error

g.var <- ggplot(variability %>% filter(temp != "all"), aes(x=g, y=var, shape = g))+
  #geom_col(color="black", size=0.5)+
  geom_point(size=4, alpha = 1)+
  #geom_errorbar(aes(ymin=var, ymax=max), width=0.5)+ #add errorbars representing +/- 1 SD
  labs(x="\nTemperature", y="Variance\n", shape="Temperature", color="Temperature")+
  scale_shape_manual(values = c(0, 1, 2, 3))+
  ylim(c(0,1.5))

g.var

#quantify w stats
glm1 <- glmmTMB(fever_change ~ (1|dpi) + (1|band_number), data=ti.f)

glm1

ti.f$resid <- resid(glm1)
hist(ti.f$resid)
lm1 <- glm(resid ~ temp*treatment, data=ti.f)
lm2 <- glm(resid ~ temp, data=ti.f)
anova(lm2)

int <- aov(resid ~ temp*treatment, data=ti.f)
summary(int)

ggplot(ti.f, aes(x=groups, y=abs(resid)))+
  geom_point()


####CVs####
#use fever_score instead of fever_change because it is just a measure of variability
cv.all <- ti.f %>% 
  group_by(groups)%>%
  summarise(bird_cv = sd(fever_score)/mean(fever_score), #this doesn't actually make sense
            bird_sd = sd(fever_score))

#generate error bar values
cv.all$max <- cv.all$bird_cv + cv.all$bird_sd
cv.all$min <- cv.all$bird_cv - cv.all$bird_sd

cv.all
#Graph of CVs
g.cv.all <- ggplot(cv.all, aes(x=groups, y=bird_cv, shape=groups))+
  geom_point(aes(shape=groups), color="black", size=2)+
  scale_shape_manual(values=c(0,1,2,3))+
  #geom_errorbar(aes(ymax=max, ymin=min), color="black", size=0.5, width=0.5)+
  labs(x="Treatment Groups", y="CV", shape="Treatment")+
  scale_fill_manual(values=c("gray85", "white", "gray65", "gray25"))+
  theme_minimal()

g.cv.all

#CV by day
cv.all.day <- ti.f %>%
  group_by(dpi, groups, treatment, temp) %>%
  summarise(
    bird_cv = sd(fever_score) / mean(fever_score),
    bird_sd = sd(fever_score),
    max = bird_cv + bird_sd,
    min = bird_cv - bird_sd
  ) %>%
  ungroup()  # Remove grouping for further operations

cv.allday <- ggplot(cv.all.day, aes(x=dpi, y=bird_cv, shape=groups, color=temp))+
  geom_line(aes(x=dpi, y=bird_cv, linetype = treatment, color=temp))+
  geom_point(size=3)+
  labs(x="Day Post Inoculation", y="CV", linetype= "Treatment Groups", shape="Treatment Groups",
       color = "Temperature")+
  scale_color_manual(values=c("blue", "red"))+
  scale_shape_manual(values=c(0, 15, 1, 16))+
  scale_linetype_manual(values=c("dashed","solid"))+
  ylim(c(0,.1))+
  theme_bw()

cv.allday



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

#Add variable indicating the phago_score mean for each treatment group
f.ti <- ti.f%>%
  group_by(groups)%>%
  mutate(groupmean = mean(fever_change))%>%
  ungroup()
f.ti$groupmean
g.ti <- f.ti %>%
  mutate(meanall = mean(fever_change)) %>%
  ungroup()
g.ti$meanall

h.ti <- ti.f %>%
  group_by(groups, dpi) %>%
  mutate(daymean = mean(fever_change))%>%
  ungroup()

#histogram showing variability in fever change across treatment groups
ggplot(f.ti, aes(x=fever_change, fill=groups))+
  geom_histogram(binwidth = 0.5, position="identity", alpha=1, color="black")+
  #geom_density()+
  geom_vline(data = f.ti, aes(xintercept = groupmean, color = "red"), linetype = "dotted",alpha = 1, show.legend=FALSE) + #mean of each group to show differences
  geom_vline(data=g.ti, aes(xintercept = meanall, color = "black"), linetype="dashed", alpha=1, show.legend=FALSE)+ #mean of all groups combined to compare
  scale_fill_manual(values=c( "gray75", "white", "gray35","black"), labels=treat_names)+
  scale_color_manual(values = c("black", "red"), guide = guide_legend(title= NULL))+
  facet_wrap(~groups, nrow=4, labeller = as_labeller(treat_names_n), scales="fixed")+
  labs(x="Fever Change", y="Count", fill="Treatment")+
  facet_wrap(~groups, nrow=4)+
  theme_bw()
  
#animation
anim_hist<- ggplot(f.ti, aes(x=fever_change, fill=groups))+
  geom_histogram(binwidth = 0.5, position="identity", alpha=1, color="black")+
  geom_vline(data = f.ti, aes(xintercept = groupmean, color = "gray"), linetype = "dashed",alpha = 1, show.legend=FALSE) + #mean of each group to show differences
  geom_vline(data = h.ti, aes(xintercept = daymean, color = "red"), linetype = "dotted",alpha = 1, show.legend=FALSE) + #mean of each group to show differences
  #geom_vline(data=g.ti, aes(xintercept = meanall, color = "black"), linetype="solid", alpha=1, show.legend=FALSE)+ #mean of all groups combined to compare
  scale_fill_manual(values=c( "gray75", "white", "gray35","black"), labels=treat_names)+
  scale_color_manual(values = c("gray", "red", "black"), guide = guide_legend(title= NULL))+
  facet_wrap(~groups, nrow=4, labeller = as_labeller(treat_names_n), scales="fixed")+
  labs(title = "Fever: DPI {closest_state}", x="Fever Change", y="Count", fill="Treatment")+
  facet_wrap(~groups, nrow=4)+
  transition_states(dpi, transition_length = 2, state_length = 2) +
  exit_fade()+
  theme_bw()

#animate(anim_hist, renderer = gifski_renderer())
#anim_save("/Users/jesse/Documents/Virginia Tech/Research/Temp + Immunity 2022/Animated Plots/fever_score_hist_grps.gif",
#          animation = last_animation(), renderer = gifski_renderer(), format = "gif")

anim_dens<- ggplot(f.ti, aes(x=fever_change, fill=groups))+
    #geom_histogram(binwidth = 0.5, position="identity", alpha=1, color="black")+
    geom_density()+
    geom_vline(data = f.ti, aes(xintercept = groupmean, color = "gray"), linetype = "dashed",alpha = 1, show.legend=FALSE) + #mean of each group to show differences
    geom_vline(data = h.ti, aes(xintercept = daymean, color = "red"), linetype = "dotted",alpha = 1, show.legend=FALSE) + #mean of each group to show differences
    #geom_vline(data=g.ti, aes(xintercept = meanall, color = "black"), linetype="solid", alpha=1, show.legend=FALSE)+ #mean of all groups combined to compare
    scale_fill_manual(values=c( "gray75", "white", "gray35","black"), labels=treat_names)+
    scale_color_manual(values = c("gray", "red", "black"), guide = guide_legend(title= NULL))+
    facet_wrap(~groups, nrow=4, labeller = as_labeller(treat_names_n), scales="fixed")+
    labs(title = "Fever: DPI {closest_state}", x="Fever Change", y="Count", fill="Treatment")+
    facet_wrap(~groups, nrow=4)+
    transition_states(dpi, transition_length = 2, state_length = 2) +
    exit_fade()+
  theme_bw()
  
  
  #animate(anim_dens, renderer = gifski_renderer())
  #anim_save("/Users/jesse/Documents/Virginia Tech/Research/Temp + Immunity 2022/Animated Plots/fever_score_dens.gif",
  #          animation = last_animation(), renderer = gifski_renderer(), format = "gif")
