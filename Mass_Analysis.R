#T+I 2022 MG Bird Mass Analysis
rm(list=ls())
setwd("/Users/jesse/Documents/Virginia Tech/Research/Temp + Immunity 2022/RAW DATA")
library(ggplot2)
library(tidyverse)

ti <- read.csv("ti_merged_data.csv")

#combine dpi 9 and 10
ti$dpi <- ifelse(ti$dpi %in% c(9, 10), "9", as.integer(ti$dpi))
ti$dpi <- as.integer(ti$dpi)
unique(ti$dpi)

ti$threshold_cutoff = 96
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

#Add column for total eye score by adding l and r eye score
ti <- ti %>%
  group_by(bird_ID) %>%
  mutate(total_eye_score=(l_eye_score+r_eye_score))

#make total eye score an integer by multiplying by two for models
ti$tes <- ti$total_eye_score*2

#Add column for fever score by finding average score
ti <- ti %>%
  group_by(bird_ID) %>%
  mutate(fever_score=((((l_min + l_max)/2) + (r_min + r_max)/2)/2))

table(ti$dpi)
ti <- ti %>%
  filter(dpi !=-28)%>%
  drop_na(mass)

#mass change from baseline
ti.m <- ti %>%
  group_by(band_number) %>%
  mutate(mass_change = mass - first(mass))


#what does mass look like
ggplot(ti.m, aes(x=dpi, y=mass_change, color=groups))+
  geom_point()+
  geom_line(aes(group=as.factor(band_number)))+
  geom_hline(yintercept = 0)+
  labs(x="Days Post Infection", y="Mass Change from Baseline", color= "Treatment Groups")+
  #scale_color_manual(values = c("blue", "red"))+
  facet_wrap(~groups)

ggplot(ti.m, aes(x=mass_change, fill=groups))+
  geom_histogram(color="black", binwidth =0.25)+
  geom_vline(xintercept = 0)+
  facet_wrap(~groups, nrow=4)

#ggplot(ti.m, aes(x=mass_change, fill=groups))+
  geom_histogram(color="black", binwidth=0.25)+
  geom_vline(xintercept = 0)+
  facet_wrap(~groups, ncol=1)+
  labs(title= "Day Post Inoculation: {closest_state}", fill = "Groups")+
  transition_states(dpi, transition_length = 2, state_length = 2) +
  enter_fade()+
  exit_fade()


lm1 <- glmmTMB(mass_change ~temp+treatment +(1|band_number),data=ti.m)
summary(lm1)
plot(allEffects(lm1))
unique(ti.m$dpi)
lm2 <- glm(mass_change~temp+treatment, data=ti.m %>% filter(dpi==21))
summary(lm2)

unique(ti.m$dpi)

ggplot(ti.m %>% filter(dpi != -12), aes(x=groups, y=mass_change, shape = groups))+
  geom_jitter(width=0.25)+
  stat_summary(aes(group=groups), fun=mean, geom="point", size=2, color="red")+
  scale_shape_manual(values = c(1, 16, 0, 15))+
  facet_wrap(~dpi, ncol=6)
 