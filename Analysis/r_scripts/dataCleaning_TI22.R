#Create Merged Dataset for T+I 2022

setwd('/Users/jesse/Documents/GitHub/Temperature-and-Immunity-2022/Analysis/')
library(ggplot2)
library(tidyverse)
library(tibble)
library(gridExtra)
library(gtsummary)
ti <- read.csv("data_frames/TI22_merged_data.csv")

ti %>%
  dplyr::select(dpi, treatment, temp, groups)%>%
  tbl_summary(
    by=dpi
  )%>%
  modify_header(
    label ~ "**All Birds**"
  )

ti_full <- ti

colnames(ti_full)

ti <- ti_full %>% 
  dplyr::select(1:4, 6:52, 61:70)

#Add column for total eye score by adding l and r eye score
ti <- ti %>%
  group_by(bird_ID) %>%
  mutate(total_eye_score=(l_eye_score+r_eye_score)) %>%
  ungroup()

colnames(ti)

#ti$dpi <- as.factor(ti$dpi)
ti$band_number <- as.factor(ti$band_number)
ti$sex <- as.factor(ti$sex)
ti$temp <- as.factor(ti$temp)
ti$treatment <- as.factor(ti$treatment)
ti$groups <- as.factor(ti$groups)

summary(ti)

ggplot(ti, aes(x=dpi, y=total_eye_score))+
  geom_jitter(width=0.2)

ggplot(ti, aes(x=dpi, y=l_max+r_max, color=temp))+
  geom_jitter(width=0.2)


ti %>%
  dplyr::select(dpi, treatment, temp, groups)%>%
  tbl_summary(
    by=dpi
  )%>%
  modify_header(
    label ~ "**All Birds**"
  )


#combine dpi 9 and 10
ti$dpi <- ifelse(ti$dpi %in% c(9, 10), "9", as.integer(ti$dpi))
ti$dpi <- as.integer(ti$dpi)
unique(ti$dpi)
ti$quantity
#add infected and seropositivity thresholds
ti.cont <- ti%>%
  filter(treatment == "Control")%>%
  drop_na(quantity)

ti.cont
range(ti.cont$quantity) #highest control quantity = 95.38
ti$quant_cutoff = 50
ti.cont$elisa_od
ti$seropos_cutoff = 0.061
ti$sympt_cutoff = 0.1

#Add column for total eye score by adding l and r eye score
ti <- ti %>%
  group_by(bird_ID) %>%
  mutate(total_eye_score=(l_eye_score+r_eye_score))

#diseased = eye score
#infected = eye score or pathology >= cutoff
ti <- ti %>%
  mutate(diseased= ifelse(total_eye_score>=sympt_cutoff, 1, 0),
         sick = ifelse(total_eye_score>=sympt_cutoff | quantity >= quant_cutoff, 1, 0),
         infected = ifelse(quantity >= quant_cutoff, 1, 0))

#new column with whether the bird is ever diseased
ti <- ti %>%
  group_by(band_number) %>%
  mutate(ever_diseased = ifelse(any(coalesce(diseased, 0) == 1), 1, 0),
         ever_infected = ifelse(any(coalesce(infected, 0) == 1), 1, 0),
         ever_sick = ifelse(any(coalesce(sick, 0) == 1), 1,0))%>%
  ungroup()

ti <- ti %>%
  mutate(sex = dplyr::recode(sex,
                      "F" = "Female",
                      "M" = "Male"))

ti <- ti %>%
  mutate(temp = fct_relevel(temp,
                            "Warm", "Cold"))

#identify seropositive birds from quarantine in dataset
Qinf <- ti %>%
  dplyr::filter(band_number %in% c(2667, 2750)) %>%
    #dpi=="-28" & elisa_od >=0.061)%>%
  dplyr::select(band_number, elisa_od, groups, dpi, elisa_od)
Qinf


#remove birds that were seropositive dpi -28 from antibody analysis
ti.ab <- ti %>%
  filter(band_number != 2667 & band_number != 2750)
