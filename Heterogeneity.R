##Heterogeneity
rm(list=ls())

library(tidyverse)
library(dplyr)
setwd("/Users/jesse/Documents/Virginia Tech/Research/Temp + Immunity 2022/RAW DATA")
ti <- read.csv("Temp_Immunity_Master_Data_RAW.csv")
hapto <- read.csv("haptoglobin_interpolated.csv")
vital <- read.csv("temp_immunity_vital_data_raw.csv")

####Haptoglobin####
ti <-left_join(ti, vital)

ti$inf_temp <- paste(ti$treatment, ti$temp, sep = "_")

ti<- mutate(ti,inf_temp=factor (inf_temp, 
                                levels=c("Sham_TN", "MG_TN", "Sham_ST", "MG_ST")))

ti <- mutate(ti,treatment=factor (treatment, 
                                  levels=c("Sham", "MG")))

ti$temp_names <- ifelse(ti$temp == "ST", "Cold",
                        ifelse(ti$temp == "TN", "Warm", NA))

ti$Groups <- ifelse(ti$inf_temp == "Sham_TN", "Warm Control",
                    ifelse(ti$inf_temp == "MG_TN", "Warm Infected",
                           ifelse(ti$inf_temp == "Sham_ST", "Cold Control",
                                  ifelse(ti$inf_temp == "MG_ST", "Cold Infected", NA))))

ti <- mutate(ti,Groups=factor (Groups, 
                               levels=c("Warm Control", "Warm Infected", "Cold Control", "Cold Infected")))
#left join ti, hapto, and vital data
ti.h<-left_join(ti, hapto, by='bird_ID')
ti.h<-left_join(ti.h, vital)
ti.h$hapto <- as.numeric(ti.h$hapto_ext)

#exclude plate 1 that was run incorrectly
ti.h<-ti.h %>%
  filter(hapto_plate > 1)%>%
  drop_na(hapto_ext, treatment)

#create dataframe
ti.hp <- ti.h[c('hapto_ext','red', 'treatment', 'dpi', 'temp', 'Groups', 'inf_temp', 'band_number', 'bird_ID', 'Sex')]

##Calculate hapto difference while maintaining bird_ID column
ti.hpb <- ti.hp %>%
  group_by(band_number)%>%
  filter(all(c(-12,2) %in% dpi)) %>%
  pull(band_number) %>%
  unique()
ti.hpb

ti.hp$band_number <- as.numeric(ti.hp$band_number)
ti.hpb$band_number <- as.numeric(ti.hpb)
all.bn <- ti.hp$band_number

#make list of all band_number that don't have hapto values on both dpi -12 and 2
#anti_join list of all band_numbers and list of band_numbers with values above filtered out
a= anti_join(ti.hp$band_number, ti.hpb)

filtered_ti.hpb <- ti.hp %>%
  filter(band_number %in% ti.hpb)
filtered_ti.hpb

filtered_ti.hpb$band_number <- as.factor(filtered_ti.hpb$band_number)
detach(package:plyr)
difference <- filtered_ti.hpb %>%
  dplyr::group_by(band_number)%>%
  dplyr::summarize(diff = diff(hapto_ext[dpi == 2]- hapto_ext[dpi == -12]))#subtract dpi 2 hapto_ext from dpi -12 hapto_ext
difference

#I can't get this to work - This code returns the differences but I need it to return an NA every other line
difference <- filtered_ti.hpb %>%
  group_by(bird_ID, band_number) %>%
  mutate(diff = hapto_ext[dpi == 2] - hapto_ext[dpi == -12])

ti.hr <- filtered_ti.hpb %>%
  drop_na(band_number, hapto_ext, dpi, treatment)

ti.hr$band_number_c = as.character(ti.hr$band_number)

unique(ti.hr$dpi)
str(ti.hr$dpi)
unique(ti.hr$band_number_c)
ti.hr$timepoint[ti.hr$dpi==(-12)]=1
ti.hr$timepoint[ti.hr$dpi==(2)]=2


ti.hr = ti.hr %>%
  dplyr::group_by(band_number_c)%>%
  dplyr::arrange(band_number, timepoint)%>%
  mutate(hapt.time.point.neg12 = dplyr::lag(hapto_ext))%>%
  mutate(hapt.diff = dplyr::lag(hapto_ext)-hapto_ext)
ti.hr$hapt.diff


ti.hp <- ti.hp %>%
  left_join(difference, by = c("bird_ID", "band_number"))


ti.hp <- ti.hr[c('hapto_ext', 'hapt.diff', 'red', 'treatment', 'dpi', 'temp', 'Groups', 'inf_temp', 'band_number', 'bird_ID', 'Sex')]


#subset into treatment groups and calculate means for each group
group.means <- ti.hp %>%
  dplyr::group_by(Groups)%>%
  dplyr::summarise(hapto.means = mean(hapt.diff))


#merge into larger dataset
ti.hp <- ti.hp%>%
  left_join(group.means, by= "Groups")
ti.hp

#Calculate residuals
ti.hp <- ti.hp %>%
  group_by(Groups)%>%
  mutate(residuals2= (hapto_ext - hapto.means)^2)

ti.hp$residuals2

g.resid <-ggplot(data=ti.hp, aes(x=Groups, 
                                    y=residuals2, color=Groups))+
  geom_histogram()
g.resid



####Fever####
fever <- read.csv("fever_assay_raw.csv")

ti.fever <- left_join(ti, fever)

ti.fever <- ti.fever %>%
  group_by(bird_ID) %>%
  mutate(fever_score=(((l_min+r_min)/2)))

ti.fever$fever_score

ti.fever<- ti.fever%>%
  drop_na(fever_score, band_number, dppi, treatment)
