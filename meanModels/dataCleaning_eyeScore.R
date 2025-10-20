setwd('/Users/jesse/Documents/GitHub/Temperature-and-Immunity-2022/meanModels')
ti <- read.csv("ti_merged_data.csv")

#Add column for total eye score by adding l and r eye score
ti <- ti %>%
  group_by(bird_ID) %>%
  mutate(total_eye_score=(l_eye_score+r_eye_score))

#make total eye score an integer by multiplying by two for models
ti$tes <- ti$total_eye_score*2

ti$quant_cutoff = 100
ti$sympt_cutoff = 0.1

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

#df with only eye score days
ti.e <- ti %>%
  filter(!is.na(total_eye_score))

#data frame with only inoculated birds
ti.mg <- ti.e%>%
  filter(treatment == "Inoculated")

hist(ti.mg$total_eye_score)
hist(ti$total_eye_score)
hist(ti.mg$tes)
ti.mg$dpi <- as.factor(ti.mg$dpi)

unique(ti.mg$dpi)

