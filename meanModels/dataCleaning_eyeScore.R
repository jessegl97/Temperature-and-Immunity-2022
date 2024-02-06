ti <- read.csv("ti_merged_data.csv")

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

ti.mg$dpi <- as.factor(ti.mg$dpi)

ti.mg.mod <- ti.mg %>% filter(dpi != "-12" & dpi != "3"& dpi != "35")
