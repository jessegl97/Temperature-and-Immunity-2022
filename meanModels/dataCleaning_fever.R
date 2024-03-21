ti <- read.csv("ti_merged_data.csv")

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

#new column with whether the bird is diseased on a particular day
ti <- ti %>%
  mutate(diseased= ifelse(total_eye_score>sympt_cutoff, 1, 0), na.rm=T)


#new column with whether the bird is ever diseased
ti <- ti %>%
  group_by(band_number) %>%
  mutate(ever_diseased = ifelse(any(coalesce(diseased, 0) == 1), 1, 0), na.rm=T) %>%
  ungroup()

#Add column for fever score by finding average score
ti <- ti %>%
  group_by(bird_ID) %>%
  mutate(fever_score=((((l_min + l_max)/2) + (r_min + r_max)/2)/2))

ti.tr <- ti %>%
  drop_na (band_number, fever_score, dpi, treatment, temp, groups, sex)

ti.tr <- ti.tr %>%
  filter(dpi !=24)

#fever change from baseline
ti.f <- ti.tr %>%
  group_by(band_number) %>%
  mutate(fever_change = fever_score - first(fever_score))

#fever change from previous score
ti.f <- ti.f %>%
  group_by(band_number) %>%
  mutate(fever_diff = fever_score - lag(fever_score))
