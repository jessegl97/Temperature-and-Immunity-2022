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

#Add column for fever score by finding average score
ti <- ti %>%
  group_by(bird_ID) %>%
  mutate(fever_score=((((l_min + l_max)/2) + (r_min + r_max)/2)/2))


ti.f <- ti %>%
  filter(!is.na(fever_score))

ti.f <- ti.f %>%
  filter(dpi !=24)

#fever change = score - baseline
ti.f <- ti.f %>%
  group_by(band_number) %>%
  mutate(fever_change = fever_score - first(fever_score))

#fever change from previous score
ti.f <- ti.f %>%
  group_by(band_number) %>%
  mutate(fever_diff = fever_score - lag(fever_score))

ti.f %>%
  dplyr::select(dpi, treatment, temp, groups)%>%
  tbl_summary(
    by="dpi"
  )%>%
  modify_header(
    label ~ "**All Birds**"
  )
