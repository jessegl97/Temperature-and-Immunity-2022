ti <- read.csv("data_frames/TI22_merged_data.csv")

library(tidyverse)
library(gtsummary)
#Add column for total eye score by adding l and r eye score
ti <- ti %>%
  group_by(bird_ID) %>%
  mutate(total_eye_score=(l_eye_score+r_eye_score))

ti <- ti %>%
  dplyr::mutate(groups = factor(groups,
                                    levels = c("Warm Control", "Warm Inoculated",
                                               "Cold Control", "Cold Inoculated")))

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
  mutate(fever_score=((((l_min + l_max)/2) + (r_min + r_max)/2)/2))%>%
  ungroup()

ti %>%
  dplyr::select(dpi, treatment, temp, groups)%>%
  tbl_summary(
    by="dpi"
  )%>%
  modify_header(
    label ~ "**All Birds**"
  )



ti.f <- ti %>%
  filter(dpi %in% c(0, 3, 7, 14, 18, 24, 28, 35) & !is.na(fever_score))

ti.fw <- ti %>%
  filter(dpi %in% c(0, 3, 7, 14, 18, 24, 28, 35))

missing.fever<-anti_join(ti.fw, ti.f, by="bird_ID") %>%
  dplyr::select(dpi, band_number, bird_ID, temp, treatment, fever_score, experiment_notes)

# ti.f <- ti.f %>%
#   filter(dpi !=24)

ti.f %>%
  dplyr::select(dpi, treatment, temp, groups)%>%
  tbl_summary(
    by="dpi"
  )%>%
  modify_header(
    label ~ "**Fever; NA Removed**"
  )


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
    label ~ "**Fever Sample Sizes**"
  )
