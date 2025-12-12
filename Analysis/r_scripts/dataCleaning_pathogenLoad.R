setwd('/Users/jesse/Documents/GitHub/Temperature-and-Immunity-2022/Analysis/data_frames/')
ti <- read.csv("TI22_merged_data.csv")

#Add column for total eye score by adding l and r eye score
ti <- ti %>%
  group_by(bird_ID) %>%
  mutate(total_eye_score=(l_eye_score+r_eye_score))

ti$quant_cutoff = 50
ti$sympt_cutoff = 0.1

#Combine DPI 9 and 10 as DPI 9 to analyze together
  #Inocualted birds sampled DPI 9, control birds sampled DPI 10 so timing should not matter
ti$dpi <- ifelse(ti$dpi %in% c(9, 10), "9", as.integer(ti$dpi))
ti$dpi <- as.integer(ti$dpi)

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


ti %>%
  dplyr::select(dpi, treatment, temp, groups)%>%
  tbl_summary(
    by="dpi"
  )%>%
  modify_header(
    label ~ "**Quantity Before Removal**"
  )

#who is missing dpi 21
who <- ti %>%
  filter(dpi %in% c(9,21, 28) & groups == "Warm Inoculated")

#2581 DPI 21 not extracted (Warm Inoculated)
who %>% filter(band_number == 2581) %>% dplyr::select(band_number, dpi, total_eye_score, quantity, experiment_notes)

#df with only pathogen load days
  #Removes 2581 DPI 21 and keeps only DPIs -12, 3, 9, 10, 14, 21, and 28
ti.pl <- ti %>%
  filter(!is.na(quantity))

ti.pl %>%
  dplyr::select(dpi, treatment, temp, groups)%>%
  tbl_summary(
    by="dpi"
  )%>%
  modify_header(
    label ~ "**Quantity After Removal**"
  )
