ti <- read.csv("ti_merged_data.csv")

ti$quant_cutoff = 100
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

#select only phagocytosis assay and make new wbc total, phago total, and phago score columns
ti.p <- ti %>%
  filter(Well.1.WBC != "NA") %>% #select only phagocytosis assay
  mutate(wbc_total=Well.1.WBC+Well.2.WBC+Well.3.WBC+Well.4.WBC)%>% #sum of wbc quadruplicates
  mutate(phago_total=Well.1.Phago+Well.2.Phago+Well.3.Phago+Well.4.Phago) %>% #sum of phago quadruplicates
  mutate(phago_score = phago_total/(wbc_total+phago_total)) #phago score = phagocytic / total

ti.p <- ti.p%>%
  filter(dpi != 2) %>% #remove dpi 2 because only a small subset were run day 2
  drop_na(phago_score)
