ti <- read.csv("ti_merged_data.csv")

#combine dpi 9 and 10
ti$dpi <- ifelse(ti$dpi %in% c(9, 10), "9", as.integer(ti$dpi))
ti$dpi <- as.integer(ti$dpi)
unique(ti$dpi)
ti$quant
#add infected and seropositivity thresholds
ti.cont <- ti%>%
  filter(treatment == "Control")%>%
  drop_na(quantity)

ti.cont
range(ti.cont$quantity) #highest control quantity = 95.38
ti$threshold_cutoff = 96
ti.cont$elisa_od
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

#identify seropositive birds from quarantine in dataset
Qinf <- ti %>%
  filter(dpi==-28)%>%
  filter(elisa_od >=0.061)%>%
  select(band_number, elisa_od, groups, dpi)
Qinf


#remove birds that were seropositive dpi -28 from antibody analysis
ti <- ti %>%
  filter(band_number != 2667 & band_number != 2750)

#remove dpi that didn't have ELISA data run
ti <- ti %>%
  filter(dpi == -28 | dpi == 9 | dpi == 28)

####Analysis####
### Overall Hypotheses: Pathogen-induced antibody response will differ between hot and cold rooms
#Ho: Antibody levels were the same between hot and cold infected groups
#Ha: Antibody levels differed between hot and cold infected groups

#Across all days post-infection (dpi 9 and 28)
ti.pi <- ti %>%
  filter(dpi > 0)
unique(ti.pi$dpi)