ti <- read.csv("ti_merged_data.csv")

ti %>%
  dplyr::select(dpi, treatment, temp, groups)%>%
  tbl_summary(
    by="dpi"
  )%>%
  modify_header(
    label ~ "**All Birds**"
  )

#combine dpi 9 and 10
ti$dpi <- ifelse(ti$dpi %in% c(9, 10), "9", as.integer(ti$dpi))
ti$dpi <- as.integer(ti$dpi)
unique(ti$dpi)
ti$quant
#add infected and seropositivity thresholds
ti.cont <- ti%>%
  filter(treatment == "Sham")%>%
  drop_na(quantity)

ti.cont
range(ti.cont$quantity) #highest control quantity = 95.38
ti$quant_cutoff = 100
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

#identify seropositive birds from quarantine in dataset
Qinf <- ti %>%
  filter(dpi==-28)%>%
  filter(elisa_od >=0.061)%>%
  dplyr::select(band_number, elisa_od, groups, dpi)
Qinf



#remove birds that were seropositive dpi -28 from antibody analysis
ti <- ti %>%
  filter(!band_number %in% c(2667, 2750))

#remove dpi that didn't have ELISA data run
ti <- ti %>%
  filter(dpi %in% c(-28, 9, 28))

ab.miss <- ti %>%
  filter(is.na(elisa_od)) %>%
  dplyr::select(dpi, band_number, treatment, temp, elisa_od, total_eye_score, groups)

ab.miss %>%
  dplyr::select(dpi, treatment, temp, groups)%>%
  tbl_summary(
    by="dpi"
  )%>%
  modify_header(
    label ~ "**All Birds**"
  )


ggplot(ab.miss, aes(x=dpi, y=total_eye_score))+
  geom_jitter(width=0.5, height=0)

#remove any bird missing at least one elisa_od
ti.rm <- ti %>%
  filter(!band_number %in% c(unique(ab.miss$band_number)))

unique(ti.rm$band_number)
#Across all days post-infection (dpi 9 and 28)
ti.pi <- ti #%>%
  #filter(dpi > 0)
unique(ti.pi$dpi)
