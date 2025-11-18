ti <- read.csv("data_frames/TI22_merged_data.csv")

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

#remove dpi that didn't have ELISA data run
ti.ab <- ti %>%
  filter(dpi %in% c(-28, 9, 28))


##Identify birds to remove from antibody analysis
ti.ab %>%
  dplyr::select(dpi, groups)%>%
  tbl_summary(
    by=dpi
  )%>%
  modify_header(
    label ~ "**before removal**"
  )

#identify seropositive birds from quarantine in dataset
Qinf <- ti %>%
  filter(dpi==-28)%>%
  filter(elisa_od >=0.061)%>%
  dplyr::select(band_number, elisa_od, groups, dpi)
Qinf

#DPI -28: 2667, 2750 to remove seropositive at baseline
ab.miss <- ti.ab %>%
  filter(is.na(elisa_od)) %>%
  dplyr::select(dpi, band_number, treatment, temp, elisa_od, total_eye_score, exclude, reason)

#Only ran a subset of controls on dpi 9/10
#Not 2684, 2893, 2847

#DPI 9/10: 2684, 2893, 2847 not run (Controls)
#DPI 28: 2684, 2893 not run (Controls)

#Remove birds with high CV or other issues
ti.ab.rm <- ti.ab %>%
  filter(exclude == 1) %>%
  select(band_number, bird_ID, groups, elisa_od, elisa_cv, reason)

#DPI 9/10: 2573 removed high CV
#DPI 28: 2921, 2847 removed high CV
#DPI 28: 2869 removed not enough plasma first run

##Total to remove:

#DPI -28:  (n=2) 2667, 2750 to remove - seropositive at baseline
#DPI 9/10: (n=4) 2684, 2893, 2847, 2573 to remove - exclude controls and high CV
#DPI 28:   (n=3) 2684, 2893, 2869 to remove - exclude controls and high CV

#Accounting for repeated measures, birds removed:
  #2667, 2750, 2684, 2893, 2847, 2573, 2869
    # n = 7 
remove.ab <- ti.ab %>% filter(band_number %in% c(2667, 2750, 2684, 2893, 2847, 2573, 2869))

remove.ab %>%
  dplyr::select(dpi, groups)%>%
  tbl_summary(
    by=dpi
  )%>%
  modify_header(
    label ~ "**birds to remove**"
  )

#Remove birds
ti.ab <- ti.ab %>%
  filter(!band_number %in% c(2667, 2750, 2684, 2893, 2847, 2573, 2869))

ti.ab %>%
  dplyr::select(dpi, groups)%>%
  tbl_summary(
    by=dpi
  )%>%
  modify_header(
    label ~ "**after removal**"
  )


unique(ti.ab$band_number)

