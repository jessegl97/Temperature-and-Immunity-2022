ti <- read.csv("ti_merged_data.csv")

#select only phagocytosis assay and make new wbc total, phago total, and phago score columns
ti <- ti %>%
  filter(Well.1.WBC != "NA") %>% #select only phagocytosis assay
  mutate(wbc_total=Well.1.WBC+Well.2.WBC+Well.3.WBC+Well.4.WBC)%>% #sum of wbc quadruplicates
  mutate(phago_total=Well.1.Phago+Well.2.Phago+Well.3.Phago+Well.4.Phago) %>% #sum of phago quadruplicates
  mutate(phago_score = phago_total/(wbc_total+phago_total)) #phago score = phagocytic / total

ti <- ti%>%
  filter(dpi != 2) %>% #remove dpi 2 because only a small subset were run day 2
  drop_na(phago_score)