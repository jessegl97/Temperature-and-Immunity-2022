####When do birds first develop eye scores and when are they first infected?####

ggplot(ti.inoc, aes(x=dpi, y=total_eye_score, color=temp))+
  geom_point(alpha=0.75)+
  geom_line(aes(group = band_number), alpha=0.5)+
  facet_wrap(~ever_diseased)




ggplot(ti.inoc, aes(x=dpi, y=total_eye_score, color=temp))+
  geom_point(alpha=0.75)+
  geom_line(aes(group = band_number), alpha=0.5)+
  facet_wrap(~ever_infected)

ggplot(ti.inoc, aes(x=dpi, y=total_eye_score, color=temp))+
  geom_point(alpha=0.75)+
  geom_line(aes(group = band_number), alpha=0.5)+
  facet_wrap(~ever_sick)

ti.e <- ti.inoc %>%
  dplyr::select(band_number, dpi, total_eye_score, quantity, log10_quantity, treatment, temp,
                diseased, ever_diseased, infected, ever_infected)

#first time infected or diseased
#grouping by ever_diseased and ever_infected means that I can pick which I want to use at later analyses. 
firsts <- ti.e %>%
  #filter(ever_diseased == 1) %>%
  group_by(band_number, ever_diseased, ever_infected) %>%
  mutate(
    # First pathology
    first_path = if (any(diseased == 1, na.rm = TRUE)) {
      min(dpi[diseased == 1], na.rm = TRUE)
    } else {
      NA
    },
    # First infection
    first_inf = if (any(infected == 1, na.rm = TRUE)) {
      min(dpi[infected == 1], na.rm = TRUE)
    } else {
      NA
    }
  ) %>%
  ungroup()

# First infection (days to first infection for inoculated birds that became infected by pathogen load)
first_inf <- ti.e %>%
  filter(ever_infected == 1) %>%
  group_by(band_number) %>%
  mutate(
    first_inf = if (any(infected ==1, na.rm = TRUE)) {
      min(dpi[infected ==1], na.rm = TRUE)
    } else {
      NA
    }
  ) %>%
  ungroup()

#first infection by pathogen load for inoculated birds that became infected by pathogen load
avg_time_inf <- first_inf %>% 
  filter(dpi == 3) %>%
  group_by(temp)%>%
  summarise(avg_to_inf = mean(first_inf, na.rm=T),
            med_to_inf = median(first_inf, na.rm=T),
            n_inf = sum(!is.na(first_inf)))

#Average time to pathology for birds that are infected and median time. These both may be
#misleading because eye score was only taken DPI 9 and 14. It is likely birds developed scores earlier
#But including time to infection helps a little bit because most of the birds showed path loads by day 9
avg_time <- firsts %>% 
  filter(dpi == 3) %>%
  group_by(temp)%>%
  summarise(avg_to_path = mean(first_path, na.rm=T),
            med_to_path = median(first_path, na.rm=T),
            n_path = sum(!is.na(first_path)),
            avg_to_inf = mean(first_inf, na.rm=T),
            med_to_inf = median(first_inf, na.rm=T),
            n_inf = sum(!is.na(first_inf)))

# Reshape the data for easier plotting
avg_time_long <- avg_time %>%
  pivot_longer(
    cols = c(avg_to_path, med_to_path, avg_to_inf, med_to_inf),
    names_to = "metric",
    values_to = "time"
  ) %>%
  mutate(
    definition = case_when(
      grepl("path", metric) ~ "Pathology",
      grepl("inf", metric) ~ "Infection"
    ),
    statistic = case_when(
      grepl("avg", metric) ~ "Average",
      grepl("med", metric) ~ "Median"
    )
  )

# Plot
ggplot(avg_time_long, aes(x = definition, y = time, fill = definition)) +
  geom_bar(stat = "identity", position = "dodge", color="black") +
  facet_wrap(~ statistic) +
  labs(
    x = "Temperature",
    y = "Time (dpi)",
    fill = "Event Type",
    title = "Average and Median Times to Pathology and Infection"
  ) +
  scale_fill_manual(values = c("Pathology" = "black", "Infection" = "white")) +
  theme_minimal()

# Reshape the count data
count_data <- avg_time %>%
  dplyr::select(temp, n_path, n_inf) %>%
  pivot_longer(
    cols = c(n_path, n_inf),
    names_to = "event",
    values_to = "count"
  ) %>%
  mutate(event = case_when(
    event == "n_path" ~ "Diseased (Eye Score)",
    event == "n_inf" ~ "Infected (Path Load)"
  ))

# Plot
ggplot(count_data, aes(x = event, y = count, fill = temp)) +
  geom_bar(stat = "identity", position = "dodge", width=0.5) +
  labs(
    x = "Definition",
    y = "Number of Birds",
    fill = "Temperature",
    title = "Number of Birds Diseased and Infected by Temperature"
  ) +
  scale_fill_manual(values = c("Warm" = "#C1121F", "Cold" = "#669BBC")) +
  theme_minimal(base_size=16)

#Still only inoculated birds with 
ti.ep <- merge(ti.e, firsts, by=c("band_number", "dpi"))

colnames(ti.ep)

colnames(ti.ep) <- gsub("\\.x$", "", colnames(ti.ep))
ti.ep <- ti.ep[, !grepl("\\.y$", colnames(ti.ep))]
colnames(ti.ep)


#how many birds were infected / diseased each sampling day? 
#includes birds that did not get sick but were inoculated
summary_counts <- ti.ep %>%
  group_by(temp, dpi) %>%
  dplyr::summarize(
    n_diseased = sum(diseased, na.rm = TRUE),
    n_infected = sum(infected, na.rm = TRUE),
    total_birds = n()
  ) %>%
  ungroup()


#Visualize the number of individuals who were diseased or infected across infection
ggplot(summary_counts, aes(x = dpi, color = temp)) +
  geom_line(aes(y = n_diseased, linetype = "Diseased"), size = 1) +
  geom_line(data=summary_counts %>% filter(dpi %in% c(3, 9, 14, 28)),
            aes(y = n_infected, linetype = "Infected"), size = 0.5) +
  geom_point(aes(y = n_diseased), size = 2, shape = 16) +
  geom_point(data=summary_counts %>% filter(dpi %in% c(3, 9, 14, 28)),
             aes(y = n_infected), size = 2, shape = 17) +
  labs(
    title = "Number of Birds Infected / Diseased",
    x = "Days Post-Infection (dpi)",
    y = "Number of Birds",
    color = "Temperature",
    linetype = "Metric"
  ) +
  scale_color_manual(values = c("Cold" = "#669BBC", "Warm" = "#C1121F"))


#survival curve
library(survival)
library(survminer)

# Prepare the data
#still includes birds that didn't get sick that were inoculated
#Kaplan-Meir survival curves incorporates censoring (birds that didn't get sick)
surv_data <- firsts %>%
  filter(dpi == 3)%>%
  mutate(
    time_to_disease = first_path,           # Time to first pathology
    event_disease = ifelse(!is.na(first_path), 1, 0),  # Event indicator for disease
    time_to_infection = first_inf,          # Time to first infection
    event_infection = ifelse(!is.na(first_inf), 1, 0)  # Event indicator for infection
  )

# Kaplan-Meier survival fit for disease
surv_fit_disease <- survfit(Surv(time_to_disease, event_disease) ~ temp, data = surv_data)

# Kaplan-Meier survival fit for infection
surv_fit_infection <- survfit(Surv(time_to_infection, event_infection) ~ temp, data = surv_data)

# Plot the survival curves for disease
ggsurvplot(
  surv_fit_disease,
  data = surv_data,
  xlab = "Days Post-Infection (dpi)",
  ylab = "Proportion Not Diseased",
  title = "Survival Curve: Time to First Disease",
  color = "temp",
  palette = c("Cold" = "#669BBC", "Warm" = "#C1121F"),
  conf.int = TRUE,
  pval= TRUE,
  risk.table = TRUE   # Add risk table below the plot
)

# Plot the survival curves for infection
ggsurvplot(
  surv_fit_infection,
  data = surv_data,
  xlab = "Days Post-Infection (dpi)",
  ylab = "Proportion Not Infected",
  title = "Survival Curve: Time to First Infection",
  color = "temp",
  palette = c("Cold" = "#669BBC", "Warm" = "#C1121F"),
  conf.int = TRUE,
  pval= TRUE,
  risk.table = TRUE   # Add risk table below the plot
)

set.seed(1)
#All Birds that were inoculated
ggplot(data=ti.e %>%filter(dpi > 0), aes(x=dpi, y=total_eye_score, color = temp))+
  geom_jitter(size=2, width = 0.2, height=0, alpha=0.5)+
  #geom_smooth(aes(x=dpi, y=total_eye_score, color = temp, group=(band_number)), method="loess", se=FALSE, span=0.5, size=0.1)+
  labs(title ="Total Eye Score of all birds that were inoculated", x="Days Post Infection", y="Eye Score", color="Temperature")+
  stat_summary(aes(group=temp, color=temp), fun=mean, geom="line", alpha=1, size=1)+
  geom_vline(data=avg_time, aes(xintercept = avg_to_path, color=temp), linetype="dashed")+
  scale_color_manual(values=c("#669BBC", "#C1121F"))+
  theme_minimal()

#All birds dpi 14
ggplot(data=ti.e %>% filter(dpi < 14), aes(x=dpi, y=total_eye_score, color = temp))+
  geom_jitter(size=2, width = 0.1, height=0, alpha=0.5)+
  #geom_smooth(aes(x=dpi, y=total_eye_score, color = temp, group=(band_number)), method="loess", se=FALSE, span=0.5, size=0.1)+
  geom_line(aes(x=dpi, y=total_eye_score, color=temp, group=band_number), alpha=0.2)+
  labs(x="Days Post Infection", y="Eye Score", color="Temperature")+
  stat_summary(aes(group=temp, color=temp), fun=mean, geom="line", alpha=1, size=1)+
  geom_vline(data=avg_time, aes(xintercept = avg_to_path, color=temp), linetype="dashed")+
  #geom_vline(data=avg_time, aes(xintercept = avg_to_inf, color=temp), linetype="solid", alpha=0.5)+
  scale_color_manual(values=c("#669BBC", "#C1121F"))+
  theme_minimal()

#Only diseased birds time to infection
ggplot(data=ti.ep %>% filter(ever_infected==1), aes(x=dpi, y=total_eye_score, color = temp))+
  geom_jitter(size=2, width = 0.2, height=0, alpha=0.5)+
  #geom_line(aes(x=dpi, y=total_eye_score, color = temp, group=(band_number)), method="loess", se=FALSE, span=0.5, size=0.1)+
  labs(title="Diseased Birds Only", x="Days Post Infection", y="Eye Score", color="Temperature")+
  stat_summary(aes(group=temp, color=temp), fun=mean, geom="line", alpha=1, size=1)+
  #geom_vline(data=avg_time, aes(xintercept = avg_to_path, color=temp), linetype="dashed", alpha=0.5)+
  #geom_vline(data=avg_time, aes(xintercept = avg_to_inf, color=temp), linetype="solid", alpha=1)+
  scale_color_manual(values=c("#669BBC", "#C1121F"))+
  theme_minimal()

#variability in time to infection?
#install.packages("remotes")
#remotes::install_github("T-Engel/CValternatives")
library(CValternatives)
# Define functions for calculating CV, PV, and V2
calculate_cv <- function(data) {
  return(sd(data) / mean(data))
}

calculate_pv <- function(data) {
  return(PV(data))
}

calculate_v2 <- function(data) {
  mean_x <- mean(data)
  sd_x <- sd(data)
  v2 <- sd_x^2 / (sd_x^2 + mean_x^2)
  return(v2)
}

n_boot = 1000

firsts.dis.var <- firsts %>%
  group_by(temp) %>%
  filter(ever_diseased == 1) %>%
  dplyr::reframe(
    temp = temp,
    band_number = band_number,
    total_eye_score = total_eye_score,
    treatment = treatment,
    
    # Calculations for first_path
    mean_t2p = mean(first_path, na.rm = TRUE),
    median_t2p = median(first_path, na.rm = TRUE),
    bird_cv_t2p = calculate_cv(first_path),
    bird_pv_t2p = calculate_pv(first_path),
    bird_v2_t2p = calculate_v2(first_path),
    
    # Calculations for first_inf
    mean_t2i = mean(first_inf, na.rm = TRUE),
    median_t2i = median(first_inf, na.rm = TRUE),
    bird_cv_t2i = calculate_cv(first_inf),
    bird_pv_t2i = calculate_pv(first_inf),
    bird_v2_t2i = calculate_v2(first_inf),
    
    # Bootstrapping for first_path
    cv_bootstrap_t2p = list(replicate(n_boot, calculate_cv(sample(first_path, replace = TRUE)))),
    pv_bootstrap_t2p = list(replicate(n_boot, calculate_pv(sample(first_path, replace = TRUE)))),
    v2_bootstrap_t2p = list(replicate(n_boot, calculate_v2(sample(first_path, replace = TRUE)))),
    
    # Bootstrapping for first_inf
    cv_bootstrap_t2i = list(replicate(n_boot, calculate_cv(sample(first_inf, replace = TRUE)))),
    pv_bootstrap_t2i = list(replicate(n_boot, calculate_pv(sample(first_inf, replace = TRUE)))),
    v2_bootstrap_t2i = list(replicate(n_boot, calculate_v2(sample(first_inf, replace = TRUE))))
  ) %>%
  mutate(
    # Confidence intervals for first_path
    cv_lower_ci_t2p = map_dbl(cv_bootstrap_t2p, ~ quantile(.x, 0.025, na.rm = TRUE)),
    cv_upper_ci_t2p = map_dbl(cv_bootstrap_t2p, ~ quantile(.x, 0.975, na.rm = TRUE)),
    pv_lower_ci_t2p = map_dbl(pv_bootstrap_t2p, ~ quantile(.x, 0.025, na.rm = TRUE)),
    pv_upper_ci_t2p = map_dbl(pv_bootstrap_t2p, ~ quantile(.x, 0.975, na.rm = TRUE)),
    v2_lower_ci_t2p = map_dbl(v2_bootstrap_t2p, ~ quantile(.x, 0.025, na.rm = TRUE)),
    v2_upper_ci_t2p = map_dbl(v2_bootstrap_t2p, ~ quantile(.x, 0.975, na.rm = TRUE)),
    
    # Confidence intervals for first_inf
    cv_lower_ci_t2i = map_dbl(cv_bootstrap_t2i, ~ quantile(.x, 0.025, na.rm = TRUE)),
    cv_upper_ci_t2i = map_dbl(cv_bootstrap_t2i, ~ quantile(.x, 0.975, na.rm = TRUE)),
    pv_lower_ci_t2i = map_dbl(pv_bootstrap_t2i, ~ quantile(.x, 0.025, na.rm = TRUE)),
    pv_upper_ci_t2i = map_dbl(pv_bootstrap_t2i, ~ quantile(.x, 0.975, na.rm = TRUE)),
    v2_lower_ci_t2i = map_dbl(v2_bootstrap_t2i, ~ quantile(.x, 0.025, na.rm = TRUE)),
    v2_upper_ci_t2i = map_dbl(v2_bootstrap_t2i, ~ quantile(.x, 0.975, na.rm = TRUE))
  ) %>%
  dplyr::select(-cv_bootstrap_t2p, -pv_bootstrap_t2p, -v2_bootstrap_t2p,
                -cv_bootstrap_t2i, -pv_bootstrap_t2i, -v2_bootstrap_t2i) %>%
  ungroup()

summary_tibble <- firsts.dis.var %>%
  group_by(temp) %>%
  reframe(
    # Summaries for first_path
    mean_bird_cv_t2p = mean(bird_cv_t2p, na.rm = TRUE),
    mean_bird_pv_t2p = mean(bird_pv_t2p, na.rm = TRUE),
    mean_bird_v2_t2p = mean(bird_v2_t2p, na.rm = TRUE),
    mean_cv_lower_ci_t2p = mean(cv_lower_ci_t2p, na.rm = TRUE),
    mean_cv_upper_ci_t2p = mean(cv_upper_ci_t2p, na.rm = TRUE),
    mean_pv_lower_ci_t2p = mean(pv_lower_ci_t2p, na.rm = TRUE),
    mean_pv_upper_ci_t2p = mean(pv_upper_ci_t2p, na.rm = TRUE),
    mean_v2_lower_ci_t2p = mean(v2_lower_ci_t2p, na.rm = TRUE),
    mean_v2_upper_ci_t2p = mean(v2_upper_ci_t2p, na.rm = TRUE),
    
    # Summaries for first_inf
    mean_bird_cv_t2i = mean(bird_cv_t2i, na.rm = TRUE),
    mean_bird_pv_t2i = mean(bird_pv_t2i, na.rm = TRUE),
    mean_bird_v2_t2i = mean(bird_v2_t2i, na.rm = TRUE),
    mean_cv_lower_ci_t2i = mean(cv_lower_ci_t2i, na.rm = TRUE),
    mean_cv_upper_ci_t2i = mean(cv_upper_ci_t2i, na.rm = TRUE),
    mean_pv_lower_ci_t2i = mean(pv_lower_ci_t2i, na.rm = TRUE),
    mean_pv_upper_ci_t2i = mean(pv_upper_ci_t2i, na.rm = TRUE),
    mean_v2_lower_ci_t2i = mean(v2_lower_ci_t2i, na.rm = TRUE),
    mean_v2_upper_ci_t2i = mean(v2_upper_ci_t2i, na.rm = TRUE)
  )

#Variability in when individuals became infected or sick
#There was basically no variability in when individuals became infected, but there was variability in pathology
#because this is using ever_diseased so only show t2p
ggplot(summary_tibble, aes(x = temp, color = temp)) +
  geom_point(aes(y = mean_bird_pv_t2p, shape = "PV (Pathology)"), size = 2) +
  geom_errorbar(aes(y = mean_bird_pv_t2p, ymin = mean_pv_lower_ci_t2p, ymax = mean_pv_upper_ci_t2p), width = 0) +
  geom_point(aes(y = mean_bird_cv_t2p, shape = "CV (Pathology)"), size = 2) +
  geom_point(aes(y = mean_bird_v2_t2p, shape = "V2 (Pathology)"), size = 2) +
  
  # Add corresponding points for first_inf
  # geom_point(aes(y = mean_bird_pv_t2i, shape = "PV (Infection)"), position = position_dodge2(width = 2.5), size = 2) +
  # geom_errorbar(aes(y = mean_bird_pv_t2i, ymin = mean_pv_lower_ci_t2i, ymax = mean_pv_upper_ci_t2i), position = position_dodge2(width = 2.5), width = 0) +
  # geom_point(aes(y = mean_bird_cv_t2i, shape = "CV (Infection)"), size = 2, position = position_dodge2(width = 2.5)) +
  # geom_point(aes(y = mean_bird_v2_t2i, shape = "V2 (Infection)"), size = 2, position = position_dodge2(width = 2.5)) +
  
  labs(
    y = "Metric Value",
    x = "Temperature",
    shape = "Metric Type",
    color = "Temperature",
    title = "Variability in Time to Pathology"
  ) +
  scale_color_manual(values = c("#669BBC", "#C1121F"))+
  scale_shape_manual(name = "Variability Metric",
                     values = c("PV (Pathology)" = 16, "CV (Pathology)" = 2, "V2 (Pathology)" = 0,
                                "PV (Infection)" = 1, "CV (Infection)" = 6, "V2 (Infection)" = 7)) +
  coord_flip()+
  theme_minimal()


summary_tibble.pv <- summary_tibble %>%
  dplyr::select(temp, mean_bird_pv_t2p, mean_bird_pv_t2i, mean_pv_lower_ci_t2i, mean_pv_upper_ci_t2i, mean_pv_lower_ci_t2p, mean_pv_upper_ci_t2p)

# Write the summarized tibble to a CSV file
#write.csv(summary_tibble.pv, "/Users/jesse/Documents/GitHub/Temperature-and-Immunity-2022/Variability/Time_To_Pathology_Diseased_Variability.csv")

#For ever_infected 
firsts.inf.var <- firsts %>%
  group_by(temp) %>%
  filter(ever_infected == 1) %>%
  dplyr::reframe(
    temp = temp,
    band_number = band_number,
    total_eye_score = total_eye_score,
    treatment = treatment,
    
    
    # Calculations for first_inf
    mean_t2i = mean(first_inf, na.rm = TRUE),
    median_t2i = median(first_inf, na.rm = TRUE),
    bird_cv_t2i = calculate_cv(first_inf),
    bird_pv_t2i = calculate_pv(first_inf),
    bird_v2_t2i = calculate_v2(first_inf),
    
    # Bootstrapping for first_inf
    cv_bootstrap_t2i = list(replicate(n_boot, calculate_cv(sample(first_inf, replace = TRUE)))),
    pv_bootstrap_t2i = list(replicate(n_boot, calculate_pv(sample(first_inf, replace = TRUE)))),
    v2_bootstrap_t2i = list(replicate(n_boot, calculate_v2(sample(first_inf, replace = TRUE))))
  ) %>%
  mutate(
    
    # Confidence intervals for first_inf
    cv_lower_ci_t2i = map_dbl(cv_bootstrap_t2i, ~ quantile(.x, 0.025, na.rm = TRUE)),
    cv_upper_ci_t2i = map_dbl(cv_bootstrap_t2i, ~ quantile(.x, 0.975, na.rm = TRUE)),
    pv_lower_ci_t2i = map_dbl(pv_bootstrap_t2i, ~ quantile(.x, 0.025, na.rm = TRUE)),
    pv_upper_ci_t2i = map_dbl(pv_bootstrap_t2i, ~ quantile(.x, 0.975, na.rm = TRUE)),
    v2_lower_ci_t2i = map_dbl(v2_bootstrap_t2i, ~ quantile(.x, 0.025, na.rm = TRUE)),
    v2_upper_ci_t2i = map_dbl(v2_bootstrap_t2i, ~ quantile(.x, 0.975, na.rm = TRUE))
  ) %>%
  dplyr::select(-cv_bootstrap_t2i, -pv_bootstrap_t2i, -v2_bootstrap_t2i) %>%
  ungroup()

summary_tibble.inf <- firsts.inf.var %>%
  group_by(temp) %>%
  reframe(
    
    # Summaries for first_inf
    mean_bird_cv_t2i = mean(bird_cv_t2i, na.rm = TRUE),
    mean_bird_pv_t2i = mean(bird_pv_t2i, na.rm = TRUE),
    mean_bird_v2_t2i = mean(bird_v2_t2i, na.rm = TRUE),
    mean_cv_lower_ci_t2i = mean(cv_lower_ci_t2i, na.rm = TRUE),
    mean_cv_upper_ci_t2i = mean(cv_upper_ci_t2i, na.rm = TRUE),
    mean_pv_lower_ci_t2i = mean(pv_lower_ci_t2i, na.rm = TRUE),
    mean_pv_upper_ci_t2i = mean(pv_upper_ci_t2i, na.rm = TRUE),
    mean_v2_lower_ci_t2i = mean(v2_lower_ci_t2i, na.rm = TRUE),
    mean_v2_upper_ci_t2i = mean(v2_upper_ci_t2i, na.rm = TRUE)
  )

#for ever_infected
ggplot(summary_tibble.inf, aes(x = fct_rev(temp), color = temp)) +
  # geom_point(aes(y = mean_bird_pv_t2p, shape = "PV (Pathology)"), size = 2) +
  # geom_errorbar(aes(y = mean_bird_pv_t2p, ymin = mean_pv_lower_ci_t2p, ymax = mean_pv_upper_ci_t2p), width = 0) +
  # geom_point(aes(y = mean_bird_cv_t2p, shape = "CV (Pathology)"), size = 2) +
  # geom_point(aes(y = mean_bird_v2_t2p, shape = "V2 (Pathology)"), size = 2) +
  
  # Add corresponding points for first_inf
  geom_point(aes(y = mean_bird_pv_t2i, shape = "PV (Infection)"), size = 2) +
  geom_errorbar(aes(y = mean_bird_pv_t2i, ymin = mean_pv_lower_ci_t2i, ymax = mean_pv_upper_ci_t2i), width = 0) +
  geom_point(aes(y = mean_bird_cv_t2i, shape = "CV (Infection)"), size = 2) +
  geom_point(aes(y = mean_bird_v2_t2i, shape = "V2 (Infection)"), size = 2) +
  
  labs(
    y = "Metric Value",
    x = "Temperature",
    shape = "Metric Type",
    color = "Temperature",
    title= "Variability in time to Infection"
  ) +
  scale_color_manual(values = c("#669BBC","#C1121F"))+
  scale_shape_manual(name = "Variability Metric",
                     values = c("PV (Pathology)" = 16, "CV (Pathology)" = 2, "V2 (Pathology)" = 0,
                                "PV (Infection)" = 1, "CV (Infection)" = 6, "V2 (Infection)" = 7)) +
  coord_flip()+
  theme_minimal()

#Questions
# 1) Does ambient temperature predict whether birds get infected upon inoculation?
# 2) Do cold temperatures affect pathology?
# 3) Do cold temperatures affect time to pathology?
# 4) Do cold temperatures affect overall infection length?

#adjust reference categories
ti.inoc$treatment <- factor(ti.inoc$treatment,levels = c("Sham", "Infected"))
ti.inoc$temp <- factor(ti.inoc$temp,levels = c("Cold", "Warm"))
#1) Binomial logistic regression
lm1 <- glm(ever_diseased ~ temp, data= ti.inoc %>% filter(dpi ==35), family=binomial)
summary(lm1)
car::Anova(lm1, type="III")
simulateResiduals(lm1, plot=T)

lm1.5 <- glm(ever_infected ~ temp, data= ti.inoc %>% filter(dpi ==35), family=binomial)
summary(lm1.5)

#Temperature does not predict whether birds that are inoculated develop pathology or become infected or not,
#but more birds developed pathology in the warm temperature.

# Filter data
data_filtered <- ti.inoc %>% filter(dpi == 35)

# Calculate proportions
proportions <- data_filtered %>%
  group_by(temp) %>%
  dplyr::summarize(
    diseased = mean(ever_diseased, na.rm = TRUE),
    diseased.n = sum(ever_diseased, na.rm=TRUE),
    infected = mean(ever_infected, na.rm=TRUE),
    infected.n = sum(ever_infected, na.rm=TRUE),
    count = n()
  )

# Plot proportion diseased
ggplot(proportions, aes(x = temp, y = diseased)) +
  geom_bar(aes(fill=temp), stat = "identity", width = 0.6) +
  #geom_text(aes(label = paste0(round(diseased * 100, 1), "%")), vjust = -0.5) +
  geom_text(aes(label = paste0(diseased.n, " / ", count)), vjust = -0.5) +
  scale_fill_manual(values=c("#669BBC", "#C1121F"))+
  labs(
    title = "Proportion of Birds Diseased by Temperature",
    x = "Temperature",
    y = "Proportion Diseased",
    fill = "Temperature"
  ) +
  theme_minimal()

ggplot(count_data, aes(x = event, y = count, fill = temp)) +
  geom_bar(stat = "identity", position = "dodge", width=0.5) +
  labs(
    x = "Definition",
    y = "Number of Birds",
    fill = "Temperature",
    title = "Number of Birds Diseased and Infected by Temperature"
  ) +
  scale_fill_manual(values = c("Warm" = "#C1121F", "Cold" = "#669BBC")) +
  theme_minimal(base_size=16)

#proportion infected
ggplot(proportions, aes(x = temp, y = infected)) +
  geom_bar(aes(fill=temp), stat = "identity", width = 0.6) +
  #geom_text(aes(label = paste0(round(infected * 100, 1), "%")), vjust = -0.5) +
  geom_text(aes(label = paste0(infected.n, " / ", count)), vjust = -0.5) +
  scale_fill_manual(values=c("#669BBC", "#C1121F"))+
  labs(
    title = "Proportion of Birds Infected by Temperature",
    x = "Temperature",
    y = "Proportion Infected",
    fill = "Temperature"
  ) +
  theme_minimal()
library(DHARMa)
library(glmmTMB)
library(emmeans)
#2) Tough one; Use a zero-inflated gamma model allowing for temperature to dictate the probability of zero
#dpi is continuous here which may be an issue, but does not converge with it as a factor
lm2 <- glmmTMB(total_eye_score ~ temp*dpi + (1|band_number),
               data=ti.inoc, 
               ziformula = ~ temp,
               family = ziGamma(link = "log"))
summary(lm2)
simulateResiduals(lm2, plot=T)
car::Anova(lm2, type = "III")
emmeans(lm2, pairwise ~ temp*dpi)
coef(summary(lm2))

#Temperature does not predict pathology

#what about max eye score
ti.inoc.tes <- ti.inoc %>%
  group_by(band_number)%>%
  reframe(max_eye_score = max(total_eye_score),
          treatment = treatment,
          temp = temp,
          band_number = band_number,
          dpi = dpi)

ggplot(ti.inoc.tes %>% filter(dpi == 28), aes(x=temp, y=max_eye_score, color=temp))+
  geom_jitter(width=0.1, height =0, alpha=0.5, size=1.5)+
  scale_color_manual(values = c("#669BBC", "#C1121F"))+
  stat_summary(aes(groups=temp, color=temp), geom="point", fun="mean", size=2, shape=15)

lm2.5 <- glm(max_eye_score ~ temp, data=ti.inoc.tes %>% filter(dpi ==28))
summary(lm2.5)

# 3) Lm testing whether time to first pathology is predicted by temperature. Subset to only individuals that became diseased.
ti.inoc.d <- ti.inoc %>%
  filter(ever_diseased == 1 & dpi == 35)

#Very few birds got sick
ti.inoc.d %>%
  dplyr::select(dpi, temp, ever_diseased)%>%
  tbl_summary(
    by=temp
  )%>%
  modify_header(
    label ~ "**Diseased Birds**"
  )

firsts$temp <- factor(firsts$temp,levels = c("Warm", "Cold"))
hist(firsts$first_path)

lm3 <- lm(first_path ~ temp, data=firsts %>% filter(ever_diseased == 1 & dpi ==3))
summary(lm3)
simulateResiduals(lm3, plot=T)

lm3.5 <- lm(first_inf ~ temp, data=firsts %>% filter(ever_infected == 1 & dpi ==3))
summary(lm3.5)

t.test(first_path ~ temp, data= firsts %>% filter(ever_diseased ==1 & dpi ==3))
t.test(first_inf ~ temp, data= firsts %>% filter(ever_infected ==1 & dpi ==3))

wilcox.test(first_path ~ temp, data = firsts %>% filter(ever_diseased == 1 & dpi ==3))

#Temperature does not predict time to first pathology, although trending towards cold temperatures delaying onset of pathology
set.seed(3)
ggplot(firsts %>% filter(dpi == 3), aes(x=temp, y=first_path, fill=temp))+
  #geom_jitter(height=0, width=0.5, alpha=0.75, size=2)+
  geom_dotplot(binaxis = "y", stackdir = "center", aes(fill=temp), dotsize=0.75)+
  stat_summary(fun=mean, geom="point", shape=17)+
  #stat_summary(fun=median, geom="point", shape=2)+
  scale_fill_manual(values=c("#C1121F","#669BBC"))+
  labs(y="Days to First Pathology", x="Temperature", fill="Temperature")+
  coord_flip()


####What about pathogen load? DPI 3, 9, 14, 28
summary(ti)

ti.q <- ti.ep %>%
  filter(dpi %in% c(3, 9, 14, 28))
ti.q$quantity1 <- ti.q$quantity+1

ggplot(ti.q, aes(x=dpi, y=quantity1, color=temp))+
  geom_point()+
  geom_line(aes(group = band_number), alpha=0.25)+
  scale_color_manual(values=c("#669BBC","#C1121F"))+
  stat_summary(aes(group=temp), fun="mean", geom="line", size =1, lty="dashed")+
  #scale_y_log10()+
  geom_hline(yintercept = 100)+
  facet_wrap(~treatment)

#max pathogen load
ti.q.max <- ti.q %>%
  group_by(band_number)%>%
  reframe(max_quantity1 = max(quantity1),
          treatment = treatment,
          temp = temp,
          band_number = band_number,
          dpi = dpi,
          ever_infected = ever_infected)

ggplot(ti.q.max %>% filter(dpi ==28 & ever_infected == 1), aes(x=temp, y=max_quantity1, color=temp))+
  geom_jitter(height=0, width=0.1, alpha=0.5)+
  scale_color_manual(values=c("#669BBC","#C1121F"))+
  stat_summary(aes(group=temp), fun="mean", geom="point", size =2)

lm4 <- glm(max_quantity1 ~ temp, data=ti.q.max)
summary(lm4)

first_sick <- ti.e %>%
  filter(ever_infected == 1) %>%
  group_by(band_number) %>%
  mutate(
    first_path = if (any(total_eye_score > 0, na.rm = TRUE)) {
      min(dpi[total_eye_score > 0], na.rm = TRUE)
    } else {
      NA
    }
  ) %>%
  ungroup()

#Average time to pathology for birds that are infected and median time. These both may be
#misleading because eye score was only taken DPI 9 and 14. It is likely birds developed scores earlier
avg_time <- firsts %>% 
  group_by(treatment, temp)%>%
  summarise(avg_to_path = mean(first_path, na.rm=T),
            med_to_path = median(first_path, na.rm=T))

ti.ep <- merge(ti.e, first_pathology, by=c("band_number", "dpi"))

colnames(ti.ep)

colnames(ti.ep) <- gsub("\\.x$", "", colnames(ti.ep))
ti.ep <- ti.ep[, !grepl("\\.y$", colnames(ti.ep))]
colnames(ti.ep)  

