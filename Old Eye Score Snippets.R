##Old eye score models for TI22
####Eye Score: Ordinal Mixed Effects Model####
#Sham birds did not get infected therefore did not have eyescores. Remove from analysis so that model converges.
#No eye score on DPI 3 so also remove. Keep DPI 0 as reference category.
ti.mod <- ti %>%
  dplyr::select(date, dpi, band_number, bird_ID, treatment, temp, total_eye_score)%>%
  filter(treatment != "Sham" & dpi %in% c(0, 7,9,14,18,21,24,28))%>%
  mutate(temp = relevel(factor(temp), ref = "Warm"),
         dpi_f  = factor(dpi, levels = c(0,7,9,14,18,21,24,28)))

#We need to define thresholds based off of the observed data. The technical bounds are [0,6], but we only have [0,4]
levels_obs <- sort(unique(ti.mod$total_eye_score))

#turn observed numeric scores into an ordered factor with factors
ti.mod$eye_score_ord <- ordered(
  ti.mod$total_eye_score,
  levels=levels_obs
)

with(ti.mod, xtabs(~ eye_score_ord + temp + dpi_f))

ti.mod <- ti.mod %>%
  mutate(eye5 = fct_collapse(eye_score_ord,
                             "0"="0",
                             "1"=c("0.5","1"),
                             "2"=c("1.5","2"),
                             "3"=c("2.5","3"),
                             "4"=c("3.5","4")))


ti.mod <- ti.mod %>%
  mutate(
    dpi_f = relevel(dpi_f, ref = "0"),
    temp  = relevel(temp,  ref = "Warm")
  )

#Model Comparisons
clm.a <- clmm(eye_score_ord ~ temp * dpi_f + (1|band_number),
              data=ti.mod,
              link = "logit")

clm.b <- clmm(eye5 ~ temp * dpi_f + (1|band_number),
              data=ti.mod,
              link = "logit")

clm.c <- clmm(eye_score_ord ~ temp + dpi_f + (1|band_number),
              data=ti.mod,
              link = "logit")

clm.d <- clmm(eye5 ~ temp + dpi_f + (1|band_number),
              data=ti.mod,
              link = "logit")

summary(clm.a)
summary(clm.b)
summary(clm.c)
summary(clm.d)

aictab(cand.set=list(clm.a, clm.c),
       modnames= c("clm.a", "clm.c"))

aictab(cand.set=list(clm.b, clm.d),
       modnames= c("clm.b", "clm.d"))

#Does collapsing eye score into 5 categories hurt model?

#can't directly compare, but lower AIC = good
AIC(clm.c, clm.d)

#Lower logLikelihood = good
logLik(clm.c)
logLik(clm.d)

#R2 same if not slightly better = good
r2_nakagawa(clm.c)
r2_nakagawa(clm.d)

#predictions look the same
emmip(clm.c, temp ~ dpi_f)
emmip(clm.d, temp ~ dpi_f)

#use eye5 in final model
clm1 <- clmm(eye5 ~ temp + dpi_f + (1|band_number),
             data=ti.mod,
             link = "logit")

summary(clm1)

ti.mod <- ti.mod %>%
  mutate(
    eye5_num = as.numeric(as.character(eye5)),               # 0..4
    dpi_num  = as.numeric(as.character(dpi_f))               # e.g., 0,7,9,14,...
  )

# 2) Model predictions: expected category (mean class) + 95% CIs
emm_mean <- emmeans(clm1, ~ temp * dpi_f,
                    type = "response", mode = "mean.class") %>%
  as.data.frame() %>%
  mutate(
    yhat = estimate - 1,
    ymin = asymp.LCL - 1,
    ymax = asymp.UCL - 1,
    dpi_num = as.numeric(as.character(dpi_f))
  )

ggplot() +
  # raw observations
  geom_jitter(
    data = ti.mod,
    aes(x = dpi_num, y = eye5_num, color = temp),
    width = 0.6, height = 0.08, alpha = 0.35, size = 1.7
  ) +
  # model mean class with CI
  geom_errorbar(
    data = emm_mean,
    aes(x = dpi_num, ymin = ymin, ymax = ymax, color = temp),
    width = 0.8, size = 0.5
  ) +
  geom_line(
    data = emm_mean,
    aes(x = dpi_num, y = yhat, color = temp, group = temp),
    size = 0.7
  ) +
  geom_point(
    data = emm_mean,
    aes(x = dpi_num, y = yhat, color = temp),
    size = 2.5
  ) +
  scale_y_continuous(breaks = 0:4, limits = c(-1, 4.1)) +
  labs(x = "Days Post Inoculation", y = "Total Eye Score",
       color = "Temperature") +
  theme_minimal(base_size = 12)


ggplot(ti.inoc, aes(x=dpi, y=total_eye_score, color=temp))+
  geom_jitter(width=0.5, height=0.1, alpha=0.25)+
  #geom_line(aes(group = band_number), alpha=0.25)+
  stat_summary(aes(x=dpi, y=total_eye_score, color = temp), geom = "point", size=3)+
  stat_summary(aes(x=dpi, y=total_eye_score, color = temp), geom = "errorbar", width=0.25)+
  stat_summary(aes(x=dpi, y=total_eye_score, color = temp), geom = "line", width=0.25, alpha=0.5)+
  scale_color_manual(values = temp_colors)+
  labs(x="Days Post Inoculation", y="Total Eye Score", color="Temperature")




#####Eye score: GAM####
mod.e <- ti.inoc %>%
  dplyr::filter(treatment == "Inoculated" & dpi > 0)%>%
  dplyr::mutate(total_eye_score1 = total_eye_score+0.001, #add small constant for gamma model
                dpi = as.numeric(dpi),
                temp = factor(temp),
                band_number = factor(band_number)
  )

#Use GAM to ask how temperature affects eye score
# Fit a GAM with smooth function of dpi, separate smooths for each temperature
#7 dpi, so need k to be less than 7

##Gamma Model - Old
# gam1 <- gam(
#   total_eye_score1 ~ temp + s(dpi, by = temp, k = 5, bs = "cr") +
#     s(band_number, bs="re"),
#   data   = mod.e,
#   family = Gamma(link = "log"),
#   method = "REML",
#   select = TRUE #Select = TRUE tells the model to add extra penalty to each smooth term; if smoothing doesn't help, the model drops it
# )
# 
# 
# summary(gam1)
# gam.check(gam1)   # look at k-index & residuals
# plot(gam1, pages = 1, shade = TRUE)
# 
# #Test whether temperature smooth improves fit - null model
# gam_null <- gam(
#   total_eye_score1 ~ temp + s(dpi, k = 5, bs = "cr") + s(band_number, bs="re"),
#   data   = mod.e,
#   family = Gamma(link="log"),
#   method = "REML",
#   select = TRUE
# )
# 
# #Does including temperature  smooth improve model fit?
# anova(gam_null, gam1, test = "Chisq")
# 
# newdat <- expand.grid(
#   dpi  = seq(min(mod.e$dpi, na.rm=TRUE), max(mod.e$dpi, na.rm=TRUE), length.out = 200),
#   temp = levels(mod.e$temp)
# )
# 
# # dummy level (any level from the fit is fine)
# newdat$band_number <- levels(mod.e$band_number)[1]
# 
# # exclude the RE smooth when predicting
# pr <- predict(gam1, newdat, type = "link", se.fit = TRUE, 
#               exclude = "s(band_number)")
# 
# pr   <- predict(gam1, newdat, type = "link", se.fit = TRUE)
# newdat$fit <- family(gam1)$linkinv(pr$fit)
# newdat$lwr <- family(gam1)$linkinv(pr$fit - 1.96*pr$se.fit)
# newdat$upr <- family(gam1)$linkinv(pr$fit + 1.96*pr$se.fit)
# 
# ggplot(newdat, aes(dpi, fit, color = temp, fill = temp)) +
#   geom_jitter(data=mod.e, aes(x=dpi, y=total_eye_score, color= temp), height=0.05, width=0)+
#   geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.15, linewidth = 0) +
#   geom_line(linewidth = 1.1) +
#   scale_fill_manual(values=c(temp_colors))+
#   scale_color_manual(values=c(temp_colors))+
#   labs(x = "Days Post Infection", y = "Total Eye Score", fill= "Temperature", color = "Temperature")+
#   theme_bw()

#tweedie GAM to handle zero inflation instead of adding small constant
gam1 <- gam(
  total_eye_score ~ temp + s(dpi, by = temp, k = 5, bs = "cr") +
    s(band_number, bs = "re"),
  data = mod.e,
  family = tw(link = "log"),
  method = "REML"
)

summary(gam1)
gam.check(gam1)   # look at k-index & residuals
plot(gam1, pages = 1, shade = TRUE)

#Test whether temperature smooth improves fit - null model
gam_null <- gam(
  total_eye_score ~ temp + s(dpi, k = 5, bs = "cr") + s(band_number, bs="re"),
  data   = mod.e,
  family = tw(link = "log"),
  method = "REML",
  select = TRUE
)

#Does including temperature  smooth improve model fit?
anova(gam_null, gam1, test = "Chisq")

newdat <- expand.grid(
  dpi  = seq(min(mod.e$dpi, na.rm=TRUE), max(mod.e$dpi, na.rm=TRUE), length.out = 200),
  temp = levels(mod.e$temp)
)

# dummy level (any level from the fit is fine)
newdat$band_number <- levels(mod.e$band_number)[1]

# exclude the RE smooth when predicting
pr <- predict(gam1, newdat, type = "link", se.fit = TRUE, 
              exclude = "s(band_number)")

pr   <- predict(gam1, newdat, type = "link", se.fit = TRUE)
newdat$fit <- family(gam1)$linkinv(pr$fit)
newdat$lwr <- family(gam1)$linkinv(pr$fit - 1.96*pr$se.fit)
newdat$upr <- family(gam1)$linkinv(pr$fit + 1.96*pr$se.fit)

ggplot(newdat, aes(dpi, fit, color = temp, fill = temp)) +
  geom_jitter(data=mod.e, aes(x=dpi, y=total_eye_score, color= temp), height=0.05, width=0.5, size=2, alpha=0.5)+
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.15, linewidth = 0) +
  geom_line(linewidth = 1.1) +
  scale_fill_manual(values=c(temp_colors))+
  scale_color_manual(values=c(temp_colors))+
  labs(x = "Days Post Infection", y = "Total Eye Score", fill= "Temperature", color = "Temperature")+
  theme_bw()

#Save statistics
s <- summary(gam1)

# Parametric (fixed) terms
param_df <- as_tibble(s$p.table, rownames = "term") %>%
  rename(
    estimate  = `Estimate`,
    std_error = `Std. Error`
  ) %>%
  mutate(
    statistic = `t value`,
    p_value   = `Pr(>|t|)`,
    model     = "Tweedie_GAM"
  ) %>%
  dplyr::select(model, term, estimate, std_error, statistic, p_value)

# Smooth terms 
smooth_df <- as_tibble(s$s.table, rownames = "smooth") %>%
  mutate(
    edf    = edf,
    ref_df = `Ref.df`,
    stat   = `F`,
    p_value = `p-value`,
    model  = "Tweedie_GAM"
  ) %>%
  dplyr::select(model, smooth, edf, ref_df, stat, p_value)

# Model-level fit metrics
fit_df <- tibble(
  model         = "Tweedie_GAM",
  family        = gam1$family$family,
  link          = gam1$family$link,
  shape_p       = if (!is.null(gam1$family$variance)) s$family$p else NA_real_,  # Tweedie p
  AIC           = AIC(gam1),
  dev_expl_frac = s$dev.expl,                      # 0–1
  dev_expl_pct  = round(100 * s$dev.expl, 1),
  r2_adj        = s$r.sq %||% NA_real_,
  scale_est     = s$scale,
  n             = s$n
)

cmp_tbl <- anova(gam_null, gam1, test = "Chisq") %>%
  as_tibble() %>%
  mutate(comparison = "Tweedie: null vs temp-specific")

#write to excel
# write_xlsx(
#   list(
#     tweedie_param   = param_df,
#     tweedie_smooth  = smooth_df,
#     tweedie_fit     = fit_df,
#     tweedie_compare = cmp_tbl
#   ),
#   path = "/Users/jesse/Documents/GitHub/Temperature-and-Immunity-2022/Results/Eye_Score_tweedie_gam_Results.xlsx"
# )



l1_poly <- glmmTMB(total_eye_score ~ temp * poly(dpi, 2) + (1|band_number),
                   data=ti.inoc, 
                   ziformula = ~ temp,
                   family = ziGamma(link = "log"))

l2_poly <- glmmTMB(
  total_eye_score ~ sex + temp * poly(dpi, 2) + (1 | band_number),
  data      = ti.inoc,
  ziformula = ~ sex,
  family    = ziGamma(link = "log")
)

l3_poly <- glmmTMB(total_eye_score ~ temp* poly(dpi, 2) + sex + (1|band_number),
                   data=ti.inoc, 
                   ziformula = ~ temp,
                   family = ziGamma(link = "log"))

l4_poly <- glmmTMB(total_eye_score ~ temp*poly(dpi, 2) * sex + (1|band_number),
                   data=ti.inoc, 
                   ziformula = ~ temp,
                   family = ziGamma(link = "log"))

l5_poly <- glmmTMB(total_eye_score ~ sex + temp * poly(dpi, 2) + (1|band_number),
                   data=ti.inoc, 
                   ziformula = ~ temp * poly(dpi, 2) + sex, #zero inflation: 
                   family = ziGamma(link = "log"))

l6_poly <- glmmTMB(total_eye_score ~ sex * temp * poly(dpi, 2) + (1|band_number),
                   data=ti.inoc, 
                   ziformula = ~ temp * poly(dpi, 2) + sex, #zero inflation: 
                   family = ziGamma(link = "log"))

l7_poly <- glmmTMB(total_eye_score ~ sex + temp * poly(dpi, 2) + (1|band_number),
                   data=ti.inoc, 
                   ziformula = ~ temp + poly(dpi, 2) + sex, #zero inflation: 
                   family = ziGamma(link = "log"))

l8_poly <- glmmTMB(total_eye_score ~ sex * temp * poly(dpi, 2) + (1|band_number),
                   data=ti.inoc, 
                   ziformula = ~ temp + poly(dpi, 2) + sex, #zero inflation: 
                   family = ziGamma(link = "log"))

l9_poly <- glmmTMB(total_eye_score ~  temp * poly(dpi, 2) + (1|band_number),
                   data=ti.inoc, 
                   ziformula = ~ temp + poly(dpi, 2) + sex, #zero inflation: 
                   family = ziGamma(link = "log"))

l10_poly <- glmmTMB(total_eye_score ~  temp + poly(dpi, 2) + (1|band_number),
                    data=ti.inoc, 
                    ziformula = ~ temp + poly(dpi, 2) + sex, #zero inflation: 
                    family = ziGamma(link = "log"))

aictab(cand.set=list(l1_poly, l2_poly, l3_poly, l4_poly, l5_poly, l6_poly, l7_poly, l8_poly, l9_poly, l10_poly, l2), 
       modnames=c("l1_poly", "l2_poly", "l3_poly", "l4_poly", "l5_poly", "l6_poly", "l7_poly", "l8_poly", "l9_poly", "l10_poly", "l2"))

lp1 <- glmmTMB(total_eye_score ~ sex + temp * poly(dpi, 2) + (1|band_number),
               data=ti.inoc, 
               ziformula = ~ temp + poly(dpi, 2) + sex, #zero inflation: 
               family = ziGamma(link = "log"))

lp2 <- glmmTMB(total_eye_score ~ sex + temp * poly(dpi, 2) + (1|band_number),
               data=ti.inoc, 
               ziformula = ~ temp + poly(dpi, 2), #zero inflation: 
               family = ziGamma(link = "log"))

lp3 <- glmmTMB(total_eye_score ~  temp * poly(dpi, 2) + (1|band_number),
               data=ti.inoc, 
               ziformula = ~ temp + poly(dpi, 2) , #zero inflation: 
               family = ziGamma(link = "log"))
lp4 <- glmmTMB(total_eye_score ~  temp * poly(dpi, 2) + (1|band_number),
               data=ti.inoc, 
               ziformula = ~ temp + poly(dpi, 2) + sex, #zero inflation: 
               family = ziGamma(link = "log"))

AIC(lp1, lp2, lp3, lp4) #keep sex

#####Max Eye Score####
ti.mod.m <- ti.inoc %>%
  filter(dpi > 0)%>%
  group_by(band_number, temp, treatment, groups, ever_infected, ever_diseased, sex) %>%
  reframe(
    max_tes = max(total_eye_score))

ggplot(ti.mod.m, aes(x=max_tes, fill=temp))+
  geom_histogram(aes(groups=temp), position = "dodge")+
  facet_wrap(~ever_infected)

#Do max eye scores differ between temperatures in 1) birds that were inoculated with MG, 2) birds that became infected, or 3) birds that developed disease?
##No for all

#Inoculated
wilcox.test(max_tes ~ temp, data = ti.mod.m %>% filter(treatment == "Inoculated"))
wilcox.test(max_tes ~ sex, data = ti.mod.m %>% filter(treatment == "Inoculated"))

#Infected
wilcox.test(max_tes ~ temp, data = ti.mod.m %>% filter(ever_infected == 1))
wilcox.test(max_tes ~ sex, data = ti.mod.m %>% filter(ever_infected == 1))

#Diseased
wilcox.test(max_tes ~ temp, data = ti.mod.m %>% filter(ever_diseased == 1))
wilcox.test(max_tes ~ sex, data = ti.mod.m %>% filter(ever_diseased == 1))

# data:  max_tes by temp
# W = 108.5, p-value = 0.6225

ggplot(ti.mod.m %>% filter(treatment == "Inoculated"), aes(x = temp, y = max_tes, color = temp, fill=temp)) +
  geom_jitter(width=0.1, height=0, size = 2.5, alpha=0.75) +
  #geom_boxplot(alpha=0.3, color="black", width=0.05)+
  stat_summary(geom="point", fun="mean", shape = 16, size=3, color="black")+
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0., color = "black") +
  # geom_point(data=dat.new, aes(x=temp, y=yhat), size=3)+
  # geom_errorbar(data=dat.new, aes(x=temp, y=yhat, ymin=Lower, ymax = Upper), width=0.05)+
  scale_color_manual(values = temp_colors) +
  scale_fill_manual(values = temp_colors) +
  labs(x = "Group", y = "Max Eyescore", color="Temperature", fill="Temperature")




#Use GAM to ask how temperature affects pathogen load
# Fit a GAM with smooth function of dpi, separate smooths for each temperature
#7 dpi, so need k to be less than 7

#Model Selection
#Use ML (maximum likelihood for model selection to compare fixed effects) to see if adding sex helps
#asking whether allowing separate smooths for temp, dpi, or sex improves model fit
g1 <- gam(
  quantity1 ~  s(dpi, k = 5, bs = "cr") +
    s(band_number, bs="re"),
  data   = dat.q,
  family = Gamma(link = "log"),
  method = "ML",
  select = TRUE #Select = TRUE tells the model to add extra penalty to each smooth term; if smoothing doesn't help, the model drops it
)

g2 <- gam(
  quantity1 ~ temp + s(dpi, by = temp, k = 5, bs = "cr") +
    s(band_number, bs="re"),
  data   = dat.q,
  family = Gamma(link = "log"),
  method = "ML",
  select = TRUE #Select = TRUE tells the model to add extra penalty to each smooth term; if smoothing doesn't help, the model drops it
)

g3 <- gam(
  quantity1 ~ temp * sex + s(dpi, by = interaction(temp, sex), k = 5, bs = "cr") +
    s(band_number, bs="re"),
  data   = dat.q,
  family = Gamma(link = "log"),
  method = "ML",
  select = TRUE #Select = TRUE tells the model to add extra penalty to each smooth term; if smoothing doesn't help, the model drops it
)

null <- gam(
  quantity1 ~ 1 + s(dpi, k = 5, bs = "cr") +
    s(band_number, bs="re"),
  data   = dat.q,
  family = Gamma(link = "log"),
  method = "ML",
  select = TRUE #Select = TRUE tells the model to add extra penalty to each smooth term; if smoothing doesn't help, the model drops it
)

aictab(cand.set = list(g1, g2, g3, null),
       modnames = c("g1", "g2", "g3", "null"))


AIC(g1, g2, g3, gam2, null)

gam2 <- gam(
  quantity1 ~ temp + s(dpi, by = temp, k = 5, bs = "cr") +
    s(band_number, bs="re"),
  data   = dat.q,
  family = Gamma(link = "log"),
  method = "REML",
  select = TRUE #Select = TRUE tells the model to add extra penalty to each smooth term; if smoothing doesn't help, the model drops it
)


summary(gam2)
summary(g3)
gam.check(gam2)   # look at k-index & residuals
plot(gam2, pages = 1, shade = TRUE)

newdat <- expand.grid(
  dpi  = seq(min(dat.q$dpi, na.rm=TRUE), max(dat.q$dpi, na.rm=TRUE), length.out = 200),
  temp = levels(dat.q$temp)
)

# dummy level (any level from the fit is fine)
newdat$band_number <- levels(dat.q$band_number)[1]

# exclude the RE smooth when predicting
pr <- predict(gam2, newdat, type = "link", se.fit = TRUE, 
              exclude = "s(band_number)")

pr   <- predict(gam2, newdat, type = "link", se.fit = TRUE)
newdat$fit <- family(gam2)$linkinv(pr$fit)
newdat$lwr <- family(gam2)$linkinv(pr$fit - 1.96*pr$se.fit)
newdat$upr <- family(gam2)$linkinv(pr$fit + 1.96*pr$se.fit)

ggplot(newdat, aes(dpi, fit, color = temp, fill = temp)) +
  geom_jitter(data=dat.q, aes(x=dpi, y=quantity1, color= temp),height=0.05, width=0.5, size=2, alpha=0.5)+
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.15, linewidth = 0) +
  geom_line(linewidth = 1.1) +
  scale_fill_manual(values=c(temp_colors))+
  scale_color_manual(values=c(temp_colors))+
  labs(x = "Days Post Infection", y = "Pathogen Load", fill= "Temperature", color = "Temperature")+
  scale_y_log10()
#Test whether temperature smooth improves fit - null model
gam_null <- gam(
  quantity1 ~ temp + s(dpi, k = 5, bs = "cr") + s(band_number, bs="re"),
  data   = dat.q,
  family = Gamma(link="log"),
  method = "REML",
  select = TRUE
)

#Does including temperature  smooth improve model fit?
anova(gam_null, gam2, test = "Chisq")


s <- summary(gam2)

# Parametric (fixed) terms
param_df <- as_tibble(s$p.table, rownames = "term") %>%
  rename(
    estimate  = `Estimate`,
    std_error = `Std. Error`
  ) %>%
  mutate(
    statistic = `t value`,
    p_value   = `Pr(>|t|)`,
    model     = "Gamma_GAM"
  ) %>%
  dplyr::select(model, term, estimate, std_error, statistic, p_value)

# Smooth terms 
smooth_df <- as_tibble(s$s.table, rownames = "smooth") %>%
  mutate(
    edf    = edf,
    ref_df = `Ref.df`,
    stat   = `F`,
    p_value = `p-value`,
    model  = "Gamma_GAM"
  ) %>%
  dplyr::select(model, smooth, edf, ref_df, stat, p_value)

# Model-level fit metrics
fit_df <- tibble(
  model         = "Gamma_GAM",
  family        = gam2$family$family,
  link          = gam2$family$link,
  shape_p       = if (!is.null(gam2$family$variance)) s$family$p else NA_real_,  # Tweedie p
  AIC           = AIC(gam2),
  dev_expl_frac = s$dev.expl,                      # 0–1
  dev_expl_pct  = round(100 * s$dev.expl, 1),
  r2_adj        = s$r.sq %||% NA_real_,
  scale_est     = s$scale,
  n             = s$n
)

cmp_tbl <- anova(gam_null, gam2, test = "Chisq") %>%
  as_tibble() %>%
  mutate(comparison = "Gamma: null vs temp-specific")

#write to excel
# write_xlsx(
#   list(
#     gamma_param   = param_df,
#     gamma_smooth  = smooth_df,
#     gamma_fit     = fit_df,
#     gamma_compare = cmp_tbl
#   ),
#   path = "/Users/jesse/Documents/GitHub/Temperature-and-Immunity-2022/Results/Path_Load_gamma_gam_Results.xlsx"
# )

#Model Selection
p1 <- glmmTMB(log10_quantity1 ~ temp + sex + dpi.f + (1|band_number), data=dat.q, family=Gamma(link="log"))
p2 <- glmmTMB(log10_quantity1 ~ temp * sex + dpi.f + (1|band_number), data=dat.q, family=Gamma(link="log"))
p3 <- glmmTMB(log10_quantity1 ~ temp * sex * dpi.f + (1|band_number), data=dat.q, family=Gamma(link="log"))
p4 <- glmmTMB(log10_quantity1 ~ temp + sex * dpi.f + (1|band_number), data=dat.q, family=Gamma(link="log"))
p4.5 <- glmmTMB(log10_quantity1 ~ temp + sex * dpi.f * ever_diseased + (1|band_number), data=dat.q, family=Gamma(link="log"))
p5 <- glmmTMB(log10_quantity1 ~ temp * dpi.f + (1|band_number), data=dat.q, family=Gamma(link="log"))
null <- glmmTMB(log10_quantity1 ~ 1 + (1|band_number), data=dat.q, family=Gamma(link="log"))

aictab(cand.set = list(p1, p2, p3, p4, p5, null),
       modnames = c("p1", "p2", "p3", "p4", "p5", "null"))

p1 <- glmmTMB(quantity1 ~ temp + sex + dpi.f + (1|band_number), data=dat.q, family=Gamma(link="log"))
p2 <- glmmTMB(quantity1 ~ temp * sex + dpi.f + (1|band_number), data=dat.q, family=Gamma(link="log"))
p3 <- glmmTMB(quantity1 ~ temp * sex * dpi.f + (1|band_number), data=dat.q, family=Gamma(link="log"))
p4 <- glmmTMB(quantity1 ~ temp + sex * dpi.f + (1|band_number), data=dat.q, family=Gamma(link="log"))
p4.5 <- glmmTMB(quantity1 ~ temp + sex * dpi.f * ever_diseased + (1|band_number), data=dat.q, family=Gamma(link="log"))
p5 <- glmmTMB(quantity1 ~ temp * dpi.f + (1|band_number), data=dat.q, family=Gamma(link="log"))
null <- glmmTMB(quantity1 ~ 1 + (1|band_number), data=dat.q, family=Gamma(link="log"))

aictab(cand.set = list(p1, p2, p3, p4, p5, null),
       modnames = c("p1", "p2", "p3", "p4", "p5", "null"))

simulateResiduals(p4, plot=T)
hist(resid(p4))
summary(p4)
car::Anova(p4, type= "III")


AIC(p4, p4.25, p4.5)
simulateResiduals(p4.25, plot=T)

simulateResiduals(p3, plot=T)
summary(p3)
car::Anova(p3, type = "III")

glm.p<- glmmTMB(log10_quantity1 ~ temp + sex * dpi.f + (1|band_number), data=dat.q, family=Gamma(link="log"))
simulateResiduals(glm.p, plot=T)
hist(resid(glm.p))
summary(glm.p)
car::Anova(glm.p, type = "III")

#To compare GLM and GAM: same structure other than smooth
m_glm <- glmmTMB(quantity1 ~ temp + sex + dpi.f + (1|band_number),
                 family=Gamma(log), data=dat.q)

m_gam1 <- gam(quantity1 ~ temp + sex + s(dpi,k=5) + s(band_number, bs = "re"),
              family=Gamma(log), data=dat.q, method="ML")


BIC(m_glm, m_gam1, glm.zi.pl, glm.p)
AIC(glm.zi.pl, glm.p)

emm_df <- emmeans(glm.p, ~  sex * dpi.f, type="response", re.form = NA) %>%
  as.data.frame()

dat.q <- dat.q %>%
  mutate(sex_band = interaction(sex, band_number, drop = TRUE))

ggplot(dat.q, aes(x = dpi.f, y = log10_quantity1, color = sex, shape = sex)) +
  # raw data
  geom_point(size = 2, alpha = 0.75, aes(group = sex), position= position_dodge(width=0.4)) +
  # model CI ribbon
  # geom_ribbon(data= emm_df,
  #   aes(x = dpi.f,
  #       ymin = asymp.LCL,
  #       ymax = asymp.UCL,
  #       group = sex,
  #       fill = sex),
  #   inherit.aes = FALSE,
  #   alpha       = 0.15,
  #   linewidth   = 0.5,
  #   linetype = "dashed"
  # ) +
  geom_errorbar(
    data= emm_df,
    aes(
      x=dpi.f,
      ymin = asymp.LCL,
      ymax = asymp.UCL,
      group=sex,
      color=sex),
    inherit.aes=FALSE,
    width=0.05)+
  # model mean line
  geom_line(
    data        = emm_df,
    aes(x = dpi.f,
        y = response,
        group = sex,
        color = sex),
    inherit.aes = FALSE,
    linewidth   = 0.5
  ) +
  # individual bird lines
  geom_line(
    data = dat.q,
    aes(x = dpi.f,
        y = log10_quantity1,
        group = as.factor(sex_band),
        color = sex),
    position = position_dodge(width = 0.4),
    inherit.aes = FALSE,
    linewidth   = 0.1
  ) +
  geom_point(data=emm_df, aes(x= dpi.f, y= response, group = sex, color = sex), size=2)+
  # scale_y_continuous(
  #   trans = "log10"       # puts both points + emmeans on log10 scale
  # ) +
  scale_color_manual(values = (sex_colors))+
  scale_fill_manual(values = (sex_colors))+
  labs(
    x     = "Days Post Infection",
    y     = "Log10(Pathogen Load + 1)",
    fill  = "Sex",
    color = "Sex",
    shape = "Sex"
  ) +
  facet_grid(~temp)+
  theme(strip.text = element_text(size=13))


dodge_width <- 0.4

ggplot(dat.q, aes(x = dpi.f, y = quantity1, color = sex, shape = sex)) +
  
  # raw points
  geom_point(
    aes(group = sex),
    size = 2, alpha = 0.75,
    position = position_dodge(width = dodge_width)
  ) +
  
  # ribbon
  # geom_ribbon(
  #   data = emm_df,
  #   aes(
  #     x = dpi.f,
  #     ymin = asymp.LCL,
  #     ymax = asymp.UCL,
  #     group = interaction(temp, sex),
  #     fill = sex
  #   ),
  #   inherit.aes = FALSE,
  #   alpha = 0.15,
  #   linewidth = 0.5,
  #   linetype = "dashed"
  # ) +
  #Or errorbar
  geom_errorbar(
    data= emm_df,
    aes(
      x=dpi.f,
      ymin = asymp.LCL,
      ymax = asymp.UCL,
      group=sex,
      color=sex),
    inherit.aes=FALSE,
    width=0.05)+
  
  # model line
  geom_line(
    data = emm_df,
    aes(x = dpi.f, y = response, group = sex, color = sex),
    inherit.aes = FALSE,
    linewidth = 0.5
  ) +
  
  # individual bird trajectories — dodged by sex
  # geom_line(
  #   data = dat.q,
  #   aes(
  #     x = dpi.f,
  #     y = quantity1,
  #     group = as.factor(sex_band),    # key change
  #     color = sex
  #   ),
  #   position = position_dodge(width = dodge_width),
  #   inherit.aes = FALSE,
  #   linewidth = 0.15
  # ) +
  
  # emmean points
  geom_point(
    data = emm_df,
    aes(x = dpi.f, y = response, color = sex),
    size = 2,
    inherit.aes = FALSE
  ) +
  
  scale_y_continuous(trans = "log10") +
  scale_color_manual(values = sex_colors) +
  scale_fill_manual(values = sex_colors) +
  labs(
    x     = "Days Post Infection",
    y     = "Log10(Pathogen Load + 1)",
    fill  = "Sex",
    color = "Sex",
    shape = "Sex"
  )


#Test whether temperature smooth improves fit - null model
gam_null <- gam(
  quantity1 ~ temp + s(dpi, k = 5, bs = "cr") + s(band_number, bs="re"),
  data   = dat.q,
  family = Gamma(link="log"),
  method = "REML",
  select = TRUE
)
