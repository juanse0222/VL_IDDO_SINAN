# ══════════════════════════════════════════════════════════════════════════════
# Project : Visceral Leishmaniasis in Brazil, 2007–2023
# Script  : model_test.R
# Purpose : Model exploration and manuscript figure generation.  Covers:
#             - unadjusted Cox model (overall and stratified by HIV status)
#             - adjusted Cox model (delay + covariates → VL mortality)
#             - logistic severity models (bleeding, jaundice, oedema) by HIV
#             - composite complication models (Figure Sx)
#             - symptom probability panel (11 symptoms × HIV, Figures 2 & 3;
#               final clean version lives in request_prabin.R)
#             - GLM logistic model (exploratory, unadjusted and adjusted)
#             - final adjusted Cox models by HIV status (also in request_prabin.R)
# Input   : data base and clean/data_6.rds  (SINAN VL notifications, curated)
# Output  : figures and tables/severe_symptoms_color_graph.png
# Author  : Juan Sebastián Hurtado Zapata
#           IDDO VL Fellowship | TDR/WHO fellow 
# Updated : 2026-07-29
# ══════════════════════════════════════════════════════════════════════════════


# ── 1. PACKAGES ───────────────────────────────────────────────────────────────
# Uncomment to install from a fresh R environment:
# install.packages(c("survival","rms","coxme","mice","Hmisc",
#                    "sjPlot","splines","broom"))

pacman::p_load(
  tidyverse, janitor, rio, here,
  survival,  # Surv(), coxph(), survfit()
  rms,       # cph(), lrm(), rcs(), Predict(), datadist()
  coxme,     # mixed-effects Cox (frailty models, exploratory)
  mice,      # multiple imputation (sensitivity analysis)
  Hmisc,     # supports rms internals
  survminer, # ggsurvplot(), ggsurvplot_facet()
  sjPlot,    # exploratory model tables
  splines,   # ns() — natural splines (alternative to rcs)
  broom,     # tidy(), glance() for model summaries
  ggpubr,    # ggarrange(), annotate_figure()
  grid,      # textGrob() for shared axis labels
  scales     # percent_format()
)


# ── 2. DATA IMPORT AND VARIABLE ENGINEERING ───────────────────────────────────
data_6 <- import(here("data base and clean", "data_6.rds"))

# IDate columns from data.table are not recognised by survival/rms; convert.
data_6 <- data_6 %>%
  dplyr::mutate(across(where(~ inherits(.x, "IDate")), as.Date))

## 2a. Severity composite indicator
# Oedema or haemorrhagic phenomena are the two cardinal complications
# used in Brazilian VL guidelines to define severe disease.
data_6 <- data_6 %>%
  mutate(severe = if_else(hemorrhagic_phenomena == 1 | edema == 1, "Yes", "No"))

## 2b. Demographic and clinical variables
data_6 <- data_6 %>%
  mutate(
    age_group_2 = case_when(
      age_years >= 0  & age_years <= 1  ~ "0-1",
      age_years > 1   & age_years <= 5  ~ "1-5",
      age_years > 5   & age_years <= 15 ~ "6-15",
      age_years > 15  & age_years <= 50 ~ "16-50",
      age_years > 50                    ~ ">50",
      TRUE ~ NA
    ),
    age_group_3 = case_when(
      age_years >= 0  & age_years <= 1  ~ "0-1",
      age_years >= 1  & age_years <= 2  ~ "1-2",
      age_years >= 2  & age_years <= 3  ~ "2-3",
      age_years >= 3  & age_years <= 5  ~ "3-5",
      age_years >= 5  & age_years <= 15 ~ "5-15",
      age_years >= 15 & age_years <= 20 ~ "15-20",
      age_years >= 20 & age_years <= 40 ~ "20-40",
      age_years >= 40 & age_years <= 60 ~ "40-60",
      age_years > 60                    ~ ">60",
      TRUE ~ NA
    ),
    year        = lubridate::year(clean_date_of_notification),
    year_period = case_when(
      str_detect(year, "2007|2008|2009")           ~ "2007-2009",
      str_detect(year, "2010|2011|2012|2013|2014") ~ "2010-2014",
      str_detect(year, "2015|2016|2017|2018")      ~ "2015-2018",
      str_detect(year, "2019|2020|2021|2022|2023") ~ "2019-2023",
      TRUE ~ NA
    ),
    race_2 = case_when(
      race == 1 ~ "White",         race == 2 ~ "Afro-descendants",
      race == 3 ~ "Yellow",        race == 4 ~ "Brown",
      race == 5 ~ "Indigenous",    race == 9 ~ "Unknown",
      TRUE ~ NA
    ),
    region = case_when(
      federative_unit %in% c(11, 13, 14, 15, 16, 17)             ~ "North",
      federative_unit %in% c(21, 22, 23, 24, 25, 26, 27, 28, 29) ~ "Northeast",
      federative_unit %in% c(31, 32, 33, 35)                     ~ "Southeast",
      federative_unit %in% c(41, 42, 43)                         ~ "South",
      federative_unit %in% c(50, 51, 52, 53)                     ~ "Center-West",
      TRUE                                                        ~ NA_character_
    ),
    age_group_3     = factor(age_group_3, levels = c(
      "0-1","1-2","2-3","3-5","5-15","15-20","20-40","40-60",">60"
    )),
    education_level = factor(education_level, levels = c(
      "Incomplete elementary school",
      "Complete elementary and incomplete middle school",
      "Complete middle school and incomplete high school",
      "Complete high school", "Incomplete university",
      "Complete university", "N/A"
    )),
    age_group_2 = factor(age_group_2, levels = c("0-1","1-5","6-15","16-50",">50")),
    fever_2 = factor(
      case_when(fever == 1 ~ "Yes", fever == 2 ~ "No", fever == 9 ~ "Unknown", TRUE ~ NA),
      levels = c("No", "Yes", "Unknown")
    ),
    edema_2 = factor(
      case_when(edema == 1 ~ "Yes", edema == 2 ~ "No", edema == 9 ~ "Unknown", TRUE ~ NA),
      levels = c("No", "Yes", "Unknown")
    ),
    hemorrhagic_phenomena_2 = factor(
      case_when(
        hemorrhagic_phenomena == 1 ~ "Yes",
        hemorrhagic_phenomena == 2 ~ "No",
        hemorrhagic_phenomena == 9 ~ "Unknown",
        TRUE ~ NA
      ),
      levels = c("No", "Yes", "Unknown")
    ),
    hiv_2 = case_when(
      hiv == 1 ~ "Yes", hiv == 2 ~ "No", hiv == 9 ~ "Unknown", TRUE ~ NA
    ),
    name_tx = case_when(
      type_of_rx == 1 ~ "Antimonial Pentavalent",
      type_of_rx == 2 ~ "Amphotericin B deoxycholate",
      type_of_rx == 3 ~ "Others",
      type_of_rx == 4 ~ "Liposomal amphotericin B",
      type_of_rx == 5 ~ "Others",
      type_of_rx == 6 ~ NA_character_,
      TRUE            ~ NA_character_
    )
  )

## 2c. Period-level summary for state-level mapping (used downstream)
data_ost <- data_6 %>%
  mutate(
    year_3 = case_when(
      year >= 2007 & year <= 2009 ~ "2007 - 2009",
      year >= 2010 & year <= 2014 ~ "2010 - 2014",
      year >= 2015 & year <= 2018 ~ "2015 - 2018",
      year >= 2019 & year <= 2023 ~ "2019 - 2023",
      TRUE ~ NA
    )
  ) %>%
  group_by(federative_unit, year_3) %>%
  summarise(
    median    = median(dif_day_tx_symp),
    n         = n(),
    mortality = sum(death_vl),
    lethality = round(mortality / n * 100, digits = 2),
    median_age = median(age_years),
    .groups = "drop"
  )

data_ost_2 <- data_ost %>% filter(n >= 100)

# States with fewer than 100 cases in a given period are excluded from
# state-level plots to avoid unstable median estimates.
state_100 <- unique(data_ost$federative_unit[data_ost$n < 100])


## 2d. Derived dataset with binary/numeric indicator variables
data_7 <- data_6 %>%
  mutate(
    edema_2 = factor(
      if_else(edema == 1, "Yes oedema", "No oedema"),
      levels = c("Yes oedema", "No oedema")
    ),
    severe = factor(
      if_else(edema == 1 | hemorrhagic_phenomena == 1, "Severe", "No severe"),
      levels = c("Severe", "No severe")
    ),
    # hiv_3: 0=No, 1=Yes, 2=Unknown — used as numeric indicator for some models
    # NOTE: replace_na(hiv_2 == "Unknown") is used here as in the original code;
    #       it returns FALSE/NA for missing values rather than a default.
    hiv_3 = replace_na(hiv_2 == "Unknown"),
    hiv_2 = factor(hiv_2, levels = c("Yes", "No"))
  )

data_8 <- data_7 %>%
  mutate(
    hiv_2    = if_else(hiv_2 == "Unknown", NA, hiv_2),
    gender_2 = factor(if_else(gender == "Male", 1, 0)),
    hiv_3    = factor(case_when(
      hiv_2 == "Yes"     ~ 1,
      hiv_2 == "No"      ~ 0,
      hiv_2 == "Unknown" ~ 2,
      TRUE ~ NA
    )),
    severe_2        = if_else(edema == 1, 1, 0),
    federative_unit = as.factor(federative_unit),
    # Numeric delay category (1–4) needed for ordinal or trend tests
    dif_treat_start_symp_cat_2 = case_when(
      dif_treat_start_symp_cat == "< 30 days"   ~ 1,
      dif_treat_start_symp_cat == "30-89 days"  ~ 2,
      dif_treat_start_symp_cat == "90-119 days" ~ 3,
      dif_treat_start_symp_cat == "≥120 days"   ~ 4,
      TRUE ~ NA
    )
  )

# Cache data_8 for downstream scripts that import it directly
export(data_8, here("data base and clean", "data_8.rds"))

data_8 <- data_8 %>%
  dplyr::mutate(across(where(~ inherits(.x, "IDate")), as.Date))


# ── 3. UNADJUSTED COX MODEL — OVERALL AND BY HIV STATUS ───────────────────────
# cph() (rms) is preferred for spline-based HR prediction via Predict().
# The unadjusted model serves as a descriptive check of the delay–mortality
# association before covariate adjustment.

## 3a. Overall (all patients, not filtered by HIV status)
data_m1 <- data_8 %>%
  filter(!is.na(death_vl)) %>%
  dplyr::select(
    dif_day_tx_symp, dif_treat_start_symp_cat_2, dif_treat_outcome,
    age_group_2, age_years, gender_2, hiv_3, hiv_2, gender,
    death_vl, name_tx, severe_2, race_2,
    federative_unit, education_level, region
  )

dd <- datadist(data_m1)
options(datadist = "dd")

# Unadjusted cph: 4-knot RCS on delay, no other covariates
m_overall <- cph(
  Surv(dif_treat_outcome, death_vl) ~ rcs(dif_day_tx_symp, 4),
  data = data_m1
)
summary(m_overall)

m_overall_pred <- Predict(m_overall, dif_day_tx_symp, fun = exp, ref.zero = FALSE)
m_overall_df   <- as.data.frame(m_overall_pred)

# NOTE: m1 is reused below as a ggplot object — here it stores the HR plot
m1 <- ggplot(m_overall_df, aes(dif_day_tx_symp, yhat)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_hline(aes(yintercept = 1), linetype = "dashed") +
  labs(
    x = "Delay to treatment initiation (days)",
    y = "Hazard ratio (Death by VL)"
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 14, color = "black", face = "bold"),
    axis.text  = element_text(size = 12, color = "black")
  ) +
  scale_x_continuous(n.breaks = 8) +
  scale_y_continuous(n.breaks = 10)
m1


## 3b. Unadjusted model excluding HIV-positive patients
# HIV co-infection is a strong competing predictor; removing HIV+ patients
# isolates the delay effect in the HIV-naive VL population.
data_m1_nohiv <- data_8 %>%
  filter(!is.na(death_vl), !hiv_2 == "Yes") %>%
  dplyr::select(
    dif_day_tx_symp, dif_treat_start_symp_cat_2, dif_treat_outcome,
    age_group_2, age_years, gender_2, hiv_3, hiv_2, gender,
    death_vl, name_tx, severe_2, race_2,
    federative_unit, education_level, region
  )

dd <- datadist(data_m1_nohiv)
options(datadist = "dd")

m1_fit_nohiv <- coxph(
  Surv(dif_treat_outcome, death_vl) ~ rcs(dif_day_tx_symp, 4),
  data = data_m1_nohiv
)
summary(m1_fit_nohiv)
termplot(m1_fit_nohiv, term = 1, se = TRUE, col.term = 1, col.se = 1)

# Prediction grid for plotting (linear predictor → HR)
newd <- data.frame(
  dif_day_tx_symp = seq(
    min(data_m1_nohiv$dif_day_tx_symp, na.rm = TRUE),
    max(data_m1_nohiv$dif_day_tx_symp, na.rm = TRUE),
    length.out = 200
  )
)
pred <- predict(m1_fit_nohiv, newdata = newd, type = "lp", se.fit = TRUE)
newd$hr    <- exp(pred$fit)
newd$lower <- exp(pred$fit - 1.96 * pred$se.fit)
newd$upper <- exp(pred$fit + 1.96 * pred$se.fit)

m1_whiv <- ggplot(newd, aes(dif_day_tx_symp, hr)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(
    x = "Delay to treatment initiation (days)",
    y = "Hazard ratio (Death by VL)"
  ) +
  theme_bw() +
  scale_x_continuous(n.breaks = 8) +
  scale_y_continuous(n.breaks = 10, limits = c(0.75, 2.8))
m1_whiv


## 3c. Unadjusted model for HIV-positive patients only
data_hiv_only <- data_8 %>%
  filter(!is.na(death_vl), !hiv_2 == "No")

m1_fit_hiv <- coxph(
  Surv(dif_treat_outcome, death_vl) ~ rcs(dif_day_tx_symp, 4),
  data = data_hiv_only
)
summary(m1_fit_hiv)

newd_hiv <- data.frame(
  dif_day_tx_symp = seq(
    min(data_hiv_only$dif_day_tx_symp, na.rm = TRUE),
    max(data_hiv_only$dif_day_tx_symp, na.rm = TRUE),
    length.out = 200
  )
)

pred_hiv       <- predict(m1_fit_hiv, newdata = newd_hiv, type = "lp", se.fit = TRUE)
newd_hiv$hr    <- exp(pred_hiv$fit)
newd_hiv$lower <- exp(pred_hiv$fit - 1.96 * pred_hiv$se.fit)
newd_hiv$upper <- exp(pred_hiv$fit + 1.96 * pred_hiv$se.fit)

m1_hiv_plot <- ggplot(newd_hiv, aes(dif_day_tx_symp, hr)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(
    x = "Delay to treatment initiation (days)",
    y = "Hazard ratio (Death by VL)"
  ) +
  theme_bw() +
  scale_x_continuous(n.breaks = 8) +
  scale_y_continuous(n.breaks = 10, limits = c(0.75, 2.8))

# Side-by-side panel: HIV-negative (left) and HIV-positive (right)
hiv_model_panel <- ggarrange(
  m1_whiv, m1_hiv_plot,
  labels    = c("a)", "b)"),
  hjust     = c(-0.7, -0.7),
  font.label = list(family = "Times New Roman")
)
hiv_model_panel

# Overlay both curves on the same axis to compare shapes
combine_model_1 <- ggplot(newd, aes(dif_day_tx_symp, hr, color = "No HIV")) +
  geom_line(linewidth = 0.5) +
  geom_ribbon(aes(ymin = lower, ymax = upper, color = "No HIV"),
              fill = "darkgreen", alpha = 0.2) +
  geom_line(
    data = newd_hiv,
    aes(dif_day_tx_symp, hr, color = "HIV"),
    linewidth = 0.5, linetype = "dashed"
  ) +
  geom_ribbon(
    data = newd_hiv,
    aes(ymin = lower, ymax = upper, color = "HIV"),
    fill = "brown2", alpha = 0.2
  ) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(
    x     = "Delay to treatment initiation (days)",
    y     = "Hazard ratio (Death by VL)",
    color = "Condition"
  ) +
  scale_x_continuous(n.breaks = 8) +
  scale_y_continuous(n.breaks = 10, limits = c(0.75, 2.8)) +
  scale_color_manual(values = c("HIV" = "brown2", "No HIV" = "darkgreen")) +
  theme_bw() +
  theme(
    axis.text    = element_text(size = 12, color = "black"),
    axis.title   = element_text(size = 14, color = "black", face = "bold"),
    legend.title = element_text(size = 14, color = "black", face = "bold"),
    legend.text  = element_text(size = 12, color = "black")
  )
combine_model_1


# ── 4. ADJUSTED COX MODEL (coxph) ─────────────────────────────────────────────
# This version uses coxph() + manual prediction grid (instead of cph + Predict)
# to verify results from the rms approach.  The two should yield identical HRs.

var_order_2 <- c("age_group_2", "gender", "region", "education_level", "race_2")

level_orders <- list(
  age_group_2     = c("0-1", "1-5", "6-15", "16-50", ">50"),
  gender          = c("Female", "Male"),
  region          = c("North", "Center-West", "Northeast", "Southeast", "South"),
  education_level = c(
    "Complete elementary and incomplete middle school",
    "Incomplete elementary school",
    "Complete middle school and incomplete high school",
    "Complete high school", "Incomplete university",
    "Complete university", "N/A"
  ),
  race_2 = c("Brown", "White", "Afro-descendants", "Yellow", "Indigenous", "Unknown")
)

data_mort <- data_m1_nohiv

for (v in var_order_2) {
  if (!is.null(level_orders[[v]])) {
    data_mort[[v]] <- factor(data_mort[[v]], levels = level_orders[[v]])
  } else {
    data_mort[[v]] <- droplevels(as.factor(data_mort[[v]]))
  }
}

# Adjusted Cox with RCS on delay — reference profile held at first levels
m2 <- coxph(
  Surv(dif_treat_outcome, death_vl) ~
    rcs(dif_day_tx_symp, k = 4) +
    age_group_2 + gender + race_2 + region + education_level,
  data = data_8
)
summary(m2)

# Prediction over the 1st–99th percentile of observed delays
delay_seq <- seq(
  quantile(data_mort$dif_day_tx_symp, 0.01, na.rm = TRUE),
  quantile(data_mort$dif_day_tx_symp, 0.99, na.rm = TRUE),
  length.out = 100
)

# Reference profile: females, North region, most common education, Brown race
newdat <- data.frame(
  dif_day_tx_symp = delay_seq,
  age_group_2     = ">50",
  gender          = "Female",
  region          = "North",
  education_level = "Complete elementary and incomplete middle school",
  race_2          = "Brown"
)
for (v in var_order_2) {
  newdat[[v]] <- factor(newdat[[v]], levels = levels(data_mort[[v]]))
}

lp     <- predict(m2, newdata = newdat, type = "lp", se.fit = TRUE)
ref_lp <- predict(m2, newdata = newdat[1, , drop = FALSE], type = "lp")

newdat <- newdat %>%
  mutate(
    hr    = exp(lp$fit - ref_lp),
    lower = exp(lp$fit - ref_lp - 1.96 * lp$se.fit),
    upper = exp(lp$fit - ref_lp + 1.96 * lp$se.fit)
  )

ggplot(newdat, aes(x = dif_day_tx_symp, y = hr)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(
    x     = "Delay between symptoms onset and treatment initiation (days)",
    y     = "Hazard ratio (Death by VL)",
    title = "Adjusted HR of VL death vs delay symptoms onset to treatment initiation"
  ) +
  scale_x_continuous(limits = c(0, 250), n.breaks = 15) +
  theme_bw() +
  theme(
    axis.title  = element_text(size = 14, color = "black", face = "bold"),
    axis.text   = element_text(size = 12, color = "black"),
    plot.title  = element_text(hjust = 0.5, color = "black", size = 16, face = "bold")
  )


# ── 5. SEVERITY SYMPTOM MODELS BY HIV STATUS ──────────────────────────────────
# Three severe complications are modelled individually (bleeding, jaundice,
# oedema) using lrm() with RCS on delay.  Models are run separately for
# HIV-negative and HIV-positive patients to assess whether HIV modifies
# the delay–complication relationship.
# fun = exp returns OR (odds ratio); ref.zero = TRUE sets reference at median delay.

## ── 5a. HIV-negative patients ────────────────────────────────────────────────
# NOTE: data_symp_no_hiv is defined HERE (HIV-negative) and reused in
# Model 1 below.  The naming reflects the subgroup, not the symptom status.
data_symp_no_hiv <- data_8 %>%
  mutate(
    jaundice_2 = case_when(jaundice == 1 ~ 1, jaundice == 2 ~ 0, TRUE ~ NA),
    bleeding   = case_when(
      hemorrhagic_phenomena == 1 ~ 1, hemorrhagic_phenomena == 2 ~ 0, TRUE ~ NA
    ),
    edema_3    = case_when(edema == 1 ~ 1, edema == 2 ~ 0, TRUE ~ NA)
  ) %>%
  filter(!hiv_2 == "Yes")   # retains HIV-negative and Unknown

### Bleeding model — no HIV
data_bleeding <- data_symp_no_hiv %>%
  dplyr::select(bleeding, dif_day_tx_symp) %>%
  filter(!is.na(bleeding))

dd <- datadist(data_bleeding)
options(datadist = "dd")

bleeding_fit    <- lrm(bleeding ~ rcs(dif_day_tx_symp, 4), data = data_bleeding)
p_bleeding      <- Predict(bleeding_fit, dif_day_tx_symp, fun = exp, ref.zero = TRUE)
p_bleeding_df   <- as.data.frame(p_bleeding)

bleeding_m1 <- ggplot(p_bleeding_df, aes(dif_day_tx_symp, yhat)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  labs(x = "Delay to treatment initiation (days)", y = "Odds (Bleeding)") +
  theme_bw()
bleeding_m1

### Jaundice model — no HIV
data_jaundice <- data_symp_no_hiv %>%
  dplyr::select(jaundice_2, dif_day_tx_symp) %>%
  dplyr::filter(!is.na(jaundice_2))

dd <- datadist(data_jaundice)
options(datadist = "dd")

jaundice_fit    <- lrm(jaundice_2 ~ rcs(dif_day_tx_symp, 4), data = data_jaundice)
p_jaundice      <- Predict(jaundice_fit, dif_day_tx_symp, fun = exp, ref.zero = TRUE)
p_jaundice_df   <- as.data.frame(p_jaundice)

jaundice_m1 <- ggplot(p_jaundice_df, aes(dif_day_tx_symp, yhat)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  labs(x = "Delay to treatment initiation (days)", y = "Odds (Jaundice)") +
  theme_bw()
jaundice_m1

### Oedema model — no HIV
data_oedema <- data_symp_no_hiv %>%
  dplyr::select(dif_day_tx_symp, edema_3) %>%
  filter(!is.na(edema_3))

dd <- datadist(data_oedema)
options(datadist = "dd")

edema_fit  <- lrm(edema_3 ~ rcs(dif_day_tx_symp, 4), data = data_oedema)
p_edema    <- Predict(edema_fit, dif_day_tx_symp, fun = exp, ref.zero = TRUE)
p_edema_df <- as.data.frame(p_edema)

edema_plot <- ggplot(p_edema_df, aes(dif_day_tx_symp, yhat)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  labs(x = "Delay to treatment initiation (days)", y = "Odds (Oedema)") +
  theme_bw()
edema_plot

### Combined panel — no HIV
symptoms_m1 <- ggplot(p_bleeding_df, aes(dif_day_tx_symp, yhat)) +
  geom_line(aes(color = "Bleeding"), linewidth = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "brown2") +
  geom_line(data = p_jaundice_df,
            aes(dif_day_tx_symp, yhat, color = "Jaundice"), linewidth = 1) +
  geom_ribbon(data = p_jaundice_df,
              aes(x = dif_day_tx_symp, ymin = lower, ymax = upper),
              fill = "orange2", alpha = 0.2, inherit.aes = FALSE) +
  geom_line(data = p_edema_df,
            aes(dif_day_tx_symp, yhat, color = "Oedema"), linewidth = 1) +
  geom_ribbon(data = p_edema_df,
              aes(x = dif_day_tx_symp, ymin = lower, ymax = upper),
              fill = "darkgreen", alpha = 0.2, inherit.aes = FALSE) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(x = "Delay to treatment initiation (days)", y = "Odds Ratio", color = "Symptoms") +
  scale_y_continuous(limits = c(0, 5)) +
  scale_color_manual(values = c(
    "Jaundice" = "orange2", "Bleeding" = "brown2", "Oedema" = "darkgreen"
  )) +
  theme_bw()
symptoms_m1


## ── 5b. HIV-positive patients ─────────────────────────────────────────────────
# data_symp_hiv_pos contains HIV-positive patients only.
# NOTE: the original script named this data_symp_no_hiv (confusing); renamed
# here for clarity.  The filter !hiv_2 == "No" retains hiv_2 == "Yes" and
# hiv_2 == "Unknown".
data_symp_hiv_pos <- data_8 %>%
  mutate(
    jaundice_2 = case_when(jaundice == 1 ~ 1, jaundice == 2 ~ 0, TRUE ~ NA),
    bleeding   = case_when(
      hemorrhagic_phenomena == 1 ~ 1, hemorrhagic_phenomena == 2 ~ 0, TRUE ~ NA
    ),
    edema_3    = case_when(edema == 1 ~ 1, edema == 2 ~ 0, TRUE ~ NA)
  ) %>%
  filter(!hiv_2 == "No")   # retains HIV-positive and Unknown

### Bleeding model — HIV+
data_bleeding_2 <- data_symp_hiv_pos %>%
  dplyr::select(bleeding, dif_day_tx_symp) %>%
  filter(!is.na(bleeding))

dd <- datadist(data_bleeding_2)
options(datadist = "dd")

bleeding_fit_2   <- lrm(bleeding ~ rcs(dif_day_tx_symp, 4), data = data_bleeding_2)
p_bleeding_2     <- Predict(bleeding_fit_2, dif_day_tx_symp, fun = exp, ref.zero = TRUE)
p_bleeding_df_2  <- as.data.frame(p_bleeding_2)

bleeding_m2 <- ggplot(p_bleeding_df_2, aes(dif_day_tx_symp, yhat)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  labs(x = "Delay to treatment initiation (days)", y = "Odds (Bleeding)") +
  theme_bw()
bleeding_m2

### Jaundice model — HIV+
data_jaundice_2 <- data_symp_hiv_pos %>%
  dplyr::select(jaundice_2, dif_day_tx_symp) %>%
  dplyr::filter(!is.na(jaundice_2))

dd <- datadist(data_jaundice_2)
options(datadist = "dd")

jaundice_fit_2   <- lrm(jaundice_2 ~ rcs(dif_day_tx_symp, 4), data = data_jaundice_2)
p_jaundice_2     <- Predict(jaundice_fit_2, dif_day_tx_symp, fun = exp, ref.zero = TRUE)
p_jaundice_df_2  <- as.data.frame(p_jaundice_2)

jaundice_m2 <- ggplot(p_jaundice_df_2, aes(dif_day_tx_symp, yhat)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  labs(x = "Delay to treatment initiation (days)", y = "Odds (Jaundice)") +
  theme_bw()
jaundice_m2

### Oedema model — HIV+
data_oedema_2 <- data_symp_hiv_pos %>%
  dplyr::select(dif_day_tx_symp, edema_3) %>%
  filter(!is.na(edema_3))

dd <- datadist(data_oedema_2)
options(datadist = "dd")

edema_fit_2   <- lrm(edema_3 ~ rcs(dif_day_tx_symp, 4), data = data_oedema_2)
p_edema_2     <- Predict(edema_fit_2, dif_day_tx_symp, fun = exp, ref.zero = TRUE)
p_edema_df_2  <- as.data.frame(p_edema_2)

edema_2 <- ggplot(p_edema_df_2, aes(dif_day_tx_symp, yhat)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  labs(x = "Delay to treatment initiation (days)", y = "Odds (Oedema)") +
  theme_bw()
edema_2

### Combined panel — HIV+
symptoms_m2 <- ggplot(p_bleeding_df_2, aes(dif_day_tx_symp, yhat)) +
  geom_line(aes(color = "Bleeding"), linewidth = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "brown2") +
  geom_line(data = p_jaundice_df_2,
            aes(dif_day_tx_symp, yhat, color = "Jaundice"), linewidth = 1) +
  geom_ribbon(data = p_jaundice_df_2,
              aes(x = dif_day_tx_symp, ymin = lower, ymax = upper),
              fill = "orange2", alpha = 0.2, inherit.aes = FALSE) +
  geom_line(data = p_edema_df_2,
            aes(dif_day_tx_symp, yhat, color = "Oedema"), linewidth = 1) +
  geom_ribbon(data = p_edema_df_2,
              aes(x = dif_day_tx_symp, ymin = lower, ymax = upper),
              fill = "darkgreen", alpha = 0.2, inherit.aes = FALSE) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(x = "Delay to treatment initiation (days)", y = "Odds Ratio", color = "Symptoms") +
  scale_y_continuous(limits = c(0, 5)) +
  scale_color_manual(values = c(
    "Jaundice" = "orange2", "Bleeding" = "brown2", "Oedema" = "darkgreen"
  )) +
  theme_bw()
symptoms_m2


## ── 5c. Overall (all HIV statuses) ───────────────────────────────────────────
data_symp_overall_hiv <- data_8 %>%
  mutate(
    jaundice_2 = case_when(jaundice == 1 ~ 1, jaundice == 2 ~ 0, TRUE ~ NA),
    bleeding   = case_when(
      hemorrhagic_phenomena == 1 ~ 1, hemorrhagic_phenomena == 2 ~ 0, TRUE ~ NA
    ),
    edema_3    = case_when(edema == 1 ~ 1, edema == 2 ~ 0, TRUE ~ NA)
  )

### Bleeding — overall
data_bleeding_3 <- data_symp_overall_hiv %>%
  dplyr::select(bleeding, dif_day_tx_symp) %>%
  filter(!is.na(bleeding))

dd <- datadist(data_bleeding_3)
options(datadist = "dd")

bleeding_fit_3  <- lrm(bleeding ~ rcs(dif_day_tx_symp, 4), data = data_bleeding_3)
p_bleeding_3    <- Predict(bleeding_fit_3, dif_day_tx_symp, fun = exp, ref.zero = TRUE)
p_bleeding_df_3 <- as.data.frame(p_bleeding_3)

bleeding_m3 <- ggplot(p_bleeding_df_3, aes(dif_day_tx_symp, yhat)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  labs(x = "Delay to treatment initiation (days)", y = "Odds (Bleeding)") +
  theme_bw()
bleeding_m3

### Jaundice — overall
data_jaundice_3 <- data_symp_overall_hiv %>%
  dplyr::select(jaundice_2, dif_day_tx_symp) %>%
  filter(!is.na(jaundice_2))

dd <- datadist(data_jaundice_3)
options(datadist = "dd")

jaundice_fit_3  <- lrm(jaundice_2 ~ rcs(dif_day_tx_symp, 4), data = data_jaundice_3)
p_jaundice_3    <- Predict(jaundice_fit_3, dif_day_tx_symp, fun = exp, ref.zero = TRUE)
p_jaundice_df_3 <- as.data.frame(p_jaundice_3)

jaundice_m3 <- ggplot(p_jaundice_df_3, aes(dif_day_tx_symp, yhat)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  labs(x = "Delay to treatment initiation (days)", y = "Odds (Jaundice)") +
  theme_bw()
jaundice_m3

### Oedema — overall
data_oedema_3 <- data_symp_overall_hiv %>%
  dplyr::select(dif_day_tx_symp, edema_3) %>%
  filter(!is.na(edema_3))

dd <- datadist(data_oedema_3)
options(datadist = "dd")

edema_fit_3  <- lrm(edema_3 ~ rcs(dif_day_tx_symp, 4), data = data_oedema_3)
p_edema_3    <- Predict(edema_fit_3, dif_day_tx_symp, fun = exp, ref.zero = TRUE)
p_edema_df_3 <- as.data.frame(p_edema_3)

edema_3 <- ggplot(p_edema_df_3, aes(dif_day_tx_symp, yhat)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  labs(x = "Delay to treatment initiation (days)", y = "Odds (Oedema)") +
  theme_bw()
edema_3

### Combined panel — overall
symptoms_m3 <- ggplot(p_bleeding_df_3, aes(dif_day_tx_symp, yhat)) +
  geom_line(aes(color = "Bleeding"), linewidth = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "brown2") +
  geom_line(data = p_jaundice_df_3,
            aes(dif_day_tx_symp, yhat, color = "Jaundice"), linewidth = 1) +
  geom_ribbon(data = p_jaundice_df_3,
              aes(x = dif_day_tx_symp, ymin = lower, ymax = upper),
              fill = "orange2", alpha = 0.2, inherit.aes = FALSE) +
  geom_line(data = p_edema_df_3,
            aes(dif_day_tx_symp, yhat, color = "Oedema"), linewidth = 1) +
  geom_ribbon(data = p_edema_df_3,
              aes(x = dif_day_tx_symp, ymin = lower, ymax = upper),
              fill = "darkgreen", alpha = 0.2, inherit.aes = FALSE) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(x = "Delay to treatment initiation (days)", y = "Odds Ratio", color = "Symptoms") +
  theme_bw() +
  scale_y_continuous(limits = c(0, 5)) +
  scale_color_manual(values = c(
    "Jaundice" = "orange2", "Bleeding" = "brown2", "Oedema" = "darkgreen"
  ))
symptoms_m3

## ── Final three-panel figure (No HIV | HIV+ | Overall) ───────────────────────
final_models <- ggarrange(
  symptoms_m1 + rremove("ylab") + rremove("xlab") + scale_y_continuous(limits = c(0.5, 5)),
  symptoms_m2 + rremove("xlab") + rremove("ylab") + scale_y_continuous(limits = c(0.5, 5)),
  symptoms_m3 + rremove("xlab") + rremove("ylab") + scale_y_continuous(limits = c(0.5, 5)),
  labels        = "AUTO",
  nrow          = 1,
  legend        = "right",
  hjust         = -0.5,
  common.legend = TRUE,
  font.label    = list(size = 14, color = "black", face = "bold",
                       family = "Times New Roman", position = "bottom")
)

symptoms <- annotate_figure(
  final_models,
  bottom = textGrob(
    "Delay between symptoms onset and treatment",
    gp = gpar(fontsize = 14, fontface = "bold", fontfamily = "Times New Roman")
  ),
  left = textGrob(
    "Odds Ratio",
    gp = gpar(fontsize = 14, fontface = "bold", fontfamily = "Times New Roman"),
    rot = 90, vjust = 1
  )
)

ggsave(
  symptoms,
  path = here("figures and tables"),
  width = 10, height = 8, dpi = 300,
  file  = "severe_symptoms_color_graph.png"
)


# ── 6. COMPOSITE COMPLICATION MODELS (FIGURE SX) ──────────────────────────────
# Three composite endpoints capturing increasing severity:
#   complication_2: any of oedema, jaundice, or haemorrhage
#   complication_3: all three simultaneously (most severe)
#   complication_4: oedema AND haemorrhage (prognostic VL phenotype)

data_complication <- data_8 %>%
  mutate(
    complication_2 = if_else(edema == 1 | jaundice == 1 | hemorrhagic_phenomena == 1, 1, 0),
    complication_3 = if_else(edema == 1 & jaundice == 1 & hemorrhagic_phenomena == 1, 1, 0),
    complication_4 = if_else(edema == 1 & hemorrhagic_phenomena == 1, 1, 0)
  )

## ── Model function for composite outcomes ────────────────────────────────────
fit_severity_model <- function(data, outcome, y_label = "Odds Ratio (complications)", knots = 4) {

  df <- data %>%
    dplyr::select(all_of(c(outcome, "dif_day_tx_symp"))) %>%
    filter(!is.na(.data[[outcome]]), !is.na(dif_day_tx_symp))

  dd <- datadist(df)
  options(datadist = "dd")

  fit <- lrm(
    as.formula(paste0(outcome, " ~ rcs(dif_day_tx_symp,", knots, ")")),
    data = df
  )

  pred    <- Predict(fit, dif_day_tx_symp, fun = exp, ref.zero = TRUE)
  pred_df <- as.data.frame(pred)

  ggplot(pred_df, aes(dif_day_tx_symp, yhat)) +
    geom_line(linewidth = 1) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
    labs(x = "Time to treatment initiation (days)", y = y_label) +
    theme_bw() +
    theme(
      axis.title  = element_text(size = 14, color = "black", face = "bold"),
      axis.text.x = element_text(size = 12, color = "black")
    )
}

severity_specs <- tibble::tribble(
  ~outcome,         ~y_label,
  "complication_2", "Any complication",
  "complication_3", "All three complications",
  "complication_4", "Oedema + haemorrhage"
)

## ── Overall (all HIV statuses)
severity_plots <- purrr::pmap(severity_specs, fit_severity_model, data = data_complication)
names(severity_plots) <- severity_specs$y_label
final_severity <- ggpubr::ggarrange(plotlist = severity_plots, nrow = 1, common.legend = FALSE)
final_severity

## ── HIV-positive patients
severity_plots_hiv <- purrr::pmap(
  severity_specs, fit_severity_model,
  data = data_complication %>% filter(hiv_2 == "Yes")
)
names(severity_plots_hiv) <- severity_specs$y_label
final_severity_hiv <- ggpubr::ggarrange(
  plotlist = severity_plots_hiv, nrow = 1, common.legend = FALSE
)
final_severity_hiv

## ── HIV-negative patients
severity_plots_nohiv <- purrr::pmap(
  severity_specs, fit_severity_model,
  data = data_complication %>% filter(hiv_2 == "No")
)
names(severity_plots_nohiv) <- severity_specs$y_label
final_severity_nohiv <- ggpubr::ggarrange(
  plotlist = severity_plots_nohiv, nrow = 1, common.legend = FALSE
)
final_severity_nohiv

## ── Stacked three-row panel (overall / HIV+ / HIV-)
panel_model_complications <- ggarrange(
  final_severity,
  final_severity_hiv,
  final_severity_nohiv,
  nrow = 3, labels = "AUTO"
)
panel_model_complications


# ── 7. SYMPTOM PROBABILITY PANEL BY HIV STATUS ────────────────────────────────
# Runs lrm with RCS on delay for 11 symptoms, stratified by HIV status.
# The final publication version (with age + gender adjustment) lives in
# request_prabin.R; this version uses delay only as the predictor.

data_symptoms <- data_8 %>%
  dplyr::select(
    dif_day_tx_symp, fever, pallor, weakness, edema, weight_loss,
    cough_and_dyspnea, increase_infection, spleen_enlargement,
    hemorrhagic_phenomena, liver_enlargement, jaundice, hiv_2
  )

fit_symptoms_model_simple <- function(data, outcome, y_label = "Complication", knots = 4) {

  df <- data %>%
    filter(
      !is.na(.data[[outcome]]),
      !is.na(dif_day_tx_symp),
      hiv_2 == "No"
    ) %>%
    dplyr::select(all_of(c(outcome, "dif_day_tx_symp")))

  dd <- datadist(df)
  options(datadist = "dd")

  fit <- lrm(
    as.formula(paste0(outcome, " ~ rcs(dif_day_tx_symp,", knots, ")")),
    data = df
  )

  df_2 <- data %>%
    filter(
      !is.na(.data[[outcome]]),
      !is.na(dif_day_tx_symp),
      hiv_2 == "Yes"
    ) %>%
    dplyr::select(all_of(c(outcome, "dif_day_tx_symp")))

  dd <<- datadist(df_2)
  options(datadist = "dd")

  fit_2 <- lrm(
    as.formula(paste0(outcome, " ~ rcs(dif_day_tx_symp,", knots, ")")),
    data = df_2
  )

  pred   <- Predict(fit,   dif_day_tx_symp, fun = plogis, ref.zero = TRUE)
  pred_2 <- Predict(fit_2, dif_day_tx_symp, fun = plogis, ref.zero = TRUE)

  pred_df   <- as.data.frame(pred)
  pred_df_2 <- as.data.frame(pred_2)

  ggplot(pred_df, aes(dif_day_tx_symp, yhat, color = "HIV(-)")) +
    geom_line(linewidth = 1) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = "CI"), alpha = 0.2, color = NA) +
    geom_line(data = pred_df_2, aes(dif_day_tx_symp, yhat, color = "HIV(+)"), linewidth = 1) +
    geom_ribbon(data = pred_df_2, aes(ymin = lower, ymax = upper, fill = "CI"),
                alpha = 0.2, color = NA) +
    labs(
      x     = "Delay to treatment initiation (days)",
      y     = bquote(bolditalic("P(x)") ~ bold("(" * .(y_label) * ")")),
      color = "Group", fill = NULL
    ) +
    theme_bw() +
    theme(
      axis.title  = element_text(size = 12, color = "black", face = "bold"),
      axis.text.x = element_text(size = 11, color = "black")
    ) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_color_manual(values = c("HIV(-)" = "darkgreen", "HIV(+)" = "brown2")) +
    scale_fill_manual(values  = c("CI" = "grey80"))
}

symptoms_specs <- tibble::tribble(
  ~outcome,                ~y_label,
  "fever",                 "Fever",
  "pallor",                "Pallor",
  "weakness",              "Weakness",
  "edema",                 "Oedema",
  "weight_loss",           "Weight loss",
  "cough_and_dyspnea",     "Cough and dyspnea",
  "spleen_enlargement",    "Splenomegaly",
  "increase_infection",    "Secondary infection",
  "hemorrhagic_phenomena", "Hemorrhage",
  "liver_enlargement",     "Hepatomegaly",
  "jaundice",              "Jaundice"
)

symptoms_plots <- purrr::pmap(symptoms_specs, fit_symptoms_model_simple, data = data_symptoms)
names(symptoms_plots) <- symptoms_specs$y_label

# Bottom row (panels 9–11) retains x-axis; all others drop it to save space
bottom_idx <- 9:11
symptoms_plots2 <- purrr::map2(
  symptoms_plots, seq_along(symptoms_plots),
  ~ if (.y %in% bottom_idx) .x else .x +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
)

library(cowplot)
legend <- cowplot::get_legend(
  symptoms_plots2[[1]] +
    theme(legend.position = "right",
          legend.text  = element_text(size = 12, color = "black"),
          legend.title = element_text(size = 16, color = "black", face = "bold"))
)

plots_no_legend <- purrr::map(symptoms_plots2, ~ .x + theme(legend.position = "none"))

main_panel <- ggpubr::ggarrange(plotlist = plots_no_legend, nrow = 4, ncol = 3)

final_severity <- ggdraw(main_panel) +
  draw_grob(legend, x = 0.78, y = 0.02, width = 0.20, height = 0.20)

ggsave(
  filename = "distribution_symptoms_final.png",
  plot     = final_severity,
  path     = here("figures and tables"),
  width = 10, height = 8, dpi = 300
)


# ── 8. GLM LOGISTIC MODEL (EXPLORATORY) ───────────────────────────────────────
# Logistic regression as an alternative to Cox to assess delay → mortality.
# Not used in the final manuscript (Cox model is the primary analysis) but
# kept here for sensitivity checking.

m_glm <- glm(
  death_vl ~ rcs(dif_day_tx_symp, 4),
  family = "binomial",
  data   = data_6 %>% filter(!hiv_3 == 1)
)
summary(m_glm)

newdat_glm <- data.frame(
  dif_day_tx_symp = seq(
    min(data_8$dif_day_tx_symp, na.rm = TRUE),
    max(data_8$dif_day_tx_symp, na.rm = TRUE),
    length.out = 100
  )
)
newdat_glm$pred <- predict(m_glm, newdata = newdat_glm, type = "response")

plot(newdat_glm$dif_day_tx_symp, newdat_glm$pred * 100, type = "l", lwd = 2,
     xlab = "Days between symptoms and treatment",
     ylab = "Predicted probability of death (%)")


# ── 9. GLM ADJUSTED LOGISTIC MODEL ────────────────────────────────────────────
# Multivariable logistic model including all covariates.
# Geometric mean ratios (GMR) from this model feed into Table 2.

var_order_2 <- c(
  "age_group_2", "gender", "region",
  "education_level", "race_2", "hiv_2", "name_tx"
)

level_orders_glm <- list(
  age_group_2     = c("0-1", "1-5", "6-15", "16-50", ">50"),
  gender          = c("Female", "Male"),
  region          = c("Center-West", "North", "Northeast", "Southeast", "South"),
  education_level = c(
    "Incomplete elementary school",
    "Complete elementary and incomplete middle school",
    "Complete middle school and incomplete high school",
    "Complete high school", "Incomplete university",
    "Complete university", "N/A"
  ),
  race_2  = c("Brown", "White", "Afro-descendants", "Yellow", "Indigenous", "Unknown"),
  hiv_2   = c("No", "Yes", "Unknown"),
  name_tx = c(
    "Antimonial Pentavalent", "Amphotericin B deoxycholate",
    "Liposomal amphotericin B", "Others"
  )
)

data_mort_glm <- data_6 %>%
  filter(!is.na(death_vl), !is.na(dif_day_tx_symp), dif_day_tx_symp > 0)

for (v in var_order_2) {
  if (!is.null(level_orders_glm[[v]])) {
    data_mort_glm[[v]] <- factor(data_mort_glm[[v]], levels = level_orders_glm[[v]])
  } else {
    data_mort_glm[[v]] <- droplevels(as.factor(data_mort_glm[[v]]))
  }
}

form_multi <- as.formula(
  paste("death_vl ~ rcs(dif_day_tx_symp, 4) +",
        paste(var_order_2, collapse = " + "))
)

m_glm_adj <- glm(form_multi, data = data_mort_glm, family = binomial())
summary(m_glm_adj)

# Reference profile for prediction (represents a typical high-risk VL patient)
delay_seq_glm <- seq(
  quantile(data_mort_glm$dif_day_tx_symp, 0.01, na.rm = TRUE),
  quantile(data_mort_glm$dif_day_tx_symp, 0.99, na.rm = TRUE),
  length.out = 100
)

newdat_adj <- data.frame(
  dif_day_tx_symp = delay_seq_glm,
  age_group_2     = "16-50",
  gender          = "Male",
  region          = "South",
  education_level = "Incomplete elementary school",
  race_2          = "Afro-descendants",
  hiv_2           = "Yes",
  name_tx         = "Antimonial Pentavalent"
)

for (v in var_order_2) {
  newdat_adj[[v]] <- factor(newdat_adj[[v]], levels = levels(data_mort_glm[[v]]))
}

pred_adj  <- predict(m_glm_adj, newdata = newdat_adj, type = "link", se.fit = TRUE)

newdat_adj <- newdat_adj %>%
  mutate(
    fit    = pred_adj$fit,
    se.fit = pred_adj$se.fit,
    prob   = plogis(fit),
    lower  = plogis(fit - 1.96 * se.fit),
    upper  = plogis(fit + 1.96 * se.fit)
  )

ggplot(newdat_adj, aes(x = dif_day_tx_symp, y = prob)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line(linewidth = 1) +
  geom_hline(aes(yintercept = 0.05), linetype = "dotted", linewidth = 0.8) +
  geom_hline(aes(yintercept = 0.07), linetype = "dotted", linewidth = 0.8) +
  geom_rug(data = data_mort_glm,
           aes(x = dif_day_tx_symp, y = 0),
           inherit.aes = FALSE, sides = "b", alpha = 0.3) +
  scale_y_continuous(
    labels  = scales::percent_format(accuracy = 1),
    name    = "Predicted probability of death (%)",
    n.breaks = 15
  ) +
  labs(
    x        = "Days between symptoms and treatment",
    title    = "Predicted death probability vs treatment delay",
    subtitle = "Adjusted for age group, gender, region, education, race, HIV and treatment"
  ) +
  theme_classic() +
  scale_x_continuous(limits = c(0, 250), n.breaks = 15)
