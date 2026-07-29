# ══════════════════════════════════════════════════════════════════════════════
# Project : Visceral Leishmaniasis in Brazil, 2007–2023
# Script  : Prabin requirements
# Purpose : Manuscript adjustments following internal review by supervisors
#           (Prabin, Philippe, Makoto).  Includes:
#             - refined delay categorisation requested by reviewers
#             - descriptive statistics by region
#             - Kaplan–Meier curves by delay category (including sensitivity
#               analyses at 6 months)
#             - symptom probability models by HIV status (Figures 2 & 3)
#             - adjusted Cox proportional-hazards models stratified by HIV
#               status (Figure 5)
#             - baseline characteristics table stratified by delay (Table 1)
# Input   : data base and clean/data_6.rds  (SINAN VL notifications, curated)
# Output  : figures and tables/distribution_symptoms_final.png
#           figures and tables/distribution_symptoms_HIV.png
#           figures and tables/cox_model_hiv.png
#           figures and tables/cox_model_nohiv.png
#           figures and tables/cox_model.png
#           figures and tables/table_1.docx
# Author  : Juan Sebastián Hurtado Zapata
#           IDDO VL Fellowship | TDR/WHO Fellow
# Updated : 2026-07-29
# ══════════════════════════════════════════════════════════════════════════════


# ── 1. PACKAGES ───────────────────────────────────────────────────────────────
# Run install.packages("pacman") once if pacman is not yet available.
# Uncomment the lines below only for a fresh R environment:
# install.packages(c("Hmisc", "rms", "cowplot"))

pacman::p_load(
  flextable, tidyverse, janitor, epikit, gtsummary, rio,
  survival, survminer, sf, here, incidence, ggstance, ggbeeswarm,
  grid, ggsurvfit, geobr, rms, patchwork, biostat3, simukde,
  Hmisc, cowplot, ggpubr, broom, officer, scales
)


# ── 2. DATA IMPORT AND VARIABLE ENGINEERING ───────────────────────────────────
# data_prabin_2 extends the base dataset with the finer delay categorisation
# (delay_2, delay_categories) and region label requested by Prabin, as well
# as a three-period year grouping (year_3) used in reviewer queries.
# data_6 is re-imported separately because it undergoes additional variable
# engineering (symptom labels, treatment names) used in Cox/GLM sections below.

data_prabin <- import(here("data base and clean", "data_6.rds"))

data_prabin_2 <- data_prabin %>%
  mutate(
    # Coarser 5-level delay grouping (Table 1 and survival curves)
    delay_2 = case_when(
      dif_day_tx_symp <= 14                                    ~ "0-14 days",
      dif_day_tx_symp > 14 & dif_day_tx_symp <= 30            ~ "15-30 days",
      TRUE                                                     ~ dif_treat_start_symp_cat
    ),
    delay_2 = factor(delay_2, levels = c(
      "0-14 days", "15-30 days", "30-89 days", "90-119 days", "≥120 days"
    )),

    # Finer 8-level delay grouping (sensitivity survival analysis)
    delay_categories = case_when(
      dif_day_tx_symp >= 0   & dif_day_tx_symp <= 14  ~ "≤14 days",
      dif_day_tx_symp >= 15  & dif_day_tx_symp <= 30  ~ "15 to ≤ 30 days",
      dif_day_tx_symp >= 31  & dif_day_tx_symp <= 60  ~ "31 to ≤ 60 days",
      dif_day_tx_symp >= 61  & dif_day_tx_symp <= 90  ~ "61 to ≤ 90 days",
      dif_day_tx_symp >= 91  & dif_day_tx_symp <= 120 ~ "91 to ≤ 120 days",
      dif_day_tx_symp >= 121 & dif_day_tx_symp <= 150 ~ "121 to ≤ 150 days",
      dif_day_tx_symp >= 151 & dif_day_tx_symp <= 180 ~ "151 to ≤ 180 days",
      dif_day_tx_symp >= 181                           ~ "≥6 months",
      TRUE ~ NA
    ),
    delay_categories = factor(delay_categories, levels = c(
      "≤14 days", "15 to ≤ 30 days", "31 to ≤ 60 days",
      "61 to ≤ 90 days", "91 to ≤ 120 days", "121 to ≤ 150 days",
      "151 to ≤ 180 days", "≥6 months"
    )),

    # Four surveillance periods for trend analysis
    year_3 = case_when(
      year_of_notification >= 2007 & year_of_notification <= 2009 ~ "2007 - 2009",
      year_of_notification >= 2010 & year_of_notification <= 2014 ~ "2010 - 2014",
      year_of_notification >= 2015 & year_of_notification <= 2018 ~ "2015 - 2018",
      year_of_notification >= 2019 & year_of_notification <= 2023 ~ "2019 - 2023",
      TRUE ~ NA
    ),

    # Brazilian macro-regions using IBGE federative unit codes
    region = case_when(
      federative_unit %in% c(11, 13, 14, 15, 16, 17)             ~ "North",
      federative_unit %in% c(21, 22, 23, 24, 25, 26, 27, 28, 29) ~ "Northeast",
      federative_unit %in% c(31, 32, 33, 35)                     ~ "Southeast",
      federative_unit %in% c(41, 42, 43)                         ~ "South",
      federative_unit %in% c(50, 51, 52, 53)                     ~ "Center-West",
      TRUE                                                        ~ NA_character_
    )
  ) %>%
  # IDate columns from data.table must be converted; rms and survival expect Date
  dplyr::mutate(across(where(~ inherits(.x, "IDate")), as.Date))


## ── data_6: base dataset with full variable engineering ──────────────────────
# This re-import applies additional label/factor engineering used in Cox and
# symptom models.  It avoids overwriting data_prabin_2's delay variables.
data_6 <- import(here("data base and clean", "data_6.rds")) %>%
  dplyr::mutate(across(where(~ inherits(.x, "IDate")), as.Date)) %>%
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
      str_detect(year, "2007|2008|2009")             ~ "2007-2009",
      str_detect(year, "2010|2011|2012|2013|2014")   ~ "2010-2014",
      str_detect(year, "2015|2016|2017|2018")        ~ "2015-2018",
      str_detect(year, "2019|2020|2021|2022|2023")   ~ "2019-2023",
      TRUE ~ NA
    ),
    race_2 = case_when(
      race == 1 ~ "White",
      race == 2 ~ "Afro-descendants",
      race == 3 ~ "Yellow",
      race == 4 ~ "Brown",
      race == 5 ~ "Indigenous",
      race == 9 ~ "Unknown",
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
    # Ordered factors keep tables and models in the intended display order
    age_group_3 = factor(age_group_3, levels = c(
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
    # Binary symptom labels (SINAN coding: 1=Yes, 2=No, 9=Unknown)
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
    # Treatment label (type_of_rx == 6 = no treatment recorded → NA)
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


# ── 3. DESCRIPTIVE STATISTICS BY REGION (REVIEWER REQUEST — PRABIN) ───────────
# Median delay (ts→tt) and IQR by Brazilian macro-region

data_prabin_2 %>%
  group_by(region) %>%
  summarise(
    n        = n(),
    median_r = median(dif_day_tx_symp),
    IQL      = quantile(dif_day_tx_symp, 0.25),
    IQS      = quantile(dif_day_tx_symp, 0.75)
  ) %>%
  adorn_totals()


# ── 4. TREATMENT OUTCOME CROSS-TABULATION ─────────────────────────────────────
# Distribution of outcomes (cure vs. death) by time from treatment start to
# outcome; early deaths (≤7 days) were flagged by co-author Makoto as potentially
# reflecting pre-treatment severity rather than treatment failure.

data_6 %>%
  filter(!is.na(name_tx), !is.na(death_vl)) %>%
  mutate(
    days_p = case_when(
      dif_treat_outcome <= 7                              ~ "≤7 days",
      dif_treat_outcome >= 8 & dif_treat_outcome <= 14   ~ "7-14 days",
      TRUE                                               ~ "Others"
    ),
    days_p = factor(days_p, levels = c("≤7 days", "7-14 days", "Others"))
  ) %>%
  group_by(days_p, evolution) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(id_cols = "days_p", names_from = "evolution", values_from = n) %>%
  adorn_totals("both") %>%
  adorn_percentages(denominator = "all") %>%
  adorn_pct_formatting() %>%
  adorn_ns(position = "front")

data_6 %>% tabyl(evolution)


# ── 5. SURVIVAL ANALYSIS (KM CURVES BY DELAY CATEGORY) ───────────────────────

## 5a. KM curve — 5-level delay grouping (delay_2)
# This categorisation splits the original "< 30 days" category into 0–14 and
# 15–30 days, which was requested because early-treatment patients (<14 days)
# may have a different risk profile.

delay <- survfit2(Surv(dif_treat_outcome, death_vl) ~ delay_2, data = data_prabin_2)

summary(delay, times = c(180), data.frame = TRUE)

gd <- ggsurvfit(delay, linewidth = 1) +
  add_risktable(
    risktable_stats = "{n.risk} ({cum.event})",
    theme           = theme_bw(),
    risktable_height = 0.2
  ) +
  add_censor_mark() +
  add_pvalue(
    caption  = "Log-rank {p.value}",
    location = "annotation", x = 200, y = 0.98
  ) +
  scale_ggsurvfit() +
  coord_cartesian(ylim = c(0.85, 1)) +
  scale_x_continuous(n.breaks = 15) +
  scale_y_continuous(n.breaks = 12, label = scales::label_percent()) +
  # Reference line at 30 days marks the threshold for the original delay grouping
  add_quantile(x_value = 30, linetype = "dotted", color = "grey30", linewidth = 0.8) +
  labs(x = "Time (days after treatment start)") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.text       = element_text(color = "black", size = 10),
    axis.title      = element_text(color = "black", size = 14, face = "bold"),
    legend.title    = element_text(face = "bold", size = 14, color = "black"),
    legend.text     = element_text(color = "black", size = 12)
  )

print(gd)


## 5b. Sensitivity: people evaluated within the first six months (Makoto query)
# Restricting follow-up to 180 days ensures comparable censoring across groups.

data_prabin_2 %>%
  filter(!is.na(death_vl), dif_treat_outcome < 181) %>%
  count()


## 5c. KM curve — comparing extreme delay categories (≤14 days vs. ≥6 months)
# Peto log-rank test (rho = 1) gives more weight to early events, which is
# appropriate given the peak in VL mortality occurs in the first few weeks
# after treatment initiation.

data_prabin_3 <- data_6 %>%
  filter(delay_categories %in% c("≤14 days", "≥6 months"))

survdiff(Surv(dif_treat_outcome, death_vl) ~ delay_categories, rho = 1,
         data = data_prabin_3)

delay_extreme <- survfit2(
  Surv(dif_treat_outcome, death_vl) ~ delay_categories,
  data = data_prabin_3
)

summary(delay_extreme, times = c(360), data.frame = TRUE)

ggsurvfit(delay_extreme, linewidth = 1) +
  add_risktable(
    risktable_stats  = "{n.risk} ({cum.event})",
    theme            = theme_bw(),
    risktable_height = 0.2
  ) +
  add_censor_mark() +
  add_pvalue(
    caption  = "Log-rank {p.value}",
    location = "annotation", x = 200, y = 0.98
  ) +
  scale_ggsurvfit() +
  coord_cartesian(ylim = c(0.85, 1)) +
  scale_x_continuous(n.breaks = 15) +
  scale_y_continuous(n.breaks = 12, label = scales::label_percent()) +
  add_quantile(x_value = 30, linetype = "dotted", color = "grey30", linewidth = 0.8) +
  labs(x = "Time (days after treatment start)") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.text       = element_text(color = "black", size = 10),
    axis.title      = element_text(color = "black", size = 14, face = "bold"),
    legend.title    = element_text(face = "bold", size = 14, color = "black"),
    legend.text     = element_text(color = "black", size = 12)
  )


## 5d. Overall fatality risk (crude KM)
fatality_risk <- survfit(Surv(dif_treat_outcome, death_vl) ~ 1, data = data_6)
summary(fatality_risk, times = c(0, 30, 60, 90, 180, 365))

gm <- ggsurvplot(
  fit          = fatality_risk,
  data         = data_6,
  ggtheme      = theme_bw(),
  tables.theme = theme_cleantable(),
  ylim         = c(0.87, 1), risk.table = TRUE,
  n.break.y    = 10,
  font.x       = c(15, "bold"), font.y = c(15, "bold"),
  legend       = "bottom",
  legend.title = "Overall VL fatality probability"
)
print(gm)


# ── 6. REVIEWER QUERIES (PHILIPPE & MAKOTO) ───────────────────────────────────

## 6a. Median delay by age group (<15 vs ≥15 years) — Makoto abstract comment
data_prabin_2 <- data_prabin_2 %>%
  mutate(
    age_15 = if_else(age_years < 15, "<15", "≥15"),
    # Fine-grained age grouping for supplementary table
    age_group_4 = case_when(
      age_years >= 0  & age_years < 1  ~ "0-1",
      age_years >= 1  & age_years < 5  ~ "1-5",
      age_years >= 5  & age_years < 15 ~ "6-15",
      age_years >= 15 & age_years < 50 ~ "16-50",
      age_years >= 50                  ~ ">50",
      TRUE ~ NA_character_
    ),
    age_group_4 = factor(age_group_4, levels = c("0-1","1-5","6-15","16-50",">50"))
  )

summary(data_prabin_2$dif_day_tx_symp[data_prabin_2$age_15 == "<15"])
summary(data_prabin_2$dif_day_tx_symp[data_prabin_2$age_15 == "≥15"])


## 6b. Mortality count by delay category
data_prabin_2 %>%
  tabyl(delay_categories, death_vl) %>%
  adorn_totals("both")

delay_cat <- survfit2(
  Surv(dif_treat_outcome, death_vl) ~ delay_categories,
  data = data_prabin_2
)

ggsurvfit(delay_cat, type = "survival") +
  add_risktable(risktable_stats = "{n.risk} ({cum.event})")


## 6c. Gender distribution
data_prabin_2 %>%
  tabyl(gender) %>%
  adorn_totals()


## 6d. Region-specific median delay for the most recent period (2019–2023)
summary(data_prabin_2$dif_day_tx_symp[data_prabin_2$year_3 == "2019 - 2023"])

data_prabin_2 %>%
  filter(year_3 == "2019 - 2023") %>%
  group_by(region) %>%
  summarise(
    n      = n(),
    median = median(dif_day_tx_symp),
    IQL    = quantile(dif_day_tx_symp, 0.25),
    IQU    = quantile(dif_day_tx_symp, 0.75)
  )

## 6e. Median delay to treatment for patients who died within 7 and 14 days
data_prabin_2 %>%
  filter(!is.na(death_vl), evolution %in% c("Cure", "Death for VL")) %>%
  mutate(
    death_7  = if_else(dif_treat_outcome <= 7,  "≤7 days",  ">7 days"),
    death_14 = if_else(dif_treat_outcome <= 14, "≤14 days", ">14 days")
  ) %>%
  group_by(death_14) %>%
  summarise(median_delay = median(dif_day_tx_symp))

data_prabin_2 %>%
  filter(evolution %in% c("Cure", "Death for VL")) %>%
  tabyl(death_vl)


# ── 7. SYMPTOM PROBABILITY MODELS BY HIV STATUS (FIGURES 2 & 3) ───────────────
# For each of 11 VL symptoms, we fit a logistic regression model (lrm from rms)
# with a restricted cubic spline (4 knots) on treatment delay, adjusted for
# age (spline) and gender.  The model is run separately for HIV(-) and HIV(+)
# patients and the two marginal probability curves are overlaid.
# fun = plogis converts the log-odds linear predictor to predicted probability.
# ref.zero = TRUE sets the reference at the median of dif_day_tx_symp.

data_7 <- data_6 %>%
  mutate(
    edema_2 = factor(
      if_else(edema == 1, "Yes oedema", "No oedema"),
      levels = c("Yes oedema", "No oedema")
    ),
    severe  = factor(
      if_else(edema == 1 | hemorrhagic_phenomena == 1, "Severe", "No severe"),
      levels = c("Severe", "No severe")
    ),
    hiv_2   = factor(hiv_2, levels = c("Yes", "No"))
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
    severe_2         = if_else(edema == 1, 1, 0),
    federative_unit  = as.factor(federative_unit),
    # Numeric delay category for ordinal modelling
    dif_treat_start_symp_cat_2 = case_when(
      dif_treat_start_symp_cat == "< 30 days"   ~ 1,
      dif_treat_start_symp_cat == "30-89 days"  ~ 2,
      dif_treat_start_symp_cat == "90-119 days" ~ 3,
      dif_treat_start_symp_cat == "≥120 days"   ~ 4,
      TRUE ~ NA
    )
  )

data_symptoms <- data_8 %>%
  dplyr::select(
    dif_day_tx_symp, fever, pallor, weakness, edema, weight_loss,
    cough_and_dyspnea, increase_infection, spleen_enlargement,
    hemorrhagic_phenomena, liver_enlargement, jaundice, hiv_2,
    age_years, gender
  )

## ── Model function ────────────────────────────────────────────────────────────
# Fits lrm(outcome ~ rcs(delay, 4) + rcs(age, 4) + gender) for HIV(-) and HIV(+)
# subgroups and returns a ggplot object with both curves overlaid.
# datadist must be re-set inside the function because each subset has different
# variable distributions (rms stores it in options, not in the fit object).

fit_symptoms_model <- function(data, outcome, y_label = "Complication", knots = 4) {

  # ── HIV(-) group
  df <- data %>%
    filter(
      !is.na(.data[[outcome]]),
      !is.na(dif_day_tx_symp),
      hiv_2 == "No"
    ) %>%
    dplyr::select(all_of(c(outcome, "dif_day_tx_symp", "age_years", "gender")))

  dd <- datadist(df)
  options(datadist = dd)

  fit <- lrm(
    as.formula(
      paste0(outcome, " ~ rcs(dif_day_tx_symp,", knots,
             ") + rcs(age_years,", knots, ") + gender")
    ),
    data = df
  )

  # ── HIV(+) group
  df_2 <- data %>%
    filter(
      !is.na(.data[[outcome]]),
      !is.na(dif_day_tx_symp),
      hiv_2 == "Yes"
    ) %>%
    dplyr::select(all_of(c(outcome, "dif_day_tx_symp", "age_years", "gender")))

  # Global assignment required because rms reads datadist from options at
  # prediction time (after fit_2 leaves this scope)
  dd <<- datadist(df_2)
  options(datadist = "dd")

  fit_2 <- lrm(
    as.formula(
      paste0(outcome, " ~ rcs(dif_day_tx_symp,", knots,
             ") + rcs(age_years,", knots, ") + gender")
    ),
    data = df_2
  )

  pred   <- Predict(fit,   dif_day_tx_symp, fun = plogis, ref.zero = TRUE)
  pred_2 <- Predict(fit_2, dif_day_tx_symp, fun = plogis, ref.zero = TRUE)

  pred_df   <- as.data.frame(pred)
  pred_df_2 <- as.data.frame(pred_2)

  ggplot(pred_df, aes(dif_day_tx_symp, yhat, color = "HIV(-)")) +
    geom_line(linewidth = 1) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = "CI"), alpha = 0.2, color = NA) +
    geom_line(
      data = pred_df_2,
      aes(dif_day_tx_symp, yhat, color = "HIV(+)"),
      linewidth = 1
    ) +
    geom_ribbon(
      data = pred_df_2,
      aes(ymin = lower, ymax = upper, fill = "CI"),
      alpha = 0.2, color = NA
    ) +
    labs(
      x     = "Time to treatment initiation (days)",
      y     = bquote(bolditalic("P") ~ bold("(" * .(y_label) * ")")),
      color = "Group", fill = NULL
    ) +
    theme_bw() +
    theme(
      axis.title   = element_text(size = 12, color = "black", face = "bold"),
      axis.text.x  = element_text(size = 11, color = "black")
    ) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_color_manual(values = c("HIV(-)" = "darkgreen", "HIV(+)" = "brown2")) +
    scale_fill_manual(values = c("CI" = "grey80"))
}

## ── Outcome specifications ────────────────────────────────────────────────────
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

## ── Run models ────────────────────────────────────────────────────────────────
# pmap iterates over rows of symptoms_specs, passing each (outcome, y_label)
# pair together with the shared data argument.
symptoms_plots <- purrr::pmap(
  symptoms_specs,
  fit_symptoms_model,
  data = data_symptoms
)
names(symptoms_plots) <- symptoms_specs$y_label

## ── Strip x-axis from all but the bottom row (rows 1–3 of a 4×3 grid) ────────
# Panels 9–11 (bottom row) keep the x-axis; the rest share it via annotation.
bottom_idx <- 9:11

symptoms_plots2 <- purrr::map2(
  symptoms_plots,
  seq_along(symptoms_plots),
  ~ if (.y %in% bottom_idx) .x else .x +
    theme(
      axis.title.x = element_blank(),
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank()
    )
)

## ── Extract shared legend from the first panel ────────────────────────────────
legend <- cowplot::get_legend(
  symptoms_plots2[[1]] +
    theme(
      legend.position = "right",
      legend.text     = element_text(size = 12, color = "black"),
      legend.title    = element_text(size = 16, color = "black", face = "bold")
    )
)

# Remove per-panel legends before arranging (legend is drawn separately below)
plots_no_legend <- purrr::map(symptoms_plots2, ~ .x + theme(legend.position = "none"))

## ── Assemble panel and place legend in empty bottom-right cell ────────────────
main_panel <- ggpubr::ggarrange(
  plotlist = plots_no_legend,
  nrow = 4, ncol = 3
)

final_severity <- ggdraw(main_panel) +
  draw_grob(legend, x = 0.78, y = 0.02, width = 0.20, height = 0.20)

ggsave(
  filename = "distribution_symptoms_final.png",
  plot     = final_severity,
  path     = here("figures and tables"),
  width = 10, height = 8, dpi = 300
)

## ── Second pmap: data_symptoms passed unfiltered so the function handles
#  the HIV(-)/HIV(+) split internally (filtering upstream would empty one arm)
symptoms_plots_2 <- purrr::pmap(
  symptoms_specs,
  fit_symptoms_model,
  data = data_symptoms
)
names(symptoms_plots_2) <- symptoms_specs$y_label

symptoms_plots3 <- purrr::map2(
  symptoms_plots_2,
  seq_along(symptoms_plots_2),
  ~ if (.y %in% bottom_idx) .x else .x +
    theme(
      axis.title.x = element_blank(),
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank()
    )
)

legend_2 <- cowplot::get_legend(
  symptoms_plots3[[1]] +
    theme(
      legend.position = "right",
      legend.text     = element_text(size = 12, color = "black"),
      legend.title    = element_text(size = 16, color = "black", face = "bold")
    )
)

plots_no_legend_2 <- purrr::map(symptoms_plots3, ~ .x + theme(legend.position = "none"))

final_severity_2 <- ggpubr::ggarrange(
  plotlist      = symptoms_plots3,
  nrow = 4, ncol = 3,
  common.legend = FALSE
)

ggsave(
  filename = "distribution_symptoms_HIV.png",
  plot     = final_severity_2,
  path     = here("figures and tables"),
  width = 10, height = 8, dpi = 300
)

# Combined figure (all symptoms, both HIV groups side by side)
ggarrange(final_severity, final_severity_2, nrow = 2)


# ── 8. ADJUSTED COX MODELS STRATIFIED BY HIV STATUS (FIGURE 5) ───────────────
# cph() from rms is preferred over coxph() here because it integrates with
# rms::Predict(), which computes adjusted HRs over the continuous delay spline
# while holding other predictors at reference levels.
# x = TRUE, y = TRUE stores design matrix and response for Predict().
# ref.zero = TRUE sets HR reference at the median of dif_day_tx_symp.

## ── Shared factor-level ordering ─────────────────────────────────────────────
# Reference category is the first level; chosen to match clinical baseline
# (adult age, female, Northeast region — highest VL burden in Brazil).
level_orders_cox <- list(
  age_group_2     = c("16-50", "0-1", "1-5", "6-15", ">50"),
  gender          = c("Female", "Male"),
  region          = c("Northeast", "Center-West", "North", "Southeast", "South"),
  education_level = c(
    "Complete elementary and incomplete middle school",
    "Incomplete elementary school",
    "Complete middle school and incomplete high school",
    "Complete high school", "Incomplete university",
    "Complete university", "N/A"
  ),
  race_2 = c("Brown", "White", "Afro-descendants", "Yellow", "Indigenous", "Unknown")
)

var_order_cox <- c("age_group_2", "gender", "region", "education_level", "race_2")

## 8a. HIV-positive patients ───────────────────────────────────────────────────
data_mort_hiv <- data_6 %>%
  filter(!is.na(death_vl), !is.na(dif_day_tx_symp)) %>%
  dplyr::select(
    dif_day_tx_symp, death_vl, dif_treat_outcome,
    age_group_2, gender, region, education_level, race_2, hiv_2
  ) %>%
  filter(hiv_2 == "Yes")

for (v in names(level_orders_cox)) {
  data_mort_hiv[[v]] <- factor(data_mort_hiv[[v]], levels = level_orders_cox[[v]])
}

dd <- datadist(data_mort_hiv)
options(datadist = "dd")

m_hiv <- cph(
  Surv(dif_treat_outcome, death_vl) ~
    rcs(dif_day_tx_symp, 4) + age_group_2 + gender + region + education_level + race_2,
  data = data_mort_hiv,
  x = TRUE, y = TRUE
)
summary(m_hiv)

pred_hiv    <- Predict(m_hiv, dif_day_tx_symp, fun = exp, ref.zero = TRUE)
pred_hiv_df <- as.data.frame(pred_hiv)

cox_hiv <- ggplot(pred_hiv_df, aes(x = dif_day_tx_symp, y = yhat)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_vline(
    aes(xintercept = median(data_mort_hiv$dif_day_tx_symp, na.rm = TRUE)),
    color = "red", linetype = "dashed"
  ) +
  labs(
    x = "Days from symptom onset to treatment",
    y = paste0(
      "Hazard ratio (reference = ",
      median(data_mort_hiv$dif_day_tx_symp[data_mort_hiv$hiv_2 == "Yes"], na.rm = TRUE),
      " days)"
    )
  ) +
  theme_bw() +
  theme(
    axis.title   = element_text(size = 14, color = "black", face = "bold"),
    axis.text.x  = element_text(size = 12, color = "black")
  ) +
  scale_x_continuous(n.breaks = 10) +
  scale_y_continuous(n.breaks = 12, limits = c(0.4, 5))

ggsave(cox_hiv,
       path = here("figures and tables"),
       width = 10, height = 8, dpi = 300,
       file = "cox_model_hiv.png")


## 8b. HIV-negative patients ───────────────────────────────────────────────────
data_mort_nohiv <- data_6 %>%
  filter(!is.na(death_vl), !is.na(dif_day_tx_symp)) %>%
  dplyr::select(
    dif_day_tx_symp, death_vl, dif_treat_outcome,
    age_group_2, gender, region, education_level, race_2, hiv_2
  ) %>%
  filter(hiv_2 == "No")

for (v in names(level_orders_cox)) {
  data_mort_nohiv[[v]] <- factor(data_mort_nohiv[[v]], levels = level_orders_cox[[v]])
}

dd <- datadist(data_mort_nohiv)
options(datadist = "dd")

m_nohiv <- cph(
  Surv(dif_treat_outcome, death_vl) ~
    rcs(dif_day_tx_symp, 4) + age_group_2 + gender + region + education_level + race_2,
  data = data_mort_nohiv,
  x = TRUE, y = TRUE
)
summary(m_nohiv)

pred_nohiv    <- Predict(m_nohiv, dif_day_tx_symp, fun = exp, ref.zero = TRUE)
pred_nohiv_df <- as.data.frame(pred_nohiv)

cox_nohiv <- ggplot(pred_nohiv_df, aes(x = dif_day_tx_symp, y = yhat)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_vline(
    aes(xintercept = median(data_mort_nohiv$dif_day_tx_symp, na.rm = TRUE)),
    color = "red", linetype = "dashed"
  ) +
  labs(
    x = "Days from symptom onset to treatment",
    y = paste0(
      "Hazard ratio (reference = ",
      median(data_mort_nohiv$dif_day_tx_symp[data_mort_nohiv$hiv_2 == "No"], na.rm = TRUE),
      " days)"
    )
  ) +
  theme_bw() +
  theme(
    axis.title   = element_text(size = 14, color = "black", face = "bold"),
    axis.text.x  = element_text(size = 12, color = "black")
  ) +
  scale_x_continuous(n.breaks = 10) +
  scale_y_continuous(n.breaks = 12, limits = c(0.4, 5))

ggsave(cox_nohiv,
       path = here("figures and tables"),
       width = 10, height = 8, dpi = 300,
       file = "cox_model_nohiv.png")

## 8c. Combined panel (Figure 5)
ggarrange(
  cox_nohiv, cox_hiv,
  labels  = c("A. No HIV", "B. HIV"),
  label.x = 0.2, label.y = 0.99
)
ggsave(
  path = here("figures and tables"),
  width = 10, height = 8, dpi = 300,
  file  = "cox_model.png"
)


# ── 9. BASELINE CHARACTERISTICS TABLE BY DELAY INTERVAL (TABLE 1) ─────────────
# Median [IQR] and range for three time intervals (ts→td, td→tt, ts→tt)
# stratified by patient characteristic.  This version uses the full delay
# range (no HIV filter) as requested for the revised Table 1.

var_order <- c(
  "age_group_2", "gender", "year_period", "region",
  "education_level", "race_2", "hiv_2", "evolution"
)

# Factor-level ordering for this table (reference = first level).
# NOTE: this differs from level_orders_cox above (different reference chosen
# for the Cox model vs. the descriptive table).
level_orders_tab1 <- list(
  age_group_2     = c("0-1", "1-5", "6-15", "16-50", ">50"),
  gender          = c("Female", "Male"),
  year_period     = c("2007-2009", "2010-2014", "2015-2018", "2019-2023"),
  region          = c("Center-West", "North", "Northeast", "Southeast", "South"),
  education_level = c(
    "Incomplete elementary school",
    "Complete elementary and incomplete middle school",
    "Complete middle school and incomplete high school",
    "Complete high school", "Incomplete university",
    "Complete university", "N/A"
  ),
  race_2    = c("Brown", "White", "Afro-descendants", "Yellow", "Indigenous", "Unknown"),
  hiv_2     = c("No", "Yes", "Unknown"),
  evolution = c("Cure", "Death for VL", "Death for any-cause", "Neglect", "Transferred", NA)
)

## ── Function: median + IQR + range for ts→td, td→tt, ts→tt ─────────────────
median_time_order <- function(vars, level_orders = list(), data = data_6) {
  map_dfr(vars, function(v) {

    d_sub <- data %>%
      filter(
        !is.na(.data[[v]]),
        !is.na(dif_day_tx_symp),
        dif_day_tx_symp > 0
      )

    if (nrow(d_sub) == 0) return(NULL)

    if (!is.null(level_orders[[v]])) {
      d_sub[[v]] <- factor(d_sub[[v]], levels = level_orders[[v]])
    } else {
      d_sub[[v]] <- factor(d_sub[[v]])
    }

    denom <- nrow(d_sub)

    df_sum <- d_sub %>%
      group_by(.data[[v]]) %>%
      summarise(
        n = n(),
        median_1 = paste0(
          round(median(dif_day_dx_symp,  na.rm = TRUE), 0), " [",
          round(quantile(dif_day_dx_symp, 0.25, na.rm = TRUE), 0), "-",
          round(quantile(dif_day_dx_symp, 0.75, na.rm = TRUE), 0), "]"
        ),
        range_1 = paste0(
          min(dif_day_dx_symp, na.rm = TRUE), "-",
          max(dif_day_dx_symp, na.rm = TRUE)
        ),
        median_2 = paste0(
          round(median(dif_day_tx_dx,  na.rm = TRUE), 0), " [",
          round(quantile(dif_day_tx_dx, 0.25, na.rm = TRUE), 0), "-",
          round(quantile(dif_day_tx_dx, 0.75, na.rm = TRUE), 0), "]"
        ),
        range_2 = paste0(
          min(dif_day_tx_dx, na.rm = TRUE), "-",
          max(dif_day_tx_dx, na.rm = TRUE)
        ),
        median_3 = paste0(
          round(median(dif_day_tx_symp,  na.rm = TRUE), 0), " [",
          round(quantile(dif_day_tx_symp, 0.25, na.rm = TRUE), 0), "-",
          round(quantile(dif_day_tx_symp, 0.75, na.rm = TRUE), 0), "]"
        ),
        range_3 = paste0(
          min(dif_day_tx_symp, na.rm = TRUE), "-",
          max(dif_day_tx_symp, na.rm = TRUE)
        ),
        .groups = "drop"
      ) %>%
      mutate(
        variable = v,
        level    = as.character(.data[[v]]),
        p = paste0(
          scales::percent(n / denom, accuracy = 0.01, suffix = "%"),
          "(", scales::comma(n, big.mark = ","), "/",
          scales::comma(denom, big.mark = ","), ")"
        )
      ) %>%
      dplyr::select(variable, level, p,
                    median_1, range_1,
                    median_2, range_2,
                    median_3, range_3)

    if (!is.null(level_orders[[v]])) {
      df_sum <- df_sum %>%
        mutate(level = factor(level, levels = level_orders[[v]]))
    } else {
      df_sum <- df_sum %>% mutate(level = factor(level))
    }

    df_sum
  }) %>%
    mutate(variable = factor(variable, levels = var_order)) %>%
    arrange(variable, level)
}

dif_table <- median_time_order(var_order, level_orders_tab1, data = data_6)

border_line <- officer::fp_border(color = "black", width = 1.38)

tab_1 <- dif_table %>%
  mutate(
    variable = recode(
      variable,
      age_group_2     = "Age group",
      gender          = "Gender",
      year_period     = "Period (year)",
      region          = "Region",
      education_level = "Education level",
      race_2          = "Race/colour",
      hiv_2           = "HIV status",
      evolution       = "Evolution"
    )
  ) %>%
  rename(
    "Characteristic"           = variable,
    "Level"                    = level,
    "% (n/N)"                  = p,
    "Median time ts-td (days)" = median_1,
    "Range ts-td (days)"       = range_1,
    "Median time td-tt (days)" = median_2,
    "Range td-tt (days)"       = range_2,
    "Median time ts-tt (days)" = median_3,
    "Range ts-tt (days)"       = range_3
  ) %>%
  flextable() %>%
  # Merge characteristic column for each variable group
  merge_at(j = 1, i = 1:5,   part = "body") %>%
  merge_at(j = 1, i = 6:7,   part = "body") %>%
  merge_at(j = 1, i = 8:11,  part = "body") %>%
  merge_at(j = 1, i = 12:16, part = "body") %>%
  merge_at(j = 1, i = 17:23, part = "body") %>%
  merge_at(j = 1, i = 24:29, part = "body") %>%
  merge_at(j = 1, i = 30:32, part = "body") %>%
  merge_at(j = 1, i = 33:37, part = "body") %>%
  hline(i = c(5, 7, 11, 16, 23, 29, 32), part = "body", border = border_line) %>%
  align(j = 3:9, part = "body",   align = "center") %>%
  align(j = 3:9, part = "header", align = "center") %>%
  bold(j = 1, part = "all") %>%
  bold(i = 1, part = "header") %>%
  autofit()

tab_1

tab_1 %>%
  save_as_docx(path = here("figures and tables", "table_1.docx"))
