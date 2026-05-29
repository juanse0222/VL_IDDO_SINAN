pacman::p_load(flextable, tidyverse, janitor, epikit, gtsummary, rio, 
               survival, survminer, sf, here, incidence, ggstance, ggbeeswarm, 
               grid, ggsurvfit, geobr, rms, patchwork, biostat3, simukde)

################# Manuscript adjusts - Requririements #####################################################################

data_prabin <- import(here("data base and clean", "data_6.rds"))

data_prabin_2 <- data_prabin %>% 
  mutate(delay_2 = case_when(dif_day_tx_symp <= 14 ~ "0-14 days",
                             dif_day_tx_symp > 14 & dif_day_tx_symp <= 30 ~ "15-30 days",
                             TRUE ~ dif_treat_start_symp_cat),
         delay_2 = factor(delay_2, levels = c("0-14 days",
                                              "15-30 days",
                                              "30-89 days",
                                              "90-119 days",
                                              "≥120 days")),
         delay_categories = case_when(dif_day_tx_symp >= 0 & dif_day_tx_symp <= 14 ~ "≤14 days",
                                             dif_day_tx_symp >= 15 & dif_day_tx_symp <= 30 ~ "15 to ≤ 30 days",
                                             dif_day_tx_symp >= 31 & dif_day_tx_symp <= 60 ~ "31 to ≤ 60 days",
                                             dif_day_tx_symp >= 61 & dif_day_tx_symp <= 90 ~ "61 to ≤ 90 days",
                                             dif_day_tx_symp >= 91 & dif_day_tx_symp <= 120 ~ "91 to ≤ 120 days",
                                             dif_day_tx_symp >= 121 & dif_day_tx_symp <= 150 ~ "121 to ≤ 150 days",
                                             dif_day_tx_symp >= 151 & dif_day_tx_symp <= 180 ~ "151 to ≤ 180 days",
                                             dif_day_tx_symp >= 181 ~ "≥6 months",
                                             TRUE ~ NA),
         delay_categories = factor(delay_categories, levels = c("≤14 days", 
                                                                       "15 to ≤ 30 days",
                                                                       "31 to ≤ 60 days",
                                                                       "61 to ≤ 90 days",
                                                                       "91 to ≤ 120 days",
                                                                       "121 to ≤ 150 days",
                                                                       "151 to ≤ 180 days",
                                                                       "≥6 months")),
         year_3 = case_when(year_of_notification >= 2007 & year_of_notification <= 2009 ~ "2007 - 2009",
                            year_of_notification >= 2010 & year_of_notification <= 2014 ~ "2010 - 2014",
                            year_of_notification >= 2015 & year_of_notification <= 2018 ~ "2015 - 2018",
                            year_of_notification >= 2019 & year_of_notification <= 2023 ~ "2019 - 2023",
                             TRUE ~ NA),
         region = case_when(
           federative_unit %in% c(11, 13, 14, 15, 16, 17) ~ "North",
           federative_unit %in% c(21, 22, 23, 24, 25, 26, 27, 28, 29) ~ "Northeast",
           federative_unit %in% c(31, 32, 33, 35)         ~ "Southeast",
           federative_unit %in% c(41, 42, 43)             ~ "South",
           federative_unit %in% c(50, 51, 52, 53)         ~ "Center-West",
           TRUE                                           ~ NA_character_))

data_prabin_2 <- data_prabin_2 %>%
  dplyr::mutate(across(where(~ inherits(.x, "IDate")), as.Date))

################ Median and descriptive values by region ##################################################################
data_6 <- import(here("data base and clean", "data_6.rds"))  ### Clean database
data_6 <- data_6 %>%
  dplyr::mutate(across(where(~ inherits(.x, "IDate")), as.Date))

data_8 <- data_8 %>%
  dplyr::mutate(across(where(~ inherits(.x, "IDate")), as.Date))


data_prabin_2 %>% 
  group_by(region) %>% 
  summarise(n = n(),
         median_r = median(dif_day_tx_symp),
         IQL = quantile(dif_day_tx_symp, 0.25),
         IQS = quantile(dif_day_tx_symp, 0.75)) %>% 
  adorn_totals()


data_6 %>% 
  filter(!is.na(name_tx)) %>% 
  filter(!is.na(death_vl)) %>% 
  #filter(evolution %in% c("Cure", "Death for VL")) %>% 
  mutate(days_p = case_when(dif_treat_outcome <= 7 ~ "≤7 days",
                            dif_treat_outcome >= 8 &  dif_treat_outcome <= 14 ~ "7-14 days",
                            TRUE ~ "Others"),
         days_p = factor(days_p, levels = c("≤7 days", "7-14 days", "Others"))) %>% 
  group_by(days_p, evolution) %>% 
  summarise(n = n(), .groups = "drop") %>% 
  pivot_wider(id_cols = "days_p", names_from = "evolution", values_from = n) %>% 
  adorn_totals("both") %>% 
  adorn_percentages(denominator = "all") %>% 
  adorn_pct_formatting() %>% 
  adorn_ns(position = "front")

data_6 %>% 
  tabyl(evolution)

############## Survival double check evaluation ###########################################################################

### Delay symptoms onset to treatment adjusted considering people under 30 days (14-30)
delay <- survfit2(Surv(dif_treat_outcome, death_vl) ~ delay_2, data = data_prabin_2) 

summary(delay, times = c(180), data.frame = T)

gd <- ggsurvfit(delay, linewidth = 1) +
  #add_confidence_interval() +
  add_risktable(risktable_stats = "{n.risk} ({cum.event})", 
                theme = theme_bw(), risktable_height = 0.2) +
  add_censor_mark() +
  add_pvalue(caption="Log-rank {p.value}",
          location = "annotation", x = 200, y = 0.98) +
  scale_ggsurvfit() + 
  coord_cartesian(ylim = c(0.85, 1)) +
  scale_x_continuous(n.breaks = 15) +
  scale_y_continuous(n.breaks = 12, label = scales::label_percent()) +
  add_quantile(x_value = 30, linetype = "dotted", 
               color = "grey30", linewidth = 0.8) +
  #geom_vline(xintercept = 30, linetype = "dashed", colour = "grey30") +
  labs(x = "Time (days after to start treatment)") +
  theme_bw() +
  theme(legend.position = "bottom", 
        axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(color = "black", size = 14, face = "bold"),
        legend.title = element_text(face = "bold", size = 14, color = "black"),
        legend.text = element_text(color = "black", size = 12),
  ) 

print(gd)  


### Makoto query # people evaluated within the first six months 

data_prabin_2 %>% 
  filter(!is.na(death_vl)) %>% 
  filter(dif_treat_outcome < 181) %>% 
  count()

data_prabin_3 <- data_6 %>% 
  filter(delay_categories %in% c("≤14 days", "≥6 months"))


survdiff(Surv(dif_treat_outcome, death_vl) ~ delay_categories, rho = 1, data = data_prabin_3)


delay <- survfit2(Surv(dif_treat_outcome, death_vl) ~ delay_categories, data = data_prabin_3)

summary(delay, times = c(360), data.frame = T)

gd <- ggsurvfit(delay, linewidth = 1) +
  #add_confidence_interval() +
  add_risktable(risktable_stats = "{n.risk} ({cum.event})", 
                theme = theme_bw(), risktable_height = 0.2) +
  add_censor_mark() +
  add_pvalue(caption="Log-rank {p.value}",
             location = "annotation", x = 200, y = 0.98) +
  scale_ggsurvfit() + 
  coord_cartesian(ylim = c(0.85, 1)) +
  scale_x_continuous(n.breaks = 15) +
  scale_y_continuous(n.breaks = 12, label = scales::label_percent()) +
  add_quantile(x_value = 30, linetype = "dotted", 
               color = "grey30", linewidth = 0.8) +
  #geom_vline(xintercept = 30, linetype = "dashed", colour = "grey30") +
  labs(x = "Time (days after to start treatment)") +
  theme_bw() +
  theme(legend.position = "bottom", 
        axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(color = "black", size = 14, face = "bold"),
        legend.title = element_text(face = "bold", size = 14, color = "black"),
        legend.text = element_text(color = "black", size = 12),
)

#### Upon mortality after treatment iniatiation within the first 7 days. 

fatality_risk <- survfit(Surv(dif_treat_outcome, death_vl) ~ 1, data = data_6) 

summary(fatality_risk, times = c(0,30,60,90,180,365))

gm <- ggsurvplot(fit = lethality_risk,
                 data = data_6,
                 ggtheme = theme_bw(),
                 tables.theme = theme_cleantable(),
                 ylim = c(0.87, 1), risk.table = T, 
                 n.break.y = 10,
                 font.x = c(15, "bold"), font.y = c(15, "bold"), 
                 # xlab = "days (delay symptoms onset to treatment)",
                 legend = c("bottom"), legend.title = "Survival regarding the delay symptoms onset to treatment") 

print(gm)

data_prabin_2 %>% 
  filter(!is.na(death_vl)) %>% 
  filter(evolution %in% c("Cure", "Death for VL")) %>% 
  mutate(death_7 = if_else(dif_treat_outcome <= 7, "7 days", "≥ 7 days"),
         death_14 = if_else(dif_treat_outcome <= 14, "14 days", "≥ 14 days")) %>% 
  group_by(death_14) %>% 
  summarise(median_delay = median(dif_day_tx_symp))



asddata_prabin_2 %>% 
  filter(evolution %in% c("Cure", "Death for VL")) %>% 
  tabyl(death_vl) %>% 

######################## Adjustment of the model ######################################################

#### COX model patients with HIV

# Prepare data
# ----------------
data_mort <- data_6 %>%
  filter(!is.na(death_vl),
         !is.na(dif_day_tx_symp)) %>% 
  dplyr::select(
    dif_day_tx_symp, death_vl, dif_treat_outcome,
    age_group_2, gender, region,
    education_level, race_2, hiv_2
  ) %>% 
  filter(hiv_2 == "Yes")

# ----------------
# Variable order
# ----------------
var_order_2 <- c(
  "age_group_2", "gender", "region",
  "education_level", "race_2"
)

# ----------------
# Explicit factor level ordering (reference = first)
# ----------------
level_orders <- list(
  age_group_2 = c("16-50", "0-1", "1-5", "6-15", ">50"),
  gender      = c("Female", "Male"),
  region      = c("Northeast", "Center-West", "North", "Southeast", "South"),
  education_level = c(
    "Complete elementary and incomplete middle school",
    "Incomplete elementary school",
    "Complete middle school and incomplete high school",
    "Complete high school",
    "Incomplete university",
    "Complete university",
    "N/A"
  ),
  race_2 = c("Brown", "White", "Afro-descendants", "Yellow", "Indigenous", "Unknown")
)

for (v in names(level_orders)) {
  data_mort[[v]] <- factor(data_mort[[v]], levels = level_orders[[v]])
}

# ----------------
# datadist (required by rms)
# ----------------
dd <- datadist(data_mort)
options(datadist = "dd")

# ----------------
# Fit Cox model with restricted cubic spline
# ----------------
m1 <- survival::coxph(
  Surv(dif_treat_outcome, death_vl) ~
    rcs(dif_day_tx_symp, 4) +
    age_group_2 + gender + region +
    education_level + race_2,
  data = data_mort,
  x = TRUE,
  y = TRUE
)

cox_table <- tidy(m1, exponentiate = TRUE, conf.int = TRUE)


tbl_regression(m1, exponentiate = TRUE)
plot(Predict(m1, dif_day_tx_symp, fun = exp))

remove.packages("rms")
install.packages("rms", dependencies = TRUE)

#=======================================================
# Comments by Philippe and Makoto 
#=======================================================
### 1. ge group <15 vs >15 abstract comment by Makoto

age_group_2 = c("0-1", "1-5", "6-15", "16-50", ">50")

data_prabin_2 <- data_prabin_2 %>%  
  mutate(
    age_15 = if_else(age_years < 15, "<15", "≥15"),
    age_group_4 = case_when(
      age_years >= 0 & age_years < 1  ~ "0-1",
      age_years >= 1 & age_years < 5  ~ "1-5",
      age_years >= 5 & age_years < 15 ~ "6-15",
      age_years >= 15 & age_years < 50 ~ "16-50",
      age_years >= 50 ~ ">50",
      TRUE ~ NA_character_
    ),
    age_group_4 = factor(age_group_4, levels = c("0-1", "1-5", "6-15", "16-50", ">50"))
  )

summary(data_prabin_2$dif_day_tx_symp[data_prabin_2$age_15 == "<15"])
summary(data_prabin_2$dif_day_tx_symp[data_prabin_2$age_15 == "≥15"])

### 2. Delay group with more mortality. 

data_prabin_2 %>% 
  tabyl(delay_categories, death_vl) %>% 
  adorn_totals("both")

delay_cat <- survfit2(Surv(dif_treat_outcome, death_vl) ~ delay_categories, data = data_prabin_2) 
times <- c(30, 60, 90, 180, 320) 

ggsurvfit(delay_cat, type = "survival") +
  #add_confidence_interval() +
  add_risktable(risktable_stats = "{n.risk} ({cum.event})") 

### 3. Age groups distribution and percentage

data_prabin_2 %>% 
  tabyl(gender) %>% 
  adorn_totals()

### 4. Region time to treatment for the last period 

summary(data_prabin_2$dif_day_tx_symp[data_prabin_2$year_3 == "2019-2023"])

data_prabin_2 %>% 
  filter(year_3 == "2019 - 2023") %>% 
  group_by(region) %>% 
  summarise(n = n(),
            median = median(dif_day_tx_symp),
            IQL = quantile(dif_day_tx_symp, 0.25),
            IQU = quantile(dif_day_tx_symp, 0.75))


################################### Model un-adjusted by symptoms - Figure 2 and 3 ##########################################################################

data_7 <- data_6 %>% 
  mutate(edema_2 = if_else(edema == 1, "Yes oedema", "No oedema"),
         edema_2 = factor(edema_2, levels = c("Yes oedema", "No oedema")), 
         severe = if_else(edema == 1 | hemorrhagic_phenomena == 1, "Severe", "No severe"),
         severe = factor(severe, levels = c("Severe", "No severe")),
         hiv_3 = replace_na(hiv_2 == "Unknown"),
         hiv_2 = factor(hiv_2, levels = c("Yes", "No")))

data_8 <- data_7 %>% 
  mutate(hiv_2 = if_else(hiv_2 == "Unknown", NA, hiv_2)) %>% 
  mutate(gender_2 = if_else(gender == "Male", 1, 0),
         hiv_3 = case_when(hiv_2 == "Yes" ~ 1, 
                           hiv_2 == "No" ~ 0,
                           hiv_2 == "Unknown" ~ 2,
                           TRUE ~ NA),
         severe_2 = if_else(edema == 1, 1, 0),
         federative_unit = factor(federative_unit),
         gender_2 = factor(gender_2),
         hiv_3 = factor(hiv_3),
         federative_unit = as.factor(federative_unit)) %>% 
  mutate(dif_treat_start_symp_cat_2 = case_when(dif_treat_start_symp_cat == "< 30 days" ~ 1,
                                                dif_treat_start_symp_cat == "30-89 days" ~ 2,
                                                dif_treat_start_symp_cat == "90-119 days" ~ 3,
                                                dif_treat_start_symp_cat == "≥120 days" ~ 4,
                                                TRUE ~ NA))

data_symptoms <- data_8 %>% 
  dplyr::select(dif_day_tx_symp, fever, pallor, weakness, edema, weight_loss, cough_and_dyspnea,
                increase_infection, spleen_enlargement, hemorrhagic_phenomena, liver_enlargement, jaundice, hiv_2)

fit_symptoms_model <- function(
    data,
    outcome,
    y_label = "Complication",
    knots = 4
) {
  
  df <- data %>%
    filter(!is.na(.data[[outcome]]),
           !is.na(dif_day_tx_symp),
           hiv_2 == "No") %>% 
    dplyr::select(all_of(c(outcome, "dif_day_tx_symp"))) 
  
  dd <- datadist(df)
  options(datadist = "dd")
  
  fit <- lrm(
    as.formula(
      paste0(outcome, " ~ rcs(dif_day_tx_symp, ", knots, ")")),
    data = df 
  )
  
  df_2 <- data %>%
    filter(!is.na(.data[[outcome]]),
           !is.na(dif_day_tx_symp),
           hiv_2 == "Yes") %>% 
    dplyr::select(all_of(c(outcome, "dif_day_tx_symp")))
  
  dd <- datadist(df_2)
  options(datadist = "dd")
  
  fit_2 <- lrm(
    as.formula(
      paste0(outcome, " ~ rcs(dif_day_tx_symp, ", knots, ")")),
    data = df_2 
  )
  
  pred <- Predict( fit, dif_day_tx_symp,
                   fun = plogis,
                   ref.zero = TRUE
  )
  
  pred_2 <- Predict( fit_2, dif_day_tx_symp,
                     fun = plogis,
                     ref.zero = TRUE
                     
  )                 
  
  pred_df <- as.data.frame(pred)
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
      alpha = 0.2,
      color = NA
    ) +
    labs(
      x = "Delay to treatment initiation (days)",
      y = bquote(bolditalic("P(x)") ~ bold("(" * .(y_label) * ")")),
      color = "Group",
      fill = NULL
    ) +
    theme_bw() +
    theme(
      axis.title = element_text(size = 12, color = "black", face = "bold"),
      axis.text.x = element_text(size = 11, color = "black")
    ) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_color_manual(values = c("HIV(-)" = "darkgreen", "HIV(+)" = "brown2")) +
    scale_fill_manual(values = c("CI" = "grey80"))
}

symptoms_specs <- tibble::tribble(
  ~outcome,               ~y_label,
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

symptoms_plots <- purrr::pmap(
  symptoms_specs,
  fit_symptoms_model,
  data = data_symptoms 
  #filter(hiv_2 == "No")
)

symptoms_plots[[1]]  # == severity_m1
symptoms_plots[[2]]  # == severity_m2
symptoms_plots[[3]]  # == severity_m3
symptoms_plots[[4]]  # == severity_m1
symptoms_plots[[5]]  # == severity_m2
symptoms_plots[[6]]  # == severity_m3
symptoms_plots[[7]]  # == severity_m1
symptoms_plots[[8]]  # == severity_m2
symptoms_plots[[9]]  # == severity_m3
symptoms_plots[[10]]  # == severity_m1
symptoms_plots[[11]]  # == severity_m2

names(symptoms_plots) <- symptoms_specs$y_label

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

library(cowplot)

# 1) Extract legend from one plot
legend <- cowplot::get_legend(
  symptoms_plots2[[1]] +
    theme(legend.position = "right",
          legend.text = element_text(size = 12, color = "black"),
          legend.title = element_text(size = 16, color = "black", face = "bold"))
)

# 3) Arrange plots
main_panel <- ggpubr::ggarrange(
  plotlist = plots_no_legend,
  nrow = 4,
  ncol = 3
)

# 4) Put legend inside the figure at exact position
final_severity <- ggdraw(main_panel) +
  draw_grob(legend, x = 0.78, y = 0.02, width = 0.20, height = 0.20)


ggsave(
  filename = "distribution_symptoms_final.png",
  plot = final_severity,
  path = here("figures and tables"),
  width = 10, height = 8, dpi = 300
)

#### Symptoms model with HIV

symptoms_plots_2 <- purrr::pmap(
  symptoms_specs,
  fit_symptoms_model,
  data = data_symptoms %>% filter(hiv_2 == "Yes")
)

symptoms_plots_2[[1]]  # == severity_m1
symptoms_plots_2[[2]]  # == severity_m2
symptoms_plots_2[[3]]  # == severity_m3
symptoms_plots_2[[4]]  # == severity_m1
symptoms_plots_2[[5]]  # == severity_m2
symptoms_plots_2[[6]]  # == severity_m3
symptoms_plots_2[[7]]  # == severity_m1
symptoms_plots_2[[8]]  # == severity_m2
symptoms_plots_2[[9]]  # == severity_m3
symptoms_plots_2[[10]]  # == severity_m1
symptoms_plots_2[[11]]  # == severity_m2

names(symptoms_plots_2) <- symptoms_specs$y_label

bottom_idx <- 9:11

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

final_severity_2 <- ggpubr::ggarrange(
  plotlist = symptoms_plots3,
  nrow = 4, ncol = 3,
  common.legend = FALSE
)

ggsave(
  filename = "distribution_symptoms_HIV.png",
  plot = final_severity_2,
  path = here("figures and tables"),
  width = 10, height = 8, dpi = 300
)

ggarrange(final_severity, final_severity_2, nrow = 2)

################################### Model un-adjusted by symptoms - Figure 2 and 3 ##########################################################################

data_7 <- data_6 %>% 
  mutate(edema_2 = if_else(edema == 1, "Yes oedema", "No oedema"),
         edema_2 = factor(edema_2, levels = c("Yes oedema", "No oedema")), 
         severe = if_else(edema == 1 | hemorrhagic_phenomena == 1, "Severe", "No severe"),
         severe = factor(severe, levels = c("Severe", "No severe")),
         hiv_3 = replace_na(hiv_2 == "Unknown"),
         hiv_2 = factor(hiv_2, levels = c("Yes", "No")))

data_8 <- data_7 %>% 
  mutate(hiv_2 = if_else(hiv_2 == "Unknown", NA, hiv_2)) %>% 
  mutate(gender_2 = if_else(gender == "Male", 1, 0),
         hiv_3 = case_when(hiv_2 == "Yes" ~ 1, 
                           hiv_2 == "No" ~ 0,
                           hiv_2 == "Unknown" ~ 2,
                           TRUE ~ NA),
         severe_2 = if_else(edema == 1, 1, 0),
         federative_unit = factor(federative_unit),
         gender_2 = factor(gender_2),
         hiv_3 = factor(hiv_3),
         federative_unit = as.factor(federative_unit)) %>% 
  mutate(dif_treat_start_symp_cat_2 = case_when(dif_treat_start_symp_cat == "< 30 days" ~ 1,
                                                dif_treat_start_symp_cat == "30-89 days" ~ 2,
                                                dif_treat_start_symp_cat == "90-119 days" ~ 3,
                                                dif_treat_start_symp_cat == "≥120 days" ~ 4,
                                                TRUE ~ NA))

data_symptoms <- data_8 %>% 
  dplyr::select(dif_day_tx_symp, fever, pallor, weakness, edema, weight_loss, cough_and_dyspnea,
                increase_infection, spleen_enlargement, hemorrhagic_phenomena, liver_enlargement,
                jaundice, hiv_2, age_years, gender)

# ── Model function ─────────────────────────────────────────────────────────────
fit_symptoms_model <- function(data, outcome, y_label = "Complication", knots = 4) {
  
  # HIV(-) group
  df <- data %>%
    filter(!is.na(.data[[outcome]]),
           !is.na(dif_day_tx_symp),
           hiv_2 == "No") %>%
    dplyr::select(all_of(c(outcome, "dif_day_tx_symp", "age_years", "gender")))
  
  # ✅ Fixed — pass the object directly, no environment lookup needed
  dd <- datadist(df)
  options(datadist = dd)
  
  # FIX 1: corrected paste0 — proper string concatenation, comma inside rcs()
  fit <- lrm(
    as.formula(
      paste0(outcome, " ~ rcs(dif_day_tx_symp, ", knots, ") + rcs(age_years, ", knots, ") + gender")
    ),
    data = df
  )
  
  # HIV(+) group
  df_2 <- data %>%
    filter(!is.na(.data[[outcome]]),
           !is.na(dif_day_tx_symp),
           hiv_2 == "Yes") %>%
    dplyr::select(all_of(c(outcome, "dif_day_tx_symp", "age_years", "gender")))
  
  dd <<- datadist(df_2)        # ← same here
  options(datadist = "dd")
  
  # FIX 1 (repeated): same formula fix for fit_2
  fit_2 <- lrm(
    as.formula(
      paste0(outcome, " ~ rcs(dif_day_tx_symp, ", knots, ") + rcs(age_years, ", knots, ") + gender")
    ),
    data = df_2
  )
  
  pred <- Predict(fit, dif_day_tx_symp, fun = plogis, ref.zero = TRUE)
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
      x = "Delay to treatment initiation (days)",
      y = bquote(bolditalic("P") ~ bold("(" * .(y_label) * ")")),
      color = "Group",
      fill = NULL
    ) +
    theme_bw() +
    theme(
      axis.title = element_text(size = 12, color = "black", face = "bold"),
      axis.text.x = element_text(size = 11, color = "black")
    ) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_color_manual(values = c("HIV(-)" = "darkgreen", "HIV(+)" = "brown2")) +
    scale_fill_manual(values = c("CI" = "grey80"))
}

# ── Outcome specs ──────────────────────────────────────────────────────────────
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

# ── Run models ─────────────────────────────────────────────────────────────────
symptoms_plots <- purrr::pmap(
  symptoms_specs,
  fit_symptoms_model,
  data = data_symptoms
)

names(symptoms_plots) <- symptoms_specs$y_label

# ── Strip x-axis from all plots except the bottom row ─────────────────────────
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

library(cowplot)

# ── Extract shared legend ──────────────────────────────────────────────────────
legend <- cowplot::get_legend(
  symptoms_plots2[[1]] +
    theme(legend.position = "right",
          legend.text  = element_text(size = 12, color = "black"),
          legend.title = element_text(size = 16, color = "black", face = "bold"))
)

# FIX 2: define plots_no_legend before passing to ggarrange
plots_no_legend <- purrr::map(symptoms_plots2, ~ .x + theme(legend.position = "none"))

# ── Arrange panels ─────────────────────────────────────────────────────────────
main_panel <- ggpubr::ggarrange(
  plotlist = plots_no_legend,
  nrow = 4,
  ncol = 3
)

final_severity <- ggdraw(main_panel) +
  draw_grob(legend, x = 0.78, y = 0.02, width = 0.20, height = 0.20)

ggsave(
  filename = "distribution_symptoms_final.png",
  plot     = final_severity,
  path     = here("figures and tables"),
  width = 10, height = 8, dpi = 300
)

# ── FIX 3: second pmap uses full data (no pre-filter); function handles the
#    HIV split internally — pre-filtering to "Yes" would empty the HIV(-) branch
symptoms_plots_2 <- purrr::pmap(
  symptoms_specs,
  fit_symptoms_model,
  data = data_symptoms        # ← removed the erroneous filter(hiv_2 == "Yes")
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

# ── Extract shared legend ──────────────────────────────────────────────────────
legend_2 <- cowplot::get_legend(
  symptoms_plots3[[1]] +
    theme(legend.position = "right",
          legend.text  = element_text(size = 12, color = "black"),
          legend.title = element_text(size = 16, color = "black", face = "bold"))
)

# FIX 2: define plots_no_legend before passing to ggarrange
plots_no_legend_2 <- purrr::map(symptoms_plots3, ~ .x + theme(legend.position = "none"))

# ── Arrange panels ─────────────────────────────────────────────────────────────
main_panel <- ggpubr::ggarrange(
  plotlist = plots_no_legend_2,
  nrow = 4,
  ncol = 3
)


final_severity_2 <- ggpubr::ggarrange(
  plotlist = symptoms_plots3,
  nrow = 4, ncol = 3,
  common.legend = FALSE
)

ggsave(
  filename = "distribution_symptoms_HIV.png",
  plot     = final_severity_2,
  path     = here("figures and tables"),
  width = 10, height = 8, dpi = 300
)

ggarrange(final_severity, final_severity_2, nrow = 2)


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

final_severity_2 <- ggpubr::ggarrange(
  plotlist = symptoms_plots3,
  nrow = 4, ncol = 3,
  common.legend = FALSE
)

ggsave(
  filename = "distribution_symptoms_HIV.png",
  plot = final_severity_2,
  path = here("figures and tables"),
  width = 10, height = 8, dpi = 300
)

ggarrange(final_severity, final_severity_2, nrow = 2)



############################################################# COX Model Adjusted - Figure 5 ############################################################################

#### COX model patients with HIV

# Prepare data
# ----------------
data_mort <- data_6 %>%
  filter(!is.na(death_vl),
         !is.na(dif_day_tx_symp)) %>% 
  dplyr::select(
    dif_day_tx_symp, death_vl, dif_treat_outcome,
    age_group_2, gender, region,
    education_level, race_2, hiv_2
  ) %>% 
  filter(hiv_2 == "Yes")

# ----------------
# Variable order
# ----------------
var_order_2 <- c(
  "age_group_2", "gender", "region",
  "education_level", "race_2"
)

# ----------------
# Explicit factor level ordering (reference = first)
# ----------------
level_orders <- list(
  age_group_2 = c("16-50", "0-1", "1-5", "6-15", ">50"),
  gender      = c("Female", "Male"),
  region      = c("Northeast", "Center-West", "North", "Southeast", "South"),
  education_level = c(
    "Complete elementary and incomplete middle school",
    "Incomplete elementary school",
    "Complete middle school and incomplete high school",
    "Complete high school",
    "Incomplete university",
    "Complete university",
    "N/A"
  ),
  race_2 = c("Brown", "White", "Afro-descendants", "Yellow", "Indigenous", "Unknown")
)

for (v in names(level_orders)) {
  data_mort[[v]] <- factor(data_mort[[v]], levels = level_orders[[v]])
}

# ----------------
# datadist (required by rms)
# ----------------
dd <- datadist(data_mort)
options(datadist = "dd")

# ----------------
# Fit Cox model with restricted cubic spline
# ----------------
m1 <- cph(
  Surv(dif_treat_outcome, death_vl) ~
    rcs(dif_day_tx_symp, 4) +
    age_group_2 + gender + region +
    education_level + race_2,
  data = data_mort,
  x = TRUE,
  y = TRUE
)

summary(m1)

# ----------------
# Predict adjusted hazard ratios
# Reference = median delay
# ----------------
pred <- Predict(
  m1,
  dif_day_tx_symp,
  fun = exp,        # HR
  ref.zero = TRUE  # reference = median of dif_day_tx_symp
)

pred_df <- as.data.frame(pred)

# ----------------
# Plot adjusted HR curve
# ----------------
cox_hiv <- ggplot(pred_df, aes(x = dif_day_tx_symp, y = yhat)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_vline(aes(xintercept = median(data_mort$dif_day_tx_symp, na.rm = T)), 
             color = "red", linetype = "dashed") +
  labs(
    x = "Days between symptoms and treatment",
    y = paste0("Hazard ratio (reference = ", 
               median(data_mort$dif_day_tx_symp[data_mort$hiv_2 == "Yes"], na.rm = T), " days",")")
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 14, color = "black", face = "bold"),
    axis.text.x = element_text(size = 12, color = "black")
  ) +
  scale_x_continuous(n.breaks = 10) +
  scale_y_continuous(n.breaks = 12, limits = c(0.4,5))

ggsave(cox_hiv, path = here("figures and tables"),
       width = 10, height = 8, dpi = 300,
       file = "cox_model_hiv.png")


#### COX model patients without HIV
# --------------
# Prepare data
# ----------------
data_mort <- data_6 %>%
  filter(!is.na(death_vl),
         !is.na(dif_day_tx_symp)) %>% 
  dplyr::select(
    dif_day_tx_symp, death_vl, dif_treat_outcome,
    age_group_2, gender, region,
    education_level, race_2, hiv_2
  ) %>% 
  filter(hiv_2 == "No")

# ----------------
# Variable order
# ----------------
var_order_2 <- c(
  "age_group_2", "gender", "region",
  "education_level", "race_2"
)

# ----------------
# Explicit factor level ordering (reference = first)
# ----------------
level_orders <- list(
  age_group_2 = c("16-50", "0-1", "1-5", "6-15", ">50"),
  gender      = c("Female", "Male"),
  region      = c("Northeast", "Center-West", "North", "Southeast", "South"),
  education_level = c(
    "Complete elementary and incomplete middle school",
    "Incomplete elementary school",
    "Complete middle school and incomplete high school",
    "Complete high school",
    "Incomplete university",
    "Complete university",
    "N/A"
  ),
  race_2 = c("Brown", "White", "Afro-descendants", "Yellow", "Indigenous", "Unknown")
)

for (v in names(level_orders)) {
  data_mort[[v]] <- factor(data_mort[[v]], levels = level_orders[[v]])
}

# ----------------
# datadist
# ----------------
dd <- datadist(data_mort)
options(datadist = "dd")

# ----------------
# Fit Cox model with restricted cubic spline
# ----------------
m1 <- cph(
  Surv(dif_treat_outcome, death_vl) ~
    rcs(dif_day_tx_symp, 4) +
    age_group_2 + gender + region +
    education_level + race_2,
  data = data_mort,
  x = TRUE,
  y = TRUE
)

# ----------------
# Predict adjusted hazard ratios
# Reference = median delay
# ----------------
pred <- Predict(
  m1,
  dif_day_tx_symp,
  fun = exp,        # HR
  ref.zero = TRUE  # reference = median of dif_day_tx_symp
)

pred_df <- as.data.frame(pred)

# ----------------
# Plot adjusted HR curve
# ----------------
cox_nohiv <- ggplot(pred_df, aes(x = dif_day_tx_symp, y = yhat)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_vline(aes(xintercept = median(data_mort$dif_day_tx_symp, na.rm = T)), 
             color = "red", linetype = "dashed") +
  labs(
    x = "Days between symptoms and treatment",
    y = paste0("Hazard ratio (reference = ", 
               median(data_mort$dif_day_tx_symp[data_mort$hiv_2 == "No"], na.rm = T), " days",")")
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 14, color = "black", face = "bold"),
    axis.text.x = element_text(size = 12, color = "black")
  ) +
  scale_x_continuous(n.breaks = 10) +
  scale_y_continuous(n.breaks = 12, limits = c(0.4, 5))

ggsave(cox_nohiv, path = here("figures and tables"),
       width = 10, height = 8, dpi = 300,
       file = "cox_model_nohiv.png")

ggarrange(cox_nohiv, 
          cox_hiv, labels = c("A. No HIV", "B. HIV"), 
          label.x = (0.2), label.y = (0.99)) 
ggsave(path = here("figures and tables"),
       width = 10, height = 8, dpi = 300,
       file = "cox_model.png")


summary(m1)


############################# Baseline characteristics by delay between symptoms onset and treatment initiation ########################################################

### Distribution and direction time between recurrences by characteristics and variables selected 
### variables in the table
var_order <- c("age_group_2", "gender", "year_period", "region",
               "education_level", "race_2", "hiv_2", "evolution")

### desired level orders WITH the reference (or desired) first
level_orders <- list(
  age_group_2 = c("0-1", "1-5", "6-15", "16-50", ">50"),
  gender      = c("Female", "Male"),
  year_period = c("2007-2009", "2010-2014", "2015-2018", "2019-2023"),
  region      = c("Center-West", "North", "Northeast", "Southeast", "South"),
  education_level = c(
    "Incomplete elementary school",
    "Complete elementary and incomplete middle school",
    "Complete middle school and incomplete high school",
    "Complete high school",
    "Incomplete university",
    "Complete university",
    "N/A"
  ),
  race_2 = c("Brown", "White", "Afro-descendants", "Yellow", "Indigenous", "Unknown"),
  hiv_2  = c("No", "Yes", "Unknown"),
  evolution = c("Cure", "Death for VL", "Death for any-cause", "Neglect", "Transferred", NA)
)

### Function to make baseline table by delay variables
median_time_order <- function(vars, level_orders = list(), data = data_6) {
  map_dfr(vars, function(v) {
    
    # subset: valid characteristic and non-missing ts-tt delay > 0
    d_sub <- data %>%
      filter(!is.na(.data[[v]]),
             !is.na(dif_day_tx_symp),
             dif_day_tx_symp > 0)
    
    if (nrow(d_sub) == 0) return(NULL)
    
    # enforce factor levels according to level_orders (this sets display/reference)
    if (!is.null(level_orders[[v]])) {
      d_sub[[v]] <- factor(d_sub[[v]], levels = level_orders[[v]])
    } else {
      d_sub[[v]] <- factor(d_sub[[v]])
    }
    
    # denominator for n/N (here: all rows in data_6, as you had it)
    denom <- nrow(d_sub)
    
    # summary: n/N, medians [IQR] and ranges for 3 time intervals
    df_sum <- d_sub %>%
      group_by(.data[[v]]) %>%
      summarise(
        n = n(),
        
        median_1 = paste0(
          round(median(dif_day_dx_symp, na.rm = TRUE), 0), " [",
          round(quantile(dif_day_dx_symp, 0.25, na.rm = TRUE), 0), "-",
          round(quantile(dif_day_dx_symp, 0.75, na.rm = TRUE), 0), "]"
        ),
        range_1 = paste0(
          min(dif_day_dx_symp, na.rm = TRUE), "-",
          max(dif_day_dx_symp, na.rm = TRUE)
        ),
        
        median_2 = paste0(
          round(median(dif_day_tx_dx, na.rm = TRUE), 0), " [",
          round(quantile(dif_day_tx_dx, 0.25, na.rm = TRUE), 0), "-",
          round(quantile(dif_day_tx_dx, 0.75, na.rm = TRUE), 0), "]"
        ),
        range_2 = paste0(
          min(dif_day_tx_dx, na.rm = TRUE), "-",
          max(dif_day_tx_dx, na.rm = TRUE)
        ),
        
        median_3 = paste0(
          round(median(dif_day_tx_symp, na.rm = TRUE), 0), " [",
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
        p = paste0(scales::percent(n/denom, accuracy = 0.01, suffix = "%"), 
                   "(", 
                   scales::comma(n,    big.mark = ","), "/", 
                   scales::comma(denom, big.mark = ","),
                   ")"
        )
      ) %>%
      dplyr::select(variable, level, p,
             median_1, range_1,
             median_2, range_2,
             median_3, range_3)
    
    # apply desired level order to the summarised table
    if (!is.null(level_orders[[v]])) {
      df_sum <- df_sum %>% mutate(level = factor(level, levels = level_orders[[v]]))
    } else {
      df_sum <- df_sum %>% mutate(level = factor(level))
    }
    
    df_sum
  }) %>%
    mutate(variable = factor(variable, levels = var_order)) %>%
    arrange(variable, level)
}

# run
vars <- var_order
dif_table <- median_time_order(vars, level_orders, data = data_6)

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
    "Characteristic"              = variable,
    "Level"                       = level,
    "% (n/N)"              = p,
    "Median time ts-td (days)"    = median_1,
    "Range ts-td (days)"          = range_1,
    "Median time td-tt (days)"    = median_2,
    "Range td-tt (days)"          = range_2,
    "Median time ts-tt (days)"    = median_3,
    "Range ts-tt (days)"          = range_3
  ) %>% 
  flextable() %>% 
  # group rows (you may need to adjust indices if your data changes)
  merge_at(j = 1, i = 1:5,  part = "body") %>%  
  merge_at(j = 1, i = 6:7,  part = "body") %>% 
  merge_at(j = 1, i = 8:11, part = "body") %>% 
  merge_at(j = 1, i = 12:16,part = "body") %>% 
  merge_at(j = 1, i = 17:23,part = "body") %>% 
  merge_at(j = 1, i = 24:29,part = "body") %>% 
  merge_at(j = 1, i = 30:32,part = "body") %>% 
  merge_at(j = 1, i = 33:37,part = "body") %>% 
  # merge_at(j = 1, i = 38:42,part = "body") %>% 
  hline(i = c(5, 7, 11, 16, 23, 29, 32), part = "body", border = border_line) %>% 
  # center numeric / formatted columns
  align(j = 3:9, part = "body",   align = "center") %>% 
  align(j = 3:9, part = "header", align = "center") %>% 
  # bold characteristic column & header
  bold(j = 1, part = "all") %>% 
  bold(i = 1, part = "header") %>% 
  autofit() 

tab_1

# save
tab_1 %>% 
  save_as_docx(path = here("figures and tables", "table_1.docx"))

#################################### Severity adjusted model #######################################################




