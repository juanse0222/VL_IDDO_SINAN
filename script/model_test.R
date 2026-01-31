#Analysis model using delay category 
install.packages(c("survival","rms","coxme","mice","Hmisc", "sjPlot", "splines"))
library(survival)
library(rms)    # for rcs / cph
library(coxme)  # mixed effects cox
library(mice)
library(Hmisc)
library(survminer)
library(sjPlot)
library(splines)

### Import work database final for running up the model
data_6 <- import(here("data base and clean", "data_6.rds")) 

### Adjusting variables 
data_6 <- data_6 %>% 
  mutate(severe = if_else(hemorrhagic_phenomena == 1 | edema == 1, "Yes", "No"))

### Adjusting variables 
data_6 <- data_6 %>%
  mutate(age_group_2 = case_when(age_years >= 0 & age_years <= 1 ~ "0-1",
                                 age_years > 1 & age_years <= 5 ~ "1-5",
                                 age_years > 5 & age_years <= 15 ~ "6-15",
                                 age_years > 15 & age_years <= 50 ~ "16-50",
                                 age_years > 50  ~ ">50",
                                 TRUE ~ NA),
         age_group_3 = case_when(age_years >= 0 & age_years <= 1 ~ "0-1",
                                 age_years >= 1 & age_years <= 2 ~ "1-2",
                                 age_years >= 2 & age_years <= 3 ~ "2-3",
                                 age_years >= 3 & age_years <= 5 ~ "3-5",
                                 age_years >= 5 & age_years <= 15 ~ "5-15",
                                 age_years >= 15 & age_years <= 20 ~ "15-20",
                                 age_years >= 20 & age_years <= 40 ~ "20-40",
                                 age_years >= 40 & age_years <= 60 ~ "40-60",
                                 age_years > 60  ~ ">60",
                                 TRUE ~ NA), 
         year = lubridate::year(clean_date_of_notification),
         year_period = case_when(str_detect(year, "2007|2008|2009") ~ "2007-2009",
                                 str_detect(year, "2010|2011|2012|2013|2014") ~ "2010-2014",
                                 str_detect(year, "2015|2016|2017|2018") ~ "2015-2018",
                                 str_detect(year, "2019|2020|2021|2022|2023") ~ "2019-2023",
                                 TRUE ~ NA),
         race_2 = case_when(race == 1 ~ "White",
                            race == 2 ~ "Afro-descendants",
                            race == 3 ~ "Yellow",
                            race == 4 ~ "Brown",
                            race == 5 ~ "Indigenous",
                            race == 9 ~ "Unknown",
                            TRUE ~ NA),
         region = case_when(
           federative_unit %in% c(11, 13, 14, 15, 16, 17) ~ "North",
           federative_unit %in% c(21, 22, 23, 24, 25, 26, 27, 28, 29) ~ "Northeast",
           federative_unit %in% c(31, 32, 33, 35)         ~ "Southeast",
           federative_unit %in% c(41, 42, 43)             ~ "South",
           federative_unit %in% c(50, 51, 52, 53)         ~ "Center-West",
           TRUE                                           ~ NA_character_),
         age_group_3 = factor(age_group_3, levels = c("0-1", "1-2", "2-3", "3-5", "5-15",
                                                      "15-20", "20-40", "40-60",
                                                      ">60")),
         education_level = factor(education_level, levels = c("Incomplete elementary school",
                                                              "Complete elementary and incomplete middle school",
                                                              "Complete middle school and incomplete high school",
                                                              "Complete high school",
                                                              "Incomplete university",
                                                              "Complete university",
                                                              "N/A")),
         age_group_2 = factor(age_group_2, levels = c("0-1", "1-5", "6-15", "16-50", ">50")),
         fever_2 = case_when(fever == 1 ~ "Yes",
                             fever == 2 ~ "No",
                             fever == 9 ~ "Unknown",
                             TRUE ~ NA),
         edema_2 = case_when(edema == 1 ~ "Yes",
                             edema == 2 ~ "No",
                             edema == 9 ~ "Unknown",
                             TRUE ~ NA),
         hemorrhagic_phenomena_2 = case_when(hemorrhagic_phenomena == 1 ~ "Yes",
                                             hemorrhagic_phenomena == 2 ~ "No",
                                             hemorrhagic_phenomena == 9 ~ "Unknown",
                                             TRUE ~ NA),
         hiv_2 = case_when(hiv == 1 ~ "Yes",
                           hiv == 2 ~ "No",
                           hiv == 9 ~ "Unknown",
                           TRUE ~ NA),
         fever_2 = factor(fever_2, levels = c("No", "Yes", "Unknown")),
         edema_2 = factor(edema_2, levels = c("No", "Yes", "Unknown")),
         hemorrhagic_phenomena_2 = factor(hemorrhagic_phenomena_2, levels = c("No", "Yes", "Unknown")),
         name_tx = case_when(
           type_of_rx == 1 ~ "Antimonial Pentavalent",
           type_of_rx == 2 ~ "Amphotericin B deoxycholate",
           type_of_rx == 3 ~ "Others",                     ### Pentamidine - Others
           type_of_rx == 4 ~ "Liposomal amphotericin B",
           type_of_rx == 5 ~ "Others",
           type_of_rx == 6 ~ NA_character_,                # Nothing - NA (excluded)??
           TRUE ~ NA_character_)) 
# mutate_at(vars(dif_day_dx_symp, dif_day_tx_dx), ~ if_else(. < 0, NA, .))

data_ost <- data_6 %>% 
  mutate(year_3 = case_when(year >= 2007 & year <= 2009 ~ "2007 - 2009",
                            year >= 2010 & year <= 2014 ~ "2010 - 2014",
                            year >= 2015 & year <= 2018 ~ "2015 - 2018",
                            year >= 2019 & year <= 2023 ~ "2019 - 2023",
                            TRUE ~ NA)) %>% 
  group_by(federative_unit, year_3) %>% 
  summarise(median = median(dif_day_tx_symp), 
            n = n(), 
            mortality = sum(death_vl),
            lethality = round(mortality/n*100, digits = 2),
            median_age = median(age_years), .groups = "drop")

data_ost_2 <- data_ost %>% 
  filter(n >= 100)

state_100 <- unique(data_ost$federative_unit[data_ost$n < 100])

######################################################################################################################

### Adjusting and creating new variables to include in the models 

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

data_m1 <- data_8 %>% 
  filter(!is.na(death_vl)) %>% 
  #filter(!hiv_2 == "Yes") %>% 
  select(dif_day_tx_symp, dif_treat_start_symp_cat_2, dif_treat_outcome, age_group_2,
         age_years, gender_2, hiv_3, hiv_2, gender, death_vl, name_tx, severe_2, race_2,
         federative_unit, education_level, region, race_2) 
#filter(!federative_unit == state_100)

export(data_8, here("data base and clean", "data_8.rds")) ### Database to work with time between treatment initiation and death by VL

###########################################################################################################################################

### Model overall - including HIV patients and unadjusted

data_m1 <- data_8 %>% 
  filter(!is.na(death_vl)) %>% 
  filter(!hiv_2 == "Yes") %>% 
  select(dif_day_tx_symp, dif_treat_start_symp_cat_2, dif_treat_outcome, age_group_2,
         age_years, gender_2, hiv_3, hiv_2, gender, death_vl, name_tx, severe_2, race_2,
         federative_unit, education_level, region, race_2) 

dd <- datadist(data_m1)
options(datadist = "dd")  ### rcs

m1 <- cph(Surv(dif_treat_outcome, death_vl) ~ 
              rcs(dif_day_tx_symp, 4), 
            data = data_m1
)

summary(m1)

m1_p <- Predict(
  m1, 
  dif_day_tx_symp,
  fun = exp,        
  ref.zero = F   # reference = median
)

m1_p_df <- as.data.frame(m1_p)

m1 <- ggplot(m1_p_df, aes(dif_day_tx_symp, yhat)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_hline(aes(yintercept = 1), linetype = "dashed") +
  #geom_hline(yintercept = 1, linetype = "dashed") +
  labs(
    x = "Delay to treatment initiation (days)",
    y = "Hazard ratio (Death by VL)"
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 14, color = "black", face = "bold"),
    axis.text = element_text(size = 12, color = "black"),
  ) +
  scale_x_continuous(n.breaks = 8) +
  scale_y_continuous(n.breaks = 10)

m1


### 1. Cox regression model  - delay and mortality without HIV

data_m1 <- data_8 %>% 
  filter(!is.na(death_vl)) %>% 
  #filter(!hiv_2 == "Yes") %>% 
  select(dif_day_tx_symp, dif_treat_start_symp_cat_2, dif_treat_outcome, age_group_2,
         age_years, gender_2, hiv_3, hiv_2, gender, death_vl, name_tx, severe_2, race_2,
         federative_unit, education_level, region, race_2) 

dd <- datadist(data_m1)
options(datadist = "dd")  ### rcs

m1 <- coxph(Surv(dif_treat_outcome, death_vl) ~ 
              rcs(dif_day_tx_symp, 4), 
            data = data_m1
)

summary(m1)
termplot(m1, term=1, se=TRUE, col.term=1, col.se=1)

newd <- data.frame(
  dif_day_tx_symp = seq(
    min(data_m1$dif_day_tx_symp, na.rm = TRUE),
    max(data_m1$dif_day_tx_symp, na.rm = TRUE),
    length.out = 200
  )
)

pred <- predict(m1, newdata = newd, type = "lp", se.fit = TRUE)

newd$hr <- exp(pred$fit)
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

###### Model with HIV and delay as a predictor. 

m1_hiv <- coxph(Surv(dif_treat_outcome, death_vl) ~ 
              rcs(dif_day_tx_symp, 4), 
            data = data_8
)

summary(m1_hiv)

newd_hiv <- data.frame(
  dif_day_tx_symp = seq(
    min(data_8$dif_day_tx_symp, na.rm = TRUE),
    max(data_8$dif_day_tx_symp, na.rm = TRUE),
    length.out = 200
  )
)

pred <- predict(m1_hiv, newdata = newd, type = "lp", se.fit = TRUE)

newd_hiv$hr <- exp(pred$fit)
newd_hiv$lower <- exp(pred$fit - 1.96 * pred$se.fit)
newd_hiv$upper <- exp(pred$fit + 1.96 * pred$se.fit)

m1_hiv <- ggplot(newd_hiv, aes(dif_day_tx_symp, hr)) +
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

hiv_model_panel <- ggarrange(m1_whiv, m1_hiv, labels = c("a)", "b)"), 
                             hjust = c(-0.7, -0.7), font.label =  list(family = "Times New Roman") )

hiv_model_panel

combine_model_1 <- ggplot(newd, aes(dif_day_tx_symp, hr, color = "No HIV")) +
  geom_line(linewidth = 0.5) +
  geom_ribbon(aes(ymin = lower, ymax = upper, color = "No HIV"),  
              fill = "darkgreen", alpha = 0.2) +
  geom_line(data =newd_hiv, aes(dif_day_tx_symp, hr, color = "HIV"), 
            linewidth = 0.5, linetype = "dashed") +
  geom_ribbon(data = newd_hiv, aes(ymin = lower, ymax = upper, color = "HIV"), 
              fill = "brown2", alpha = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(
    x = "Delay to treatment initiation (days)",
    y = "Hazard ratio (Death by VL)"
  ) +
  labs(color = "Condition") +
  scale_x_continuous(n.breaks = 8) +
  scale_y_continuous(n.breaks = 10, limits = c(0.75, 2.8)) +
  scale_color_manual(values = c("HIV" = "brown2", "No HIV" = "darkgreen")) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, color = "black", face = "bold"),
    legend.title = element_text(size = 14, color = "black", face = "bold"),
    legend.text = element_text(size = 12, color = "black")
  ) 

combine_model_1
############################################################################################################

###############################################################################
## 1) Variables to adjust for and factor level ordering
###############################################################################

var_order_2 <- c(
  "age_group_2", "gender", "region",
  "education_level", "race_2"
  #"hiv_2"
)

level_orders <- list(
  age_group_2 = c("0-1", "1-5", "6-15", "16-50", ">50"),
  gender      = c("Female", "Male"),
  region      = c("North", "Center-West", "Northeast", "Southeast", "South"),
  education_level = c(
    "Complete elementary and incomplete middle school",
    "Incomplete elementary school",
    "Complete middle school and incomplete high school",
    "Complete high school",
    "Incomplete university",
    "Complete university",
    "N/A"
  ),
  #hiv_2 = c("Yes", "No"),
  race_2 = c("Brown", "White", "Afro-descendants", "Yellow", "Indigenous", "Unknown")
)

###############################################################################
## 2) Modelling dataset
###############################################################################

data_mort <- data_m1

###############################################################################
## 3) Enforce factor level orders (reference = first level)
###############################################################################

for (v in var_order_2) {
  if (!is.null(level_orders[[v]])) {
    data_mort[[v]] <- factor(data_mort[[v]], levels = level_orders[[v]])
  } else {
    data_mort[[v]] <- droplevels(as.factor(data_mort[[v]]))
  }
}

###############################################################################
## 4) Cox proportional hazards model with spline on treatment delay
###############################################################################

m2 <- coxph(
  Surv(dif_treat_outcome, death_vl) ~
    rcs(dif_day_tx_symp, k = 4) +
    age_group_2 + gender + race_2 + region +
    education_level,
  data = data_8
)

str(m2)

summary(m2)

###############################################################################
## 5) Prediction grid over treatment delay (reference profile fixed)
###############################################################################

delay_seq <- seq(
  quantile(data_mort$dif_day_tx_symp, 0.01, na.rm = TRUE),
  quantile(data_mort$dif_day_tx_symp, 0.99, na.rm = TRUE),
  length.out = 100
)

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

###############################################################################
## 6) Adjusted hazard ratios (reference = shortest delay)
###############################################################################

lp <- predict(m2, newdata = newdat, type = "lp", se.fit = TRUE)

ref_lp <- predict(
  m2,
  newdata = newdat[1, , drop = FALSE],
  type = "lp"
)

newdat <- newdat %>%
  mutate(
    hr    = exp(lp$fit - ref_lp),
    lower = exp(lp$fit - ref_lp - 1.96 * lp$se.fit),
    upper = exp(lp$fit - ref_lp + 1.96 * lp$se.fit)
  )

###############################################################################
## 7) Publication-ready plot: hazard ratio vs treatment delay
###############################################################################

ggplot(newdat, aes(x = dif_day_tx_symp, y = hr)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  #geom_rug(
   # data = data_mort,
  #  aes(x = dif_day_tx_symp),
  #  inherit.aes = FALSE,
  #  sides = "b",
  #  alpha = 0.3
  #) +
  labs(
    x = "Delay between symptoms onset and treatment initiation (days)",
    y = "Hazard ratio (Death by VL)",
    title = "Adjusted Hazard Ratio (HR) of death vs Delay symptoms onset and treatmet initiation",
  ) +
  scale_x_continuous(limits = c(0, 250), n.breaks = 15) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 14, color = "black", face = "bold"),
    axis.text = element_text(size = 12, color = "black"),
    plot.title = element_text(hjust = 0.5, color = "black", size = 16, face = "bold")
  )
  
  
####################################################################################################################################

### Logistic regression model based on symptoms
data_symp_no_hiv <- data_8 %>% 
  mutate(jaundice_2 = case_when(jaundice == 1 ~ 1,
                                jaundice == 2 ~ 0,
                                TRUE ~ NA),
         bleeding = case_when(hemorrhagic_phenomena == 1 ~ 1,
                                             hemorrhagic_phenomena == 2 ~ 0,
                                             TRUE ~ NA),
         edema_3 = case_when(edema == 1 ~ 1,
                             edema == 2 ~ 0,
                             TRUE ~NA)) %>% 
  filter(!hiv_2 == "Yes")


### Other way to approach of the model using Predic instead of predict 

### Bleeding model option 2

data_bleeding <- data_symp_no_hiv %>% 
  select(bleeding, dif_day_tx_symp) %>% 
  filter(!is.na(bleeding))

dd <- datadist(data_bleeding)
options(datadist = "dd")

bleeding_fit <- lrm(bleeding ~ rcs(dif_day_tx_symp, 4), data = data_bleeding)

p_bleeding <- Predict(
  bleeding_fit,
  dif_day_tx_symp,
  fun = exp,        # odds ratios
  ref.zero = TRUE   # reference = median
)

p_bleeding_df <- as.data.frame(p_bleeding)

bleeding_m1 <- ggplot(p_bleeding_df, aes(dif_day_tx_symp, yhat)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  #geom_hline(yintercept = 1, linetype = "dashed") +
  labs(
    x = "Delay to treatment initiation (days)",
    y = "Odds (Bleeding)"
  ) +
  theme_bw()
  #scale_x_continuous(n.breaks = 8) +
  #scale_y_continuous(n.breaks = 10, limits = c(0.75, 2.8))

bleeding_m1

### Jaundice
data_jaundice <- data_symp_no_hiv %>% 
  select(jaundice_2, dif_day_tx_symp) %>% 
  filter(!is.na(jaundice_2))

dd <- datadist(data_jaundice)
options(datadist = "dd")

jaundice_fit <- lrm(jaundice_2 ~ rcs(dif_day_tx_symp, 4), data = data_jaundice)

p_jaundice <- Predict(
  jaundice_fit,
  dif_day_tx_symp,
  fun = exp,        # odds ratios
  ref.zero = TRUE   # reference = median
)

p_jaundice_df <- as.data.frame(p_jaundice)

jaundice_m1 <- ggplot(p_jaundice_df, aes(dif_day_tx_symp, yhat)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  #geom_hline(yintercept = 1, linetype = "dashed") +
  labs(
    x = "Delay to treatment initiation (days)",
    y = "Odds (Jaundice)"
  ) +
  theme_bw()
#scale_x_continuous(n.breaks = 8) +
#scale_y_continuous(n.breaks = 10, limits = c(0.75, 2.8))

jaundice_m1

### oedema

data_oedema <- data_symp_no_hiv %>%
  select(dif_day_tx_symp, edema_3) %>% 
  filter(!is.na(edema_3))

#oedema_m <- glm(formula = edema_3 ~ rcs(dif_day_tx_symp, 4), family = "binomial",
 #                 data = data_oedema)

edema_fit <- lrm(edema_3 ~ rcs(dif_day_tx_symp, 4), data = data_oedema)

p_edema <- Predict(
  edema_fit,
  dif_day_tx_symp,
  fun = exp,        # odds ratios
  ref.zero = TRUE   # reference = median
)

p_edema_df <- as.data.frame(p_edema)

edema <- ggplot(p_edema_df, aes(dif_day_tx_symp, yhat)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  #geom_hline(yintercept = 1, linetype = "dashed") +
  labs(
    x = "Delay to treatment initiation (days)",
    y = "Odds (Oedema)"
  ) +
  theme_bw()
#scale_x_continuous(n.breaks = 8) +
#scale_y_continuous(n.breaks = 10, limits = c(0.75, 2.8))

edema

symptoms_m1 <- ggplot(
  data = p_bleeding_df,
  aes(dif_day_tx_symp, yhat)
) +
  geom_line(aes(color = "Bleeding"), linewidth = 1) +
  geom_ribbon(
    aes(ymin = lower, ymax = upper),
    alpha = 0.2, fill = "brown2"
  ) +
  
  geom_line(
    data = p_jaundice_df,
    aes(dif_day_tx_symp, yhat, color = "Jaundice"),
    linewidth = 1,
  ) +
  geom_ribbon(
    data = p_jaundice_df,
    aes(
      x = dif_day_tx_symp,
      ymin = lower,
      ymax = upper,
    ), fill = "orange2",
    alpha = 0.2,
    inherit.aes = FALSE
  ) +
  
  geom_line(
    data = p_edema_df,
    aes(dif_day_tx_symp, yhat, color = "Oedema"),
    linewidth = 1,
  ) +
  geom_ribbon(
    data = p_edema_df,
    aes(
      x = dif_day_tx_symp,
      ymin = lower,
      ymax = upper
    ), fill = "darkgreen",
    alpha = 0.2,
    inherit.aes = FALSE
  ) +
  
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(
    x = "Delay to treatment initiation (days)",
    y = "Odds Ratio", 
    color = "Symptoms"
  ) +
  scale_y_continuous(limits = c(0, 5)) +
  theme_bw() +
  scale_color_manual(values = c("Jaundice" = "orange2", 
                                "Bleeding" = "brown2", 
                                "Oedema" = "darkgreen")) 
  
symptoms_m1
  

data_symp_hiv <- data_8 %>% 
  mutate(jaundice_2 = case_when(jaundice == 1 ~ 1,
                                jaundice == 2 ~ 0,
                                TRUE ~ NA),
         bleeding = case_when(hemorrhagic_phenomena == 1 ~ 1,
                              hemorrhagic_phenomena == 2 ~ 0,
                              TRUE ~ NA),
         edema_3 = case_when(edema == 1 ~ 1,
                             edema == 2 ~ 0,
                             TRUE ~NA)) %>% 
  filter(hiv_2 == "Yes")


### Other way to approach of the model using Predic instead of predict 

### Bleeding model option 2

data_bleeding <- data_symp_hiv %>% 
  select(bleeding, dif_day_tx_symp) %>% 
  filter(!is.na(bleeding))

dd <- datadist(data_bleeding)
options(datadist = "dd")

bleeding_fit <- lrm(bleeding ~ rcs(dif_day_tx_symp, 4), data = data_bleeding)

p_bleeding <- Predict(
  bleeding_fit,
  dif_day_tx_symp,
  fun = exp,        # odds ratios
  ref.zero = TRUE   # reference = median
)

p_bleeding_df <- as.data.frame(p_bleeding)

bleeding_m1 <- ggplot(p_bleeding_df, aes(dif_day_tx_symp, yhat)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  #geom_hline(yintercept = 1, linetype = "dashed") +
  labs(
    x = "Delay to treatment initiation (days)",
    y = "Odds (Bleeding)"
  ) +
  theme_bw()
#scale_x_continuous(n.breaks = 8) +
#scale_y_continuous(n.breaks = 10, limits = c(0.75, 2.8))

bleeding_m1

### Jaundice
data_jaundice <- data_symp_hiv %>% 
  select(jaundice_2, dif_day_tx_symp) %>% 
  filter(!is.na(jaundice_2))

dd <- datadist(data_jaundice)
options(datadist = "dd")

jaundice_fit <- lrm(jaundice_2 ~ rcs(dif_day_tx_symp, 4), data = data_jaundice)

p_jaundice <- Predict(
  jaundice_fit,
  dif_day_tx_symp,
  fun = exp,        # odds ratios
  ref.zero = TRUE   # reference = median
)

p_jaundice_df <- as.data.frame(p_jaundice)

jaundice_m1 <- ggplot(p_jaundice_df, aes(dif_day_tx_symp, yhat)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  #geom_hline(yintercept = 1, linetype = "dashed") +
  labs(
    x = "Delay to treatment initiation (days)",
    y = "Odds (Jaundice)"
  ) +
  theme_bw()
#scale_x_continuous(n.breaks = 8) +
#scale_y_continuous(n.breaks = 10, limits = c(0.75, 2.8))

jaundice_m1

### oedema

data_oedema <- data_symp_hiv %>%
  select(dif_day_tx_symp, edema_3) %>% 
  filter(!is.na(edema_3))

#oedema_m <- glm(formula = edema_3 ~ rcs(dif_day_tx_symp, 4), family = "binomial",
#                 data = data_oedema)

edema_fit <- lrm(edema_3 ~ rcs(dif_day_tx_symp, 4), data = data_oedema)

p_edema <- Predict(
  edema_fit,
  dif_day_tx_symp,
  fun = exp,        # odds ratios
  ref.zero = TRUE   # reference = median
)

p_edema_df <- as.data.frame(p_edema)

edema <- ggplot(p_edema_df, aes(dif_day_tx_symp, yhat)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  #geom_hline(yintercept = 1, linetype = "dashed") +
  labs(
    x = "Delay to treatment initiation (days)",
    y = "Odds (Oedema)"
  ) +
  theme_bw()
#scale_x_continuous(n.breaks = 8) +
#scale_y_continuous(n.breaks = 10, limits = c(0.75, 2.8))

edema

symptoms_m2 <- ggplot(
  data = p_bleeding_df,
  aes(dif_day_tx_symp, yhat)
) +
  geom_line(aes(color = "Bleeding"), linewidth = 1) +
  geom_ribbon(
    aes(ymin = lower, ymax = upper),
    alpha = 0.2, fill = "brown2"
  ) +
  
  geom_line(
    data = p_jaundice_df,
    aes(dif_day_tx_symp, yhat, color = "Jaundice"),
    linewidth = 1,
  ) +
  geom_ribbon(
    data = p_jaundice_df,
    aes(
      x = dif_day_tx_symp,
      ymin = lower,
      ymax = upper,
    ), fill = "orange2",
    alpha = 0.2,
    inherit.aes = FALSE
  ) +
  
  geom_line(
    data = p_edema_df,
    aes(dif_day_tx_symp, yhat, color = "Oedema"),
    linewidth = 1,
  ) +
  geom_ribbon(
    data = p_edema_df,
    aes(
      x = dif_day_tx_symp,
      ymin = lower,
      ymax = upper
    ), fill = "darkgreen",
    alpha = 0.2,
    inherit.aes = FALSE
  ) +
  
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(
    x = "Delay to treatment initiation (days)",
    y = "Odds Ratio", 
    color = "Symptoms"
  ) +
  theme_bw() +
  scale_y_continuous(limits = c(0, 5)) +
  scale_color_manual(values = c("Jaundice" = "orange2", 
                                "Bleeding" = "brown2", 
                                "Oedema" = "darkgreen")) 

symptoms_m2

data_symp_overall <- data_8 %>% 
  mutate(jaundice_2 = case_when(jaundice == 1 ~ 1,
                                jaundice == 2 ~ 0,
                                TRUE ~ NA),
         bleeding = case_when(hemorrhagic_phenomena == 1 ~ 1,
                              hemorrhagic_phenomena == 2 ~ 0,
                              TRUE ~ NA),
         edema_3 = case_when(edema == 1 ~ 1,
                             edema == 2 ~ 0,
                             TRUE ~NA)) 


### Other way to approach of the model using Predic instead of predict 

### Bleeding model option 2

data_bleeding <- data_symp_overall %>% 
  select(bleeding, dif_day_tx_symp) %>% 
  filter(!is.na(bleeding))

dd <- datadist(data_bleeding)
options(datadist = "dd")

bleeding_fit <- lrm(bleeding ~ rcs(dif_day_tx_symp, 4), data = data_bleeding)

p_bleeding <- Predict(
  bleeding_fit,
  dif_day_tx_symp,
  fun = exp,        # odds ratios
  ref.zero = TRUE   # reference = median
)

p_bleeding_df <- as.data.frame(p_bleeding)

bleeding_m1 <- ggplot(p_bleeding_df, aes(dif_day_tx_symp, yhat)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  #geom_hline(yintercept = 1, linetype = "dashed") +
  labs(
    x = "Delay to treatment initiation (days)",
    y = "Odds (Bleeding)",
  ) +
  theme_bw()
#scale_x_continuous(n.breaks = 8) +
#scale_y_continuous(n.breaks = 10, limits = c(0.75, 2.8))

bleeding_m1

### Jaundice
data_jaundice <- data_symp_overall %>% 
  select(jaundice_2, dif_day_tx_symp) %>% 
  filter(!is.na(jaundice_2))

dd <- datadist(data_jaundice)
options(datadist = "dd")

jaundice_fit <- lrm(jaundice_2 ~ rcs(dif_day_tx_symp, 4), data = data_jaundice)

p_jaundice <- Predict(
  jaundice_fit,
  dif_day_tx_symp,
  fun = exp,        # odds ratios
  ref.zero = TRUE   # reference = median
)

p_jaundice_df <- as.data.frame(p_jaundice)

jaundice_m1 <- ggplot(p_jaundice_df, aes(dif_day_tx_symp, yhat)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  #geom_hline(yintercept = 1, linetype = "dashed") +
  labs(
    x = "Delay to treatment initiation (days)",
    y = "Odds (Jaundice)"
  ) +
  theme_bw()
#scale_x_continuous(n.breaks = 8) +
#scale_y_continuous(n.breaks = 10, limits = c(0.75, 2.8))

jaundice_m1

### oedema

data_oedema <- data_symp_overall %>%
  select(dif_day_tx_symp, edema_3) %>% 
  filter(!is.na(edema_3))

#oedema_m <- glm(formula = edema_3 ~ rcs(dif_day_tx_symp, 4), family = "binomial",
#                 data = data_oedema)

edema_fit <- lrm(edema_3 ~ rcs(dif_day_tx_symp, 4), data = data_oedema)

p_edema <- Predict(
  edema_fit,
  dif_day_tx_symp,
  fun = exp,        # odds ratios
  ref.zero = TRUE   # reference = median
)

p_edema_df <- as.data.frame(p_edema)

edema <- ggplot(p_edema_df, aes(dif_day_tx_symp, yhat)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  #geom_hline(yintercept = 1, linetype = "dashed") +
  labs(
    x = "Delay to treatment initiation (days)",
    y = "Odds (Oedema)"
  ) +
  theme_bw()
#scale_x_continuous(n.breaks = 8) +
#scale_y_continuous(n.breaks = 10, limits = c(0.75, 2.8))

edema

symptoms_m3 <- ggplot(
  data = p_bleeding_df,
  aes(dif_day_tx_symp, yhat)
) +
  geom_line(aes(color = "Bleeding"), linewidth = 1) +
  geom_ribbon(
    aes(ymin = lower, ymax = upper),
    alpha = 0.2, fill = "brown2"
  ) +
  
  geom_line(
    data = p_jaundice_df,
    aes(dif_day_tx_symp, yhat, color = "Jaundice"),
    linewidth = 1,
  ) +
  geom_ribbon(
    data = p_jaundice_df,
    aes(
      x = dif_day_tx_symp,
      ymin = lower,
      ymax = upper,
    ), fill = "orange2",
    alpha = 0.2,
    inherit.aes = FALSE
  ) +
  
  geom_line(
    data = p_edema_df,
    aes(dif_day_tx_symp, yhat, color = "Oedema"),
    linewidth = 1,
  ) +
  geom_ribbon(
    data = p_edema_df,
    aes(
      x = dif_day_tx_symp,
      ymin = lower,
      ymax = upper
    ), fill = "darkgreen",
    alpha = 0.2,
    inherit.aes = FALSE
  ) +
  
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(
    x = "Delay to treatment initiation (days)",
    y = "Odds Ratio", 
    color = "Symptoms"
  ) +
  theme_bw() +
  scale_y_continuous(limits = c(0, 5)) +
  scale_color_manual(values = c("Jaundice" = "orange2", 
                                "Bleeding" = "brown2", 
                                "Oedema" = "darkgreen")) 

symptoms_m3

final_models <- ggarrange(symptoms_m1 + rremove("xlab") + rremove("ylab") + 
                            scale_y_continuous(limits = c(0.5, 5)),
                          symptoms_m2 + rremove("xlab") + rremove("ylab") +
                            scale_y_continuous(limits = c(0.5, 5)), 
                          symptoms_m3 + rremove("xlab") + rremove("ylab") +
                            scale_y_continuous(limits = c(0.5, 5)), 
                          common.legend = T, labels = list("a)", "b)", "c)"),
                          nrow = 1, legend = "right", hjust = -0.5,
                          font.label = list(size = 14, color = "black",
                                            face = "bold", family = "Times New Roman",
                                            position = "bottom")) 

annotate_figure(final_models, 
                bottom = textGrob("Delay between symptoms onset and treatment",
                                  gp = gpar(fontsize = 14,
                                            fontface = "bold",
                                            fontfamily = "Times New Roman")), 
                #fig.lab = "Figure 1. Odds ratio of severe VL according with the delay between symptoms 
                #onset to treatment initiation \a) No HIV, b) HIV and c) Overall",
                #fig.lab.pos = "bottom",
                left = textGrob("Odds Ratio", 
                                gp = gpar(fontsize = 14,
                                          fontface = "bold",
                                          fontfamily = "Times New Roman"), 
                                rot = 90,vjust = 1), 
                )

#######################################################################################################################################

#### Model based on sever and complications taken cardinal symptoms. 

data_complication <- data_8 %>% 
  mutate(complication_2 = if_else(edema == 1|jaundice == 1|hemorrhagic_phenomena == 1, 1, 0),
         complication_3 = if_else(edema == 1 & jaundice == 1 & hemorrhagic_phenomena == 1, 1, 0),
         complication_4 = if_else(edema == 1 & hemorrhagic_phenomena == 1, 1, 0)
           )


fit_severity_model <- function(
    data,
    outcome,
    y_label = "Odds Ratio (complications)",
    knots = 4) 
  {
  
  df <- data %>%
    select(all_of(c(outcome, "dif_day_tx_symp"))) %>%
    filter(!is.na(.data[[outcome]]),
           !is.na(dif_day_tx_symp))
  
  dd <- datadist(df)
  options(datadist = "dd")
  
  fit <- lrm(
    as.formula(
      paste0(outcome, " ~ rcs(dif_day_tx_symp, ", knots, ")")),
    data = df
  )
  
  pred <- Predict( fit,dif_day_tx_symp,
    fun = exp,
    ref.zero = TRUE
  )
  
  pred_df <- as.data.frame(pred)
  
  ggplot(pred_df, aes(dif_day_tx_symp, yhat)) +
    geom_line(linewidth = 1) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
    labs(
      x = "Delay to treatment initiation (days)",
      y = y_label
    ) +
    theme_bw() +
    theme(
      axis.title = element_text(size = 14, color = "black", face = "bold"),
      axis.text.x = element_text(size = 12, color = "black")
    )
}

severity_specs <- tibble::tribble(
  ~outcome,           ~y_label,
  "complication_2",   "Any complication",
  "complication_3",   "All three complications",
  "complication_4",   "Edema + hemorrhage"
)

severity_plots <- purrr::pmap(
  severity_specs,
  fit_severity_model,
  data = data_complication
)

severity_plots[[1]]  # == severity_m1
severity_plots[[2]]  # == severity_m2
severity_plots[[3]]  # == severity_m3

names(severity_plots) <- severity_specs$y_label

final_severity <- ggpubr::ggarrange(
  plotlist = severity_plots,
  nrow = 1,
  common.legend = FALSE
  )



### Model patients with HIV

severity_plots <- purrr::pmap(
  severity_specs,
  fit_severity_model,
  data = data_complication %>% filter(hiv_2 == "Yes")
)

severity_plots[[1]]  # == severity_m1
severity_plots[[2]]  # == severity_m2
severity_plots[[3]]  # == severity_m3

names(severity_plots) <- severity_specs$y_label

final_severity_2 <- ggpubr::ggarrange(
  plotlist = severity_plots,
  nrow = 1,
  common.legend = FALSE
)

final_severity_2


### Model patients without HIV 

severity_plots <- purrr::pmap(
  severity_specs,
  fit_severity_model,
  data = data_complication %>% filter(hiv_2 == "No")
)

severity_plots[[1]]  # == severity_m1
severity_plots[[2]]  # == severity_m2
severity_plots[[3]]  # == severity_m3

names(severity_plots) <- severity_specs$y_label

final_severity_3 <- ggpubr::ggarrange(
  plotlist = severity_plots,
  nrow = 1,
  common.legend = FALSE
)

final_severity_3

panel_model_complications <- ggarrange(final_severity, 
                                       final_severity_2, 
                                       final_severity_3, 
                                       nrow = 3, labels = "AUTO"
                                       )

panel_model_complications












## --------------------------------------------------
## 2. Function: fit spline model + predict + plot
##    (HIV used ONLY as a filter)
## --------------------------------------------------

fit_severity_model <- function(
    data,
    outcome,
    y_label,
    hiv_status = NULL,
    knots = 4,
    min_events = 5
) {
  
  df <- data %>%
    select(all_of(c(outcome, "dif_day_tx_symp", "hiv_3"))) %>%
    filter(
      !is.na(.data[[outcome]]),
      !is.na(dif_day_tx_symp)
    )
  
  # Apply HIV filter
  if (!is.null(hiv_status)) {
    df <- df %>% filter(hiv_3 == hiv_status)
  }
  
  # ---- SAFETY CHECKS ----
  n_total  <- nrow(df)
  n_events <- sum(df[[outcome]] == 1, na.rm = TRUE)
  
  if (n_total < 10 || n_events < min_events) {
    message(
      paste(
        "Skipping:",
        outcome,
        ifelse(is.null(hiv_status), "(All)",
               ifelse(hiv_status == 1, "(HIV+)", "(HIV−)")),
        "- insufficient data (n =", n_total,
        ", events =", n_events, ")"
      )
    )
    return(NULL)
  }
  
  # ---- rms workflow ----
  dd <- datadist(df)
  options(datadist = "dd")
  
  fit <- lrm(
    as.formula(
      paste0(outcome, " ~ rcs(dif_day_tx_symp, ", knots, ")")
    ),
    data = df
  )
  
  pred <- Predict(
    fit,
    dif_day_tx_symp,
    fun = exp,
    ref.zero = TRUE
  )
  
  pred_df <- as.data.frame(pred)
  
  hiv_label <- if (is.null(hiv_status)) {
    "All"
  } else if (hiv_status == 1) {
    "HIV+"
  } else {
    "HIV−"
  }
  
  ggplot(pred_df, aes(dif_day_tx_symp, yhat)) +
    geom_line(linewidth = 1) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
    labs(
      x = "Delay to treatment initiation (days)",
      y = "Odds Ratio",
      title = paste(y_label, "-", hiv_label)
    ) +
    theme_bw()
}

## --------------------------------------------------
## 3. Specification table: outcomes × HIV strata
## --------------------------------------------------

severity_specs <- crossing(
  tibble::tribble(
    ~outcome,           ~y_label,
    "complication_2",   "Any complication",
    "complication_3",   "All three complications",
    "complication_4",   "Edema + hemorrhage"
  ),
  tibble::tribble(
    ~hiv_status,
    NA_integer_,   # All
    0L,            # HIV−
    1L             # HIV+
  )
)

## --------------------------------------------------
## 4. Run ALL models automatically
## --------------------------------------------------

severity_plots <- pmap(
  severity_specs,
  fit_severity_model,
  data = data_complication
)

## Name plots for easy selection
names(severity_plots) <- with(
  severity_specs,
  paste(y_label,
        ifelse(is.na(hiv_status), "All",
               ifelse(hiv_status == 1, "HIV+", "HIV−")),
        sep = " – ")
)

## --------------------------------------------------
## 5. Examples of arranging plots
## --------------------------------------------------

## Example 1: One outcome, stratified by HIV
plots_any_complication <- severity_plots[
  grepl("^Any complication", names(severity_plots))
]

ggarrange(
  plotlist = plots_any_complication,
  nrow = 1
)

## Example 2: HIV− only, all outcomes
plots_hiv_neg <- severity_plots[grepl("HIV−$", names(severity_plots))]

ggarrange(
  plotlist = plots_hiv_neg,
  nrow = 1
)

## Example 3: Full grid (3 outcomes × 3 HIV strata)
ggarrange(
  plotlist = severity_plots,
  nrow = 3,
  ncol = 3
)






  
fit_symptom_rms <- function(
    data,
    outcome,
    outcome_label,
    hiv_filter = NULL   # NULL = all, 0 = HIV-, 1 = HIV+
) {
  
  df <- data_symp %>%
    select(all_of(c(outcome, "dif_day_tx_symp", "hiv_3"))) %>%
    filter(!is.na(.data[[outcome]]),
           !is.na(dif_day_tx_symp))
  
  # Apply HIV filter if requested
  if (!is.null(hiv_filter)) {
    df <- df %>% filter(hiv_3 == hiv_filter)
  }
  
  dd <- datadist(df)
  options(datadist = "dd")
  
  fit <- lrm(
    formula = as.formula(
      paste0(outcome, " ~ rcs(dif_day_tx_symp, 4)")
    ),
    data = df
  )
  
  pred <- Predict(
    fit,
    dif_day_tx_symp,
    fun = exp,
    ref.zero = TRUE
  )
  
  as.data.frame(pred) %>%
    mutate(
      symptom = outcome_label,
      hiv_group = case_when(
        is.null(hiv_filter) ~ "All",
        hiv_filter == 0     ~ "HIV−",
        hiv_filter == 1     ~ "HIV+"
      )
    )
}

symptom_specs <- tibble::tribble(
  ~outcome,      ~outcome_label,
  "bleeding",    "Bleeding",
  "jaundice_2",  "Jaundice",
  "edema_3",     "Oedema"
)

pred_all <- purrr::pmap_dfr(
  symptom_specs,
  fit_symptom_rms,
  data = data_symp,
  hiv_filter = NULL
)

pred_hiv_neg <- purrr::pmap_dfr(
  symptom_specs,
  fit_symptom_rms,
  data = data_symp,
  hiv_filter = 0
)

pred_hiv_pos <- purrr::pmap_dfr(
  symptom_specs,
  fit_symptom_rms,
  data = data_symp,
  hiv_filter = 1
)

plot_all <- ggplot(
  pred_all,
  aes(dif_day_tx_symp, yhat, color = symptom, fill = symptom)
) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(
    x = "Delay to treatment initiation (days)",
    y = "Odds Ratio",
    title = "All patients"
  ) +
  theme_bw() +
  scale_color_manual(values = c(
    "Bleeding" = "brown2",
    "Jaundice" = "orange2",
    "Oedema"   = "darkgreen"
  )) +
  scale_fill_manual(values = c(
    "Bleeding" = "brown2",
    "Jaundice" = "orange2",
    "Oedema"   = "darkgreen"
  ))

plot_all
  
  
  
  
  
  
  
  
  ######################################### GLM model ###################################################################

var_order_2 <- c("age_group_2", "gender",  "region",
                 "education_level", "race_2", "hiv_2", "name_tx")


m1 <- glm(death_vl ~ rcs(dif_day_tx_symp, 4) ,
          family = "binomial", data_6 %>% filter(!hiv_3 == 1))

summary(m1)

newdat <- data.frame(
  dif_day_tx_symp = seq(min(data_8$dif_day_tx_symp, na.rm = TRUE),
                        max(data_8$dif_day_tx_symp, na.rm = TRUE),
                        length.out = 100)
  #age_years = mean(data_8$age_years, na.rm = TRUE),
  #hiv_3 = 1  # choose one category (0 or 1)
)

newdat$pred <- predict(m1, newdata = newdat, type = "response")

plot(newdat$dif_day_tx_symp, newdat$pred*100, type = "l", lwd = 2,
     xlab = "Days between Treatment and Symptoms",
     ylab = "Predicted Probability of Death (VL)",
     main = "Predicted Death Probability vs. Treatment Delay")


plot(m1)

summary(m1)
############################################ GLM Adjusted #######################################################

library(dplyr)
library(rms)

## 1) Variables to adjust for in the model
var_order_2 <- c("age_group_2", "gender", "region",
                 "education_level", "race_2", "hiv_2", "name_tx")

## 2) Desired level orders (reference = first level)
level_orders <- list(
  age_group_2 = c("0-1", "1-5", "6-15", "16-50", ">50"),
  gender      = c("Female", "Male"),
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
  name_tx = c(
    "Antimonial Pentavalent",
    "Amphotericin B deoxycholate",
    "Liposomal amphotericin B",
    "Others"
  )
)

## 3) Modelling dataset: non-missing outcome and delay, positive delay
data_mort <- data_6 %>%
  filter(!is.na(death_vl),
         !is.na(dif_day_tx_symp),
         dif_day_tx_symp > 0)

## 4) Enforce factor level orders so references are as desired
for (v in var_order_2) {
  if (!is.null(level_orders[[v]])) {
    data_mort[[v]] <- factor(data_mort[[v]], levels = level_orders[[v]])
  } else {
    data_mort[[v]] <- droplevels(as.factor(data_mort[[v]]))
  }
}

## 5) Multivariable logistic model with restricted cubic spline on delay
form_multi <- as.formula(
  paste(
    "death_vl ~ rcs(dif_day_tx_symp, 4) +",
    paste(var_order_2, collapse = " + ")
  )
)

m1 <- glm(form_multi,
          data   = data_mort,
          family = binomial())

summary(m1)

## 6) Prediction grid over delay for a reference profile
# Reference profile:
#   age_group_2     = 0-1
#   gender          = Female
#   region          = Center-West
#   education_level = Incomplete elementary school
#   race_2          = Brown
#   hiv_2           = No
#   name_tx         = Antimonial Pentavalent

delay_seq <- seq(
  quantile(data_mort$dif_day_tx_symp, 0.01, na.rm = TRUE),
  quantile(data_mort$dif_day_tx_symp, 0.99, na.rm = TRUE),
  length.out = 100
)

newdat <- data.frame(
  dif_day_tx_symp = delay_seq,
  age_group_2     = "16-50",
  gender          = "Male",
  region          = "South",
  education_level = "Incomplete elementary school",
  race_2          = "Afro-descendants",
  hiv_2           = "Yes",
  name_tx         = "Antimonial Pentavalent"
)

# make sure factor levels in newdat match data_mort
for (v in var_order_2) {
  newdat[[v]] <- factor(newdat[[v]], levels = levels(data_mort[[v]]))
}

## 7) Predicted probabilities and 95% CI
pred <- predict(m1, newdata = newdat, type = "link", se.fit = TRUE)

newdat$fit    <- pred$fit
newdat$se.fit <- pred$se.fit

newdat <- newdat %>%
  mutate(
    prob  = plogis(fit),
    lower = plogis(fit - 1.96 * se.fit),
    upper = plogis(fit + 1.96 * se.fit)
  )

## 8) Plot: predicted probability of death vs treatment delay
plot(newdat$dif_day_tx_symp, newdat$prob * 100,
     type = "l", lwd = 2,
     xlab = "Days between symptoms and treatment",
     ylab = "Predicted probability of death (%)",
     main = "Predicted death probability vs treatment delay\n(adjusted for covariates)")

lines(newdat$dif_day_tx_symp, newdat$lower * 100, lty = 2)
lines(newdat$dif_day_tx_symp, newdat$upper * 100, lty = 2)

# Optional: rug to show distribution of observed delays
rug(data_mort$dif_day_tx_symp, ticksize = 0.02)


ggplot(newdat, aes(x = dif_day_tx_symp, y = prob)) +
  # 95% CI band
  geom_ribbon(aes(ymin = lower, ymax = upper),
              alpha = 0.2) +
  # predicted probability curve
  geom_line(size = 1) +
  geom_hline(aes(yintercept = 0.05), linetype = "dotted", linewidth = 0.8) +
  geom_hline(aes(yintercept = 0.07), linetype = "dotted", linewidth = 0.8) +
  # optional: rug showing observed delays
  geom_rug(data = data_mort,
           aes(x = dif_day_tx_symp, y = 0),
           inherit.aes = FALSE,
           sides = "b",
           alpha = 0.3) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    name   = "Predicted probability of death (%)",
    n.breaks = 15
  ) +
  labs(
    x     = "Days between symptoms and treatment",
    title = "Predicted death probability vs treatment delay",
    subtitle = "Adjusted for age group, gender, region, education, race, HIV status and treatment"
  ) +
  theme_classic() +
  scale_x_continuous(limits = c(0, 250), n.breaks = 15)



data_mort <- data_m1 %>%
  filter(
    !is.na(death_vl)) %>% 
    filter(!hiv_2 == "Yes")

for (v in var_order_2) {
  if (!is.null(level_orders[[v]])) {
    data_mort[[v]] <- factor(data_mort[[v]], levels = level_orders[[v]])
  } else {
    data_mort[[v]] <- droplevels(as.factor(data_mort[[v]]))
  }
}


############################################################# COX Model Adjusted ##############################################################

## 1) Variables to adjust for in the model
var_order_2 <- c("age_group_2", "gender", "region",
                 "education_level", "race_2")

## 2) Desired level orders (reference = first level)
level_orders <- list(
  age_group_2 = c("0-1", "1-5", "6-15", "16-50", ">50"),
  gender      = c("Female", "Male"),
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
  race_2 = c("Brown", "White", "Afro-descendants", "Yellow", "Indigenous", "Unknown")
  )

dd <- datadist(data_mort)
options(datadist = "dd")

form_multi <- as.formula(
  paste(
    "Surv(dif_treat_outcome, death_vl) ~ rcs(dif_day_tx_symp, 4) +",
    paste(var_order_2, collapse = " + ")
  )
)

m1 <- coxph(form_multi, data = data_mort, x = TRUE)

delay_seq <- seq(
  quantile(data_mort$dif_day_tx_symp, 0.01, na.rm = TRUE),
  quantile(data_mort$dif_day_tx_symp, 0.99, na.rm = TRUE),
  length.out = 100
)

newdat <- data.frame(
  dif_day_tx_symp = delay_seq,
  age_group_2     = "0-1",
  gender          = "Female",
  region          = "Northeast",
  education_level = "Incomplete elementary school",
  race_2          = "Brown"
)

for (v in var_order_2) {
  newdat[[v]] <- factor(newdat[[v]], levels = levels(data_mort[[v]]))
}

pred <- predict(m1, newdata = newdat, type = "lp", se.fit = TRUE)

# Reference = median delay
ref_lp <- predict(
  m1,
  newdata = 
    transform(newdat, dif_day_tx_symp = median(data_mort$dif_day_tx_symp)),
  type = "lp"
)

newdat <- newdat %>%
  mutate(
    hr    = exp(pred$fit - ref_lp),
    lower = exp(pred$fit - ref_lp - 1.96 * pred$se.fit),
    upper = exp(pred$fit - ref_lp + 1.96 * pred$se.fit)
  )

ggplot(newdat, aes(x = dif_day_tx_symp, y = hr)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(
    x = "Days between symptoms and treatment",
    y = "Hazard ratio for death (median delay 29 days)",
    #title = "Hazard of death by treatment delay",
    #subtitle = "Adjusted for age, gender, region, education, race, HIV status, and treatment"
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 14, color = "black", face = "bold"),
    axis.text.x = element_text(size = 12, color = "black")
    ) +
  scale_x_continuous(n.breaks = 10) +
  scale_y_continuous(n.breaks = 12)
  

data_6 %>% 
  summarise(media = median(age_years, na.rm = T), 
            IQU = quantile(age_years, na.rm = T), 
            IQL = quantile(age_years, na.rm = T), 
            rangeU = max(age_years, na.rm = T), 
            rangeL = min(age_years, na.rm = T))


summary(data_6$age_years)




dd <- datadist(data_mort)
options(datadist = "dd")

form_multi <- as.formula(
  paste(
    "Surv(dif_treat_outcome, death_vl) ~ rcs(dif_day_tx_symp, 4) +",
    paste(var_order_2, collapse = " + ")
  )
)

m1 <- coxph(form_multi, data = data_mort, x = TRUE)

delay_seq <- seq(
  quantile(data_mort$dif_day_tx_symp, 0.01, na.rm = TRUE),
  quantile(data_mort$dif_day_tx_symp, 0.99, na.rm = TRUE),
  length.out = 100
)

newdat <- data.frame(
  dif_day_tx_symp = delay_seq,
  age_group_2     = "0-1",
  gender          = "Female",
  region          = "Northeast",
  education_level = "Incomplete elementary school",
  race_2          = "Brown"
)

for (v in var_order_2) {
  newdat[[v]] <- factor(newdat[[v]], levels = levels(data_mort[[v]]))
}

pred <- predict(m1, newdata = newdat, type = "lp", se.fit = TRUE)

# Reference = median delay
ref_lp <- predict(
  m1,
  newdata = 
    transform(newdat, dif_day_tx_symp = median(data_mort$dif_day_tx_symp)),
  type = "lp"
)

newdat <- newdat %>%
  mutate(
    hr    = exp(pred$fit - ref_lp),
    lower = exp(pred$fit - ref_lp - 1.96 * pred$se.fit),
    upper = exp(pred$fit - ref_lp + 1.96 * pred$se.fit)
  )

ggplot(newdat, aes(x = dif_day_tx_symp, y = hr)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(
    x = "Days between symptoms and treatment",
    y = "Hazard ratio for death (median delay 29 days)",
    #title = "Hazard of death by treatment delay",
    #subtitle = "Adjusted for age, gender, region, education, race, HIV status, and treatment"
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 14, color = "black", face = "bold"),
    axis.text.x = element_text(size = 12, color = "black")
  ) +
  scale_x_continuous(n.breaks = 10) +
  scale_y_continuous(n.breaks = 12)
