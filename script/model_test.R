#Analysis model using delay category 
install.packages(c("survival","rms","coxme","mice","Hmisc", "sjPlot", "splines", "broom"))
library(survival)
library(rms)    # for rcs / cph
library(coxme)  # mixed effects cox
library(mice)
library(Hmisc)
library(survminer)
library(sjPlot)
library(splines)

pacman::p_load(broom)

######################################## Adjust variables and set up database ###################################################################
### Import work database final for running up the model
data_6 <- import(here("data base and clean", "data_6.rds")) 

data_6 <- data_6 %>%
  dplyr::mutate(across(where(~ inherits(.x, "IDate")), as.Date))

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

######################################## Adjust variables and set up database opt 2 ##############################################################################

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
  dplyr::select(dif_day_tx_symp, dif_treat_start_symp_cat_2, dif_treat_outcome, age_group_2,
         age_years, gender_2, hiv_3, hiv_2, gender, death_vl, name_tx, severe_2, race_2,
         federative_unit, education_level, , race_2) 
#filter(!federative_unit == state_100)

export(data_8, here("data base and clean", "data_8.rds")) ### Database to work with time between treatment initiation and death by VL

data_8 <- data_8 %>%
  dplyr::mutate(across(where(~ inherits(.x, "IDate")), as.Date))

######################################## Non-adjust cox regression model by HIV status and Overall ###################################################################################################

### Model overall - including HIV patients and unadjusted

data_m1 <- data_8 %>% 
  #filter(evolution %in% c("Cure", "Death for VL")) %>% 
  filter(!is.na(death_vl)) %>% 
  filter(!hiv_2 == "Yes") %>% 
  dplyr::select(dif_day_tx_symp, dif_treat_start_symp_cat_2, dif_treat_outcome, age_group_2,
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
  #filter(evolution %in% c("Cure", "Death for VL")) %>% 
  #filter(!hiv_2 == "Yes") %>% 
  dplyr::select(dif_day_tx_symp, dif_treat_start_symp_cat_2, dif_treat_outcome, age_group_2,
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

data_8.2 <- data_8 %>% 
  #filter(evolution %in% c("Cure", "Death for VL")) %>% 
  filter(!is.na(death_vl)) %>% 
  filter(!hiv_2 == "Yes") 
  

m1_hiv <- coxph(Surv(dif_treat_outcome, death_vl) ~ 
              rcs(dif_day_tx_symp, 4), 
            data = data_8.2
)

summary(m1_hiv)

newd_hiv <- data.frame(
  dif_day_tx_symp = seq(
    min(data_8.2$dif_day_tx_symp, na.rm = TRUE),
    max(data_8.2$dif_day_tx_symp, na.rm = TRUE),
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
######################################## Adjusted cox regression model using coxph and predict functions#######################################
#### Model created to evaluate direction and distribution of delay treatment initiation adjusted by 
#### age, gender, region, education level and race

## 1) Variables to adjust for and factor level ordering
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

## 2) Modelling dataset

data_mort <- data_m1

## 3) Enforce factor level orders (reference = first level)

for (v in var_order_2) {
  if (!is.null(level_orders[[v]])) {
    data_mort[[v]] <- factor(data_mort[[v]], levels = level_orders[[v]])
  } else {
    data_mort[[v]] <- droplevels(as.factor(data_mort[[v]]))
  }
}

## 4) Cox proportional hazards model with spline on treatment delay

m2 <- coxph(
  Surv(dif_treat_outcome, death_vl) ~
    rcs(dif_day_tx_symp, k = 4) +
    age_group_2 + gender + race_2 + region +
    education_level,
  data = data_8
)

str(m2)

summary(m2)

## 5) Prediction grid over treatment delay (reference profile fixed)

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

## 6) Adjusted hazard ratios (reference = shortest delay)

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

## 7) Publication-ready plot: hazard ratio vs treatment delay

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
  
  
######################################## Logistic regression model by severe symptoms ############################################################################################

### Selecting symptoms associated with severity and presentation 
### 1. Patients HIV positive

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
  filter(!hiv_2 == "Yes")

### Bleeding model 1

data_bleeding <- data_symp_no_hiv %>% 
  dplyr::select(bleeding, dif_day_tx_symp) %>% 
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

### Jaundice model 1
data_jaundice <- data_symp_no_hiv %>% 
  dplyr::select(jaundice_2, dif_day_tx_symp) %>% 
  dplyr::filter(!is.na(jaundice_2))

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

### oedema model 1

data_oedema <- data_symp_no_hiv %>%
  dplyr::select(dif_day_tx_symp, edema_3) %>% 
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
  


### 1. Patients not HIV positive

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
  filter(!hiv_2 == "No")

### Bleeding model 2

data_bleeding_2 <- data_symp_no_hiv %>% 
  dplyr::select(bleeding, dif_day_tx_symp) %>% 
  filter(!is.na(bleeding))

dd <- datadist(data_bleeding_2)
options(datadist = "dd")

bleeding_fit <- lrm(bleeding ~ rcs(dif_day_tx_symp, 4), data = data_bleeding_2)

p_bleeding_2 <- Predict(
  bleeding_fit,
  dif_day_tx_symp,
  fun = exp,        # odds ratios
  ref.zero = TRUE   # reference = median
)

p_bleeding_df_2 <- as.data.frame(p_bleeding_2)

bleeding_m2 <- ggplot(p_bleeding_df_2, aes(dif_day_tx_symp, yhat)) +
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

bleeding_m2

### Jaundice model 1
data_jaundice_2 <- data_symp_no_hiv %>% 
  dplyr::select(jaundice_2, dif_day_tx_symp) %>% 
  dplyr::filter(!is.na(jaundice_2))

dd <- datadist(data_jaundice_2)
options(datadist = "dd")

jaundice_fit_2 <- lrm(jaundice_2 ~ rcs(dif_day_tx_symp, 4), data = data_jaundice_2)

p_jaundice_2 <- Predict(
  jaundice_fit_2,
  dif_day_tx_symp,
  fun = exp,        # odds ratios
  ref.zero = TRUE   # reference = median
)

p_jaundice_df_2 <- as.data.frame(p_jaundice_2)

jaundice_m2 <- ggplot(p_jaundice_df_2, aes(dif_day_tx_symp, yhat)) +
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

jaundice_m2

### oedema model 1

data_oedema_2 <- data_symp_no_hiv %>%
  dplyr::select(dif_day_tx_symp, edema_3) %>% 
  filter(!is.na(edema_3))

#oedema_m <- glm(formula = edema_3 ~ rcs(dif_day_tx_symp, 4), family = "binomial",
#                 data = data_oedema)

edema_fit_2 <- lrm(edema_3 ~ rcs(dif_day_tx_symp, 4), data = data_oedema_2)

p_edema_2 <- Predict(
  edema_fit_2,
  dif_day_tx_symp,
  fun = exp,        # odds ratios
  ref.zero = TRUE   # reference = median
)

p_edema_df_2 <- as.data.frame(p_edema_2)

edema_2 <- ggplot(p_edema_df_2, aes(dif_day_tx_symp, yhat)) +
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

edema_2

symptoms_m2 <- ggplot(
  data = p_bleeding_df_2,
  aes(dif_day_tx_symp, yhat)
) +
  geom_line(aes(color = "Bleeding"), linewidth = 1) +
  geom_ribbon(
    aes(ymin = lower, ymax = upper),
    alpha = 0.2, fill = "brown2"
  ) +
  
  geom_line(
    data = p_jaundice_df_2,
    aes(dif_day_tx_symp, yhat, color = "Jaundice"),
    linewidth = 1,
  ) +
  geom_ribbon(
    data = p_jaundice_df_2,
    aes(
      x = dif_day_tx_symp,
      ymin = lower,
      ymax = upper,
    ), fill = "orange2",
    alpha = 0.2,
    inherit.aes = FALSE
  ) +
  
  geom_line(
    data = p_edema_df_2,
    aes(dif_day_tx_symp, yhat, color = "Oedema"),
    linewidth = 1,
  ) +
  geom_ribbon(
    data = p_edema_df_2,
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

symptoms_m2


#### 3. Overall model (HIV, No HIV and Unknown HIV)

data_symp_overall_hiv <- data_8 %>% 
  mutate(jaundice_2 = case_when(jaundice == 1 ~ 1,
                                jaundice == 2 ~ 0,
                                TRUE ~ NA),
         bleeding = case_when(hemorrhagic_phenomena == 1 ~ 1,
                              hemorrhagic_phenomena == 2 ~ 0,
                              TRUE ~ NA),
         edema_3 = case_when(edema == 1 ~ 1,
                             edema == 2 ~ 0,
                             TRUE ~NA)) 

### Bleeding model option 3

data_bleeding_3 <- data_symp_overall_hiv %>% 
  dplyr::select(bleeding, dif_day_tx_symp) %>% 
  filter(!is.na(bleeding))

dd <- datadist(data_bleeding_3)
options(datadist = "dd")

bleeding_fit_3 <- lrm(bleeding ~ rcs(dif_day_tx_symp, 4), data = data_bleeding_3)

p_bleeding_3 <- Predict(
  bleeding_fit_3,
  dif_day_tx_symp,
  fun = exp,        # odds ratios
  ref.zero = TRUE   # reference = median
)

p_bleeding_df_3 <- as.data.frame(p_bleeding_3)

bleeding_m3 <- ggplot(p_bleeding_df_3, aes(dif_day_tx_symp, yhat)) +
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

bleeding_m3

### Jaundice
data_jaundice_3 <- data_symp_overall_hiv %>% 
  dplyr::select(jaundice_2, dif_day_tx_symp) %>% 
  filter(!is.na(jaundice_2))

dd <- datadist(data_jaundice_3)
options(datadist = "dd")

jaundice_fit_3 <- lrm(jaundice_2 ~ rcs(dif_day_tx_symp, 4), data = data_jaundice_3)

p_jaundice_3 <- Predict(
  jaundice_fit_3,
  dif_day_tx_symp,
  fun = exp,        # odds ratios
  ref.zero = TRUE   # reference = median
)

p_jaundice_df_3 <- as.data.frame(p_jaundice_3)

jaundice_m3 <- ggplot(p_jaundice_df_3, aes(dif_day_tx_symp, yhat)) +
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

jaundice_m3

### oedema

data_oedema_3 <- data_symp_overall_hiv %>%
  dplyr::select(dif_day_tx_symp, edema_3) %>% 
  filter(!is.na(edema_3))

#oedema_m <- glm(formula = edema_3 ~ rcs(dif_day_tx_symp, 4), family = "binomial",
#                 data = data_oedema)

edema_fit_3 <- lrm(edema_3 ~ rcs(dif_day_tx_symp, 4), data = data_oedema_3)

p_edema_3 <- Predict(
  edema_fit_3,
  dif_day_tx_symp,
  fun = exp,        # odds ratios
  ref.zero = TRUE   # reference = median
)

p_edema_df_3 <- as.data.frame(p_edema_3)

edema_3 <- ggplot(p_edema_df_3, aes(dif_day_tx_symp, yhat)) +
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

edema_3

symptoms_m3 <- ggplot(
  data = p_bleeding_df_3,
  aes(dif_day_tx_symp, yhat)
) +
  geom_line(aes(color = "Bleeding"), linewidth = 1) +
  geom_ribbon(
    aes(ymin = lower, ymax = upper),
    alpha = 0.2, fill = "brown2"
  ) +
  
  geom_line(
    data = p_jaundice_df_3,
    aes(dif_day_tx_symp, yhat, color = "Jaundice"),
    linewidth = 1,
  ) +
  geom_ribbon(
    data = p_jaundice_df_3,
    aes(
      x = dif_day_tx_symp,
      ymin = lower,
      ymax = upper,
    ), fill = "orange2",
    alpha = 0.2,
    inherit.aes = FALSE
  ) +
  
  geom_line(
    data = p_edema_df_3,
    aes(dif_day_tx_symp, yhat, color = "Oedema"),
    linewidth = 1,
  ) +
  geom_ribbon(
    data = p_edema_df_3,
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

final_models <- ggarrange(symptoms_m1 + rremove("ylab") + rremove("xlab") + 
                            scale_y_continuous(limits = c(0.5, 5)),
                          symptoms_m2 + rremove("xlab") + rremove("ylab") +
                            scale_y_continuous(limits = c(0.5, 5)), 
                          symptoms_m3 + rremove("xlab") + rremove("ylab") +
                            scale_y_continuous(limits = c(0.5, 5)), labels = "AUTO",
                          nrow = 1, legend = "right", hjust = -0.5,
                          common.legend = T,
                          font.label = list(size = 14, color = "black",
                                            face = "bold", family = "Times New Roman",
                                            position = "bottom")) 

symptoms <- annotate_figure(final_models, 
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
ggsave(symptoms, path = here("figures and tables"),
       width = 10, height = 8, dpi = 300,
       file = "severe_symptoms_color_graph.png")
######################################## Model 2 (figure S2) based on sever and complications taken cardinal symptoms ###############################################################

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
    dplyr::select(all_of(c(outcome, "dif_day_tx_symp"))) %>%
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

######################################## Model by symptoms and HIV status #############################################################################

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
  
################################################# GLM model ##################################################################################################

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
############################################ GLM Adjusted ###########################################################################################

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


############################################################# COX Model Adjusted ############################################################################

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
