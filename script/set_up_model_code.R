
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
  mutate(complication = if_else(hemorrhagic_phenomena == 1 | edema == 1, "Yes", "No"))

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
                                                TRUE ~ NA)) %>% 
  select(dif_day_tx_symp, dif_treat_start_symp_cat_2, dif_treat_outcome, 
         age_years, gender_2, hiv_3, hiv_2, gender, death_vl, name_tx, severe_2, federative_unit) 
  #filter(!federative_unit == state_100)
##########################################################################################################################################################

### Cox regression model 

dd <- datadist(data_8)
options(datadist = "dd")  ### rcs

### different models applied for get to know wich one has better adjustment. 
fit_rms <- coxph(Surv(dif_treat_outcome, death_vl) ~ rcs(dif_day_tx_symp, 4) + rcs(age_years, 4) + 
                   hiv_3 + gender_2 + severe_2,
                 data = data_8)


fit_rms_2 <- coxph(Surv(dif_treat_outcome, death_vl) ~ rcs(dif_day_tx_symp, 4) + rcs(age_years, 4) + hiv_3,
                 data = data_8)

fit_rms_3 <- coxph(Surv(dif_treat_outcome, death_vl) ~ ns(dif_day_tx_symp, df = 4) + ns(age_years, df = 4) + hiv_3,
                   data = data_8)

fit_rms_4 <- coxph(Surv(dif_treat_outcome, death_vl) ~ ns(dif_day_tx_symp, df = 4) + ns(age_years, df = 4) + hiv_3,
                   data = data_8, cluster = federative_unit)


summary(fit_rms)


plot(survfit(fit_rms_2), ylim=c(0.9, 1), xlab="Weeks",
     ylab="Proportion of mortality by VL")


fit_spline_cluster <- coxph(
  Surv(dif_treat_outcome, death_vl) ~ (ns(dif_day_tx_symp, df = 4) + ns(age_years, df = 4) +
                                         hiv_2 + gender), data = data_8)

plot(survfit(fit_rms_2)$cumhaz)

sfit <- survfit(fit_rms)
cumhaz.upper <- -log(sfit$upper)
cumhaz.lower <- -log(sfit$lower)
cumhaz <- sfit$cumhaz # same as -log(sfit$surv)

plot(cumhaz, xlab="Delay (Days between symptoms onset to treatment)", ylab="cumulative hazard",
     ylim=c(min(cumhaz.lower), 0.09))
lines(cumhaz.lower)
lines(cumhaz.upper)

ggforest(fit_rms)


basehaz(fit_spline_cluster)

baseline_hazard <-  basehaz(fit_spline_cluster)

baseline_hazard$survival <- exp(-baseline_hazard$hazard)



nd <- predict()

fit_rms_2 <- coxph(Surv(dif_treat_outcome, death_vl) ~ rcs(dif_day_tx_symp, 4)
                   + rcs(age_years, 4) + hiv_3,
                   data = data_8)


# grid for delay
delay_seq <- seq(min(data_8$dif_day_tx_symp, na.rm=TRUE),
                 max(data_8$dif_day_tx_symp, na.rm=TRUE),
                 length.out = 365)

# choose typical values for other covariates: median age, mode for categorical
ref_age <- median(data_8$age_years, na.rm=TRUE)
mode_val <- function(x) { ux <- unique(na.omit(x)); ux[which.max(tabulate(match(x, ux)))] }
#ref_hiv <- mode_val(data_8$hiv_3)    # note you used hiv_2 in model
ref_hiv <- "1"
#ref_gender <- mode_val(data_8$gender)

# construct newdata (100+ rows, each col same length)
newd <- data.frame(
  dif_day_tx_symp = delay_seq,
  age_years = rep(ref_age, length(delay_seq)),
  hiv_3 = if (is.factor(data_8$hiv_3)) factor(rep(ref_hiv, length(delay_seq)), levels = levels(data_8$hiv_3)) else rep(ref_hiv, length(delay_seq))
  #gender = if (is.factor(data_8$gender)) factor(rep(ref_gender, length(delay_seq)), levels = levels(data_8$gender)) else rep(ref_gender, length(delay_seq))
)

# predict linear predictor and SE from coxph
pred_lp <- predict(fit_rms_2, newdata = newd, type = "lp", se.fit = TRUE)

newd$lp <- pred_lp$fit
newd$se_lp <- pred_lp$se.fit

# convert to HR (and CI)
newd$HR <- exp(newd$lp)
newd$HR_low <- exp(newd$lp 
                   - 1.96 * newd$se_lp)
newd$HR_high <- exp(newd$lp + 1.96 * newd$se_lp)

# normalize relative to a reference day (e.g. 30 days). Change ref_day if desired.
ref_day <- 30
idx_ref <- which.min(abs(delay_seq - ref_day))
ref_hr <- newd$HR[idx_ref]

newd$HR_rel <- newd$HR / ref_hr
newd$HR_low_rel <- newd$HR_low / ref_hr
newd$HR_high_rel <- newd$HR_high / ref_hr

# ggplot
ggplot(newd, aes(x = dif_day_tx_symp, y = HR_rel)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = HR_low_rel, ymax = HR_high_rel), alpha = 0.2) +
  geom_vline(xintercept = c(16, 30, 60), linetype = "dashed", color = "black") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  labs(x = "Delay (days from symptom onset to treatment)",
       y = paste0("Adjusted hazard ratio (ref = ", delay_seq[idx_ref], " days)"),
       title = "Adjusted HR vs Delay (spline model)") +
  theme_bw(base_size = 14) +
  scale_x_continuous(n.breaks = 15) +
  scale_y_continuous(n.breaks = 10) +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, color = "black", face = "bold"),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
  )


fdd <- datadist(data_8); options(datadist = "dd")
fit_rms_cph <- cph(Surv(dif_treat_outcome, death_vl) ~ rcs(dif_day_tx_symp, 4) + rcs(age_years,4) + 
                     hiv_3, data = data_8, x = TRUE, y = TRUE, surv = TRUE)

p_hiv3_1 <- Predict(fit_rms_cph, dif_day_tx_symp, hiv_3 = factor("1", levels = levels(data_8$hiv_3)), fun = exp)
plot(p_hiv3_1, ylab = "Hazard ratio (hiv_3 = 1)")

summary(fit_rms_cph)

# then you can do Predict + plot
plot(Predict(fit_rms_cph, dif_day_tx_symp, fun = exp),
     xlab = "Delay days",
     ylab = "Hazard ratio",
     rug = TRUE,          # adds tick marks for data points
     ticknumber = 10)

plot(p_hiv3_1)

######################################### GLM model ###################################################################

var_order_2 <- c("age_group_2", "gender",  "region",
                 "education_level", "race_2", "hiv_2")


m1 <- glm(death_vl ~ rcs(dif_day_tx_symp, 4) ,
          family = "binomial", data_6)

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



























fit_coxme <- coxme(Surv(dif_treat_outcome, death_vl) ~ rcs(dif_day_tx_symp,4) + rcs(age_years, 4) +  
                        hiv_3 + gender_2 + (1|federative_unit),
                   data = data_8)


summary(fit_coxme)


fit_coxme_strat <- coxme(
  Surv(dif_treat_outcome, death_vl) ~ rcs(age_years,4) + hiv_3 +
    strata(dif_treat_start_symp_cat_2) + (1 | federative_unit),
  data = data_8
)
summary(fit_coxme_strat)


days_of_interest <- c(30, 60, 90, 120, 200)
nd_oi <- data.frame(
  dif_day_tx_symp = days_of_interest,
  age_years = rep(ref_age, length(days_of_interest)),
  hiv_2 = if (is.factor(data_8$hiv_2)) factor(rep(ref_hiv, length(days_of_interest)), levels = levels(data_8$hiv_2)) else rep(ref_hiv, length(days_of_interest)),
  gender = if (is.factor(data_8$gender)) factor(rep(ref_gender, length(days_of_interest)), levels = levels(data_8$gender)) else rep(ref_gender, length(days_of_interest))
)
pred_oi <- predict(fit_spline_cluster, newdata = nd_oi, type = "lp", se.fit = TRUE)
HRs <- exp(pred_oi$fit)
low <- exp(pred_oi$fit - 1.96*pred_oi$se.fit)
high <- exp(pred_oi$fit + 1.96*pred_oi$se.fit)
result_table <- data.frame(day = days_of_interest, HR = HRs, CI_low = low, CI_high = high)
result_table
# To show relative to 30 days:
ref_val <- result_table$HR[result_table$day==30]
result_table$HR_rel <- result_table$HR / ref_val
result_table$CI_low_rel <- result_table$CI_low / ref_val
result_table$CI_high_rel <- result_table$CI_high / ref_val
result_table
###################################################################################################################################

data_6 %>% 
  filter(!is.na(name_tx)) %>%
  filter(!is.na(evolution)) %>% 
  filter(evolution %in% c("Cure", "Death for VL")) %>% 
  mutate(time = if_else(dif_treat_outcome <= 7, "Yes", "No")) %>% 
  tabyl(evolution, time) %>% 
  adorn_totals("both") %>% 
  adorn_percentages(denominator = "col") %>% 
  adorn_pct_formatting() %>% 
  adorn_ns(position = "front")

summary(data_6$dif_day_tx_symp)
