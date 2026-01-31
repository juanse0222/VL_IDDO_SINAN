### Import dataset by SINAN 
data <- import("data base and clean/data_raw/combined_dataset_after_drop.csv") 

### 1. Classification criteria: Only patients that have parasitological or molecular confirmation 

data %>% 
  filter(!type_of_cases_at_baseline == 1) %>% 
  tabyl(type_of_cases_at_baseline) %>% 
  adorn_totals()

data_3 <- data %>% 
  clean_names() %>% 
  filter(final_calssification == 1 & confirmation_criteria == 1) %>% 
  filter(type_of_cases_at_baseline == 1) 

data_3 <- data_3 %>% 
  mutate(name_tx = case_when(type_of_rx == 1 ~ "Antimonial Pentavalent",
                             type_of_rx == 2 ~ "Amphotericin B deoxycholate",
                             type_of_rx == 3 ~ "Pentamidine",
                             type_of_rx == 4 ~ "Liposomal amphotericin B",
                             type_of_rx == 5 ~ "Others",
                             type_of_rx == 6 ~ "Nothing",
                             TRUE ~ NA))

#### 1.1 Review data based on name of treatments and dates to identify gaps or errors
data_3 %>% 
  summarise(n = sum(name_tx %in% c(NA, "Nothing")))

data_3 %>%
  mutate(equal = if_else(date_symptoms_started == date_of_rx_started, 1L, 0L)) %>% 
  summarise(n = sum(equal, na.rm = T))

#### 2. Apply the rules to choose only cases fulfilling with it. 
data_4 <- data_3 %>% 
  filter(!name_tx %in% c(NA, "Nothing")) %>% ### Delete registers without treatment's name 
  filter(!is.na(date_of_rx_started)) %>% #### To make sure that the patients chosen everyone have date of treatment initiation
  filter(!is.na(date_symptoms_started)) %>%  ### Date of symptom onset for each register selected
  filter(date_symptoms_started != clean_date_of_rx_started) 
#### Don't have the same date of the treatment initiation regard to symptoms onset

#### 3. 

data_4 <- data_4 %>% 
mutate(death = if_else(!is.na(date_of_death), 1, 0), ### We need to evaluate if this information is correct. 
       cure = if_else(evolution_of_cases == 1, 1, 0),
       death_vl = if_else(evolution_of_cases == 3, 1, 0),
       pregnant = case_when(gestational_age == 1 ~ "1 Trimester",
                            gestational_age == 2 ~ "2 Trimester",
                            gestational_age == 3 ~ "3 Trimester",
                            gestational_age == 4 ~ "Unknown",
                            gestational_age == 5 ~ "No",
                            gestational_age %in% c(6,9) ~ "N/A",
                            TRUE ~ NA),
       race_2 = case_when(race == 1 ~ "White",
                          race == 2 ~ "Afro-descendants",
                          race == 3 ~ "Yellow",
                          race == 4 ~ "Brown",
                          race == 5 ~ "Indigenous",
                          race == 9 ~ "Unknown",
                          TRUE ~ NA),
       education_level = case_when(level_of_education == 1 ~ "Incomplete elementary school",
                                   level_of_education %in% c(2,3) ~ "Complete elementary and incomplete middle school",
                                   level_of_education %in% c(4,5) ~ "Complete middle school and incomplete high school",
                                   level_of_education == 6 ~ "Complete high school",
                                   level_of_education == 7 ~ "Incomplete university",
                                   level_of_education == 8 ~ "Complete university",
                                   # level_of_education == 9 ~ "Unknown",
                                   level_of_education == 10 ~ "N/A",
                                   TRUE ~ NA),
       identification_case = case_when(type_of_cases_at_baseline == 1 ~ "New case",
                                       type_of_cases_at_baseline == 2 ~ "Releapse",
                                       type_of_cases_at_baseline == 3 ~ "Transferred",
                                       # type_of_cases_at_baseline == 9 ~ "No data",
                                       TRUE ~ NA),
       death_any_causes = if_else(death == 1, T, F),
       death_vl_2 = if_else(death_vl == 1, T, F),
       evolution = case_when(evolution_of_cases == 1 ~ "Cure",
                             evolution_of_cases == 2 ~ "Neglect",
                             evolution_of_cases == 4 ~ "Death for any-cause",
                             evolution_of_cases == 3 ~ "Death for VL",
                             evolution_of_cases == 5 ~ "Transferred",
                             # is.na(evolution_of_cases) ~ "Unknown",
                             TRUE ~ NA)) 

### 4. Adjust dates and evaluation/review differences of days between symptoms, diagnosis and treatment

data_4.1 <- data_4 %>%
  mutate(across(c(year_of_notification, date_of_notification, date_symptoms_started,
                  date_of_investigation, date_of_rx_started, date_of_death), as.Date)) %>%
  
  # 2) Keep your treatment labels
  mutate(
    death    = if_else(!is.na(date_of_death), 1L, 0L),
    cure     = if_else(evolution_of_cases == 1, 1L, 0L),
    death_vl = if_else(evolution_of_cases == 3, 1L, 0L)
    ) %>%
  
  # 3) Compute raw diffs once
  mutate(
    dif_day_tx_symp = as.numeric(date_of_rx_started - date_symptoms_started),
    dif_day_dx_symp = as.numeric(date_of_investigation - date_symptoms_started),
    dif_day_tx_dx   = as.numeric(date_of_rx_started - date_of_investigation)) %>%
  
  # 4) Fix investigation date (your swap + validity rules) for only “problem rows”
  mutate(
    date_not = date_of_notification,
    date_symp = date_symptoms_started,
    date_inv = date_of_investigation,
    date_tx = date_of_rx_started,
    inv_swapped = make_date(year(date_of_investigation),
                            day(date_of_investigation),
                            month(date_of_investigation)),
    
    inv_candidate = case_when(
      (dif_day_dx_symp < 0) ~ inv_swapped,          # your earlier correction
      (dif_day_tx_dx < 0)   ~ inv_swapped,          # tx before dx → try swap
      TRUE                  ~ date_of_investigation
    ),
    
    inv_candidate = coalesce(inv_candidate, date_of_notification),
    tx_dx   = as.numeric(date_of_rx_started - inv_candidate),
    dx_symp = as.numeric(inv_candidate - date_symptoms_started),
    inv_candidate = if_else(dx_symp < 0, date_of_notification, inv_candidate),
    inv_candidate = if_else(tx_dx < 0, date_of_notification, inv_candidate),
    dx_symp = as.numeric(inv_candidate - date_symptoms_started),
    tx_dx = as.numeric(date_of_rx_started - inv_candidate),
    
    # rule: invalid if tx_dx <0, tx_dx > 50, or dx_symp <0 → use notification date
    inv_candidate = if_else(tx_dx > 50, date_of_notification,
                            inv_candidate),
    tx_dx   = as.numeric(date_of_rx_started - inv_candidate),
    
    # 5) Final diffs based on corrected dx date
    date_of_investigation = inv_candidate,
    dif_day_tx_dx   = tx_dx,
    dif_day_dx_symp = dx_symp
  ) %>%
  # 6) Clean helper columns if you don’t want them
  select(-inv_swapped, -inv_candidate, -tx_dx, -dx_symp)

quantile(data_4.1$dif_day_dx_symp, 0.99)
quantile(data_4.1$dif_day_tx_dx, 0.99)
quantile(data_4.1$dif_day_tx_symp, 0.99)

data_4.1 %>% 
  ggplot() +
  geom_histogram(aes(x = dif_day_tx_symp)) +
  scale_x_log10()

### Identify 99% distribution
quantile(as.numeric(data_4.1$dif_day_tx_symp), probs = 0.99, na.rm = TRUE)

### adjusting database using 99% distribution
data_5 <- data_4.1 %>% 
  #mutate(code_id = paste(date_of_birth, gender, code_of_the_municipality_of_residence, regional_health_area, 
   #                      race, municipality_code_where_a_person_get_infected, code_of_the_municipality, 
    #                     code_of_the_health_unit, health_region_of_residence)) %>% 
  #distinct(clean_date_of_birth, clean_date_of_notification, clean_date_of_investigation, 
  #        gender, age_as_pt_repored, .keep_all = T) %>%   ### Duplicated registers using unique code
  filter(dif_day_tx_symp < quantile(as.numeric(data_4.1$dif_day_tx_symp), probs = 0.99, na.rm = TRUE)) %>% 
  filter(!dif_day_tx_dx < 0)
#filter(!is.na(evolution_of_cases)) #### Only it will be used for the model

data_5 %>% 
  tabyl(dif_day_tx_dx)

tabyl(data_4$evolution) ### Distribution of evolution variable


### Categorising time between symptom onset to treatment initiation, and creating a variable between treatment to outcome (death or alive)
data_6 <- data_5 %>% 
  mutate(dif_treat_start_symp_cat = case_when(dif_day_tx_symp < 30 ~ "< 30 days",
                                              dif_day_tx_symp >= 30 & dif_day_tx_symp < 90 ~ "30-89 days",
                                              dif_day_tx_symp >=90  & dif_day_tx_symp < 120 ~ "90-119 days",
                                              dif_day_tx_symp >= 120  ~ "≥120 days",
                                              TRUE ~ NA),
         dif_treat_start_symp_cat = factor(dif_treat_start_symp_cat, levels = c("< 30 days", "30-89 days", "90-119 days", "≥120 days")),
         dif_treat_death = as.numeric(date_of_death) - as.numeric(date_of_rx_started),
         day_of_outcome = if_else(is.na(date_of_death), as.Date(date_of_close_out), date_of_death),
         dif_treat_outcome = as.numeric(day_of_outcome) - as.numeric(date_of_rx_started),
         dif_day_tx_symp = as.numeric(dif_day_tx_symp)
  ) %>% 
  filter(!dif_treat_outcome <0 & dif_treat_outcome < 366) 


### Recodification, and standardize age groups, fever and HIV.
data_6 <- data_6 %>% 
  mutate(age_group_2 = case_when(age_years >= 0 & age_years <= 1 ~ "0-1",
                                 age_years > 1 & age_years <= 5 ~ "1-5",
                                 age_years > 5 & age_years <= 15 ~ "5-15",
                                 age_years > 15 & age_years <= 45 ~ "15-45",
                                 age_years > 45  ~ ">45",
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
         age_group_3 = factor(age_group_3, levels = c("0-1", "1-2", "2-3", "3-5", "5-15",
                                                      "15-20", "20-40", "40-60",
                                                      ">60")),
         age_group_2 = factor(age_group_2, levels = c("0-1", "1-5", "5-15", "15-45", ">45")),
         fever_2 = case_when(fever == 1 ~ "Yes",
                             fever == 2 ~ "No",
                             fever == 9 ~ "Unknown",
                             TRUE ~ NA),
         hiv_2 = case_when(hiv == 1 ~ "Yes",
                           hiv == 2 ~ "No",
                           hiv == 9 ~ "Unknown",
                           TRUE ~ NA),
         fever_2 = factor(fever_2, levels = c("No", "Yes", "Unknown")),
         hiv_2 = factor(hiv_2, levels = c("No", "Yes")))

data_7 <- data_6 %>%
  mutate(tx_label = case_when(
    name_tx == "Amphotericin B deoxycholate" ~ "Amphotericin B Deoxycholate",
    name_tx == "Antimonial Pentavalent" ~ "Antimonial Pentavalent",
    name_tx == "Liposomal amphotericin B" ~ "Liposomal Amphotericin B",
    TRUE ~ name_tx),
    delay_cat_label = case_when(
      dif_treat_start_symp_cat == "< 30 days" ~ "<30 days",
      dif_treat_start_symp_cat == "30-89 days" ~ "30–89 days",
      dif_treat_start_symp_cat == "90-119 days" ~ "90–119 days",
      dif_treat_start_symp_cat == "> 120 days" ~ "≥120 days",
      TRUE ~ dif_treat_start_symp_cat),
    strata_combined = paste(tx_label, "-", delay_cat_label))


### Finala database 
rio::export(data_6, here("data base and clean", "data_6.rds"))
rio::export(data_6, here("data base and clean", "data_6.csv"))
###############################################################################################################

