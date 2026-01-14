pacman::p_load(flextable, tidyverse, janitor, epikit, gtsummary, rio, 
               survival, survminer, sf, here, incidence, scales)
############################################################################################################################

### Import dataset by SINAN 
data <- import("data base and clean/data_raw/combined_dataset_after_drop.csv") 
data <- data %>% 
  #filter(!is.na(age_years)) %>% 
  mutate(year_of_notification = as.Date(year_of_notification)) %>% 
  mutate(age_group = age_categories(x = age_years, lower = min(age_years, na.rm = T),
                                    upper = 90, by = 10)) 


data %>% summarise(n = sum(is.na(age_years))) ### Number of registers without age. 
summary(data$date_of_notification) ### Evalution about date


### Evaluation dataset
colnames(data)
summary(data)

### Adjust names treatment according with SINA
data <- data %>% 
  mutate(name_tx = case_when(type_of_rx == 1 ~ "Antimonial Pentavalent",
                             type_of_rx == 2 ~ "Amphotericin B deoxycholate",
                             type_of_rx == 3 ~ "Pentamidine",
                             type_of_rx == 4 ~ "Liposomal amphotericin B",
                             type_of_rx == 5 ~ "Others",
                             type_of_rx == 6 ~ "Nothing",
                             TRUE ~ NA))

### Reviewing dates and days between onset symptoms and start treatment. 
data_2 <- data %>%
  mutate(dif_day_tx_symp = as.Date(date_of_rx_started) - as.Date(date_symptoms_started),
         dif_day_dx_symp = as.Date(date_of_investigation) - as.Date(date_symptoms_started),
         dif_day_tx_dx = as.Date(date_of_rx_started) - as.Date(date_of_investigation),
         death = if_else(!is.na(date_of_death), 1, 0),
         cure = if_else(evolution_of_cases == 1, 1, 0),
         death_vl = if_else(evolution_of_cases == 3, 1, 0)) %>%
  mutate(date_of_investigation = if_else(dif_day_dx_symp <0,
                                         make_date(year(date_of_investigation), 
                                                   day(date_of_investigation), 
                                                   month(date_of_investigation)),
                                         as.Date(date_of_investigation))) ## Adjust date of investigation for working forward with it.


data_2 %>% 
  tabyl(dif_day_dx_symp)

### Classification criteria: Only patients that have parasitological or molecular confirmation 

data_2 %>% 
  filter(!type_of_cases_at_baseline == 1) %>% 
  tabyl(type_of_cases_at_baseline) %>% 
  adorn_totals()

data_3 <- data_2 %>% 
  filter(final_calssification == 1 & confirmation_criteria == 1) %>% 
  filter(type_of_cases_at_baseline == 1)

data_3 %>% 
  summarise(n = sum(name_tx %in% c(NA, "Nothing")))

data_3 %>%
  mutate(equal = if_else(date_symptoms_started == date_of_rx_started,
                         1L, 0L)) %>% 
  summarise(n = sum(equal, na.rm = T))

data_3 %>% 
  tabyl(dif_day_dx_symp)

data_4 <- data_3 %>% 
  filter(!name_tx %in% c(NA, "Nothing")) %>% ### Delete registers without treatment's name 
  filter(!is.na(clean_date_of_rx_started)) %>% #### To make sure that the patients chosen everyone have date of treatment initiation
  filter(!is.na(clean_date_symptoms_started)) %>%  ### Date of symptom onset for each register selected
  filter(clean_date_symptoms_started != clean_date_of_rx_started) 
#### Don't have the same date of the treatment initiation regard to symptoms onset
  
 data_4 <- data_4 %>% 
  #filter(clean_date_of_close_out > "2014-01-01" & clean_date_of_close_out < "2025-01-01") %>% 
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
         treatment_failure = case_when(`treatment failure` == 1 ~ "Amphotericin B",
                                       `treatment failure` == 2 ~ "Amphotericin B liposomal",
                                       `treatment failure` == 3 ~ "Other",
                                       `treatment failure` == 4 ~ "N/A",
                                       TRUE ~ NA),
         tx_failure = if_else(!is.na(`treatment failure`), 1, 0),
         tx_failure_2 = if_else(!is.na(`treatment failure`), TRUE, FALSE),
         death_any_causes = if_else(death == 1, T, F),
         death_vl_2 = if_else(death_vl == 1, T, F),
         evolution = case_when(evolution_of_cases == 1 ~ "Cure",
                               evolution_of_cases == 2 ~ "Neglect",
                               evolution_of_cases == 4 ~ "Death for any-cause",
                               evolution_of_cases == 3 ~ "Death for VL",
                               evolution_of_cases == 5 ~ "Transferred",
                               # is.na(evolution_of_cases) ~ "Unknown",
                               TRUE ~ NA)) %>% 
  filter(!dif_day_tx_symp < 0)

median(data_4$dif_day_tx_symp) ### Median between symptoms onset to treatment 
median(data_4$dif_day_tx_symp[death_vl = 1]) ### Median of distribution symptoms to treatment by VL death

table(data_4$evolution_of_cases)

ggplot(data = data_4, aes(x = as.numeric(dif_day_tx_symp))) +
  geom_histogram(na.rm = T) +
  geom_vline(xintercept = median(as.numeric(data_4$dif_day_tx_symp)), linetype = "dashed") +
  geom_vline(xintercept = median(as.numeric(data_4$dif_day_tx_symp[death_vl = 1])), linetype = "dashed", color = "brown2") +
  scale_x_continuous(transform = "log10") +
  facet_wrap(~death_vl == 0) ### The difference between death by VL and cure were 39:29


### Rename variables and recode levels
data_4 <- data_4 %>% 
  mutate(education_level = factor(education_level, levels = c("Incomplete elementary school",
                                                              "Complete elementary and incomplete middle school",
                                                              "Complete middle school and incomplete high school",
                                                              "Complete high school",
                                                              "Incomplete university",
                                                              "Complete university",
                                                              "N/A")),
         year = year(clean_date_symptoms_started), 
         hiv_2 = case_when(hiv == 1 ~ "Yes",
                           hiv == 2 ~ "No",
                           hiv == 9 ~ "Unknown",
                           TRUE ~ NA),
         hiv_2 = factor(hiv_2, levels = c("Yes", "No", "Unknown"))) %>% 
  mutate_at(vars(dif_day_tx_symp, dif_day_dx_symp), as.numeric)

##############################################################################################################################
### Transforming variables class date to number. 
data_5 <- data_4 %>%
  mutate(dif_day_tx_symp = as.numeric(dif_day_tx_symp), ### N° days between symptoms to treatment initiation
         dif_day_dx_symp = as.numeric(dif_day_dx_symp), ### N° days between symptoms to diagnosis 
         dif_day_tx_dx = as.numeric(dif_day_tx_dx)) ### N° days between diagnosis to treatment initiation


### Distribution of delay to treatment initiation among federative states
boxplot(data_5$dif_day_tx_symp ~ data_5$federative_unit)
abline(h=c(14,30,60,90,120,180), col="red")


### Identify 99% distribution
quantile(as.numeric(data_4$dif_day_tx_symp), probs = 0.99, na.rm = TRUE)

### adjusting database using 99% distribution
data_5 <- data_4 %>% 
  mutate(code_id = paste(clean_date_of_birth, gender, code_of_the_municipality_of_residence, regional_health_area, 
                         race, municipality_code_where_a_person_get_infected, code_of_the_municipality, 
                         code_of_the_health_unit, health_region_of_residence)) %>% 
  #distinct(clean_date_of_birth, clean_date_of_notification, clean_date_of_investigation, 
   #        gender, age_as_pt_repored, .keep_all = T) %>%   ### Duplicated registers using unique code
  filter(dif_day_tx_symp < quantile(as.numeric(data_4$dif_day_tx_symp), probs = 0.99, na.rm = TRUE)) 
  #filter(!is.na(evolution_of_cases)) #### Only it will be used for the model

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
         day_of_outcome = if_else(is.na(date_of_death), date_of_close_out, date_of_death),
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
rio::export(data_6, here("data base and clean", "data_6.rda"))
rio::export(data_6, here("data base and clean", "data_6.csv"))

##################################################### Relapse patient #####################################################################
data_5 <- data_5 %>% 
  arrange(code_id, clean_date_of_notification) %>% 
  group_by(code_id) %>% 
  mutate(consult_num = dense_rank(clean_date_of_notification),
         criterial_case_previous = lag(identification_case),
         cure_case_previous = lag(cure),
         day_last_treat = lag(clean_date_of_rx_started),
         day_last_not = lag(clean_date_of_notification),
         dif_day_last_treat = as.Date(clean_date_of_rx_started) - as.Date(day_last_treat),
         dif_day_last_not = as.Date(clean_date_of_notification) - as.Date(day_last_not),
         dif_treat_start_symp = as.Date(clean_date_of_rx_started) - as.Date(clean_date_symptoms_started),
         dif_not_start_symp = as.Date(clean_date_of_notification) - as.Date(clean_date_symptoms_started),
         prev_death = lag(evolution),
         relapse = if_else(criterial_case_previous == "New case" & cure_case_previous == 1 & !type_of_cases_at_baseline == 3 &
                             (as.numeric(dif_day_last_treat) > 28 | as.numeric(dif_day_last_not) > 28) & 
                             (as.numeric(dif_treat_start_symp) > 0 | as.numeric(dif_not_start_symp) > 0) & 
                             !str_detect(prev_death, "Death"), 1, 0)) %>% 
  ungroup() 

data_5 <- data_5 %>% 
  mutate(time_y = round(dif_day_last_not / 365.25, digits = 1),
         time_y = case_when(is.na(dif_day_last_not) & !is.na(date_of_death) ~ round(as.numeric(as.Date("2023-12-31") - as.Date(date_of_death)) / 365.25, digits = 1),
                            is.na(dif_day_last_not) & is.na(date_of_death) ~ round(as.numeric(as.Date("2023-12-31") - as.Date(clean_date_of_notification)) / 365.25, digits = 1),
                            !is.na(dif_day_last_not) ~ as.numeric(time_y),
                            TRUE ~ NA_real_),
         relapse = replace_na(relapse, 0), 
         death_vl = replace_na(death_vl, 0),
         age_group_2 = case_when(age_years >= 0 & age_years <= 1 ~ "0-1",
                                 age_years >= 1 & age_years <= 2 ~ "1-2",
                                 age_years >= 2 & age_years <= 3 ~ "2-3",
                                 age_years >= 3 & age_years <= 5 ~ "3-5",
                                 age_years >= 5 & age_years <= 15 ~ "5-15",
                                 age_years >= 15 & age_years <= 20 ~ "15-20",
                                 age_years >= 20 & age_years <= 40 ~ "20-40",
                                 age_years >= 40 & age_years <= 60 ~ "40-60",
                                 age_years > 60  ~ ">60",
                                 TRUE ~ NA), 
         age_group_2 = factor(age_group_2, levels = c("0-1", "1-2", "2-3", "3-5", "5-15",
                                                      "15-20", "20-40", "40-60",
                                                      ">60")))

quantile(data_4$age_years, 0.99)

data_rev <- data_4 %>% 
  select(59, 60, 63:65, 57, year, 91:98, 87, 88)

data_4 %>% 
  tabyl(relapse, type_of_cases_at_baseline) %>% 
  adorn_totals("col")
#######################################################################################################################################################
### Database to work
rio::export(data_5, here("data base and clean", "data_4.rda"))

####################################################################################################################################################
### Table 1 characteristics and generalities 
data_4 <- import(here("data base and clean", "data_4.rda"))

theme_gtsummary_journal("lancet")

data_4 %>% 
  mutate(dif_day_tx_symp = as.numeric(dif_day_tx_symp),
         relapse_2 = if_else(relapse == 1, T, F)) %>%
  select(age_years, gender, education_level, race_2, pregnant, identification_case, 
         hiv_2, dif_day_tx_symp, name_tx, treatment_failure, relapse_2, evolution) %>% 
  tbl_summary(by = evolution, 
              #missing_stat = "{n} ({p}%)",
              missing = "no",
              statistic = list(all_categorical() ~ "{n} ({p}%)",
                               all_continuous() ~ "{median} ({p25},{p75})"), 
              label = list(age_years = "Age years",
                           # age_group = "Age group",
                           gender = "Gender",
                           education_level = "Education level",
                           race_2 = "Race/Color",
                           pregnant = "Gestational age",
                           identification_case = "Case identification",
                           dif_day_tx_symp = "Difference days treatment and onset symptoms",                           
                           name_tx = "Treatment",
                           hiv_2 = "HIV status",
                           treatment_failure = "Treatment failure",
                           relapse_2 = "Relapse")) %>% 
  add_overall(last = T) %>% 
  bold_labels() %>% 
  as_flex_table() %>% 
  autofit()
#save_as_docx(path = here("figures and tables", "table_1.docx"))

################################# Endemic channel to visceral leishmaniasis in Brazil last 5 years ############################################

### database with cases in Brazil for each year. 
year <- c(2013:2023) ### years was analyzed in the endemic channel

### Example of endemic channel using one year
ec_vl_brazil <- data_4 %>% 
  mutate(fis = as.Date(clean_date_symptoms_started)) %>%
  # group_by(fis) %>% 
  #  summarise(cases = n()) %>% 
  #  complete(fis = seq.Date(as.Date("2007-01-01"), as.Date("2023-12-01"), by = "1 day")) %>% 
  #  mutate(cases = replace_na(cases, 0)) %>% 
  mutate(ew = floor_date(fis, unit = "week"),
         ew_2 = epiweek(fis),
         year = year(fis)) %>% 
  filter(!ew_2 == 53)

historic <- ec_vl_brazil %>% 
  filter(year >= 2008 & year < 2013)

observated <- ec_vl_brazil %>%
  filter(year == 2013) %>% 
  group_by(ew_2) %>%
  summarise(cases = n(), .groups = "drop") %>% 
  filter(!ew_2 == 53) 

incidence_historic <- incidence(historic$fis, interval = "1 epiweek")
incidence_observated <- observated$cases

canal <- epiCo::endemic_channel(
  incidence_historic = incidence_historic,
  observations = incidence_observated,
  method = "geometric",
  ci = 0.95,
  plot = TRUE
)


# I make the incidence object for historic data
endemic_channel <- function(year) {
  results <- list()  # List of results. 
  
  for (year_eval in year) {
    cat("➤ Year observed:", year_eval, "\n")
    
    plots_list <- list()
    
    # History from last 5 years
    historic <- ec_vl_brazil %>% filter(year >= (year - 5) & year < year_eval)
    observated <- ec_vl_brazil %>% filter(year == year_eval)
    
    # Crear serie histórica de incidencia
    incidence_historic <- incidence(historic$fis, interval = "1 epiweek")
    
    # Agregar casos observados por semana
    obs_week <- observated %>%
      group_by(ew_2) %>%
      summarise(cases = n(), .groups = "drop") 
    
    obs_week_2 <- obs_week$cases
    
    
    # Canal endémico
    chanell <- epiCo::endemic_channel(
      incidence_historic = incidence_historic,
      observations = obs_week_2,
      method = "geometric",
      ci = 0.975,
      plot = TRUE
    )
    
    # Crear gráfico base
    p <- chanell$plot +
      labs(
        title = paste(year),
        subtitle = paste("Observated year:", year_eval),
        x = "Epidemiological week", y = "Cases"
      ) +
      theme_minimal(base_size = 10)
    
    # Guardar resultado
    plots_list[[as.character(year_eval)]] <- p
    results[[as.character(year_eval)]] <- chanell
  }
  
  return(results)
}

ec_vl <- endemic_channel(year = year)

#### I make the epidemiological calendar for using within endemic channel 
ec_2013 <- as.data.frame(epiCo::epi_calendar(year = 2013), optional = T)
ec_2013$week <- seq_len(length.out = nrow(ec_2013))
ec_2013$year <- rep(2013, length(nrow(ec_2013)))
ec_2014 <- as.data.frame(epiCo::epi_calendar(year = 2014), optional = T)
ec_2014$week <- seq_len(length.out = nrow(ec_2014))
ec_2014$year <- rep(2014, length(nrow(ec_2014)))
ec_2015 <- as.data.frame(epiCo::epi_calendar(year = 2015), optional = T)
ec_2015$week <- seq_len(length.out = nrow(ec_2015))
ec_2015$year <- rep(2015, length(nrow(ec_2015)))
ec_2016 <- as.data.frame(epiCo::epi_calendar(year = 2016), optional = T)
ec_2016$week <- seq_len(length.out = nrow(ec_2016))
ec_2016$year <- rep(2016, length(nrow(ec_2016)))
ec_2017 <- as.data.frame(epiCo::epi_calendar(year = 2017), optional = T)
ec_2017$week <- seq_len(length.out = nrow(ec_2017))
ec_2017$year <- rep(2017, length(nrow(ec_2017)))
ec_2018 <- as.data.frame(epiCo::epi_calendar(year = 2018), optional = T)
ec_2018$week <- seq_len(length.out = nrow(ec_2018))
ec_2018$year <- rep(2018, length(nrow(ec_2018)))
ec_2019 <- as.data.frame(epiCo::epi_calendar(year = 2019), optional = T)
ec_2019$week <- seq_len(length.out = nrow(ec_2019))
ec_2019$year <- rep(2019, length(nrow(ec_2019)))
ec_2020 <- as.data.frame(epiCo::epi_calendar(year = 2020), optional = T)
ec_2020$week <- seq_len(length.out = nrow(ec_2020))
ec_2020$year <- rep(2020, length(nrow(ec_2020)))
ec_2021 <- as.data.frame(epiCo::epi_calendar(year = 2021), optional = T)
ec_2021$week <- seq_len(length.out = nrow(ec_2021))
ec_2021$year <- rep(2021, length(nrow(ec_2021)))
ec_2022 <- as.data.frame(epiCo::epi_calendar(year = 2022), optional = T)
ec_2022$week <- seq_len(length.out = nrow(ec_2022))
ec_2022$year <- rep(2022, length(nrow(ec_2022)))
ec_2023 <- as.data.frame(epiCo::epi_calendar(year = 2023), optional = T)
ec_2023$week <- seq_len(length.out = nrow(ec_2023))
ec_2023$year <- rep(2023, length(nrow(ec_2023)))

### Join endemic calendar for each year
ec <- as.data.frame(rbind(ec_2013, ec_2014, ec_2015, ec_2016, ec_2017, ec_2018, ec_2019,
                          ec_2020, ec_2021, ec_2022, ec_2023))


ec <- ec %>%  
  rename("date" = c(1)) %>% 
  mutate(year_week = paste(year, week, sep = "-")) %>% 
  filter(week < 53)

ec_vl_2013 <- ec_vl$'2013'$data
ec_vl_2014 <- ec_vl$'2014'$data
ec_vl_2015 <- ec_vl$'2015'$data
ec_vl_2016 <- ec_vl$'2016'$data
ec_vl_2017 <- ec_vl$'2017'$data
ec_vl_2018 <- ec_vl$'2018'$data
ec_vl_2019 <- ec_vl$'2019'$data
ec_vl_2020 <- ec_vl$'2020'$data
ec_vl_2021 <- ec_vl$'2021'$data
ec_vl_2022 <- ec_vl$'2022'$data
ec_vl_2023 <- ec_vl$'2023'$data


ec_vl <- rbind(ec_vl_2013, ec_vl_2014, ec_vl_2015, ec_vl_2016, ec_vl_2017, ec_vl_2018,
               ec_vl_2019, ec_vl_2020, ec_vl_2021, ec_vl_2022, ec_vl_2023)

ec_vl$week <- rep(seq(from = 1, to = 52, by = 1), times = 11)
ec_vl$year <- rep(2013:2023, each = 52) 

ec_vl_2 <- ec_vl %>% 
  mutate(year_week = paste(year, week, sep = "-"))

ec_vl_2 <- left_join(ec_vl_2, select(ec, date, year_week), by = "year_week")

ec_vl_3 <- ec_vl_2 %>% 
  mutate(month = floor_date(date, unit = "month")) %>% 
  group_by(month) %>% 
  summarise(obs = sum(obs), 
            low_lim = sum(low_lim),
            central = sum(central),
            up_lim = sum(up_lim), .groups = "drop")

rio::export(ec_vl_2, here("data base and clean", "database_endemic_channel.rda")) #### Export database for endemic channel

ec_vl_b <- ggplot(data = ec_vl_2, mapping = aes(x = date, y = obs)) +
  geom_line(aes(y = low_lim, color = "Low limit")) +
  geom_line(aes(y = central, color = "Central")) +
  geom_line(aes(y = up_lim, color = "Upper limit")) +
  geom_ribbon(aes(ymin = central, ymax = up_lim), alpha = 0.3, fill = "darkred") +
  geom_ribbon(aes(ymin = low_lim, ymax = central), alpha = 0.3, fill = "orange") +
  geom_ribbon(aes(ymin = min(low_lim), ymax = low_lim), alpha = 0.3, fill = "darkgreen") +
  geom_point() +
  geom_line() +
  labs(x = "Date", y = "Visceral Leishmaniasis Cases") + 
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, size = 12, color = "black"), 
    legend.position = "no",
    axis.title = element_text(face = "bold", size = 14, color = "black"),
  ) +
  scale_x_date(date_breaks = "6 months", date_labels = "%b/%y") +
  scale_y_continuous(n.breaks = 10, limits = c(0, 130)) +
  scale_color_manual(values = c("Low limit" = "darkgreen", "Upper limit" = "darkred", "Central" = "orange")) 

print(ec_vl_b)


ec_vl_b_2 <- ggplot(data = ec_vl_3 %>% filter(month >= "2013-01-01"), 
                    mapping = aes(x = month, y = obs)) +
  geom_line(aes(y = low_lim, color = "Safety Zone")) +
  geom_line(aes(y = central, color = "Alert Zone")) +
  geom_line(aes(y = up_lim, color = "Alarm Zone")) +
  geom_ribbon(aes(ymin = central, ymax = up_lim), alpha = 0.3, fill = "darkred") +
  geom_ribbon(aes(ymin = low_lim, ymax = central), alpha = 0.3, fill = "orange") +
  geom_ribbon(aes(ymin = 0, ymax = low_lim), alpha = 0.3, fill = "darkgreen") +
  geom_point() +
  geom_line() +
  labs(x = "Date", y = "Visceral Leishmaniasis Cases",
       title = "Endemic channel Visceral Leishmaniasis in Brazil, 2013-2023", 
       color = "Zones") + 
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, size = 12, color = "black"), 
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 14, color = "black"),
    legend.text = element_text(size = 12, color = "black"),
    axis.title = element_text(face = "bold", size = 14, color = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold", color = "black", size = 16),
  ) +
  scale_x_date(date_breaks = "3 months", date_labels = "%b/%y") +
  scale_y_continuous(n.breaks = 30, limits = c(0, 450)) +
  scale_color_manual(values = c("Safety Zone" = "darkgreen", "Alarm Zone" = "darkred", "Alert Zone" = "orange")) 

ec_vl_b_2 

ggplotly(ec_vl_b_2)

##################################################################################################################################
## Import work's data base 
data_6 <- import(here("data base and clean", "data_6.rda"))

data_6 <- data_6 %>% 
  mutate(relapse_2 = if_else(relapse == 1, "Yes", "No"),
         dif_day_last_treat_2 = if_else(is.na(dif_day_last_treat), as.numeric(dif_treat_outcome), 
                                        as.numeric(dif_day_last_treat))) %>% 
  mutate(complication = if_else(hemorrhagic_phenomena == 1 | liver_enlargement == 1 |
                                  spleen_enlargement == 1| edema == 1, "Yes", "No")) %>% 
  filter(dif_day_tx_symp > 14) ### Attention with it, filter after fourteen days onset symptoms to treatment. 



### Distributiona de federative units according with the media between symptoms onset to treatmnet. 

map_uf_shp <- read_sf(here("BR_UF_2023", "BR_UF_2023.shp")) %>% 
  mutate(CD_UF = as.integer(CD_UF))

data_7 <- left_join(data_6, select(map_uf_shp, CD_UF, NM_UF), 
                    by = c("federative_unit" = "CD_UF"))

tabla_fu <- data_7 %>% 
  mutate(year_2 = as.Date(year, format = "%Y")) %>% 
  group_by(NM_UF, year) %>% 
  summarise(median = median(dif_day_tx_symp),
            cases = n(), .groups = "drop") 


ggplot(data = tabla_fu, aes(x = year, y = median)) +
  geom_line(aes(color = NM_UF))


ggplot(data = tabla_fu, aes(x = year, y = median)) +
  geom_boxplot(aes(group = NM_UF)) +
  geom_jitter(aes(colour = NM_UF)) +
  theme_classic() +
  theme(legend.position = "bottom") +
  scale_x_date(date_breaks = "year", date_labels = "%Y")



library(gganimate)
ggplot(data = tabla_fu, aes(x = cases, y = median, size = median, color = NM_UF)) +
  geom_point() +
  theme_classic() +
  theme(legend.position = "bottom") +
  transition_time(year) +
  labs(title = 'Year: {frame_time}', x = 'Cases', y = 'Median') 

anim_save(
  filename = "271-ggplot2-animated-gif-chart-with-gganimate1.gif",
  animation = last_animation(),
  width = 1024,     # width in pixels
  height = 800,    # height in pixels
  res = 600        # DPI (dots per inch)
)

#################################################################################################################################################

### Kaplan Meier analysis: Lethality risk regarding treatment / relapse
### Lethally comparing time between onset symptoms and start tretment 

## Survival analysis: Evaluating relapse and lethally

### Behaviour lethality according with delay symptom onset to tretament

lethality_risk <- survfit(Surv(dif_treat_outcome, death_vl) ~ 1, data = data_6) 

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

lethaly_risk_symp_tx <- survfit(Surv(dif_treat_outcome, death_vl) ~ dif_treat_start_symp_cat , data = data_6
                                #  filter(str_detect(name_tx, "Antimonial|Amphote|Lipos"))
)

gstx <- ggsurvplot(fit = lethaly_risk_symp_tx,
                   data = data_6, 
                   pval = T, pval.coord = c(330, .95),
                   #ggtheme = theme_bw(),
                   #risk.table = T,
                   tables.theme = theme_cleantable(),
                   ylim = c(0.85, 1), 
                   font.x = c(15, "bold"), font.y = c(15, "bold"),
                   xlab = "days",
                   risk.table = T,
                   legend = c("bottom"), 
                   legend.title = "Delay",
                   legend.labs = c("<= 30 days", "30–89 days", "90–119 days", ">= 120 days"),
                   legend.nrow = 2) 

gstx <- gstx$plot + theme_bw() +
  theme(legend.title = element_text(face = "bold", size = 14),
        legend.text =  element_text(size = 12), 
        axis.title = element_text(face = "bold"), 
        legend.position = "bottom") +
  
  scale_x_continuous(n.breaks = 15)

print(gstx)


### Death comparing time between onset symptoms and start tretment  and hiv 
death_risk_day_hiv <- survfit(Surv(dif_treat_outcome, death_vl) ~ dif_treat_start_symp_cat + hiv_2 , 
                              data = data_6)

gstx_hiv <- ggsurvplot_facet(fit = death_risk_day_hiv,
                             data = data_6,
                             facet.by = "dif_treat_start_symp_cat",
                             palette = "lancet",
                             pval = T, pval.coord = c(330, .95),
                             #ggtheme = theme_bw(),
                             risk.table = T,
                             ylim = c(0.6, 1),
                             font.x = c(15, "bold"), font.y = c(15, "bold"),
                             xlab = "days",
                             legend.title = "HIV", 
                             legend = "bottom",
                             font.legend = list(size = 14, color = "black", face = "bold"),
                             legend.labs = c("No", "Yes", "Unknown"),  
                             panel.labs.font = list(size = 12, face = "bold"),
                             short.panel.labs = T) +
  scale_x_continuous(n.breaks = 15)


gstx_hiv

data_6 %>% 
  tabyl(hiv_2, dif_treat_start_symp_cat)

### Death comparing time between symptoms onset to start treatment and and type of the treatment

data_6 <- data_6 %>% 
  mutate(edema_2 = if_else(edema == 1, "Yes oedema", "No oedema"),
         edema_2 = factor(edema_2, levels = c("Yes oedema", "No oedema")), 
         severe = if_else(edema == 1 | hemorrhagic_phenomena == 1, "Severe", "No severe"),
         severe = factor(severe, levels = c("Severe", "No severe")))

centr = theme_bw() + theme(plot.title = element_text(hjust = 0.5, face = "bold"))

death_risk_day_tx <- survfit(Surv(dif_treat_outcome, death_vl) ~ dif_treat_start_symp_cat + name_tx + severe, 
                             data = data_6 %>% 
                               filter(str_detect(name_tx, "Antimonial|Amphote|Lipos|Pentamidine"))) 

death_risk_day_tx_2 <- survfit(Surv(dif_treat_outcome, death_vl) ~ dif_treat_start_symp_cat + name_tx + hiv_2, 
                               data = data_6%>% 
                                 filter(str_detect(name_tx, "Antimonial|Amphote|Lipos|Pentamidine")) %>% 
                                 filter(dif_day_tx_symp > 14)) 

data_6 %>% 
  tabyl(name_tx)

### ### Death comparing time between symptoms onset to start treatment and and type of the treatment based on complications
gstx_tx <- ggsurvplot_facet(fit = death_risk_day_tx,
                            data = data_6,
                            facet.by = c("dif_treat_start_symp_cat", "severe"),
                            palette = "lancet",
                            pval = T, pval.coord = c(330, .99),
                            #title = "Time to death according with the delay symptoms onset to treatment and 
                            #type of treatment used and presence of complications",
                            font.title = list(size = 14, face ="bold"), 
                            #ggtheme=centr,
                            ggtheme = theme_bw(),
                            #risk.table = T,
                            ylim = c(0.4, 1),
                            font.x = c(15, "bold"), font.y = c(15, "bold"),
                            xlab = "Time (days after to start treatment)",
                            legend.title = "Treatment", 
                            legend = "bottom",
                            font.legend = list(size = 14, color = "black", face = "bold"),
                            legend.labs = c("Amphotericin B deoxycholate",
                                            "Antimonial Pentavalent",
                                            "Liposomal amphotericin B",
                                            "Pentamidine"),  
                            panel.labs.font = list(size = 12, face = "bold"),
                            short.panel.labs = T) +
  scale_x_continuous(n.breaks = 15)


gstx_tx


gstx_tx <- ggsurvplot_facet(fit = death_risk_day_tx_2,
                            data = data_6,
                            facet.by = c("dif_treat_start_symp_cat", "hiv_2"),
                            palette = "lancet",
                            pval = T, pval.coord = c(330, .99),
                            title = "Time to death according with the delay symptoms onset to treatment and 
                            type of treatment used and presence of complications",
                            font.title = list(size = 14, face ="bold"),
                            ggtheme=centr,
                            #ggtheme = theme_bw(),
                            #risk.table = T,
                            ylim = c(0.6, 1),
                            font.x = c(15, "bold"), font.y = c(15, "bold"),
                            xlab = "days",
                            legend.title = "Treatment", 
                            legend = "bottom",
                            font.legend = list(size = 14, color = "black", face = "bold"),
                            legend.labs = c("Amphotericin B deoxycholate",
                                            "Antimonial Pentavalent",
                                            "Liposomal amphotericin B",
                                            "Pentamidine"),  
                            panel.labs.font = list(size = 12, face = "bold"),
                            short.panel.labs = T) +
  scale_x_continuous(n.breaks = 15)


gstx_tx

data_6 %>% 
  filter(dif_day_tx_symp > 14) %>% 
  group_by(name_tx, dif_treat_start_symp_cat) %>% 
  summarise(n = n()) %>% 
  ggplot(aes(x = name_tx, y = n)) +
  geom_col(fill = "skyblue2", colour = "black") +
  labs(y = "Cases", x = "Treatment") +
  facet_grid(~dif_treat_start_symp_cat) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90),
    axis.title = element_text(color = "black", size = 14, face = "bold"),
    axis.text = element_text(color = "black", size = 12),
    strip.text = element_text(color = "black",size = 12, face = "bold")
  ) 

## Time to relapse according with delay treatment from symptoms onset


death_risk_realpse <- survfit(Surv(dif_treat_outcome, relapse) ~ dif_treat_start_symp_cat, 
                              data = data_6 %>% filter(dif_day_last_treat_2 < 366))

centr = theme_bw() + theme(plot.title = element_text(hjust = 0.5, face = "bold"))

gstx_relapse <- ggsurvplot(fit = death_risk_realpse,
                           data = data_6,
                           pval = T, pval.coord = c(330, .98), 
                           legend.labs = c("<= 30 days", "30–89 days", "90–119 days", ">= 120 days"),
                           title = "Time to relapse according with delay symptoms onset to treatment",
                           font.title = list(size = 14, face ="bold"),
                           ggtheme=centr,
                           legend.title = "Delay symptoms onset to treatment",
                           legend = "bottom",
                           font.legend = list(size = 14, color = "black", face = "bold"),
                           #ggtheme = theme_bw(),
                           #risk.table = T,
                           ylim = c(0.94, 1),
                           font.x = c(15, "bold"), font.y = c(15, "bold"),
                           xlab = "days")


gstx_relapse                      


## Time to relapse according with delay symptoms onset to treatment and treatment
data_6_tx <- data_6 %>% 
  filter((str_detect(name_tx, "Antimonial|Amphote|Lipos")))

death_risk_realpse <- survfit(Surv(dif_treat_outcome, relapse) ~ dif_treat_start_symp_cat + name_tx, 
                              data = data_6_tx %>% filter(dif_day_last_treat_2 < 366 ))


gstx_relapse <- ggsurvplot_facet(fit = death_risk_realpse,
                                 data = data_6_tx,
                                 facet.by = "name_tx",
                                 short.panel.labs = T,
                                 pval = T, pval.coord = c(330, .98),
                                 legend.labs = c("<= 30 days", "30–89 days", "90–119 days", ">= 120 days"),
                                 legend.title = "Delay symptoms onset to treatment",
                                 panel.labs.font = list(size = 12, face = "bold"),
                                 legend = "bottom",
                                 font.legend = list(size = 14, color = "black", face = "bold"),
                                 #ggtheme = theme_bw(),
                                 #risk.table = T,
                                 ylim = c(0.85, 1),
                                 font.x = c(15, "bold"), font.y = c(15, "bold"),
                                 xlab = "days") +
  scale_x_continuous(n.breaks = 12)


gstx_relapse    



#####
death_risk_realpse_2 <- survfit(Surv(dif_treat_outcome, relapse) ~  name_tx + dif_treat_start_symp_cat, 
                                data = data_6_tx %>% filter(dif_day_last_treat_2 < 366 ))

tabyl(data_6_tx$name_tx)

gstx_relapse_2 <- ggsurvplot_facet(fit = death_risk_realpse_2,
                                   data = data_6_tx,
                                   facet.by = "dif_treat_start_symp_cat",
                                   short.panel.labs = T,
                                   pval = T, pval.coord = c(330, .98),
                                   legend.labs = c("Amphotericin B deoxycholate",
                                                   "Antimonial Pentavalent",
                                                   "Liposomal amphotericin B"),
                                   legend.title = "Treatment",
                                   panel.labs.font = list(size = 12, face = "bold"),
                                   legend = "bottom",
                                   font.legend = list(size = 14, color = "black", face = "bold"),
                                   #ggtheme = theme_bw(),
                                   #risk.table = T,
                                   ylim = c(0.85, 1),
                                   font.x = c(15, "bold"), font.y = c(15, "bold"),
                                   xlab = "days") +
  scale_x_continuous(n.breaks = 12)


gstx_relapse_2


###############################################################################################################

data_6 %>% 
  tabyl(dif_day_dx_symp, dif_day_tx_dx)

min(data_6$dif_day_dx_symp, na.rm = T)
data_6 %>% 
  mutate(rev = if_else(dif_day_tx_dx < 0, 1, 0)) %>% 
  summarise(n = sum(rev, na.rm = T)) %>% 
  ungroup()

data_6 %>% 
  tabyl(evolution)
