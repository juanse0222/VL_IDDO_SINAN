pacman::p_load(flextable, tidyverse, janitor, epikit, gtsummary, rio, 
               survival, survminer, sf, here, incidence, ggstance, ggbeeswarm, 
               grid, ggsurvfit, geobr, rms)
##########################################################################################################################
### Descriptive analysis. 
data_6 <- import(here("data base and clean", "data_6.rda")) 

## Initial distribution evaluation about delay treatment 
ggplot(data = data_6, aes(x = dif_day_tx_symp)) +
  geom_histogram(colour = "black", fill = "lightblue") +
  geom_vline(aes(xintercept = median(dif_day_tx_symp[death_vl == 1], na.rm = T)), 
             linetype = "dashed", color = "brown2") + 
  geom_vline(aes(xintercept = median(dif_day_tx_symp[death_vl == 0], na.rm = T)), 
             linetype = "dashed", color = "black") + 
  scale_x_continuous(transform = "log", n.breaks = 7, 
                     labels = scales::number_format(accuracy = 2)) +
  labs(x = "Delay between symptom onset to treatment", y = "Cases") +
  scale_y_continuous(n.breaks = 10) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    axis.text = element_text(color = "black", size = 12),
    axis.title = element_text(color = "black", size = 14, face = "bold"),
    legend.title = element_text(face = "bold", size = 14, color = "black"),
    legend.text = element_text(color = "black", size = 12))


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




### Distribution and direction time between recurrences by characteristics and some variables chosen 

## variables in the table
var_order <- c("age_group_2", "gender", "year_period", "region",
               "education_level", "race_2", "hiv_2", "name_tx", "evolution")

## desired level orders WITH the reference as first element
level_orders <- list(
  age_group_2 = c( "16-50", "0-1", "1-5", "6-15", ">50"),
  gender      = c("Female", "Male"),          # adjust if your coding is "F"/"M"
  year_period = c("2007-2009", "2010-2014", "2015-2018", "2019-2023"),
  region      = c("Center-West", "North", "Northeast", "Southeast", "South"),
  education_level = c(
    "Incomplete elementary school",
    "Complete elementary and incomplete middle school",
    "Complete middle school and incomplete high school",
    "Complete high school",
    "Incomplete university",
    "Complete university",
    "N/A"),
  race_2 = c("Brown", "White", "Afro-descendants", "Yellow", "Indigenous", "Unknown"),
  hiv_2  = c("No", "Yes", "Unknown"),
  name_tx = c(
    "Antimonial Pentavalent",
    "Amphotericin B deoxycholate",
    "Liposomal amphotericin B",
    "Others")
  #evolution = c("Cure", "Death", "Unknown")   # adjust to your actual labels
)


data_6 %>% 
  tabyl(dif_day_tx_symp)

### Function to make the table 1 - Baseline characteristics based on delay treatment initiation
median_time_glm_ordered <- function(vars, level_orders = list(), data = data_6) {
  map_dfr(vars, function(v) {
    
    # subset: valid x and dif_day_tx_symp > 0 for log
    d_sub <- data %>%
      filter(!is.na(.data[[v]]),
             !is.na(dif_day_tx_symp),
             dif_day_tx_symp > 0) %>% 
      mutate(log_delay = log(dif_day_tx_symp))
    
    if (nrow(d_sub) == 0) return(NULL)
    
    # enforce factor levels according to level_orders (this sets the reference)
    if (!is.null(level_orders[[v]])) {
      d_sub[[v]] <- factor(d_sub[[v]], levels = level_orders[[v]])
    } else {
      d_sub[[v]] <- factor(d_sub[[v]])
    }
    
    # summary: n, median, range
    df_sum <- d_sub %>%
      group_by(.data[[v]]) %>%
      summarise(
        n      = n(),
        median = paste0(round(median(dif_day_tx_symp, na.rm = TRUE), 0)," ", "[",
                        round(quantile(dif_day_tx_symp, 0.25, na.rm = TRUE), 0), "-",
                        round(quantile(dif_day_tx_symp, 0.75, na.rm = TRUE), 0), "]"),
        range  = paste0(min(dif_day_tx_symp, na.rm = TRUE), "-",
                        max(dif_day_tx_symp, na.rm = TRUE)),
        .groups = "drop"
      ) %>%
      mutate(
        variable = v,
        level    = as.character(.data[[v]])
      ) %>%
      select(variable, level, n, median, range)
    
    # GLM: log(dif_day_tx_symp) ~ v
    form <- as.formula(paste("log_delay ~", v))
    fit  <- glm(form, data = d_sub, family = gaussian())
    
    sm <- summary(fit)$coefficients
    coef_df <- as.data.frame(sm)
    coef_df$term <- rownames(coef_df)
    
    # reference = first level (given by level_orders)
    base_level <- levels(d_sub[[v]])[1]
    
    base_row <- data.frame(
      level    = base_level,
      gm_ratio = 1,
      gm_lower = 1,
      gm_upper = 1
    )
    
    others <- subset(coef_df, term != "(Intercept)")
    
    if (nrow(others) > 0) {
      others_df <- data.frame(
        level    = sub(paste0("^", v), "", others$term),
        gm_ratio = round(exp(others$Estimate), 3),
        gm_lower = round(exp(others$Estimate - 1.96 * others$`Std. Error`), 3),
        gm_upper = round(exp(others$Estimate + 1.96 * others$`Std. Error`), 3)
      )
      glm_tbl <- rbind(base_row, others_df)
    } else {
      glm_tbl <- base_row
 }
    
    # join summary + GLM CI
    df <- df_sum %>%
      left_join(glm_tbl, by = "level")
    
    # order levels for display too
    if (!is.null(level_orders[[v]])) {
      df <- df %>% mutate(level = factor(level, levels = level_orders[[v]]))
    } else {
      df <- df %>% mutate(level = factor(level))
    }
    
    df
  }) %>%
    mutate(variable = factor(variable, levels = var_order)) %>%
    arrange(variable, level)
}


vars <- var_order

dif_table <- median_time_glm_ordered(vars, level_orders, data = data_6)

dif_table <- dif_table %>% 
  mutate(
    is_ref = gm_ratio == 1 & gm_lower == 1 & gm_upper == 1,
    gm_ratio_chr = ifelse(is_ref, "Reference",
                          sprintf("%.2f", gm_ratio)),
    gm_lower_chr = ifelse(is_ref, "-",
                          sprintf("%.2f", gm_lower)),
    gm_upper_chr = ifelse(is_ref, "-",
                          sprintf("%.2f", gm_upper))
  ) %>% 
  mutate(gm_ratio_chr = paste0(gm_ratio_chr, " ", "[", gm_lower_chr,
                           "-", gm_upper_chr, "]")) %>% 
  select(-gm_lower_chr, -gm_upper_chr, -gm_lower, -gm_upper, 
         -is_ref, -gm_ratio) 

border_line = officer::fp_border(color = "black", width = 1.38) 
flextable(dif_table)

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
      name_tx         = "Treatment",
      evolution       = "Evolution"
    )
  ) %>% 
  rename(
    "Characteristic"          = variable,
    "Level"                   = level,
    "N"                       = n,
    "Median time (days)"      = median,
    "Range (days)"            = range,
    "GLM (geometric mean), [95% CI]"    = gm_ratio_chr
  ) %>% 
  flextable() %>% 
  merge_at(j = 1, i = 1:5, part = "body") %>% 
  merge_at(j = 1, i = 6:7, part = "body") %>% 
  merge_at(j = 1, i = 8:12, part = "body") %>% 
  merge_at(j = 1, i = 13:19, part = "body") %>% 
  merge_at(j = 1, i = 20:25, part = "body") %>% 
  merge_at(j = 1, i = 26:28, part = "body") %>% 
  merge_at(j = 1, i = 29:31, part = "body") %>% 
  merge_at(j = 1, i = 32:36, part = "body") %>% 
  hline(i = c(5, 7, 12, 19, 25, 28, 31), part = "body", border = border_line) %>% 
  align(j = 3:6, part = "body", align = "center") %>% 
  align(j = 3:6, part = "header", align = "center") %>% 
  bold(j = 1, part = "all") %>% 
  bold(i = 1, part = "header") %>% 
  align(j = 3:5, part = "all", align = "center") %>% 
  autofit() 

tab_1 %>% 
  save_as_docx(path = here("figures and tables", "table_1.docx"))

######################################################## GAM and logsitic model for mortality evaluation ###################################################

var_order_2 <- c("age_group_2", "gender",  "region",
                 "education_level", "race_2", "hiv_2", "name_tx")

# keep only rows with non-missing outcome and delay
data_mort <- data_6 %>%
  filter(!is.na(death_vl),
         !is.na(dif_day_tx_symp))

# IMPORTANT: enforce factor levels so that the first level in level_orders
# is the reference used by glm()
for (v in var_order_2) {
  if (!is.null(level_orders[[v]])) {
    data_mort[[v]] <- factor(data_mort[[v]], levels = level_orders[[v]])
  } else {
    data_mort[[v]] <- factor(data_mort[[v]])
  }
}

# multivariable with spline
form_multi <- as.formula(
  paste(
    "death_vl ~ rcs(dif_day_tx_symp, 4) +",
    paste(var_order_2, collapse = " + "))
)

fit_multi <- glm(form_multi, data = data_mort, family = binomial())


get_or_from_fit <- function(fit, v, level_orders) {
  td <- broom::tidy(fit) %>%
    filter(term != "(Intercept)",
           !grepl("^rcs\\(", term)) %>%          # drop spline basis terms
    filter(grepl(paste0("^", v), term))          # keep only this variable
  
  if (nrow(td) == 0) {
    return(NULL)
  }
  
  td %>%
    mutate(
      variable = v,
      level    = sub(paste0("^", v), "", term),
      OR       = exp(estimate),
      lower    = exp(estimate - 1.96 * std.error),
      upper    = exp(estimate + 1.96 * std.error)
    ) %>%
    select(variable, level, OR, lower, upper)
}


death_vl ~ rcs(dif_day_tx_symp, 4) + v


# list to collect univariate results
univ_list <- map_dfr(var_order_2, function(v) {
  form_uni <- as.formula(paste("death_vl ~ rcs(dif_day_tx_symp, 4) +", v))
  
  fit_uni <- glm(form_uni, data = data_mort, family = binomial())
  
  get_or_from_fit(fit_uni, v, level_orders = NULL) %>%
    rename(
      OR_uni    = OR,
      lower_uni = lower,
      upper_uni = upper
    )
})


multi_list <- map_dfr(var_order_2, function(v) {
  get_or_from_fit(fit_multi, v, level_orders = NULL) %>%
    rename(
      OR_multi    = OR,
      lower_multi = lower,
      upper_multi = upper
    )
})

# assumes you still have `level_orders` defined as before
# (with the reference category as the first level for each variable)

# function to build full table for one variable
build_full_for_var <- function(v) {
  lvls <- level_orders[[v]]
  if (is.null(lvls)) {
    # fall back to levels from data if not in level_orders
    lvls <- levels(droplevels(data_mort[[v]]))
  }
  
  base <- tibble(
    variable = v,
    level    = lvls
  )
  
  u <- univ_list %>% filter(variable == v)
  m <- multi_list %>% filter(variable == v)
  
  df <- base %>%
    left_join(u, by = c("variable", "level")) %>%
    left_join(m, by = c("variable", "level"))
  
  # first level is reference
  df <- df %>%
    mutate(
      OR_uni    = ifelse(row_number() == 1 & is.na(OR_uni),    1, OR_uni),
      lower_uni = ifelse(row_number() == 1 & is.na(lower_uni), 1, lower_uni),
      upper_uni = ifelse(row_number() == 1 & is.na(upper_uni), 1, upper_uni),
      
      OR_multi    = ifelse(row_number() == 1 & is.na(OR_multi),    1, OR_multi),
      lower_multi = ifelse(row_number() == 1 & is.na(lower_multi), 1, lower_multi),
      upper_multi = ifelse(row_number() == 1 & is.na(upper_multi), 1, upper_multi)
    )
}

mort_table <- map_dfr(var_order_2, build_full_for_var)

print(mort_table, n = 100)

mort_table <- mort_table %>%
  mutate(
    # flags for reference in univariate and multivariable
    is_ref_uni   = OR_uni    == 1 & lower_uni == 1 & upper_uni == 1,
    is_ref_multi = OR_multi  == 1 & lower_multi == 1 & upper_multi == 1,
    
    OR_uni_chr    = ifelse(is_ref_uni,   "Reference", sprintf("%.2f", OR_uni)),
    lower_uni_chr = ifelse(is_ref_uni,   "-",         sprintf("%.2f", lower_uni)),
    upper_uni_chr = ifelse(is_ref_uni,   "-",         sprintf("%.2f", upper_uni)),
    
    OR_multi_chr    = ifelse(is_ref_multi, "Reference", sprintf("%.2f", OR_multi)),
    lower_multi_chr = ifelse(is_ref_multi, "-",         sprintf("%.2f", lower_multi)),
    upper_multi_chr = ifelse(is_ref_multi, "-",         sprintf("%.2f", upper_multi))
  ) %>%
  select(variable, level,
         OR_uni_chr, lower_uni_chr, upper_uni_chr,
         OR_multi_chr, lower_multi_chr, upper_multi_chr) %>%
  mutate(OR_uni_chr = paste0(OR_uni_chr, " ", "[", lower_uni_chr,
                             "-", upper_uni_chr, "]")) %>%
  mutate(OR_multi_chr = paste0(OR_multi_chr, " ", "[", lower_multi_chr,
                               "-", upper_multi_chr, "]")) %>% 
  select(-lower_uni_chr, -lower_multi_chr, -upper_uni_chr, -upper_multi_chr) %>% 
  rename(
    OR_uni       = OR_uni_chr,
    OR_multi     = OR_multi_chr,
  )


table_2 <- mort_table %>%
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
      name_tx         = "Treatment"
    )
  ) %>%
  rename(
    Characteristic         = variable,
    Level                  = level,
    "Univariable odds ratio [95% CI]" = OR_uni,
    "Multivariable odds ratio [95% CI]" = OR_multi
  ) %>%
  flextable() %>%
  merge_at(j = 1, i = 1:5, part = "body") %>% 
  merge_at(j = 1, i = 6:7, part = "body") %>% 
  merge_at(j = 1, i = 8:12, part = "body") %>% 
  merge_at(j = 1, i = 13:19, part = "body") %>% 
  merge_at(j = 1, i = 20:25, part = "body") %>% 
  merge_at(j = 1, i = 26:28, part = "body") %>% 
  merge_at(j = 1, i = 29:32, part = "body") %>% 
  hline(i = c(5, 7, 12, 19, 25, 28, 32), part = "body", border = border_line) %>% 
  align(j = 2:4, part = "body", align = "center") %>% 
  align(j = 2:4, part = "header", align = "center") %>% 
  bold(j = 1, part = "all") %>% 
  bold(i = 1, part = "header") %>% 
  align(j = 3:4, part = "all", align = "center") %>% 
  autofit() 

table_2 %>% 
  save_as_docx(path = here("figures and tables", "table_2.docx"))
###################################################################################################################################################
### Behavior difference start treatment and onset symptoms. 

q25 <- quantile(as.numeric(data_6$dif_treat_start_symp), 0.25, na.rm = TRUE)
q50 <- quantile(as.numeric(data_6$dif_treat_start_symp), 0.50, na.rm = TRUE)
q75 <- quantile(as.numeric(data_6$dif_treat_start_symp), 0.75, na.rm = TRUE)
q95 <- quantile(as.numeric(data_6$dif_treat_start_symp), 0.95, na.rm = TRUE)

quant_df <- data.frame(
  Quantile = c("Q25 = 12", "Q50 = 27", "Q75 = 60", "Q95 = 165"),
  value = c(q25, q50, q75, q95))

data_6 %>% 
  filter(!is.na(dif_day_tx_dx)) %>% 
  summarise(n = quantile(dif_day_tx_dx, 0.5))


### Table 1. General characteristics in database

data_6$year_period

data_6 %>% 
  filter(evolution %in% c("Death for VL", "Cure")) %>% 
  mutate(dif_day_tx_symp = as.numeric(dif_day_tx_symp)) %>%
  select(age_group_2, gender, race_2, year_period, region,
         hiv_2, dif_day_tx_symp, name_tx, evolution,
        ) %>% 
  tbl_summary(by = evolution, 
              #missing_stat = "{n} ({p}%)", 
              missing = "no",
              statistic = list(all_categorical() ~ "{n} ({p}%)",
                               all_continuous() ~ "{median} ({p25},{p75})"), 
              label = list(age_group_2 = "Age years",
                           # age_group = "Age group",
                           gender = "Gender",
                           race_2 = "Race/Color",
                           year_period = "Period",
                           region = "Region",
                           hiv_2 = "HIV",
                           dif_day_tx_symp = "Days symptoms to treatment",
                           name_tx = "Treatment"
                           )) %>% 
  add_overall(last = T) %>%
  add_p() %>% 
  bold_labels() %>% 
  as_flex_table() %>% 
  autofit() %>% 
  save_as_docx(path = here("figures and tables", "table_1.0.docx"))

data_6 %>% 
  tabyl(evolution)


data_6 %>% 
  mutate(dif_day_tx_symp = as.numeric(dif_day_tx_symp)) %>%
  select(age_years, gender, education_level, race_2, pregnant, evolution, dif_day_dx_symp, dif_day_tx_symp, 
         dif_treat_start_symp_cat) %>% 
  tbl_summary(by = dif_treat_start_symp_cat, 
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
                           dif_day_tx_symp = "Days symptoms to treatment",         
                           dif_day_dx_symp = "Days symptoms to diagnosis")) %>% 
  add_overall(last = T) %>% 
  add_p() %>% 
  bold_labels() %>% 
  as_flex_table() %>% 
  autofit()  %>% 
  save_as_docx(path = here("figures and tables", "table_1.docx"))

### Table 2. Clinic and treatment features in database by delay symptoms onset to treatment, 

data_6 %>% 
  select(dif_treat_start_symp_cat, fever_2, edema_2, 
         hemorrhagic_phenomena_2, evolution, hiv_2, name_tx) %>% 
  tbl_summary(by = dif_treat_start_symp_cat, 
              #missing_stat = "{n} ({p}%)",
              missing = "no", 
              statistic = list(all_categorical() ~ "{n} ({p}%)",
                               all_continuous() ~ "{median} ({p25},{p75})"), 
              label = list(dif_treat_start_symp_cat = "Delay symptoms to treatment",                      
                           name_tx = "Treatment",
                           hiv_2 = "HIV status",
                           edema_2 = "Oedema",
                           fever_2 = "Fever",
                           hemorrhagic_phenomena_2 = "Bleeding")) %>% 
  add_overall(last = T) %>% 
  modify_table_body(filter, !(row_type == "level" &
        variable %in% c("fever_2", "edema_2", "hemorrhagic_phenomena_2") &
        label == "Unknown")) %>% 
  add_p() %>% 
  bold_labels() %>% 
  as_flex_table() %>% 
  autofit() %>% 
  save_as_docx(path = here("figures and tables", "table_2.docx"))

data_6 %>% 
  tabyl(hemorrhagic_phenomena)


#################### Box-plot delay distribution by year distribution ##########################################

ggplot(data = data_6 %>% filter(year > 2006 & evolution %in% c("Cure", "Death for VL")),
       mapping = aes(x = as.factor(year), y = dif_day_tx_symp, fill = evolution)) +
  geom_boxplot(position = position_dodge()) +
  scale_fill_manual(
    values = c(
      "Cure" = "#00BFC4",        
      "Death for VL" = "#F8766D")) +
  scale_y_continuous(n.breaks = 15) +
  labs(
    x = "Year", 
    y = "Delay between symptoms onset to treatment iniation",
    fill = "Condition") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.text = element_text(color = "black", size = 12),
    axis.title = element_text(color = "black", size = 14, face = "bold"),
    legend.title = element_text(face = "bold", size = 14, color = "black"),
    legend.text = element_text(color = "black", size = 12))
  

ggplot(data = data_6 %>% filter(year > 2006 & evolution %in% c("Cure", "Death for VL") & hiv_2 %in% c("Yes", "No")) %>% 
         mutate(hiv_3 = if_else(hiv_2 == "Yes", "HIV", "Not HIV")), 
       mapping = aes(x = as.factor(year), y = dif_day_tx_symp, fill = evolution)) +
  geom_boxplot(position = position_dodge()) +
  scale_y_continuous(n.breaks = 15) + 
  facet_wrap(~ hiv_3)


### Box-plot distribution delay between symptom onset to treatment initation by symptom 
ss_table <- data_6 %>% 
  select(fever, weakness, edema, weight_loss, 
         spleen_enlargement, increase_infection, 
         hemorrhagic_phenomena, liver_enlargement, 
         jaundice, dif_day_tx_symp) %>% 
  rename("Fever" = fever,
         "Weakness" = weakness,
         "Edema" = edema,
         "Weight loss" = weight_loss,
         "Spleen enlargment" = spleen_enlargement,
         "Infection" = increase_infection,
         "Hemorrhagic" = hemorrhagic_phenomena,
         "Liver enlargement" = liver_enlargement,
         "Jaundice" = jaundice) %>% 
  mutate(across(Fever:Jaundice, ~ na_if(., 9)),
         across(1:9, ~ (if_else(. == 1, "Yes", "No")))) %>% 
  pivot_longer(
    cols = Fever:Jaundice,
    names_to  = "ss",
    values_to = "ss_value")


ggplot(ss_table %>% filter(!is.na(ss_value)), 
       mapping = aes(x = ss, y = dif_day_tx_symp, fill = ss_value)) +
  geom_boxplot(position = position_dodge()) +
  scale_y_continuous(n.breaks = 15) +
  scale_fill_manual(
    values = c(
      "No" = "darkgreen",        
      "Yes" = "brown2")) +
  labs(
    x = "Symptom", 
    y = "Delay between symptoms onset to treatment iniation",
    fill = "Condition") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.text = element_text(color = "black", size = 12),
    axis.text.x = element_text(angle = 90),
    axis.title = element_text(color = "black", size = 14, face = "bold"),
    legend.title = element_text(face = "bold", size = 14, color = "black"),
    legend.text = element_text(color = "black", size = 12)) 


### Histogram difference days between symptoms onset 

distribution_histogram <- ggplot(data = data_6 %>% filter(!dif_day_tx_dx < 0 & dif_day_tx_symp < 366), aes(x = as.numeric(dif_day_dx_symp))) +
  geom_histogram(position = "identity", aes(y = after_stat(density), fill = "Symptoms to diagnosis"), color = "black") +
  geom_histogram(aes(x = as.numeric(dif_day_tx_dx), y = after_stat(density), fill = "Diagnosis to treatment"), alpha = 0.3, color = "black") +
  scale_x_continuous(transform = "log", n.breaks = 10, labels = label_number(accuracy = 1)) +
  scale_y_continuous(n.breaks = 10) +
  scale_fill_manual(values = c("Symptoms to diagnosis" = "skyblue", 
                               "Diagnosis to treatment" = "brown2")) +
  labs(fill = "Delay time", x = "Time (days)", y = "Probability") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.text = element_text(color = "black", size = 12),
    axis.title = element_text(color = "black", size = 14, face = "bold"),
    legend.title = element_text(face = "bold", size = 14, color = "black"),
    legend.text = element_text(color = "black", size = 12),
  )

distribution_histogram

### General cases of visceral leishmaniasis delay symptoms to treatment
distribution_histogram_2 <- ggplot(data = data_6 %>%  filter(!dif_treat_start_symp < 0 & dif_treat_start_symp < 366), aes(x = as.numeric(dif_treat_start_symp))) +
  geom_histogram(position = "identity", aes(y = after_stat(density), fill = "Symptoms to treatment"), color = "black", alpha = 0.2) +
  geom_vline(xintercept = c(q25, q50, q75, q95), linetype = "dashed", colour = c("darkgreen", "gold3", "brown2", "brown")) +
  scale_x_continuous(transform = "log", n.breaks = 10, labels = label_number(accuracy = 1)) +
  scale_y_continuous(n.breaks = 10) +
  scale_fill_manual(values = c("Symptoms to treatment" = "brown2")) +
  labs(fill = "Delay time", x = "Time (days)", y = "Probability") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.text = element_text(color = "black", size = 12),
    axis.title = element_text(color = "black", size = 14, face = "bold"),
    legend.title = element_text(face = "bold", size = 14, color = "black"),
    legend.text = element_text(color = "black", size = 12),
  ) +
  annotate(geom = "text", x = c(13, 25, 51, 138), y = c(0.5), 
           label = c("q25", "q50", "q75", "q95"), fontface = "italic",
           fontface = "bold",
           angle = 90)
  

distribution_histogram_2

## HIV histogram
h_hiv <- ggplot(data = data_6, aes(x = as.numeric(dif_day_tx_symp))) +
  geom_histogram(aes(y = after_stat(density), fill = hiv_2), position = "identity", color = "black", alpha = 0.3) +
  geom_density(aes(color = hiv_2), lwd = 1) +
  labs(x = "Difference days (Onset symptoms to treatment)", y = "Density", fill = "HIV", color = "HIV") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold")) +
  scale_x_continuous(n.breaks = 15, transform = "log", 
                     labels = function(x) format(round(x), nsmall = 0))
  
h_hiv


h_hiv <- ggplot(data = data_6, aes(x = as.numeric(dif_day_tx_symp))) +
  geom_histogram(aes( y = after_stat(density),
                      fill = hiv_2), 
                 position = "identity",
                 color = "black", alpha = 0.3) +
  geom_density(aes(color = hiv_2), lwd = 1) +
  labs(x = "Difference days (Onset symptoms to treatment)", y = "Probability", fill = "HIV", color = "HIV") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold")) +
  scale_x_continuous(n.breaks = 15, transform = "log10", 
                     labels = function(x) format(round(x), nsmall = 0))

h_hiv


#### Age distribution 

gender <- c("Male" = "#00BFC4", "Female" = "#F8766D")

age <- ggplot(data = data_6 %>% filter(!is.na(gender)), 
              aes(x = age_years, y = as.numeric(dif_day_tx_symp), color = gender)) +
  geom_point() +
  geom_smooth(method = "gam") +  # by gender
  geom_smooth(aes(x = age_years, y = as.numeric(dif_day_tx_symp), color = "General"), 
              method = "gam", se = FALSE, inherit.aes = FALSE, 
              linetype = "dashed") +  # global line
  theme_bw() +
  labs(x = "Age (years)", , 
       y = "Delay between symptoms onset to treatment",
       color = "Gender/condtion") +   # change legend title to be more general
  theme(
    legend.title = element_text(size = 14, color = "black", face = "bold"),
    axis.title = element_text(size = 14, color = "black", face = "bold"),
    axis.text = element_text(size = 12, color = "black")
  ) +
  scale_x_continuous(n.breaks = 15) +
  scale_y_continuous(n.breaks = 15) +
  scale_color_manual(values = c(gender, "General" = "black")) 
 #facet_grid(~gender)

age

### Type of cases
id_type <- ggplot(data = data_6, aes(x = as.numeric(dif_day_tx_symp))) +
  geom_histogram(aes(y = after_stat(density),fill = identification_case), 
                 position = "identity" , color = "black", alpha = 0.3) +
  geom_density(aes(color = identification_case), lwd = 1) +
  labs(x = "Difference days (Onset symptoms to treatment)", y = "Probability", 
       fill = "Type case", color = "Type case") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold")) +
  scale_x_continuous(n.breaks = 15, transform = "log10", 
                     labels = function(x) format(round(x), nsmall = 0))

id_type


### Geneder (sex)
gender <- ggplot(data = data_6 %>% filter(!is.na(gender)), aes(x = as.numeric(dif_day_tx_symp), y = -0.5)) +
  geom_boxploth (aes(fill = gender), color = "black", alpha = 0.3) +
  geom_density(aes(x = as.numeric(dif_day_tx_symp), y = after_stat(density), fill = gender), inherit.aes = FALSE) +
  #geom_point(aes(x = as.numeric(dif_day_tx_symp, y = -0.3))) +
  facet_grid(gender ~ .) +
  scale_fill_discrete() +
  #coord_flip() +
  labs(x = "Difference days (Onset symptoms to treatment)", y = "Probability", 
       fill = "Gender") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold")) +
  scale_x_continuous(n.breaks = 15, transform = "log10", 
                     labels = function(x) format(round(x), nsmall = 0))

gender


data_6 %>% 
  tabyl(dif_treat_start_symp_cat)

### age and sex pyramid 
p_1 <- apyramid::age_pyramid(data = data_6 %>% filter(dif_treat_start_symp_cat == "< 30 days"), age_group = age_group_2, 
                           split_by = gender, proportional = T, show_midpoint = F) +
  theme_bw() +
  xlab("Age group") +
  ylab("Proportion") +
  labs(fill = "Gender") +
  theme(legend.position = "bottom", 
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))

p_2 <- apyramid::age_pyramid(data = data_6 %>% filter(dif_treat_start_symp_cat == "30-89 days"), age_group = age_group_2, 
                          split_by = gender, proportional = T, show_midpoint = F) +
  theme_bw() +
  xlab("Age group") +
  ylab("Proportion") +
  labs(fill = "Gender") +
  theme(legend.position = "bottom", 
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))

p_3 <- apyramid::age_pyramid(data = data_6 %>% filter(dif_treat_start_symp_cat == "90-119 days"), age_group = age_group_2, 
                             split_by = gender, proportional = T, show_midpoint = F) +
  theme_bw() +
  xlab("Age group") +
  ylab("Proportion") +
  labs(fill = "Gender") +
  theme(legend.position = "bottom", 
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))


p_4 <- apyramid::age_pyramid(data = data_6 %>% filter(dif_treat_start_symp_cat == "> 120 days"), age_group = age_group_2, 
                             split_by = gender, proportional = T, show_midpoint = F) +
  theme_bw() +
  xlab("Age group") +
  ylab("Proportion") +
  labs(fill = "Gender") +
  theme(legend.position = "bottom", 
        axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))


ap <- ggarrange(p_1 + rremove("ylab") + rremove("xlab") ,
                p_2 + rremove("ylab") + rremove("xlab") + theme(axis.text.y = element_blank()), 
                p_3 + rremove("ylab") + rremove("xlab"), 
                p_4 + rremove("ylab") + rremove("xlab") + theme(axis.text.y = element_blank()), 
                ncol = 2, nrow = 2, common.legend = T, 
                labels = c("<30 days","30–89 days", "90–119 days", "≥120 days"),
                legend = "right", align = "hv",
                label.x = 0.05,    
                label.y = 0.98, 
                font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "lef"))

ap_annotated <- annotate_figure(ap, bottom = textGrob("Proportion", gp = gpar(cex = 1.3, fontface = "bold")),
  left = textGrob("Age group", rot = 90, vjust = 0.3, gp = gpar(cex = 1.3, fontface = "bold")))

ap_annotated 

data_6$delay_cat_label %>% tabyl()


### Fever
fever <- ggplot(data = data_6, aes(x = as.numeric(dif_day_tx_symp))) +
  geom_histogram(aes(y = after_stat(density),fill = fever_2), 
                 position = "identity" , color = "black", alpha = 0.3) +
  geom_density(aes(color = fever_2), lwd = 1) +
  labs(x = "Difference days (Onset symptoms to treatment)", y = "Probability", 
       fill = "Fever", color = "Fever") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold")) +
  scale_x_continuous(n.breaks = 15, transform = "log10", 
                     labels = function(x) format(round(x), nsmall = 0))

fever


g2 <- ggarrange(h_hiv + rremove("ylab") + rremove("xlab") + theme(axis.text.x = element_blank()), 
          fever + rremove("ylab") + rremove("xlab") + theme(axis.text.x = element_blank()), 
          gender + rremove("ylab") + rremove("xlab"),
          ncol = 1, nrow = 3)


g2_annotated <- annotate_figure(g2, bottom = textGrob("Difference days (Onset symptoms to treatment)", gp = gpar(cex = 1.3, fontface = "bold")),
                                left = textGrob("Probability", rot = 90, vjust = 0.3, gp = gpar(cex = 1.3, fontface = "bold")))


g2_annotated

gender
########################################################################################################################################################
### Distribution by federative state

map_brazil <- read_municipality(year = 2022) #### Brazil map level municipality
map_brazil <- map_brazil %>% 
  mutate(code_muni_2 = substr(code_muni, 1, 6)) %>% 
  mutate(code_muni_2 = as.double(code_muni_2))

map_ufh <- read_health_region(year = 2013) ### Health region in Brazil 

code_ufh <- map_ufh %>%  ### Selection specific variables for doing it 
  select(code_health_region, name_health_region) %>% 
  mutate(code_health_region = as.integer(code_health_region)) %>% 
  as.data.frame() %>% 
  select(-c(geom)) %>% 
  unique()

code_uf <- map_brazil %>% ### selection specific variables in federative units
  select(code_state, name_state) %>% 
  as.data.frame() %>% 
  select(-c(geom)) %>% 
  unique()


#### Attention here you need use data_7
### Join UF's names with work's database. 
data_6 <- left_join(data_6, select(code_uf, code_state, name_state), ### Federative units 
                          by = c("federative_unit" = "code_state"))

data_6 <- left_join(data_6, select(code_ufh, code_health_region, name_health_region), ### Health region 
                    by = c("regional_health_area" = "code_health_region"))

q25 <- quantile(as.numeric(data_6$dif_day_tx_symp), 0.25, na.rm = TRUE)
q50 <- quantile(as.numeric(data_6$dif_day_tx_symp), 0.50, na.rm = TRUE)
q75 <- quantile(as.numeric(data_6$dif_day_tx_symp), 0.75, na.rm = TRUE)
q95 <- quantile(as.numeric(data_6$dif_day_tx_symp), 0.95, na.rm = TRUE)

data_6.1 <- data_6 %>% 
  mutate(more_100 = if_else(federative_unit == state_100, "Endemic", "No endemic"))

data_6 %>% tabyl(name)

quant_df <- data.frame(
  Quantile = c("Q25 = 12", "Q50 = 27", "Q75 = 60", "Q95 = 165"),
  value = c(q25, q50, q75, q95))

ggplot(data = data_6.1, aes(x = name_state, y = dif_day_tx_symp)) + 
  geom_boxplot() +
  geom_hline(data = quant_df, aes(yintercept = value, color = Quantile), linetype = "dashed", size = 1) +
  theme_bw() +
  theme(
    axis.text = element_text(angle = 90)) +
  labs(x = "Name of state", y = "Delay (Symptoms onset to treatment)", 
       color = "Quantile") +
  theme(
    legend.position = "bottom",
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(color = "black", size = 14, face = "bold"),
    legend.title = element_text(face = "bold", size = 14, color = "black"),
    legend.text = element_text(color = "black", size = 12),
  ) +
  scale_y_continuous(n.breaks = 20) +
  scale_color_manual(values = c("Q25 = 12" = "seagreen", "Q50 = 27" = "dodgerblue3", 
                                "Q75 = 60" = "brown3", "Q95 = 165" = "orange3")) 


ggplot(data = data_6, aes(x = name_health_region, y = dif_day_tx_symp)) + 
  geom_boxplot() +
  geom_hline(data = quant_df, aes(yintercept = value, color = Quantile), linetype = "dashed", size = 1) +
  theme_bw() +
  theme(
    axis.text = element_text(angle = 90)) +
  labs(x = "Name of state", y = "Difference days (Symptoms onset to treatment)", 
       color = "Quantile") +
  theme(
    legend.position = "bottom",
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(color = "black", size = 14, face = "bold"),
    legend.title = element_text(face = "bold", size = 14, color = "black"),
    legend.text = element_text(color = "black", size = 12),
  ) +
  scale_y_continuous(n.breaks = 20) +
  scale_color_manual(values = c("Q25 = 12" = "seagreen", "Q50 = 27" = "dodgerblue3", 
                                "Q75 = 60" = "brown3", "Q95 = 165" = "orange3")) 


#####################################################################################################################################
### MAP delay days between onset symptoms to treatment. 
map_uf_shp <- read_sf(here("BR_UF_2023", "BR_UF_2023.shp")) %>% 
  clean_names()
map_uf <- as.data.frame(read_state(year = 2020))

mortality <- data_6
  

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

state_100 <- unique(data_ost$federative_unit[data_ost$n > 100])

data_6 %>% 
  tabyl()

years <- c("2007 - 2009", "2010 - 2014", "2015 - 2018", "2019 - 2023")

map_uf_shp_2 <- map_uf_shp %>% 
  slice(rep(1:n(), each = length(years))) %>%
  mutate(year = rep(years, times = nrow(map_uf_shp))) %>% 
  rename("year_3" = year,
         "federative_unit" = cd_uf) %>% 
  mutate(federative_unit = as.double(federative_unit))


map_uf_shp_3 <- full_join(map_uf_shp_2, select(data_ost_2, federative_unit, median, lethality, median_age, year_3)) %>% 
  mutate(median = replace_na(median, 0)) %>% 
  mutate(lethality = replace_na(lethality, 0)) %>% 
  mutate(category = factor(case_when(median > 0  & median < 20 ~ "1 to 19",
                              median >= 20 & median <40 ~ "20 to 39",
                              median >= 40 & median <59 ~ "40 to 59",
                              median >= 60 & median <79 ~ "60 to 79",
                              median >= 80 ~ "> 80",
                              TRUE ~ "Less than 100 cases")), 
         nm_uf = na_if(str_trim(nm_uf), "Distrito Federal")) %>%
  st_as_sf()

breaks <- c(10, 20, 40, 60, 80, Inf)

ggplot(map_uf_shp_3) +
  geom_sf(aes(fill = median), color = "white") +
  labs(fill = "Median delay") +
  ggtitle("Visceral Leishmaniasis Brazil, federative states with more than 100 cases") +
  scale_fill_viridis_c(limits = c(1, 80), 
                       breaks = breaks,  option = "viridis", 
                       na.value = "lightgray",) +
  theme_minimal() +
  theme(legend.position = "right",
        legend.title.position = "bottom",
        plot.title.position = element_text(hjust = 0.5)) +
  facet_wrap(~year_3, nrow = 2) 


pal_blue_red <- c(
  "Less than 100 cases" = "gray",  # light blue
  "1 to 19"  = "#A9CCF7",  # medium soft blue
  "20 to 39" = "#6FA6E0",  # deeper blue
  "40 to 59" = "#F5A3A0",  # soft red
  "60 to 79" = "#CC4C46",  # coral
  ">= 80"    = "brown" # strong warm red
)

pal_green_red <- c(
  "0"        = "#EAF7EF",
  "1 to 19"  = "#D3EFDB",
  "20 to 39" = "#B9E4C3",
  "40 to 59" = "#FFEDE2",
  "60 to 79" = "#F9C8B5",
  "80 to 99" = "#EE9E90",
  "> 100"    = "#E67C73")

  
  # get centroids of each federative unit
  map_centroids <- st_centroid(map_uf_shp_3 %>% filter(lethality > 0 & !is.na(median_age)))
  map_points <- st_point_on_surface(map_uf_shp_3)

  ggplot(map_uf_shp_3) +
    geom_sf(aes(fill = category), color = "white") +
    geom_sf(data = map_centroids, 
            aes(size = lethality), 
            color = "pink", alpha = 0.7) +
    labs(fill = "Median delay",
         size = "Lethality rate") +
    ggtitle("Visceral Leishmaniasis Brazil\nFederative states with more than 100 cases") +
    scale_fill_manual(
      values = pal_blue_red, drop = FALSE,
      name = "Delay (sympt → treatment)"
    ) +
    scale_size_continuous(range = c(1, 12)) +  # adjust dot size range
    
    theme_void() +
    theme(
      legend.position = "right",
      legend.title.position = "bottom",
      plot.title = element_text(hjust = 0.5)
    ) +
    facet_wrap(~year_3, nrow = 2)

#### Final map to identify conditions related lethality and distribution of cases in Brazil 
  ggplot(map_uf_shp_3) +
    geom_sf(aes(fill = category), color = "white") +
   # geom_sf( data = map_centroids, aes(
    #  color = lethality),
    #  alpha = 0.5,
    #  size = 3) +
    geom_sf_label(data = map_points, inherit.aes = FALSE,
      aes(label = nm_uf, geometry = geometry, vjust = -0.3), fill = NA,
      color = "black", label.size = 0, size = 2, fontface = "bold") +
    labs(fill = "Median delay", color = "Lethality rate") +
    labs(fill = "Median delay", color = "Lethality rate") +
    ggtitle("Visceral Leishmaniasis Brazil\nFederative states with more than 100 cases") +
    scale_fill_manual(values = pal_blue_red, drop = FALSE,
                      name = "Delay (sympt → treatment)") +
    scale_size_area(max_size = 12) +
    scale_color_viridis_c(option = "plasma") +
    theme_void() +
    theme(
      legend.position = "right",
      legend.title.position = "bottom",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14, color = "black"),
      strip.text = element_text(size = 12)) +
    facet_wrap(~year_3, nrow = 2) 
  
####################################################################################################################################
  
## Survival analysis 
  
  library(ggsurvfit)
  
  ### Behaviour lethality according with delay symptom onset to tretament
  
mortality_risk <- survfit2(Surv(dif_treat_outcome, death_vl) ~ 1, data = data_6) 
  
  gm <- ggsurvfit(mortality_risk, linewidth = 1, color = "pink2") +
    add_confidence_interval() +
    add_risktable(risktable_stats = "{n.risk} ({cum.event})") +
    add_censor_mark(color = "pink3") +
    #add_pvalue(caption="Log-rank {p.value}",
      #         location = "annotation", x = 50, y = 0.68) +
    scale_ggsurvfit() + 
    coord_cartesian(ylim = c(0.88, 1)) +
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
  
  print(gm)  
  
  
  ### Behaviour lethality according with delay symptom onset to tretament
  
  delay <- survfit2(Surv(dif_treat_outcome, death_vl) ~ dif_treat_start_symp_cat, data = data_6) 
  
  gd <- ggsurvfit(delay, linewidth = 1) +
    #add_confidence_interval() +
    add_risktable(risktable_stats = "{n.risk} ({cum.event})", 
                  theme = theme_bw(), risktable_height = 0.2) +
    add_censor_mark() +
    #add_pvalue(caption="Log-rank {p.value}",
    #         location = "annotation", x = 50, y = 0.68) +
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
  
  
  ### groups of age
  
  group_age <- survfit2(Surv(dif_treat_outcome, death_vl) ~ age_group_2, data = data_6) 
  
  ga <- ggsurvfit(group_age, linewidth = 1) +
    add_confidence_interval() +
    add_risktable(risktable_stats = "{n.risk} ({cum.event})", 
                  theme = theme_bw()) +
    add_censor_mark() +
    add_pvalue(caption="Log-rank {p.value}",
               location = "annotation", x = 50, y = 0.68) +
    add_risktable_strata_symbol(symbol = "\U25CF", size = 10) +
    scale_ggsurvfit() + 
    coord_cartesian(ylim = c(0.65, 1)) +
    scale_x_continuous(n.breaks = 15) +
    scale_y_continuous(n.breaks = 12, label = scales::label_percent()) +
    add_quantile(x_value = 30, linetype = "dotted", 
                 color = "grey30", linewidth = 0.8) +
    #geom_vline(xintercept = 30, linetype = "dashed", colour = "grey30") +
    labs(x = "Time (days after to start treatment)", color = "Age group", 
         fill = "Age group") +
    theme_bw() +
    theme(legend.position = "bottom", 
          axis.text = element_text(color = "black", size = 10),
          axis.title = element_text(color = "black", size = 14, face = "bold"),
          legend.title = element_text(face = "bold", size = 14, color = "black"),
          legend.text = element_text(color = "black", size = 12),
    ) 
    
  print(ga)
  
  
  ### HIV survival curves
  
  hiv <- survfit2(Surv(dif_treat_outcome, death_vl) ~ hiv_2, data = data_6) 
  
  ghiv <- ggsurvfit(hiv, linewidth = 1) +
    add_confidence_interval() +
    add_risktable(risktable_stats = "{n.risk} ({cum.event})", 
                  theme = theme_bw(), risktable_height = 0.2) +
    add_censor_mark() +
    add_pvalue(caption="Log-rank {p.value}",
               location = "annotation", x = 60, y = 0.7) +
    add_risktable_strata_symbol(symbol = "\U25CF", size = 10) +
    #scale_ggsurvfit() + 
    coord_cartesian(ylim = c(0.68, 1)) +
    scale_x_continuous(n.breaks = 15) +
    scale_y_continuous(n.breaks = 12, label = scales::label_percent()) +
    add_quantile(x_value = 30, linetype = "dotted", 
                 color = "grey30", linewidth = 0.8) +
    #geom_vline(xintercept = 30, linetype = "dashed", colour = "grey30") +
    labs(x = "Time (days after to start treatment)", color = "HIV", 
         fill = "HIV") +
    theme_bw() +
    theme(legend.position = "bottom", 
          axis.text = element_text(color = "black", size = 10),
          axis.title = element_text(color = "black", size = 14, face = "bold"),
          legend.title = element_text(face = "bold", size = 14, color = "black"),
          legend.text = element_text(color = "black", size = 10)) +
    scale_color_manual(values = c(
      "Yes"      = "#F8766D",   # red / salmon
      "No"     = "#00BA38",   # green
      "Unknown" = "#619CFF"    # blue
    )) +
    scale_fill_manual(values = c(
      "Yes"      = alpha("#F8766D", 0.18),
      "No"     = alpha("#00BA38", 0.18),
      "Unknown" = alpha("#619CFF", 0.18)
    ))
  
  
  print(ghiv)
  
  
  ### HIV survival treatment ----- be careful
  
  tto <- survfit2(Surv(dif_treat_outcome, death_vl) ~ name_tx + hiv_2, data = data_6) 
  
  gtto <- ggsurvfit(tto, linewidth = 1) +
    #add_confidence_interval() +
    add_risktable(risktable_stats = "{n.risk} ({cum.event})", 
                  theme = theme_bw(), risktable_height = 0.2) +
    add_censor_mark() +
    add_pvalue(caption="Log-rank {p.value}",
               location = "annotation", x = 60, y = 0.7) +
    add_risktable_strata_symbol(symbol = "\U25CF", size = 10) +
    #scale_ggsurvfit() + 
    coord_cartesian(ylim = c(0.68, 1)) +
    scale_x_continuous(n.breaks = 15) +
    scale_y_continuous(n.breaks = 12, label = scales::label_percent()) +
    add_quantile(x_value = 30, linetype = "dotted", 
                 color = "grey30", linewidth = 0.8) +
    #geom_vline(xintercept = 30, linetype = "dashed", colour = "grey30") +
    labs(x = "Time (days after to start treatment)", color = "Treatment", 
         fill = "Treatment") +
    theme_bw() +
    theme(legend.position = "bottom", 
          axis.text = element_text(color = "black", size = 10),
          axis.title = element_text(color = "black", size = 14, face = "bold"),
          legend.title = element_text(face = "bold", size = 14, color = "black"),
          legend.text = element_text(color = "black", size = 10)) 
  
  
  gtto +
    facet_wrap(~ hiv_2)
  

  ### Mortality comparing time between symptoms onset to start treatment and and type of the treatment
  
  data_6 <- data_6 %>% 
    mutate(edema_2 = if_else(edema == 1, "Yes oedema", "No oedema"),
           edema_2 = factor(edema_2, levels = c("Yes oedema", "No oedema")), 
           severe = if_else(edema == 1 | hemorrhagic_phenomena == 1, "Severe", "No severe"),
           severe = factor(severe, levels = c("Severe", "No severe")),
           hiv_3 = replace_na(hiv_2 == "Unknown"),
           hiv_2 = factor(hiv_2, levels = c("Yes", "No")))
  
  data_7 <- data_6 %>% 
    filter(!is.na(hiv_2))
  
  centr = theme_bw() + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  death_risk_day_tx <- survfit(Surv(dif_treat_outcome, death_vl) ~ dif_treat_start_symp_cat + name_tx + severe, 
                               data = data_6 %>% 
                                 filter(str_detect(name_tx, "Antimonial|Amphote|Lipos|Pentamidine"))) 
  
  death_risk_day_tx_2 <- survfit2(Surv(dif_treat_outcome, death_vl) ~ dif_treat_start_symp_cat + name_tx + hiv_2, 
                                 data = data_7%>% 
                                   filter(str_detect(name_tx, "Antimonial|Amphote|Lipos|Pentamidine") & !is.na(hiv_2)))
                                
  
  death_risk_day_hiv <- survfit2(Surv(dif_treat_outcome, death_vl) ~ name_tx + hiv_2 , 
                                data = data_6 %>% 
                                filter(str_detect(name_tx, "Antimonial|Amphote|Lipos") & !hiv_2 == "Unknown"))
  
   lethaly_risk_symp_tx <- survfit(Surv(dif_treat_outcome, death_vl) ~ dif_treat_start_symp_cat , data = data_6 %>% 
                                    filter(str_detect(name_tx, "Antimonial|Amphote|Lipos")))
  
  gstx <- ggsurvplot(fit = lethaly_risk_symp_tx,
                     data = data_6, 
                     pval = T, pval.coord = c(330, .95),
                     #ggtheme = theme_bw(),
                     #risk.table = T,
                     tables.theme = theme_cleantable(),
                     ylim = c(0.8, 1), 
                     font.x = c(15, "bold"), font.y = c(15, "bold"),
                     xlab = "days",
                     legend = c("bottom"), 
                     legend.title = "Delay days symptoms onset and treatment",
                     legend.labs = c("<= 30 days", "30–89 days", "90–119 days", ">= 120 days"),
                     legend.nrow = 2) 
  
  gstx <- gstx$plot + theme_bw() +
    theme(legend.title = element_text(face = "bold", size = 14),
          legend.text =  element_text(size = 12), 
          axis.title = element_text(face = "bold"), 
          legend.position = "bottom") +
    
    scale_x_continuous(n.breaks = 15)
  
  print(gstx)
  
  
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
  
  
##### Treatment by HIV condition - survival analysis. 
  
  gstx_tx <- ggsurvplot_facet(fit = death_risk_day_tx_2,
                              data = data_7,
                              facet.by = c("dif_treat_start_symp_cat", "hiv_2"),
                              palette = "lancet",
                              pval = T, pval.coord = c(330, .99),
                              #title = "Time to death according with the delay symptoms onset to treatment and 
                              #type of treatment used and presence of complications",
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
                                              "Pentamidine"
                                              ),  
                              panel.labs.font = list(size = 12, face = "bold"),
                              short.panel.labs = T) +
    scale_x_continuous(n.breaks = 15)
  
  
  gstx_tx
  
  ##################################################################################################################
  
  
  
