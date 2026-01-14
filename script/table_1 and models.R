pacman::p_load(rms)

## variables in the table
var_order <- c("age_group_2", "gender", "year_period", "region",
  "education_level", "race_2", "hiv_2", "name_tx", "evolution")

## desired level orders WITH the reference as first element
level_orders <- list(
  age_group_2 = c("0-1", "1-5", "6-15", "16-50", ">50"),
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
        median = median(dif_day_tx_symp, na.rm = TRUE),
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
  select(-gm_ratio, -gm_lower, -gm_upper, -is_ref) %>% 
  rename(
    gm_ratio = gm_ratio_chr,
    gm_lower = gm_lower_chr,
    gm_upper = gm_upper_chr
  )

border_line = officer::fp_border(color = "black", width = 1.38) 

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
    "Odds Ratio"              = gm_ratio,
    "Lower 95% CI"            = gm_lower,
    "Upper 95% CI"            = gm_upper
  ) %>% 
  flextable() %>% 
  merge_at(j = 1, i = 1:5, part = "body") %>% 
  merge_at(j = 1, i = 6:7, part = "body") %>% 
  merge_at(j = 1, i = 8:11, part = "body") %>% 
  merge_at(j = 1, i = 12:16, part = "body") %>% 
  merge_at(j = 1, i = 17:23, part = "body") %>% 
  merge_at(j = 1, i = 24:29, part = "body") %>% 
  merge_at(j = 1, i = 30:32, part = "body") %>% 
  merge_at(j = 1, i = 33:36, part = "body") %>% 
  merge_at(j = 1, i = 37:41, part = "body") %>% 
  hline(i = c(5, 7, 11, 16, 23, 29, 32, 36), part = "body", border = border_line) %>% 
  align(j = 3:8, part = "body", align = "center") %>% 
  align(j = 3:8, part = "header", align = "center") %>% 
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
         !is.na(dif_day_tx_symp)) %>% 
  mutate(delay5 = dif_day_tx_symp / 7)

fit_delay_uni <- glm(death_vl ~ delay5,
                     data = data_mort,
                     family = binomial())

coef_uni <- summary(fit_delay_uni)$coefficients["delay5", ]

OR_delay_uni    <- exp(coef_uni["Estimate"])
L95_delay_uni   <- exp(coef_uni["Estimate"] - 1.96 * coef_uni["Std. Error"])
U95_delay_uni   <- exp(coef_uni["Estimate"] + 1.96 * coef_uni["Std. Error"])

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

