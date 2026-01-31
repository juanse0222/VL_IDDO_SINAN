library(read.dbc)

# I have downloded the .dbc files into the local folder "~/data/arquivo/"
data_folder <- "arquivo/"
all_files <- list.files(data_folder, pattern = "\\.dbc$", full.names = TRUE)

# exclude database from 2005, 2006 and 2024 (2005 and 2006 are .net system, and 2024 data are preliminary)
filtered_files <- all_files[!grepl("(05|06|24)\\.dbc$", all_files)]

# import files using read.dbc function from read.dbc package
data_list <- list()

for (i in seq_along(filtered_files)) {
  data_list[[i]] <- read.dbc(filtered_files[i])
}

# 151,604 rows in the binded datset for 2007-2023
combined_data <- do.call(rbind, data_list)
combined_data %>% names()

saveRDS(combined_data, "data base and clean/publicSINAN.rds")

combined_data %>% count(SG_UF_NOT, ID_MUNICIP)
######################################################################################################################################

### Review and evaluation dataset with cases from SINAN

df <- combined_data %>% 
  clean_names()

df_2 <- df %>% 
  mutate(dif_dx_symp = dt_invest - dt_sin_pri)

summary(df_2$dif_dx_symp)

df_2.1 <- df_2 %>% 
  filter(dif_dx_symp %in% ())

df_3 <- df_2 %>% filter(dif_dx_symp < 0)

df_3 <- df_3 %>% 
  mutate(dt_invest_2 = paste0(year(dt_invest),"-",
                                    day(dt_invest), "-",
                                    month(dt_invest)),
         dt_invest_2 = as.Date(dt_invest_2)) %>% 
  mutate(dif_dx_symp_2 = dt_invest_2 - dt_sin_pri,
         dt_invest = dt_invest_2,
         dif_dx_symp = dt_invest - dt_sin_pri) %>% 
  select(-c(dif_dx_symp_2, dt_invest_2))

df_4 <- rbind(df_2.1, df_3)         
