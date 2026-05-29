library(paletteer)

### Import dataset with population data of Brasil between 2007 and 2023
data_pop <- import(here("data base and clean", "population_brasil.xlsx"))
data_6 <- import(here("data base and clean", "data_6.rds"))
map <- import(here("data base and clean", "map_uf_shp_df.rds"))

### 1. Clean and curated dataset and set new features (variables) - Brasil population
data_pop_2 <- data_pop %>% 
  clean_names() %>% 
  group_by(cod_uf, uf, year) %>% 
  summarise(pop = sum(population, na.rm = T), .groups = "drop") %>% 
  mutate(year_3 = case_when(year >= 2007 & year <= 2009 ~ "2007 - 2009",
                            year >= 2010 & year <= 2014 ~ "2010 - 2014",
                            year >= 2015 & year <= 2018 ~ "2015 - 2018",
                            year >= 2019 & year <= 2023 ~ "2019 - 2023",
                            TRUE ~ NA)) %>% 
  filter(!is.na(year_3)) %>% 
  filter(year %in% c(2008, 2012, 2017, 2021))

### 2. Clean and curated dataset and set new features (variables) - case of VL (SINAN dataset)
data_cases <- data_6 %>% 
  mutate(
    year    = lubridate::year(as.Date(date_of_notification)),  # ← as.Date() repairs the IDate corruption
    year_3  = case_when(
      year >= 2007 & year <= 2009 ~ "2007 - 2009",
      year >= 2010 & year <= 2014 ~ "2010 - 2014",
      year >= 2015 & year <= 2018 ~ "2015 - 2018",
      year >= 2019 & year <= 2023 ~ "2019 - 2023",
      TRUE ~ NA_character_                                     # ← keep NA_character_ to match string type
    )
  ) %>% 
  group_by(federative_unit, year_3) %>% 
  summarise(cases = n(), .groups = "drop")

### 3. 

map <- as.data.frame(map_uf_shp_3)

map_2 <- left_join(map, dplyr::select(data_pop_2, cod_uf, year_3, pop),
                   by = c("year_3", "federative_unit" = "cod_uf"))

map_3 <- left_join(map_2, dplyr::select(data_cases, federative_unit, year_3, cases),
                   by = c("year_3", "federative_unit"))

export(map_3, here("data base and clean", "map_3.rds")) ### export incidence map format rds

map_4 <- map_3 %>% 
  mutate(incidence = round(cases/pop*100000, digits = 2)) %>% 
  st_as_sf()

### 5. Plot incidence map 
ggplot(map_4) +
  geom_sf(aes(fill = incidence), color = "white") +
  geom_sf_label(data = map_points, inherit.aes = FALSE,
                aes(label = nm_uf, geometry = geometry, vjust = -0.3), fill = NA,
                color = "black", label.size = 0, size = 2, fontface = "bold") +
  ggtitle("Incidence of VL cases by range of years in Brasil") +
  scale_fill_gradientn(
    colours = paletteer_c("ggthemes::Temperature Diverging", 30),
    name = "Incidence x 100,000 habitants"
  ) +
  #scale_color_viridis_c(option = "plasma") +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title.position = "bottom",
    legend.title = element_text(size = 12, face = "bold.italic", color = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14, color = "black"),
    strip.text = element_text(size = 12)) +
  facet_wrap(~year_3, nrow = 2) 
