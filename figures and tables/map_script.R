pacman::p_load(geobr, terra, sf, gganimate, tidyverse, flextable, rio, here, gtsummary, janitor, epikit)
install.packages("geobr")
library(geobr)
####################################################################################################
#### Import Brazil map from package geobr
map_brazil_hr <- read_health_region(year = 2013) ### importing Brazil healthcare region map
map_brazil_region <- read_state(year = 2020)   ### importing Brazil region map
map_brazil <- read_municipality(year = 2022) #### Brazil map level municipality
map_brazil <- map_brazil %>% 
  mutate(code_muni_2 = substr(code_muni, 1, 6)) %>% 
  mutate(code_muni_2 = as.double(code_muni_2)) ### Ajusting municipalitiy code because VL database 

# Define the years you want to assign
years <- 2007:2023

# Repeat the data and assign years
map_brazil_2 <- map_brazil %>% 
  slice(rep(1:n(), each = length(years))) %>%
  mutate(year = rep(years, times = nrow(map_brazil)))


map_brazil <- map_brazil %>%   
  mutate(code_state = as.character(code_state),
         code_name = paste(name_muni, abbrev_state, sep = ""))

#################################################################################################################################
### Import database with cases
data <- import(here("data base and clean", "data_4.rda"))
data_p <- import(here("data base and clean", "population_brasil.xlsx")) ### population map 

#### Cases per municipality for including in the map and showing incidence over time (2007-2023)
map_cases <- data %>% 
  group_by(code_of_the_municipality, year) %>% 
  summarise(n = n(), .groups = "drop") %>% ### Cases per municipality including year and geom 
  rename("code_muni_2" = code_of_the_municipality)


map_brazil_3 <- left_join(map_brazil_2, select(map_cases, code_muni_2, year, n),
                          by = c("code_muni_2", "year")) ## Joining cases with space information about Brazil. 

### Export joined database for using in Shiny. 
rio::export(map_brazil_3, here("data base and clean", "map_brazil.rda"))
map_brazil_3 <- import(here("data base and clean", "map_brazil.rda"))

### Replace NA values in n (counts) to be able reading when I make the map. 
map_brazil_3 <- map_brazil_3 %>% 
  mutate(n = replace_na(n, 0)) %>% 
  st_as_sf()

### I review database and map database after the joined have same number of municipalities with cases. 
map_brazil_3 %>% 
  filter(!is.na(n)) %>% 
  nrow()

### Single map for eaching year.
animate_map <- ggplot(map_brazil_3) +
  geom_sf(aes(fill = n), lwd = 0.01) + 
  labs(fill = "Cases") +
  ggtitle("Visceral Leishmaniasis Brazil,
          {frame_time}") +
  scale_fill_viridis_c(limits = c(1, 300), option = "plasma", 
                       na.value = 0, ) +
  theme_minimal() +
  transition_time(year) +
  shadow_mark()

# Number of frames = number of years from 2007 to 2023
num_years <- 2023 - 2007 + 1

# Create animation object
anim_obj <- animate(animate_map, nframes = num_years, fps = 2,
                    width = 1000,
                    height = 800,
                    res = 200,
                    renderer = gifski_renderer(loop = TRUE))

# Save animation
anim_save("vl_Animation.gif", animation = anim_obj)


############################################################################################################################################################

### Joining brazil map and database VL
#data_4 <- left_join(data_4, select(map_brazil, geom, code_muni_2, name_muni),
 #                   by = c("code_of_the_municipality" = "code_muni_2") 
                    
map_cases_uf <- data_6 %>% 
  group_by(federative_unit, geom, name_state, year) %>% 
  summarise(n = n(), .groups = "drop") %>% ### Cases per municipality including year and geom 


map_cases %>% 
  tabyl(year)

data_p_4 <- full_join(data_p_3, select(map_cases, n, year, code_muni_2)) %>% 
  mutate(n = replace_na(n, 0)) %>% 
  mutate(year = as.numeric(year)) %>% 
  st_as_sf() ### Final database for making map 


### Brazil map cases 
plot(map_brazil_2$geom)

animate_map <- ggplot(data_p_4) +
  geom_sf(aes(fill = n), color = "black") + 
  labs(fill = "Cases") +
  theme_minimal() +
  ggtitle('Year: {year}') +
  scale_fill_viridis_c(limits = c(1, 300)) +
  transition_time(year) +
  shadow_mark()

# Number of frames = number of years from 2008 to 2023
num_years <- 2023 - 2008 + 1

# Create animation object
anim_obj <- animate(animate_map, nframes = num_years, fps = 2, 
                    width = 1000,     # increase width
                    height = 800,     # increase height
                    res = 200,        # higher DPI,
                    renderer = gifski_renderer(loop = T))

# Save animation
anim_save("vl_Animation.gif", animation = anim_obj)


#### Brazil population 
data_p <- import(here("data base and clean", "population_brasil.xlsx"))

### Clean database population brazil
data_p <- data_p %>% 
  clean_names() %>% 
  select(cod_r, population, year) %>% 
  rename("code_of_the_municipality" = cod_r) 

### group data for calculating cases per municipality 
rate_map <- data_4 %>% 
  group_by(year, code_of_the_municipality, name_muni, geom) %>% 
  summarise(n = n(), .groups = "drop") %>% 
  mutate(code_of_the_municipality = as.character(code_of_the_municipality))

#### Calculate of rate for eaching municipalty by year   
rate_map <- full_join(rate_map, data_p) 
rate_map_2 <- rate_map %>% 
  filter(year %in% c(2008, 2009, 2011, 2012, 2013, 2014,
                     2015, 2016, 2017, 2018, 2019, 2020,
                     2021)) %>% 
  #filter(!is.na(population)) %>% 
  #filter(!is.na(n)) %>% 
  mutate(rate = round(n/population*100000, digits = 2)) %>% 
  mutate(rate = replace_na(rate, 0)) 
#st_as_sf()

#### rate adjusted whole country 
rate_map_3 <- rate_map_2 %>% 
  group_by(code_of_the_municipality, name_muni) %>% 
  summarise(median = median(n),
            population = mean(population),
            rate = round(median/population*10000, digits = 3), .groups = "drop") 

rate_map_4 <- rate_map_2 %>% 
  mutate(uf = substr(code_of_the_municipality, 1, 2)) %>% 
  mutate_at(vars(n, population), ~ replace_na(., 0)) %>% 
  group_by(uf) %>% 
  summarise(cases = sum(n),
            population = sum(population),
            rate = round(cases/population*100000, digits = 3), .groups = "drop")

map_brazil <-  map_brazil %>% 
  mutate(code_muni_2 = as.character(code_muni_2))

####  Region map Brazil consolid
rate_map_4 <- left_join(rate_map_4, select(map_brazil_region, code_state, name_state, geom),
                        by = c("uf" = "code_state")) %>% 
  st_as_sf()

map_brazil_3 <- left_join(map_brazil, select(rate_map_2, rate, year, code_of_the_municipality),
                          by = c("code_muni_2" = "code_of_the_municipality")) %>% 
  mutate(rate = replace_na(rate, 0)) 

my_breaks <- c(0, 2.5, 5, 7.5, 10, 12.5, 15, Inf)

ggplot(rate_map_4) +
  geom_sf(aes(fill = rate)) +
  labs(fill = "Rate x 100.000 people") + 
  theme_minimal() +
  #geom_text(aes(label = name_state)) +
  #facet_wrap(~year) +
  scale_fill_viridis_c(breaks = my_breaks, labels = my_breaks, limits = c(1, 17))