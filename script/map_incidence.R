# ══════════════════════════════════════════════════════════════════════════════
# Project : Visceral Leishmaniasis in Brazil, 2007–2023
# Script  : map_incidence.R
# Purpose : Compute VL incidence per 100,000 inhabitants by federative state
#           and four-year time period; produce choropleth maps faceted by period.
# Input   : data base and clean/data_6.rds           (SINAN VL notifications)
#           data base and clean/population_brasil.xlsx (IBGE population estimates)
#           data base and clean/map_uf_shp_df.rds     (sf polygon object with
#             period-specific delay/lethality summaries; produced by
#             descriptive_analysis.R and cached as .rds)
# Output  : data base and clean/map_3.rds
#           figures and tables/incidence_map.png
# Depends : map_uf_shp_3 and map_points must exist in the R environment
#           (created in descriptive_analysis.R) or map_uf_shp_df.rds must be
#           present in "data base and clean/".
#           If running this script standalone, rebuild map_points with:
#             map_points <- sf::st_point_on_surface(map_4)
# Author  : Juan Sebastián Hurtado Zapata
#           IDDO VL Fellowship | CIDEIM, Cali, Colombia
# Updated : 2026-07-29
# ══════════════════════════════════════════════════════════════════════════════


# ── 1. PACKAGES ───────────────────────────────────────────────────────────────
# pacman::p_load installs missing packages and loads them in one call.
# To install pacman itself: install.packages("pacman")
pacman::p_load(
  tidyverse,   # data manipulation and ggplot2
  janitor,     # clean_names()
  rio,         # import/export
  here,        # project-relative paths
  sf,          # spatial operations (st_as_sf, st_point_on_surface)
  paletteer,   # extended colour palettes (paletteer_c)
  ggpubr       # ggsave helpers and figure annotation
)


# ── 2. DATA IMPORT ────────────────────────────────────────────────────────────

# Population denominators by state and year (IBGE census projections)
data_pop <- import(here("data base and clean", "population_brasil.xlsx"))

# SINAN VL notifications (curated; produced by curated_sinan_data.R)
data_6 <- import(here("data base and clean", "data_6.rds"))

# Spatial polygon object with period-aggregated delay/lethality metrics.
# Loaded from cache so this script can run without re-running
# the full descriptive_analysis.R pipeline.
map <- import(here("data base and clean", "map_uf_shp_df.rds"))


# ── 3. POPULATION DENOMINATORS BY STATE AND PERIOD ───────────────────────────
# Aggregate total population by state for each four-year period.
# We select the mid-year of each period as the representative denominator:
#   2007–2009 → 2008 | 2010–2014 → 2012 | 2015–2018 → 2017 | 2019–2023 → 2021
data_pop_2 <- data_pop %>%
  clean_names() %>%
  group_by(cod_uf, uf, year) %>%
  summarise(pop = sum(population, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    year_3 = case_when(
      year >= 2007 & year <= 2009 ~ "2007 - 2009",
      year >= 2010 & year <= 2014 ~ "2010 - 2014",
      year >= 2015 & year <= 2018 ~ "2015 - 2018",
      year >= 2019 & year <= 2023 ~ "2019 - 2023",
      TRUE ~ NA
    )
  ) %>%
  filter(!is.na(year_3)) %>%
  filter(year %in% c(2008, 2012, 2017, 2021))  # mid-period denominators


# ── 4. CASE COUNTS BY STATE AND PERIOD ───────────────────────────────────────
data_cases <- data_6 %>%
  mutate(
    year   = lubridate::year(as.Date(date_of_notification)),
    year_3 = case_when(
      year >= 2007 & year <= 2009 ~ "2007 - 2009",
      year >= 2010 & year <= 2014 ~ "2010 - 2014",
      year >= 2015 & year <= 2018 ~ "2015 - 2018",
      year >= 2019 & year <= 2023 ~ "2019 - 2023",
      TRUE ~ NA_character_
    )
  ) %>%
  group_by(federative_unit, year_3) %>%
  summarise(cases = n(), .groups = "drop")


# ── 5. JOIN POPULATION AND CASES TO SPATIAL OBJECT ───────────────────────────
# Convert to plain data frame first; spatial attributes are re-added at step 6.
map <- as.data.frame(map)

map_2 <- left_join(
  map,
  dplyr::select(data_pop_2, cod_uf, year_3, pop),
  by = c("year_3", "federative_unit" = "cod_uf")
)

map_3 <- left_join(
  map_2,
  dplyr::select(data_cases, federative_unit, year_3, cases),
  by = c("year_3", "federative_unit")
)

# Cache joined object to avoid re-running this pipeline downstream
export(map_3, here("data base and clean", "map_3.rds"))


# ── 6. COMPUTE INCIDENCE RATE ─────────────────────────────────────────────────
# Incidence = VL cases / mid-period population × 100,000 person-years
map_4 <- map_3 %>%
  mutate(incidence = round(cases / pop * 100000, digits = 2)) %>%
  st_as_sf()


# ── 7. INCIDENCE CHOROPLETH MAP ───────────────────────────────────────────────
# Diverging temperature palette highlights the gradient from low- to
# high-burden states. State name labels come from map_points (centroids),
# produced in descriptive_analysis.R via st_point_on_surface(map_uf_shp_3).
# If running standalone, rebuild with: map_points <- st_point_on_surface(map_4)

incidence_map <- ggplot(map_4) +
  geom_sf(aes(fill = incidence), color = "white") +
  geom_sf_label(
    data        = map_points, inherit.aes = FALSE,
    aes(label   = nm_uf, geometry = geometry, vjust = -0.3),
    fill        = NA, color = "black",
    label.size  = 0, size = 2, fontface = "bold"
  ) +
  ggtitle("Incidence of VL cases by time period — Brazil") +
  scale_fill_gradientn(
    colours = paletteer_c("ggthemes::Temperature Diverging", 30),
    name    = "Incidence per 100,000 inhabitants"
  ) +
  theme_void() +
  theme(
    legend.position       = "right",
    legend.title.position = "bottom",
    legend.title          = element_text(size = 12, face = "bold.italic", color = "black"),
    plot.title            = element_text(hjust = 0.5, face = "bold", size = 14, color = "black"),
    strip.text            = element_text(size = 12, face = "bold")
  ) +
  facet_wrap(~year_3, nrow = 2)

print(incidence_map)

ggsave(
  plot     = incidence_map,
  path     = here("figures and tables"),
  filename = "incidence_map.png",
  width    = 10, height = 8, dpi = 300
)
