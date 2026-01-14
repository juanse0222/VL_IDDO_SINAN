stacked_plot <- data_7 %>% 
  filter(!evolution == "Neglect" & !is.na(hiv_2)) %>% 
  group_by(year, evolution, hiv_2) %>% 
  summarise(n = n(), .groups = "drop") %>% 
  mutate(freq = n / sum(n),
         per = scales::label_percent()(freq),
         hiv_2 = case_when(hiv_2 == "Yes" ~ "HIV",
                           hiv_2 == "No" ~ "No HIV", 
                           TRUE ~ NA))

ggplot(stacked_plot, 
       aes(x = factor(year), y = freq, fill = evolution)) +
  geom_bar(stat = "identity", position = "fill", color = "black") +
  labs(y = "Percentage", x = "Year", fill = "Evolution") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 14, face = "bold", color = "black"),
    legend.title = element_text(size = 14, face = "bold", color = "black"),
    legend.text = element_text(size = 12, color = "black"),
    axis.text = element_text(size = 12, color = "black"),
    axis.text.x = element_text(angle = 90),
    strip.text = element_text(size = 14, color = "black", face = "bold")
  ) +
  scale_y_continuous(labels = scales::label_percent(), n.breaks = 10) +
  facet_wrap(~hiv_2)
  #geom_text(aes(label = per), 
          #  position = position_fill(vjust = 0.5))



#################### Treatment according with delay ######################################################

stacked_plot <- data_7 %>% 
  filter(!evolution == "Neglect" & !is.na(hiv_2)) %>% 
  mutate(dif_delay_cat_2 = factor(case_when(dif_day_tx_symp > 6 & dif_day_tx_symp < 16 ~ "7 to 15 days",
                                     dif_day_tx_symp > 15 & dif_day_tx_symp < 31 ~ "16 to 30 days",
                                     dif_day_tx_symp > 30 & dif_day_tx_symp < 61 ~ "31 to 60 days",
                                     dif_day_tx_symp > 60 & dif_day_tx_symp < 91 ~ "61 to 90 days",
                                     dif_day_tx_symp > 90 & dif_day_tx_symp < 181 ~ "90 to 180 days",
                                     dif_day_tx_symp > 180 & dif_day_tx_symp < 366 ~ "181 to 365 days",
                                     TRUE ~ NA))) %>% 
  group_by(dif_delay_cat_2, name_tx, evolution, hiv_2) %>% 
  summarise(n = n(), .groups = "drop") %>% 
  mutate(freq = n / sum(n),
         per = scales::label_percent()(freq),
         dif_delay_cat_2 = factor(dif_delay_cat_2, levels = c("7 to 15 days", "16 to 30 days",
                                                              "31 to 60 days", "61 to 90 days",
                                                              "90 to 180 days", "181 to 365 days")),
         hiv_2 = case_when(hiv_2 == "Yes" ~ "HIV",
                           hiv_2 == "No" ~ "No HIV", 
                           TRUE ~ NA)) %>% 
  filter(!is.na(dif_delay_cat_2) & evolution %in% c("Death for VL", "Cure")) %>% 
  mutate(evolution = if_else(evolution == "Cure", "Cure", "Death by VL"))

ggplot(stacked_plot, 
       aes(x = dif_delay_cat_2, y = freq, fill = hiv_2)) +
  geom_bar(stat = "identity", position = "fill", color = "black") +
  labs(y = "Percentage", x = "Delay symptom onset to treatment", fill = "Treatment") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 14, face = "bold", color = "black"),
    legend.title = element_text(size = 14, face = "bold", color = "black"),
    legend.text = element_text(size = 12, color = "black"),
    axis.text = element_text(size = 12, color = "black"),
    axis.text.x = element_text(angle = 90),
    strip.text = element_text(size = 14, color = "black", face = "bold")
  ) +
  scale_y_continuous(labels = scales::label_percent(), n.breaks = 10) +
  facet_wrap(~evolution)


data_x <- data_7 %>% 
  filter(dif_day_tx_symp > 7 & !is.na(hiv_2)) %>% 
  group_by(dif_day_tx_symp, evolution, hiv_2, age_group_2) %>% 
  summarise(n = n(), .groups = "drop") %>% 
  ggplot(aes( x= dif_day_tx_symp, y = n)) +
  geom_line(aes(color = evolution)) +
  facet_wrap(~hiv_2 + age_group_2)


data_x_2 <- data_6 %>% 
  filter(dif_day_tx_symp > 7 & !is.na(hiv_2)) %>% 
  group_by(name_tx, evolution, hiv_2, age_group_2) %>% 
  summarise(n = n(), .groups = "drop") 
  #mutate(freq = n / sum(n),
   #      per = scales::label_percent()(freq),
    #     hiv_2 = case_when(hiv_2 == "Yes" ~ "HIV",
     #                      hiv_2 == "No" ~ "No HIV", 
      #                     TRUE ~ NA))
  
plot_data <- data_x_2 %>%
  filter(hiv_2 == "Yes")

ggplot(data_x_2, aes(x = age_group_2, y = n, fill = name_tx)) +
  geom_col(position = position_dodge(width = 0.9), color = "black", width = 0.8) +
  # show case counts above bars, aligned with dodge
  geom_text(aes(label = n),
            position = position_dodge(width = 0.9),
            vjust = -0.3, size = 3.5) +
  scale_y_continuous(n.breaks = 10, expand = expansion(mult = c(0, 0.05))) +
  labs(y = "Cases", x = "Age (yrs)", fill = "Treatment") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 14, face = "bold", color = "black"),
    legend.title = element_text(size = 14, face = "bold", color = "black"),
    legend.text = element_text(size = 12, color = "black"),
    axis.text = element_text(size = 12, color = "black"),
    axis.text.x = element_text(angle = 0, vjust = 0.5),
    strip.text = element_text(size = 14, color = "black", face = "bold")
  ) +
  facet_wrap(~ evolution + hiv_2)   # or facet_wrap(~ evolution + hiv_2) if


  labs(y = "Percentage", x = "Years", fill = "Treatment") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 14, face = "bold", color = "black"),
    legend.title = element_text(size = 14, face = "bold", color = "black"),
    legend.text = element_text(size = 12, color = "black"),
    axis.text = element_text(size = 12, color = "black"),
    axis.text.x = element_text(angle = 0),
    strip.text = element_text(size = 14, color = "black", face = "bold")
  ) +
  scale_y_continuous(n.breaks = 10) +
  facet_wrap(~evolution + hiv_2)
  

#geom_text(aes(label = per), 
#  position = position_fill(vjust = 0.5))


