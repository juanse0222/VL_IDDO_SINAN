
color_distribution <- c("log-normal" = "red",
                        "gamma" = "orange")

# 2. Visualize to choose candidates
m1 <- simukde::find_best_fit(x = data_6$dif_day_tx_symp, positive = TRUE)

# 3. Fit different distributions (e.g., Normal, Gamma, Lognormal)
fit_norm <- fitdist(data_6$dif_day_tx_symp, "norm")
fit_gamma <- fitdist(data_6$dif_day_tx_symp, "gamma")
fit_lnorm <- fitdist(data_6$dif_day_tx_symp, "lnorm")
fit_cauchy <- fitdist(data_6$dif_day_tx_symp, "cauchy")


# 4. Compare AIC values
summary(fit_norm)$aic
summary(fit_gamma)$aic
summary(fit_lnorm)$aic
summary(fit_cauchy)$aic

summary(fit_lnorm)

# 5. Comparative goodness-of-fit statistics
###

media_d0 <- log10(median(data_6$dif_day_tx_symp, na.rm = T))

d0 <- ggplot(
  data_6,
  aes(x = log10(dif_day_tx_symp))) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 30,
    fill = "purple2",
    alpha = 0.3,
    color = "black") +
  stat_function(
    fun = function(z)
      dgamma(10^z,
             shape = fit_gamma$estimate["shape"],
             rate  = fit_gamma$estimate["rate"]) * log(10) * 10^z,
    aes(color = "gamma"),
    linewidth = 1) +
  stat_function(
    fun = function(z)
      dlnorm(10^z,
             meanlog = fit_lnorm$estimate["meanlog"],
             sdlog   = fit_lnorm$estimate["sdlog"]) *
      log(10) * 10^z,
    aes(color = "log-normal"),
    linewidth = 1
  ) +
  geom_vline(xintercept = media_d0, linetype = "dashed", 
             linewidth = 1, color = "black") +
  labs(
    x = expression(bold(log[10]~"(Days since symptom onset)")),
    y = "Density",
    title = "Overall", 
    color = "Distribution"
  ) +
  coord_cartesian(xlim = c(0, log10(500))) +
  scale_x_continuous(
    breaks = log10(break_vals),
    labels = break_vals,
    limits = c(0, log10(500))   # keep or change as needed
  ) +
  scale_y_continuous(limits = c(0, 1.5),
                     n.breaks = 7) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, color = "black", size = 16, face = "bold"),
    axis.title = element_text(size = 14, color = "black", face = "bold"),
    axis.text = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 14, color = "black", face = "bold"),
    legend.text = element_text(size = 12, color = "black"),
    legend.position = c(0.9, 0.9)
  ) +
  scale_color_manual(values = color_distribution)



##### Regions distribution ################
break_vals <- c(1, 3, 5, 10, 30, 50, 100, 300)


### Center-West
# 2. Visualize to choose candidates
ne <- simukde::find_best_fit(x = data_6$dif_day_tx_symp[data_6$region == "Center-West"], positive = TRUE)

# 3. Fit different distributions (e.g., Normal, Gamma, Lognormal)
fit_norm_1 <- fitdist(data_6$dif_day_tx_symp[data_6$region == "Center-West"], "norm")
fit_gamma_1 <- fitdist(data_6$dif_day_tx_symp[data_6$region == "Center-West"], "gamma")
fit_lnorm_1 <- fitdist(data_6$dif_day_tx_symp[data_6$region == "Center-West"], "lnorm")
fit_cauchy_1 <- fitdist(data_6$dif_day_tx_symp[data_6$region == "Center-West"], "cauchy")


# 4. Compare AIC values
summary(fit_norm_1)$aic
summary(fit_gamma_1)$aic
summary(fit_lnorm_1)$aic
summary(fit_cauchy_1)$aic

summary(fit_lnorm)

media_d1 <- log10(median(data_6$dif_day_tx_symp[data_6$region == "Center-West"], na.rm = T))
10^(media_d1)

d1 <- ggplot(
  data_6 %>% filter(region == "Center-West"),
  aes(x = log10(dif_day_tx_symp))) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 30,
    fill = "grey80",
    color = "black") +
  stat_function(
    fun = function(z)
      dgamma(10^z,
             shape = fit_gamma_1$estimate["shape"],
             rate  = fit_gamma_1$estimate["rate"]) * log(10) * 10^z,
    aes(color = "gamma"),
    linewidth = 1) +
  stat_function(
    fun = function(z)
      dlnorm(10^z,
             meanlog = fit_lnorm_1$estimate["meanlog"],
             sdlog   = fit_lnorm_1$estimate["sdlog"]) *
      log(10) * 10^z,
    aes(color = "log-normal"),
    linewidth = 1
  ) +
  geom_vline(xintercept = media_d1, linetype = "dashed", 
             linewidth = 1, color = "black") +
  labs(
    x = expression(bold(log[10]~"(Days since symptom onset)")),
    y = "Density",
    title = "Center-West",
    color = "Distribution"
  ) +
  coord_cartesian(xlim = c(0, log10(500))) +
  scale_x_continuous(
    breaks = log10(break_vals),
    labels = break_vals,
    limits = c(0, log10(500))   # keep or change as needed
  ) +
  scale_y_continuous(limits = c(0, 1.5),
                     n.breaks = 7) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, color = "black", size = 16, face = "bold"),
    axis.title = element_text(size = 14, color = "black", face = "bold"),
    axis.text = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 14, color = "black", face = "bold"),
    legend.text = element_text(size = 12, color = "black"),
    legend.position = c(0.9, 0.9)
  ) +
  scale_color_manual(values = color_distribution)


### Northeast
# 2. Visualize to choose candidates
ne <- simukde::find_best_fit(x = data_6$dif_day_tx_symp[data_6$region == "Northeast"], positive = TRUE)

# 3. Fit different distributions (e.g., Normal, Gamma, Lognormal)
fit_norm_2 <- fitdist(data_6$dif_day_tx_symp[data_6$region == "Northeast"], "norm")
fit_gamma_2 <- fitdist(data_6$dif_day_tx_symp[data_6$region == "Northeast"], "gamma")
fit_lnorm_2 <- fitdist(data_6$dif_day_tx_symp[data_6$region == "Northeast"], "lnorm")
fit_cauchy_2 <- fitdist(data_6$dif_day_tx_symp[data_6$region == "Northeast"], "cauchy")


# 4. Compare AIC values
summary(fit_norm_2)$aic
summary(fit_gamma_2)$aic
summary(fit_lnorm_2)$aic
summary(fit_cauchy_2)$aic

summary(fit_gamma)
plot(fit_gamma)
plot(fit_lnorm)

media_d2 <- log10(median(data_6$dif_day_tx_symp[data_6$region == "Northeast"], na.rm = T))
10^(media_d2)

d2 <- ggplot(
  data_6 %>% filter(region == "Northeast"),
  aes(x = log10(dif_day_tx_symp))) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 30,
    fill = "grey80",
    color = "black") +
  stat_function(
    fun = function(z)
      dgamma(10^z,
             shape = fit_gamma_2$estimate["shape"],
             rate  = fit_gamma_2$estimate["rate"]) * log(10) * 10^z,
    aes(color = "gamma"),
    linewidth = 1) +
  stat_function(
    fun = function(z)
      dlnorm(10^z,
             meanlog = fit_lnorm_2$estimate["meanlog"],
             sdlog   = fit_lnorm_2$estimate["sdlog"]) *
      log(10) * 10^z,
    aes(color = "log-normal"),
    linewidth = 1
  ) +
  geom_vline(xintercept = media_d2, linetype = "dashed", 
             linewidth = 1, color = "black") +
  labs(
    x = expression(bold(log[10]~"(Days since symptom onset)")),
    y = "Density",
    title = "Notheast",
    color = "Distribution",
  ) +
  coord_cartesian(xlim = c(0, log10(500))) +
  scale_x_continuous(
    breaks = log10(break_vals),
    labels = break_vals,
    limits = c(0, log10(500))   # keep or change as needed
  ) +
  scale_y_continuous(limits = c(0, 1.5),
                     n.breaks = 7) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, color = "black", size = 16, face = "bold"),
    axis.title = element_text(size = 14, color = "black", face = "bold"),
    axis.text = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 14, color = "black", face = "bold"),
    legend.text = element_text(size = 12, color = "black"),
    legend.position = c(0.9, 0.9)
  ) +
  scale_color_manual(values = color_distribution)


### North
# 2. Visualize to choose candidates
n <- simukde::find_best_fit(x = data_6$dif_day_tx_symp[data_6$region == "North"], positive = TRUE)


# 3. Fit different distributions (e.g., Normal, Gamma, Lognormal)
fit_norm_3 <- fitdist(data_6$dif_day_tx_symp[data_6$region == "North"], "norm")
fit_gamma_3 <- fitdist(data_6$dif_day_tx_symp[data_6$region == "North"], "gamma")
fit_lnorm_3 <- fitdist(data_6$dif_day_tx_symp[data_6$region == "North"], "lnorm")
fit_cauchy_3 <- fitdist(data_6$dif_day_tx_symp[data_6$region == "North"], "cauchy")


# 4. Compare AIC values
summary(fit_norm_3)$aic
summary(fit_gamma_3)$aic
summary(fit_lnorm_3)$aic
summary(fit_cauchy_3)$aic

summary(fit_lnorm)

# 5. Comparative goodness-of-fit statistics
media_d3 <- log10(median(data_6$dif_day_tx_symp[data_6$region == "North"], na.rm = T))
10^(media_d3)

d3 <- ggplot(
  data_6 %>% filter(region == "North"),
  aes(x = log10(dif_day_tx_symp))
) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 30,
    fill = "grey80",
    color = "black"
  ) +
  stat_function(
    fun = function(z)
      dgamma(10^z,
             shape = fit_gamma_3$estimate["shape"],
             rate  = fit_gamma_3$estimate["rate"]) * log(10) * 10^z,
    aes(color = "gamma"),
    linewidth = 1) +
  stat_function(
    fun = function(z)
      dlnorm(10^z,
             meanlog = fit_lnorm_3$estimate["meanlog"],
             sdlog   = fit_lnorm_3$estimate["sdlog"]) *
      log(10) * 10^z,
    aes(color = "log-normal"),
    linewidth = 1
  ) +
  geom_vline(xintercept = media_d3, linetype = "dashed", 
             linewidth = 1, color = "black") +
  labs(
    x = expression(bold(log[10]~"(Days since symptom onset)")),
    y = "Density",
    title = "North",
    color = "Distribution",
  ) +
  coord_cartesian(xlim = c(0, log10(500))) +
  scale_x_continuous(
    breaks = log10(break_vals),
    labels = break_vals,
    limits = c(0, log10(500))   # keep or change as needed
  ) +
  scale_y_continuous(limits = c(0, 1.5),
                     n.breaks = 7) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, color = "black", size = 16, face = "bold"),
    axis.title = element_text(size = 14, color = "black", face = "bold"),
    axis.text = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 14, color = "black", face = "bold"),
    legend.text = element_text(size = 12, color = "black"),
    legend.position = c(0.9, 0.9)
  ) +
  scale_color_manual(values = color_distribution)

### Southeast
# 2. Visualize to choose candidates
ne <- simukde::find_best_fit(x = data_6$dif_day_tx_symp[data_6$region == "Southeast"], positive = TRUE)


# 3. Fit different distributions (e.g., Normal, Gamma, Lognormal)
fit_norm_4 <- fitdist(data_6$dif_day_tx_symp[data_6$region == "Southeast"], "norm")
fit_gamma_4 <- fitdist(data_6$dif_day_tx_symp[data_6$region == "Southeast"], "gamma")
fit_lnorm_4 <- fitdist(data_6$dif_day_tx_symp[data_6$region == "Southeast"], "lnorm")
fit_cauchy_4 <- fitdist(data_6$dif_day_tx_symp[data_6$region == "Southeast"], "cauchy")


# 4. Compare AIC values
summary(fit_norm_4)$aic
summary(fit_gamma_4)$aic
summary(fit_lnorm_4)$aic
summary(fit_cauchy_4)$aic

summary(fit_lnorm)
summary(fit_gamma)

# 5. Comparative goodness-of-fit statistics
media_d4 <- log10(median(data_6$dif_day_tx_symp[data_6$region == "Southeast"], na.rm = T))
10^(media_d4)

d4 <- ggplot(
  data_6 %>% filter(region == "Southeast"),
  aes(x = log10(dif_day_tx_symp))
) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 30,
    fill = "grey80",
    color = "black"
  ) +
  stat_function(
    fun = function(z)
      dgamma(10^z,
             shape = fit_gamma_4$estimate["shape"],
             rate  = fit_gamma_4$estimate["rate"]) * log(10) * 10^z,
    aes(color = "gamma"),
    linewidth = 1) +
  stat_function(
    fun = function(z)
      dlnorm(10^z,
             meanlog = fit_lnorm_4$estimate["meanlog"],
             sdlog   = fit_lnorm_4$estimate["sdlog"]) *
      log(10) * 10^z,
    aes(color = "log-normal"),
    linewidth = 1
  ) +
  geom_vline(xintercept = media_d4, linetype = "dashed", 
             linewidth = 1, color = "black") +
  labs(
    x = expression(bold(log[10]~"(Days since symptom onset)")),
    y = "Density",
    title = "Southeast",
    color = "Distribution"
  ) +
  coord_cartesian(xlim = c(0, log10(500))) +
  scale_x_continuous(
    breaks = log10(break_vals),
    labels = break_vals,
    limits = c(0, log10(500))   # keep or change as needed
  ) +
  scale_y_continuous(limits = c(0, 1.5),
                     n.breaks = 7) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, color = "black", size = 16, face = "bold"),
    axis.title = element_text(size = 14, color = "black", face = "bold"),
    axis.text = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 14, color = "black", face = "bold"),
    legend.text = element_text(size = 12, color = "black"),
    legend.position = c(0.9, 0.9)
  ) +
  scale_color_manual(values = color_distribution)


### South
# 2. Visualize to choose candidates
ne <- simukde::find_best_fit(x = data_6$dif_day_tx_symp[data_6$region == "South"], positive = TRUE)


# 3. Fit different distributions (e.g., Normal, Gamma, Lognormal)
fit_norm_5 <- fitdist(data_6$dif_day_tx_symp[data_6$region == "South"], "norm")
fit_gamma_5 <- fitdist(data_6$dif_day_tx_symp[data_6$region == "South"], "gamma")
fit_lnorm_5 <- fitdist(data_6$dif_day_tx_symp[data_6$region == "South"], "lnorm")
fit_cauchy_5 <- fitdist(data_6$dif_day_tx_symp[data_6$region == "South"], "cauchy")


# 4. Compare AIC values
summary(fit_norm_5)$aic
summary(fit_gamma_5)$aic
summary(fit_lnorm_5)$aic
summary(fit_cauchy_5)$aic

summary(fit_lnorm)
summary(fit_gamma)

# 5. Comparative goodness-of-fit statistics
media_d5 <- log10(median(data_6$dif_day_tx_symp[data_6$region == "South"], na.rm = T))
10^(media_d5)


d5 <- ggplot(
  data_6 %>% filter(region == "South"),
  aes(x = log10(dif_day_tx_symp))
) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 30,
    fill = "grey80",
    color = "black"
  ) +
  stat_function(
    fun = function(z)
      dgamma(10^z,
             shape = fit_gamma_5$estimate["shape"],
             rate  = fit_gamma_5$estimate["rate"]) *
      log(10) * 10^z,
    aes(color = "gamma"),
    linewidth = 1
  ) +
  stat_function(
    fun = function(z)
      dlnorm(10^z,
             meanlog = fit_lnorm_5$estimate["meanlog"],
             sdlog   = fit_lnorm_5$estimate["sdlog"]) *
      log(10) * 10^z,
    aes(color = "log-normal"),
    linewidth = 1
  ) +
  geom_vline(xintercept = media_d5, linetype = "dashed", 
             linewidth = 1, color = "black") +
  labs(
    x = expression(bold(log[10]~"(Days since symptom onset)")),
    y = "Density",
    title = "South",
    color = "Distribution",
  ) +
  coord_cartesian(xlim = c(0, log10(500))) +
  scale_x_continuous(
    breaks = log10(break_vals),
    labels = break_vals,
    limits = c(0, log10(500))   # keep or change as needed
  ) +
  
  scale_y_continuous(limits = c(0, 1.5),
                     n.breaks = 7) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, color = "black", size = 16, face = "bold"),
    axis.title = element_text(size = 14, color = "black", face = "bold"),
    axis.text = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 14, color = "black", face = "bold"),
    legend.text = element_text(size = 12, color = "black"),
    legend.position = c(0.9, 0.9)
  ) +
  scale_color_manual(values = color_distribution)

ggarrange(d0 + rremove("xlab"),
          d1 + rremove("xlab") + rremove("ylab"),
          d3 + rremove("xlab") + rremove("ylab"),
          d2,
          d5 + rremove("ylab"),
          d4 + rremove("ylab"), 
          common.legend = T, 
          legend = "bottom")
ggsave(path = here("figures and tables"),
       width = 12, height = 8, dpi = 300,
       file = "regional_delay_distribution.png")

Region <- c("Overall", "Center-West", "North", "Northeast", "Sout", "Southeast")
Normal <- c(fit_norm$aic, fit_norm_1$aic, fit_norm_3$aic, fit_norm_2$aic, fit_norm_5$aic, fit_norm_4$aic)
log_normal <- c(fit_lnorm$aic, fit_lnorm_1$aic, fit_lnorm_3$aic, fit_lnorm_2$aic, fit_lnorm_5$aic, fit_lnorm_4$aic)
gamma <- c(fit_gamma$aic, fit_gamma_1$aic, fit_gamma_3$aic, fit_gamma_2$aic, fit_gamma_5$aic, fit_gamma_4$aic)
cauchy <- c(fit_cauchy$aic, fit_cauchy_1$aic, fit_cauchy_3$aic, fit_cauchy_2$aic, fit_cauchy_5$aic, fit_cauchy_4$aic)
lmedia <- c(fit_lnorm$estimate["meanlog"], fit_lnorm_1$estimate["meanlog"], 
            fit_lnorm_3$estimate["meanlog"], fit_lnorm_2$estimate["meanlog"], 
            fit_lnorm_5$estimate["meanlog"], fit_lnorm_4$estimate["meanlog"])
lsd <- c(fit_lnorm$estimate["sdlog"], fit_lnorm_1$estimate["sdlog"], 
            fit_lnorm_3$estimate["sdlog"], fit_lnorm_2$estimate["sdlog"], 
            fit_lnorm_5$estimate["sdlog"], fit_lnorm_4$estimate["sdlog"])
gamma_shape <- c(fit_gamma$estimate["shape"], fit_gamma_1$estimate["shape"], 
                 fit_gamma_3$estimate["shape"], fit_gamma_2$estimate["shape"], 
                 fit_gamma_5$estimate["shape"], fit_gamma_4$estimate["shape"])
gamma_rate <- c(fit_gamma$estimate["rate"], fit_gamma_1$estimate["rate"], 
                 fit_gamma_3$estimate["rate"], fit_gamma_2$estimate["rate"], 
                 fit_gamma_5$estimate["rate"], fit_gamma_4$estimate["rate"])


table_aic <- data.frame(Region = Region,
                        Normal = Normal,
                        `log-normal` = log_normal,
                        Gamma = gamma,
                        Cauchy = cauchy, 
                        lnorm_media = lmedia,
                        lnorm_sd = lsd,
                        gamma_s = gamma_shape,
                        gamma_r = gamma_rate) 

# ---- ensure the log-normal column name is exactly "log-normal" ----
if ("log.normal" %in% names(table_aic)) {
  names(table_aic)[names(table_aic) == "log.normal"] <- "log-normal"
  names(table_aic)[names(table_aic) == "lnorm_media"] <- "log media"
  names(table_aic)[names(table_aic) == "lnorm_sd"] <- "log sd"
  names(table_aic)[names(table_aic) == "gamma_s"] <- "Gamma shape"
  names(table_aic)[names(table_aic) == "gamma_r"] <- "Gamma rate"
}
# (if you already have `log-normal` as column name, this does nothing)

# ---- prepare a display copy (round numeric values for nicer printing) ----
table_disp <- table_aic %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

# ---- build flextable ----
ft <- flextable(table_disp)

# numeric columns (exclude Region)
num_cols <- setdiff(names(table_disp), "Region")

# define the allowed column index range (2 through 5)
# NOTE: these are positions in the data.frame (1 = Region)
start_col <- 2
end_col   <- 5

# ensure we don't exceed available columns
ncols <- ncol(table_disp)
allowed_positions <- seq(start_col, min(end_col, ncols))
allowed_cols <- names(table_disp)[allowed_positions]

# loop rows and bold the minimum(s) only if they lie in allowed_cols
for (r in seq_len(nrow(table_disp))) {
  # row values for numeric columns (in the same order as num_cols)
  row_vals <- as.numeric(table_disp[r, num_cols])
  if (all(is.na(row_vals))) next
  min_val <- min(row_vals, na.rm = TRUE)
  # candidate cols (names) where value equals the min (handles ties)
  candidate_cols <- num_cols[which(!is.na(row_vals) & row_vals == min_val)]
  # restrict to allowed columns (intersection)
  cols_to_bold <- intersect(candidate_cols, allowed_cols)
  # apply bold if any
  if (length(cols_to_bold) > 0) {
    ft <- bold(ft, i = r, j = cols_to_bold, bold = TRUE)
  }
}

# final formatting and display
ft <- autofit(ft)

# show it
ft %>% 
  save_as_docx(path = here("figures and tables", "table_aic_distribution.docx"))
#############################################################################################################################


#This is to summarise the Area Under the Curve of the best fitting distribution: Log-Normal. 

#-------------------------------------------------------------------------------------
# Using the best fitting Log Normal to describe/summarise the delay characteristics
#-------------------------------------------------------------------------------------

# Parameters

mu <- 3.270
sigma <- 1.130

# Time grid (days)

t <- seq(0.01, 400, length.out = 1000)

# Log-normal density

dens <- dlnorm(t, meanlog = mu, sdlog = sigma)

# Plot

plot(
  t, dens, type = "l", lwd = 2,
  xlab = "Time (days)",
  ylab = "Density",
  main = "Log-normal Distribution (meanlog = 3.27, sdlog = 1.13)"
)


## These are P(T > x), i.e. the survival function:

p_gt_30 <- 1 - plnorm(30, meanlog = mu, sdlog = sigma)
p_gt_60 <- 1 - plnorm(60, meanlog = mu, sdlog = sigma)
p_gt_90 <- 1 - plnorm(90, meanlog = mu, sdlog = sigma)
p_gt_120 <- 1 - plnorm(120, meanlog = mu, sdlog = sigma)
p_gt_150 <- 1 - plnorm(150, meanlog = mu, sdlog = sigma)
p_gt_180 <- 1 - plnorm(180, meanlog = mu, sdlog = sigma)
p_gt_30
p_gt_60
p_gt_90
p_gt_120
p_gt_150
p_gt_180

# Quantiles
# Days corresponding to cumulative AUC levels
# These are quantiles, i.e. values t such that
# P(T ≤ t) = 0.80, 0.90, 0.95

q_50 <- qlnorm(0.50, meanlog = mu, sdlog = sigma)
q_75 <- qlnorm(0.75, meanlog = mu, sdlog = sigma)
q_80 <- qlnorm(0.80, meanlog = mu, sdlog = sigma)
q_90 <- qlnorm(0.90, meanlog = mu, sdlog = sigma)
q_95 <- qlnorm(0.95, meanlog = mu, sdlog = sigma)
q_50
q_75
q_80
q_90
q_95

# END


