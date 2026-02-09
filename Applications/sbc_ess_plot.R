library(dplyr)
library(ggplot2)
library(tidyr)

N_simulations <- 1000
leps_vec <- c(-53, -32, -16, -8, -4, -2)

ess_summary <- lapply(leps_vec, function(leps_exp) {
  
  file_name <- sprintf(
    "sbc_N%i_leps%.0f.csv",
    N_simulations,
    leps_exp
  )
  
  df <- read.csv(file_name)
  
  tibble(
    leps = leps_exp,
    very_bad = sum(df$ess_mean < 20, na.rm = TRUE),
    bad      = sum(df$ess_mean >= 20 & df$ess_mean < 40, na.rm = TRUE),
    ok       = sum(df$ess_mean >= 40 & df$ess_mean < 100, na.rm = TRUE),
    good     = sum(df$ess_mean >= 100, na.rm = TRUE)
  )
}) %>%
  bind_rows()

ess_long <- ess_long %>%
  mutate(
    category = factor(
      category,
      levels = c("very_bad", "bad", "ok", "good")
    )
  )

ggplot(ess_long, aes(x = leps, y = count, color = category)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  scale_color_manual(
    values = c(
      very_bad = "red",
      bad      = "orange",
      ok       = "darkgreen",
      good     = "blue"
    ),
    labels = c(
      very_bad = "Very bad (ESS < 20)",
      bad      = "Bad (20 ≤ ESS < 40)",
      ok       = "OK (40 ≤ ESS < 100)",
      good     = "Good (100 ≤ ESS)"
    )
  ) +
  scale_x_continuous(breaks = leps_vec) +
  labs(
    title = "ESS quality vs log(error)",
    x = "log(error)",
    y = "Number of simulations",
    color = "ESS category"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "top"
  )
