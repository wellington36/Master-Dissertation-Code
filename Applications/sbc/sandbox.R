library(dplyr)
library(ggplot2)
library(purrr)

N_simulations <- 1000
leps_vec <- c(-53, -32, -16, -8, -4, -2)

# Read and combine all ESS values
ess_all <- map_dfr(leps_vec, function(leps_exp) {
  
  file_name <- sprintf(
    "results/sbc_N%i_leps%.0f.csv",
    N_simulations,
    leps_exp
  )
  
  df <- read.csv(file_name)
  
  df %>%
    transmute(
      leps = factor(leps_exp, levels = leps_vec),
      ess_per_sec = ess_per_sec
    )
})

# Sanity check plot
ggplot(ess_all, aes(x = leps, y = ess_per_sec)) +
  geom_boxplot(
    outlier.shape = NA,
    fill = "skyblue",
    alpha = 0.5
  ) +
  geom_jitter(
    width = 0.15,
    alpha = 0.3,
    size = 1
  ) +
  labs(
    title = "Sanity check: ESS per second distribution per log(error)",
    x = "log(error)",
    y = "ESS per sec"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold")
  )
