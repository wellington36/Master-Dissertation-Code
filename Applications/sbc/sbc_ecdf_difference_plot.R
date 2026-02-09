library(dplyr)
library(ggplot2)
library(patchwork)

N_simulations <- 1000
L_SAMPLES     <- 200

plot_ecdf_difference <- function(ranks, L_SAMPLES, leps_exp) {
  
  N_RANKS <- length(ranks)
  
  # Normalised ranks
  rank_norm <- sort(ranks / L_SAMPLES)
  
  ecdf_val <- ecdf(rank_norm)(rank_norm)
  
  ecdf_diff_df <- data.frame(
    rank_norm = rank_norm,
    diff = ecdf_val - rank_norm
  )
  
  # Elliptical confidence region
  z_crit <- qnorm(0.975)
  
  band_df <- data.frame(
    rank_norm = seq(0, 1, length.out = 200)
  ) %>%
    mutate(
      std_dev = sqrt(rank_norm * (1 - rank_norm) / N_RANKS),
      y_upper =  z_crit * std_dev,
      y_lower = -z_crit * std_dev
    )
  
  ggplot(ecdf_diff_df, aes(x = rank_norm, y = diff)) +
    geom_ribbon(
      data = band_df,
      aes(
        x = rank_norm,
        ymin = y_lower,
        ymax = y_upper
      ),
      fill = "skyblue",
      alpha = 0.25,
      inherit.aes = FALSE
    ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      color = "black",
      linewidth = 0.8
    ) +
    geom_step(
      color = "blue",
      linewidth = 0.8
    ) +
    labs(
      title = sprintf("ECDF Difference (log(error) = %d)", leps_exp),
      x = "Normalised Rank",
      y = "ECDF − Uniform"
    ) +
    coord_cartesian(
      xlim = c(0, 1),
      ylim = c(
        -1.5 * max(band_df$y_upper),
        1.5 * max(band_df$y_upper)
      )
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
}

plot_ecdf_differences_for_N <- function(N_simulations,
                                        L_SAMPLES = 200,
                                        leps_vec = c(-53, -32, -16, -8, -4, -2)) {
  
  plots <- lapply(leps_vec, function(leps_exp) {
    
    file_name <- sprintf(
      "results/sbc_N%i_leps%.0f.csv",
      N_simulations,
      leps_exp
    )
    
    sbc_data <- read.csv(file_name)
    
    plot_ecdf_difference(
      ranks = sbc_data$rank_mu,
      L_SAMPLES = L_SAMPLES,
      leps_exp = leps_exp
    )
  })
  
  wrap_plots(plots, ncol = 2) +
    plot_annotation()
}

final_plot <- plot_ecdf_differences_for_N(
  N_simulations = N_simulations,
  L_SAMPLES = L_SAMPLES
)

print(final_plot)
