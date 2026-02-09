library(dplyr)
library(ggplot2)
library(patchwork)


# select file
leps_exp      <- -53
N_simulations <- 1000

file_name = sprintf("results/sbc_N%i_leps%.0f.csv", N_simulations, leps_exp)

sbc_data <- read.csv(file_name)


# --- Plotting the SBC Diagnostics (3-Panel Figure) ---
plot_sbc_diagnostics <- function(ranks, L_SAMPLES) {
  
  N_RANKS <- length(ranks)
  N_BINS <- 20
  
  # --- Rank Histogram (Plot A) ---
  expected_count <- N_RANKS / N_BINS
  alpha <- 0.05
  ci_hist <- qbinom(c(alpha/2, 1 - alpha/2), size = N_RANKS, prob = 1/N_BINS)
  
  rank_df <- data.frame(rank = ranks)
  
  p_hist <- ggplot(rank_df, aes(x = rank)) +
    geom_histogram(bins = N_BINS, color = "black", fill = "gray40", boundary = 0) +
    annotate("rect", xmin = 0, xmax = 200,
             ymin = ci_hist[1], ymax = ci_hist[2],
             fill = "skyblue", alpha = 0.2) +
    geom_hline(yintercept = expected_count, color = "blue", linewidth = 0.8) +
    labs(
      title = "a) Rank Histogram",
      x = "Rank Statistic",
      y = "Count"
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
  
  # --- ECDF Plot (Plot B) ---
  ecdf_data <- ecdf(ranks / L_SAMPLES)
  ecdf_df <- data.frame(
    rank_norm = sort(ranks / L_SAMPLES),
    ecdf_val = ecdf_data(sort(ranks / L_SAMPLES))
  )
  
  # Kolmogorov-Smirnov critical value for 95% CI
  ks_crit <- 1.36 / sqrt(N_RANKS)
  
  # Data for the confidence band
  band_data <- data.frame(
    rank_norm = seq(0, 1, length.out = 100),
    y_lower = pmax(0, seq(0, 1, length.out = 100) - ks_crit),
    y_upper = pmin(1, seq(0, 1, length.out = 100) + ks_crit)
  )
  
  p_ecdf <- ggplot(ecdf_df, aes(x = rank_norm, y = ecdf_val)) +
    geom_ribbon(data = band_data, aes(x = rank_norm, ymin = y_lower, ymax = y_upper), 
                fill = "skyblue", alpha = 0.3, inherit.aes = FALSE) + # KS Confidence Band
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", linewidth = 0.8) +
    geom_step(color = "blue", linewidth = 0.8) +
    labs(
      title = "b) ECDF",
      x = "Normalized Rank",
      y = "Cumulative Probability"
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
  
  # --- ECDF Difference Plot (Plot C) ---
  ecdf_diff_df <- ecdf_df %>%
    mutate(diff = ecdf_val - rank_norm)
  
  # --- Elliptical Confidence Region (Plot C) ---
  z_crit <- qnorm(0.975)
  
  # Calculate the boundary of the elliptical region
  elliptical_band_data <- data.frame(
    rank_norm = seq(0, 1, length.out = 200)
  ) %>%
    mutate(
      # The standard deviation of the difference is sqrt(x(1-x) / N_RANKS)
      std_dev = sqrt(rank_norm * (1 - rank_norm) / N_RANKS),
      y_upper = z_crit * std_dev,
      y_lower = -z_crit * std_dev
    )
  
  p_diff <- ggplot(ecdf_diff_df, aes(x = rank_norm, y = diff)) +
    geom_ribbon(data = elliptical_band_data, aes(x = rank_norm, ymin = y_lower, ymax = y_upper), 
                fill = "skyblue", alpha = 0.2, inherit.aes = FALSE) + # Elliptical Confidence Band
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
    geom_step(color = "blue", linewidth = 0.8) +
    labs(
      title = "c) ECDF Difference",
      x = "Normalized Rank",
      y = "Difference"
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(-max(elliptical_band_data$y_upper) * 1.5, max(elliptical_band_data$y_upper) * 1.5)) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
  
  # --- Combine Plots ---
  final_plot <- 
    p_hist / (p_ecdf + p_diff) +
    plot_annotation(
      theme = theme(plot.title = element_text(size = 16, face = "plain"))
    )
  
  return(final_plot)
}

print(plot_sbc_diagnostics(sbc_data$rank_mu, 200))
