# --- Libraries ---
library(cmdstanr)
library(posterior)
library(ggplot2)
library(dplyr)

# --- SBC Parameters ---
STAN_FILE <- "normal_normal_sbc.stan"
N_SIMS <- 1000       # Number of SBC runs (simulated datasets)
N_OBS <- 10          # Number of observations per dataset
L_SAMPLES <- 200     # Number of posterior samples (L) for rank calculation
N_CHAINS <- 1        # Number of chains per run
N_WARMUP <- 1000      # Warmup iterations

# Prior parameters for the simulation and the Stan model
PRIOR_MU_MEAN <- 0
PRIOR_MU_SD <- 1
SIGMA_Y <- 1         # Known data-level standard deviation

# --- Compile Stan Model ---
tryCatch({
  mod <- cmdstan_model(STAN_FILE)
}, error = function(e) {
  message("Error compiling Stan model. Ensure cmdstanr is installed and configured, and the Stan file exists.")
  stop(e)
})

# --- SBC Simulation Function ---
run_sbc_normal_normal <- function(mod, N_SIMS, N_OBS, L_SAMPLES, N_CHAINS, N_WARMUP) {
  
  ranks <- numeric(N_SIMS)
  
  for (s in 1:N_SIMS) {
    # Simulate true parameter (mu) from the prior
    mu_true <- rnorm(1, PRIOR_MU_MEAN, PRIOR_MU_SD)
    
    # Simulate data (y) from the likelihood
    y <- rnorm(N_OBS, mu_true, SIGMA_Y)
    
    data_list <- list(
      N = N_OBS, 
      y = y, 
      sigma_y = SIGMA_Y,
      prior_mu_mean = PRIOR_MU_MEAN, 
      prior_mu_sd = PRIOR_MU_SD
    )
    
    # Run inference
    fit <- mod$sample(
      data = data_list, 
      chains = N_CHAINS,
      iter_sampling = L_SAMPLES, 
      iter_warmup = N_WARMUP, 
      refresh = 0, 
      show_messages = FALSE,
      sig_figs = 10,
      save_warmup = FALSE
    )
    
    # Check for sampling errors
    if (fit$return_codes() != 0) {
      warning(paste("Sampling failed for simulation", s))
      ranks[s] <- NA
      next
    }
    
    # Calculate rank statistic
    mu_draws <- as_draws_matrix(fit$draws(variable = "mu"))[, 1]
    ranks[s] <- sum(mu_draws < mu_true)
  }
  
  return(ranks[!is.na(ranks)])
}

# --- Run SBC ---
cat(paste("Starting SBC run with", N_SIMS, "simulations...\n"))
ranks <- run_sbc_normal_normal(mod, N_SIMS, N_OBS, L_SAMPLES, N_CHAINS, N_WARMUP)
cat(paste("SBC run complete. Successfully ran", length(ranks), "simulations.\n"))

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
  final_plot <- (p_hist) / (p_ecdf + p_diff) + 
    plot_annotation() & theme(plot.title = element_text(size = 16, face = "plain"))
  
  return(final_plot)
}

print(plot_sbc_diagnostics(ranks, L_SAMPLES))
