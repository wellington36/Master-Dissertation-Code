# --- Libraries ---
library(cmdstanr)
library(ggplot2)
library(patchwork)
library(dplyr)

# Compile model --------------------------------------------------------------
mod <- cmdstan_model("SBC/simple_normal_model.stan")

set.seed(123)
n_sims <- 500        # number of simulated datasets
n_obs  <- 5          # number of observations per dataset
sigma  <- 1

# Four SBC scenarios (as in Mendes Fig. 6)
scenarios <- list(
  list(prior_mean = 0, prior_sd = 1,  label = "(a) Prior N(0,1)  — Correct"),
  list(prior_mean = 0, prior_sd = 0.2,label = "(b) Prior N(0,0.2) — Under-dispersed"),
  list(prior_mean = 0, prior_sd = 20,  label = "(c) Prior N(0,5)  — Over-dispersed"),
  list(prior_mean = 5, prior_sd = 1,  label = "(d) Prior N(2,1)  — Biased")
)

# --------------------------------------------------------------------------- #
# Helper: run one SBC scenario
# --------------------------------------------------------------------------- #
run_sbc <- function(prior_mean, prior_sd, n_sims, n_obs, sigma = 1) {
  ranks <- numeric(n_sims)
  all_draws <- numeric()
  
  for (s in 1:n_sims) {
    # --- Simulate true parameter and data
    theta_true <- rnorm(1, 0, 1)   # true always from N(0,1)
    y <- rnorm(n_obs, theta_true, sigma)
    data_list <- list(
      N = n_obs, y = y, sigma = sigma,
      prior_mean = prior_mean, prior_sd = prior_sd
    )
    
    # --- Posterior inference
    fit <- mod$sample(
      data = data_list, chains = 1, iter_warmup = 300,
      iter_sampling = 200, refresh = 0, show_messages = FALSE
    )
    
    draws <- fit$draws(variable = "theta", format = "matrix")[, 1]
    ranks[s] <- sum(draws < theta_true)
    all_draws <- c(all_draws, draws)
  }
  
  list(ranks = ranks, post_draws = all_draws)
}

# --------------------------------------------------------------------------- #
# Plot builder for one scenario
# --------------------------------------------------------------------------- #
plot_sbc_case <- function(res, prior_mean, prior_sd, label, n_sims) {
  ranks <- res$ranks
  post_draws <- res$post_draws
  
  # ----- Prior vs Replicate-averaged Posterior (RAP)
  theta_vals <- seq(-5, 5, length.out = 400)
  df_prior <- data.frame(theta = theta_vals,
                         density = dnorm(theta_vals, 0, 1),
                         Type = "True Prior N(0,1)")
  df_post <- data.frame(theta = post_draws, Type = "Replicate-averaged Posterior")
  
  p_prior <- ggplot() +
    geom_density(data = df_post, aes(x = theta, fill = Type),
                 color = "black", alpha = 0.5, adjust = 1.2) +
    geom_line(data = df_prior, aes(x = theta, y = density),
              color = "skyblue4", linewidth = 1.1) +
    coord_cartesian(xlim = c(-5, 5), ylim = c(0, 0.5)) +
    scale_fill_manual(values = c("Replicate-averaged Posterior" = "navy")) +
    theme_minimal(base_size = 12) +
    labs(x = expression(theta), y = "Density", title = label, fill = NULL) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "top")
  
  # ----- Rank histogram
  B <- 20
  expected <- n_sims / B
  sd_count <- sqrt(n_sims * (1 / B) * (1 - 1 / B))
  ci_low <- expected - 1.96 * sd_count
  ci_high <- expected + 1.96 * sd_count
  rank_df <- data.frame(rank = ranks)
  
  p_hist <- ggplot(rank_df, aes(x = rank)) +
    geom_histogram(bins = B, fill = "gray40", color = "black", boundary = 0) +
    geom_hline(yintercept = expected, color = "black", linewidth = 0.8) +
    annotate("rect", xmin = -Inf, xmax = Inf,
             ymin = ci_low, ymax = ci_high,
             fill = "skyblue", alpha = 0.2) +
    theme_minimal(base_size = 12) +
    labs(x = "Rank", y = "Count", title = "Rank Histogram") +
    theme(plot.title = element_text(hjust = 0.5))
  
  # ----- ECDF of normalized ranks
  p_ecdf <- ggplot(rank_df, aes(rank / max(rank))) +
    stat_ecdf(geom = "step", color = "black") +
    geom_abline(slope = 1, intercept = 0, color = "skyblue3", linewidth = 1) +
    theme_minimal(base_size = 12) +
    labs(x = "Normalized Rank", y = "Cumulative Probability",
         title = "ECDF Plot") +
    theme(plot.title = element_text(hjust = 0.5))
  
  p_prior | p_hist | p_ecdf
}

# --------------------------------------------------------------------------- #
# Run all four scenarios
# --------------------------------------------------------------------------- #
plot_list <- list()
for (sc in scenarios) {
  cat("Running", sc$label, "...\n")
  res <- run_sbc(sc$prior_mean, sc$prior_sd, n_sims, n_obs, sigma)
  plot_list[[sc$label]] <-
    plot_sbc_case(res, sc$prior_mean, sc$prior_sd, sc$label, n_sims)
}

# Combine in a single 4×3 figure (one row per scenario)
final_plot <- wrap_plots(plot_list, ncol = 1)
final_plot
