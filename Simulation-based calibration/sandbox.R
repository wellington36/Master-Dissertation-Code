# --- Libraries ---
library(ggplot2)
library(patchwork)
library(dplyr)

set.seed(123)

# --------------------------------------------------------------------------- #
# 1. Define 4 conceptual scenarios (matching Mendes et al. Fig. 6)
# --------------------------------------------------------------------------- #
scenarios <- list(
  list(name = "(a) Correct",        prior_mean = 0, prior_sd = 1,
       post_mean = 0, post_sd = 1,   rank_pattern = "uniform"),
  list(name = "(b) Under-dispersed",prior_mean = 0, prior_sd = 1,
       post_mean = 0, post_sd = 0.3, rank_pattern = "under"),
  list(name = "(c) Over-dispersed", prior_mean = 0, prior_sd = 1,
       post_mean = 0, post_sd = 2,   rank_pattern = "over"),
  list(name = "(d) Biased",         prior_mean = 0, prior_sd = 1,
       post_mean = 2, post_sd = 1,   rank_pattern = "bias")
)

# --------------------------------------------------------------------------- #
# 2. Simulate rank patterns (approximate shapes of SBC results)
# --------------------------------------------------------------------------- #
simulate_ranks <- function(n = 1000, pattern = c("uniform", "under", "over", "bias")) {
  pattern <- match.arg(pattern)
  if (pattern == "uniform") {
    runif(n)
  } else if (pattern == "under") {
    # under-dispersed → more ranks near 0 or 1 (horns)
    c(rbeta(n/2, 0.4, 2), rbeta(n/2, 2, 0.4))
  } else if (pattern == "over") {
    # over-dispersed → central bulge
    rbeta(n, 2, 2)
  } else if (pattern == "bias") {
    # biased → skewed to one side
    rbeta(n, 0.5, 2)
  }
}

# --------------------------------------------------------------------------- #
# 3. Function to build one row of 3 plots
# --------------------------------------------------------------------------- #
make_row <- function(scen) {
  # --- Prior & Posterior ---
  theta_vals <- seq(-5, 5, length.out = 400)
  df_prior <- data.frame(
    theta = theta_vals,
    density = dnorm(theta_vals, scen$prior_mean, scen$prior_sd),
    Type = "Prior"
  )
  df_post <- data.frame(
    theta = theta_vals,
    density = dnorm(theta_vals, scen$post_mean, scen$post_sd),
    Type = "Posterior"
  )
  p_prior <- ggplot() +
    geom_area(data = df_prior, aes(theta, density, fill = Type),
              alpha = 0.4, color = "black") +
    geom_area(data = df_post, aes(theta, density, fill = Type),
              alpha = 0.6, color = "black") +
    scale_fill_manual(values = c("Prior" = "skyblue", "Posterior" = "navy")) +
    coord_cartesian(xlim = c(-5, 5), ylim = c(0, 0.5)) +
    theme_minimal(base_size = 12) +
    labs(x = expression(theta), y = "Density", title = scen$name, fill = NULL) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "top")
  
  # --- Simulate ranks and histogram ---
  ranks <- simulate_ranks(1000, scen$rank_pattern)
  rank_df <- data.frame(rank = ranks)
  B <- 20
  expected <- length(ranks) / B
  sd_count <- sqrt(length(ranks) * (1/B) * (1 - 1/B))
  ci_low <- expected - 1.96 * sd_count
  ci_high <- expected + 1.96 * sd_count
  
  p_hist <- ggplot(rank_df, aes(x = rank)) +
    geom_histogram(bins = B, fill = "gray40", color = "black", boundary = 0) +
    geom_hline(yintercept = expected, color = "black", linewidth = 0.8) +
    annotate("rect", xmin = -Inf, xmax = Inf,
             ymin = ci_low, ymax = ci_high,
             fill = "skyblue", alpha = 0.2) +
    theme_minimal(base_size = 12) +
    labs(x = "Rank", y = "Count", title = "Rank Histogram") +
    theme(plot.title = element_text(hjust = 0.5))
  
  # --- ECDF ---
  p_ecdf <- ggplot(rank_df, aes(rank)) +
    stat_ecdf(geom = "step", color = "black") +
    geom_abline(slope = 1, intercept = 0, color = "skyblue3", linewidth = 1) +
    theme_minimal(base_size = 12) +
    labs(x = "Normalized Rank", y = "Cumulative Probability",
         title = "ECDF Plot") +
    theme(plot.title = element_text(hjust = 0.5))
  
  # Combine horizontally
  p_prior | p_hist | p_ecdf
}

# --------------------------------------------------------------------------- #
# 4. Combine all four rows into a 4×3 figure
# --------------------------------------------------------------------------- #
rows <- lapply(scenarios, make_row)
final_plot <- wrap_plots(rows, ncol = 1)
final_plot
