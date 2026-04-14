library(dplyr)
library(ggplot2)
library(patchwork)

N_simulations <- 1000

plot_rank_histogram <- function(ranks, leps_exp, N_bins = 20) {
  
  N_RANKS <- length(ranks)
  expected_count <- N_RANKS / N_bins
  alpha <- 0.05
  
  ci_hist <- qbinom(
    c(alpha / 2, 1 - alpha / 2),
    size = N_RANKS,
    prob = 1 / N_bins
  )
  
  rank_df <- data.frame(rank = ranks)
  
  ggplot(rank_df, aes(x = rank)) +
    geom_histogram(
      bins = N_bins,
      color = "black",
      fill = "gray40",
      boundary = 0
    ) +
    annotate(
      "rect",
      xmin = 0,
      xmax = 200,
      ymin = ci_hist[1],
      ymax = ci_hist[2],
      fill = "skyblue",
      alpha = 0.25
    ) +
    geom_hline(
      yintercept = expected_count,
      color = "blue",
      linewidth = 0.8
    ) +
    labs(
      title = sprintf("ECDF Difference (log(error) = %d)", leps_exp),
      x = "Rank Statistic",
      y = "Count"
    ) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
}

plot_rank_histograms_for_N <- function(N_simulations,
                                       leps_vec = c(-53, -32, -16, -8, -4, -2),
                                       L_SAMPLES = 200) {
  
  all_data <- lapply(leps_vec, function(leps_exp) {
    file_name <- sprintf(
      "results/sbc_N%i_leps%.0f.csv",
      N_simulations,
      leps_exp
    )
    read.csv(file_name)
  })
  
  # Compute global max count across all histograms
  max_count <- max(sapply(all_data, function(df) {
    hist(df$rank_mu, breaks = 20, plot = FALSE)$counts
  }))
  
  plots <- mapply(function(sbc_data, leps_exp) {
    
    plot_rank_histogram(
      ranks = sbc_data$rank_mu,
      leps_exp = leps_exp
    ) +
      ylim(0, max_count * 1.10)  # same y-axis for all
    
  }, all_data, leps_vec, SIMPLIFY = FALSE)
  
  wrap_plots(plots, ncol = 2) +
    plot_annotation()
}

final_plot <- plot_rank_histograms_for_N(N_simulations)
print(final_plot)
