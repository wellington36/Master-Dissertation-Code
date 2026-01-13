library(dplyr)
library(ggplot2)
library(patchwork)

source("com_poisson_pmf.R")
source("rcomp_exact.R")
source("rcomp_rou.R")
source("rcomp_rejection_benson.R")

# Dummy sampler
rcomp_dummy <- function(n, lambda, nu) {
  rpois(n, lambda)
}

# --- Parameters ---
lambda <- 2
n_sample <- 10000

# --- Global style maps ---
method_colors <- c(
  "True PMF" = "black",
  "Rejection RoU" = "blue",
  "Rejection Benson" = "orange",
  "Exact" = "green",
  "Poisson" = "red"
)

method_linetypes <- c(
  "True PMF" = "solid",
  "Rejection RoU" = "dashed",
  "Rejection Benson" = "dashed",
  "Exact" = "dashed",
  "Poisson" = "dotted"
)

make_plot_pair <- function(nu, methods_to_use) {
  
  log_lambda <- log(lambda)
  
  # --- Sampling ---
  samples_list <- list()
  if ("rou" %in% methods_to_use)
    samples_list$rou <- rcomp_rou(n_sample, lambda, nu)
  if ("rej_ben" %in% methods_to_use)
    samples_list$rej_ben <- rcomp_rejection_benson(n_sample, lambda, nu)
  if ("exa" %in% methods_to_use)
    samples_list$exa <- rcomp_exact(n_sample, lambda, nu)
  if ("dummy" %in% methods_to_use)
    samples_list$dummy <- rcomp_dummy(n_sample, lambda, nu)
  
  # --- True PMF ---
  K <- max(unlist(samples_list)) + 10
  support <- 0:K
  log_p <- sapply(
    support,
    function(y) com_poisson_log_lpmf(y, log_lambda, nu, 1e-12)
  )
  p_true <- exp(log_p) / sum(exp(log_p))
  
  df_true <- data.frame(
    x = support,
    y = p_true,
    Method = "True PMF"
  )
  
  # --- Empirical frequencies ---
  make_freq_df <- function(samples, name) {
    tbl <- as.data.frame(table(x = samples))
    data.frame(
      x = as.numeric(as.character(tbl$x)),
      y = tbl$Freq / sum(tbl$Freq),
      Method = name
    )
  }
  
  df_plot <- list(df_true)
  if ("rou" %in% methods_to_use)
    df_plot <- c(df_plot, list(make_freq_df(samples_list$rou, "Rejection RoU")))
  if ("rej_ben" %in% methods_to_use)
    df_plot <- c(df_plot, list(make_freq_df(samples_list$rej_ben, "Rejection Benson")))
  if ("exa" %in% methods_to_use)
    df_plot <- c(df_plot, list(make_freq_df(samples_list$exa, "Exact")))
  if ("dummy" %in% methods_to_use)
    df_plot <- c(df_plot, list(make_freq_df(samples_list$dummy, "Poisson")))
  
  df_plot <- do.call(rbind, df_plot)
  
  # --- Deviations ---
  df_diff <- df_plot %>%
    left_join(
      df_true %>% select(x, y_true = y),
      by = "x"
    ) %>%
    mutate(diff = y - y_true) %>%
    filter(Method != "True PMF")
  
  
  if (nu == 1) {
    # --- PMF plot ---
    p1 <- ggplot(df_plot, aes(x, y, color = Method, linetype = Method)) +
      geom_point(size = 1.5) +
      geom_line(linewidth = 0.5) +
      guides(
        color = guide_legend(nrow = 1),
        linetype = guide_legend(nrow = 1)
      ) + 
      scale_color_manual(values = method_colors) +
      scale_linetype_manual(values = method_linetypes) +
      labs(
        title = paste0("λ = ", lambda, ", ν = ", nu, " (PMF)"),
        x = "y",
        y = "Probability / Rel. Frequency"
      ) +
      theme_minimal() + 
      theme(legend.title = element_blank())
  } else {
    # --- PMF plot ---
    p1 <- ggplot(df_plot, aes(x, y, color = Method, linetype = Method)) +
      geom_point(size = 1.5) +
      geom_line(linewidth = 0.5) +
      guides(
        color = guide_legend(nrow = 1),
        linetype = guide_legend(nrow = 1)
      ) + 
      scale_color_manual(values = method_colors) +
      scale_linetype_manual(values = method_linetypes) +
      labs(
        title = paste0("λ = ", lambda, ", ν = ", nu, " (PMF)"),
        x = "y",
        y = "Probability / Rel. Frequency"
      ) +
      theme_minimal() +
      theme(legend.position = "none")
  }
    
  # --- Difference plot ---
  p2 <- ggplot(df_diff, aes(x, diff, color = Method, linetype = Method)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_point(size = 1.5) +
    geom_line(linewidth = 0.5) +
    scale_color_manual(values = method_colors) +
    scale_linetype_manual(values = method_linetypes) +
    labs(
      title = paste0("λ = ", lambda, ", ν = ", nu, " (Deviation)"),
      x = "y",
      y = "Empirical − True"
    ) +
    theme_minimal() +
    theme(legend.position = "none")
  
  list(p1 = p1, p2 = p2)
}

# --- Build grid ---
plots <- list(
  plots_05 <- make_plot_pair(0.5, c("rou", "rej_ben", "exa")),
  plots_1  <- make_plot_pair(1,   c("rou", "rej_ben", "exa", "dummy")),
  plots_2  <- make_plot_pair(2,   c("rou", "rej_ben", "exa"))
)

# --- Final combined plot ---
final_plot <-
  (plots_05$p1 | plots_05$p2) /
  (plots_1$p1  | plots_1$p2 ) /
  (plots_2$p1  | plots_2$p2 ) +
  plot_layout(guides = "collect") +
  plot_annotation(
    theme = theme(
      legend.position = "top",
      legend.direction = "horizontal",
      legend.title = element_blank()
    )
  )

final_plot
