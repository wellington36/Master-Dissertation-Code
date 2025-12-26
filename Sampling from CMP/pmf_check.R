library(dplyr)
library(ggplot2)
library(patchwork)


#set.seed(121)
source("com_poisson_pmf.R")
source("rcomp_exact.R")       # exact sampler
source("rcomp_rou.R")   # RoU sampler
source("rcomp_rejection_benson.R") # rejection sampler

# Dummy sampler for test
rcomp_dummy <- function(n, lambda, nu) {
  rpois(n, lambda)  # just Poisson as placeholder
}

# --- Choose methods to use ---
methods_to_use <- c("rou", "rej_ben", "exa")  
# options: "rou", "rej_ben", "exa", "dummy"

# --- Parameters ---
#lambda <- runif(1, 0.5, 10)
#nu     <- runif(1, 0.5, 2)
lambda <- 2
nu <- 0.5
log_lambda <- log(lambda)
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
  if ("rou" %in% methods_to_use) samples_list$rou <- rcomp_rou(n_sample, lambda, nu)
  if ("rej_ben" %in% methods_to_use) samples_list$rej_ben <- rcomp_rejection_benson(n_sample, lambda, nu)
  if ("exa" %in% methods_to_use) samples_list$exa <- rcomp_exact(n_sample, lambda, nu)
  if ("dummy" %in% methods_to_use) samples_list$dummy <- rcomp_dummy(n_sample, lambda, nu)
  
  # --- True pmf ---
  K <- max(unlist(samples_list)) + 10
  support <- 0:K
  log_p <- sapply(support, function(y) com_poisson_log_lpmf(y, log_lambda, nu, 1e-12))
  p_true <- exp(log_p) / sum(exp(log_p))
  df_true <- data.frame(x = support, y = p_true, Method = "True PMF")
  
  # --- Empirical histograms ---
  make_freq_df <- function(samples, name) {
    tbl <- as.data.frame(table(x = samples))
    data.frame(
      x = as.numeric(as.character(tbl$x)),
      y = tbl$Freq / sum(tbl$Freq),
      Method = name
    )
  }
  
  df_plot <- list(df_true)
  if ("rou" %in% methods_to_use) df_plot <- c(df_plot, list(make_freq_df(samples_list$rou, "Rejection RoU")))
  if ("rej_ben" %in% methods_to_use) df_plot <- c(df_plot, list(make_freq_df(samples_list$rej_ben, "Rejection Benson")))
  if ("exa" %in% methods_to_use) df_plot <- c(df_plot, list(make_freq_df(samples_list$exa, "Exact")))
  if ("dummy" %in% methods_to_use) df_plot <- c(df_plot, list(make_freq_df(samples_list$dummy, "Poisson")))
  df_plot <- do.call(rbind, df_plot)
  
  # --- PMF plot ---
  p1 <- ggplot(df_plot, aes(x = x, y = y, color = Method, linetype = Method)) +
    geom_point(size = 2) +
    geom_line(size = 0.5) +
    labs(title = paste0("λ=", lambda, ", ν=", nu, " (PMF)"),
         x = "y", y = "Probability / Rel. Frequency") +
    theme_minimal()
  
  # --- Deviation plot ---
  df_diff <- df_plot %>%
    left_join(df_true %>% select(x, y_true = y), by = "x") %>%
    mutate(diff = y - y_true) %>%
    filter(Method != "True PMF")
  
  p2 <- ggplot(df_diff, aes(x = x, y = diff, color = Method, linetype = Method)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_point(size = 2) +
    geom_line() +
    labs(title = paste0("λ=", lambda, ", ν=", nu, " (Deviation)"),
         x = "y", y = "Empirical - True") +
    theme_minimal()
  
  p1 + p2
}

# --- Build grid of plots ---
plots <- list(
  make_plot_pair(0.5, c("rou", "rej_ben", "exa")),
  make_plot_pair(1,   c("rou", "rej_ben", "exa", "dummy")),
  make_plot_pair(2,   c("rou", "rej_ben", "exa"))
)

final_plot <- (plots[[1]] / plots[[2]] / plots[[3]]) +
  plot_layout(guides = "collect") &  # collect all legends into one
  scale_color_manual(
    values = method_colors,
    breaks = c("True PMF", "Rejection RoU", "Rejection Benson", "Exact", "Poisson")
  ) &
  scale_linetype_manual(
    values = method_linetypes,
    breaks = c("True PMF", "Rejection RoU", "Rejection Benson", "Exact", "Poisson")
  ) &
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    legend.title = element_blank(),
    legend.justification = "center"
  )

final_plot
