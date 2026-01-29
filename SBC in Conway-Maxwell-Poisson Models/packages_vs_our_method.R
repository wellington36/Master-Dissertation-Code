source("rcomp_benson.R")
source("rcomp.R")

library(COMPoissonReg)
library(ggplot2)
library(dplyr)

method_colors <- c(
  "compoisson"    = "#1b9e77",
  "CompGLM"       = "#d95f02",
  "COMPoissonReg" = "#e7298a",
  "rcomp_benson"  = "#7570b3",
  "rcomp (exact)" = "lightskyblue",
  "Brute force"   = "black"
)

set.seed(123)

# Functions from discontinued packages)
com_log_sum <- function(x, y) {
  if (x == -Inf) return(y) else if (y == -Inf) return(x)
  if (x > y) return(x + log(1 + exp(y - x))) else return(y + log(1 + exp(x - y)))
}

# Implementation from 'compoisson' package
rcom_legacy <- function(n, lambda, nu) {
  log_z <- -Inf; z_last <- 0; j <- 0
  while (abs(log_z - z_last) > 0.001) {
    z_last <- log_z
    log_z <- com_log_sum(log_z, j * log(lambda) - nu * lfactorial(j))
    j <- j + 1
  }
  r <- numeric(n)
  for (i in 1:n) {
    log_u <- log(runif(1))
    j <- 0; curr_log_cum <- -Inf
    while (TRUE) {
      log_pj <- (j * log(lambda) - nu * lfactorial(j)) - log_z
      if (log_u < log(exp(curr_log_cum) + exp(log_pj))) break
      curr_log_cum <- com_log_sum(curr_log_cum, log_pj)
      j <- j + 1
      if (j > 1000) break
    }
    r[i] <- j
  }
  return(r)
}

# Implementation from 'CompGLM' package
rcomp_legacy <- function(n, lambda, nu, sumTo = 100) {
  j <- 0:sumTo
  log_probs <- j * log(lambda) - nu * lfactorial(j)
  probs <- exp(log_probs - max(log_probs))
  probs <- probs / sum(probs)
  return(sample(0:sumTo, size = n, replace = TRUE, prob = probs))
}

# Brute force with high number of terms
brute_force <- function(lambda, nu, sumTo = 100000) {
  j <- 0:sumTo
  log_z <- j * log(lambda) - nu * lfactorial(j)
  p <- exp(log_z - max(log_z))
  p <- p / sum(p)
  m <- sum(j * p)
  v <- sum((j^2) * p) - (m^2)
  return(list(mean = m, var = v))
}

brute_force_pmf <- function(lambda, nu, sumTo = 100000) {
  j <- 0:sumTo
  log_p <- j * log(lambda) - nu * lfactorial(j)
  p <- exp(log_p - max(log_p))
  p <- p / sum(p)
  data.frame(x = j, p = p)
}

empirical_pmf <- function(samples, max_x = NULL) {
  if (is.null(max_x)) max_x <- max(samples)
  tab <- table(factor(samples, levels = 0:max_x))
  data.frame(
    x = 0:max_x,
    p = as.numeric(tab) / sum(tab)
  )
}


# Evaluation
evaluate_benckmarks <- function(n = 5000, lambda, nu) {
  brute_force <- brute_force(lambda, nu, sumTo = n)
  
  # Measure time and generate samples
  cat("Running benchmarks...\n")
  
  t1 <- system.time(s1 <- rcom_legacy(n, lambda, nu))["elapsed"]
  t2 <- system.time(s2 <- rcomp_legacy(n, lambda, nu))["elapsed"]
  
  # COMPoissonReg check
  if(requireNamespace("COMPoissonReg", quietly = TRUE)) {
    t3 <- system.time(s3 <- COMPoissonReg::rcmp(n, lambda, nu))["elapsed"]
  } else {
    s3 <- rep(NA, n); t3 <- NA
  }
  
  # Benson function check
  if(exists("rcomp_benson")) {
    mu <- lambda^(1/nu)
    t4 <- system.time(s4 <- rcomp_benson(n, mu, nu))["elapsed"]
  } else {
    s4 <- rep(NA, n); t4 <- NA
  }
  
  # Our function check
  if(exists("rcomp")) {
    t5 <- system.time(s5 <- rcomp(n, lambda, nu))["elapsed"]
  } else {
    s5 <- rep(NA, n); t5 <- NA
  }
  
  # Create results table
  results <- data.frame(
    Method = c("Brute force", "compoisson", "CompGLM", "COMPoissonReg", "rcomp_benson", "rcomp (exact)"),
    Mean = c(brute_force$mean, mean(s1), mean(s2), mean(s3, na.rm=T), mean(s4, na.rm=T), mean(s5, na.rm=T)),
    Var  = c(brute_force$var,  var(s1),  var(s2),  var(s3, na.rm=T),  var(s4, na.rm=T), var(s5, na.rm=T)),
    Time_Sec = c(0, t1, t2, t3, t4, t5)
  )
  
  # Add error columns
  results$Mean_Err <- abs(results$Mean - brute_force$mean)
  results$Var_Err  <- abs(results$Var - brute_force$var)
  
  # Reorder columns for better reading
  results <- results[, c("Method", "Mean", "Mean_Err", "Var", "Var_Err", "Time_Sec")]
  
  print(results, digits = 4)
}


evaluate_plot_data <- function(n, lambda, nu, max_x_plot) {
  
  # Brute force PMF
  bf <- brute_force_pmf(lambda, nu, sumTo = max_x_plot) |>
    mutate(Method = "Brute force")
  
  # Samples
  s1 <- rcom_legacy(n, lambda, nu)
  s2 <- rcomp_legacy(n, lambda, nu)
  s3 <- COMPoissonReg::rcmp(n, lambda, nu)
  
  s4 <- {
    mu <- lambda^(1 / nu)
    rcomp_benson(n, mu, nu)
  }
  
  s5 <- rcomp(n, lambda, nu)
  
  pmfs <- bind_rows(
    empirical_pmf(s1, max_x_plot) |> mutate(Method = "compoisson"),
    empirical_pmf(s2, max_x_plot) |> mutate(Method = "CompGLM"),
    empirical_pmf(s3, max_x_plot) |> mutate(Method = "COMPoissonReg"),
    #empirical_pmf(s4, max_x_plot) |> mutate(Method = "rcomp_benson"),
    empirical_pmf(s5, max_x_plot) |> mutate(Method = "rcomp (exact)"),
    bf
  )
  
  pmfs |>
    mutate(
      scenario = sprintf("λ=%.2f, ν=%.2f, n=%i", lambda, nu, n)
    )
}

#plot_data <- bind_rows(
#  evaluate_plot_data(n = 10000, lambda = 2, nu = 0.2, max_x_plot = 80),
#  evaluate_plot_data(n = 500,   lambda = 2, nu = 1,   max_x_plot = 5),
#  evaluate_plot_data(n = 500,   lambda = 2, nu = 2,   max_x_plot = 5)
#)


plot_data <- bind_rows(
  evaluate_plot_data(n = 10000, lambda = 2, nu = 0.15, max_x_plot = 200),
  evaluate_plot_data(n = 10000, lambda = 2, nu = 0.1, max_x_plot = 2000)
)

plot_data <- plot_data |>
  mutate(
    scenario = factor(
      scenario,
      levels = c(
        "λ=2.00, ν=0.15, n=10000",
        "λ=2.00, ν=0.10, n=10000"
      )
    )
  )

# Create labels (a), (b), (c)
scenario_levels <- unique(plot_data$scenario)

panel_labels <- setNames(
  paste0("(", letters[seq_along(scenario_levels)], ") ", scenario_levels),
  scenario_levels
)


ggplot(plot_data, aes(x = x, y = p)) +
  
# Other methods (draw first)
geom_ribbon(
  data = subset(plot_data, Method != "Brute force"),
  aes(ymin = 0, ymax = p, fill = Method),
  alpha = 0.05,
  linewidth = 0
) +
  geom_line(
    data = subset(plot_data, Method != "Brute force"),
    aes(color = Method),
    linewidth = 0.8
  ) +
  
# Brute force (draw last)
geom_ribbon(
  data = subset(plot_data, Method == "Brute force"),
  aes(ymin = 0, ymax = p, fill = Method),
  alpha = 0.05,
  linewidth = 0
) +
  geom_line(
    data = subset(plot_data, Method == "Brute force"),
    aes(color = Method),
    linewidth = 0.8
  ) +
  
  facet_wrap(
    ~ scenario,
    scales = "free",
    ncol = 1,
    labeller = labeller(scenario = panel_labels)
  ) +
  
  scale_color_manual(values = method_colors) +
  scale_fill_manual(values = method_colors) +
  
  labs(
    x = "x",
    y = "Probability"
  ) +
  
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

evaluate_benckmarks(n = 1000, lambda = 0.1, nu = 2.3)