source("rcomp_benson.R")

# --- Conditional Moment Approximations ---
conditional_mean_default <- function(mu, nu) {
  #' Approximation for E[Y | mu, nu] for the CMP distribution.
  return(mu - (nu - 1) / (2 * nu))
}

conditional_variance_default <- function(mu, nu) {
  #' Approximation for Var[Y | mu, nu] for the CMP distribution.
  return(mu / nu)
}

# Gaunt-1term approximation
conditional_mean_gaunt_1term <- function(mu, nu) {
  #' Approximation for E[Y | mu, nu] for the CMP distribution.
  
  term1 <- mu - (nu - 1) / (2 * nu)
  term2 <- (nu^2 - 1) / (24 * nu^2 * mu)
  
  return(term1 - term2)
}

conditional_variance_gaunt_1term <- function(mu, nu) {
  #' Approximation for Var[Y | mu, nu] for the CMP distribution.
  
  factor <- mu / nu
  term_in_paren <- 1 + (nu^2 - 1) / (24 * nu^2 * mu^2)
  
  return(factor * term_in_paren)
}

# Gaunt-2terms approximation
conditional_mean_gaunt_2terms <- function(mu, nu) {
  #' Approximation for E[Y | mu, nu] for the CMP distribution.

  term1 <- mu - (nu - 1) / (2 * nu)
  term2 <- (nu^2 - 1) / (24 * nu^2 * mu)
  term3 <- (nu^2 - 1) / (24 * nu^3 * mu^2)
  
  return(term1 - term2 - term3)
}

conditional_variance_gaunt_2terms <- function(mu, nu) {
  #' Approximation for Var[Y | mu, nu] for the CMP distribution.

  factor <- mu / nu
  term_in_paren <- 1 + (nu^2 - 1) / (24 * nu^2 * mu^2) + (nu^2 - 1) / (12 * nu^3 * mu^3)
  
  return(factor * term_in_paren)
}

# Threshold approximation ---
threshold_sum <- function(f, mu, nu, eps) {
  k <- 1
  a_k <- f(mu, nu, k)
  logsum <- a_k
  
  terms <- numeric(10)
  while(a_k > eps || k < mu+1) {
    for(i in 1:10) {
      k=k+1
      terms[i] <- f(mu, nu, k)
    }
    
    a_k <- terms[10]
    logsum <- logsumexp(c(logsum, terms))
    terms <- numeric(10)
  }
  
  return(exp(logsum))
}

log_z_cmp <- function(mu, nu, i) {
  return(nu * ((i-1) * log(mu) - lgamma(i)))
}

log_z_moment_1 <- function(mu, nu, i) {
  return(log(i-1) + nu * ((i-1) * log(mu) - lgamma(i)))
}

log_z_moment_2 <- function(mu, nu, i) {
  return(2*log(i-1) + nu * ((i-1) * log(mu) - lgamma(i)))
}

conditional_mean_threshold <- function(mu, nu) {
  eps <- 1e-6
  
  Z_approx = threshold_sum(log_z_cmp, mu, nu, eps)
  First_moment_approx = threshold_sum(log_z_moment_1, mu, nu, eps)
  
  return(First_moment_approx/Z_approx)
}

conditional_variance_threshold <- function(mu, nu) {
  eps <- 1e-6
  
  Z_approx = threshold_sum(log_z_cmp, mu, nu, eps)
  First_moment_approx = threshold_sum(log_z_moment_1, mu, nu, eps)
  Second_moment_approx = threshold_sum(log_z_moment_2, mu, nu, eps)
  
  return(Second_moment_approx/Z_approx - (First_moment_approx/Z_approx)^2)
}

# --- Simulation Parameters ---
# Grid of (mu, nu) parameters to test
mu_values <- c(0.5, 1, 5, 10)
nu_values <- c(0.5, 1, 2, 5) # nu=1 is Poisson, nu<1 is under-dispersed, nu>1 is over-dispersed

# Number of samples for each (mu, nu) combination
N_SAMPLES <- 1000000

# Data frame to store all results
results_df <- expand.grid(mu = mu_values, nu = nu_values)
results_list <- list()

# --- Simulation Function ---
run_simulation <- function(mu, nu, n_samples) {
  
  samples <- rcomp_benson(n_samples, mu, nu)
  
  emp_mean <- mean(samples)
  emp_variance <- var(samples)
  
  # Default
  approx_mean_def <- conditional_mean_default(mu, nu)
  approx_var_def <- conditional_variance_default(mu, nu)
  
  # Gaunt 1-term
  approx_mean_gaunt1t <- conditional_mean_gaunt_1term(mu, nu)
  approx_var_gaunt1t <- conditional_variance_gaunt_1term(mu, nu)
  
  # Gaunt 2-terms
  approx_mean_gaunt2t <- conditional_mean_gaunt_2terms(mu, nu)
  approx_var_gaunt2t <- conditional_variance_gaunt_2terms(mu, nu)
  
  ## Threshold
  #approx_mean_threshold <- conditional_mean_threshold(mu, nu)
  #approx_var_threshold <- conditional_variance_threshold(mu, nu)
  
  data.frame(
    mu = mu,
    nu = nu,
    N = n_samples,
    
    # Empirical
    Empirical_Mean = emp_mean,
    Empirical_Variance = emp_variance,
    
    # Default Approximation
    Approx_Mean_Default = approx_mean_def,
    Approx_Var_Default = approx_var_def,
    
    # Gaunt-1term Approximation
    Approx_Mean_Gaunt_1term = approx_mean_gaunt1t,
    Approx_Var_Gaunt_1term = approx_var_gaunt1t,
    
    # Gaunt-2term Approximation
    Approx_Mean_Gaunt_2terms = approx_mean_gaunt2t,
    Approx_Var_Gaunt_2terms = approx_var_gaunt2t
    
    
    ## Threshold Approximation
    #Approx_Mean_Threshold = approx_mean_threshold,
    #Approx_Var_Threshold = approx_var_threshold
  )
}

cat("Starting simulation for", nrow(results_df), "parameter combinations...\n")

for (i in 1:nrow(results_df)) {
  mu_i <- results_df$mu[i]
  nu_i <- results_df$nu[i]
  
  cat(sprintf("Running for mu = %.2f, nu = %.2f...\n", mu_i, nu_i))
  
  result <- run_simulation(mu_i, nu_i, N_SAMPLES)
  results_list[[i]] <- result
}

final_results <- bind_rows(results_list)

final_results <- final_results %>%
  mutate(
    # Mean Errors
    Error_Mean_Default = Approx_Mean_Default - Empirical_Mean,
    Error_Mean_Gaunt_1term = Approx_Mean_Gaunt_1term - Empirical_Mean,
    Error_Mean_Gaunt_2terms = Approx_Mean_Gaunt_2terms - Empirical_Mean,
    #Error_Mean_Threshold = Approx_Mean_Threshold - Empirical_Mean,
    Rel_Error_Mean_Default = Error_Mean_Default / Empirical_Mean,
    Rel_Error_Mean_Gaunt_1term = Error_Mean_Gaunt_1term / Empirical_Mean,
    Rel_Error_Mean_Gaunt_2terms = Error_Mean_Gaunt_2terms / Empirical_Mean,
    #Rel_Error_Mean_Threshold = Error_Mean_Threshold / Empirical_Mean,
    
    # Variance Errors
    Error_Var_Default = Approx_Var_Default - Empirical_Variance,
    Error_Var_Gaunt_1term = Approx_Var_Gaunt_1term - Empirical_Variance,
    Error_Var_Gaunt_2terms = Approx_Var_Gaunt_2terms - Empirical_Variance,
    #Error_Var_Threshold = Approx_Var_Threshold - Empirical_Variance,
    Rel_Error_Var_Default = Error_Var_Default / Empirical_Variance,
    Rel_Error_Var_Gaunt_1term = Error_Var_Gaunt_1term / Empirical_Variance,
    Rel_Error_Var_Gaunt_2terms = Error_Var_Gaunt_2terms / Empirical_Variance
    #Rel_Error_Var_Threshold = Error_Var_Threshold / Empirical_Variance
  )

# Print the final table of results and errors
print(final_results)

plot_data_all <- final_results %>%
  select(mu, nu,
         starts_with("Rel_Error_Mean"),
         starts_with("Rel_Error_Var")) %>%
  pivot_longer(
    cols = starts_with("Rel_Error"),
    names_to = "Quantity",
    values_to = "Relative_Error"
  ) %>%
  mutate(
    Metric = ifelse(grepl("Mean", Quantity), "Mean", "Variance"),
    Approximation = Quantity %>%
      gsub("Rel_Error_(Mean|Var)_", "", .)
  
  )

ggplot(
  plot_data_all,
  aes(x = factor(mu),
      y = Relative_Error,
      fill = Approximation)
) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  facet_grid(
    rows = vars(Metric),
    cols = vars(nu),
    labeller = label_both
  ) +
  scale_fill_brewer() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(
    title = "Relative Error of COM-Poisson Moment Approximations",
    x = expression(mu),
    y = "Relative Error ((Approx − Empirical) / Empirical)",
    fill = "Approximation"
  ) +
  theme_minimal() +
  theme(
    panel.spacing = unit(1, "lines"),
    strip.text = element_text(face = "bold")
  )


# Prepare data for plotting (Mean)
plot_data_mean <- final_results %>%
  select(mu, nu, starts_with("Rel_Error_Mean")) %>%
  pivot_longer(
    cols = starts_with("Rel_Error_Mean"),
    names_to = "Approximation",
    values_to = "Relative_Error"
  ) %>%
  mutate(
    Approximation = gsub("Rel_Error_Mean_", "", Approximation)
  )

# Plot Relative Error for Mean
plot_mean <- ggplot(plot_data_mean, aes(x = factor(mu), y = Relative_Error, fill = Approximation)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap(~ factor(nu), labeller = label_both) +
  labs(
    title = "Relative Error of Mean Approximations for COM-Poisson",
    x = expression(mu),
    y = "Relative Error ((Approx - Empirical) / Empirical)",
    fill = "Approximation Type"
  ) +
  theme_minimal() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red")

plot_mean

plot_data_var <- final_results %>%
  select(mu, nu, starts_with("Rel_Error_Var")) %>%
  pivot_longer(
    cols = starts_with("Rel_Error_Var"),
    names_to = "Approximation",
    values_to = "Relative_Error"
  ) %>%
  mutate(
    Approximation = gsub("Rel_Error_Var_", "", Approximation)
  )

# Plot Relative Error for Variance
plot_variance <- ggplot(plot_data_var, aes(x = factor(mu), y = Relative_Error, fill = Approximation)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap(~ factor(nu), labeller = label_both) +
  labs(
    title = "Relative Error of Variance Approximations for COM-Poisson",
    x = expression(mu),
    y = "Relative Error ((Approx - Empirical) / Empirical)",
    fill = "Approximation Type"
  ) +
  theme_minimal() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red")

plot_variance
