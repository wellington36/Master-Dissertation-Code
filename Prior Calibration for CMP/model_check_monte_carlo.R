# --- Conditional Moment Approximations ---

conditional_mean <- function(mu, nu) {
  #' Approximation for E[Y | mu, nu] for the CMP distribution.

  return(mu + (nu - 1) / (2 * nu))
}

conditional_variance <- function(mu, nu) {
  #' Approximation for Var[Y | mu, nu] for the CMP distribution.
  
  return(mu / nu)
}

# --- Monte Carlo Moment Estimation ---

estimate_total_moments <- function(hyperparameters, N_samples = 10000) {
  #' Estimates E[Y] and Var[Y] using Monte Carlo simulation and the Law of Total
  #' Expectation and Variance with conditional moment approximations.
  #'
  #' @param hyperparameters A vector: [alpha_mu, beta_mu, alpha_nu, beta_nu]
  #' @param N_samples Number of Monte Carlo samples.
  #' @return A list with estimated E_Y_hat and Var_Y_hat.
  
  alpha_mu <- hyperparameters[1]
  beta_mu <- hyperparameters[2]
  alpha_nu <- hyperparameters[3]
  beta_nu <- hyperparameters[4]

  
  # Sample mu and nu from their Gamma priors
  mu_samples <- rgamma(N_samples, shape = alpha_mu, rate = beta_mu)
  nu_samples <- rgamma(N_samples, shape = alpha_nu, rate = beta_nu)
  
  # Calculate conditional moments for each sample
  M_samples <- conditional_mean(mu_samples, nu_samples)
  V_samples <- conditional_variance(mu_samples, nu_samples)
  
  # Estimate Total Expectation: E[Y] = E[E[Y | mu, nu]]
  E_Y_hat <- mean(M_samples)
  
  # Estimate Total Variance: Var[Y] = Var[E[Y | mu, nu]] + E[Var[Y | mu, nu]]
  Var_E_Y_hat <- var(M_samples)
  E_Var_Y_hat <- mean(V_samples)
  
  Var_Y_hat <- Var_E_Y_hat + E_Var_Y_hat
  
  return(list(E_Y_hat = E_Y_hat, Var_Y_hat = Var_Y_hat))
}


hyperparameters <- c(0.01, 0.1, 1, 1)
print(estimate_total_moments(hyperparameters))