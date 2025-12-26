source("rcomp.R")
source("rcomp_benson.R")

random_variables_cmp_mu_nu <- function(N, mu_alpha_prior, mu_beta_prior,
                                       nu_alpha_prior, nu_beta_prior) {
  out <- numeric(N)
  out_mu <- numeric(N)
  out_nu <- numeric(N)
  
  for (i in seq_len(N)) {
    mu <- rgamma(1, shape = mu_alpha_prior, rate = mu_beta_prior)
    nu <- rgamma(1, shape = nu_alpha_prior, rate = nu_beta_prior)
    lambda <- mu^nu
    
    print(sprintf("CMP Mean %.16f, Dispersion %.4f", mu, nu))
    out[i] <- rcomp_benson(1, lambda, nu)
    out_mu[i] <- mu
    out_nu[i] <- nu
  }
  
  data_mean_mu <- mean(out_mu)
  data_var_mu  <- var(out_mu)
  data_mean_nu <- mean(out_nu)
  data_var_nu  <- var(out_nu)
  
  cat("Empirical mean mu:", data_mean_mu, "\n")
  cat("Empirical variance mu:", data_var_mu, "\n")
  cat("Empirical mean nu:", data_mean_nu, "\n")
  cat("Empirical variance nu:", data_var_nu, "\n")
  
  return(out)
}

# -------------------------------
# Parameters
# -------------------------------
N <- 10^5

mu_alpha_prior <- 0.1
mu_beta_prior <- 0.01
nu_alpha_prior <- 1
nu_beta_prior <- 1


mu_alpha_prior <- 1
mu_beta_prior <- 1
nu_alpha_prior <- 0.0625
nu_beta_prior <- 0.25

data <- random_variables_cmp_mu_nu(N, mu_alpha_prior, mu_beta_prior,
                                   nu_alpha_prior, nu_beta_prior)

# -------------------------------
# Empirical statistics
# -------------------------------
data_mean <- mean(data)
data_var  <- var(data)

cat("Empirical mean:", data_mean, "\n")
cat("Empirical variance:", data_var, "\n")

# -------------------------------
# Plotting
# -------------------------------
par(mfrow = c(1, 3))

# 1. Histogram
hist(
  data,
  main = "Histogram of CMP Samples",
  xlab = "Value",
  col = "lightblue",
  border = "white"
)
mtext(sprintf("Mean = %.2f, Var = %.2f", data_mean, data_var), side = 3, line = 0.3)

# 2. ECDF
plot(
  ecdf(data),
  main = "ECDF of CMP Samples",
  xlab = "Value",
  ylab = "F(x)"
)
mtext(sprintf("Mean = %.2f", data_mean), side = 3, line = 0.3)

# 3. Scatter (index vs value)
plot(
  data,
  pch = 19,
  col = "blue",
  main = "Scatter Plot of CMP Samples",
  xlab = "Index",
  ylab = "Value"
)
mtext(sprintf("Var = %.2f", data_var), side = 3, line = 0.3)
