source("rcomp_benson.R")

random_variables_cmp_mu_nu <- function(N, mu_alpha_prior, mu_beta_prior,
                                       nu_alpha_prior, nu_beta_prior) {
  out <- numeric(N)
  
  for (i in seq_len(N)) {
    mu <- rgamma(1, shape = mu_alpha_prior, rate = mu_beta_prior)
    nu <- rgamma(1, shape = nu_alpha_prior, rate = nu_beta_prior)
    lambda <- mu^nu
    
    print(sprintf("CMP Mean %.16f, Dispersion %.4f", mu, nu))
    out[i] <- rcomp_benson(1, lambda, nu)
  }
  
  return(out)
}

# -------------------------------
# Parameters
# -------------------------------
N <- 10000

mu_alpha_prior <- 1
mu_beta_prior <- 0.01
nu_alpha_prior <- 1
nu_beta_prior <- 1

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
mtext(sprintf("Mean = %.2e, Var = %.2e", data_mean, data_var), side = 3, line = 0.3)

# 2. ECDF
plot(
  ecdf(data),
  main = "ECDF of CMP Samples",
  xlab = "Value",
  ylab = "F(x)"
)
mtext(sprintf("Mean = %.2e", data_mean), side = 3, line = 0.3)

# 3. Scatter (index vs value)
plot(
  data,
  pch = 19,
  col = "blue",
  main = "Scatter Plot of CMP Samples",
  xlab = "Index",
  ylab = "Value"
)
mtext(sprintf("Var = %.2e", data_var), side = 3, line = 0.3)
