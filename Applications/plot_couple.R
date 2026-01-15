library(COMPoissonReg)
source("rcomp.R")

# Data
data(couple)

tab <- table(couple$UPB)

x_obs <- as.numeric(names(tab))
y_obs <- as.numeric(tab)
N     <- sum(y_obs)
K     <- length(x_obs)

# Estimated parameters (30k iterations 20k warmup eps = 2^-32)
mu_hat     <- 0.001178342
nu_hat     <- 0.1886298
lambda_hat <- mu_hat^nu_hat

lambda_poi_hat <- exp(0.824999)

# Monte Carlo settings
B <- 10000   # number of simulated datasets

# Storage
freq_cmp <- matrix(0, nrow = B, ncol = K)
freq_poi <- matrix(0, nrow = B, ncol = K)

# Monte Carlo simulation
set.seed(123)

for (b in 1:B) {
  
  # COM-Poisson
  sim_cmp <- rcomp(
    n      = N,
    lambda = lambda_hat,
    nu     = nu_hat
  )
  
  tab_cmp <- table(factor(sim_cmp, levels = x_obs))
  freq_cmp[b, ] <- as.numeric(tab_cmp)
  
  # Poisson
  sim_poi <- rpois(
    n      = N,
    lambda = lambda_poi_hat
  )
  
  tab_poi <- table(factor(sim_poi, levels = x_obs))
  freq_poi[b, ] <- as.numeric(tab_poi)
}

# Expected fitted frequencies (Monte Carlo mean)
y_cmp <- colMeans(freq_cmp)
y_poi <- colMeans(freq_poi)

# Errors
err_cmp <- y_obs - y_cmp
err_poi <- y_obs - y_poi

# Plot
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))

# Frequencies
plot(
  x_obs,
  y_obs,
  type = "h",
  lwd  = 3,
  col  = "black",
  xlab = "Unwanted Pursuit Behaviour (UPB)",
  ylab = "Frequency",
  main = "Observed vs Fitted (Simulation)"
)

lines(x_obs, y_cmp, lwd = 3, col = "steelblue")
lines(x_obs, y_poi, lwd = 3, col = "firebrick")

legend(
  "topright",
  legend = c("Observed", "COM-Poisson", "Poisson"),
  col    = c("black", "steelblue", "firebrick"),
  lwd    = 3,
  bty    = "n"
)

# Errors
plot(
  x_obs,
  err_poi,
  type = "l",
  lwd  = 3,
  col  = "firebrick",
  xlab = "UPB",
  ylab = "Observed - Fitted",
  main = "Fitting Error",
  ylim = range(c(err_cmp, err_poi))
)

abline(h = 0, lty = 2, col = "gray50")
lines(x_obs, err_cmp, lwd = 3, col = "steelblue")

legend(
  "topright",
  legend = c("COM-Poisson Error", "Poisson Error"),
  col    = c("steelblue", "firebrick"),
  lwd    = 3,
  bty    = "n"
)

par(mfrow = c(1, 1))
