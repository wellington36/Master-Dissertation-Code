library(cmdstanr)
library(posterior)
library(MASS)

# data
data(epil)

y_obs <- epil$y
N <- length(y_obs)

stan_data <- list(
  N    = N,
  y    = y_obs,
  leps = - 32 * log(2)
)

iter_warmup <- 7000
iter_sampling <- 3000
chains <- 4


# compile and sample stan 
mod_cmp <- cmdstan_model("compoisson.stan")
mod_poi <- cmdstan_model("poisson.stan")
mod_neb <- cmdstan_model("nb.stan")


fit_cmp <- mod_cmp$sample(
  data = stan_data,
  chains = chains,
  iter_warmup = iter_warmup,
  iter_sampling = iter_sampling
)

fit_poi <- mod_poi$sample(
  data = stan_data,
  chains = chains,
  iter_warmup = iter_warmup,
  iter_sampling = iter_sampling
)

fit_neb <- mod_neb$sample(
  data = stan_data,
  chains = chains,
  iter_warmup = iter_warmup,
  iter_sampling = iter_sampling
)

draws_df_cmp <- as_draws_df(fit_cmp)
draws_df_poi <- as_draws_df(fit_poi)
draws_df_neg <- as_draws_df(fit_neb)

mu_hat_cmp   <- mean(draws_df_cmp$mu)
nu_hat_cmp   <- mean(draws_df_cmp$nu)

lambda_hat_poi <-  mean(draws_df_poi$lambda)

mu_hat_neb <-  mean(draws_df_neg$mu)
phi_hat_neb <-  mean(draws_df_neg$phi)



# COM-Poisson pmf truncated (ONLY FOR THE PLOT)
dcom_trunc <- function(x, mu, nu, K) {
  logZ <- log(sum(exp(
    sapply(0:K, function(k)
      nu * (k * log(mu) - lgamma(k + 1)))
  )))
  
  exp(nu * (x * log(mu) - lgamma(x + 1)) - logZ)
}

x_vals <- 0:max(y_obs)
p_hat_cmp <- sapply(
  x_vals,
  dcom_trunc,
  mu = mu_hat_cmp,
  nu = nu_hat_cmp,
  K  = max(y_obs) * 20
)

p_hat_poi <- sapply(
  x_vals,
  dpois,
  lambda = lambda_hat_poi
)

p_hat_neg <- sapply(
  x_vals,
  dnbinom,
  size = phi_hat_neb,
  mu = mu_hat_neb
)


# observed empirical distribution
obs_tab <- table(y_obs)
obs_p   <- as.numeric(obs_tab) / sum(obs_tab)
obs_x   <- as.numeric(names(obs_tab))

# plot
par(mar = c(5, 5, 4, 2) + 0.1)
plot(
  NA, NA,
  xlim = range(x_vals),
  ylim = c(0, max(c(obs_p, p_hat_cmp, p_hat_poi, p_hat_neg)) * 1.1),
  xlab = "Count of seizures",
  ylab = "Probability",
  main = "Observed vs Posterior Predictive Distributions",
  cex.lab = 1.2,
  cex.main = 1.3
)

# grid
grid(nx = NA, ny = NULL, col = "gray85", lty = "dotted")

# observed distribution as bars
bar_width <- 0.4
rect(
  obs_x - bar_width / 2,
  0,
  obs_x + bar_width / 2,
  obs_p,
  col = gray(0.7),
  border = gray(0.4)
)

# cmp curve
lines(
  x_vals, p_hat_cmp,
  lwd = 2.5,
  col = "steelblue"
)

# poisson curve
lines(
  x_vals, p_hat_poi,
  lwd = 2.5,
  col = "firebrick"
)

# negative binomial curve
lines(
  x_vals, p_hat_neg,
  lwd = 2.5,
  col = "darkgreen"
)

legend(
  "topright",
  legend = c(
    "Observed (empirical)",
    "COM-Poisson (posterior mean)",
    "Poisson (posterior mean)",
    "Neg-Binomial (posterior mean)"
  ),
  pch = c(15, NA, NA, NA),
  pt.cex = c(1.6, NA, NA, NA),
  col = c(gray(0.7), "steelblue", "firebrick", "darkgreen"),
  lwd = c(NA, 2.5, 2.5, 2.5),
  lty = c(NA, 1, 1, 1),
  bty = "n",
  cex = 1.05
)
