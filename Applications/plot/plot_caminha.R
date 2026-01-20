library(cmdstanr)
library(posterior)

# data
caminha <- read.csv("https://wellington36.github.io/word_count_caminha.csv")

y_obs <- caminha$count
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

# posterior means
draws_df_cmp <- as_draws_df(fit_cmp)
draws_df_poi <- as_draws_df(fit_poi)
draws_df_neg <- as_draws_df(fit_neb)

mu_hat_cmp   <- mean(draws_df_cmp$mu)
nu_hat_cmp   <- mean(draws_df_cmp$nu)

lambda_hat_poi <- mean(draws_df_poi$lambda)

mu_hat_neb  <- mean(draws_df_neg$mu)
phi_hat_neb <- mean(draws_df_neg$phi)

# COM-Poisson pmf truncated (ONLY FOR THE PLOT)
dcom_trunc <- function(x, mu, nu, K) {
  logZ <- log(sum(exp(
    sapply(0:K, function(k)
      nu * (k * log(mu) - lgamma(k + 1)))
  )))
  
  exp(nu * (x * log(mu) - lgamma(x + 1)) - logZ)
}

# X support with 40+ grouping
x_max_plot <- 40
x_max_full <- max(y_obs) * 20

x_vals_full <- 0:x_max_full
x_vals_plot <- 0:x_max_plot

# predictive distributions (FULL support first)
p_full_cmp <- sapply(
  x_vals_full,
  dcom_trunc,
  mu = mu_hat_cmp,
  nu = nu_hat_cmp,
  K  = x_max_full
)

p_full_poi <- dpois(
  x_vals_full,
  lambda = lambda_hat_poi
)

p_full_neg <- dnbinom(
  x_vals_full,
  size = phi_hat_neb,
  mu   = mu_hat_neb
)

# GROUP tail at 40+
p_hat_cmp <- c(
  p_full_cmp[1:x_max_plot],
  sum(p_full_cmp[(x_max_plot + 1):length(p_full_cmp)])
)

p_hat_poi <- c(
  p_full_poi[1:x_max_plot],
  sum(p_full_poi[(x_max_plot + 1):length(p_full_poi)])
)

p_hat_neg <- c(
  p_full_neg[1:x_max_plot],
  sum(p_full_neg[(x_max_plot + 1):length(p_full_neg)])
)

# observed empirical distribution
obs_tab <- table(y_obs)

obs_le_40 <- obs_tab[as.numeric(names(obs_tab)) < x_max_plot]
obs_ge_40 <- sum(obs_tab[as.numeric(names(obs_tab)) >= x_max_plot])

obs_tab_grp <- c(obs_le_40, "40+" = obs_ge_40)

obs_p <- as.numeric(obs_tab_grp) / sum(obs_tab_grp)

obs_x <- c(
  as.numeric(names(obs_le_40)),
  x_max_plot
)

# plot
par(mar = c(5, 5, 4, 2) + 0.1)

plot(
  NA, NA,
  xlim = c(0, x_max_plot),
  ylim = c(0, max(c(obs_p, p_hat_cmp, p_hat_poi, p_hat_neg)) * 1.1),
  xlab = "Count of word frequency",
  ylab = "Probability",
  main = "Observed vs Posterior Predictive Distributions",
  cex.lab = 1.2,
  cex.main = 1.3,
  xaxt = "n"
)
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
  x_vals_plot, p_hat_cmp,
  lwd = 2.5,
  col = "steelblue"
)

# poisson curve
lines(
  x_vals_plot, p_hat_poi,
  lwd = 2.5,
  col = "firebrick"
)

# negative binomial
lines(
  x_vals_plot, p_hat_neg,
  lwd = 2.5,
  col = "darkgreen"
)

axis(1, at = x_vals_plot, labels = c(0:(x_max_plot - 1), "40+"))

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
