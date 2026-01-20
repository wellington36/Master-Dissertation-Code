library(cmdstanr)
library(pbapply)
library(posterior)
library(tidyr)


# thread control (important for benchmarking)
Sys.setenv(
  OMP_NUM_THREADS = 1,
  MKL_NUM_THREADS = 1,
  OPENBLAS_NUM_THREADS = 1,
  VECLIB_MAXIMUM_THREADS = 1,
  NUMEXPR_NUM_THREADS = 1
)

# mcmc parameters
leps          <- - 1 * log(2)
N_simulations <- 100
chains        <- 1
iter_warmup   <- 1000
iter_sampling <- 500
core_number   <- 100      # 3 for my machine 34, 50 or 100 for virtual machine

# data
house_prices <- read.csv("https://wellington36.github.io/house_prices.csv")
house_bathroom <- drop_na(house_prices, "Bathroom")$Bathroom

y_obs <- house_bathroom
N <- length(y_obs)

stan_data <- list(
  N    = N,
  y    = y_obs,
  leps = leps
)


# compile stan
mod_cmp <- cmdstan_model("compoisson.stan")

one_mcmc <- function(i) {
  t0 <- Sys.time()
  
  fit_cmp <- mod_cmp$sample(
    data = stan_data,
    chains = chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling
  )
  
  t1 <- Sys.time()
  
  draws_df_cmp <- as_draws_df(fit_cmp)
  sum_cmp <- fit_cmp$summary()
  
  elapsed_sec <- as.numeric(difftime(t1, t0, units = "secs"))
  
  # keep only parameters
  pars <- c("mu", "nu")
  idx  <- sum_cmp$variable %in% pars
  
  ess_mean  <- mean(sum_cmp$ess_bulk[idx])
  rhat_max  <- max(sum_cmp$rhat[idx])
  
  data.frame(
    mu_hat_cmp  = mean(draws_df_cmp$mu),
    nu_hat_cmp  = mean(draws_df_cmp$nu),
    ess_mean    = ess_mean,
    rhat_max    = rhat_max,
    time_sec    = elapsed_sec,
    ess_per_sec = ess_mean / elapsed_sec
  )
}


# run in parallel
results_list <- pblapply(1:N_simulations,
                         FUN = function(i) {
                           try({
                             one_mcmc(i)
                           }, silent = TRUE)
                         },
                         cl = core_number)
results <- do.call(rbind, results_list)

cat("\n=== Overall diagnostics ===\n")
cat("Mean mu:", mean(results$mu_hat_cmp), "\n")
cat("Mean nu:", mean(results$nu_hat_cmp), "\n")
cat("Mean ESS:", mean(results$ess_mean), "\n")
cat("Mean Rhat:", mean(results$rhat_max), "\n")
cat("Mean time (s):", mean(results$time_sec), "\n")
cat("Mean ESS/sec:", mean(results$ess_mean)/mean(results$time_sec), "\n")
