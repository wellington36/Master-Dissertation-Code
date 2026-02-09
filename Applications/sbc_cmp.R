library(pbapply)
library(cmdstanr)
library(posterior)

source("rcomp.R")

# thread control (important for benchmarking)
Sys.setenv(
  OMP_NUM_THREADS = 1,
  MKL_NUM_THREADS = 1,
  OPENBLAS_NUM_THREADS = 1,
  VECLIB_MAXIMUM_THREADS = 1,
  NUMEXPR_NUM_THREADS = 1
)

# sbc parameters
leps_exp      <- -32
N_simulations <- 1000
J             <- 100
chains        <- 1
iter_warmup   <- 1800
iter_sampling <- 200
core_number   <- 100     # 3 for my machine 50 or 100 for virtual machine

leps <- leps_exp * log(2)

# compile stan
mod_cmp <- cmdstan_model("compoisson.stan")

# sbc prior parameters list
theta_list <- lapply(seq_len(N_simulations), function(i) {
  c(
    mu = rgamma(1, 0.1, 0.01),
    nu = rgamma(1, 1, 1)
  )
})

# order the dificult of the problem
nu_vals <- sapply(theta_list, function(x) x["nu"])
theta_list <- theta_list[order(nu_vals)]


# one sbc run
sbc_run <- function(i) {
  t0 <- Sys.time()
  
  # prior draw
  mu_sim <- unname(theta_list[[i]]["mu"])
  nu_sim <- unname(theta_list[[i]]["nu"])
  
  # simulate data
  data <- rcomp(J, mu_sim, nu_sim)
  
  # stan fit
  stan_data <- list(
    N    = J,
    y    = data,
    leps = leps
  )
  
  fit_cmp <- mod_cmp$sample(
    data = stan_data,
    chains = chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    adapt_delta = 0.8,
    max_treedepth = 8
  )
  
  # extract draws and diagnostics
  t1 <- Sys.time()
  
  draws_df_cmp <- as_draws_df(fit_cmp)
  sum_cmp <- fit_cmp$summary()
  
  mu_mean <- mean(draws_df_cmp$mu)
  nu_mean <- mean(draws_df_cmp$nu)
  
  elapsed_sec <- as.numeric(difftime(t1, t0, units = "secs"))
  
  pars <- c("mu", "nu")
  idx  <- sum_cmp$variable %in% pars
  
  ess_mean  <- mean(sum_cmp$ess_bulk[idx])
  rhat_max  <- max(sum_cmp$rhat[idx])
  
  list(
    rank_mu     = sum(draws_df_cmp$mu < mu_sim),
    rank_nu     = sum(draws_df_cmp$nu < nu_sim),
    mu_mean     = mu_mean,
    nu_mean     = nu_mean,
    mu_sim      = mu_sim,
    nu_sim      = nu_sim,
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
                             sbc_run(i)
                           }, silent = TRUE)
                         },
                         cl = core_number)
results <- do.call(rbind, results_list)

# save
file_name = sprintf("sbc_N%i_leps%.0f.csv", N_simulations, leps_exp)
write.csv(results, file_name, row.names = FALSE)
