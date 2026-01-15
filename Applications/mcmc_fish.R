library(pbapply)
library(brms)

# Thread control (important for benchmarking)
Sys.setenv(
  OMP_NUM_THREADS = 1,
  MKL_NUM_THREADS = 1,
  OPENBLAS_NUM_THREADS = 1,
  VECLIB_MAXIMUM_THREADS = 1,
  NUMEXPR_NUM_THREADS = 1
)

zinb <- read.csv("https://paul-buerkner.github.io/data/fish.csv")


# MCMC run with brms
leps          <- - 32 * log(2)
N_simulations <- 1
stan_chains   <- 1
stan_iter     <- 30000     # reduce while testing
stan_warmup   <- 20000
core_number   <- 3      # 3 for my machine 50 for virtual


one_mcmc <- function(i) {
  t0 <- Sys.time()
  
  fit <- update(base_fit,
                newdata = zinb, 
                chains = stan_chains,
                iter = stan_iter,
                warmup = stan_warmup,
                recompile = FALSE)   # reuse compiled model
  
  t1 <- Sys.time()
  
  # --- Extract draws ---
  draws <- as_draws_df(fit)
  mu_draws <- draws$b_Intercept
  nu_draws <- draws$shape
  
  # --- Diagnostics ---
  coef_tab <- summary(fit)$fixed
  rhat_intercept <- coef_tab["Intercept", "Rhat"]
  estimate_intercept <- coef_tab["Intercept", "Estimate"]
  ess_bulk_intercept <- coef_tab["Intercept", "Bulk_ESS"]
  
  elapsed_sec <- as.numeric(difftime(t1, t0, units = "secs"))
  
  list(
    mu_mean     = mean(mu_draws),
    nu_mean     = mean(nu_draws),
    rhat        = rhat_intercept,
    ess         = ess_bulk_intercept,
    time_sec    = elapsed_sec,
    ess_per_sec = ess_bulk_intercept / elapsed_sec
  )
}

# Precompile once !!! NO WORK FOR VARIABLE DATASET !!!
stan_leps <- stanvar(
  scode = sprintf(
    "real leps_custom() { return %f; }",
    leps
  ),
  block = "functions"
)

base_fit <- brm(
  count ~ 1,
  data = zinb[1:5, ],   # tiny dummy dataset
  chains = 0,           # just compile, no sampling
  prior = prior("gamma(0.1,0.01)", class = "Intercept", lb = 0) +
    prior("gamma(1, 1)", class = "shape", lb = 0),
  stanvars = stan_leps,
  backend = "cmdstanr",
  family = "com_poisson"
)

# --- Run in parallel ---
results_list <- pblapply(1:N_simulations,
                         FUN = function(i) {
                           try({
                             one_mcmc(i)   # check for inconsistenses
                           }, silent = TRUE)
                         },
                         cl = core_number)
results <- do.call(rbind, lapply(results_list, as.data.frame))

cat("\n=== Overall diagnostics ===\n")
cat("Mean ESS:", mean(results$ess), "\n")
cat("Mean Rhat:", mean(results$rhat), "\n")
cat("Mean time (s):", mean(results$time_sec), "\n")
cat("Mean ESS/sec:", mean(results$ess_per_sec), "\n")


