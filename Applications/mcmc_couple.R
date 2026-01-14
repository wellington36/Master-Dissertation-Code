library(COMPoissonReg)
library(pbapply)
library(brms)
library(posterior)

# Thread control (important for benchmarking)
Sys.setenv(
  OMP_NUM_THREADS = 1,
  MKL_NUM_THREADS = 1,
  OPENBLAS_NUM_THREADS = 1,
  VECLIB_MAXIMUM_THREADS = 1,
  NUMEXPR_NUM_THREADS = 1
)

# Custom COM-Poisson family
com_poisson <- custom_family(
  name  = "com_poisson",
  dpars = c("mu", "nu"),
  links = c("log", "log"),
  lb    = c(0, 0),
  type  = "int"
)

# Stan functions
stan_code <- paste(readLines("comp_pmf.stan"), collapse = "\n")

stan_lik <- stanvar(
  scode  = stan_code,
  block = "functions"
)

leps <- -32 * log(2)

stan_leps <- stanvar(
  scode = sprintf(
    "real leps_custom() { return %f; }",
    leps
  ),
  block = "functions"
)

all_stanvars <- c(stan_lik, stan_leps)

# Data
data(couple)

# MCMC configuration
N_simulations <- 1
stan_chains   <- 1
stan_iter     <- 1000
stan_warmup   <- 800
core_number   <- 1

# Compile ONCE (dummy data)
base_fit <- brm(
  UPB ~ EDUCATION + ANXIETY,
  data = couple[1:5, ],   # tiny dataset → compile only
  chains = 0,
  family = com_poisson, # if running with Poisson comment
  prior =
    prior("gamma(0.1,0.01)", class = "Intercept", lb = 0) +
    prior("gamma(1, 1)", class = "nu", lb = 0),
  stanvars = all_stanvars,
  #family = "poisson",
  backend  = "cmdstanr"
)

# Single benchmark run
one_mcmc <- function(i) {
  
  t0 <- Sys.time()
  
  fit <- update(
    base_fit,
    newdata   = couple,
    chains    = stan_chains,
    iter      = stan_iter,
    warmup   = stan_warmup,
    recompile = FALSE
  )
  
  t1 <- Sys.time()
  
  # --- Extract draws ---
  draws <- as_draws_df(fit)
  mu_draws <- draws$b_Intercept
  nu_draws <- draws$nu  # if running with Poisson comment
  
  # --- Diagnostics ---
  coef_tab <- summary(fit)$fixed
  
  elapsed_sec <- as.numeric(difftime(t1, t0, units = "secs"))
  
  list(
    mu_mean     = mean(mu_draws),
    nu_mean     = mean(nu_draws), # if running with Poisson comment
    rhat        = coef_tab["Intercept", "Rhat"],
    ess         = coef_tab["Intercept", "Bulk_ESS"],
    time_sec    = elapsed_sec,
    ess_per_sec = coef_tab["Intercept", "Bulk_ESS"] / elapsed_sec
  )
}

# Parallel runs
results_list <- pblapply(
  1:N_simulations,
  FUN = function(i) one_mcmc(i),
  cl  = core_number
)

results <- do.call(rbind, lapply(results_list, as.data.frame))

# Summary
cat("\n=== Overall diagnostics ===\n")
cat("Mean ESS:", mean(results$ess), "\n")
cat("Mean Rhat:", mean(results$rhat), "\n")
cat("Mean time (s):", mean(results$time_sec), "\n")
cat("Mean ESS/sec:", mean(results$ess_per_sec), "\n")
