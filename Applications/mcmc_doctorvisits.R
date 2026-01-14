library(AER)
library(pbapply)
library(brms)

source("rcomp.R")

data("DoctorVisits")


# MCMC run with brms
leps <- - 32 * log(2)
N_simulations <- 1   # increase later
stan_chains <- 1
stan_iter <- 30000     # reduce while testing
stan_warmup <- 20000
core_number <- 4      # 4 for my machine 34 for virtual



one_test <- function(i) {
  t0 <- Sys.time()
  
  fit <- update(base_fit,
                newdata = DoctorVisits, 
                chains = stan_chains,
                iter = stan_iter,
                warmup = stan_warmup,
                recompile = FALSE)   # reuse compiled model
  
  t1 <- Sys.time()
  
  # --- Extract draws ---
  draws <- as_draws_df(fit)
  mu_draws <- draws$b_Intercept
  nu_draws <- draws$shape
  
  # --- Diagnostics with posterior package ---
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
scode_string <- sprintf("real leps_custom() { return %f; }", leps)
custom_stanvars <- stanvar(scode = scode_string, block = "functions")

base_fit <- brm(
  visits ~ age + income + illness,
  data = DoctorVisits[1:5, ],   # tiny dummy dataset
  chains = 0,           # just compile, no sampling
  #prior = prior("gamma(0.1,0.01)", class = "Intercept", lb = 0) +
  #  prior("gamma(1, 1)", class = "shape", lb = 0),
  #stanvars = custom_stanvars,
  backend = "cmdstanr",
  family = "poisson"
)

# --- Run in parallel ---
results_list <- pblapply(1:N_simulations,
                         FUN = function(i) {
                           try({
                             one_test(i)   # check for inconsistenses
                           }, silent = TRUE)
                         },
                         cl = core_number)
results <- do.call(rbind, lapply(results_list, as.data.frame))

cat("\n=== Overall diagnostics ===\n")
cat("Mean ESS:", mean(results$ess), "\n")
cat("Mean time (s):", mean(results$time_sec), "\n")
cat("Mean ESS/sec:", mean(results$ess_per_sec), "\n")
cat("MCMCs converged:", sum(results$converged), "of", N_simulations, "\n")
