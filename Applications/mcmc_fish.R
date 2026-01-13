library(pbapply)
library(brms)

source("rcomp.R")

zinb <- read.csv("https://paul-buerkner.github.io/data/fish.csv")


# MCMC run with brms
leps <- - 16 * log(2)
intercept_gold <- 0.000228946
N_simulations <- 1   # increase later
stan_chains <- 1
stan_iter <- 1000     # reduce while testing
stan_warmup <- 800
core_number <- 3      # 3 for my machine 34 for virtual


check_convergency <- function(rhat_intercept, estimate_intercept, intercept_gold, 
                              rel_tol = 0.01, abs_tol = 0.01) {
  # Calculate the absolute error
  abs_error <- abs(estimate_intercept - intercept_gold)
  
  # Dynamic threshold: handles both large scales (relative) and near-zero (absolute)
  threshold <- max(rel_tol * abs(intercept_gold), abs_tol)
  
  # Convergence criteria
  if (rhat_intercept < 1.01 && abs_error < threshold) {
    return(1)
  } else if (rhat_intercept < 1.05 && abs_error < (threshold * 0.5)) {
    # Slightly higher R-hat requires better accuracy
    return(1)
  } else if (rhat_intercept < 1.10 && abs_error < (threshold * 0.1)) {
    # Even higher R-hat requires much better accuracy
    return(1)
  } else {
    return(0)
  }
}


one_test <- function(i) {
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
    ess_per_sec = ess_bulk_intercept / elapsed_sec,
    converged   = check_convergency(rhat_intercept, 
                                    estimate_intercept, 
                                    intercept_gold)
  )
}

# Precompile once !!! NO WORK FOR VARIABLE DATASET !!!
scode_string <- sprintf("real leps_custom() { return %f; }", leps)
custom_stanvars <- stanvar(scode = scode_string, block = "functions")

base_fit <- brm(
  count ~ 1,
  data = zinb[1:5, ],   # tiny dummy dataset
  chains = 0,           # just compile, no sampling
  prior = prior("gamma(0.1,0.01)", class = "Intercept", lb = 0) +
    prior("gamma(1, 1)", class = "shape", lb = 0),
  stanvars = custom_stanvars,
  backend = "cmdstanr",
  family = "com_poisson"
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


