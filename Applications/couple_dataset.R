library(COMPoissonReg)
library(parallel)
library(pbapply)
library(brms)
data(couple)

tab <- table(couple$UPB)
intercept_gold <- -0.35

# Check data
plot(
  as.numeric(names(tab)),
  as.numeric(tab),
  type = "b",
  lwd = 2,
  col = "steelblue",
  xlab = "Unwanted Pursuit Behaviour (UPB)",
  ylab = "Frequency"
)

# MCMC run with brms
leps <- - 16 * log(2)
N_simulations <- 1   # increase later
stan_chains <- 1
stan_iter <- 10000     # reduce while testing
stan_warmup <- 8000

check_convergency <- function(rhat_intercept, estimate_intercept, intercept_gold) {
  if (rhat_intercept < 1.01 & 
      (estimate_intercept > intercept_gold - 0.01 & 
       estimate_intercept < intercept_gold + 0.01)) {
    return(1)
  } else if (rhat_intercept < 1.05 & 
             (estimate_intercept > intercept_gold - 0.005 & 
              estimate_intercept < intercept_gold + 0.005)) {
    return(1)
  } else if (rhat_intercept < 1.20 & 
             (estimate_intercept > intercept_gold - 0.0002 & 
              estimate_intercept < intercept_gold + 0.0002)) {
    return(1)
  } else {
    return(0)
  }
}

one_test <- function(i) {
  t0 <- Sys.time()
  
  fit <- update(base_fit,
                newdata = couple, 
                chains = stan_chains,
                iter = stan_iter,
                warmup = stan_warmup,
                recompile = FALSE)   # reuse compiled model
  
  t1 <- Sys.time()
  
  # --- Extract draws ---
  draws <- as_draws_df(fit)
  lambda_draws <- exp(draws$b_Intercept)
  nu_draws <- draws$shape
  
  # --- Diagnostics with posterior package ---
  coef_tab <- summary(fit)$fixed
  rhat_intercept <- coef_tab["Intercept", "Rhat"]
  estimate_intercept <- coef_tab["Intercept", "Estimate"]
  ess_bulk_intercept <- coef_tab["Intercept", "Bulk_ESS"]
  
  elapsed_sec <- as.numeric(difftime(t1, t0, units = "secs"))
  
  list(
    lambda_mean = log(mean(lambda_draws)),
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
  UPB ~ 1,
  data = couple[1:5, ],   # tiny dummy dataset
  chains = 0,           # just compile, no sampling
  prior = prior(normal(0, 1), class = "Intercept") +
    prior(lognormal(0, 1), class = "shape"),
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
                         cl = 3)
results <- do.call(rbind, lapply(results_list, as.data.frame))

cat("\n=== Overall diagnostics ===\n")
cat("Mean ESS:", mean(results$ess), "\n")
cat("Mean time (s):", mean(results$time_sec), "\n")
cat("Mean ESS/sec:", mean(results$ess_per_sec), "\n")
cat("MCMCs converged:", sum(results$converged), "of", N_simulations, "\n")


