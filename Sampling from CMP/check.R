#set.seed(121)
source("com_poisson_pmf.R")
source("rcomp_exact.R")       # exact sampler
source("rcomp_rou.R")   # RoU sampler
source("rcomp_rejection_benson.R") # rejection sampler

# Dummy sampler for test
rcomp_dummy <- function(n, lambda, nu) {
  rpois(n, lambda)  # just Poisson as placeholder
}

# --- Register methods in a list ---
methods <- list(
  rou        = rcomp_rou,
  rejection_benson = rcomp_rejection_benson,
  exact            = rcomp_exact,
  dummy            = rcomp_dummy
)

# --- Benchmark wrapper: return multiple metrics ---
benchmark_method <- function(method_fun, n_iter = 200, n_sample = 10000, eps = 1e-8) {
  times  <- numeric(n_iter)
  pvals  <- numeric(n_iter)
  diffs  <- numeric(n_iter)
  
  for (i in seq_len(n_iter)) {
    lambda_i <- runif(1, 0.5, 20)
    nu_i     <- runif(1, 0.5, 5)
    #nu_i <- 1  # Poisson
    log_lambda <- log(lambda_i)
    
    # --- timing + sampling ---
    t0 <- Sys.time()
    samples <- method_fun(n = n_sample, lambda = lambda_i, nu = nu_i)
    t1 <- Sys.time()
    times[i] <- as.numeric(difftime(t1, t0, units = "secs"))
    
    # --- sample mean ---
    sample_mean <- mean(samples)
    
    # --- true distribution (truncate at max obs + margin) ---
    K <- max(samples) + 50
    support <- 0:K
    log_p <- sapply(support, function(y)
      com_poisson_log_lpmf(y, log_lambda, nu_i, eps))
    p <- exp(log_p); p <- p / sum(p)
    
    true_mean <- sum(support * p)
    diffs[i] <- abs(sample_mean - true_mean)
    
    # --- chi-square p-value ---
    obs <- table(factor(samples, levels = support))
    pvals[i] <- tryCatch({
      suppressWarnings(chisq.test(obs, p = p)$p.value)
    }, error = function(e) NA)
  }
  
  # Return named metrics
  c(
    Time = mean(times, na.rm = TRUE),
    Chi2Pval = mean(pvals, na.rm = TRUE),
    MeanDiff = mean(diffs, na.rm = TRUE)
  )
}


# --- Run all methods ---
eval <- function() {
  results_list <- lapply(methods, benchmark_method)
  results <- do.call(cbind, results_list)
  colnames(results) <- names(methods)

  # --- Format table for display ---
  tab <- rbind(
    c("", colnames(results)),
    c("Time(s)", sprintf("%.4f", results["Time", ])),
    c("Chi2 p-val", sprintf("%.4f", results["Chi2Pval", ])),
    c("|Mean - TrueMean|", sprintf("%.4f", results["MeanDiff", ]))
  )

  print(tab, quote = FALSE)
}

eval()