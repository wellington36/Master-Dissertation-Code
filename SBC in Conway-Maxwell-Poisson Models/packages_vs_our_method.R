source("rcomp_benson.R")
source("rcomp.R")

library(COMPoissonReg)


# Functions from discontinued packages)
com_log_sum <- function(x, y) {
  if (x == -Inf) return(y) else if (y == -Inf) return(x)
  if (x > y) return(x + log(1 + exp(y - x))) else return(y + log(1 + exp(x - y)))
}

# Implementation from 'compoisson' package
rcom_legacy <- function(n, lambda, nu) {
  log_z <- -Inf; z_last <- 0; j <- 0
  while (abs(log_z - z_last) > 0.001) {
    z_last <- log_z
    log_z <- com_log_sum(log_z, j * log(lambda) - nu * lfactorial(j))
    j <- j + 1
  }
  r <- numeric(n)
  for (i in 1:n) {
    log_u <- log(runif(1))
    j <- 0; curr_log_cum <- -Inf
    while (TRUE) {
      log_pj <- (j * log(lambda) - nu * lfactorial(j)) - log_z
      if (log_u < log(exp(curr_log_cum) + exp(log_pj))) break
      curr_log_cum <- com_log_sum(curr_log_cum, log_pj)
      j <- j + 1
      if (j > 1000) break
    }
    r[i] <- j
  }
  return(r)
}

# Implementation from 'CompGLM' package
rcomp_legacy <- function(n, lambda, nu, sumTo = 100) {
  j <- 0:sumTo
  log_probs <- j * log(lambda) - nu * lfactorial(j)
  probs <- exp(log_probs - max(log_probs))
  probs <- probs / sum(probs)
  return(sample(0:sumTo, size = n, replace = TRUE, prob = probs))
}

# Threshold with high number of terms
threshold <- function(lambda, nu, sumTo = 100000) {
  j <- 0:sumTo
  log_z <- j * log(lambda) - nu * lfactorial(j)
  p <- exp(log_z - max(log_z))
  p <- p / sum(p)
  m <- sum(j * p)
  v <- sum((j^2) * p) - (m^2)
  return(list(mean = m, var = v))
}

# Evaluation
evaluate <- function(n = 5000, lambda, nu) {
  threshold <- threshold(lambda, nu)
  
  # Measure time and generate samples
  cat("Running benchmarks...\n")
  
  t1 <- system.time(s1 <- rcom_legacy(n, lambda, nu))["elapsed"]
  t2 <- system.time(s2 <- rcomp_legacy(n, lambda, nu))["elapsed"]
  
  # COMPoissonReg check
  if(requireNamespace("COMPoissonReg", quietly = TRUE)) {
    t3 <- system.time(s3 <- COMPoissonReg::rcmp(n, lambda, nu))["elapsed"]
  } else {
    s3 <- rep(NA, n); t3 <- NA
  }
  
  # Benson function check
  if(exists("rcomp_benson")) {
    mu <- lambda^(1/nu)
    t4 <- system.time(s4 <- rcomp_benson(n, mu, nu))["elapsed"]
  } else {
    s4 <- rep(NA, n); t4 <- NA
  }

  # Our function check
  if(exists("rcomp")) {
    t5 <- system.time(s5 <- rcomp(n, lambda, nu))["elapsed"]
  } else {
    s5 <- rep(NA, n); t5 <- NA
  }
  
  # Create results table
  results <- data.frame(
    Method = c("Threshold", "compoisson", "CompGLM", "COMPoissonReg", "rcomp_benson", "rcomp"),
    Mean = c(threshold$mean, mean(s1), mean(s2), mean(s3, na.rm=T), mean(s4, na.rm=T), mean(s5, na.rm=T)),
    Var  = c(threshold$var,  var(s1),  var(s2),  var(s3, na.rm=T),  var(s4, na.rm=T), var(s5, na.rm=T)),
    Time_Sec = c(0, t1, t2, t3, t4, t5)
  )
  
  # Add error columns
  results$Mean_Err <- abs(results$Mean - threshold$mean)
  results$Var_Err  <- abs(results$Var - threshold$var)
  
  # Reorder columns for better reading
  results <- results[, c("Method", "Mean", "Mean_Err", "Var", "Var_Err", "Time_Sec")]
  
  print(results, digits = 4)
}


evaluate(n = 1000000, lambda = 2, nu = 0.15)