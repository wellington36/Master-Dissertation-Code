source("com_poisson_pmf.R")

log1mexp <- function(d) {
  # computes log(1 - exp(d)) for d <= 0
  stopifnot(d <= 0)
  
  if (d <= -log(2)) {
    # safe because exp(d) is small enough
    return(log1p(-exp(d)))
  } else {
    # use expm1 for precision near 0
    return(log(-expm1(d)))
  }
}

# Exact sampler for n samples
rcomp_exact <- function(n, lambda, nu) {
  log_lambda <- log(lambda)
  
  # eval blindly until the maximum
  k <- 1
  a <- c(log_k_term(log_lambda, nu, k))
  
  while ((k < (floor(lambda^(1/nu)) + 2))) {
    k <- k + 1
    a <- append(a, log_k_term(log_lambda, nu, k))
  }
  
  out <- numeric(n)
  i <- 1
  while (i <= n) {
    z <- runif(1)
    logz <- log(z)
    j <- 1

    log_F_j <- logsumexp(a[1:j])
    log_F_k <- logsumexp(a[1:k])

    log_L_k <- log_F_k
    log_U_k <- logsumexp(c(log_F_k, a[k] - log1mexp(a[k] - a[k-1])))

 
    log_F_j_k_lower <- log_F_j - log_U_k
    log_F_j_k_upper <- log_F_j - log_L_k

    while (log_F_j_k_lower < logz) {
      if (log_F_j_k_lower < logz && logz < log_F_j_k_upper) {
        k <- k + 1
      
        a <- append(a, log_k_term(log_lambda, nu, k))
        log_F_k <- logsumexp(a[1:k])
      
        log_L_k <- log_F_k
        log_U_k <- logsumexp(c(log_F_k, a[k] - log1mexp(a[k] - a[k-1])))
      } else if (logz > log_F_j_k_upper) {
        j <- j + 1
      
        log_F_j <- logsumexp(a[1:j])
      }
      
      log_F_j_k_lower <- log_F_j - log_U_k
      log_F_j_k_upper <- log_F_j - log_L_k
    }
    
    out[i] <- j - 1 # -1 because R starts with 1 but we need to start with 0
    i <- i + 1
  }
  return(out)
}
