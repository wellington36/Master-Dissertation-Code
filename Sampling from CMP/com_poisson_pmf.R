# Helper function for logsumexp, as it's not native in base R
logsumexp <- function(x) {
  xmax <- max(x)
  xmax + log(sum(exp(x - xmax)))
}

log_k_term <- function(log_lambda, nu, k) {
  return((k - 1) * log_lambda - nu * lgamma(k))
}

bound_remainder <- function(k_current_term, k_previous_term) {
  return(k_current_term - log(-expm1(k_current_term - k_previous_term)))
}

stopping_criterio_bucket <- function(k_current_term, k_previous_term, k, leps, amean) {
  if (k %% max(50, amean) == 0) {
    return(bound_remainder(k_current_term, k_previous_term) >= leps)
  }
  return(TRUE)
}

log_Z_com_poisson <- function(log_lambda, nu, eps) {
  if (nu == 1) {
    return(log_lambda)
  }
  if (nu <= 0) {
    stop("nu must be positive")
  }
  if (is.infinite(nu)) {
    stop("nu must be finite")
  }
  
  M <- 10
  amean <- ceiling(exp(log_lambda)^(1/nu) - (nu - 1)/(2*nu))
  leps <- log(eps)
  
  # direct computation of the truncated series
  # check if the Mth term of the series pass in the stopping criteria
  while ((log_k_term(log_lambda, nu, M) > log_k_term(log_lambda, nu, M - 1)) ||
         (bound_remainder(log_k_term(log_lambda, nu, M),
                          log_k_term(log_lambda, nu, M - 1)) >= leps)) {
    M = 10 * M
  }
  
  log_Z_terms <- numeric(M)
  log_Z_terms[1] <- log_k_term(log_lambda, nu, 1)
  log_Z_terms[2] <- log_k_term(log_lambda, nu, 2)
  
  k <- 2
  while (((log_Z_terms[k] >= log_Z_terms[k-1]) ||
          (stopping_criterio_bucket(log_Z_terms[k],
                                    log_Z_terms[k-1],
                                    k,
                                    leps,
                                    amean))) &&
         k < M) {
    k <- k + 1
    log_Z_terms[k] <- log_k_term(log_lambda, nu, k)
  }
  log_Z <- logsumexp(log_Z_terms[1:k])
  
  return(log_Z)
}

com_poisson_log_lpmf <- function(y, log_lambda, nu, eps) {
  if (nu == 1) {
    return(dpois(y, exp(log_lambda), log = TRUE))
  }
  return(y * log_lambda - nu * lgamma(y + 1) - log_Z_com_poisson(log_lambda, nu, eps))
}