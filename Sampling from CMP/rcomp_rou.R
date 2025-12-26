# --- corrected Rou / RoU sampler for COM-Poisson(lambda, nu) ---

log_pi_tilde <- function(x, lambda, nu) {
  x * log(lambda) - nu * lgamma(x + 1)
}

# b_x (ratio) in normal space, defined for x >= 0
b_x <- function(x, lambda, nu) {
  if (x == 0L) {
    return(exp(0.5 * (log_pi_tilde(1, lambda, nu) - log_pi_tilde(0, lambda, nu))))
  }
  # compute in log-space for stability
  log_num <- log(x + 1L) + 0.5 * log_pi_tilde(x + 1L, lambda, nu)
  log_den <- log(x)      + 0.5 * log_pi_tilde(x, lambda, nu)
  exp(log_num - log_den)
}

# find mode x0: smallest x such that \tilde\pi(x+1) <= \tilde\pi(x)
find_mode <- function(lambda, nu, start = NULL) {
  if (is.null(start)) start <- max(0L, floor(lambda^(1 / nu)))
  x <- as.integer(start)
  # increase until mass stops increasing
  while (TRUE) {
    logp_x   <- log_pi_tilde(x,   lambda, nu)
    logp_xp1 <- log_pi_tilde(x+1L, lambda, nu)
    if (logp_xp1 <= logp_x) break
    x <- x + 1L
  }
  # step backward if we overshot (shouldn't happen given the while)
  while (x > 0L) {
    logp_prev <- log_pi_tilde(x-1L, lambda, nu)
    if (logp_prev > logp_x) {
      x <- x - 1L
      logp_x <- logp_prev
    } else break
  }
  x
}

# find x* = minimal x >= 0 such that b_x <= 1
find_x_star <- function(lambda, nu, start = NULL) {
  if (is.null(start)) start <- max(1L, floor(lambda^(1 / nu)))
  x <- as.integer(start)
  
  while (b_x(x, lambda, nu) > 1) {
    x <- x + 1L

    if (x > 1e8) stop("find_x_star: exceeded search limit")
  }
  
  while (x > 0L && b_x(x - 1L, lambda, nu) <= 1) {
    x <- x - 1L
  }
  
  x
}

# Rou / RoU rejection sampler for n samples
rcomp_rou <- function(n, lambda, nu) {
  if (n <= 0) return(integer(0))
  
  # --- setup ---
  x0 <- find_mode(lambda, nu)
  A <- sqrt(exp(log_pi_tilde(x0, lambda, nu)))   # A = sqrt(pi~(x0))
  
  x_star <- find_x_star(lambda, nu)

  x_star_for_B <- max(1L, x_star)
  B <- x_star_for_B * sqrt(exp(log_pi_tilde(x_star_for_B, lambda, nu)))
  
  # ensure A and B are positive numeric
  if (!is.finite(A) || A <= 0) stop("Bad A (non-finite or <=0)")
  if (!is.finite(B) || B <= 0) stop("Bad B (non-finite or <=0)")
  
  out <- integer(n)
  i <- 1L
  
  while (i <= n) {
    U <- runif(1, 0, A)
    V <- runif(1, 0, B)
    
    # Candidate by floor, not round
    candidate <- as.integer(floor(V / U))
    
    # candidate must be >= 0 (floor can give -Inf if U==0, but runif(0, A) won't be exactly 0 with prob 1)
    if (candidate < 0L) next
    
    # accept if U <= sqrt(pi~(candidate))
    if (U <= sqrt(exp(log_pi_tilde(candidate, lambda, nu)))) {
      out[i] <- candidate
      i <- i + 1L
    }
  }
  
  out
}
