# COM–Poisson term
a_n <- function(n, mu, nu) {
  (mu^n / factorial(n))^nu
}

# Partial sum Z_n
Z_n <- function(n, mu, nu) {
  i <- 0:n
  sum((mu^i / factorial(i))^nu)
}

# "True" value Z (high-accuracy reference)
# adaptive truncation until terms are negligible
Z_true <- function(mu, nu, tol = 1e-16) {
  s <- 0
  n <- 0
  repeat {
    term <- a_n(n, mu, nu)
    s <- s + term
    if (term < tol && n > mu) break
    n <- n + 1
    if (n > 500) break  # safety
  }
  s
}

# Find truncation index n via sum-to-threshold
# (simple and safe for COM–Poisson)
find_n <- function(mu, nu, eps) {
  n <- 0
  repeat {
    term <- a_n(n, mu, nu)
    if (term < eps) return(n)
    n <- n + 1
    if (n > 300) return(NA)
  }
}

# Bound for |1/Z - 1/Z_n|
# derived from ratio bounding pair
bound_inverse_Z <- function(n, mu, nu, Zn) {
  tail_bound <- mu^nu / ((n + 1)^nu - mu^nu)
  denom <- Zn^2 - tail_bound^2
  if (denom <= 0) return(NA)
  tail_bound / denom
}

# TEST GRID
mu_values <- c(0.5, 1, 2, 5, 10)
nu_values <- c(0.7, 1.0, 1.3)
eps_values <- c(1e-6, 1e-10)

cat("Checking bound for |1/Z - 1/Z_n|\n\n")

results <- data.frame()

for (mu in mu_values) {
  for (nu in nu_values) {
    
    Zinf <- Z_true(mu, nu)
    
    for (eps in eps_values) {
      
      n <- find_n(mu, nu, eps)
      if (is.na(n)) next
      
      Zn <- Z_n(n, mu, nu)
      
      actual_err <- abs(1 / Zinf - 1 / Zn)
      bound_val <- bound_inverse_Z(n, mu, nu, Zn)
      
      valid <- actual_err <= bound_val
      
      cat(sprintf(
        "mu=%4.1f  nu=%4.1f  eps=%1.0e  n=%3d  actual=%9.2e  bound=%9.2e  OK=%s\n",
        mu, nu, eps, n, actual_err, bound_val, valid
      ))
      
      results <- rbind(
        results,
        data.frame(mu, nu, eps, n, actual_err, bound_val, valid)
      )
    }
  }
}

# SUMMARY
cat("\nSummary:\n")
cat(sprintf(
  "Passed: %d / %d tests\n",
  sum(results$valid, na.rm = TRUE),
  nrow(results)
))
