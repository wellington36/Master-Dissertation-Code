source("rcomp_benson.R")

cmp_hyper_check <- function(N,
                            alpha_mu, beta_mu,
                            alpha_nu, beta_nu) {
  
  Y <- numeric(N)
  
  for (i in seq_len(N)) {
    mu_i <- rgamma(1, shape = alpha_mu, rate = beta_mu)
    nu_i <- rgamma(1, shape = alpha_nu, rate = beta_nu)
    lambda_i <- mu_i^nu_i
    Y[i] <- rcomp_benson(lambda = lambda_i, nu = nu_i, n = 1)
  }
  
  emp_mean <- mean(Y)
  emp_var  <- var(Y)
  
  # ============================
  # Layout: 1 row, 2 columns
  # ============================
  par(mfrow = c(1, 2))
  
  # ============================
  # Empirical PDF
  # ============================
  pdf_tab <- table(Y) / length(Y)
  y_vals  <- as.numeric(names(pdf_tab))
  pdf_vals <- as.numeric(pdf_tab)
  
  plot(y_vals, pdf_vals,
       type = "h",
       lwd = 3,
       col = "darkred",
       main = "Empirical PDF of Y",
       xlab = "y",
       ylab = "P(Y = y)")
  points(y_vals, pdf_vals,
         pch = 19,
         col = "darkred")
  grid()
  
  # ============================
  # Empirical CDF
  # ============================
  plot(ecdf(Y),
       main = "Empirical CDF of Y",
       xlab = "y",
       ylab = "F_Y(y)",
       col = "blue",
       lwd = 2)
  grid()
  
  # Reset layout
  par(mfrow = c(1, 1))
  
  # Print results
  cat("\n----- Results -----\n")
  cat("Empirical mean     =", emp_mean, "\n")
  cat("Empirical variance =", emp_var, "\n")
  cat("-------------------\n")
  
  invisible(list(
    Y = Y,
    pdf = pdf_tab,
    cdf = ecdf(Y),
    emp_mean = emp_mean,
    emp_var = emp_var
  ))
}


out <- cmp_hyper_check(
  N = 5000,
  alpha_mu = 0.01, beta_mu = 0.1,
  alpha_nu = 1, beta_nu = 1
)
