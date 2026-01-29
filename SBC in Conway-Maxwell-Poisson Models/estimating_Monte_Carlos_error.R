source("rcomp.R")

N <- 10^6
mu_values <- c(0.5, 1, 2)
nu_values <- c(0.5, 1, 2, 5) # nu=1 is Poisson, nu<1 is under-dispersed, nu>1 is over-dispersed

max_var <- Inf
max_mu <- Inf
max_nu <- Inf

for (mu in mu_values) {
  for (nu in nu_values) {
    X_list <- rcomp(N, mu^nu, nu)
    
    var_x <- var(X_list)
    
    if (max_var > var_x) {
      max_var <- var_x
      max_mu <- mu
      max_nu <- nu
    }
  }
}

max_var
max_mu
max_nu