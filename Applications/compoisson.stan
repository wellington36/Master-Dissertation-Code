functions {
  real log_k_term(real log_mu, real nu, int k) {
    return nu * (k * log_mu - lgamma(k + 1));
  }
  
  real bound_remainder(real k_current_term, real k_previous_term) {
    return k_current_term - log(- expm1(k_current_term - k_previous_term));
  }

  int stopping_criterio_bucket(real k_current_term, real k_previous_term, int k, real leps) {
    if (k % 5 == 0) {
      return (bound_remainder(k_current_term, k_previous_term) >= leps);
    }
    return (1e300 >= leps); // Int > leps
  }

  real log_Z_com_poisson(real log_mu, real nu, real leps) {
    real log_Z;
    int k = 2;
    int M = 10000;
    vector[M] log_Z_terms;

    if (nu == 1) {
      return exp(log_mu);
    }
    // nu == 0 or Inf will fail in this parameterization
    if (nu <= 0) {
      reject("nu must be positive");
    }
    if (nu == positive_infinity()) {
      reject("nu must be finite");
    }

    // direct computation of the truncated series
    // check if the Mth term of the series pass in the stopping criteria
    if (bound_remainder(log_k_term(log_mu, nu, M),
                        log_k_term(log_mu, nu, M - 1)) >= leps) {
      reject("nu is too close to zero.");
    }

    // first 2 terms of the series
    log_Z_terms[1] = log_k_term(log_mu, nu, 0);
    log_Z_terms[2] = log_k_term(log_mu, nu, 1);

    while (((log_Z_terms[k] > log_Z_terms[k-1]) ||
      (stopping_criterio_bucket(log_Z_terms[k], log_Z_terms[k-1], k, leps))) &&
      k < M) {
      k += 1;
      log_Z_terms[k] = log_k_term(log_mu, nu, k-1);
    }
    log_Z = log_sum_exp(log_Z_terms[1:(k-1)]);

    return log_Z;
  }

  real compoisson_lpmf(int y, real mu, real nu, real leps) {
    real logZ = log_Z_com_poisson(log(mu), nu, leps);
    return nu * (y * log(mu) - lgamma(y + 1)) - logZ;
  }
}

data {
  int<lower=1> N;          // number of observations
  array[N] int<lower=0> y;       // count data
  real leps;          // truncation level
}

parameters {
  real<lower=0> mu;
  real<lower=1e-2> nu;
}

model {
  // weakly informative priors
  mu     ~ gamma(0.1, 0.01);
  nu     ~ gamma(1, 1);

  for (n in 1:N)
    y[n] ~ compoisson(mu, nu, leps);
}
