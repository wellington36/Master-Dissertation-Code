  // log of kth term of the normalizing series of the COM Poisson distribution
  // Args:
  //   log_mu: log location parameter
  //   shape: positive shape parameter
  //   k: k-th term
  real log_k_term(real log_mu, real nu, int k) {
    return nu * ((k - 1) * log_mu - lgamma(k));
  }

  // bound for the remainder of the normalizing series of the COM Poisson
  // distribution given the last two terms in log-scale
  // Args:
  //   k_current_term: the log of a_k term
  //   k_previous_term: the log of a_(k-1) term
  real bound_remainder(real k_current_term, real k_previous_term) {
    return k_current_term - log(- expm1(k_current_term - k_previous_term));
  }

  // stopping criterio with bucket
  // Args:
  //   k_current_term: the log of a_k term
  //   k_previous_term: the log of a_(k-1) term
  //   k: k term of the series
  //   leps: log(eps)
  int stopping_criterio_bucket(real k_current_term, real k_previous_term, int k, real leps) {
    if (k % 10 == 0) {
      return (bound_remainder(k_current_term, k_previous_term) >= leps);
    }
    return (1e300 >= leps); // Int > leps
  }

  // log normalizing constant of the COM Poisson distribution
  // implementation inspired by code of Ben Goodrich
  // improved following suggestions of Sebastian Weber (#892)
  // Args:
  //   log_mu: log location parameter
  //   shape: positive shape parameter
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
    log_Z_terms[1] = log_k_term(log_mu, nu, 1);
    log_Z_terms[2] = log_k_term(log_mu, nu, 2);

    while (((k < exp(log_mu) + 2) ||
      (stopping_criterio_bucket(log_Z_terms[k], log_Z_terms[k-1], k, leps))) &&
      k < M) {
      k += 1;
      log_Z_terms[k] = log_k_term(log_mu, nu, k);
    }
    log_Z = log_sum_exp(log_Z_terms[1:k]);

    return log_Z;
  }
  // COM Poisson log-PMF for a single response (log parameterization)
  // Args:
  //   y: the response value
  //   log_mu: log location parameter
  //   shape: positive shape parameter
  real com_poisson_log_lpmf(int y, real log_mu, real nu) {
    real current_leps = leps_custom();
    if (nu == 1) return poisson_log_lpmf(y | log_mu);
    return nu * (y * log_mu - lgamma(y + 1)) - log_Z_com_poisson(log_mu, nu, current_leps);
  }
  // COM Poisson log-PMF for a single response
  real com_poisson_lpmf(int y, real mu, real nu) {
    if (nu == 1) return poisson_lpmf(y | mu);
    return com_poisson_log_lpmf(y | log(mu), nu);
  }
  // COM Poisson log-CDF for a single response
  real com_poisson_lcdf(int y, real mu, real nu) {
    real current_leps = leps_custom();
    real log_mu;
    real log_Z;  // log denominator
    vector[y] log_num_terms;  // terms of the log numerator
    if (nu == 1) {
      return poisson_lcdf(y | mu);
    }
    // nu == 0 or Inf will fail in this parameterization
    if (nu <= 0) {
      reject("nu must be positive");
    }
    if (nu == positive_infinity()) {
      reject("nu must be finite");
    }
    if (y > 10000) {
      reject("cannot handle y > 10000");
    }
    log_mu = log(mu);
    if (nu * (y * log_mu - lgamma(y + 1)) <= -36.0) {
      // y is large enough for the CDF to be very close to 1;
      return 0;
    }
    log_Z = log_Z_com_poisson(log_mu, nu, current_leps);
    if (y == 0) {
      return -log_Z;
    }
    // first 2 terms of the series
    log_num_terms[1] = log1p_exp(nu * log_mu);
    // remaining terms of the series until y
    for (k in 2:y) {
      log_num_terms[k] = nu * (k * log_mu - lgamma(k + 1));
    }
    return log_sum_exp(log_num_terms) - log_Z;
  }
  // COM Poisson log-CCDF for a single response
  real com_poisson_lccdf(int y, real mu, real nu) {
    return log1m_exp(com_poisson_lcdf(y | mu, nu));
  }
