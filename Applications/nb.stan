data {
  int<lower=1> N;                 // number of observations
  array[N] int<lower=0> y;        // count data
  real leps;                      // kept for compatibility
}

parameters {
  real<lower=0> mu;               // mean
  real<lower=0> phi;              // dispersion (larger = closer to Poisson)
}

model {
  // weakly informative priors
  mu  ~ gamma(0.1, 0.01);
  phi ~ gamma(0.1, 0.01);

  for (n in 1:N)
    y[n] ~ neg_binomial_2(mu, phi);
}
