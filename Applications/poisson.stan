data {
  int<lower=1> N;          // number of observations
  array[N] int<lower=0> y;       // count data
  real leps;           // just to match with compoisson data
}

parameters {
  real<lower=0> lambda;
}

model {
  // weakly informative priors
  lambda     ~ gamma(0.1, 0.01);

  for (n in 1:N)
    y[n] ~ poisson(lambda);
}
