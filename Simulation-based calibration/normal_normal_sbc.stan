data {
  int<lower=0> N;
  vector[N] y;
  
  real<lower=0> sigma_y;
  real prior_mu_mean;
  real<lower=0> prior_mu_sd;
}

parameters {
  real mu;
}

model {
  mu ~ normal(prior_mu_mean, prior_mu_sd);
  y ~ normal(mu, sigma_y);
}


generated quantities {
}
