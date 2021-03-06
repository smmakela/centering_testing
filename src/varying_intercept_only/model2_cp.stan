# centered, varying intercept only, with group-level predictor
# yi ~ N(beta0j[i], sigma_y)
# beta0j ~ N(mu0 + gamma0 * w, tau0)
data {
  int n;
  int J;
  int grp_id[n];
  vector[n] y;
  vector[J] w;
}
parameters {
  vector[J] beta0;
  real mu0;
  real gamma0;
  real<lower=0> tau0;
  real<lower=0> sigma_y;
}
transformed parameters {
  vector[n] ymean;
  for (i in 1:n) {
    ymean[i] = beta0[grp_id[i]];
  }
}
model {
  mu0 ~ normal(0, 1);
  gamma0 ~ normal(0, 1);
  tau0 ~ cauchy(0, 2.5);
  beta0 ~ normal(mu0 + gamma0 * w, tau0);
  sigma_y ~ cauchy(0, 2.5);
  y ~ normal(ymean, sigma_y);
}
