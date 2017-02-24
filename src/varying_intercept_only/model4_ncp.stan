# noncentered, varying intercept with unit-level predictor, with group-level predictor
# yi ~ N(beta0j[i] + beta1 * x, sigma_y)
# beta0j ~ N(mu0 + gamma0 * w, tau0)
data {
  int n;
  int J;
  int grp_id[n];
  vector[n] y;
  vector[n] x;
  vector[J] w;
}
parameters {
  vector[J] eta0;
  real mu0;
  real gamma0;
  real<lower=0> tau0;
  real beta1;
  real<lower=0> sigma_y;
}
transformed parameters {
  vector[J] beta0;
  vector[n] ymean;
  beta0 = mu0 + gamma0 * w + eta0 * tau0;
  for (i in 1:n) {
    ymean[i] = beta0[grp_id[i]] + beta1 * x[i];
  }
}
model {
  mu0 ~ normal(0, 1);
  gamma0 ~ normal(0, 1);
  tau0 ~ cauchy(0, 2.5);
  eta0 ~ normal(0, 1);
  beta1 ~ normal(0, 1);
  sigma_y ~ cauchy(0, 2.5);
  y ~ normal(ymean, sigma_y);
}
