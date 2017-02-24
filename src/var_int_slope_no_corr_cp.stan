# centered param, varying intercepts and slopes, independent
data {
  int n;
  int J;
  vector[n] x;
  vector[n] y;
  int grp_id[n];
  real sigma_y;
}
parameters {
  real mu;
  real eta;
  real<lower=0> tau;
  real<lower=0> omega;
  vector[J] beta0;
  vector[J] beta1;
}
transformed parameters {
  vector[n] theta;
  for (i in 1:n) {
    theta[i] = beta0[grp_id[i]] + beta1[grp_id[i]] * x[i];
  }
}
model {
  mu ~ normal(0, 1);
  tau ~ cauchy(0, 2.5);
  beta0 ~ normal(mu, tau);
  eta ~ normal(0, 1);
  omega ~ cauchy(0, 2.5);
  beta1 ~ normal(eta, omega);
  y ~ normal(theta, sigma_y);
}
