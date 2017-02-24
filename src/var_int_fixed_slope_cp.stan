# centered param, varying intercepts only
data {
  int n;
  int J;
  vector[n] x;
  vector[n] y;
  int grp_id[n];
  real sigma_y;
}
parameters {
  vector[J] beta0;
  real beta1;
  real mu;
  real<lower=0> tau;
}
transformed parameters {
  vector[n] theta;
  for (i in 1:n) {
    theta[i] = beta0[grp_id[i]] + beta1 * x[i];
  }
}
model {
  mu ~ normal(0, 1);
  tau ~ cauchy(0, 2.5);
  beta0 ~ normal(mu, tau);
  beta1 ~ normal(0, 1);
  y ~ normal(theta, sigma_y);
}
