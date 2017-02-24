# noncentered param, varying intercept only
data {
  int n;
  int J;
  vector[n] x;
  vector[n] y;
  int grp_id[n];
  real sigma_y;
}
parameters {
  real beta1;
  real mu;
  real<lower=0> tau;
  vector[J] eta;
}
transformed parameters {
  vector[J] beta0;
  vector[n] theta;
  beta0 = mu + tau * eta;
  for (i in 1:n) {
    theta[i] = beta0[grp_id[i]] + beta1 * x[i];
  }
}
model {
  mu ~ normal(0, 1);
  tau ~ cauchy(0, 2.5);
  eta ~ normal(0, 1);
  beta1 ~ normal(0, 1);
  y ~ normal(theta, sigma_y);
}
