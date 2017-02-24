# noncentered param, varying intercepts and slopes, independent
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
  vector[J] epsilon0;
  vector[J] epsilon1;
}
transformed parameters {
  vector[J] beta0;
  vector[J] beta1;
  vector[n] theta;
  beta0 = mu + tau * epsilon0;
  beta1 = eta + omega * epsilon1;
  for (i in 1:n) {
    theta[i] = beta0[grp_id[i]] + beta1[grp_id[i]] * x[i];
  }
}
model {
  mu ~ normal(0, 1);
  eta ~ normal(0, 1);
  tau ~ cauchy(0, 2.5);
  omega ~ cauchy(0, 2.5);
  epsilon0 ~ normal(0, 1);
  epsilon1 ~ normal(0, 1);
  y ~ normal(theta, sigma_y);
}
