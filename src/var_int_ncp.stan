# no matt trick, varying intercept only
data {
  int n;
  int J;
  //vector[n] x;
  vector[n] y;
  int grp_id[n];
  real sigma_y;
}
parameters {
  //real beta1;
  vector[J] eta;
  real mu;
  real<lower=0> tau;
  //real<lower=0> sigma_y;
}
transformed parameters {
  vector[J] beta0;
  vector[n] theta;
  beta0 = mu + tau * eta;
  for (i in 1:n) {
    theta[i] = beta0[grp_id[i]]; // + beta1 * x[i];
  }
}
model {
  //beta1 ~ normal(0, 1);
  mu ~ normal(0, 1);
  tau ~ cauchy(0, 2.5);
  //sigma_y ~ cauchy(0, 2.5);
  //for (i in 1:n) {
  //  y[i] ~ normal(beta0[grp_id[i]] + beta1 * x[i], sigma_y);
  //}
  //target += normal_lpdf(eta | 0, 1);
  //target += normal_lpdf(y | theta, sigma_y);
  eta ~ normal(0, 1);
  y ~ normal(theta, sigma_y);
}
