// noncentered param, varying intercepts and slopes, with correlation
// beta ~ multi_normal(mu, Tau),
// where Tau = diag_matrix(tau) * L * L' * diag_matrix(tau) and
// L is the cholesky factor of the correlation matrix
data {
  int n;
  int J;
  vector[n] x;
  vector[n] y;
  int grp_id[n];
  real sigma_y;
}
parameters {
  matrix[J, 2] beta;
  vector[2] mu;
  vector<lower=0> tau;
  cholesky_factor_corr[2] L; // cholesky factor of corr matrix
}
transformed parameters {
  vector[n] theta;
  for (i in 1:n) {
    theta[i] = beta[grp_id[i], 1] + beta[grp_id[i], 2] * x[i];
  }
}
model {
  mu ~ normal(0, 1);
  tau ~ cauchy(0, 2.5);
  L ~ lkj_corr_cholesky(2);
  beta ~ multi_normal(mu, diag_matrix(tau) * L * L' * diag_matrix(tau));
  y ~ normal(theta, sigma_y);
}
generated quantities {
  matrix[2, 2] Rho;
  Rho = multiply_lower_tri_self_transpose(L);
}
