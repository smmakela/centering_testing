// noncentered param, varying intercepts and slopes, with correlation
// beta ~ multi_normal(mu, Tau),
// where Tau = diag_matrix(tau) * L * L' * diag_matrix(tau) and
// L is the cholesky factor of the correlation matrix
// we do the noncentered parameterization by doing beta = mu + Lz, we just
// have to replicate the transpose of mu so that it's J rows
// (rep_matrix(mu', J)), pre-multiply tau and L so that it's
// actually diag(tau) * L * L' * diag(tau) (diag_pre_multiply(tau, L)), then
// multiply that by z (diag_pre_multiply(tau, L) * z), and take the transpose so
// that it's the right shape ((diag_pre_multiply(tau, L) * z)')
data {
  int n;
  int J;
  vector[n] x;
  vector[n] y;
  int grp_id[n];
  real sigma_y;
}
parameters {
  vector[2] mu;
  matrix[2, J] z;
  vector<lower=0> tau;
  cholesky_factor_corr[2] L;
}
transformed parameters {
  matrix[J, 2] beta;
  vector[n] theta;
  beta = rep_matrix(mu', J) + (diag_pre_multiply(tau, L) * z)';
  for (i in 1:n) {
    theta[i] = beta[grp_id[i], 1] + beta[grp_id[i], 2] * x[i];
  }
}
model {
  mu ~ normal(0, 1);
  to_vector(z) ~ normal(0, 1);
  tau ~ cauchy(0, 2.5);
  L ~ lkj_corr_cholesky(2);
  y ~ normal(theta, sigma_y);
}
generated quantities {
  matrix[2, 2] Rho;
  Rho = multiply_lower_tri_self_transpose(L);
}
