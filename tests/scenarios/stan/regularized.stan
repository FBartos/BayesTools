data {
  int<lower=1> N;
  int<lower=1> K;
  matrix[N, K] X;
  vector[N] y;
  int<lower=1, upper=4> model_id;  // ridge, lasso, horseshoe, finite-slab horseshoe
}
transformed data {
  int L_dim = model_id >= 3 ? K : 0;
  matrix[K, K] gram = crossprod(X);
}
parameters {
  real alpha;
  real<lower=0> tau;
  vector[model_id <= 2 ? K : 0] z;
  vector<lower=0>[model_id >= 3 ? K : 0] lambda;
  vector<lower=0>[model_id == 4 ? 1 : 0] slab_aux;
}
transformed parameters {
  vector[L_dim] local_sd = tau * lambda;
  matrix[L_dim, L_dim] L;
  vector[L_dim] proj;
  if (model_id >= 3) {
    if (model_id == 4) {
      // Piironen & Vehtari (2017), Appendix C.1: c^2 = slab_scale^2 * caux.
      real c = 2 * sqrt(slab_aux[1]);
      for (k in 1:K) local_sd[k] = c * local_sd[k] / hypot(c, local_sd[k]);
    }
    L = cholesky_decompose(add_diag(gram, inv_square(local_sd)));
    proj = mdivide_left_tri_low(L, X' * (y - alpha));
  }
}
model {
  alpha ~ normal(0, 1);
  tau ~ cauchy(0, 0.2);
  if (model_id == 2) {
    z ~ double_exponential(0, 1);
  } else {
    z ~ std_normal();
  }
  lambda ~ cauchy(0, 1);
  slab_aux ~ inv_gamma(2, 2);
  if (model_id <= 2) {
    y ~ normal(alpha + X * (tau * z), 1);
  } else {
    // Integrate Gaussian coefficients out conditional on the shrinkage scales.
    target += -sum(log(local_sd)) - sum(log(diagonal(L)))
              + 0.5 * (dot_self(proj) - dot_self(y - alpha));
  }
}
generated quantities {
  vector[K] beta;
  if (model_id <= 2) {
    beta = tau * z;
  } else {
    // Draw coefficients from the exact conditional posterior.
    vector[K] noise;
    for (k in 1:K) noise[k] = normal_rng(0, 1);
    beta = mdivide_right_tri_low((proj + noise)', L)';
  }
}
