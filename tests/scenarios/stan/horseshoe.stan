// beta[k] ~ normal(0, slope_sd[k]); y ~ normal(alpha + X * beta, 1).
// Integrate the slopes out when sampling the heavy-tailed scales, then draw
// beta from its conditional Gaussian posterior in generated quantities.
data {
  int<lower=1> N;
  int<lower=1> K;
  matrix[N, K] X;
  vector[N] y;
  int<lower=1> J;                // number of distinct local scales
  array[K] int<lower=1, upper=J> unit;
}
transformed data {
  matrix[K, K] gram = crossprod(X);
}
parameters {
  real alpha;
  real<lower=0> tau;
  vector<lower=0>[J] lambda;
}
transformed parameters {
  vector[K] slope_sd = tau * lambda[unit];
  matrix[K, K] L;
  vector[K] projected_y;

  L = cholesky_decompose(add_diag(gram, inv_square(slope_sd)));
  projected_y = mdivide_left_tri_low(L, X' * (y - alpha));
}
model {
  target += normal_lpdf(alpha | 0, 1);
  target += cauchy_lpdf(tau | 0, 0.2) + log(2);
  target += cauchy_lpdf(lambda | 0, 1) + J * log(2);

  // Normalized Gaussian likelihood integrated over beta.
  target += -0.5 * N * log(2 * pi()) - sum(log(slope_sd)) - sum(log(diagonal(L)))
            + 0.5 * (dot_self(projected_y) - dot_self(y - alpha));
}
generated quantities {
  vector[K] noise;
  vector[K] beta;
  for (k in 1:K) noise[k] = normal_rng(0, 1);
  beta = mdivide_right_tri_low((projected_y + noise)', L)';
}
