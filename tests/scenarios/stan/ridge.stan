data {
  int<lower=1> N;
  int<lower=1> K;
  matrix[N, K] X;
  vector[N] y;
}
parameters {
  real alpha;
  real<lower=0> tau;
  vector[K] z;
}
transformed parameters {
  vector[K] beta = tau * z;
}
model {
  // Keep the normalizing constants for marginal likelihoods.
  target += normal_lpdf(alpha | 0, 1);
  target += cauchy_lpdf(tau | 0, 0.2) + log(2);
  target += std_normal_lpdf(z);
  target += normal_lpdf(y | alpha + X * beta, 1);
}
