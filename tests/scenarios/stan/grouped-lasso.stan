data {
  int<lower=1> N;
  int<lower=2> K;
  matrix[N, K] X;
  vector[N] y;
}
parameters {
  real alpha;
  real<lower=0> tau;
  vector[K] beta;
}
model {
  alpha ~ normal(0, 1);
  tau ~ cauchy(0, 0.2);

  // The first K-1 columns code month; the last is standardized wind.
  target += -(K - 1) * log(tau) - norm2(head(beta, K - 1)) / tau;
  beta[K] ~ double_exponential(0, tau);
  y ~ normal(alpha + X * beta, 1);
}
