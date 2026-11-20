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
  alpha ~ normal(0, 1);
  tau ~ cauchy(0, 0.2);
  z ~ double_exponential(0, 1);
  y ~ normal(alpha + X * beta, 1);
}
