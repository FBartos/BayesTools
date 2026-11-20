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
  target += normal_lpdf(alpha | 0, 1);
  target += cauchy_lpdf(tau | 0, 0.2) + log(2);
  target += double_exponential_lpdf(z | 0, 1);
  target += normal_lpdf(y | alpha + X * beta, 1);
}
