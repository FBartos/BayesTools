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
  target += normal_lpdf(alpha | 0, 1);
  target += cauchy_lpdf(tau | 0, 0.2) + log(2);

  // The first K-1 columns code month; the last is standardized wind.
  // Divide by sphere area times the integral of r^(d-1) exp(-r/tau).
  {
    int d = K - 1;
    target += lgamma(0.5 * d) - log(2) - 0.5 * d * log(pi()) - lgamma(d)
              - d * log(tau) - norm2(head(beta, d)) / tau;
  }
  target += double_exponential_lpdf(beta[K] | 0, tau);
  target += normal_lpdf(y | alpha + X * beta, 1);
}
