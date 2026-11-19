data {
  int<lower=1> N;
  int<lower=3> K;  // K-1 centered factor coefficients, followed by one numeric slope
  matrix[N, K] X;
  vector[N] y;
  int<lower=1, upper=4> family;
}
transformed data {
  int gaussian = family != 2;
  int n_gaussian = gaussian ? K : 0;
  matrix[K, K] gram = crossprod(X);
}
parameters {
  real alpha;
  real<lower=0> tau;
  vector[gaussian ? 0 : K] laplace_beta;
  vector<lower=0>[family >= 3 ? 2 : 0] lambda;
  vector<lower=0>[family == 4 ? 1 : 0] slab_aux;
}
transformed parameters {
  matrix[n_gaussian, n_gaussian] L;
  vector[n_gaussian] proj;
  vector[n_gaussian] s = rep_vector(tau, n_gaussian);
  if (gaussian) {
    if (family >= 3) {
      s[1:(K-1)] *= lambda[1];
      s[K] *= lambda[2];
      if (family == 4) {
        real c = 2 * sqrt(slab_aux[1]);
        for (k in 1:K) s[k] = c * s[k] / hypot(c, s[k]);
      }
    }
    L = cholesky_decompose(add_diag(gram, inv_square(s)));
    proj = mdivide_left_tri_low(L, X' * (y - alpha));
  }
}
model {
  alpha ~ normal(0, 1);
  tau ~ cauchy(0, 0.2);
  lambda ~ cauchy(0, 1);
  slab_aux ~ inv_gamma(2, 2);
  if (gaussian) {
    // Integrate out Gaussian coefficients conditional on their scales.
    target += -sum(log(s)) - sum(log(diagonal(L)))
              + 0.5 * (dot_self(proj) - dot_self(y - alpha));
  } else {
    // The spherical density on K-1 contrasts has scale factor tau^(-(K-1)).
    target += -(K-1) * log(tau) - norm2(head(laplace_beta, K-1)) / tau;
    laplace_beta[K] ~ double_exponential(0, tau);
    y ~ normal(alpha + X * laplace_beta, 1);
  }
}
generated quantities {
  vector[K] beta;
  if (gaussian) {
    vector[K] noise;
    for (k in 1:K) noise[k] = normal_rng(0, 1);
    beta = mdivide_right_tri_low((proj + noise)', L)';
  } else {
    beta = laplace_beta;
  }
}
