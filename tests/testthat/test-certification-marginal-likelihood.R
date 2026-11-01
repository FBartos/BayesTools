skip_if_not_test_profile("certification")
skip_if_not_installed("rjags")
skip_if_not_installed("runjags")
skip_if_not_installed("bridgesampling")

# Bridge-sampling log marginal likelihoods of Gaussian linear models against
# exact values. With N(0, 1) priors on the coefficients (the coefficients of
# the mean-difference factor columns included), the coefficients and the
# random intercepts integrate out analytically,
# y ~ N(0, X X' + tau^2 Z Z' + sigma^2 I); the SDs, with half-normal(0, 1)
# priors, are integrated by adaptive quadrature to a relative error of 1e-10.
# The bridge estimates of 4 x 5000 draws are held to 0.05 of the exact values
# and to a Monte Carlo SD below 0.01 across 10 bridge seeds of each fit. In
# development the fixed-effect estimates were within 0.006 of exact, and the
# random-intercept estimates (13 bridge coordinates, or 7 with the intercepts
# integrated out) 0.005 to 0.025 below it across four fit seeds, a small bias
# that the tolerance allows twice over; the exact values agree with an
# independent grid sum (step 0.002) to 1e-5.

certification_data <- local({
  set.seed(1)
  n_groups <- 6L
  data <- data.frame(
    x = stats::rnorm(5L * n_groups),
    f = factor(rep(c("a", "b", "c"), length.out = 5L * n_groups)),
    g = factor(rep(seq_len(n_groups), each = 5L))
  )
  group_effect <- stats::rnorm(n_groups, 0, 0.5)
  data$y <- 0.3 + 0.5 * data$x + c(-0.2, 0, 0.3)[as.integer(data$f)] +
    group_effect[as.integer(data$g)] + stats::rnorm(nrow(data), 0, 0.7)
  data
})

certification_design <- cbind(
  1,
  certification_data$x,
  contr.meandif(levels(certification_data$f))[as.integer(certification_data$f), ]
)
certification_same_group <- outer(certification_data$g, certification_data$g, "==") * 1
certification_tolerance  <- 0.05

gaussian_log_density <- function(y, covariance) {
  factor <- chol(covariance)
  z <- backsolve(factor, y, transpose = TRUE)
  -0.5 * (length(y) * log(2 * pi) + sum(z^2)) - sum(log(diag(factor)))
}

# log of the integral of exp(log_f(s)) times the half-normal(0, 1) density,
# split at the mode of the integrand; the prior mass beyond 10 is below 1e-22
log_integrate_half_normal <- function(log_f) {
  log_g <- function(s) {
    vapply(s, function(value) log_f(value) + log(2) + stats::dnorm(value, log = TRUE), numeric(1))
  }
  mode <- stats::optimize(log_g, c(1e-8, 10), maximum = TRUE)
  g    <- function(s) exp(log_g(s) - mode$objective)
  value <- stats::integrate(g, 0, mode$maximum, rel.tol = 1e-10)$value +
    stats::integrate(g, mode$maximum, 10, rel.tol = 1e-10)$value
  mode$objective + log(value)
}

fixed_covariance <- tcrossprod(certification_design)

gaussian_likelihood <- "model{
  for(i in 1:N){
    y[i] ~ dnorm(mu[i], pow(sigma, -2))
  }
}
"
known_sigma_likelihood <- "model{
  for(i in 1:N){
    y[i] ~ dnorm(mu[i], pow(0.7, -2))
  }
}
"
certification_jags_data <- list(y = certification_data$y, N = nrow(certification_data))
unit_prior    <- prior("normal", list(mean = 0, sd = 1))
unit_sd_prior <- prior("normal", list(mean = 0, sd = 1), truncation = list(0, Inf))
fixed_priors  <- list(mu = list(
  intercept = unit_prior,
  x         = unit_prior,
  f         = prior_factor("mnormal", list(mean = 0, sd = 1), contrast = "meandif")
))

fit_certification_model <- function(likelihood, formula, prior_list = NULL, random = NULL) {
  JAGS_fit(
    likelihood, certification_jags_data, prior_list = prior_list,
    formula_list = list(mu = formula), formula_data_list = list(mu = certification_data),
    formula_prior_list = fixed_priors, formula_random_prior_list = random,
    chains = 4, adapt = 1000, burnin = 1000, sample = 5000, seed = 1
  )
}

log_posterior_known_sigma <- function(parameters, data) {
  sum(stats::dnorm(data$y, parameters$mu, 0.7, log = TRUE))
}
log_posterior_sigma <- function(parameters, data) {
  sum(stats::dnorm(data$y, parameters$mu, parameters$sigma, log = TRUE))
}

# the bridge estimates of a fit under several seeds
bridge_logml <- function(fit, log_posterior, prior_list = NULL, seeds = 1:10, ...) {
  vapply(seeds, function(seed) {
    JAGS_bridgesampling(fit, log_posterior = log_posterior, data = certification_jags_data,
                        prior_list = prior_list, seed = seed, ...)$logml
  }, numeric(1))
}

expect_bridge_exact <- function(logml, exact) {
  expect_lt(stats::sd(logml), 0.01)
  expect_equal_each(logml, rep(exact, length(logml)), tolerance = certification_tolerance)
}

test_that("bridge sampling recovers the exact log marginal likelihood with known sigma", {

  fit   <- fit_certification_model(known_sigma_likelihood, ~ x + f)
  exact <- gaussian_log_density(certification_data$y, fixed_covariance + 0.7^2 * diag(nrow(certification_data)))

  expect_bridge_exact(bridge_logml(fit, log_posterior_known_sigma), exact)
})

test_that("bridge sampling recovers the exact log marginal likelihood of a fixed-effects model", {

  fit   <- fit_certification_model(gaussian_likelihood, ~ x + f, prior_list = list(sigma = unit_sd_prior))
  exact <- log_integrate_half_normal(function(sigma) {
    gaussian_log_density(certification_data$y, fixed_covariance + sigma^2 * diag(nrow(certification_data)))
  })

  expect_bridge_exact(bridge_logml(fit, log_posterior_sigma, prior_list = list(sigma = unit_sd_prior)), exact)
})

test_that("bridge sampling recovers the exact log marginal likelihood of a random-intercept model", {

  random <- list(mu = prior_random(g = random_block(sd = unit_sd_prior)))
  fit    <- fit_certification_model(gaussian_likelihood, ~ x + f + (1 | g),
                                    prior_list = list(sigma = unit_sd_prior), random = random)
  exact  <- log_integrate_half_normal(function(tau) {
    log_integrate_half_normal(function(sigma) {
      gaussian_log_density(
        certification_data$y,
        fixed_covariance + tau^2 * certification_same_group + sigma^2 * diag(nrow(certification_data))
      )
    })
  })

  # the latent intercepts as bridge coordinates
  expect_bridge_exact(bridge_logml(fit, log_posterior_sigma, prior_list = list(sigma = unit_sd_prior)), exact)

  # the intercepts integrated out of the bridge target
  log_posterior_marginal <- function(parameters, data, bridge_context) {
    covariance <- bridge_context$marginalized_random$mu$covariance + parameters$sigma^2 * diag(length(data$y))
    gaussian_log_density(data$y - parameters$mu, covariance)
  }
  expect_bridge_exact(
    bridge_logml(fit, log_posterior_marginal, prior_list = list(sigma = unit_sd_prior),
                 bridge_context = "marginal", formula_random_effects_marginalize_list = list(mu = "g")),
    exact
  )
})
