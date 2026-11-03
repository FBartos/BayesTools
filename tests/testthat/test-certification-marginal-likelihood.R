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
# independent grid sum (step 0.002) to 1e-5. A model with 50 random intercepts
# (54 bridge coordinates) is held to the same tolerance under 3 bridge seeds;
# in development the errors did not grow with the number of groups (6, 20,
# and 50 groups: at most 0.014 with the latent intercepts and 0.007 with them
# integrated out, across 3 bridge seeds).

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
# the bridge target with the random intercepts integrated out
log_posterior_marginal <- function(parameters, data, bridge_context) {
  covariance <- bridge_context$marginalized_random$mu$covariance + parameters$sigma^2 * diag(length(data$y))
  gaussian_log_density(data$y - parameters$mu, covariance)
}

# the bridge estimates of a fit under several seeds
bridge_logml <- function(fit, log_posterior, prior_list = NULL, seeds = 1:10,
                         data = certification_jags_data, ...) {
  vapply(seeds, function(seed) {
    JAGS_bridgesampling(fit, log_posterior = log_posterior, data = data,
                        prior_list = prior_list, seed = seed, ...)$logml
  }, numeric(1))
}

expect_bridge_exact <- function(logml, exact) {
  expect_lt(stats::sd(logml), 0.01)
  expect_lt(max(abs(logml - exact)), certification_tolerance)
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
  expect_bridge_exact(
    bridge_logml(fit, log_posterior_marginal, prior_list = list(sigma = unit_sd_prior),
                 bridge_context = "marginal", formula_random_effects_marginalize_list = list(mu = "g")),
    exact
  )
})

test_that("bridge sampling stays accurate with 50 random intercepts", {

  set.seed(50)
  data <- data.frame(x = stats::rnorm(250L), g = factor(rep(seq_len(50L), each = 5L)))
  data$y <- 0.3 + 0.5 * data$x + stats::rnorm(50L, 0, 0.5)[as.integer(data$g)] + stats::rnorm(250L, 0, 0.7)
  jags_data <- list(y = data$y, N = nrow(data))

  # for each tau, X X' + tau^2 Z Z' + sigma^2 I has the eigenvectors of
  # X X' + tau^2 Z Z' and its eigenvalues plus sigma^2
  same_group <- outer(data$g, data$g, "==") * 1
  exact <- log_integrate_half_normal(function(tau) {
    decomposition <- eigen(tcrossprod(cbind(1, data$x)) + tau^2 * same_group, symmetric = TRUE)
    squares <- drop(crossprod(decomposition$vectors, data$y))^2
    log_integrate_half_normal(function(sigma) {
      variances <- decomposition$values + sigma^2
      -0.5 * (length(data$y) * log(2 * pi) + sum(log(variances)) + sum(squares / variances))
    })
  })

  fit <- JAGS_fit(
    gaussian_likelihood, jags_data, prior_list = list(sigma = unit_sd_prior),
    formula_list = list(mu = ~ x + (1 | g)), formula_data_list = list(mu = data),
    formula_prior_list = list(mu = list(intercept = unit_prior, x = unit_prior)),
    formula_random_prior_list = list(mu = prior_random(g = random_block(sd = unit_sd_prior))),
    chains = 4, adapt = 1000, burnin = 1000, sample = 5000, seed = 1
  )

  latent   <- bridge_logml(fit, log_posterior_sigma, prior_list = list(sigma = unit_sd_prior),
                           seeds = 1:3, data = jags_data)
  marginal <- bridge_logml(fit, log_posterior_marginal, prior_list = list(sigma = unit_sd_prior),
                           seeds = 1L, data = jags_data, bridge_context = "marginal",
                           formula_random_effects_marginalize_list = list(mu = "g"))
  expect_lt(max(abs(latent - exact)), certification_tolerance)
  expect_lt(abs(marginal - exact), certification_tolerance)
})

test_that("bridge sampling recovers the exact log marginal likelihood of mean-centered blocks", {

  # The mean-centered parameterization translates the group effects around the
  # fixed terms they absorb; the bridge target (standardized latents) and the
  # exact value are those of the noncentered model. Mean centering suits data
  # that inform the group effects (with weakly informed effects it mixes worse
  # than noncentering), so these models have four replicates per cell and a
  # residual SD of 0.25. Each model is integrated over its two scale
  # parameters. The estimates of 4 x 10000 draws are held to 0.05 of the
  # exact values and to a Monte Carlo SD below 0.025 across 10 bridge seeds:
  # with these 25 to 30 bridge coordinates the bridge estimates spread more
  # than those of the smaller models above. In development, across three fit
  # seeds, the largest errors were 0.013, 0.027, and 0.042 and the SDs at most
  # 0.005, 0.011, and 0.018 (noncentered fits of the same models, 4 x 5000
  # draws: errors up to 0.042, 0.070, and 0.153, SDs up to 0.017, 0.023, and
  # 0.046).
  data <- local({
    set.seed(7)
    data <- expand.grid(rep = 1:4, f = factor(c("a", "b", "c")), g = factor(1:6))
    data$x <- stats::rnorm(nrow(data))
    data$w <- c(-0.8, -0.3, 0, 0.2, 0.5, 0.9)[as.integer(data$g)]
    cell <- (as.integer(data$g) - 1L) * 3L + as.integer(data$f)
    data$y <- 0.3 + 0.5 * data$x + c(-0.2, 0, 0.3)[as.integer(data$f)] + 0.4 * data$w +
      stats::rnorm(6L, 0, 0.4)[as.integer(data$g)] + stats::rnorm(18L, 0, 0.3)[cell] +
      stats::rnorm(nrow(data), 0, 0.25)
    data
  })
  n <- nrow(data)
  jags_data  <- list(y = data$y, N = n)
  same_group <- outer(data$g, data$g, "==") * 1
  same_cell  <- same_group * outer(data$f, data$f, "==")
  design     <- cbind(1, data$x, contr.meandif(levels(data$f))[as.integer(data$f), ])
  centered_block <- function(...) random_block(sd = unit_sd_prior, parameterization = "mean_centered", ...)
  fit_mean_centered <- function(likelihood, formula, random, prior_list = list(sigma = unit_sd_prior),
                                formula_prior_list = fixed_priors) {
    JAGS_fit(
      likelihood, jags_data, prior_list = prior_list,
      formula_list = list(mu = formula), formula_data_list = list(mu = data),
      formula_prior_list = formula_prior_list, formula_random_prior_list = list(mu = random),
      chains = 4, adapt = 1000, burnin = 1000, sample = 10000, seed = 1
    )
  }
  expect_bridge_close <- function(logml, exact) {
    expect_lt(stats::sd(logml), 0.025)
    expect_lt(max(abs(logml - exact)), certification_tolerance)
  }
  absorbed <- function(fit) {
    lapply(JAGS_formula_design(fit)$mu$random_effects, function(block) block$mean_translation$absorbed_terms)
  }

  # one SD for the level effects of f within g: the block spans the intercept
  # and the mean-difference columns of f
  fit <- fit_mean_centered(gaussian_likelihood, ~ x + f + id(0 + f | g),
                           prior_random(g = centered_block(contrasts = c(f = "independent"))))
  expect_identical(absorbed(fit), list(c("intercept", "f")))
  exact <- log_integrate_half_normal(function(tau) {
    log_integrate_half_normal(function(sigma) {
      gaussian_log_density(data$y, tcrossprod(design) + tau^2 * same_cell + sigma^2 * diag(n))
    })
  })
  expect_bridge_close(bridge_logml(fit, log_posterior_sigma, prior_list = list(sigma = unit_sd_prior),
                                   data = jags_data), exact)

  # a group-level moderator enters the location of the random intercept
  fit <- fit_mean_centered(gaussian_likelihood, ~ x + w + (1 | g), prior_random(g = centered_block()),
                           formula_prior_list = list(mu = list(intercept = unit_prior, x = unit_prior, w = unit_prior)))
  expect_identical(absorbed(fit), list(c("intercept", "w")))
  exact <- log_integrate_half_normal(function(tau) {
    log_integrate_half_normal(function(sigma) {
      gaussian_log_density(data$y, tcrossprod(cbind(1, data$x, data$w)) + tau^2 * same_group + sigma^2 * diag(n))
    })
  })
  expect_bridge_close(bridge_logml(fit, log_posterior_sigma, prior_list = list(sigma = unit_sd_prior),
                                   data = jags_data), exact)

  # nested blocks with the residual SD known: the cells of g and f take f and
  # the locations of g, which takes the intercept
  known_sigma <- "model{
    for(i in 1:N){
      y[i] ~ dnorm(mu[i], pow(0.25, -2))
    }
  }
  "
  fit <- fit_mean_centered(known_sigma, ~ x + f + (1 | g/f),
                           prior_random(sd = unit_sd_prior, g = centered_block(), f_g = centered_block()),
                           prior_list = NULL)
  expect_identical(absorbed(fit), list("f", "intercept"))
  exact <- log_integrate_half_normal(function(tau_g) {
    log_integrate_half_normal(function(tau_cell) {
      gaussian_log_density(data$y, tcrossprod(design) + tau_g^2 * same_group + tau_cell^2 * same_cell +
                             0.25^2 * diag(n))
    })
  })
  log_posterior <- function(parameters, data) sum(stats::dnorm(data$y, parameters$mu, 0.25, log = TRUE))
  expect_bridge_close(bridge_logml(fit, log_posterior, data = jags_data), exact)
})
