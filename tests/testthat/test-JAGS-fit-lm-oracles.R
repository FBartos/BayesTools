skip_if_not_test_profile("fit")

# ============================================================================ #
# TEST FILE: JAGS Fit LM Oracles
# ============================================================================ #
#
# PURPOSE:
#   Fit-profile semantic oracle tests for Gaussian JAGS formula models. These
#   compare automatic formula scaling and manual standardization against lm()
#   predictions on the original data scale, and compare a known-sigma Gaussian
#   model against its closed-form posterior.
#
# TAGS: @fit, @JAGS, @formula, @standardization
# ============================================================================ #

skip_on_cran()
skip_if_not_installed("rjags")
skip_if_not_installed("runjags")

.fit_gaussian_formula_oracle <- function(data, formula_data, formula_scale_list = NULL, seed = 1L) {
  model_syntax <- paste0(
    "model{\n",
    "for(i in 1:N){\n",
    "  y[i] ~ dnorm(mu[i], pow(sigma, -2))\n",
    "}\n",
    "}"
  )

  JAGS_fit(
    model_syntax = model_syntax,
    data = list(y = data$y, N = nrow(data)),
    prior_list = list(
      sigma = prior("lognormal", list(0, 1))
    ),
    formula_list = list(mu = ~ x1 * x2),
    formula_data_list = list(mu = formula_data),
    formula_prior_list = list(
      mu = list(
        "intercept" = prior("normal", list(0, 10)),
        "x1" = prior("normal", list(0, 5)),
        "x2" = prior("normal", list(0, 5)),
        "x1:x2" = prior("normal", list(0, 5))
      )
    ),
    formula_scale_list = formula_scale_list,
    chains = 2,
    adapt = 250,
    burnin = 500,
    sample = 1200,
    seed = seed,
    silent = TRUE
  )
}

.fit_known_sigma_gaussian_formula_oracle <- function(data, sigma_known, seed = 1L) {
  model_syntax <- paste0(
    "model{\n",
    "for(i in 1:N){\n",
    "  y[i] ~ dnorm(mu[i], pow(sigma_known, -2))\n",
    "}\n",
    "}"
  )

  JAGS_fit(
    model_syntax = model_syntax,
    data = list(y = data$y, N = nrow(data), sigma_known = sigma_known),
    prior_list = NULL,
    formula_list = list(mu = ~ x1 * x2),
    formula_data_list = list(mu = data[c("x1", "x2")]),
    formula_prior_list = list(
      mu = list(
        "intercept" = prior("normal", list(0, 10)),
        "x1" = prior("normal", list(0, 5)),
        "x2" = prior("normal", list(0, 5)),
        "x1:x2" = prior("normal", list(0, 5))
      )
    ),
    formula_scale_list = list(mu = list(x1 = TRUE, x2 = TRUE)),
    chains = 2,
    adapt = 500,
    burnin = 1000,
    sample = 3000,
    seed = seed,
    silent = TRUE
  )
}

test_that("Gaussian JAGS formula fit agrees with lm oracle after scaling", {
  data <- bayestools_oracle_gaussian_regression_data()
  manual_formula_data <- bayestools_manual_scaled_data(data[c("x1", "x2")], c("x1", "x2"))
  manual_scale <- attr(manual_formula_data, "manual_scale")

  fit_manual <- .fit_gaussian_formula_oracle(
    data = data,
    formula_data = manual_formula_data,
    formula_scale_list = NULL,
    seed = 7701L
  )
  attr(fit_manual, "manual_scale") <- manual_scale

  fit_auto <- .fit_gaussian_formula_oracle(
    data = data,
    formula_data = data[c("x1", "x2")],
    formula_scale_list = list(mu = list(x1 = TRUE, x2 = TRUE)),
    seed = 7701L
  )

  expect_formula_scale_equal(attr(fit_auto, "formula_scale")$mu, manual_scale)

  scaled_parameters <- c("mu_intercept", "mu_x1", "mu_x2", "mu_x1__xXx__x2")
  posterior_manual <- as.matrix(suppressWarnings(coda::as.mcmc(fit_manual)))
  posterior_auto <- as.matrix(suppressWarnings(coda::as.mcmc(fit_auto)))

  expect_posterior_summary_close(
    posterior_auto[, scaled_parameters, drop = FALSE],
    colMeans(posterior_manual[, scaled_parameters, drop = FALSE]),
    tolerance = 0.05
  )

  posterior_original <- transform_scale_samples(
    posterior_auto[, scaled_parameters, drop = FALSE],
    attr(fit_auto, "formula_scale")
  )
  posterior_manual_original <- transform_scale_samples(
    posterior_manual[, scaled_parameters, drop = FALSE],
    list(mu = manual_scale)
  )

  lm_fit <- stats::lm(y ~ x1 * x2, data = data)
  lm_coef <- stats::coef(lm_fit)
  names(lm_coef) <- bayestools_lm_coef_to_jags_names(names(lm_coef))

  expect_posterior_summary_close(
    posterior_original[, names(lm_coef), drop = FALSE],
    lm_coef,
    tolerance = 0.30
  )

  newdata <- data.frame(
    x1 = stats::quantile(data$x1, probs = c(0.10, 0.50, 0.90), names = FALSE),
    x2 = stats::quantile(data$x2, probs = c(0.20, 0.60, 0.80), names = FALSE)
  )
  prediction_design <- stats::model.matrix(~ x1 * x2, data = newdata)
  prediction_names <- bayestools_lm_coef_to_jags_names(colnames(prediction_design))
  prediction_auto <- as.numeric(prediction_design %*% colMeans(
    posterior_original[, prediction_names, drop = FALSE]
  ))
  prediction_manual <- as.numeric(prediction_design %*% colMeans(
    posterior_manual_original[, prediction_names, drop = FALSE]
  ))
  expect_equal(prediction_auto, prediction_manual, tolerance = 0.15)

  expect_lm_predictions_equal(
    lm_fit = lm_fit,
    posterior = posterior_original,
    newdata = newdata,
    formula = ~ x1 * x2,
    tolerance = 0.30
  )

  sigma_samples <- posterior_auto[, "sigma"]
  expect_equal(mean(sigma_samples), summary(lm_fit)$sigma, tolerance = 0.25)

  posterior_intervals <- apply(
    posterior_original[, names(lm_coef), drop = FALSE],
    2,
    stats::quantile,
    probs = c(0.005, 0.995)
  )
  expect_true(all(lm_coef >= posterior_intervals[1, ] & lm_coef <= posterior_intervals[2, ]))
})

test_that("known-sigma Gaussian JAGS formula fit matches closed-form posterior", {
  data <- bayestools_oracle_gaussian_regression_data()
  sigma_known <- 0.6
  parameters <- c("mu_intercept", "mu_x1", "mu_x2", "mu_x1__xXx__x2")

  fit <- .fit_known_sigma_gaussian_formula_oracle(
    data = data,
    sigma_known = sigma_known,
    seed = 7711L
  )

  posterior <- as.matrix(suppressWarnings(coda::as.mcmc(fit)))
  posterior <- posterior[, parameters, drop = FALSE]

  scale_info <- attr(fit, "formula_scale")$mu
  expected_scale <- bayestools_manual_scaled_data(data[c("x1", "x2")], c("x1", "x2"))
  expect_formula_scale_equal(scale_info, attr(expected_scale, "manual_scale"))

  design <- stats::model.matrix(~ x1 * x2, data = expected_scale)
  jags_design <- JAGS_formula_design(fit, "mu")$model_matrix
  expect_equal(colnames(design), c("(Intercept)", "x1", "x2", "x1:x2"))
  expect_equal(colnames(jags_design), c("(Intercept)", "x1", "x2", "x1__xXx__x2"))
  expect_equal(unname(jags_design), unname(design), tolerance = 1e-12)
  expect_setequal(colnames(posterior), parameters)

  oracle <- bayestools_gaussian_posterior_oracle(
    X = design,
    y = data$y,
    sigma = sigma_known,
    prior_mean = c(0, 0, 0, 0),
    prior_sd = c(10, 5, 5, 5)
  )
  names(oracle$mean) <- parameters
  names(oracle$sd) <- parameters

  expect_equal(colMeans(posterior), oracle$mean, tolerance = 0.07)
  expect_equal(apply(posterior, 2, stats::sd), oracle$sd, tolerance = 0.04)
  expect_equal(unname(stats::cov(posterior)), unname(oracle$cov), tolerance = 0.02)

  posterior_intervals <- apply(posterior, 2, stats::quantile, probs = c(0.01, 0.99))
  expect_true(all(oracle$mean >= posterior_intervals[1, ] & oracle$mean <= posterior_intervals[2, ]))
})
