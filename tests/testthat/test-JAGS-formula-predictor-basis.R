skip_if_not_test_profile("unit")

.formula_predictor_basis_fit <- function(formula_result, samples,
                                         extra_priors = list()) {

  samples <- coda::mcmc(as.matrix(samples))
  class(samples) <- c("BayesTools_fit", class(samples))
  parameter <- formula_result$formula_design$parameter
  attr(samples, "prior_list") <- c(
    formula_result$prior_list,
    extra_priors
  )
  attr(samples, "formula_design") <- stats::setNames(
    list(formula_result$formula_design),
    parameter
  )
  if (!is.null(formula_result$formula_scale)) {
    attr(samples, "formula_scale") <- stats::setNames(
      list(formula_result$formula_scale),
      parameter
    )
  }
  attach_test_parameter_map(samples)
}

.formula_predictor_basis_multi_fit <- function(formula_results, samples,
                                               extra_priors = list()) {

  samples <- coda::mcmc(as.matrix(samples))
  class(samples) <- c("BayesTools_fit", class(samples))
  attr(samples, "prior_list") <- c(
    do.call(c, unname(lapply(formula_results, `[[`, "prior_list"))),
    extra_priors
  )
  attr(samples, "formula_design") <- lapply(
    formula_results,
    `[[`,
    "formula_design"
  )
  attach_test_parameter_map(samples)
}


test_that("formula predictor bases preserve exact factor columns and ordering", {

  data <- data.frame(
    group = factor(
      c("sensitivity", "specificity", "sensitivity", "specificity"),
      levels = c("sensitivity", "specificity")
    )
  )
  formula_result <- JAGS_formula(
    formula = ~ 0 + group,
    parameter = "mu",
    data = data,
    prior_list = list(
      group = prior_factor(
        "normal",
        list(mean = 0, sd = 1),
        contrast = "independent"
      )
    )
  )
  coordinate_names <- c("mu_group[1]", "mu_group[2]")
  samples <- matrix(
    c(-0.2, 0.5, 0.1, 0.8),
    nrow = 2L,
    byrow = TRUE,
    dimnames = list(NULL, coordinate_names)
  )
  fit <- .formula_predictor_basis_fit(formula_result, samples)
  directions <- matrix(
    c(1, 0, 0.5, -1),
    nrow = 2L,
    byrow = TRUE,
    dimnames = list(NULL, coordinate_names)
  )

  basis <- JAGS_formula_predictor_basis(
    fit,
    directions = directions,
    posterior_samples = samples
  )
  expected <- directions %*% t(
    formula_result$formula_design$model_matrix[, c("group1", "group2")]
  )

  expect_s3_class(basis, "BayesTools_formula_predictor_basis")
  expect_identical(basis$status, "affine")
  expect_identical(basis$parameter, "mu")
  expect_identical(basis$coordinates, coordinate_names)
  expect_identical(unname(basis$basis), unname(expected))
})


test_that("formula predictor bases apply persisted scaling and multipliers", {

  x_prior <- prior("normal", list(mean = 0, sd = 1))
  attr(x_prior, "multiply_by") <- "sigma"
  formula_result <- JAGS_formula(
    formula = ~ 1 + x,
    parameter = "mu",
    data = data.frame(x = c(2, 4, 6)),
    prior_list = list(
      intercept = prior("normal", list(mean = 0, sd = 1)),
      x = x_prior
    ),
    formula_scale = TRUE
  )
  sample_names <- c("mu_intercept", "mu_x", "sigma")
  samples <- matrix(
    c(0.1, -0.2, 2, -0.1, 0.3, 4),
    nrow = 2L,
    byrow = TRUE,
    dimnames = list(NULL, sample_names)
  )
  fit <- .formula_predictor_basis_fit(
    formula_result,
    samples,
    extra_priors = list(sigma = prior("gamma", list(shape = 2, rate = 1)))
  )
  directions <- matrix(
    c(1, 0.5),
    ncol = 1L,
    dimnames = list(NULL, "mu_x")
  )

  basis <- JAGS_formula_predictor_basis(
    fit,
    directions = directions,
    posterior_samples = samples
  )
  x_column <- formula_result$formula_design$model_matrix[, "x"]
  expected <- outer(directions[, 1L] * samples[, "sigma"], x_column)

  expect_identical(basis$status, "affine")
  expect_equal(basis$basis, expected, tolerance = 0)
})


test_that("formula predictor bases expose metadata-declared fallbacks", {

  log_formula <- ~ 1
  attr(log_formula, "log(intercept)") <- TRUE
  log_result <- JAGS_formula(
    formula = log_formula,
    parameter = "mu",
    data = data.frame(row = 1:2),
    prior_list = list(
      intercept = prior("gamma", list(shape = 2, rate = 1))
    )
  )
  log_samples <- matrix(
    c(0.5, 1.5),
    ncol = 1L,
    dimnames = list(NULL, "mu_intercept")
  )
  log_fit <- .formula_predictor_basis_fit(log_result, log_samples)
  log_basis <- JAGS_formula_predictor_basis(
    log_fit,
    directions = c(mu_intercept = 1),
    posterior_samples = log_samples
  )

  expect_identical(log_basis$status, "non_affine")
  expect_match(log_basis$reason, "log\\(\\)")
  expect_null(log_basis$basis)

  z_prior <- prior("normal", list(mean = 0, sd = 1))
  attr(z_prior, "multiply_by") <- "mu_x"
  multiplier_result <- JAGS_formula(
    formula = ~ 0 + x + z,
    parameter = "mu",
    data = data.frame(x = c(-1, 1), z = c(1, -1)),
    prior_list = list(
      x = prior("normal", list(mean = 0, sd = 1)),
      z = z_prior
    )
  )
  multiplier_samples <- matrix(
    c(0.2, 0.3, -0.1, 0.4),
    nrow = 2L,
    byrow = TRUE,
    dimnames = list(NULL, c("mu_x", "mu_z"))
  )
  multiplier_fit <- .formula_predictor_basis_fit(
    multiplier_result,
    multiplier_samples
  )
  multiplier_basis <- JAGS_formula_predictor_basis(
    multiplier_fit,
    directions = c(mu_x = 1),
    posterior_samples = multiplier_samples[1L, , drop = FALSE]
  )

  expect_identical(multiplier_basis$status, "unsupported")
  expect_match(multiplier_basis$reason, "multipliers")
  expect_null(multiplier_basis$basis)
})


test_that("formula dependencies cover every persisted predictor", {

  mu_result <- JAGS_formula(
    formula = ~ 0 + x,
    parameter = "mu",
    data = data.frame(x = c(-1, 1)),
    prior_list = list(x = prior("normal", list(mean = 0, sd = 1)))
  )
  z_prior <- prior("normal", list(mean = 0, sd = 1))
  attr(z_prior, "multiply_by") <- "mu_x"
  scale_result <- JAGS_formula(
    formula = ~ 0 + z,
    parameter = "log_tau",
    data = data.frame(z = c(1, 2)),
    prior_list = list(z = z_prior)
  )
  samples <- matrix(
    c(0.2, -0.1, 0.4, 0.3),
    nrow = 2L,
    byrow = TRUE,
    dimnames = list(NULL, c("mu_x", "log_tau_z"))
  )
  fit <- .formula_predictor_basis_multi_fit(
    list(mu = mu_result, log_tau = scale_result),
    samples
  )

  dependencies <- JAGS_formula_coordinate_dependencies(fit, "mu_x")
  expect_equal(
    dependencies,
    data.frame(
      coordinate_name = c("mu_x", "mu_x"),
      formula_parameter = c("mu", "log_tau"),
      dependency_type = c("coefficient", "multiplier"),
      stringsAsFactors = FALSE
    )
  )

  basis <- JAGS_formula_predictor_basis(
    fit,
    directions = c(mu_x = 1),
    posterior_samples = samples[1L, , drop = FALSE]
  )
  expect_identical(basis$status, "unsupported")
  expect_match(basis$reason, "multipliers")
  expect_null(basis$basis)
})


test_that("formula dependencies include external random SD sources", {

  result <- JAGS_formula(
    formula = ~ 1 + diag(1 | study),
    parameter = "mu",
    data = data.frame(study = factor(c("a", "a", "b"))),
    prior_list = list(
      intercept = prior("normal", list(mean = 0, sd = 1))
    ),
    prior_random = prior_random(
      study = random_block(sd_source = random_sd_source("tau"))
    )
  )
  samples <- matrix(
    c(0.1, 0.5, -0.2, 0.8),
    nrow = 2L,
    byrow = TRUE,
    dimnames = list(NULL, c("mu_intercept", "tau"))
  )
  fit <- .formula_predictor_basis_fit(
    result,
    samples,
    extra_priors = list(tau = prior("gamma", list(shape = 2, rate = 1)))
  )

  dependencies <- JAGS_formula_coordinate_dependencies(fit, "tau")
  expect_equal(dependencies$coordinate_name, "tau")
  expect_equal(dependencies$formula_parameter, "mu")
  expect_equal(dependencies$dependency_type, "random_sd_source")
})
