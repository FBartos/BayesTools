skip_if_not_test_profile("fixture")

# ============================================================================ #
# TEST FILE: Model Averaging Functions
# ============================================================================ #
#
# PURPOSE:
#   Tests for compute_inference, ensemble_inference, mix_posteriors, models_inference,
#   inclusion_BF, weightfunctions_mapping, and related Bayesian model averaging
#   functions in R/model-averaging.R
#
# DEPENDENCIES:
#   - bridgesampling: Required for marginal likelihood computation
#   - rjags: For tests using pre-fitted models
#   - common-functions.R: Test helpers
#
# SKIP CONDITIONS:
#   - skip_if_not_installed("bridgesampling")
#   - skip_if_not_installed("rjags")
#   - skip_if_no_fits()
#
# MODELS/FIXTURES:
#   - Uses pre-computed marginal likelihoods and pre-fitted models
#
# TAGS: @evaluation, @model-averaging
# ============================================================================ #

source(testthat::test_path("common-functions.R"))


# Deterministic model-averaging algebra and validation tests live in the
# unit-profile semantic suite: test-model-averaging-edge-cases.R.

.current_model_probabilities <- function(models) {

  prior_weights <- vapply(models, function(model) model$prior_weights, numeric(1))
  prior_probs <- prior_weights / sum(prior_weights)
  log_margliks <- vapply(models, function(model) model$marglik$logml, numeric(1))
  log_post_weights <- log(prior_probs) + log_margliks
  log_post_weights <- log_post_weights - max(log_post_weights)
  post_probs <- exp(log_post_weights)
  post_probs <- post_probs / sum(post_probs)

  list(
    prior_probs = prior_probs,
    log_margliks = log_margliks,
    post_probs = post_probs
  )
}

.current_inclusion_bf <- function(prior_probs, post_probs, is_null) {

  prior_null <- sum(prior_probs[is_null])
  prior_alt <- sum(prior_probs[!is_null])
  post_null <- sum(post_probs[is_null])
  post_alt <- sum(post_probs[!is_null])

  if (prior_null == 0 || prior_alt == 0 || post_null == 0 || post_alt == 0) {
    return(NA_real_)
  }

  (post_alt / post_null) / (prior_alt / prior_null)
}

.current_fit_samples <- function(fit) {

  samples <- BayesTools:::.extract_posterior_samples(
    fit,
    as_list = FALSE
  )
  if (!is.matrix(samples)) {
    samples <- matrix(samples, ncol = 1)
    colnames(samples) <- fit$monitor
  }

  as.matrix(samples)
}

# ============================================================================ #
# SECTION 1: mix_posteriors tests
# ============================================================================ #
test_that("mix_posteriors handles various prior types correctly", {

  skip_on_cran()
  skip_if_not_installed("rjags")
  skip_if_no_fits()

  # Load fits with margliks
  fit_simple_normal <- readRDS(file.path(temp_fits_dir, "fit_simple_normal.RDS"))
  marglik_simple_normal <- readRDS(file.path(temp_marglik_dir, "fit_simple_normal.RDS"))

  fit_simple_spike <- readRDS(file.path(temp_fits_dir, "fit_simple_spike.RDS"))
  marglik_simple_spike <- readRDS(file.path(temp_marglik_dir, "fit_simple_spike.RDS"))

  # Create model list for simple priors
  models_simple <- list(
    list(fit = fit_simple_normal, marglik = marglik_simple_normal, prior_weights = 1),
    list(fit = fit_simple_spike, marglik = marglik_simple_spike, prior_weights = 1)
  )

  # Test mix_posteriors with simple priors
  mixed <- mix_posteriors(
    model_list = models_simple,
    parameters = c("m", "s"),
    is_null_list = list("m" = c(FALSE, TRUE), "s" = c(FALSE, FALSE)),
    seed = 1,
    n_samples = 1000
  )

  expect_true(inherits(mixed, "mixed_posteriors"))
  expect_named(mixed, c("m", "s"))
  expect_s3_class(mixed$m, "mixed_posteriors.simple")
  expect_s3_class(mixed$s, "mixed_posteriors.simple")
  expect_equal(length(mixed$m), 1000)
  expect_equal(length(mixed$s), 1000)
  expect_equal(attr(mixed$m, "models_ind"), attr(mixed$s, "models_ind"))
  expect_equal(attr(mixed$m, "sample_ind"), attr(mixed$s, "sample_ind"))

  probabilities_simple <- .current_model_probabilities(models_simple)
  expect_equal(
    vapply(attr(mixed$m, "prior_list"), function(x) x$prior_weights, numeric(1)),
    probabilities_simple$prior_probs,
    tolerance = 1e-12
  )
  expect_equal(
    sum(tabulate(attr(mixed$m, "models_ind"), nbins = length(models_simple))),
    1000
  )

  # Test with conditional = TRUE
  mixed_conditional <- mix_posteriors(
    model_list = models_simple,
    parameters = c("m"),
    is_null_list = list("m" = c(FALSE, TRUE)),
    conditional = TRUE,
    seed = 1,
    n_samples = 1000
  )

  expect_true(inherits(mixed_conditional, "mixed_posteriors"))
  expect_named(mixed_conditional, "m")
  expect_equal(length(mixed_conditional$m), 1000)
  expect_true(all(attr(mixed_conditional$m, "models_ind") == 1))
  expect_equal(
    vapply(attr(mixed_conditional$m, "prior_list"), function(x) x$prior_weights, numeric(1)),
    c(1, 0),
    tolerance = 1e-12
  )
})


test_that("mix_posteriors handles weightfunction priors", {

  skip_on_cran()
  skip_if_not_installed("rjags")
  skip_if_no_fits()

  # Load summary models which have weightfunction priors
  fit_summary0 <- readRDS(file.path(temp_fits_dir, "fit_summary0.RDS"))
  marglik_summary0 <- readRDS(file.path(temp_marglik_dir, "fit_summary0.RDS"))

  fit_summary1 <- readRDS(file.path(temp_fits_dir, "fit_summary1.RDS"))
  marglik_summary1 <- readRDS(file.path(temp_marglik_dir, "fit_summary1.RDS"))

  fit_summary2 <- readRDS(file.path(temp_fits_dir, "fit_summary2.RDS"))
  marglik_summary2 <- readRDS(file.path(temp_marglik_dir, "fit_summary2.RDS"))

  models_wf <- list(
    list(fit = fit_summary0, marglik = marglik_summary0, prior_weights = 1),
    list(fit = fit_summary1, marglik = marglik_summary1, prior_weights = 1),
    list(fit = fit_summary2, marglik = marglik_summary2, prior_weights = 1)
  )

  mixed_wf <- mix_posteriors(
    model_list = models_wf,
    parameters = c("m", "omega"),
    is_null_list = list("m" = c(FALSE, FALSE, FALSE), "omega" = c(TRUE, FALSE, FALSE)),
    seed = 1,
    n_samples = 1000
  )

  expect_true(inherits(mixed_wf, "mixed_posteriors"))
  expect_named(mixed_wf, c("m", "omega"))

  omega_samples <- mixed_wf$omega
  omega_priors <- lapply(models_wf, function(model) {
    prior <- attr(model$fit, "prior_list")[["omega"]]
    if (is.null(prior)) prior_none() else prior
  })
  omega_mapping <- weightfunctions_mapping(omega_priors)
  omega_cuts <- weightfunctions_mapping(omega_priors, cuts_only = TRUE)
  omega_names <- vapply(
    seq_len(length(omega_cuts) - 1L),
    function(i) paste0("omega[", omega_cuts[i], ",", omega_cuts[i + 1L], "]"),
    character(1)
  )

  expect_s3_class(omega_samples, "mixed_posteriors.weightfunction")
  expect_true(is.matrix(omega_samples))
  expect_identical(dim(omega_samples), c(1000L, length(omega_names)))
  expect_identical(colnames(omega_samples), omega_names)

  models_ind <- attr(omega_samples, "models_ind")
  sample_ind <- attr(omega_samples, "sample_ind")
  expect_length(models_ind, nrow(omega_samples))
  expect_length(sample_ind, nrow(omega_samples))
  expect_true(all(models_ind %in% seq_along(models_wf)))

  for (model_i in seq_along(models_wf)) {
    rows <- models_ind == model_i
    if (!any(rows)) {
      next
    }

    fit_samples <- .current_fit_samples(models_wf[[model_i]]$fit)
    expect_true(all(sample_ind[rows] %in% seq_len(nrow(fit_samples))))

    if (is.prior.weightfunction(omega_priors[[model_i]])) {
      source_names <- paste0("omega[", omega_mapping[[model_i]], "]")
      expected <- fit_samples[sample_ind[rows], source_names, drop = FALSE]
    } else {
      expected <- matrix(1, nrow = sum(rows), ncol = length(omega_names))
    }

    expect_equal(
      as.numeric(omega_samples[rows, , drop = FALSE]),
      as.numeric(expected)
    )
  }
})


test_that("mix_posteriors handles factor priors", {

  skip_on_cran()
  skip_if_not_installed("rjags")
  skip_if_no_fits()

  # Load the orthonormal factor models (have both factor priors and marginal likelihoods)
  fit_orthonormal_0 <- readRDS(file.path(temp_fits_dir, "fit_orthonormal_0.RDS"))
  marglik_orthonormal_0 <- readRDS(file.path(temp_marglik_dir, "fit_orthonormal_0.RDS"))

  fit_orthonormal_1 <- readRDS(file.path(temp_fits_dir, "fit_orthonormal_1.RDS"))
  marglik_orthonormal_1 <- readRDS(file.path(temp_marglik_dir, "fit_orthonormal_1.RDS"))

  # Create model list with two different models
  models_factor <- list(
    list(fit = fit_orthonormal_0, marglik = marglik_orthonormal_0, prior_weights = 1),
    list(fit = fit_orthonormal_1, marglik = marglik_orthonormal_1, prior_weights = 1)
  )

  # Get the parameters from the model
  prior_list <- attr(fit_orthonormal_1, "prior_list")
  param_names <- names(prior_list)

  # Filter to factor parameters only
  factor_params <- param_names[sapply(prior_list, is.prior.factor)]

  mixed_factor <- mix_posteriors(
    model_list = models_factor,
    parameters = factor_params[1],  # Just test one
    is_null_list = setNames(list(c(TRUE, FALSE)), factor_params[1]),
    seed = 1,
    n_samples = 1000
  )

  expect_true(inherits(mixed_factor, "mixed_posteriors"))
  expect_named(mixed_factor, factor_params[1])

  factor_parameter <- factor_params[1]
  factor_samples <- mixed_factor[[factor_parameter]]
  factor_levels <- as.integer(BayesTools:::.get_prior_factor_levels(
    prior_list[[factor_parameter]]
  ))
  factor_names <- paste0(factor_parameter, "[", seq_len(factor_levels), "]")

  expect_s3_class(factor_samples, "mixed_posteriors.factor")
  expect_s3_class(factor_samples, "mixed_posteriors.vector")
  expect_true(is.matrix(factor_samples))
  expect_identical(dim(factor_samples), c(1000L, factor_levels))
  expect_identical(colnames(factor_samples), factor_names)

  models_ind <- attr(factor_samples, "models_ind")
  sample_ind <- attr(factor_samples, "sample_ind")
  effective_priors <- attr(factor_samples, "prior_list")
  expect_length(models_ind, nrow(factor_samples))
  expect_length(sample_ind, nrow(factor_samples))
  expect_true(all(models_ind %in% seq_along(models_factor)))

  for (model_i in seq_along(models_factor)) {
    rows <- models_ind == model_i
    if (!any(rows)) {
      next
    }

    fit_samples <- .current_fit_samples(models_factor[[model_i]]$fit)
    expect_true(all(sample_ind[rows] %in% seq_len(nrow(fit_samples))))

    if (is.prior.point(effective_priors[[model_i]])) {
      location <- effective_priors[[model_i]]$parameters[["location"]]
      expected <- matrix(
        rep(location, times = sum(rows)),
        nrow = sum(rows),
        ncol = factor_levels,
        byrow = TRUE
      )
    } else {
      expected <- fit_samples[sample_ind[rows], factor_names, drop = FALSE]
    }

    expect_equal(
      as.numeric(factor_samples[rows, , drop = FALSE]),
      as.numeric(expected)
    )
  }
})


test_that("mix_posteriors handles vector priors", {

  skip_on_cran()
  skip_if_not_installed("rjags")
  skip_if_no_fits()

  # Load vector prior models
  fit_vector_mnormal <- readRDS(file.path(temp_fits_dir, "fit_vector_mnormal.RDS"))

  # Create a mock marglik for testing (we only need the structure)
  mock_marglik <- bridgesampling_object(-100)

  models_vector <- list(
    list(fit = fit_vector_mnormal, marglik = mock_marglik, prior_weights = 1),
    list(fit = fit_vector_mnormal, marglik = mock_marglik, prior_weights = 1)
  )

  prior_list <- attr(fit_vector_mnormal, "prior_list")
  vector_params <- names(prior_list)[sapply(prior_list, is.prior.vector)]

  mixed_vector <- mix_posteriors(
    model_list = models_vector,
    parameters = vector_params[1],
    is_null_list = setNames(list(c(FALSE, FALSE)), vector_params[1]),
    seed = 1,
    n_samples = 1000
  )

  expect_true(inherits(mixed_vector, "mixed_posteriors"))
  expect_named(mixed_vector, vector_params[1])

  vector_parameter <- vector_params[1]
  vector_samples <- mixed_vector[[vector_parameter]]
  vector_length <- as.integer(
    prior_list[[vector_parameter]]$parameters[["K"]]
  )
  vector_names <- paste0(vector_parameter, "[", seq_len(vector_length), "]")

  expect_s3_class(vector_samples, "mixed_posteriors.vector")
  expect_true(is.matrix(vector_samples))
  expect_identical(dim(vector_samples), c(1000L, vector_length))
  expect_identical(colnames(vector_samples), vector_names)

  models_ind <- attr(vector_samples, "models_ind")
  sample_ind <- attr(vector_samples, "sample_ind")
  expect_length(models_ind, nrow(vector_samples))
  expect_length(sample_ind, nrow(vector_samples))
  expect_true(all(models_ind %in% seq_along(models_vector)))

  for (model_i in seq_along(models_vector)) {
    rows <- models_ind == model_i
    if (!any(rows)) {
      next
    }

    fit_samples <- .current_fit_samples(models_vector[[model_i]]$fit)
    expect_true(all(sample_ind[rows] %in% seq_len(nrow(fit_samples))))
    expected <- fit_samples[sample_ind[rows], vector_names, drop = FALSE]

    expect_equal(
      as.numeric(vector_samples[rows, , drop = FALSE]),
      as.numeric(expected)
    )
  }
})


test_that("mix_posteriors respects seed across prior types", {

  skip_on_cran()
  skip_if_not_installed("rjags")
  skip_if_no_fits()

  expect_seed_changes_draws <- function(draw_factory, parameter) {
    same_seed_a <- draw_factory(11L)
    same_seed_b <- draw_factory(11L)
    diff_seed   <- draw_factory(12L)

    expect_equal(
      attr(same_seed_a[[parameter]], "sample_ind"),
      attr(same_seed_b[[parameter]], "sample_ind")
    )
    expect_false(identical(
      attr(same_seed_a[[parameter]], "sample_ind"),
      attr(diff_seed[[parameter]], "sample_ind")
    ))
  }

  fit_simple_normal <- readRDS(file.path(temp_fits_dir, "fit_simple_normal.RDS"))
  marglik_simple_normal <- readRDS(file.path(temp_marglik_dir, "fit_simple_normal.RDS"))
  fit_simple_spike <- readRDS(file.path(temp_fits_dir, "fit_simple_spike.RDS"))
  marglik_simple_spike <- readRDS(file.path(temp_marglik_dir, "fit_simple_spike.RDS"))

  models_simple <- list(
    list(fit = fit_simple_normal, marglik = marglik_simple_normal, prior_weights = 1),
    list(fit = fit_simple_spike, marglik = marglik_simple_spike, prior_weights = 1)
  )

  expect_seed_changes_draws(
    draw_factory = function(seed) {
      mix_posteriors(
        model_list = models_simple,
        parameters = c("m", "s"),
        is_null_list = list("m" = c(FALSE, TRUE), "s" = c(FALSE, FALSE)),
        seed = seed,
        n_samples = 500
      )
    },
    parameter = "m"
  )

  fit_summary0 <- readRDS(file.path(temp_fits_dir, "fit_summary0.RDS"))
  marglik_summary0 <- readRDS(file.path(temp_marglik_dir, "fit_summary0.RDS"))
  fit_summary1 <- readRDS(file.path(temp_fits_dir, "fit_summary1.RDS"))
  marglik_summary1 <- readRDS(file.path(temp_marglik_dir, "fit_summary1.RDS"))
  fit_summary2 <- readRDS(file.path(temp_fits_dir, "fit_summary2.RDS"))
  marglik_summary2 <- readRDS(file.path(temp_marglik_dir, "fit_summary2.RDS"))

  models_wf <- list(
    list(fit = fit_summary0, marglik = marglik_summary0, prior_weights = 1),
    list(fit = fit_summary1, marglik = marglik_summary1, prior_weights = 1),
    list(fit = fit_summary2, marglik = marglik_summary2, prior_weights = 1)
  )

  expect_seed_changes_draws(
    draw_factory = function(seed) {
      mix_posteriors(
        model_list = models_wf,
        parameters = c("m", "omega"),
        is_null_list = list("m" = c(FALSE, FALSE, FALSE), "omega" = c(TRUE, FALSE, FALSE)),
        seed = seed,
        n_samples = 500
      )
    },
    parameter = "omega"
  )

  fit_orthonormal_0 <- readRDS(file.path(temp_fits_dir, "fit_orthonormal_0.RDS"))
  marglik_orthonormal_0 <- readRDS(file.path(temp_marglik_dir, "fit_orthonormal_0.RDS"))
  fit_orthonormal_1 <- readRDS(file.path(temp_fits_dir, "fit_orthonormal_1.RDS"))
  marglik_orthonormal_1 <- readRDS(file.path(temp_marglik_dir, "fit_orthonormal_1.RDS"))

  models_factor <- list(
    list(fit = fit_orthonormal_0, marglik = marglik_orthonormal_0, prior_weights = 1),
    list(fit = fit_orthonormal_1, marglik = marglik_orthonormal_1, prior_weights = 1)
  )

  factor_prior_list <- attr(fit_orthonormal_1, "prior_list")
  factor_parameter <- names(factor_prior_list)[sapply(factor_prior_list, is.prior.factor)][1]

  expect_seed_changes_draws(
    draw_factory = function(seed) {
      mix_posteriors(
        model_list = models_factor,
        parameters = factor_parameter,
        is_null_list = setNames(list(c(TRUE, FALSE)), factor_parameter),
        seed = seed,
        n_samples = 500
      )
    },
    parameter = factor_parameter
  )

  fit_vector_mnormal <- readRDS(file.path(temp_fits_dir, "fit_vector_mnormal.RDS"))
  mock_marglik <- bridgesampling_object(-100)
  models_vector <- list(
    list(fit = fit_vector_mnormal, marglik = mock_marglik, prior_weights = 1),
    list(fit = fit_vector_mnormal, marglik = mock_marglik, prior_weights = 1)
  )

  vector_prior_list <- attr(fit_vector_mnormal, "prior_list")
  vector_parameter <- names(vector_prior_list)[sapply(vector_prior_list, is.prior.vector)][1]

  expect_seed_changes_draws(
    draw_factory = function(seed) {
      mix_posteriors(
        model_list = models_vector,
        parameters = vector_parameter,
        is_null_list = setNames(list(c(FALSE, FALSE)), vector_parameter),
        seed = seed,
        n_samples = 500
      )
    },
    parameter = vector_parameter
  )
})


# ============================================================================ #
# SECTION 6: ensemble_inference tests
# ============================================================================ #
test_that("ensemble_inference handles different configurations", {

  skip_on_cran()
  skip_if_not_installed("rjags")
  skip_if_no_fits()

  # Load fits with margliks
  fit_simple_normal <- readRDS(file.path(temp_fits_dir, "fit_simple_normal.RDS"))
  marglik_simple_normal <- readRDS(file.path(temp_marglik_dir, "fit_simple_normal.RDS"))

  fit_simple_spike <- readRDS(file.path(temp_fits_dir, "fit_simple_spike.RDS"))
  marglik_simple_spike <- readRDS(file.path(temp_marglik_dir, "fit_simple_spike.RDS"))

  models <- list(
    list(fit = fit_simple_normal, marglik = marglik_simple_normal, prior_weights = 1),
    list(fit = fit_simple_spike, marglik = marglik_simple_spike, prior_weights = 1)
  )

  # Test with integer is_null specification
  inference_int <- ensemble_inference(
    model_list = models,
    parameters = "m",
    is_null_list = list("m" = 2)  # Second model is null
  )

  expect_true(inherits(inference_int$m, "inference"))
  current_probabilities <- .current_model_probabilities(models)
  expected_is_null <- c(FALSE, TRUE)
  expect_identical(attr(inference_int$m, "is_null"), expected_is_null)
  expect_equal(
    inference_int$m$prior_probs,
    current_probabilities$prior_probs,
    tolerance = 1e-12
  )
  expect_equal(
    inference_int$m$post_probs,
    current_probabilities$post_probs,
    tolerance = 1e-12
  )
  expect_equal(sum(inference_int$m$prior_probs), 1, tolerance = 1e-12)
  expect_equal(sum(inference_int$m$post_probs), 1, tolerance = 1e-12)
  expect_equal(
    inference_int$m$BF,
    .current_inclusion_bf(
      current_probabilities$prior_probs,
      current_probabilities$post_probs,
      expected_is_null
    ),
    tolerance = 1e-12
  )

  # Test conditional inference
  inference_cond <- ensemble_inference(
    model_list = models,
    parameters = "m",
    is_null_list = list("m" = c(FALSE, TRUE)),
    conditional = TRUE
  )

  expect_true(attr(inference_cond, "conditional"))
  expect_equal(inference_cond$m$prior_probs, c(1, 0), tolerance = 1e-12)
  expect_equal(inference_cond$m$post_probs, c(1, 0), tolerance = 1e-12)
  expect_equal(inference_cond$m$BF, inference_int$m$BF, tolerance = 1e-12)
  expect_equal(
    sum(inference_cond$m$prior_probs),
    1,
    tolerance = 1e-12
  )
  expect_equal(
    sum(inference_cond$m$post_probs),
    1,
    tolerance = 1e-12
  )

})


# ============================================================================ #
# SECTION 7: models_inference tests
# ============================================================================ #
test_that("models_inference computes correctly", {

  skip_on_cran()
  skip_if_not_installed("rjags")
  skip_if_no_fits()

  # Load fits with margliks
  fit_simple_normal <- readRDS(file.path(temp_fits_dir, "fit_simple_normal.RDS"))
  marglik_simple_normal <- readRDS(file.path(temp_marglik_dir, "fit_simple_normal.RDS"))

  fit_simple_spike <- readRDS(file.path(temp_fits_dir, "fit_simple_spike.RDS"))
  marglik_simple_spike <- readRDS(file.path(temp_marglik_dir, "fit_simple_spike.RDS"))

  models <- list(
    list(fit = fit_simple_normal, marglik = marglik_simple_normal, prior_weights = 1),
    list(fit = fit_simple_spike, marglik = marglik_simple_spike, prior_weights = 2)
  )

  models_with_inference <- models_inference(models)

  # Check that inference was added to each model
  expect_true("inference" %in% names(models_with_inference[[1]]))
  expect_true("inference" %in% names(models_with_inference[[2]]))

  current_probabilities <- .current_model_probabilities(models)
  expect_equal(
    vapply(models_with_inference, function(model) model$inference$m_number, numeric(1)),
    seq_along(models)
  )
  expect_equal(
    vapply(models_with_inference, function(model) model$inference$marglik, numeric(1)),
    current_probabilities$log_margliks,
    tolerance = 1e-12
  )
  expect_equal(
    vapply(models_with_inference, function(model) model$inference$prior_prob, numeric(1)),
    current_probabilities$prior_probs,
    tolerance = 1e-12
  )
  expect_equal(
    vapply(models_with_inference, function(model) model$inference$post_prob, numeric(1)),
    current_probabilities$post_probs,
    tolerance = 1e-12
  )

  expected_inclusion_bf <- vapply(seq_along(models), function(model_i) {
    is_null <- rep(TRUE, length(models))
    is_null[model_i] <- FALSE
    .current_inclusion_bf(
      current_probabilities$prior_probs,
      current_probabilities$post_probs,
      is_null
    )
  }, numeric(1))
  expect_equal(
    vapply(models_with_inference, function(model) model$inference$inclusion_BF, numeric(1)),
    expected_inclusion_bf,
    tolerance = 1e-12
  )

  # Check prior probs reflect weights (1:2 ratio)
  expect_equal(models_with_inference[[1]]$inference$prior_prob, 1/3, tolerance = 1e-10)
  expect_equal(models_with_inference[[2]]$inference$prior_prob, 2/3, tolerance = 1e-10)

  # Check posterior probs sum to 1
  total_post_prob <- sum(sapply(models_with_inference, function(m) m$inference$post_prob))
  expect_equal(total_post_prob, 1, tolerance = 1e-10)
  expect_equal(
    sum(sapply(models_with_inference, function(m) m$inference$prior_prob)),
    1,
    tolerance = 1e-10
  )

})


# ============================================================================ #
# SECTION 8: as_mixed_posteriors and as_marginal_inference tests
# ============================================================================ #
test_that("as_mixed_posteriors works correctly with BayesTools_fit objects", {

  skip_on_cran()
  skip_if_not_installed("rjags")
  skip_if_no_fits()

  # Load a fitted model
  fit_simple_normal <- readRDS(file.path(temp_fits_dir, "fit_simple_normal.RDS"))

  # as_mixed_posteriors needs a BayesTools_fit object
  mixed <- as_mixed_posteriors(fit_simple_normal, parameters = c("m", "s"))

  expect_true(inherits(mixed, "mixed_posteriors"))
  expect_named(mixed, c("m", "s"))

  current_samples <- suppressWarnings(coda::as.mcmc(fit_simple_normal))
  if (!is.matrix(current_samples)) {
    current_samples <- matrix(current_samples, ncol = 1)
    colnames(current_samples) <- fit_simple_normal$monitor
  }

  expect_equal(as.numeric(mixed$m), as.numeric(current_samples[, "m"]))
  expect_equal(as.numeric(mixed$s), as.numeric(current_samples[, "s"]))
  expect_equal(attr(mixed$m, "models_ind"), rep(1, nrow(current_samples)))
  expect_equal(attr(mixed$s, "models_ind"), rep(1, nrow(current_samples)))
  expect_identical(attr(mixed, "prior_list"), attr(fit_simple_normal, "prior_list"))
})
