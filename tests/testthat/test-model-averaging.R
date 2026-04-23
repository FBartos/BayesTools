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

# Reference directory for text output comparisons
REFERENCE_DIR <<- testthat::test_path("..", "results", "model-averaging")

source(testthat::test_path("common-functions.R"))


# ============================================================================ #
# SECTION 1: compute_inference tests
# ============================================================================ #
test_that("compute_inference works correctly", {

  skip_if_not_installed("bridgesampling")

  # Test basic inference computation
  prior_weights <- c(1, 1)
  margliks <- c(-10, -11)

  result <- compute_inference(prior_weights, margliks)

  expect_equal(length(result$prior_probs), 2)
  expect_equal(sum(result$prior_probs), 1)
  expect_equal(length(result$post_probs), 2)
  expect_true(sum(result$post_probs) > 0.99)  # Should be close to 1
  expect_true(is.numeric(result$BF))

  # Test with is_null as logical vector
  result2 <- compute_inference(prior_weights, margliks, is_null = c(TRUE, FALSE))
  expect_true(is.numeric(result2$BF))
  expect_true(attr(result2, "is_null")[1])
  expect_false(attr(result2, "is_null")[2])

  # Test with is_null as integer vector
  result3 <- compute_inference(prior_weights, margliks, is_null = 1)
  expect_true(attr(result3, "is_null")[1])
  expect_false(attr(result3, "is_null")[2])

  # Test conditional inference
  result4 <- compute_inference(prior_weights, margliks, is_null = c(TRUE, FALSE), conditional = TRUE)
  expect_equal(result4$prior_probs[1], 0)  # Null model should have 0 prior prob in conditional
  expect_equal(result4$prior_probs[2], 1)  # Alternative should have all weight
  expect_true(attr(result4, "conditional"))

  # Test with unequal prior weights
  result5 <- compute_inference(c(1, 3), margliks)
  expect_equal(result5$prior_probs[1], 0.25)
  expect_equal(result5$prior_probs[2], 0.75)

})


test_that("compute_inference input validation works", {

  skip_if_not_installed("bridgesampling")

  # Wrong is_null type
  expect_error(
    compute_inference(c(1, 1), c(-10, -11), is_null = "TRUE"),
    "must be either logical vector, integer vector, or NULL"
  )

  # is_null wrong length
  expect_error(
    compute_inference(c(1, 1), c(-10, -11), is_null = c(TRUE, FALSE, TRUE)),
    "must have length"
  )

  # mismatched lengths
  expect_error(
    compute_inference(c(1, 1, 1), c(-10, -11)),
    "must have length"
  )

})


test_that("inclusion_BF works correctly", {

  # Test with posterior probabilities
  prior_probs <- c(0.5, 0.5)
  post_probs <- c(0.8, 0.2)
  is_null <- c(TRUE, FALSE)

  BF <- inclusion_BF(prior_probs = prior_probs, post_probs = post_probs, is_null = is_null)
  expect_true(is.numeric(BF))
  # With equal prior, BF should be 0.2/0.8 = 0.25 (for alternative vs null)
  expect_equal(BF, 0.25)

  # Test with marginal likelihoods
  margliks <- c(-10, -10)  # Equal margliks
  BF2 <- inclusion_BF(prior_probs = prior_probs, margliks = margliks, is_null = is_null)
  expect_equal(BF2, 1)  # Should be 1 with equal margliks and equal priors

  # Test with integer is_null
  BF3 <- inclusion_BF(prior_probs = prior_probs, post_probs = post_probs, is_null = 1)
  expect_equal(BF3, BF)

  # Test all null scenario
  BF_all_null <- inclusion_BF(prior_probs = prior_probs, post_probs = post_probs, is_null = c(TRUE, TRUE))
  expect_equal(BF_all_null, 0)

  # Test all alternative scenario
  BF_all_alt <- inclusion_BF(prior_probs = prior_probs, post_probs = post_probs, is_null = c(FALSE, FALSE))
  expect_equal(BF_all_alt, Inf)

})


test_that("inclusion_BF input validation works", {

  # Wrong is_null type
  expect_error(
    inclusion_BF(prior_probs = c(0.5, 0.5), post_probs = c(0.8, 0.2), is_null = "TRUE"),
    "must be either logical vector, integer vector, or NULL"
  )

  # Missing arguments
  expect_error(
    inclusion_BF(prior_probs = c(0.5, 0.5), is_null = c(TRUE, FALSE)),
    "'prior_probs' and either 'post_probs' or 'marglik' must be specified"
  )

})


test_that(".inclusion_BF.probs edge cases work", {

  # Test when posterior is fully concentrated on alternative
  prior_probs <- c(0.5, 0.5)
  post_probs <- c(0, 1)
  is_null <- c(TRUE, FALSE)

  BF <- BayesTools:::.inclusion_BF.probs(prior_probs, post_probs, is_null)
  expect_equal(BF, Inf)

  # Test when posterior is fully concentrated on null
  post_probs2 <- c(1, 0)
  BF2 <- BayesTools:::.inclusion_BF.probs(prior_probs, post_probs2, is_null)
  expect_equal(BF2, 0)

})


test_that("weightfunctions_mapping works correctly", {

  # Create weightfunction priors
  wf_onesided <- prior_weightfunction("one.sided", list(c(0.05), c(1, 1)))
  wf_twosided <- prior_weightfunction("two.sided", list(c(0.05), c(1, 1)))

  # Test with single weightfunction
  mapping1 <- weightfunctions_mapping(list(wf_onesided))
  expect_true(is.list(mapping1))

  # Test with cuts_only = TRUE
  cuts <- weightfunctions_mapping(list(wf_onesided), cuts_only = TRUE)
  expect_true(is.numeric(cuts))
  expect_true(0 %in% cuts)
  expect_true(1 %in% cuts)

  # Test with one_sided = TRUE
  mapping2 <- weightfunctions_mapping(list(wf_twosided), one_sided = TRUE)
  expect_true(is.list(mapping2))

  # Test with point prior (should be handled gracefully)
  p_point <- prior("point", list(1))
  mapping3 <- weightfunctions_mapping(list(wf_onesided, p_point))
  expect_true(is.list(mapping3))

})


test_that("weightfunctions_mapping input validation works", {

  # Non-weightfunction priors should fail
  p_normal <- prior("normal", list(0, 1))
  expect_error(
    weightfunctions_mapping(list(p_normal)),
    "must be a list of weightfunction priors"
  )

})

# ============================================================================ #
# SECTION 5: mix_posteriors tests
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
  # Capture a summary of the mixed posteriors structure for reference
  mixed_info <- paste0(
    "Class: ", paste(class(mixed), collapse = ", "), "\n",
    "Parameters: ", paste(names(mixed), collapse = ", "), "\n",
    "Sample size m: ", length(mixed$m), "\n",
    "Sample size s: ", length(mixed$s)
  )
  test_reference_text(mixed_info, "mix_posteriors_simple_info.txt")
  expect_equal(length(mixed$m), 1000)
  expect_equal(length(mixed$s), 1000)

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
})


test_that("mix_posteriors handles vector priors", {

  skip_on_cran()
  skip_if_not_installed("rjags")
  skip_if_no_fits()

  # Load vector prior models
  fit_vector_mnormal <- readRDS(file.path(temp_fits_dir, "fit_vector_mnormal.RDS"))

  # Create a mock marglik for testing (we only need the structure)
  mock_marglik <- structure(
    list(logml = -100, niter = 1000, method = "warp3"),
    class = "bridge"
  )

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
  mock_marglik <- structure(
    list(logml = -100, niter = 1000, method = "warp3"),
    class = "bridge"
  )
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
  inference_int_info <- paste0(
    "BF: ", round(inference_int$m$BF, 4), "\n",
    "is_null: ", paste(attr(inference_int$m, "is_null"), collapse = ", "), "\n",
    "prior_probs: ", paste(round(inference_int$m$prior_probs, 4), collapse = ", "), "\n",
    "post_probs: ", paste(round(inference_int$m$post_probs, 4), collapse = ", ")
  )
  test_reference_text(inference_int_info, "ensemble_inference_int_spec.txt")

  # Test conditional inference
  inference_cond <- ensemble_inference(
    model_list = models,
    parameters = "m",
    is_null_list = list("m" = c(FALSE, TRUE)),
    conditional = TRUE
  )

  expect_true(attr(inference_cond, "conditional"))
  inference_cond_info <- paste0(
    "Conditional: ", attr(inference_cond, "conditional"), "\n",
    "BF: ", round(inference_cond$m$BF, 4)
  )
  test_reference_text(inference_cond_info, "ensemble_inference_conditional.txt")

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

  # Create reference output for models_inference structure
  models_inf_info <- paste0(
    "Model 1 inference:\n",
    "  m_number: ", models_with_inference[[1]]$inference$m_number, "\n",
    "  prior_prob: ", round(models_with_inference[[1]]$inference$prior_prob, 6), "\n",
    "  post_prob: ", round(models_with_inference[[1]]$inference$post_prob, 6), "\n",
    "Model 2 inference:\n",
    "  m_number: ", models_with_inference[[2]]$inference$m_number, "\n",
    "  prior_prob: ", round(models_with_inference[[2]]$inference$prior_prob, 6), "\n",
    "  post_prob: ", round(models_with_inference[[2]]$inference$post_prob, 6), "\n",
    "Total post_prob: ", round(sum(sapply(models_with_inference, function(m) m$inference$post_prob)), 6)
  )
  test_reference_text(models_inf_info, "models_inference_output.txt")

  # Check prior probs reflect weights (1:2 ratio)
  expect_equal(models_with_inference[[1]]$inference$prior_prob, 1/3, tolerance = 1e-10)
  expect_equal(models_with_inference[[2]]$inference$prior_prob, 2/3, tolerance = 1e-10)

  # Check posterior probs sum to 1
  total_post_prob <- sum(sapply(models_with_inference, function(m) m$inference$post_prob))
  expect_equal(total_post_prob, 1, tolerance = 1e-10)

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
})
