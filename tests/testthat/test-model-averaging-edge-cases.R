context("Model averaging edge cases and comprehensive tests")

# Reference directory for text output comparisons
REFERENCE_DIR <<- testthat::test_path("..", "results", "model-averaging-edge-cases")

source(testthat::test_path("common-functions.R"))

# ============================================================================ #
# SECTION 1: mix_posteriors edge cases
# ============================================================================ #
test_that("mix_posteriors handles various prior types correctly", {

  skip_on_cran()
  skip_if_not_installed("rjags")
  skip_if_no_fits()

  # Load fits with margliks
  fit_simple_normal <- readRDS(file.path(temp_fits_dir, "fit_simple_normal.RDS"))
  marglik_simple_normal <- readRDS(file.path(temp_fits_dir, "fit_simple_normal_marglik.RDS"))

  fit_simple_spike <- readRDS(file.path(temp_fits_dir, "fit_simple_spike.RDS"))
  marglik_simple_spike <- readRDS(file.path(temp_fits_dir, "fit_simple_spike_marglik.RDS"))

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
  marglik_summary0 <- readRDS(file.path(temp_fits_dir, "fit_summary0_marglik.RDS"))

  fit_summary1 <- readRDS(file.path(temp_fits_dir, "fit_summary1.RDS"))
  marglik_summary1 <- readRDS(file.path(temp_fits_dir, "fit_summary1_marglik.RDS"))

  fit_summary2 <- readRDS(file.path(temp_fits_dir, "fit_summary2.RDS"))
  marglik_summary2 <- readRDS(file.path(temp_fits_dir, "fit_summary2_marglik.RDS"))

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
  marglik_orthonormal_0 <- readRDS(file.path(temp_fits_dir, "fit_orthonormal_0_marglik.RDS"))

  fit_orthonormal_1 <- readRDS(file.path(temp_fits_dir, "fit_orthonormal_1.RDS"))
  marglik_orthonormal_1 <- readRDS(file.path(temp_fits_dir, "fit_orthonormal_1_marglik.RDS"))

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


# ============================================================================ #
# SECTION 2: ensemble_inference edge cases
# ============================================================================ #
test_that("ensemble_inference handles different configurations", {

  skip_on_cran()
  skip_if_not_installed("rjags")
  skip_if_no_fits()

  # Load fits with margliks
  fit_simple_normal <- readRDS(file.path(temp_fits_dir, "fit_simple_normal.RDS"))
  marglik_simple_normal <- readRDS(file.path(temp_fits_dir, "fit_simple_normal_marglik.RDS"))

  fit_simple_spike <- readRDS(file.path(temp_fits_dir, "fit_simple_spike.RDS"))
  marglik_simple_spike <- readRDS(file.path(temp_fits_dir, "fit_simple_spike_marglik.RDS"))

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


test_that("models_inference computes correctly", {

  skip_on_cran()
  skip_if_not_installed("rjags")
  skip_if_no_fits()

  # Load fits with margliks
  fit_simple_normal <- readRDS(file.path(temp_fits_dir, "fit_simple_normal.RDS"))
  marglik_simple_normal <- readRDS(file.path(temp_fits_dir, "fit_simple_normal_marglik.RDS"))

  fit_simple_spike <- readRDS(file.path(temp_fits_dir, "fit_simple_spike.RDS"))
  marglik_simple_spike <- readRDS(file.path(temp_fits_dir, "fit_simple_spike_marglik.RDS"))

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
# SECTION 3: as_mixed_posteriors and as_marginal_inference
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


test_that("as_marginal_inference works correctly", {

  skip_on_cran()
  skip_if_not_installed("rjags")
  skip_if_no_fits()

  # as_marginal_inference requires a BayesTools_fit object - load one
  fit_simple_normal <- readRDS(file.path(temp_fits_dir, "fit_simple_normal.RDS"))

  # Get parameter names from the fitted model
  prior_list <- attr(fit_simple_normal, "prior_list")
  params <- names(prior_list)

  # Test basic as_marginal_inference call
  # This requires the model to have spike_and_slab or mixture priors
  # For now, just test that it errors correctly with non-matching parameters
  expect_error(
    as_marginal_inference(fit_simple_normal, marginal_parameters = "nonexistent"),
    regexp = NULL  # Any error is expected
  )

})


# ============================================================================ #
# SECTION 4: Inclusion BF edge cases
# ============================================================================ #
test_that("inclusion_BF handles edge cases", {

  # All null models - should return 0
  prior_probs <- c(0.5, 0.5)
  post_probs <- c(0.5, 0.5)
  is_null <- c(TRUE, TRUE)

  BF <- inclusion_BF(prior_probs = prior_probs, post_probs = post_probs, is_null = is_null)
  expect_equal(BF, 0)

  # All alternative models - should return Inf
  is_null <- c(FALSE, FALSE)
  BF <- inclusion_BF(prior_probs = prior_probs, post_probs = post_probs, is_null = is_null)
  expect_equal(BF, Inf)

  # Single model
  prior_probs <- 1
  post_probs <- 1
  is_null <- FALSE
  BF <- inclusion_BF(prior_probs = prior_probs, post_probs = post_probs, is_null = is_null)
  expect_equal(BF, Inf)

  # Test with marginal likelihoods only
  prior_probs <- c(0.5, 0.5)
  margliks <- c(-10, -10)  # Equal margliks
  is_null <- c(TRUE, FALSE)
  BF <- inclusion_BF(prior_probs = prior_probs, margliks = margliks, is_null = is_null)
  expect_equal(BF, 1)
})


# ============================================================================ #
# SECTION 5: weightfunctions_mapping edge cases
# ============================================================================ #
test_that("weightfunctions_mapping handles various configurations", {

  # Create one-sided weightfunction prior
  wf_onesided <- prior_weightfunction("one.sided", list(c(0.05), c(1, 1)))

  # Create two-sided weightfunction prior
  wf_twosided <- prior_weightfunction("two.sided", list(c(0.05), c(1, 1)))

  # Test with one-sided prior list
  mapping <- weightfunctions_mapping(list(wf_onesided))
  mapping_info <- paste0(
    "One-sided mapping:\n",
    "  Length: ", length(mapping), "\n",
    "  Is list: ", is.list(mapping)
  )

  # Test with two-sided prior list
  mapping_two <- weightfunctions_mapping(list(wf_twosided))

  # Test with one_sided = TRUE
  mapping_one <- weightfunctions_mapping(list(wf_twosided), one_sided = TRUE)

  # Test cuts_only = TRUE
  cuts <- weightfunctions_mapping(list(wf_onesided), cuts_only = TRUE)

  wf_mapping_info <- paste0(
    "One-sided mapping length: ", length(mapping), "\n",
    "Two-sided mapping length: ", length(mapping_two), "\n",
    "Two-sided with one_sided=TRUE length: ", length(mapping_one), "\n",
    "Cuts: ", paste(cuts, collapse = ", ")
  )
  test_reference_text(wf_mapping_info, "weightfunctions_mapping_info.txt")

})

