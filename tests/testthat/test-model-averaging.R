# ============================================================================ #
# TEST FILE: Model Averaging Functions
# ============================================================================ #
#
# PURPOSE:
#   Tests for compute_inference, ensemble_inference, and related Bayesian
#   model averaging functions in R/model-averaging.R
#
# DEPENDENCIES:
#   - bridgesampling: Required for marginal likelihood computation
#   - common-functions.R: Test helpers
#
# SKIP CONDITIONS:
#   - skip_if_not_installed("bridgesampling")
#
# MODELS/FIXTURES:
#   - Uses pre-computed marginal likelihoods, not pre-fitted models
#
# TAGS: @evaluation, @model-averaging
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
