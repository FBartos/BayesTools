# ============================================================================ #
# TEST FILE: Model Averaging Edge Cases
# ============================================================================ #
#
# PURPOSE:
#   Edge case tests for model averaging functions including input validation,
#   boundary conditions for Bayes factors, and weightfunction mapping edge cases.
#
# DEPENDENCIES:
#   - common-functions.R: test_reference_text
#
# SKIP CONDITIONS:
#   - None (these are simple edge case tests that don't require fitted models)
#
# MODELS/FIXTURES:
#   - None required
#
# TAGS: @edge-cases, @model-averaging, @input-validation
# ============================================================================ #

# Reference directory for text output comparisons
REFERENCE_DIR <<- testthat::test_path("..", "results", "model-averaging-edge-cases")

source(testthat::test_path("common-functions.R"))


# ============================================================================ #
# SECTION 1: inclusion_BF boundary conditions
# ============================================================================ #
test_that("inclusion_BF handles all-null models", {

  # All null models - should return 0
  prior_probs <- c(0.5, 0.5)
  post_probs <- c(0.5, 0.5)
  is_null <- c(TRUE, TRUE)

  BF <- inclusion_BF(prior_probs = prior_probs, post_probs = post_probs, is_null = is_null)
  expect_equal(BF, 0)

})


test_that("inclusion_BF handles all-alternative models", {

  # All alternative models - should return Inf
  prior_probs <- c(0.5, 0.5)
  post_probs <- c(0.5, 0.5)
  is_null <- c(FALSE, FALSE)

  BF <- inclusion_BF(prior_probs = prior_probs, post_probs = post_probs, is_null = is_null)
  expect_equal(BF, Inf)

})


test_that("inclusion_BF handles single model case", {

  prior_probs <- 1
  post_probs <- 1
  is_null <- FALSE

  BF <- inclusion_BF(prior_probs = prior_probs, post_probs = post_probs, is_null = is_null)
  expect_equal(BF, Inf)

  # Single null model
  is_null <- TRUE
  BF <- inclusion_BF(prior_probs = prior_probs, post_probs = post_probs, is_null = is_null)
  expect_equal(BF, 0)

})


test_that("inclusion_BF works with marginal likelihoods only", {

  # Test with marginal likelihoods instead of posterior probs
  prior_probs <- c(0.5, 0.5)
  margliks <- c(-10, -10)  # Equal margliks
  is_null <- c(TRUE, FALSE)

  BF <- inclusion_BF(prior_probs = prior_probs, margliks = margliks, is_null = is_null)
  expect_equal(BF, 1)

  # Unequal margliks - alternative has higher marglik
  margliks <- c(-10, -8)  # Alternative model is better
  BF <- inclusion_BF(prior_probs = prior_probs, margliks = margliks, is_null = is_null)
  expect_true(BF > 1)

  # Unequal margliks - null has higher marglik
  margliks <- c(-8, -10)  # Null model is better
  BF <- inclusion_BF(prior_probs = prior_probs, margliks = margliks, is_null = is_null)
  expect_true(BF < 1)

})


test_that("inclusion_BF ignores zero-prior marginal likelihoods in stability algebra", {
  expect_equal(
    inclusion_BF(
      prior_probs = c(0, 0.5, 0.5),
      margliks    = c(1000, 0, 0),
      is_null     = c(FALSE, TRUE, FALSE)
    ),
    1
  )

  expect_equal(
    inclusion_BF(
      prior_probs = c(0, 1),
      margliks    = c(1000, 0),
      is_null     = c(FALSE, TRUE)
    ),
    0
  )

  expect_equal(
    inclusion_BF(
      prior_probs = c(1, 0),
      margliks    = c(0, 1000),
      is_null     = c(FALSE, TRUE)
    ),
    Inf
  )

  expect_equal(
    inclusion_BF(
      prior_probs = c(0, 0.5, 0.5),
      margliks    = c(Inf, 0, 0),
      is_null     = c(FALSE, TRUE, FALSE)
    ),
    1
  )
})


# ============================================================================ #
# SECTION 2: weightfunctions_mapping edge cases
# ============================================================================ #
test_that("weightfunctions_mapping handles one-sided priors", {

  # Create one-sided weightfunction prior
  wf_onesided <- prior_weightfunction("one-sided", c(0.05), wf_cumulative(c(1, 1)))

  mapping <- weightfunctions_mapping(list(wf_onesided))

  expect_true(is.list(mapping))
  expect_equal(mapping[[1]], c(1, 2))

})


test_that("weightfunctions_mapping handles two-sided priors", {

  # Create two-sided weightfunction prior
  wf_twosided <- prior_weightfunction("two-sided", c(0.05), wf_cumulative(c(1, 1)))

  mapping <- weightfunctions_mapping(list(wf_twosided))

  expect_true(is.list(mapping))
  expect_equal(mapping[[1]], c(1, 2))

})


test_that("weightfunctions_mapping handles one_sided argument", {

  # Create two-sided weightfunction prior
  wf_twosided <- prior_weightfunction("two-sided", c(0.05), wf_cumulative(c(1, 1)))

  # Test with one_sided = TRUE
  mapping_one <- weightfunctions_mapping(list(wf_twosided), one_sided = TRUE)

  expect_true(is.list(mapping_one))
  expect_equal(mapping_one[[1]], c(1, 2, 1))

})


test_that("weightfunctions_mapping cuts_only option works", {

  # Create one-sided weightfunction prior
  wf_onesided <- prior_weightfunction("one-sided", c(0.05), wf_cumulative(c(1, 1)))

  # Test cuts_only = TRUE
  cuts <- weightfunctions_mapping(list(wf_onesided), cuts_only = TRUE)

  expect_equal(cuts, c(0.00, 0.05, 1.00))
})


test_that("weightfunctions_mapping handles mixed prior list", {

  # Multiple weightfunction priors with different configurations
  wf_onesided <- prior_weightfunction("one-sided", c(0.05), wf_cumulative(c(1, 1)))
  wf_twosided <- prior_weightfunction("two-sided", c(0.05, 0.10), wf_cumulative(c(1, 1, 1)))

  mapping <- weightfunctions_mapping(list(wf_onesided, wf_twosided))

  expect_true(is.list(mapping))

  wf_mapping_info <- paste0(
    "Mixed mapping length: ", length(mapping), "\n",
    "Inx 1: ", paste0(mapping[[1]], collapse = ","), "\n",
    "Inx 2: ", paste0(mapping[[2]], collapse = ",")
  )
  test_reference_text(wf_mapping_info, "weightfunctions_mapping_info.txt")

})
