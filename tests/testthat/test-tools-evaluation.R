skip_if_not_test_profile("unit")

# ============================================================================ #
# TEST FILE: Utility Functions Evaluation Tests
# ============================================================================ #
#
# PURPOSE:
#   Tests for utility functions behavior (not input validation).
#   Includes .is.wholenumber, transformation checks, stan extraction, etc.
#
# DEPENDENCIES:
#   - No external packages required beyond testthat
#
# SKIP CONDITIONS:
#   - None (fast, pure R tests)
#
# MODELS/FIXTURES:
#   - None required
#
# TAGS: @evaluation, @fast
# ============================================================================ #


test_that(".is.wholenumber works correctly", {

  # Positive cases

  expect_true(BayesTools:::.is.wholenumber(0))
  expect_true(BayesTools:::.is.wholenumber(5))
  expect_true(BayesTools:::.is.wholenumber(-3))
  expect_true(BayesTools:::.is.wholenumber(1e10))

  # Negative cases
  expect_false(BayesTools:::.is.wholenumber(0.5))
  expect_false(BayesTools:::.is.wholenumber(1.1))
  expect_false(BayesTools:::.is.wholenumber(-3.5))

  # NA handling
  expect_true(is.na(BayesTools:::.is.wholenumber(NA)))
  expect_equal(BayesTools:::.is.wholenumber(NA, na.rm = TRUE), logical(0))

  # Vector input
  expect_equal(BayesTools:::.is.wholenumber(c(1, 2, 3.5)), c(TRUE, TRUE, FALSE))
  expect_equal(BayesTools:::.is.wholenumber(c(1, NA, 3.5)), c(TRUE, NA, FALSE))
})


test_that("transformation input validation works", {

  # Valid transformation
  expect_null(.check_transformation_input(transformation = list(
    "fun" = function(x) exp(x),
    "inv" = function(x) log(x),
    "jac" = function(x) exp(x)
  ), NULL, FALSE))

  # Missing 'jac' component
  expect_error(.check_transformation_input(transformation = list(
    "fun" = function(x) exp(x),
    "inv" = function(x) log(x),
    "err" = function(x) exp(x)
  ), NULL, FALSE), "The 'jac' objects are missing in the 'transformation' argument.")

  # Invalid format
  expect_error(.check_transformation_input(transformation = 1, NULL, FALSE),
               "Uknown format of the 'transformation' argument.")
})


test_that("stan extraction requires rstan fit", {
  expect_error(.extract_stan(NULL), "'fit' must be an rstan fit")
})


test_that("depreciation warnings work", {
  expect_warning(.depreciate.transform_orthonormal(TRUE, FALSE),
                 "'transform_orthonormal' argument will be depreciated in favor of 'transform_factors' argument.")
})
