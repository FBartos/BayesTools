# ============================================================================ #
# TEST FILE: Model Averaging Plots Edge Cases
# ============================================================================ #
#
# PURPOSE:
#   Edge case tests for plot functions including input validation and
#   error handling for invalid prior configurations.
#
# DEPENDENCIES:
#   - None (pure R testing)
#
# SKIP CONDITIONS:
#   - None (can run on CRAN)
#
# MODELS/FIXTURES:
#   - None required
#
# TAGS: @edge-cases, @plots, @input-validation
# ============================================================================ #


# ============================================================================ #
# SECTION 1: plot_prior_list input validation
# ============================================================================ #
test_that("plot_prior_list rejects non-list input", {

  expect_error(
    plot_prior_list(prior("normal", list(0, 1))),
    "must be a list of priors"
  )

})


test_that("plot_prior_list rejects PET-PEESE without prior_list_mu", {

  pet_list <- list(
    p1 = prior_PET("normal", list(0, 1))
  )
  expect_error(
    plot_prior_list(pet_list),
    "prior_list_mu"
  )

})


test_that("plot_prior_list rejects prior_list_mu when not needed", {

  simple_list <- list(
    p1 = prior("normal", list(0, 1))
  )
  expect_error(
    plot_prior_list(simple_list, prior_list_mu = list(prior("spike", list(0)))),
    "prior_list_mu"
  )

})


