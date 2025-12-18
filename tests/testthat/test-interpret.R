# ============================================================================ #
# TEST FILE: Interpret Functions
# ============================================================================ #
#
# PURPOSE:
#   Tests for interpret and interpret2 functions that generate human-readable
#   summaries of Bayesian inference results.
#
# DEPENDENCIES:
#   - common-functions.R: test_reference_text, REFERENCE_DIR
#
# SKIP CONDITIONS:
#   - None (can run on CRAN - pure R with reference file testing)
#
# TAGS: @evaluation, @interpret, @output
# ============================================================================ #

REFERENCE_DIR <<- testthat::test_path("..", "results", "interpret")
source(testthat::test_path("common-functions.R"))


test_that("interpret2 function works", {

  set.seed(1)

  # Test basic interpret2 with all fields
  info1 <- list(
    list(
      inference_name        = "Effect",
      inference_BF_name     = "BF10",
      inference_BF          = 3.5,
      estimate_name         = "mu",
      estimate_samples      = rnorm(1000, 0.3, 0.15),
      estimate_units        = "kg",
      estimate_conditional  = FALSE
    )
  )

  result1 <- interpret2(info1, "RoBMA")
  test_reference_text(result1, "interpret2_basic.txt")
  expect_match(result1, "RoBMA found moderate evidence in favor of the Effect")
  expect_match(result1, "BF10 = 3.50")
  expect_match(result1, "model-averaged")
  expect_match(result1, "kg")

  # Test with conditional = TRUE
  info2 <- list(
    list(
      inference_name        = "Effect",
      inference_BF_name     = "BF10",
      inference_BF          = 15,
      estimate_name         = "mu",
      estimate_samples      = rnorm(1000, 0.5, 0.1),
      estimate_units        = NULL,
      estimate_conditional  = TRUE
    )
  )

  result2 <- interpret2(info2, "Test")
  test_reference_text(result2, "interpret2_conditional.txt")
  expect_match(result2, "strong evidence in favor")
  expect_match(result2, "conditional")
  expect_false(grepl("model-averaged", result2))

  # Test evidence against (BF < 1)
  info3 <- list(
    list(
      inference_name        = "Effect",
      inference_BF_name     = "BF01",
      inference_BF          = 0.1,
      estimate_name         = "mu",
      estimate_samples      = rnorm(1000, 0, 0.05),
      estimate_units        = NULL,
      estimate_conditional  = NULL
    )
  )

  result3 <- interpret2(info3, "Method")
  test_reference_text(result3, "interpret2_evidence_against.txt")
  expect_match(result3, "moderate evidence against the Effect")
  expect_match(result3, "BF01 = 0.100")

  # Test weak evidence
  info4 <- list(
    list(
      inference_name        = "Effect",
      inference_BF_name     = "BF",
      inference_BF          = 1.5,
      estimate_name         = "delta",
      estimate_samples      = rnorm(1000, 0.1, 0.1),
      estimate_units        = NULL,
      estimate_conditional  = FALSE
    )
  )

  result4 <- interpret2(info4, "Test")
  test_reference_text(result4, "interpret2_weak_evidence.txt")
  expect_match(result4, "weak evidence in favor")

  # Test without estimate samples (inference only)
  info5 <- list(
    list(
      inference_name        = "Bias",
      inference_BF_name     = "BF_pb",
      inference_BF          = 5
    )
  )

  result5 <- interpret2(info5, "RoBMA")
  test_reference_text(result5, "interpret2_inference_only.txt")
  expect_match(result5, "RoBMA found moderate evidence in favor of the Bias")
  expect_match(result5, "BF_pb = 5.00")
  expect_false(grepl("estimate", result5))

  # Test multiple specifications
  info6 <- list(
    list(
      inference_name        = "Effect",
      inference_BF_name     = "BF10",
      inference_BF          = 10,
      estimate_name         = "mu",
      estimate_samples      = rnorm(1000, 0.3, 0.1),
      estimate_units        = NULL,
      estimate_conditional  = FALSE
    ),
    list(
      inference_name        = "Bias",
      inference_BF_name     = "BF_pb",
      inference_BF          = 0.5
    )
  )

  result6 <- interpret2(info6, "Test")
  test_reference_text(result6, "interpret2_multiple.txt")
  expect_match(result6, "Effect")
  expect_match(result6, "Bias")

  # Test without method
  info7 <- list(
    list(
      inference_name        = "Effect",
      inference_BF_name     = "BF",
      inference_BF          = 2
    )
  )

  result7 <- interpret2(info7, NULL)
  test_reference_text(result7, "interpret2_no_method.txt")
  expect_match(result7, "found weak evidence")

})


test_that(".interpret.BF helper function works", {

  # Strong evidence in favor (BF > 10)
  result_strong_favor <- BayesTools:::.interpret.BF(15, "effect", "BF10")
  test_reference_text(result_strong_favor, "interpret_BF_strong_favor.txt")
  expect_match(result_strong_favor, "strong evidence in favor of the effect")
  expect_match(result_strong_favor, "BF10 = 15.00")

  # Moderate evidence in favor (3 < BF < 10)
  result_moderate_favor <- BayesTools:::.interpret.BF(5, "effect", NULL)
  test_reference_text(result_moderate_favor, "interpret_BF_moderate_favor.txt")
  expect_match(result_moderate_favor, "moderate evidence in favor")
  expect_match(result_moderate_favor, "BF = 5.00")

  # Weak evidence in favor (1 < BF < 3)
  result_weak_favor <- BayesTools:::.interpret.BF(1.5, "effect", "BF")
  test_reference_text(result_weak_favor, "interpret_BF_weak_favor.txt")
  expect_match(result_weak_favor, "weak evidence in favor")

  # Strong evidence against (BF < 0.1)
  result_strong_against <- BayesTools:::.interpret.BF(0.05, "effect", "BF01")
  test_reference_text(result_strong_against, "interpret_BF_strong_against.txt")
  expect_match(result_strong_against, "strong evidence against the effect")
  expect_match(result_strong_against, "BF01 = 0.050")

  # Moderate evidence against (0.1 <= BF < 1/3)
  result_moderate_against1 <- BayesTools:::.interpret.BF(0.1, "effect", NULL)
  test_reference_text(result_moderate_against1, "interpret_BF_moderate_against1.txt")
  expect_match(result_moderate_against1, "moderate evidence against")
  
  result_moderate_against2 <- BayesTools:::.interpret.BF(0.2, "effect", NULL)
  test_reference_text(result_moderate_against2, "interpret_BF_moderate_against2.txt")
  expect_match(result_moderate_against2, "moderate evidence against")

  # Weak evidence against (1/3 < BF < 1)
  result_weak_against <- BayesTools:::.interpret.BF(0.5, "effect", NULL)
  test_reference_text(result_weak_against, "interpret_BF_weak_against.txt")
  expect_match(result_weak_against, "weak evidence against")

})


test_that(".interpret.par helper function works", {

  set.seed(42)
  samples <- rnorm(10000, 0.5, 0.1)

  # Test model-averaged (conditional = FALSE)
  result1 <- BayesTools:::.interpret.par(samples, "mu", NULL, FALSE)
  test_reference_text(result1, "interpret_par_model_averaged.txt")
  expect_match(result1, "model-averaged estimate mu")
  expect_match(result1, "95% CI")

  # Test model-averaged (conditional = NULL)
  result2 <- BayesTools:::.interpret.par(samples, "delta", NULL, NULL)
  test_reference_text(result2, "interpret_par_model_averaged_null.txt")
  expect_match(result2, "model-averaged")

  # Test conditional
  result3 <- BayesTools:::.interpret.par(samples, "mu", NULL, TRUE)
  test_reference_text(result3, "interpret_par_conditional.txt")
  expect_match(result3, "conditional estimate mu")
  expect_false(grepl("model-averaged", result3))

  # Test with units
  result4 <- BayesTools:::.interpret.par(samples, "weight", "kg", FALSE)
  test_reference_text(result4, "interpret_par_with_units.txt")
  expect_match(result4, "kg")

})


test_that("interpret function input validation works", {

  # Test specification validation
  expect_error(interpret(list(), list(), "not a list", "Test"))

  # Test invalid specification elements
  expect_error(interpret(list(), list(), list(list(inference = 1)), "Test"))

})
