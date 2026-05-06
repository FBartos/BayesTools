skip_if_not_test_profile("unit")

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

  result_none <- BayesTools:::.interpret.BF(1, "effect", NULL)
  expect_match(result_none, "no evidence for or against the effect")
  expect_match(result_none, "BF = 1.00")

})


test_that(".interpret.BF preserves threshold boundary semantics", {

  boundary_cases <- list(
    list(BF = 1 / 10, expected = "moderate evidence against the effect, BF = 0.100"),
    list(BF = 1 / 3,  expected = "weak evidence against the effect, BF = 0.333"),
    list(BF = 1,      expected = "no evidence for or against the effect, BF = 1.00"),
    list(BF = 3,      expected = "weak evidence in favor of the effect, BF = 3.00"),
    list(BF = 10,     expected = "moderate evidence in favor of the effect, BF = 10.00")
  )

  for(case in boundary_cases){
    expect_identical(
      BayesTools:::.interpret.BF(case$BF, "effect", NULL),
      case$expected
    )
  }

})


test_that(".interpret.BF reports finite-sample BF bounds", {

  inclusion_BF <- 9
  attr(inclusion_BF, "bound_operator") <- ">"
  expect_identical(
    BayesTools:::.interpret.BF(inclusion_BF, "effect", "Inclusion BF"),
    "at least moderate evidence in favor of the effect, Inclusion BF > 9.00"
  )

  bounded_raw <- c(9, 1 / 9)
  attr(bounded_raw, "bound_operator") <- c(">", "<")
  bounded_column <- format_BF(bounded_raw, inclusion = TRUE)
  expect_identical(
    BayesTools:::.interpret.BF(bounded_column[1], "effect", attr(bounded_column, "name")),
    "at least moderate evidence in favor of the effect, Inclusion BF > 9.00"
  )

  exclusion_BF <- 1 / 9
  attr(exclusion_BF, "bound_operator") <- "<"
  expect_identical(
    BayesTools:::.interpret.BF(exclusion_BF, "effect", "Inclusion BF"),
    "at least moderate evidence against the effect, Inclusion BF < 0.111"
  )

  expect_identical(
    interpret2(
      list(list(
        inference_name = "effect",
        inference_BF = 9,
        inference_BF_name = "Inclusion BF",
        inference_BF_bound_operator = ">"
      )),
      "Method"
    ),
    "Method found at least moderate evidence in favor of the effect, Inclusion BF > 9.00."
  )

  expect_identical(
    interpret(
      list(effect = list(BF = 9, BF_bound_operator = ">")),
      list(dummy = 1),
      list(list(inference = "effect", inference_name = "effect", inference_BF_name = "Inclusion BF")),
      "Method"
    ),
    "Method found at least moderate evidence in favor of the effect, Inclusion BF > 9.00."
  )

})


test_that(".interpret.BF rejects invalid Bayes factors before formatting", {

  invalid_BFs <- list(0, Inf, NA_real_, NaN, -1)

  for(BF in invalid_BFs){
    expect_error(
      BayesTools:::.interpret.BF(BF, "effect", NULL),
      "inference_BF"
    )
  }

  expect_error(
    BayesTools:::.interpret.BF("2", "effect", NULL),
    "numeric vector"
  )
  expect_error(
    interpret2(list(list(inference_name = "effect", inference_BF = "2")), "Method"),
    "numeric vector"
  )
  expect_error(
    interpret(
      list(effect = list(BF = "2")),
      list(dummy = 1),
      list(list(inference = "effect", inference_name = "effect")),
      "Method"
    ),
    "numeric vector"
  )

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


test_that(".interpret.par rejects empty and malformed estimate samples", {

  expect_error(
    BayesTools:::.interpret.par(numeric(0), "mu", NULL, FALSE),
    "estimate_samples"
  )
  expect_error(
    BayesTools:::.interpret.par(c(NA_real_, NA_real_), "mu", NULL, FALSE),
    "estimate_samples"
  )
  expect_error(
    BayesTools:::.interpret.par(c(1, Inf), "mu", NULL, FALSE),
    "estimate_samples"
  )
  expect_error(
    BayesTools:::.interpret.par(c(1, NaN), "mu", NULL, FALSE),
    "estimate_samples"
  )
  expect_error(
    BayesTools:::.interpret.par("not numeric", "mu", NULL, FALSE),
    "estimate_samples"
  )
  expect_error(
    BayesTools:::.interpret.par(c(1, 2), NULL, NULL, FALSE),
    "estimate_name"
  )

})


test_that("interpret2 maps wrappers to core helpers and validates missing fields", {

  info <- list(
    list(
      inference_name        = "Effect",
      inference_BF_name     = "BF10",
      inference_BF          = 3,
      estimate_name         = "mu",
      estimate_samples      = c(1, 2, 3),
      estimate_units        = "kg",
      estimate_conditional  = TRUE
    )
  )

  expected <- paste0(
    "Method found ",
    BayesTools:::.interpret.BF(3, "Effect", "BF10"),
    ", ",
    BayesTools:::.interpret.par(c(1, 2, 3), "mu", "kg", TRUE),
    "."
  )
  expect_identical(interpret2(info, "Method"), expected)

  inference_only <- list(list(inference_name = "Effect", inference_BF = 2))
  expect_identical(
    interpret2(inference_only, NULL),
    paste0(" found ", BayesTools:::.interpret.BF(2, "Effect", NULL), ".")
  )

  expect_error(
    interpret2(list(list(inference_name = "Effect")), "Method"),
    "inference_BF"
  )
  expect_error(
    interpret2(list(list(inference_name = "Effect", inference_BF = 2, estimate_samples = c(1, 2))), "Method"),
    "estimate_name"
  )

})


test_that("interpret wrapper maps named inference and samples to core helpers", {

  inference <- list(effect = list(BF = 10))
  samples <- list(theta = c(1, 2, 3))
  specification <- list(
    list(
      inference              = "effect",
      inference_name         = "Effect",
      inference_BF_name      = "BF10",
      samples                = "theta",
      samples_name           = "mu",
      samples_units          = NULL,
      samples_conditional    = FALSE
    )
  )

  expected <- paste0(
    "Method found ",
    BayesTools:::.interpret.BF(10, "Effect", "BF10"),
    ", ",
    BayesTools:::.interpret.par(c(1, 2, 3), "mu", NULL, FALSE),
    "."
  )
  expect_identical(interpret(inference, samples, specification, "Method"), expected)

  expect_error(
    interpret(list(effect = list()), list(dummy = 1), list(list(inference = "effect")), "Method"),
    "inference_BF"
  )
  expect_error(
    interpret(list(effect = list(BF = 2)), list(dummy = 1), list(list(inference = "effect", samples = "theta")), "Method"),
    "theta.*missing|missing.*theta"
  )

})


test_that("interpret_records normalizes ordered table and direct-record sources", {

  effect_BF <- 9
  attr(effect_BF, "bound_operator") <- ">"
  effect_BF <- format_BF(effect_BF, logBF = TRUE, BF01 = TRUE, inclusion = TRUE)

  component_tests <- data.frame(
    prior_prob = 0.5,
    post_prob = 1,
    check.names = FALSE
  )
  component_tests[["inclusion_BF"]] <- effect_BF
  rownames(component_tests) <- "Effect"
  class(component_tests) <- c("BayesTools_table", "data.frame")
  attr(component_tests, "type") <- c("prior_prob", "post_prob", "inclusion_BF")

  moderator_BF <- format_BF(c(2, 1 / 4), inclusion = TRUE)
  moderator_tests <- data.frame(
    prior_prob = c(0.5, 0.5),
    post_prob = c(2 / 3, 0.2),
    check.names = FALSE
  )
  moderator_tests[["inclusion_BF"]] <- moderator_BF
  rownames(moderator_tests) <- c("x1", "x2")
  class(moderator_tests) <- c("BayesTools_table", "data.frame")
  attr(moderator_tests, "type") <- c("prior_prob", "post_prob", "inclusion_BF")

  moderator_estimates <- data.frame(
    Mean = c(0.20, -0.10),
    "0.025" = c(0.05, -0.30),
    "0.975" = c(0.35, 0.10),
    check.names = FALSE
  )
  rownames(moderator_estimates) <- c("x1", "x2")

  sources <- list(
    component_tests = component_tests,
    pooled_effect = list(
      type = "record",
      data = list(
        kind = "estimate",
        parameter = "effect odds ratio",
        central_name = "mode",
        central_value = 1.25,
        lower_value = 1.05,
        upper_value = 1.50,
        lower_prob = 0.025,
        upper_prob = 0.975,
        interval_level = 0.95,
        conditioning = "conditional on effect inclusion"
      )
    ),
    moderator_tests = moderator_tests,
    moderator_estimates = list(
      data = moderator_estimates,
      schema = list(
        central = "Mean",
        lower = "0.025",
        upper = "0.975",
        units = "d",
        conditioning = "model-averaged"
      )
    )
  )

  plan <- list(
    list(kind = "header", section = "model", item_id = "header", order = 0, text = "RoBMA model."),
    list(
      kind = "pair",
      section = "primary",
      item_id = "effect",
      order = 10,
      evidence = list(source = "component_tests", row = "Effect", label = "the effect"),
      estimate = list(source = "pooled_effect", label = "pooled effect")
    ),
    list(
      kind = "for_each",
      section = "moderators",
      item_id = "moderator",
      order = 100,
      source = "moderator_tests",
      pair_with = "moderator_estimates",
      rows = "source_order"
    )
  )

  records <- interpret_records(sources, plan)

  expect_s3_class(records, "BayesTools_interpret_records")
  expect_equal(records$kind, c("header", "evidence", "estimate", "evidence", "estimate", "evidence", "estimate"))
  expect_equal(records$record_id[1:3], c("model.header.header", "primary.effect.evidence", "primary.effect.estimate"))
  expect_equal(records$source[2:3], c("component_tests", "pooled_effect"))
  expect_equal(records$row[4:7], c("x1", "x1", "x2", "x2"))

  effect <- records[records$record_id == "primary.effect.evidence", ]
  expect_equal(effect$BF_value, log(1 / 9), tolerance = 1e-12)
  expect_equal(effect$BF_scale, "log")
  expect_equal(effect$BF_orientation, "exclusion_over_inclusion")
  expect_equal(effect$BF_bound_operator, "<")
  expect_equal(effect$BF_canonical_value, 9, tolerance = 1e-12)
  expect_equal(effect$BF_canonical_bound_operator, ">")

  estimate <- records[records$record_id == "primary.effect.estimate", ]
  expect_equal(estimate$central_name, "mode")
  expect_equal(estimate$central_value, 1.25)
  expect_equal(estimate$conditioning, "conditional on effect inclusion")

  moderator_estimate <- records[records$record_id == "moderators.moderator.x1.estimate", ]
  expect_equal(moderator_estimate$central_name, "mean")
  expect_equal(moderator_estimate$lower_prob, 0.025)
  expect_equal(moderator_estimate$upper_prob, 0.975)
  expect_equal(moderator_estimate$interval_level, 0.95)
  expect_equal(moderator_estimate$units, "d")

  text <- interpret_records(sources, plan, output = "text")
  expect_match(paste(text, collapse = "\n"), "Inclusion BF > 9.00", fixed = TRUE)
  expect_match(paste(text, collapse = "\n"), "conditional on effect inclusion", fixed = TRUE)

})


test_that("interpret_records selects direct records by record_id and row columns", {

  records_source <- data.frame(
    record_id = c("bf.record", "estimate.record"),
    row = c("bf-row", "estimate-row"),
    kind = c("evidence", "estimate"),
    label = c("effect", "theta"),
    BF_value = c(4, NA),
    BF_orientation = c("alternative_over_null", NA),
    central_name = c(NA, "mean"),
    central_value = c(NA, 1.5),
    stringsAsFactors = FALSE
  )
  rownames(records_source) <- c("not-bf", "not-estimate")

  sources <- list(records = list(type = "records", data = records_source))
  plan <- list(
    list(kind = "evidence", source = "records", row = "bf.record", section = "s", item_id = "bf"),
    list(kind = "estimate", source = "records", row = "estimate-row", section = "s", item_id = "est")
  )

  out <- interpret_records(sources, plan)

  expect_equal(out$record_id, c("bf.record", "estimate.record"))
  expect_equal(out$row, c("bf.record", "estimate-row"))
  expect_equal(out$BF_canonical_value[1], 4)
  expect_equal(out$central_value[2], 1.5)
})


test_that("interpret_records handles ambiguous multi-row sources according to missing policy", {

  table <- data.frame(BF = c(2, 3), check.names = FALSE)
  rownames(table) <- c("a", "b")

  sources <- list(tests = table)
  plan <- list(list(kind = "evidence", source = "tests", section = "s", item_id = "missing-row"))

  expect_error(
    interpret_records(sources, plan),
    "multiple rows; specify 'row'",
    fixed = TRUE
  )

  skipped <- interpret_records(sources, plan, missing = "skip")
  expect_s3_class(skipped, "BayesTools_interpret_records")
  expect_equal(nrow(skipped), 0)

  expect_warning(
    warned <- interpret_records(sources, plan, missing = "warn"),
    "multiple rows; specify 'row'",
    fixed = TRUE
  )
  expect_equal(nrow(warned), 0)
})


test_that("interpret_records matches explicitly requested padded probability columns", {

  estimates <- data.frame(
    Mean = 0,
    "0.100" = -2,
    "0.200" = -1,
    "0.800" = 1,
    "0.900" = 2,
    check.names = FALSE
  )
  rownames(estimates) <- "theta"

  out <- interpret_records(
    sources = list(estimates = estimates),
    plan = list(list(
      kind = "estimate",
      source = "estimates",
      row = "theta",
      lower_prob = .2,
      upper_prob = .8
    ))
  )

  expect_equal(out$central_name, "mean")
  expect_equal(out$lower_value, -1)
  expect_equal(out$upper_value, 1)
  expect_equal(out$lower_prob, .2)
  expect_equal(out$upper_prob, .8)
  expect_equal(out$interval_level, .6)
})


test_that("interpret_tables aliases interpret_records and supports optional missing entries", {

  table <- data.frame(
    BF = 4,
    check.names = FALSE
  )
  rownames(table) <- "joint"

  sources <- list(joint = list(
    data = table,
    schema = list(
      BF = "BF",
      BF_orientation = "alternative_over_null",
      BF_scale = "linear"
    )
  ))
  spec <- list(
    list(kind = "evidence", source = "joint", row = "joint", section = "moderators", item_id = "joint", label = "moderators"),
    list(kind = "evidence", source = "missing_joint", optional = TRUE, section = "moderators", item_id = "optional")
  )

  records <- interpret_tables(sources, spec)

  expect_equal(nrow(records), 1)
  expect_equal(records$record_id, "moderators.joint.evidence")
  expect_equal(records$BF_canonical_value, 4)

})


test_that("interpret function input validation works", {

  # Test specification validation
  expect_error(interpret(list(), list(), "not a list", "Test"))

  # Test invalid specification elements
  expect_error(interpret(list(), list(), list(list(inference = 1)), "Test"))

})
