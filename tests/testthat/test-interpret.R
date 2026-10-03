skip_if_not_test_profile("unit")

# ============================================================================ #
# TEST FILE: Interpret Functions
# ============================================================================ #
#
# PURPOSE:
#   Tests for interpret and interpret_records functions that generate
#   human-readable summaries of Bayesian inference results.
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
  expect_equal(estimate$interval_level, 0.95)

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


test_that("interpret_records labels fallback intervals with their own probabilities", {

  estimates <- data.frame(
    Mean = 0.3,
    "0.025" = 0.1,
    "0.5" = 0.3,
    "0.975" = 0.5,
    check.names = FALSE
  )
  rownames(estimates) <- "mu"
  source <- list(
    type = "table",
    data = estimates,
    schema = list(lower_prob = 0.05, upper_prob = 0.95)
  )
  plan <- list(list(kind = "estimate", source = "est", row = "mu"))

  out <- interpret_records(sources = list(est = source), plan = plan)
  expect_equal(out$lower_value, 0.1)
  expect_equal(out$upper_value, 0.5)
  expect_equal(out$lower_prob, 0.025)
  expect_equal(out$upper_prob, 0.975)
  expect_equal(out$interval_level, 0.95)

  text <- interpret_records(sources = list(est = source), plan = plan, output = "text")
  expect_match(text, "95%", fixed = TRUE)
  expect_false(grepl("90%", text, fixed = TRUE))
})


test_that("interpret_records supports central-only estimate tables", {

  estimates <- ensemble_estimates_table(
    samples = list(theta = c(-1, 0, 2)),
    parameters = "theta",
    probs = NULL
  )
  plan <- list(list(
    kind = "estimate",
    source = "estimates",
    row = "theta"
  ))

  out <- interpret_records(
    sources = list(estimates = estimates),
    plan = plan
  )
  text <- interpret_records(
    sources = list(estimates = estimates),
    plan = plan,
    output = "text"
  )

  expect_equal(out$central_name, "mean")
  expect_equal(out$central_value, mean(c(-1, 0, 2)))
  expect_true(is.na(out$lower_value))
  expect_true(is.na(out$upper_value))
  expect_true(is.na(out$interval_level))
  expect_false(grepl("interval", text, fixed = TRUE))
})

test_that("interpret_records derives interval levels from endpoint probabilities", {
  estimate <- list(
    kind = "estimate",
    parameter = "theta",
    central_name = "mean",
    central_value = 0,
    lower_value = -1,
    upper_value = 1,
    lower_prob = 0.025,
    upper_prob = 0.975
  )
  out <- interpret_records(
    sources = list(estimate = estimate),
    plan = list(list(kind = "estimate", source = "estimate"))
  )
  expect_equal(out$interval_level, 0.95)

  arbitrary <- data.frame(
    Mean = 0,
    lower = -1,
    upper = 1,
    check.names = FALSE
  )
  arbitrary_out <- interpret_records(
    sources = list(arbitrary = list(
      data = arbitrary,
      schema = list(
        central = "Mean",
        lower = "lower",
        upper = "upper"
      )
    )),
    plan = list(list(
      kind = "estimate",
      source = "arbitrary",
      row = 1
    ))
  )
  arbitrary_text <- interpret_records(
    sources = list(arbitrary = list(
      data = arbitrary,
      schema = list(
        central = "Mean",
        lower = "lower",
        upper = "upper"
      )
    )),
    plan = list(list(
      kind = "estimate",
      source = "arbitrary",
      row = 1
    )),
    output = "text"
  )
  expect_true(is.na(arbitrary_out$interval_level))
  expect_match(arbitrary_text, "uncertainty interval", fixed = TRUE)

  invalid_inputs <- list(
    schema = function(){
      interpret_records(
        sources = list(estimates = list(
          data = arbitrary,
          schema = list(
            central = "Mean",
            lower = "lower",
            upper = "upper",
            interval_level = 0.95
          )
        )),
        plan = list(list(kind = "estimate", source = "estimates"))
      )
    },
    plan = function(){
      interpret_records(
        sources = list(estimates = arbitrary),
        plan = list(list(
          kind = "estimate",
          source = "estimates",
          interval_level = 0.95
        ))
      )
    },
    record = function(){
      bad_record <- estimate
      bad_record$interval_level <- 0.95
      interpret_records(
        sources = list(estimate = bad_record),
        plan = list(list(kind = "estimate", source = "estimate"))
      )
    },
    records = function(){
      bad_records <- as.data.frame(estimate, check.names = FALSE)
      bad_records$interval_level <- 0.95
      interpret_records(
        sources = list(estimates = list(
          type = "records",
          data = bad_records
        )),
        plan = list(list(kind = "estimate", source = "estimates"))
      )
    }
  )
  for(invalid_input in invalid_inputs){
    expect_error(
      invalid_input(),
      "'interval_level' is derived output and cannot be supplied",
      fixed = TRUE
    )
  }
})


test_that("interpret_records supports optional missing entries", {

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

  records <- interpret_records(sources, spec)

  expect_equal(nrow(records), 1)
  expect_equal(records$record_id, "moderators.joint.evidence")
  expect_equal(records$BF_canonical_value, 4)

})


test_that("the removed interpret2 and interpret_tables are not exported", {

  exports <- getNamespaceExports("BayesTools")
  expect_true("interpret" %in% exports)
  expect_true("interpret_records" %in% exports)
  expect_false(any(c("interpret2", "interpret_tables") %in% exports))
})


test_that("interpret function input validation works", {

  # Test specification validation
  expect_error(interpret(list(), list(), "not a list", "Test"))

  # Test invalid specification elements
  expect_error(interpret(list(), list(), list(list(inference = 1)), "Test"))

})

test_that("N22 fallback endpoints retain their actual probability labels", {

  source <- data.frame(Mean = .3, "0.025" = .1, "0.975" = .5,
                       check.names = FALSE, row.names = "mu")
  plan <- list(list(kind = "estimate", source = "est", row = "mu",
                    lower_prob = .05, upper_prob = .95, label = "chosen label"))
  records <- interpret_records(list(est = source), plan)
  expect_equal(unlist(records[c("lower_value", "upper_value")], use.names = FALSE), c(.1, .5))
  expect_equal(unlist(records[c("lower_prob", "upper_prob", "interval_level")],
                      use.names = FALSE), c(.025, .975, .95), tolerance = 1e-15)
  expect_identical(records$label, "chosen label")
  text <- interpret_records(list(est = source), plan, output = "text")
  expect_match(paste(text, collapse = "\n"), "95%", fixed = TRUE)
  expect_false(grepl("90%", paste(text, collapse = "\n"), fixed = TRUE))
})

test_that("N22 exact requested and inferred probabilities keep the source interval", {

  source <- data.frame(Mean = .3, "0.025" = .1, "0.975" = .5,
                       check.names = FALSE, row.names = "mu")
  for(probabilities in list(list(lower_prob = .025, upper_prob = .975), list())){
    plan <- list(c(list(kind = "estimate", source = "est", row = "mu"), probabilities))
    records <- interpret_records(list(est = source), plan)
    expect_equal(unlist(records[c("lower_prob", "upper_prob", "interval_level")],
                        use.names = FALSE), c(.025, .975, .95), tolerance = 1e-15)
  }
})

test_that("N73 expanded estimates stay together before the next plan item", {

  source <- data.frame(Mean = c(.3, .4), "0.025" = c(.1, .2), "0.975" = c(.5, .6),
                       check.names = FALSE, row.names = c("a", "b"))
  plan <- list(list(kind = "for_each", source = "est", rows = c("a", "b"),
                    template = "estimate"), list(kind = "note", text = "after"))
  records <- interpret_records(list(est = source), plan)
  expect_identical(records$kind, c("estimate", "estimate", "note"))
  expect_identical(records$row[1:2], c("a", "b"))
  expect_equal(records$order, c(1, 1, 2), tolerance = 0)
  text <- interpret_records(list(est = source), plan, output = "text")
  expect_identical(tail(text, 1L), "after")
})

test_that("N73 paired child records stay adjacent in row order", {

  estimates <- data.frame(Mean = c(.3, .4), "0.025" = c(.1, .2), "0.975" = c(.5, .6),
                          check.names = FALSE, row.names = c("a", "b"))
  evidence <- data.frame(prior_prob = c(.5, .5), post_prob = c(2 / 3, .2),
                         inclusion_BF = c(2, 1 / 4), row.names = c("a", "b"))
  plan <- list(list(kind = "for_each", source = "tests", pair_with = "est",
                    rows = c("a", "b")), list(kind = "note", text = "after"))
  records <- interpret_records(list(tests = evidence, est = estimates), plan)
  expect_identical(records$kind, c("evidence", "estimate", "evidence", "estimate", "note"))
  expect_identical(records$row[1:4], c("a", "a", "b", "b"))
  expect_equal(records$order, c(1, 1, 1, 1, 2), tolerance = 0)
})

test_that("N73 ordinary explicitly ordered plan items remain stable", {

  plan <- list(list(kind = "note", order = 3, text = "last"),
               list(kind = "note", order = 1, text = "first"),
               list(kind = "note", order = 1, text = "second"))
  expect_identical(as.character(interpret_records(list(est = data.frame(Mean = .3, row.names = "mu")), plan, output = "text")),
                   c("first", "second", "last"))
})

test_that("N73 explicit pair reference order overrides retain precedence", {

  estimates <- data.frame(Mean = .3, row.names = "a")
  evidence <- data.frame(prior_prob = .5, post_prob = 2 / 3, inclusion_BF = 2, row.names = "a")
  plan <- list(list(kind = "pair", order = 2,
                    evidence = list(source = "tests", row = "a", order = 5),
                    estimate = list(source = "est", row = "a", order = 1)),
               list(kind = "note", order = 3, text = "between"))
  records <- interpret_records(list(tests = evidence, est = estimates), plan)
  expect_identical(records$kind, c("estimate", "note", "evidence"))
  expect_equal(records$order, c(1, 3, 5), tolerance = 0)
})
