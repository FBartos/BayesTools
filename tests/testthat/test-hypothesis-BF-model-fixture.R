skip_if_not_test_profile("fixture")

# ============================================================================ #
# TEST FILE: Fit-level hypothesis Bayes factors on cached fits
# ============================================================================ #
#
# PURPOSE:
#   The one-call form hypothesis_BF(fit, "(mu) cor = 0.9") equals the
#   three-step extraction (parameter_catalog_resolve(), parameter_mixed_posterior(),
#   hypothesis_BF() of the draws), also on the original and fitted scale of
#   scaled predictors, for spike-and-slab coefficients, and for the quantities
#   of fits with several formula parameters.
#
# MODELS/FIXTURES (test-00-model-fits.R):
#   fit_hypothesis_hcs, fit_hypothesis_car, fit_hypothesis_spike_ar1,
#   fit_label_random_slope, fit_re_summary_point_components,
#   fit_dual_param_regression, fit_formula_treatment
#
# TAGS: @fixture, @hypothesis
# ============================================================================ #

source(testthat::test_path("common-functions.R"))

skip_if_no_fits()

.hypothesis_model_fixture_fit <- function(name){

  skip_if_missing_fits(name)
  readRDS(file.path(temp_fits_dir, paste0(name, ".RDS")))
}

# The three-step extraction of a catalog quantity.
.hypothesis_model_three_step <- function(fit, name, statement, ...){

  draws <- parameter_mixed_posterior(
    fit, parameter_catalog_resolve(parameter_catalog(fit), name), ...
  )
  hypothesis_BF(draws, hypothesis = statement)
}


test_that("structured correlations equal the three-step extraction", {

  for(case in list(
    list(fit = "fit_hypothesis_hcs", statement = "(mu) cor = 0.9"),
    list(fit = "fit_hypothesis_car", statement = "(mu) cor = 0.5")
  )){
    fit <- .hypothesis_model_fixture_fit(case$fit)
    one_call <- hypothesis_BF(fit, case$statement, BF01 = TRUE, seed = 1)
    quoted   <- sub("^\\(mu\\) cor", "`(mu) cor`", case$statement)
    manual   <- hypothesis_BF(
      parameter_mixed_posterior(
        fit, parameter_catalog_resolve(parameter_catalog(fit), "(mu) cor")
      ),
      hypothesis = quoted, BF01 = TRUE, seed = 1
    )
    expect_identical(attr(one_call, "raw_log_BF"), attr(manual, "raw_log_BF"), info = case$fit)
    expect_identical(as.numeric(one_call$BF), as.numeric(manual$BF), info = case$fit)
    expect_identical(rownames(one_call), rownames(manual))
    expect_true(is.finite(attr(one_call, "raw_log_BF")))
    expect_true(attr(one_call, "BF01"))

    # region statements and vs statements of the same quantity
    both <- hypothesis_BF(fit, c(case$statement, "(mu) cor > 0.3 vs (mu) cor < 0.3"),
                          seed = 1)
    expect_identical(rownames(both), c("(mu) cor (1)", "(mu) cor (2)"))
    expect_identical(attr(both, "raw_log_BF")[[1L]],
                     attr(hypothesis_BF(fit, case$statement, seed = 1), "raw_log_BF"))
  }

})


test_that("transform_scaled tests the fitted or the original scale", {

  fit   <- .hypothesis_model_fixture_fit("fit_label_random_slope")
  route <- JAGS_formula_coefficient_route(fit, "mu", "mu_x")
  expect_identical(route$type, "affine")
  scale <- unname(route$weights[["mu_x"]])

  # a coefficient: the value on the original scale is the fitted value times
  # the weight, and the Savage-Dickey ratio does not depend on the scale
  original <- hypothesis_BF(fit, "x = 0.05", transform_scaled = TRUE, seed = 1)
  fitted   <- hypothesis_BF(fit, paste0("x = ", format(0.05 / scale, digits = 17)), seed = 1)
  expect_equal(attr(original, "raw_log_BF"), attr(fitted, "raw_log_BF"), tolerance = 1e-6)

  # a random-effect SD: the draws of the two scales differ by the same weight,
  # the default is the fitted scale, and the one-call form equals the
  # three-step extraction on either
  selection <- parameter_catalog_resolve(parameter_catalog(fit), "(mu) sd(x)")
  on_fitted   <- parameter_mixed_posterior(fit, selection, transform_scaled = FALSE)
  on_original <- parameter_mixed_posterior(fit, selection)
  expect_identical(as.numeric(on_original),
                   as.numeric(parameter_mixed_posterior(fit, selection, transform_scaled = TRUE)))
  expect_equal(as.numeric(on_original) / as.numeric(on_fitted),
               rep(scale, length(on_fitted)), tolerance = 1e-10)
  estimates <- function(transform_scaled){
    JAGS_estimates_table(fit, transform_scaled = transform_scaled)["(mu) sd(x)", "Mean"]
  }
  expect_equal(mean(as.numeric(on_fitted)), as.numeric(estimates(FALSE)), tolerance = 1e-10)
  expect_equal(mean(as.numeric(on_original)), as.numeric(estimates(TRUE)), tolerance = 1e-10)

  for(transform_scaled in c(FALSE, TRUE)){
    one_call <- hypothesis_BF(fit, "(mu) sd(x) = 0.2", transform_scaled = transform_scaled, seed = 1)
    manual   <- hypothesis_BF(
      parameter_mixed_posterior(fit, selection, transform_scaled = transform_scaled),
      hypothesis = "`(mu) sd(x)` = 0.2", seed = 1
    )
    expect_identical(attr(one_call, "raw_log_BF"), attr(manual, "raw_log_BF"),
                     info = transform_scaled)
  }
  fitted_default <- hypothesis_BF(fit, "(mu) sd(x) = 0.2", seed = 1)
  expect_identical(attr(fitted_default, "raw_log_BF"),
                   attr(hypothesis_BF(fit, "(mu) sd(x) = 0.2", transform_scaled = FALSE, seed = 1),
                        "raw_log_BF"))
  expect_false(isTRUE(all.equal(
    attr(fitted_default, "raw_log_BF"),
    attr(hypothesis_BF(fit, "(mu) sd(x) = 0.2", transform_scaled = TRUE, seed = 1), "raw_log_BF")
  )))

  # a variance is tested through its standard deviation
  variance <- hypothesis_BF(fit, "(mu) var(x) = 0.04", seed = 1)
  expect_equal(attr(variance, "raw_log_BF"),
               attr(hypothesis_BF(fit, "(mu) sd(x) = 0.2", seed = 1), "raw_log_BF"),
               tolerance = 1e-8)

})


test_that("spike-and-slab coefficients are tested within their inclusion event", {

  fit <- .hypothesis_model_fixture_fit("fit_re_summary_point_components")

  expect_error(hypothesis_BF(fit, "x = 0"), class = "BayesTools_point_mass_at_null")
  expect_error(hypothesis_BF(fit, "(mu) sd(intercept) = 0"),
               class = "BayesTools_point_mass_at_null")
  indicator <- tryCatch(hypothesis_BF(fit, "mu_x_indicator = 0"), error = identity)
  expect_s3_class(indicator, "BayesTools_hypothesis_target")
  expect_s3_class(indicator, "BayesTools_hypothesis_unavailable")

  # conditional = TRUE: Savage-Dickey within the inclusion event
  samples <- as_mixed_posteriors(fit, "mu_x", conditional = "mu_x",
                                 conditional_rule = "AND", n_prior_samples = 10000)
  posterior <- marginal_posterior(samples, "mu_x", prior_samples = TRUE,
                                  use_formula = FALSE, n_samples = 10000)
  conditional <- hypothesis_BF(fit, c("x = 0", "x = 0.1"), conditional = TRUE, seed = 1)
  manual <- hypothesis_BF.default(posterior, hypothesis = c("mu_x = 0", "mu_x = 0.1"),
                                  parameter = "mu_x", seed = 1)
  expect_identical(attr(conditional, "raw_log_BF"), attr(manual, "raw_log_BF"))
  expect_identical(attr(conditional, "type"), attr(manual, "type"))

  # without conditional, the test is within the full ensemble and warns once
  expect_warning(
    full <- hypothesis_BF(fit, "x = 0.1", seed = 1),
    "Model-averaged parameter test: this hypothesis uses the full ensemble for mu_x"
  )
  expect_no_warning(hypothesis_BF(fit, "x = 0.1", conditional = FALSE, seed = 1))
  expect_false(isTRUE(all.equal(attr(full, "raw_log_BF")[[1L]],
                                attr(conditional, "raw_log_BF")[[2L]])))
  region <- suppressWarnings(hypothesis_BF(fit, "x > 0", seed = 1))
  samples <- as_mixed_posteriors(fit, "mu_x", n_prior_samples = 10000)
  posterior <- marginal_posterior(samples, "mu_x", prior_samples = TRUE,
                                  use_formula = FALSE, n_samples = 10000)
  expect_identical(
    attr(region, "raw_log_BF"),
    attr(hypothesis_BF.default(posterior, hypothesis = "mu_x > 0", parameter = "mu_x",
                               seed = 1), "raw_log_BF")
  )

})


test_that("spike-and-slab random-effect SDs are tested within their inclusion event", {

  fit <- .hypothesis_model_fixture_fit("fit_re_summary_point_components")
  catalog <- parameter_catalog(fit)
  three_step <- function(name, statement, ...){
    draws <- parameter_mixed_posterior(
      fit, parameter_catalog_resolve(catalog, name), ...
    )
    hypothesis_BF(draws, hypothesis = statement, seed = 1)
  }

  # conditional = TRUE equals the three-step extraction exactly: point values,
  # regions, and the variance (tested through its SD, or on the variance draws)
  for(case in list(
    list(statement = "(mu) sd(intercept) = 0.5", name = "(mu) sd(intercept)",
         manual = "`(mu) sd(intercept)` = 0.5"),
    list(statement = "(mu) sd(intercept) > 0.1", name = "(mu) sd(intercept)",
         manual = "`(mu) sd(intercept)` > 0.1"),
    list(statement = "(mu) var(intercept) = 0.25", name = "(mu) sd(intercept)",
         manual = "`(mu) sd(intercept)` = 0.5"),
    list(statement = "(mu) var(intercept) > 0.04", name = "(mu) var(intercept)",
         manual = "`(mu) var(intercept)` > 0.04")
  )){
    one_call <- hypothesis_BF(fit, case$statement, conditional = TRUE, seed = 1)
    manual   <- three_step(case$name, case$manual, conditional = TRUE)
    expect_identical(attr(one_call, "raw_log_BF"), attr(manual, "raw_log_BF"),
                     info = case$statement)
    expect_true(is.finite(attr(one_call, "raw_log_BF")), info = case$statement)
  }

  # a point null at the spike stops without conditional = TRUE; region tests keep
  # working on the full ensemble
  expect_error(hypothesis_BF(fit, "(mu) sd(intercept) = 0"),
               "or 'conditional = TRUE' to test the quantity within its inclusion event",
               fixed = TRUE, class = "BayesTools_point_mass_at_null")
  expect_error(hypothesis_BF(fit, "(mu) var(intercept) = 0", conditional = FALSE),
               class = "BayesTools_point_mass_at_null")
  region <- hypothesis_BF(fit, "(mu) sd(intercept) > 0.1", conditional = FALSE, seed = 1)
  expect_identical(
    attr(region, "raw_log_BF"),
    attr(three_step("(mu) sd(intercept)", "`(mu) sd(intercept)` > 0.1"), "raw_log_BF")
  )
  expect_false(isTRUE(all.equal(
    attr(region, "raw_log_BF"),
    attr(hypothesis_BF(fit, "(mu) sd(intercept) > 0.1", conditional = TRUE, seed = 1),
         "raw_log_BF")
  )))
})


test_that("factor levels and several formula parameters resolve by their names", {

  fit <- .hypothesis_model_fixture_fit("fit_formula_treatment")
  samples <- as_mixed_posteriors(fit, "mu_x_fac2t", n_prior_samples = 10000)
  levels  <- marginal_posterior(samples, "mu_x_fac2t", prior_samples = TRUE,
                                use_formula = FALSE, n_samples = 10000)
  engine <- hypothesis_BF(fit, "x_fac2t[B] = 0.1", seed = 1)
  manual <- hypothesis_BF.default(levels, hypothesis = "mu_x_fac2t[B] = 0.1",
                                  parameter = "mu_x_fac2t", seed = 1)
  expect_identical(attr(engine, "raw_log_BF"), attr(manual, "raw_log_BF"))
  expect_error(hypothesis_BF(fit, "x_fac2t[A] = 0"), class = "BayesTools_hypothesis_fixed")
  expect_error(hypothesis_BF(fit, "x_fac2t[Z] = 0"), class = "BayesTools_parameter_not_found")

  fit <- .hypothesis_model_fixture_fit("fit_dual_param_regression")
  ambiguous <- tryCatch(hypothesis_BF(fit, "intercept > 0.5"), error = identity)
  expect_s3_class(ambiguous, "BayesTools_hypothesis_ambiguous")
  expect_s3_class(ambiguous, "BayesTools_hypothesis_statement")
  expect_match(conditionMessage(ambiguous), "Set 'namespace' to 'mu', 'log_sigma'", fixed = TRUE)
  for(namespace in c("mu", "log_sigma")){
    selected <- suppressWarnings(
      hypothesis_BF(fit, "intercept > 0.5", namespace = namespace, seed = 1)
    )
    explicit <- suppressWarnings(
      hypothesis_BF(fit, paste0(namespace, "_intercept > 0.5"), seed = 1)
    )
    expect_identical(attr(selected, "raw_log_BF"), attr(explicit, "raw_log_BF"), info = namespace)
    expect_identical(rownames(selected), rownames(explicit), info = namespace)
  }
  expect_error(hypothesis_BF(fit, "mu_intercept > 0.5", namespace = "log_sigma"),
               class = "BayesTools_hypothesis_component_mismatch")
  expect_error(hypothesis_BF(fit, "mu_intercept > 0.5", namespace = "nothing"), "should be one of")
})


test_that("spike-and-slab structured correlations of a cached fit are tested within their inclusion event", {

  fit <- .hypothesis_model_fixture_fit("fit_hypothesis_spike_ar1")
  three_step <- function(name, statement, ...){
    draws <- parameter_mixed_posterior(
      fit, parameter_catalog_resolve(parameter_catalog(fit), name), ...
    )
    hypothesis_BF(draws, hypothesis = statement, seed = 1)
  }
  for(case in list(
    list(statement = "(mu) cor = 0.2", name = "(mu) cor", manual = "`(mu) cor` = 0.2"),
    list(statement = "(mu) cor > 0.1", name = "(mu) cor", manual = "`(mu) cor` > 0.1"),
    list(statement = "(mu) sd = 0.5", name = "(mu) sd", manual = "`(mu) sd` = 0.5"),
    list(statement = "(mu) var > 0.04", name = "(mu) var", manual = "`(mu) var` > 0.04")
  )){
    one_call <- hypothesis_BF(fit, case$statement, conditional = TRUE, seed = 1)
    manual   <- three_step(case$name, case$manual, conditional = TRUE)
    expect_identical(attr(one_call, "raw_log_BF"), attr(manual, "raw_log_BF"),
                     info = case$statement)
    expect_true(is.finite(attr(one_call, "raw_log_BF")), info = case$statement)
  }

  # a point null at a spike stops without conditional = TRUE; region tests keep
  # working on the full ensemble
  expect_error(hypothesis_BF(fit, "(mu) cor = 0", conditional = FALSE),
               "or 'conditional = TRUE' to test the quantity within its inclusion event",
               fixed = TRUE, class = "BayesTools_point_mass_at_null")
  expect_error(hypothesis_BF(fit, "(mu) sd = 0", conditional = FALSE),
               class = "BayesTools_point_mass_at_null")
  region <- hypothesis_BF(fit, "(mu) cor > 0.1", conditional = FALSE, seed = 1)
  expect_identical(attr(region, "raw_log_BF"),
                   attr(three_step("(mu) cor", "`(mu) cor` > 0.1"), "raw_log_BF"))
})


test_that("the random-effect callbacks of the default provider read the quantity on the requested scale", {

  # the fitted-scale quantity of a standardized block is a direct coordinate,
  # the original-scale quantity of the intercept SD is a composite
  fit <- .hypothesis_model_fixture_fit("fit_label_random_slope")
  provider <- hypothesis_provider(fit)
  seen <- list()
  local_mocked_bindings(
    .hypothesis_model_default_sd_pair = function(fit, selection){
      seen[[length(seen) + 1L]] <<- list(
        scales = length(attr(fit, "formula_scale", exact = TRUE)),
        source = selection$quantities$source_type
      )
      NULL
    },
    .package = "BayesTools"
  )
  for(transform_scaled in c(FALSE, TRUE)){
    selected <- provider$random$select("(mu) sd(intercept)", transform_scaled)
    provider$random$sd_pair(selected)
  }
  expect_identical(seen[[1L]]$source, "identity")
  expect_identical(seen[[1L]]$scales, 0L)
  expect_identical(seen[[2L]]$source, "composite")
  expect_gt(seen[[2L]]$scales, 0L)

  # the eligible draws of a conditional quantity are those in the inclusion
  # event of the gates of the fitted-scale quantity
  fit <- .hypothesis_model_fixture_fit("fit_re_summary_allocated")
  provider <- hypothesis_provider(fit)
  parameter <- "(mu) s: sd(intercept)"
  selected  <- provider$random$select(parameter, FALSE)
  expect_null(parameter_gate_states(fit, selected$entry$selection)$event)
  conditioned <- parameter_mixed_posterior(
    fit, selected$entry$selection, conditional = TRUE, transform_scaled = FALSE
  )
  partial <- as.numeric(conditioned)[-1L]
  posterior_metadata(partial, "undefined_draws") <- stats::setNames("test", parameter)
  expect_equal(
    provider$random$defined_share(partial, selected, conditional = TRUE),
    list(defined = length(conditioned) - 1L, eligible = length(conditioned))
  )

})


test_that("a selection of the fitted-scale catalog is extracted with the gates of the fitted scale", {

  fit  <- .hypothesis_model_fixture_fit("fit_re_summary_allocated")
  name <- "(mu) s: sd(intercept)"
  original <- parameter_catalog_resolve(parameter_catalog(fit), name)
  fitted   <- parameter_catalog_resolve(parameter_catalog(fit, transform_scaled = FALSE), name)

  # the original-scale SD of the standardized intercept column combines
  # several fitted coordinates: its support has no declared atoms and no gates
  original_support <- posterior_metadata(parameter_draws(fit, original), "support")[[1L]]
  expect_false(original_support$exact)
  expect_null(parameter_gate_states(fit, original)$event)
  # the fitted-scale SD is the allocation-derived coordinate, which is 0 where its gates are off
  fitted_support <- posterior_metadata(parameter_draws(fit, fitted), "support")[[1L]]
  expect_true(fitted_support$exact)
  expect_identical(fitted_support$type, "mixed")
  expect_identical(fitted_support$points, 0)
  gates <- parameter_gate_states(fit, fitted)
  expect_true(any(gates$event))
  expect_false(all(gates$event))

  # the conditional mixed posterior of the fitted scale is the draws in the event of those gates
  conditioned <- parameter_mixed_posterior(fit, fitted, conditional = TRUE, transform_scaled = FALSE)
  expect_identical(length(conditioned), sum(gates$event))
  expect_identical(
    conditioned,
    parameter_mixed_posterior(fit, original, conditional = TRUE, transform_scaled = FALSE)
  )
  expect_error(
    parameter_mixed_posterior(fit, fitted, conditional = TRUE, transform_scaled = TRUE),
    "belongs to the catalog of the fitted scale", fixed = TRUE
  )

})
