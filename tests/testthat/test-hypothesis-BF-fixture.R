skip_if_not_test_profile("fixture")

# ============================================================================ #
# TEST FILE: Fit-level hypothesis Bayes factors on cached fits
# ============================================================================ #
#
# PURPOSE:
#   hypothesis_BF() of cached fits: the one-call form against the three-step extraction,
#   contracts of the fit-level engine, and agreement with bridge-sampling model Bayes
#   factors. The synthetic-fit counterparts are in test-hypothesis-BF-model.R.
# ============================================================================ #

source(testthat::test_path("common-functions.R"))

skip_if_no_fits()

# ============================================================================ #
# SECTION: Fit-level hypothesis Bayes factors on cached fits
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


test_that("the fitted-scale catalog of a cached fit keeps its aliases", {

  fit <- .hypothesis_model_fixture_fit("fit_dual_param_regression")
  original <- parameter_catalog(fit)
  fitted   <- parameter_catalog(fit, transform_scaled = FALSE)
  expect_identical(fitted$aliases, original$aliases)
  expect_true(all(c("exp(intercept)", "(log_sigma) exp(intercept)") %in% fitted$aliases$alias))
  for(catalog in list(original, fitted)){
    expect_identical(
      parameter_catalog_resolve(catalog, "exp(intercept)", namespace = "log_sigma")$quantity_id,
      parameter_catalog_resolve(catalog, "(log_sigma) exp(intercept)")$quantity_id
    )
  }
  expect_identical(
    parameter_catalog_resolve(fitted, "exp(intercept)", namespace = "log_sigma")$quantity_id,
    parameter_catalog_resolve(original, "exp(intercept)", namespace = "log_sigma")$quantity_id
  )

})

# ============================================================================ #
# SECTION: Fit-level hypothesis Bayes factors, contracts on cached fits
# ============================================================================ #
#
# PURPOSE:
#   The level token (`component`) of a factor term, the quantities table of a
#   fit without tested entries, the type and cause of the refusal of a
#   random-effect density target, and the request (with the declared
#   conditioning of the draws) that every callback of a density backend
#   receives, on the cached fits of test-00-model-fits.R. The synthetic-fit
#   counterparts are in test-hypothesis-BF-model.R.
#
# MODELS/FIXTURES (test-00-model-fits.R):
#   fit_formula_orthonormal, fit_formula_treatment, fit_vector_mnormal,
#   fit_label_random_slope, fit_re_summary_point_components
#
# TAGS: @fixture, @hypothesis
# ============================================================================ #


.hypothesis_review_fit <- function(name){

  skip_if_missing_fits(name)
  readRDS(file.path(temp_fits_dir, paste0(name, ".RDS")))
}

# An ordinate backend that attaches the exact normal ordinate of the draws with
# the conditioning of the request, and records what its callbacks receive.
.hypothesis_review_backend <- function(state, ...){

  state$requests <- list()
  state$warnings <- list()
  state$diagnostics <- list()
  hypothesis_density_backend(
    name = "reviewNormal", type = "ordinate",
    attach = function(posterior, request, control, cache){
      state$requests[[length(state$requests) + 1L]] <- request
      sample <- if(identical(request$kind, "level")) posterior[[request$level]] else posterior
      draws  <- as.numeric(sample)
      entry  <- do.call(posterior_ordinate_attribute, c(
        list(value = request$values,
             ordinate = stats::dnorm(request$values, mean(draws), stats::sd(draws)),
             method = "review", density_method = "reviewNormal",
             parameter = request$parameter),
        request$condition
      ))
      posterior_metadata(sample, "posterior_ordinate") <- posterior_ordinate_append(
        posterior_metadata(sample, "posterior_ordinate"), entry
      )
      if(identical(request$kind, "level")){
        posterior[[request$level]] <- sample
        posterior
      }else{
        sample
      }
    },
    warnings = function(posterior, request){
      state$warnings[[length(state$warnings) + 1L]] <- request
      BayesTools:::.hypothesis_model_ordinate_warning_records(posterior)
    },
    diagnostics = function(posterior, request){
      state$diagnostics[[length(state$diagnostics) + 1L]] <- request
      NULL
    },
    ...
  )
}


test_that("a level token selects one level of a cached factor term", {

  for(case in list(
    list(fit = "fit_formula_orthonormal", term = "x_fac3o", levels = c("A", "B", "C")),
    list(fit = "fit_formula_treatment", term = "x_fac2t", levels = c("B"))
  )){
    fit <- .hypothesis_review_fit(case$fit)
    for(level in case$levels){
      selected <- hypothesis_BF(fit, paste0(case$term, " = 0.1"), component = level, seed = 1)
      explicit <- hypothesis_BF(fit, paste0(case$term, "[", level, "] = 0.1"), seed = 1)
      expect_identical(nrow(selected), 1L, info = paste(case$fit, level))
      expect_identical(attr(selected, "raw_log_BF"), attr(explicit, "raw_log_BF"),
                       info = paste(case$fit, level))
    }
  }

  # the fixed level of a treatment factor is refused as a statement on it,
  # while the other levels are tested
  fit <- .hypothesis_review_fit("fit_formula_treatment")
  expect_error(hypothesis_BF(fit, "x_fac2t = 0.1", component = "A"),
               class = "BayesTools_hypothesis_fixed")
  expect_error(hypothesis_BF(fit, "x_fac2t = 0.1"), class = "BayesTools_hypothesis_fixed")

  # three levels of a term: the whole term gives three rows, the token one
  fit <- .hypothesis_review_fit("fit_formula_orthonormal")
  expect_identical(nrow(hypothesis_BF(fit, "x_fac3o = 0.1", seed = 1)), 3L)

})


test_that("a fit without tested entries has an empty quantities table", {

  fit <- .hypothesis_review_fit("fit_vector_mnormal")
  expect_identical(nrow(hypothesis_provider(fit)$metadata()$entries), 0L)
  expect_no_warning(quantities <- hypothesis_quantities(fit))
  expect_s3_class(quantities, "data.frame")
  expect_identical(nrow(quantities), 0L)
  expect_identical(names(quantities), names(hypothesis_quantities(
    .hypothesis_review_fit("fit_formula_treatment")
  )))

})


test_that("the refusal of a random-effect density target keeps its type and cause", {

  fit <- .hypothesis_review_fit("fit_label_random_slope")
  model <- structure(list(fit = fit), class = "bt_hypothesis_review_random")
  registerS3method("hypothesis_provider", "bt_hypothesis_review_random", function(object, ...){
    provider <- hypothesis_provider(object$fit)
    provider$object <- object
    provider$random$density_target <- function(parameter, transform_scaled, operation){
      list(reason = paste0("The allocated SD '", parameter, "' has no original-scale ordinate."),
           type = "target", cause = "original_scale")
    }
    provider
  }, envir = asNamespace("BayesTools"))
  backend <- hypothesis_density_backend(
    "reviewRefusal", type = "ordinate",
    attach = function(posterior, request, control, cache) stop("unexpected attach"),
    refusal_classes = function(cause){
      c(paste0("review_density_", cause), "review_density_unavailable")
    }
  )
  register_hypothesis_density_backend(backend)
  withr::defer(unregister_hypothesis_density_backend("reviewRefusal"))

  # the classes of the engine's refusal, then the backend's
  expected <- c("BayesTools_hypothesis_target", "BayesTools_hypothesis_unavailable",
                "review_density_original_scale", "review_density_unavailable")
  for(statement in c("(mu) sd(x) = 0.2", "(mu) var(x) = 0.04")){
    refused <- tryCatch(
      hypothesis_BF_model(model, statement, density_method = "reviewRefusal"),
      error = identity
    )
    expect_identical(class(refused)[seq_along(expected)], expected, info = statement)
    expect_match(conditionMessage(refused), "has no original-scale ordinate", info = statement)
  }
  # the SD of a variance names the refusal: the reason is the variance's own
  variance <- tryCatch(
    hypothesis_BF_model(model, "(mu) var(x) = 0.04", density_method = "reviewRefusal"),
    error = identity
  )
  expect_match(conditionMessage(variance), "'(mu) var(x)'", fixed = TRUE)
  # region statements need no density method
  expect_s3_class(
    hypothesis_BF_model(model, "(mu) sd(x) > 0.1", density_method = "reviewRefusal", seed = 1),
    "BayesTools_hypothesis_BF"
  )

})


test_that("the callbacks of a backend receive the declared conditioning of the draws", {

  state   <- new.env()
  backend <- .hypothesis_review_backend(state)
  register_hypothesis_density_backend(backend)
  withr::defer(unregister_hypothesis_density_backend("reviewNormal"))

  # conditional draws: the request carries the declared conditioning
  fit <- .hypothesis_review_fit("fit_re_summary_point_components")
  result <- hypothesis_BF(fit, "x = 0.1", conditional = TRUE,
                          density_method = "reviewNormal", seed = 1)
  samples <- as_mixed_posteriors(fit, "mu_x", conditional = "mu_x",
                                 conditional_rule = "AND", n_prior_samples = 100)
  declared <- posterior_metadata(
    marginal_posterior(samples, "mu_x", prior_samples = TRUE, use_formula = FALSE,
                       n_samples = 100),
    "condition"
  )
  request <- state$requests[[1L]]
  expect_identical(request$conditional, "mu_x")
  expect_identical(request$condition$conditional, "mu_x")
  expect_identical(request$condition$conditional_rule, "AND")
  expect_identical(request$condition$condition_key, declared$condition_key)
  expect_s3_class(request$condition$condition_event, "BayesTools_condition_event")
  expect_s3_class(request$condition$resolved_condition_event, "BayesTools_condition_event")
  expect_identical(request$condition$condition_key, "AND\r1\rmu_x")
  # the same request for the warnings and the diagnostics of the group
  expect_identical(state$warnings[[1L]], request)
  expect_identical(state$diagnostics[[1L]], request)
  expect_s3_class(result, "BayesTools_hypothesis_BF")

  # the draws that carry a stored ordinate and their conditioning reach reuse()
  posterior <- marginal_posterior(samples, "mu_x", prior_samples = TRUE,
                                  use_formula = FALSE, n_samples = 100)
  posterior_metadata(posterior, "posterior_ordinate") <- do.call(
    posterior_ordinate_attribute,
    c(list(value = 0.1, ordinate = 1, method = "stored", density_method = "reviewNormal",
           parameter = "mu_x"), request$condition)
  )
  seen <- NULL
  stored <- hypothesis_density_backend(
    "reviewNormal", type = "ordinate",
    attach = function(posterior, request, control, cache) stop("unexpected attach"),
    reuse = function(entry, request, control, posterior){
      seen <<- list(entry = entry, request = request, posterior = posterior)
      identical(entry$condition_key, request$condition$condition_key) &&
        identical(posterior_metadata(posterior, "condition")$condition_key,
                  request$condition$condition_key)
    }
  )
  plan  <- list(kind = "scalar")
  stored_request <- BayesTools:::.hypothesis_model_request(
    hypothesis_provider(fit), stored, plan, kind = "scalar", parameter = "mu_x",
    label = "x", values = 0.1, sample = posterior, conditional = "mu_x"
  )
  attached <- BayesTools:::.hypothesis_model_attach_request(
    posterior, stored, stored_request, NULL, new.env()
  )
  expect_identical(seen$request$condition, request$condition)
  expect_identical(seen$posterior, posterior)
  expect_identical(
    BayesTools:::.hypothesis_model_ordinate_entries(
      posterior_metadata(attached, "posterior_ordinate"))[[1L]]$method,
    "stored"
  )

  # random-effect quantities: the SD, and the variance through its SD
  fit <- .hypothesis_review_fit("fit_label_random_slope")
  for(case in list(c("(mu) sd(x) = 0.2", "(mu) sd(x)"), c("(mu) var(x) = 0.04", "(mu) sd(x)"))){
    state$requests <- list(); state$warnings <- list(); state$diagnostics <- list()
    hypothesis_BF(fit, case[[1L]], density_method = "reviewNormal", seed = 1)
    request <- state$requests[[1L]]
    expect_identical(request$kind, "scalar", info = case[[1L]])
    expect_identical(request$route, "random", info = case[[1L]])
    expect_identical(request$parameter, case[[2L]], info = case[[1L]])
    expect_true(is.list(request$condition), info = case[[1L]])
    expect_identical(state$warnings[[1L]], request, info = case[[1L]])
    expect_identical(state$diagnostics[[1L]], request, info = case[[1L]])
  }

})

# ============================================================================ #
# SECTION: Hypothesis Bayes Factor Bridge Comparisons
# ============================================================================ #


.hypothesis_treatment_level_marginal_for_test <- function(fit){

  model <- coda::as.mcmc(as.matrix(fit[["mcmc"]]))
  class(model) <- c("BayesTools_fit", class(model))
  attr(model, "prior_list") <- attr(fit, "prior_list", exact = TRUE)
  model <- attach_test_parameter_map(model)

  mixed <- as_mixed_posteriors(
    model,
    parameters      = "mu_x_fac2t",
    n_prior_samples = 20000
  )

  marginal_posterior(
    samples       = mixed,
    parameter     = "mu_x_fac2t",
    prior_samples = TRUE,
    use_formula   = FALSE,
    n_samples     = 20000
  )
}


# Log posterior odds of `parameter > 0` in the fitted draws and their
# delta-method Monte Carlo standard error, 1 / sqrt(ESS p (1 - p)), from the
# effective sample size of the region indicator across chains. For the
# treatment contrast, level B minus level A is the coefficient itself.
.hypothesis_draws_log_odds_for_test <- function(fit, parameter){

  indicator <- coda::mcmc.list(lapply(fit[["mcmc"]], function(chain){
    coda::mcmc(as.numeric(chain[, parameter] > 0))
  }))
  p   <- mean(unlist(indicator))
  ess <- unname(coda::effectiveSize(indicator))

  list(
    log_odds = log(p / (1 - p)),
    mcse     = 1 / sqrt(ess * p * (1 - p))
  )
}


# Checks an inequality hypothesis BF against the ratio of the truncated-prior
# bridge marginal likelihoods (the encompassing-prior identity; the symmetric
# N(0, 1) treatment prior has prior odds 1). The fixture's log BF is small
# (0.22), so the direction is asserted explicitly: a reciprocal BF or BF = 1
# fails. The magnitude must agree within 3 Monte Carlo standard errors of the
# difference: the bridge MCSEs of both marginal likelihoods and the posterior
# draw error of the encompassing fit (together about 0.09 on the log scale).
# The hypothesis BF must also reproduce the log posterior odds of the same
# draws within 0.08, four standard errors of resampling 20000 posterior and
# 20000 prior draws (sqrt(1 / (20000 p (1 - p)) + 1 / (20000 / 4)) = 0.02).
.expect_inequality_BF_matches_bridge <- function(bf_hyp, fit, ml_positive, ml_negative){

  log_BF_hypothesis <- log(attr(bf_hyp, "raw_BF"))
  log_BF_bridge     <- ml_positive[["logml"]] - ml_negative[["logml"]]
  draws             <- .hypothesis_draws_log_odds_for_test(fit, "mu_x_fac2t")
  bridge_mcse       <- sqrt(sum(c(
    ml_positive[["repetitions"]][["mcse"]],
    ml_negative[["repetitions"]][["mcse"]]
  )^2))
  difference_mcse   <- sqrt(draws$mcse^2 + bridge_mcse^2)

  expect_gt(log_BF_bridge, 0)
  expect_gt(log_BF_hypothesis, 0)
  expect_lt(abs(log_BF_hypothesis - log_BF_bridge), 3 * difference_mcse)
  expect_lt(abs(log_BF_hypothesis - draws$log_odds), 0.08)
}


test_that("hypothesis_BF point-null agrees with bridge-sampling model BF", {

  fit_alt <- readRDS(file.path(temp_fits_dir, "fit_simple_normal.RDS"))
  fit_null <- readRDS(file.path(temp_fits_dir, "fit_simple_spike.RDS"))
  ml_alt <- readRDS(file.path(temp_marglik_dir, "fit_simple_normal.RDS"))
  ml_null <- readRDS(file.path(temp_marglik_dir, "fit_simple_spike.RDS"))

  post <- as.matrix(fit_alt[["mcmc"]])[, "m"]
  prior <- attr(fit_alt, "prior_list")[["m"]]

  bf_hyp <- hypothesis_BF(
    posterior  = post,
    prior      = prior,
    hypothesis = "m == 0"
  )
  bf_bridge <- exp(ml_alt[["logml"]] - ml_null[["logml"]])

  expect_equal(log(attr(bf_hyp, "raw_BF")), log(bf_bridge), tolerance = 0.25)
  expect_s3_class(fit_null, "BayesTools_fit")
})


test_that("hypothesis_BF explicit level point-null agrees with bridge-sampling model BF", {

  fit_alt <- readRDS(file.path(temp_fits_dir, "fit_formula_treatment.RDS"))
  ml_alt <- readRDS(file.path(temp_marglik_dir, "fit_formula_treatment.RDS"))
  ml_null <- readRDS(file.path(temp_marglik_dir, "fit_formula_simple.RDS"))

  posterior <- .hypothesis_treatment_level_marginal_for_test(fit_alt)

  bf_hyp <- hypothesis_BF(
    posterior  = posterior,
    hypothesis = "mu_x_fac2t[B] = 0",
    seed       = 101
  )
  bf_bridge <- exp(ml_alt[["logml"]] - ml_null[["logml"]])

  expect_equal(log(attr(bf_hyp, "raw_BF")), log(bf_bridge), tolerance = 0.35)
})


test_that("hypothesis_BF treatment-level point contrast agrees with bridge-sampling model BF", {

  fit_alt <- readRDS(file.path(temp_fits_dir, "fit_formula_treatment.RDS"))
  ml_alt <- readRDS(file.path(temp_marglik_dir, "fit_formula_treatment.RDS"))
  ml_null <- readRDS(file.path(temp_marglik_dir, "fit_formula_simple.RDS"))

  posterior <- .hypothesis_treatment_level_marginal_for_test(fit_alt)

  bf_hyp <- hypothesis_BF(
    posterior  = posterior,
    hypothesis = "mu_x_fac2t[B] - mu_x_fac2t[A] = 0",
    seed       = 102
  )
  bf_bridge <- exp(ml_alt[["logml"]] - ml_null[["logml"]])

  expect_equal(log(attr(bf_hyp, "raw_BF")), log(bf_bridge), tolerance = 0.45)
})


test_that("hypothesis_BF treatment-level inequality agrees with truncated-prior bridge BF", {

  fit_alt <- readRDS(file.path(temp_fits_dir, "fit_formula_treatment.RDS"))
  ml_positive <- readRDS(file.path(temp_marglik_dir, "fit_formula_treatment_positive.RDS"))
  ml_negative <- readRDS(file.path(temp_marglik_dir, "fit_formula_treatment_negative.RDS"))

  posterior <- .hypothesis_treatment_level_marginal_for_test(fit_alt)

  bf_hyp <- hypothesis_BF(
    posterior  = posterior,
    hypothesis = "mu_x_fac2t[B] > mu_x_fac2t[A] vs mu_x_fac2t[B] < mu_x_fac2t[A]",
    seed       = 103
  )

  .expect_inequality_BF_matches_bridge(bf_hyp, fit_alt, ml_positive, ml_negative)
})


test_that("hypothesis_BF transformed level inequality agrees with truncated-prior bridge BF", {

  fit_alt <- readRDS(file.path(temp_fits_dir, "fit_formula_treatment.RDS"))
  ml_positive <- readRDS(file.path(temp_marglik_dir, "fit_formula_treatment_positive.RDS"))
  ml_negative <- readRDS(file.path(temp_marglik_dir, "fit_formula_treatment_negative.RDS"))

  posterior <- .hypothesis_treatment_level_marginal_for_test(fit_alt)

  bf_hyp <- hypothesis_BF(
    posterior  = posterior,
    hypothesis = paste(
      "exp(mu_x_fac2t[B]) > exp(mu_x_fac2t[A])",
      "vs",
      "exp(mu_x_fac2t[B]) < exp(mu_x_fac2t[A])"
    ),
    seed       = 104
  )

  .expect_inequality_BF_matches_bridge(bf_hyp, fit_alt, ml_positive, ml_negative)
})


test_that("hypothesis_BF point-vs-level-region agrees with truncated-prior bridge BF", {

  fit_alt <- readRDS(file.path(temp_fits_dir, "fit_formula_treatment.RDS"))
  ml_null <- readRDS(file.path(temp_marglik_dir, "fit_formula_simple.RDS"))
  ml_positive <- readRDS(file.path(temp_marglik_dir, "fit_formula_treatment_positive.RDS"))

  posterior <- .hypothesis_treatment_level_marginal_for_test(fit_alt)

  bf_hyp <- hypothesis_BF(
    posterior  = posterior,
    hypothesis = paste(
      "mu_x_fac2t[B] - mu_x_fac2t[A] = 0",
      "vs",
      "mu_x_fac2t[B] - mu_x_fac2t[A] > 0"
    ),
    seed       = 105
  )
  bf_bridge <- exp(ml_null[["logml"]] - ml_positive[["logml"]])

  expect_equal(log(attr(bf_hyp, "raw_BF")), log(bf_bridge), tolerance = 0.55)
})


test_that("hypothesis_BF nested coefficient agrees with bridge-sampling model BF", {

  fit_alt <- readRDS(file.path(temp_fits_dir, "fit_formula_treatment.RDS"))
  ml_alt <- readRDS(file.path(temp_marglik_dir, "fit_formula_treatment.RDS"))
  ml_null <- readRDS(file.path(temp_marglik_dir, "fit_formula_simple.RDS"))

  post <- as.matrix(fit_alt[["mcmc"]])[, "mu_x_fac2t"]
  prior <- attr(fit_alt, "prior_list")[["mu_x_fac2t"]]

  bf_hyp <- hypothesis_BF(
    posterior  = post,
    prior      = prior,
    hypothesis = "mu_x_fac2t == 0"
  )
  bf_bridge <- exp(ml_alt[["logml"]] - ml_null[["logml"]])

  expect_equal(log(attr(bf_hyp, "raw_BF")), log(bf_bridge), tolerance = 0.35)
})
