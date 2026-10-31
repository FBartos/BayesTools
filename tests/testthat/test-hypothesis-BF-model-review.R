skip_if_not_test_profile("fixture")

# ============================================================================ #
# TEST FILE: Fit-level hypothesis Bayes factors, contracts on cached fits
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

source(testthat::test_path("common-functions.R"))

skip_if_no_fits()

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
