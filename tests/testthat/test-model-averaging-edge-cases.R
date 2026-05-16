skip_if_not_test_profile("unit")

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

.expected_post_probs <- function(prior_weights, margliks) {
  prior_probs <- prior_weights / sum(prior_weights)
  log_weights <- log(prior_probs) + margliks
  log_weights[prior_probs == 0 | is.na(log_weights)] <- -Inf
  log_weights <- log_weights - max(log_weights[is.finite(log_weights)])
  weights <- exp(log_weights)
  weights / sum(weights)
}

.expected_inclusion_bf <- function(prior_probs, post_probs, is_null) {
  (sum(post_probs[!is_null]) / sum(post_probs[is_null])) /
    (sum(prior_probs[!is_null]) / sum(prior_probs[is_null]))
}

.mock_bridge <- function(logml) {
  structure(list(logml = logml), class = "bridge")
}

.mock_runjags_fit_for_mixing <- function(samples, prior_list) {
  samples <- coda::mcmc(as.matrix(samples))
  fit <- structure(
    list(
      mcmc = coda::mcmc.list(samples),
      sample = nrow(samples),
      summary.pars = list(mutate = NULL),
      monitor = colnames(samples)
    ),
    class = c("runjags", "BayesTools_fit", "list")
  )
  attr(fit, "prior_list") <- prior_list
  fit
}

.mock_mixing_model <- function(offset, logml, prior_weight = 1) {
  posterior <- cbind(
    theta = offset + seq_len(20),
    beta  = offset + 100 + seq_len(20)
  )
  list(
    fit = .mock_runjags_fit_for_mixing(
      posterior,
      prior_list = list(
        theta = prior("normal", list(0, 1)),
        beta  = prior("normal", list(0, 1))
      )
    ),
    marglik = .mock_bridge(logml),
    prior_weights = prior_weight
  )
}


# ============================================================================ #
# SECTION 1: inclusion_BF boundary conditions
# ============================================================================ #
test_that("compute_inference matches exact posterior probability algebra", {

  prior_weights <- c(null = 2, alt_wide = 3, alt_narrow = 5)
  margliks <- log(c(null = 0.25, alt_wide = 1.5, alt_narrow = 0.5))
  is_null <- c(TRUE, FALSE, FALSE)

  expected_prior <- prior_weights / sum(prior_weights)
  expected_post <- .expected_post_probs(prior_weights, margliks)
  expected_BF <- .expected_inclusion_bf(expected_prior, expected_post, is_null)

  inference <- compute_inference(prior_weights, margliks, is_null = is_null)

  expect_equal(unname(inference$prior_probs), unname(expected_prior), tolerance = 1e-12)
  expect_equal(unname(inference$post_probs), unname(expected_post), tolerance = 1e-12)
  expect_equal(inference$BF, expected_BF, tolerance = 1e-12)
  expect_equal(sum(inference$prior_probs), 1, tolerance = 1e-12)
  expect_equal(sum(inference$post_probs), 1, tolerance = 1e-12)
  expect_identical(attr(inference, "is_null"), unname(is_null))
  expect_false(attr(inference, "conditional"))
})

test_that("model averaging normalizes explicit unnormalized prior probabilities", {

  prior_weights <- c(null = 2, alt_strong = 3, alt_regular = 5)
  margliks <- log(c(null = 1 / 2, alt_strong = 2, alt_regular = 1))
  is_null <- c(TRUE, FALSE, FALSE)

  expected_prior <- c(null = 1 / 5, alt_strong = 3 / 10, alt_regular = 1 / 2)
  expected_post <- c(null = 1 / 12, alt_strong = 1 / 2, alt_regular = 5 / 12)
  expected_BF <- 11 / 4

  inference <- compute_inference(prior_weights, margliks, is_null = is_null)

  expect_equal(inference$prior_probs, expected_prior, tolerance = 1e-12)
  expect_equal(unname(inference$post_probs), unname(expected_post), tolerance = 1e-12)
  expect_equal(inference$BF, expected_BF, tolerance = 1e-12)

  model_list <- lapply(seq_along(prior_weights), function(i) {
    list(marglik = .mock_bridge(margliks[i]), prior_weights = prior_weights[i])
  })
  model_summary <- models_inference(model_list)

  expect_equal(
    vapply(model_summary, function(model) model$inference$prior_prob, numeric(1)),
    unname(expected_prior),
    tolerance = 1e-12
  )
  expect_equal(
    vapply(model_summary, function(model) model$inference$post_prob, numeric(1)),
    unname(expected_post),
    tolerance = 1e-12
  )
})

test_that("compute_inference is invariant to model order", {

  prior_weights <- c(null = 2, alt_wide = 3, alt_narrow = 5)
  margliks <- log(c(null = 0.25, alt_wide = 1.5, alt_narrow = 0.5))
  is_null <- c(TRUE, FALSE, FALSE)
  order <- c(3, 1, 2)

  reference <- compute_inference(prior_weights, margliks, is_null = is_null)
  reordered <- compute_inference(
    prior_weights[order],
    margliks[order],
    is_null = is_null[order]
  )

  expect_equal(reordered$prior_probs[order(order)], reference$prior_probs, tolerance = 1e-12)
  expect_equal(reordered$post_probs[order(order)], reference$post_probs, tolerance = 1e-12)
  expect_equal(reordered$BF, reference$BF, tolerance = 1e-12)
})

test_that("conditional compute_inference renormalizes alternatives but keeps full inclusion BF", {

  prior_weights <- c(null = 2, alt_wide = 3, alt_narrow = 5)
  margliks <- log(c(null = 0.25, alt_wide = 1.5, alt_narrow = 0.5))
  is_null <- c(TRUE, FALSE, FALSE)

  unconditional <- compute_inference(prior_weights, margliks, is_null = is_null)
  conditional <- compute_inference(prior_weights, margliks, is_null = is_null, conditional = TRUE)
  conditional_weights <- ifelse(is_null, 0, prior_weights)

  expect_equal(conditional$prior_probs, conditional_weights / sum(conditional_weights), tolerance = 1e-12)
  expect_equal(unname(conditional$post_probs), unname(.expected_post_probs(conditional_weights, margliks)), tolerance = 1e-12)
  expect_equal(conditional$prior_probs[is_null], 0)
  expect_equal(conditional$post_probs[is_null], 0)
  expect_equal(conditional$BF, unconditional$BF, tolerance = 1e-12)
  expect_true(attr(conditional, "conditional"))
  expect_error(
    compute_inference(c(1, 2), c(0, 1), is_null = c(TRUE, TRUE), conditional = TRUE),
    "Conditional inference requires at least one non-null model.",
    fixed = TRUE
  )
})

test_that("compute_inference rejects invalid model-index indicators", {
  expect_error(
    compute_inference(c(1, 1), c(0, 0), is_null = 0),
    "equal or higher than 1"
  )
  expect_error(
    compute_inference(c(1, 1), c(0, 0), is_null = c(TRUE, NA)),
    "cannot contain NA"
  )
  expect_error(
    compute_inference(c(1, 1), c(0, 0), is_null = 1, conditional = NA),
    "cannot contain NA"
  )
})

test_that("ensemble_inference applies exact algebra independently per parameter", {

  prior_weights <- c(1, 2, 3)
  margliks <- log(c(2, 1, 4))
  model_list <- lapply(seq_along(prior_weights), function(i) {
    list(marglik = .mock_bridge(margliks[i]), prior_weights = prior_weights[i])
  })
  is_null_list <- list(
    theta = c(TRUE, FALSE, FALSE),
    beta  = c(FALSE, TRUE, FALSE)
  )

  inference <- ensemble_inference(
    model_list   = model_list,
    parameters   = c("theta", "beta"),
    is_null_list = is_null_list
  )

  expected_prior <- prior_weights / sum(prior_weights)
  expected_post <- .expected_post_probs(prior_weights, margliks)

  for (parameter in names(is_null_list)) {
    is_null <- is_null_list[[parameter]]
    expect_equal(inference[[parameter]]$prior_probs, expected_prior, tolerance = 1e-12)
    expect_equal(inference[[parameter]]$post_probs, expected_post, tolerance = 1e-12)
    expect_equal(
      inference[[parameter]]$BF,
      .expected_inclusion_bf(expected_prior, expected_post, is_null),
      tolerance = 1e-12
    )
    expect_identical(attr(inference[[parameter]], "is_null"), is_null)
    expect_equal(attr(inference[[parameter]], "parameter_name"), parameter)
    expect_false(attr(inference[[parameter]], "conditional"))
  }
  expect_false(attr(inference, "conditional"))

  conditional <- ensemble_inference(
    model_list   = model_list,
    parameters   = c("theta", "beta"),
    is_null_list = is_null_list,
    conditional  = TRUE
  )

  for (parameter in names(is_null_list)) {
    is_null <- is_null_list[[parameter]]
    conditional_weights <- ifelse(is_null, 0, prior_weights)

    expect_equal(
      conditional[[parameter]]$prior_probs,
      conditional_weights / sum(conditional_weights),
      tolerance = 1e-12
    )
    expect_equal(
      conditional[[parameter]]$post_probs,
      .expected_post_probs(conditional_weights, margliks),
      tolerance = 1e-12
    )
    expect_equal(conditional[[parameter]]$BF, inference[[parameter]]$BF, tolerance = 1e-12)
    expect_true(attr(conditional[[parameter]], "conditional"))
  }
  expect_true(attr(conditional, "conditional"))
})

test_that("models_inference preserves model indexes and exact probabilities", {

  prior_weights <- c(1, 1, 2, 2)
  margliks <- log(c(1, 4, 2, 8))
  model_list <- lapply(seq_along(prior_weights), function(i) {
    list(marglik = .mock_bridge(margliks[i]), prior_weights = prior_weights[i])
  })

  out <- models_inference(model_list)
  expected_prior <- prior_weights / sum(prior_weights)
  expected_post <- .expected_post_probs(prior_weights, margliks)
  expected_inclusion_BF <- c(5 / 24, 20 / 21, 8 / 21, 32 / 9)

  expect_equal(vapply(out, function(model) model$inference$m_number, integer(1)), seq_along(model_list))
  expect_equal(vapply(out, function(model) model$inference$marglik, numeric(1)), margliks)
  expect_equal(vapply(out, function(model) model$inference$prior_prob, numeric(1)), expected_prior, tolerance = 1e-12)
  expect_equal(vapply(out, function(model) model$inference$post_prob, numeric(1)), expected_post, tolerance = 1e-12)
  expect_equal(
    vapply(out, function(model) model$inference$inclusion_BF, numeric(1)),
    expected_inclusion_BF,
    tolerance = 1e-12
  )

  perm <- c(4, 1, 3, 2)
  reordered <- models_inference(model_list[perm])
  expect_equal(
    vapply(reordered, function(model) model$inference$post_prob, numeric(1)),
    expected_post[perm],
    tolerance = 1e-12
  )
})

test_that("compute_inference handles zero and tiny prior weights semantically", {

  zero <- compute_inference(
    prior_weights = c(0, 1, 1),
    margliks = log(c(1e100, 1, 4)),
    is_null = c(FALSE, TRUE, FALSE)
  )

  expect_equal(unname(zero$prior_probs), c(0, .5, .5), tolerance = 1e-12)
  expect_equal(unname(zero$post_probs), c(0, .2, .8), tolerance = 1e-12)
  expect_equal(zero$BF, 4, tolerance = 1e-12)

  tiny <- compute_inference(
    prior_weights = c(1e-300, 1, 1),
    margliks = log(c(2, 1, 3)),
    is_null = c(FALSE, TRUE, FALSE)
  )

  expect_gt(tiny$post_probs[1], 0)
  expect_lt(tiny$post_probs[1], 1e-299)
  expect_equal(unname(tiny$post_probs[-1]), c(.25, .75), tolerance = 1e-12)
  expect_equal(tiny$BF, 3, tolerance = 1e-12)
})

test_that("model averaging rejects invalid or unavailable positive-prior marginal likelihoods", {

  expect_error(
    compute_inference(c(1, 1), c(Inf, 0), is_null = c(TRUE, FALSE)),
    "Infinite positive marginal likelihoods are not supported.",
    fixed = TRUE
  )
  expect_error(
    compute_inference(c(1, 0), c(NA_real_, 1000), is_null = c(TRUE, FALSE)),
    "No finite marginal likelihoods are available for models with positive prior probability.",
    fixed = TRUE
  )

  inference <- compute_inference(
    prior_weights = c(1, 1, 0),
    margliks = c(NA_real_, 0, 1000),
    is_null = c(TRUE, FALSE, FALSE)
  )
  expect_equal(inference$post_probs, c(0, 1, 0), tolerance = 1e-12)
})

test_that("model averaging rejects malformed prior weights at each public entry point", {

  expect_error(compute_inference(c(0, 0), c(0, 1)), "At least one prior model weight")
  expect_error(compute_inference(c(-1, 1), c(0, 1)))
  expect_error(compute_inference(c(NA_real_, 1), c(0, 1)))
  expect_error(compute_inference(c(Inf, 1), c(0, 1)), "finite")

  expect_error(
    models_inference(list(
      list(marglik = .mock_bridge(0), prior_weights = 0),
      list(marglik = .mock_bridge(1), prior_weights = 0)
    )),
    "At least one prior model weight"
  )
  expect_error(
    models_inference(list(
      list(marglik = .mock_bridge(0), prior_weights = -1),
      list(marglik = .mock_bridge(1), prior_weights = 1)
    ))
  )

  all_zero_models <- list(
    .mock_mixing_model(offset = 100, logml = 0, prior_weight = 0),
    .mock_mixing_model(offset = 200, logml = 1, prior_weight = 0)
  )
  expect_error(
    mix_posteriors(
      all_zero_models,
      parameters   = "theta",
      is_null_list = list(theta = c(TRUE, FALSE)),
      n_samples    = 4
    ),
    "At least one prior model weight"
  )

  bad_model <- .mock_mixing_model(offset = 100, logml = 0, prior_weight = NA_real_)
  expect_error(
    mix_posteriors(
      list(bad_model),
      parameters   = "theta",
      is_null_list = list(theta = FALSE),
      n_samples    = 4
    )
  )
})

test_that("mixture sample counts are deterministic, exact length, and retain positive components", {
  counts <- BayesTools:::.posterior_mixture_sample_counts(c(.999, .001), 1000)

  expect_equal(sum(counts), 1000)
  expect_equal(counts, c(999L, 1L))
  expect_equal(BayesTools:::.posterior_mixture_sample_counts(c(1, 3), 400), c(100L, 300L))
  expect_equal(BayesTools:::.posterior_mixture_sample_counts(c(.9999, .0001), 1000), c(999L, 1L))
})

test_that("mix_posteriors preserves model and sample alignment across parameters", {

  model_list <- list(
    .mock_mixing_model(offset = 100, logml = log(1)),
    .mock_mixing_model(offset = 200, logml = log(2)),
    .mock_mixing_model(offset = 300, logml = log(1))
  )
  is_null_list <- list(
    theta = c(TRUE, FALSE, FALSE),
    beta  = c(FALSE, TRUE, FALSE)
  )

  mixed <- mix_posteriors(
    model_list   = model_list,
    parameters   = c("theta", "beta"),
    is_null_list = is_null_list,
    seed         = 20260504,
    n_samples    = 12
  )

  expected_models <- c(rep(1L, 3L), rep(2L, 6L), rep(3L, 3L))

  expect_equal(attr(mixed$theta, "models_ind"), expected_models)
  expect_equal(attr(mixed$beta, "models_ind"), expected_models)
  expect_equal(attr(mixed$theta, "sample_ind"), attr(mixed$beta, "sample_ind"))
  expect_equal(as.numeric(mixed$beta - mixed$theta), rep(100, 12))
  expect_equal(as.integer(table(factor(attr(mixed$theta, "models_ind"), levels = 1:3))), c(3L, 6L, 3L))

  sample_ind <- attr(mixed$theta, "sample_ind")
  for(model_i in seq_along(model_list)){
    model_rows <- attr(mixed$theta, "models_ind") == model_i
    expect_true(all(sample_ind[model_rows] >= 1L & sample_ind[model_rows] <= 20L))
    expect_equal(
      unname(mixed$theta[model_rows]),
      c(100, 200, 300)[[model_i]] + sample_ind[model_rows]
    )
  }
})

test_that("conditional mix_posteriors excludes null models from samples and prior weights", {

  model_list <- list(
    .mock_mixing_model(offset = 100, logml = log(1)),
    .mock_mixing_model(offset = 200, logml = log(2)),
    .mock_mixing_model(offset = 300, logml = log(1))
  )

  mixed <- mix_posteriors(
    model_list   = model_list,
    parameters   = "theta",
    is_null_list = list(theta = c(TRUE, FALSE, FALSE)),
    conditional  = TRUE,
    seed         = 20260504,
    n_samples    = 8
  )

  expect_equal(as.integer(table(factor(attr(mixed$theta, "models_ind"), levels = 1:3))), c(0L, 5L, 3L))
  expect_false(any(attr(mixed$theta, "models_ind") == 1L))
  expect_equal(attr(mixed$theta, "sample_ind"), attr(mixed$theta, "sample_ind")[attr(mixed$theta, "models_ind") != 1L])

  mixed_priors <- attr(mixed$theta, "prior_list")
  expect_equal(
    vapply(mixed_priors, `[[`, numeric(1), "prior_weights"),
    c(0, 0.5, 0.5),
    tolerance = 1e-12
  )
})

test_that("inclusion_BF handles all-null models", {

  # All null models assign zero prior mass to the alternative comparison.
  prior_probs <- c(0.5, 0.5)
  post_probs <- c(0.5, 0.5)
  is_null <- c(TRUE, TRUE)

  BF <- inclusion_BF(prior_probs = prior_probs, post_probs = post_probs, is_null = is_null)
  expect_true(is.na(BF))

})


test_that("inclusion_BF handles all-alternative models", {

  # All alternative models assign zero prior mass to the null comparison.
  prior_probs <- c(0.5, 0.5)
  post_probs <- c(0.5, 0.5)
  is_null <- c(FALSE, FALSE)

  BF <- inclusion_BF(prior_probs = prior_probs, post_probs = post_probs, is_null = is_null)
  expect_true(is.na(BF))

})


test_that("inclusion_BF handles single model case", {

  prior_probs <- 1
  post_probs <- 1
  is_null <- FALSE

  BF <- inclusion_BF(prior_probs = prior_probs, post_probs = post_probs, is_null = is_null)
  expect_true(is.na(BF))

  # Single null model
  is_null <- TRUE
  BF <- inclusion_BF(prior_probs = prior_probs, post_probs = post_probs, is_null = is_null)
  expect_true(is.na(BF))

})


test_that("inclusion_BF treats prior and posterior boundaries differently", {

  is_null <- c(TRUE, FALSE)

  expect_equal(
    inclusion_BF(
      prior_probs = c(null = 0.5, alternative = 0.5),
      post_probs  = c(null = 0,   alternative = 1),
      is_null     = is_null
    ),
    Inf
  )
  expect_equal(
    inclusion_BF(
      prior_probs = c(null = 0.5, alternative = 0.5),
      post_probs  = c(null = 1,   alternative = 0),
      is_null     = is_null
    ),
    0
  )
  expect_true(is.na(inclusion_BF(
    prior_probs = c(null = 1, alternative = 0),
    post_probs  = c(null = 0.5, alternative = 0.5),
    is_null     = is_null
  )))
  expect_true(is.na(inclusion_BF(
    prior_probs = c(null = 0, alternative = 1),
    post_probs  = c(null = 0.5, alternative = 0.5),
    is_null     = is_null
  )))

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
    NA_real_
  )

  expect_equal(
    inclusion_BF(
      prior_probs = c(1, 0),
      margliks    = c(0, 1000),
      is_null     = c(FALSE, TRUE)
    ),
    NA_real_
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
