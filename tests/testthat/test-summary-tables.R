skip_if_not_test_profile("unit")

# ============================================================================ #
# TEST FILE: Summary Tables
# ============================================================================ #
#
# PURPOSE:
#   Tests for summary table functions including ensemble_estimates_table,
#   ensemble_inference_table, ensemble_summary_table, ensemble_diagnostics_table,
#   model_summary_table, and print methods.
#
# DEPENDENCIES:
#   - common-functions.R: deterministic registry and reference helpers
#
# SKIP CONDITIONS:
#   - runjags: One in-memory runjags summary transformation test
#
# MODELS/FIXTURES:
#   - Fixed in-memory posterior, inference, prior, and diagnostics objects
#
# TAGS: @evaluation, @summary-tables
# ============================================================================ #

REFERENCE_DIR <<- testthat::test_path("..", "results", "summary-tables")
source(testthat::test_path("common-functions.R"))

.summary_table_fixed_samples <- function(target_mean, probabilities,
                                         quantiles, bounds){

  # With 201 ordered values, the probabilities used below are exact order
  # statistics. Mixing linear and step interpolation changes the mean without
  # moving those quantiles or introducing random input.
  n_samples <- 201L
  knots <- c(0, probabilities, 1)
  knot_indices <- 1L + as.integer(knots * (n_samples - 1L))
  knot_values <- c(bounds[1L], quantiles, bounds[2L])
  linear <- stats::approx(
    x = knot_indices,
    y = knot_values,
    xout = seq_len(n_samples)
  )$y
  lower_step <- vapply(seq_len(n_samples), function(i){
    knot_values[max(which(knot_indices <= i))]
  }, numeric(1))
  upper_step <- vapply(seq_len(n_samples), function(i){
    knot_values[min(which(knot_indices >= i))]
  }, numeric(1))
  step <- if(target_mean < mean(linear)) lower_step else upper_step
  weight <- (target_mean - mean(linear)) / (mean(step) - mean(linear))
  if(!is.finite(weight) || weight < 0 || weight > 1){
    stop("Fixed summary samples cannot attain the requested mean.", call. = FALSE)
  }

  linear + weight * (step - linear)
}

.summary_table_ensemble_inference_for_test <- function(){

  models <- list(
    list(
      marglik = BayesTools:::.bt_marglik_manual_result(0),
      prior_weights = 1
    ),
    list(
      marglik = BayesTools:::.bt_marglik_manual_result(log(1.68)),
      prior_weights = 1
    )
  )
  ensemble_inference(
    model_list = models,
    parameters = c("m", "omega"),
    is_null_list = list(m = c(FALSE, FALSE), omega = c(TRUE, FALSE))
  )
}

.summary_table_fit_summary_for_test <- function(MCMC_error, MCMC_SD_error,
                                                ESS, R_hat){

  out <- data.frame(
    MCMC_error = MCMC_error,
    MCMC_SD_error = MCMC_SD_error,
    ESS = ESS,
    R_hat = R_hat
  )
  class(out) <- c(
    "BayesTools_table",
    "BayesTools_runjags_summary",
    class(out)
  )

  out
}

.summary_table_models_for_test <- function(){

  s_prior <- prior(
    "normal",
    list(0, 1),
    truncation = list(lower = 0, upper = Inf)
  )
  prior_lists <- list(
    list(m = prior("normal", list(0, 1)), s = s_prior),
    list(m = prior("point", list(0)), s = s_prior)
  )
  fits <- lapply(prior_lists, function(prior_list){
    fit <- list()
    attr(fit, "prior_list") <- prior_list
    fit
  })
  marglik_ratio <- 11.633
  margliks <- c(-29.482 - log(marglik_ratio), -29.482)
  margliks <- lapply(
    margliks,
    BayesTools:::.bt_marglik_manual_result
  )

  models_inference(list(
    list(
      fit = fits[[1L]],
      marglik = margliks[[1L]],
      prior_weights = 1,
      fit_summary = .summary_table_fit_summary_for_test(
        MCMC_error = 0.00211,
        MCMC_SD_error = 0.045,
        ESS = 491,
        R_hat = 1.010
      )
    ),
    list(
      fit = fits[[2L]],
      marglik = margliks[[2L]],
      prior_weights = 1,
      fit_summary = .summary_table_fit_summary_for_test(
        MCMC_error = 0.00207,
        MCMC_SD_error = 0.047,
        ESS = 455,
        R_hat = 1.002
      )
    )
  ))
}

test_that("Stan transformed summaries are recomputed from transformed draws", {

  transformed_samples <- c(1, 2, 4, 8)
  transformed_summary <- BayesTools:::.stan_transformed_summary(transformed_samples, SSeff = 4)

  expect_equal(unname(transformed_summary["Mean"]), mean(transformed_samples))
  expect_equal(unname(transformed_summary["SD"]), stats::sd(transformed_samples))
  expect_equal(
    unname(transformed_summary[c("Lower95", "Median", "Upper95")]),
    unname(stats::quantile(transformed_samples, probs = c(0.025, 0.5, 0.975), names = FALSE))
  )
  expect_equal(unname(transformed_summary["MCerr"]), stats::sd(transformed_samples) / sqrt(4))
})


test_that("ensemble estimates use equal-tailed 95 percent defaults", {
  estimates <- ensemble_estimates_table(
    list(theta = c(1, 2, 3, 4)),
    parameters = "theta"
  )

  expect_equal(colnames(estimates), c("Mean", "Median", "0.025", "0.975"))
})


test_that("update.BayesTools_table remove_parameters removes matching rows", {
  table <- data.frame(value = c("1.000", "2.000"), row.names = c("keep", "drop"))
  class(table) <- c("BayesTools_table", class(table))
  attr(table, "type") <- "estimate"

  updated <- update(table, remove_parameters = "drop")
  expect_equal(rownames(updated), "keep")
})


# ============================================================================ #
# SECTION 1: ensemble_estimates_table tests
# ============================================================================ #
test_that("ensemble_estimates_table handles matrix posteriors", {

  probabilities <- c(0.025, 0.1, 0.5, 0.9, 0.975)
  m_samples <- .summary_table_fixed_samples(
    target_mean = 0.172,
    probabilities = probabilities,
    quantiles = c(-0.221, -0.107, 0.179, 0.427, 0.587),
    bounds = c(-0.3, 0.65)
  )
  omega_samples <- cbind(
    "omega[0,0.05]" = rep(1, length(m_samples)),
    "omega[0.05,1]" = .summary_table_fixed_samples(
      target_mean = 0.677,
      probabilities = probabilities,
      quantiles = c(0.022, 0.134, 0.804, 1, 1),
      bounds = c(0, 1)
    )
  )
  mixed_posteriors <- list(m = m_samples, omega = omega_samples)
  class(mixed_posteriors) <- "mixed_posteriors"

  # Test basic table creation
  estimates_table <- ensemble_estimates_table(
    mixed_posteriors,
    parameters = c("m", "omega")
  )

  test_reference_table(estimates_table, "ensemble_estimates_basic.txt")

  # Test with custom probs
  estimates_table_probs <- ensemble_estimates_table(
    mixed_posteriors,
    parameters = c("m", "omega"),
    probs = c(0.10, 0.50, 0.90)
  )

  test_reference_table(estimates_table_probs, "ensemble_estimates_custom_probs.txt")

})


test_that("ensemble_estimates_table handles transformed factor posteriors", {

  probabilities <- c(0.025, 0.5, 0.975)
  factor_samples <- cbind(
    "mu_x_fac3o[dif: A]" = .summary_table_fixed_samples(
      target_mean = 0.023,
      probabilities = probabilities,
      quantiles = c(-0.185, 0.019, 0.220),
      bounds = c(-0.25, 0.30)
    ),
    "mu_x_fac3o[dif: B]" = .summary_table_fixed_samples(
      target_mean = -0.305,
      probabilities = probabilities,
      quantiles = c(-0.520, -0.320, 0),
      bounds = c(-0.6, 0)
    ),
    "mu_x_fac3o[dif: C]" = .summary_table_fixed_samples(
      target_mean = 0.282,
      probabilities = probabilities,
      quantiles = c(0, 0.289, 0.509),
      bounds = c(0, 0.6)
    )
  )
  class(factor_samples) <- c(
    "mixed_posteriors",
    "mixed_posteriors.factor",
    "mixed_posteriors.vector",
    "mixed_posteriors.formula",
    "mixed_posteriors.orthonormal_transformed",
    class(factor_samples)
  )
  attr(factor_samples, "formula_parameter") <- "mu"
  mixed_posteriors <- list(mu_x_fac3o = factor_samples)
  class(mixed_posteriors) <- "mixed_posteriors"

  # Test with transform_factors = TRUE
  estimates_table_transform <- ensemble_estimates_table(
    mixed_posteriors,
    parameters = "mu_x_fac3o",
    transform_factors = TRUE
  )

  test_reference_table(estimates_table_transform, "ensemble_estimates_transform_factors.txt")

})

test_that("ensemble_estimates_table handles multi-factor transformed interactions", {

  df <- expand.grid(
    a = factor(c("a1", "a2"), levels = c("a1", "a2")),
    b = factor(c("b1", "b2", "b3"), levels = c("b1", "b2", "b3"))
  )
  formula_result <- JAGS_formula(
    formula = ~ a * b,
    parameter = "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      a         = prior_factor("mnormal", list(0, 1), contrast = "orthonormal"),
      b         = prior_factor("mnormal", list(0, 1), contrast = "meandif"),
      "a:b"     = prior_factor("mnormal", list(0, 1), contrast = "orthonormal")
    )
  )
  interaction_prior <- formula_result$prior_list$mu_a__xXx__b

  interaction_samples <- matrix(seq_len(20), nrow = 10, ncol = 2)
  colnames(interaction_samples) <- paste0("mu_a__xXx__b[", 1:2, "]")
  class(interaction_samples) <- c("mixed_posteriors", "mixed_posteriors.factor", "mixed_posteriors.vector", "mixed_posteriors.formula")
  attr(interaction_samples, "levels")            <- BayesTools:::.get_prior_factor_levels(interaction_prior)
  attr(interaction_samples, "level_names")       <- attr(interaction_prior, "level_names")
  attr(interaction_samples, "interaction")       <- TRUE
  attr(interaction_samples, "interaction_terms") <- attr(interaction_prior, "interaction_terms")
  attr(interaction_samples, "term_components")   <- attr(interaction_prior, "term_components")
  attr(interaction_samples, "factor_terms")      <- attr(interaction_prior, "factor_terms")
  attr(interaction_samples, "factor_contrasts")  <- attr(interaction_prior, "factor_contrasts")
  attr(interaction_samples, "factor_design")     <- attr(interaction_prior, "factor_design")
  attr(interaction_samples, "factor_cell_names") <- attr(interaction_prior, "factor_cell_names")
  attr(interaction_samples, "orthonormal")       <- TRUE
  attr(interaction_samples, "meandif")           <- FALSE
  attr(interaction_samples, "treatment")         <- FALSE
  attr(interaction_samples, "independent")       <- FALSE
  attr(interaction_samples, "formula_parameter") <- "mu"

  samples <- list(mu_a__xXx__b = interaction_samples)
  class(samples) <- "mixed_posteriors"

  estimates_table <- ensemble_estimates_table(
    samples = samples,
    parameters = "mu_a__xXx__b",
    transform_factors = TRUE
  )

  expect_equal(
    rownames(estimates_table),
    paste0(
      "(mu) a[dif: ",
      rep(c("a1", "a2"), times = 3),
      "]:b[dif: ",
      rep(c("b1", "b2", "b3"), each = 2),
      "]"
    )
  )
  expect_equal(anyDuplicated(rownames(estimates_table)), 0L)
})

test_that("runjags_estimates_table unscales before parameter filtering", {

  skip_if_not_installed("runjags")

  posterior <- matrix(
    c(
      1.0, 1.2,
      0.4, 0.5,
      2.0, 2.2,
      0.8, 1.0
    ),
    nrow = 2,
    byrow = FALSE
  )
  colnames(posterior) <- c("mu_intercept", "mu_x", "mu_a", "mu_x__xXx__a[1]")

  factor_prior <- prior_factor("normal", list(0, 1), contrast = "treatment")
  attr(factor_prior, "levels") <- 2
  attr(factor_prior, "level_names") <- c("A", "B")
  attr(factor_prior, "parameter") <- "mu"

  interaction_prior <- prior_factor("normal", list(0, 1), contrast = "treatment")
  attr(interaction_prior, "levels") <- 2
  attr(interaction_prior, "level_names") <- list(a = c("A", "B"))
  attr(interaction_prior, "interaction") <- TRUE
  attr(interaction_prior, "interaction_terms") <- c("x", "a")
  attr(interaction_prior, "term_components") <- c("x", "a")
  attr(interaction_prior, "factor_terms") <- "a"
  attr(interaction_prior, "factor_contrasts") <- c(a = "contr.treatment")
  attr(interaction_prior, "factor_design") <- stats::contr.treatment(c("A", "B"))
  attr(interaction_prior, "factor_cell_names") <- c("A", "B")
  attr(interaction_prior, "parameter") <- "mu"

  prior_list <- list(
    mu_intercept = prior("normal", list(0, 1)),
    mu_x = prior("normal", list(0, 1)),
    mu_a = factor_prior,
    mu_x__xXx__a = interaction_prior
  )
  for(parameter in c("mu_intercept", "mu_x")){
    attr(prior_list[[parameter]], "parameter") <- "mu"
  }

  fit <- list(
    mcmc = coda::mcmc.list(coda::mcmc(posterior)),
    summary.pars = list(mutate = NULL),
    monitor = colnames(posterior)
  )
  class(fit) <- c("runjags", "BayesTools_fit")
  attr(fit, "prior_list") <- prior_list
  attr(fit, "formula_scale") <- list(mu = list(mu_x = list(mean = 10, sd = 2)))
  fit <- attach_test_parameter_registry(
    fit,
    monitor_names = colnames(posterior)
  )

  samples <- suppressWarnings(runjags_estimates_table(
    fit,
    keep_parameters = "mu_a",
    transform_scaled = TRUE,
    formula_prefix = FALSE,
    return_samples = TRUE
  ))
  expected <- posterior[, "mu_a"] - (posterior[, "mu_x__xXx__a[1]"] / 2) * 10

  expect_equal(as.numeric(samples[, "a[B]"]), expected, tolerance = 1e-10)
  expect_false("x:a[B]" %in% colnames(samples))
})


test_that("ensemble_estimates_table handles formula posteriors", {

  probabilities <- c(0.025, 0.5, 0.975)
  intercept_samples <- .summary_table_fixed_samples(
    target_mean = 0.514,
    probabilities = probabilities,
    quantiles = c(0.355, 0.512, 0.678),
    bounds = c(0.3, 0.75)
  )
  class(intercept_samples) <- c(
    "mixed_posteriors",
    "mixed_posteriors.vector",
    "mixed_posteriors.formula",
    class(intercept_samples)
  )
  attr(intercept_samples, "formula_parameter") <- "mu"
  sigma_samples <- .summary_table_fixed_samples(
    target_mean = 0.887,
    probabilities = probabilities,
    quantiles = c(0.783, 0.885, 1.004),
    bounds = c(0.7, 1.1)
  )
  mixed_posteriors <- list(
    mu_intercept = intercept_samples,
    sigma = sigma_samples
  )
  class(mixed_posteriors) <- "mixed_posteriors"
  params <- names(mixed_posteriors)

  # Test with formula_prefix = TRUE
  estimates_prefix_true <- ensemble_estimates_table(
    mixed_posteriors,
    parameters = params,
    formula_prefix = TRUE
  )

  # Test with formula_prefix = FALSE
  estimates_prefix_false <- ensemble_estimates_table(
    mixed_posteriors,
    parameters = params,
    formula_prefix = FALSE
  )

  test_reference_table(estimates_prefix_true, "ensemble_estimates_formula_prefix_true.txt")
  test_reference_table(estimates_prefix_false, "ensemble_estimates_formula_prefix_false.txt")

})


# ============================================================================ #
# SECTION 2: ensemble_inference_table tests
# ============================================================================ #
test_that("ensemble_inference_table handles multiple parameters", {

  inference <- .summary_table_ensemble_inference_for_test()

  # Basic table
  inference_table <- ensemble_inference_table(inference, names(inference))
  test_reference_table(inference_table, "ensemble_inference_basic.txt")

  # With logBF
  inference_table_log <- ensemble_inference_table(inference, names(inference), logBF = TRUE)
  test_reference_table(inference_table_log, "ensemble_inference_logBF.txt")

  # With BF01
  inference_table_bf01 <- ensemble_inference_table(inference, names(inference), BF01 = TRUE)
  test_reference_table(inference_table_bf01, "ensemble_inference_BF01.txt")

  # With both
  inference_table_both <- ensemble_inference_table(inference, names(inference), logBF = TRUE, BF01 = TRUE)
  test_reference_table(inference_table_both, "ensemble_inference_both.txt")

})


# ============================================================================ #
# SECTION 3: ensemble_summary_table and ensemble_diagnostics_table tests
# ============================================================================ #
test_that("ensemble_summary_table handles different model configurations", {

  models <- .summary_table_models_for_test()

  # Test summary table
  summary_table <- ensemble_summary_table(models, c("m", "s"))
  test_reference_table(summary_table, "ensemble_summary_basic.txt")

  # Test with short_name
  summary_table_short <- ensemble_summary_table(models, c("m", "s"), short_name = TRUE)
  test_reference_table(summary_table_short, "ensemble_summary_short_name.txt")

  # Test with logBF and BF01
  summary_table_bf <- ensemble_summary_table(models, c("m", "s"), logBF = TRUE, BF01 = TRUE)
  test_reference_table(summary_table_bf, "ensemble_summary_bf_options.txt")

  # Test with remove_spike_0
  summary_table_no_spike <- ensemble_summary_table(models, c("m", "s"), remove_spike_0 = FALSE)
  test_reference_table(summary_table_no_spike, "ensemble_summary_no_spike.txt")

})


test_that("ensemble_summary_table handles parameters as list", {

  models <- .summary_table_models_for_test()

  # Test with parameters supplied as a list
  pars <- list("m" = "m", "renamed 2" = "s")
  summary_table_list <- ensemble_summary_table(models, pars)
  test_reference_table(summary_table_list, "ensemble_summary_params_list.txt")

})


test_that("ensemble_diagnostics_table handles different configurations", {

  models <- .summary_table_models_for_test()

  # Test diagnostics table
  diagnostics_table <- ensemble_diagnostics_table(models, c("m", "s"))
  test_reference_table(diagnostics_table, "ensemble_diagnostics_basic.txt")

  # Test with short_name
  diagnostics_short <- ensemble_diagnostics_table(models, c("m", "s"), short_name = TRUE)
  test_reference_table(diagnostics_short, "ensemble_diagnostics_short_name.txt")

  # Test with remove_spike_0
  diagnostics_no_spike <- ensemble_diagnostics_table(models, c("m", "s"), remove_spike_0 = FALSE)
  test_reference_table(diagnostics_no_spike, "ensemble_diagnostics_no_spike.txt")

})


# ============================================================================ #
# SECTION 4: marginal_estimates_table tests
# ============================================================================ #
test_that("marginal_estimates_table handles various inputs", {

  # Use a fixed quantile grid rather than a seeded stochastic realization so
  # the exact presentation references remain deterministic numerical oracles.
  probability_grid <- (seq_len(1001L) - 0.5) / 1001
  samples <- list(
    mu = stats::qnorm(probability_grid)
  )

  inference <- list(
    mu = structure(list(
      BF = 2.5,
      prior_probs = c(0.5, 0.5),
      post_probs = c(0.4, 0.6)
    ), class = c("list", "marginal_inference"))
  )

  attr(inference$mu, "is_null") <- c(TRUE, FALSE)
  attr(inference$mu, "prior_list") <- list(
    prior("spike", list(0)),
    prior("normal", list(0, 1))
  )

  marginal_table <- marginal_estimates_table(
    samples = samples,
    inference = inference,
    parameters = "mu"
  )

  test_reference_table(marginal_table, "marginal_estimates_basic.txt")

  # With logBF
  marginal_table_log <- marginal_estimates_table(
    samples = samples,
    inference = inference,
    parameters = "mu",
    logBF = TRUE
  )
  test_reference_table(marginal_table_log, "marginal_estimates_logBF.txt")

  # With BF01
  marginal_table_bf01 <- marginal_estimates_table(
    samples = samples,
    inference = inference,
    parameters = "mu",
    BF01 = TRUE
  )
  test_reference_table(marginal_table_bf01, "marginal_estimates_BF01.txt")

})

test_that("marginal_estimates_table reports Savage-Dickey BF error attributes", {

  set.seed(1)
  samples <- list(
    mu = list(
      a = rnorm(100),
      b = rnorm(100, mean = 1)
    )
  )
  BF_a <- 2.5
  BF_b <- 4
  attr(BF_a, "BF_error_percent") <- 6.25
  inference <- list(
    mu = list(
      a = BF_a,
      b = BF_b
    )
  )

  marginal_table <- marginal_estimates_table(
    samples    = samples,
    inference  = inference,
    parameters = "mu"
  )

  expect_true("BF_error_percent" %in% colnames(marginal_table))
  expect_equal(attr(marginal_table, "type"), c("estimate", "estimate", "estimate", "estimate", "estimate", "inclusion_BF", "BF_error"))
  expect_equal(as.numeric(marginal_table[["BF_error_percent"]]), c(6.25, NA_real_))
  expect_equal(attr(marginal_table[["BF_error_percent"]], "name"), "error%(Inclusion BF)")

  marginal_table_BF01 <- update(marginal_table, BF01 = TRUE)
  expect_equal(as.numeric(marginal_table_BF01[["BF_error_percent"]]), c(6.25, NA_real_))
  expect_equal(attr(marginal_table_BF01[["BF_error_percent"]], "name"), "error%(Inclusion BF)")
})


# ============================================================================ #
# SECTION 5: model_summary_table tests
# ============================================================================ #
test_that("model_summary_table handles various configurations", {

  model <- .summary_table_models_for_test()[[2L]]
  model <- models_inference(list(model))[[1L]]

  # Basic model summary
  summary_table <- model_summary_table(model)
  test_reference_table(summary_table, "model_summary_basic.txt")

  # With short_name
  summary_short <- model_summary_table(model, short_name = TRUE)
  test_reference_table(summary_short, "model_summary_short_name.txt")

  # With remove_spike_0 (should remove 'm' which has spike at zero)
  summary_no_spike <- model_summary_table(model, remove_spike_0 = TRUE)
  test_reference_table(summary_no_spike, "model_summary_no_spike.txt")

})


# ============================================================================ #
# SECTION 6: update.BayesTools_table tests
# ============================================================================ #
test_that("update.BayesTools_table works correctly", {

  inference <- .summary_table_ensemble_inference_for_test()

  # Create inference table
  inference_table <- ensemble_inference_table(inference, names(inference))

  # Update with new title
  updated_table <- update(inference_table, title = "Updated Title")
  test_reference_table(updated_table, "update_table_new_title.txt")

  # Update with footnotes
  updated_footnotes <- update(inference_table, footnotes = "This is a footnote")
  test_reference_table(updated_footnotes, "update_table_footnotes.txt")

  # Update with warnings
  updated_warnings <- update(inference_table, warnings = "This is a warning")
  test_reference_table(updated_warnings, "update_table_warnings.txt")

  # Update with logBF
  updated_logbf <- update(inference_table, logBF = TRUE)
  test_reference_table(updated_logbf, "update_table_logBF.txt")

  # Update with BF01
  updated_bf01 <- update(inference_table, BF01 = TRUE)
  test_reference_table(updated_bf01, "update_table_BF01.txt")

})
