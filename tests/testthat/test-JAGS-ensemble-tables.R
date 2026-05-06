skip_if_not_test_profile("fixture")

# ============================================================================ #
# TEST FILE: JAGS Ensemble Tables Functions
# ============================================================================ #
#
# PURPOSE:
#   Tests for ensemble table generation functions (ensemble_estimates_table,
#   ensemble_inference_table, ensemble_diagnostics_table).
#
# DEPENDENCIES:
#   - rjags, bridgesampling: For tests using pre-fitted models
#   - common-functions.R: temp_fits_dir, skip_if_no_fits, test_reference_table
#
# SKIP CONDITIONS:
#   - First section (empty tables): Can run on CRAN (pure R)
#   - Second section (advanced features): skip_if_no_fits()
#
# MODELS/FIXTURES:
#   - fit_summary0/1/2, fit_formula_interaction_mix/fac
#
# TAGS: @evaluation, @JAGS, @model-averaging, @tables
# ============================================================================ #

REFERENCE_DIR <<- testthat::test_path("..", "results", "JAGS-ensemble-tables")
source(testthat::test_path("common-functions.R"))

.ensemble_estimate_row_names_for_test <- function(samples, parameter, formula_prefix = TRUE) {

  if (is.matrix(samples[[parameter]])) {
    if (inherits(samples[[parameter]], "mixed_posteriors.formula")) {
      return(format_parameter_names(
        colnames(samples[[parameter]]),
        formula_parameters = attr(samples[[parameter]], "formula_parameter"),
        formula_prefix = formula_prefix
      ))
    }
    return(colnames(samples[[parameter]]))
  }

  if (inherits(samples[[parameter]], "mixed_posteriors.formula")) {
    parameter_name <- gsub(
      paste0(attr(samples[[parameter]], "formula_parameter"), "_"),
      if (formula_prefix) paste0("(", attr(samples[[parameter]], "formula_parameter"), ") ") else "",
      parameter
    )
    return(gsub("__xXx__", ":", parameter_name))
  }

  parameter
}

.expect_ensemble_estimates_semantics <- function(table, samples, parameters, probs) {

  expected <- NULL
  for (parameter in parameters) {
    if (is.matrix(samples[[parameter]])) {
      par_expected <- cbind(
        "Mean"   = apply(samples[[parameter]], 2, mean),
        "Median" = apply(samples[[parameter]], 2, stats::median)
      )
      for (prob in probs) {
        par_expected <- cbind(par_expected, apply(samples[[parameter]], 2, stats::quantile, probs = prob))
        colnames(par_expected)[ncol(par_expected)] <- prob
      }
      rownames(par_expected) <- .ensemble_estimate_row_names_for_test(samples, parameter)
    } else {
      par_expected <- c(
        "Mean"   = mean(samples[[parameter]]),
        "Median" = stats::median(samples[[parameter]])
      )
      for (prob in probs) {
        par_expected <- c(par_expected, stats::quantile(samples[[parameter]], probs = prob))
        names(par_expected)[length(par_expected)] <- prob
      }
      par_expected <- matrix(par_expected, nrow = 1, dimnames = list(parameter, names(par_expected)))
      rownames(par_expected) <- .ensemble_estimate_row_names_for_test(samples, parameter)
    }
    expected <- rbind(expected, par_expected)
  }

  expect_equal(rownames(table), rownames(expected))
  expect_equal(colnames(table), colnames(expected))
  expect_equal(unname(as.matrix(table)), unname(expected), tolerance = 1e-12)
  expect_equal(attr(table, "type"), rep("estimate", ncol(table)))
  expect_true(attr(table, "rownames"))
}

.expect_ensemble_inference_semantics <- function(table, inference, parameters, is_null_list) {

  expect_equal(rownames(table), unname(vapply(inference[parameters], attr, character(1), "parameter_name")))
  expect_equal(table$models, unname(vapply(inference[parameters], function(x) sum(!attr(x, "is_null")), numeric(1))))
  expect_equal(attr(table, "n_models"), unname(vapply(inference[parameters], function(x) length(attr(x, "is_null")), integer(1))))

  for (parameter in parameters) {
    is_null <- attr(inference[[parameter]], "is_null")
    parameter_name <- attr(inference[[parameter]], "parameter_name")
    expect_equal(sum(inference[[parameter]]$prior_probs), 1, tolerance = 1e-12)
    expect_equal(sum(inference[[parameter]]$post_probs), 1, tolerance = 1e-12)
    expect_equal(table[parameter_name, "prior_prob"], sum(inference[[parameter]]$prior_probs[!is_null]), tolerance = 1e-12)
    expect_equal(table[parameter_name, "post_prob"], sum(inference[[parameter]]$post_probs[!is_null]), tolerance = 1e-12)
    expect_equal(as.numeric(table[parameter_name, "inclusion_BF"]), inference[[parameter]]$BF, tolerance = 1e-12)
  }

  expect_equal(attr(table, "type"), c("n_models", "prior_prob", "post_prob", "inclusion_BF"))
  expect_true(attr(table, "rownames"))
}

.expect_ensemble_model_table_semantics <- function(summary_table, diagnostics_table, models, parameters) {

  max_or_na <- function(x) {
    if (all(is.na(x))) {
      NA_real_
    } else {
      max(x, na.rm = TRUE)
    }
  }

  min_or_na <- function(x) {
    if (all(is.na(x))) {
      NA_real_
    } else {
      min(x, na.rm = TRUE)
    }
  }

  expect_equal(summary_table$Model, seq_along(models))
  expect_equal(diagnostics_table$Model, seq_along(models))
  expect_equal(summary_table$prior_prob, vapply(models, function(model) model$inference$prior_prob, numeric(1)), tolerance = 1e-12)
  expect_equal(summary_table$marglik, vapply(models, function(model) model$inference$marglik, numeric(1)), tolerance = 1e-12)
  expect_equal(summary_table$post_prob, vapply(models, function(model) model$inference$post_prob, numeric(1)), tolerance = 1e-12)
  expect_equal(as.numeric(summary_table$inclusion_BF), vapply(models, function(model) model$inference$inclusion_BF, numeric(1)), tolerance = 1e-12)
  expect_equal(sum(summary_table$prior_prob), 1, tolerance = 1e-12)
  expect_equal(sum(summary_table$post_prob), 1, tolerance = 1e-12)

  expected_diagnostics <- data.frame(
    max_MCMC_error = vapply(models, function(model) max_or_na(model$fit_summary[, "MCMC_error"]), numeric(1)),
    max_MCMC_SD_error = vapply(models, function(model) max_or_na(model$fit_summary[, "MCMC_SD_error"]), numeric(1)),
    min_ESS = vapply(models, function(model) min_or_na(model$fit_summary[, "ESS"]), numeric(1)),
    max_R_hat = vapply(models, function(model) max_or_na(model$fit_summary[, "R_hat"]), numeric(1))
  )
  actual_diagnostics <- diagnostics_table[, names(expected_diagnostics)]
  actual_diagnostics[] <- lapply(actual_diagnostics, as.numeric)
  expect_equal(actual_diagnostics, expected_diagnostics, tolerance = 1e-12, ignore_attr = TRUE)

  expect_equal(
    attr(summary_table, "type"),
    c("integer", rep("prior", length(parameters)), "prior_prob", "marglik", "post_prob", "inclusion_BF")
  )
  expect_equal(
    attr(diagnostics_table, "type"),
    c("integer", rep("prior", length(parameters)), "max_MCMC_error", "max_MCMC_SD_error", "min_ESS", "max_R_hat")
  )
  expect_false(attr(summary_table, "rownames"))
  expect_false(attr(diagnostics_table, "rownames"))
}

# ============================================================================ #
# SECTION 1: Test Empty Tables (Can run on CRAN - pure R)
# ============================================================================ #
test_that("Empty summary tables work correctly", {

  ensemble_estimates_empty <- ensemble_estimates_empty_table()
  ensemble_inference_empty <- ensemble_inference_empty_table()
  ensemble_diagnostics_empty <- ensemble_diagnostics_empty_table()

  expect_equal(nrow(ensemble_estimates_empty), 0, ignore_attr = TRUE)
  expect_equal(nrow(ensemble_inference_empty), 0, ignore_attr = TRUE)
  expect_equal(nrow(ensemble_diagnostics_empty), 0, ignore_attr = TRUE)
  expect_equal(colnames(ensemble_estimates_empty), c("Mean", "Median", "0.025", "0.975"))
  expect_equal(colnames(ensemble_inference_empty), c("models", "prior_prob", "post_prob", "inclusion_BF"))
  expect_equal(colnames(ensemble_diagnostics_empty), c("Model", "max_MCMC_error", "max_MCMC_SD_error", "min_ESS", "max_R_hat"))
  expect_equal(attr(ensemble_estimates_empty, "type"), rep("estimate", 4))
  expect_equal(attr(ensemble_inference_empty, "type"), c("n_models", "prior_prob", "post_prob", "inclusion_BF"))
  expect_equal(attr(ensemble_diagnostics_empty, "type"), c("integer", "max_MCMC_error", "max_MCMC_SD_error", "min_ESS", "max_R_hat"))

  # Test that empty tables have correct structure
  expect_s3_class(ensemble_estimates_empty, "BayesTools_table")
  expect_s3_class(ensemble_inference_empty, "BayesTools_table")
  expect_s3_class(ensemble_diagnostics_empty, "BayesTools_table")

  test_reference_table(ensemble_estimates_empty, "empty_ensemble_estimates.txt", "Empty ensemble_estimates table mismatch")
  test_reference_table(ensemble_inference_empty, "empty_ensemble_inference.txt", "Empty ensemble_inference table mismatch")
  test_reference_table(ensemble_diagnostics_empty, "empty_ensemble_diagnostics.txt", "Empty ensemble_diagnostics table mismatch")
})

# ============================================================================ #
# SECTION 2: Test Advanced Features (Requires pre-fitted models)
# ============================================================================ #
test_that("Summary table advanced features work correctly", {

  skip_if_no_fits()
  skip_if_not_installed("rjags")
  skip_if_not_installed("bridgesampling")

  # 1. Simple models (m, omega)
  # -------------------------------------------------------------- #
  fit_summary0 <- readRDS(file.path(temp_fits_dir, "fit_summary0.RDS"))
  marglik_summary0 <- readRDS(file.path(temp_marglik_dir, "fit_summary0.RDS"))

  fit_summary1 <- readRDS(file.path(temp_fits_dir, "fit_summary1.RDS"))
  marglik_summary1 <- readRDS(file.path(temp_marglik_dir, "fit_summary1.RDS"))

  fit_summary2 <- readRDS(file.path(temp_fits_dir, "fit_summary2.RDS"))
  marglik_summary2 <- readRDS(file.path(temp_marglik_dir, "fit_summary2.RDS"))

  models <- list(
    list(fit = fit_summary0, marglik = marglik_summary0, prior_weights = 1, fit_summary = runjags_estimates_table(fit_summary0)),
    list(fit = fit_summary1, marglik = marglik_summary1, prior_weights = 1, fit_summary = runjags_estimates_table(fit_summary1)),
    list(fit = fit_summary2, marglik = marglik_summary2, prior_weights = 1, fit_summary = runjags_estimates_table(fit_summary2))
  )
  models <- models_inference(models)

  # Create inference and mixed posteriors
  inference <- ensemble_inference(model_list = models, parameters = c("m", "omega"), is_null_list = list("m" = c(F,F,F), "omega" = c(T,F,F)), conditional = FALSE)
  mixed_posteriors <- mix_posteriors(model_list = models, parameters = c("m", "omega"), is_null_list = list("m" = c(F,F,F), "omega" = c(T,F,F)), seed = 1)

  # Test tables
  estimates_table <- ensemble_estimates_table(mixed_posteriors, parameters = c("m", "omega"), probs = c(.025, 0.95))
  inference_table <- ensemble_inference_table(inference, names(inference))
  summary_table   <- ensemble_summary_table(models, c("m", "omega"))
  diagnostics_table <- ensemble_diagnostics_table(models, c("m", "omega"))

  # Check structure
  expect_s3_class(estimates_table, "BayesTools_table")
  expect_s3_class(inference_table, "BayesTools_table")
  expect_s3_class(summary_table, "BayesTools_table")
  expect_s3_class(diagnostics_table, "BayesTools_table")

  # Check semantic table content directly, not just the printed table snapshots
  .expect_ensemble_estimates_semantics(estimates_table, mixed_posteriors, c("m", "omega"), c(.025, 0.95))
  .expect_ensemble_inference_semantics(inference_table, inference, c("m", "omega"), list("m" = c(FALSE, FALSE, FALSE), "omega" = c(TRUE, FALSE, FALSE)))
  .expect_ensemble_model_table_semantics(summary_table, diagnostics_table, models, c("m", "omega"))

  inference_conditional <- ensemble_inference(
    model_list = models,
    parameters = c("m", "omega"),
    is_null_list = list("m" = c(FALSE, FALSE, FALSE), "omega" = c(TRUE, FALSE, FALSE)),
    conditional = TRUE
  )
  expect_true(attr(inference_conditional, "conditional"))
  expect_equal(inference_conditional$m$prior_probs, inference$m$prior_probs, tolerance = 1e-12)
  expect_equal(inference_conditional$m$post_probs, inference$m$post_probs, tolerance = 1e-12)
  expect_equal(inference_conditional$omega$prior_probs[1], 0, tolerance = 1e-12)
  expect_equal(inference_conditional$omega$post_probs[1], 0, tolerance = 1e-12)
  expect_equal(sum(inference_conditional$omega$prior_probs[-1]), 1, tolerance = 1e-12)
  expect_equal(sum(inference_conditional$omega$post_probs[-1]), 1, tolerance = 1e-12)
  expect_equal(inference_conditional$omega$BF, inference$omega$BF, tolerance = 1e-12)
  expect_error(
    ensemble_inference_table(inference_conditional, names(inference_conditional)),
    "cannot be 'conditional'",
    fixed = TRUE
  )

  mixed_posteriors_conditional <- mix_posteriors(
    model_list = models,
    parameters = "omega",
    is_null_list = list("omega" = c(TRUE, FALSE, FALSE)),
    conditional = TRUE,
    seed = 1,
    n_samples = 1000
  )
  estimates_conditional <- ensemble_estimates_table(mixed_posteriors_conditional, parameters = "omega", probs = c(.025, 0.95))
  expect_equal(rownames(estimates_conditional), rownames(estimates_table)[grepl("^omega\\[", rownames(estimates_table))])
  expect_false(any(attr(mixed_posteriors_conditional$omega, "models_ind") == 1))
  expect_true(any(attr(mixed_posteriors$omega, "models_ind") == 1))

  malformed_model <- models[[1]]
  malformed_model$fit_summary <- data.frame(MCMC_error = 0)
  expect_error(
    ensemble_diagnostics_table(list(malformed_model), c("m", "omega")),
    "fit_summary",
    fixed = TRUE
  )
  expect_error(
    ensemble_summary_table(list(list(fit = fit_summary0)), c("m", "omega")),
    "inference",
    fixed = TRUE
  )

  # Check content with reference files
  test_reference_table(estimates_table, "simple_ensemble_estimates.txt")
  test_reference_table(inference_table, "simple_ensemble_inference.txt")
  test_reference_table(summary_table, "simple_ensemble_summary.txt")
  test_reference_table(diagnostics_table, "simple_ensemble_diagnostics.txt")

  # Test remove_column on diagnostics table
  diagnostics_table.trimmed <- remove_column(diagnostics_table, 2)
  diagnostics_table.trimmed <- remove_column(diagnostics_table.trimmed, 2)
  test_reference_table(diagnostics_table.trimmed, "simple_ensemble_diagnostics_trimmed.txt")

  # Test that trimmed diagnostics table matches empty table structure
  ensemble_diagnostics_empty <- ensemble_diagnostics_empty_table()
  expect_equal(colnames(ensemble_diagnostics_empty), colnames(diagnostics_table.trimmed))
  expect_equal(capture_output_lines(ensemble_diagnostics_empty, width = 150)[1], capture_output_lines(diagnostics_table.trimmed, width = 150)[1])

  # # Test interpret
  inference_for_interpret <- inference
  inference_for_interpret[["m"]][["BF"]] <- 100
  interpretation_samples <- list(m = as.numeric(mixed_posteriors[["m"]]))
  interpretation <- interpret(inference_for_interpret, interpretation_samples, list(
    list(
      inference         = "m",
      samples           = "m",
      inference_name    = "effect",
      inference_BF_name = "BF_10",
      samples_name      = "y",
      samples_units     = NULL
    )
  ), "Test")

  test_reference_text(interpretation, "simple_interpretation.txt")

  # Test interpret 2 (modified inference)
  inference[["m"]][["BF"]] <- 1/5
  interpretation2 <- interpret(inference, interpretation_samples, list(
    list(
      inference           = "m",
      samples             = "m",
      inference_name      = "effect",
      inference_BF_name   = "BF_10",
      samples_name        = "y",
      samples_units       = "mm",
      samples_conditional = TRUE
    ),
    list(
      inference           = "omega",
      inference_name      = "bias",
      inference_BF_name   = "BF_pb"
    )
  ), "Test2")

  test_reference_text(interpretation2, "simple_interpretation2.txt")


  # 2. Complex models (Formula)
  # -------------------------------------------------------------- #
  fit_formula_simple <- readRDS(file.path(temp_fits_dir, "fit_formula_simple.RDS"))
  marglik_formula_simple <- readRDS(file.path(temp_marglik_dir, "fit_formula_simple.RDS"))

  fit_formula_treatment <- readRDS(file.path(temp_fits_dir, "fit_formula_treatment.RDS"))
  marglik_formula_treatment <- readRDS(file.path(temp_marglik_dir, "fit_formula_treatment.RDS"))

  fit_formula_orthonormal <- readRDS(file.path(temp_fits_dir, "fit_formula_orthonormal.RDS"))
  marglik_formula_orthonormal <- readRDS(file.path(temp_marglik_dir, "fit_formula_orthonormal.RDS"))

  models_complex <- list(
    list(fit = fit_formula_simple, marglik = marglik_formula_simple, prior_weights = 1, fit_summary = runjags_estimates_table(fit_formula_simple)),
    list(fit = fit_formula_treatment, marglik = marglik_formula_treatment, prior_weights = 1, fit_summary = runjags_estimates_table(fit_formula_treatment)),
    list(fit = fit_formula_orthonormal, marglik = marglik_formula_orthonormal, prior_weights = 1, fit_summary = runjags_estimates_table(fit_formula_orthonormal))
  )
  models_complex <- models_inference(models_complex)

  parameters_complex <- c("mu_x_cont1", "mu_x_fac2t", "mu_x_fac3o")
  is_null_list_complex <- list(
    "mu_x_cont1" = c(FALSE, FALSE, FALSE),
    "mu_x_fac2t" = c(TRUE,  FALSE, TRUE),
    "mu_x_fac3o" = c(TRUE,  TRUE,  FALSE)
  )

  inference_complex <- ensemble_inference(
    model_list   = models_complex,
    parameters   = parameters_complex,
    is_null_list = is_null_list_complex,
    conditional  = FALSE
  )

  mixed_posteriors_complex <- mix_posteriors(
    model_list   = models_complex,
    parameters   = parameters_complex,
    is_null_list = is_null_list_complex,
    seed = 1, n_samples = 10000
  )

  # Tables
  estimates_table_complex <- ensemble_estimates_table(mixed_posteriors_complex, parameters = parameters_complex, probs = c(.025, 0.95))
  inference_table_complex <- ensemble_inference_table(inference_complex, names(inference_complex))
  summary_table_complex   <- ensemble_summary_table(models_complex, parameters_complex)
  diagnostics_table_complex <- ensemble_diagnostics_table(models_complex, parameters_complex)

  .expect_ensemble_estimates_semantics(estimates_table_complex, mixed_posteriors_complex, parameters_complex, c(.025, 0.95))
  .expect_ensemble_inference_semantics(inference_table_complex, inference_complex, parameters_complex, is_null_list_complex)
  .expect_ensemble_model_table_semantics(summary_table_complex, diagnostics_table_complex, models_complex, parameters_complex)

  test_reference_table(estimates_table_complex, "complex_ensemble_estimates.txt")
  test_reference_table(inference_table_complex, "complex_ensemble_inference.txt")
  test_reference_table(summary_table_complex, "complex_ensemble_summary.txt")
  test_reference_table(diagnostics_table_complex, "complex_ensemble_diagnostics.txt")

  # 3. Simple Spike vs Normal (Model Averaging)
  # -------------------------------------------------------------- #
  fit_simple_spike <- readRDS(file.path(temp_fits_dir, "fit_simple_spike.RDS"))
  marglik_simple_spike <- readRDS(file.path(temp_marglik_dir, "fit_simple_spike.RDS"))

  fit_simple_normal <- readRDS(file.path(temp_fits_dir, "fit_simple_normal.RDS"))
  marglik_simple_normal <- readRDS(file.path(temp_marglik_dir, "fit_simple_normal.RDS"))

  models_simple_ma <- list(
    list(fit = fit_simple_spike, marglik = marglik_simple_spike, prior_weights = 1, fit_summary = runjags_estimates_table(fit_simple_spike)),
    list(fit = fit_simple_normal, marglik = marglik_simple_normal, prior_weights = 1, fit_summary = runjags_estimates_table(fit_simple_normal))
  )
  models_simple_ma <- models_inference(models_simple_ma)

  inference_simple_ma <- ensemble_inference(
    model_list = models_simple_ma,
    parameters = c("m", "s"),
    is_null_list = list("m" = 1, "s" = 0), # m is spike in model 1 (null), s is never null
    conditional = FALSE
  )

  mixed_posteriors_simple_ma <- mix_posteriors(
    model_list = models_simple_ma,
    parameters = c("m", "s"),
    is_null_list = list("m" = 1, "s" = 0),
    seed = 1
  )

  estimates_simple_ma <- ensemble_estimates_table(mixed_posteriors_simple_ma, parameters = c("m", "s"))
  inference_simple_ma_table <- ensemble_inference_table(inference_simple_ma, names(inference_simple_ma))

  .expect_ensemble_estimates_semantics(estimates_simple_ma, mixed_posteriors_simple_ma, c("m", "s"), c(.025, .975))
  .expect_ensemble_inference_semantics(inference_simple_ma_table, inference_simple_ma, c("m", "s"), list("m" = 1, "s" = 0))

  test_reference_table(estimates_simple_ma, "simple_ma_estimates.txt")
  test_reference_table(inference_simple_ma_table, "simple_ma_inference.txt")


  # 4. Fixed Weightfunctions
  # -------------------------------------------------------------- #
  # Re-using summary models 0-2 and adding a fixed weightfunction model
  fit_summary3 <- readRDS(file.path(temp_fits_dir, "fit_summary3.RDS"))
  marglik_summary3 <- readRDS(file.path(temp_marglik_dir, "fit_summary3.RDS"))

  models_fixed_wf <- list(
    list(fit = fit_summary0, marglik = marglik_summary0, prior_weights = 1, fit_summary = runjags_estimates_table(fit_summary0)),
    list(fit = fit_summary1, marglik = marglik_summary1, prior_weights = 1, fit_summary = runjags_estimates_table(fit_summary1)),
    list(fit = fit_summary2, marglik = marglik_summary2, prior_weights = 1, fit_summary = runjags_estimates_table(fit_summary2)),
    list(fit = fit_summary3, marglik = marglik_summary3, prior_weights = 1, fit_summary = runjags_estimates_table(fit_summary3))
  )
  models_fixed_wf <- models_inference(models_fixed_wf)

  inference_fixed_wf <- ensemble_inference(
    model_list = models_fixed_wf,
    parameters = c("m", "omega"),
    is_null_list = list("m" = 0, "omega" = 1),
    conditional = FALSE
  )

  mixed_posteriors_fixed_wf <- mix_posteriors(
    model_list = models_fixed_wf,
    parameters = c("m", "omega"),
    is_null_list = list("m" = 0, "omega" = 1),
    seed = 1
  )

  estimates_fixed_wf <- ensemble_estimates_table(mixed_posteriors_fixed_wf, parameters = c("m", "omega"))
  inference_fixed_wf_table <- ensemble_inference_table(inference_fixed_wf, names(inference_fixed_wf))

  .expect_ensemble_estimates_semantics(estimates_fixed_wf, mixed_posteriors_fixed_wf, c("m", "omega"), c(.025, .975))
  .expect_ensemble_inference_semantics(inference_fixed_wf_table, inference_fixed_wf, c("m", "omega"), list("m" = 0, "omega" = 1))

  test_reference_table(estimates_fixed_wf, "fixed_wf_estimates.txt")
  test_reference_table(inference_fixed_wf_table, "fixed_wf_inference.txt")

  # 5. Interactions
  # -------------------------------------------------------------- #
  fit_formula_interaction_mix <- readRDS(file.path(temp_fits_dir, "fit_formula_interaction_mix.RDS"))
  marglik_formula_interaction_mix <- structure(list(logml = -20), class = "bridge")

  fit_formula_interaction_mix_main <- readRDS(file.path(temp_fits_dir, "fit_formula_interaction_mix_main.RDS"))
  marglik_formula_interaction_mix_main <- structure(list(logml = -22), class = "bridge")

  models_interaction <- list(
    list(fit = fit_formula_interaction_mix_main, marglik = marglik_formula_interaction_mix_main, prior_weights = 1, fit_summary = runjags_estimates_table(fit_formula_interaction_mix_main)),
    list(fit = fit_formula_interaction_mix, marglik = marglik_formula_interaction_mix, prior_weights = 1, fit_summary = runjags_estimates_table(fit_formula_interaction_mix))
  )
  models_interaction <- models_inference(models_interaction)

  parameters_int <- c("mu_x_cont1", "mu_x_fac3o", "mu_x_cont1__xXx__x_fac3o")
  is_null_list_int <- list(
    "mu_x_cont1" = c(FALSE, FALSE),
    "mu_x_fac3o" = c(FALSE, FALSE),
    "mu_x_cont1:x_fac3o" = c(TRUE, FALSE)
  )

  inference_interaction <- ensemble_inference(
    model_list = models_interaction,
    parameters = parameters_int,
    is_null_list = is_null_list_int,
    conditional = FALSE
  )

  mixed_posteriors_interaction <- mix_posteriors(
    model_list = models_interaction,
    parameters = parameters_int,
    is_null_list = is_null_list_int,
    seed = 1
  )

  estimates_interaction <- ensemble_estimates_table(mixed_posteriors_interaction, parameters = parameters_int)
  inference_interaction_table <- ensemble_inference_table(inference_interaction, names(inference_interaction))
  summary_interaction_table <- ensemble_summary_table(models_interaction, parameters_int)
  diagnostics_interaction_table <- ensemble_diagnostics_table(models_interaction, parameters_int)

  .expect_ensemble_estimates_semantics(estimates_interaction, mixed_posteriors_interaction, parameters_int, c(.025, .975))
  .expect_ensemble_inference_semantics(inference_interaction_table, inference_interaction, parameters_int, is_null_list_int)
  .expect_ensemble_model_table_semantics(summary_interaction_table, diagnostics_interaction_table, models_interaction, parameters_int)

  test_reference_table(estimates_interaction, "interaction_ensemble_estimates.txt")
  test_reference_table(inference_interaction_table, "interaction_ensemble_inference.txt")
  test_reference_table(summary_interaction_table, "interaction_ensemble_summary.txt")

  # 6. Spike Factors (using marginal distribution models)
  # -------------------------------------------------------------- #
  # Using fit_marginal_0 (spike) and fit_marginal_1 (normal) which have meandif factors
  fit_spike_factors_null <- readRDS(file.path(temp_fits_dir, "fit_marginal_0.RDS"))
  marglik_spike_factors_null <- readRDS(file.path(temp_marglik_dir, "fit_marginal_0.RDS"))

  fit_spike_factors_alt <- readRDS(file.path(temp_fits_dir, "fit_marginal_1.RDS"))
  marglik_spike_factors_alt <- readRDS(file.path(temp_marglik_dir, "fit_marginal_1.RDS"))

  models_spike_factors <- list(
    list(fit = fit_spike_factors_null, marglik = marglik_spike_factors_null, prior_weights = 1, fit_summary = runjags_estimates_table(fit_spike_factors_null)),
    list(fit = fit_spike_factors_alt, marglik = marglik_spike_factors_alt, prior_weights = 1, fit_summary = runjags_estimates_table(fit_spike_factors_alt))
  )
  models_spike_factors <- models_inference(models_spike_factors)

  inference_spike_factors <- ensemble_inference(
    model_list = models_spike_factors,
    parameters = c("mu_x_fac3md"),
    is_null_list = list("mu_x_fac3md" = c(TRUE, FALSE)),
    conditional = FALSE
  )

  mixed_posteriors_spike_factors <- mix_posteriors(
    model_list = models_spike_factors,
    parameters = c("mu_x_fac3md"),
    is_null_list = list("mu_x_fac3md" = c(TRUE, FALSE)),
    seed = 1
  )

  estimates_spike_factors <- ensemble_estimates_table(mixed_posteriors_spike_factors, parameters = c("mu_x_fac3md"))
  inference_spike_factors_table <- ensemble_inference_table(inference_spike_factors, names(inference_spike_factors))

  .expect_ensemble_estimates_semantics(estimates_spike_factors, mixed_posteriors_spike_factors, c("mu_x_fac3md"), c(.025, .975))
  .expect_ensemble_inference_semantics(inference_spike_factors_table, inference_spike_factors, c("mu_x_fac3md"), list("mu_x_fac3md" = c(TRUE, FALSE)))

  test_reference_table(estimates_spike_factors, "spike_factors_estimates.txt")
  test_reference_table(inference_spike_factors_table, "spike_factors_inference.txt")

})


test_that("Simplified interpret2 function", {

  set.seed(1)
  information <- list(
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

  expect_equal(
    interpret2(information, "RoBMA"),
    "RoBMA found moderate evidence in favor of the Effect, BF10 = 3.50, with mean model-averaged estimate mu = 0.298 kg, 95% CI [-0.020,  0.601]."
  )

})

test_that("as_mixed_posteriors works with ensemble tables", {

  skip_if_no_fits()
  skip_if_not_installed("RoBMA")
  skip_if_missing_fits(c("fit_complex_mixed", "fit_simple_formula_mixed"))
  skip_if_not_installed("rjags")
  skip_if_not_installed("bridgesampling")

  # 1. Complex Mixed Model
  fit_complex_mixed <- readRDS(file.path(temp_fits_dir, "fit_complex_mixed.RDS"))

  mixed_posteriors_complex <- as_mixed_posteriors(
    mode       = fit_complex_mixed,
    parameters = names(attr(fit_complex_mixed, "prior_list"))
  )

  # Generate estimates table
  estimates_table_complex <- ensemble_estimates_table(
    mixed_posteriors_complex,
    parameters = names(attr(fit_complex_mixed, "prior_list")),
    probs = c(.025, 0.95)
  )

  test_reference_table(estimates_table_complex, "as_mixed_posteriors_complex_estimates.txt")


  # 2. Simple Formula Mixed Model
  fit_simple_formula_mixed <- readRDS(file.path(temp_fits_dir, "fit_simple_formula_mixed.RDS"))

  mixed_posteriors_simple <- as_mixed_posteriors(
    mode       = fit_simple_formula_mixed,
    parameters = names(attr(fit_simple_formula_mixed, "prior_list"))
  )

  # Generate estimates table
  estimates_table_simple <- ensemble_estimates_table(
    mixed_posteriors_simple,
    parameters = names(attr(fit_simple_formula_mixed, "prior_list")),
    probs = c(.025, 0.95)
  )

  test_reference_table(estimates_table_simple, "as_mixed_posteriors_simple_estimates.txt")

})

