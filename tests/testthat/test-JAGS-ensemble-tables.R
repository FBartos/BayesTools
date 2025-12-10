context("JAGS ensemble tables functions")

REFERENCE_DIR <<- testthat::test_path("..", "results", "JAGS-ensemble-tables")
source(testthat::test_path("common-functions.R"))

# ============================================================================ #
# SECTION 1: Test Empty Tables
# ============================================================================ #
test_that("Empty summary tables work correctly", {

  ensemble_estimates_empty <- ensemble_estimates_empty_table()
  ensemble_inference_empty <- ensemble_inference_empty_table()
  ensemble_diagnostics_empty <- ensemble_diagnostics_empty_table()

  expect_equivalent(nrow(ensemble_estimates_empty), 0)
  expect_equivalent(nrow(ensemble_inference_empty), 0)
  expect_equivalent(nrow(ensemble_diagnostics_empty), 0)

  # Test that empty tables have correct structure
  expect_s3_class(ensemble_estimates_empty, "BayesTools_table")
  expect_s3_class(ensemble_inference_empty, "BayesTools_table")
  expect_s3_class(ensemble_diagnostics_empty, "BayesTools_table")

  test_reference_table(ensemble_estimates_empty, "empty_ensemble_estimates.txt", "Empty ensemble_estimates table mismatch")
  test_reference_table(ensemble_inference_empty, "empty_ensemble_inference.txt", "Empty ensemble_inference table mismatch")
  test_reference_table(ensemble_diagnostics_empty, "empty_ensemble_diagnostics.txt", "Empty ensemble_diagnostics table mismatch")
})

# ============================================================================ #
# SECTION 2: Test Advanced Features (Transformations, Formula Handling, etc.)
# ============================================================================ #
test_that("Summary table advanced features work correctly", {

  skip_if_not_installed("rjags")
  skip_if_not_installed("bridgesampling")

  # 1. Simple models (m, omega)
  # -------------------------------------------------------------- #
  fit_summary0 <- readRDS(file.path(temp_fits_dir, "fit_summary0.RDS"))
  marglik_summary0 <- readRDS(file.path(temp_fits_dir, "fit_summary0_marglik.RDS"))

  fit_summary1 <- readRDS(file.path(temp_fits_dir, "fit_summary1.RDS"))
  marglik_summary1 <- readRDS(file.path(temp_fits_dir, "fit_summary1_marglik.RDS"))

  fit_summary2 <- readRDS(file.path(temp_fits_dir, "fit_summary2.RDS"))
  marglik_summary2 <- readRDS(file.path(temp_fits_dir, "fit_summary2_marglik.RDS"))

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

  # Test interpret
  interpretation <- interpret(inference, mixed_posteriors, list(
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
  interpretation2 <- interpret(inference, mixed_posteriors, list(
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
  marglik_formula_simple <- readRDS(file.path(temp_fits_dir, "fit_formula_simple_marglik.RDS"))

  fit_formula_treatment <- readRDS(file.path(temp_fits_dir, "fit_formula_treatment.RDS"))
  marglik_formula_treatment <- readRDS(file.path(temp_fits_dir, "fit_formula_treatment_marglik.RDS"))

  fit_formula_orthonormal <- readRDS(file.path(temp_fits_dir, "fit_formula_orthonormal.RDS"))
  marglik_formula_orthonormal <- readRDS(file.path(temp_fits_dir, "fit_formula_orthonormal_marglik.RDS"))

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

  test_reference_table(estimates_table_complex, "complex_ensemble_estimates.txt")
  test_reference_table(inference_table_complex, "complex_ensemble_inference.txt")
  test_reference_table(summary_table_complex, "complex_ensemble_summary.txt")
  test_reference_table(diagnostics_table_complex, "complex_ensemble_diagnostics.txt")

  # 3. Simple Spike vs Normal (Model Averaging)
  # -------------------------------------------------------------- #
  fit_simple_spike <- readRDS(file.path(temp_fits_dir, "fit_simple_spike.RDS"))
  marglik_simple_spike <- readRDS(file.path(temp_fits_dir, "fit_simple_spike_marglik.RDS"))

  fit_simple_normal <- readRDS(file.path(temp_fits_dir, "fit_simple_normal.RDS"))
  marglik_simple_normal <- readRDS(file.path(temp_fits_dir, "fit_simple_normal_marglik.RDS"))

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

  test_reference_table(estimates_simple_ma, "simple_ma_estimates.txt")
  test_reference_table(inference_simple_ma_table, "simple_ma_inference.txt")


  # 4. Fixed Weightfunctions
  # -------------------------------------------------------------- #
  # Re-using summary models 0-2 and adding a fixed weightfunction model
  fit_summary3 <- readRDS(file.path(temp_fits_dir, "fit_summary3.RDS"))
  marglik_summary3 <- readRDS(file.path(temp_fits_dir, "fit_summary3_marglik.RDS"))

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

  test_reference_table(estimates_fixed_wf, "fixed_wf_estimates.txt")
  test_reference_table(inference_fixed_wf_table, "fixed_wf_inference.txt")

  # 5. Interactions
  # -------------------------------------------------------------- #
  fit_formula_interaction_mix <- readRDS(file.path(temp_fits_dir, "fit_formula_interaction_mix.RDS"))
  marglik_formula_interaction_mix <- readRDS(file.path(temp_fits_dir, "fit_formula_interaction_mix_marglik.RDS"))

  fit_formula_interaction_mix_main <- readRDS(file.path(temp_fits_dir, "fit_formula_interaction_mix_main.RDS"))
  marglik_formula_interaction_mix_main <- readRDS(file.path(temp_fits_dir, "fit_formula_interaction_mix_main_marglik.RDS"))

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

  test_reference_table(estimates_interaction, "interaction_ensemble_estimates.txt")
  test_reference_table(inference_interaction_table, "interaction_ensemble_inference.txt")
  test_reference_table(summary_interaction_table, "interaction_ensemble_summary.txt")

  # 6. Spike Factors (using marginal distribution models)
  # -------------------------------------------------------------- #
  # Using fit_marginal_0 (spike) and fit_marginal_1 (normal) which have meandif factors
  fit_spike_factors_null <- readRDS(file.path(temp_fits_dir, "fit_marginal_0.RDS"))
  marglik_spike_factors_null <- readRDS(file.path(temp_fits_dir, "fit_marginal_0_marglik.RDS"))

  fit_spike_factors_alt <- readRDS(file.path(temp_fits_dir, "fit_marginal_1.RDS"))
  marglik_spike_factors_alt <- readRDS(file.path(temp_fits_dir, "fit_marginal_1_marglik.RDS"))

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

