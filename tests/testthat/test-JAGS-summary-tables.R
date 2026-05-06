skip_if_not_test_profile("fixture")

# ============================================================================ #
# TEST FILE: JAGS Summary Tables Functions
# ============================================================================ #
#
# PURPOSE:
#   Tests for runjags_estimates_table, runjags_inference_table, and related
#   summary table generation functions.
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
#   - fit_formula_interaction_cont, fit_factor_treatment, fit_spike_slab_factor
#   - fit_factor_orthonormal
#
# TAGS: @evaluation, @JAGS, @summary-tables
# ============================================================================ #

REFERENCE_DIR <<- testthat::test_path("..", "results", "JAGS-summary-tables")
source(testthat::test_path("common-functions.R"))

make_adversarial_indicator_fit <- function(indicator = c(rep(1, 9900), rep(0, 100))){

  samples <- coda::mcmc(
    cbind(theta = indicator, theta_indicator = indicator),
    start = 1,
    end   = length(indicator),
    thin  = 1
  )
  fit <- structure(
    list(mcmc = coda::mcmc.list(samples), sample = length(indicator)),
    class = c("runjags", "BayesTools_fit", "list")
  )
  attr(fit, "prior_list") <- list(
    theta = prior_spike_and_slab(
      prior("normal", list(0, 1)),
      prior_inclusion = prior("point", list(0.5))
    )
  )

  fit
}

expect_summary_row_order <- function(table, rows) {
  expect_equal(rownames(table), rows)
}

expect_summary_values <- function(table, expected, tolerance = 5e-3) {
  missing_rows <- setdiff(rownames(expected), rownames(table))
  missing_cols <- setdiff(colnames(expected), colnames(table))

  expect_equal(missing_rows, character())
  expect_equal(missing_cols, character())

  actual_values <- as.matrix(table[rownames(expected), colnames(expected), drop = FALSE])
  storage.mode(actual_values) <- "numeric"

  expect_equal(unname(actual_values), unname(expected), tolerance = tolerance)
}

summary_values <- function(..., columns = c("Mean", "SD", "0.025", "0.5", "0.975")) {
  values <- rbind(...)
  colnames(values) <- columns
  values
}

# ============================================================================ #
# SECTION 1: Test Empty Tables
# ============================================================================ #
test_that("Empty summary tables work correctly", {

  runjags_summary_empty <- runjags_estimates_empty_table()

  expect_equal(nrow(runjags_summary_empty), 0, ignore_attr = TRUE)

  # Test that empty tables have correct structure
  expect_s3_class(runjags_summary_empty, "BayesTools_table")

  test_reference_table(runjags_summary_empty, "empty_runjags_estimates.txt", "Empty runjags_estimates table mismatch")
})


test_that("runjags_inference_table reports known BF error percent in adversarial near-boundary case", {

  fit <- make_adversarial_indicator_fit()

  known_post_prob <- 0.99
  known_mcse      <- 0.0025
  expected_BF     <- (known_post_prob / (1 - known_post_prob)) / (0.5 / (1 - 0.5))
  expected_error  <- 100 * known_mcse / (known_post_prob * (1 - known_post_prob))

  testthat::local_mocked_bindings(
    .indicator_MCMC_summary = function(indicator_list){
      list(
        post_prob     = known_post_prob,
        MCMC_error    = known_mcse,
        MCMC_SD_error = NA_real_,
        ESS           = 1600,
        visits        = sum(unlist(indicator_list)),
        n_samples     = length(unlist(indicator_list))
      )
    },
    .package = "BayesTools"
  )

  inference <- runjags_inference_table(fit, BF_diagnostics = TRUE)

  expect_equal(as.numeric(inference["theta", "prior_prob"]), 0.5)
  expect_equal(as.numeric(inference["theta", "post_prob"]), known_post_prob)
  expect_equal(as.numeric(inference["theta", "inclusion_BF"]), expected_BF)
  expect_equal(as.numeric(inference["theta", "MCMC_error"]), known_mcse)
  expect_equal(as.numeric(inference["theta", "BF_error_percent"]), expected_error)
  expect_gt(as.numeric(inference["theta", "BF_error_percent"]), 25)

  inference_BF01 <- runjags_inference_table(fit, BF_diagnostics = TRUE, BF01 = TRUE)
  expect_equal(as.numeric(inference_BF01["theta", "inclusion_BF"]), 1 / expected_BF)
  expect_equal(as.numeric(inference_BF01["theta", "BF_error_percent"]), expected_error)
  expect_equal(attr(inference_BF01[["BF_error_percent"]], "name"), "error%(Inclusion BF)")

  test_reference_table(inference, "runjags_inference_adversarial_diagnostics.txt", "Adversarial BF diagnostics table mismatch")
  test_reference_table(inference_BF01, "runjags_inference_adversarial_diagnostics_BF01.txt", "Adversarial BF01 diagnostics table mismatch")
})


test_that("runjags_inference_table supports selected BF diagnostics from arguments and options", {

  fit <- make_adversarial_indicator_fit()

  inference_ess <- runjags_inference_table(fit, BF_diagnostic_columns = "ESS")
  expect_equal(colnames(inference_ess), c("prior_prob", "post_prob", "inclusion_BF", "ESS"))
  expect_equal(attr(inference_ess, "type"), c("prior_prob", "post_prob", "inclusion_BF", "ESS"))
  expect_null(attr(inference_ess, "footnotes"))
  expect_true(is.finite(as.numeric(inference_ess["theta", "ESS"])))

  inference_none <- runjags_inference_table(fit, BF_diagnostic_columns = FALSE)
  expect_equal(colnames(inference_none), c("prior_prob", "post_prob", "inclusion_BF"))
  expect_equal(attr(inference_none, "type"), c("prior_prob", "post_prob", "inclusion_BF"))

  old_options <- options(BayesTools.JAGS_BF_diagnostic_columns = c("MCMC_error", "BF_error_percent"))
  on.exit(options(old_options), add = TRUE)

  inference_option <- runjags_inference_table(fit)
  expect_equal(colnames(inference_option), c("prior_prob", "post_prob", "inclusion_BF", "MCMC_error", "BF_error_percent"))
  expect_equal(attr(inference_option, "type"), c("prior_prob", "post_prob", "inclusion_BF", "MCMC_error", "BF_error"))
  expect_equal(attr(inference_option[["MCMC_error"]], "name"), "error(Post. prob.)")
  expect_equal(attr(inference_option[["BF_error_percent"]], "name"), "error%(Inclusion BF)")
  expect_null(attr(inference_option, "footnotes"))
})


test_that("runjags_inference_table propagates real indicator MCSE to BF error percent", {

  fit       <- make_adversarial_indicator_fit()
  inference <- runjags_inference_table(fit, BF_diagnostics = TRUE)

  post_prob      <- as.numeric(inference["theta", "post_prob"])
  prior_prob     <- as.numeric(inference["theta", "prior_prob"])
  expected_BF    <- (post_prob / (1 - post_prob)) / (prior_prob / (1 - prior_prob))
  expected_error <- 100 * as.numeric(inference["theta", "MCMC_error"]) / (post_prob * (1 - post_prob))

  expect_equal(post_prob, 0.99)
  expect_equal(prior_prob, 0.5)
  expect_equal(as.numeric(inference["theta", "inclusion_BF"]), expected_BF)
  expect_true(is.finite(as.numeric(inference["theta", "MCMC_error"])))
  expect_equal(as.numeric(inference["theta", "BF_error_percent"]), expected_error, tolerance = 1e-12)
})


# ============================================================================ #
# SECTION 2: Test Advanced Features (Transformations, Formula Handling, etc.)
# ============================================================================ #
test_that("Summary table advanced features work correctly", {

  skip_if_not_installed("rjags")
  skip_if_not_installed("bridgesampling")
  skip_if_no_fits()

  # Use fit_formula_interaction_cont for testing advanced features
  # This model has continuous interactions and formulas
  fit_complex <- readRDS(file.path(temp_fits_dir, "fit_formula_interaction_cont.RDS"))

  # Test 1: Parameter transformations
  runjags_summary_transform <- runjags_estimates_table(
    fit_complex,
    transformations = list("mu_intercept" = list(fun = exp))
  )

  # Test 2: Formula handling with prefix
  runjags_summary_prefix_true <- runjags_estimates_table(fit_complex, formula_prefix = TRUE)
  runjags_summary_prefix_false <- runjags_estimates_table(fit_complex, formula_prefix = FALSE)

  # Test 3: Conditional vs unconditional
  runjags_summary_conditional <- runjags_estimates_table(fit_complex, conditional = TRUE)
  runjags_summary_unconditional <- runjags_estimates_table(fit_complex, conditional = FALSE)

  # Test 4: Factor transformations (use fit_factor_treatment for this)
  fit_factor <- readRDS(file.path(temp_fits_dir, "fit_factor_treatment.RDS"))

  runjags_summary_factor <- runjags_estimates_table(fit_factor)

  # Test 5: Use fit with spike and slab
  fit_spike <- readRDS(file.path(temp_fits_dir, "fit_spike_slab_factor.RDS"))

  runjags_summary_spike <- runjags_estimates_table(fit_spike)
  runjags_inference_spike <- runjags_inference_table(fit_spike)
  runjags_inference_spike_diagnostics <- runjags_inference_table(fit_spike, BF_diagnostics = TRUE)
  inclusion_rows <- grepl("(inclusion)", rownames(runjags_summary_spike), fixed = TRUE)
  expect_true(any(inclusion_rows))
  expect_true(any(is.finite(runjags_summary_spike[inclusion_rows, "MCMC_error"])))
  expect_true(any(is.finite(runjags_summary_spike[inclusion_rows, "ESS"])))
  expect_equal(colnames(runjags_inference_spike_diagnostics), c("prior_prob", "post_prob", "inclusion_BF", "ESS", "MCMC_error", "BF_error_percent"))
  expect_equal(attr(runjags_inference_spike_diagnostics[["BF_error_percent"]], "name"), "error%(Inclusion BF)")
  expect_null(attr(runjags_inference_spike_diagnostics, "footnotes"))

  runjags_inference_spike_log_BF01 <- runjags_inference_table(fit_spike, BF_diagnostics = TRUE, logBF = TRUE, BF01 = TRUE)
  expect_equal(as.numeric(runjags_inference_spike_log_BF01[["inclusion_BF"]]), log(1 / as.numeric(runjags_inference_spike_diagnostics[["inclusion_BF"]])), tolerance = 1e-10)
  expect_equal(as.numeric(runjags_inference_spike_log_BF01[["BF_error_percent"]]), as.numeric(runjags_inference_spike_diagnostics[["BF_error_percent"]]), tolerance = 1e-10)
  expect_equal(attr(runjags_inference_spike_log_BF01[["BF_error_percent"]], "name"), "error%(Inclusion BF)")

  # Test 6: Orthonormal contrast transformations to differences from the mean
  fit_orthonormal <- readRDS(file.path(temp_fits_dir, "fit_factor_orthonormal.RDS"))

  runjags_summary_orthonormal  <- suppressMessages(runjags_estimates_table(fit_orthonormal, transform_factors = TRUE))
  runjags_summary_orthonormal2 <- suppressMessages(runjags_estimates_table(fit_orthonormal, transform_factors = TRUE,
                                                                          transformations = list("p1" = list(fun = exp))))

  # Test 7: Custom transformations with transform_factors = FALSE
  # Use a model with factor parameters for transformation testing
  runjags_summary_custom_transform <- suppressMessages(runjags_estimates_table(
    fit_factor,
    transform_factors = FALSE,
    transformations = list("mu_x_fac3t[2]" = list(fun = exp))
  ))

  # Test 8: Conditional with remove_inclusion
  runjags_summary_remove_inclusion <- suppressMessages(runjags_estimates_table(
    fit_spike,
    transform_factors = TRUE,
    conditional = TRUE,
    remove_inclusion = TRUE
  ))

  # Test 9: Conditional estimates with mixture priors
  complex_optional_fit_files <- file.path(
    temp_fits_dir,
    c("fit_complex_bias.RDS", "fit_complex_mixed.RDS")
  )
  has_complex_optional_fits <- all(file.exists(complex_optional_fit_files))
  if (has_complex_optional_fits) {
    fit_complex_bias  <- readRDS(complex_optional_fit_files[[1]])
    fit_complex_mixed <- readRDS(complex_optional_fit_files[[2]])
    runjags_summary_complex2 <- runjags_estimates_table(fit_complex_bias, conditional = TRUE)
    runjags_summary_complex3 <- runjags_estimates_table(fit_complex_mixed, conditional = TRUE)
  }

  # Test basic properties
  expect_s3_class(runjags_summary_transform, "BayesTools_table")
  expect_s3_class(runjags_summary_prefix_true, "BayesTools_table")
  expect_s3_class(runjags_summary_prefix_false, "BayesTools_table")
  expect_s3_class(runjags_summary_conditional, "BayesTools_table")
  expect_s3_class(runjags_summary_unconditional, "BayesTools_table")
  expect_s3_class(runjags_summary_factor, "BayesTools_table")
  expect_s3_class(runjags_summary_spike, "BayesTools_table")
  expect_s3_class(runjags_inference_spike, "BayesTools_table")
  expect_s3_class(runjags_summary_orthonormal, "BayesTools_table")
  expect_s3_class(runjags_summary_orthonormal2, "BayesTools_table")
  expect_s3_class(runjags_summary_custom_transform, "BayesTools_table")
  expect_s3_class(runjags_summary_remove_inclusion, "BayesTools_table")
  if (has_complex_optional_fits) {
    expect_s3_class(runjags_summary_complex2, "BayesTools_table")
    expect_s3_class(runjags_summary_complex3, "BayesTools_table")
  }

  # Test that row names differ with different formula_prefix settings
  expect_false(identical(rownames(runjags_summary_prefix_true),
                        rownames(runjags_summary_prefix_false)))

  # Test that remove_inclusion reduces the number of rows
  expect_true(nrow(runjags_summary_remove_inclusion) <= nrow(runjags_summary_spike))

  # Semantic companions for the presentation snapshots below.
  expected_formula_rows_prefixed <- c(
    "(mu) intercept",
    "(mu) x_cont1",
    "(mu) x_cont2",
    "(mu) x_cont1:x_cont2",
    "sigma"
  )
  expected_formula_rows_unprefixed <- c(
    "intercept",
    "x_cont1",
    "x_cont2",
    "x_cont1:x_cont2",
    "sigma"
  )
  expected_formula_values <- summary_values(
    "(mu) intercept"        = c( 0.035, 0.104, -0.178,  0.040,  0.243),
    "(mu) x_cont1"          = c( 0.361, 0.123,  0.124,  0.365,  0.587),
    "(mu) x_cont2"          = c(-0.029, 0.109, -0.241, -0.028,  0.190),
    "(mu) x_cont1:x_cont2"  = c(-0.391, 0.149, -0.685, -0.390, -0.104),
    "sigma"                 = c( 1.040, 0.077,  0.905,  1.033,  1.203)
  )
  expected_transform_values <- expected_formula_values
  expected_transform_values["(mu) intercept", ] <- c(1.041, 0.109, 0.837, 1.040, 1.276)

  expect_summary_row_order(runjags_summary_prefix_true, expected_formula_rows_prefixed)
  expect_summary_row_order(runjags_summary_prefix_false, expected_formula_rows_unprefixed)
  expect_summary_row_order(runjags_summary_transform, expected_formula_rows_prefixed)
  expect_summary_values(runjags_summary_prefix_true, expected_formula_values)
  expect_summary_values(runjags_summary_transform, expected_transform_values)
  expect_equal(attr(runjags_summary_prefix_true, "parameters"),
               c("mu_intercept", "mu_x_cont1", "mu_x_cont2", "mu_x_cont1__xXx__x_cont2", "sigma"))

  expect_equal(colnames(runjags_summary_conditional), c("Mean", "SD", "0.025", "0.5", "0.975"))
  expect_equal(attr(runjags_summary_conditional, "type"), rep("estimate", 5))
  expect_summary_values(runjags_summary_conditional, expected_formula_values)
  expect_true(all(c("MCMC_error", "MCMC_SD_error", "ESS", "R_hat") %in% colnames(runjags_summary_unconditional)))
  expect_equal(attr(runjags_summary_unconditional, "type"),
               c(rep("estimate", 5), "MCMC_error", "MCMC_SD_error", "ESS", "R_hat"))

  expect_summary_row_order(runjags_summary_factor, "p1[2]")
  expect_summary_values(runjags_summary_factor, summary_values("p1[2]" = c(0.497, 0.291, 0.022, 0.491, 0.973)))

  expect_summary_row_order(runjags_summary_spike, c("beta (inclusion)", "beta[1]", "beta[2]"))
  expect_summary_values(
    runjags_summary_spike,
    summary_values(
      "beta[1]" = c(0.034, 0.747, -1.569, 0.000, 1.759),
      "beta[2]" = c(0.005, 0.721, -1.680, 0.000, 1.639)
    )
  )
  expect_equal(as.numeric(runjags_summary_spike["beta (inclusion)", "Mean"]), 0.527, tolerance = 5e-3)
  expect_equal(rownames(runjags_inference_spike), "beta")
  expect_equal(as.numeric(runjags_inference_spike["beta", "prior_prob"]), 0.5)
  expect_equal(as.numeric(runjags_inference_spike["beta", "post_prob"]),
               as.numeric(runjags_summary_spike["beta (inclusion)", "Mean"]), tolerance = 5e-3)
  expect_equal(as.numeric(runjags_inference_spike["beta", "inclusion_BF"]), 1.114, tolerance = 5e-3)

  expect_summary_row_order(runjags_summary_orthonormal, c("p1[dif: 1]", "p1[dif: 2]", "p1[dif: 3]"))
  expect_summary_values(
    runjags_summary_orthonormal,
    summary_values(
      "p1[dif: 1]" = c( 0.041, 0.818, -1.631,  0.060, 1.605),
      "p1[dif: 2]" = c(-0.033, 0.796, -1.612, -0.029, 1.527),
      "p1[dif: 3]" = c(-0.008, 0.811, -1.550, -0.009, 1.564)
    )
  )
  expect_equal(sum(runjags_summary_orthonormal[,"Mean"]), 0, tolerance = 1e-12)

  expect_summary_row_order(runjags_summary_remove_inclusion, c("beta[dif: 1]", "beta[dif: 2]", "beta[dif: 3]"))
  expect_false(any(grepl("inclusion", rownames(runjags_summary_remove_inclusion), fixed = TRUE)))
  expect_equal(colnames(runjags_summary_remove_inclusion), c("Mean", "SD", "0.025", "0.5", "0.975"))

  test_reference_table(runjags_summary_transform, "advanced_transform.txt", "Transform table mismatch")
  test_reference_table(runjags_summary_prefix_true, "advanced_formula_prefix_true.txt", "Formula prefix true table mismatch")
  test_reference_table(runjags_summary_prefix_false, "advanced_formula_prefix_false.txt", "Formula prefix false table mismatch")
  test_reference_table(runjags_summary_conditional, "advanced_conditional.txt", "Conditional table mismatch")
  test_reference_table(runjags_summary_unconditional, "advanced_unconditional.txt", "Unconditional table mismatch")
  test_reference_table(runjags_summary_factor, "advanced_factor_treatment.txt", "Factor treatment table mismatch")
  test_reference_table(runjags_summary_spike, "advanced_spike_slab_estimates.txt", "Spike slab estimates table mismatch")
  test_reference_table(runjags_inference_spike, "advanced_spike_slab_inference.txt", "Spike slab inference table mismatch")
  test_reference_table(runjags_summary_orthonormal, "advanced_orthonormal_transform.txt", "Orthonormal transform table mismatch")
  test_reference_table(runjags_summary_orthonormal2, "advanced_orthonormal_transform2.txt", "Orthonormal transform2 table mismatch")
  test_reference_table(runjags_summary_custom_transform, "advanced_custom_transform.txt", "Custom transform table mismatch")
  test_reference_table(runjags_summary_remove_inclusion, "advanced_remove_inclusion.txt", "Remove inclusion table mismatch")
  if (has_complex_optional_fits) {
    expect_summary_row_order(
      runjags_summary_complex2,
      c("mu", "bias (inclusion)", "omega[0,0.025]", "omega[0.025,0.05]",
        "omega[0.05,0.975]", "omega[0.975,1]", "PET", "PEESE")
    )
    expect_summary_values(
      runjags_summary_complex2,
      summary_values(
        "mu"                = c(1.024, 0.631, 0.208, 0.894, 2.639),
        "omega[0,0.025]"   = c(1.000, 0.000, 1.000, 1.000, 1.000),
        "PET"               = c(0.816, 0.659, 0.063, 0.583, 2.295),
        "PEESE"             = c(1.434, 1.318, 0.031, 1.043, 4.889)
      )
    )
    expect_equal(as.numeric(runjags_summary_complex2["bias (inclusion)", "Mean"]), 0.526, tolerance = 5e-3)
    expect_true(any(grepl("Conditional summary for PET is based on 79 samples.",
                          attr(runjags_summary_complex2, "warnings"), fixed = TRUE)))

    expect_summary_row_order(
      runjags_summary_complex3,
      c("(mu) intercept (inclusion)", "(mu) intercept",
        "(mu) x_cont1 (inclusion)", "(mu) x_cont1",
        "(mu) x_fac2t (inclusion)", "(mu) x_fac2t",
        "(mu) x_fac3t (inclusion)", "(mu) x_fac3t[1]", "(mu) x_fac3t[2]",
        "sigma (inclusion: normal)", "sigma (inclusion: lognormal)",
        "sigma[normal]", "sigma[lognormal]", "bias (inclusion)",
        "omega[0,0.025]", "omega[0.025,0.05]", "omega[0.05,0.975]",
        "omega[0.975,1]", "PET")
    )
    expect_equal(as.numeric(runjags_summary_complex3["(mu) x_cont1 (inclusion)", "Mean"]), 1)
    expect_equal(
      sum(as.numeric(runjags_summary_complex3[c("sigma (inclusion: normal)", "sigma (inclusion: lognormal)"), "Mean"])),
      1,
      tolerance = 5e-3
    )

    test_reference_table(runjags_summary_complex2, "runjags_summary_complex2.txt", "Custom probs table mismatch")
    test_reference_table(runjags_summary_complex3, "runjags_summary_complex3.txt", "Custom probs table mismatch")
  }

  # Removal of formula and parameter names
  fit_dual_param  <- readRDS(file.path(temp_fits_dir, "fit_dual_param_regression.RDS"))

  runjags_summary_removal_01 <- JAGS_estimates_table(fit_dual_param)
  runjags_summary_removal_02 <- JAGS_estimates_table(fit_dual_param, remove_formulas = "mu")
  runjags_summary_removal_03 <- JAGS_estimates_table(fit_dual_param, keep_formulas = "log_sigma")

  expected_dual_rows <- c("(mu) intercept", "(mu) x_mu", "(log_sigma) intercept", "(log_sigma) x_sigma")
  expected_log_sigma_rows <- c("(log_sigma) intercept", "(log_sigma) x_sigma")
  expected_dual_values <- summary_values(
    "(mu) intercept"        = c( 2.498, 0.008,  2.482,  2.498,  2.514),
    "(mu) x_mu"             = c( 0.631, 0.008,  0.617,  0.631,  0.646),
    "(log_sigma) intercept" = c( 0.285, 0.006,  0.273,  0.285,  0.298),
    "(log_sigma) x_sigma"   = c(-0.324, 0.024, -0.370, -0.325, -0.280)
  )
  expect_summary_row_order(runjags_summary_removal_01, expected_dual_rows)
  expect_summary_row_order(runjags_summary_removal_02, expected_log_sigma_rows)
  expect_summary_row_order(runjags_summary_removal_03, expected_log_sigma_rows)
  expect_summary_values(runjags_summary_removal_01, expected_dual_values)
  expect_summary_values(runjags_summary_removal_02, expected_dual_values[expected_log_sigma_rows,, drop = FALSE])
  expect_equal(runjags_summary_removal_02, runjags_summary_removal_03, ignore_attr = TRUE)

  test_reference_table(runjags_summary_removal_01, "summary_parameter_or_formula_removal01.txt", "Parameter/formula removal")
  test_reference_table(runjags_summary_removal_02, "summary_parameter_or_formula_removal02.txt", "Parameter/formula removal")
  test_reference_table(runjags_summary_removal_03, "summary_parameter_or_formula_removal03.txt", "Parameter/formula removal")
  if (has_complex_optional_fits) {
    runjags_summary_removal_04 <- JAGS_estimates_table(fit_complex_mixed)
    runjags_summary_removal_05 <- JAGS_estimates_table(fit_complex_mixed, remove_parameters = TRUE)
    runjags_summary_removal_06 <- JAGS_estimates_table(fit_complex_mixed, remove_parameters = TRUE, remove_inclusion = TRUE)
    runjags_summary_removal_07 <- JAGS_estimates_table(fit_complex_mixed, remove_parameters = "bias")
    runjags_summary_removal_08 <- JAGS_estimates_table(fit_complex_mixed, remove_parameters = "sigma")
    runjags_summary_removal_09 <- JAGS_estimates_table(fit_complex_mixed, remove_formulas = "mu")

    expect_true(any(grepl("inclusion", rownames(runjags_summary_removal_04), fixed = TRUE)))
    expect_false(any(grepl("inclusion", rownames(runjags_summary_removal_06), fixed = TRUE)))
    expect_false(any(grepl("^bias", rownames(runjags_summary_removal_07))))
    expect_false(any(grepl("^sigma", rownames(runjags_summary_removal_08))))
    expect_false(any(grepl("^\\(mu\\)", rownames(runjags_summary_removal_09))))

    test_reference_table(runjags_summary_removal_04, "summary_parameter_or_formula_removal04.txt", "Parameter/formula removal")
    test_reference_table(runjags_summary_removal_05, "summary_parameter_or_formula_removal05.txt", "Parameter/formula removal")
    test_reference_table(runjags_summary_removal_06, "summary_parameter_or_formula_removal06.txt", "Parameter/formula removal")
    test_reference_table(runjags_summary_removal_07, "summary_parameter_or_formula_removal07.txt", "Parameter/formula removal")
    test_reference_table(runjags_summary_removal_08, "summary_parameter_or_formula_removal08.txt", "Parameter/formula removal")
    test_reference_table(runjags_summary_removal_09, "summary_parameter_or_formula_removal09.txt", "Parameter/formula removal")
  }

  # Custom probs
  runjags_summary_probs_01 <- JAGS_estimates_table(fit_dual_param)
  runjags_summary_probs_02 <- JAGS_estimates_table(fit_dual_param, probs = c(0.5))
  runjags_summary_probs_03 <- JAGS_estimates_table(fit_dual_param, probs = c(0.25, 0.20, 0.99))

  expect_summary_row_order(runjags_summary_probs_01, expected_dual_rows)
  expect_summary_row_order(runjags_summary_probs_02, expected_dual_rows)
  expect_summary_row_order(runjags_summary_probs_03, expected_dual_rows)
  expect_equal(colnames(runjags_summary_probs_02),
               c("Mean", "SD", "0.5", "MCMC_error", "MCMC_SD_error", "ESS", "R_hat"))
  expect_equal(colnames(runjags_summary_probs_03),
               c("Mean", "SD", "0.25", "0.2", "0.99", "MCMC_error", "MCMC_SD_error", "ESS", "R_hat"))
  expect_summary_values(
    runjags_summary_probs_02,
    summary_values(
      "(mu) intercept"        = c( 2.498, 0.008,  2.498),
      "(mu) x_mu"             = c( 0.631, 0.008,  0.631),
      "(log_sigma) intercept" = c( 0.285, 0.006,  0.285),
      "(log_sigma) x_sigma"   = c(-0.324, 0.024, -0.325),
      columns = c("Mean", "SD", "0.5")
    )
  )

  test_reference_table(runjags_summary_probs_01, "summary_parameter_probs1.txt", "Parameter/formula removal")
  test_reference_table(runjags_summary_probs_02, "summary_parameter_probs2.txt", "Parameter/formula removal")
  test_reference_table(runjags_summary_probs_03, "summary_parameter_probs3.txt", "Parameter/formula removal")

  # Remove diagnostics
  runjags_remove_diagnostics <- JAGS_estimates_table(fit_dual_param, remove_diagnostics = TRUE)
  expect_summary_row_order(runjags_remove_diagnostics, expected_dual_rows)
  expect_equal(colnames(runjags_remove_diagnostics), c("Mean", "SD", "0.025", "0.5", "0.975"))
  expect_summary_values(runjags_remove_diagnostics, expected_dual_values)
  test_reference_table(runjags_remove_diagnostics, "runjags_remove_diagnostics.txt", "Diagnostics removal")

  # Selected diagnostics
  runjags_selected_diagnostics <- JAGS_estimates_table(fit_dual_param, diagnostic_columns = c("ESS", "R_hat"))
  expect_summary_row_order(runjags_selected_diagnostics, expected_dual_rows)
  expect_equal(colnames(runjags_selected_diagnostics), c("Mean", "SD", "0.025", "0.5", "0.975", "ESS", "R_hat"))
  expect_equal(attr(runjags_selected_diagnostics, "type"), c(rep("estimate", 5), "ESS", "R_hat"))
  expect_summary_values(runjags_selected_diagnostics, expected_dual_values)

  old_options <- options(BayesTools.JAGS_estimates_diagnostic_columns = "ESS")
  on.exit(options(old_options), add = TRUE)

  runjags_option_diagnostics <- JAGS_estimates_table(fit_dual_param)
  expect_equal(colnames(runjags_option_diagnostics), c("Mean", "SD", "0.025", "0.5", "0.975", "ESS"))
  expect_equal(attr(runjags_option_diagnostics, "type"), c(rep("estimate", 5), "ESS"))
})

# ============================================================================ #
# SECTION 3: Test Summary Tables for All Saved Models
# ============================================================================ #
test_that("Summary tables for all saved models", {

  skip_if_not_installed("rjags")
  skip_if_not_installed("bridgesampling")
  skip_if_no_fits()

  runjags::runjags.options(silent.jags = TRUE, silent.runjags = TRUE)

  # Load model registry to get list of all fitted models
  registry_file <- file.path(test_files_dir, "model_registry.RDS")

  model_registry <- readRDS(registry_file)
  model_names <- model_registry$model_name

  print_dir <- testthat::test_path("..", "results", "print")

  for (model_name in model_names) {
    fit_file <- file.path(temp_fits_dir, paste0(model_name, ".RDS"))
    marglik_file <- file.path(temp_marglik_dir, paste0(model_name, ".RDS"))

    fit <- readRDS(fit_file)
    has_marglik <- file.exists(marglik_file)

    if (has_marglik) {
      marglik <- readRDS(marglik_file)
    }

    # Process model summary table
    if (has_marglik) {
      model_list <- list(
        list(fit = fit, marglik = marglik, prior_weights = 1,
             fit_summary = runjags_estimates_table(fit))
      )
      model_list <- models_inference(model_list)
      model_summary <- model_summary_table(model_list[[1]])
      test_reference_table(model_summary, paste0(model_name, "_model_summary.txt"),
                       paste0("Model summary mismatch for ", model_name))
    }

    # Process runjags estimates table
    runjags_summary <- runjags_estimates_table(fit)
    test_reference_table_numeric(
      runjags_summary,
      paste0(model_name, "_runjags_estimates.txt"),
      info_msg = paste0("Runjags estimates mismatch for ", model_name)
    )

  }
})


# ============================================================================ #
# SECTION 4: Test runjags_estimates_table with conditional=TRUE on various priors
# ============================================================================ #
test_that("runjags_estimates_table with conditional=TRUE on various prior types", {

  skip_if_not_installed("rjags")
  skip_if_not_installed("bridgesampling")
  skip_on_cran()
  skip_if_no_fits()

  # Test with publication bias priors
  fit_pub_bias <- readRDS(file.path(temp_fits_dir, "fit_simple_pub_bias.RDS"))
  runjags_pub_bias_conditional <- runjags_estimates_table(fit_pub_bias, conditional = TRUE)
  test_reference_table(runjags_pub_bias_conditional, "runjags_pub_bias_conditional.txt")

  # Test with factor priors
  fit_factor <- readRDS(file.path(temp_fits_dir, "fit_factor_orthonormal.RDS"))
  runjags_factor_conditional <- runjags_estimates_table(fit_factor, conditional = TRUE)
  test_reference_table(runjags_factor_conditional, "runjags_factor_conditional.txt")

  runjags_factor_conditional_transformed <- runjags_estimates_table(fit_factor, conditional = TRUE, transform_factors = TRUE)
  test_reference_table(runjags_factor_conditional_transformed, "runjags_factor_conditional_transformed.txt")

  # Test with mixture priors
  fit_mixture <- readRDS(file.path(temp_fits_dir, "fit_mixture_simple.RDS"))
  runjags_mixture_conditional <- runjags_estimates_table(fit_mixture, conditional = TRUE)
  test_reference_table(runjags_mixture_conditional, "runjags_mixture_conditional.txt")

  # Test with spike and slab priors
  fit_spike_slab <- readRDS(file.path(temp_fits_dir, "fit_spike_slab_simple.RDS"))
  runjags_spike_slab_conditional <- runjags_estimates_table(fit_spike_slab, conditional = TRUE)
  test_reference_table(runjags_spike_slab_conditional, "runjags_spike_slab_conditional.txt")

})


# ============================================================================ #
# SECTION 5: Test runjags_inference_table with mixture and formula priors
# ============================================================================ #
test_that("runjags_inference_table with mixture priors", {

  skip_if_not_installed("rjags")
  skip_if_not_installed("bridgesampling")
  skip_on_cran()
  skip_if_no_fits()

  # Test with mixture priors
  fit_mixture <- readRDS(file.path(temp_fits_dir, "fit_mixture_simple.RDS"))
  runjags_mixture_inference <- runjags_inference_table(fit_mixture)
  test_reference_table(runjags_mixture_inference, "runjags_mixture_inference.txt")

  runjags_mixture_inference_diagnostics <- runjags_inference_table(fit_mixture, BF_diagnostics = TRUE)
  expect_equal(attr(runjags_mixture_inference_diagnostics, "type"), c("prior_prob", "post_prob", "inclusion_BF", "ESS", "MCMC_error", "BF_error"))
  expect_equal(colnames(runjags_mixture_inference_diagnostics), c("prior_prob", "post_prob", "inclusion_BF", "ESS", "MCMC_error", "BF_error_percent"))
  expect_equal(attr(runjags_mixture_inference_diagnostics[["MCMC_error"]], "name"), "error(Post. prob.)")
  expect_null(attr(runjags_mixture_inference_diagnostics, "footnotes"))
  expect_true(all(c("ESS", "MCMC_error", "BF_error_percent") %in% colnames(runjags_mixture_inference_diagnostics)))

  runjags_mixture_inference_log <- runjags_inference_table(fit_mixture, BF_diagnostics = TRUE, logBF = TRUE)
  expect_equal(as.numeric(runjags_mixture_inference_log[["inclusion_BF"]]), log(as.numeric(runjags_mixture_inference_diagnostics[["inclusion_BF"]])), tolerance = 1e-10)
  expect_equal(as.numeric(runjags_mixture_inference_log[["BF_error_percent"]]), as.numeric(runjags_mixture_inference_diagnostics[["BF_error_percent"]]), tolerance = 1e-10)

  runjags_mixture_inference_BF01 <- runjags_inference_table(fit_mixture, BF_diagnostics = TRUE, BF01 = TRUE)
  expect_equal(as.numeric(runjags_mixture_inference_BF01[["inclusion_BF"]]), 1 / as.numeric(runjags_mixture_inference_diagnostics[["inclusion_BF"]]), tolerance = 1e-10)
  expect_equal(as.numeric(runjags_mixture_inference_BF01[["BF_error_percent"]]), as.numeric(runjags_mixture_inference_diagnostics[["BF_error_percent"]]), tolerance = 1e-10)
  expect_equal(attr(runjags_mixture_inference_BF01[["BF_error_percent"]], "name"), "error%(Inclusion BF)")

  runjags_mixture_inference_log_BF01 <- runjags_inference_table(fit_mixture, BF_diagnostics = TRUE, logBF = TRUE, BF01 = TRUE)
  expect_equal(as.numeric(runjags_mixture_inference_log_BF01[["inclusion_BF"]]), log(1 / as.numeric(runjags_mixture_inference_diagnostics[["inclusion_BF"]])), tolerance = 1e-10)
  expect_equal(as.numeric(runjags_mixture_inference_log_BF01[["BF_error_percent"]]), as.numeric(runjags_mixture_inference_diagnostics[["BF_error_percent"]]), tolerance = 1e-10)
  expect_equal(attr(runjags_mixture_inference_log_BF01[["BF_error_percent"]], "name"), "error%(Inclusion BF)")

  # Test with mixture containing spike
  fit_mixture_spike <- readRDS(file.path(temp_fits_dir, "fit_mixture_spike.RDS"))
  runjags_mixture_spike_inference <- runjags_inference_table(fit_mixture_spike)
  test_reference_table(runjags_mixture_spike_inference, "runjags_mixture_spike_inference.txt")

  # Test component-specific product-space diagnostics
  fit_components <- readRDS(file.path(temp_fits_dir, "fit_mixture_components.RDS"))
  runjags_components_inference <- runjags_inference_table(fit_components, BF_diagnostics = TRUE)
  components <- attr(attr(fit_components, "prior_list")[["beta"]], "components")
  samples    <- suppressWarnings(coda::as.mcmc(fit_components))
  expected_post <- sapply(unique(components), function(component){
    mean(samples[,"beta_indicator"] %in% which(components == component))
  })
  names(expected_post) <- paste0("beta [", names(expected_post), "]")
  expect_equal(as.numeric(runjags_components_inference[names(expected_post), "post_prob"]), as.numeric(expected_post), tolerance = 1e-10)
  expect_true(all(is.finite(runjags_components_inference[,"MCMC_error"])))
  expect_true(all(is.finite(runjags_components_inference[,"BF_error_percent"])))

})


test_that("runjags_inference_table with formula priors", {

  skip_if_not_installed("rjags")
  skip_if_not_installed("bridgesampling")
  skip_on_cran()
  skip_if_no_fits()

  # Test with formula + mixture priors (mixture on factor predictor)
  fit_formula_mixture <- readRDS(file.path(temp_fits_dir, "fit_formula_factor_mixture.RDS"))
  runjags_formula_mixture_inference <- runjags_inference_table(fit_formula_mixture)
  test_reference_table(runjags_formula_mixture_inference, "runjags_formula_mixture_inference.txt")

  # Test with joint complex model (formula + mixture + spike-and-slab)
  fit_joint <- readRDS(file.path(temp_fits_dir, "fit_joint_complex.RDS"))
  runjags_joint_inference <- runjags_inference_table(fit_joint)
  test_reference_table(runjags_joint_inference, "runjags_joint_complex_inference.txt")

})


# ============================================================================ #
# SECTION 6: Test empty tables
# ============================================================================ #
test_that("model_summary_empty_table works correctly", {

  empty_table <- model_summary_empty_table()
  expect_s3_class(empty_table, "BayesTools_table")
  expect_true(nrow(empty_table) > 0)
  expect_true(ncol(empty_table) > 0)
  test_reference_table(empty_table, "model_summary_empty.txt")

})


test_that("runjags_inference_empty_table works correctly", {

  empty_table <- runjags_inference_empty_table()
  expect_s3_class(empty_table, "BayesTools_table")
  expect_equal(nrow(empty_table), 0)
  expect_true(ncol(empty_table) > 0)
  test_reference_table(empty_table, "runjags_inference_empty.txt")

  empty_diagnostics <- runjags_inference_empty_table(BF_diagnostics = TRUE, logBF = TRUE, BF01 = TRUE)
  expect_s3_class(empty_diagnostics, "BayesTools_table")
  expect_equal(colnames(empty_diagnostics), c("prior_prob", "post_prob", "inclusion_BF", "ESS", "MCMC_error", "BF_error_percent"))
  expect_equal(attr(empty_diagnostics[["inclusion_BF"]], "name"), "log(Exclusion BF)")
  expect_equal(attr(empty_diagnostics[["BF_error_percent"]], "name"), "error%(Inclusion BF)")
  expect_null(attr(empty_diagnostics, "footnotes"))

})


# ============================================================================ #
# SECTION 7: Test stan_estimates_table with stored RDS file
# ============================================================================ #
test_that("stan_estimates_table works with stored fit", {

  skip_if_not_installed("rstan")

  # Load stored stan fit from tests/results/fits
  stan_fit_file <- testthat::test_path("..", "results", "fits", "fit_RoBTT.RDS")

  fit_stan <- readRDS(stan_fit_file)

  # Test basic stan_estimates_table
  stan_summary <- stan_estimates_table(fit_stan)
  test_reference_table(stan_summary, "stan_estimates_basic.txt")

  stan_summary2 <- stan_estimates_table(fit_stan, transformations = list("mu" = list(fun = exp)))
  test_reference_table(stan_summary2, "stan_estimates_basic2.txt")

})
