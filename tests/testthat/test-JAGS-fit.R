skip_if_not_test_profile("fixture")

# ============================================================================ #
# TEST FILE: JAGS Fit Functions
# ============================================================================ #
#
# PURPOSE:
#   Tests for JAGS fitting functions including JAGS_add_priors, JAGS_get_inits,
#   JAGS_to_monitor, JAGS_check_convergence, JAGS_extend, and related utilities.
#
# DEPENDENCIES:
#   - rjags: For JAGS model syntax generation and testing
#   - common-functions.R: REFERENCE_DIR, test_reference_text, skip_if_no_fits
#
# SKIP CONDITIONS:
#   - skip_if_not_installed("rjags"): For all tests
#   - skip_if_no_fits(): For tests using pre-fitted models
#
# MODELS/FIXTURES:
#   - Some tests use pre-fitted models from test-00-model-fits.R
#
# TAGS: @evaluation, @JAGS
# ============================================================================ #

# Reference directory for text output comparisons
REFERENCE_DIR <<- testthat::test_path("..", "results", "JAGS-fit")

source(testthat::test_path("common-functions.R"))


# ============================================================================ #
# SECTION 1: JAGS_add_priors tests
# ============================================================================ #
test_that("JAGS_add_priors handles various prior types", {

  skip_if_not_installed("rjags")

  # Test with simple priors
  syntax_simple <- "model{}"
  priors_simple <- list(
    mu = prior("normal", list(0, 1)),
    sigma = prior("gamma", list(2, 1))
  )

  result_simple <- JAGS_add_priors(syntax_simple, priors_simple)
  test_reference_text(result_simple, "JAGS_add_priors_simple.txt")

  # Test with truncated priors
  priors_truncated <- list(
    mu = prior("normal", list(0, 1), list(0, Inf))
  )

  result_truncated <- JAGS_add_priors(syntax_simple, priors_truncated)
  test_reference_text(result_truncated, "JAGS_add_priors_truncated.txt")

  # Test with point prior
  priors_point <- list(
    mu = prior("point", list(0))
  )

  result_point <- JAGS_add_priors(syntax_simple, priors_point)
  test_reference_text(result_point, "JAGS_add_priors_point.txt")

  # Test with factor priors
  priors_factor <- list(
    p1 = prior_factor("mnorm", list(mean = 0, sd = 1), contrast = "orthonormal")
  )
  attr(priors_factor[[1]], "levels") <- 3

  result_factor <- JAGS_add_priors(syntax_simple, priors_factor)
  test_reference_text(result_factor, "JAGS_add_priors_factor.txt")

  # Test with weightfunction priors
  priors_wf <- list(
    omega = prior_weightfunction("one-sided", c(0.05), wf_cumulative(c(1, 1)))
  )

  result_wf <- JAGS_add_priors(syntax_simple, priors_wf)
  test_reference_text(result_wf, "JAGS_add_priors_weightfunction.txt")

})


test_that("JAGS_add_priors handles spike_and_slab priors", {

  skip_if_not_installed("rjags")

  priors_sas <- list(
    mu = prior_spike_and_slab(
      prior("normal", list(0, 1)),
      prior_inclusion = prior("beta", list(1, 1))
    )
  )

  result <- JAGS_add_priors("model{}", priors_sas)
  expect_true(grepl("mu_variable", result))
  expect_true(grepl("mu_inclusion", result))
  expect_true(grepl("mu_indicator", result))

  # Test inits
  inits <- JAGS_get_inits(priors_sas, chains = 2, seed = 1)
  expect_true("mu_variable" %in% names(inits[[1]]) || "mu_inclusion" %in% names(inits[[1]]))

  # Test monitor
  monitor <- JAGS_to_monitor(priors_sas)
  expect_true("mu_indicator" %in% monitor)

})


test_that("JAGS_add_priors handles standard prior_mixture (non-bias)", {

  skip_if_not_installed("rjags")

  # Standard mixture (not bias mixture)
  mix <- prior_mixture(list(
    prior("normal", list(0, 0.5)),
    prior("normal", list(0, 1))
  ), is_null = c(TRUE, FALSE))

  priors_mix <- list(mu = mix)

  result <- JAGS_add_priors("model{}", priors_mix)
  expect_true(grepl("mu_indicator", result))
  expect_true(grepl("mu_component_1", result))
  expect_true(grepl("mu_component_2", result))

  # Test inits
  inits <- JAGS_get_inits(priors_mix, chains = 2, seed = 1)
  expect_true("mu_indicator" %in% names(inits[[1]]))

  # Test monitor
  monitor <- JAGS_to_monitor(priors_mix)
  expect_true("mu_indicator" %in% monitor)
  expect_true("mu" %in% monitor)

})


test_that("JAGS_add_priors handles mixture with PEESE prior", {

  skip_if_not_installed("rjags")

  # Create a bias mixture with PEESE prior
  bias_mixture <- prior_mixture(list(
    prior_none(prior_weights = 1),
    prior_PEESE("normal", list(0, 1), prior_weights = 1)
  ))

  priors_peese <- list(
    bias = bias_mixture
  )

  result_peese <- JAGS_add_priors("model{}", priors_peese)
  test_reference_text(result_peese, "JAGS_add_priors_peese_mixture.txt")

})


test_that("JAGS_add_priors handles mixture with PET prior", {

  skip_if_not_installed("rjags")

  # Create a bias mixture with PET prior
  bias_mixture <- prior_mixture(list(
    prior_none(prior_weights = 1),
    prior_PET("normal", list(0, 1), prior_weights = 1)
  ))

  priors_pet <- list(
    bias = bias_mixture
  )

  result_pet <- JAGS_add_priors("model{}", priors_pet)
  test_reference_text(result_pet, "JAGS_add_priors_pet_mixture.txt")

})


# ============================================================================ #
# SECTION 2: JAGS_get_inits tests
# ============================================================================ #
test_that("JAGS_get_inits handles various prior types", {

  skip_if_not_installed("rjags")

  # Test with simple priors
  priors_simple <- list(
    mu = prior("normal", list(0, 1)),
    sigma = prior("gamma", list(2, 1))
  )

  inits1 <- JAGS_get_inits(priors_simple, chains = 2, seed = 1)
  expect_equal(length(inits1), 2)
  expect_true("mu" %in% names(inits1[[1]]))
  expect_true("sigma" %in% names(inits1[[1]]))

  # Same seed should give same results
  inits2 <- JAGS_get_inits(priors_simple, chains = 2, seed = 1)
  expect_equal(inits1, inits2)

  # Different seeds should give different results
  inits3 <- JAGS_get_inits(priors_simple, chains = 2, seed = 123)
  expect_false(isTRUE(all.equal(inits1, inits3)))

  # Test with truncated priors
  priors_truncated <- list(
    mu = prior("normal", list(0, 1), list(0, Inf))
  )

  inits_truncated <- JAGS_get_inits(priors_truncated, chains = 2, seed = 1)
  expect_true(all(sapply(inits_truncated, function(i) i$mu >= 0)))

  # Test with point prior
  priors_point <- list(
    mu = prior("point", list(5))
  )

  inits_point <- JAGS_get_inits(priors_point, chains = 2, seed = 1)
  # Point priors should not generate inits (they're fixed)
  expect_true(!("mu" %in% names(inits_point[[1]])) || all(sapply(inits_point, function(i) i$mu == 5)))

  # Test with factor priors
  priors_factor <- list(
    p1 = prior_factor("mnorm", list(mean = 0, sd = 1), contrast = "orthonormal")
  )
  attr(priors_factor[[1]], "levels") <- 3

  inits_factor <- JAGS_get_inits(priors_factor, chains = 2, seed = 1)
  expect_true("p1" %in% names(inits_factor[[1]]))

})


# ============================================================================ #
# SECTION 3: JAGS_check_convergence tests
# ============================================================================ #
test_that("JAGS_check_convergence works with fitted models", {

  skip_if_not_installed("rjags")
  skip_if_no_fits()

  fit_simple <- readRDS(file.path(temp_fits_dir, "fit_simple_normal.RDS"))
  prior_list <- attr(fit_simple, "prior_list")

  # Test convergence check with prior_list
  convergence <- JAGS_check_convergence(fit_simple, prior_list = prior_list)
  expect_true(is.logical(convergence) || is.list(convergence))

  # Test with NULL prior_list
  convergence_null <- JAGS_check_convergence(fit_simple, prior_list = NULL)
  expect_true(is.logical(convergence_null) || is.list(convergence_null))

})

.fit_convergence_roles <- function(fit){

  coordinates <- parameter_coordinates(fit)
  stats::setNames(coordinates$convergence_role, coordinates$coordinate_name)
}

test_that("JAGS_fit declares monitored fully observed data structural", {

  skip_if_not_installed("rjags")
  skip_if_missing_fits("fit_convergence_observed_data")

  # 'N' and 'x' are data that the syntax reads but never defines. JAGS_fit()
  # passes the fully observed data: partly observed 'y' stays sampled.
  real <- readRDS(file.path(temp_fits_dir, "fit_convergence_observed_data.RDS"))
  real_roles <- .fit_convergence_roles(real)
  expect_true(all(real_roles[c("N", paste0("x[", 1:10, "]"))] == "structural"))
  expect_true(all(real_roles[paste0("y[", 1:10, "]")] == "sampled"))
})

test_that("JAGS LKJ correlation diagonals are exactly 1 and structural in convergence checks", {

  skip_if_not_installed("rjags")
  skip_if_missing_fits(c("fit_lkj_diagonal_K2", "fit_lkj_diagonal_K3"))

  # bt_lkj_corr() returned the diagonal as rounded row sums of squares
  # (1 +/- 2.2e-16), so the monitored R[k,k] had undefined ESS, were "not
  # assessable", and every fit with an unstructured block failed its check.
  fits <- list(
    K2 = readRDS(file.path(temp_fits_dir, "fit_lkj_diagonal_K2.RDS")),
    K3 = readRDS(file.path(temp_fits_dir, "fit_lkj_diagonal_K3.RDS"))
  )

  for(K in 2:3){
    fit <- fits[[paste0("K", K)]]
    samples <- as.matrix(fit$mcmc)
    stem <- attr(fit, "formula_design")$mu$random_effects[[1L]]$parameter_stem
    cell <- function(matrix_name, row, column){
      samples[, paste0(stem, "_xRE_CORx_", matrix_name, "[", row, ",", column, "]")]
    }
    R_diagonal <- paste0(stem, "_xRE_CORx_R[", seq_len(K), ",", seq_len(K), "]")

    # Every monitored diagonal entry is exactly 1 in every draw; the
    # off-diagonal entries are the correlations of the Cholesky factor.
    expect_true(all(samples[, R_diagonal] == 1))
    for(row in 2:K){
      for(column in seq_len(row - 1L)){
        from_L <- rowSums(vapply(seq_len(column), function(m){
          cell("L", row, m) * cell("L", column, m)
        }, numeric(nrow(samples))))
        expect_equal(cell("R", row, column), from_L, tolerance = 1e-12)
        expect_identical(cell("R", row, column), cell("R", column, row))
      }
    }

    check <- JAGS_check_convergence(
      fit, attr(fit, "prior_list"),
      max_Rhat = 2, min_ESS = 1, max_error = NULL, max_SD_error = NULL
    )
    diagnostics <- attr(check, "diagnostics")
    expect_equal(
      diagnostics$state[match(R_diagonal, diagnostics$parameter)],
      rep("structural_constant", K)
    )
    expect_false(any(grepl("_xRE_CORx_R[", attr(check, "errors"), fixed = TRUE)))
    expect_true(check)
  }
})


# ============================================================================ #
# SECTION 4: JAGS_to_monitor tests
# ============================================================================ #
test_that("JAGS_to_monitor generates correct monitor strings", {

  skip_if_not_installed("rjags")

  # Test with simple priors
  priors_simple <- list(
    mu = prior("normal", list(0, 1)),
    sigma = prior("gamma", list(2, 1))
  )

  monitor <- JAGS_to_monitor(priors_simple)
  test_reference_text(paste(sort(monitor), collapse = ","), "JAGS_to_monitor_simple.txt")

  # Test with point prior
  priors_with_point <- list(
    mu = prior("normal", list(0, 1)),
    fixed = prior("point", list(0))
  )

  monitor_point <- JAGS_to_monitor(priors_with_point)
  expect_equal(sort(monitor_point), c("fixed", "mu"))
  test_reference_text(paste(sort(monitor_point), collapse = ", "), "JAGS_to_monitor_point.txt")

  monitor_point_only <- JAGS_to_monitor(list(fixed = prior("point", list(0))))
  expect_equal(monitor_point_only, "fixed")

  monitor_mpoint <- JAGS_to_monitor(list(
    fixed_vector = prior("mpoint", list(1, 2))
  ))
  expect_equal(monitor_mpoint, "fixed_vector")

  monitor_factor_point <- JAGS_to_monitor(list(
    fixed_factor = prior_factor("point", list(0), contrast = "treatment")
  ))
  expect_equal(monitor_factor_point, "fixed_factor")

  # Test with factor priors
  priors_factor <- list(
    p1 = prior_factor("mnorm", list(mean = 0, sd = 1), contrast = "orthonormal")
  )
  attr(priors_factor[[1]], "levels") <- 3

  monitor_factor <- JAGS_to_monitor(priors_factor)
  test_reference_text(paste(sort(monitor_factor), collapse = ","), "JAGS_to_monitor_factor.txt")

})


# ============================================================================ #
# SECTION 5: JAGS_fit attribute preservation
# ============================================================================ #
test_that("JAGS_fit preserves attributes", {

  skip_if_not_installed("rjags")
  skip_on_cran()
  skip_if_no_fits()

  fit_simple <- readRDS(file.path(temp_fits_dir, "fit_simple_normal.RDS"))

  # Check that prior_list attribute is preserved
  prior_list <- attr(fit_simple, "prior_list")
  expect_true(!is.null(prior_list))
  expect_true(is.list(prior_list))

  # Check class
  expect_true(inherits(fit_simple, "BayesTools_fit") || inherits(fit_simple, "runjags"))

})

test_that("functions reading fitted metadata refuse fits without the current contract", {

  skip_if_not_installed("rjags")
  skip_if_missing_fits("fit_simple_thin")
  fit <- readRDS(file.path(temp_fits_dir, "fit_simple_thin.RDS"))
  models <- function(object){
    list(list(fit = object, marglik = bridgesampling_object(0), prior_weights = 1))
  }
  calls <- list(
    JAGS_check_convergence = function(object){
      JAGS_check_convergence(object)
    },
    JAGS_diagnostics = function(object){
      JAGS_diagnostics(object, "mu", type = "trace", plot_type = "ggplot")
    },
    as_mixed_posteriors = function(object){
      as_mixed_posteriors(object, "mu")
    },
    mix_posteriors = function(object){
      mix_posteriors(models(object), "mu", list(mu = FALSE))
    },
    JAGS_extend = function(object){
      JAGS_extend(object, autofit_control = list(max_extend = 1, sample_extend = 10))
    }
  )

  # The metadata state of a BayesTools 0.3.0 fit: no parameter map, contract,
  # or draw geometry. The message is that of the summary tables.
  stripped <- fit
  for(name in c("parameter_map", "fit_contract", "draw_geometry")){
    attr(stripped, name) <- NULL
  }
  missing_map <- paste0(
    "The fitted object does not contain parameter-map metadata. ",
    "Refit the model with the current BayesTools version."
  )
  expect_error(runjags_estimates_table(stripped), missing_map, fixed = TRUE,
               class = "BayesTools_refit_required")

  without_contract <- fit
  attr(without_contract, "fit_contract") <- NULL
  missing_contract <- paste0(
    "The fitted object does not contain a supported schema contract. ",
    "Refit the model with this version of BayesTools."
  )

  previous_map <- fit
  map <- attr(previous_map, "parameter_map", exact = TRUE)
  map$schema_version <- map$schema_version - 1L
  attr(previous_map, "parameter_map") <- map
  unsupported_map <- paste0(
    "Parameter-map metadata are missing, malformed, or unsupported. ",
    "Refit the model with the current BayesTools version."
  )

  for(name in names(calls)){
    expect_error(calls[[name]](stripped), missing_map, fixed = TRUE,
                 class = "BayesTools_refit_required", info = name)
    expect_error(calls[[name]](without_contract), missing_contract, fixed = TRUE,
                 class = "BayesTools_refit_required", info = name)
    expect_error(calls[[name]](previous_map), unsupported_map, fixed = TRUE,
                 class = "BayesTools_refit_required", info = name)
  }

  # Plain runjags objects carry no fitted metadata at all.
  plain <- fit
  class(plain) <- "runjags"
  expect_error(
    JAGS_check_convergence(plain),
    "'fit' must be a 'BayesTools_fit' created by JAGS_fit(). Refit the model with this version of BayesTools.",
    fixed = TRUE,
    class = "BayesTools_refit_required"
  )
  expect_error(
    mix_posteriors(models(plain), "mu", list(mu = FALSE)),
    "'model_list:fit' must be a 'BayesTools_fit' created by JAGS_fit(). Refit the model with this version of BayesTools.",
    fixed = TRUE,
    class = "BayesTools_refit_required"
  )

  # Fits of this version pass.
  expect_true(JAGS_check_convergence(
    fit, max_Rhat = 2, min_ESS = 1, max_error = NULL, max_SD_error = NULL
  ))
  expect_s3_class(as_mixed_posteriors(fit, "mu")$mu, "mixed_posteriors")
})


# ============================================================================ #
# SECTION 6: runjags_estimates_table tests
# ============================================================================ #
test_that("runjags_estimates_table works with fitted models", {

  skip_if_not_installed("rjags")
  skip_on_cran()
  skip_if_no_fits()

  fit_simple <- readRDS(file.path(temp_fits_dir, "fit_simple_normal.RDS"))
  fit_draws <- do.call(rbind, lapply(fit_simple$mcmc, as.matrix))
  expect_estimates_from_fit <- function(table){
    parameters <- attr(table, "parameters")
    probs <- c(0.025, 0.5, 0.975)
    estimate_names <- c("Mean", "SD", as.character(probs))
    expected <- t(vapply(parameters, function(parameter){
      draws <- fit_draws[, parameter]
      c(
        Mean = mean(draws, na.rm = TRUE),
        SD = stats::sd(draws, na.rm = TRUE),
        vapply(
          probs,
          function(prob) unname(stats::quantile(
            draws,
            probs = prob,
            na.rm = TRUE
          )),
          numeric(1)
        )
      )
    }, numeric(length(estimate_names))))
    colnames(expected) <- estimate_names

    expect_equal(
      unname(as.matrix(table[, estimate_names, drop = FALSE])),
      unname(expected),
      tolerance = 1e-12
    )
  }

  # Test basic estimates table
  estimates_table <- runjags_estimates_table(fit_simple)
  test_reference_table_stochastic(
    estimates_table,
    "runjags_estimates_simple.txt"
  )
  expect_identical(attr(estimates_table, "parameters"), c("m", "s"))
  expect_estimates_from_fit(estimates_table)

  # Test without specific parameters
  estimates_table_param <- runjags_estimates_table(fit_simple, remove_parameters = "m")
  test_reference_table_stochastic(
    estimates_table_param,
    "runjags_estimates_param_m.txt"
  )
  expect_identical(attr(estimates_table_param, "parameters"), "s")
  expect_estimates_from_fit(estimates_table_param)

})


# ============================================================================ #
# SECTION 7: JAGS_extend tests
# ============================================================================ #


# ============================================================================ #
# SECTION 8: JAGS handles specific prior types
# ============================================================================ #
test_that("JAGS handles invgamma prior", {

  skip_if_not_installed("rjags")

  priors_inv <- list(tau = prior("invgamma", list(3, 2)))

  # Test syntax
  result <- JAGS_add_priors("model{}", priors_inv)
  expect_true(grepl("tau ~ dbt_invgamma(3,2)", result, fixed = TRUE))
  expect_false(grepl("inv_tau", result, fixed = TRUE))
  expect_false(grepl("pow(inv_tau", result, fixed = TRUE))

  # Test inits
  inits <- JAGS_get_inits(priors_inv, chains = 2, seed = 1)
  expect_true("tau" %in% names(inits[[1]]))
  expect_false("inv_tau" %in% names(inits[[1]]))

  # Test monitor
  monitor <- JAGS_to_monitor(priors_inv)
  expect_equal(monitor, "tau")

})


test_that("JAGS handles independent weightfunction priors", {

  skip_if_not_installed("rjags")

  priors_wf2 <- list(omega = prior_weightfunction("one-sided", c(0.05, 0.60), wf_independent(prior("beta", list(1, 1)))))

  # Test syntax
  result <- JAGS_add_priors("model{}", priors_wf2)
  expect_true(grepl("omega\\[2\\] ~ dbeta", result))
  expect_true(grepl("omega\\[3\\] ~ dbeta", result))

  # Test inits
  inits <- JAGS_get_inits(priors_wf2, chains = 2, seed = 1)
  expect_false("eta" %in% names(inits[[1]]))

  # Test monitor
  monitor <- JAGS_to_monitor(priors_wf2)
  expect_true("omega" %in% monitor)

})


test_that("JAGS handles weightfunction fixed prior", {

  skip_if_not_installed("rjags")

  priors_wf_fixed <- list(omega = prior_weightfunction("one-sided", c(0.05), wf_fixed(c(1, 0.5))))

  # Test syntax - fixed weightfunction has no eta parameters to sample
  result <- JAGS_add_priors("model{}", priors_wf_fixed)
  expect_true(grepl("omega", result))

  # Test inits - fixed weightfunction should return empty inits for eta
  inits <- JAGS_get_inits(priors_wf_fixed, chains = 2, seed = 1)
  # Should not have eta since it's fixed
  expect_true(!("eta" %in% names(inits[[1]])))

  # Test monitor
  monitor <- JAGS_to_monitor(priors_wf_fixed)
  expect_true("omega" %in% monitor)

})


test_that("JAGS handles factor treatment/independent priors", {

  skip_if_not_installed("rjags")

  # Treatment contrast
  prior_treat <- prior_factor("normal", list(0, 1), contrast = "treatment")
  attr(prior_treat, "levels") <- 3

  priors_treat <- list(fac = prior_treat)
  result_treat <- JAGS_add_priors("model{}", priors_treat)
  expect_true(grepl("fac\\[i\\]", result_treat))

  # Independent contrast
  prior_indep <- prior_factor("gamma", list(2, 1), contrast = "independent")
  attr(prior_indep, "levels") <- 2

  priors_indep <- list(fac = prior_indep)
  result_indep <- JAGS_add_priors("model{}", priors_indep)
  expect_true(grepl("dgamma", result_indep))

})


test_that("JAGS handles vector mt prior", {

  skip_if_not_installed("rjags")

  prior_mt <- prior("mt", list(location = 0, scale = 1, df = 5, K = 2))
  priors_mt <- list(p = prior_mt)

  # Test syntax
  result <- JAGS_add_priors("model{}", priors_mt)
  expect_true(grepl("prior_par_s_p", result))
  expect_true(grepl("prior_par_z_p", result))

  # Test inits
  inits <- JAGS_get_inits(priors_mt, chains = 2, seed = 1)
  expect_true("prior_par_s_p" %in% names(inits[[1]]))
  expect_true("prior_par_z_p" %in% names(inits[[1]]))

})


test_that("JAGS handles bias mixture with weightfunction", {

  skip_if_not_installed("rjags")

  bias_mix_wf <- prior_mixture(list(
    prior_none(prior_weights = 1),
    prior_weightfunction("one-sided", c(0.05), wf_cumulative(c(1, 1)), prior_weights = 1)
  ))

  priors_bias_wf <- list(bias = bias_mix_wf)

  result <- JAGS_add_priors("model{}", priors_bias_wf)
  expect_true(grepl("bias_indicator", result))
  expect_true(grepl("omega", result))
  expect_true(grepl("omega_ratio_component_2", result, fixed = TRUE))

  # Test inits
  inits <- JAGS_get_inits(priors_bias_wf, chains = 2, seed = 1)
  expect_true("bias_indicator" %in% names(inits[[1]]))

  # Test monitor
  monitor <- JAGS_to_monitor(priors_bias_wf)
  expect_true("bias_indicator" %in% monitor)
  expect_true("omega" %in% monitor)

})


# ============================================================================ #
# SECTION 9: JAGS_check_and_list_autofit_settings
# ============================================================================ #
test_that("JAGS_check_and_list_autofit_settings validates all parameters", {

  # Valid settings
  valid_settings <- list(
    max_Rhat = 1.05,
    min_ESS = 500,
    max_error = 0.01,
    max_SD_error = 0.05,
    max_time = list(time = 1, unit = "mins"),
    sample_extend = 100,
    restarts = 3,
    max_extend = 10
  )
  expect_silent(JAGS_check_and_list_autofit_settings(valid_settings))

  # max_time without names - should auto-assign
  unnamed_time <- list(
    max_Rhat = 1.05, min_ESS = 500, max_error = 0.01, max_SD_error = 0.05,
    max_time = list(1, "mins"), sample_extend = 100
  )
  expect_silent(JAGS_check_and_list_autofit_settings(unnamed_time))

})
