skip_if_not_test_profile("fit")

# ============================================================================ #
# TEST FILE: JAGS Fit Edge Cases
# ============================================================================ #
#
# PURPOSE:
#   Edge case tests for JAGS fitting functions including input validation,
#   error handling, and boundary conditions.
#
# DEPENDENCIES:
#   - rjags: For JAGS model syntax generation and testing
#   - common-functions.R: REFERENCE_DIR, test_reference_text, skip_if_no_fits
#
# SKIP CONDITIONS:
#   - skip_if_not_installed("rjags"): For all tests
#
# MODELS/FIXTURES:
#   - Some tests use pre-fitted models from test-00-model-fits.R
#
# TAGS: @edge-cases, @JAGS, @input-validation
# ============================================================================ #

# Reference directory for text output comparisons
REFERENCE_DIR <<- testthat::test_path("..", "results", "JAGS-fit-edge-cases")

source(testthat::test_path("common-functions.R"))


# ============================================================================ #
# SECTION 1: Input validation tests
# ============================================================================ #
test_that("JAGS_add_priors input validation works", {

  # Empty prior_list returns original syntax
  expect_equal(JAGS_add_priors("model{}", list()), "model{}")

  # prior_list must be a list of priors
  expect_error(JAGS_add_priors("model{}", list(x = 1)), "'prior_list' must be a list of priors.")
  expect_error(JAGS_add_priors("model{}", prior("normal", list(0, 1))), "'prior_list' must be a list of priors.")

})


test_that("JAGS_get_inits input validation works", {

  # Empty prior_list returns empty list
  expect_equal(JAGS_get_inits(list(), chains = 2, seed = 1), list())

  # Input validation
  expect_error(JAGS_get_inits(list(x = 1), chains = 2, seed = 1), "'prior_list' must be a list of priors.")
  expect_error(JAGS_get_inits(prior("normal", list(0, 1)), chains = 2, seed = 1), "'prior_list' must be a list of priors.")

})


test_that("JAGS_to_monitor input validation works", {

  # Empty prior_list returns empty string
  expect_equal(JAGS_to_monitor(list()), "")

  # Input validation
  expect_error(JAGS_to_monitor(list(x = 1)), "'prior_list' must be a list of priors.")
  expect_error(JAGS_to_monitor(prior("normal", list(0, 1))), "'prior_list' must be a list of priors.")

})


test_that(".check_JAGS_syntax validates syntax correctly", {

  # Test with valid syntax
  expect_silent(JAGS_add_priors("model{}", list(mu = prior("normal", list(0, 1)))))

  # Test with missing "model" keyword
  expect_error(
    JAGS_add_priors("invalid{}", list(mu = prior("normal", list(0, 1)))),
    "syntax must be a JAGS model syntax"
  )

  # Test with missing opening brace
  expect_error(
    JAGS_add_priors("model}", list(mu = prior("normal", list(0, 1)))),
    "syntax must be a JAGS model syntax"
  )

  # Test with missing closing brace
  expect_error(
    JAGS_add_priors("model{", list(mu = prior("normal", list(0, 1)))),
    "syntax must be a JAGS model syntax"
  )

  # Test with non-character input
  expect_error(
    JAGS_add_priors(123, list(mu = prior("normal", list(0, 1)))),
    "must be a character"
  )

})


test_that("JAGS_extend error handling", {

  skip_if_not_installed("rjags")
  skip_on_cran()

  # Test error when fit is not a BayesTools_fit
  expect_error(
    JAGS_extend(list(), autofit_control = list()),
    "'fit' must be a 'BayesTools_fit'"
  )

})


test_that("required packages are checked locally and on parallel workers", {

  expect_silent({
    local_loaded <- .JAGS_require_packages("stats")
  })
  expect_equal(unname(local_loaded), TRUE)
  expect_error(
    .JAGS_require_packages("BayesToolsMissingPackageForTest"),
    "Required packages are not available: 'BayesToolsMissingPackageForTest'.",
    fixed = TRUE
  )
  expect_error(
    JAGS_fit("model{}", required_packages = NA_character_),
    "The 'required_packages' argument cannot contain NA/NaN values.",
    fixed = TRUE
  )

  cl <- parallel::makePSOCKcluster(2)
  on.exit(parallel::stopCluster(cl), add = TRUE)

  expect_silent({
    worker_loaded <- .JAGS_require_packages("stats", cl)
  })
  expect_equal(unname(worker_loaded), TRUE)
  expect_error(
    .JAGS_require_packages(c("stats", "BayesToolsMissingPackageForTest"), cl),
    "Required packages are not available: 'BayesToolsMissingPackageForTest'.",
    fixed = TRUE
  )

})

test_that("JAGS_fit stores formula design metadata on fitted formula models", {

  skip_if_not_installed("runjags")
  skip_if_not_installed("rjags")
  skip_on_cran()

  df <- data.frame(
    y = c(-0.2, 0.1, 0.4, 0.7),
    x = c(-1, 0, 1, 2)
  )

  fit <- JAGS_fit(
    model_syntax = "model{
      for(i in 1:N_mu){
        y[i] ~ dnorm(mu[i], 1)
      }
    }",
    data = list(y = df$y),
    formula_list = list(mu = ~ x),
    formula_data_list = list(mu = df),
    formula_prior_list = list(mu = list(
      intercept = prior("normal", list(0, 1)),
      x         = prior("normal", list(0, 1))
    )),
    chains = 1,
    adapt = 50,
    burnin = 50,
    sample = 100,
    silent = TRUE,
    seed = 123
  )

  design <- JAGS_formula_design(fit, "mu")

  expect_s3_class(fit, "BayesTools_fit")
  expect_s3_class(design, "BayesTools_formula_design")
  expect_equal(colnames(design$model_matrix), c("(Intercept)", "x"))
  expect_equal(unname(design$model_matrix[, "x"]), df$x)
  expect_equal(design$jags_data_names$x, "mu_data_x")
  expect_identical(JAGS_formula_design(fit), list(mu = design))
})

test_that("JAGS_fit runs dummy structured random-effect formula models", {

  skip_if_not_installed("runjags")
  skip_if_not_installed("rjags")
  skip_on_cran()

  sd_prior <- prior("normal", list(0, 0.5), truncation = list(lower = 0, upper = Inf))
  rho_prior <- prior("normal", list(0, 0.5))
  model_syntax <- "model{
    for(i in 1:N_mu){
      y[i] ~ dnorm(mu[i], 4)
    }
  }"

  df_cs <- data.frame(
    y = c(-0.2, 0.1, 0.4, 0.7, -0.1, 0.3, 0.5, 0.9),
    x = c(-1, 0, 1, 2, -2, 0.5, 1.5, 2.5),
    idx = factor(rep(c("t1", "t2"), 4), levels = c("t1", "t2")),
    id = factor(rep(c("g1", "g2", "g3", "g4"), each = 2))
  )

  fit_cs <- suppressWarnings(JAGS_fit(
    model_syntax = model_syntax,
    data = list(y = df_cs$y),
    formula_list = list(mu = ~ 1 + x + cs(idx | id)),
    formula_data_list = list(mu = df_cs),
    formula_prior_list = list(mu = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    )),
    formula_random_prior_list = list(mu = prior_random(
      id = random_block(
        sd = sd_prior,
        rho = rho_prior,
        monitor = random_monitor(latent = FALSE, coefficients = FALSE, correlation = TRUE)
      )
    )),
    chains = 1,
    adapt = 50,
    burnin = 50,
    sample = 100,
    silent = TRUE,
    seed = 101
  ))

  expect_s3_class(fit_cs, "BayesTools_fit")
  expect_equal(attr(fit_cs, "jags_modules"), character())
  expect_equal(JAGS_formula_design(fit_cs, "mu")$random_effects[[1]]$structure, "cs")
  expect_true("mu__xREx__id_rho" %in% colnames(as.matrix(fit_cs$mcmc)))

  df_ar1 <- data.frame(
    y = c(-0.2, 0.1, 0.4, 0.7, -0.1, 0.3, 0.5, 0.9, 0.2, 0.6, 0.8, 1.0),
    f = factor(rep(c("a", "b", "c"), 4), levels = c("a", "b", "c")),
    id = factor(rep(c("g1", "g2", "g3", "g4"), each = 3))
  )

  fit_ar1 <- suppressWarnings(JAGS_fit(
    model_syntax = model_syntax,
    data = list(y = df_ar1$y),
    formula_list = list(mu = ~ 1 + ar1(f | id)),
    formula_data_list = list(mu = df_ar1),
    formula_prior_list = list(mu = list(
      intercept = prior("normal", list(0, 1))
    )),
    formula_random_prior_list = list(mu = prior_random(
      id = random_block(
        sd = sd_prior,
        rho = rho_prior,
        monitor = random_monitor(latent = FALSE, coefficients = FALSE, correlation = TRUE)
      )
    )),
    chains = 1,
    adapt = 50,
    burnin = 50,
    sample = 100,
    silent = TRUE,
    seed = 102
  ))

  expect_s3_class(fit_ar1, "BayesTools_fit")
  ar1_design <- JAGS_formula_design(fit_ar1, "mu")$random_effects[[1]]
  expect_equal(ar1_design$structure, "ar1")
  expect_equal(ar1_design$column_names, c("fa", "fb", "fc"))
  expect_true("mu__xREx__id_rho" %in% colnames(as.matrix(fit_ar1$mcmc)))

  df_us <- data.frame(
    y = c(-0.2, 0.1, 0.4, 0.7, -0.1, 0.3),
    x = c(-1, 0, 1, 2, -2, 0.5),
    id = factor(c("g1", "g1", "g2", "g2", "g3", "g3"))
  )

  fit_us <- suppressWarnings(JAGS_fit(
    model_syntax = model_syntax,
    data = list(y = df_us$y),
    formula_list = list(mu = ~ 1 + x + (1 + x | id)),
    formula_data_list = list(mu = df_us),
    formula_prior_list = list(mu = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    )),
    formula_random_prior_list = list(mu = prior_random(
      id = random_block(
        sd = sd_prior,
        cor = prior_lkj(eta = 1, include_correlation = FALSE)
      )
    )),
    chains = 1,
    adapt = 50,
    burnin = 50,
    sample = 100,
    silent = TRUE,
    seed = 103
  ))

  expect_s3_class(fit_us, "BayesTools_fit")
  expect_equal(attr(fit_us, "jags_modules"), "BayesTools")
  expect_equal(JAGS_formula_design(fit_us, "mu")$random_effects[[1]]$structure, "us")
  expect_true(any(grepl("mu__xREx__id_xRE_CORx_L", colnames(as.matrix(fit_us$mcmc)), fixed = TRUE)))
})

test_that("JAGS_fit monitors random coefficients for observed-level prediction", {

  skip_if_not_installed("runjags")
  skip_if_not_installed("rjags")
  skip_on_cran()

  df <- data.frame(
    y = c(-0.2, 0.1, 0.4, 0.7, -0.1, 0.3),
    id = factor(c("g1", "g1", "g2", "g2", "g3", "g3"))
  )

  fit <- suppressWarnings(JAGS_fit(
    model_syntax = "model{
      for(i in 1:N_mu){
        y[i] ~ dnorm(mu[i], 4)
      }
    }",
    data = list(y = df$y),
    formula_list = list(mu = ~ 1 + diag(1 | id)),
    formula_data_list = list(mu = df),
    formula_prior_list = list(mu = list(
      intercept = prior("normal", list(0, 1))
    )),
    formula_random_prior_list = list(mu = prior_random(
      id = random_block(
        sd = prior("gamma", list(2, 2)),
        monitor = random_monitor(coefficients = TRUE, correlation = FALSE)
      )
    )),
    chains = 1,
    adapt = 50,
    burnin = 50,
    sample = 100,
    silent = TRUE,
    seed = 103
  ))

  posterior <- as.matrix(fit$mcmc)
  expect_true(any(grepl("mu__xREx__id_xRE_COEFx", colnames(posterior), fixed = TRUE)))

  prediction <- JAGS_evaluate_formula(
    fit = fit,
    formula = ~ 1 + diag(1 | id),
    parameter = "mu",
    data = df,
    prior_list = attr(fit, "prior_list")
  )

  expect_equal(dim(prediction), c(nrow(df), nrow(posterior)))
  expect_true(all(is.finite(prediction)))
})

test_that("JAGS_fit predicts observed random effects from latent monitors", {

  skip_if_not_installed("runjags")
  skip_if_not_installed("rjags")
  skip_on_cran()

  df <- data.frame(
    y = c(-0.2, 0.1, 0.4, 0.7, -0.1, 0.3),
    id = factor(c("g1", "g1", "g2", "g2", "g3", "g3"))
  )

  fit <- suppressWarnings(JAGS_fit(
    model_syntax = "model{
      for(i in 1:N_mu){
        y[i] ~ dnorm(mu[i], 4)
      }
    }",
    data = list(y = df$y),
    formula_list = list(mu = ~ 1 + diag(1 | id)),
    formula_data_list = list(mu = df),
    formula_prior_list = list(mu = list(
      intercept = prior("normal", list(0, 1))
    )),
    formula_random_prior_list = list(mu = prior_random(
      id = random_block(
        sd = prior("gamma", list(2, 2)),
        monitor = random_monitor(coefficients = FALSE, correlation = FALSE)
      )
    )),
    chains = 1,
    adapt = 50,
    burnin = 50,
    sample = 100,
    silent = TRUE,
    seed = 104
  ))

  posterior <- as.matrix(fit$mcmc)
  expect_true(any(grepl("mu__xREx__id_xRE_Zx", colnames(posterior), fixed = TRUE)))
  expect_false(any(grepl("mu__xREx__id_xRE_COEFx", colnames(posterior), fixed = TRUE)))

  prediction <- JAGS_evaluate_formula(
    fit = fit,
    formula = ~ 1 + diag(1 | id),
    parameter = "mu",
    data = df,
    prior_list = attr(fit, "prior_list")
  )

  expect_equal(dim(prediction), c(nrow(df), nrow(posterior)))
  expect_true(all(is.finite(prediction)))
})

test_that("JAGS_fit stores formula design metadata on failed sampling objects", {

  skip_if_not_installed("runjags")
  skip_if_not_installed("rjags")
  skip_on_cran()

  df <- data.frame(x = c(-1, 0, 1))

  fit <- suppressWarnings(JAGS_fit(
    model_syntax = "model{
      broken_node <- missing_node
    }",
    formula_list = list(mu = ~ x),
    formula_data_list = list(mu = df),
    formula_prior_list = list(mu = list(
      intercept = prior("normal", list(0, 1)),
      x         = prior("normal", list(0, 1))
    )),
    chains = 1,
    adapt = 50,
    burnin = 50,
    sample = 100,
    silent = TRUE,
    seed = 123
  ))

  design <- JAGS_formula_design(fit, "mu")

  expect_s3_class(fit, "BayesTools_fit")
  expect_s3_class(fit, "error")
  expect_s3_class(design, "BayesTools_formula_design")
  expect_equal(unname(design$model_matrix[, "x"]), df$x)
})


# ============================================================================ #
# SECTION 2: Convergence edge cases
# ============================================================================ #
test_that("autofit settings keep indicator checks off by default", {

  settings <- JAGS_check_and_list_autofit_settings(list(
    max_Rhat = 1.05,
    min_ESS = 500,
    max_error = 0.01,
    max_SD_error = 0.05,
    max_time = list(time = 60, unit = "mins"),
    sample_extend = 1000,
    restarts = 10,
    max_extend = 10
  ))
  expect_false(settings$check_indicators)

  settings <- JAGS_check_and_list_autofit_settings(list(
    max_Rhat = 1.05,
    min_ESS = 500,
    max_error = 0.01,
    max_SD_error = 0.05,
    max_time = list(time = 60, unit = "mins"),
    sample_extend = 1000,
    restarts = 10,
    max_extend = 10,
    check_indicators = TRUE
  ))
  expect_true(settings$check_indicators)

})


test_that("JAGS_check_convergence ignores indicator variables unless requested", {

  set.seed(1)
  chain_1 <- cbind(mu = rnorm(100), mu_indicator = rep(0, 100))
  chain_2 <- cbind(mu = rnorm(100), mu_indicator = rep(1, 100))
  fit <- list(
    mcmc         = coda::mcmc.list(coda::mcmc(chain_1), coda::mcmc(chain_2)),
    summary.pars = list(mutate = NULL)
  )
  class(fit) <- "runjags"

  prior_list <- list(mu = prior("normal", list(0, 1)))

  expect_true(JAGS_check_convergence(
    fit,
    prior_list       = prior_list,
    max_Rhat        = 1.05,
    min_ESS         = NULL,
    max_error       = NULL,
    max_SD_error    = NULL,
    check_indicators = FALSE
  ))

  with_indicators <- JAGS_check_convergence(
    fit,
    prior_list       = prior_list,
    max_Rhat        = 1.05,
    min_ESS         = NULL,
    max_error       = NULL,
    max_SD_error    = NULL,
    check_indicators = TRUE
  )
  expect_false(with_indicators)
  expect_match(attr(with_indicators, "errors"), "R-hat")

})

test_that("JAGS_check_convergence ignores add_parameters without priors", {

  set.seed(2)
  mu_values <- rnorm(100)
  chain_1 <- cbind(mu = mu_values, "aux[1]" = rep(0, 100))
  chain_2 <- cbind(mu = mu_values, "aux[1]" = rep(1, 100))
  fit <- list(
    mcmc         = coda::mcmc.list(coda::mcmc(chain_1), coda::mcmc(chain_2)),
    summary.pars = list(mutate = NULL)
  )
  class(fit) <- "runjags"

  prior_list <- list(mu = prior("normal", list(0, 1)))

  without_aux <- JAGS_check_convergence(
    fit,
    prior_list    = prior_list,
    max_Rhat     = 1.05,
    min_ESS      = NULL,
    max_error    = NULL,
    max_SD_error = NULL
  )
  expect_false(without_aux)

  expect_true(JAGS_check_convergence(
    fit,
    prior_list     = prior_list,
    max_Rhat      = 1.05,
    min_ESS       = NULL,
    max_error     = NULL,
    max_SD_error  = NULL,
    add_parameters = "aux"
  ))
})


test_that("JAGS_check_convergence handles single chain (R-hat warning)", {

  skip_if_not_installed("rjags")
  skip_on_cran()

  prior_list <- list(mu = prior("normal", list(0, 1)))
  model_syntax <- JAGS_add_priors("model{}", prior_list)
  old_silent.runjags <- runjags::runjags.getOption("silent.runjags")
  on.exit(runjags::runjags.options(silent.runjags = old_silent.runjags), add = TRUE)
  runjags::runjags.options(silent.runjags = TRUE)

  set.seed(1)
  fit <- suppressWarnings(runjags::run.jags(
    model = model_syntax,
    monitor = "mu",
    n.chains = 1,  # Single chain - R-hat cannot be computed
    adapt = 50,
    burnin = 50,
    sample = 100,
    silent.jags = TRUE
  ))

  # Should warn about single chain R-hat
  expect_warning(
    JAGS_check_convergence(fit, prior_list = prior_list, max_Rhat = 1.05),
    "Only one chain was run"
  )

})


test_that("JAGS_check_convergence handles ESS and error checks", {

  skip_if_not_installed("rjags")
  skip_on_cran()

  prior_list <- list(mu = prior("normal", list(0, 1)))
  model_syntax <- JAGS_add_priors("model{}", prior_list)
  old_silent.runjags <- runjags::runjags.getOption("silent.runjags")
  on.exit(runjags::runjags.options(silent.runjags = old_silent.runjags), add = TRUE)
  runjags::runjags.options(silent.runjags = TRUE)

  set.seed(1)
  fit <- suppressWarnings(runjags::run.jags(
    model = model_syntax,
    monitor = "mu",
    n.chains = 2,
    adapt = 50,
    burnin = 50,
    sample = 50,  # Small sample for testing convergence failures
    silent.jags = TRUE
  ))

  # Test with very strict ESS requirement (should fail)
  result_ess <- JAGS_check_convergence(fit, prior_list = prior_list, max_Rhat = NULL, min_ESS = 10000, max_error = NULL, max_SD_error = NULL, fail_fast = FALSE)
  expect_false(result_ess)
  expect_true(!is.null(attr(result_ess, "errors")))

  # Test with very strict error requirement
  result_err <- JAGS_check_convergence(fit, prior_list = prior_list, max_Rhat = NULL, min_ESS = NULL, max_error = 0.00001, max_SD_error = NULL, fail_fast = FALSE)
  expect_false(result_err)

  # Test with very strict SD error requirement
  result_sd <- JAGS_check_convergence(fit, prior_list = prior_list, max_Rhat = NULL, min_ESS = NULL, max_error = NULL, max_SD_error = 0.00001, fail_fast = FALSE)
  expect_false(result_sd)

})


# ============================================================================ #
# SECTION 3: JAGS_fit with is_JASP mode
# ============================================================================ #
test_that("JAGS_fit works with is_JASP mode", {

  skip_if_not_installed("rjags")
  skip_on_cran()

  # Simple model for testing is_JASP mode
  set.seed(1)
  data <- list(
    y = rnorm(20, 0.5, 1),
    N = 20
  )

  prior_list <- list(
    mu    = prior("normal", list(0, 1)),
    sigma = prior("normal", list(0, 1), list(0, Inf))
  )

  model_syntax <- "model{
    for(i in 1:N){
      y[i] ~ dnorm(mu, 1/pow(sigma, 2))
    }
  }"

  # Mock JASP progress bar functions (they should be skipped if not available)
  # The is_JASP mode should work but simply skip progress bars if functions don't exist
  fit_jasp <- capture.output(tryCatch({
    suppressWarnings(JAGS_fit(
      model_syntax = model_syntax,
      data = data,
      prior_list = prior_list,
      chains = 1,
      adapt = 50,
      burnin = 50,
      sample = 100,
      seed = 1,
      silent = TRUE,
      is_JASP = TRUE,
      is_JASP_prefix = "Test"
    ))
  }, error = function(e) {
    # If JASP functions don't exist, this should still produce a fit
    # or fail gracefully
    if (grepl("JASP", e$message)) {
      skip("JASP progress bar functions not available")
    }
    stop(e)
  }))

  test_reference_text(paste0(fit_jasp, collapse = ","), "fit_jasp.txt")

})
