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
  expect_equal(JAGS_add_priors(NULL, list()), "model{}")
  expect_match(JAGS_add_priors(NULL, list(mu = prior("normal", list(0, 1)))), "^model\\{")

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

test_that("JAGS_fit rejects unrecognized formula_scale_list entries", {

  skip_if_not_installed("runjags")

  data <- data.frame(x = c(-1, 0, 1))
  expect_error(
    JAGS_fit(
      model_syntax = "model{}",
      formula_list = list(mu = ~ x),
      formula_data_list = list(mu = data),
      formula_prior_list = list(mu = list(
        intercept = prior("normal", list(0, 1)),
        x = prior("normal", list(0, 1))
      )),
      formula_scale_list = list(muu = list(x = TRUE))
    ),
    "not recognized"
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

test_that("JAGS_fit does not emit unused sampled data for marginalized random effects", {

  skip_if_not_installed("runjags")
  skip_if_not_installed("rjags")
  skip_on_cran()

  df <- data.frame(
    y = c(-0.2, 0.1, 0.3, -0.1),
    study = factor(c("s1", "s1", "s2", "s2")),
    estimate = factor(c("e1", "e2", "e3", "e4"))
  )
  sd_prior <- prior(
    "normal",
    list(mean = 0, sd = 1),
    truncation = list(lower = 0, upper = Inf)
  )
  prior_random <- prior_random(
    study = random_block(sd = sd_prior),
    estimate = random_block(sd = sd_prior)
  )
  warnings <- character()

  fit <- withCallingHandlers(
    JAGS_fit(
      model_syntax = "model{
        for(i in 1:N_mu){
          y[i] ~ dnorm(mu[i], 1)
        }
      }",
      data = list(y = df$y),
      formula_list = list(
        mu = ~ 1 +
          random(1 | study, name = "study", covariance = "diag") +
          random(1 | estimate, name = "estimate", covariance = "diag")
      ),
      formula_data_list = list(mu = df),
      formula_prior_list = list(mu = list(
        intercept = prior("normal", list(mean = 0, sd = 1))
      )),
      formula_random_prior_list = list(mu = prior_random),
      formula_random_effects_compile_list = list(
        mu = random_effects_compile(marginalized = "estimate")
      ),
      chains = 1,
      adapt = 100,
      burnin = 100,
      sample = 100,
      silent = FALSE,
      seed = 123
    ),
    warning = function(w){
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  design <- JAGS_formula_design(fit, "mu")
  estimate_term <- design$random_effects[[which(vapply(
    design$random_effects,
    function(random_term) identical(random_term$block_name, "estimate"),
    logical(1)
  ))]]

  expect_equal(warnings, character())
  expect_equal(estimate_term$compile_mode, "marginalized")
  expect_equal(estimate_term$jags_data_names, character())
  expect_false(grepl("mu__xREx__estimate_xRE_DATAx", attr(fit, "model"), fixed = TRUE))
  expect_false(grepl("mu__xREx__estimate_xRE_MAPx", attr(fit, "model"), fixed = TRUE))
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

test_that("centered and noncentered fits preserve substantive output schemas", {

  skip_if_not_installed("runjags")
  skip_if_not_installed("rjags")
  skip_on_cran()

  df <- data.frame(
    y = c(-0.4, -0.1, 0.2, 0.5, 0.8, 0.1, 0.3, 0.6, 0.9, 1.1,
          -0.3, 0, 0.1, 0.4, 0.7, 0.2, 0.5, 0.7, 1, 1.2),
    index = factor(rep(c("i1", "i2"), 10L)),
    id = factor(rep(paste0("g", seq_len(4L)), each = 5L))
  )
  fit_parameterization <- function(parameterization, seed){

    suppressWarnings(JAGS_fit(
      model_syntax = "model{
        for(i in 1:N_mu){
          y[i] ~ dnorm(mu[i], 4)
        }
      }",
      data = list(y = df$y),
      formula_list = list(mu = ~ 1 + cs(index | id)),
      formula_data_list = list(mu = df),
      formula_prior_list = list(mu = list(
        intercept = prior("normal", list(0, 1))
      )),
      formula_random_prior_list = list(mu = prior_random(
        id = random_block(
          sd = prior("gamma", list(2, 2)),
          rho = prior("normal", list(0, 0.5)),
          parameterization = parameterization
        )
      )),
      chains = 1,
      adapt = 100,
      burnin = 200,
      sample = 400,
      silent = TRUE,
      seed = seed
    ))
  }

  noncentered <- fit_parameterization("noncentered", 207)
  centered    <- fit_parameterization("centered", 208)
  noncentered_draws <- as.matrix(noncentered$mcmc)
  centered_draws    <- as.matrix(centered$mcmc)
  centered_term <- JAGS_formula_design(centered, "mu")$random_effects[[1L]]
  sd_name       <- centered_term$sd_parameter_names[[1L]]
  rho_name      <- centered_term$correlation$rho_name

  expect_setequal(colnames(centered_draws), colnames(noncentered_draws))
  expect_equal(
    mean(centered_draws[, "mu_intercept"]),
    mean(noncentered_draws[, "mu_intercept"]),
    tolerance = 0.2
  )
  expect_equal(
    mean(centered_draws[, sd_name]),
    mean(noncentered_draws[, sd_name]),
    tolerance = 0.25
  )
  expect_equal(
    mean(centered_draws[, rho_name]),
    mean(noncentered_draws[, rho_name]),
    tolerance = 0.25
  )
  expect_identical(centered_term$parameterization_resolved, "centered")
})

test_that("group-local structured fits monitor only active latent cells", {

  skip_if_not_installed("runjags")
  skip_if_not_installed("rjags")
  skip_on_cran()

  K <- 40L
  df <- data.frame(
    y = seq(-0.5, 0.5, length.out = K),
    index = factor(
      paste0("level_", seq_len(K)),
      levels = paste0("level_", seq_len(K))
    ),
    id = factor(rep(paste0("group_", seq_len(8L)), length.out = K))
  )
  fit <- suppressWarnings(JAGS_fit(
    model_syntax = "model{
      for(i in 1:N_mu){
        y[i] ~ dnorm(mu[i], 4)
      }
    }",
    data = list(y = df$y),
    formula_list = list(mu = ~ 1 + cs(index | id)),
    formula_data_list = list(mu = df),
    formula_prior_list = list(mu = list(
      intercept = prior("normal", list(0, 1))
    )),
    formula_random_prior_list = list(mu = prior_random(
      id = random_block(
        sd = prior("point", list(location = 1)),
        rho = prior("point", list(location = 0.2))
      )
    )),
    chains = 1,
    adapt = 50,
    burnin = 50,
    sample = 100,
    silent = TRUE,
    seed = 209
  ))
  term <- JAGS_formula_design(fit, "mu")$random_effects[[1L]]
  posterior <- as.matrix(fit$mcmc)
  z_names <- grep("_xRE_Zx[", colnames(posterior), fixed = TRUE, value = TRUE)

  expect_s3_class(
    term$latent_layout,
    "BayesTools_random_effect_structured_local_layout"
  )
  expect_length(z_names, K)
  expect_lt(length(z_names), term$n_groups * term$n_columns)
  expect_false(any(grepl("_xRE_CORx_L[", colnames(posterior), fixed = TRUE)))
  expect_false(any(grepl("_xRE_CORx_R[", colnames(posterior), fixed = TRUE)))

  prediction <- JAGS_evaluate_formula(
    fit = fit,
    formula = ~ 1 + cs(index | id),
    parameter = "mu",
    data = df,
    prior_list = attr(fit, "prior_list")
  )
  expect_equal(dim(prediction), c(K, nrow(posterior)))
  expect_true(all(is.finite(prediction)))
})

test_that("JAGS transformed scalar rho remains representably inside support", {

  skip_if_not_installed("runjags")
  skip_if_not_installed("rjags")
  skip_on_cran()

  df <- data.frame(
    y = c(-0.2, 0.1, 0.4, -0.1, 0.2, 0.5),
    index = factor(rep(c("i1", "i2", "i3"), 2L)),
    id = factor(rep(c("g1", "g2"), each = 3L))
  )
  fit <- suppressWarnings(JAGS_fit(
    model_syntax = "model{
      for(i in 1:N_mu){
        y[i] ~ dnorm(mu[i], 4)
      }
    }",
    data = list(y = df$y),
    formula_list = list(mu = ~ 1 + ar1(index | id)),
    formula_data_list = list(mu = df),
    formula_prior_list = list(mu = list(
      intercept = prior("normal", list(0, 1))
    )),
    formula_random_prior_list = list(mu = prior_random(
      id = random_block(
        sd = prior("point", list(location = 1)),
        rho = prior("point", list(location = 1e300))
      )
    )),
    chains = 1,
    adapt = 50,
    burnin = 50,
    sample = 100,
    silent = TRUE,
    seed = 210
  ))
  random_term <- JAGS_formula_design(fit, "mu")$random_effects[[1L]]
  interior <- .bt_random_effect_representable_rho_bounds(
    random_term$correlation$bounds,
    random_term$structure
  )
  rho <- as.matrix(fit$mcmc)[, random_term$correlation$rho_name]

  expect_true(all(is.finite(rho)))
  expect_equal(unname(rho), rep(interior[["upper"]], length(rho)), tolerance = 0)
  expect_lt(max(rho), random_term$correlation$bounds[["upper"]])
})

test_that("JAGS_fit predicts row-indexed external SD random effects from latent monitors", {

  skip_if_not_installed("runjags")
  skip_if_not_installed("rjags")
  skip_on_cran()

  df <- data.frame(
    y = c(-0.2, 0.1, 0.4, 0.7),
    study = factor(c("s1", "s1", "s2", "s2")),
    drug = factor(c("a", "b", "a", "b")),
    tau_data = c(0.5, 0.75, 1.0, 1.25)
  )
  formula <- ~ 1 +
    random(1 | study, name = "study", covariance = "diag") +
    random(1 | drug, name = "drug", covariance = "diag")

  fit <- suppressWarnings(JAGS_fit(
    model_syntax = "model{
      for(i in 1:N_mu){
        y[i] ~ dnorm(mu[i], 4)
        tau[i] <- tau_data[i]
      }
    }",
    data = list(y = df$y, tau_data = df$tau_data),
    formula_list = list(mu = formula),
    formula_data_list = list(mu = df),
    formula_prior_list = list(mu = list(
      intercept = prior("normal", list(0, 1))
    )),
    formula_random_prior_list = list(mu = prior_random(
      allocation = random_variance_allocation(
        sd_source = random_sd_source("tau", shape = "row"),
        weights = prior("dirichlet", list(alpha = c(2, 3)))
      )
    )),
    add_parameters = c("tau", "mu"),
    chains = 1,
    adapt = 50,
    burnin = 50,
    sample = 100,
    silent = TRUE,
    seed = 105
  ))

  posterior <- as.matrix(fit$mcmc)
  tau_names <- paste0("tau[", seq_len(nrow(df)), "]")
  mu_names <- paste0("mu[", seq_len(nrow(df)), "]")

  expect_true(all(tau_names %in% colnames(posterior)))
  expect_true(all(mu_names %in% colnames(posterior)))
  expect_true(any(grepl("mu__xREx__study_xRE_Zx", colnames(posterior), fixed = TRUE)))
  expect_true(any(grepl("mu__xREx__drug_xRE_Zx", colnames(posterior), fixed = TRUE)))
  expect_false(any(grepl("_xRE_COEFx", colnames(posterior), fixed = TRUE)))
  expect_equal(
    unname(posterior[, tau_names, drop = FALSE]),
    unname(matrix(
      rep(df$tau_data, each = nrow(posterior)),
      nrow = nrow(posterior),
      ncol = nrow(df)
    )),
    tolerance = 1e-12
  )

  prediction <- JAGS_evaluate_formula(
    fit = fit,
    formula = formula,
    parameter = "mu",
    data = df,
    fitted_rows = seq_len(nrow(df)),
    prior_list = attr(fit, "prior_list")
  )

  expect_equal(dim(prediction), c(nrow(df), nrow(posterior)))
  expect_true(all(is.finite(prediction)))
  expect_equal(
    unname(t(prediction)),
    unname(posterior[, mu_names, drop = FALSE]),
    tolerance = 1e-8
  )
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
