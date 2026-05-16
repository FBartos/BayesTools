skip_if_not_test_profile("fixture")

# ============================================================================ #
# TEST FILE: JAGS Formula Handling
# ============================================================================ #
#
# PURPOSE:
#   Tests for JAGS formula parsing, parameter naming, and prediction functions
#   in R/JAGS-formula.R. Includes JAGS_evaluate_formula and helper utilities.
#
# DEPENDENCIES:
#   - rjags: Required for JAGS model evaluation
#   - common-functions.R: Test helpers and pre-fitted model access
#
# SKIP CONDITIONS:
#   - First section (parameter name tools): Can run on CRAN (pure R)
#   - Second section (JAGS evaluation): skip_if_not_installed("rjags")
#   - skip_on_os(): Multivariate sampling consistency (meandif priors)
#
# MODELS/FIXTURES:
#   - Uses pre-fitted models from test-00-model-fits.R via temp_fits_dir
#
# TAGS: @evaluation, @JAGS, @formula
# ============================================================================ #

# Load common test helpers
source(testthat::test_path("common-functions.R"))

test_that("JAGS formula tools work", {

  # additional tools work
  expect_equal(
    format_parameter_names(c("mu_x_cont", "mu_x_fac3t", "mu_x_fac3t__xXx__x_cont")),
    c("mu_x_cont", "mu_x_fac3t", "mu_x_fac3t:x_cont")
  )
  expect_equal(
    format_parameter_names(c("mu_x_cont", "mu_x_fac3t", "mu_x_fac3t__xXx__x_cont"), formula_parameters = "mu"),
    c("(mu) x_cont", "(mu) x_fac3t", "(mu) x_fac3t:x_cont")
  )
  expect_equal(
    format_parameter_names(c("mu_x_cont", "mu_x_fac3t", "mu_x_fac3t__xXx__x_cont"), formula_parameters = "mu", formula_prefix = FALSE),
    c("x_cont", "x_fac3t", "x_fac3t:x_cont")
  )

  expect_equal(
    JAGS_parameter_names(c("x_cont", "x_fac3t", "x_fac3t:x_cont")),
    c("x_cont", "x_fac3t", "x_fac3t__xXx__x_cont")
  )
  expect_equal(
    JAGS_parameter_names(c("x_cont", "x_fac3t", "x_fac3t:x_cont"), formula_parameter = "mu"),
    c("mu_x_cont", "mu_x_fac3t", "mu_x_fac3t__xXx__x_cont")
  )

})

test_that("formula_add_intercept repairs only top-level no-intercept terms", {

  env <- new.env(parent = emptyenv())
  formula <- stats::as.formula("~ I(x - 1) + offset(z - 1) - 1", env = env)
  attr(formula, "custom") <- "metadata"

  out <- formula_add_intercept(formula)

  expect_identical(environment(out), env)
  expect_identical(attr(out, "custom"), "metadata")
  expect_equal(attr(stats::terms(out), "intercept"), 1)
  expect_equal(attr(stats::terms(out), "term.labels"), "I(x - 1)")
  expect_equal(attr(stats::terms(out), "offset"), 2)

  expect_equal(formula_add_intercept(~ x - 1), ~ x, ignore_formula_env = TRUE)
  expect_equal(formula_add_intercept(~ 0 + x), ~ x, ignore_formula_env = TRUE)
  expect_equal(formula_add_intercept(~ x + 0), ~ x, ignore_formula_env = TRUE)
  expect_equal(formula_add_intercept(~ x + -1), ~ x, ignore_formula_env = TRUE)
  expect_equal(formula_add_intercept(y ~ 0 + x), y ~ x, ignore_formula_env = TRUE)
  expect_equal(formula_add_intercept(~ -1), ~ 1, ignore_formula_env = TRUE)

  intercept_formula <- ~ x + y
  expect_identical(formula_add_intercept(intercept_formula), intercept_formula)

  expect_error(formula_add_intercept("~ x - 1"), "'formula' must be a formula.", fixed = TRUE)
})

test_that("JAGS_formula stores exact fitted formula design metadata", {

  df <- data.frame(
    x  = c(-1, 0, 1, 2),
    x2 = c(-1, 0, 1, 2),
    f  = factor(c("a", "b", "a", "b"), levels = c("a", "b"))
  )

  result <- JAGS_formula(
    formula = ~ x + x2 + f,
    parameter = "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x         = prior("normal", list(0, 1)),
      x2        = prior("normal", list(0, 1)),
      f         = prior_factor("normal", list(0, 1), contrast = "treatment")
    ),
    formula_scale = list(x = TRUE)
  )

  design <- result$formula_design

  expect_s3_class(design, "BayesTools_formula_design")
  expect_identical(design$parameter, "mu")
  expect_equal(design$formula, ~ x + x2 + f, ignore_formula_env = TRUE)
  expect_equal(nrow(design$model_frame), nrow(df))
  expect_equal(colnames(design$model_matrix), c("(Intercept)", "x", "x2", "fb"))
  expect_equal(design$column_names, colnames(design$model_matrix))
  expect_equal(design$raw_column_names, c("(Intercept)", "x", "x2", "fb"))
  expect_equal(design$assign, attr(design$model_matrix, "assign"))
  expect_equal(attr(design$terms, "term.labels"), c("x", "x2", "f"))
  expect_equal(design$contrasts$f, "contr.treatment")
  expect_equal(design$xlevels$f, c("a", "b"))
  expect_equal(design$predictors, c("x", "x2", "f"))
  expect_equal(design$predictor_types, c(x = "continuous", x2 = "continuous", f = "factor"))
  expect_equal(design$model_terms, c("intercept", "x", "x2", "f"))
  expect_equal(design$model_terms_type, c(intercept = "continuous", x = "continuous", x2 = "continuous", f = "factor"))

  expect_equal(unname(design$model_matrix[, "x"]), unname(result$data$mu_data_x))
  expect_equal(unname(design$model_matrix[, "x2"]), unname(result$data$mu_data_x2))
  expect_equal(unname(design$model_matrix[, "fb", drop = FALSE]), unname(result$data$mu_data_f))
  expect_equal(design$jags_data_names$x, "mu_data_x")
  expect_equal(design$jags_data_names$x2, "mu_data_x2")
  expect_equal(design$jags_data_names$f, "mu_data_f")
  expect_equal(design$formula_scale$mu_x$mean, mean(df$x))
  expect_equal(design$formula_scale$mu_x$sd, stats::sd(df$x))
  expect_equal(names(design$prior_list), names(result$prior_list))

  expect_lt(design$rank, ncol(design$model_matrix))
  expect_named(design$aliased, colnames(design$model_matrix))
  expect_true(any(design$aliased))
})

test_that("JAGS_formula_design accessor returns stored designs predictably", {

  design_mu <- structure(list(parameter = "mu"), class = c("BayesTools_formula_design", "list"))
  design_tau <- structure(list(parameter = "log_tau"), class = c("BayesTools_formula_design", "list"))
  fit <- structure(list(), class = "BayesTools_fit")
  attr(fit, "formula_design") <- list(mu = design_mu, log_tau = design_tau)

  expect_identical(JAGS_formula_design(fit), list(mu = design_mu, log_tau = design_tau))
  expect_identical(JAGS_formula_design(fit, "mu"), design_mu)
  expect_identical(JAGS_formula_design(fit, "log_tau"), design_tau)
  expect_null(JAGS_formula_design(structure(list(), class = "BayesTools_fit")))
  expect_error(JAGS_formula_design(fit, "sigma"), "Formula design for parameter 'sigma' was not found.", fixed = TRUE)
})

test_that("JAGS_formula stores multi-factor contrast metadata", {

  df <- expand.grid(
    a = factor(c("a1", "a2"), levels = c("a1", "a2")),
    b = factor(c("b1", "b2", "b3"), levels = c("b1", "b2", "b3")),
    x = c(-1, 1)
  )

  formula_result <- JAGS_formula(
    formula = ~ x * a * b,
    parameter = "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1)),
      a = prior_factor("mnormal", list(0, 1), contrast = "orthonormal"),
      b = prior_factor("mnormal", list(0, 1), contrast = "meandif"),
      "x:a" = prior_factor("mnormal", list(0, 1), contrast = "orthonormal"),
      "x:b" = prior_factor("mnormal", list(0, 1), contrast = "meandif"),
      "a:b" = prior_factor("mnormal", list(0, 1), contrast = "orthonormal"),
      "x:a:b" = prior_factor("mnormal", list(0, 1), contrast = "orthonormal")
    )
  )

  interaction_prior <- formula_result$prior_list$mu_x__xXx__a__xXx__b
  level_grid <- expand.grid(
    a = c("a1", "a2"),
    b = c("b1", "b2", "b3"),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  a_contrast <- contr.orthonormal(c("a1", "a2"))
  b_contrast <- contr.meandif(c("b1", "b2", "b3"))
  a_index <- match(level_grid$a, c("a1", "a2"))
  b_index <- match(level_grid$b, c("b1", "b2", "b3"))
  expected_design <- cbind(
    a_contrast[a_index, 1] * b_contrast[b_index, 1],
    a_contrast[a_index, 1] * b_contrast[b_index, 2]
  )

  expect_equal(attr(interaction_prior, "term_components"), c("x", "a", "b"))
  expect_equal(attr(interaction_prior, "factor_terms"), c("a", "b"))
  expect_equal(
    attr(interaction_prior, "factor_contrasts"),
    c(a = "contr.orthonormal", b = "contr.meandif")
  )
  expect_equal(dim(attr(interaction_prior, "factor_design")), c(6L, 2L))
  expect_equal(
    attr(interaction_prior, "factor_cell_names"),
    c("a=a1, b=b1", "a=a2, b=b1", "a=a1, b=b2", "a=a2, b=b2", "a=a1, b=b3", "a=a2, b=b3")
  )
  expect_equal(unname(attr(interaction_prior, "factor_design")), unname(expected_design))
})

test_that(".factor_term_design_from_formula handles no-intercept factor interactions", {

  df <- expand.grid(
    a = factor(c("a1", "a2"), levels = c("a1", "a2")),
    b = factor(c("b1", "b2", "b3"), levels = c("b1", "b2", "b3"))
  )
  stats::contrasts(df$a) <- "contr.orthonormal"
  stats::contrasts(df$b) <- "contr.meandif"

  formula <- ~ a + b + a:b - 1
  predictors <- as.character(attr(stats::terms(formula), "variables"))[-1]
  model_matrix <- stats::model.matrix(stats::model.frame(formula, data = df), formula = formula, data = df)
  term_index <- which(attr(stats::terms(formula), "term.labels") == "a:b")

  design_info <- BayesTools:::.factor_term_design_from_formula(
    formula = formula,
    data = df,
    predictors = predictors,
    predictors_type = c(a = "factor", b = "factor"),
    term_index = term_index,
    term_components = c("a", "b"),
    factor_terms = c("a", "b"),
    has_intercept = FALSE
  )

  expect_equal(
    unname(design_info$design),
    unname(model_matrix[, attr(model_matrix, "assign") == term_index, drop = FALSE])
  )
  expect_equal(
    design_info$cell_names,
    c("a=a1, b=b1", "a=a2, b=b1", "a=a1, b=b2", "a=a2, b=b2", "a=a1, b=b3", "a=a2, b=b3")
  )
})

# ============================================================================ #
# SECTION: Tests requiring JAGS (skip conditions per test)
# ============================================================================ #

test_that("JAGS evaluate formula works with spike priors", {

  # Test JAGS_evaluate_formula with spike prior distributions using pre-fitted model
  skip_on_os(c("mac", "linux", "solaris"))
  skip_on_cran()
  skip_if_no_fits()
  skip_if_not_installed("rjags")

  # Load pre-fitted model with spike factor priors (all 4 contrast types)
  fit_spike <- readRDS(file.path(temp_fits_dir, "fit_spike_factors.RDS"))

  # New data for prediction
  new_data <- data.frame(
    x_fac2i  = factor(c("A", "B", "A"), levels = c("A", "B")),
    x_fac3o  = factor(c("A", "A", "B"), levels = c("A", "B", "C")),
    x_fac3t  = factor(c("A", "B", "C"), levels = c("A", "B", "C")),
    x_fac3md = factor(c("B", "B", "C"), levels = c("A", "B", "C"))
  )

  # Note: fit_spike_factors uses formula ~ x_fac2i + x_fac3o + x_fac3t + x_fac3md - 1
  # with spike priors: independent(1), orthonormal(0), treatment(2), meandif(0)
  prior_list <- attr(fit_spike, "prior_list")
  new_samples <- JAGS_evaluate_formula(fit_spike, ~ x_fac2i + x_fac3o + x_fac3t + x_fac3md - 1, "mu", new_data, prior_list)
  new_samples_mean <- apply(new_samples, 1, mean)

  # Verify spike values are correctly applied:
  # - x_fac2i independent spike(1): each level gets value 1
  # - x_fac3o orthonormal spike(0): contrast coefficients are 0
  # - x_fac3t treatment spike(2): non-reference levels get value 2
  # - x_fac3md meandif spike(0): differences from mean are 0
  # Row 1: A(1) + A(0) + A(ref=0) + B(0) = 1
  # Row 2: B(1) + A(0) + B(2) + B(0) = 3
  # Row 3: A(1) + B(0) + C(2) + C(0) = 3
  expect_equal(new_samples_mean[1], 1, tolerance = 0.01, ignore_attr = TRUE)
  expect_equal(new_samples_mean[2], 3, tolerance = 0.01, ignore_attr = TRUE)
  expect_equal(new_samples_mean[3], 3, tolerance = 0.01, ignore_attr = TRUE)
})

test_that("JAGS evaluate formula works with spike-and-slab and mixture priors", {

  # Test JAGS_evaluate_formula with spike-and-slab and mixture priors using pre-fitted model
  skip_on_os(c("mac", "linux", "solaris"))
  skip_on_cran()
  skip_if_no_fits()
  skip_if_not_installed("rjags")

  # Load pre-fitted joint complex model (mixture intercept, spike-and-slab continuous, spike-and-slab factor)
  fit_joint <- readRDS(file.path(temp_fits_dir, "fit_joint_complex.RDS"))

  # New data for prediction
  new_data <- data.frame(
    x_cont1 = c(0, 1, -1),
    x_fac3t = factor(c("A", "B", "C"), levels = c("A", "B", "C"))
  )

  # fit_joint_complex uses formula ~ x_cont1 + x_fac3t
  prior_list <- attr(fit_joint, "prior_list")
  new_samples <- JAGS_evaluate_formula(fit_joint, ~ x_cont1 + x_fac3t, "mu", new_data, prior_list)

  # Should return samples for 3 new data points x number of posterior samples
  expect_equal(nrow(new_samples), 3)
  expect_equal(ncol(new_samples), 1000)
})

test_that("Expression handling functions work", {

  f1 <- formula(y ~ 1)
  f2 <- formula(y ~ z)
  f3 <- formula(y ~ expression(x))
  f4 <- formula(y ~ z + expression(x))
  f5 <- formula(y ~ expression(x) + z)
  f6 <- formula(y ~ expression(x) + z + expression(b))

  expect_true(!.has_expression(f1))
  expect_true(!.has_expression(f2))
  expect_true(.has_expression(f3))
  expect_true(.has_expression(f4))
  expect_true(.has_expression(f5))
  expect_true(.has_expression(f6))

  expect_equal(.extract_expressions(f3), list("x"))
  expect_equal(.extract_expressions(f4), list("x"))
  expect_equal(.extract_expressions(f5), list("x"))
  expect_equal(.extract_expressions(f6), list("x", "b"))

  expect_equal(.remove_expressions(f1), formula(y ~ 1), ignore_formula_env = TRUE)
  expect_equal(.remove_expressions(f2), formula(y ~ z), ignore_formula_env = TRUE)
  expect_equal(.remove_expressions(f3), formula(y ~ 1), ignore_formula_env = TRUE)
  expect_equal(.remove_expressions(f4), formula(y ~ z), ignore_formula_env = TRUE)
  expect_equal(.remove_expressions(f5), formula(y ~ z), ignore_formula_env = TRUE)
  expect_equal(.remove_expressions(f6), formula(y ~ z), ignore_formula_env = TRUE)
})

test_that("Random effects handling functions work", {

  f1 <- formula( ~ 1)
  f2 <- formula( ~ x_cont1)
  f3 <- formula( ~ (1 | id))
  f4 <- formula( ~ (1 + x_cont1 | id))
  f5 <- formula( ~ (1 + x_cont1 | id) + x_cont1)
  f6 <- formula( ~ x_cont1 + (1 + x_cont1 | id))
  f7 <- formula( ~ x_cont1 + (  x_cont1 | id   ) +   x_cont2 + (0   + x_cont2 ||  group))

  expect_true(!.has_random_effects(f1))
  expect_true(!.has_random_effects(f2))
  expect_true(.has_random_effects(f3))
  expect_true(.has_random_effects(f4))
  expect_true(.has_random_effects(f5))
  expect_true(.has_random_effects(f6))
  expect_true(.has_random_effects(f7))

  parsed_f3 <- .bt_parse_random_effects(f3)$terms
  parsed_f4 <- .bt_parse_random_effects(f4)$terms
  parsed_f7 <- .bt_parse_random_effects(f7)$terms

  expect_length(parsed_f3, 1L)
  expect_equal(parsed_f3[[1]]$term_formula, ~ 1, ignore_formula_env = TRUE)
  expect_equal(parsed_f3[[1]]$group_label, "id")
  expect_false(isTRUE(parsed_f3[[1]]$independent))

  expect_length(parsed_f4, 1L)
  expect_equal(parsed_f4[[1]]$term_formula, ~ 1 + x_cont1, ignore_formula_env = TRUE)
  expect_equal(parsed_f4[[1]]$group_label, "id")
  expect_false(isTRUE(parsed_f4[[1]]$independent))

  expect_length(parsed_f7, 2L)
  expect_equal(parsed_f7[[1]]$term_formula, ~ x_cont1, ignore_formula_env = TRUE)
  expect_equal(parsed_f7[[1]]$group_label, "id")
  expect_false(isTRUE(parsed_f7[[1]]$independent))
  expect_equal(parsed_f7[[2]]$term_formula, ~ 0 + x_cont2, ignore_formula_env = TRUE)
  expect_equal(parsed_f7[[2]]$group_label, "group")
  expect_true(isTRUE(parsed_f7[[2]]$independent))

  expect_equal(.remove_random_effects(f1), formula( ~ 1), ignore_formula_env = TRUE)
  expect_equal(.remove_random_effects(f2), formula( ~ x_cont1), ignore_formula_env = TRUE)
  expect_equal(.remove_random_effects(f3), formula( ~ 1), ignore_formula_env = TRUE)
  expect_equal(.remove_random_effects(f4), formula( ~ 1), ignore_formula_env = TRUE)
  expect_equal(.remove_random_effects(f5), formula( ~ x_cont1), ignore_formula_env = TRUE)
  expect_equal(.remove_random_effects(f6), formula( ~ x_cont1), ignore_formula_env = TRUE)
  expect_equal(.remove_random_effects(f7), formula( ~ x_cont1 + x_cont2), ignore_formula_env = TRUE)

})

test_that("-1 (no intercept) formula handling works correctly", {

  # setup test data
  set.seed(1)
  df_test <- data.frame(
    x_fac3md = factor(rep(c("A", "B", "C"), 20), levels = c("A", "B", "C")),
    x_fac3i  = factor(rep(c("A", "B", "C"), 20), levels = c("A", "B", "C")),
    x_cont   = rnorm(60)
  )

  # Test 1: Basic -1 formula functionality
  prior_list_basic <- list(
    "x_fac3md" = prior_factor("mnormal", contrast = "meandif", list(0, 1))
  )
  result_basic <- JAGS_formula(~ x_fac3md - 1, parameter = "mu",
                              data = df_test[, "x_fac3md", drop = FALSE],
                              prior_list = prior_list_basic)

  # The -1 should automatically add spike(0) intercept
  expect_true("mu_intercept" %in% names(result_basic$prior_list))
  expect_true(is.prior.point(result_basic$prior_list$mu_intercept))
  expect_equal(result_basic$prior_list$mu_intercept$parameters$location, 0)
  expect_true(grepl("mu_intercept", result_basic$formula_syntax))

  # Test 2: Helper function test
  expect_equal(.add_intercept_to_formula(~ x - 1), ~ x, ignore_formula_env = TRUE)
  expect_equal(.add_intercept_to_formula(~ x + y - 1), ~ x + y, ignore_formula_env = TRUE)
  expect_equal(.add_intercept_to_formula(~ - 1), ~ 1, ignore_formula_env = TRUE)

  expect_equal(.add_intercept_to_formula(~ x + 0), ~ x, ignore_formula_env = TRUE)
  expect_equal(.add_intercept_to_formula(~ x + y + 0), ~ x + y, ignore_formula_env = TRUE)
  expect_equal(.add_intercept_to_formula(~ 0), ~ 1, ignore_formula_env = TRUE)

  skip_if_not_installed("coda")

  prior_list_continuous <- list(
    "x_cont" = prior("normal", list(0, 1))
  )
  result_continuous <- JAGS_formula(~ x_cont - 1, parameter = "mu",
                                    data = df_test[, "x_cont", drop = FALSE],
                                    prior_list = prior_list_continuous)
  posterior <- matrix(c(1, 2), nrow = 2)
  colnames(posterior) <- "mu_x_cont"
  posterior <- coda::as.mcmc(posterior)

  expect_false("mu_intercept" %in% colnames(posterior))
  expect_equal(
    JAGS_evaluate_formula(
      posterior,
      ~ x_cont - 1,
      "mu",
      data.frame(x_cont = c(3, 4)),
      result_continuous$prior_list
    ),
    cbind(c(3, 4), c(6, 8)),
    ignore_attr = TRUE
  )

})

test_that("log(intercept) attribute works for specifying log(int) + sum(beta_i * x_i) models", {

  # this is helpful for specifying models for e.g., standard deviation where the output must be positive,
  # but we want the intercept to be specified on the original scale - we can take exp() of the whole formula output

  # setup test data
  set.seed(1)
  df_test <- data.frame(
    x_fac3md = factor(rep(c("A", "B", "C"), 20), levels = c("A", "B", "C")),
    x_fac3i  = factor(rep(c("A", "B", "C"), 20), levels = c("A", "B", "C")),
    x_cont   = rnorm(60)
  )

  # Test 1: Basic -1 formula functionality
  prior_list_basic <- list(
    "intercept" = prior("normal", list(0, 1)),
    "x_fac3md"  = prior_factor("mnormal", contrast = "meandif", list(0, 1))
  )

  # no log intercept
  result_basic <- JAGS_formula(~ 1 + x_fac3md, parameter = "mu",
                               data = df_test[, "x_fac3md", drop = FALSE],
                               prior_list = prior_list_basic)

  # log intercept
  formula <- ~ 1 + x_fac3md
  attr(formula, "log(intercept)") <- TRUE
  result_log   <- JAGS_formula(formula, parameter = "mu",
                               data = df_test[, "x_fac3md", drop = FALSE],
                               prior_list = prior_list_basic)

  # generates normal intercept
  expect_equal(
    result_basic[["formula_syntax"]],
    "for(i in 1:N_mu){\n  mu[i] = mu_intercept + inprod(mu_x_fac3md, mu_data_x_fac3md[i,])\n}\n"
  )

  # generates log intercept
  expect_equal(
    result_log[["formula_syntax"]],
    "for(i in 1:N_mu){\n  mu[i] = log(mu_intercept) + inprod(mu_x_fac3md, mu_data_x_fac3md[i,])\n}\n"
  )

  # everything else should match
  result_basic[["formula_syntax"]] <- NULL
  result_log[["formula_syntax"]]   <- NULL
  result_basic[["formula"]] <- NULL
  result_log[["formula"]]   <- NULL
  expect_equal(result_basic, result_log)
})

test_that("JAGS_evaluate_formula works with log(intercept) attribute", {

  # Test that JAGS_evaluate_formula correctly applies log() transformation to intercept
  # when the formula has the log(intercept) attribute set

  skip_if_not_installed("coda")

  # Setup: simple data for testing
  set.seed(1)
  df_test <- data.frame(
    x_cont = rnorm(10)
  )

  # Create prior list with gamma prior for intercept (must be positive for log)
  prior_list <- list(
    "intercept" = prior("gamma", list(2, 1)),
    "x_cont"    = prior("normal", list(0, 1))
  )

  # Process formula to get prior_list with parameter names
  formula_result <- JAGS_formula(~ x_cont, parameter = "mu", data = df_test, prior_list = prior_list)
  prior_list_processed <- formula_result$prior_list


  # Create mock samples: intercept = 2, x_cont = 0.5
  samples <- matrix(c(2, 0.5), nrow = 1)
  colnames(samples) <- c("mu_intercept", "mu_x_cont")
  samples <- coda::as.mcmc.list(coda::as.mcmc(samples))

  # New data for prediction
  new_data <- data.frame(x_cont = c(0, 1, -1))

  # Test without log(intercept): result = intercept + x_cont * data
  # For x_cont =  0: result = 2 + 0.5 * 0 = 2
  # For x_cont =  1: result = 2 + 0.5 * 1 = 2.5
  # For x_cont = -1: result = 2 + 0.5 * (-1) = 1.5
  formula_no_log <- ~ x_cont
  result_no_log <- JAGS_evaluate_formula(samples, formula_no_log, "mu", new_data, prior_list_processed)
  expect_equal(as.vector(result_no_log[,1]), c(2, 2.5, 1.5), tolerance = 1e-10)

  # Test with log(intercept): result = log(intercept) + x_cont * data
  # For x_cont =  0: result = log(2) + 0.5 * 0 = log(2)
  # For x_cont =  1: result = log(2) + 0.5 * 1 = log(2) + 0.5
  # For x_cont = -1: result = log(2) + 0.5 * (-1) = log(2) - 0.5
  formula_log <- ~ x_cont
  attr(formula_log, "log(intercept)") <- TRUE
  result_log <- JAGS_evaluate_formula(samples, formula_log, "mu", new_data, prior_list_processed)
  expect_equal(as.vector(result_log[,1]), c(log(2), log(2) + 0.5, log(2) - 0.5), tolerance = 1e-10)
})

test_that("Default priors (__default_factor and __default_continuous) work correctly", {

  # setup test data
  set.seed(1)
  df_test <- data.frame(
    x_cont1 = rnorm(60),
    x_cont2 = rnorm(60),
    x_fac3  = factor(rep(c("A", "B", "C"), 20), levels = c("A", "B", "C")),
    x_fac2  = factor(rep(c("X", "Y"), 30), levels = c("X", "Y"))
  )

  # Test 1: Only __default_continuous - applies to intercept and continuous predictors
  prior_list_cont_default <- list(
    "__default_continuous" = prior("normal", list(0, 1))
  )
  result1 <- JAGS_formula(~ x_cont1 + x_cont2, parameter = "mu",
                          data = df_test, prior_list = prior_list_cont_default)

  # Check that intercept and both continuous predictors got the default prior
  expect_true("mu_intercept" %in% names(result1$prior_list))
  expect_true("mu_x_cont1" %in% names(result1$prior_list))
  expect_true("mu_x_cont2" %in% names(result1$prior_list))
  expect_equal(result1$prior_list$mu_intercept$distribution, "normal")
  expect_equal(result1$prior_list$mu_x_cont1$distribution, "normal")
  expect_equal(result1$prior_list$mu_x_cont2$distribution, "normal")

  # Test 2: Only __default_factor - continuous predictors must still be specified
  prior_list_fac_default <- list(
    "intercept"            = prior("normal", list(0, 5)),
    "x_cont1"              = prior("cauchy", list(0, 1)),
    "__default_factor"     = prior_factor("normal", list(0, 0.5), contrast = "treatment")
  )
  result2 <- JAGS_formula(~ x_cont1 + x_fac3 + x_fac2, parameter = "mu",
                          data = df_test, prior_list = prior_list_fac_default)

  # Check that factors got the default prior
  expect_true("mu_x_fac3" %in% names(result2$prior_list))
  expect_true("mu_x_fac2" %in% names(result2$prior_list))
  expect_equal(result2$prior_list$mu_x_fac3$distribution, "normal")
  expect_equal(result2$prior_list$mu_x_fac2$distribution, "normal")
  # Check that explicit priors are preserved
  expect_equal(result2$prior_list$mu_intercept$distribution, "normal")
  expect_equal(result2$prior_list$mu_intercept$parameters$mean, 0)
  expect_equal(result2$prior_list$mu_intercept$parameters$sd, 5)
  expect_equal(result2$prior_list$mu_x_cont1$distribution, "t")  # cauchy is internally stored as t

  # Test 3: Both defaults - all terms get assigned correctly
  prior_list_both_defaults <- list(
    "__default_continuous" = prior("normal", list(0, 2)),
    "__default_factor"     = prior_factor("normal", list(0, 1), contrast = "treatment")
  )
  result3 <- JAGS_formula(~ x_cont1 + x_fac3, parameter = "mu",
                          data = df_test, prior_list = prior_list_both_defaults)

  expect_equal(result3$prior_list$mu_intercept$distribution, "normal")
  expect_equal(result3$prior_list$mu_intercept$parameters$sd, 2)  # from continuous default
  expect_equal(result3$prior_list$mu_x_cont1$parameters$sd, 2)    # from continuous default
  expect_equal(result3$prior_list$mu_x_fac3$parameters$sd, 1)     # from factor default

  # Test 4: Explicit priors override defaults
  prior_list_override <- list(
    "intercept"            = prior("cauchy", list(0, 10)),  # explicit override
    "__default_continuous" = prior("normal", list(0, 1)),
    "__default_factor"     = prior_factor("normal", list(0, 0.5), contrast = "treatment"),
    "x_fac3"               = prior_factor("mnormal", list(0, 2), contrast = "orthonormal")  # explicit override
  )
  result4 <- JAGS_formula(~ x_cont1 + x_fac3 + x_fac2, parameter = "mu",
                          data = df_test, prior_list = prior_list_override)

  # Explicit priors should be used
  expect_equal(result4$prior_list$mu_intercept$distribution, "t")  # cauchy is internally stored as t
  expect_equal(result4$prior_list$mu_x_fac3$distribution, "mnormal")
  expect_equal(result4$prior_list$mu_x_fac3$parameters$sd, 2)
  # Default priors for non-specified terms
  expect_equal(result4$prior_list$mu_x_cont1$distribution, "normal")
  expect_equal(result4$prior_list$mu_x_fac2$distribution, "normal")
  expect_equal(result4$prior_list$mu_x_fac2$parameters$sd, 0.5)

  # Test 5: Interactions use factor default when they involve factors
  prior_list_interaction <- list(
    "__default_continuous" = prior("normal", list(0, 1)),
    "__default_factor"     = prior_factor("mnormal", list(0, 0.5), contrast = "orthonormal")
  )
  result5 <- JAGS_formula(~ x_cont1 * x_fac3, parameter = "mu",
                          data = df_test, prior_list = prior_list_interaction)

  # x_cont1:x_fac3 interaction involves a factor, so should get factor default
  expect_true("mu_x_cont1__xXx__x_fac3" %in% names(result5$prior_list))
  expect_equal(result5$prior_list[["mu_x_cont1__xXx__x_fac3"]]$distribution, "mnormal")

  # Test 6: Error when term is missing and no appropriate default
  prior_list_missing <- list(
    "__default_continuous" = prior("normal", list(0, 1))
    # no __default_factor, and x_fac3 not specified
 )
  expect_error(
    JAGS_formula(~ x_cont1 + x_fac3, parameter = "mu",
                 data = df_test, prior_list = prior_list_missing),
    "missing"
  )

  # Test 7: Reserved names cannot be used as variable names in data
  df_bad <- data.frame(
    `__default_factor` = rnorm(10),
    x = rnorm(10),
    check.names = FALSE
  )
  expect_error(
    JAGS_formula(~ x, parameter = "mu", data = df_bad,
                 prior_list = list("intercept" = prior("normal", list(0, 1)),
                                   "x" = prior("normal", list(0, 1)))),
    "__default_factor"
  )

})
