context("JAGS formula")

# This file tests the JAGS formula functions
# - Helper functions for parameter naming and formula handling
# - JAGS_evaluate_formula function for prediction
# Uses pre-fitted models from test-00-model-fits.R per testing guidelines

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

test_that("JAGS evaluate formula works", {

  # Test JAGS_evaluate_formula by comparing against lm() predictions using ML estimates.
  # This test constructs samples manually (from ML estimates) - no pre-fitted JAGS model needed.
  skip_on_os(c("mac", "linux", "solaris")) # multivariate sampling does not exactly match across OSes
  skip_on_cran()
  skip_if_not_installed("rjags")

  # Setup: complex formula including scaling
  set.seed(1)
  df_all <- data.frame(
    x_cont1 = rnorm(60),
    x_cont2 = rnorm(60),
    x_fac2t = factor(rep(c("A", "B"), 30), levels = c("A", "B")),
    x_fac3o = factor(rep(c("A", "B", "C"), 20), levels = c("A", "B", "C"))
  )
  df_all$y <- rnorm(60, 0.1, 0.5) + 0.30 * df_all$x_cont1 - 0.15 * df_all$x_cont1 * df_all$x_cont2 +
    ifelse(df_all$x_fac3o == "A", 0.2, ifelse(df_all$x_fac3o == "B", -0.2, 0))

  prior_list_all <- list(
    "intercept"       = prior("normal", list(0, 5)),
    "x_cont1"         = prior("normal", list(0, 1)),
    "x_cont2"         = prior("normal", list(0, 1)),
    "x_fac2t"         = prior_factor("normal",  contrast = "treatment",   list(0, 1)),
    "x_fac3o"         = prior_factor("mnormal", contrast = "orthonormal", list(0, 1)),
    "x_cont1:x_fac3o" = prior_factor("mnormal", contrast = "orthonormal", list(0,  2))
  )
  prior_list2  <- list(
    "sigma" = prior("cauchy", list(0, 1), list(0, 1))
  )

  # Use JAGS_formula to process formula and get prior_list with correct parameter names
  formula_result <- JAGS_formula(~ x_fac2t + x_cont2 + x_cont1 * x_fac3o, parameter = "mu", data = df_all, prior_list = prior_list_all)
  prior_list     <- c(formula_result$prior_list, prior_list2)

  # Define expected column names for samples (must match what JAGS_formula produces)
  col_names <- c("mu_intercept", "mu_x_cont1", "mu_x_cont1__xXx__x_fac3o[1]", "mu_x_cont1__xXx__x_fac3o[2]",
                 "mu_x_cont2", "mu_x_fac2t", "mu_x_fac3o[1]", "mu_x_fac3o[2]", "sigma")

  new_data <-  data.frame(
    x_cont1 = c(0, 0, 1, 1),
    x_cont2 = c(0, 1, 0, 1),
    x_fac2t = factor(c("A", "B", "A", "B"), levels = c("A", "B")),
    x_fac3o = factor(c("A", "B", "C", "A"), levels = c("A", "B", "C"))
  )

  # Test the results against the lm function (by passing the ML estimates)
  contrasts(df_all$x_fac3o) <- contr.orthonormal(levels(df_all$x_fac3o))
  fit_lm <- stats::lm(y ~ x_fac2t + x_cont2 + x_cont1 * x_fac3o, data = df_all)

  # Create mock samples from ML estimates
  samples_new <- c(coef(fit_lm), sigma = sigma(fit_lm))[c("(Intercept)","x_cont1","x_cont1:x_fac3o1","x_cont1:x_fac3o2","x_cont2","x_fac2tB","x_fac3o1","x_fac3o2","sigma")]
  samples_new <- matrix(samples_new, nrow = 1)
  colnames(samples_new) <- col_names
  samples_new <- coda::as.mcmc.list(coda::as.mcmc(samples_new))

  expect_equal(predict(fit_lm, newdata = new_data), JAGS_evaluate_formula(samples_new, ~ x_fac2t + x_cont2 + x_cont1 * x_fac3o, "mu", new_data, prior_list)[,1])

  # For a posterior samples matrix (multiple rows)
  samples_new <- c(coef(fit_lm), sigma = sigma(fit_lm))[c("(Intercept)","x_cont1","x_cont1:x_fac3o1","x_cont1:x_fac3o2","x_cont2","x_fac2tB","x_fac3o1","x_fac3o2","sigma")]
  samples_new <- matrix(samples_new, nrow = 5, ncol = length(samples_new), byrow = TRUE)
  colnames(samples_new) <- col_names
  samples_new <- coda::as.mcmc.list(coda::as.mcmc(samples_new))

  expect_equal(matrix(predict(fit_lm, newdata = new_data), nrow = 4, ncol = 5), unname(JAGS_evaluate_formula(samples_new, ~ x_fac2t + x_cont2 + x_cont1 * x_fac3o, "mu", new_data, prior_list)))

  # Check filling in missing or miss-ordered factor levels
  samples_new <- c(coef(fit_lm), sigma = sigma(fit_lm))[c("(Intercept)","x_cont1","x_cont1:x_fac3o1","x_cont1:x_fac3o2","x_cont2","x_fac2tB","x_fac3o1","x_fac3o2","sigma")]
  samples_new <- matrix(samples_new, nrow = 1)
  colnames(samples_new) <- col_names
  samples_new <- coda::as.mcmc.list(coda::as.mcmc(samples_new))

  new_data2         <- new_data
  new_data2$x_fac2t <- factor(as.character(new_data2$x_fac2t), levels = c("B", "A"))
  expect_equal(predict(fit_lm, newdata = new_data), JAGS_evaluate_formula(samples_new, ~ x_fac2t + x_cont2 + x_cont1 * x_fac3o, "mu", new_data2, prior_list)[,1])

  new_data3         <- new_data
  new_data3$x_fac3o <- factor(c("A", "B", "A", "B"), levels = c("B", "A"))
  expect_equal(predict(fit_lm, newdata = new_data3), JAGS_evaluate_formula(samples_new, ~ x_fac2t + x_cont2 + x_cont1 * x_fac3o, "mu", new_data3, prior_list)[,1])

  new_data4         <- new_data
  new_data4$x_fac3o <- c("A", "B", "A", "B")
  expect_equal(predict(fit_lm, newdata = new_data3), JAGS_evaluate_formula(samples_new, ~ x_fac2t + x_cont2 + x_cont1 * x_fac3o, "mu", new_data4, prior_list)[,1])

  # Check scaling works (by multiplying by zero)
  prior_list_scaled <- prior_list
  attr(prior_list_scaled$mu_x_cont2, "multiply_by") <- 0
  attr(prior_list_scaled$mu_x_fac2t, "multiply_by") <- 0

  samples_new2 <- c(coef(fit_lm), sigma = sigma(fit_lm))[c("(Intercept)","x_cont1","x_cont1:x_fac3o1","x_cont1:x_fac3o2","x_cont2","x_fac2tB","x_fac3o1","x_fac3o2","sigma")]
  samples_new2 <- matrix(samples_new2, nrow = 1)
  colnames(samples_new2) <- col_names
  samples_new2[,"mu_x_cont2"] <- 0
  samples_new2[,"mu_x_fac2t"] <- 0
  samples_new2 <- coda::as.mcmc.list(coda::as.mcmc(samples_new2))

  expect_equal(JAGS_evaluate_formula(samples_new, ~ x_fac2t + x_cont2 + x_cont1 * x_fac3o, "mu", new_data, prior_list_scaled)[,1],
               JAGS_evaluate_formula(samples_new2, ~ x_fac2t + x_cont2 + x_cont1 * x_fac3o, "mu", new_data, prior_list)[,1])

  ### Test input validation
  expect_error(JAGS_evaluate_formula(samples_new, ~ x_fac2t + x_cont2 + x_cont1 * x_fac3o, "mu", new_data[,1:3], prior_list),
               "The 'x_fac3o' predictor variable is missing in the data.")
  expect_error(JAGS_evaluate_formula(samples_new, ~ x_fac2t + x_cont2 + x_cont1 * x_fac3o, "mu", new_data, prior_list[-1]),
               "The prior distribution for the 'x_fac2t' term is missing in the prior_list.")

  bad_data         <- new_data
  bad_data$x_fac2t <- factor(c("C", "B", "C", "B"), levels = c("B", "C"))

  expect_error(JAGS_evaluate_formula(samples_new, ~ x_fac2t + x_cont2 + x_cont1 * x_fac3o, "mu", bad_data, prior_list),
               "Levels specified in the 'x_fac2t' factor variable do not match the levels used for model specification.")

  bad_data2         <- new_data
  bad_data2$x_fac2t <- c("C", "B", "C", "B")

  expect_error(JAGS_evaluate_formula(samples_new, ~ x_fac2t + x_cont2 + x_cont1 * x_fac3o, "mu", bad_data2, prior_list),
               "Levels specified in the 'x_fac2t' factor variable do not match the levels used for model specification.")
})


test_that("JAGS evaluate formula works with spike priors", {

  # Test JAGS_evaluate_formula with spike prior distributions using pre-fitted model
  skip_on_os(c("mac", "linux", "solaris"))
  skip_on_cran()
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
  expect_equivalent(new_samples_mean[1], 1, tolerance = 0.01)
  expect_equivalent(new_samples_mean[2], 3, tolerance = 0.01)
  expect_equivalent(new_samples_mean[3], 3, tolerance = 0.01)
})


test_that("JAGS evaluate formula works with spike-and-slab and mixture priors", {

  # Test JAGS_evaluate_formula with spike-and-slab and mixture priors using pre-fitted model
  skip_on_os(c("mac", "linux", "solaris"))
  skip_on_cran()
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

  expect_equal(.remove_expressions(f1), formula(y ~ 1))
  expect_equal(.remove_expressions(f2), formula(y ~ z))
  expect_equal(.remove_expressions(f3), formula(y ~ 1))
  expect_equal(.remove_expressions(f4), formula(y ~ z))
  expect_equal(.remove_expressions(f5), formula(y ~ z))
  expect_equal(.remove_expressions(f6), formula(y ~ z))
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

  t1 <- list("1 | id")
  t2 <- list("1 + x_cont1 | id")
  t3 <- list("x_cont1 | id", "0 + x_cont2 || group")
  attr(t1[[1]], "grouping_factor") <- "id"
  attr(t2[[1]], "grouping_factor") <- "id"
  attr(t3[[1]], "grouping_factor") <- "id"
  attr(t3[[2]], "grouping_factor") <- "group"
  attr(t1[[1]], "independent") <- FALSE
  attr(t2[[1]], "independent") <- FALSE
  attr(t3[[1]], "independent") <- FALSE
  attr(t3[[2]], "independent") <- TRUE

  expect_equal(.extract_random_effects(f1), list())
  expect_equal(.extract_random_effects(f2), list())
  expect_equal(.extract_random_effects(f3), t1)
  expect_equal(.extract_random_effects(f4), t2)
  expect_equal(.extract_random_effects(f5), t2)
  expect_equal(.extract_random_effects(f6), t2)
  expect_equal(.extract_random_effects(f7), t3)

  expect_equal(.remove_random_effects(f1), formula( ~ 1))
  expect_equal(.remove_random_effects(f2), formula( ~ x_cont1))
  expect_equal(.remove_random_effects(f3), formula( ~ 1))
  expect_equal(.remove_random_effects(f4), formula( ~ 1))
  expect_equal(.remove_random_effects(f5), formula( ~ x_cont1))
  expect_equal(.remove_random_effects(f6), formula( ~ x_cont1))
  expect_equal(.remove_random_effects(f7), formula( ~ x_cont1 + x_cont2))

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
  expect_equal(.add_intercept_to_formula(~ x - 1), ~ x)
  expect_equal(.add_intercept_to_formula(~ x + y - 1), ~ x + y)
  expect_equal(.add_intercept_to_formula(~ - 1), ~ 1)

  expect_equal(.add_intercept_to_formula(~ x + 0), ~ x)
  expect_equal(.add_intercept_to_formula(~ x + y + 0), ~ x + y)
  expect_equal(.add_intercept_to_formula(~ 0), ~ 1)

})
