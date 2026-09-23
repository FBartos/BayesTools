skip_if_not_test_profile("unit")

# ============================================================================ #
# TEST FILE: Utility Functions Evaluation Tests
# ============================================================================ #
#
# PURPOSE:
#   Tests for utility functions behavior (not input validation).
#   Includes .is.wholenumber, transformation checks, stan extraction, etc.
#
# DEPENDENCIES:
#   - No external packages required beyond testthat
#
# SKIP CONDITIONS:
#   - None (fast, pure R tests)
#
# MODELS/FIXTURES:
#   - None required
#
# TAGS: @evaluation, @fast
# ============================================================================ #


test_that(".is.wholenumber works correctly", {

  # Positive cases

  expect_true(BayesTools:::.is.wholenumber(0))
  expect_true(BayesTools:::.is.wholenumber(5))
  expect_true(BayesTools:::.is.wholenumber(-3))
  expect_true(BayesTools:::.is.wholenumber(1e10))

  # Negative cases
  expect_false(BayesTools:::.is.wholenumber(0.5))
  expect_false(BayesTools:::.is.wholenumber(1.1))
  expect_false(BayesTools:::.is.wholenumber(-3.5))

  # NA handling
  expect_true(is.na(BayesTools:::.is.wholenumber(NA)))
  expect_equal(BayesTools:::.is.wholenumber(NA, na.rm = TRUE), logical(0))

  # Vector input
  expect_equal(BayesTools:::.is.wholenumber(c(1, 2, 3.5)), c(TRUE, TRUE, FALSE))
  expect_equal(BayesTools:::.is.wholenumber(c(1, NA, 3.5)), c(TRUE, NA, FALSE))
})


test_that("transformation input validation works", {

  # Valid transformation
  expect_null(.check_transformation_input(transformation = list(
    "fun" = function(x) exp(x),
    "inv" = function(x) log(x),
    "jac" = function(x) exp(x)
  ), NULL, FALSE))

  # Missing 'jac' component
  expect_error(.check_transformation_input(transformation = list(
    "fun" = function(x) exp(x),
    "inv" = function(x) log(x),
    "err" = function(x) exp(x)
  ), NULL, FALSE), "The 'jac' objects are missing in the 'transformation' argument.")

  # Invalid format
  expect_error(.check_transformation_input(transformation = 1, NULL, FALSE),
               "Uknown format of the 'transformation' argument.")
})


test_that("stan extraction requires rstan fit", {
  expect_error(.extract_stan(NULL), "'fit' must be an rstan fit")
})

.stan_csv_fit_for_test <- function(chains){

  paths <- file.path(tempfile("stan_chain_"), paste0("chain_", seq_along(chains), ".csv"))
  dir.create(dirname(paths[[1L]]))
  for(i in seq_along(chains)){
    header <- c(
      "# model = mock_model",
      "# method = sample (Default)",
      "#   sample",
      paste0("#     num_samples = ", nrow(chains[[i]])),
      "#     num_warmup = 0",
      "#     save_warmup = 0 (Default)",
      "#     thin = 1 (Default)",
      paste0("#   id = ", i),
      "#   random",
      "#     seed = 1"
    )
    footer <- c(
      "#  Elapsed Time: 0.01 seconds (Warm-up)",
      "#                0.02 seconds (Sampling)",
      "#                0.03 seconds (Total)"
    )
    body <- utils::capture.output(utils::write.csv(
      chains[[i]], row.names = FALSE, quote = FALSE
    ))
    writeLines(c(header, body, footer), paths[[i]])
  }
  rstan::read_stan_csv(paths)
}

test_that("stan extraction keeps matrix-valued parameters and draw counts", {

  skip_if_not_installed("rstan")
  set.seed(91)
  make_chain <- function(n){
    omega <- round(stats::runif(n, -0.9, 0.9), 6)
    cbind(
      lp__ = round(stats::rnorm(n), 6),
      accept_stat__ = round(stats::runif(n), 6),
      mu = round(stats::rnorm(n), 6),
      Omega.1.1 = 1,
      Omega.2.1 = omega,
      Omega.1.2 = omega,
      Omega.2.2 = 1,
      beta.1 = round(stats::rnorm(n), 6),
      sigma = round(stats::rexp(n), 6)
    )
  }
  chains <- list(make_chain(50L), make_chain(50L))
  fit <- .stan_csv_fit_for_test(chains)
  stacked <- do.call(rbind, chains)
  element_names <- c(
    "mu", "Omega[1,1]", "Omega[2,1]", "Omega[1,2]", "Omega[2,2]",
    "beta[1]", "sigma", "lp__"
  )
  csv_names <- c(
    "mu", "Omega.1.1", "Omega.2.1", "Omega.1.2", "Omega.2.2",
    "beta.1", "sigma", "lp__"
  )

  # Every retained draw once, in chain order, one column per element.
  flat <- .extract_stan(fit, drop = FALSE)
  expect_identical(dim(flat), c(100L, 8L))
  expect_identical(colnames(flat), element_names)
  expect_equal(unname(flat), unname(stacked[, csv_names]), tolerance = 0)

  # Single-element arrays drop their index only when requested.
  dropped <- .extract_stan(fit)
  expect_identical(
    colnames(dropped),
    replace(element_names, element_names == "beta[1]", "beta")
  )
  expect_equal(unname(dropped), unname(flat), tolerance = 0)

  # Transformations reach matrix-valued elements of the summary table.
  attr(fit, "prior_list") <- list(mu = prior("normal", list(0, 1)))
  transformed <- stan_estimates_table(
    fit,
    transformations = list("Omega[2,1]" = list(fun = function(x) 2 * x + 1))
  )
  expected <- 2 * stacked[, "Omega.2.1"] + 1
  expect_equal(transformed["Omega[2,1]", "Mean"], mean(expected), tolerance = 1e-12)
  expect_equal(transformed["Omega[2,1]", "SD"], stats::sd(expected), tolerance = 1e-12)
  expect_equal(
    transformed["Omega[2,1]", "0.5"],
    stats::median(expected),
    tolerance = 1e-12
  )
})


test_that("depreciation warnings work", {
  expect_warning(.depreciate.transform_orthonormal(TRUE, FALSE),
                 "'transform_orthonormal' argument will be depreciated in favor of 'transform_factors' argument.")
})
