skip_if_not_test_profile("fixture")

source(testthat::test_path("common-functions.R"))
skip_if_missing_fits(c("fit_formula_draws_scaled_random",
  "fit_formula_draws_expression_data"))

test_that("JAGS_formula_draws() evaluates draws without a fit through the design of JAGS_fit()", {

  skip_if_not_installed("rjags")
  skip_if_not_installed("runjags")
  fit <- readRDS(file.path(temp_fits_dir, "fit_formula_draws_scaled_random.RDS"))
  inputs <- attr(fit, "formula_draws_inputs", exact = TRUE)
  data <- inputs$data
  prior_list <- inputs$prior_list
  formula <- inputs$formula
  y <- inputs$y
  random_priors <- inputs$random_priors
  posterior <- as.matrix(BayesTools:::.fit_to_posterior(fit))

  # the draws of the fit, without the fit, with the design built from the
  # formula, data, and priors of the fit
  draws <- JAGS_formula_draws(
    posterior, formula = formula, parameter = "mu", data = data,
    prior_list = prior_list, formula_scale = list(x = TRUE),
    prior_random = random_priors
  )
  expect_s3_class(draws, "mcmc")
  expect_identical(JAGS_formula_design(draws, "mu"), JAGS_formula_design(fit, "mu"))
  expect_identical(attr(draws, "formula_scale")$mu, attr(fit, "formula_scale")$mu)

  # the evaluation equals the evaluation of the fit on its data and on new
  # data (standardized and factor-coded as fitted), for every target
  newdata <- data.frame(
    x = c(-1, 4, 9),
    d = factor(c("c", "a", "b")),
    g = factor(c("g2", "g5", "g1"))
  )
  for(prediction_data in list(NULL, newdata)){
    for(formula_target in c("fixed", "conditional")){
      expect_equal(
        JAGS_evaluate_formula(draws, parameter = "mu", data = prediction_data,
                              formula_target = formula_target),
        JAGS_evaluate_formula(fit, parameter = "mu", data = prediction_data,
                              formula_target = formula_target),
        tolerance = 1e-14
      )
    }
  }
  expect_equal(
    JAGS_predict_formula(draws, "mu", formula_target = "marginal")[c("value", "vcov")],
    JAGS_predict_formula(fit, "mu", formula_target = "marginal")[c("value", "vcov")],
    tolerance = 1e-14
  )

  # draws without the design are not evaluated from the priors
  expect_error(
    JAGS_evaluate_formula(coda::mcmc(posterior), formula, "mu", data, prior_list),
    paste0(
      "JAGS_evaluate_formula() needs the fitted formula design of parameter ",
      "'mu': pass a fit from JAGS_fit() with a formula for 'mu', or posterior ",
      "draws with the design built by JAGS_formula_draws(). Refit the model ",
      "with the current BayesTools version if it was fitted by BayesTools 0.3.0."
    ),
    fixed = TRUE,
    class = "BayesTools_refit_required"
  )
  # draws that lack a coefficient of the design
  expect_error(
    JAGS_evaluate_formula(
      JAGS_formula_draws(
        posterior[, colnames(posterior) != "mu_x", drop = FALSE],
        formula = formula, parameter = "mu", data = data,
        prior_list = prior_list, formula_scale = list(x = TRUE),
        prior_random = random_priors
      ),
      parameter = "mu", formula_target = "fixed"
    ),
    paste0(
      "JAGS_evaluate_formula() needs the posterior draws of the coefficient(s) ",
      "'mu_x' of parameter 'mu', which the draws do not contain."
    ),
    fixed = TRUE
  )
  expect_error(
    JAGS_formula_draws(fit, formula, "mu", data, prior_list),
    "'draws' must be draws without a fit",
    fixed = TRUE
  )
  expect_error(
    JAGS_formula_draws(posterior[, c(1L, 1L)], formula, "mu", data, prior_list),
    "'draws' must have unique, non-empty column names",
    fixed = TRUE
  )
})

test_that("JAGS_formula_draws() rebuilds expression terms that read JAGS model data", {

  skip_if_not_installed("rjags")
  skip_if_not_installed("runjags")
  fit <- readRDS(file.path(temp_fits_dir, "fit_formula_draws_expression_data.RDS"))
  inputs <- attr(fit, "formula_draws_inputs", exact = TRUE)
  data <- inputs$data
  prior_list <- inputs$prior_list
  formula <- inputs$formula
  y <- inputs$y
  v <- inputs$v
  posterior <- as.matrix(BayesTools:::.fit_to_posterior(fit))

  # with the model data of the fit, the design is the fitted design
  draws <- JAGS_formula_draws(
    posterior, formula, "mu", data, prior_list,
    model_data = list(y = y, v = v)
  )
  expect_identical(JAGS_formula_design(draws, "mu"), JAGS_formula_design(fit, "mu"))
  expect_equal(
    JAGS_evaluate_formula(draws, parameter = "mu"),
    JAGS_evaluate_formula(fit, parameter = "mu"),
    tolerance = 1e-14
  )
  newdata <- data.frame(x = c(0, 1), v = c(1, 3))
  expect_equal(
    JAGS_evaluate_formula(draws, parameter = "mu", data = newdata),
    JAGS_evaluate_formula(fit, parameter = "mu", data = newdata),
    tolerance = 1e-14
  )

  # without them the expression cannot be replayed
  expect_error(
    JAGS_formula_draws(posterior, formula, "mu", data, prior_list),
    "expression() term 'v[i]' is not replayable: unknown replay dependency 'v'.",
    fixed = TRUE
  )
})
