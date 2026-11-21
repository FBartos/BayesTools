skip_if_not_test_profile("unit")

# ============================================================================ #
# TEST FILE: JAGS Formula - Fixed-effect formula evaluation
# ============================================================================ #
#
# PURPOSE:
#   JAGS_evaluate_formula() on fixed-effect formulas: interaction-only predictors and
#   agreement with lm() predictions.
# ============================================================================ #

.jags_formula_posterior_from_lm <- function(formula_result, lm_fit) {
  design <- formula_result$formula_design
  coefficients <- stats::coef(lm_fit)
  posterior <- numeric()

  intercept_column <- which(design$assign == 0L)
  if (length(intercept_column) == 1L) {
    posterior[JAGS_parameter_names("intercept", formula_parameter = design$parameter)] <-
      coefficients[[design$raw_column_names[[intercept_column]]]]
  }

  term_labels <- attr(design$terms, "term.labels")
  for (term_index in seq_along(term_labels)) {
    raw_term <- term_labels[[term_index]]
    jags_term <- gsub(":", "__xXx__", raw_term, fixed = TRUE)
    columns <- which(design$assign == term_index)
    values <- coefficients[design$raw_column_names[columns]]
    base_name <- JAGS_parameter_names(jags_term, formula_parameter = design$parameter)

    if (design$model_terms_type[[jags_term]] == "factor" && length(values) > 1L) {
      names(values) <- paste0(base_name, "[", seq_along(values), "]")
    } else {
      names(values) <- base_name
    }
    posterior <- c(posterior, values)
  }

  posterior <- coda::mcmc(matrix(posterior, nrow = 1, dimnames = list(NULL, names(posterior))))
  attr(posterior, "formula_design") <- stats::setNames(list(design), design$parameter)
  posterior
}

test_that("JAGS_evaluate_formula resolves interaction-only continuous predictors", {
  fitted_data <- data.frame(
    x = c(-2, -1, 1, 2),
    z = c(1, 3, 2, 4)
  )
  formula_result <- JAGS_formula(
    ~ x:z,
    "mu",
    fitted_data,
    list(
      intercept = prior("normal", list(0, 1)),
      "x:z" = prior("normal", list(0, 1))
    )
  )
  posterior <- coda::mcmc(matrix(
    c(1.5, -0.75),
    nrow = 1,
    dimnames = list(NULL, c("mu_intercept", "mu_x__xXx__z"))
  ))
  attr(posterior, "formula_design") <- list(mu = formula_result$formula_design)

  newdata <- data.frame(x = c(-3, 0.5, 4), z = c(2, -1, 3))
  expected <- drop(
    stats::model.matrix(~ x:z, data = newdata) %*% c(1.5, -0.75)
  )
  prediction <- JAGS_evaluate_formula(
    posterior,
    ~ x:z,
    "mu",
    newdata,
    formula_result$prior_list
  )

  expect_equal(unname(drop(prediction)), unname(expected), tolerance = 1e-12)

  # draws without the fitted design are not evaluated from the priors
  attr(posterior, "formula_design") <- NULL
  expect_error(
    JAGS_evaluate_formula(
      posterior,
      ~ x:z,
      "mu",
      newdata,
      formula_result$prior_list
    ),
    "JAGS_evaluate_formula() needs the fitted formula design of parameter 'mu': pass a fit from JAGS_fit() with a formula for 'mu', or posterior draws with the design built by JAGS_formula_draws(). Refit the model with the current BayesTools version if it was fitted by BayesTools 0.3.0.",
    fixed = TRUE,
    class = "BayesTools_refit_required"
  )
})

test_that("JAGS_evaluate_formula replays interaction-only factor metadata", {
  fitted_data <- data.frame(
    x = c(-2, -1, 1, 2, 3, 4),
    g = factor(c("a", "b", "c", "a", "b", "c"), levels = c("c", "a", "b"))
  )
  formula_result <- JAGS_formula(
    ~ x:g,
    "mu",
    fitted_data,
    list(
      intercept = prior("normal", list(0, 1)),
      "x:g" = prior_factor("normal", list(0, 1), contrast = "treatment")
    )
  )
  posterior <- coda::mcmc(matrix(
    c(1.5, -0.75, 0.25, 1.25),
    nrow = 1,
    dimnames = list(
      NULL,
      c(
        "mu_intercept",
        "mu_x__xXx__g[1]",
        "mu_x__xXx__g[2]",
        "mu_x__xXx__g[3]"
      )
    )
  ))
  attr(posterior, "formula_design") <- list(mu = formula_result$formula_design)

  newdata <- data.frame(
    x = c(-3, 0.5, 4),
    g = factor(c("b", "c", "a"), levels = c("unused", "a", "b", "c"))
  )
  canonical_newdata <- newdata
  canonical_newdata$g <- factor(
    canonical_newdata$g,
    levels = formula_result$formula_design$xlevels$g
  )
  stats::contrasts(canonical_newdata$g) <- "contr.treatment"
  expected <- drop(
    stats::model.matrix(~ x:g, data = canonical_newdata) %*%
      c(1.5, -0.75, 0.25, 1.25)
  )
  prediction <- JAGS_evaluate_formula(
    posterior,
    ~ x:g,
    "mu",
    newdata,
    formula_result$prior_list
  )

  expect_equal(unname(drop(prediction)), unname(expected), tolerance = 1e-12)

  # draws without the fitted design are not evaluated from the priors
  attr(posterior, "formula_design") <- NULL
  expect_error(
    JAGS_evaluate_formula(
      posterior,
      ~ x:g,
      "mu",
      newdata,
      formula_result$prior_list
    ),
    "JAGS_evaluate_formula() needs the fitted formula design of parameter 'mu': pass a fit from JAGS_fit() with a formula for 'mu', or posterior draws with the design built by JAGS_formula_draws(). Refit the model with the current BayesTools version if it was fitted by BayesTools 0.3.0.",
    fixed = TRUE,
    class = "BayesTools_refit_required"
  )
})

test_that("JAGS_evaluate_formula matches lm predictions with automatic scaling", {

  data <- bayestools_oracle_gaussian_regression_data()
  formula_data <- data[c("x1", "x2")]
  prior_list <- list(
    intercept = prior("normal", list(0, 10)),
    x1 = prior("normal", list(0, 5)),
    x2 = prior("normal", list(0, 5)),
    "x1:x2" = prior("normal", list(0, 5))
  )

  formula_result <- JAGS_formula(
    formula = ~ x1 * x2,
    parameter = "mu",
    data = formula_data,
    prior_list = prior_list,
    formula_scale = list(x1 = TRUE, x2 = TRUE)
  )

  scaled_formula_data <- bayestools_manual_scaled_data(formula_data, c("x1", "x2"))
  expect_formula_scale_equal(formula_result$formula_scale, attr(scaled_formula_data, "manual_scale"))

  lm_fit <- stats::lm(y ~ x1 * x2, data = cbind(y = data$y, scaled_formula_data))
  scaled_coefficients <- stats::coef(lm_fit)
  names(scaled_coefficients) <- bayestools_lm_coef_to_jags_names(names(scaled_coefficients))

  fit <- JAGS_formula_draws(
    matrix(
      scaled_coefficients,
      nrow = 1,
      dimnames = list(NULL, names(scaled_coefficients))
    ),
    formula = ~ x1 * x2,
    parameter = "mu",
    data = formula_data,
    prior_list = prior_list,
    formula_scale = list(x1 = TRUE, x2 = TRUE)
  )

  fitted_values <- JAGS_evaluate_formula(
    fit = fit,
    formula = ~ x1 * x2,
    parameter = "mu",
    data = formula_data,
    prior_list = formula_result$prior_list
  )
  expect_equal(unname(drop(fitted_values)), unname(stats::fitted(lm_fit)), tolerance = 1e-12)

  newdata <- data.frame(
    x1 = stats::quantile(data$x1, probs = c(0.10, 0.50, 0.90), names = FALSE),
    x2 = stats::quantile(data$x2, probs = c(0.20, 0.60, 0.80), names = FALSE)
  )
  scaled_newdata <- newdata
  for (variable in c("x1", "x2")) {
    scale_info <- formula_result$formula_scale[[paste0("mu_", variable)]]
    scaled_newdata[[variable]] <- (newdata[[variable]] - scale_info$mean) / scale_info$sd
  }

  predicted <- JAGS_evaluate_formula(
    fit = fit,
    formula = ~ x1 * x2,
    parameter = "mu",
    data = newdata,
    prior_list = formula_result$prior_list
  )
  expect_equal(
    unname(drop(predicted)),
    unname(stats::predict(lm_fit, newdata = scaled_newdata)),
    tolerance = 1e-12
  )
})

test_that("JAGS_evaluate_formula matches lm predictions for factors and no-intercept formulas", {

  factor_data <- data.frame(
    y = c(3.0, 3.8, 4.6, 5.4, 1.0, 1.7, 2.4, 3.1, 4.8, 5.7, 6.6, 7.5),
    x = rep(c(-1, 0, 1, 2), 3),
    group = factor(rep(c("b", "a", "c"), each = 4), levels = c("c", "a", "b"))
  )
  factor_prior <- list(
    intercept = prior("normal", list(0, 10)),
    x = prior("normal", list(0, 5)),
    group = prior_factor("normal", list(0, 5), contrast = "treatment"),
    "x:group" = prior_factor("normal", list(0, 5), contrast = "treatment")
  )
  factor_formula <- ~ x * group
  factor_result <- JAGS_formula(
    formula = factor_formula,
    parameter = "mu",
    data = factor_data[c("x", "group")],
    prior_list = factor_prior
  )
  factor_lm <- stats::lm(y ~ x * group, data = factor_data)
  factor_fit <- .jags_formula_posterior_from_lm(factor_result, factor_lm)

  factor_newdata <- data.frame(
    x = c(-.5, .5, 1.5),
    group = factor(c("b", "c", "a"), levels = c("b", "a", "c"))
  )
  factor_prediction <- JAGS_evaluate_formula(
    fit = factor_fit,
    formula = factor_formula,
    parameter = "mu",
    data = factor_newdata,
    prior_list = factor_result$prior_list
  )
  expect_equal(
    unname(drop(factor_prediction)),
    unname(stats::predict(factor_lm, newdata = factor_newdata)),
    tolerance = 1e-12
  )

  no_intercept_data <- data.frame(y = c(-2, -1, 1, 2, 4), x = c(-2, -1, 1, 2, 4))
  no_intercept_result <- JAGS_formula(
    formula = ~ x - 1,
    parameter = "mu",
    data = no_intercept_data["x"],
    prior_list = list(x = prior("normal", list(0, 5)))
  )
  no_intercept_lm <- stats::lm(y ~ x - 1, data = no_intercept_data)
  no_intercept_fit <- coda::mcmc(
    matrix(
      stats::coef(no_intercept_lm),
      nrow = 1,
      dimnames = list(NULL, JAGS_parameter_names("x", formula_parameter = "mu"))
    )
  )
  attr(no_intercept_fit, "formula_design") <- list(mu = no_intercept_result$formula_design)
  no_intercept_newdata <- data.frame(x = c(-3, 0, 3))
  no_intercept_prediction <- JAGS_evaluate_formula(
    fit = no_intercept_fit,
    formula = ~ x - 1,
    parameter = "mu",
    data = no_intercept_newdata,
    prior_list = no_intercept_result$prior_list
  )
  expect_equal(
    unname(drop(no_intercept_prediction)),
    unname(stats::predict(no_intercept_lm, newdata = no_intercept_newdata)),
    tolerance = 1e-12
  )

  scaled_no_intercept_result <- JAGS_formula(
    formula = ~ x - 1,
    parameter = "mu",
    data = no_intercept_data["x"],
    prior_list = list(x = prior("normal", list(0, 5))),
    formula_scale = list(x = TRUE)
  )
  scale_info <- scaled_no_intercept_result$formula_scale$mu_x
  scaled_slope <- 3
  scaled_posterior <- matrix(
    scaled_slope,
    nrow = 1,
    dimnames = list(NULL, "mu_x")
  )
  original_posterior <- BayesTools:::.bt_transform_scale_posterior(
    scaled_posterior,
    list(mu = scaled_no_intercept_result$formula_scale)
  )
  expect_true("mu_intercept" %in% colnames(original_posterior))
  expect_equal(
    unname(original_posterior[1, "mu_x"]),
    scaled_slope / scale_info$sd,
    tolerance = 1e-12
  )
  expect_equal(
    unname(original_posterior[1, "mu_intercept"]),
    -scaled_slope * scale_info$mean / scale_info$sd,
    tolerance = 1e-12
  )
  expect_equal(
    unname(original_posterior[1, "mu_intercept"] +
      original_posterior[1, "mu_x"] * no_intercept_newdata$x),
    scaled_slope * unname((no_intercept_newdata$x - scale_info$mean) / scale_info$sd),
    tolerance = 1e-12
  )
})

test_that("concrete full-rank factor contrasts survive neutral intercept replay", {

  data <- data.frame(
    group = factor(c("a", "b", "a"), levels = c("a", "b"))
  )
  attr(data$group, "contrasts") <- contr.independent(2L)
  formula <- ~ group
  model_frame <- stats::model.frame(formula, data = data)

  model_matrix <- BayesTools:::.bt_model_matrix(
    model_frame = model_frame,
    formula = formula,
    data = data
  )

  expect_equal(
    as.vector(model_matrix),
    as.vector(cbind(
      1,
      contr.independent(2L)[c(1L, 2L, 1L), , drop = FALSE]
    ))
  )
  expect_identical(dim(model_matrix), c(3L, 3L))
  expect_equal(
    attr(model_matrix, "assign"),
    c(0L, 1L, 1L)
  )
})

test_that("JAGS_evaluate_formula has stable semantics for aliased rank-deficient designs", {

  data <- bayestools_oracle_formula_design_data()
  formula_result <- JAGS_formula(
    formula = ~ x + x_alias,
    parameter = "mu",
    data = data,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1)),
      x_alias = prior("normal", list(0, 1))
    )
  )

  expect_true(any(formula_result$formula_design$aliased))
  expect_equal(formula_result$formula_design$rank, 2L)

  coefficients <- c(mu_intercept = 0.75, mu_x = 1.25, mu_x_alias = -0.5)
  fit <- coda::mcmc(
    matrix(coefficients, nrow = 1, dimnames = list(NULL, names(coefficients)))
  )
  attr(fit, "formula_design") <- list(mu = formula_result$formula_design)

  newdata <- data.frame(
    x = c(-2, 0, 3),
    x_alias = c(-4, 0, 6)
  )
  model_matrix <- stats::model.matrix(~ x + x_alias, data = newdata)

  expect_equal(
    unname(drop(JAGS_evaluate_formula(
      fit = fit,
      formula = ~ x + x_alias,
      parameter = "mu",
      data = newdata,
      prior_list = formula_result$prior_list
    ))),
    unname(drop(model_matrix %*% coefficients)),
    tolerance = 1e-12
  )

  off_alias_newdata <- data.frame(
    x = c(-2, 0, 3),
    x_alias = c(5, 6, 7)
  )
  off_alias_matrix <- stats::model.matrix(~ x + x_alias, data = off_alias_newdata)

  expect_equal(
    unname(drop(JAGS_evaluate_formula(
      fit = fit,
      formula = ~ x + x_alias,
      parameter = "mu",
      data = off_alias_newdata,
      prior_list = formula_result$prior_list
    ))),
    unname(drop(off_alias_matrix %*% coefficients)),
    tolerance = 1e-12
  )
})

test_that("JAGS_evaluate_formula preserves factor metadata and validates edge cases", {

  factor_data <- data.frame(
    y = c(3.0, 3.8, 4.6, 5.4, 1.0, 1.7, 2.4, 3.1, 4.8, 5.7, 6.6, 7.5),
    x = rep(c(-1, 0, 1, 2), 3),
    group = factor(rep(c("b", "a", "c"), each = 4), levels = c("c", "a", "b"))
  )
  prior_list <- list(
    intercept = prior("normal", list(0, 10)),
    x = prior("normal", list(0, 5)),
    group = prior_factor("normal", list(0, 5), contrast = "treatment")
  )
  formula <- ~ x + group

  formula_result <- JAGS_formula(
    formula = formula,
    parameter = "mu",
    data = factor_data[c("x", "group")],
    prior_list = prior_list
  )
  lm_fit <- stats::lm(y ~ x + group, data = factor_data)
  fit <- .jags_formula_posterior_from_lm(formula_result, lm_fit)

  reordered_newdata <- data.frame(
    x = c(-.5, .5, 1.5),
    group = factor(c("b", "c", "a"), levels = c("b", "a", "c"))
  )
  expected <- stats::predict(lm_fit, newdata = reordered_newdata)

  expect_equal(
    unname(drop(JAGS_evaluate_formula(
      fit = fit,
      formula = formula,
      parameter = "mu",
      data = reordered_newdata,
      prior_list = formula_result$prior_list
    ))),
    unname(expected),
    tolerance = 1e-12
  )

  character_newdata <- reordered_newdata
  character_newdata$group <- as.character(character_newdata$group)
  expect_equal(
    unname(drop(JAGS_evaluate_formula(
      fit = fit,
      formula = formula,
      parameter = "mu",
      data = character_newdata,
      prior_list = formula_result$prior_list
    ))),
    unname(expected),
    tolerance = 1e-12
  )

  prior_list_scaled <- formula_result$prior_list
  attr(prior_list_scaled$mu_x, "multiply_by") <- 0
  zero_x_samples <- as.matrix(fit)
  zero_x_samples[, "mu_x"] <- 0
  zero_x_fit <- coda::mcmc(zero_x_samples)
  attr(zero_x_fit, "formula_design") <- attr(fit, "formula_design")

  expect_equal(
    unname(JAGS_evaluate_formula(
      fit = fit,
      formula = formula,
      parameter = "mu",
      data = reordered_newdata,
      prior_list = prior_list_scaled
    )),
    unname(JAGS_evaluate_formula(
      fit = zero_x_fit,
      formula = formula,
      parameter = "mu",
      data = reordered_newdata,
      prior_list = formula_result$prior_list
    )),
    tolerance = 1e-12
  )

  expect_error(
    JAGS_evaluate_formula(
      fit = fit,
      formula = formula,
      parameter = "mu",
      data = reordered_newdata["x"],
      prior_list = formula_result$prior_list
    ),
    "predictor variable is missing"
  )

  missing_group_prior <- formula_result$prior_list[
    setdiff(names(formula_result$prior_list), "mu_group")
  ]
  expect_error(
    JAGS_evaluate_formula(
      fit = fit,
      formula = formula,
      parameter = "mu",
      data = reordered_newdata,
      prior_list = missing_group_prior
    ),
    "prior distribution for the 'group' term is missing"
  )

  missing_intercept_prior <- formula_result$prior_list[
    setdiff(names(formula_result$prior_list), "mu_intercept")
  ]
  expect_error(
    JAGS_evaluate_formula(
      fit = fit,
      formula = formula,
      parameter = "mu",
      data = reordered_newdata,
      prior_list = missing_intercept_prior
    ),
    "prior distribution for the 'intercept' term is missing"
  )

  bad_factor_newdata <- reordered_newdata
  bad_factor_newdata$group <- factor(c("d", "b", "d"), levels = c("b", "d"))
  expect_error(
    JAGS_evaluate_formula(
      fit = fit,
      formula = formula,
      parameter = "mu",
      data = bad_factor_newdata,
      prior_list = formula_result$prior_list
    ),
    "Levels specified"
  )

  bad_character_newdata <- reordered_newdata
  bad_character_newdata$group <- c("d", "b", "d")
  expect_error(
    JAGS_evaluate_formula(
      fit = fit,
      formula = formula,
      parameter = "mu",
      data = bad_character_newdata,
      prior_list = formula_result$prior_list
    ),
    "Levels specified"
  )
})
