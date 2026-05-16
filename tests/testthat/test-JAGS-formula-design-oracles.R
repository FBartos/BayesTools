skip_if_not_test_profile("unit")

# ============================================================================ #
# TEST FILE: JAGS Formula Design Oracles
# ============================================================================ #
#
# PURPOSE:
#   Deterministic model.matrix() oracles for JAGS_formula() design metadata.
#
# TAGS: @evaluation, @formula, @model-matrix
# ============================================================================ #

.jags_formula_oracle_expected_data <- function(data, factor_contrasts = list(),
                                               formula_scale = NULL) {
  out <- data

  for (factor_name in names(factor_contrasts)) {
    stats::contrasts(out[[factor_name]]) <- factor_contrasts[[factor_name]]
  }

  continuous <- names(out)[vapply(out, is.numeric, logical(1))]
  for (variable in continuous) {
    should_scale <- isTRUE(formula_scale)
    if (is.list(formula_scale) && !is.null(formula_scale[[variable]])) {
      should_scale <- isTRUE(formula_scale[[variable]])
    }
    if (should_scale) {
      out[[variable]] <- (out[[variable]] - mean(out[[variable]])) / stats::sd(out[[variable]])
    }
  }

  out
}

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

  posterior <- matrix(posterior, nrow = 1, dimnames = list(NULL, names(posterior)))
  coda::mcmc(posterior)
}

test_that("JAGS_formula design metadata matches stats::model.matrix for fixed effects", {

  df <- bayestools_oracle_formula_design_data()

  cases <- list(
    continuous_main_effects = list(
      formula = ~ x + z,
      prior_list = list(
        intercept = prior("normal", list(0, 1)),
        x = prior("normal", list(0, 1)),
        z = prior("normal", list(0, 1))
      ),
      formula_scale = NULL,
      factor_contrasts = list()
    ),
    continuous_interaction_scaled = list(
      formula = ~ x * z,
      prior_list = list(
        intercept = prior("normal", list(0, 1)),
        x = prior("normal", list(0, 1)),
        z = prior("normal", list(0, 1)),
        "x:z" = prior("normal", list(0, 1))
      ),
      formula_scale = list(x = TRUE),
      factor_contrasts = list()
    ),
    treatment_factor_interaction_unused_level = list(
      formula = ~ x * group,
      prior_list = list(
        intercept = prior("normal", list(0, 1)),
        x = prior("normal", list(0, 1)),
        group = prior_factor("normal", list(0, 1), contrast = "treatment"),
        "x:group" = prior_factor("normal", list(0, 1), contrast = "treatment")
      ),
      formula_scale = NULL,
      factor_contrasts = list(group = "contr.treatment")
    ),
    reordered_treatment_factor = list(
      formula = ~ group_reordered,
      prior_list = list(
        intercept = prior("normal", list(0, 1)),
        group_reordered = prior_factor("normal", list(0, 1), contrast = "treatment")
      ),
      formula_scale = NULL,
      factor_contrasts = list(group_reordered = "contr.treatment")
    ),
    independent_factor_no_intercept = list(
      formula = ~ group - 1,
      expected_formula = formula_add_intercept(~ group - 1),
      prior_list = list(
        group = prior_factor("normal", list(0, 1), contrast = "independent")
      ),
      formula_scale = NULL,
      factor_contrasts = list(group = "contr.independent")
    ),
    factor_by_factor_orthonormal_meandif = list(
      formula = ~ a * b,
      prior_list = list(
        intercept = prior("normal", list(0, 1)),
        a = prior_factor("mnormal", list(0, 1), contrast = "orthonormal"),
        b = prior_factor("mnormal", list(0, 1), contrast = "meandif"),
        "a:b" = prior_factor("mnormal", list(0, 1), contrast = "orthonormal")
      ),
      formula_scale = NULL,
      factor_contrasts = list(a = "contr.orthonormal", b = "contr.meandif")
    ),
    aliased_continuous_design = list(
      formula = ~ x + x_alias,
      prior_list = list(
        intercept = prior("normal", list(0, 1)),
        x = prior("normal", list(0, 1)),
        x_alias = prior("normal", list(0, 1))
      ),
      formula_scale = NULL,
      factor_contrasts = list()
    )
  )

  for (case_name in names(cases)) {
    case <- cases[[case_name]]
    result <- JAGS_formula(
      formula = case$formula,
      parameter = "mu",
      data = df,
      prior_list = case$prior_list,
      formula_scale = case$formula_scale
    )

    expected_formula <- if (is.null(case$expected_formula)) case$formula else case$expected_formula
    expected_data <- .jags_formula_oracle_expected_data(
      data = df,
      factor_contrasts = case$factor_contrasts,
      formula_scale = case$formula_scale
    )

    expect_jags_formula_design_matches_model_matrix(
      result$formula_design,
      formula = expected_formula,
      data = expected_data
    )
    expect_equal(result$formula_design$formula, expected_formula, ignore_formula_env = TRUE)

    expected_matrix <- stats::model.matrix(expected_formula, data = expected_data)
    expected_rank <- qr(expected_matrix)$rank
    expected_aliased <- rep(FALSE, ncol(expected_matrix))
    if (expected_rank < ncol(expected_matrix)) {
      expected_aliased[qr(expected_matrix)$pivot[(expected_rank + 1L):ncol(expected_matrix)]] <- TRUE
    }
    names(expected_aliased) <- gsub(":", "__xXx__", colnames(expected_matrix), fixed = TRUE)

    expect_equal(result$formula_design$rank, expected_rank)
    expect_equal(result$formula_design$aliased, expected_aliased)
  }
})

test_that("JAGS_formula rejects non-syntactic formula names explicitly", {
  data <- bayestools_oracle_nonsyntactic_formula_data()
  prior_list <- setNames(
    list(prior("normal", list(0, 1)), prior("normal", list(0, 1))),
    c("intercept", "x weird")
  )
  nonsyntactic_formula <- stats::as.formula(
    paste0("~ ", intToUtf8(96), "x weird", intToUtf8(96))
  )

  expect_error(
    JAGS_formula(
      formula = nonsyntactic_formula,
      parameter = "mu",
      data = data,
      prior_list = prior_list
    ),
    "predictor variable is missing"
  )
})

test_that("JAGS_formula rejects malformed formula_scale input", {
  data <- data.frame(
    x = c(-2, -1, 0, 1, 2),
    f = factor(c("a", "b", "a", "b", "a"))
  )
  prior_list <- list(
    intercept = prior("normal", list(0, 1)),
    x = prior("normal", list(0, 1)),
    f = prior_factor("normal", list(0, 1), contrast = "treatment")
  )

  expect_error(
    JAGS_formula(~ x + f, "mu", data, prior_list, formula_scale = c(TRUE, FALSE)),
    "length '1'"
  )
  expect_error(
    JAGS_formula(~ x + f, "mu", data, prior_list, formula_scale = list(TRUE)),
    "named list"
  )
  expect_error(
    JAGS_formula(~ x + f, "mu", data, prior_list, formula_scale = list(x = TRUE, x = FALSE)),
    "unique"
  )
  expect_error(
    JAGS_formula(~ x + f, "mu", data, prior_list, formula_scale = list(xx = TRUE)),
    "not predictor variables"
  )
  expect_error(
    JAGS_formula(~ x + f, "mu", data, prior_list, formula_scale = list(f = TRUE)),
    "Only continuous predictors"
  )
  expect_error(
    JAGS_formula(~ x + f, "mu", data, prior_list, formula_scale = list(x = NA)),
    "cannot contain NA"
  )
  expect_error(
    JAGS_formula(~ x + f, "mu", data, prior_list, formula_scale = list(x = 1)),
    "logical vector"
  )

  expect_error(
    JAGS_formula(~ x + f, "mu", data, prior_list, formula_scale = list(x = TRUE, f = FALSE)),
    NA
  )
})

test_that("JAGS_formula rejects missing fixed-effect predictors before row dropping", {
  data <- data.frame(x = c(1, NA_real_, 3))
  prior_list <- list(
    intercept = prior("normal", list(0, 1)),
    x = prior("normal", list(0, 1))
  )

  expect_error(
    JAGS_formula(~ x, "mu", data, prior_list),
    "Formula predictors contain missing values.",
    fixed = TRUE
  )

  formula_result <- JAGS_formula(~ x, "mu", data.frame(x = c(1, 2, 3)), prior_list)
  posterior <- coda::mcmc(
    matrix(c(mu_intercept = 0, mu_x = 1), nrow = 1, dimnames = list(NULL, c("mu_intercept", "mu_x")))
  )
  expect_error(
    JAGS_evaluate_formula(posterior, ~ x, "mu", data, formula_result$prior_list),
    "Formula predictors contain missing values.",
    fixed = TRUE
  )
})

test_that("JAGS_formula handles character and interaction-only factor predictors", {
  character_data <- data.frame(g = c("a", "b", "a", "c"))
  factor_data <- data.frame(g = factor(character_data$g))
  prior_list <- list(
    intercept = prior("normal", list(0, 1)),
    g = prior_factor("normal", list(0, 1), contrast = "treatment")
  )

  character_result <- JAGS_formula(~ g, "mu", character_data, prior_list)
  factor_result <- JAGS_formula(~ g, "mu", factor_data, prior_list)

  expect_equal(character_result$formula_design$model_matrix, factor_result$formula_design$model_matrix)
  expect_equal(character_result$formula_design$xlevels, factor_result$formula_design$xlevels)

  interaction_data <- data.frame(
    x = c(-1, 0, 1, 2),
    g = factor(c("a", "b", "a", "b"))
  )
  interaction_prior <- list(
    intercept = prior("normal", list(0, 1)),
    "x:g" = prior_factor("normal", list(0, 1), contrast = "treatment")
  )
  interaction_result <- JAGS_formula(~ x:g, "mu", interaction_data, interaction_prior)
  expect_jags_formula_design_matches_model_matrix(
    interaction_result$formula_design,
    formula = ~ x:g,
    data = interaction_data
  )
})

test_that("formula expression terms are parsed structurally", {
  expect_equal(.extract_expressions(~ expression(log(x))), list("log(x)"))
  expect_equal(.extract_expressions(y ~ z + expression(log(x)) + expression(exp(b))), list("log(x)", "exp(b)"))
  expect_equal(.remove_expressions(y ~ expression(log(x))), formula(y ~ 1), ignore_formula_env = TRUE)
  expect_equal(.remove_expressions(~ z + expression(log(x))), formula(~ z), ignore_formula_env = TRUE)

  expression_result <- JAGS_formula(
    y ~ expression(log(x)),
    "mu",
    data.frame(x = c(1, 2, 3)),
    list(intercept = prior("normal", list(0, 1)))
  )
  expect_equal(expression_result$formula_design$transformed_terms, list("log(x)"))
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

  fit <- coda::mcmc(
    matrix(
      scaled_coefficients,
      nrow = 1,
      dimnames = list(NULL, names(scaled_coefficients))
    )
  )
  attr(fit, "formula_scale") <- list(mu = formula_result$formula_scale)

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
    "prior distribution.*missing"
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

test_that("JAGS_formula records scaling and JAGS data columns by model term", {

  df <- bayestools_oracle_formula_design_data()
  result <- JAGS_formula(
    formula = ~ x * group,
    parameter = "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1)),
      group = prior_factor("normal", list(0, 1), contrast = "treatment"),
      "x:group" = prior_factor("normal", list(0, 1), contrast = "treatment")
    ),
    formula_scale = list(x = TRUE)
  )
  design <- result$formula_design

  expected_data <- .jags_formula_oracle_expected_data(
    df,
    factor_contrasts = list(group = "contr.treatment"),
    formula_scale = list(x = TRUE)
  )
  expected <- stats::model.matrix(~ x * group, data = expected_data)

  expect_equal(result$formula_scale$mu_x$mean, mean(df$x))
  expect_equal(result$formula_scale$mu_x$sd, stats::sd(df$x))
  expect_equal(unname(result$data$mu_data_x), unname(expected[, "x"]))
  expect_equal(unname(result$data$mu_data_group), unname(expected[, c("groupb", "groupc", "groupunused")]))
  expect_equal(
    unname(result$data$mu_data_x__xXx__group),
    unname(expected[, c("x:groupb", "x:groupc", "x:groupunused")])
  )
  expect_equal(
    unlist(design$jags_data_names, use.names = FALSE),
    c("mu_data_x", "mu_data_group", "mu_data_x__xXx__group")
  )
})

test_that("random-effect design exposes grouping maps and correlated syntax", {

  df <- data.frame(
    x = c(-1, 0, 1, 2, -2, 3),
    idx = factor(rep(c("t1", "t2"), 3), levels = c("t1", "t2")),
    id = factor(c("b", "a", "b", "c", "a", "c"), levels = c("a", "b", "c"))
  )
  random_term <- BayesTools:::.bt_parse_random_effects(~ 1 + x + (1 + x || id))$terms[[1]]

  expect_s3_class(random_term, "BayesTools_random_effect_term")
  expect_equal(random_term$group_label, "id")
  expect_equal(random_term$structure, "diag")
  expect_false("covariance" %in% names(random_term))
  expect_equal(attr(random_term, "structure"), "diag")
  expect_null(attr(random_term, "covariance"))
  expect_true(isTRUE(random_term$independent))
  expect_equal(random_term$term_formula, ~ 1 + x, ignore_formula_env = TRUE)
  expect_equal(attr(random_term, "grouping_factor"), "id")
  expect_true(isTRUE(attr(random_term, "independent")))

  reordered_wrapper <- BayesTools:::.bt_parse_random_effects(
    ~ 1 + x + random(name = "study", covariance = "diag", 1 + x | id)
  )$terms[[1]]
  expect_equal(reordered_wrapper$block_name, "study")
  expect_equal(reordered_wrapper$structure, "diag")
  expect_equal(reordered_wrapper$term_formula, ~ 1 + x, ignore_formula_env = TRUE)

  reordered_special <- BayesTools:::.bt_parse_random_effects(
    ~ 1 + x + diag(name = "study", 1 + x | id)
  )$terms[[1]]
  expect_equal(reordered_special$block_name, "study")
  expect_equal(reordered_special$structure, "diag")

  explicit_diag_double_bar <- BayesTools:::.bt_parse_random_effects(
    ~ 1 + x + random(1 + x || id, name = "study", covariance = "diag")
  )$terms[[1]]
  expect_equal(explicit_diag_double_bar$block_name, "study")
  expect_equal(explicit_diag_double_bar$structure, "diag")

  expect_error(
    BayesTools:::.bt_parse_random_effects(
      ~ 1 + x + random(1 + x || id, name = "study", covariance = "us")
    ),
    "cannot be combined with '||' syntax",
    fixed = TRUE
  )
  expect_error(
    BayesTools:::.bt_parse_random_effects(
      ~ 1 + x +
        random(1 | drug, name = "drug", covariance = "diag") +
        us(1 + x || id)
    ),
    "cannot be combined with '||' syntax",
    fixed = TRUE
  )

  bare_and_wrapper <- BayesTools:::.bt_parse_random_effects(
    ~ 1 + (1 | id) + random(1 | id, name = "study", covariance = "diag")
  )
  expect_equal(
    vapply(bare_and_wrapper$terms, function(term) term$block_name, character(1)),
    c("id", "study")
  )
  expect_equal(
    vapply(bare_and_wrapper$terms, function(term) term$structure, character(1)),
    c("us", "diag")
  )

  bare_and_special <- BayesTools:::.bt_parse_random_effects(
    ~ 1 + (1 | id) + diag(name = "study", 1 | id)
  )
  expect_equal(
    vapply(bare_and_special$terms, function(term) term$block_name, character(1)),
    c("id", "study")
  )
  expect_equal(
    vapply(bare_and_special$terms, function(term) term$structure, character(1)),
    c("us", "diag")
  )

  ordered_wrapper_terms <- BayesTools:::.bt_parse_random_effects(
    ~ 1 +
      random(1 | id, name = "study", covariance = "diag") +
      diag(1 | drug, name = "drug")
  )
  expect_equal(
    vapply(ordered_wrapper_terms$terms, function(term) term$block_name, character(1)),
    c("study", "drug")
  )

  reversed_wrapper_terms <- BayesTools:::.bt_parse_random_effects(
    ~ 1 +
      diag(1 | drug, name = "drug") +
      random(1 | id, name = "study", covariance = "diag")
  )
  expect_equal(
    vapply(reversed_wrapper_terms$terms, function(term) term$block_name, character(1)),
    c("drug", "study")
  )

  formula_env <- new.env(parent = baseenv())
  formula_env$make_group_for_random_effect_test <- function(x) paste0("g_", x)
  grouped_formula <- ~ 1 + random(1 | make_group_for_random_effect_test(id), name = "derived", covariance = "diag")
  environment(grouped_formula) <- formula_env
  grouped_term <- BayesTools:::.bt_parse_random_effects(grouped_formula)$terms[[1]]
  expect_equal(
    BayesTools:::.bt_random_group_values(grouped_term, data.frame(id = c("a", "b"))),
    c("g_a", "g_b")
  )

  expect_error(
    BayesTools:::.bt_parse_random_effects(~ 1 + random(1 | id, TRUE, name = "id")),
    "exactly one unnamed term",
    fixed = TRUE
  )
  expect_error(
    BayesTools:::.bt_parse_random_effects(~ 1 + diag(1 | id, TRUE)),
    "exactly one unnamed term",
    fixed = TRUE
  )

  correlated <- BayesTools:::.JAGS_random_effect_formula(
    formula = BayesTools:::.bt_parse_random_effects(~ 1 + x + (1 + x | id))$terms[[1]],
    parameter = "mu",
    data = df,
    prior_random = prior_random(
      id = random_block(
        sd = prior("gamma", list(2, 2)),
        cor = prior_lkj(eta = 1),
        monitor = random_monitor(latent = FALSE, coefficients = FALSE, correlation = TRUE)
      )
    )
  )

  expect_equal(correlated$random_effect$structure, "us")
  expect_match(paste(correlated$random_syntax, collapse = "\n"), "dbt_lkj_cpc", fixed = TRUE)
  expect_equal(correlated$jags_modules, "BayesTools")
  expect_equal(
    correlated$add_parameters,
    c(
      "mu__xREx__id_xRE_CORx_L",
      "mu__xREx__id_xRE_CORx_R",
      "mu__xREx__id_xRE_CORx_lkj_u[1]"
    )
  )
})

test_that("random-effect terms use structure as canonical internal metadata", {

  df <- data.frame(
    x = c(-1, 0, 1, 2, -2, 3),
    id = factor(c("b", "a", "b", "c", "a", "c"), levels = c("a", "b", "c"))
  )

  parsed <- BayesTools:::.bt_parse_random_effects(
    ~ 1 + x +
      (1 + x | id) +
      random(1 + x || id, name = "independent", covariance = "diag") +
      ar1(x | id, name = "ordered")
  )

  formula_result <- JAGS_formula(
    formula = ~ 1 + x + random(1 + x | id, name = "id", covariance = "us"),
    parameter = "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    ),
    prior_random = prior_random(
      id = random_block(
        sd = prior("gamma", list(2, 2)),
        cor = prior_lkj(eta = 1)
      )
    )
  )

  terms <- c(parsed$terms, formula_result$formula_design$random_effects)
  expect_true(length(terms) > 0L)

  for (random_term in terms) {
    expect_s3_class(random_term, "BayesTools_random_effect_term")
    expect_true("structure" %in% names(random_term))
    expect_false("covariance" %in% names(random_term))
    expect_true(nzchar(random_term$structure))
    expect_equal(attr(random_term, "structure"), random_term$structure)
    expect_null(attr(random_term, "covariance"))
  }

  malformed_term <- terms[[1]]
  malformed_term$covariance <- malformed_term$structure
  malformed_term$structure <- NULL
  attr(malformed_term, "covariance") <- malformed_term$covariance
  attr(malformed_term, "structure") <- NULL
  expect_error(
    BayesTools:::.bt_random_effect_structure(malformed_term),
    "missing canonical 'random_term\\$structure'"
  )
})

test_that("JAGS_formula random-effect design exposes grouping maps and public variable names", {

  df <- data.frame(
    x = c(-1, 0, 1, 2, -2, 3),
    id = factor(c("b", "a", "b", "c", "a", "c"), levels = c("c", "a", "b"))
  )
  prior_list <- list(
    intercept = prior("normal", list(0, 1)),
    x = prior("normal", list(0, 1))
  )
  random_prior_diag <- prior_random(
    id = random_block(sd = prior("gamma", list(2, 2)))
  )
  random_prior_us <- prior_random(
    id = random_block(
      sd = prior("gamma", list(2, 2)),
      cor = prior_lkj(eta = 1)
    )
  )

  result <- JAGS_formula(
    formula = ~ 1 + x + (1 + x || id),
    parameter = "mu",
    data = df,
    prior_list = prior_list,
    prior_random = random_prior_diag
  )

  expected_map <- match(as.character(df$id), levels(df$id))

  expect_equal(result$data$mu__xREx__id_xRE_MAPx, expected_map)
  expect_equal(dim(result$data$mu__xREx__id_xRE_DATAx), c(6L, 2L))
  expect_equal(colnames(result$data$mu__xREx__id_xRE_DATAx), c("(Intercept)", "x"))
  expect_equal(unname(result$data$mu__xREx__id_xRE_DATAx[, "x"]), df$x)

  expect_equal(
    result$formula_design$jags_data_names[["__xREx__id"]],
    c("mu__xREx__id_xRE_DATAx", "mu__xREx__id_xRE_MAPx")
  )
  expect_true(isTRUE(attr(result$formula_design$random_effects[[1]], "independent")))
  expect_equal(attr(result$formula_design$random_effects[[1]], "grouping_factor"), "id")
  expect_equal(result$formula_design$random_effects[[1]]$group_map, expected_map)
  expect_equal(result$formula_design$random_effects[[1]]$group_levels, levels(df$id))
  expect_equal(result$formula_design$random_effects[[1]]$column_names, c("(Intercept)", "x"))
  expect_equal(result$formula_design$random_effects[[1]]$jags_data_names, c("mu__xREx__id_xRE_DATAx", "mu__xREx__id_xRE_MAPx"))
  expect_s3_class(
    result$formula_design$random_effects[[1]]$sd_leaves,
    "BayesTools_random_effect_sd_leaves"
  )
  expect_equal(
    result$formula_design$random_effects[[1]]$sd_leaves$leaf_names_by_column,
    result$formula_design$random_effects[[1]]$sd_parameter_names
  )
  expect_equal(
    result$formula_design$random_effects[[1]]$sd_leaves$leaf_terms_by_column,
    c("intercept", "x")
  )
  expect_equal(
    names(result$prior_list),
    c("mu_intercept", "mu_x", "mu__xREx__id_intercept", "mu__xREx__id_x")
  )
  expect_true(isTRUE(attr(result$prior_list$mu__xREx__id_x, "random_sd")))
  expect_equal(attr(result$prior_list$mu__xREx__id_x, "random_factor"), "id")
  expect_match(result$formula_syntax, "for\\(i in 1:3\\)")
  expect_match(result$formula_syntax, "dmnorm\\(rep\\(0, 2\\)")

  diag_result <- JAGS_formula(
    formula = ~ 1 + x + diag(1 + x | id),
    parameter = "mu",
    data = df,
    prior_list = prior_list,
    prior_random = random_prior_diag
  )

  expect_equal(
    diag_result$data$mu__xREx__id_xRE_MAPx,
    result$data$mu__xREx__id_xRE_MAPx
  )
  expect_equal(
    unname(diag_result$data$mu__xREx__id_xRE_DATAx),
    unname(result$data$mu__xREx__id_xRE_DATAx)
  )
  expect_equal(diag_result$formula_design$random_effects[[1]]$structure, "diag")

  us_result <- JAGS_formula(
    formula = ~ 1 + x + (1 + x | id),
    parameter = "mu",
    data = df,
    prior_list = prior_list,
    prior_random = random_prior_us
  )
  expect_equal(us_result$formula_design$random_effects[[1]]$structure, "us")
  expect_equal(us_result$jags_modules, "BayesTools")
  expect_true(any(grepl("xRE_CORx_R", us_result$add_parameters, fixed = TRUE)))

  explicit_us_result <- JAGS_formula(
    formula = ~ 1 + x + us(1 + x | id),
    parameter = "mu",
    data = df,
    prior_list = prior_list,
    prior_random = random_prior_us
  )
  expect_equal(explicit_us_result$formula_design$random_effects[[1]]$structure, "us")

  single_us_result <- JAGS_formula(
    formula = ~ 1 + (1 | id),
    parameter = "mu",
    data = df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(id = random_block(sd = prior("gamma", list(2, 2))))
  )
  expect_equal(single_us_result$jags_modules, character())
  expect_null(single_us_result$formula_design$random_effects[[1]]$correlation)
  expect_false(grepl("xRE_CORx", single_us_result$formula_syntax, fixed = TRUE))
  expect_error(
    JAGS_formula(
      formula = ~ 1 + (1 | id),
      parameter = "mu",
      data = df,
      prior_list = list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(id = random_block(
        sd = prior("gamma", list(2, 2)),
        cor = prior_lkj(eta = 1)
      ))
    ),
    "Single-column random-effect structure 'us' has no correlation parameter",
    fixed = TRUE
  )

  expect_error(
    JAGS_formula(
      formula = ~ 1 + x + diag(1 + x | id, unknown = TRUE),
      parameter = "mu",
      data = df,
      prior_list = prior_list,
      prior_random = random_prior_diag
    ),
    "The 'diag' random-effect covariance structure does not support extra arguments yet.",
    fixed = TRUE
  )
  expect_error(
    JAGS_formula(
      formula = ~ 1 + x + random(1 | id, name = "xRE_Zx_block"),
      parameter = "mu",
      data = df,
      prior_list = prior_list,
      prior_random = prior_random(
        xRE_Zx_block = random_block(sd = prior("gamma", list(2, 2)))
      )
    ),
    "internally used",
    fixed = TRUE
  )
  expect_error(
    JAGS_formula(
      formula = ~ 1 + x + cs(1 + x | id, hom = FALSE),
      parameter = "mu",
      data = df,
      prior_list = list(
        intercept = prior("normal", list(0, 1)),
        x = prior("normal", list(0, 1))
      ),
      prior_random = prior_random(id = random_block(sd = prior("gamma", list(2, 2))))
    ),
    "Use 'hcs' for heteroscedastic standard deviations",
    fixed = TRUE
  )
  expect_error(
    JAGS_formula(
      formula = ~ 1 + x + ar1(1 + x | id, hom = FALSE),
      parameter = "mu",
      data = df,
      prior_list = list(
        intercept = prior("normal", list(0, 1)),
        x = prior("normal", list(0, 1))
      ),
      prior_random = prior_random(id = random_block(sd = prior("gamma", list(2, 2))))
    ),
    "Use 'har' for heteroscedastic standard deviations",
    fixed = TRUE
  )
  expect_error(
    JAGS_formula(
      formula = ~ 1 + x + diag(1 + x | id, hom = NA),
      parameter = "mu",
      data = df,
      prior_list = list(
        intercept = prior("normal", list(0, 1)),
        x = prior("normal", list(0, 1))
      ),
      prior_random = prior_random(id = random_block(sd = prior("gamma", list(2, 2))))
    ),
    "'hom' must be TRUE or FALSE",
    fixed = TRUE
  )
  expect_error(
    JAGS_formula(
      formula = ~ 1 + x + hcs(1 + x | id, hom = TRUE),
      parameter = "mu",
      data = df,
      prior_list = list(
        intercept = prior("normal", list(0, 1)),
        x = prior("normal", list(0, 1))
      ),
      prior_random = prior_random(id = random_block(sd = prior("gamma", list(2, 2))))
    ),
    "already heteroscedastic",
    fixed = TRUE
  )
  expect_error(
    JAGS_formula(
      formula = ~ 1 + x + us(1 + x | id, hom = TRUE),
      parameter = "mu",
      data = df,
      prior_list = list(
        intercept = prior("normal", list(0, 1)),
        x = prior("normal", list(0, 1))
      ),
      prior_random = prior_random(id = random_block(sd = prior("gamma", list(2, 2))))
    ),
    "Homogeneous random-effect standard deviations are not supported yet",
    fixed = TRUE
  )
})

test_that("prior_random rejects unsupported and ignored production settings", {

  sd_prior <- prior("normal", list(0, 1), truncation = list(lower = 0, upper = Inf))
  rho_prior <- prior("normal", list(0, 0.5))
  lkj_prior <- prior_lkj(eta = 2)
  df <- data.frame(
    x = c(-1, 0, 1, 2, -2, 3),
    idx = factor(rep(c("t1", "t2"), 3), levels = c("t1", "t2")),
    id = factor(c("b", "a", "b", "c", "a", "c"), levels = c("a", "b", "c"))
  )
  fixed_priors <- list(
    intercept = prior("normal", list(0, 1)),
    x = prior("normal", list(0, 1))
  )
  group_as_data_frame <- function(x){
    data.frame(group = x)
  }

  expect_error(
    JAGS_formula(
      formula = ~ 1 + diag(1 | id),
      parameter = "mu",
      data = transform(df, id = factor(c("b", NA, "b", "c", "a", "c"))),
      prior_list = list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(id = random_block(sd = sd_prior))
    ),
    "must not contain missing values",
    fixed = TRUE
  )
  expect_error(
    JAGS_formula(
      formula = ~ 1 + random(1 | group_as_data_frame(id), name = "bad_group", covariance = "diag"),
      parameter = "mu",
      data = df,
      prior_list = list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(bad_group = random_block(sd = sd_prior))
    ),
    "must evaluate to one atomic value per row of data",
    fixed = TRUE
  )
  expect_error(
    JAGS_formula(
      formula = ~ 1 + random(1 | id[1], name = "bad_group", covariance = "diag"),
      parameter = "mu",
      data = df,
      prior_list = list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(bad_group = random_block(sd = sd_prior))
    ),
    "must evaluate to one value per row of data",
    fixed = TRUE
  )
  expect_error(
    JAGS_formula(
      formula = ~ 1 + random(1 | unknown_group_function(id), name = "bad_group", covariance = "diag"),
      parameter = "mu",
      data = df,
      prior_list = list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(bad_group = random_block(sd = sd_prior))
    ),
    "Could not evaluate random-effect grouping expression",
    fixed = TRUE
  )
  expect_error(
    prior_random(random_block(sd = sd_prior)),
    "must be named",
    fixed = TRUE
  )
  expect_error(
    prior_random(
      study = random_block(sd = sd_prior),
      study = random_block(sd = sd_prior)
    ),
    "must be unique",
    fixed = TRUE
  )
  expect_error(
    prior_random(study = sd_prior),
    "must be created with random_block",
    fixed = TRUE
  )
  expect_error(
    JAGS_formula(
      formula = ~ 1 + diag(1 | id),
      parameter = "mu",
      data = df,
      prior_list = list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(id = random_block(sd = prior("point", list(-1))))
    ),
    "point mass must be nonnegative",
    fixed = TRUE
  )
  expect_error(
    JAGS_formula(
      formula = ~ 1 + diag(0 + missing_x | id),
      parameter = "mu",
      data = df["id"],
      prior_list = list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(id = random_block(sd = sd_prior))
    ),
    "predictor variable is missing",
    fixed = TRUE
  )
  expect_error(
    JAGS_formula(
      formula = ~ 1 + diag(0 + x | id),
      parameter = "mu",
      data = transform(df, x = c(-1, NA, 1, 2, -2, 3)),
      prior_list = list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(id = random_block(sd = sd_prior))
    ),
    "missing predictor values",
    fixed = TRUE
  )

  allocation_prior <- random_variance_allocation(
    terms = c("study", "drug"),
    sd = sd_prior,
    allocation = prior("dirichlet", list(alpha = c(2, 3)))
  )
  expect_s3_class(allocation_prior, "random_variance_allocation")
  expect_s3_class(prior_random(allocation = allocation_prior), "prior_random")
  named_dot_allocation <- prior_random(total = allocation_prior)
  expect_equal(
    names(BayesTools:::.bt_random_allocation_list(named_dot_allocation$allocation)),
    "total"
  )
  expect_error(
    prior_random(`bad-name` = allocation_prior),
    "letters, numbers, and underscores",
    fixed = TRUE
  )
  expect_error(
    prior_random(allocation = stats::setNames(list(allocation_prior), "bad-name")),
    "letters, numbers, and underscores",
    fixed = TRUE
  )
  expect_error(
    random_variance_allocation(
      terms = c("study", "study"),
      sd = sd_prior,
      allocation = prior("dirichlet", list(alpha = c(1, 1)))
    ),
    "must be unique",
    fixed = TRUE
  )
  expect_error(
    random_variance_allocation(
      terms = "study",
      sd = sd_prior
    ),
    "at least two 'terms'",
    fixed = TRUE
  )
  expect_error(
    random_variance_allocation(
      terms = c("study", "drug"),
      sd = sd_prior,
      allocation = prior("beta", list(1, 1))
    ),
    "Dirichlet simplex priors",
    fixed = TRUE
  )
  expect_error(
    random_variance_allocation(
      terms = c("study", "drug"),
      sd = sd_prior,
      allocation = prior("dirichlet", list(alpha = c(1, 1, 1)))
    ),
    "dimension must match",
    fixed = TRUE
  )
  expect_error(
    random_variance_allocation(
      terms = c("study", "drug"),
      sd = sd_prior,
      name = "bad-name"
    ),
    "letters, numbers, and underscores",
    fixed = TRUE
  )
  expect_error(
    prior_random(sd = sd_prior, allocation = list(total = sd_prior)),
    "entries must be created with random_variance_allocation",
    fixed = TRUE
  )
  expect_error(
    random_block(allocation = list(total = sd_prior)),
    "Block-local variance allocation priors are not implemented",
    fixed = TRUE
  )
  expect_error(
    prior_random(sd = sd_prior, new_levels = random_new_levels(allow = TRUE)),
    "New-level random-effect prediction is not implemented yet",
    fixed = TRUE
  )
  expect_error(
    random_covariance(rho = prior_factor("normal", list(0, 1), contrast = "treatment")),
    "ordinary scalar prior",
    fixed = TRUE
  )
  expect_error(
    random_covariance(rho = prior("mnormal", list(mean = 0, sd = 1, K = 2))),
    "ordinary scalar prior",
    fixed = TRUE
  )
  expect_error(
    random_covariance(rho = prior_none(), rho_scale = "logit"),
    "cannot use prior_none",
    fixed = TRUE
  )
  expect_error(
    random_covariance(
      rho = prior_mixture(list(
        prior_factor("normal", list(0, 1), contrast = "treatment"),
        prior("point", list(0))
      ))
    ),
    "ordinary scalar prior",
    fixed = TRUE
  )
  expect_error(
    random_covariance(structure = "diag", eta = 2),
    "structure 'diag' has no correlation parameter",
    fixed = TRUE
  )
  expect_error(
    random_covariance(eta = 2, cor = lkj_prior),
    "explicit 'cor' or 'rho'",
    fixed = TRUE
  )
  expect_error(
    random_covariance(eta = 2, rho = rho_prior),
    "explicit 'cor' or 'rho'",
    fixed = TRUE
  )
  expect_s3_class(random_covariance(structure = "car"), "random_covariance")
  expect_error(
    random_covariance(structure = "car", cor = prior_lkj(eta = 1)),
    "uses a scalar correlation prior",
    fixed = TRUE
  )
  expect_error(
    JAGS_formula(
      formula = ~ 1 + x + car(1 | id),
      parameter = "mu",
      data = df,
      prior_list = fixed_priors,
      prior_random = prior_random(id = random_block(sd = sd_prior))
    ),
    "CAR random-effect terms require exactly one untransformed time variable"
  )
  expect_error(
    JAGS_formula(
      formula = ~ 1 + x + (1 + x | id),
      parameter = "mu",
      data = df,
      prior_list = fixed_priors,
      prior_random = prior_random(
        id = random_block(
          sd = sd_prior,
          covariance = random_covariance(structure = "diag")
        )
      )
    ),
    "The formula owns the covariance structure",
    fixed = TRUE
  )
  expect_error(
    JAGS_formula(
      formula = ~ 1 + x + (1 + x | id),
      parameter = "mu",
      data = df,
      prior_list = fixed_priors,
      prior_random = prior_random(id = random_block(sd = sd_prior))
    ),
    "requires an LKJ correlation prior",
    fixed = TRUE
  )
  expect_error(
    JAGS_formula(
      formula = ~ 1 + x + diag(1 + x | id),
      parameter = "mu",
      data = df,
      prior_list = fixed_priors,
      prior_random = prior_random(
        id = random_block(
          sd = sd_prior,
          covariance = random_covariance(sd = sd_prior)
        )
      )
    ),
    "SD prior was supplied both",
    fixed = TRUE
  )
  expect_error(
    JAGS_formula(
      formula = ~ 1 + diag(0 + x | id),
      parameter = "mu",
      data = df,
      prior_list = list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(
        id = random_block(sd = prior_factor("gamma", list(2, 2), contrast = "treatment"))
      )
    ),
    "ordinary scalar prior",
    fixed = TRUE
  )
  expect_error(
    JAGS_formula(
      formula = ~ 1 + diag(1 | id),
      parameter = "mu",
      data = df,
      prior_list = list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(
        id = random_block(sd = prior_mixture(list(prior("point", list(-1)), prior("gamma", list(2, 2)))))
      )
    ),
    "point mass must be nonnegative",
    fixed = TRUE
  )
  zero_sd <- JAGS_formula(
    formula = ~ 1 + diag(1 | id),
    parameter = "mu",
    data = df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(id = random_block(sd = prior("point", list(0))))
  )
  expect_s3_class(zero_sd$prior_list$mu__xREx__id_intercept, "prior.point")

  expect_error(
    JAGS_formula(
      formula = ~ 1 + x + diag(1 + x | id),
      parameter = "mu",
      data = df,
      prior_list = fixed_priors,
      prior_random = prior_random(id = random_block(sd = sd_prior, cor = lkj_prior))
    ),
    "structure 'diag' has no correlation parameter",
    fixed = TRUE
  )
  expect_error(
    JAGS_formula(
      formula = ~ 1 + x + (1 + x | id),
      parameter = "mu",
      data = df,
      prior_list = fixed_priors,
      prior_random = prior_random(id = random_block(sd = sd_prior, rho = rho_prior))
    ),
    "structure 'us' uses an LKJ correlation prior",
    fixed = TRUE
  )
  expect_error(
    JAGS_formula(
      formula = ~ 1 + x + ar1(idx | id),
      parameter = "mu",
      data = df,
      prior_list = fixed_priors,
      prior_random = prior_random(id = random_block(sd = sd_prior, cor = lkj_prior))
    ),
    "structure 'ar1' uses a scalar correlation prior",
    fixed = TRUE
  )
  expect_error(
    JAGS_formula(
      formula = ~ 1 + x + id(1 + x | id),
      parameter = "mu",
      data = df,
      prior_list = fixed_priors,
      prior_random = prior_random(
        covariance = random_covariance(structure = "id"),
        id = random_block(sd = sd_prior, rho = rho_prior)
      )
    ),
    "Block covariance override supplies a scalar correlation prior, but structure 'id' has no correlation parameter",
    fixed = TRUE
  )
})

test_that("variance allocation priors generate shared total SD and Dirichlet allocation syntax", {

  df <- data.frame(
    study = factor(c("s1", "s1", "s2", "s2", "s3", "s3")),
    drug = factor(c("a", "b", "a", "b", "a", "b")),
    x = c(-1, 0, 1, 2, -2, 3)
  )
  fixed_priors <- list(intercept = prior("normal", list(0, 1)))
  random_prior <- prior_random(
    allocation = random_variance_allocation(
      sd = prior("gamma", list(2, 2)),
      allocation = prior("dirichlet", list(alpha = c(2, 3)))
    )
  )

  result <- JAGS_formula(
    formula = ~ 1 +
      random(1 | study, name = "study", covariance = "diag") +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = df,
    prior_list = fixed_priors,
    prior_random = random_prior
  )

  expect_equal(
    names(result$prior_list),
    c(
      "mu_intercept",
      "mu__xRE_ALLOCx_allocation_total_sd",
      "mu__xRE_ALLOCx_allocation_weight"
    )
  )
  expect_s3_class(result$prior_list$mu__xRE_ALLOCx_allocation_weight, "prior.simplex")
  expect_match(
    result$formula_syntax,
    "mu__xREx__study_intercept = mu__xRE_ALLOCx_allocation_total_sd * sqrt(mu__xRE_ALLOCx_allocation_weight[1])",
    fixed = TRUE
  )
  expect_match(
    result$formula_syntax,
    "mu__xREx__drug_intercept = mu__xRE_ALLOCx_allocation_total_sd * sqrt(mu__xRE_ALLOCx_allocation_weight[2])",
    fixed = TRUE
  )
  expect_equal(
    result$add_parameters,
    c(
      "mu__xREx__study_intercept",
      "mu__xREx__study_xRE_Zx",
      "mu__xREx__drug_intercept",
      "mu__xREx__drug_xRE_Zx"
    )
  )
  expect_equal(
    result$formula_design$random_effects[[1]]$allocation$weight_name,
    "mu__xRE_ALLOCx_allocation_weight"
  )
  expect_equal(result$formula_design$random_effects[[1]]$allocation$index, 1L)
  expect_equal(result$formula_design$random_effects[[2]]$allocation$index, 2L)
  expect_equal(
    result$formula_design$random_effects[[1]]$sd_parameter_names,
    "mu__xREx__study_intercept"
  )

  reversed_result <- JAGS_formula(
    formula = ~ 1 +
      random(1 | study, name = "study", covariance = "diag") +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = df,
    prior_list = fixed_priors,
    prior_random = prior_random(
      allocation = random_variance_allocation(
        terms = c("drug", "study"),
        sd = prior("gamma", list(2, 2)),
        allocation = prior("dirichlet", list(alpha = c(2, 3)))
      )
    )
  )
  expect_equal(reversed_result$formula_design$random_effects[[1]]$allocation$index, 2L)
  expect_equal(reversed_result$formula_design$random_effects[[2]]$allocation$index, 1L)
  expect_match(
    reversed_result$formula_syntax,
    "mu__xREx__study_intercept = mu__xRE_ALLOCx_allocation_total_sd * sqrt(mu__xRE_ALLOCx_allocation_weight[2])",
    fixed = TRUE
  )
  expect_match(
    reversed_result$formula_syntax,
    "mu__xREx__drug_intercept = mu__xRE_ALLOCx_allocation_total_sd * sqrt(mu__xRE_ALLOCx_allocation_weight[1])",
    fixed = TRUE
  )

  full_syntax <- JAGS_add_priors("model{}", result$prior_list)
  expect_match(
    full_syntax,
    "prior_par_eta_mu__xRE_ALLOCx_allocation_weight[1] ~ dgamma(2, 1)",
    fixed = TRUE
  )
  expect_match(
    full_syntax,
    "mu__xRE_ALLOCx_allocation_weight[2] <- prior_par_eta_mu__xRE_ALLOCx_allocation_weight[2] / sum(prior_par_eta_mu__xRE_ALLOCx_allocation_weight[1:2])",
    fixed = TRUE
  )

  posterior <- matrix(
    c(
      2, 1, 3,
      4, 2, 2
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      NULL,
      c(
        "mu__xRE_ALLOCx_allocation_total_sd",
        "prior_par_eta_mu__xRE_ALLOCx_allocation_weight[1]",
        "prior_par_eta_mu__xRE_ALLOCx_allocation_weight[2]"
      )
    )
  )
  study_sd <- BayesTools:::.bt_random_effect_sd_draws(
    random_term = result$formula_design$random_effects[[1]],
    n_columns = 1,
    posterior = posterior,
    prior_list = result$prior_list
  )
  drug_sd <- BayesTools:::.bt_random_effect_sd_draws(
    random_term = result$formula_design$random_effects[[2]],
    n_columns = 1,
    posterior = posterior,
    prior_list = result$prior_list
  )
  expect_equal(study_sd[, 1], c(2 * sqrt(1 / 4), 4 * sqrt(2 / 4)))
  expect_equal(drug_sd[, 1], c(2 * sqrt(3 / 4), 4 * sqrt(2 / 4)))
  reversed_study_sd <- BayesTools:::.bt_random_effect_sd_draws(
    random_term = reversed_result$formula_design$random_effects[[1]],
    n_columns = 1,
    posterior = posterior,
    prior_list = reversed_result$prior_list
  )
  reversed_drug_sd <- BayesTools:::.bt_random_effect_sd_draws(
    random_term = reversed_result$formula_design$random_effects[[2]],
    n_columns = 1,
    posterior = posterior,
    prior_list = reversed_result$prior_list
  )
  expect_equal(reversed_study_sd[, 1], c(2 * sqrt(3 / 4), 4 * sqrt(2 / 4)))
  expect_equal(reversed_drug_sd[, 1], c(2 * sqrt(1 / 4), 4 * sqrt(2 / 4)))

  bridge_row <- posterior[1, ]
  expect_equal(
    BayesTools:::.bt_JAGS_marglik_random_effect_sd_values(
      samples = bridge_row,
      random_term = result$formula_design$random_effects[[1]],
      prior_list = result$prior_list
    ),
    2 * sqrt(1 / 4)
  )

  single_factor_df <- data.frame(
    f = factor(c("a", "b", "a", "b"), levels = c("a", "b")),
    study = factor(c("s1", "s1", "s2", "s2")),
    drug = factor(c("d1", "d2", "d1", "d2"))
  )
  single_factor_result <- JAGS_formula(
    formula = ~ 1 +
      random(0 + f | study, name = "study", covariance = "diag") +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = single_factor_df,
    prior_list = fixed_priors,
    prior_random = random_prior
  )
  study_random <- single_factor_result$formula_design$random_effects[[1]]
  expect_equal(study_random$n_columns, 1L)
  expect_equal(study_random$sd_parameter_names, "mu__xREx__study_f[1]")
  expect_equal(study_random$sd_leaves$leaf_names_by_column, study_random$sd_parameter_names)
  expect_match(
    single_factor_result$formula_syntax,
    "mu__xREx__study_f[1] = mu__xRE_ALLOCx_allocation_total_sd * sqrt(mu__xRE_ALLOCx_allocation_weight[1])",
    fixed = TRUE
  )

  random_bridge <- BayesTools:::.bt_JAGS_formula_random_bridge_parameters(
    list(mu = result$formula_design)
  )
  expect_true(all(grepl("_xRE_Zx", random_bridge$parameters, fixed = TRUE)))
  expect_false(any(grepl("_intercept", random_bridge$parameters, fixed = TRUE)))

  bridge_samples <- c(
    "mu_intercept" = 10,
    "mu__xRE_ALLOCx_allocation_total_sd" = 2,
    "prior_par_eta_mu__xRE_ALLOCx_allocation_weight[1]" = 1,
    "prior_par_eta_mu__xRE_ALLOCx_allocation_weight[2]" = 3,
    "mu__xREx__study_xRE_Zx[1,1]" = 0.1,
    "mu__xREx__study_xRE_Zx[2,1]" = 0.2,
    "mu__xREx__study_xRE_Zx[3,1]" = 0.3,
    "mu__xREx__drug_xRE_Zx[1,1]" = 1,
    "mu__xREx__drug_xRE_Zx[2,1]" = 2
  )
  fixed_formula_priors <- BayesTools:::.bt_JAGS_marglik_formula_fixed_priors(
    result$prior_list,
    "mu"
  )
  expect_equal(names(fixed_formula_priors), "mu_intercept")

  reconstructed <- JAGS_marglik_parameters_formula(
    samples = bridge_samples,
    formula_list = list(mu = result$formula),
    formula_data_list = list(mu = result$data),
    formula_prior_list = list(mu = result$prior_list),
    prior_list_parameters = list(),
    formula_design_list = list(mu = result$formula_design)
  )
  expected_mu <- 10 +
    c(0.1, 0.1, 0.2, 0.2, 0.3, 0.3) +
    sqrt(3) * c(1, 2, 1, 2, 1, 2)
  expect_equal(reconstructed$mu, expected_mu, tolerance = 1e-12)

  expect_equal(
    JAGS_marglik_priors_formula(
      samples = bridge_samples,
      formula_prior_list = list(mu = result$prior_list)
    ),
    stats::dnorm(10, 0, 1, log = TRUE) +
      stats::dgamma(2, shape = 2, rate = 2, log = TRUE) +
      stats::dgamma(1, shape = 2, rate = 1, log = TRUE) +
    stats::dgamma(3, shape = 3, rate = 1, log = TRUE),
    tolerance = 1e-12
  )
  edge_bridge_samples <- bridge_samples
  edge_bridge_samples[["prior_par_eta_mu__xRE_ALLOCx_allocation_weight[1]"]] <- 0
  expect_equal(
    JAGS_marglik_priors_formula(
      samples = edge_bridge_samples,
      formula_prior_list = list(mu = result$prior_list)
    ),
    -Inf
  )
  expect_error(
    JAGS_marglik_parameters_formula(
      samples = edge_bridge_samples,
      formula_list = list(mu = result$formula),
      formula_data_list = list(mu = result$data),
      formula_prior_list = list(mu = result$prior_list),
      prior_list_parameters = list(),
      formula_design_list = list(mu = result$formula_design)
    ),
    "Dirichlet allocation auxiliary samples must be positive"
  )
  expect_equal(
    BayesTools:::.bt_JAGS_marglik_priors_formula_random(
      bridge_samples,
      list(mu = result$formula_design)
    ),
    sum(stats::dnorm(bridge_samples[grepl("_xRE_Zx", names(bridge_samples), fixed = TRUE)], log = TRUE)),
    tolerance = 1e-12
  )
  expect_error(
    BayesTools:::.bt_JAGS_marglik_random_effect_sd_values(
      samples = bridge_samples[!grepl("prior_par_eta_", names(bridge_samples), fixed = TRUE)],
      random_term = result$formula_design$random_effects[[1]],
      prior_list = result$prior_list
    ),
    "missing Dirichlet allocation coordinates",
    fixed = TRUE
  )
  expect_error(
    BayesTools:::.bt_JAGS_marglik_random_effect_sd_values(
      samples = bridge_samples[names(bridge_samples) != "mu__xRE_ALLOCx_allocation_total_sd"],
      random_term = result$formula_design$random_effects[[1]],
      prior_list = result$prior_list
    ),
    "does not contain all monitored formula prior parameters",
    fixed = TRUE
  )

  expect_error(
    JAGS_formula(
      formula = ~ 1 +
        random(1 | study, name = "study", covariance = "diag") +
        random(1 | drug, name = "drug", covariance = "diag"),
      parameter = "mu",
      data = df,
      prior_list = fixed_priors,
      prior_random = prior_random(
        allocation = random_variance_allocation(
          terms = c("study", "missing"),
          sd = prior("gamma", list(2, 2))
        )
      )
    ),
    "unknown random-effect block",
    fixed = TRUE
  )
  expect_error(
    JAGS_formula(
      formula = ~ 1 +
        random(1 | study, name = "study", covariance = "diag") +
        random(1 | drug, name = "drug", covariance = "diag"),
      parameter = "mu",
      data = df,
      prior_list = fixed_priors,
      prior_random = prior_random(
        study = random_block(sd = prior("gamma", list(4, 1))),
        allocation = random_variance_allocation(
          terms = c("study", "drug"),
          sd = prior("gamma", list(2, 2))
        )
      )
    ),
    "cannot supply block-specific SD",
    fixed = TRUE
  )
  expect_error(
    JAGS_formula(
      formula = ~ 1 +
        random(1 + x | study, name = "study", covariance = "diag") +
        random(1 | drug, name = "drug", covariance = "diag"),
      parameter = "mu",
      data = df,
      prior_list = fixed_priors,
      prior_random = random_prior
    ),
    "one SD component",
    fixed = TRUE
  )
})

test_that("bridge auxiliary positive support is enforced before reconstruction", {

  dirichlet_prior <- prior("dirichlet", list(alpha = c(1, 2)))
  dirichlet_samples <- c(
    "prior_par_eta_w[1]" = 0,
    "prior_par_eta_w[2]" = 2
  )
  expect_equal(JAGS_marglik_priors(dirichlet_samples, list(w = dirichlet_prior)), -Inf)
  expect_error(
    JAGS_marglik_parameters(dirichlet_samples, list(w = dirichlet_prior)),
    "out-of-support positive auxiliary coordinate"
  )

  invgamma_prior <- prior("invgamma", list(1, 2))
  invgamma_samples <- c("inv_sigma" = 0)
  expect_equal(JAGS_marglik_priors(invgamma_samples, list(sigma = invgamma_prior)), -Inf)
  expect_error(
    JAGS_marglik_parameters(invgamma_samples, list(sigma = invgamma_prior)),
    "out-of-support positive auxiliary coordinate"
  )

  formula_prior_list <- list(
    mu = list(mu_intercept = invgamma_prior)
  )
  expect_equal(
    JAGS_marglik_priors_formula(
      samples = c("inv_mu_intercept" = 0),
      formula_prior_list = formula_prior_list
    ),
    -Inf
  )
  expect_error(
    JAGS_marglik_parameters_formula(
      samples = c("inv_mu_intercept" = 0),
      formula_list = list(mu = ~ 1),
      formula_data_list = list(mu = list(N_mu = 1)),
      formula_prior_list = formula_prior_list,
      prior_list_parameters = list()
    ),
    "out-of-support positive auxiliary coordinate"
  )
})

test_that("variance allocation graph supports child block and SD-leaf allocations", {

  sd_prior <- prior("gamma", list(2, 2))
  fixed_priors <- list(intercept = prior("normal", list(0, 1)))

  df_nested <- data.frame(
    study = factor(c("s1", "s1", "s2", "s2")),
    paper = factor(c("p1", "p2", "p1", "p2")),
    drug = factor(c("a", "b", "a", "b"))
  )
  nested_prior <- prior_random(
    random_variance_allocation(
      name = "total_re",
      terms = c(nested = "nested", drug = "drug"),
      sd = sd_prior,
      allocation = prior("dirichlet", list(alpha = c(1, 1)))
    ),
    random_variance_allocation(
      name = "nested_split",
      parent = allocation_ref("total_re", "nested"),
      terms = c(study = "study", paper = "paper"),
      allocation = prior("dirichlet", list(alpha = c(3, 2)))
    )
  )
  nested_result <- JAGS_formula(
    formula = ~ 1 +
      random(1 | study, name = "study", covariance = "diag") +
      random(1 | paper, name = "paper", covariance = "diag") +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = df_nested,
    prior_list = fixed_priors,
    prior_random = nested_prior
  )

  expect_equal(
    names(nested_result$prior_list),
    c(
      "mu_intercept",
      "mu__xRE_ALLOCx_total_re_total_sd",
      "mu__xRE_ALLOCx_total_re_weight",
      "mu__xRE_ALLOCx_nested_split_weight"
    )
  )
  expect_match(
    nested_result$formula_syntax,
    "mu__xRE_ALLOCx_total_re_nested_sd = mu__xRE_ALLOCx_total_re_total_sd * sqrt(mu__xRE_ALLOCx_total_re_weight[1])",
    fixed = TRUE
  )
  expect_match(
    nested_result$formula_syntax,
    "mu__xREx__study_intercept = mu__xRE_ALLOCx_total_re_nested_sd * sqrt(mu__xRE_ALLOCx_nested_split_weight[1])",
    fixed = TRUE
  )
  expect_match(
    nested_result$formula_syntax,
    "mu__xREx__drug_intercept = mu__xRE_ALLOCx_total_re_total_sd * sqrt(mu__xRE_ALLOCx_total_re_weight[2])",
    fixed = TRUE
  )
  expect_equal(length(nested_result$formula_design$random_effects[[1]]$allocation$factors), 2L)
  expect_equal(length(nested_result$formula_design$random_effects[[3]]$allocation$factors), 1L)

  nested_posterior <- matrix(
    c(4, 1, 3, 3, 2),
    nrow = 1,
    dimnames = list(
      NULL,
      c(
        "mu__xRE_ALLOCx_total_re_total_sd",
        "prior_par_eta_mu__xRE_ALLOCx_total_re_weight[1]",
        "prior_par_eta_mu__xRE_ALLOCx_total_re_weight[2]",
        "prior_par_eta_mu__xRE_ALLOCx_nested_split_weight[1]",
        "prior_par_eta_mu__xRE_ALLOCx_nested_split_weight[2]"
      )
    )
  )
  expect_equal(
    BayesTools:::.bt_random_effect_sd_draws(
      random_term = nested_result$formula_design$random_effects[[1]],
      n_columns = 1,
      posterior = nested_posterior,
      prior_list = nested_result$prior_list
    )[, 1],
    4 * sqrt(1 / 4) * sqrt(3 / 5),
    tolerance = 1e-12
  )
  expect_equal(
    BayesTools:::.bt_random_effect_sd_draws(
      random_term = nested_result$formula_design$random_effects[[3]],
      n_columns = 1,
      posterior = nested_posterior,
      prior_list = nested_result$prior_list
    )[, 1],
    4 * sqrt(3 / 4),
    tolerance = 1e-12
  )

  df_sd <- data.frame(
    study = factor(c("s1", "s1", "s2", "s2")),
    drug = factor(c("a", "b", "a", "b")),
    x = c(-1, 0, 1, 2)
  )
  sd_leaf_prior <- prior_random(
    random_variance_allocation(
      name = "total_re",
      terms = c(het = "het", simple = "drug"),
      sd = sd_prior,
      allocation = prior("dirichlet", list(alpha = c(1, 1)))
    ),
    random_variance_allocation(
      name = "het_sd",
      parent = allocation_ref("total_re", "het"),
      terms = "het",
      components = "sd",
      scale = "mean_variance",
      allocation = prior("dirichlet", list(alpha = c(2, 2)))
    )
  )
  sd_leaf_result <- JAGS_formula(
    formula = ~ 1 +
      random(1 + x | study, name = "het", covariance = "diag") +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = df_sd,
    prior_list = fixed_priors,
    prior_random = sd_leaf_prior
  )

  expect_match(
    sd_leaf_result$formula_syntax,
    "mu__xREx__het_intercept = mu__xRE_ALLOCx_total_re_het_sd * sqrt(2 * mu__xRE_ALLOCx_het_sd_weight[1])",
    fixed = TRUE
  )
  expect_match(
    sd_leaf_result$formula_syntax,
    "mu__xREx__het_x = mu__xRE_ALLOCx_total_re_het_sd * sqrt(2 * mu__xRE_ALLOCx_het_sd_weight[2])",
    fixed = TRUE
  )
  expect_equal(
    sd_leaf_result$formula_design$random_effects[[1]]$allocation$components,
    "sd"
  )
  expect_equal(
    sd_leaf_result$formula_design$random_effects[[1]]$allocation$leaf_terms,
    c(mu__xREx__het_intercept = "intercept", mu__xREx__het_x = "x")
  )
  expect_equal(
    sd_leaf_result$formula_design$random_effects[[1]]$sd_leaves$leaf_names,
    names(sd_leaf_result$formula_design$random_effects[[1]]$allocation$leaf_terms)
  )
  expect_equal(
    sd_leaf_result$formula_design$random_effects[[1]]$sd_leaves$leaf_index_by_column,
    c(1L, 2L)
  )

  sd_leaf_samples <- c(
    "mu__xRE_ALLOCx_total_re_total_sd" = 2,
    "prior_par_eta_mu__xRE_ALLOCx_total_re_weight[1]" = 1,
    "prior_par_eta_mu__xRE_ALLOCx_total_re_weight[2]" = 3,
    "prior_par_eta_mu__xRE_ALLOCx_het_sd_weight[1]" = 1,
    "prior_par_eta_mu__xRE_ALLOCx_het_sd_weight[2]" = 3
  )
  expect_equal(
    BayesTools:::.bt_JAGS_marglik_random_effect_sd_values(
      samples = sd_leaf_samples,
      random_term = sd_leaf_result$formula_design$random_effects[[1]],
      prior_list = sd_leaf_result$prior_list
    ),
    c(
      2 * sqrt(1 / 4) * sqrt(2 * 1 / 4),
      2 * sqrt(1 / 4) * sqrt(2 * 3 / 4)
    ),
    tolerance = 1e-12
  )
})

test_that("variance allocation graph rejects ambiguous or conflicting specifications", {

  sd_prior <- prior("gamma", list(2, 2))
  fixed_priors <- list(intercept = prior("normal", list(0, 1)))
  df <- data.frame(
    study = factor(c("s1", "s1", "s2", "s2")),
    drug = factor(c("a", "b", "a", "b")),
    x = c(-1, 0, 1, 2)
  )

  expect_error(
    random_variance_allocation(
      name = "bad_child",
      parent = allocation_ref("total_re", "study"),
      terms = "study",
      sd = sd_prior,
      allocation = prior("dirichlet", list(alpha = c(1, 1)))
    ),
    "must not specify 'sd'",
    fixed = TRUE
  )
  expect_error(
    random_variance_allocation(
      name = "bad_mean",
      terms = c("study", "drug"),
      sd = sd_prior,
      scale = "mean_variance"
    ),
    "components = \"sd\"",
    fixed = TRUE
  )
  expect_error(
    random_variance_allocation(
      name = "bad_sd_terms",
      terms = c("study", "drug"),
      sd = sd_prior,
      components = "sd"
    ),
    "exactly one",
    fixed = TRUE
  )
  expect_error(
    random_variance_allocation(
      terms = c(study = "study", "drug"),
      sd = sd_prior,
      allocation = prior("dirichlet", list(alpha = c(1, 1)))
    ),
    "either all named or all unnamed",
    fixed = TRUE
  )
  expect_error(
    JAGS_formula(
      formula = ~ 1 +
        random(1 | study, name = "study", covariance = "diag") +
        random(1 | drug, name = "drug", covariance = "diag"),
      parameter = "mu",
      data = df,
      prior_list = fixed_priors,
      prior_random = prior_random(
        random_variance_allocation(
          sd = sd_prior,
          allocation = prior("dirichlet", list(alpha = c(1, 1)))
        ),
        random_variance_allocation(
          name = "second",
          terms = c("study", "drug"),
          sd = sd_prior,
          allocation = prior("dirichlet", list(alpha = c(1, 1)))
        )
      )
    ),
    "without explicit 'terms'",
    fixed = TRUE
  )
  expect_error(
    JAGS_formula(
      formula = ~ 1 +
        random(1 | study, name = "study", covariance = "diag") +
        random(1 | drug, name = "drug", covariance = "diag"),
      parameter = "mu",
      data = df,
      prior_list = fixed_priors,
      prior_random = prior_random(
        random_variance_allocation(
          name = "duplicate",
          terms = c("study", "drug"),
          sd = sd_prior,
          allocation = prior("dirichlet", list(alpha = c(1, 1)))
        ),
        random_variance_allocation(
          name = "duplicate",
          terms = c("study", "drug"),
          sd = sd_prior,
          allocation = prior("dirichlet", list(alpha = c(1, 1)))
        )
      )
    ),
    "labels must be unique",
    fixed = TRUE
  )
  expect_error(
    JAGS_formula(
      formula = ~ 1 +
        random(1 | study, name = "study", covariance = "diag") +
        random(1 | drug, name = "drug", covariance = "diag"),
      parameter = "mu",
      data = df,
      prior_list = fixed_priors,
      prior_random = prior_random(
        random_variance_allocation(
          name = "total_re",
          terms = c(symbolic = "symbolic", drug = "drug"),
          sd = sd_prior,
          allocation = prior("dirichlet", list(alpha = c(1, 1)))
        )
      )
    ),
    "Unknown targets are allowed only when consumed by a child allocation",
    fixed = TRUE
  )
  expect_error(
    JAGS_formula(
      formula = ~ 1 +
        random(1 | study, name = "study", covariance = "diag") +
        random(1 | drug, name = "drug", covariance = "diag"),
      parameter = "mu",
      data = df,
      prior_list = fixed_priors,
      prior_random = prior_random(
        random_variance_allocation(
          name = "child_first",
          parent = allocation_ref("total_re", "study"),
          terms = "study",
          components = "sd",
          scale = "mean_variance",
          allocation = prior("dirichlet", list(alpha = c(1, 1)))
        ),
        random_variance_allocation(
          name = "total_re",
          terms = c("study", "drug"),
          sd = sd_prior,
          allocation = prior("dirichlet", list(alpha = c(1, 1)))
        )
      )
    ),
    "must be defined before child allocation",
    fixed = TRUE
  )
  expect_error(
    JAGS_formula(
      formula = ~ 1 +
        random(1 | study, name = "study", covariance = "diag") +
        random(1 | drug, name = "drug", covariance = "diag"),
      parameter = "mu",
      data = df,
      prior_list = fixed_priors,
      prior_random = prior_random(
        random_variance_allocation(
          name = "total_re",
          terms = c("study", "drug"),
          sd = sd_prior,
          allocation = prior("dirichlet", list(alpha = c(1, 1)))
        ),
        random_variance_allocation(
          name = "bad_ref",
          parent = allocation_ref("total_re", "missing"),
          terms = "study",
          components = "sd",
          scale = "mean_variance",
          allocation = prior("dirichlet", list(alpha = c(1, 1)))
        )
      )
    ),
    "does not contain component 'missing'",
    fixed = TRUE
  )
  expect_error(
    JAGS_formula(
      formula = ~ 1 +
        random(1 | study, name = "study", covariance = "diag") +
        random(1 | drug, name = "drug", covariance = "diag"),
      parameter = "mu",
      data = df,
      prior_list = fixed_priors,
      prior_random = prior_random(
        random_variance_allocation(
          name = "total_re",
          terms = c(study = "study", drug = "drug"),
          sd = sd_prior,
          allocation = prior("dirichlet", list(alpha = c(1, 1)))
        ),
        random_variance_allocation(
          name = "study_split_1",
          parent = allocation_ref("total_re", "study"),
          terms = "study",
          components = "sd",
          allocation = prior("dirichlet", list(alpha = c(1, 1)))
        ),
        random_variance_allocation(
          name = "study_split_2",
          parent = allocation_ref("total_re", "study"),
          terms = "drug",
          components = "sd",
          allocation = prior("dirichlet", list(alpha = c(1, 1)))
        )
      )
    ),
    "can be consumed by only one child",
    fixed = TRUE
  )
  expect_error(
    JAGS_formula(
      formula = ~ 1 +
        random(1 | study, name = "study", covariance = "diag") +
        random(1 | drug, name = "drug", covariance = "diag"),
      parameter = "mu",
      data = df,
      prior_list = fixed_priors,
      prior_random = prior_random(
        random_variance_allocation(
          name = "total_re",
          terms = c(study = "study", drug = "drug"),
          sd = sd_prior,
          allocation = prior("dirichlet", list(alpha = c(1, 1)))
        ),
        random_variance_allocation(
          name = "study_split",
          parent = allocation_ref("total_re", "study"),
          terms = c("study", "drug"),
          allocation = prior("dirichlet", list(alpha = c(1, 1)))
        )
      )
    ),
    "also names a random-effect block",
    fixed = TRUE
  )
  expect_error(
    BayesTools:::.bt_random_variance_allocation_component_labels(c("a-b", "a_b")),
    "unique after sanitization",
    fixed = TRUE
  )
  expect_error(
    JAGS_formula(
      formula = ~ 1 + random(1 + x | study, name = "study", covariance = "diag"),
      parameter = "mu",
      data = df,
      prior_list = fixed_priors,
      prior_random = prior_random(
        random_variance_allocation(
          name = "study_sd",
          terms = "study",
          sd = sd_prior,
          components = "sd",
          allocation = prior("dirichlet", list(alpha = c(1, 1, 1)))
        )
      )
    ),
    "dimension must match",
    fixed = TRUE
  )
  expect_error(
    JAGS_formula(
      formula = ~ 1 + random(1 | study, name = "study", covariance = "diag"),
      parameter = "mu",
      data = df,
      prior_list = fixed_priors,
      prior_random = prior_random(
        random_variance_allocation(
          name = "study_sd",
          terms = "study",
          sd = sd_prior,
          components = "sd",
          allocation = prior("dirichlet", list(alpha = c(1, 1)))
        )
      )
    ),
    "at least two resolved SD components",
    fixed = TRUE
  )
})

test_that("random-effect summary samples expose semantic SD, rho, and allocation filters", {

  sd_prior <- prior("gamma", list(2, 2))
  fixed_priors <- list(intercept = prior("normal", list(0, 1)))

  df_nested <- data.frame(
    study = factor(c("s1", "s1", "s2", "s2")),
    paper = factor(c("p1", "p2", "p1", "p2")),
    drug = factor(c("a", "b", "a", "b"))
  )
  nested_result <- JAGS_formula(
    formula = ~ 1 +
      random(1 | study, name = "study", covariance = "diag") +
      random(1 | paper, name = "paper", covariance = "diag") +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = df_nested,
    prior_list = fixed_priors,
    prior_random = prior_random(
      random_variance_allocation(
        name = "total_re",
        terms = c(nested = "nested", drug = "drug"),
        sd = sd_prior,
        allocation = prior("dirichlet", list(alpha = c(1, 1)))
      ),
      random_variance_allocation(
        name = "nested_split",
        parent = allocation_ref("total_re", "nested"),
        terms = c(study = "study", paper = "paper"),
        allocation = prior("dirichlet", list(alpha = c(3, 2)))
      )
    )
  )
  nested_posterior <- matrix(
    c(
      4, 1, 3, 3, 2,
      5, 2, 2, 1, 3
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      NULL,
      c(
        "mu__xRE_ALLOCx_total_re_total_sd",
        "prior_par_eta_mu__xRE_ALLOCx_total_re_weight[1]",
        "prior_par_eta_mu__xRE_ALLOCx_total_re_weight[2]",
        "prior_par_eta_mu__xRE_ALLOCx_nested_split_weight[1]",
        "prior_par_eta_mu__xRE_ALLOCx_nested_split_weight[2]"
      )
    )
  )

  nested_summary <- BayesTools:::.bt_random_effect_summary_samples(
    model_samples = nested_posterior,
    prior_list = nested_result$prior_list,
    formula_design = list(mu = nested_result$formula_design),
    mode = "standard"
  )
  expect_false(any(grepl("__xRE_ALLOCx", colnames(nested_summary$model_samples), fixed = TRUE)))
  expect_true(all(c(
    "mu__xRE_SUMMARY__sd_total__total_re",
    "mu__xRE_SUMMARY__var_frac__total_re__nested",
    "mu__xRE_SUMMARY__var_frac__nested_split__study",
    "mu__xRE_SUMMARY__sd__study__intercept",
    "mu__xRE_SUMMARY__sd__drug__intercept"
  ) %in% colnames(nested_summary$model_samples)))
  expect_equal(
    nested_summary$model_samples[, "mu__xRE_SUMMARY__var_frac__total_re__nested"],
    c(1 / 4, 2 / 4),
    tolerance = 1e-12
  )
  expect_equal(
    unname(nested_summary$model_samples[1, "mu__xRE_SUMMARY__sd__study__intercept"]),
    4 * sqrt(1 / 4) * sqrt(3 / 5),
    tolerance = 1e-12
  )

  display_names <- BayesTools:::.bt_random_effect_summary_display_names(
    names = colnames(nested_summary$model_samples),
    raw_names = colnames(nested_summary$model_samples),
    prior_list = nested_summary$prior_list,
    formula_prefix = TRUE
  )
  expect_true("(mu) sd_total(total_re)" %in% display_names)
  expect_true("(mu) var_frac(total_re: nested)" %in% display_names)
  expect_true("(mu) sd(intercept | study)" %in% display_names)

  missing_nested_posterior <- nested_posterior[
    ,
    colnames(nested_posterior) != "prior_par_eta_mu__xRE_ALLOCx_nested_split_weight[2]",
    drop = FALSE
  ]
  expect_error(
    BayesTools:::.bt_random_effect_summary_samples(
      model_samples = missing_nested_posterior,
      prior_list = nested_result$prior_list,
      formula_design = list(mu = nested_result$formula_design),
      mode = "standard"
    ),
    "missing canonical SD coordinates"
  )
  edge_nested_posterior <- nested_posterior
  edge_nested_posterior[1, "prior_par_eta_mu__xRE_ALLOCx_total_re_weight[1]"] <- 0
  expect_error(
    BayesTools:::.bt_random_effect_summary_samples(
      model_samples = edge_nested_posterior,
      prior_list = nested_result$prior_list,
      formula_design = list(mu = nested_result$formula_design),
      mode = "standard"
    ),
    "Dirichlet allocation auxiliary samples must be positive"
  )

  direct_sd_names <- unique(unlist(lapply(
    nested_result$formula_design$random_effects,
    function(random_term) random_term$sd_parameter_names
  )))
  allocation_missing_weight_posterior <- matrix(
    c(4, 1, 2, 3),
    nrow = 1,
    dimnames = list(
      NULL,
      c("mu__xRE_ALLOCx_total_re_total_sd", direct_sd_names)
    )
  )
  expect_error(
    BayesTools:::.bt_random_effect_summary_samples(
      model_samples = allocation_missing_weight_posterior,
      prior_list = nested_result$prior_list,
      formula_design = list(mu = nested_result$formula_design),
      mode = "standard"
    ),
    "missing Dirichlet allocation coordinates",
    fixed = TRUE
  )

  df_slash <- data.frame(
    batch = factor(c("b1", "b1", "b2", "b2")),
    cask = factor(c("c1", "c2", "c1", "c2"))
  )
  slash_result <- JAGS_formula(
    formula = ~ 1 + (1 | batch/cask),
    parameter = "mu",
    data = df_slash,
    prior_list = fixed_priors,
    prior_random = prior_random(sd = sd_prior)
  )
  slash_summary <- BayesTools:::.bt_random_effect_summary_samples(
    model_samples = matrix(
      c(1, 2, 3, 4),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(
        NULL,
        c("mu__xREx__cask_batch_intercept", "mu__xREx__batch_intercept")
      )
    ),
    prior_list = slash_result$prior_list,
    formula_design = list(mu = slash_result$formula_design),
    mode = "standard"
  )
  expect_true("mu__xRE_SUMMARY__sd__cask_batch__intercept" %in% colnames(slash_summary$model_samples))
  slash_display_names <- BayesTools:::.bt_random_effect_summary_display_names(
    names = colnames(slash_summary$model_samples),
    raw_names = colnames(slash_summary$model_samples),
    prior_list = slash_summary$prior_list,
    formula_prefix = TRUE
  )
  expect_true("(mu) sd(intercept | cask:batch)" %in% slash_display_names)
  expect_false("(mu) sd(intercept | cask_batch)" %in% slash_display_names)

  slash_raw_summary <- BayesTools:::.bt_random_effect_summary_samples(
    model_samples = matrix(
      c(1, 2, 3, 4),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(
        NULL,
        c("mu__xREx__cask_batch_intercept", "mu__xREx__batch_intercept")
      )
    ),
    prior_list = slash_result$prior_list,
    formula_design = list(mu = slash_result$formula_design),
    mode = "raw"
  )
  slash_raw_samples <- BayesTools:::.rename_factor_levels(
    slash_raw_summary$model_samples,
    slash_raw_summary$prior_list
  )
  slash_raw_names <- colnames(slash_raw_samples)
  slash_raw_display_names <- BayesTools:::.bt_random_effect_summary_display_names(
    names = format_parameter_names(
      parameters = slash_raw_names,
      formula_parameters = unique(unlist(lapply(slash_raw_summary$prior_list, attr, which = "parameter"))),
      formula_random = unique(unlist(lapply(slash_raw_summary$prior_list, attr, which = "random_factor"))),
      formula_prefix = TRUE
    ),
    raw_names = slash_raw_names,
    prior_list = slash_raw_summary$prior_list,
    formula_prefix = TRUE,
    formula_design = list(mu = slash_result$formula_design)
  )
  expect_true("(mu) sd(intercept | cask:batch)" %in% slash_raw_display_names)
  expect_false(any(grepl("cask_batch", slash_raw_display_names, fixed = TRUE)))
  expect_false(any(grepl("sd((mu)", slash_raw_display_names, fixed = TRUE)))

  slash_metadata <- BayesTools:::.bt_random_effect_summary_metadata_table(
    parameter_names = colnames(slash_summary$model_samples),
    prior_list = slash_summary$prior_list,
    formula_design = list(mu = slash_result$formula_design)
  )
  slash_cask_row <- match(
    "mu__xRE_SUMMARY__sd__cask_batch__intercept",
    colnames(slash_summary$model_samples)
  )
  expect_equal(slash_metadata[slash_cask_row, "Random name"], "cask:batch")
  expect_equal(slash_metadata[slash_cask_row, "Random grouping"], "cask:batch")
  expect_equal(
    slash_metadata[slash_cask_row, "Random structure"],
    slash_result$formula_design$random_effects[[1]]$structure
  )

  named_result <- JAGS_formula(
    formula = ~ 1 + diag(1 | study, name = "study_re"),
    parameter = "mu",
    data = df_nested,
    prior_list = fixed_priors,
    prior_random = prior_random(study_re = random_block(sd = sd_prior))
  )
  named_summary <- BayesTools:::.bt_random_effect_summary_samples(
    model_samples = matrix(
      c(1, 2),
      nrow = 2,
      dimnames = list(NULL, "mu__xREx__study_re_intercept")
    ),
    prior_list = named_result$prior_list,
    formula_design = list(mu = named_result$formula_design),
    mode = "standard"
  )
  named_metadata <- BayesTools:::.bt_random_effect_summary_metadata_table(
    parameter_names = colnames(named_summary$model_samples),
    prior_list = named_summary$prior_list,
    formula_design = list(mu = named_result$formula_design)
  )
  expect_equal(named_metadata[1, "Random name"], "study_re")
  expect_equal(named_metadata[1, "Random grouping"], "study")
  expect_equal(named_metadata[1, "Random structure"], "diag")
  expect_error(
    BayesTools:::.bt_random_effect_summary_samples(
      model_samples = matrix(
        1,
        nrow = 1,
        dimnames = list(NULL, "unrelated")
      ),
      prior_list = named_result$prior_list,
      formula_design = list(mu = named_result$formula_design),
      mode = "standard"
    ),
    "missing canonical SD coordinates"
  )

  malformed_named <- named_result$formula_design
  malformed_named$random_effects[[1]]$structure <- NULL
  malformed_named$random_effects[[1]]$covariance <- "diag"
  expect_error(
    BayesTools:::.bt_random_effect_summary_samples(
      model_samples = matrix(
        c(1, 2),
        nrow = 2,
        dimnames = list(NULL, "mu__xREx__study_re_intercept")
      ),
      prior_list = named_result$prior_list,
      formula_design = list(mu = malformed_named),
      mode = "standard"
    ),
    "missing canonical 'random_term\\$structure'"
  )
  expect_error(
    BayesTools:::.bt_random_effect_summary_metadata_table(
      parameter_names = "mu__xREx__study_re_xRE_Zx[1,1]",
      prior_list = named_result$prior_list,
      formula_design = list(mu = malformed_named)
    ),
    "missing canonical 'random_term\\$structure'"
  )

  mixed_prior_list <- c(list(intercept = fixed_priors$intercept), named_summary$prior_list)
  remove_for_named <- BayesTools:::.filter_parameters(
    mixed_prior_list,
    keep_random_effects = "study_re",
    remove_spike_0 = FALSE
  )
  kept_named <- setdiff(names(mixed_prior_list), remove_for_named)
  expect_true("intercept" %in% kept_named)
  expect_true("mu_intercept" %in% kept_named)
  expect_true("mu__xRE_SUMMARY__sd__study_re__intercept" %in% kept_named)

  remove_for_other_name <- BayesTools:::.filter_parameters(
    mixed_prior_list,
    keep_random_effects = "other_re",
    remove_spike_0 = FALSE
  )
  kept_other_name <- setdiff(names(mixed_prior_list), remove_for_other_name)
  expect_equal(kept_other_name, c("intercept", "mu_intercept"))

  remove_for_structure <- BayesTools:::.filter_parameters(
    mixed_prior_list,
    keep_random_structures = "diag",
    remove_spike_0 = FALSE
  )
  kept_structure <- setdiff(names(mixed_prior_list), remove_for_structure)
  expect_true("intercept" %in% kept_structure)
  expect_true("mu_intercept" %in% kept_structure)
  expect_true("mu__xRE_SUMMARY__sd__study_re__intercept" %in% kept_structure)

  remove_diag_structure <- BayesTools:::.filter_parameters(
    mixed_prior_list,
    remove_random_structures = "diag",
    remove_spike_0 = FALSE
  )
  kept_after_remove_diag <- setdiff(names(mixed_prior_list), remove_diag_structure)
  expect_equal(kept_after_remove_diag, c("intercept", "mu_intercept"))

  remove_for_sd <- BayesTools:::.filter_parameters(
    nested_summary$prior_list,
    keep_parameters = "random_sd",
    remove_spike_0 = FALSE
  )
  kept_sd <- setdiff(names(nested_summary$prior_list), remove_for_sd)
  expect_true("mu__xRE_SUMMARY__sd__study__intercept" %in% kept_sd)
  expect_true("mu__xRE_SUMMARY__sd_total__total_re" %in% kept_sd)
  expect_false(any(grepl("__var_frac__", kept_sd, fixed = TRUE)))

  remove_for_fraction <- BayesTools:::.filter_parameters(
    nested_summary$prior_list,
    keep_parameters = "random_variance_fraction",
    remove_spike_0 = FALSE
  )
  kept_fraction <- setdiff(names(nested_summary$prior_list), remove_for_fraction)
  expect_true(all(grepl("__var_frac__", kept_fraction, fixed = TRUE)))

  df_consumed <- data.frame(
    study = factor(c("s1", "s1", "s2", "s2")),
    paper = factor(c("p1", "p2", "p1", "p2")),
    country = factor(c("c1", "c1", "c2", "c2")),
    source = factor(c("a", "b", "a", "b"))
  )
  consumed_result <- JAGS_formula(
    formula = ~ 1 +
      random(1 | study, name = "study", covariance = "diag") +
      random(1 | paper, name = "paper", covariance = "diag") +
      random(1 | country, name = "country", covariance = "diag") +
      random(1 | source, name = "source", covariance = "diag"),
    parameter = "mu",
    data = df_consumed,
    prior_list = fixed_priors,
    prior_random = prior_random(
      random_variance_allocation(
        name = "total_re",
        terms = c(site = "site", lab = "lab"),
        sd = sd_prior,
        allocation = prior("dirichlet", list(alpha = c(1, 1)))
      ),
      random_variance_allocation(
        name = "site_split",
        parent = allocation_ref("total_re", "site"),
        terms = c(study = "study", paper = "paper"),
        allocation = prior("dirichlet", list(alpha = c(1, 1)))
      ),
      random_variance_allocation(
        name = "lab_split",
        parent = allocation_ref("total_re", "lab"),
        terms = c(country = "country", source = "source"),
        allocation = prior("dirichlet", list(alpha = c(1, 1)))
      )
    )
  )
  consumed_summary <- BayesTools:::.bt_random_effect_summary_samples(
    model_samples = matrix(
      c(2, 1, 3, 2, 2, 3, 1),
      nrow = 1,
      dimnames = list(
        NULL,
        c(
          "mu__xRE_ALLOCx_total_re_total_sd",
          "prior_par_eta_mu__xRE_ALLOCx_total_re_weight[1]",
          "prior_par_eta_mu__xRE_ALLOCx_total_re_weight[2]",
          "prior_par_eta_mu__xRE_ALLOCx_site_split_weight[1]",
          "prior_par_eta_mu__xRE_ALLOCx_site_split_weight[2]",
          "prior_par_eta_mu__xRE_ALLOCx_lab_split_weight[1]",
          "prior_par_eta_mu__xRE_ALLOCx_lab_split_weight[2]"
        )
      )
    ),
    prior_list = consumed_result$prior_list,
    formula_design = list(mu = consumed_result$formula_design),
    mode = "standard"
  )
  expect_true(all(c(
    "mu__xRE_SUMMARY__var_frac__total_re__site",
    "mu__xRE_SUMMARY__var_frac__total_re__lab",
    "mu__xRE_SUMMARY__var_frac__site_split__study",
    "mu__xRE_SUMMARY__var_frac__lab_split__source"
  ) %in% colnames(consumed_summary$model_samples)))
  expect_equal(
    unname(consumed_summary$model_samples[1, "mu__xRE_SUMMARY__var_frac__total_re__site"]),
    1 / 4,
    tolerance = 1e-12
  )

  df_sd <- data.frame(
    study = factor(c("s1", "s1", "s2", "s2")),
    drug = factor(c("a", "b", "a", "b")),
    x = c(-1, 0, 1, 2)
  )
  sd_leaf_result <- JAGS_formula(
    formula = ~ 1 +
      random(1 + x | study, name = "het", covariance = "diag") +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = df_sd,
    prior_list = fixed_priors,
    prior_random = prior_random(
      random_variance_allocation(
        name = "total_re",
        terms = c(het = "het", simple = "drug"),
        sd = sd_prior,
        allocation = prior("dirichlet", list(alpha = c(1, 1)))
      ),
      random_variance_allocation(
        name = "het_sd",
        parent = allocation_ref("total_re", "het"),
        terms = "het",
        components = "sd",
        scale = "mean_variance",
        allocation = prior("dirichlet", list(alpha = c(2, 2)))
      )
    )
  )
  sd_leaf_summary <- BayesTools:::.bt_random_effect_summary_samples(
    model_samples = matrix(
      c(2, 1, 3, 1, 3),
      nrow = 1,
      dimnames = list(
        NULL,
        c(
          "mu__xRE_ALLOCx_total_re_total_sd",
          "prior_par_eta_mu__xRE_ALLOCx_total_re_weight[1]",
          "prior_par_eta_mu__xRE_ALLOCx_total_re_weight[2]",
          "prior_par_eta_mu__xRE_ALLOCx_het_sd_weight[1]",
          "prior_par_eta_mu__xRE_ALLOCx_het_sd_weight[2]"
        )
      )
    ),
    prior_list = sd_leaf_result$prior_list,
    formula_design = list(mu = sd_leaf_result$formula_design),
    mode = "full"
  )
  expect_equal(
    unname(sd_leaf_summary$model_samples[1, "mu__xRE_SUMMARY__sd_mult__het_sd__x"]),
    sqrt(2 * 3 / 4),
    tolerance = 1e-12
  )

  homogeneous_result <- JAGS_formula(
    formula = ~ 1 + id(1 + x | study),
    parameter = "mu",
    data = df_sd,
    prior_list = fixed_priors,
    prior_random = prior_random(
      study = random_block(sd = sd_prior)
    )
  )
  homogeneous_summary <- BayesTools:::.bt_random_effect_summary_samples(
    model_samples = matrix(
      2,
      nrow = 1,
      dimnames = list(NULL, "mu__xREx__study_sd")
    ),
    prior_list = homogeneous_result$prior_list,
    formula_design = list(mu = homogeneous_result$formula_design),
    mode = "standard"
  )
  expect_true("mu__xRE_SUMMARY__sd__study__shared" %in%
    colnames(homogeneous_summary$model_samples))
  expect_false("mu__xRE_SUMMARY__sd__study__sd" %in%
    colnames(homogeneous_summary$model_samples))

  scaled_allocation_result <- JAGS_formula(
    formula = ~ 1 +
      random(0 + x | study, name = "slope", covariance = "diag") +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = df_sd,
    prior_list = fixed_priors,
    formula_scale = list(x = TRUE),
    prior_random = prior_random(
      random_variance_allocation(
        name = "total_re",
        terms = c(slope = "slope", drug = "drug"),
        sd = sd_prior,
        allocation = prior("dirichlet", list(alpha = c(1, 1)))
      )
    )
  )
  scaled_allocation_summary <- BayesTools:::.bt_random_effect_summary_samples(
    model_samples = matrix(
      c(3, 1, 3),
      nrow = 1,
      dimnames = list(
        NULL,
        c(
          "mu__xRE_ALLOCx_total_re_total_sd",
          "prior_par_eta_mu__xRE_ALLOCx_total_re_weight[1]",
          "prior_par_eta_mu__xRE_ALLOCx_total_re_weight[2]"
        )
      )
    ),
    prior_list = scaled_allocation_result$prior_list,
    formula_design = list(mu = scaled_allocation_result$formula_design),
    mode = "standard",
    formula_scale = list(mu = scaled_allocation_result$formula_scale)
  )
  expect_equal(
    unname(scaled_allocation_summary$model_samples[1, "mu__xRE_SUMMARY__sd__slope__x"]),
    3 * sqrt(1 / 4) / scaled_allocation_result$formula_scale$mu_x$sd,
    tolerance = 1e-12
  )
  allocation_internal_samples <- matrix(
    c(3, 0.4, 0.6, 1.5),
    nrow = 1,
    dimnames = list(
      NULL,
      c(
        "mu__xRE_ALLOCx_total_re_total_sd",
        "mu__xRE_ALLOCx_total_re_weight[1]",
        "mu__xRE_ALLOCx_total_re_weight[2]",
        "mu__xRE_SUMMARY__sd__slope__x"
      )
    )
  )
  allocation_internal_warnings <- character()
  transformed_allocation_internal <- withCallingHandlers(
    transform_scale_samples(
      allocation_internal_samples,
      formula_scale = list(mu = scaled_allocation_result$formula_scale)
    ),
    warning = function(w) {
      allocation_internal_warnings <<- c(
        allocation_internal_warnings,
        conditionMessage(w)
      )
      invokeRestart("muffleWarning")
    }
  )
  expect_equal(allocation_internal_warnings, character())
  expect_equal(transformed_allocation_internal, allocation_internal_samples)

  df_rho <- data.frame(
    id = factor(c("a", "a", "b", "b")),
    idx = factor(c("t1", "t2", "t1", "t2"), levels = c("t1", "t2")),
    x = c(-1, 0, 1, 2)
  )
  ar1_result <- JAGS_formula(
    formula = ~ 1 + ar1(idx | id),
    parameter = "mu",
    data = df_rho,
    prior_list = fixed_priors,
    prior_random = prior_random(
      id = random_block(
        sd = sd_prior,
        rho = prior("normal", list(0, 0.5))
      )
    )
  )
  rho_posterior <- matrix(
    c(2, 0.25, 99, 3, -0.5, 88),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      NULL,
      c("mu__xREx__id_sd", "mu__xREx__id_rho_z", "mu__xREx__id_xRE_Zx[1,1]")
    )
  )
  rho_summary <- BayesTools:::.bt_random_effect_summary_samples(
    model_samples = rho_posterior,
    prior_list = ar1_result$prior_list,
    formula_design = list(mu = ar1_result$formula_design),
    mode = "standard"
  )
  expect_false(any(grepl("_xRE_Zx", colnames(rho_summary$model_samples), fixed = TRUE)))
  expect_equal(
    rho_summary$model_samples[, "mu__xRE_SUMMARY__rho__id"],
    tanh(c(0.25, -0.5)),
    tolerance = 1e-12
  )
  remove_for_rho <- BayesTools:::.filter_parameters(
    rho_summary$prior_list,
    keep_parameters = "random_rho",
    remove_spike_0 = FALSE
  )
  kept_rho <- setdiff(names(rho_summary$prior_list), remove_for_rho)
  expect_equal(kept_rho, "mu__xRE_SUMMARY__rho__id")

  us_result <- JAGS_formula(
    formula = ~ 1 + us(1 + x | id),
    parameter = "mu",
    data = df_rho,
    prior_list = fixed_priors,
    prior_random = prior_random(
      id = random_block(
        sd = sd_prior,
        cor = prior_lkj(eta = 1)
      )
    )
  )
  us_random_term <- us_result$formula_design$random_effects[[1]]
  L_names <- BayesTools:::.bt_random_effect_cholesky_names(us_random_term, 2L)
  us_posterior <- matrix(
    c(2, 3, 1, 0.4, 0, sqrt(1 - 0.4^2)),
    nrow = 1,
    dimnames = list(NULL, c(
      "mu__xREx__id_intercept",
      "mu__xREx__id_x",
      as.vector(L_names)
    ))
  )
  us_summary <- BayesTools:::.bt_random_effect_summary_samples(
    model_samples = us_posterior,
    prior_list = us_result$prior_list,
    formula_design = list(mu = us_result$formula_design),
    mode = "standard"
  )
  expect_equal(
    unname(us_summary$model_samples[1, "mu__xRE_SUMMARY__cor__id__intercept__x"]),
    0.4,
    tolerance = 1e-12
  )

  u_names <- BayesTools:::.bt_random_effect_lkj_primitive_names(us_random_term, 2L)
  us_primitive_posterior <- matrix(
    c(2, 3, 0.7),
    nrow = 1,
    dimnames = list(NULL, c(
      "mu__xREx__id_intercept",
      "mu__xREx__id_x",
      u_names
    ))
  )
  us_primitive_summary <- BayesTools:::.bt_random_effect_summary_samples(
    model_samples = us_primitive_posterior,
    prior_list = us_result$prior_list,
    formula_design = list(mu = us_result$formula_design),
    mode = "standard"
  )
  expect_equal(
    unname(us_primitive_summary$model_samples[1, "mu__xRE_SUMMARY__cor__id__intercept__x"]),
    0.4,
    tolerance = 1e-12
  )

  us_missing_correlation <- us_primitive_posterior[
    ,
    c("mu__xREx__id_intercept", "mu__xREx__id_x"),
    drop = FALSE
  ]
  expect_error(
    BayesTools:::.bt_random_effect_summary_samples(
      model_samples = us_missing_correlation,
      prior_list = us_result$prior_list,
      formula_design = list(mu = us_result$formula_design),
      mode = "standard"
    ),
    "missing or invalid canonical correlation coordinates",
    fixed = TRUE
  )

  us_edge_posterior <- us_primitive_posterior
  us_edge_posterior[1, u_names] <- 0
  expect_error(
    BayesTools:::.bt_random_effect_summary_samples(
      model_samples = us_edge_posterior,
      prior_list = us_result$prior_list,
      formula_design = list(mu = us_result$formula_design),
      mode = "standard"
    ),
    "strictly between 0 and 1",
    fixed = TRUE
  )

  malformed_us_design <- us_result$formula_design
  malformed_us_design$random_effects[[1]]$correlation$primitive_names <- character()
  expect_error(
    BayesTools:::.bt_random_effect_summary_samples(
      model_samples = us_primitive_posterior,
      prior_list = us_result$prior_list,
      formula_design = list(mu = malformed_us_design),
      mode = "standard"
    ),
    "random_term\\$correlation\\$primitive_names"
  )

  partial_scaled_result <- JAGS_formula(
    formula = ~ 1 + us(1 + x | id),
    parameter = "mu",
    data = df_rho,
    prior_list = fixed_priors,
    formula_scale = list(x = TRUE),
    prior_random = prior_random(
      id = random_block(
        sd = sd_prior,
        terms = list(intercept = prior("point", list(2))),
        cor = prior_lkj(eta = 1)
      )
    )
  )
  partial_term <- partial_scaled_result$formula_design$random_effects[[1]]
  partial_L_names <- BayesTools:::.bt_random_effect_cholesky_names(partial_term, 2L)
  partial_source_sd <- c(2, 3)
  partial_source_cor <- matrix(c(1, 0.4, 0.4, 1), 2, 2)
  partial_source_L <- t(chol(partial_source_cor))
  partial_posterior <- matrix(
    c(partial_source_sd[2], as.vector(partial_source_L)),
    nrow = 1,
    dimnames = list(NULL, c("mu__xREx__id_x", as.vector(partial_L_names)))
  )
  partial_summary <- BayesTools:::.bt_random_effect_summary_samples(
    model_samples = partial_posterior,
    prior_list = partial_scaled_result$prior_list,
    formula_design = list(mu = partial_scaled_result$formula_design),
    mode = "standard",
    formula_scale = list(mu = partial_scaled_result$formula_scale)
  )
  partial_scale <- partial_scaled_result$formula_scale$mu_x
  partial_M <- matrix(
    c(1, -partial_scale$mean / partial_scale$sd, 0, 1 / partial_scale$sd),
    nrow = 2,
    byrow = TRUE
  )
  partial_source_cov <- diag(partial_source_sd, nrow = 2) %*% partial_source_cor %*%
    diag(partial_source_sd, nrow = 2)
  partial_expected_cov <- partial_M %*% partial_source_cov %*% t(partial_M)
  partial_expected_sd <- sqrt(diag(partial_expected_cov))
  partial_expected_cor <- partial_expected_cov[1, 2] / prod(partial_expected_sd)

  expect_equal(
    unname(partial_summary$model_samples[1, "mu__xRE_SUMMARY__sd__id__intercept"]),
    partial_expected_sd[1],
    tolerance = 1e-12
  )
  expect_equal(
    unname(partial_summary$model_samples[1, "mu__xRE_SUMMARY__sd__id__x"]),
    partial_expected_sd[2],
    tolerance = 1e-12
  )
  expect_equal(
    unname(partial_summary$model_samples[1, "mu__xRE_SUMMARY__cor__id__intercept__x"]),
    partial_expected_cor,
    tolerance = 1e-12
  )
})

test_that("random-effect formulas are guarded in fixed-only downstream evaluators", {

  df <- data.frame(
    x = c(-1, 0, 1, 2),
    id = factor(c("a", "a", "b", "b"))
  )
  fit <- coda::mcmc(matrix(
    c(0, 0),
    ncol = 1,
    dimnames = list(NULL, "mu_intercept")
  ))
  fixed_prior <- prior("normal", list(0, 1))
  attr(fixed_prior, "parameter") <- "mu"

  expect_error(
    JAGS_evaluate_formula(
      fit = fit,
      formula = ~ 1 + (1 || id),
      parameter = "mu",
      data = df,
      prior_list = list(mu_intercept = fixed_prior)
    ),
    "needs fitted formula design metadata",
    fixed = TRUE
  )

  formula_result <- JAGS_formula(
    formula = ~ 1 + diag(1 | id),
    parameter = "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1))
    ),
    prior_random = prior_random(
      id = random_block(sd = prior("gamma", list(2, 2)))
    )
  )
  attr(fit, "formula_design") <- list(mu = formula_result$formula_design)

  expect_error(
    JAGS_evaluate_formula(
      fit = fit,
      formula = ~ 1 + diag(1 | id),
      parameter = "mu",
      data = df,
      prior_list = formula_result$prior_list
    ),
    "cannot be reconstructed",
    fixed = TRUE
  )

  expect_error(
    JAGS_evaluate_formula(
      fit = fit,
      formula = ~ 1,
      parameter = "mu",
      data = df,
      prior_list = formula_result$prior_list
    ),
    "cannot currently evaluate random-effect fits without silently dropping group-level contributions",
    fixed = TRUE
  )

  expect_error(
    JAGS_bridgesampling(
      fit = fit,
      log_posterior = function(parameters, data) 0,
      data = list(),
      prior_list = list(),
      formula_list = list(mu = ~ 1 + diag(1 | id)),
      formula_data_list = list(mu = df),
      formula_prior_list = list(mu = list(
        intercept = prior("normal", list(0, 1))
      ))
    ),
    "requires 'formula_random_prior_list'",
    fixed = TRUE
  )

  bridge_incomplete_fit <- coda::mcmc(matrix(
    c(0, 1),
    nrow = 1,
    dimnames = list(NULL, c("mu_intercept", "mu__xREx__id_sd"))
  ))
  expect_error(
    JAGS_bridgesampling(
      fit = bridge_incomplete_fit,
      log_posterior = function(parameters, data) 0,
      data = list(),
      prior_list = list(),
      formula_list = list(mu = ~ 1 + id(1 | id)),
      formula_data_list = list(mu = df),
      formula_prior_list = list(mu = list(
        intercept = prior("normal", list(0, 1))
      )),
      formula_random_prior_list = list(mu = prior_random(
        id = random_block(sd = prior("gamma", list(2, 2)))
      ))
    ),
    "requires posterior samples of standardized latent random effects",
    fixed = TRUE
  )
  expect_error(
    JAGS_bridgesampling(
      fit = fit,
      log_posterior = function(parameters, data) 0,
      data = list(),
      prior_list = formula_result$prior_list
    ),
    "requires 'formula_list'",
    fixed = TRUE
  )
})

test_that("formula random effects require prior_random with no legacy fallback", {

  df <- data.frame(
    id = factor(c("a", "a", "b", "b"))
  )

  expect_error(
    JAGS_formula(
      formula = ~ 1 + diag(1 | id),
      parameter = "mu",
      data = df,
      prior_list = list(
        intercept = prior("normal", list(0, 1)),
        "intercept|id" = prior("gamma", list(2, 2))
      )
    ),
    "Formula random effects require 'prior_random'.",
    fixed = TRUE
  )
  expect_error(
    JAGS_fit(
      model_syntax = "model{}",
      formula_list = list(mu = ~ 1 + diag(1 | id)),
      formula_data_list = list(mu = df),
      formula_prior_list = list(mu = list(intercept = prior("normal", list(0, 1))))
    ),
    "JAGS_fit() requires 'formula_random_prior_list'",
    fixed = TRUE
  )

  expect_error(
    JAGS_formula(
      formula = ~ 1 + diag(1 | id),
      parameter = "mu",
      data = df,
      prior_list = list(
        intercept = prior("normal", list(0, 1)),
        "intercept|id" = prior("gamma", list(2, 2))
      ),
      prior_random = prior_random(
        id = random_block(sd = prior("gamma", list(2, 2)))
      )
    ),
    "no longer supported",
    fixed = TRUE
  )
})

test_that("formula random effects expose bridge-ready stochastic coordinates", {

  df <- data.frame(
    x = c(0, 1, 2, 3),
    idx = factor(c("t1", "t2", "t1", "t2"), levels = c("t1", "t2")),
    id = factor(c("a", "a", "b", "b"), levels = c("a", "b"))
  )

  formula_result <- JAGS_formula(
    formula = ~ 1 + x + diag(1 + x | id),
    parameter = "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    ),
    prior_random = prior_random(
      id = random_block(sd = prior("gamma", list(2, 2)))
    )
  )

  bridge <- BayesTools:::.bt_JAGS_formula_random_bridge_parameters(
    list(mu = formula_result$formula_design)
  )
  expected_z <- c(
    "mu__xREx__id_xRE_Zx[1,1]",
    "mu__xREx__id_xRE_Zx[2,1]",
    "mu__xREx__id_xRE_Zx[1,2]",
    "mu__xREx__id_xRE_Zx[2,2]"
  )
  expect_equal(bridge$parameters, expected_z)
  expect_equal(bridge$bounds$lb, stats::setNames(rep(-Inf, 4), expected_z))
  expect_equal(bridge$bounds$ub, stats::setNames(rep( Inf, 4), expected_z))
  expect_error(
    BayesTools:::.bt_JAGS_bridge_merge_add_parameters(
      add_parameters = expected_z[1],
      add_bounds = list(
        lb = stats::setNames(0, expected_z[1]),
        ub = stats::setNames(Inf, expected_z[1])
      ),
      bridge_parameters = bridge$parameters,
      bridge_bounds = bridge$bounds
    ),
    "conflict",
    fixed = TRUE
  )
  expect_error(
    BayesTools:::.bt_JAGS_bridge_merge_add_parameters(
      add_parameters = expected_z[1],
      add_bounds = list(
        lb = stats::setNames(-Inf, "wrong_name"),
        ub = stats::setNames( Inf, expected_z[1])
      ),
      bridge_parameters = bridge$parameters,
      bridge_bounds = bridge$bounds
    ),
    "names must match",
    fixed = TRUE
  )
  expect_error(
    BayesTools:::.bt_JAGS_bridge_merge_add_parameters(
      add_parameters = c(expected_z[1], expected_z[2]),
      add_bounds = list(
        lb = stats::setNames(c(-Inf, -Inf), c(expected_z[1], expected_z[1])),
        ub = stats::setNames(c( Inf,  Inf), c(expected_z[1], expected_z[2]))
      ),
      bridge_parameters = bridge$parameters,
      bridge_bounds = bridge$bounds
    ),
    "unique",
    fixed = TRUE
  )
  merged_bridge <- BayesTools:::.bt_JAGS_bridge_merge_add_parameters(
    add_parameters = expected_z[1],
    add_bounds = list(
      lb = stats::setNames(-Inf, expected_z[1]),
      ub = stats::setNames( Inf, expected_z[1])
    ),
    bridge_parameters = bridge$parameters,
    bridge_bounds = bridge$bounds
  )
  expect_equal(merged_bridge$add_parameters, expected_z)

  samples <- c(
    mu_intercept = 10,
    mu_x = 1,
    mu__xREx__id_intercept = 0.5,
    mu__xREx__id_x = 0.25,
    "mu__xREx__id_xRE_Zx[1,1]" = 1,
    "mu__xREx__id_xRE_Zx[1,2]" = 2,
    "mu__xREx__id_xRE_Zx[2,1]" = 3,
    "mu__xREx__id_xRE_Zx[2,2]" = 4
  )
  parameters <- JAGS_marglik_parameters_formula(
    samples = samples,
    formula_list = list(mu = formula_result$formula),
    formula_data_list = list(mu = formula_result$data),
    formula_prior_list = list(mu = formula_result$prior_list),
    prior_list_parameters = list(),
    formula_design_list = list(mu = formula_result$formula_design)
  )

  expect_equal(parameters$mu, c(10.5, 12, 15.5, 17.5), tolerance = 1e-12)
  expect_equal(
    BayesTools:::.bt_JAGS_marglik_priors_formula_random(
      samples,
      list(mu = formula_result$formula_design)
    ),
    sum(stats::dnorm(samples[expected_z], log = TRUE)),
    tolerance = 1e-12
  )
})

test_that("formula random-effect bridge helpers handle LKJ and scalar rho blocks", {

  df <- data.frame(
    x = c(0, 1, 2, 3),
    idx = factor(c("t1", "t2", "t1", "t2"), levels = c("t1", "t2")),
    id = factor(c("a", "a", "b", "b"), levels = c("a", "b"))
  )

  lkj_result <- JAGS_formula(
    formula = ~ 1 + x + us(1 + x | id),
    parameter = "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    ),
    prior_random = prior_random(
      id = random_block(
        sd = prior("gamma", list(2, 2)),
        cor = prior_lkj(eta = 1.5)
      )
    )
  )
  lkj_u <- "mu__xREx__id_xRE_CORx_lkj_u[1]"
  expect_true(lkj_u %in% lkj_result$add_parameters)
  lkj_bridge <- BayesTools:::.bt_JAGS_formula_random_bridge_parameters(
    list(mu = lkj_result$formula_design)
  )
  expect_true(lkj_u %in% lkj_bridge$parameters)
  expect_equal(lkj_bridge$bounds$lb[[lkj_u]], 0)
  expect_equal(lkj_bridge$bounds$ub[[lkj_u]], 1)

  bridge_fitted <- list(mu = lkj_result$formula_design)
  bridge_changed <- bridge_fitted
  bridge_changed$mu$random_effects[[1]]$homogeneous_sd <- NULL
  expect_error(
    BayesTools:::.bt_JAGS_bridge_validate_formula_random_designs(
      bridge_fitted,
      bridge_changed
    ),
    "missing canonical 'random_term\\$homogeneous_sd'"
  )

  bridge_changed <- bridge_fitted
  bridge_changed$mu$random_effects[[1]]$correlation <- NULL
  expect_error(
    BayesTools:::.bt_JAGS_bridge_validate_formula_random_designs(
      bridge_fitted,
      bridge_changed
    ),
    "missing canonical 'random_term\\$correlation'"
  )

  malformed_lkj <- lkj_result$formula_design
  malformed_lkj$random_effects[[1]]$structure <- NULL
  malformed_lkj$random_effects[[1]]$covariance <- "us"
  expect_error(
    BayesTools:::.bt_JAGS_formula_random_bridge_parameters(
      list(mu = malformed_lkj)
    ),
    "missing canonical 'random_term\\$structure'"
  )

  lkj_samples <- c(
    mu_intercept = 0,
    mu_x = 0,
    mu__xREx__id_intercept = 0.5,
    mu__xREx__id_x = 0.25,
    "mu__xREx__id_xRE_Zx[1,1]" = 1,
    "mu__xREx__id_xRE_Zx[1,2]" = 2,
    "mu__xREx__id_xRE_Zx[2,1]" = 3,
    "mu__xREx__id_xRE_Zx[2,2]" = 4,
    "mu__xREx__id_xRE_CORx_lkj_u[1]" = 0.75
  )
  lkj_parameters <- JAGS_marglik_parameters_formula(
    samples = lkj_samples,
    formula_list = list(mu = lkj_result$formula),
    formula_data_list = list(mu = lkj_result$data),
    formula_prior_list = list(mu = lkj_result$prior_list),
    prior_list_parameters = list(),
    formula_design_list = list(mu = lkj_result$formula_design)
  )
  L <- BayesTools:::.bt_lkj_cholesky_cpc_u_to_L(0.75, K = 2)
  z <- matrix(c(1, 3, 2, 4), nrow = 2, ncol = 2)
  coef <- sweep(z %*% t(L), 2, c(0.5, 0.25), "*")
  expected_lkj <- rowSums(coef[lkj_result$formula_design$random_effects[[1]]$group_map, ] *
    lkj_result$formula_design$random_effects[[1]]$model_matrix)
  expect_equal(lkj_parameters$mu, unname(expected_lkj), tolerance = 1e-12)
  expect_equal(
    BayesTools:::.bt_JAGS_marglik_priors_formula_random(
      lkj_samples,
      list(mu = lkj_result$formula_design)
    ),
    sum(stats::dnorm(lkj_samples[grep("xRE_Zx", names(lkj_samples), value = TRUE)], log = TRUE)) +
      BayesTools:::.bt_lkj_cholesky_cpc_u_log_prior(0.75, K = 2, eta = 1.5),
    tolerance = 1e-12
  )
  expect_error(
    BayesTools:::.bt_JAGS_marglik_priors_formula_random(
      lkj_samples,
      list(mu = malformed_lkj)
    ),
    "missing canonical 'random_term\\$structure'"
  )

  missing_correlation_lkj <- lkj_result$formula_design
  missing_correlation_lkj$random_effects[[1]]$correlation <- NULL
  expect_error(
    BayesTools:::.bt_JAGS_marglik_priors_formula_random(
      lkj_samples,
      list(mu = missing_correlation_lkj)
    ),
    "missing canonical 'random_term\\$correlation'"
  )
  expect_error(
    BayesTools:::.bt_JAGS_marglik_random_effect_cholesky(
      lkj_samples,
      missing_correlation_lkj$random_effects[[1]]
    ),
    "missing canonical 'random_term\\$correlation'"
  )
  expect_error(
    BayesTools:::.bt_random_effect_summary_samples(
      model_samples = matrix(
        unname(lkj_samples),
        nrow = 1,
        dimnames = list(NULL, names(lkj_samples))
      ),
      prior_list = lkj_result$prior_list,
      formula_design = list(mu = missing_correlation_lkj),
      mode = "standard"
    ),
    "missing canonical 'random_term\\$correlation'"
  )

  missing_eta_lkj <- lkj_result$formula_design
  missing_eta_lkj$random_effects[[1]]$correlation$eta <- NULL
  expect_error(
    BayesTools:::.bt_JAGS_marglik_priors_formula_random(
      lkj_samples,
      list(mu = missing_eta_lkj)
    ),
    "missing canonical 'random_term\\$correlation\\$eta'"
  )
  expect_error(
    BayesTools:::.bt_JAGS_marglik_random_effect_cholesky(
      lkj_samples,
      malformed_lkj$random_effects[[1]]
    ),
    "missing canonical 'random_term\\$structure'"
  )
  expect_error(
    BayesTools:::.bt_JAGS_marglik_random_effect_cholesky(
      lkj_samples[names(lkj_samples) != lkj_u],
      lkj_result$formula_design$random_effects[[1]]
    ),
    "missing random-effect correlation coordinates",
    fixed = TRUE
  )

  ar_result <- JAGS_formula(
    formula = ~ 1 + x + ar1(idx | id),
    parameter = "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    ),
    prior_random = prior_random(
      id = random_block(
        sd = prior("gamma", list(2, 2)),
        rho = prior("normal", list(0, 0.5))
      )
    )
  )
  expect_equal(ar_result$formula_design$random_effects[[1]]$correlation$sample_name, "mu__xREx__id_rho_z")

  ar_rho_bridge <- BayesTools:::.bt_JAGS_formula_random_scalar_rho_bridge_parameters(
    list(mu = ar_result$formula_design)
  )
  expect_equal(ar_rho_bridge$parameters, "mu__xREx__id_rho_z")
  expect_equal(ar_rho_bridge$bounds$lb[["mu__xREx__id_rho_z"]], -Inf)
  expect_equal(ar_rho_bridge$bounds$ub[["mu__xREx__id_rho_z"]],  Inf)

  ar_samples <- c(
    mu_intercept = 0,
    mu_x = 0,
    mu__xREx__id_sd = 0.5,
    mu__xREx__id_rho_z = atanh(0.25),
    "mu__xREx__id_xRE_Zx[1,1]" = 1,
    "mu__xREx__id_xRE_Zx[1,2]" = 2,
    "mu__xREx__id_xRE_Zx[2,1]" = 3,
    "mu__xREx__id_xRE_Zx[2,2]" = 4
  )
  ar_parameters <- JAGS_marglik_parameters_formula(
    samples = ar_samples,
    formula_list = list(mu = ar_result$formula),
    formula_data_list = list(mu = ar_result$data),
    formula_prior_list = list(mu = ar_result$prior_list),
    prior_list_parameters = list(),
    formula_design_list = list(mu = ar_result$formula_design)
  )
  ar_L <- t(chol(matrix(c(1, 0.25, 0.25, 1), nrow = 2)))
  ar_coef <- 0.5 * (z %*% t(ar_L))
  expected_ar <- rowSums(ar_coef[ar_result$formula_design$random_effects[[1]]$group_map, ] *
    ar_result$formula_design$random_effects[[1]]$model_matrix)
  expect_equal(ar_parameters$mu, unname(expected_ar), tolerance = 1e-12)
  expect_equal(
    BayesTools:::.bt_JAGS_marglik_random_effect_rho(
      c(ar_samples, mu__xREx__id_rho = 0.99),
      ar_result$formula_design$random_effects[[1]]
    ),
    0.25,
    tolerance = 1e-12
  )

  missing_rho_scale <- ar_result$formula_design
  missing_rho_scale$random_effects[[1]]$correlation$rho_scale <- NULL
  expect_error(
    BayesTools:::.bt_JAGS_marglik_random_effect_rho(
      ar_samples,
      missing_rho_scale$random_effects[[1]]
    ),
    "missing canonical 'random_term\\$correlation\\$rho_scale'"
  )

  missing_rho_bounds <- ar_result$formula_design
  missing_rho_bounds$random_effects[[1]]$correlation$bounds <- NULL
  expect_error(
    BayesTools:::.bt_JAGS_marglik_random_effect_rho(
      ar_samples,
      missing_rho_bounds$random_effects[[1]]
    ),
    "missing canonical 'random_term\\$correlation\\$bounds'"
  )

  car_df <- data.frame(
    time = c(0, 0.5, 2, 0, 0.5, 2),
    id = factor(c("a", "a", "a", "b", "b", "b"), levels = c("a", "b"))
  )
  car_result <- JAGS_formula(
    formula = ~ 1 + car(time | id),
    parameter = "mu",
    data = car_df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(
        sd = prior("gamma", list(2, 2)),
        rho = prior("normal", list(0, 0.5))
      )
    )
  )
  expect_equal(car_result$formula_design$random_effects[[1]]$correlation$sample_name, "mu__xREx__id_rho_z")
  expect_equal(car_result$formula_design$random_effects[[1]]$correlation$time_values, c(0, 0.5, 2))

  car_samples <- c(
    mu_intercept = 0,
    mu__xREx__id_sd = 0.5,
    mu__xREx__id_rho_z = atanh(0.25),
    "mu__xREx__id_xRE_Zx[1,1]" = 1,
    "mu__xREx__id_xRE_Zx[1,2]" = 2,
    "mu__xREx__id_xRE_Zx[1,3]" = 3,
    "mu__xREx__id_xRE_Zx[2,1]" = 4,
    "mu__xREx__id_xRE_Zx[2,2]" = 5,
    "mu__xREx__id_xRE_Zx[2,3]" = 6
  )
  car_parameters <- JAGS_marglik_parameters_formula(
    samples = car_samples,
    formula_list = list(mu = car_result$formula),
    formula_data_list = list(mu = car_result$data),
    formula_prior_list = list(mu = car_result$prior_list),
    prior_list_parameters = list(),
    formula_design_list = list(mu = car_result$formula_design)
  )
  car_R <- 0.25^abs(outer(c(0, 0.5, 2), c(0, 0.5, 2), "-"))
  car_L <- t(chol(car_R))
  car_z <- matrix(c(1, 4, 2, 5, 3, 6), nrow = 2, ncol = 3)
  car_coef <- 0.5 * (car_z %*% t(car_L))
  expected_car <- rowSums(car_coef[car_result$formula_design$random_effects[[1]]$group_map, ] *
    car_result$formula_design$random_effects[[1]]$model_matrix)
  expect_equal(car_parameters$mu, unname(expected_car), tolerance = 1e-12)

  car_logit <- car_result
  car_logit$formula_design$random_effects[[1]]$correlation$rho_scale <- "logit"
  car_logit$formula_design$random_effects[[1]]$correlation$sample_name <- "mu__xREx__id_rho_logit"
  expect_equal(
    BayesTools:::.bt_JAGS_marglik_random_effect_rho(
      c(car_samples[names(car_samples) != "mu__xREx__id_rho_z"], mu__xREx__id_rho_logit = stats::qlogis(0.25)),
      car_logit$formula_design$random_effects[[1]]
    ),
    0.25,
    tolerance = 1e-12
  )

  hcs_df <- data.frame(
    f = factor(rep(c("a", "b", "c"), 2), levels = c("a", "b", "c")),
    id = factor(rep(c("g1", "g2"), each = 3), levels = c("g1", "g2"))
  )
  hcs_result <- JAGS_formula(
    formula = ~ 1 + hcs(f | id),
    parameter = "mu",
    data = hcs_df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(
        sd = prior("gamma", list(2, 2)),
        rho = prior("normal", list(0, 0.5))
      )
    )
  )
  hcs_random_term <- hcs_result$formula_design$random_effects[[1]]
  hcs_rho_bridge <- BayesTools:::.bt_JAGS_formula_random_scalar_rho_bridge_parameters(
    list(mu = hcs_result$formula_design)
  )
  expect_equal(hcs_rho_bridge$parameters, "mu__xREx__id_rho_z")
  expect_equal(
    hcs_rho_bridge$bounds$lb[["mu__xREx__id_rho_z"]],
    atanh(-1 / 2),
    tolerance = 1e-12
  )
  expect_equal(hcs_rho_bridge$bounds$ub[["mu__xREx__id_rho_z"]], Inf)

  hcs_z_names <- as.vector(BayesTools:::.bt_random_effect_latent_names(
    random_term = hcs_random_term,
    n_groups = hcs_random_term$n_groups,
    n_columns = hcs_random_term$n_columns
  ))
  hcs_boundary_samples <- c(
    stats::setNames(rep(0, length(hcs_z_names)), hcs_z_names),
    mu__xREx__id_rho_z = atanh(-1 / 2)
  )
  expect_equal(
    BayesTools:::.bt_JAGS_marglik_priors_formula_random(
      hcs_boundary_samples,
      list(mu = hcs_result$formula_design)
    ),
    -Inf
  )
  hcs_interior_samples <- hcs_boundary_samples
  hcs_interior_samples[["mu__xREx__id_rho_z"]] <- 0
  expect_true(is.finite(
    BayesTools:::.bt_JAGS_marglik_priors_formula_random(
      hcs_interior_samples,
      list(mu = hcs_result$formula_design)
    )
  ))

  hcs_fixed_rho <- JAGS_formula(
    formula = ~ 1 + hcs(f | id),
    parameter = "mu",
    data = hcs_df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(
        sd = prior("gamma", list(2, 2)),
        rho = prior("point", list(location = 0.2))
      )
    )
  )
  hcs_fixed_term <- hcs_fixed_rho$formula_design$random_effects[[1]]
  expect_null(BayesTools:::.bt_JAGS_bridge_scalar_rho_parameter(hcs_fixed_term))
  expect_equal(hcs_fixed_term$correlation$sample_fixed, 0.2)
  hcs_fixed_z_names <- as.vector(BayesTools:::.bt_random_effect_latent_names(
    random_term = hcs_fixed_term,
    n_groups = hcs_fixed_term$n_groups,
    n_columns = hcs_fixed_term$n_columns
  ))
  hcs_fixed_samples <- stats::setNames(rep(0, length(hcs_fixed_z_names)), hcs_fixed_z_names)
  expect_true(is.finite(
    BayesTools:::.bt_JAGS_marglik_priors_formula_random(
      hcs_fixed_samples,
      list(mu = hcs_fixed_rho$formula_design)
    )
  ))
  expect_equal(
    BayesTools:::.bt_JAGS_marglik_random_effect_rho(
      hcs_fixed_samples,
      hcs_fixed_term
    ),
    tanh(0.2),
    tolerance = 1e-12
  )
})

test_that("JAGS_evaluate_formula evaluates monitored observed random effects", {

  df <- data.frame(
    x = c(-1, 0, 1, 2),
    id = factor(c("a", "a", "b", "b"), levels = c("a", "b"))
  )
  formula_result <- JAGS_formula(
    formula = ~ 1 + x + diag(1 + x | id),
    parameter = "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    ),
    prior_random = prior_random(
      id = random_block(
        sd = prior("gamma", list(2, 2)),
        monitor = random_monitor(latent = FALSE, coefficients = TRUE, correlation = FALSE)
      )
    )
  )
  posterior <- matrix(
    c(
      10, 1,  0.5, 0.1, -0.5, -0.2,
      20, 2,  0.6, 0.2, -0.6, -0.3
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(NULL, c(
      "mu_intercept",
      "mu_x",
      "mu__xREx__id_xRE_COEFx[1,1]",
      "mu__xREx__id_xRE_COEFx[1,2]",
      "mu__xREx__id_xRE_COEFx[2,1]",
      "mu__xREx__id_xRE_COEFx[2,2]"
    ))
  )
  fit <- coda::mcmc(posterior)
  attr(fit, "formula_design") <- list(mu = formula_result$formula_design)

  prediction <- JAGS_evaluate_formula(
    fit = fit,
    formula = ~ 1 + x + diag(1 + x | id),
    parameter = "mu",
    data = df,
    prior_list = formula_result$prior_list
  )

  expected <- cbind(
    c(9.4, 10.5, 10.3, 11.1),
    c(18.4, 20.6, 21.1, 22.8)
  )
  expect_equal(unname(prediction), expected)

  expect_error(
    JAGS_evaluate_formula(
      fit = fit,
      formula = ~ 1 + x + diag(1 + x | id),
      parameter = "mu",
      data = data.frame(x = 0, id = factor("c", levels = "c")),
      prior_list = formula_result$prior_list
    ),
    "New random-effect level",
    fixed = TRUE
  )
  expect_error(
    JAGS_evaluate_formula(
      fit = fit,
      formula = ~ 1 + x +
        diag(1 + x | id) +
        random(1 | drug, name = "drug", covariance = "diag"),
      parameter = "mu",
      data = transform(df, drug = factor(c("d1", "d1", "d2", "d2"))),
      prior_list = formula_result$prior_list
    ),
    "do not match the fitted formula",
    fixed = TRUE
  )
  expect_error(
    JAGS_evaluate_formula(
      fit = fit,
      formula = ~ 1 + x + random(1 + x | drug, name = "id", covariance = "diag"),
      parameter = "mu",
      data = transform(df, drug = id),
      prior_list = formula_result$prior_list
    ),
    "does not match the fitted formula",
    fixed = TRUE
  )
  expect_error(
    JAGS_evaluate_formula(
      fit = fit,
      formula = ~ 1 + x + us(1 + x | id),
      parameter = "mu",
      data = df,
      prior_list = formula_result$prior_list
    ),
    "does not match the fitted formula",
    fixed = TRUE
  )
  expect_error(
    JAGS_evaluate_formula(
      fit = fit,
      formula = ~ 1 + x + diag(1 + x | id, hom = TRUE),
      parameter = "mu",
      data = df,
      prior_list = formula_result$prior_list
    ),
    "does not match the fitted formula",
    fixed = TRUE
  )

  slope_only_result <- JAGS_formula(
    formula = ~ 1 + diag(0 + x | id),
    parameter = "mu",
    data = df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(
        sd = prior("gamma", list(2, 2)),
        monitor = random_monitor(latent = FALSE, coefficients = TRUE, correlation = FALSE)
      )
    )
  )
  slope_only_fit <- coda::mcmc(matrix(
    c(0, 1, 2),
    nrow = 1,
    dimnames = list(NULL, c(
      "mu_intercept",
      "mu__xREx__id_xRE_COEFx[1,1]",
      "mu__xREx__id_xRE_COEFx[2,1]"
    ))
  ))
  attr(slope_only_fit, "formula_design") <- list(mu = slope_only_result$formula_design)
  expect_error(
    JAGS_evaluate_formula(
      fit = slope_only_fit,
      formula = ~ 1 + diag(0 + x | id),
      parameter = "mu",
      data = transform(df, x = c(-1, NA, 1, 2)),
      prior_list = slope_only_result$prior_list
    ),
    "missing predictor values",
    fixed = TRUE
  )
})

test_that("JAGS_evaluate_formula reconstructs observed random effects from latent samples", {

  df <- data.frame(
    x = c(-1, 0, 1, 2),
    id = factor(c("a", "a", "b", "b"), levels = c("a", "b"))
  )
  sd_prior <- prior("gamma", list(2, 2))
  formula_result <- JAGS_formula(
    formula = ~ 1 + x + diag(1 + x | id),
    parameter = "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    ),
    prior_random = prior_random(
      id = random_block(
        sd = sd_prior,
        monitor = random_monitor(latent = TRUE, coefficients = FALSE, correlation = FALSE)
      )
    )
  )
  expect_true("mu__xREx__id_xRE_Zx" %in% formula_result$add_parameters)
  expect_false("mu__xREx__id_xRE_COEFx" %in% formula_result$add_parameters)

  posterior <- matrix(
    c(
      10, 1, 2, 3,  0.5, 0.1, -0.5, -0.2,
      20, 2, 4, 5,  0.6, 0.2, -0.6, -0.3
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(NULL, c(
      "mu_intercept",
      "mu_x",
      "mu__xREx__id_intercept",
      "mu__xREx__id_x",
      "mu__xREx__id_xRE_Zx[1,1]",
      "mu__xREx__id_xRE_Zx[1,2]",
      "mu__xREx__id_xRE_Zx[2,1]",
      "mu__xREx__id_xRE_Zx[2,2]"
    ))
  )
  fit <- coda::mcmc(posterior)
  attr(fit, "formula_design") <- list(mu = formula_result$formula_design)

  prediction <- JAGS_evaluate_formula(
    fit = fit,
    formula = ~ 1 + x + diag(1 + x | id),
    parameter = "mu",
    data = df,
    prior_list = formula_result$prior_list
  )

  expected <- cbind(
    c(9.7, 11.0, 9.4, 9.8),
    c(19.4, 22.4, 18.1, 18.6)
  )
  expect_equal(unname(prediction), expected)

  incomplete_fit <- coda::mcmc(posterior[, colnames(posterior) != "mu__xREx__id_xRE_Zx[2,2]", drop = FALSE])
  attr(incomplete_fit, "formula_design") <- list(mu = formula_result$formula_design)
  expect_error(
    JAGS_evaluate_formula(
      fit = incomplete_fit,
      formula = ~ 1 + x + diag(1 + x | id),
      parameter = "mu",
      data = df,
      prior_list = formula_result$prior_list
    ),
    "cannot be reconstructed from the posterior samples",
    fixed = TRUE
  )
})

test_that("JAGS_evaluate_formula reconstructs allocated random effects from simplex coordinates", {

  df <- data.frame(
    study = factor(c("s1", "s1", "s2", "s2", "s3", "s3")),
    drug = factor(c("a", "b", "a", "b", "a", "b"))
  )
  formula_result <- JAGS_formula(
    formula = ~ 1 +
      random(1 | study, name = "study", covariance = "diag") +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      allocation = random_variance_allocation(
        sd = prior("gamma", list(2, 2)),
        allocation = prior("dirichlet", list(alpha = c(2, 3)))
      )
    )
  )

  posterior <- matrix(
    c(
      10, 2, 1, 3,  0.1, 0.2, 0.3,  1, 2,
      20, 4, 2, 2,  0.5, 0.6, 0.7, -1, 1
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(NULL, c(
      "mu_intercept",
      "mu__xRE_ALLOCx_allocation_total_sd",
      "prior_par_eta_mu__xRE_ALLOCx_allocation_weight[1]",
      "prior_par_eta_mu__xRE_ALLOCx_allocation_weight[2]",
      "mu__xREx__study_xRE_Zx[1,1]",
      "mu__xREx__study_xRE_Zx[2,1]",
      "mu__xREx__study_xRE_Zx[3,1]",
      "mu__xREx__drug_xRE_Zx[1,1]",
      "mu__xREx__drug_xRE_Zx[2,1]"
    ))
  )
  fit <- coda::mcmc(posterior)
  attr(fit, "formula_design") <- list(mu = formula_result$formula_design)

  prediction <- JAGS_evaluate_formula(
    fit = fit,
    formula = ~ 1 +
      random(1 | study, name = "study", covariance = "diag") +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = df,
    prior_list = formula_result$prior_list
  )

  first_draw <- 10 +
    c(0.1, 0.1, 0.2, 0.2, 0.3, 0.3) +
    sqrt(3) * c(1, 2, 1, 2, 1, 2)
  second_sd <- 4 * sqrt(0.5)
  second_draw <- 20 +
    second_sd * c(0.5, 0.5, 0.6, 0.6, 0.7, 0.7) +
    second_sd * c(-1, 1, -1, 1, -1, 1)
  expect_equal(unname(prediction), unname(cbind(first_draw, second_draw)), tolerance = 1e-12)

  public_simplex <- matrix(
    c(
      5, 2, 0.25, 0.75,  0.1, 0.2, 0.3,  1, 2
    ),
    nrow = 1,
    dimnames = list(NULL, c(
      "mu_intercept",
      "mu__xRE_ALLOCx_allocation_total_sd",
      "mu__xRE_ALLOCx_allocation_weight[1]",
      "mu__xRE_ALLOCx_allocation_weight[2]",
      "mu__xREx__study_xRE_Zx[1,1]",
      "mu__xREx__study_xRE_Zx[2,1]",
      "mu__xREx__study_xRE_Zx[3,1]",
      "mu__xREx__drug_xRE_Zx[1,1]",
      "mu__xREx__drug_xRE_Zx[2,1]"
    ))
  )
  public_fit <- coda::mcmc(public_simplex)
  attr(public_fit, "formula_design") <- list(mu = formula_result$formula_design)
  public_prediction <- JAGS_evaluate_formula(
    fit = public_fit,
    formula = ~ 1 +
      random(1 | study, name = "study", covariance = "diag") +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = df,
    prior_list = formula_result$prior_list
  )
  expect_equal(
    unname(drop(public_prediction)),
    5 +
      c(0.1, 0.1, 0.2, 0.2, 0.3, 0.3) +
      sqrt(3) * c(1, 2, 1, 2, 1, 2),
    tolerance = 1e-12
  )
})

test_that("JAGS_evaluate_formula reconstructs correlated latent random effects", {

  df <- data.frame(
    x = c(0, 1, 2, -1),
    id = factor(c("a", "a", "b", "b"), levels = c("a", "b"))
  )
  formula_result <- JAGS_formula(
    formula = ~ 1 + x + random(1 + x | id, name = "study"),
    parameter = "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    ),
    prior_random = prior_random(
      study = random_block(
        sd = prior("gamma", list(2, 2)),
        cor = prior_lkj(eta = 2),
        monitor = random_monitor(latent = TRUE, coefficients = FALSE, correlation = FALSE)
      )
    )
  )
  expect_true("mu__xREx__study_xRE_Zx" %in% formula_result$add_parameters)
  expect_true("mu__xREx__study_xRE_CORx_L" %in% formula_result$add_parameters)
  expect_false("mu__xREx__study_xRE_COEFx" %in% formula_result$add_parameters)

  L22 <- sqrt(1 - 0.5^2)
  posterior <- matrix(
    c(
      0, 0,
      2, 3,
      1, 0, 0.5, L22,
      1, 2, -1, 0.25
    ),
    nrow = 1,
    dimnames = list(NULL, c(
      "mu_intercept",
      "mu_x",
      "mu__xREx__study_intercept",
      "mu__xREx__study_x",
      "mu__xREx__study_xRE_CORx_L[1,1]",
      "mu__xREx__study_xRE_CORx_L[1,2]",
      "mu__xREx__study_xRE_CORx_L[2,1]",
      "mu__xREx__study_xRE_CORx_L[2,2]",
      "mu__xREx__study_xRE_Zx[1,1]",
      "mu__xREx__study_xRE_Zx[1,2]",
      "mu__xREx__study_xRE_Zx[2,1]",
      "mu__xREx__study_xRE_Zx[2,2]"
    ))
  )
  fit <- coda::mcmc(posterior)
  attr(fit, "formula_design") <- list(mu = formula_result$formula_design)

  prediction <- JAGS_evaluate_formula(
    fit = fit,
    formula = ~ 1 + x + random(1 + x | id, name = "study"),
    parameter = "mu",
    data = df,
    prior_list = formula_result$prior_list
  )

  coef_a_x <- 3 * (0.5 * 1 + L22 * 2)
  coef_b_x <- 3 * (0.5 * -1 + L22 * 0.25)
  expect_equal(
    unname(drop(prediction)),
    c(2, 2 + coef_a_x, -2 + 2 * coef_b_x, -2 - coef_b_x),
    tolerance = 1e-12
  )
})

test_that("JAGS_evaluate_formula reconstructs point-SD latent random effects", {

  df <- data.frame(id = factor(c("a", "b"), levels = c("a", "b")))
  formula_result <- JAGS_formula(
    formula = ~ 1 + diag(1 | id),
    parameter = "mu",
    data = df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(
        sd = prior("point", list(2)),
        monitor = random_monitor(latent = TRUE, coefficients = FALSE, correlation = FALSE)
      )
    )
  )
  posterior <- matrix(
    c(1, 0.25, -0.5),
    nrow = 1,
    dimnames = list(NULL, c(
      "mu_intercept",
      "mu__xREx__id_xRE_Zx[1,1]",
      "mu__xREx__id_xRE_Zx[2,1]"
    ))
  )
  fit <- coda::mcmc(posterior)
  attr(fit, "formula_design") <- list(mu = formula_result$formula_design)

  prediction <- JAGS_evaluate_formula(
    fit = fit,
    formula = ~ 1 + diag(1 | id),
    parameter = "mu",
    data = df,
    prior_list = formula_result$prior_list
  )

  expect_equal(unname(drop(prediction)), c(1.5, 0))
})

test_that("JAGS_evaluate_formula scales random-slope prediction data", {

  df <- data.frame(
    x = c(1, 2, 3, 4),
    id = factor(c("a", "a", "b", "b"), levels = c("a", "b"))
  )
  formula_result <- JAGS_formula(
    formula = ~ 1 + x + diag(0 + x | id),
    parameter = "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    ),
    formula_scale = list(x = TRUE),
    prior_random = prior_random(
      id = random_block(
        sd = prior("gamma", list(2, 2)),
        monitor = random_monitor(latent = TRUE, coefficients = FALSE, correlation = FALSE)
      )
    )
  )
  scale_info <- formula_result$formula_scale$mu_x
  posterior <- matrix(
    c(0, 0, 2, 1, -1),
    nrow = 1,
    dimnames = list(NULL, c(
      "mu_intercept",
      "mu_x",
      "mu__xREx__id_x",
      "mu__xREx__id_xRE_Zx[1,1]",
      "mu__xREx__id_xRE_Zx[2,1]"
    ))
  )
  fit <- coda::mcmc(posterior)
  attr(fit, "formula_design") <- list(mu = formula_result$formula_design)
  attr(fit, "formula_scale") <- list(mu = formula_result$formula_scale)

  newdata <- data.frame(
    x = c(2, 4),
    id = factor(c("a", "b"), levels = c("a", "b"))
  )
  prediction <- JAGS_evaluate_formula(
    fit = fit,
    formula = ~ 1 + x + diag(0 + x | id),
    parameter = "mu",
    data = newdata,
    prior_list = formula_result$prior_list
  )

  expected_x <- (newdata$x - scale_info$mean) / scale_info$sd
  expect_equal(unname(drop(prediction)), expected_x * c(2, -2))
})

test_that("structured random-effect prediction indexes raw data when fixed predictors are scaled", {

  df <- data.frame(
    time = rep(c(1, 2, 3), 2),
    id = factor(rep(c("a", "b"), each = 3), levels = c("a", "b"))
  )
  sd_prior <- prior("normal", list(0, 1), truncation = list(lower = 0, upper = Inf))
  rho_prior <- prior("normal", list(0, 0.5))

  for (structure in c("cs", "hcs", "ar1", "car", "har")) {
    structured_formula <- stats::as.formula(
      if (identical(structure, "car")) {
        "~ 1 + time + car(0 + time | id)"
      } else {
        paste0("~ 1 + time + ", structure, "(time | id)")
      }
    )
    formula_result <- JAGS_formula(
      formula = structured_formula,
      parameter = "mu",
      data = df,
      prior_list = list(
        intercept = prior("normal", list(0, 1)),
        time = prior("normal", list(0, 1))
      ),
      formula_scale = list(time = TRUE),
      prior_random = prior_random(
        id = random_block(
          sd = sd_prior,
          rho = rho_prior,
          monitor = random_monitor(latent = FALSE, coefficients = TRUE, correlation = FALSE)
        )
      )
    )
    random_term <- formula_result$formula_design$random_effects[[1]]
    coefficient_names <- BayesTools:::.bt_random_effect_coefficient_names(
      random_term = random_term,
      n_groups = length(random_term$group_levels),
      n_columns = random_term$n_columns
    )
    target_column_name <- if (identical(structure, "car")) "time_2" else "time2"
    target_column <- match(target_column_name, random_term$column_names)
    expect_false(is.na(target_column), info = structure)

    posterior_names <- c("mu_intercept", "mu_time", as.vector(coefficient_names))
    posterior <- matrix(
      0,
      nrow = 1,
      ncol = length(posterior_names),
      dimnames = list(NULL, posterior_names)
    )
    posterior[, coefficient_names[match("b", random_term$group_levels), target_column]] <- 42
    fit <- coda::mcmc(posterior)
    attr(fit, "formula_design") <- list(mu = formula_result$formula_design)
    attr(fit, "formula_scale") <- list(mu = formula_result$formula_scale)

    prediction <- JAGS_evaluate_formula(
      fit = fit,
      formula = structured_formula,
      parameter = "mu",
      data = data.frame(time = 2, id = factor("b", levels = c("a", "b"))),
      prior_list = formula_result$prior_list
    )
    expect_equal(unname(drop(prediction)), 42, info = structure)
  }
})

test_that("JAGS_evaluate_formula rejects duplicate requested random-effect block names", {

  df <- data.frame(id = factor(c("a", "b"), levels = c("a", "b")))
  formula_result <- JAGS_formula(
    formula = ~ 1 + diag(1 | id),
    parameter = "mu",
    data = df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(sd = prior("gamma", list(2, 2)))
    )
  )
  fit <- coda::mcmc(matrix(
    0,
    nrow = 1,
    dimnames = list(NULL, "mu_intercept")
  ))
  attr(fit, "formula_design") <- list(mu = formula_result$formula_design)

  expect_error(
    JAGS_evaluate_formula(
      fit = fit,
      formula = ~ 1 + diag(1 | id) + random(1 | id, name = "id", covariance = "diag"),
      parameter = "mu",
      data = df,
      prior_list = formula_result$prior_list
    ),
    "Random-effect block names in the supplied formula must be unique.",
    fixed = TRUE
  )
})

test_that("JAGS_evaluate_formula rejects malformed fitted random-effect metadata without structure", {

  df <- data.frame(id = factor(c("a", "b"), levels = c("a", "b")))
  formula_result <- JAGS_formula(
    formula = ~ 1 + diag(1 | id),
    parameter = "mu",
    data = df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(sd = prior("gamma", list(2, 2)))
    )
  )
  malformed_design <- formula_result$formula_design
  malformed_design$random_effects[[1]]$structure <- NULL
  malformed_design$random_effects[[1]]$covariance <- "diag"

  fit <- coda::mcmc(matrix(
    0,
    nrow = 1,
    dimnames = list(NULL, "mu_intercept")
  ))
  attr(fit, "formula_design") <- list(mu = malformed_design)

  expect_error(
    JAGS_evaluate_formula(
      fit = fit,
      formula = ~ 1 + diag(1 | id),
      parameter = "mu",
      data = df,
      prior_list = formula_result$prior_list
    ),
    "missing canonical 'random_term\\$structure'"
  )

  malformed_design <- formula_result$formula_design
  malformed_design$random_effects[[1]]$homogeneous_sd <- NULL
  attr(fit, "formula_design") <- list(mu = malformed_design)

  expect_error(
    JAGS_evaluate_formula(
      fit = fit,
      formula = ~ 1 + diag(1 | id),
      parameter = "mu",
      data = df,
      prior_list = formula_result$prior_list
    ),
    "missing canonical 'random_term\\$homogeneous_sd'"
  )
})

test_that("JAGS_formula scales predictors used only in random effects", {

  df <- data.frame(
    x = c(1, 2, 3, 4),
    id = factor(c("a", "a", "b", "b"), levels = c("a", "b"))
  )
  formula_result <- JAGS_formula(
    formula = ~ 1 + diag(0 + x | id),
    parameter = "mu",
    data = df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    formula_scale = list(x = TRUE),
    prior_random = prior_random(
      id = random_block(
        sd = prior("gamma", list(2, 2)),
        monitor = random_monitor(latent = TRUE, coefficients = FALSE, correlation = FALSE)
      )
    )
  )
  scale_info <- formula_result$formula_scale$mu_x
  expect_equal(scale_info$mean, mean(df$x))
  expect_equal(scale_info$sd, stats::sd(df$x))
  expect_equal(
    unname(formula_result$data$mu__xREx__id_xRE_DATAx[, "x"]),
    unname((df$x - scale_info$mean) / scale_info$sd)
  )
  scale_leaves <- attr(formula_result$formula_scale, "random_effect_sd_leaves")
  expect_s3_class(scale_leaves[["__xREx__id"]], "BayesTools_random_effect_sd_leaves")
  expect_equal(scale_leaves[["__xREx__id"]]$leaf_terms, c(mu__xREx__id_x = "x"))

  posterior <- matrix(
    c(0, 2, 1, -1),
    nrow = 1,
    dimnames = list(NULL, c(
      "mu_intercept",
      "mu__xREx__id_x",
      "mu__xREx__id_xRE_Zx[1,1]",
      "mu__xREx__id_xRE_Zx[2,1]"
    ))
  )
  fit <- coda::mcmc(posterior)
  attr(fit, "formula_design") <- list(mu = formula_result$formula_design)
  attr(fit, "formula_scale") <- list(mu = formula_result$formula_scale)

  newdata <- data.frame(
    x = c(2, 4),
    id = factor(c("a", "b"), levels = c("a", "b"))
  )
  prediction <- JAGS_evaluate_formula(
    fit = fit,
    formula = ~ 1 + diag(0 + x | id),
    parameter = "mu",
    data = newdata,
    prior_list = formula_result$prior_list
  )

  expected_x <- (newdata$x - scale_info$mean) / scale_info$sd
  expect_equal(unname(drop(prediction)), expected_x * c(2, -2))
})

test_that("random-effect grouping expressions use raw data when predictors are scaled", {

  df <- data.frame(
    x = c(-1, 0.1, 0.2, 100),
    id = factor(c("a", "a", "b", "b"))
  )
  formula_result <- JAGS_formula(
    formula = ~ 1 + x + random(1 | x > 0, name = "xpos", covariance = "diag"),
    parameter = "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    ),
    formula_scale = list(x = TRUE),
    prior_random = prior_random(
      xpos = random_block(
        sd = prior("gamma", list(2, 2)),
        monitor = random_monitor(latent = FALSE, coefficients = TRUE, correlation = FALSE)
      )
    )
  )

  raw_group_map <- as.numeric(factor(df$x > 0, levels = levels(as.factor(df$x > 0))))
  expect_equal(formula_result$formula_design$random_effects[[1]]$group_map, raw_group_map)

  coefficient_names <- BayesTools:::.bt_random_effect_coefficient_names(
    random_term = formula_result$formula_design$random_effects[[1]],
    n_groups = 2,
    n_columns = 1
  )
  posterior <- matrix(
    c(0, 0, 10, 20),
    nrow = 1,
    dimnames = list(NULL, c("mu_intercept", "mu_x", as.vector(coefficient_names)))
  )
  fit <- coda::mcmc(posterior)
  attr(fit, "formula_design") <- list(mu = formula_result$formula_design)
  attr(fit, "formula_scale") <- list(mu = formula_result$formula_scale)

  prediction <- JAGS_evaluate_formula(
    fit = fit,
    formula = ~ 1 + x + random(1 | x > 0, name = "xpos", covariance = "diag"),
    parameter = "mu",
    data = data.frame(x = 0.1),
    prior_list = formula_result$prior_list
  )
  expect_equal(unname(drop(prediction)), 20)
})

test_that("random-effect grouping factors preserve unused training levels", {

  df <- data.frame(
    id = factor(c("a", "a", "b", "b"), levels = c("a", "b", "unused"))
  )
  formula_result <- JAGS_formula(
    formula = ~ 1 + diag(1 | id),
    parameter = "mu",
    data = df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(
        sd = prior("gamma", list(2, 2)),
        monitor = random_monitor(latent = FALSE, coefficients = TRUE, correlation = FALSE)
      )
    )
  )
  random_term <- formula_result$formula_design$random_effects[[1]]
  expect_equal(random_term$group_levels, c("a", "b", "unused"))

  coefficient_names <- BayesTools:::.bt_random_effect_coefficient_names(
    random_term = random_term,
    n_groups = length(random_term$group_levels),
    n_columns = random_term$n_columns
  )
  posterior <- matrix(
    c(0, 10, 20, 30),
    nrow = 1,
    dimnames = list(NULL, c("mu_intercept", as.vector(coefficient_names)))
  )
  fit <- coda::mcmc(posterior)
  attr(fit, "formula_design") <- list(mu = formula_result$formula_design)

  prediction <- JAGS_evaluate_formula(
    fit = fit,
    formula = ~ 1 + diag(1 | id),
    parameter = "mu",
    data = data.frame(id = factor("unused", levels = c("a", "b", "unused"))),
    prior_list = formula_result$prior_list
  )
  expect_equal(unname(drop(prediction)), 30)
})

test_that("random effects evaluate on non-mu formula parameters", {

  df <- data.frame(
    id = factor(c("a", "a", "b", "b"))
  )
  formula_result <- JAGS_formula(
    formula = ~ 1 + diag(1 | id),
    parameter = "sigma",
    data = df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(
        sd = prior("gamma", list(2, 2)),
        monitor = random_monitor(latent = FALSE, coefficients = TRUE, correlation = FALSE)
      )
    )
  )
  expect_equal(
    names(formula_result$prior_list),
    c("sigma_intercept", "sigma__xREx__id_intercept")
  )

  random_term <- formula_result$formula_design$random_effects[[1]]
  coefficient_names <- BayesTools:::.bt_random_effect_coefficient_names(
    random_term = random_term,
    n_groups = length(random_term$group_levels),
    n_columns = random_term$n_columns
  )
  posterior <- matrix(
    c(1, 10, 20),
    nrow = 1,
    dimnames = list(NULL, c("sigma_intercept", as.vector(coefficient_names)))
  )
  fit <- coda::mcmc(posterior)
  attr(fit, "formula_design") <- list(sigma = formula_result$formula_design)

  prediction <- JAGS_evaluate_formula(
    fit = fit,
    formula = ~ 1 + diag(1 | id),
    parameter = "sigma",
    data = data.frame(id = factor(c("a", "b"))),
    prior_list = formula_result$prior_list
  )
  expect_equal(unname(drop(prediction)), c(11, 21))
})

test_that("character random-factor predictors keep fitted levels for prediction", {

  df <- data.frame(
    f = c("b", "a", "c", "a", "b", "c"),
    id = factor(c("g1", "g1", "g1", "g2", "g2", "g2"))
  )
  formula_result <- JAGS_formula(
    formula = ~ 1 + diag(0 + f | id),
    parameter = "mu",
    data = df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(
        sd = prior("gamma", list(2, 2)),
        monitor = random_monitor(latent = FALSE, coefficients = TRUE, correlation = FALSE)
      )
    )
  )
  random_term <- formula_result$formula_design$random_effects[[1]]
  expect_equal(random_term$xlevels$f, c("a", "b", "c"))

  coefficient_names <- BayesTools:::.bt_random_effect_coefficient_names(
    random_term = random_term,
    n_groups = length(random_term$group_levels),
    n_columns = random_term$n_columns
  )
  posterior <- matrix(
    c(0, seq_along(coefficient_names)),
    nrow = 1,
    dimnames = list(NULL, c("mu_intercept", as.vector(coefficient_names)))
  )
  fit <- coda::mcmc(posterior)
  attr(fit, "formula_design") <- list(mu = formula_result$formula_design)

  prediction <- JAGS_evaluate_formula(
    fit = fit,
    formula = ~ 1 + diag(0 + f | id),
    parameter = "mu",
    data = data.frame(f = c("c", "a"), id = factor(c("g1", "g2"))),
    prior_list = formula_result$prior_list
  )

  expected <- c(
    posterior[1, coefficient_names[1, 2]],
    0
  )
  expect_equal(unname(drop(prediction)), unname(expected))

  prediction_extra_level <- JAGS_evaluate_formula(
    fit = fit,
    formula = ~ 1 + diag(0 + f | id),
    parameter = "mu",
    data = data.frame(
      f = factor(c("c", "a"), levels = c("a", "b", "c", "unused")),
      id = factor(c("g1", "g2"))
    ),
    prior_list = formula_result$prior_list
  )
  expect_equal(unname(drop(prediction_extra_level)), unname(expected))
  expect_error(
    JAGS_evaluate_formula(
      fit = fit,
      formula = ~ 1 + diag(0 + f | id),
      parameter = "mu",
      data = data.frame(
        f = factor(c("unused", "a"), levels = c("a", "b", "c", "unused")),
        id = factor(c("g1", "g2"))
      ),
      prior_list = formula_result$prior_list
    ),
    "do not match the levels",
    fixed = TRUE
  )
})

test_that("transform_scale_samples unscales correlated random-effect SDs with covariance", {

  posterior <- matrix(
    c(1, 2, 1, 0.8, 0.8, 1, 1, 0.8, 0, 0.6),
    nrow = 1,
    dimnames = list(NULL, c(
      "mu__xREx__id_intercept",
      "mu__xREx__id_x",
      "mu__xREx__id_xRE_CORx_R[1,1]",
      "mu__xREx__id_xRE_CORx_R[1,2]",
      "mu__xREx__id_xRE_CORx_R[2,1]",
      "mu__xREx__id_xRE_CORx_R[2,2]",
      "mu__xREx__id_xRE_CORx_L[1,1]",
      "mu__xREx__id_xRE_CORx_L[2,1]",
      "mu__xREx__id_xRE_CORx_L[1,2]",
      "mu__xREx__id_xRE_CORx_L[2,2]"
    ))
  )
  formula_scale <- list(
    mu = list(mu_x = list(mean = 5, sd = 1))
  )
  attr(formula_scale$mu, "random_effect_terms") <- c(
    "mu__xREx__id_intercept" = "intercept",
    "mu__xREx__id_x" = "x"
  )

  transformed <- transform_scale_samples(posterior, formula_scale)

  expect_equal(
    unname(transformed[1, "mu__xREx__id_intercept"]),
    sqrt(1^2 + 5^2 * 2^2 - 2 * 5 * 0.8 * 1 * 2),
    tolerance = 1e-12
  )
  expect_equal(unname(transformed[1, "mu__xREx__id_x"]), 2, tolerance = 1e-12)
  expected_cor <- (0.8 * 1 * 2 - 5 * 2^2) /
    (sqrt(1^2 + 5^2 * 2^2 - 2 * 5 * 0.8 * 1 * 2) * 2)
  expect_equal(
    unname(transformed[1, "mu__xREx__id_xRE_CORx_R[1,2]"]),
    expected_cor,
    tolerance = 1e-12
  )
  expect_equal(
    unname(transformed[1, "mu__xREx__id_xRE_CORx_R[2,1]"]),
    expected_cor,
    tolerance = 1e-12
  )
  expected_L <- t(chol(matrix(c(1, expected_cor, expected_cor, 1), 2, 2)))
  expect_equal(
    unname(transformed[1, "mu__xREx__id_xRE_CORx_L[2,1]"]),
    expected_L[2, 1],
    tolerance = 1e-12
  )
  expect_equal(
    unname(transformed[1, "mu__xREx__id_xRE_CORx_L[2,2]"]),
    expected_L[2, 2],
    tolerance = 1e-12
  )
})

test_that("transform_scale_samples requires correlations for scaled correlated random-effect SDs", {

  df <- data.frame(
    x = c(1, 2, 3, 4),
    id = factor(c("a", "a", "b", "b"))
  )
  formula_result <- JAGS_formula(
    formula = ~ 1 + us(1 + x | id),
    parameter = "mu",
    data = df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    formula_scale = list(x = TRUE),
    prior_random = prior_random(
      id = random_block(
        sd = prior("gamma", list(2, 2)),
        cor = prior_lkj(eta = 1)
      )
    )
  )
  random_term <- formula_result$formula_design$random_effects[[1]]
  sd_cols <- random_term$sd_parameter_names
  posterior <- matrix(
    c(1, 2),
    nrow = 1,
    dimnames = list(NULL, sd_cols)
  )

  expect_error(
    transform_scale_samples(posterior, list(mu = formula_result$formula_scale)),
    "requires random-effect correlation samples",
    fixed = TRUE
  )

  R_names <- outer(
    seq_len(2),
    seq_len(2),
    Vectorize(function(row, column){
      paste0(random_term$parameter_stem, "_xRE_CORx_R[", row, ",", column, "]")
    })
  )
  partial_R_posterior <- cbind(posterior, 1)
  colnames(partial_R_posterior)[ncol(partial_R_posterior)] <- R_names[1, 1]
  expect_error(
    transform_scale_samples(partial_R_posterior, list(mu = formula_result$formula_scale)),
    "correlation samples are incomplete",
    fixed = TRUE
  )

  expect_error(
    BayesTools:::.bt_random_effect_summary_samples(
      model_samples = posterior,
      prior_list = formula_result$prior_list,
      formula_design = list(mu = formula_result$formula_design),
      mode = "standard",
      formula_scale = list(mu = formula_result$formula_scale)
    ),
    "missing or invalid canonical correlation coordinates",
    fixed = TRUE
  )
})

test_that("JAGS_estimates_table backtransforms random-effect correlations", {

  testthat::skip_if_not_installed("runjags")

  df <- data.frame(
    x = c(4, 5, 6),
    id = factor(c("a", "a", "b"))
  )
  fixed_priors <- list(
    intercept = prior("normal", list(0, 1)),
    x = prior("normal", list(0, 1))
  )
  formula_result <- JAGS_formula(
    formula = ~ 1 + x + random(1 + x | id, name = "id", covariance = "us"),
    parameter = "mu",
    data = df,
    prior_list = fixed_priors,
    formula_scale = list(x = TRUE),
    prior_random = prior_random(
      id = random_block(
        sd = prior("normal", list(0, 1), truncation = list(0, Inf)),
        cor = prior_lkj(eta = 1),
        monitor = random_monitor(
          latent = FALSE,
          coefficients = FALSE,
          correlation = TRUE
        )
      )
    )
  )

  random_term <- formula_result$formula_design$random_effects[[1]]
  sd_names <- random_term$sd_parameter_names
  R_names <- outer(
    seq_len(2),
    seq_len(2),
    Vectorize(function(row, column) {
      paste0(random_term$parameter_stem, "_xRE_CORx_R[", row, ",", column, "]")
    })
  )
  L_names <- BayesTools:::.bt_random_effect_cholesky_names(random_term, 2L)
  source_sd <- c(1, 2)
  source_cor <- matrix(c(1, 0.8, 0.8, 1), 2, 2)
  source_L <- t(chol(source_cor))
  posterior <- matrix(
    rep(c(source_sd, as.vector(source_cor), as.vector(source_L)), 2),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(NULL, c(sd_names, as.vector(R_names), as.vector(L_names)))
  )

  fit <- list(
    mcmc = coda::mcmc.list(coda::mcmc(posterior)),
    summary.pars = list(mutate = NULL),
    monitor = colnames(posterior),
    sample = nrow(posterior)
  )
  class(fit) <- c("runjags", "BayesTools_fit")
  attr(fit, "prior_list") <- formula_result$prior_list
  attr(fit, "formula_design") <- list(mu = formula_result$formula_design)
  attr(fit, "formula_scale") <- list(mu = formula_result$formula_scale)

  scaled_samples <- JAGS_estimates_table(
    fit,
    transform_scaled = FALSE,
    random_effects_summary = "standard",
    remove_diagnostics = TRUE,
    return_samples = TRUE
  )
  original_samples <- JAGS_estimates_table(
    fit,
    transform_scaled = TRUE,
    random_effects_summary = "standard",
    remove_diagnostics = TRUE,
    return_samples = TRUE
  )
  original_table <- JAGS_estimates_table(
    fit,
    transform_scaled = TRUE,
    random_effects_summary = "standard",
    remove_diagnostics = TRUE
  )

  scale_info <- formula_result$formula_scale$mu_x
  M <- matrix(
    c(1, -scale_info$mean / scale_info$sd, 0, 1 / scale_info$sd),
    nrow = 2,
    byrow = TRUE
  )
  source_cov <- diag(source_sd, nrow = 2) %*% source_cor %*%
    diag(source_sd, nrow = 2)
  expected_cov <- M %*% source_cov %*% t(M)
  expected_sd <- sqrt(diag(expected_cov))
  expected_cor <- expected_cov[1, 2] / prod(expected_sd)

  expect_equal(
    unname(scaled_samples[, "(mu) cor(intercept,x | id)"]),
    rep(0.8, 2),
    tolerance = 1e-12
  )
  expect_equal(
    unname(original_samples[, "(mu) sd(intercept | id)"]),
    rep(expected_sd[1], 2),
    tolerance = 1e-12
  )
  expect_equal(
    unname(original_samples[, "(mu) sd(x | id)"]),
    rep(expected_sd[2], 2),
    tolerance = 1e-12
  )
  expect_equal(
    unname(original_samples[, "(mu) cor(intercept,x | id)"]),
    rep(expected_cor, 2),
    tolerance = 1e-12
  )
  expect_equal(
    unname(original_table["(mu) cor(intercept,x | id)", "Mean"]),
    expected_cor,
    tolerance = 1e-12
  )
  expect_lt(expected_cor, 0)
})

test_that("transform_scale_samples updates valid random-effect correlations draw-wise", {

  posterior <- matrix(
    c(
      1, 2, 1, 0.8, 0.8, 1, 1, 0.8, 0, 0.6,
      1, 2, 1, 0.8, 0.8, 1, 1, 0.8, 0, 0.6
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(NULL, c(
      "mu__xREx__id_intercept",
      "mu__xREx__id_x",
      "mu__xREx__id_xRE_CORx_R[1,1]",
      "mu__xREx__id_xRE_CORx_R[1,2]",
      "mu__xREx__id_xRE_CORx_R[2,1]",
      "mu__xREx__id_xRE_CORx_R[2,2]",
      "mu__xREx__id_xRE_CORx_L[1,1]",
      "mu__xREx__id_xRE_CORx_L[2,1]",
      "mu__xREx__id_xRE_CORx_L[1,2]",
      "mu__xREx__id_xRE_CORx_L[2,2]"
    ))
  )
  formula_scale <- list(
    mu = list(mu_x = list(mean = 5, sd = 1))
  )
  attr(formula_scale$mu, "random_effect_terms") <- c(
    "mu__xREx__id_intercept" = "intercept",
    "mu__xREx__id_x" = "x"
  )

  transformed <- transform_scale_samples(posterior, formula_scale)

  expected_cor <- (0.8 * 1 * 2 - 5 * 2^2) /
    (sqrt(1^2 + 5^2 * 2^2 - 2 * 5 * 0.8 * 1 * 2) * 2)
  expected_L <- t(chol(matrix(c(1, expected_cor, expected_cor, 1), 2, 2)))
  expect_equal(
    unname(transformed[1, "mu__xREx__id_xRE_CORx_R[1,2]"]),
    expected_cor,
    tolerance = 1e-12
  )
  expect_equal(
    unname(transformed[2, "mu__xREx__id_xRE_CORx_R[1,2]"]),
    expected_cor,
    tolerance = 1e-12
  )
  expect_equal(
    unname(transformed[2, "mu__xREx__id_xRE_CORx_L[2,1]"]),
    expected_L[2, 1],
    tolerance = 1e-12
  )
})

test_that("transform_scale_samples clears invalid transformed random-effect correlations", {

  posterior <- matrix(
    c(
      0, 0, 1, 0.8, 0.8, 1, 1, 0.8, 0, 0.6,
      1, 2, 1, 0.8, 0.8, 1, 1, 0.8, 0, 0.6
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(NULL, c(
      "mu__xREx__id_intercept",
      "mu__xREx__id_x",
      "mu__xREx__id_xRE_CORx_R[1,1]",
      "mu__xREx__id_xRE_CORx_R[1,2]",
      "mu__xREx__id_xRE_CORx_R[2,1]",
      "mu__xREx__id_xRE_CORx_R[2,2]",
      "mu__xREx__id_xRE_CORx_L[1,1]",
      "mu__xREx__id_xRE_CORx_L[2,1]",
      "mu__xREx__id_xRE_CORx_L[1,2]",
      "mu__xREx__id_xRE_CORx_L[2,2]"
    ))
  )
  formula_scale <- list(
    mu = list(mu_x = list(mean = 5, sd = 1))
  )
  attr(formula_scale$mu, "random_effect_terms") <- c(
    "mu__xREx__id_intercept" = "intercept",
    "mu__xREx__id_x" = "x"
  )

  transformed <- transform_scale_samples(posterior, formula_scale)

  expect_true(is.na(transformed[1, "mu__xREx__id_xRE_CORx_R[1,2]"]))
  expect_true(is.na(transformed[1, "mu__xREx__id_xRE_CORx_L[2,1]"]))
  expect_false(is.na(transformed[2, "mu__xREx__id_xRE_CORx_R[1,2]"]))

  df <- data.frame(
    x = c(1, 2, 3, 4),
    id = factor(c("a", "a", "b", "b"))
  )
  formula_result <- JAGS_formula(
    formula = ~ 1 + us(1 + x | id),
    parameter = "mu",
    data = df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    formula_scale = list(x = TRUE),
    prior_random = prior_random(
      id = random_block(
        sd = prior("gamma", list(2, 2)),
        cor = prior_lkj(eta = 1)
      )
    )
  )
  transformed_summary <- BayesTools:::.bt_random_effect_summary_samples(
    model_samples = transformed,
    prior_list = formula_result$prior_list,
    formula_design = list(mu = formula_result$formula_design),
    mode = "standard"
  )
  summary_name <- "mu__xRE_SUMMARY__cor__id__intercept__x"
  expect_true(summary_name %in% colnames(transformed_summary$model_samples))
  expect_true(is.na(transformed_summary$model_samples[1, summary_name]))
  expect_false(is.na(transformed_summary$model_samples[2, summary_name]))
})

test_that("transform_scale_samples keeps indexed random-effect correlations in one block", {

  df <- data.frame(
    x = c(1, 2, 3, 4, 5, 6),
    f = factor(c("a", "b", "a", "b", "a", "b")),
    id = factor(c("g1", "g1", "g2", "g2", "g3", "g3"))
  )
  formula_result <- JAGS_formula(
    formula = ~ 1 + us(0 + f + x:f | id),
    parameter = "mu",
    data = df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    formula_scale = list(x = TRUE),
    prior_random = prior_random(
      id = random_block(
        sd = prior("gamma", list(2, 2)),
        cor = prior_lkj(eta = 1)
      )
    )
  )

  random_term <- formula_result$formula_design$random_effects[[1]]
  sd_cols <- random_term$sd_parameter_names
  expect_equal(
    unname(random_term$sd_leaves$leaf_terms[sd_cols]),
    c("f[b]", "f[a]__xXx__x", "f[b]__xXx__x")
  )

  source_sd <- c(1.1, 0.8, 1.2)
  source_cor <- matrix(0.25, nrow = 3, ncol = 3)
  diag(source_cor) <- 1
  cor_names <- outer(
    seq_along(sd_cols),
    seq_along(sd_cols),
    Vectorize(function(row, column){
      paste0("mu__xREx__id_xRE_CORx_R[", row, ",", column, "]")
    })
  )
  posterior <- matrix(
    c(source_sd, as.vector(source_cor)),
    nrow = 1,
    dimnames = list(NULL, c(sd_cols, as.vector(cor_names)))
  )

  transformed <- transform_scale_samples(
    posterior,
    list(mu = formula_result$formula_scale)
  )

  scale_info <- formula_result$formula_scale$mu_x
  M <- diag(3)
  M[1, 3] <- -scale_info$mean / scale_info$sd
  M[2, 2] <- 1 / scale_info$sd
  M[3, 3] <- 1 / scale_info$sd
  expected_cov <- M %*% diag(source_sd, nrow = 3) %*% source_cor %*%
    diag(source_sd, nrow = 3) %*% t(M)
  expected_sd <- sqrt(diag(expected_cov))
  expected_cor <- expected_cov / tcrossprod(expected_sd)

  expect_equal(unname(transformed[1, sd_cols]), expected_sd, tolerance = 1e-12)
  expect_equal(
    unname(transformed[1, as.vector(cor_names)]),
    as.vector(expected_cor),
    tolerance = 1e-12
  )
})

test_that("transform_scale_samples guards homogeneous random-effect SD scaling", {

  df <- data.frame(
    x = c(1, 2, 3, 4),
    id = factor(c("a", "a", "b", "b"))
  )
  slope_result <- JAGS_formula(
    formula = ~ 1 + id(0 + x | id),
    parameter = "mu",
    data = df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    formula_scale = list(x = TRUE),
    prior_random = prior_random(id = random_block(sd = prior("gamma", list(2, 2))))
  )
  slope_posterior <- matrix(
    2,
    nrow = 1,
    dimnames = list(NULL, "mu__xREx__id_sd")
  )
  slope_transformed <- transform_scale_samples(
    slope_posterior,
    list(mu = slope_result$formula_scale)
  )
  expect_equal(
    unname(slope_transformed[1, "mu__xREx__id_sd"]),
    2 / slope_result$formula_scale$mu_x$sd,
    tolerance = 1e-12
  )

  intercept_slope_result <- JAGS_formula(
    formula = ~ 1 + id(1 + x | id),
    parameter = "mu",
    data = df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    formula_scale = list(x = TRUE),
    prior_random = prior_random(id = random_block(sd = prior("gamma", list(2, 2))))
  )
  expect_error(
    transform_scale_samples(
      slope_posterior,
      list(mu = intercept_slope_result$formula_scale)
    ),
    "Cannot unscale homogeneous random-effect SD",
    fixed = TRUE
  )
})

test_that("JAGS_evaluate_formula sums multiple monitored random-effect blocks", {

  df <- data.frame(
    x = c(0, 1, 2, 3),
    id = factor(c("a", "a", "b", "b"), levels = c("a", "b")),
    drug = factor(c("d1", "d2", "d1", "d2"), levels = c("d1", "d2"))
  )
  formula_result <- JAGS_formula(
    formula = ~ 1 + x +
      random(1 | id, name = "study", covariance = "diag") +
      random(0 + x | drug, name = "drug_slope", covariance = "diag"),
    parameter = "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    ),
    prior_random = prior_random(
      study = random_block(sd = prior("gamma", list(2, 2))),
      drug_slope = random_block(sd = prior("gamma", list(2, 2)))
    )
  )
  posterior <- matrix(
    c(
      1, 2,
      0.5, -0.5,
      0.1, 0.2
    ),
    nrow = 1,
    dimnames = list(NULL, c(
      "mu_intercept",
      "mu_x",
      "mu__xREx__study_xRE_COEFx[1,1]",
      "mu__xREx__study_xRE_COEFx[2,1]",
      "mu__xREx__drug_slope_xRE_COEFx[1,1]",
      "mu__xREx__drug_slope_xRE_COEFx[2,1]"
    ))
  )
  fit <- coda::mcmc(posterior)
  attr(fit, "formula_design") <- list(mu = formula_result$formula_design)

  prediction <- JAGS_evaluate_formula(
    fit = fit,
    formula = ~ 1 + x +
      random(1 | id, name = "study", covariance = "diag") +
      random(0 + x | drug, name = "drug_slope", covariance = "diag"),
    parameter = "mu",
    data = df,
    prior_list = formula_result$prior_list
  )

  expect_equal(unname(drop(prediction)), c(1.5, 3.7, 4.7, 7.1))
})

test_that("JAGS_evaluate_formula reconstructs factor random-slope designs", {

  df <- data.frame(
    f = factor(c("a", "b", "c", "a"), levels = c("a", "b", "c")),
    id = factor(c("g1", "g1", "g2", "g2"), levels = c("g1", "g2"))
  )
  formula_result <- JAGS_formula(
    formula = ~ 1 + diag(0 + f | id),
    parameter = "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1))
    ),
    prior_random = prior_random(
      id = random_block(
        sd = prior_factor("gamma", list(2, 2), contrast = "treatment"),
        monitor = random_monitor(latent = FALSE, coefficients = TRUE, correlation = FALSE)
      )
    )
  )
  posterior <- matrix(
    c(1, 10, 20, 30, 40),
    nrow = 1,
    dimnames = list(NULL, c(
      "mu_intercept",
      "mu__xREx__id_xRE_COEFx[1,1]",
      "mu__xREx__id_xRE_COEFx[1,2]",
      "mu__xREx__id_xRE_COEFx[2,1]",
      "mu__xREx__id_xRE_COEFx[2,2]"
    ))
  )
  fit <- coda::mcmc(posterior)
  attr(fit, "formula_design") <- list(mu = formula_result$formula_design)

  prediction <- JAGS_evaluate_formula(
    fit = fit,
    formula = ~ 1 + diag(0 + f | id),
    parameter = "mu",
    data = data.frame(
      f = c("c", "b"),
      id = factor(c("g1", "g2"), levels = c("g2", "g1"))
    ),
    prior_list = formula_result$prior_list
  )

  expect_equal(unname(drop(prediction)), c(21, 31))
})

test_that("prior_random maps to explicitly named random-effect blocks", {

  df <- data.frame(
    x = c(-1, 0, 1, 2, -2, 3),
    id = factor(c("b", "a", "b", "c", "a", "c"), levels = c("c", "a", "b")),
    drug = factor(c("d1", "d1", "d2", "d2", "d3", "d3"))
  )
  sd_prior <- prior("normal", list(0, 1), truncation = list(lower = 0, upper = Inf))

  result <- JAGS_formula(
    formula = ~ 1 + x +
      random(1 + x | id, name = "study") +
      diag(0 + x | drug, name = "drug_slope"),
    parameter = "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    ),
    prior_random = prior_random(
      study = random_block(sd = sd_prior, cor = prior_lkj(eta = 2, backend = "syntax")),
      drug_slope = random_block(sd = sd_prior)
    )
  )

  expect_equal(
    names(result$prior_list),
    c(
      "mu_intercept",
      "mu_x",
      "mu__xREx__study_intercept",
      "mu__xREx__study_x",
      "mu__xREx__drug_slope_x"
    )
  )
  expect_equal(
    names(result$formula_design$jags_data_names),
    c("x", "__xREx__study", "__xREx__drug_slope")
  )
  expect_equal(result$formula_design$random_effects[[1]]$block_name, "study")
  expect_equal(result$formula_design$random_effects[[1]]$structure, "us")
  expect_equal(result$formula_design$random_effects[[2]]$block_name, "drug_slope")
  expect_equal(result$formula_design$random_effects[[2]]$structure, "diag")
  expect_equal(result$jags_modules, character())
  expect_true(any(grepl("dbeta", result$formula_syntax, fixed = TRUE)))

  inherited_covariance <- JAGS_formula(
    formula = ~ 1 + x + random(1 + x | id, name = "study"),
    parameter = "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    ),
    prior_random = prior_random(
      cor = prior_lkj(eta = 3, backend = "syntax"),
      study = random_block(sd = sd_prior)
    )
  )
  expect_equal(inherited_covariance$jags_modules, character())
  expect_match(inherited_covariance$formula_syntax, "dbeta(3, 3)", fixed = TRUE)

  backend_override <- JAGS_formula(
    formula = ~ 1 + x + random(1 + x | id, name = "study"),
    parameter = "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    ),
    prior_random = prior_random(
      cor = prior_lkj(eta = 3, backend = "module"),
      study = random_block(
        sd = sd_prior,
        covariance = random_covariance(backend = "syntax")
      )
    )
  )
  expect_equal(backend_override$jags_modules, character())
  expect_match(backend_override$formula_syntax, "dbeta(3, 3)", fixed = TRUE)

  direct_backend_override <- JAGS_formula(
    formula = ~ 1 + x + random(1 + x | id, name = "study"),
    parameter = "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    ),
    prior_random = prior_random(
      covariance = random_covariance(
        cor = prior_lkj(eta = 3),
        backend = "syntax"
      ),
      study = random_block(sd = sd_prior)
    )
  )
  expect_equal(direct_backend_override$jags_modules, character())
  expect_match(direct_backend_override$formula_syntax, "dbeta(3, 3)", fixed = TRUE)

  inherited_backend_block_cor <- JAGS_formula(
    formula = ~ 1 + x + random(1 + x | id, name = "study"),
    parameter = "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    ),
    prior_random = prior_random(
      covariance = random_covariance(backend = "syntax"),
      study = random_block(sd = sd_prior, cor = prior_lkj(eta = 3))
    )
  )
  expect_equal(inherited_backend_block_cor$jags_modules, character())
  expect_match(inherited_backend_block_cor$formula_syntax, "dbeta(3, 3)", fixed = TRUE)

  explicit_block_backend <- JAGS_formula(
    formula = ~ 1 + x + random(1 + x | id, name = "study"),
    parameter = "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    ),
    prior_random = prior_random(
      covariance = random_covariance(backend = "syntax"),
      study = random_block(
        sd = sd_prior,
        cor = prior_lkj(eta = 3, backend = "module")
      )
    )
  )
  expect_equal(explicit_block_backend$jags_modules, "BayesTools")
  expect_false(grepl("dbeta(3, 3)", explicit_block_backend$formula_syntax, fixed = TRUE))

  mixed_covariance <- JAGS_formula(
    formula = ~ 1 + x +
      random(1 + x | id, name = "study") +
      random(0 + x | drug, name = "drug_slope", covariance = "diag"),
    parameter = "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    ),
    prior_random = prior_random(
      sd = sd_prior,
      cor = prior_lkj(eta = 3, backend = "syntax"),
      drug_slope = random_block(covariance = random_covariance(structure = "diag"))
    )
  )
  expect_equal(mixed_covariance$formula_design$random_effects[[1]]$structure, "us")
  expect_equal(mixed_covariance$formula_design$random_effects[[2]]$structure, "diag")
  expect_match(mixed_covariance$formula_syntax, "dbeta(3, 3)", fixed = TRUE)
  expect_false(any(grepl("drug_slope_xRE_CORx", mixed_covariance$formula_syntax, fixed = TRUE)))

  wrapper_covariance <- JAGS_formula(
    formula = ~ 1 + x + random(0 + x | id, name = "id_slope", covariance = "diag"),
    parameter = "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    ),
    prior_random = prior_random(
      id_slope = random_block(sd = sd_prior)
    )
  )
  expect_equal(wrapper_covariance$formula_design$random_effects[[1]]$block_name, "id_slope")
  expect_equal(wrapper_covariance$formula_design$random_effects[[1]]$structure, "diag")
  expect_equal(names(wrapper_covariance$prior_list), c("mu_intercept", "mu_x", "mu__xREx__id_slope_x"))

  unnamed_blocks <- JAGS_formula(
    formula = ~ 1 + x + (1 | id) + (0 + x | drug),
    parameter = "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    ),
    prior_random = prior_random(sd = sd_prior)
  )
  expect_equal(
    vapply(unnamed_blocks$formula_design$random_effects, `[[`, character(1), "block_name"),
    c("id", "drug")
  )
  expect_equal(
    names(unnamed_blocks$prior_list),
    c("mu_intercept", "mu_x", "mu__xREx__id_intercept", "mu__xREx__drug_x")
  )
  expect_error(
    JAGS_formula(
      formula = ~ 1 + x +
        random(1 | id, name = "study") +
        random(0 + x | drug, name = "study"),
      parameter = "mu",
      data = df,
      prior_list = list(
        intercept = prior("normal", list(0, 1)),
        x = prior("normal", list(0, 1))
      ),
      prior_random = prior_random(study = random_block(sd = sd_prior))
    ),
    "Random-effect block names must be unique",
    fixed = TRUE
  )
  expect_error(
    JAGS_formula(
      formula = ~ 1 + x +
        random(1 | id, name = "study-id") +
        random(0 + x | drug, name = "study_id"),
      parameter = "mu",
      data = df,
      prior_list = list(
        intercept = prior("normal", list(0, 1)),
        x = prior("normal", list(0, 1))
      ),
      prior_random = prior_random(study_id = random_block(sd = sd_prior))
    ),
    "Random-effect block names must be unique",
    fixed = TRUE
  )
  expect_error(
    JAGS_formula(
      formula = ~ 1 + x + random(1 | id, name = "!!!"),
      parameter = "mu",
      data = df,
      prior_list = list(
        intercept = prior("normal", list(0, 1)),
        x = prior("normal", list(0, 1))
      ),
      prior_random = prior_random(sd = sd_prior)
    ),
    "must contain at least one letter",
    fixed = TRUE
  )
})

test_that("nested grouping formulas expand to evaluable random-effect blocks", {

  df <- data.frame(
    x = 1:6,
    study = factor(c("s1", "s1", "s1", "s2", "s2", "s2")),
    paper = factor(c("p1", "p1", "p2", "p1", "p2", "p2"))
  )

  result <- JAGS_formula(
    formula = ~ x +
      random(1 | paper:study, name = "paper_study", covariance = "diag") +
      random(1 | study, name = "study", covariance = "diag"),
    parameter = "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    ),
    prior_random = prior_random(
      paper_study = random_block(sd = prior("gamma", list(2, 2))),
      study = random_block(sd = prior("gamma", list(2, 2)))
    )
  )

  expect_equal(
    names(result$formula_design$jags_data_names),
    c("x", "__xREx__paper_study", "__xREx__study")
  )
  expect_equal(
    names(result$prior_list),
    c(
      "mu_intercept",
      "mu_x",
      "mu__xREx__paper_study_intercept",
      "mu__xREx__study_intercept"
    )
  )
  expect_equal(result$formula_design$random_effects[[1]]$group_label, "paper:study")
  expect_equal(result$formula_design$random_effects[[1]]$block_name, "paper_study")
  expect_equal(result$formula_design$random_effects[[2]]$group_label, "study")
  expect_equal(result$formula_design$random_effects[[1]]$n_groups, 4L)
  expect_equal(result$formula_design$random_effects[[2]]$n_groups, 2L)

  char_df <- transform(
    df,
    study = as.character(study),
    paper = as.character(paper)
  )
  char_result <- JAGS_formula(
    formula = ~ x +
      random(1 | paper:study, name = "paper_study", covariance = "diag") +
      random(1 | study, name = "study", covariance = "diag"),
    parameter = "mu",
    data = char_df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    ),
    prior_random = prior_random(
      paper_study = random_block(sd = prior("gamma", list(2, 2))),
      study = random_block(sd = prior("gamma", list(2, 2)))
    )
  )
  expect_equal(char_result$formula_design$random_effects[[1]]$group_label, "paper:study")
  expect_equal(char_result$formula_design$random_effects[[1]]$n_groups, 4L)
  expect_equal(char_result$formula_design$random_effects[[2]]$n_groups, 2L)
})

test_that("diag random-effect syntax handles slope-only and factor-slope designs", {

  df <- data.frame(
    x = c(-1, 0, 1, 2, -2, 3),
    id = factor(c("b", "a", "b", "c", "a", "c"), levels = c("c", "a", "b"))
  )

  slope_result <- JAGS_formula(
    formula = ~ 1 + x + diag(0 + x | id),
    parameter = "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    ),
    prior_random = prior_random(
      id = random_block(sd = prior("gamma", list(2, 2)))
    )
  )

  expect_equal(.remove_random_effects(~ z + diag(0 + x | id)), ~ z, ignore_formula_env = TRUE)
  expect_equal(slope_result$formula_design$random_effects[[1]]$term_formula, ~ 0 + x, ignore_formula_env = TRUE)
  expect_equal(slope_result$formula_design$random_effects[[1]]$structure, "diag")
  expect_true(isTRUE(slope_result$formula_design$random_effects[[1]]$independent))
  expect_equal(slope_result$formula_design$random_effects[[1]]$n_columns, 1L)
  expect_equal(colnames(slope_result$data$mu__xREx__id_xRE_DATAx), "x")
  expect_equal(unname(slope_result$data$mu__xREx__id_xRE_DATAx[, "x"]), df$x)
  expect_equal(
    names(slope_result$prior_list),
    c("mu_intercept", "mu_x", "mu__xREx__id_x")
  )

  factor_df <- data.frame(
    f = factor(rep(c("a", "b", "c"), 2), levels = c("a", "b", "c")),
    id = factor(rep(c("g1", "g2", "g3"), each = 2))
  )

  factor_result <- JAGS_formula(
    formula = ~ 1 + f + diag(0 + f | id),
    parameter = "mu",
    data = factor_df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      f = prior_factor("normal", list(0, 1), contrast = "treatment")
    ),
    prior_random = prior_random(
      id = random_block(sd = prior("gamma", list(2, 2)))
    )
  )

  expected_factor_df <- factor_df
  stats::contrasts(expected_factor_df$f) <- "contr.treatment"
  expected_factor_design <- stats::model.matrix(~ f, data = expected_factor_df)[, -1, drop = FALSE]

  expect_equal(colnames(factor_result$data$mu__xREx__id_xRE_DATAx), c("fb", "fc"))
  expect_equal(
    as.vector(factor_result$data$mu__xREx__id_xRE_DATAx),
    as.vector(expected_factor_design)
  )
  expect_equal(factor_result$formula_design$random_effects[[1]]$model_terms_type[["f"]], "factor")
  expect_equal(factor_result$formula_design$random_effects[[1]]$n_columns, 2L)
  expect_equal(
    names(factor_result$prior_list),
    c("mu_intercept", "mu_f", "mu__xREx__id_f")
  )
  expect_true(is.prior.treatment(factor_result$prior_list$mu__xREx__id_f))
  expect_equal(attr(factor_result$prior_list$mu__xREx__id_f, "random_factor"), "id")

  ordered_df <- factor_df
  ordered_df$f <- ordered(ordered_df$f, levels = c("a", "b", "c"))
  old_contrasts <- options("contrasts")
  on.exit(options(old_contrasts), add = TRUE)
  options(contrasts = c("contr.sum", "contr.poly"))
  ordered_result <- JAGS_formula(
    formula = ~ 1 + diag(0 + f | id),
    parameter = "mu",
    data = ordered_df,
    prior_list = list(
      intercept = prior("normal", list(0, 1))
    ),
    prior_random = prior_random(
      id = random_block(sd = prior("gamma", list(2, 2)))
    )
  )
  expected_ordered_df <- ordered_df
  stats::contrasts(expected_ordered_df$f) <- "contr.treatment"
  expected_ordered <- stats::model.matrix(~ f, data = expected_ordered_df)[, -1, drop = FALSE]
  expect_equal(colnames(ordered_result$data$mu__xREx__id_xRE_DATAx), c("fb", "fc"))
  expect_equal(
    as.vector(ordered_result$data$mu__xREx__id_xRE_DATAx),
    as.vector(expected_ordered)
  )
  expect_true(is.prior.treatment(ordered_result$prior_list$mu__xREx__id_f))

  random_only_orthonormal <- suppressWarnings(
    JAGS_formula(
      formula = ~ 1 + diag(0 + f | id),
      parameter = "mu",
      data = factor_df,
      prior_list = list(
        intercept = prior("normal", list(0, 1))
      ),
      prior_random = prior_random(
        id = random_block(sd = prior_factor("mnormal", list(0, 1), contrast = "orthonormal"))
      )
    )
  )
  orthonormal_df <- factor_df
  stats::contrasts(orthonormal_df$f) <- "contr.orthonormal"
  expected_orthonormal <- stats::model.matrix(~ f, data = orthonormal_df)[, -1, drop = FALSE]

  expect_equal(
    as.vector(random_only_orthonormal$data$mu__xREx__id_xRE_DATAx),
    as.vector(expected_orthonormal),
    tolerance = 1e-12
  )
  expect_true(is.prior.orthonormal(random_only_orthonormal$prior_list$mu__xREx__id_f))

  random_only_meandif <- suppressWarnings(
    JAGS_formula(
      formula = ~ 1 + diag(0 + f | id),
      parameter = "mu",
      data = factor_df,
      prior_list = list(
        intercept = prior("normal", list(0, 1))
      ),
      prior_random = prior_random(
        id = random_block(sd = prior_factor("mnormal", list(0, 1), contrast = "meandif"))
      )
    )
  )
  meandif_df <- factor_df
  stats::contrasts(meandif_df$f) <- "contr.meandif"
  expected_meandif <- stats::model.matrix(~ f, data = meandif_df)[, -1, drop = FALSE]

  expect_equal(
    as.vector(random_only_meandif$data$mu__xREx__id_xRE_DATAx),
    as.vector(expected_meandif),
    tolerance = 1e-12
  )
  expect_true(is.prior.meandif(random_only_meandif$prior_list$mu__xREx__id_f))

  interaction_df <- data.frame(
    f = factor(rep(c("a", "b", "c"), 4), levels = c("a", "b", "c")),
    g = factor(rep(c("u", "v"), each = 6), levels = c("u", "v")),
    id = factor(rep(c("s1", "s2", "s3", "s4"), each = 3))
  )
  interaction_result <- JAGS_formula(
    formula = ~ 1 + f * g + diag(0 + f:g | id),
    parameter = "mu",
    data = interaction_df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      f = prior_factor("normal", list(0, 1), contrast = "treatment"),
      g = prior_factor("normal", list(0, 1), contrast = "treatment"),
      "f:g" = prior_factor("normal", list(0, 1), contrast = "treatment")
    ),
    prior_random = prior_random(
      id = random_block(sd = prior("gamma", list(2, 2)))
    )
  )
  expected_interaction_df <- interaction_df
  stats::contrasts(expected_interaction_df$f) <- "contr.treatment"
  stats::contrasts(expected_interaction_df$g) <- "contr.treatment"
  expected_interaction <- stats::model.matrix(~ f:g, data = expected_interaction_df)[, -1, drop = FALSE]
  colnames(expected_interaction) <- gsub(":", "__xXx__", colnames(expected_interaction), fixed = TRUE)

  expect_equal(
    colnames(interaction_result$data$mu__xREx__id_xRE_DATAx),
    colnames(expected_interaction)
  )
  expect_equal(
    as.vector(interaction_result$data$mu__xREx__id_xRE_DATAx),
    as.vector(expected_interaction),
    tolerance = 1e-12
  )
  expect_equal(
    attr(interaction_result$prior_list$mu__xREx__id_f__xXx__g, "level_names"),
    list(f = c("a", "b", "c"), g = c("u", "v"))
  )

  mixed_contrast_df <- interaction_df
  stats::contrasts(mixed_contrast_df$f) <- "contr.treatment"
  stats::contrasts(mixed_contrast_df$g) <- "contr.orthonormal"
  expect_error(
    JAGS_formula(
      formula = ~ 1 + diag(0 + f:g | id),
      parameter = "mu",
      data = mixed_contrast_df,
      prior_list = list(
        intercept = prior("normal", list(0, 1))
      ),
      prior_random = prior_random(
        id = random_block(sd = prior("gamma", list(2, 2)))
      )
    ),
    "mixed factor contrast families",
    fixed = TRUE
  )
})

test_that("structured random-effect terms use level-indexed factor columns and scalar rho priors", {

  factor_df <- data.frame(
    f = factor(rep(c("a", "b", "c"), 4), levels = c("a", "b", "c")),
    g = factor(rep(c("u", "v"), 6), levels = c("u", "v")),
    id = factor(rep(c("g1", "g2", "g3", "g4"), each = 3))
  )
  continuous_df <- data.frame(
    x = c(-1, 0, 1, 2, -2, 0.5),
    id = factor(c("a", "a", "b", "b", "c", "c"), levels = c("a", "b", "c"))
  )
  car_df <- data.frame(
    time = c(0, 0.5, 2, 0, 0.5, 2),
    id = factor(c("a", "a", "a", "b", "b", "b"), levels = c("a", "b"))
  )
  sd_prior <- prior("normal", list(0, 1), truncation = list(lower = 0, upper = Inf))

  expect_error(
    JAGS_formula(
      formula = ~ 1 + cs(f | id),
      parameter = "mu",
      data = factor_df,
      prior_list = list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(id = random_block(sd = sd_prior))
    ),
    "requires a scalar correlation prior",
    fixed = TRUE
  )

  default_cs <- JAGS_formula(
    formula = ~ 1 + cs(f | id),
    parameter = "mu",
    data = factor_df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(sd = sd_prior, rho = prior("normal", list(0, 0.5)))
    )
  )

  expect_equal(names(default_cs$prior_list), c("mu_intercept", "mu__xREx__id_sd", "mu__xREx__id_rho_z"))
  expect_equal(default_cs$prior_list$mu__xREx__id_rho_z$distribution, "normal")
  expect_equal(default_cs$prior_list$mu__xREx__id_rho_z$parameters, list(mean = 0, sd = 0.5))
  expect_true(default_cs$formula_design$random_effects[[1]]$homogeneous_sd)
  expect_equal(colnames(default_cs$data$mu__xREx__id_xRE_DATAx), c("fa", "fb", "fc"))
  expect_equal(
    default_cs$formula_design$random_effects[[1]]$sd_parameter_names,
    rep("mu__xREx__id_sd", 3)
  )

  expect_error(
    JAGS_formula(
      formula = ~ 1 + cs(1 | id),
      parameter = "mu",
      data = continuous_df,
      prior_list = list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(
        id = random_block(sd = sd_prior, rho = prior("normal", list(0, 0.5)))
      )
    ),
    "does not support explicit '1', '0', or '-1' terms",
    fixed = TRUE
  )

  id_result <- JAGS_formula(
    formula = ~ 1 + x + id(1 + x | id),
    parameter = "mu",
    data = continuous_df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    ),
    prior_random = prior_random(id = random_block(sd = sd_prior))
  )
  expect_equal(names(id_result$prior_list), c("mu_intercept", "mu_x", "mu__xREx__id_sd"))
  expect_true(id_result$formula_design$random_effects[[1]]$homogeneous_sd)
  expect_equal(
    id_result$formula_design$random_effects[[1]]$sd_parameter_names,
    rep("mu__xREx__id_sd", 2)
  )

  diag_heterogeneous <- JAGS_formula(
    formula = ~ 1 + x + diag(1 + x | id, hom = FALSE),
    parameter = "mu",
    data = continuous_df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    ),
    prior_random = prior_random(id = random_block(sd = sd_prior))
  )
  expect_equal(names(diag_heterogeneous$prior_list), c("mu_intercept", "mu_x", "mu__xREx__id_intercept", "mu__xREx__id_x"))
  expect_false(diag_heterogeneous$formula_design$random_effects[[1]]$homogeneous_sd)
  expect_equal(
    diag_heterogeneous$formula_design$random_effects[[1]]$sd_parameter_names,
    c("mu__xREx__id_intercept", "mu__xREx__id_x")
  )

  hom_env <- new.env(parent = globalenv())
  hom_env$hom_flag <- FALSE
  hom_formula <- stats::as.formula("~ 1 + x + diag(1 + x | id, hom = hom_flag)", env = hom_env)
  hom_from_env <- JAGS_formula(
    formula = hom_formula,
    parameter = "mu",
    data = continuous_df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    ),
    prior_random = prior_random(id = random_block(sd = sd_prior))
  )
  expect_false(hom_from_env$formula_design$random_effects[[1]]$homogeneous_sd)

  hcs_index <- JAGS_formula(
    formula = ~ 1 + hcs(f | id),
    parameter = "mu",
    data = factor_df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(sd = sd_prior, rho = prior("normal", list(0, 0.5)))
    )
  )
  expect_equal(
    names(hcs_index$prior_list),
    c("mu_intercept", "mu__xREx__id_f", "mu__xREx__id_rho_z")
  )
  expect_false(hcs_index$formula_design$random_effects[[1]]$homogeneous_sd)
  expect_equal(
    hcs_index$formula_design$random_effects[[1]]$sd_parameter_names,
    paste0("mu__xREx__id_f[", 1:3, "]")
  )

  hcs_composite <- JAGS_formula(
    formula = ~ 1 + hcs(f + g | id),
    parameter = "mu",
    data = factor_df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(sd = sd_prior, rho = prior("normal", list(0, 0.5)))
    )
  )
  expect_equal(hcs_composite$formula_design$random_effects[[1]]$structured_index$variables, c("f", "g"))
  expect_equal(hcs_composite$formula_design$random_effects[[1]]$structured_index$label, "f:g")
  expect_equal(ncol(hcs_composite$data$mu__xREx__id_xRE_DATAx), 6L)
  expect_equal(
    hcs_composite$formula_design$random_effects[[1]]$sd_leaves$leaf_terms_by_column,
    paste0("f_g[", c("a.u", "a.v", "b.u", "b.v", "c.u", "c.v"), "]")
  )
  expect_equal(
    BayesTools:::.bt_random_effect_summary_sd_components(
      hcs_composite$formula_design$random_effects[[1]],
      paste0("mu__xREx__id_f_g[", 1:6, "]")
    ),
    paste0("f:g[", c("a.u", "a.v", "b.u", "b.v", "c.u", "c.v"), "]")
  )
  hcs_composite_raw <- BayesTools:::.bt_random_effect_summary_samples(
    model_samples = matrix(
      seq_len(7),
      nrow = 1,
      dimnames = list(
        NULL,
        c(paste0("mu__xREx__id_f_g[", 1:6, "]"), "mu__xREx__id_rho")
      )
    ),
    prior_list = hcs_composite$prior_list,
    formula_design = list(mu = hcs_composite$formula_design),
    mode = "raw"
  )
  hcs_composite_raw_samples <- BayesTools:::.rename_factor_levels(
    hcs_composite_raw$model_samples,
    hcs_composite_raw$prior_list
  )
  hcs_composite_raw_names <- colnames(hcs_composite_raw_samples)
  hcs_composite_raw_display <- BayesTools:::.bt_random_effect_summary_display_names(
    names = format_parameter_names(
      parameters = hcs_composite_raw_names,
      formula_parameters = unique(unlist(lapply(hcs_composite_raw$prior_list, attr, which = "parameter"))),
      formula_random = unique(unlist(lapply(hcs_composite_raw$prior_list, attr, which = "random_factor"))),
      formula_prefix = TRUE
    ),
    raw_names = hcs_composite_raw_names,
    prior_list = hcs_composite_raw$prior_list,
    formula_prefix = TRUE,
    formula_design = list(mu = hcs_composite$formula_design)
  )
  expect_true("(mu) sd(f:g[a.u] | id)" %in% hcs_composite_raw_display)
  expect_true("(mu) rho(id)" %in% hcs_composite_raw_display)
  expect_false(any(grepl("f_g", hcs_composite_raw_display, fixed = TRUE)))
  expect_false(any(grepl("sd((mu)", hcs_composite_raw_display, fixed = TRUE)))

  hcs_composite_allocation <- JAGS_formula(
    formula = ~ 1 + hcs(f + g | id),
    parameter = "mu",
    data = factor_df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      random_variance_allocation(
        name = "leaf_alloc",
        terms = "id",
        components = "sd",
        sd = sd_prior,
        allocation = prior("dirichlet", list(alpha = rep(1, 6)))
      ),
      id = random_block(rho = prior("normal", list(0, 0.5)))
    )
  )
  hcs_composite_allocation_summary <- BayesTools:::.bt_random_effect_summary_samples(
    model_samples = matrix(
      rep(1, 7),
      nrow = 1,
      dimnames = list(
        NULL,
        c(
          "mu__xRE_ALLOCx_leaf_alloc_total_sd",
          paste0("prior_par_eta_mu__xRE_ALLOCx_leaf_alloc_weight[", 1:6, "]")
        )
      )
    ),
    prior_list = hcs_composite_allocation$prior_list,
    formula_design = list(mu = hcs_composite_allocation$formula_design),
    mode = "full"
  )
  hcs_composite_allocation_display <- BayesTools:::.bt_random_effect_summary_display_names(
    names = colnames(hcs_composite_allocation_summary$model_samples),
    raw_names = colnames(hcs_composite_allocation_summary$model_samples),
    prior_list = hcs_composite_allocation_summary$prior_list
  )
  expect_true("(mu) var_frac(leaf_alloc: f:g[a.u])" %in% hcs_composite_allocation_display)
  expect_true("(mu) sd_mult(leaf_alloc: f:g[a.u])" %in% hcs_composite_allocation_display)
  expect_false(any(grepl("leaf_alloc: f_g", hcs_composite_allocation_display, fixed = TRUE)))

  malformed_allocation_design <- hcs_composite_allocation$formula_design
  malformed_allocation_design$random_effects[[1]]$allocation$target <- "sd"
  malformed_allocation_design$random_effects[[1]]$allocation$components <- NULL
  expect_error(
    BayesTools:::.bt_random_effect_summary_samples(
      model_samples = matrix(
        rep(1, 7),
        nrow = 1,
        dimnames = list(
          NULL,
          c(
            "mu__xRE_ALLOCx_leaf_alloc_total_sd",
            paste0("prior_par_eta_mu__xRE_ALLOCx_leaf_alloc_weight[", 1:6, "]")
          )
        )
      ),
      prior_list = hcs_composite_allocation$prior_list,
      formula_design = list(mu = malformed_allocation_design),
      mode = "full"
    ),
    "missing canonical 'allocation\\$components'"
  )

  malformed_allocation_design <- hcs_composite_allocation$formula_design
  malformed_allocation_design$random_effects[[1]]$allocation$components <- "variance"
  expect_error(
    BayesTools:::.bt_random_effect_summary_samples(
      model_samples = matrix(
        rep(1, 7),
        nrow = 1,
        dimnames = list(
          NULL,
          c(
            "mu__xRE_ALLOCx_leaf_alloc_total_sd",
            paste0("prior_par_eta_mu__xRE_ALLOCx_leaf_alloc_weight[", 1:6, "]")
          )
        )
      ),
      prior_list = hcs_composite_allocation$prior_list,
      formula_design = list(mu = malformed_allocation_design),
      mode = "full"
    ),
    "missing canonical 'allocation\\$components'"
  )

  malformed_allocation_design <- hcs_composite_allocation$formula_design
  malformed_allocation_design$random_effects[[1]]$allocation$scale <- NULL
  expect_error(
    BayesTools:::.bt_random_effect_summary_samples(
      model_samples = matrix(
        rep(1, 7),
        nrow = 1,
        dimnames = list(
          NULL,
          c(
            "mu__xRE_ALLOCx_leaf_alloc_total_sd",
            paste0("prior_par_eta_mu__xRE_ALLOCx_leaf_alloc_weight[", 1:6, "]")
          )
        )
      ),
      prior_list = hcs_composite_allocation$prior_list,
      formula_design = list(mu = malformed_allocation_design),
      mode = "full"
    ),
    "missing canonical 'allocation\\$scale'"
  )

  malformed_allocation_design <- hcs_composite_allocation$formula_design
  malformed_allocation_design$random_effects[[1]]$allocation$factors <- NULL
  expect_error(
    BayesTools:::.bt_random_effect_summary_samples(
      model_samples = matrix(
        rep(1, 7),
        nrow = 1,
        dimnames = list(
          NULL,
          c(
            "mu__xRE_ALLOCx_leaf_alloc_total_sd",
            paste0("prior_par_eta_mu__xRE_ALLOCx_leaf_alloc_weight[", 1:6, "]")
          )
        )
      ),
      prior_list = hcs_composite_allocation$prior_list,
      formula_design = list(mu = malformed_allocation_design),
      mode = "full"
    ),
    "missing canonical 'allocation\\$factors'"
  )

  ar1_result <- JAGS_formula(
    formula = ~ 1 + ar1(f | id),
    parameter = "mu",
    data = factor_df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(sd = sd_prior, rho = prior("normal", list(0, 0.5)))
    )
  )

  expect_equal(colnames(ar1_result$data$mu__xREx__id_xRE_DATAx), c("fa", "fb", "fc"))
  expect_equal(names(ar1_result$prior_list), c("mu_intercept", "mu__xREx__id_sd", "mu__xREx__id_rho_z"))
  expect_true(ar1_result$formula_design$random_effects[[1]]$homogeneous_sd)
  expect_equal(ar1_result$formula_design$random_effects[[1]]$structure, "ar1")
  expect_equal(
    ar1_result$formula_design$random_effects[[1]]$sd_parameter_names,
    rep("mu__xREx__id_sd", 3)
  )
  expect_match(ar1_result$formula_syntax, "pow(mu__xREx__id_rho, 2)", fixed = TRUE)

  ar1_implicit_levels <- JAGS_formula(
    formula = ~ 1 + ar1(f | id),
    parameter = "mu",
    data = factor_df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(sd = sd_prior, rho = prior("normal", list(0, 0.5)))
    )
  )
  expect_equal(
    ar1_implicit_levels$data$mu__xREx__id_xRE_DATAx,
    ar1_result$data$mu__xREx__id_xRE_DATAx
  )
  expect_equal(
    ar1_implicit_levels$formula_design$random_effects[[1]]$term_formula,
    ~ f - 1,
    ignore_formula_env = TRUE
  )

  hcs_result <- JAGS_formula(
    formula = ~ 1 + hcs(f | id),
    parameter = "mu",
    data = factor_df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(sd = sd_prior, rho = prior("normal", list(0, 0.5)))
    )
  )

  expect_equal(colnames(hcs_result$data$mu__xREx__id_xRE_DATAx), c("fa", "fb", "fc"))
  expect_equal(names(hcs_result$prior_list), c("mu_intercept", "mu__xREx__id_f", "mu__xREx__id_rho_z"))
  expect_false(hcs_result$formula_design$random_effects[[1]]$homogeneous_sd)
  expect_equal(
    hcs_result$formula_design$random_effects[[1]]$sd_parameter_names,
    paste0("mu__xREx__id_f[", 1:3, "]")
  )
  expect_equal(hcs_result$prior_list$mu__xREx__id_rho_z$truncation$lower, atanh(-1 / 2), tolerance = 1e-12)
  expect_s3_class(hcs_result$prior_list$mu__xREx__id_f, "prior.independent")
  expect_equal(attr(hcs_result$prior_list$mu__xREx__id_f, "levels"), 3L)
  expect_equal(attr(hcs_result$prior_list$mu__xREx__id_f, "level_names"), c("a", "b", "c"))
  expect_equal(
    hcs_result$formula_design$random_effects[[1]]$sd_leaves$leaf_terms_by_column,
    c("f[a]", "f[b]", "f[c]")
  )
  hcs_summary <- BayesTools:::.bt_random_effect_summary_samples(
    model_samples = matrix(
      c(1, 2, 3),
      nrow = 1,
      dimnames = list(NULL, paste0("mu__xREx__id_f[", 1:3, "]"))
    ),
    prior_list = hcs_result$prior_list,
    formula_design = list(mu = hcs_result$formula_design),
    mode = "standard"
  )
  hcs_display_names <- BayesTools:::.bt_random_effect_summary_display_names(
    names = colnames(hcs_summary$model_samples),
    raw_names = colnames(hcs_summary$model_samples),
    prior_list = hcs_summary$prior_list
  )
  expect_true("(mu) sd(f[a] | id)" %in% hcs_display_names)
  expect_false("(mu) sd(f[1] | id)" %in% hcs_display_names)
  hcs_implicit_levels <- JAGS_formula(
    formula = ~ 1 + hcs(f | id),
    parameter = "mu",
    data = factor_df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(sd = sd_prior, rho = prior("normal", list(0, 0.5)))
    )
  )
  expect_equal(
    hcs_implicit_levels$data$mu__xREx__id_xRE_DATAx,
    hcs_result$data$mu__xREx__id_xRE_DATAx
  )

  har_implicit_levels <- JAGS_formula(
    formula = ~ 1 + har(f | id),
    parameter = "mu",
    data = factor_df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(sd = sd_prior, rho = prior("normal", list(0, 0.5)))
    )
  )
  expect_equal(colnames(har_implicit_levels$data$mu__xREx__id_xRE_DATAx), c("fa", "fb", "fc"))
  expect_false(har_implicit_levels$formula_design$random_effects[[1]]$homogeneous_sd)

  expect_error(
    JAGS_formula(
      formula = ~ 1 + hcs(f | id),
      parameter = "mu",
      data = factor_df,
      prior_list = list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(
        id = random_block(
          sd = prior_factor("mnormal", list(0, 1), contrast = "orthonormal"),
          rho = prior("normal", list(0, 0.5))
        )
      )
    ),
    "ordinary scalar prior",
    fixed = TRUE
  )
  hcs_mixture_sd <- JAGS_formula(
    formula = ~ 1 + hcs(f | id),
    parameter = "mu",
    data = factor_df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(
        sd = prior_mixture(list(prior("gamma", list(2, 2)), prior("point", list(0)))),
        rho = prior("normal", list(0, 0.5))
      )
    )
  )
  expect_s3_class(hcs_mixture_sd$prior_list$mu__xREx__id_f, "prior.factor_mixture")

  raw_result <- NULL
  expect_warning(
    raw_result <- JAGS_formula(
      formula = ~ 1 + har(f | id),
      parameter = "mu",
      data = factor_df,
      prior_list = list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(
        id = random_block(
          sd = sd_prior,
          covariance = random_covariance(
            rho = prior("normal", list(0, 0.5)),
            rho_scale = "rho"
          )
        )
      )
    ),
    "valid correlation range",
    fixed = TRUE
  )

  expect_equal(names(raw_result$prior_list), c("mu_intercept", "mu__xREx__id_f", "mu__xREx__id_rho"))
  expect_equal(raw_result$prior_list$mu__xREx__id_rho$truncation, list(lower = -1, upper = 1))
  expect_equal(
    raw_result$formula_design$random_effects[[1]]$sd_parameter_names,
    paste0("mu__xREx__id_f[", 1:3, "]")
  )
  expect_false(grepl("rho_z", raw_result$formula_syntax, fixed = TRUE))
  expect_match(raw_result$formula_syntax, "pow(mu__xREx__id_rho, 2)", fixed = TRUE)

  raw_inherited <- NULL
  expect_warning(
    raw_inherited <- JAGS_formula(
      formula = ~ 1 + har(f | id),
      parameter = "mu",
      data = factor_df,
      prior_list = list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(
        covariance = random_covariance(
          rho = prior("normal", list(0, 0.5)),
          rho_scale = "rho"
        ),
        id = random_block(
          sd = sd_prior,
          covariance = random_covariance(structure = "har")
        )
      )
    ),
    "valid correlation range",
    fixed = TRUE
  )
  expect_equal(names(raw_inherited$prior_list), c("mu_intercept", "mu__xREx__id_f", "mu__xREx__id_rho"))
  expect_false(grepl("rho_z", raw_inherited$formula_syntax, fixed = TRUE))

  raw_block_override <- NULL
  expect_warning(
    raw_block_override <- JAGS_formula(
      formula = ~ 1 + har(f | id),
      parameter = "mu",
      data = factor_df,
      prior_list = list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(
        covariance = random_covariance(rho = prior("normal", list(0, 0.5))),
        id = random_block(
          sd = sd_prior,
          covariance = random_covariance(structure = "har", rho_scale = "rho")
        )
      )
    ),
    "valid correlation range",
    fixed = TRUE
  )
  expect_equal(names(raw_block_override$prior_list), c("mu_intercept", "mu__xREx__id_f", "mu__xREx__id_rho"))
  expect_false(grepl("rho_z", raw_block_override$formula_syntax, fixed = TRUE))

  raw_top_scale_block_rho <- NULL
  expect_warning(
    raw_top_scale_block_rho <- JAGS_formula(
      formula = ~ 1 + har(f | id),
      parameter = "mu",
      data = factor_df,
      prior_list = list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(
        covariance = random_covariance(rho_scale = "rho"),
        id = random_block(
          sd = sd_prior,
          rho = prior("normal", list(0, 0.5))
        )
      )
    ),
    "valid correlation range",
    fixed = TRUE
  )
  expect_equal(names(raw_top_scale_block_rho$prior_list), c("mu_intercept", "mu__xREx__id_f", "mu__xREx__id_rho"))
  expect_false(grepl("rho_z", raw_top_scale_block_rho$formula_syntax, fixed = TRUE))

  car_result <- JAGS_formula(
    formula = ~ 1 + car(0 + time | id),
    parameter = "mu",
    data = car_df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(sd = sd_prior, rho = prior("normal", list(0, 0.5)))
    )
  )
  car_term <- car_result$formula_design$random_effects[[1]]
  car_distance <- abs(outer(c(0, 0.5, 2), c(0, 0.5, 2), "-"))

  expect_equal(colnames(car_result$data$mu__xREx__id_xRE_DATAx), c("time_0", "time_0p5", "time_2"))
  expect_equal(names(car_result$prior_list), c("mu_intercept", "mu__xREx__id_sd", "mu__xREx__id_rho_z"))
  expect_true(car_term$homogeneous_sd)
  expect_equal(car_term$structure, "car")
  expect_equal(car_term$car$time_variable, "time")
  expect_equal(car_term$car$time_values, c(0, 0.5, 2))
  expect_equal(car_term$car$distance_matrix, car_distance)
  expect_equal(car_term$correlation$bounds, c(lower = 0, upper = 1))
  expect_equal(car_term$correlation$distance_matrix, car_distance)
  expect_equal(car_result$prior_list$mu__xREx__id_rho_z$truncation$lower, 0)
  expect_equal(car_result$prior_list$mu__xREx__id_rho_z$truncation$upper, Inf)
  expect_match(car_result$formula_syntax, "mu__xREx__id_rho <- 2 * ilogit(2 * mu__xREx__id_rho_z) - 1", fixed = TRUE)
  expect_match(car_result$formula_syntax, "pow(mu__xREx__id_rho, 1.5)", fixed = TRUE)
  expect_match(car_result$formula_syntax, "pow(mu__xREx__id_rho, 2)", fixed = TRUE)

  car_explicit_no_intercept <- JAGS_formula(
    formula = ~ 1 + car(0 + time | id),
    parameter = "mu",
    data = car_df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(sd = sd_prior, rho = prior("normal", list(0, 0.5)))
    )
  )
  expect_equal(
    car_explicit_no_intercept$data$mu__xREx__id_xRE_DATAx,
    car_result$data$mu__xREx__id_xRE_DATAx
  )
  expect_error(
    BayesTools:::.bt_random_effect_prediction_data(
      car_term,
      data.frame(time = 1, id = factor("a", levels = c("a", "b")))
    ),
    "New CAR time coordinate",
    fixed = TRUE
  )

  car_scale_request <- JAGS_formula(
    formula = ~ 1 + car(0 + time | id),
    parameter = "mu",
    data = car_df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    formula_scale = TRUE,
    prior_random = prior_random(
      id = random_block(sd = sd_prior, rho = prior("normal", list(0, 0.5)))
    )
  )
  expect_false("formula_scale" %in% names(car_scale_request))
  expect_equal(car_scale_request$formula_design$random_effects[[1]]$car$time_values, c(0, 0.5, 2))

  fixed_car_scale <- JAGS_formula(
    formula = ~ 1 + time + car(0 + time | id),
    parameter = "mu",
    data = car_df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      time = prior("normal", list(0, 1))
    ),
    formula_scale = TRUE,
    prior_random = prior_random(
      id = random_block(sd = sd_prior, rho = prior("normal", list(0, 0.5)))
    )
  )
  expect_equal(fixed_car_scale$formula_design$random_effects[[1]]$car$time_values, c(0, 0.5, 2))
  expect_equal(fixed_car_scale$formula_design$random_effects[[1]]$car$distance_matrix, car_distance)
  expect_equal(
    fixed_car_scale$data$mu__xREx__id_xRE_DATAx,
    car_result$data$mu__xREx__id_xRE_DATAx
  )
  expect_equal(
    unname(fixed_car_scale$data$mu_data_time),
    unname(as.numeric(scale(car_df$time)))
  )
  expect_equal(names(fixed_car_scale$formula_scale), "mu_time")

  car_prediction_term <- fixed_car_scale$formula_design$random_effects[[1]]
  car_coefficient_names <- as.vector(BayesTools:::.bt_random_effect_coefficient_names(
    random_term = car_prediction_term,
    n_groups = length(car_prediction_term$group_levels),
    n_columns = car_prediction_term$n_columns
  ))
  car_prediction_posterior <- matrix(
    0,
    nrow = 1,
    ncol = 2L + length(car_coefficient_names),
    dimnames = list(NULL, c("mu_intercept", "mu_time", car_coefficient_names))
  )
  car_prediction_posterior[, "mu_intercept"] <- 1
  car_prediction_posterior[, "mu_time"] <- 0.25
  car_prediction_fit <- coda::mcmc(car_prediction_posterior)
  attr(car_prediction_fit, "formula_design") <- list(mu = fixed_car_scale$formula_design)
  attr(car_prediction_fit, "formula_scale") <- list(mu = fixed_car_scale$formula_scale)
  expect_equal(
    unname(drop(JAGS_evaluate_formula(
      fit = car_prediction_fit,
      formula = ~ 1 + time + car(0 + time | id),
      parameter = "mu",
      data = car_df,
      prior_list = fixed_car_scale$prior_list
    ))),
    1 + 0.25 * unname(as.numeric(scale(car_df$time))),
    tolerance = 1e-12
  )

  expect_error(
    JAGS_formula(
      formula = ~ 1 + car(0 + time | id),
      parameter = "mu",
      data = transform(car_df, time = c(0, NA, 2, 0, 0.5, 2)),
      prior_list = list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(
        id = random_block(sd = sd_prior, rho = prior("normal", list(0, 0.5)))
      )
    ),
    "must contain only finite values",
    fixed = TRUE
  )
  bad_ordered_car_df <- car_df
  bad_ordered_car_df$time <- ordered(
    c("early", "mid", "late", "early", "mid", "late"),
    levels = c("early", "mid", "late")
  )
  expect_error(
    JAGS_formula(
      formula = ~ 1 + car(0 + time | id),
      parameter = "mu",
      data = bad_ordered_car_df,
      prior_list = list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(
        id = random_block(sd = sd_prior, rho = prior("normal", list(0, 0.5)))
      )
    ),
    "numeric level labels",
    fixed = TRUE
  )

  ordered_car_df <- car_df
  ordered_car_df$time <- ordered(
    as.character(ordered_car_df$time),
    levels = c("0", "0.5", "2")
  )
  ordered_car_result <- JAGS_formula(
    formula = ~ 1 + car(0 + time | id),
    parameter = "mu",
    data = ordered_car_df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(sd = sd_prior, rho = prior("normal", list(0, 0.5)))
    )
  )
  ordered_car_term <- ordered_car_result$formula_design$random_effects[[1]]
  expect_equal(ordered_car_term$car$time_values, c(0, 0.5, 2))
  expect_equal(
    ordered_car_result$data$mu__xREx__id_xRE_DATAx,
    car_result$data$mu__xREx__id_xRE_DATAx
  )

  car_logit <- JAGS_formula(
    formula = ~ 1 + car(0 + time | id),
    parameter = "mu",
    data = car_df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(
        sd = sd_prior,
        covariance = random_covariance(
          rho = prior("normal", list(0, 0.5)),
          rho_scale = "logit"
        )
      )
    )
  )
  expect_equal(names(car_logit$prior_list), c("mu_intercept", "mu__xREx__id_sd", "mu__xREx__id_rho_logit"))
  expect_match(car_logit$formula_syntax, "mu__xREx__id_rho <- 0 + 1 * ilogit(mu__xREx__id_rho_logit)", fixed = TRUE)

  ar_logit <- JAGS_formula(
    formula = ~ 1 + ar1(f | id),
    parameter = "mu",
    data = factor_df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(
        sd = sd_prior,
        covariance = random_covariance(
          rho = prior("normal", list(0, 0.5)),
          rho_scale = "logit"
        )
      )
    )
  )
  expect_equal(names(ar_logit$prior_list), c("mu_intercept", "mu__xREx__id_sd", "mu__xREx__id_rho_logit"))
  expect_match(ar_logit$formula_syntax, "mu__xREx__id_rho <- -1 + 2 * ilogit(mu__xREx__id_rho_logit)", fixed = TRUE)

  expect_error(
    JAGS_formula(
      formula = ~ 1 + car(0 + x + time | id),
      parameter = "mu",
      data = transform(car_df, x = seq_len(nrow(car_df))),
      prior_list = list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(
        id = random_block(sd = sd_prior, rho = prior("normal", list(0, 0.5)))
      )
    ),
    "exactly one untransformed time variable",
    fixed = TRUE
  )
  expect_error(
    JAGS_formula(
      formula = ~ 1 + car(0 + time | id),
      parameter = "mu",
      data = transform(car_df, time = factor(time)),
      prior_list = list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(
        id = random_block(sd = sd_prior, rho = prior("normal", list(0, 0.5)))
      )
    ),
    "must be numeric or an ordered factor",
    fixed = TRUE
  )
})

test_that("structured correlation Cholesky syntax exposes intended covariance patterns", {

  sd_prior <- prior("normal", list(0, 1), truncation = list(lower = 0, upper = Inf))
  rho_prior <- prior("normal", list(0, 0.5))
  block_prior <- random_block(sd = sd_prior, rho = rho_prior)

  cs_module <- BayesTools:::.bt_JAGS_structured_corr_cholesky(
    node_prefix = "mu__xREx__id",
    prior_prefix = "_xREx__id",
    K = 3,
    structure = "cs",
    block_prior = block_prior,
    include_correlation = TRUE
  )
  expect_match(cs_module$syntax, "mu__xREx__id_xRE_CORx_R[1,2] <- mu__xREx__id_rho", fixed = TRUE)
  expect_match(cs_module$syntax, "mu__xREx__id_xRE_CORx_L[1,2] <- 0", fixed = TRUE)
  expect_match(cs_module$syntax, "mu__xREx__id_xRE_CORx_L[2,2] <- sqrt", fixed = TRUE)

  ar_module <- BayesTools:::.bt_JAGS_structured_corr_cholesky(
    node_prefix = "mu__xREx__id",
    prior_prefix = "_xREx__id",
    K = 4,
    structure = "ar1",
    block_prior = block_prior,
    include_correlation = TRUE
  )
  expect_match(ar_module$syntax, "mu__xREx__id_xRE_CORx_R[1,3] <- pow(mu__xREx__id_rho, 2)", fixed = TRUE)
  expect_match(ar_module$syntax, "mu__xREx__id_xRE_CORx_R[1,4] <- pow(mu__xREx__id_rho, 3)", fixed = TRUE)
  expect_match(ar_module$syntax, "mu__xREx__id_xRE_CORx_L[1,4] <- 0", fixed = TRUE)
  expect_match(ar_module$syntax, "mu__xREx__id_xRE_CORx_L[4,4] <- sqrt", fixed = TRUE)

  car_distance <- abs(outer(c(0, 1.5, 3), c(0, 1.5, 3), "-"))
  car_module <- BayesTools:::.bt_JAGS_structured_corr_cholesky(
    node_prefix = "mu__xREx__id",
    prior_prefix = "_xREx__id",
    K = 3,
    structure = "car",
    block_prior = block_prior,
    include_correlation = TRUE,
    distance_matrix = car_distance
  )
  expect_match(car_module$syntax, "mu__xREx__id_xRE_CORx_R[1,2] <- pow(mu__xREx__id_rho, 1.5)", fixed = TRUE)
  expect_match(car_module$syntax, "mu__xREx__id_xRE_CORx_R[1,3] <- pow(mu__xREx__id_rho, 3)", fixed = TRUE)
  expect_equal(car_module$prior_list$`_xREx__id_rho_z`$truncation$lower, 0)
  expect_equal(car_module$bridge$bounds, c(lower = 0, upper = 1))
  expect_equal(car_module$bridge$distance_matrix, car_distance)

  logit_module <- BayesTools:::.bt_JAGS_structured_corr_cholesky(
    node_prefix = "mu__xREx__id",
    prior_prefix = "_xREx__id",
    K = 3,
    structure = "cs",
    block_prior = random_block(
      sd = sd_prior,
      covariance = random_covariance(
        rho = rho_prior,
        rho_scale = "logit"
      )
    ),
    include_correlation = TRUE
  )
  expect_match(logit_module$syntax, "mu__xREx__id_rho <- -0.5 + 1.5 * ilogit(mu__xREx__id_rho_logit)", fixed = TRUE)
  expect_true("_xREx__id_rho_logit" %in% names(logit_module$prior_list))
})

.jags_lme4_random_oracle_data <- function() {
  data.frame(
    y = c(1.5, 0.2, 2.1, 3.4, -0.7, 4.2),
    x = c(-1, 0, 1, 2, -2, 3),
    id = factor(c("b", "a", "b", "c", "a", "c"), levels = c("c", "a", "b"))
  )
}

.lme4_lFormula_or_fallback <- function(formula, fallback_formula, data) {
  result <- try(lme4::lFormula(formula, data = data), silent = TRUE)
  if (inherits(result, "try-error")) {
    result <- lme4::lFormula(fallback_formula, data = data)
  }
  result
}

.expect_jags_random_design_matches_lme4 <- function(jags_result, lme4_result,
                                                    expected_re_columns) {
  bayes_random_data <- jags_result$data$mu__xREx__id_xRE_DATAx
  bayes_random_map <- jags_result$data$mu__xREx__id_xRE_MAPx
  lme4_re_terms <- lme4_result$reTrms
  lme4_group <- lme4_re_terms$flist$id
  lme4_Z <- t(as.matrix(lme4_re_terms$Zt))

  expect_equal(dim(jags_result$formula_design$model_matrix), dim(lme4_result$X))
  expect_equal(
    as.vector(jags_result$formula_design$model_matrix),
    as.vector(lme4_result$X)
  )
  expect_equal(colnames(jags_result$formula_design$model_matrix), colnames(lme4_result$X))

  expect_equal(bayes_random_map, unname(as.integer(lme4_group)))
  expect_equal(levels(lme4_group), levels(.jags_lme4_random_oracle_data()$id))
  expect_equal(attr(jags_result$formula_design$random_effects[[1]], "grouping_factor"), "id")
  expect_true(isTRUE(attr(jags_result$formula_design$random_effects[[1]], "independent")))
  expect_s3_class(jags_result$formula_design$random_effects[[1]], "BayesTools_random_effect_term")
  expect_equal(jags_result$formula_design$random_effects[[1]]$structure, "diag")
  expect_equal(
    jags_result$formula_design$jags_data_names[["__xREx__id"]],
    c("mu__xREx__id_xRE_DATAx", "mu__xREx__id_xRE_MAPx")
  )

  expect_equal(dim(bayes_random_data), c(nrow(lme4_result$fr), length(expected_re_columns)))
  expect_equal(colnames(bayes_random_data), expected_re_columns)
  expect_equal(lme4_re_terms$cnms$id, expected_re_columns)
  expect_equal(dim(lme4_Z), c(nrow(lme4_result$fr), length(expected_re_columns) * nlevels(lme4_group)))

  for (column_i in seq_along(expected_re_columns)) {
    selected <- bayes_random_map + (column_i - 1L) * nlevels(lme4_group)
    expect_equal(unname(bayes_random_data[, column_i]), unname(lme4_Z[cbind(seq_len(nrow(lme4_Z)), selected)]))
  }
}

test_that("JAGS_formula independent random-effect design matches lme4 lFormula oracles", {

  skip_if_not_installed("lme4")

  df <- .jags_lme4_random_oracle_data()

  intercept_result <- JAGS_formula(
    formula = ~ x + (1 || id),
    parameter = "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    ),
    prior_random = prior_random(
      id = random_block(sd = prior("gamma", list(2, 2)))
    )
  )
  intercept_lme4 <- .lme4_lFormula_or_fallback(
    y ~ x + (1 || id),
    fallback_formula = y ~ x + (1 | id),
    data = df
  )

  .expect_jags_random_design_matches_lme4(
    jags_result = intercept_result,
    lme4_result = intercept_lme4,
    expected_re_columns = "(Intercept)"
  )

  slope_result <- JAGS_formula(
    formula = ~ x + (0 + x || id),
    parameter = "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    ),
    prior_random = prior_random(
      id = random_block(sd = prior("gamma", list(2, 2)))
    )
  )
  slope_lme4 <- lme4::lFormula(y ~ x + (0 + x || id), data = df)

  .expect_jags_random_design_matches_lme4(
    jags_result = slope_result,
    lme4_result = slope_lme4,
    expected_re_columns = "x"
  )

  df_us <- data.frame(
    y = seq_len(12) / 10,
    x = c(-1, 0, 1, 2, -2, 3, -0.5, 0.5, 1.5, -1.5, 2.5, -2.5),
    id = factor(rep(c("c", "a", "b"), 4), levels = c("c", "a", "b"))
  )
  us_result <- JAGS_formula(
    formula = ~ x + (1 + x | id),
    parameter = "mu",
    data = df_us,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    ),
    prior_random = prior_random(
      id = random_block(
        sd = prior("gamma", list(2, 2)),
        cor = prior_lkj(eta = 1)
      )
    )
  )
  us_lme4 <- lme4::lFormula(y ~ x + (1 + x | id), data = df_us)
  us_map <- unname(as.integer(us_lme4$reTrms$flist$id))

  expect_equal(us_result$formula_design$random_effects[[1]]$structure, "us")
  expect_equal(
    as.vector(us_result$formula_design$model_matrix),
    as.vector(us_lme4$X)
  )
  expect_equal(us_result$formula_design$random_effects[[1]]$column_names, us_lme4$reTrms$cnms$id)
  expect_equal(us_result$data$mu__xREx__id_xRE_MAPx, us_map)
  expect_equal(colnames(us_result$data$mu__xREx__id_xRE_DATAx), c("(Intercept)", "x"))
  expect_equal(unname(us_result$data$mu__xREx__id_xRE_DATAx[, "(Intercept)"]), rep(1, nrow(df_us)))
  expect_equal(unname(us_result$data$mu__xREx__id_xRE_DATAx[, "x"]), df_us$x)
  expect_match(us_result$formula_syntax, "dbt_lkj_cpc", fixed = TRUE)

  df_ar1 <- data.frame(
    y = seq_len(12) / 10,
    f = factor(rep(c("a", "b", "c"), 4), levels = c("a", "b", "c")),
    id = factor(rep(c("g1", "g2", "g3", "g4"), each = 3))
  )
  ar1_result <- JAGS_formula(
    formula = ~ 1 + ar1(f | id),
    parameter = "mu",
    data = df_ar1,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(
        sd = prior("normal", list(0, 1), truncation = list(lower = 0, upper = Inf)),
        rho = prior("normal", list(0, 0.5))
      )
    )
  )
  ar1_lme4 <- lme4::lFormula(y ~ 1 + ar1(0 + f | id), data = df_ar1)

  expect_equal(ar1_result$formula_design$random_effects[[1]]$column_names, ar1_lme4$reTrms$cnms$id)
  expect_equal(ar1_result$data$mu__xREx__id_xRE_MAPx, unname(as.integer(ar1_lme4$reTrms$flist$id)))
  expect_equal(colnames(ar1_result$data$mu__xREx__id_xRE_DATAx), c("fa", "fb", "fc"))
  expect_equal(ar1_result$formula_design$random_effects[[1]]$structure, "ar1")
  expect_true(ar1_result$formula_design$random_effects[[1]]$homogeneous_sd)
})
