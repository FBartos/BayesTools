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

test_that("independent random-effect design exposes grouping maps and rejects correlated syntax", {

  df <- data.frame(
    x = c(-1, 0, 1, 2, -2, 3),
    id = factor(c("b", "a", "b", "c", "a", "c"), levels = c("a", "b", "c"))
  )
  prior_list <- list(
    intercept = prior("gamma", list(2, 2)),
    x = prior("gamma", list(2, 2))
  )

  random <- BayesTools:::.JAGS_random_effect_formula(
    formula = .extract_random_effects(~ 1 + x + (1 + x || id))[[1]],
    parameter = "mu",
    data = df,
    prior_list = prior_list
  )

  expect_equal(random$data$mu__xREx__id_xRE_MAPx, c(2L, 1L, 2L, 3L, 1L, 3L))
  expect_equal(unname(random$data$mu__xREx__id_xRE_DATAx[, "x"]), df$x)
  expect_equal(colnames(random$data$mu__xREx__id_xRE_DATAx), c("(Intercept)", "x"))
  expect_equal(
    names(random$prior_list),
    c("_xREx__id_intercept", "_xREx__id_x")
  )
  expect_equal(
    unname(vapply(random$prior_list, attr, character(1), which = "random_factor")),
    c("id", "id")
  )
  expect_true(isTRUE(attr(random$prior_list[["_xREx__id_x"]], "random_sd")))
  expect_equal(random$random_scale_terms, c("_xREx__id_intercept" = "intercept", "_xREx__id_x" = "x"))
  expect_identical(random$formula_term, "mu__xREx__id[i]")

  expect_error(
    BayesTools:::.JAGS_random_effect_formula(
      formula = .extract_random_effects(~ 1 + x + (1 + x | id))[[1]],
      parameter = "mu",
      data = df,
      prior_list = prior_list
    ),
    "Only independent random effects are supported yet.",
    fixed = TRUE
  )
})

test_that("JAGS_formula random-effect design exposes grouping maps and public variable names", {

  df <- data.frame(
    x = c(-1, 0, 1, 2, -2, 3),
    id = factor(c("b", "a", "b", "c", "a", "c"), levels = c("c", "a", "b"))
  )
  prior_list <- list(
    intercept = prior("normal", list(0, 1)),
    x = prior("normal", list(0, 1)),
    "intercept|id" = prior("gamma", list(2, 2)),
    "x|id" = prior("gamma", list(2, 2))
  )

  result <- JAGS_formula(
    formula = ~ 1 + x + (1 + x || id),
    parameter = "mu",
    data = df,
    prior_list = prior_list
  )

  expected_map <- match(as.character(df$id), levels(df$id))

  expect_equal(result$data$mu__xREx__id_xRE_MAPx, expected_map)
  expect_equal(result$data$mu__xREx__id_xRE_MAPx, c(3L, 2L, 3L, 1L, 2L, 1L))
  expect_equal(dim(result$data$mu__xREx__id_xRE_DATAx), c(6L, 2L))
  expect_equal(colnames(result$data$mu__xREx__id_xRE_DATAx), c("(Intercept)", "x"))
  expect_equal(unname(result$data$mu__xREx__id_xRE_DATAx[, "x"]), df$x)

  expect_equal(
    result$formula_design$jags_data_names[["__xREx__id"]],
    c("mu__xREx__id_xRE_DATAx", "mu__xREx__id_xRE_MAPx")
  )
  expect_true(isTRUE(attr(result$formula_design$random_effects[[1]], "independent")))
  expect_equal(attr(result$formula_design$random_effects[[1]], "grouping_factor"), "id")
  expect_equal(
    names(result$prior_list),
    c("mu_intercept", "mu_x", "mu__xREx__id_intercept", "mu__xREx__id_x")
  )
  expect_true(isTRUE(attr(result$prior_list$mu__xREx__id_x, "random_sd")))
  expect_equal(attr(result$prior_list$mu__xREx__id_x, "random_factor"), "id")
  expect_match(result$formula_syntax, "for\\(i in 1:3\\)")
  expect_match(result$formula_syntax, "dmnorm\\(rep\\(0, 2\\)")

  expect_error(
    JAGS_formula(
      formula = ~ 1 + x + (1 + x | id),
      parameter = "mu",
      data = df,
      prior_list = prior_list
    ),
    "Only independent random effects are supported yet.",
    fixed = TRUE
  )
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
      x = prior("normal", list(0, 1)),
      "intercept|id" = prior("gamma", list(2, 2))
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
      x = prior("normal", list(0, 1)),
      "x|id" = prior("gamma", list(2, 2))
    )
  )
  slope_lme4 <- lme4::lFormula(y ~ x + (0 + x || id), data = df)

  .expect_jags_random_design_matches_lme4(
    jags_result = slope_result,
    lme4_result = slope_lme4,
    expected_re_columns = "x"
  )
})
