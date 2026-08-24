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

test_that("BayesTools contrasts resolve without search-path lookup", {

  contrast_names <- c(
    "contr.treatment",
    "contr.independent",
    "contr.orthonormal",
    "contr.meandif",
    "contr.ordered_cumulative",
    "contr.ordered_cumulative_levels"
  )

  for(contrast_name in contrast_names){
    factor_data <- data.frame(
      group = factor(c("a", "b", "c", "a"), levels = c("a", "b", "c"))
    )
    if(grepl("ordered", contrast_name, fixed = TRUE)){
      factor_data$group <- ordered(factor_data$group, levels = c("a", "b", "c"))
    }
    attr(factor_data$group, "contrasts") <- contrast_name

    model_frame <- stats::model.frame(~ group, data = factor_data)
    observed <- .bt_model_matrix(
      model_frame,
      formula = ~ group,
      data    = factor_data
    )
    expected <- stats::model.matrix(
      ~ group,
      data = factor_data,
      contrasts.arg = list(
        group = .factor_contrast_matrix(levels(factor_data$group), contrast_name)
      )
    )
    attr(expected, "contrasts")[["group"]] <- contrast_name

    expect_equal(observed, expected, info = contrast_name)
    expect_identical(
      attr(observed, "contrasts", exact = TRUE)[["group"]],
      contrast_name,
      info = contrast_name
    )

    random_design <- .bt_random_effect_design_matrix(
      ~ 0 + group,
      data = factor_data
    )
    expect_identical(
      attr(random_design$model_matrix, "contrasts", exact = TRUE)[["group"]],
      contrast_name,
      info = paste(contrast_name, "random-effect metadata")
    )
  }

  interaction_data <- expand.grid(
    group_a = factor(c("a", "b", "c"), levels = c("a", "b", "c")),
    group_b = factor(c("u", "v"), levels = c("u", "v"))
  )
  attr(interaction_data$group_a, "contrasts") <- "contr.orthonormal"
  attr(interaction_data$group_b, "contrasts") <- "contr.meandif"
  interaction_design <- .bt_random_effect_design_matrix(
    ~ 0 + group_a:group_b,
    data = interaction_data
  )

  expect_identical(
    attr(interaction_design$model_matrix, "contrasts", exact = TRUE),
    list(
      group_a = "contr.orthonormal",
      group_b = "contr.meandif"
    )
  )

  sum_data <- interaction_data
  attr(sum_data$group_a, "contrasts") <- "contr.sum"
  attr(sum_data$group_b, "contrasts") <- "contr.orthonormal"
  sum_frame <- stats::model.frame(~ group_a * group_b, data = sum_data)
  sum_observed <- .bt_model_matrix(
    sum_frame,
    formula = ~ group_a * group_b,
    data    = sum_data
  )
  sum_expected <- stats::model.matrix(
    ~ group_a * group_b,
    data = sum_data,
    contrasts.arg = list(
      group_a = stats::contr.sum(levels(sum_data$group_a)),
      group_b = contr.orthonormal(levels(sum_data$group_b))
    )
  )
  attr(sum_expected, "contrasts") <- list(
    group_a = "contr.sum",
    group_b = "contr.orthonormal"
  )

  expect_equal(sum_observed, sum_expected)

  prediction_data <- data.frame(
    group = factor(rep(c("a", "b", "c"), 2), levels = c("a", "b", "c")),
    id    = factor(rep(c("g1", "g2", "g3"), each = 2))
  )
  random_result <- suppressWarnings(JAGS_formula(
    formula   = ~ 1 + diag(0 + group | id),
    parameter = "mu",
    data      = prediction_data,
    prior_list = list(
      intercept = prior("normal", list(0, 1))
    ),
    prior_random = prior_random(
      id = random_block(
        sd = prior_factor("mnormal", list(0, 1), contrast = "orthonormal")
      )
    )
  ))
  random_term <- random_result$formula_design$random_effects[[1L]]
  prediction <- .bt_random_effect_prediction_data(
    random_term,
    data = prediction_data
  )

  expect_equal(prediction$model_matrix, random_term$model_matrix)
  expect_equal(prediction$group_map, random_term$group_map)
  expect_identical(prediction$group_levels, random_term$group_levels)
})

test_that("BayesTools contrasts accept scalar level counts and validate inputs", {
  contrast_functions <- list(
    orthonormal = contr.orthonormal,
    meandif = contr.meandif,
    independent = contr.independent,
    ordered_cumulative = contr.ordered_cumulative,
    ordered_cumulative_levels = contr.ordered_cumulative_levels
  )
  expected_columns <- c(
    orthonormal = 2L,
    meandif = 2L,
    independent = 3L,
    ordered_cumulative = 2L,
    ordered_cumulative_levels = 3L
  )

  for(name in names(contrast_functions)){
    contrast_function <- contrast_functions[[name]]
    scalar_result <- contrast_function(3)
    expect_true(is.matrix(scalar_result), info = name)
    expect_equal(dim(scalar_result), c(3L, expected_columns[[name]]), info = name)
    expect_equal(contrast_function(3, contrasts = FALSE), diag(3), info = name)
    expect_error(contrast_function(3, contrasts = NA), "cannot contain NA", info = name)
    expect_error(contrast_function(3, contrasts = c(TRUE, FALSE)), "length '1'", info = name)
  }

  expect_equal(contr.independent(1), matrix(1, 1, 1))
  expect_equal(contr.ordered_cumulative_levels(1), matrix(1, 1, 1))
  expect_equal(contr.independent("level"), matrix(1, 1, 1))
  expect_equal(contr.ordered_cumulative_levels("level"), matrix(1, 1, 1))
  expect_error(contr.orthonormal(1), "Not enough degrees of freedom")
  expect_error(contr.meandif(1), "Not enough degrees of freedom")
  expect_error(contr.ordered_cumulative(1), "Not enough degrees of freedom")
  expect_error(contr.orthonormal("level"), "Not enough degrees of freedom")
  expect_error(contr.meandif("level"), "Not enough degrees of freedom")
  expect_error(contr.ordered_cumulative("level"), "Not enough degrees of freedom")

  expect_error(contr.orthonormal(NA_real_), "cannot contain NA")
  expect_error(contr.orthonormal(Inf), "finite values")
  expect_error(contr.orthonormal(2.5), "integer vector")
})

test_that("factor contrast priors agree across main effects and interactions", {
  data <- data.frame(
    x = rep(c(-1, 1), 3),
    group = factor(
      rep(c("a", "b", "c"), each = 2),
      levels = c("a", "b", "c")
    )
  )

  expect_error(
    JAGS_formula(
      ~ x * group,
      parameter = "mu",
      data = data,
      prior_list = list(
        intercept = prior("normal", list(0, 1)),
        x = prior("normal", list(0, 1)),
        group = prior_factor("normal", list(0, 1), contrast = "treatment"),
        "x:group" = prior_factor("normal", list(0, 1), contrast = "independent")
      )
    ),
    "conflicting contrast priors across formula terms",
    fixed = TRUE
  )

  expect_error(
    JAGS_formula(
      ~ x * group,
      parameter = "mu",
      data = data,
      prior_list = list(
        intercept = prior("normal", list(0, 1)),
        x = prior("normal", list(0, 1)),
        group = prior_factor("normal", list(0, 1), contrast = "independent"),
        "x:group" = prior_factor("normal", list(0, 1), contrast = "treatment")
      )
    ),
    "conflicting contrast priors across formula terms",
    fixed = TRUE
  )
})

test_that("formula_add_intercept handles grouped and unary no-intercept terms", {
  cases <- list(
    list(~ (x - 1), ~ x),
    list(~ ((0 + x)), ~ x),
    list(~ (x + z - 1), ~ x + z),
    list(~ x + (-1), ~ x),
    list(~ x + ((-1)), ~ x),
    list(~ x + (+0), ~ x),
    list(~ x - (+1), ~ x)
  )

  for(case in cases){
    expect_equal(
      formula_add_intercept(case[[1L]]),
      case[[2L]],
      ignore_formula_env = TRUE
    )
  }

  expect_equal(
    formula_add_intercept(~ I(x - 1) - 1),
    ~ I(x - 1),
    ignore_formula_env = TRUE
  )
  expect_equal(
    formula_add_intercept(~ offset(x - 1) - 1),
    ~ offset(x - 1),
    ignore_formula_env = TRUE
  )
})

test_that("JAGS_fit rejects prior-backed add_parameters before fitting", {
  prior_list <- list(mu = prior("normal", list(0, 1)))

  expect_error(
    JAGS_fit(
      model_syntax = "model{}",
      prior_list = prior_list,
      add_parameters = "mu"
    ),
    "already monitored through 'prior_list'",
    fixed = TRUE
  )
  expect_error(
    JAGS_fit(
      model_syntax = "model{}",
      prior_list = prior_list,
      add_parameters = "mu[1]"
    ),
    "already monitored through 'prior_list'",
    fixed = TRUE
  )
  expect_error(
    JAGS_fit(
      model_syntax = "model{}",
      prior_list = prior_list,
      add_parameters = ""
    ),
    "cannot contain empty parameter names",
    fixed = TRUE
  )
  expect_error(
    JAGS_fit(
      model_syntax = "model{}",
      prior_list = prior_list,
      add_parameters = NA_character_
    ),
    "cannot contain NA",
    fixed = TRUE
  )
})

test_that("JAGS_fit requires uniquely named formula-indexed lists", {
  formula <- ~ 1
  formula_data <- data.frame(row = 1)
  formula_priors <- list(intercept = prior("normal", list(0, 1)))

  expect_error(
    JAGS_fit(
      model_syntax = "model{}",
      formula_list = list(formula),
      formula_data_list = list(mu = formula_data),
      formula_prior_list = list(mu = formula_priors)
    ),
    "'formula_list' argument must be a fully named list",
    fixed = TRUE
  )
  expect_error(
    JAGS_fit(
      model_syntax = "model{}",
      formula_list = structure(
        list(formula, formula),
        names = c("mu", "mu")
      ),
      formula_data_list = list(mu = formula_data),
      formula_prior_list = list(mu = formula_priors)
    ),
    "'formula_list' argument must not contain duplicate names",
    fixed = TRUE
  )
  expect_error(
    JAGS_fit(
      model_syntax = "model{}",
      formula_list = list(mu = formula),
      formula_data_list = list(formula_data),
      formula_prior_list = list(mu = formula_priors)
    ),
    "'formula_data_list' argument must be a fully named list",
    fixed = TRUE
  )
  expect_error(
    JAGS_fit(
      model_syntax = "model{}",
      formula_list = list(mu = formula),
      formula_data_list = list(mu = formula_data),
      formula_prior_list = list(mu = formula_priors),
      formula_scale_list = structure(
        list(NULL, NULL),
        names = c("mu", "mu")
      )
    ),
    "'formula_scale_list' argument must not contain duplicate names",
    fixed = TRUE
  )
  expect_error(
    JAGS_fit(
      model_syntax = "model{}",
      formula_data_list = list(mu = formula_data)
    ),
    "'formula_data_list' argument cannot be supplied without 'formula_list'",
    fixed = TRUE
  )
  expect_error(
    JAGS_fit(
      model_syntax = "model{}",
      formula_random_effects_compile_list = list(
        mu = random_effects_compile(marginalized = "block")
      )
    ),
    "'formula_random_effects_compile_list' argument cannot be supplied without 'formula_list'",
    fixed = TRUE
  )
})

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

test_that("JAGS_formula rejects matrix-valued continuous predictors", {
  data <- data.frame(row = seq_len(3))
  data$m <- I(matrix(seq_len(6), nrow = 3))

  expect_error(
    JAGS_formula(
      ~ m,
      "mu",
      data,
      list(
        intercept = prior("normal", list(0, 1)),
        m = prior("normal", list(0, 1))
      )
    ),
    "Matrix-valued predictor 'm' is not supported"
  )

  data$g <- factor(c("a", "b", "c"))
  expect_error(
    JAGS_formula(
      ~ g + m:g,
      "mu",
      data,
      list(
        intercept = prior("normal", list(0, 1)),
        g = prior_factor("normal", list(0, 1), contrast = "treatment"),
        "m:g" = prior_factor("normal", list(0, 1), contrast = "treatment")
      )
    ),
    "Matrix-valued predictor 'm' is not supported"
  )
})

test_that("formula design matrices reject non-finite predictor values", {
  prior_list <- list(
    intercept = prior("normal", list(0, 1)),
    x = prior("normal", list(0, 1))
  )
  expect_error(
    JAGS_formula(
      ~ x,
      "mu",
      data.frame(x = c(1, Inf, 2)),
      prior_list
    ),
    "Formula design matrix contains non-finite values",
    fixed = TRUE
  )

  fitted <- JAGS_formula(
    ~ x,
    "mu",
    data.frame(x = 1:3),
    prior_list
  )
  posterior <- coda::mcmc(matrix(
    c(0, 1),
    nrow = 1,
    dimnames = list(NULL, c("mu_intercept", "mu_x"))
  ))
  attr(posterior, "formula_design") <- list(mu = fitted$formula_design)
  expect_error(
    JAGS_evaluate_formula(
      posterior,
      ~ x,
      "mu",
      data.frame(x = Inf),
      fitted$prior_list
    ),
    "Formula design matrix contains non-finite values",
    fixed = TRUE
  )

  positive_prior <- prior(
    "normal",
    list(0, 1),
    truncation = list(lower = 0, upper = Inf)
  )
  expect_error(
    JAGS_formula(
      ~ 1 + diag(0 + x | id),
      "mu",
      data.frame(x = c(1, Inf), id = factor(c("a", "b"))),
      list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(
        id = random_block(sd = positive_prior)
      )
    ),
    "Random-effect block 'id' design matrix contains non-finite values",
    fixed = TRUE
  )
})

test_that("ordinary random slopes reject invalid values under na.pass", {
  withr::local_options(list(na.action = "na.pass"))
  positive_prior <- prior(
    "normal",
    list(0, 1),
    truncation = list(lower = 0, upper = Inf)
  )
  random_priors <- prior_random(
    id = random_block(sd = positive_prior)
  )

  for(missing_value in list(NA_real_, NaN)){
    expect_error(
      JAGS_formula(
        ~ 1 + diag(0 + x | id),
        "mu",
        data.frame(
          x = c(1, missing_value),
          id = factor(c("a", "b"))
        ),
        list(intercept = prior("normal", list(0, 1))),
        prior_random = random_priors
      ),
      paste0(
        "Random-effect block 'id' contains missing predictor values; ",
        "random-effect design matrices must have one row per data row."
      ),
      fixed = TRUE
    )
  }

  for(infinite_value in list(Inf, -Inf)){
    expect_error(
      JAGS_formula(
        ~ 1 + diag(0 + x | id),
        "mu",
        data.frame(
          x = c(1, infinite_value),
          id = factor(c("a", "b"))
        ),
        list(intercept = prior("normal", list(0, 1))),
        prior_random = random_priors
      ),
      "Random-effect block 'id' design matrix contains non-finite values.",
      fixed = TRUE
    )
  }
})

test_that("JAGS_formula validates JAGS parameter names", {
  expect_no_error(
    JAGS_formula(
      ~ x,
      "mu.x",
      data.frame(x = c(-1, 0, 1)),
      list(
        intercept = prior("normal", list(0, 1)),
        x = prior("normal", list(0, 1))
      )
    )
  )
  expect_error(
    JAGS_formula(
      ~ x,
      "bad-name",
      data.frame(x = c(-1, 0, 1)),
      list(
        intercept = prior("normal", list(0, 1)),
        x = prior("normal", list(0, 1))
      )
    ),
    "'parameter' must start with a letter"
  )
  expect_error(
    JAGS_formula(
      ~ x,
      "model",
      data.frame(x = c(-1, 0, 1)),
      list(
        intercept = prior("normal", list(0, 1)),
        x = prior("normal", list(0, 1))
      )
    ),
    "reserved JAGS keyword 'model'",
    fixed = TRUE
  )
  expect_error(
    JAGS_formula(
      ~ x,
      NA_character_,
      data.frame(x = c(-1, 0, 1)),
      list(
        intercept = prior("normal", list(0, 1)),
        x = prior("normal", list(0, 1))
      )
    ),
    "'parameter' argument cannot contain NA/NaN values",
    fixed = TRUE
  )
})

test_that("formula interfaces reject reserved tokens in categorical levels", {
  factor_prior <- prior_factor(
    "normal",
    list(0, 1),
    contrast = "treatment"
  )
  fixed_priors <- list(
    intercept = prior("normal", list(0, 1)),
    f = factor_prior
  )

  for(f in list(
    factor(c("a", "b__xXx__c")),
    c("a", "b__xXx__c"),
    factor("a", levels = c("a", "unused__xXx__level"))
  )){
    expect_error(
      JAGS_formula(
        ~ f,
        "mu",
        data.frame(f = f),
        fixed_priors
      ),
      "Factor predictor 'f' contains the internally reserved token '__xXx__'",
      fixed = TRUE
    )
  }

  valid_fixed <- JAGS_formula(
    ~ f,
    "mu",
    data.frame(f = factor(c("a", "b"))),
    fixed_priors
  )
  fixed_posterior <- coda::mcmc(matrix(
    c(0, 0),
    nrow = 1L,
    dimnames = list(NULL, c("mu_intercept", "mu_f"))
  ))
  expect_error(
    JAGS_evaluate_formula(
      fixed_posterior,
      ~ f,
      "mu",
      data.frame(f = "bad__xXx__level"),
      valid_fixed$prior_list
    ),
    "Factor predictor 'f' contains the internally reserved token '__xXx__'",
    fixed = TRUE
  )

  sd_prior <- prior(
    "normal",
    list(0, 1),
    truncation = list(lower = 0, upper = Inf)
  )
  expect_error(
    JAGS_formula(
      ~ 1 + diag(0 + f | id),
      "mu",
      data.frame(
        f = factor(c("a", "b__xRE_SUMMARY__c")),
        id = factor(c("one", "two"))
      ),
      list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(
        id = random_block(sd = sd_prior)
      )
    ),
    paste0(
      "Random-effect factor predictor 'f' contains the internally reserved ",
      "token '__xRE_SUMMARY__'"
    ),
    fixed = TRUE
  )

  valid_random <- JAGS_formula(
    ~ 1 + diag(0 + f | id),
    "mu",
    data.frame(
      f = factor(c("a", "b")),
      id = factor(c("one", "two"))
    ),
    list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(sd = sd_prior)
    )
  )
  expect_error(
    BayesTools:::.bt_random_effect_prediction_data(
      valid_random$formula_design$random_effects[[1L]],
      data.frame(
        f = "bad__xXx__level",
        id = "one"
      ),
      allow_new_groups = TRUE
    ),
    paste0(
      "Random-effect factor predictor 'f' contains the internally reserved ",
      "token '__xXx__'"
    ),
    fixed = TRUE
  )

  expect_error(
    JAGS_formula(
      ~ 1 + diag(1 | id),
      "mu",
      data.frame(id = factor(c("one", "two__xREx__group"))),
      list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(
        id = random_block(sd = sd_prior)
      )
    ),
    paste0(
      "Random-effect grouping variable 'id' contains the internally reserved ",
      "token '__xREx__' in level 'two__xREx__group'"
    ),
    fixed = TRUE
  )
})

test_that("formula interfaces reject bare language objects", {
  prior_list <- list(
    intercept = prior("normal", list(0, 1)),
    x = prior("normal", list(0, 1))
  )
  data <- data.frame(x = c(-1, 0, 1))

  for(invalid_formula in list(quote(x), quote(x + 1))){
    expect_error(
      JAGS_formula(invalid_formula, "mu", data, prior_list),
      "'formula' must be a formula",
      fixed = TRUE
    )
  }

  posterior <- coda::mcmc(matrix(
    c(0, 1),
    nrow = 1,
    dimnames = list(NULL, c("mu_intercept", "mu_x"))
  ))
  expect_error(
    JAGS_evaluate_formula(
      posterior,
      quote(x + 1),
      "mu",
      data,
      prior_list
    ),
    "'formula' must be a formula",
    fixed = TRUE
  )
})

test_that("JAGS_formula validates prior classes for every model term", {
  data <- data.frame(
    x = c(-1, 0, 1, 2),
    z = c(2, 1, 0, -1),
    g = factor(c("a", "b", "a", "b"))
  )

  expect_error(
    JAGS_formula(
      ~ x:g,
      "mu",
      data,
      list(
        intercept = prior("normal", list(0, 1)),
        "x:g" = prior("normal", list(0, 1))
      )
    ),
    "Unsupported prior distribution defined for 'x:g' factor variable",
    fixed = TRUE
  )
  expect_error(
    JAGS_formula(
      ~ x:z,
      "mu",
      data,
      list(
        intercept = prior("normal", list(0, 1)),
        "x:z" = prior_factor("normal", list(0, 1), contrast = "treatment")
      )
    ),
    "Unsupported prior distribution defined for 'x:z' continuous variable",
    fixed = TRUE
  )
  expect_error(
    JAGS_formula(
      ~ 1,
      "mu",
      data,
      list(
        intercept = prior_factor("normal", list(0, 1), contrast = "treatment")
      )
    ),
    "Unsupported prior distribution defined for 'intercept' continuous variable",
    fixed = TRUE
  )
})

test_that("JAGS_formula revalidates centered factor-prior metadata", {
  invalid_prior <- prior_factor(
    "mnormal",
    list(mean = 0, sd = 1),
    contrast = "orthonormal"
  )
  invalid_prior$parameters$mean <- 0.25

  expect_error(
    JAGS_formula(
      ~ g,
      "mu",
      data.frame(g = factor(c("a", "b", "a"))),
      list(
        intercept = prior("normal", list(0, 1)),
        g = invalid_prior
      )
    ),
    "prior_list[[\"g\"]]' mean-difference or orthonormal factor prior must be centered exactly at zero",
    fixed = TRUE
  )
})

test_that("JAGS_formula reserves the intercept predictor name", {
  expect_error(
    JAGS_formula(
      ~ intercept,
      "mu",
      data.frame(intercept = c(-1, 0, 1)),
      list(
        intercept = prior("normal", list(0, 1))
      )
    ),
    "predictor name 'intercept' is reserved"
  )
})

test_that("JAGS_formula uses a neutral point for log intercepts without an intercept", {
  formula <- ~ x - 1
  attr(formula, "log(intercept)") <- TRUE

  result <- JAGS_formula(
    formula,
    "mu",
    data.frame(x = c(-1, 0, 1)),
    list(x = prior("normal", list(0, 1)))
  )

  expect_true(is.prior.point(result$prior_list$mu_intercept))
  expect_equal(result$prior_list$mu_intercept$parameters$location, 1)
  expect_match(result$formula_syntax, "log\\(mu_intercept\\)", fixed = FALSE)
})

test_that("JAGS_formula canonicalizes prior_none to a zero point prior", {
  data <- data.frame(x = c(-1, 0, 1))
  result <- JAGS_formula(
    ~ x,
    "mu",
    data,
    list(
      intercept = prior_none(prior_weights = 2),
      x = prior_none(prior_weights = 3)
    )
  )

  expect_true(is.prior.point(result$prior_list$mu_intercept))
  expect_true(is.prior.point(result$prior_list$mu_x))
  expect_identical(result$prior_list$mu_intercept$parameters$location, 0)
  expect_identical(result$prior_list$mu_x$parameters$location, 0)
  expect_identical(.prior_model_weight(result$prior_list$mu_intercept), 2)
  expect_identical(.prior_model_weight(result$prior_list$mu_x), 3)

  default_result <- JAGS_formula(
    ~ x,
    "mu",
    data,
    list("__default_continuous" = prior_none())
  )
  expect_true(is.prior.point(default_result$prior_list$mu_intercept))
  expect_true(is.prior.point(default_result$prior_list$mu_x))

  log_formula <- ~ 1
  attr(log_formula, "log(intercept)") <- TRUE
  expect_error(
    JAGS_formula(
      log_formula,
      "mu",
      data,
      list(intercept = prior_none())
    ),
    "must have strictly positive support",
    fixed = TRUE
  )
})

test_that("fixed formulas reject unsupported calls before data lookup", {
  data <- data.frame(x = c(-1, 0, 1), z = c(1, 0, -1))
  priors <- list(
    intercept = prior("normal", list(0, 1)),
    x = prior("normal", list(0, 1)),
    z = prior("normal", list(0, 1)),
    "x:z" = prior("normal", list(0, 1))
  )
  invalid_formulas <- list(
    ~ .,
    ~ I(x^2),
    ~ stats::poly(x, 2),
    ~ scale(x),
    ~ offset(z),
    ~ custom_transform(x)
  )
  offending <- c(
    ".",
    "I(x^2)",
    "stats::poly(x, 2)",
    "scale(x)",
    "offset(z)",
    "custom_transform(x)"
  )

  for(i in seq_along(invalid_formulas)){
    expect_error(
      JAGS_formula(
        invalid_formulas[[i]],
        "mu",
        data,
        priors
      ),
      offending[[i]],
      fixed = TRUE
    )
  }

  expect_no_error(
    JAGS_formula(
      ~ x * z + expression(2),
      "mu",
      data,
      priors
    )
  )

  valid_result <- JAGS_formula(
    ~ x,
    "mu",
    data,
    priors[c("intercept", "x")]
  )
  fit <- coda::mcmc(matrix(
    c(0, 1),
    nrow = 1,
    dimnames = list(NULL, c("mu_intercept", "mu_x"))
  ))
  expect_error(
    JAGS_evaluate_formula(
      fit,
      ~ I(x^2),
      "mu",
      data,
      valid_result$prior_list
    ),
    "Unsupported fixed-formula call 'I(x^2)'",
    fixed = TRUE
  )
})

test_that("random formulas reject transformations and persist tuple ordering", {
  for(random_formula in list(
    ~ 1 + diag(1 + I(x^2) | g1),
    ~ 1 + diag(1 + scale(x) | g1),
    ~ 1 + diag(1 + custom_transform(x) | g1)
  )){
    expect_error(
      BayesTools:::.bt_parse_random_effects(random_formula),
      "Create the transformed predictor as an explicit data column",
      fixed = TRUE
    )
  }

  data <- data.frame(
    x = rep(c(-1, 1), 4),
    g1 = factor(rep(c("b", "a", "b", "a"), 2), levels = c("a", "b")),
    g2 = factor(rep(c("y", "y", "x", "x"), 2), levels = c("x", "y"))
  )
  sd_prior <- prior(
    "normal",
    list(0, 1),
    truncation = list(lower = 0, upper = Inf)
  )
  result <- JAGS_formula(
    ~ 1 + diag(1 + x | g1:g2),
    "mu",
    data,
    list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      g1_g2 = random_block(sd = sd_prior)
    )
  )
  random_term <- result$formula_design$random_effects[[1L]]
  expected_group <- interaction(
    data$g1,
    data$g2,
    drop = TRUE,
    sep = ":",
    lex.order = TRUE
  )

  expect_identical(random_term$group_levels, levels(expected_group))
  expect_identical(random_term$group_map, as.integer(expected_group))
  expect_identical(random_term$group_components, c("g1", "g2"))
  expect_identical(
    random_term$group_tuples,
    matrix(
      c("a", "x", "a", "y", "b", "x", "b", "y"),
      ncol = 2,
      byrow = TRUE,
      dimnames = list(NULL, c("g1", "g2"))
    )
  )
  expect_identical(
    unname(random_term$group_tuple_index[random_term$group_tuple_keys]),
    seq_along(random_term$group_tuple_keys)
  )

  shuffled <- data[c(4, 1, 3, 2), , drop = FALSE]
  prediction <- BayesTools:::.bt_random_effect_prediction_data(
    random_term,
    shuffled,
    group_data = shuffled
  )
  expected_shuffled <- interaction(
    shuffled$g1,
    shuffled$g2,
    drop = TRUE,
    sep = ":",
    lex.order = TRUE
  )
  expect_identical(
    prediction$group_map,
    match(as.character(expected_shuffled), random_term$group_levels)
  )

  collision_data <- data.frame(
    g1 = factor(c("a", "a:b"), levels = c("a", "a:b")),
    g2 = factor(c("b:c", "c"), levels = c("b:c", "c"))
  )
  collision_result <- JAGS_formula(
    ~ 1 + diag(1 | g1:g2),
    "mu",
    collision_data,
    list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      g1_g2 = random_block(sd = sd_prior)
    )
  )
  collision_term <- collision_result$formula_design$random_effects[[1L]]
  expect_identical(collision_term$n_groups, 2L)
  expect_true(anyDuplicated(collision_term$group_labels) > 0L)
  expect_false(anyDuplicated(collision_term$group_tuple_keys) > 0L)
  expect_false(anyDuplicated(collision_term$group_levels) > 0L)
  expect_identical(
    BayesTools:::.bt_random_group_unique_labels(
      "a:b:c",
      "3:a:b|1:c",
      existing = "a:b:c"
    ),
    "a:b:c [3:a:b|1:c]"
  )
})

test_that("log-intercept formulas require recursively positive prior support", {
  data <- data.frame(x = c(-1, 0, 1))
  log_formula <- ~ x
  attr(log_formula, "log(intercept)") <- TRUE

  invalid_priors <- list(
    prior("normal", list(0, 1)),
    prior("point", list(0)),
    prior("bernoulli", list(0.5)),
    prior_mixture(list(
      prior("gamma", list(2, 1)),
      prior("point", list(0))
    )),
    prior_spike_and_slab(prior("gamma", list(2, 1)))
  )
  for(intercept_prior in invalid_priors){
    expect_error(
      JAGS_formula(
        log_formula,
        "mu",
        data,
        list(
          intercept = intercept_prior,
          x = prior("normal", list(0, 1))
        )
      ),
      "must have strictly positive support",
      fixed = TRUE
    )
  }

  valid_priors <- list(
    prior("gamma", list(2, 1)),
    prior("normal", list(0, 1), truncation = list(lower = 0, upper = Inf)),
    prior("point", list(1)),
    prior_mixture(list(
      prior("gamma", list(2, 1)),
      prior("lognormal", list(0, 1))
    )),
    prior_spike_and_slab(
      prior("gamma", list(2, 1)),
      prior_inclusion = prior("point", list(1))
    )
  )
  for(intercept_prior in valid_priors){
    expect_no_error(
      JAGS_formula(
        log_formula,
        "mu",
        data,
        list(
          intercept = intercept_prior,
          x = prior("normal", list(0, 1))
        )
      )
    )
  }

  fitted <- coda::mcmc(matrix(
    c(1, 0),
    nrow = 1,
    dimnames = list(NULL, c("mu_intercept", "mu_x"))
  ))
  invalid_replay_priors <- JAGS_formula(
    log_formula,
    "mu",
    data,
    list(
      intercept = prior("point", list(1)),
      x = prior("normal", list(0, 1))
    )
  )$prior_list
  invalid_replay_priors$mu_intercept <- prior("point", list(0))
  attr(invalid_replay_priors$mu_intercept, "parameter") <- "mu"
  expect_error(
    JAGS_evaluate_formula(
      fitted,
      log_formula,
      "mu",
      data,
      invalid_replay_priors
    ),
    "must have strictly positive support",
    fixed = TRUE
  )
})

test_that("transform_prior_samples validates counts and seeds before sampling", {
  expect_error(
    transform_prior_samples(NULL, n_samples = NA_real_),
    "'n_samples' argument cannot contain NA/NaN values",
    fixed = TRUE
  )
  expect_error(
    transform_prior_samples(NULL, seed = NA_real_),
    "'seed' argument cannot contain NA/NaN values",
    fixed = TRUE
  )
  expect_error(
    transform_prior_samples(NULL, seed = .Machine$integer.max + 1),
    "'seed' must be equal or lower than",
    fixed = TRUE
  )
})

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

  attr(posterior, "formula_design") <- NULL
  legacy_prediction <- JAGS_evaluate_formula(
    posterior,
    ~ x:z,
    "mu",
    newdata,
    formula_result$prior_list
  )
  expect_equal(
    unname(drop(legacy_prediction)),
    unname(expected),
    tolerance = 1e-12
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

  attr(posterior, "formula_design") <- NULL
  legacy_prediction <- JAGS_evaluate_formula(
    posterior,
    ~ x:g,
    "mu",
    newdata,
    formula_result$prior_list
  )
  expect_equal(
    unname(drop(legacy_prediction)),
    unname(expected),
    tolerance = 1e-12
  )
})

test_that("formula expression terms are parsed structurally", {
  expect_equal(.extract_expressions(~ expression(log(x))), list("log(x)"))
  expect_equal(.extract_expressions(y ~ z + expression(log(x)) + expression(exp(b))), list("log(x)", "exp(b)"))
  expect_equal(.remove_expressions(y ~ expression(log(x))), formula(y ~ 1), ignore_formula_env = TRUE)
  expect_equal(.remove_expressions(~ z + expression(log(x))), formula(~ z), ignore_formula_env = TRUE)

  expression_result <- JAGS_formula(
    y ~ expression(sqrt(abs(log(x[i]))) + exp(0)),
    "mu",
    data.frame(x = c(1, 2, 3)),
    list(intercept = prior("normal", list(0, 1)))
  )
  expect_equal(
    expression_result$formula_design$transformed_terms,
    list("sqrt(abs(log(x[i]))) + exp(0)")
  )

  expect_error(
    JAGS_formula(
      ~ expression(step(x[i])),
      "mu",
      data.frame(x = c(-1, 1)),
      list(intercept = prior("normal", list(0, 1)))
    ),
    "is not replayable: unsupported call 'step'",
    fixed = TRUE
  )
  parameter_expression <- JAGS_formula(
    ~ expression(theta),
    "mu",
    data.frame(x = c(-1, 1)),
    list(intercept = prior("normal", list(0, 1)))
  )
  expect_identical(
    parameter_expression$formula_design$expression_specs[[1L]]$unresolved_dependencies,
    "theta"
  )
  expect_error(
    BayesTools:::.bt_formula_expression_finalize_design(
      design = parameter_expression$formula_design,
      formula_data = data.frame(x = c(-1, 1)),
      model_data = list(),
      parameter_names = character(),
      forbidden_parameters = "mu",
      context = "Test expression"
    ),
    "unknown replay dependency 'theta'",
    fixed = TRUE
  )
  expect_error(
    BayesTools:::.bt_formula_expression_specs(
      "theta[i, 1]",
      parameter_names = "theta"
    ),
    "sampled parameter 'theta' must use exactly one index",
    fixed = TRUE
  )
  expect_error(
    JAGS_formula(
      ~ expression(log(i)),
      "mu",
      data.frame(i = c(1, 2)),
      list(intercept = prior("normal", list(0, 1)))
    ),
    "cannot contain a column named 'i'",
    fixed = TRUE
  )

  expect_error(
    JAGS_formula(
      ~ x - expression(z[i]),
      "mu",
      data.frame(x = c(1, 2, 3), z = c(3, 2, 1)),
      list(
        intercept = prior("normal", list(0, 1)),
        x = prior("normal", list(0, 1))
      )
    ),
    "expression() terms must be additive formula terms",
    fixed = TRUE
  )
  expect_error(
    JAGS_formula(
      ~ expression(a[i]) - expression(b[i]),
      "mu",
      data.frame(a = c(1, 2, 3), b = c(3, 2, 1)),
      list(intercept = prior("normal", list(0, 1)))
    ),
    "expression() terms must be additive formula terms",
    fixed = TRUE
  )
})

test_that("JAGS_evaluate_formula includes literal expression contributions", {
  formula_result <- JAGS_formula(
    ~ x + expression(z[i]),
    "mu",
    data.frame(x = c(1, 2), z = c(10, 20)),
    list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    )
  )
  posterior <- coda::mcmc(matrix(
    c(1, 2),
    nrow = 1,
    dimnames = list(NULL, c("mu_intercept", "mu_x"))
  ))
  attr(posterior, "formula_design") <- list(mu = formula_result$formula_design)

  expect_equal(
    unname(drop(JAGS_evaluate_formula(posterior, formula = NULL, parameter = "mu"))),
    c(1 + 2 * 1 + 10, 1 + 2 * 2 + 20)
  )
  expect_equal(
    unname(drop(JAGS_evaluate_formula(
      posterior,
      ~ x + expression(z[i]),
      "mu",
      data.frame(x = 3, z = 30),
      formula_result$prior_list
    ))),
    1 + 2 * 3 + 30
  )
  expect_equal(
    unname(JAGS_evaluate_formula(
      posterior,
      ~ x,
      "mu",
      data.frame(x = 3),
      formula_result$prior_list
    )),
    matrix(7, nrow = 1)
  )
})

test_that("formula design reconstruction includes expression offsets for marglik", {
  formula_result <- JAGS_formula(
    ~ x + expression(z[i]),
    "mu",
    data.frame(x = c(1, 2), z = c(10, 20)),
    list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    )
  )
  samples <- list(
    mu_intercept = 1,
    mu_x = 2
  )
  reconstructed <- BayesTools:::.bt_JAGS_marglik_parameters_formula_design(
    samples = samples,
    design = formula_result$formula_design,
    formula_prior_list = formula_result$prior_list,
    prior_list_parameters = list(),
    log_intercept = FALSE
  )
  expect_equal(unname(reconstructed), c(1 + 2 * 1 + 10, 1 + 2 * 2 + 20))

  bridge_plan <- BayesTools:::.bt_JAGS_bridge_compile_formula_design_plan(
    design = formula_result$formula_design,
    formula_prior_list = formula_result$prior_list,
    log_intercept = FALSE
  )
  expect_equal(
    unname(bridge_plan$value(samples = samples, prior_list_parameters = list())),
    c(1 + 2 * 1 + 10, 1 + 2 * 2 + 20)
  )
})

test_that("formula expression helpers reject non-finite or wrong-length results", {
  expect_error(
    suppressWarnings(BayesTools:::.bt_formula_expression_row_values(
      expressions = list("log(-1)"),
      data = data.frame(x = 1),
      n_rows = 1L
    )),
    "non-finite"
  )
  expect_error(
    BayesTools:::.bt_formula_expression_row_values(
      expressions = list("x"),
      data = data.frame(x = c(1, 2)),
      n_rows = 1L
    ),
    "must be scalar or have 1 rows"
  )
  expect_error(
    BayesTools:::.bt_formula_expression_row_values(
      expressions = list("system('echo unsafe')"),
      data = data.frame(x = 1),
      n_rows = 1L
    ),
    "is not replayable: unsupported call 'system'",
    fixed = TRUE
  )
})

test_that("formula expressions replay sampled indexed parameters", {
  data <- data.frame(
    x = c(1, 2, 3),
    mapping = c(1L, 2L, 1L)
  )
  formula_result <- JAGS_formula(
    ~ x + expression(theta[mapping[i]]),
    "mu",
    data,
    list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    )
  )
  formula_result$formula_design <-
    BayesTools:::.bt_formula_expression_finalize_design(
      design = formula_result$formula_design,
      formula_data = data,
      model_data = list(),
      parameter_names = "theta",
      forbidden_parameters = "mu",
      context = "Test expression"
    )
  expect_identical(
    formula_result$formula_design$expression_specs[[1L]]$data_dependencies,
    "mapping"
  )
  expect_identical(
    formula_result$formula_design$expression_specs[[1L]]$parameter_dependencies,
    "theta"
  )
  expect_equal(formula_result$data$mapping, data$mapping)

  posterior <- coda::mcmc(matrix(
    c(
      1, 2, 10, 20,
      -1, 0.5, 30, 40
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      NULL,
      c("mu_intercept", "mu_x", "theta[1]", "theta[2]")
    )
  ))
  attr(posterior, "formula_design") <- list(
    mu = formula_result$formula_design
  )
  prediction <- JAGS_evaluate_formula(
    posterior,
    formula = NULL,
    parameter = "mu"
  )
  expect_equal(
    unname(prediction),
    matrix(
      c(13, 25, 17, 29.5, 40, 30.5),
      nrow = 3,
      ncol = 2
    )
  )
  new_prediction <- JAGS_evaluate_formula(
    posterior,
    formula = ~ x + expression(theta[mapping[i]]),
    parameter = "mu",
    data = data.frame(x = c(4, 5), mapping = c(2L, 1L)),
    prior_list = formula_result$prior_list
  )
  expect_equal(
    unname(new_prediction),
    matrix(c(29, 21, 41, 31.5), nrow = 2, ncol = 2)
  )
  expect_error(
    JAGS_evaluate_formula(
      posterior,
      formula = ~ x + expression(theta[mapping[i]]),
      parameter = "mu",
      data = data.frame(x = c(4, 5)),
      prior_list = formula_result$prior_list
    ),
    "unknown replay dependency 'mapping'",
    fixed = TRUE
  )

  reconstructed <- BayesTools:::.bt_JAGS_marglik_parameters_formula_design(
    samples = list(
      mu_intercept = 1,
      mu_x = 2,
      `theta[1]` = 10,
      `theta[2]` = 20
    ),
    design = formula_result$formula_design,
    formula_prior_list = formula_result$prior_list,
    prior_list_parameters = list(theta = c(10, 20)),
    log_intercept = FALSE
  )
  expect_equal(unname(reconstructed), c(13, 25, 17))

  bridge_plan <- BayesTools:::.bt_JAGS_bridge_compile_formula_design_plan(
    design = formula_result$formula_design,
    formula_prior_list = formula_result$prior_list,
    log_intercept = FALSE
  )
  expect_equal(
    unname(bridge_plan$value(
      samples = list(
        mu_intercept = 1,
        mu_x = 2,
        `theta[1]` = 10,
        `theta[2]` = 20
      ),
      prior_list_parameters = list(theta = c(10, 20))
    )),
    c(13, 25, 17)
  )
  expect_error(
    bridge_plan$value(
      samples = list(mu_intercept = 1, mu_x = 2),
      prior_list_parameters = list()
    ),
    "cannot reconstruct expression parameter 'theta'",
    fixed = TRUE
  )

  sparse_specs <- BayesTools:::.bt_formula_expression_specs(
    "theta[2]",
    parameter_names = "theta"
  )
  expect_equal(
    BayesTools:::.bt_formula_expression_contribution_matrix(
      expressions = sparse_specs,
      data = list(),
      n_rows = 3L,
      n_draws = 2L,
      samples = matrix(
        c(10, 20),
        ncol = 1L,
        dimnames = list(NULL, "theta[2]")
      )
    ),
    matrix(c(10, 10, 10, 20, 20, 20), nrow = 3L)
  )
  expect_error(
    BayesTools:::.bt_formula_expression_contribution_matrix(
      expressions = BayesTools:::.bt_formula_expression_specs(
        "theta[1]",
        parameter_names = "theta"
      ),
      data = list(),
      n_rows = 1L,
      n_draws = 1L,
      samples = matrix(
        20,
        ncol = 1L,
        dimnames = list(NULL, "theta[2]")
      )
    ),
    "produced non-finite values",
    fixed = TRUE
  )
})

test_that("sampled parameter expressions coexist with random-effect syntax", {
  sd_prior <- prior(
    "normal",
    list(0, 1),
    truncation = list(lower = 0, upper = Inf)
  )
  formula_result <- JAGS_formula(
    ~ expression(theta[mapping[i]]) + diag(1 | id),
    "mu",
    data.frame(
      mapping = c(1L, 2L, 1L),
      id = factor(c("a", "a", "b"))
    ),
    list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(sd = sd_prior)
    )
  )
  expect_length(formula_result$formula_design$random_effects, 1L)
  expect_identical(
    formula_result$formula_design$expression_specs[[1L]]$unresolved_dependencies,
    "theta"
  )
  expect_match(
    formula_result$formula_syntax,
    "theta[mapping[i]]",
    fixed = TRUE
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
  original_posterior <- transform_scale_samples(
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

test_that("nested grouping expands consistently across covariance specials and wrappers", {

  block_labels <- function(formula){
    terms <- BayesTools:::.bt_parse_random_effects(formula)$terms
    list(
      labels     = vapply(terms, function(t) t$group_label, character(1)),
      structures = vapply(terms, function(t) t$structure, character(1)),
      blocks     = vapply(terms, function(t) t$block_name, character(1))
    )
  }

  # The plain-bar reference expansion (already supported) defines the canonical
  # order/interaction convention every wrapped structure must match.
  reference <- block_labels(~ 1 + (1 | g1 / g2))
  expect_equal(reference$labels, c("g2:g1", "g1"))
  expect_equal(reference$blocks, c("g2_g1", "g1"))
  parenthesized_reference <- block_labels(~ 1 + random(1 | (g1 / g2)))
  expect_equal(parenthesized_reference$labels, c("g2:g1", "g1"))
  expect_equal(parenthesized_reference$blocks, c("g2_g1", "g1"))
  parenthesized_special <- block_labels(~ 1 + diag(1 | (g1 / g2)))
  expect_equal(parenthesized_special$labels, c("g2:g1", "g1"))
  expect_equal(parenthesized_special$blocks, c("g2_g1", "g1"))

  for(structure in c("diag", "id", "us", "un", "cs", "hcs", "ar1", "ar", "har", "car")){
    parsed <- block_labels(stats::as.formula(
      paste0("~ 1 + x + ", structure, "(1 + x | g1 / g2)")
    ))
    expect_equal(parsed$labels, c("g2:g1", "g1"),
                 info = paste("labels for", structure))
    expect_equal(parsed$blocks, c("g2_g1", "g1"),
                 info = paste("blocks for", structure))
    expect_equal(length(unique(parsed$structures)), 1L,
                 info = paste("single structure for", structure))
  }

  for(wrapper in c("random", "re")){
    parsed <- block_labels(stats::as.formula(
      paste0("~ 1 + ", wrapper, "(1 | g1 / g2)")
    ))
    expect_equal(parsed$labels, c("g2:g1", "g1"))
    expect_equal(parsed$blocks, c("g2_g1", "g1"))
    expect_equal(parsed$structures, c("us", "us"))
  }

  # Three-level nesting expands to one block per cumulative level.
  three <- block_labels(~ 1 + diag(1 | g1 / g2 / g3))
  expect_equal(three$labels, c("g3:g2:g1", "g2:g1", "g1"))
  expect_equal(three$blocks, c("g3_g2_g1", "g2_g1", "g1"))

  # Formula order is preserved when wrapped and plain terms are mixed.
  mixed <- block_labels(~ 1 + x + diag(1 | g1 / g2) + random(1 | h))
  expect_equal(mixed$blocks, c("g2_g1", "g1", "h"))
  expect_equal(mixed$structures, c("diag", "diag", "us"))

  # Named covariance / hom arguments propagate to every expanded sub-term.
  propagated <- block_labels(~ 1 + x + random(1 + x | g1 / g2, covariance = "cs"))
  expect_equal(propagated$structures, c("cs", "cs"))
  hom_propagated <- BayesTools:::.bt_parse_random_effects(~ 1 + cs(1 + x | g1 / g2, hom = TRUE))$terms
  expect_true(all(vapply(hom_propagated, function(t) isTRUE(t$hom), logical(1))))
  wrapper_hcs <- BayesTools:::.bt_parse_random_effects(
    ~ 1 + random(1 + x | g1, covariance = "cs", hom = FALSE)
  )$terms[[1]]
  expect_equal(wrapper_hcs$structure, "hcs")
  expect_null(wrapper_hcs$hom)
  wrapper_har <- BayesTools:::.bt_parse_random_effects(
    ~ 1 + re(1 + x | g1, covariance = "ar1", hom = FALSE)
  )$terms[[1]]
  expect_equal(wrapper_har$structure, "har")
  expect_null(wrapper_har$hom)
  expect_error(
    BayesTools:::.bt_parse_random_effects(~ 1 + random(cs(1 + x | g1))),
    "Do not wrap covariance-special random-effect calls",
    fixed = TRUE
  )
  expect_error(
    BayesTools:::.bt_parse_random_effects(~ 1 + re(ar1(1 + x | g1), name = "time")),
    "Do not wrap covariance-special random-effect calls",
    fixed = TRUE
  )
  expect_error(
    random_effects_formula(~ foo() + random(1 | g1)),
    "only random-effect terms",
    fixed = TRUE
  )

  # '||' under a wrapper expands and still resolves to diagonal blocks, matching
  # the plain-bar '||' nesting expansion.
  double_bar <- block_labels(~ 1 + x + random(1 + x || g1 / g2))
  expect_equal(double_bar$structures, c("diag", "diag"))
  expect_equal(double_bar$blocks, c("g2_g1", "g1"))

  # An explicit block name becomes a per-level prefix on each expanded block,
  # preserving the user's naming intent while keeping the blocks distinct.
  named_two <- block_labels(~ 1 + random(1 | g1 / g2, name = "foo"))
  expect_equal(named_two$blocks, c("foo_g2_g1", "foo_g1"))
  named_special <- block_labels(~ 1 + diag(1 | g1 / g2, name = "foo"))
  expect_equal(named_special$blocks, c("foo_g2_g1", "foo_g1"))
  named_three <- block_labels(~ 1 + random(1 | g1 / g2 / g3, name = "spatial"))
  expect_equal(named_three$blocks, c("spatial_g3_g2_g1", "spatial_g2_g1", "spatial_g1"))

  # End-to-end: numeric grouping IDs must build two nested blocks, not a single
  # block keyed on the literal 'g1 / g2' quotient.
  df_numeric <- data.frame(
    g1 = rep(1:4, each = 10),               # 4 outer groups
    g2 = rep(rep(1:2, each = 5), 4),        # 2 inner-of-g1 -> 8 (g2:g1) groups
    g3 = rep(1:5, 8)                        # 5 inner-of-(g2:g1) -> 40 (g3:g2:g1) groups
  )
  sd_prior <- prior("normal", list(0, 1), truncation = list(lower = 0, upper = Inf))
  map_level_counts <- function(formula){
    built <- JAGS_formula(
      formula,
      parameter    = "mu",
      data         = df_numeric,
      prior_list   = list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(sd = sd_prior)
    )
    map_columns <- grep("xRE_MAPx", names(built$data), value = TRUE)
    sort(vapply(map_columns, function(m) length(unique(built$data[[m]])), integer(1)))
  }

  # Two-level nesting: inner (g2:g1) has 8 groups, outer (g1) has 4.
  for(formula in list(
    ~ 1 + (1 | g1 / g2),
    ~ 1 + diag(1 | g1 / g2),
    ~ 1 + random(1 | g1 / g2),
    ~ 1 + diag(1 | (g1 / g2)),
    ~ 1 + random(1 | (g1 / g2))
  )){
    expect_equal(unname(map_level_counts(formula)), c(4L, 8L))
  }

  # Three-level nesting: 40 / 8 / 4 across the cumulative levels.
  for(formula in list(~ 1 + (1 | g1 / g2 / g3), ~ 1 + diag(1 | g1 / g2 / g3), ~ 1 + random(1 | g1 / g2 / g3))){
    expect_equal(unname(map_level_counts(formula)), c(4L, 8L, 40L))
  }
})

test_that("covariance specials preserve parenthesized bar terms", {

  diag_parsed <- BayesTools:::.bt_parse_random_effects(~ diag((1 | g)))
  expect_length(diag_parsed$terms, 1L)
  expect_equal(diag_parsed$terms[[1L]]$structure, "diag")
  expect_equal(diag_parsed$terms[[1L]]$term_formula, ~ 1, ignore_formula_env = TRUE)
  expect_equal(diag_parsed$fixed_formula, ~ 1, ignore_formula_env = TRUE)

  cs_parsed <- BayesTools:::.bt_parse_random_effects(~ cs((x | g)))
  expect_length(cs_parsed$terms, 1L)
  expect_equal(cs_parsed$terms[[1L]]$structure, "cs")
  expect_equal(cs_parsed$terms[[1L]]$term_formula, ~ x, ignore_formula_env = TRUE)
  expect_equal(cs_parsed$fixed_formula, ~ 1, ignore_formula_env = TRUE)

  named_diag <- BayesTools:::.bt_parse_random_effects(
    ~ 1 + diag((1 | g), name = "diag_block")
  )
  expect_equal(named_diag$terms[[1L]]$block_name, "diag_block")
  expect_equal(named_diag$terms[[1L]]$structure, "diag")
  expect_equal(named_diag$fixed_formula, ~ 1, ignore_formula_env = TRUE)

  named_cs <- BayesTools:::.bt_parse_random_effects(
    ~ 1 + cs((x | g), name = "cs_block", hom = TRUE)
  )
  expect_equal(named_cs$terms[[1L]]$block_name, "cs_block")
  expect_equal(named_cs$terms[[1L]]$structure, "cs")
  expect_true(named_cs$terms[[1L]]$hom)
  expect_equal(named_cs$fixed_formula, ~ 1, ignore_formula_env = TRUE)
})

test_that("random-effect formula lists retain component hierarchy", {

  df <- data.frame(
    district = factor(c("d1", "d1", "d2", "d2")),
    school   = factor(c("s1", "s2", "s1", "s2")),
    id       = factor(c("i1", "i1", "i2", "i2")),
    study    = factor(c("a", "a", "b", "b"))
  )
  sd_prior <- prior(
    "normal",
    list(mean = 0, sd = 1),
    truncation = list(lower = 0, upper = Inf)
  )

  random_effects <- random_effects_formula(
    list(
      nested = ~ 1 | district / school,
      study  = ~ 1 | study
    )
  )
  expect_equal(
    vapply(random_effects$terms, `[[`, character(1), "block_name"),
    c("nested_school_district", "nested_district", "study")
  )
  expect_equal(
    random_effects$components,
    list(nested = c("nested_school_district", "nested_district"), study = "study")
  )

  expect_equal(
    attr(random_effects$formula, "random_components"),
    random_effects$components
  )
  expect_equal(
    vapply(attr(random_effects$formula, "random_terms"), `[[`, character(1), "block_name"),
    c("nested_school_district", "nested_district", "study")
  )

  explicit_named_component <- random_effects_formula(
    list(study = ~ random(1 | id, name = "id"))
  )
  expect_equal(
    vapply(explicit_named_component$terms, `[[`, character(1), "block_name"),
    "id"
  )
  expect_equal(explicit_named_component$components, list(study = "id"))
  expect_equal(explicit_named_component$terms[[1]]$component_label, "study")
  expect_true(explicit_named_component$terms[[1]]$component_visible)

  explicit_named_result <- JAGS_formula(
    formula = explicit_named_component$formula,
    parameter = "mu",
    data = df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(id = random_block(sd = sd_prior))
  )
  expect_equal(explicit_named_result$formula_design$random_effects[[1]]$block_name, "id")
  expect_error(
    JAGS_formula(
      formula = explicit_named_component$formula,
      parameter = "mu",
      data = df,
      prior_list = list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(study = random_block(sd = sd_prior))
    ),
    "block override names were not found",
    fixed = TRUE
  )

  random_prior <- prior_random(
    random_variance_allocation(
      name  = "random_total",
      terms = c(nested = "nested", study = "study"),
      sd    = sd_prior
    ),
    random_variance_allocation(
      name   = "nested_split",
      terms  = c(school_district = "nested_school_district", district = "nested_district"),
      parent = allocation_ref("random_total", "nested")
    )
  )
  allocations <- BayesTools:::.bt_random_allocation_list(random_prior$allocation)
  expect_length(allocations, 2L)
  expect_equal(allocations[[1]]$terms, c(nested = "nested", study = "study"))
  expect_equal(allocations[[2]]$parent$allocation, "random_total")
  expect_equal(allocations[[2]]$parent$component, "nested")
  expect_equal(
    allocations[[2]]$terms,
    c(school_district = "nested_school_district", district = "nested_district")
  )

  result <- JAGS_formula(
    formula      = random_effects$formula,
    parameter    = "mu",
    data         = df,
    prior_list   = list(intercept = prior("normal", list(0, 1))),
    prior_random = random_prior
  )
  expect_match(result$formula_syntax, "random_total__component_nested_sd", fixed = TRUE)
  expect_equal(length(result$formula_design$random_effects[[1]]$sd_binding$allocations[[1L]]$factors), 2L)
  expect_equal(length(result$formula_design$random_effects[[3]]$sd_binding$allocations[[1L]]$factors), 1L)

  expect_error(
    random_effects_formula(
      list(
        `study-component` = ~ 1 | study,
        study_component   = ~ 1 | district
      )
    ),
    "unique after sanitization"
  )
  expect_error(
    random_effects_formula(
      ~ us(1 | study, name = "shared") + us(1 | district, name = "shared")
    ),
    "block names must be unique"
  )
  expect_error(
    random_effects_formula(~ x + (1 | study)),
    "only random-effect terms"
  )
})

test_that("random group covariance constructor validates and scales kernels", {

  K <- matrix(
    c(4, 1, 1,
      1, 9, 2,
      1, 2, 16),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("b", "a", "extra"), c("b", "a", "extra"))
  )
  none <- random_group_covariance(K, scale = "none")
  expect_s3_class(none, "random_group_covariance")
  expect_equal(none$covariance, K, ignore_attr = TRUE)
  expect_equal(none$scale, "none")
  expect_identical(none$triangle_source, "complete")

  lower <- matrix(
    c(2, .4, 0, 3),
    nrow = 2,
    dimnames = list(c("a", "b"), c("a", "b"))
  )
  lower_covariance <- random_group_covariance(lower, scale = "none")
  expect_identical(lower_covariance$triangle_source, "lower")
  expect_equal(
    lower_covariance$covariance,
    matrix(
      c(2, .4, .4, 3),
      nrow = 2,
      dimnames = dimnames(lower)
    ),
    ignore_attr = TRUE
  )

  upper <- t(lower)
  upper_covariance <- random_group_covariance(upper, scale = "none")
  expect_identical(upper_covariance$triangle_source, "upper")
  expect_equal(
    upper_covariance$covariance,
    lower_covariance$covariance,
    ignore_attr = TRUE
  )

  asymmetric <- matrix(
    c(2, .4, .4 + .Machine$double.eps, 3),
    nrow = 2,
    dimnames = dimnames(lower)
  )
  expect_error(
    random_group_covariance(asymmetric, scale = "none"),
    "must be exactly symmetric",
    fixed = TRUE
  )

  prepared_none <- BayesTools:::.bt_prepare_group_covariance_kernel(
    none,
    group_levels = c("a", "b"),
    block_name = "id"
  )
  expect_equal(prepared_none$levels, c("a", "b"))
  expect_equal(prepared_none$dropped_levels, "extra")
  expect_equal(prepared_none$kernel, K[c("a", "b"), c("a", "b")])
  expect_true(all(is.finite(prepared_none$precision)))
  expect_true(is.finite(prepared_none$log_det))

  K2 <- matrix(
    c(4, 1, 1, 9),
    nrow = 2,
    dimnames = list(c("a", "b"), c("a", "b"))
  )
  prepared_cor <- BayesTools:::.bt_prepare_group_covariance_kernel(
    random_group_covariance(K2, scale = "cor"),
    group_levels = c("a", "b"),
    block_name = "id"
  )
  expect_equal(prepared_cor$kernel, stats::cov2cor(K2), tolerance = 1e-12)

  prepared_cor0 <- BayesTools:::.bt_prepare_group_covariance_kernel(
    random_group_covariance(K2, scale = "cor0"),
    group_levels = c("a", "b"),
    block_name = "id"
  )
  R <- stats::cov2cor(K2)
  expect_equal(prepared_cor0$kernel, (R - min(R)) / (1 - min(R)), tolerance = 1e-12)

  prepared_cov0 <- BayesTools:::.bt_prepare_group_covariance_kernel(
    random_group_covariance(K2, scale = "cov0"),
    group_levels = c("a", "b"),
    block_name = "id"
  )
  expect_equal(prepared_cov0$kernel, K2 - min(K2), tolerance = 1e-12)

  expect_error(
    BayesTools:::.bt_prepare_group_covariance_kernel(
      random_group_covariance(K2, scale = "none"),
      group_levels = c("a", "missing"),
      block_name = "id"
    ),
    "missing fitted level",
    fixed = TRUE
  )
  singular <- matrix(
    c(1, 1, 1, 1),
    nrow = 2,
    dimnames = list(c("a", "b"), c("a", "b"))
  )
  expect_error(
    BayesTools:::.bt_prepare_group_covariance_kernel(
      random_group_covariance(singular, scale = "none"),
      group_levels = c("a", "b"),
      block_name = "id"
    ),
    "positive definite",
    fixed = TRUE
  )

  nonfinite_precision <- diag(c(1e-310, 1))
  dimnames(nonfinite_precision) <- list(c("a", "b"), c("a", "b"))
  expect_error(
    BayesTools:::.bt_prepare_group_covariance_kernel(
      random_group_covariance(nonfinite_precision, scale = "none"),
      group_levels = c("a", "b"),
      block_name = "id"
    ),
    "finite precision and log-determinant metadata",
    fixed = TRUE
  )

  expect_error(
    random_group_covariance(matrix(1, 1, 1), scale = "none"),
    "row and column names",
    fixed = TRUE
  )
  bad <- K2
  dimnames(bad) <- list(c("a", "a"), c("a", "b"))
  expect_error(
    random_group_covariance(bad, scale = "none"),
    "row names must be unique",
    fixed = TRUE
  )
})

test_that("random_effects_formula attaches known group covariance by block names", {

  study_kernel <- diag(2)
  dimnames(study_kernel) <- list(c("a", "b"), c("a", "b"))
  district_kernel <- diag(2)
  dimnames(district_kernel) <- list(c("d1", "d2"), c("d1", "d2"))

  by_block <- random_effects_formula(
    ~ random(1 | study, name = "study_block") +
      random(1 | district, name = "district_block"),
    group_covariance = list(
      study_block = random_group_covariance(study_kernel, scale = "none"),
      district_block = district_kernel
    )
  )
  expect_equal(
    vapply(by_block$terms, function(term) term$group_covariance$scale, character(1)),
    c("none", "cor")
  )
  expect_equal(names(by_block$group_covariance), c("study_block", "district_block"))

  by_group <- random_effects_formula(
    ~ random(1 | study, name = "study_block"),
    group_covariance = list(study = random_group_covariance(study_kernel))
  )
  expect_s3_class(by_group$terms[[1]]$group_covariance, "random_group_covariance")

  expect_error(
    random_effects_formula(
      ~ random(1 | study, name = "a") + random(1 | district, name = "b"),
      group_covariance = random_group_covariance(study_kernel)
    ),
    "single random-effect block",
    fixed = TRUE
  )
  expect_error(
    random_effects_formula(
      ~ random(1 | study, name = "a"),
      group_covariance = list(missing = random_group_covariance(study_kernel))
    ),
    "does not match",
    fixed = TRUE
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
  expect_error(
    BayesTools:::.bt_parse_random_effects(grouped_formula),
    "grouping expressions must be variables",
    fixed = TRUE
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
  expect_match(result$formula_syntax, "xRE_Zx\\[i,j\\] ~ dnorm\\(0, 1\\)")

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
  explicit_na_id <- factor(
    c("b", NA, "b", "c", "a", "c"),
    levels = c("a", "b", "c", NA),
    exclude = NULL
  )
  expect_false(anyNA(explicit_na_id))
  expect_true(anyNA(levels(explicit_na_id)))
  expect_error(
    JAGS_formula(
      formula = ~ 1 + diag(1 | id),
      parameter = "mu",
      data = transform(df, id = explicit_na_id),
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
    "grouping expressions must be variables",
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
    "grouping expressions must be variables",
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
    "grouping expressions must be variables",
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
  expect_s3_class(
    random_block(terms = list(intercept = sd_prior, x = random_block(sd = sd_prior))),
    "random_block"
  )
  expect_error(
    random_block(terms = list(sd_prior)),
    "must be named",
    fixed = TRUE
  )
  expect_error(
    random_block(terms = list(intercept = sd_prior, intercept = sd_prior)),
    "must be unique",
    fixed = TRUE
  )
  expect_error(
    random_block(terms = list(intercept = 1)),
    "must be a prior or random_block",
    fixed = TRUE
  )
  malformed_terms_block <- random_block(sd = sd_prior)
  malformed_terms_block$terms <- list(intercept = 1)
  expect_error(
    prior_random(id = malformed_terms_block),
    "must be a prior or random_block",
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

  allocation_prior <- random_variance_allocation(name = "allocation",
    terms = c("study", "drug"),
    sd = sd_prior,
    weights = prior("dirichlet", list(alpha = c(2, 3)))
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
    prior_random(`xRE_ALLOCx` = allocation_prior),
    "internally used",
    fixed = TRUE
  )
  expect_error(
    prior_random(allocation = stats::setNames(list(allocation_prior), "__xRE_SUMMARY__")),
    "internally used",
    fixed = TRUE
  )
  expect_error(
    random_variance_allocation(name = "allocation",
      terms = c("study", "study"),
      sd = sd_prior,
      weights = prior("dirichlet", list(alpha = c(1, 1)))
    ),
    "must be unique",
    fixed = TRUE
  )
  expect_error(
    random_variance_allocation(name = "allocation",
      terms = "study",
      sd = sd_prior
    ),
    "at least two resolved random-effect blocks",
    fixed = TRUE
  )
  expect_error(
    random_variance_allocation(name = "allocation",
      terms = c("study", "drug"),
      sd = sd_prior,
      weights = prior("beta", list(1, 1))
    ),
    "Dirichlet simplex priors",
    fixed = TRUE
  )
  expect_error(
    random_variance_allocation(name = "allocation",
      terms = c(study = "study", drug = "drug"),
      sd = sd_prior,
      inclusion = prior("spike", list(location = 0.5))
    ),
    "named list of scalar probability priors",
    fixed = TRUE
  )
  expect_error(
    random_variance_allocation(name = "allocation",
      terms = c(study = "study", drug = "drug"),
      sd = sd_prior,
      inclusion = list(study = 0.5)
    ),
    "must be a prior object",
    fixed = TRUE
  )
  expect_error(
    random_variance_allocation(name = "allocation",
      terms = c(study = "study", drug = "drug"),
      sd = sd_prior,
      inclusion = list(`bad-name` = prior("spike", list(location = 0.5)))
    ),
    "letters, numbers, and underscores",
    fixed = TRUE
  )
  expect_error(
    random_variance_allocation(name = "allocation",
      terms = c(study = "study", drug = "drug"),
      sd = sd_prior,
      inclusion = list(study = prior("spike", list(location = 1.5)))
    ),
    "bounded within 0 and 1",
    fixed = TRUE
  )
  expect_error(
    random_variance_allocation(name = "allocation",
      terms = c("study", "drug"),
      sd = sd_prior,
      weights = prior("dirichlet", list(alpha = c(1, 1, 1)))
    ),
    "dimension must match",
    fixed = TRUE
  )
  positional_allocation <- random_variance_allocation(
    "allocation",
    c("study", "drug"),
    sd_prior,
    prior("dirichlet", list(alpha = c(1, 1)))
  )

  unnamed_single <- random_effects_formula(list(~ 1 | study))
  expect_equal(
    vapply(unnamed_single$terms, `[[`, character(1), "block_name"),
    "study"
  )
  expect_equal(unnamed_single$components, list(component_1 = "study"))
  expect_equal(unnamed_single$terms[[1L]]$component, "component 1")
  expect_false(unnamed_single$terms[[1L]]$component_visible)

  bare_single <- random_effects_formula(~ 1 | study)
  expect_equal(bare_single$components, list(component_1 = "study"))
  expect_equal(bare_single$terms[[1L]]$component, "component 1")
  expect_false(bare_single$terms[[1L]]$component_visible)

  unnamed_multiple <- random_effects_formula(
    list(~ 1 | study, ~ 1 | id)
  )
  expect_equal(
    vapply(unnamed_multiple$terms, `[[`, character(1), "block_name"),
    c("component_1", "component_2")
  )
  expect_equal(
    vapply(unnamed_multiple$terms, `[[`, character(1), "component"),
    c("component 1", "component 2")
  )
  expect_true(all(vapply(
    unnamed_multiple$terms,
    `[[`,
    logical(1),
    "component_visible"
  )))
  expect_s3_class(positional_allocation$weights, "prior.simplex")
  positional_named_allocation <- random_variance_allocation(
    "total_re",
    c("study", "drug"),
    sd_prior,
    prior("dirichlet", list(alpha = c(1, 1)))
  )
  expect_equal(positional_named_allocation$name, "total_re")
  positional_child_allocation <- random_variance_allocation(
    "nested_split",
    c("paper", "estimate"),
    NULL,
    prior("dirichlet", list(alpha = c(1, 1))),
    allocation_ref("total_re", "study")
  )
  expect_equal(positional_child_allocation$name, "nested_split")
  expect_s3_class(positional_child_allocation$parent, "random_allocation_ref")
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
    random_variance_allocation(
      terms = c("study", "drug"),
      sd = sd_prior,
      name = "xRE_ALLOCx"
    ),
    "internally used",
    fixed = TRUE
  )
  expect_error(
    random_variance_allocation(
      terms = c("study", "drug"),
      sd = sd_prior,
      name = "__xRE_SUMMARY__"
    ),
    "internally used",
    fixed = TRUE
  )
  expect_error(
    random_variance_allocation(name = "allocation",
      terms = c(xRE_ALLOCx = "study", drug = "drug"),
      sd = sd_prior
    ),
    "internally used",
    fixed = TRUE
  )
  expect_error(
    allocation_ref("xRE_ALLOCx", "study"),
    "internally used",
    fixed = TRUE
  )
  expect_error(
    allocation_ref("total_re", "__xRE_SUMMARY__"),
    "internally used",
    fixed = TRUE
  )
  corrupted_ref <- allocation_ref("total_re", "study")
  corrupted_ref$component <- "__xRE_SUMMARY__"
  expect_error(
    BayesTools:::.bt_check_random_allocation_ref(corrupted_ref),
    "internally used",
    fixed = TRUE
  )
  corrupted_ref$component <- "bad-name"
  expect_error(
    BayesTools:::.bt_check_random_allocation_ref(corrupted_ref),
    "letters, numbers, and underscores",
    fixed = TRUE
  )
  expect_error(
    prior_random(sd = sd_prior, allocation = list(total = sd_prior)),
    "Variance allocation 'total' is invalid",
    fixed = TRUE
  )
  expect_error(
    prior_random(sd = sd_prior, allocation = list(allocation_prior, sd_prior)),
    "Variance allocation #2 is invalid",
    fixed = TRUE
  )
  expect_error(
    prior_random(sd = sd_prior, allocation = list(total = sd_prior)),
    "entries must be created with random_variance_allocation",
    fixed = TRUE
  )
  expect_false("allocation" %in% names(formals(random_block)))
  expect_false("allocation" %in% names(random_block()))
  expect_false("allocation" %in% names(
    BayesTools:::.bt_random_prior_for_block(prior_random(), "study")
  ))
  expect_s3_class(
    prior_random(sd = sd_prior, new_levels = random_new_levels(method = "zero")),
    "prior_random"
  )
  expect_error(
    random_covariance(cor = prior_factor("normal", list(0, 1), contrast = "treatment")),
    "ordinary scalar prior",
    fixed = TRUE
  )
  expect_error(
    random_covariance(cor = prior("mnormal", list(mean = 0, sd = 1, K = 2))),
    "ordinary scalar prior",
    fixed = TRUE
  )
  expect_error(
    random_covariance(cor = prior_none(), cor_scale = "logit"),
    "cannot use prior_none",
    fixed = TRUE
  )
  expect_error(
    random_covariance(
      cor = prior_mixture(list(
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
    "explicit 'cor'",
    fixed = TRUE
  )
  expect_error(
    random_covariance(eta = 2, cor = rho_prior),
    "explicit 'cor'",
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
  default_us <- JAGS_formula(
    formula = ~ 1 + x + (1 + x | id),
    parameter = "mu",
    data = df,
    prior_list = fixed_priors,
    prior_random = prior_random(id = random_block(sd = sd_prior))
  )
  expect_equal(default_us$formula_design$random_effects[[1]]$correlation$eta, 1)
  expect_true("mu__xREx__id_xRE_CORx_lkj_u[1]" %in% default_us$add_parameters)
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
      prior_random = prior_random(id = random_block(sd = sd_prior, cor = rho_prior))
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
        id = random_block(sd = sd_prior, cor = rho_prior)
      )
    ),
    "Block covariance override supplies a correlation prior, but structure 'id' has no correlation parameter",
    fixed = TRUE
  )
})

test_that("parameter and random SD sources validate simple external references", {

  sd_prior <- prior("gamma", list(2, 2))

  source <- parameter_source("tau", shape = "row")
  expect_s3_class(source, "parameter_source")
  expect_equal(source$shape, "row")
  expect_false(any(c("index", "expression") %in% names(source)))
  expect_equal(
    BayesTools:::.bt_parameter_source_jags_expression(source, row_index = "i"),
    "tau[i]"
  )
  expect_equal(
    BayesTools:::.bt_parameter_source_label(parameter_source("tau", shape = "row")),
    "tau[row]"
  )
  expect_equal(
    BayesTools:::.bt_parameter_source_row_names(source, 2),
    c("tau[1]", "tau[2]")
  )
  expect_error(
    BayesTools:::.bt_parameter_source_row_names(source, 1.5),
    "integer vector",
    fixed = TRUE
  )
  expect_error(
    BayesTools:::.bt_parameter_source_row_names(parameter_source("tau"), 2),
    "not row-shaped",
    fixed = TRUE
  )
  source_values <- parameter_source(
    "tau",
    shape = "row",
    values = function(parameters, data, n_rows){
      parameters$tau * data$scale[seq_len(n_rows)]
    }
  )
  expect_true(is.function(source_values$values))
  ellipsis_values <- parameter_source(
    "tau",
    shape = "row",
    values = function(...){
      args <- list(...)
      rep(args$parameters$tau, args$n_rows)
    }
  )
  expect_true(is.function(ellipsis_values$values))
  expect_error(
    parameter_source("tau", shape = "row", values = function(x)x),
    "must accept named arguments",
    fixed = TRUE
  )
  broken_values <- parameter_source("tau", shape = "row")
  broken_values$values <- function(x)x
  expect_error(
    BayesTools:::.bt_check_parameter_source(broken_values),
    "must accept named arguments",
    fixed = TRUE
  )

  sd_source <- random_sd_source(source)
  expect_s3_class(sd_source, "random_sd_source")
  expect_equal(sd_source$shape, "row")
  expect_false(any(c("index", "expression", "row_indexed", "values") %in% names(sd_source)))
  expect_true(BayesTools:::.bt_random_sd_source_is_row_indexed(sd_source))
  expect_equal(
    BayesTools:::.bt_random_sd_source_expression(sd_source, row_index = "i"),
    "tau[i]"
  )
  expect_error(
    BayesTools:::.bt_random_sd_source_expression(sd_source),
    "requires a row index",
    fixed = TRUE
  )
  expect_equal(
    BayesTools:::.bt_random_sd_source_expression(random_sd_source("tau")),
    "tau"
  )
  expect_true(is.function(BayesTools:::.bt_parameter_source_values_function(
    random_sd_source(source_values)
  )))
  source_posterior <- matrix(
    c(2, 4),
    ncol = 1,
    dimnames = list(NULL, "tau")
  )
  expect_equal(
    BayesTools:::.bt_parameter_source_value_draws(
      source = source_values,
      n_rows = 2,
      posterior = source_posterior,
      data = list(scale = c(1, 3))
    ),
    matrix(
      c(2, 6, 4, 12),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(NULL, c("tau[1]", "tau[2]"))
    ),
    tolerance = 1e-12
  )
  dependent_parameters <- BayesTools:::.bt_parameter_source_forbid_formula_parameters(
    list(theta = 2, deterministic_formula = 3, log_sigma = 4),
    "log_sigma"
  )
  allowed_source <- parameter_source(
    "allowed",
    shape = "row",
    values = function(parameters, data, n_rows){
      rep(parameters$theta + parameters[["deterministic_formula"]], n_rows)
    }
  )
  expect_equal(
    BayesTools:::.bt_parameter_source_value_draws(
      source = allowed_source,
      n_rows = 2,
      posterior = source_posterior,
      parameters = dependent_parameters
    ),
    matrix(
      rep(5, 4),
      nrow = 2,
      dimnames = list(NULL, c("allowed[1]", "allowed[2]"))
    )
  )
  forbidden_dollar_source <- parameter_source(
    "forbidden_dollar",
    shape = "row",
    values = function(parameters, data, n_rows){
      rep(parameters$log_sigma, n_rows)
    }
  )
  expect_error(
    BayesTools:::.bt_parameter_source_value_draws(
      source = forbidden_dollar_source,
      n_rows = 2,
      posterior = source_posterior,
      parameters = dependent_parameters
    ),
    "source 'forbidden_dollar\\[row\\]'.*formula parameter 'log_sigma'",
    perl = TRUE
  )
  forbidden_bracket_source <- parameter_source(
    "forbidden_bracket",
    shape = "row",
    values = function(parameters, data, n_rows){
      rep(parameters[["log_sigma"]], n_rows)
    }
  )
  expect_error(
    BayesTools:::.bt_parameter_source_value_draws(
      source = forbidden_bracket_source,
      n_rows = 2,
      posterior = source_posterior,
      parameters = dependent_parameters
    ),
    "source 'forbidden_bracket\\[row\\]'.*formula parameter 'log_sigma'",
    perl = TRUE
  )
  forbidden_subset_source <- parameter_source(
    "forbidden_subset",
    shape = "row",
    values = function(parameters, data, n_rows){
      rep(parameters["log_sigma"][[1]], n_rows)
    }
  )
  expect_error(
    BayesTools:::.bt_parameter_source_value_draws(
      source = forbidden_subset_source,
      n_rows = 2,
      posterior = source_posterior,
      parameters = dependent_parameters
    ),
    "source 'forbidden_subset\\[row\\]'.*formula parameter 'log_sigma'",
    perl = TRUE
  )
  expect_error(
    BayesTools:::.bt_parameter_source_value_draws(
      source = parameter_source("tau", shape = "row", values = function(parameters, data, n_rows) "bad"),
      n_rows = 2,
      posterior = source_posterior,
      context = "Test source"
    ),
    "must return a numeric vector",
    fixed = TRUE
  )
  expect_error(
    BayesTools:::.bt_parameter_source_value_draws(
      source = parameter_source("tau", shape = "row", values = function(parameters, data, n_rows) 1),
      n_rows = 2,
      posterior = source_posterior,
      context = "Test source"
    ),
    "must return a numeric vector of length 2",
    fixed = TRUE
  )
  expect_error(
    BayesTools:::.bt_parameter_source_value_draws(
      source = parameter_source("tau", shape = "row", values = function(parameters, data, n_rows) c(1, NA_real_)),
      n_rows = 2,
      posterior = source_posterior,
      context = "Test source"
    ),
    "returned missing values",
    fixed = TRUE
  )
  expect_error(
    BayesTools:::.bt_parameter_source_value_draws(
      source = parameter_source("tau", shape = "row", values = function(parameters, data, n_rows) stop("boom")),
      n_rows = 2,
      posterior = source_posterior,
      context = "Test source"
    ),
    "failed: boom",
    fixed = TRUE
  )
  source_data <- BayesTools:::.bt_JAGS_marglik_parameter_source_data(
    model_data = list(raw_x = c(1, 2)),
    formula_data = list(raw_x = c(1, 2), formula_only = c(3, 4)),
    design = structure(
      list(
        model_frame = data.frame(raw_x = c(-1, 1), design_only = c(5, 6)),
        source_data = data.frame(raw_x = c(1, 2), source_only = c(5, 6))
      ),
      class = "BayesTools_formula_design"
    )
  )
  expect_equal(source_data$raw_x, c(1, 2))
  expect_equal(source_data$source_only, c(5, 6))
  expect_false("formula_only" %in% names(source_data))
  expect_false("design_only" %in% names(source_data))
  expect_error(
    BayesTools:::.bt_JAGS_marglik_parameter_source_data(
      model_data = list(raw_x = c(9, 9)),
      formula_data = NULL,
      design = structure(
        list(source_data = data.frame(raw_x = c(1, 2))),
        class = "BayesTools_formula_design"
      )
    ),
    "conflict with the fitted source snapshot",
    fixed = TRUE
  )
  expect_error(
    BayesTools:::.bt_JAGS_marglik_parameter_source_data(
      model_data = list(raw_x = c(9, 9)),
      formula_data = list(raw_x = c(1, 2)),
      design = NULL
    ),
    "conflicting data for variable 'raw_x'",
    fixed = TRUE
  )
  scaled_source_result <- JAGS_formula(
    formula = ~ 1 + x + random(1 | id, name = "id", covariance = "diag"),
    parameter = "mu",
    data = data.frame(
      id = factor(c("a", "a", "b", "b")),
      x = c(10, 20, 30, 40)
    ),
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    ),
    formula_scale = list(x = TRUE),
    prior_random = prior_random(
      id = random_block(
        sd_source = random_sd_source(parameter_source(
          "tau",
          shape = "row",
          values = function(parameters, data, n_rows){
            data$x[seq_len(n_rows)]
          }
        ))
      )
    )
  )
  expect_equal(scaled_source_result$formula_design$source_data$x, c(10, 20, 30, 40))
  expect_false(isTRUE(all.equal(
    scaled_source_result$formula_design$model_frame$x,
    scaled_source_result$formula_design$source_data$x
  )))

  expect_error(
    parameter_source("1tau"),
    "must start with a letter",
    fixed = TRUE
  )
  expect_error(
    parameter_source("tau", shape = "rows"),
    "should be one of"
  )
  expect_error(
    parameter_source("tau", index = "i"),
    "unused argument",
    fixed = TRUE
  )
  expect_error(
    random_sd_source("tau", shape = "rows"),
    "should be one of"
  )
  expect_error(
    random_sd_source("tau", index = "i"),
    "unused argument",
    fixed = TRUE
  )
  expect_error(
    parameter_source("tau", values = 1),
    "'values' must be NULL or a function",
    fixed = TRUE
  )
  expect_error(
    parameter_source("tau", values = function(parameters, data, n_rows) 1),
    "only for row-shaped parameter sources",
    fixed = TRUE
  )
  expect_error(
    parameter_source("mu__xREx__study_intercept"),
    "internally used",
    fixed = TRUE
  )
  expect_error(
    random_sd_source("mu__xRE_ALLOCx_allocation__allocation_sd"),
    "internally used",
    fixed = TRUE
  )
  expect_error(
    random_sd_source("mu__xRE_SUMMARY__sd__study__intercept"),
    "internally used",
    fixed = TRUE
  )
  expect_error(
    random_sd_source(parameter_source("tau"), shape = "row"),
    "must not be supplied",
    fixed = TRUE
  )
  expect_error(
    BayesTools:::.bt_random_sd_binding(
      source = NULL,
      application = "block"
    ),
    "require a shared 'source'",
    fixed = TRUE
  )
  expect_error(
    BayesTools:::.bt_random_sd_binding(
      source = NULL,
      application = "column"
    ),
    "require shared or per-column sources",
    fixed = TRUE
  )
  expect_error(
    BayesTools:::.bt_random_sd_binding(
      sources_by_column = list(random_sd_source("tau", shape = "row")),
      application = "column"
    ),
    "Row-indexed per-column SD sources are not supported",
    fixed = TRUE
  )
  broken_sd_source <- BayesTools:::.bt_prior_owned_sd_source("tau")
  broken_sd_source$owned <- FALSE
  expect_error(
    BayesTools:::.bt_check_random_sd_binding_source(broken_sd_source),
    "inconsistent",
    fixed = TRUE
  )
  one_target_factor <- BayesTools:::.bt_random_variance_allocation_factor(
    weight_name = "w",
    index = 1L,
    scale = "total_variance",
    n_targets = 1L
  )
  expect_error(
    BayesTools:::.bt_check_random_variance_allocation_factor(one_target_factor),
    "n_targets",
    fixed = TRUE
  )
  broken_source <- parameter_source("tau")
  broken_source$shape <- "rows"
  expect_error(
    BayesTools:::.bt_check_parameter_source(broken_source),
    "metadata are inconsistent",
    fixed = TRUE
  )
  expect_error(
    random_variance_allocation(name = "allocation",
      terms = c("study", "drug"),
      sd = sd_prior,
      sd_source = random_sd_source("tau")
    ),
    "exactly one of 'sd' or 'sd_source'",
    fixed = TRUE
  )
  expect_error(
    random_variance_allocation(name = "allocation",
      terms = c("study", "drug")
    ),
    "exactly one of 'sd' or 'sd_source'",
    fixed = TRUE
  )
  expect_error(
    random_variance_allocation(name = "allocation",
      parent = allocation_ref("total", "study"),
      terms = c("paper", "estimate"),
      sd_source = random_sd_source("tau")
    ),
    "must not specify 'sd_source'",
    fixed = TRUE
  )
})

test_that("variance allocation priors generate shared total SD and Dirichlet allocation syntax", {

  df <- data.frame(
    study = factor(c("s1", "s1", "s2", "s2", "s3", "s3")),
    paper = factor(c("p1", "p2", "p1", "p2", "p3", "p3")),
    drug = factor(c("a", "b", "a", "b", "a", "b")),
    x = c(-1, 0, 1, 2, -2, 3)
  )
  fixed_priors <- list(intercept = prior("normal", list(0, 1)))
  random_prior <- prior_random(
    allocation = random_variance_allocation(name = "allocation",
      sd = prior("gamma", list(2, 2)),
      weights = prior("dirichlet", list(alpha = c(2, 3)))
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
      "mu__xRE_ALLOCx_allocation__allocation_sd",
      "mu__xRE_ALLOCx_allocation__weight"
    )
  )
  expect_s3_class(result$prior_list$mu__xRE_ALLOCx_allocation__weight, "prior.simplex")
  expect_match(
    result$formula_syntax,
    "mu__xREx__study_intercept = mu__xRE_ALLOCx_allocation__allocation_sd * sqrt(mu__xRE_ALLOCx_allocation__weight[1])",
    fixed = TRUE
  )
  expect_match(
    result$formula_syntax,
    "mu__xREx__drug_intercept = mu__xRE_ALLOCx_allocation__allocation_sd * sqrt(mu__xRE_ALLOCx_allocation__weight[2])",
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
    result$formula_design$random_effects[[1]]$sd_binding$allocations[[1L]]$weight_name,
    "mu__xRE_ALLOCx_allocation__weight"
  )
  expect_equal(result$formula_design$random_effects[[1]]$sd_binding$allocations[[1L]]$index, 1L)
  expect_equal(result$formula_design$random_effects[[2]]$sd_binding$allocations[[1L]]$index, 2L)
  expect_equal(
    result$formula_design$random_effects[[1]]$sd_parameter_names,
    "mu__xREx__study_intercept"
  )

  allocation_inclusion_result <- JAGS_formula(
    formula = ~ 1 +
      random(1 | study, name = "study", covariance = "diag") +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = df,
    prior_list = fixed_priors,
    prior_random = prior_random(
      allocation = random_variance_allocation(
        name = "total_re",
        terms = c(study = "study", drug = "drug"),
        sd = prior("gamma", list(2, 2)),
        weights = prior("dirichlet", list(alpha = c(2, 3))),
        inclusion = list(study = prior("spike", list(location = 0.5)))
      )
    )
  )
  expect_true("mu__xRE_ALLOCx_total_re__include_study_prob" %in%
                names(allocation_inclusion_result$prior_list))
  expect_match(
    allocation_inclusion_result$formula_syntax,
    "mu__xRE_ALLOCx_total_re__include_study_indicator ~ dbern(mu__xRE_ALLOCx_total_re__include_study_prob)",
    fixed = TRUE
  )
  expect_match(
    allocation_inclusion_result$formula_syntax,
    "mu__xREx__study_intercept = mu__xRE_ALLOCx_total_re__allocation_sd * mu__xRE_ALLOCx_total_re__include_study_indicator * sqrt(mu__xRE_ALLOCx_total_re__weight[1])",
    fixed = TRUE
  )
  expect_match(
    allocation_inclusion_result$formula_syntax,
    "mu__xREx__drug_intercept = mu__xRE_ALLOCx_total_re__allocation_sd * sqrt(mu__xRE_ALLOCx_total_re__weight[2])",
    fixed = TRUE
  )
  expect_false(grepl(
    "include_drug",
    allocation_inclusion_result$formula_syntax,
    fixed = TRUE
  ))
  expect_true("mu__xRE_ALLOCx_total_re__include_study_indicator" %in%
                allocation_inclusion_result$add_parameters)
  allocation_inclusion_syntax <- JAGS_add_priors(
    "model{}",
    allocation_inclusion_result$prior_list
  )
  expect_match(
    allocation_inclusion_syntax,
    "mu__xRE_ALLOCx_total_re__include_study_prob = 0.5",
    fixed = TRUE
  )

  allocation_inclusion_posterior <- matrix(
    c(
      2, 1, 3, 1,
      2, 1, 3, 0
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      NULL,
      c(
        "mu__xRE_ALLOCx_total_re__allocation_sd",
        "prior_par_eta_mu__xRE_ALLOCx_total_re__weight[1]",
        "prior_par_eta_mu__xRE_ALLOCx_total_re__weight[2]",
        "mu__xRE_ALLOCx_total_re__include_study_indicator"
      )
    )
  )
  inclusion_study_sd <- BayesTools:::.bt_random_effect_sd_draws(
    random_term = allocation_inclusion_result$formula_design$random_effects[[1]],
    n_columns = 1,
    posterior = allocation_inclusion_posterior,
    prior_list = allocation_inclusion_result$prior_list
  )
  inclusion_drug_sd <- BayesTools:::.bt_random_effect_sd_draws(
    random_term = allocation_inclusion_result$formula_design$random_effects[[2]],
    n_columns = 1,
    posterior = allocation_inclusion_posterior,
    prior_list = allocation_inclusion_result$prior_list
  )
  expect_equal(inclusion_study_sd[, 1], c(2 * sqrt(1 / 4), 0))
  expect_equal(inclusion_drug_sd[, 1], c(2 * sqrt(3 / 4), 2 * sqrt(3 / 4)))
  expect_error(
    BayesTools:::.bt_random_effect_sd_draws(
      random_term = allocation_inclusion_result$formula_design$random_effects[[1]],
      n_columns = 1,
      posterior = allocation_inclusion_posterior[
        ,
        colnames(allocation_inclusion_posterior) !=
          "mu__xRE_ALLOCx_total_re__include_study_indicator",
        drop = FALSE
      ],
      prior_list = allocation_inclusion_result$prior_list
    ),
    "missing Bernoulli indicator",
    fixed = TRUE
  )
  expect_equal(
    BayesTools:::.bt_JAGS_bridge_formula_allocation_inclusion_names(
      list(mu = allocation_inclusion_result$formula_design)
    ),
    "mu__xRE_ALLOCx_total_re__include_study_indicator"
  )
  expect_equal(
    BayesTools:::.bt_JAGS_bridge_allocation_factor_metadata(
      allocation_inclusion_result$formula_design$random_effects[[1]]$
        sd_binding$allocations[[1L]]$factors[[1L]]
    )$inclusion_name,
    "mu__xRE_ALLOCx_total_re__include_study_indicator"
  )
  expect_equal(
    BayesTools:::.bt_JAGS_bridge_allocation_metadata(
      allocation_inclusion_result$formula_design$random_effects[[1]]$
        sd_binding$allocations[[1L]]
    )$inclusion$study$indicator_name,
    "mu__xRE_ALLOCx_total_re__include_study_indicator"
  )
  expect_error(
    BayesTools:::.bt_JAGS_bridge_compile_formula_random_prior_evaluator(
      list(mu = allocation_inclusion_result$formula_design)
    ),
    "Bridge sampling for variance allocation inclusion gates",
    fixed = TRUE
  )
  expect_error(
    bayestools_reference_formula_random_log_prior(
      samples = allocation_inclusion_posterior[1, ],
      formula_design_list = list(mu = allocation_inclusion_result$formula_design)
    ),
    "Bridge sampling for variance allocation inclusion gates",
    fixed = TRUE
  )

  total_component_result <- JAGS_formula(
    formula = ~ 1 +
      random(1 | study, name = "study", covariance = "diag") +
      random(1 | paper, name = "paper", covariance = "diag") +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = df,
    prior_list = fixed_priors,
    prior_random = prior_random(
      random_variance_allocation(
        name = "a",
        terms = c(total = "nested", drug = "drug"),
        sd = prior("gamma", list(2, 2)),
        weights = prior("dirichlet", list(alpha = c(2, 3)))
      ),
      random_variance_allocation(
        name = "split",
        terms = c("study", "paper"),
        parent = allocation_ref("a", "total"),
        weights = prior("dirichlet", list(alpha = c(1, 1)))
      )
    )
  )
  expect_match(
    total_component_result$formula_syntax,
    "mu__xRE_ALLOCx_a__component_total_sd = mu__xRE_ALLOCx_a__allocation_sd * sqrt(mu__xRE_ALLOCx_a__weight[1])",
    fixed = TRUE
  )
  expect_false(grepl(
    "mu__xRE_ALLOCx_a__allocation_sd = mu__xRE_ALLOCx_a__allocation_sd",
    total_component_result$formula_syntax,
    fixed = TRUE
  ))

  malformed_binding <- result$formula_design$random_effects[[1]]$sd_binding
  malformed_binding$factors <- list()
  expect_error(
    BayesTools:::.bt_check_random_sd_binding(malformed_binding),
    "binding$factors",
    fixed = TRUE
  )

  reversed_result <- JAGS_formula(
    formula = ~ 1 +
      random(1 | study, name = "study", covariance = "diag") +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = df,
    prior_list = fixed_priors,
    prior_random = prior_random(
      allocation = random_variance_allocation(name = "allocation",
        terms = c("drug", "study"),
        sd = prior("gamma", list(2, 2)),
        weights = prior("dirichlet", list(alpha = c(2, 3)))
      )
    )
  )
  expect_equal(reversed_result$formula_design$random_effects[[1]]$sd_binding$allocations[[1L]]$index, 2L)
  expect_equal(reversed_result$formula_design$random_effects[[2]]$sd_binding$allocations[[1L]]$index, 1L)
  expect_match(
    reversed_result$formula_syntax,
    "mu__xREx__study_intercept = mu__xRE_ALLOCx_allocation__allocation_sd * sqrt(mu__xRE_ALLOCx_allocation__weight[2])",
    fixed = TRUE
  )
  expect_match(
    reversed_result$formula_syntax,
    "mu__xREx__drug_intercept = mu__xRE_ALLOCx_allocation__allocation_sd * sqrt(mu__xRE_ALLOCx_allocation__weight[1])",
    fixed = TRUE
  )

  full_syntax <- JAGS_add_priors("model{}", result$prior_list)
  expect_match(
    full_syntax,
    "prior_par_eta_mu__xRE_ALLOCx_allocation__weight[1] ~ dgamma(2, 1)",
    fixed = TRUE
  )
  expect_match(
    full_syntax,
    "mu__xRE_ALLOCx_allocation__weight[2] <- prior_par_eta_mu__xRE_ALLOCx_allocation__weight[2] / sum(prior_par_eta_mu__xRE_ALLOCx_allocation__weight[1:2])",
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
        "mu__xRE_ALLOCx_allocation__allocation_sd",
        "prior_par_eta_mu__xRE_ALLOCx_allocation__weight[1]",
        "prior_par_eta_mu__xRE_ALLOCx_allocation__weight[2]"
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
    "mu__xREx__study_f[1] = mu__xRE_ALLOCx_allocation__allocation_sd * sqrt(mu__xRE_ALLOCx_allocation__weight[1])",
    fixed = TRUE
  )

  random_bridge <- BayesTools:::.bt_JAGS_formula_random_bridge_parameters(
    list(mu = result$formula_design)
  )
  expect_true(all(grepl("_xRE_Zx", random_bridge$parameters, fixed = TRUE)))
  expect_false(any(grepl("_intercept", random_bridge$parameters, fixed = TRUE)))

  bridge_samples <- c(
    "mu_intercept" = 10,
    "mu__xRE_ALLOCx_allocation__allocation_sd" = 2,
    "prior_par_eta_mu__xRE_ALLOCx_allocation__weight[1]" = 1,
    "prior_par_eta_mu__xRE_ALLOCx_allocation__weight[2]" = 3,
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
  edge_bridge_samples[["prior_par_eta_mu__xRE_ALLOCx_allocation__weight[1]"]] <- 0
  expect_equal(
    JAGS_marglik_priors_formula(
      samples = edge_bridge_samples,
      formula_prior_list = list(mu = result$prior_list)
    ),
    -Inf
  )
  expect_error(
    JAGS_marglik_parameters_formula(
      samples                 = edge_bridge_samples,
      formula_list            = list(mu = result$formula),
      formula_data_list       = list(mu = result$data),
      formula_prior_list      = list(mu = result$prior_list),
      prior_list_parameters   = list(),
      formula_design_list     = list(mu = result$formula_design)
    ),
    "must be finite and positive",
    fixed = TRUE
  )
  expect_equal(
    bayestools_reference_formula_random_log_prior(
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
      samples = bridge_samples[names(bridge_samples) != "mu__xRE_ALLOCx_allocation__allocation_sd"],
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
        allocation = random_variance_allocation(name = "allocation",
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
        allocation = random_variance_allocation(name = "allocation",
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

test_that("bridge positive support is enforced before reconstruction", {

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
  invgamma_samples <- c("sigma" = 0)
  expect_equal(JAGS_marglik_priors(invgamma_samples, list(sigma = invgamma_prior)), -Inf)
  expect_error(
    JAGS_marglik_parameters(invgamma_samples, list(sigma = invgamma_prior)),
    "out-of-support inverse-gamma coordinate"
  )

  formula_prior_list <- list(
    mu = list(mu_intercept = invgamma_prior)
  )
  expect_equal(
    JAGS_marglik_priors_formula(
      samples = c("mu_intercept" = 0),
      formula_prior_list = formula_prior_list
    ),
    -Inf
  )
  expect_error(
    JAGS_marglik_parameters_formula(
      samples = c("mu_intercept" = 0),
      formula_list = list(mu = ~ 1),
      formula_data_list = list(mu = list(N_mu = 1)),
      formula_prior_list = formula_prior_list,
      prior_list_parameters = list()
    ),
    "out-of-support inverse-gamma coordinate"
  )

  # TODO(BayesTools 0.4.0): remove legacy inv_<parameter> inverse-gamma test.
  expect_error(
    JAGS_marglik_parameters(
      c("inv_sigma" = 0),
      list(sigma = invgamma_prior)
    ),
    "out-of-support positive auxiliary coordinate"
  )
})

test_that("external variance allocation sources generate scalar and row-indexed syntax", {

  df <- data.frame(
    study = factor(c("s1", "s1", "s2", "s2")),
    drug = factor(c("a", "b", "a", "b"))
  )
  fixed_priors <- list(intercept = prior("normal", list(0, 1)))

  direct_scalar_result <- JAGS_formula(
    formula = ~ 1 +
      random(1 | study, name = "study", covariance = "diag"),
    parameter = "mu",
    data = df,
    prior_list = fixed_priors,
    prior_random = prior_random(
      sd = prior("gamma", list(2, 2)),
      study = random_block(sd_source = random_sd_source("tau"))
    )
  )
  expect_equal(names(direct_scalar_result$prior_list), "mu_intercept")
  expect_match(
    direct_scalar_result$formula_syntax,
    "mu__xREx__study_xRE_STDx[1] = tau",
    fixed = TRUE
  )
  expect_equal(
    direct_scalar_result$formula_design$random_effects[[1]]$sd_binding$source$name,
    "tau"
  )
  expect_false(direct_scalar_result$formula_design$random_effects[[1]]$sd_binding$true_allocation)
  expect_true("tau" %in% direct_scalar_result$add_parameters)
  direct_scalar_coordinates <- BayesTools:::.bt_build_parameter_coordinates(
    columns = "tau",
    monitor_names = direct_scalar_result$add_parameters,
    prior_list = direct_scalar_result$prior_list,
    formula_design = list(mu = direct_scalar_result$formula_design)
  )
  expect_identical(direct_scalar_coordinates$role, "parameter")
  expect_identical(direct_scalar_coordinates$display_label, "tau")
  direct_scalar_summary <- .parameter_catalog_random_summary_samples(
    model_samples = matrix(
      2,
      nrow = 1,
      dimnames = list(NULL, "tau")
    ),
    prior_list = direct_scalar_result$prior_list,
    formula_design = list(mu = direct_scalar_result$formula_design)
  )
  expect_false(any(grepl("__xRE_SUMMARY__sd", colnames(direct_scalar_summary$model_samples), fixed = TRUE)))

  direct_row_result <- JAGS_formula(
    formula = ~ 1 +
      random(1 | study, name = "study", covariance = "diag"),
    parameter = "mu",
    data = df,
    prior_list = fixed_priors,
    prior_random = prior_random(
      study = random_block(sd_source = random_sd_source("tau", shape = "row"))
    )
  )
  expect_true(all(is.na(direct_row_result$formula_design$random_effects[[1]]$sd_parameter_names)))
  expect_match(
    direct_row_result$formula_syntax,
    "mu__xREx__study[i] = tau[i] * 1 * inprod(mu__xREx__study_xRE_UNIT_COEFx",
    fixed = TRUE
  )
  expect_equal(
    direct_row_result$formula_design$random_effects[[1]]$sd_binding$source$shape,
    "row"
  )

  expect_error(
    random_block(
      sd = prior("gamma", list(2, 2)),
      sd_source = random_sd_source("tau")
    ),
    "'sd_source' cannot be supplied together with block-local 'sd'",
    fixed = TRUE
  )

  scalar_result <- JAGS_formula(
    formula = ~ 1 +
      random(1 | study, name = "study", covariance = "diag") +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = df,
    prior_list = fixed_priors,
    prior_random = prior_random(
      allocation = random_variance_allocation(name = "allocation",
        sd_source = random_sd_source("tau"),
        weights = prior("dirichlet", list(alpha = c(2, 3)))
      )
    )
  )

  expect_equal(
    names(scalar_result$prior_list),
    c("mu_intercept", "mu__xRE_ALLOCx_allocation__weight")
  )
  expect_match(
    scalar_result$formula_syntax,
    "mu__xREx__study_intercept = tau * sqrt(mu__xRE_ALLOCx_allocation__weight[1])",
    fixed = TRUE
  )
  expect_match(
    scalar_result$formula_syntax,
    "mu__xREx__drug_intercept = tau * sqrt(mu__xRE_ALLOCx_allocation__weight[2])",
    fixed = TRUE
  )
  expect_equal(
    scalar_result$formula_design$random_effects[[1]]$sd_binding$allocations[[1L]]$source$name,
    "tau"
  )
  expect_equal(
    scalar_result$formula_design$random_effects[[1]]$sd_binding$source$shape,
    "scalar"
  )

  scalar_posterior <- matrix(
    c(
      2, 1, 3,
      4, 2, 2
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      NULL,
      c(
        "tau",
        "prior_par_eta_mu__xRE_ALLOCx_allocation__weight[1]",
        "prior_par_eta_mu__xRE_ALLOCx_allocation__weight[2]"
      )
    )
  )
  scalar_study_sd <- BayesTools:::.bt_random_effect_sd_draws(
    random_term = scalar_result$formula_design$random_effects[[1]],
    n_columns = 1,
    posterior = scalar_posterior,
    prior_list = scalar_result$prior_list
  )
  expect_equal(scalar_study_sd[, 1], c(2 * sqrt(1 / 4), 4 * sqrt(2 / 4)))
  scalar_monitored_leaf_posterior <- matrix(
    c(
      0.5, 1.5, 1, 3,
      0.75, 1.25, 2, 2
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      NULL,
      c(
        "mu__xREx__study_intercept",
        "mu__xREx__drug_intercept",
        "prior_par_eta_mu__xRE_ALLOCx_allocation__weight[1]",
        "prior_par_eta_mu__xRE_ALLOCx_allocation__weight[2]"
      )
    )
  )
  scalar_study_leaf_sd <- BayesTools:::.bt_random_effect_sd_draws(
    random_term = scalar_result$formula_design$random_effects[[1]],
    n_columns = 1,
    posterior = scalar_monitored_leaf_posterior,
    prior_list = scalar_result$prior_list
  )
  expect_equal(scalar_study_leaf_sd[, 1], c(0.5, 0.75))
  scalar_bridge_samples <- c(
    "tau" = 2,
    "prior_par_eta_mu__xRE_ALLOCx_allocation__weight[1]" = 1,
    "prior_par_eta_mu__xRE_ALLOCx_allocation__weight[2]" = 3,
    "mu__xREx__study_xRE_Zx[1,1]" = 0.1,
    "mu__xREx__study_xRE_Zx[2,1]" = 0.2
  )
  expect_equal(
    BayesTools:::.bt_JAGS_marglik_random_effect_value(
      samples = scalar_bridge_samples,
      random_term = scalar_result$formula_design$random_effects[[1]],
      prior_list = scalar_result$prior_list
    ),
    c(0.1, 0.1, 0.2, 0.2),
    tolerance = 1e-12
  )
  expect_equal(
    BayesTools:::.bt_JAGS_bridge_scale_metadata(scalar_result$formula_design$random_effects[[1]])$sd_binding$source$shape,
    "scalar"
  )

  scalar_sd_leaf_result <- JAGS_formula(
    formula = ~ 1 +
      random(1 + x | study, name = "study", covariance = "diag"),
    parameter = "mu",
    data = transform(df, x = c(-1, 0, 1, 2)),
    prior_list = fixed_priors,
    prior_random = prior_random(
      allocation = random_variance_allocation(name = "allocation",
        terms = "study",
        target = "sd_component",
        sd_source = random_sd_source("tau"),
        weights = prior("dirichlet", list(alpha = c(2, 3)))
      )
    )
  )
  expect_equal(
    names(scalar_sd_leaf_result$prior_list),
    c("mu_intercept", "mu__xRE_ALLOCx_allocation__weight")
  )
  expect_match(
    scalar_sd_leaf_result$formula_syntax,
    "mu__xREx__study_intercept = tau * sqrt(mu__xRE_ALLOCx_allocation__weight[1])",
    fixed = TRUE
  )
  expect_match(
    scalar_sd_leaf_result$formula_syntax,
    "mu__xREx__study_x = tau * sqrt(mu__xRE_ALLOCx_allocation__weight[2])",
    fixed = TRUE
  )
  scalar_sd_leaf_posterior <- matrix(
    c(2, 1, 3),
    nrow = 1,
    dimnames = list(
      NULL,
      c(
        "tau",
        "prior_par_eta_mu__xRE_ALLOCx_allocation__weight[1]",
        "prior_par_eta_mu__xRE_ALLOCx_allocation__weight[2]"
      )
    )
  )
  expect_equal(
    BayesTools:::.bt_random_effect_sd_draws(
      random_term = scalar_sd_leaf_result$formula_design$random_effects[[1]],
      n_columns = 2,
      posterior = scalar_sd_leaf_posterior,
      prior_list = scalar_sd_leaf_result$prior_list
    )[1, ],
    c(2 * sqrt(1 / 4), 2 * sqrt(3 / 4)),
    tolerance = 1e-12
  )

  row_result <- JAGS_formula(
    formula = ~ 1 +
      random(1 | study, name = "study", covariance = "diag") +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = df,
    prior_list = fixed_priors,
    prior_random = prior_random(
      allocation = random_variance_allocation(name = "allocation",
        sd_source = random_sd_source("tau", shape = "row"),
        weights = prior("dirichlet", list(alpha = c(2, 3)))
      )
    )
  )

  expect_equal(
    names(row_result$prior_list),
    c("mu_intercept", "mu__xRE_ALLOCx_allocation__weight")
  )
  expect_true(all(is.na(row_result$formula_design$random_effects[[1]]$sd_parameter_names)))
  expect_equal(
    row_result$formula_design$random_effects[[1]]$sd_binding$source$shape,
    "row"
  )
  expect_equal(
    BayesTools:::.bt_JAGS_bridge_scale_metadata(row_result$formula_design$random_effects[[1]])$sd_binding$source$shape,
    "row"
  )
  row_bridge_metadata <- BayesTools:::.bt_JAGS_bridge_scale_metadata(
    row_result$formula_design$random_effects[[1]]
  )
  expect_equal(row_bridge_metadata$sd_binding$source$shape, "row")
  expect_equal(row_bridge_metadata$sd_binding$source$source$shape, "row")
  expect_false(any(c("index", "expression", "row_indexed", "values") %in% names(row_bridge_metadata$sd_binding$source)))
  expect_equal(row_bridge_metadata$sd_binding$allocations[[1L]]$source$shape, "row")
  expect_false(any(c("source_index", "source_row_indexed") %in% names(row_bridge_metadata$sd_binding$allocations[[1L]])))
  expect_false(grepl("mu__xREx__study_intercept = tau[i]", row_result$formula_syntax, fixed = TRUE))
  expect_false(grepl("mu__xREx__study_xRE_STDx", row_result$formula_syntax, fixed = TRUE))
  expect_match(
    row_result$formula_syntax,
    "mu__xREx__study[i] = tau[i] * sqrt(mu__xRE_ALLOCx_allocation__weight[1]) * inprod(mu__xREx__study_xRE_UNIT_COEFx",
    fixed = TRUE
  )
  expect_equal(
    row_result$add_parameters,
    c("mu__xREx__study_xRE_Zx", "mu__xREx__drug_xRE_Zx")
  )
  row_id_result <- JAGS_formula(
    formula = ~ 1 +
      id(1 + x | study) +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = transform(df, x = c(-1, 0, 1, 2)),
    prior_list = fixed_priors,
    prior_random = prior_random(
      allocation = random_variance_allocation(name = "allocation",
        sd_source = random_sd_source("tau", shape = "row"),
        weights = prior("dirichlet", list(alpha = c(2, 3)))
      )
    )
  )
  expect_false(grepl("mu__xREx__study_xRE_STDx", row_id_result$formula_syntax, fixed = TRUE))
  expect_match(
    row_id_result$formula_syntax,
    "mu__xREx__study_xRE_UNIT_COEFx[1:2,i] = mu__xREx__study_xRE_Zx[1:2,i]",
    fixed = TRUE
  )
  expect_match(
    row_id_result$formula_syntax,
    "mu__xREx__study[i] = tau[i] * sqrt(mu__xRE_ALLOCx_allocation__weight[1]) * inprod(mu__xREx__study_xRE_UNIT_COEFx[mu__xREx__study_xRE_MAPx[i], 1:2]",
    fixed = TRUE
  )
  row_cs_result <- JAGS_formula(
    formula = ~ 1 +
      cs(f | study) +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = transform(df, f = factor(c("a", "b", "a", "b"))),
    prior_list = fixed_priors,
    prior_random = prior_random(
      allocation = random_variance_allocation(name = "allocation",
        sd_source = random_sd_source("tau", shape = "row"),
        weights = prior("dirichlet", list(alpha = c(2, 3)))
      ),
      study = random_block(cor = prior("normal", list(0, 0.5)))
    )
  )
  expect_false(grepl("mu__xREx__study_xRE_STDx", row_cs_result$formula_syntax, fixed = TRUE))
  expect_match(
    row_cs_result$formula_syntax,
    "mu__xREx__study_xRE_UNIT_COEFx[g,i] <- mu__xREx__study_xRE_CS_PREFIXx[g,i - 1] + mu__xREx__study_xRE_CS_DIAGx[i] * mu__xREx__study_xRE_Zx[g,i]",
    fixed = TRUE
  )
  expect_match(
    row_cs_result$formula_syntax,
    paste0(
      "mu__xREx__study[i] = tau[i] * ",
      "sqrt(mu__xRE_ALLOCx_allocation__weight[1]) * ",
      "mu__xREx__study_xRE_UNIT_COEFx[mu__xREx__study_xRE_MAPx[i],",
      "mu__xREx__study_xRE_COLx[i]]"
    ),
    fixed = TRUE
  )
  expect_null(row_cs_result$data$mu__xREx__study_xRE_DATAx)

  expect_error(
    BayesTools:::.bt_JAGS_marglik_random_effect_sd_values(
      samples = scalar_posterior[1, -1],
      random_term = row_result$formula_design$random_effects[[1]],
      prior_list = row_result$prior_list
    ),
    "scalar random-effect SD values are not defined for row-indexed sources",
    fixed = TRUE
  )
  row_bridge_samples <- c(
    "mu_intercept" = 10,
    "tau[1]" = 2,
    "tau[2]" = 4,
    "tau[3]" = 6,
    "tau[4]" = 8,
    "prior_par_eta_mu__xRE_ALLOCx_allocation__weight[1]" = 1,
    "prior_par_eta_mu__xRE_ALLOCx_allocation__weight[2]" = 3,
    "mu__xREx__study_xRE_Zx[1,1]" = 1,
    "mu__xREx__study_xRE_Zx[2,1]" = 2,
    "mu__xREx__drug_xRE_Zx[1,1]" = 3,
    "mu__xREx__drug_xRE_Zx[2,1]" = 4
  )
  row_bridge_posterior <- matrix(
    row_bridge_samples,
    nrow = 1,
    dimnames = list(NULL, names(row_bridge_samples))
  )
  attr(row_bridge_posterior, "lb") <- stats::setNames(
    rep(-Inf, ncol(row_bridge_posterior)),
    colnames(row_bridge_posterior)
  )
  attr(row_bridge_posterior, "ub") <- stats::setNames(
    rep(Inf, ncol(row_bridge_posterior)),
    colnames(row_bridge_posterior)
  )
  attr(row_bridge_posterior, "lb")[paste0("tau[", 1:4, "]")] <- 0
  expect_silent(BayesTools:::.bt_JAGS_bridge_check_row_indexed_external_sd_sources(
    formula_design_list = list(mu = row_result$formula_design),
    bridgesampling_posterior = row_bridge_posterior
  ))
  bad_row_bridge_posterior <- row_bridge_posterior
  attr(bad_row_bridge_posterior, "lb")["tau[2]"] <- -Inf
  expect_error(
    BayesTools:::.bt_JAGS_bridge_check_row_indexed_external_sd_sources(
      formula_design_list = list(mu = row_result$formula_design),
      bridgesampling_posterior = bad_row_bridge_posterior
    ),
    "non-negative lower bounds",
    fixed = TRUE
  )
  expect_equal(
    BayesTools:::.bt_JAGS_marglik_random_effect_value(
      samples = row_bridge_samples,
      random_term = row_result$formula_design$random_effects[[1]],
      prior_list = row_result$prior_list
    ),
    c(
      2 * sqrt(1 / 4) * 1,
      4 * sqrt(1 / 4) * 1,
      6 * sqrt(1 / 4) * 2,
      8 * sqrt(1 / 4) * 2
    ),
    tolerance = 1e-12
  )
  expect_equal(
    JAGS_marglik_parameters_formula(
      samples = row_bridge_samples,
      formula_list = list(mu = row_result$formula),
      formula_data_list = list(mu = row_result$data),
      formula_prior_list = list(mu = row_result$prior_list),
      prior_list_parameters = list(),
      formula_design_list = list(mu = row_result$formula_design)
    )$mu,
    10 + c(
      2 * (sqrt(1 / 4) * 1 + sqrt(3 / 4) * 3),
      4 * (sqrt(1 / 4) * 1 + sqrt(3 / 4) * 4),
      6 * (sqrt(1 / 4) * 2 + sqrt(3 / 4) * 3),
      8 * (sqrt(1 / 4) * 2 + sqrt(3 / 4) * 4)
    ),
    tolerance = 1e-12
  )
  value_source <- parameter_source(
    "tau",
    shape = "row",
    values = function(parameters, data, n_rows){
      parameters$total_tau * data$tau_factor[seq_len(n_rows)]
    }
  )
  value_prior_random <- prior_random(
    allocation = random_variance_allocation(name = "allocation",
      sd_source = random_sd_source(value_source),
      weights = prior("dirichlet", list(alpha = c(2, 3)))
    )
  )
  row_values_data <- transform(
    df,
    tau_factor = c(1, 2, 3, 4)
  )
  row_values_result <- JAGS_formula(
    formula = ~ 1 +
      random(1 | study, name = "study", covariance = "diag") +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = row_values_data,
    prior_list = fixed_priors,
    prior_random = value_prior_random
  )
  row_value_samples <- row_bridge_samples[
    !grepl("^tau\\[", names(row_bridge_samples))
  ]
  row_expected_value <- 10 + c(
    2 * (sqrt(1 / 4) * 1 + sqrt(3 / 4) * 3),
    4 * (sqrt(1 / 4) * 1 + sqrt(3 / 4) * 4),
    6 * (sqrt(1 / 4) * 2 + sqrt(3 / 4) * 3),
    8 * (sqrt(1 / 4) * 2 + sqrt(3 / 4) * 4)
  )
  expect_equal(
    JAGS_marglik_parameters_formula(
      samples = row_value_samples,
      formula_list = list(mu = row_values_result$formula),
      formula_data_list = list(mu = row_values_result$data),
      formula_prior_list = list(mu = row_values_result$prior_list),
      prior_list_parameters = list(total_tau = 2),
      formula_design_list = list(mu = row_values_result$formula_design),
      model_data = list(tau_factor = c(1, 2, 3, 4))
    )$mu,
    row_expected_value,
    tolerance = 1e-12
  )
  total_tau_result <- JAGS_formula(
    formula = ~ 1,
    parameter = "total_tau",
    data = df,
    prior_list = fixed_priors
  )
  expect_equal(
    JAGS_marglik_parameters_formula(
      samples = c(row_value_samples, "total_tau_intercept" = 2),
      formula_list = list(
        mu = row_values_result$formula,
        total_tau = total_tau_result$formula
      ),
      formula_data_list = list(
        mu = row_values_result$data,
        total_tau = total_tau_result$data
      ),
      formula_prior_list = list(
        mu = row_values_result$prior_list,
        total_tau = total_tau_result$prior_list
      ),
      prior_list_parameters = list(),
      formula_design_list = list(
        mu = row_values_result$formula_design,
        total_tau = total_tau_result$formula_design
      ),
      model_data = list(tau_factor = c(1, 2, 3, 4))
    )$mu,
    row_expected_value,
    tolerance = 1e-12
  )
  bridge_preflight_fit <- matrix(
    row_value_samples,
    nrow = 1,
    dimnames = list(NULL, names(row_value_samples))
  )
  bridge_preflight_fit <- bridge_preflight_fit[, names(row_value_samples) != "mu_intercept", drop = FALSE]
  bridge_preflight_fit <- cbind(
    "mu_intercept" = 0,
    bridge_preflight_fit
  )
  bridge_preflight_fit <- coda::mcmc(bridge_preflight_fit)
  attr(bridge_preflight_fit, "formula_design") <- list(
    mu = row_values_result$formula_design
  )
  bridge_value_context <- BayesTools:::.bt_JAGS_bridge_formula_context(
    fit = bridge_preflight_fit,
    formula_list = list(mu = ~ 1 +
      random(1 | study, name = "study", covariance = "diag") +
      random(1 | drug, name = "drug", covariance = "diag")),
    formula_data_list = list(mu = row_values_data),
    formula_prior_list = list(mu = fixed_priors),
    formula_scale_list = NULL,
    formula_random_prior_list = list(mu = value_prior_random)
  )
  expect_true(BayesTools:::.bt_parameter_source_has_values(
    BayesTools:::.bt_random_effect_row_indexed_source(
      bridge_value_context$formula_design_list$mu$random_effects[[1]]
    )
  ))
  expect_null(bridge_value_context$formula_data_list)
  expect_equal(
    JAGS_marglik_parameters_formula(
      samples = row_value_samples,
      formula_list = bridge_value_context$formula_list,
      formula_data_list = NULL,
      formula_prior_list = bridge_value_context$formula_prior_list,
      prior_list_parameters = list(total_tau = 2),
      formula_design_list = bridge_value_context$formula_design_list
    )$mu,
    row_expected_value,
    tolerance = 1e-12
  )
  ambiguous_row_value_samples <- c(
    row_value_samples,
    "tau[1]" = 999,
    "tau[2]" = 999,
    "tau[3]" = 999,
    "tau[4]" = 999
  )
  expect_error(
    BayesTools:::.bt_JAGS_marglik_random_effect_value(
      samples = ambiguous_row_value_samples,
      random_term = row_values_result$formula_design$random_effects[[1]],
      prior_list = row_values_result$prior_list
    ),
    "Use one row-source reconstruction path only",
    fixed = TRUE
  )
  ambiguous_bridge_posterior <- matrix(
    ambiguous_row_value_samples,
    nrow = 1,
    dimnames = list(NULL, names(ambiguous_row_value_samples))
  )
  expect_error(
    BayesTools:::.bt_JAGS_bridge_check_row_indexed_external_sd_sources(
      formula_design_list = list(mu = row_values_result$formula_design),
      bridgesampling_posterior = ambiguous_bridge_posterior
    ),
    "also appears as bridge parameter column",
    fixed = TRUE
  )
  bad_value_source <- parameter_source(
    "tau",
    shape = "row",
    values = function(parameters, data, n_rows){
      c(1, -1, 1, 1)[seq_len(n_rows)]
    }
  )
  bad_values_result <- JAGS_formula(
    formula = ~ 1 +
      random(1 | study, name = "study", covariance = "diag") +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = df,
    prior_list = fixed_priors,
    prior_random = prior_random(
      allocation = random_variance_allocation(name = "allocation",
        sd_source = random_sd_source(bad_value_source),
        weights = prior("dirichlet", list(alpha = c(2, 3)))
      )
    )
  )
  row_value_posterior <- matrix(
    row_value_samples,
    nrow = 1,
    dimnames = list(NULL, names(row_value_samples))
  )
  expect_error(
    BayesTools:::.bt_random_effect_row_indexed_contribution_from_latent(
      random_term = bad_values_result$formula_design$random_effects[[1]],
      model_matrix = bad_values_result$formula_design$random_effects[[1]]$model_matrix,
      group_map = bad_values_result$formula_design$random_effects[[1]]$group_map,
      posterior = row_value_posterior,
      prior_list = bad_values_result$prior_list,
      context = "Prediction"
    ),
    "non-finite or negative source values",
    fixed = TRUE
  )
  negative_marglik_error <- tryCatch(
    BayesTools:::.bt_JAGS_marglik_random_effect_value(
      samples = row_value_samples,
      random_term = bad_values_result$formula_design$random_effects[[1]],
      prior_list = bad_values_result$prior_list
    ),
    error = function(e)e
  )
  expect_s3_class(negative_marglik_error, "BayesTools_marglik_out_of_support")
  nan_value_source <- parameter_source(
    "tau",
    shape = "row",
    values = function(parameters, data, n_rows){
      c(1, NaN, 1, 1)[seq_len(n_rows)]
    }
  )
  nan_values_result <- JAGS_formula(
    formula = ~ 1 +
      random(1 | study, name = "study", covariance = "diag") +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = df,
    prior_list = fixed_priors,
    prior_random = prior_random(
      allocation = random_variance_allocation(name = "allocation",
        sd_source = random_sd_source(nan_value_source),
        weights = prior("dirichlet", list(alpha = c(2, 3)))
      )
    )
  )
  nan_marglik_error <- tryCatch(
    BayesTools:::.bt_JAGS_marglik_random_effect_value(
      samples = row_value_samples,
      random_term = nan_values_result$formula_design$random_effects[[1]],
      prior_list = nan_values_result$prior_list
    ),
    error = function(e)e
  )
  expect_s3_class(nan_marglik_error, "BayesTools_marglik_out_of_support")
  short_value_source <- parameter_source(
    "tau",
    shape = "row",
    values = function(parameters, data, n_rows){
      1
    }
  )
  short_values_result <- JAGS_formula(
    formula = ~ 1 +
      random(1 | study, name = "study", covariance = "diag") +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = df,
    prior_list = fixed_priors,
    prior_random = prior_random(
      allocation = random_variance_allocation(name = "allocation",
        sd_source = random_sd_source(short_value_source),
        weights = prior("dirichlet", list(alpha = c(2, 3)))
      )
    )
  )
  expect_error(
    BayesTools:::.bt_JAGS_marglik_random_effect_value(
      samples = row_value_samples,
      random_term = short_values_result$formula_design$random_effects[[1]],
      prior_list = short_values_result$prior_list
    ),
    "must return a numeric vector of length 4",
    fixed = TRUE
  )
  missing_source_fit <- bridge_preflight_fit
  attr(missing_source_fit, "formula_design") <- list(
    mu = row_result$formula_design
  )
  expect_error(
    JAGS_bridgesampling(
      fit = missing_source_fit,
      log_posterior = function(parameters, data) 0,
      data = list()
    ),
    "cannot reconstruct row-indexed external SD source 'tau[row]'",
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
        allocation = random_variance_allocation(name = "allocation",
          sd_source = random_sd_source("tau", shape = "row"),
          weights = prior("dirichlet", list(alpha = c(2, 3)))
        ),
        monitor = random_monitor(coefficients = TRUE)
      )
    ),
    "Group-level coefficient monitoring is not supported",
    fixed = TRUE
  )

  row_nested_result <- JAGS_formula(
    formula = ~ 1 +
      random(1 | study, name = "study", covariance = "diag") +
      random(1 | paper, name = "paper", covariance = "diag") +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = transform(df, paper = factor(c("p1", "p2", "p1", "p2"))),
    prior_list = fixed_priors,
    prior_random = prior_random(
      random_variance_allocation(
        name = "total_re",
        terms = c(nested = "nested", drug = "drug"),
        sd_source = random_sd_source("tau", shape = "row"),
        weights = prior("dirichlet", list(alpha = c(1, 3)))
      ),
      random_variance_allocation(
        name = "nested_split",
        parent = allocation_ref("total_re", "nested"),
        terms = c(study = "study", paper = "paper"),
        weights = prior("dirichlet", list(alpha = c(3, 2)))
      )
    )
  )
  expect_false(grepl(
    "mu__xRE_ALLOCx_total_re__component_nested_sd = tau[i]",
    row_nested_result$formula_syntax,
    fixed = TRUE
  ))
  expect_match(
    row_nested_result$formula_syntax,
    "mu__xREx__study[i] = tau[i] * sqrt(mu__xRE_ALLOCx_total_re__weight[1]) * sqrt(mu__xRE_ALLOCx_nested_split__weight[1]) * inprod(mu__xREx__study_xRE_UNIT_COEFx",
    fixed = TRUE
  )
  expect_match(
    row_nested_result$formula_syntax,
    "mu__xREx__drug[i] = tau[i] * sqrt(mu__xRE_ALLOCx_total_re__weight[2]) * inprod(mu__xREx__drug_xRE_UNIT_COEFx",
    fixed = TRUE
  )
  expect_equal(length(row_nested_result$formula_design$random_effects[[1]]$sd_binding$factors), 2L)
  expect_equal(length(row_nested_result$formula_design$random_effects[[3]]$sd_binding$factors), 1L)
  row_nested_samples <- c(
    "tau[1]" = 2,
    "tau[2]" = 4,
    "tau[3]" = 6,
    "tau[4]" = 8,
    "prior_par_eta_mu__xRE_ALLOCx_total_re__weight[1]" = 1,
    "prior_par_eta_mu__xRE_ALLOCx_total_re__weight[2]" = 3,
    "prior_par_eta_mu__xRE_ALLOCx_nested_split__weight[1]" = 3,
    "prior_par_eta_mu__xRE_ALLOCx_nested_split__weight[2]" = 2,
    "mu__xREx__study_xRE_Zx[1,1]" = 1,
    "mu__xREx__study_xRE_Zx[2,1]" = 2,
    "mu__xREx__paper_xRE_Zx[1,1]" = 3,
    "mu__xREx__paper_xRE_Zx[2,1]" = 4,
    "mu__xREx__drug_xRE_Zx[1,1]" = 5,
    "mu__xREx__drug_xRE_Zx[2,1]" = 6
  )
  expect_equal(
    BayesTools:::.bt_JAGS_marglik_random_effect_value(
      samples = row_nested_samples,
      random_term = row_nested_result$formula_design$random_effects[[1]],
      prior_list = row_nested_result$prior_list
    ),
    c(
      2 * sqrt(1 / 4) * sqrt(3 / 5) * 1,
      4 * sqrt(1 / 4) * sqrt(3 / 5) * 1,
      6 * sqrt(1 / 4) * sqrt(3 / 5) * 2,
      8 * sqrt(1 / 4) * sqrt(3 / 5) * 2
    ),
    tolerance = 1e-12
  )
  expect_equal(
    BayesTools:::.bt_JAGS_marglik_random_effect_value(
      samples = row_nested_samples,
      random_term = row_nested_result$formula_design$random_effects[[2]],
      prior_list = row_nested_result$prior_list
    ),
    c(
      2 * sqrt(1 / 4) * sqrt(2 / 5) * 3,
      4 * sqrt(1 / 4) * sqrt(2 / 5) * 4,
      6 * sqrt(1 / 4) * sqrt(2 / 5) * 3,
      8 * sqrt(1 / 4) * sqrt(2 / 5) * 4
    ),
    tolerance = 1e-12
  )
  expect_equal(
    BayesTools:::.bt_JAGS_marglik_random_effect_value(
      samples = row_nested_samples,
      random_term = row_nested_result$formula_design$random_effects[[3]],
      prior_list = row_nested_result$prior_list
    ),
    c(
      2 * sqrt(3 / 4) * 5,
      4 * sqrt(3 / 4) * 6,
      6 * sqrt(3 / 4) * 5,
      8 * sqrt(3 / 4) * 6
    ),
    tolerance = 1e-12
  )

  row_sd_leaf_result <- JAGS_formula(
    formula = ~ 1 +
      random(1 + x | study, name = "study", covariance = "diag"),
    parameter = "mu",
    data = transform(df, x = c(-1, 0, 1, 2)),
    prior_list = fixed_priors,
    prior_random = prior_random(
      allocation = random_variance_allocation(name = "allocation",
        terms = "study",
        target = "sd_component",
        sd_source = random_sd_source("tau", shape = "row"),
        weights = prior("dirichlet", list(alpha = c(1, 1)))
      )
    )
  )
  expect_true(all(is.na(row_sd_leaf_result$formula_design$random_effects[[1]]$sd_parameter_names)))
  expect_equal(row_sd_leaf_result$formula_design$random_effects[[1]]$sd_binding$application, "column")
  expect_equal(length(row_sd_leaf_result$formula_design$random_effects[[1]]$sd_binding$factors_by_column), 2L)
  expect_false(grepl("mu__xREx__study_intercept = tau[i]", row_sd_leaf_result$formula_syntax, fixed = TRUE))
  expect_match(
    row_sd_leaf_result$formula_syntax,
    "mu__xREx__study[i] = tau[i] * (sqrt(mu__xRE_ALLOCx_allocation__weight[1]) * mu__xREx__study_xRE_UNIT_COEFx[mu__xREx__study_xRE_MAPx[i],1] * mu__xREx__study_xRE_DATAx[i,1] + sqrt(mu__xRE_ALLOCx_allocation__weight[2]) * mu__xREx__study_xRE_UNIT_COEFx[mu__xREx__study_xRE_MAPx[i],2] * mu__xREx__study_xRE_DATAx[i,2])",
    fixed = TRUE
  )
  row_sd_leaf_forced_latent_result <- JAGS_formula(
    formula = ~ 1 +
      random(1 + x | study, name = "study", covariance = "diag"),
    parameter = "mu",
    data = transform(df, x = c(-1, 0, 1, 2)),
    prior_list = fixed_priors,
    prior_random = prior_random(
      allocation = random_variance_allocation(name = "allocation",
        terms = "study",
        target = "sd_component",
        sd_source = random_sd_source("tau", shape = "row"),
        weights = prior("dirichlet", list(alpha = c(1, 1)))
      ),
      study = random_block(
        monitor = random_monitor(latent = FALSE, coefficients = FALSE)
      )
    )
  )
  expect_true("mu__xREx__study_xRE_Zx" %in%
                row_sd_leaf_forced_latent_result$add_parameters)
  expect_false(any(grepl(
    "mu__xREx__study_xRE_COEFx",
    row_sd_leaf_forced_latent_result$add_parameters,
    fixed = TRUE
  )))
  row_sd_leaf_samples <- c(
    "tau[1]" = 2,
    "tau[2]" = 4,
    "tau[3]" = 6,
    "tau[4]" = 8,
    "prior_par_eta_mu__xRE_ALLOCx_allocation__weight[1]" = 1,
    "prior_par_eta_mu__xRE_ALLOCx_allocation__weight[2]" = 3,
    "mu__xREx__study_xRE_Zx[1,1]" = 1,
    "mu__xREx__study_xRE_Zx[2,1]" = 2,
    "mu__xREx__study_xRE_Zx[1,2]" = 3,
    "mu__xREx__study_xRE_Zx[2,2]" = 4
  )
  row_sd_leaf_expected <- c(
    2 * (sqrt(1 / 4) * 1 + sqrt(3 / 4) * 3 * -1),
    4 * (sqrt(1 / 4) * 1 + sqrt(3 / 4) * 3 * 0),
    6 * (sqrt(1 / 4) * 2 + sqrt(3 / 4) * 4 * 1),
    8 * (sqrt(1 / 4) * 2 + sqrt(3 / 4) * 4 * 2)
  )
  expect_equal(
    BayesTools:::.bt_JAGS_marglik_random_effect_value(
      samples = row_sd_leaf_samples,
      random_term = row_sd_leaf_result$formula_design$random_effects[[1]],
      prior_list = row_sd_leaf_result$prior_list
    ),
    row_sd_leaf_expected,
    tolerance = 1e-12
  )
  row_sd_leaf_ambiguous_samples <- c(
    row_sd_leaf_samples,
    "mu__xRE_ALLOCx_allocation__weight[1]" = 0.25,
    "mu__xRE_ALLOCx_allocation__weight[2]" = 0.75
  )
  expect_error(
    BayesTools:::.bt_JAGS_marglik_random_effect_value(
      samples = row_sd_leaf_ambiguous_samples,
      random_term = row_sd_leaf_result$formula_design$random_effects[[1]],
      prior_list = row_sd_leaf_result$prior_list
    ),
    "both normalized Dirichlet allocation coordinates",
    fixed = TRUE
  )
  expect_equal(
    unname(BayesTools:::.bt_random_effect_row_indexed_contribution_from_latent(
      random_term = row_sd_leaf_result$formula_design$random_effects[[1]],
      model_matrix = row_sd_leaf_result$formula_design$random_effects[[1]]$model_matrix,
      group_map = row_sd_leaf_result$formula_design$random_effects[[1]]$group_map,
      posterior = matrix(row_sd_leaf_samples, nrow = 1, dimnames = list(NULL, names(row_sd_leaf_samples))),
      prior_list = row_sd_leaf_result$prior_list
    )[, 1]),
    row_sd_leaf_expected,
    tolerance = 1e-12
  )
  row_sd_leaf_bad_eta <- row_sd_leaf_samples
  row_sd_leaf_bad_eta[["prior_par_eta_mu__xRE_ALLOCx_allocation__weight[1]"]] <- 0
  row_sd_leaf_bad_eta[["prior_par_eta_mu__xRE_ALLOCx_allocation__weight[2]"]] <- 0
  row_sd_leaf_bad_eta_error <- tryCatch(
    BayesTools:::.bt_JAGS_marglik_random_effect_value(
      samples = row_sd_leaf_bad_eta,
      random_term = row_sd_leaf_result$formula_design$random_effects[[1]],
      prior_list = row_sd_leaf_result$prior_list
    ),
    error = function(e)e
  )
  expect_s3_class(row_sd_leaf_bad_eta_error, "BayesTools_marglik_out_of_support")
  valid_weight_draws <- matrix(
    c(0.25, 0.75),
    nrow = 1,
    dimnames = list(NULL, c("w[1]", "w[2]"))
  )
  expect_equal(
    BayesTools:::.bt_random_effect_dirichlet_draws(
      parameter_name = "w",
      posterior = valid_weight_draws,
      prior_list = list(w = prior("dirichlet", list(alpha = c(1, 1))))
    ),
    valid_weight_draws
  )
  invalid_weight_error <- tryCatch(
    BayesTools:::.bt_random_effect_dirichlet_draws(
      parameter_name = "w",
      posterior = matrix(
        c(0.25, 0.90),
        nrow = 1,
        dimnames = list(NULL, c("w[1]", "w[2]"))
      ),
      prior_list = list(w = prior("dirichlet", list(alpha = c(1, 1))))
    ),
    error = function(e)e
  )
  expect_s3_class(
    invalid_weight_error,
    "BayesTools_random_effect_allocation_out_of_support"
  )
  bad_row_sd_leaf_target <- row_sd_leaf_result$formula_design$random_effects[[1]]
  bad_row_sd_leaf_target$sd_binding$allocations[[1L]]$target <- "block"
  expect_error(
    BayesTools:::.bt_random_effect_row_indexed_column_allocation_draws(
      random_term = bad_row_sd_leaf_target,
      posterior = matrix(
        row_sd_leaf_samples,
        nrow = 1,
        dimnames = list(NULL, names(row_sd_leaf_samples))
      ),
      prior_list = row_sd_leaf_result$prior_list,
      n_columns = 2L
    ),
    "missing 'allocation\\$factors'"
  )
  missing_allocation_term <- row_sd_leaf_result$formula_design$random_effects[[1]]
  missing_allocation_term$sd_binding$allocations <- list()
  expect_error(
    BayesTools:::.bt_random_effect_row_indexed_column_allocation_draws(
      random_term = missing_allocation_term,
      posterior = matrix(
        row_sd_leaf_samples,
        nrow = 1,
        dimnames = list(NULL, names(row_sd_leaf_samples))
      ),
      prior_list = row_sd_leaf_result$prior_list,
      n_columns = 2L
    ),
    "exactly one allocation record",
    fixed = TRUE
  )
  extra_allocation_binding <- row_sd_leaf_result$formula_design$random_effects[[1]]$sd_binding
  extra_allocation_binding$allocations <- c(
    extra_allocation_binding$allocations,
    extra_allocation_binding$allocations
  )
  expect_error(
    BayesTools:::.bt_check_random_sd_binding(extra_allocation_binding),
    "exactly one allocation record",
    fixed = TRUE
  )
  stale_allocation_binding <- row_sd_leaf_result$formula_design$random_effects[[1]]$sd_binding
  stale_allocation_binding$true_allocation <- FALSE
  expect_error(
    BayesTools:::.bt_check_random_sd_binding(stale_allocation_binding),
    "must not contain allocation records",
    fixed = TRUE
  )
  bad_row_sd_leaf_term <- row_sd_leaf_result$formula_design$random_effects[[1]]
  bad_row_sd_leaf_term$sd_binding$application <- "block"
  expect_error(
    BayesTools:::.bt_random_effect_row_indexed_column_allocation_draws(
      random_term = bad_row_sd_leaf_term,
      posterior = matrix(
        row_sd_leaf_samples,
        nrow = 1,
        dimnames = list(NULL, names(row_sd_leaf_samples))
      ),
      prior_list = row_sd_leaf_result$prior_list,
      n_columns = 2L
    ),
    "Block SD bindings must not use 'factors_by_column'",
    fixed = TRUE
  )
  empty_column_chain_term <- row_sd_leaf_result$formula_design$random_effects[[1]]
  empty_column_chain_term$sd_binding$factors_by_column <- list()
  expect_error(
    BayesTools:::.bt_random_effect_row_indexed_column_allocation_draws(
      random_term = empty_column_chain_term,
      posterior = matrix(
        row_sd_leaf_samples,
        nrow = 1,
        dimnames = list(NULL, names(row_sd_leaf_samples))
      ),
      prior_list = row_sd_leaf_result$prior_list,
      n_columns = 2L
    ),
    "missing canonical 'binding\\$factors_by_column'"
  )
  expect_error(
    BayesTools:::.bt_random_effect_row_indexed_contribution_from_latent(
      random_term = empty_column_chain_term,
      model_matrix = empty_column_chain_term$model_matrix,
      group_map = empty_column_chain_term$group_map,
      posterior = matrix(
        row_sd_leaf_samples,
        nrow = 1,
        dimnames = list(NULL, names(row_sd_leaf_samples))
      ),
      prior_list = row_sd_leaf_result$prior_list
    ),
    "missing canonical 'binding\\$factors_by_column'"
  )
  mismatched_column_chain_term <- row_sd_leaf_result$formula_design$random_effects[[1]]
  mismatched_column_chain_term$sd_binding$factors_by_column[[1L]][[1L]]$index <- 2L
  expect_error(
    BayesTools:::.bt_random_effect_row_indexed_column_allocation_draws(
      random_term = mismatched_column_chain_term,
      posterior = matrix(
        row_sd_leaf_samples,
        nrow = 1,
        dimnames = list(NULL, names(row_sd_leaf_samples))
      ),
      prior_list = row_sd_leaf_result$prior_list,
      n_columns = 2L
    ),
    "column factor chains do not match",
    fixed = TRUE
  )
  row_sd_leaf_mean_result <- JAGS_formula(
    formula = ~ 1 +
      random(1 + x | study, name = "study", covariance = "diag"),
    parameter = "mu",
    data = transform(df, x = c(-1, 0, 1, 2)),
    prior_list = fixed_priors,
    prior_random = prior_random(
      allocation = random_variance_allocation(name = "allocation",
        terms = "study",
        target = "sd_component",
        scale = "mean_variance",
        sd_source = random_sd_source("tau", shape = "row"),
        weights = prior("dirichlet", list(alpha = c(1, 1)))
      )
    )
  )
  expect_equal(
    BayesTools:::.bt_JAGS_marglik_random_effect_value(
      samples = row_sd_leaf_samples,
      random_term = row_sd_leaf_mean_result$formula_design$random_effects[[1]],
      prior_list = row_sd_leaf_mean_result$prior_list
    ),
    c(
      2 * (sqrt(2 * 1 / 4) * 1 + sqrt(2 * 3 / 4) * 3 * -1),
      4 * (sqrt(2 * 1 / 4) * 1 + sqrt(2 * 3 / 4) * 3 * 0),
      6 * (sqrt(2 * 1 / 4) * 2 + sqrt(2 * 3 / 4) * 4 * 1),
      8 * (sqrt(2 * 1 / 4) * 2 + sqrt(2 * 3 / 4) * 4 * 2)
    ),
    tolerance = 1e-12
  )

  row_child_sd_leaf_result <- JAGS_formula(
    formula = ~ 1 +
      random(1 + x | study, name = "het", covariance = "diag") +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = transform(df, x = c(-1, 0, 1, 2)),
    prior_list = fixed_priors,
    prior_random = prior_random(
      random_variance_allocation(
        name = "total_re",
        terms = c(het = "het", simple = "drug"),
        sd_source = random_sd_source("tau", shape = "row"),
        weights = prior("dirichlet", list(alpha = c(1, 3)))
      ),
      random_variance_allocation(
        name = "het_sd",
        parent = allocation_ref("total_re", "het"),
        terms = "het",
        target = "sd_component",
        weights = prior("dirichlet", list(alpha = c(1, 3)))
      )
    )
  )
  expect_false(grepl("mu__xRE_ALLOCx_total_re__component_het_sd = tau[i]", row_child_sd_leaf_result$formula_syntax, fixed = TRUE))
  expect_match(
    row_child_sd_leaf_result$formula_syntax,
    "mu__xREx__het[i] = tau[i] * (sqrt(mu__xRE_ALLOCx_total_re__weight[1]) * sqrt(mu__xRE_ALLOCx_het_sd__weight[1]) * mu__xREx__het_xRE_UNIT_COEFx[mu__xREx__het_xRE_MAPx[i],1]",
    fixed = TRUE
  )
  expect_equal(length(row_child_sd_leaf_result$formula_design$random_effects[[1]]$sd_binding$factors_by_column[[1L]]), 2L)
  row_child_sd_leaf_samples <- c(
    "tau[1]" = 2,
    "tau[2]" = 4,
    "tau[3]" = 6,
    "tau[4]" = 8,
    "prior_par_eta_mu__xRE_ALLOCx_total_re__weight[1]" = 1,
    "prior_par_eta_mu__xRE_ALLOCx_total_re__weight[2]" = 3,
    "prior_par_eta_mu__xRE_ALLOCx_het_sd__weight[1]" = 1,
    "prior_par_eta_mu__xRE_ALLOCx_het_sd__weight[2]" = 3,
    "mu__xREx__het_xRE_Zx[1,1]" = 1,
    "mu__xREx__het_xRE_Zx[2,1]" = 2,
    "mu__xREx__het_xRE_Zx[1,2]" = 3,
    "mu__xREx__het_xRE_Zx[2,2]" = 4,
    "mu__xREx__drug_xRE_Zx[1,1]" = 5,
    "mu__xREx__drug_xRE_Zx[2,1]" = 6
  )
  row_child_sd_leaf_expected <- c(
    2 * sqrt(1 / 4) * (sqrt(1 / 4) * 1 + sqrt(3 / 4) * 3 * -1),
    4 * sqrt(1 / 4) * (sqrt(1 / 4) * 1 + sqrt(3 / 4) * 3 * 0),
    6 * sqrt(1 / 4) * (sqrt(1 / 4) * 2 + sqrt(3 / 4) * 4 * 1),
    8 * sqrt(1 / 4) * (sqrt(1 / 4) * 2 + sqrt(3 / 4) * 4 * 2)
  )
  expect_equal(
    BayesTools:::.bt_JAGS_marglik_random_effect_value(
      samples = row_child_sd_leaf_samples,
      random_term = row_child_sd_leaf_result$formula_design$random_effects[[1]],
      prior_list = row_child_sd_leaf_result$prior_list
    ),
    row_child_sd_leaf_expected,
    tolerance = 1e-12
  )
  expect_equal(
    BayesTools:::.bt_JAGS_marglik_random_effect_value(
      samples = row_child_sd_leaf_samples,
      random_term = row_child_sd_leaf_result$formula_design$random_effects[[2]],
      prior_list = row_child_sd_leaf_result$prior_list
    ),
    c(
      2 * sqrt(3 / 4) * 5,
      4 * sqrt(3 / 4) * 6,
      6 * sqrt(3 / 4) * 5,
      8 * sqrt(3 / 4) * 6
    ),
    tolerance = 1e-12
  )
  row_child_bridge_design <- list(mu = row_child_sd_leaf_result$formula_design)
  row_child_bridge_changed <- row_child_bridge_design
  row_child_bridge_changed_term <- row_child_bridge_changed$mu$random_effects[[1]]
  row_child_bridge_changed_allocation <- row_child_bridge_changed_term$sd_binding$allocations[[1L]]
  row_child_bridge_changed_allocation$parent_factors[[1L]]$index <- 2L
  row_child_bridge_changed_term$sd_binding$allocations[[1L]] <- row_child_bridge_changed_allocation
  row_child_bridge_changed$mu$random_effects[[1]] <- row_child_bridge_changed_term
  expect_error(
    BayesTools:::.bt_JAGS_bridge_validate_formula_random_designs(
      row_child_bridge_design,
      row_child_bridge_changed
    ),
    "binding$factors",
    fixed = TRUE
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
      weights = prior("dirichlet", list(alpha = c(1, 1)))
    ),
    random_variance_allocation(
      name = "nested_split",
      parent = allocation_ref("total_re", "nested"),
      terms = c(study = "study", paper = "paper"),
      weights = prior("dirichlet", list(alpha = c(3, 2)))
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
      "mu__xRE_ALLOCx_total_re__allocation_sd",
      "mu__xRE_ALLOCx_total_re__weight",
      "mu__xRE_ALLOCx_nested_split__weight"
    )
  )
  expect_match(
    nested_result$formula_syntax,
    "mu__xRE_ALLOCx_total_re__component_nested_sd = mu__xRE_ALLOCx_total_re__allocation_sd * sqrt(mu__xRE_ALLOCx_total_re__weight[1])",
    fixed = TRUE
  )
  expect_match(
    nested_result$formula_syntax,
    "mu__xREx__study_intercept = mu__xRE_ALLOCx_total_re__allocation_sd * sqrt(mu__xRE_ALLOCx_total_re__weight[1]) * sqrt(mu__xRE_ALLOCx_nested_split__weight[1])",
    fixed = TRUE
  )
  expect_match(
    nested_result$formula_syntax,
    "mu__xREx__drug_intercept = mu__xRE_ALLOCx_total_re__allocation_sd * sqrt(mu__xRE_ALLOCx_total_re__weight[2])",
    fixed = TRUE
  )
  expect_equal(length(nested_result$formula_design$random_effects[[1]]$sd_binding$allocations[[1L]]$factors), 2L)
  expect_equal(length(nested_result$formula_design$random_effects[[3]]$sd_binding$allocations[[1L]]$factors), 1L)

  nested_gate_prior <- prior_random(
    random_variance_allocation(
      name = "total_re",
      terms = c(nested = "nested", drug = "drug"),
      sd = sd_prior,
      weights = prior("dirichlet", list(alpha = c(1, 1))),
      inclusion = list(nested = prior("spike", list(location = 0.5)))
    ),
    random_variance_allocation(
      name = "nested_split",
      parent = allocation_ref("total_re", "nested"),
      terms = c(study = "study", paper = "paper"),
      weights = prior("dirichlet", list(alpha = c(3, 2)))
    )
  )
  nested_gate_result <- JAGS_formula(
    formula = ~ 1 +
      random(1 | study, name = "study", covariance = "diag") +
      random(1 | paper, name = "paper", covariance = "diag") +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = df_nested,
    prior_list = fixed_priors,
    prior_random = nested_gate_prior
  )
  expect_match(
    nested_gate_result$formula_syntax,
    "mu__xRE_ALLOCx_total_re__component_nested_sd = mu__xRE_ALLOCx_total_re__allocation_sd * mu__xRE_ALLOCx_total_re__include_nested_indicator * sqrt(mu__xRE_ALLOCx_total_re__weight[1])",
    fixed = TRUE
  )
  expect_match(
    nested_gate_result$formula_syntax,
    "mu__xREx__study_intercept = mu__xRE_ALLOCx_total_re__allocation_sd * mu__xRE_ALLOCx_total_re__include_nested_indicator * sqrt(mu__xRE_ALLOCx_total_re__weight[1]) * sqrt(mu__xRE_ALLOCx_nested_split__weight[1])",
    fixed = TRUE
  )
  expect_true("mu__xRE_ALLOCx_total_re__include_nested_indicator" %in%
                nested_gate_result$add_parameters)
  nested_gate_posterior <- matrix(
    c(
      4, 1, 3, 3, 2, 1,
      4, 1, 3, 3, 2, 0
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      NULL,
      c(
        "mu__xRE_ALLOCx_total_re__allocation_sd",
        "prior_par_eta_mu__xRE_ALLOCx_total_re__weight[1]",
        "prior_par_eta_mu__xRE_ALLOCx_total_re__weight[2]",
        "prior_par_eta_mu__xRE_ALLOCx_nested_split__weight[1]",
        "prior_par_eta_mu__xRE_ALLOCx_nested_split__weight[2]",
        "mu__xRE_ALLOCx_total_re__include_nested_indicator"
      )
    )
  )
  expect_equal(
    BayesTools:::.bt_random_effect_sd_draws(
      random_term = nested_gate_result$formula_design$random_effects[[1]],
      n_columns = 1,
      posterior = nested_gate_posterior,
      prior_list = nested_gate_result$prior_list
    )[, 1],
    c(4 * sqrt(1 / 4) * sqrt(3 / 5), 0),
    tolerance = 1e-12
  )
  expect_equal(
    BayesTools:::.bt_random_effect_sd_draws(
      random_term = nested_gate_result$formula_design$random_effects[[3]],
      n_columns = 1,
      posterior = nested_gate_posterior,
      prior_list = nested_gate_result$prior_list
    )[, 1],
    c(4 * sqrt(3 / 4), 4 * sqrt(3 / 4)),
    tolerance = 1e-12
  )

  nested_external_prior <- prior_random(
    random_variance_allocation(
      name = "total_re",
      terms = c(nested = "nested", drug = "drug"),
      sd_source = random_sd_source("tau"),
      weights = prior("dirichlet", list(alpha = c(1, 1)))
    ),
    random_variance_allocation(
      name = "nested_split",
      parent = allocation_ref("total_re", "nested"),
      terms = c(study = "study", paper = "paper"),
      weights = prior("dirichlet", list(alpha = c(3, 2)))
    )
  )
  nested_external_result <- JAGS_formula(
    formula = ~ 1 +
      random(1 | study, name = "study", covariance = "diag") +
      random(1 | paper, name = "paper", covariance = "diag") +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = df_nested,
    prior_list = fixed_priors,
    prior_random = nested_external_prior
  )
  expect_equal(
    names(nested_external_result$prior_list),
    c(
      "mu_intercept",
      "mu__xRE_ALLOCx_total_re__weight",
      "mu__xRE_ALLOCx_nested_split__weight"
    )
  )
  expect_match(
    nested_external_result$formula_syntax,
    "mu__xRE_ALLOCx_total_re__component_nested_sd = tau * sqrt(mu__xRE_ALLOCx_total_re__weight[1])",
    fixed = TRUE
  )
  expect_match(
    nested_external_result$formula_syntax,
    "mu__xREx__study_intercept = tau * sqrt(mu__xRE_ALLOCx_total_re__weight[1]) * sqrt(mu__xRE_ALLOCx_nested_split__weight[1])",
    fixed = TRUE
  )
  expect_equal(
    nested_external_result$formula_design$random_effects[[1]]$sd_binding$allocations[[1L]]$source$kind,
    "external"
  )
  expect_equal(length(nested_external_result$formula_design$random_effects[[1]]$sd_binding$allocations[[1L]]$factors), 2L)

  nested_posterior <- matrix(
    c(4, 1, 3, 3, 2),
    nrow = 1,
    dimnames = list(
      NULL,
      c(
        "mu__xRE_ALLOCx_total_re__allocation_sd",
        "prior_par_eta_mu__xRE_ALLOCx_total_re__weight[1]",
        "prior_par_eta_mu__xRE_ALLOCx_total_re__weight[2]",
        "prior_par_eta_mu__xRE_ALLOCx_nested_split__weight[1]",
        "prior_par_eta_mu__xRE_ALLOCx_nested_split__weight[2]"
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
      weights = prior("dirichlet", list(alpha = c(1, 1)))
    ),
    random_variance_allocation(
      name = "het_sd",
      parent = allocation_ref("total_re", "het"),
      terms = "het",
      target = "sd_component",
      scale = "mean_variance",
      weights = prior("dirichlet", list(alpha = c(2, 2)))
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
    "mu__xREx__het_intercept = mu__xRE_ALLOCx_total_re__component_het_sd * sqrt(2 * mu__xRE_ALLOCx_het_sd__weight[1])",
    fixed = TRUE
  )
  expect_match(
    sd_leaf_result$formula_syntax,
    "mu__xREx__het_x = mu__xRE_ALLOCx_total_re__component_het_sd * sqrt(2 * mu__xRE_ALLOCx_het_sd__weight[2])",
    fixed = TRUE
  )
  expect_equal(
    sd_leaf_result$formula_design$random_effects[[1]]$sd_binding$allocations[[1L]]$target,
    "sd_component"
  )
  expect_equal(
    sd_leaf_result$formula_design$random_effects[[1]]$sd_binding$allocations[[1L]]$leaf_terms,
    c(mu__xREx__het_intercept = "intercept", mu__xREx__het_x = "x")
  )
  expect_equal(
    sd_leaf_result$formula_design$random_effects[[1]]$sd_leaves$leaf_names,
    names(sd_leaf_result$formula_design$random_effects[[1]]$sd_binding$allocations[[1L]]$leaf_terms)
  )
  expect_equal(
    sd_leaf_result$formula_design$random_effects[[1]]$sd_leaves$leaf_index_by_column,
    c(1L, 2L)
  )

  sd_leaf_samples <- c(
    "mu__xRE_ALLOCx_total_re__allocation_sd" = 2,
    "prior_par_eta_mu__xRE_ALLOCx_total_re__weight[1]" = 1,
    "prior_par_eta_mu__xRE_ALLOCx_total_re__weight[2]" = 3,
    "prior_par_eta_mu__xRE_ALLOCx_het_sd__weight[1]" = 1,
    "prior_par_eta_mu__xRE_ALLOCx_het_sd__weight[2]" = 3
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
      weights = prior("dirichlet", list(alpha = c(1, 1)))
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
    "target = \"sd_component\"",
    fixed = TRUE
  )
  expect_error(
    random_variance_allocation(
      name = "bad_sd_terms",
      terms = c("study", "drug"),
      sd = sd_prior,
      target = "sd_component"
    ),
    "exactly one",
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
          weights = prior("dirichlet", list(alpha = c(1, 1))),
          inclusion = list(missing = prior("spike", list(location = 0.5)))
        )
      )
    ),
    "unknown component label",
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
          target = "sd_component",
          weights = prior("dirichlet", list(alpha = c(1, 1))),
          inclusion = list(study = prior("spike", list(location = 0.5)))
        )
      )
    ),
    "supports target = \"block\"",
    fixed = TRUE
  )
  expect_error(
    random_variance_allocation(name = "allocation",
      terms = c(study = "study", "drug"),
      sd = sd_prior,
      weights = prior("dirichlet", list(alpha = c(1, 1)))
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
        random_variance_allocation(name = "allocation",
          sd = sd_prior,
          weights = prior("dirichlet", list(alpha = c(1, 1)))
        ),
        random_variance_allocation(
          name = "second",
          terms = c("study", "drug"),
          sd = sd_prior,
          weights = prior("dirichlet", list(alpha = c(1, 1)))
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
          weights = prior("dirichlet", list(alpha = c(1, 1)))
        ),
        random_variance_allocation(
          name = "duplicate",
          terms = c("study", "drug"),
          sd = sd_prior,
          weights = prior("dirichlet", list(alpha = c(1, 1)))
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
          weights = prior("dirichlet", list(alpha = c(1, 1)))
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
          target = "sd_component",
          scale = "mean_variance",
          weights = prior("dirichlet", list(alpha = c(1, 1)))
        ),
        random_variance_allocation(
          name = "total_re",
          terms = c("study", "drug"),
          sd = sd_prior,
          weights = prior("dirichlet", list(alpha = c(1, 1)))
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
          weights = prior("dirichlet", list(alpha = c(1, 1)))
        ),
        random_variance_allocation(
          name = "bad_ref",
          parent = allocation_ref("total_re", "missing"),
          terms = "study",
          target = "sd_component",
          scale = "mean_variance",
          weights = prior("dirichlet", list(alpha = c(1, 1)))
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
          weights = prior("dirichlet", list(alpha = c(1, 1)))
        ),
        random_variance_allocation(
          name = "study_split_1",
          parent = allocation_ref("total_re", "study"),
          terms = "study",
          target = "sd_component",
          weights = prior("dirichlet", list(alpha = c(1, 1)))
        ),
        random_variance_allocation(
          name = "study_split_2",
          parent = allocation_ref("total_re", "study"),
          terms = "drug",
          target = "sd_component",
          weights = prior("dirichlet", list(alpha = c(1, 1)))
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
          weights = prior("dirichlet", list(alpha = c(1, 1)))
        ),
        random_variance_allocation(
          name = "study_split",
          parent = allocation_ref("total_re", "study"),
          terms = c("study", "drug"),
          weights = prior("dirichlet", list(alpha = c(1, 1)))
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
          target = "sd_component",
          weights = prior("dirichlet", list(alpha = c(1, 1, 1)))
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
          target = "sd_component",
          weights = prior("dirichlet", list(alpha = c(1, 1)))
        )
      )
    ),
    "at least two resolved SD components",
    fixed = TRUE
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
    "supplied formula-related inputs do not fully match",
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
    "supplied formula inputs cannot replace fitted replay metadata",
    fixed = TRUE
  )
  expect_error(
    JAGS_bridgesampling(
      fit = fit,
      log_posterior = function(parameters, data) 0,
      data = list(),
      prior_list = NULL
    ),
    "requires posterior samples of standardized latent random effects",
    fixed = TRUE
  )
})

test_that("JAGS bridgesampling can reconstruct formula parameters from fitted design metadata", {

  df <- data.frame(x = c(10, 20, 30))
  prior_list <- list(
    intercept = prior("normal", list(0, 1)),
    x = prior("normal", list(0, 1))
  )
  formula_result <- JAGS_formula(
    formula = ~ 1 + x,
    parameter = "mu",
    data = df,
    prior_list = prior_list,
    formula_scale = list(x = TRUE)
  )
  samples <- c(mu_intercept = 10, mu_x = 2)

  reconstructed <- JAGS_marglik_parameters_formula(
    samples = samples,
    formula_list = NULL,
    formula_data_list = NULL,
    formula_prior_list = list(mu = formula_result$prior_list),
    prior_list_parameters = list(),
    formula_design_list = list(mu = formula_result$formula_design)
  )
  expected_x <- as.vector(formula_result$formula_design$model_matrix[, "x"])
  expect_equal(reconstructed$mu, 10 + 2 * expected_x, tolerance = 1e-12)

  log_formula <- ~ 1 + x
  attr(log_formula, "log(intercept)") <- TRUE
  log_result <- JAGS_formula(
    formula = log_formula,
    parameter = "mu",
    data = df,
    prior_list = list(
      intercept = prior("gamma", list(2, 2)),
      x = prior("normal", list(0, 1))
    )
  )
  log_reconstructed <- JAGS_marglik_parameters_formula(
    samples = c(mu_intercept = exp(1), mu_x = 2),
    formula_list = NULL,
    formula_data_list = NULL,
    formula_prior_list = list(mu = log_result$prior_list),
    prior_list_parameters = list(),
    formula_design_list = list(mu = log_result$formula_design)
  )
  expect_equal(log_reconstructed$mu, 1 + 2 * df$x, tolerance = 1e-12)
})

test_that("JAGS bridgesampling uses fitted formula metadata and errors on supplied mismatches", {

  df <- data.frame(x = c(10, 20, 30))
  prior_list <- list(
    intercept = prior("normal", list(0, 1)),
    x = prior("normal", list(0, 1))
  )
  formula_result <- JAGS_formula(
    formula = ~ 1 + x,
    parameter = "mu",
    data = df,
    prior_list = prior_list,
    formula_scale = list(x = TRUE)
  )
  fit <- coda::mcmc(matrix(
    0,
    nrow = 1,
    ncol = 1,
    dimnames = list(NULL, "mu_intercept")
  ))
  attr(fit, "prior_list") <- formula_result$prior_list
  attr(fit, "formula_design") <- list(mu = formula_result$formula_design)

  expect_error(
    JAGS_bridgesampling(
      fit = fit,
      log_posterior = function(parameters, data) 0,
      data = list()
    ),
    "posterior' does not contain",
    fixed = TRUE
  )

  expect_error(
    JAGS_bridgesampling(
      fit = fit,
      log_posterior = function(parameters, data) 0,
      data = list(),
      formula_list = list(mu = ~ 1 + x),
      formula_data_list = list(mu = df),
      formula_prior_list = list(mu = prior_list),
      formula_scale_list = list(wrong_name = list(x = TRUE))
    ),
    "supplied formula-related inputs do not fully match",
    fixed = TRUE
  )

  expect_error(
    JAGS_bridgesampling(
      fit = fit,
      log_posterior = function(parameters, data) 0,
      data = list(),
      formula_list = list(mu = ~ 1 + x),
      formula_data_list = list(mu = df),
      formula_prior_list = list(mu = prior_list),
      formula_scale_list = TRUE
    ),
    "supplied formula-related inputs do not fully match",
    fixed = TRUE
  )

  expect_warning(
    expect_error(
      JAGS_bridgesampling(
        fit = fit,
        log_posterior = function(parameters, data) 0,
        data = list(),
        prior_list = formula_result$prior_list
      ),
      "posterior' does not contain",
      fixed = TRUE
    ),
    "received formula priors in 'prior_list'",
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
      add_parameters = NULL,
      add_bounds = list(lb = -Inf, ub = Inf),
      bridge_parameters = bridge$parameters,
      bridge_bounds = bridge$bounds
    ),
    "requires at least one 'add_parameters'",
    fixed = TRUE
  )
  expect_error(
    BayesTools:::.bt_JAGS_bridge_merge_add_parameters(
      add_parameters = expected_z[1],
      add_bounds = list(lb = -Inf, ub = Inf),
      bridge_parameters = bridge$parameters,
      bridge_bounds = bridge$bounds
    ),
    "names must be unique and match",
    fixed = TRUE
  )
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
        lb = stats::setNames(-Inf, expected_z[1]),
        ub = stats::setNames(1 - 1e-12, expected_z[1])
      ),
      bridge_parameters = expected_z[1],
      bridge_bounds = list(
        lb = stats::setNames(-Inf, expected_z[1]),
        ub = stats::setNames(1, expected_z[1])
      )
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
    bayestools_reference_formula_random_log_prior(
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
    bayestools_reference_formula_random_log_prior(
      lkj_samples,
      list(mu = lkj_result$formula_design)
    ),
    sum(stats::dnorm(lkj_samples[grep("xRE_Zx", names(lkj_samples), value = TRUE)], log = TRUE)) +
      BayesTools:::.bt_lkj_cholesky_cpc_u_log_prior(0.75, K = 2, eta = 1.5),
    tolerance = 1e-12
  )
  expect_error(
    bayestools_reference_formula_random_log_prior(
      lkj_samples,
      list(mu = malformed_lkj)
    ),
    "missing canonical 'random_term\\$structure'"
  )

  missing_correlation_lkj <- lkj_result$formula_design
  missing_correlation_lkj$random_effects[[1]]$correlation <- NULL
  expect_error(
    bayestools_reference_formula_random_log_prior(
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
    .parameter_catalog_random_summary_samples(
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
    bayestools_reference_formula_random_log_prior(
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
        cor = prior("normal", list(0, 0.5))
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

  ar_saturated_samples <- ar_samples
  ar_saturated_samples[["mu__xREx__id_rho_z"]] <- 20
  expect_equal(
    bayestools_reference_random_effect_scalar_rho_support(
      ar_saturated_samples,
      ar_result$formula_design$random_effects[[1]]
    ),
    -Inf
  )
  expect_error(
    BayesTools:::.bt_JAGS_marglik_random_effect_rho(
      ar_saturated_samples,
      ar_result$formula_design$random_effects[[1]]
    ),
    "missing or invalid scalar correlation coordinates",
    fixed = TRUE
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
        cor = prior("normal", list(0, 0.5))
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
  car_logit_saturated_samples <- c(
    car_samples[names(car_samples) != "mu__xREx__id_rho_z"],
    mu__xREx__id_rho_logit = 1000
  )
  expect_equal(
    bayestools_reference_random_effect_scalar_rho_support(
      car_logit_saturated_samples,
      car_logit$formula_design$random_effects[[1]]
    ),
    -Inf
  )
  expect_error(
    BayesTools:::.bt_JAGS_marglik_random_effect_rho(
      car_logit_saturated_samples,
      car_logit$formula_design$random_effects[[1]]
    ),
    "missing or invalid scalar correlation coordinates",
    fixed = TRUE
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
        cor = prior("normal", list(0, 0.5))
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
    bayestools_reference_formula_random_log_prior(
      hcs_boundary_samples,
      list(mu = hcs_result$formula_design)
    ),
    -Inf
  )
  hcs_interior_samples <- hcs_boundary_samples
  hcs_interior_samples[["mu__xREx__id_rho_z"]] <- 0
  expect_true(is.finite(
    bayestools_reference_formula_random_log_prior(
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
        cor = prior("point", list(location = 0.2))
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
    bayestools_reference_formula_random_log_prior(
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

test_that("JAGS_evaluate_formula rejects new levels for known group covariance", {

  df <- data.frame(
    id = factor(c("a", "b", "a", "c"), levels = c("a", "b", "c"))
  )
  K <- matrix(
    c(2, .4, .2,
      .4, 3, .5,
      .2, .5, 4),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("a", "b", "c"), c("a", "b", "c"))
  )
  random_effects <- random_effects_formula(
    ~ 1 | id,
    group_covariance = random_group_covariance(K, scale = "none")
  )
  formula_result <- JAGS_formula(
    formula = random_effects,
    parameter = "mu",
    data = df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(
        sd = prior("gamma", list(2, 1)),
        monitor = random_monitor(latent = TRUE, coefficients = FALSE)
      )
    )
  )
  random_term <- formula_result$formula_design$random_effects[[1]]
  z_names <- as.vector(BayesTools:::.bt_random_effect_latent_names(
    random_term = random_term,
    n_groups = random_term$n_groups,
    n_columns = random_term$n_columns
  ))
  posterior <- matrix(
    c(10, 2, .1, -.2, .3),
    nrow = 1,
    dimnames = list(NULL, c("mu_intercept", random_term$sd_parameter_names, z_names))
  )
  fit <- coda::mcmc(posterior)
  attr(fit, "formula_design") <- list(mu = formula_result$formula_design)

  prediction <- JAGS_evaluate_formula(
    fit = fit,
    formula = random_effects$formula,
    parameter = "mu",
    data = df,
    prior_list = formula_result$prior_list
  )
  expect_equal(
    unname(drop(prediction)),
    10 + 2 * c(.1, -.2, .1, .3),
    tolerance = 1e-12
  )
  for(new_levels in c("zero", "sample")){
    conditional_prediction <- JAGS_evaluate_formula(
      fit = fit,
      formula = random_effects$formula,
      parameter = "mu",
      data = df,
      prior_list = formula_result$prior_list,
      formula_target = "conditional",
      new_levels = new_levels
    )
    expect_equal(conditional_prediction, prediction, tolerance = 1e-12)
  }

  new_data <- data.frame(
    id = factor(c("a", "d"), levels = c("a", "b", "c", "d"))
  )
  expect_error(
    JAGS_evaluate_formula(
      fit = fit,
      formula = random_effects$formula,
      parameter = "mu",
      data = new_data,
      prior_list = formula_result$prior_list
    ),
    "known group covariance",
    fixed = TRUE
  )
  expect_error(
    JAGS_evaluate_formula(
      fit = fit,
      formula = random_effects$formula,
      parameter = "mu",
      data = new_data,
      prior_list = formula_result$prior_list,
      formula_target = "conditional",
      new_levels = "sample"
    ),
    "known group covariance",
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
      allocation = random_variance_allocation(name = "allocation",
        sd = prior("gamma", list(2, 2)),
        weights = prior("dirichlet", list(alpha = c(2, 3)))
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
      "mu__xRE_ALLOCx_allocation__allocation_sd",
      "prior_par_eta_mu__xRE_ALLOCx_allocation__weight[1]",
      "prior_par_eta_mu__xRE_ALLOCx_allocation__weight[2]",
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
      "mu__xRE_ALLOCx_allocation__allocation_sd",
      "mu__xRE_ALLOCx_allocation__weight[1]",
      "mu__xRE_ALLOCx_allocation__weight[2]",
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

test_that("JAGS_evaluate_formula reconstructs row-indexed external SD random effects", {

  df <- data.frame(
    study = factor(c("s1", "s1", "s2", "s2")),
    drug = factor(c("a", "b", "a", "b"))
  )
  fixed_priors <- list(intercept = prior("normal", list(0, 1)))
  formula_result <- JAGS_formula(
    formula = ~ 1 +
      random(1 | study, name = "study", covariance = "diag") +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = df,
    prior_list = fixed_priors,
    prior_random = prior_random(
      allocation = random_variance_allocation(name = "allocation",
        sd_source = random_sd_source("tau", shape = "row"),
        weights = prior("dirichlet", list(alpha = c(2, 3)))
      )
    )
  )

  posterior <- matrix(
    c(
      10,
      2, 4, 6, 8,
      1, 3,
      1, 2,
      3, 4
    ),
    nrow = 1,
    dimnames = list(NULL, c(
      "mu_intercept",
      "tau[1]",
      "tau[2]",
      "tau[3]",
      "tau[4]",
      "prior_par_eta_mu__xRE_ALLOCx_allocation__weight[1]",
      "prior_par_eta_mu__xRE_ALLOCx_allocation__weight[2]",
      "mu__xREx__study_xRE_Zx[1,1]",
      "mu__xREx__study_xRE_Zx[2,1]",
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
    prior_list = formula_result$prior_list,
    fitted_rows = seq_len(nrow(df))
  )
  expect_equal(
    unname(drop(prediction)),
    10 + c(
      2 * (sqrt(1 / 4) * 1 + sqrt(3 / 4) * 3),
      4 * (sqrt(1 / 4) * 1 + sqrt(3 / 4) * 4),
      6 * (sqrt(1 / 4) * 2 + sqrt(3 / 4) * 3),
      8 * (sqrt(1 / 4) * 2 + sqrt(3 / 4) * 4)
    ),
    tolerance = 1e-12
  )

  multi_posterior <- matrix(
    c(
      10,
      2, 4, 6, 8,
      1, 3,
      1, 2,
      3, 4,

      -5,
      1, 10, 100, 1000,
      4, 1,
      -2, 5,
      7, -11
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(NULL, colnames(posterior))
  )
  multi_fit <- coda::mcmc(multi_posterior)
  attr(multi_fit, "formula_design") <- list(mu = formula_result$formula_design)
  multi_prediction <- JAGS_evaluate_formula(
    fit = multi_fit,
    formula = ~ 1 +
      random(1 | study, name = "study", covariance = "diag") +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = df,
    prior_list = formula_result$prior_list,
    fitted_rows = seq_len(nrow(df))
  )
  expect_equal(dim(multi_prediction), c(4L, 2L))
  expect_equal(
    unname(multi_prediction),
    unname(cbind(
      10 + c(
        2 * (sqrt(1 / 4) * 1 + sqrt(3 / 4) * 3),
        4 * (sqrt(1 / 4) * 1 + sqrt(3 / 4) * 4),
        6 * (sqrt(1 / 4) * 2 + sqrt(3 / 4) * 3),
        8 * (sqrt(1 / 4) * 2 + sqrt(3 / 4) * 4)
      ),
      -5 + c(
        1    * (sqrt(4 / 5) * -2 + sqrt(1 / 5) *   7),
        10   * (sqrt(4 / 5) * -2 + sqrt(1 / 5) * -11),
        100  * (sqrt(4 / 5) *  5 + sqrt(1 / 5) *   7),
        1000 * (sqrt(4 / 5) *  5 + sqrt(1 / 5) * -11)
      )
    )),
    tolerance = 1e-12
  )

  newdata <- data.frame(
    study = factor(c("s2", "s1"), levels = c("s1", "s2")),
    drug = factor(c("b", "a"), levels = c("a", "b"))
  )
  new_posterior <- posterior[, c(
    "mu_intercept",
    "tau[1]",
    "tau[4]",
    "prior_par_eta_mu__xRE_ALLOCx_allocation__weight[1]",
    "prior_par_eta_mu__xRE_ALLOCx_allocation__weight[2]",
    "mu__xREx__study_xRE_Zx[1,1]",
    "mu__xREx__study_xRE_Zx[2,1]",
    "mu__xREx__drug_xRE_Zx[1,1]",
    "mu__xREx__drug_xRE_Zx[2,1]"
  ), drop = FALSE]
  new_posterior[, "tau[1]"] <- 20
  new_posterior[, "tau[4]"] <- 10
  new_fit <- coda::mcmc(new_posterior)
  attr(new_fit, "formula_design") <- list(mu = formula_result$formula_design)
  new_prediction <- JAGS_evaluate_formula(
    fit = new_fit,
    formula = ~ 1 +
      random(1 | study, name = "study", covariance = "diag") +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = newdata,
    prior_list = formula_result$prior_list,
    fitted_rows = c(4L, 1L)
  )
  expect_equal(
    unname(drop(new_prediction)),
    10 + c(
      10 * (sqrt(1 / 4) * 2 + sqrt(3 / 4) * 4),
      20 * (sqrt(1 / 4) * 1 + sqrt(3 / 4) * 3)
    ),
    tolerance = 1e-12
  )

  sd_component_df <- data.frame(
    x = c(-1, 0, 1, 2),
    study = factor(c("s1", "s1", "s2", "s2"))
  )
  sd_component_result <- JAGS_formula(
    formula = ~ 1 +
      random(1 + x | study, name = "study", covariance = "diag"),
    parameter = "mu",
    data = sd_component_df,
    prior_list = fixed_priors,
    prior_random = prior_random(
      allocation = random_variance_allocation(name = "allocation",
        terms = "study",
        target = "sd_component",
        sd_source = random_sd_source("tau", shape = "row"),
        weights = prior("dirichlet", list(alpha = c(1, 1)))
      )
    )
  )
  sd_component_posterior <- matrix(
    c(
      5,
      2, 4, 6, 8,
      1, 3,
      1, 2,
      3, 4
    ),
    nrow = 1,
    dimnames = list(NULL, c(
      "mu_intercept",
      "tau[1]",
      "tau[2]",
      "tau[3]",
      "tau[4]",
      "prior_par_eta_mu__xRE_ALLOCx_allocation__weight[1]",
      "prior_par_eta_mu__xRE_ALLOCx_allocation__weight[2]",
      "mu__xREx__study_xRE_Zx[1,1]",
      "mu__xREx__study_xRE_Zx[2,1]",
      "mu__xREx__study_xRE_Zx[1,2]",
      "mu__xREx__study_xRE_Zx[2,2]"
    ))
  )
  sd_component_fit <- coda::mcmc(sd_component_posterior)
  attr(sd_component_fit, "formula_design") <- list(mu = sd_component_result$formula_design)
  sd_component_prediction <- JAGS_evaluate_formula(
    fit = sd_component_fit,
    formula = ~ 1 +
      random(1 + x | study, name = "study", covariance = "diag"),
    parameter = "mu",
    data = sd_component_df,
    prior_list = sd_component_result$prior_list,
    fitted_rows = seq_len(nrow(sd_component_df))
  )
  expect_equal(
    unname(drop(sd_component_prediction)),
    5 + c(
      2 * (sqrt(1 / 4) * 1 + sqrt(3 / 4) * 3 * -1),
      4 * (sqrt(1 / 4) * 1 + sqrt(3 / 4) * 3 * 0),
      6 * (sqrt(1 / 4) * 2 + sqrt(3 / 4) * 4 * 1),
      8 * (sqrt(1 / 4) * 2 + sqrt(3 / 4) * 4 * 2)
    ),
    tolerance = 1e-12
  )

  missing_source_fit <- coda::mcmc(
    posterior[, colnames(posterior) != "tau[2]", drop = FALSE]
  )
  attr(missing_source_fit, "formula_design") <- list(mu = formula_result$formula_design)
  expect_error(
    JAGS_evaluate_formula(
      fit = missing_source_fit,
      formula = ~ 1 +
        random(1 | study, name = "study", covariance = "diag") +
        random(1 | drug, name = "drug", covariance = "diag"),
      parameter = "mu",
      data = df[1:2, , drop = FALSE],
      prior_list = formula_result$prior_list,
      fitted_rows = 1:2
    ),
    "missing values for prediction row(s): 2",
    fixed = TRUE
  )

  scalar_source_posterior <- posterior[, c(
    "mu_intercept",
    "tau[1]",
    "prior_par_eta_mu__xRE_ALLOCx_allocation__weight[1]",
    "prior_par_eta_mu__xRE_ALLOCx_allocation__weight[2]",
    "mu__xREx__study_xRE_Zx[1,1]",
    "mu__xREx__study_xRE_Zx[2,1]",
    "mu__xREx__drug_xRE_Zx[1,1]",
    "mu__xREx__drug_xRE_Zx[2,1]"
  ), drop = FALSE]
  colnames(scalar_source_posterior)[colnames(scalar_source_posterior) == "tau[1]"] <- "tau"
  scalar_source_fit <- coda::mcmc(scalar_source_posterior)
  attr(scalar_source_fit, "formula_design") <- list(mu = formula_result$formula_design)
  expect_error(
    JAGS_evaluate_formula(
      fit = scalar_source_fit,
      formula = ~ 1 +
        random(1 | study, name = "study", covariance = "diag") +
        random(1 | drug, name = "drug", covariance = "diag"),
      parameter = "mu",
      data = df[1, , drop = FALSE],
      prior_list = formula_result$prior_list,
      fitted_rows = 1L
    ),
    "missing values for prediction row(s): 1",
    fixed = TRUE
  )

  coefficient_names <- as.vector(BayesTools:::.bt_random_effect_coefficient_names(
    random_term = formula_result$formula_design$random_effects[[1]],
    n_groups = formula_result$formula_design$random_effects[[1]]$n_groups,
    n_columns = formula_result$formula_design$random_effects[[1]]$n_columns
  ))
  coefficient_only_posterior <- cbind(
    posterior[, c(
      "mu_intercept",
      "tau[1]",
      "tau[2]",
      "tau[3]",
      "tau[4]",
      "prior_par_eta_mu__xRE_ALLOCx_allocation__weight[1]",
      "prior_par_eta_mu__xRE_ALLOCx_allocation__weight[2]"
    ), drop = FALSE],
    matrix(
      seq_along(coefficient_names),
      nrow = 1,
      dimnames = list(NULL, coefficient_names)
    )
  )
  coefficient_only_fit <- coda::mcmc(coefficient_only_posterior)
  attr(coefficient_only_fit, "formula_design") <- list(mu = formula_result$formula_design)
  expect_error(
    JAGS_evaluate_formula(
      fit = coefficient_only_fit,
      formula = ~ 1 +
        random(1 | study, name = "study", covariance = "diag") +
        random(1 | drug, name = "drug", covariance = "diag"),
      parameter = "mu",
      data = df,
      prior_list = formula_result$prior_list,
      fitted_rows = seq_len(nrow(df))
    ),
    "cannot be reconstructed from the posterior samples",
    fixed = TRUE
  )

  expect_error(
    JAGS_evaluate_formula(
      fit = fit,
      formula = ~ 1 +
        random(1 | study, name = "study", covariance = "diag") +
        random(1 | drug, name = "drug", covariance = "diag"),
      parameter = "mu",
      data = data.frame(
        study = factor("s3", levels = "s3"),
        drug = factor("a", levels = c("a", "b"))
      ),
      prior_list = formula_result$prior_list
    ),
    "New random-effect level",
    fixed = TRUE
  )

  nested_df <- data.frame(
    study = factor(c("s1", "s1", "s2", "s2")),
    estimate = factor(c("e1", "e2", "e3", "e4"))
  )
  nested_result <- JAGS_formula(
    formula = ~ 1 +
      random(1 | estimate, name = "estimate", covariance = "diag") +
      random(1 | study, name = "study", covariance = "diag"),
    parameter = "mu",
    data = nested_df,
    prior_list = fixed_priors,
    prior_random = prior_random(
      allocation = random_variance_allocation(name = "allocation",
        sd_source = random_sd_source("tau", shape = "row"),
        weights = prior("dirichlet", list(alpha = c(1, 1)))
      )
    )
  )
  nested_posterior <- matrix(
    c(
      10,
      0.25, 0.75,
      2, 3, 4, 5,
      0.1, 0.2, -0.1, 0.4,
      1, -0.5
    ),
    nrow = 1,
    dimnames = list(NULL, c(
      "mu_intercept",
      "mu__xRE_ALLOCx_allocation__weight[1]",
      "mu__xRE_ALLOCx_allocation__weight[2]",
      "tau[1]",
      "tau[2]",
      "tau[3]",
      "tau[4]",
      "mu__xREx__estimate_xRE_Zx[1,1]",
      "mu__xREx__estimate_xRE_Zx[2,1]",
      "mu__xREx__estimate_xRE_Zx[3,1]",
      "mu__xREx__estimate_xRE_Zx[4,1]",
      "mu__xREx__study_xRE_Zx[1,1]",
      "mu__xREx__study_xRE_Zx[2,1]"
    ))
  )
  nested_fit <- coda::mcmc(nested_posterior)
  attr(nested_fit, "formula_design") <- list(mu = nested_result$formula_design)
  nested_prediction <- JAGS_evaluate_formula(
    fit = nested_fit,
    formula = ~ 1 +
      random(1 | estimate, name = "estimate", covariance = "diag") +
      random(1 | study, name = "study", covariance = "diag"),
    parameter = "mu",
    data = nested_df,
    prior_list = nested_result$prior_list,
    fitted_rows = seq_len(nrow(nested_df))
  )
  expect_equal(
    unname(drop(nested_prediction)),
    10 + c(
      2 * (sqrt(0.25) *  0.1 + sqrt(0.75) *  1),
      3 * (sqrt(0.25) *  0.2 + sqrt(0.75) *  1),
      4 * (sqrt(0.25) * -0.1 + sqrt(0.75) * -0.5),
      5 * (sqrt(0.25) *  0.4 + sqrt(0.75) * -0.5)
    ),
    tolerance = 1e-12
  )

  slope_df <- data.frame(
    x = c(1, 2, 3, 4),
    study = factor(c("s1", "s1", "s2", "s2")),
    drug = factor(c("a", "b", "a", "b"))
  )
  slope_result <- JAGS_formula(
    formula = ~ 1 + x +
      id(1 + x | study) +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = slope_df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    ),
    formula_scale = list(x = TRUE),
    prior_random = prior_random(
      allocation = random_variance_allocation(name = "allocation",
        sd_source = random_sd_source("tau", shape = "row"),
        weights = prior("dirichlet", list(alpha = c(1, 1)))
      )
    )
  )
  slope_newdata <- data.frame(
    x = c(2, 4),
    study = factor(c("s1", "s2"), levels = c("s1", "s2")),
    drug = factor(c("a", "b"), levels = c("a", "b"))
  )
  slope_posterior <- matrix(
    c(
      0, 0,
      2, 3,
      0.25, 0.75,
      1, 2,
      -1, 0.5,
      0, 0
    ),
    nrow = 1,
    dimnames = list(NULL, c(
      "mu_intercept",
      "mu_x",
      "tau[2]",
      "tau[4]",
      "mu__xRE_ALLOCx_allocation__weight[1]",
      "mu__xRE_ALLOCx_allocation__weight[2]",
      "mu__xREx__study_xRE_Zx[1,1]",
      "mu__xREx__study_xRE_Zx[1,2]",
      "mu__xREx__study_xRE_Zx[2,1]",
      "mu__xREx__study_xRE_Zx[2,2]",
      "mu__xREx__drug_xRE_Zx[1,1]",
      "mu__xREx__drug_xRE_Zx[2,1]"
    ))
  )
  slope_fit <- coda::mcmc(slope_posterior)
  attr(slope_fit, "formula_design") <- list(mu = slope_result$formula_design)
  attr(slope_fit, "formula_scale") <- list(mu = slope_result$formula_scale)
  slope_prediction <- JAGS_evaluate_formula(
    fit = slope_fit,
    formula = ~ 1 + x +
      id(1 + x | study) +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = slope_newdata,
    prior_list = slope_result$prior_list,
    fitted_rows = c(2L, 4L)
  )
  scaled_x <- (slope_newdata$x - slope_result$formula_scale$mu_x$mean) /
    slope_result$formula_scale$mu_x$sd
  expect_equal(
    unname(drop(slope_prediction)),
    c(
      2 * sqrt(0.25) * (1 + scaled_x[1] * 2),
      3 * sqrt(0.25) * (-1 + scaled_x[2] * 0.5)
    ),
    tolerance = 1e-12
  )
  slope_value_source <- parameter_source(
    "tau",
    shape = "row",
    values = function(parameters, data, n_rows){
      parameters$base_tau * data$x[seq_len(n_rows)]
    }
  )
  slope_values_result <- JAGS_formula(
    formula = ~ 1 + x +
      id(1 + x | study) +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = slope_df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    ),
    formula_scale = list(x = TRUE),
    prior_random = prior_random(
      allocation = random_variance_allocation(name = "allocation",
        sd_source = random_sd_source(slope_value_source),
        weights = prior("dirichlet", list(alpha = c(1, 1)))
      )
    )
  )
  slope_values_posterior <- cbind(
    slope_posterior[, !colnames(slope_posterior) %in% c("tau[1]", "tau[2]"), drop = FALSE],
    "base_tau" = 1
  )
  slope_values_fit <- coda::mcmc(slope_values_posterior)
  attr(slope_values_fit, "formula_design") <- list(mu = slope_values_result$formula_design)
  attr(slope_values_fit, "formula_scale") <- list(mu = slope_values_result$formula_scale)
  slope_values_prediction <- JAGS_evaluate_formula(
    fit = slope_values_fit,
    formula = ~ 1 + x +
      id(1 + x | study) +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = slope_newdata,
    prior_list = slope_values_result$prior_list
  )
  expect_equal(
    unname(drop(slope_values_prediction)),
    c(
      2 * sqrt(0.25) * (1 + scaled_x[1] * 2),
      4 * sqrt(0.25) * (-1 + scaled_x[2] * 0.5)
    ),
    tolerance = 1e-12
  )
  slope_values_stale_fit <- coda::mcmc(cbind(
    slope_values_posterior,
    "tau[1]" = 999,
    "tau[2]" = 999
  ))
  attr(slope_values_stale_fit, "formula_design") <- list(mu = slope_values_result$formula_design)
  attr(slope_values_stale_fit, "formula_scale") <- list(mu = slope_values_result$formula_scale)
  slope_values_stale_prediction <- JAGS_evaluate_formula(
    fit = slope_values_stale_fit,
    formula = ~ 1 + x +
      id(1 + x | study) +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = slope_newdata,
    prior_list = slope_values_result$prior_list
  )
  expect_equal(
    unname(drop(slope_values_stale_prediction)),
    unname(drop(slope_values_prediction)),
    tolerance = 1e-12
  )

  cs_df <- data.frame(
    f = factor(c("a", "b", "a", "b")),
    study = factor(c("s1", "s1", "s2", "s2")),
    drug = factor(c("a", "b", "a", "b"))
  )
  cs_result <- JAGS_formula(
    formula = ~ 1 +
      cs(f | study) +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = cs_df,
    prior_list = fixed_priors,
    prior_random = prior_random(
      allocation = random_variance_allocation(name = "allocation",
        sd_source = random_sd_source("tau", shape = "row"),
        weights = prior("dirichlet", list(alpha = c(1, 1)))
      ),
      study = random_block(cor = prior("normal", list(0, 0.5)))
    )
  )
  cs_posterior <- matrix(
    c(
      0,
      2, 3,
      0.25, 0.75,
      0,
      1, 2,
      -1, 0.5,
      0, 0
    ),
    nrow = 1,
    dimnames = list(NULL, c(
      "mu_intercept",
      "tau[1]",
      "tau[2]",
      "mu__xRE_ALLOCx_allocation__weight[1]",
      "mu__xRE_ALLOCx_allocation__weight[2]",
      "mu__xREx__study_rho_z",
      "mu__xREx__study_xRE_Zx[1,1]",
      "mu__xREx__study_xRE_Zx[1,2]",
      "mu__xREx__study_xRE_Zx[2,1]",
      "mu__xREx__study_xRE_Zx[2,2]",
      "mu__xREx__drug_xRE_Zx[1,1]",
      "mu__xREx__drug_xRE_Zx[2,1]"
    ))
  )
  cs_fit <- coda::mcmc(cs_posterior)
  attr(cs_fit, "formula_design") <- list(mu = cs_result$formula_design)
  cs_prediction <- JAGS_evaluate_formula(
    fit = cs_fit,
    formula = ~ 1 +
      cs(f | study) +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = cs_df[1:2, , drop = FALSE],
    prior_list = cs_result$prior_list,
    fitted_rows = 1:2
  )
  expect_equal(
    unname(drop(cs_prediction)),
    c(
      2 * sqrt(0.25) * 1,
      3 * sqrt(0.25) * 2
    ),
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
          cor = rho_prior,
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

test_that("transform_scale_samples leaves structured random-effect SDs on the fitted scale", {

  df <- data.frame(
    time = rep(c(1, 2, 3), 2),
    id = factor(rep(c("a", "b"), each = 3), levels = c("a", "b"))
  )
  sd_prior <- prior("normal", list(0, 1), truncation = list(lower = 0, upper = Inf))
  rho_prior <- prior("normal", list(0, 0.5))

  for (structure in c("hcs", "har")) {
    formula_result <- JAGS_formula(
      formula = stats::as.formula(paste0("~ 1 + time + ", structure, "(time | id)")),
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
          cor = rho_prior
        )
      )
    )
    random_term <- formula_result$formula_design$random_effects[[1]]
    sd_cols <- unique(random_term$sd_parameter_names)
    posterior <- matrix(
      seq_along(sd_cols),
      nrow = 1,
      dimnames = list(NULL, sd_cols)
    )
    transformed <- transform_scale_samples(
      posterior,
      list(mu = formula_result$formula_scale)
    )
    expect_equal(
      unname(transformed[1, sd_cols]),
      seq_along(sd_cols),
      info = structure
    )
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

test_that("random-effect grouping variables use raw data when predictors are scaled", {

  df <- data.frame(
    x = c(-1, 0.1, 0.2, 100),
    xpos = c(FALSE, TRUE, TRUE, TRUE),
    id = factor(c("a", "a", "b", "b"))
  )
  formula_result <- JAGS_formula(
    formula = ~ 1 + x + random(1 | xpos, name = "xpos", covariance = "diag"),
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

  raw_group_map <- as.numeric(factor(df$xpos, levels = levels(as.factor(df$xpos))))
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
    formula = ~ 1 + x + random(1 | xpos, name = "xpos", covariance = "diag"),
    parameter = "mu",
    data = data.frame(x = 0.1, xpos = TRUE),
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

.canonical_us_sd_leaves <- function(){

  leaf_names <- c(
    "mu__xREx__id_intercept",
    "mu__xREx__id_x"
  )
  out <- list(
    random_structure = "us",
    leaf_names_by_column = leaf_names,
    leaf_terms_by_column = c("intercept", "x"),
    leaf_names = leaf_names,
    leaf_terms = stats::setNames(c("intercept", "x"), leaf_names),
    column_names = c("(Intercept)", "x")
  )
  class(out) <- c("BayesTools_random_effect_sd_leaves", "list")
  list(`__xREx__id` = out)
}


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
  attr(formula_scale$mu, "random_effect_sd_leaves") <-
    .canonical_us_sd_leaves()

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
    .parameter_catalog_random_summary_samples(
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
    rep(c(0, 0, source_sd, as.vector(source_cor), as.vector(source_L)), 2),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      NULL,
      c(
        "mu_intercept",
        "mu_x",
        sd_names,
        as.vector(R_names),
        as.vector(L_names)
      )
    )
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
  fit <- attach_test_parameter_map(fit)

  semantic_samples <- JAGS_estimates_table(
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
    unname(semantic_samples[, "(mu) sd(intercept)"]),
    rep(expected_sd[1], 2),
    tolerance = 1e-12
  )
  expect_equal(
    unname(semantic_samples[, "(mu) sd(x)"]),
    rep(expected_sd[2], 2),
    tolerance = 1e-12
  )
  expect_equal(
    unname(semantic_samples[, "(mu) cor(intercept,x)"]),
    rep(expected_cor, 2),
    tolerance = 1e-12
  )
  expect_equal(
    unname(original_samples[, "(mu) sd(intercept)"]),
    rep(expected_sd[1], 2),
    tolerance = 1e-12
  )
  expect_equal(
    unname(original_samples[, "(mu) sd(x)"]),
    rep(expected_sd[2], 2),
    tolerance = 1e-12
  )
  expect_equal(
    unname(original_samples[, "(mu) cor(intercept,x)"]),
    rep(expected_cor, 2),
    tolerance = 1e-12
  )
  expect_equal(
    unname(original_table["(mu) cor(intercept,x)", "Mean"]),
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
  attr(formula_scale$mu, "random_effect_sd_leaves") <-
    .canonical_us_sd_leaves()

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
  attr(formula_scale$mu, "random_effect_sd_leaves") <-
    .canonical_us_sd_leaves()

  transformed <- transform_scale_samples(posterior, formula_scale)

  expect_true(is.na(transformed[1, "mu__xREx__id_xRE_CORx_R[1,2]"]))
  expect_true(is.na(transformed[1, "mu__xREx__id_xRE_CORx_L[2,1]"]))
  expect_false(is.na(transformed[2, "mu__xREx__id_xRE_CORx_R[1,2]"]))

  all_invalid <- transform_scale_samples(
    posterior[1, , drop = FALSE],
    formula_scale
  )
  all_invalid_cor_cols <- grep("_xRE_CORx_[RL]", colnames(all_invalid), value = TRUE)
  expect_true(all(is.na(all_invalid[, all_invalid_cor_cols, drop = FALSE])))

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

test_that("transform_scale_samples unscales random-factor correlations in column space", {

  df <- data.frame(
    x = 1:9,
    f = factor(rep(c("a", "b", "c"), 3), levels = c("a", "b", "c")),
    id = factor(rep(c("g1", "g2", "g3"), each = 3))
  )
  formula_result <- suppressWarnings(
    JAGS_formula(
      formula = ~ 1 + us(1 + x + f | id),
      parameter = "mu",
      data = df,
      prior_list = list(intercept = prior("normal", list(0, 1))),
      formula_scale = list(x = TRUE),
      prior_random = prior_random(
        id = random_block(
          sd = prior("gamma", list(2, 2)),
          terms = list(
            f = prior_factor("mnormal", list(0, 1), contrast = "orthonormal")
          ),
          cor = prior_lkj(eta = 1)
        )
      )
    )
  )

  random_term <- formula_result$formula_design$random_effects[[1]]
  expect_equal(
    random_term$sd_parameter_names[3:4],
    paste0("mu__xREx__id_f[", 1:2, "]")
  )
  sd_cols <- unique(random_term$sd_parameter_names)
  R_names <- outer(
    seq_len(random_term$n_columns),
    seq_len(random_term$n_columns),
    Vectorize(function(row, column){
      paste0("mu__xREx__id_xRE_CORx_R[", row, ",", column, "]")
    })
  )
  L_names <- outer(
    seq_len(random_term$n_columns),
    seq_len(random_term$n_columns),
    Vectorize(function(row, column){
      paste0("mu__xREx__id_xRE_CORx_L[", row, ",", column, "]")
    })
  )
  source_sd <- c(1.5, 2, 3, 3)
  source_cor <- matrix(
    c(
      1, 0.2, 0.3, 0.4,
      0.2, 1, 0.5, 0.6,
      0.3, 0.5, 1, 0.7,
      0.4, 0.6, 0.7, 1
    ),
    nrow = 4,
    byrow = TRUE
  )
  source_L <- t(chol(source_cor))
  posterior <- matrix(
    c(source_sd, as.vector(source_cor), as.vector(source_L)),
    nrow = 1,
    dimnames = list(NULL, c(sd_cols, as.vector(R_names), as.vector(L_names)))
  )

  transformed <- transform_scale_samples(
    posterior,
    list(mu = formula_result$formula_scale)
  )

  scale_info <- formula_result$formula_scale$mu_x
  M <- diag(random_term$n_columns)
  M[1, 2] <- -scale_info$mean / scale_info$sd
  M[2, 2] <- 1 / scale_info$sd
  expected_cov <- M %*% diag(source_sd, nrow = 4) %*% source_cor %*%
    diag(source_sd, nrow = 4) %*% t(M)
  expected_sd <- sqrt(diag(expected_cov))
  expected_cor <- expected_cov / tcrossprod(expected_sd)
  expected_L <- t(chol(expected_cor))

  expect_equal(
    unname(transformed[1, sd_cols]),
    expected_sd,
    tolerance = 1e-12
  )
  expect_equal(
    unname(transformed[1, as.vector(R_names)]),
    as.vector(expected_cor),
    tolerance = 1e-12
  )
  expect_equal(
    unname(transformed[1, as.vector(L_names)]),
    as.vector(expected_L),
    tolerance = 1e-12
  )
  expect_false(isTRUE(all.equal(source_cor[1, 4], expected_cor[1, 4])))
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

  study_prediction <- JAGS_evaluate_formula(
    fit = fit,
    formula = ~ 1 + x +
      random(1 | id, name = "study", covariance = "diag"),
    parameter = "mu",
    data = df,
    prior_list = formula_result$prior_list
  )

  expect_equal(unname(drop(study_prediction)), c(1.5, 3.5, 4.5, 6.5))
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

test_that("factor random slopes compile identically from data frames and tibbles", {

  skip_if_not_installed("tibble")

  df <- data.frame(
    f = factor(rep(c("a", "b", "c"), 2), levels = c("a", "b", "c")),
    id = factor(rep(c("g1", "g2"), each = 3), levels = c("g1", "g2"))
  )
  compile_formula <- function(data){
    JAGS_formula(
      formula = ~ 1 + diag(0 + f | id),
      parameter = "mu",
      data = data,
      prior_list = list(
        intercept = prior("normal", list(0, 1))
      ),
      formula_scale = TRUE,
      prior_random = prior_random(
        id = random_block(sd = prior("gamma", list(2, 2)))
      )
    )
  }

  data_frame_result <- compile_formula(df)
  tibble_result <- compile_formula(tibble::as_tibble(df))
  data_frame_term <- data_frame_result$formula_design$random_effects[[1L]]
  tibble_term <- tibble_result$formula_design$random_effects[[1L]]

  expect_identical(tibble_term$model_terms_type, data_frame_term$model_terms_type)
  expect_identical(tibble_term$model_terms_type[["f"]], "factor")
  expect_identical(tibble_term$sd_parameter_names, data_frame_term$sd_parameter_names)
  expect_equal(
    tibble_result$data$mu__xREx__id_xRE_DATAx,
    data_frame_result$data$mu__xREx__id_xRE_DATAx
  )
  expect_identical(names(tibble_result$prior_list), names(data_frame_result$prior_list))
  expect_null(tibble_result$formula_scale)
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
      study = random_block(sd = sd_prior, cor = prior_lkj(eta = 2)),
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
  expect_equal(result$jags_modules, "BayesTools")
  expect_match(result$formula_syntax, "dbt_lkj_cpc", fixed = TRUE)
  expect_false("backend" %in% names(result$formula_design$random_effects[[1]]$correlation))
  expect_false(any(grepl("dbeta", result$formula_syntax, fixed = TRUE)))

  inherited_covariance <- JAGS_formula(
    formula = ~ 1 + x + random(1 + x | id, name = "study"),
    parameter = "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    ),
    prior_random = prior_random(
      cor = prior_lkj(eta = 3),
      study = random_block(sd = sd_prior)
    )
  )
  expect_equal(inherited_covariance$jags_modules, "BayesTools")
  expect_match(inherited_covariance$formula_syntax, "_lkj_alpha[1] <- 3", fixed = TRUE)
  expect_match(inherited_covariance$formula_syntax, "dbt_lkj_cpc", fixed = TRUE)

  eta_override <- JAGS_formula(
    formula = ~ 1 + x + random(1 + x | id, name = "study"),
    parameter = "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    ),
    prior_random = prior_random(
      cor = prior_lkj(eta = 3),
      study = random_block(
        sd = sd_prior,
        covariance = random_covariance(eta = 4)
      )
    )
  )
  expect_equal(eta_override$jags_modules, "BayesTools")
  expect_match(eta_override$formula_syntax, "_lkj_alpha[1] <- 4", fixed = TRUE)

  expect_error(
    JAGS_formula(
      formula = ~ 1 + x + random(x | id, name = "study", covariance = "cs"),
      parameter = "mu",
      data = df,
      prior_list = list(
        intercept = prior("normal", list(0, 1)),
        x = prior("normal", list(0, 1))
      ),
      prior_random = prior_random(
        covariance = random_covariance(cor = prior("normal", list(0, 0.5))),
        study = random_block(sd = sd_prior, covariance = random_covariance(eta = 2))
      )
    ),
    "uses a scalar correlation prior",
    fixed = TRUE
  )

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
      cor = prior_lkj(eta = 3),
      drug_slope = random_block(covariance = random_covariance(structure = "diag"))
    )
  )
  expect_equal(mixed_covariance$formula_design$random_effects[[1]]$structure, "us")
  expect_equal(mixed_covariance$formula_design$random_effects[[2]]$structure, "diag")
  expect_equal(mixed_covariance$jags_modules, "BayesTools")
  expect_match(mixed_covariance$formula_syntax, "dbt_lkj_cpc", fixed = TRUE)
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
  expect_no_error(
    JAGS_formula(
      formula = ~ 1 + x,
      parameter = "mu",
      data = df,
      prior_list = list(
        intercept = prior("normal", list(0, 1)),
        x = prior("normal", list(0, 1))
      ),
      prior_random = prior_random(sd = sd_prior)
    )
  )
  expect_error(
    JAGS_formula(
      formula = ~ 1 + x,
      parameter = "mu",
      data = df,
      prior_list = list(
        intercept = prior("normal", list(0, 1)),
        x = prior("normal", list(0, 1))
      ),
      prior_random = prior_random(study = random_block(sd = sd_prior))
    ),
    "block override names were not found",
    fixed = TRUE
  )
  expect_error(
    JAGS_formula(
      formula = ~ 1 + x,
      parameter = "mu",
      data = df,
      prior_list = list(
        intercept = prior("normal", list(0, 1)),
        x = prior("normal", list(0, 1))
      ),
      prior_random = prior_random(
        allocation = random_variance_allocation(name = "allocation",
          terms = c("study", "drug"),
          sd = sd_prior
        )
      )
    ),
    "Variance allocation priors require formula random-effect terms",
    fixed = TRUE
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
  mixed_contrast_result <- JAGS_formula(
    formula = ~ 1 + diag(0 + f:g | id),
    parameter = "mu",
    data = mixed_contrast_df,
    prior_list = list(
      intercept = prior("normal", list(0, 1))
    ),
    prior_random = prior_random(
      id = random_block(sd = prior("gamma", list(2, 2)))
    )
  )
  expected_mixed <- stats::model.matrix(
    ~ f:g,
    data = mixed_contrast_df
  )[, -1, drop = FALSE]
  mixed_term <- mixed_contrast_result$formula_design$random_effects[[1L]]
  expect_equal(
    as.vector(mixed_term$model_matrix),
    as.vector(expected_mixed),
    tolerance = 1e-12
  )
  expect_equal(
    mixed_term$contrast_matrices$f,
    stats::contrasts(mixed_contrast_df$f)
  )
  expect_equal(
    mixed_term$contrast_matrices$g,
    stats::contrasts(mixed_contrast_df$g)
  )
})

test_that("fixed and random blocks own independent concrete factor bases", {

  data <- data.frame(
    f = factor(rep(c("a", "b", "c"), 4)),
    id_treatment = factor(rep(c("s1", "s2", "s3", "s4"), each = 3)),
    id_meandif = factor(rep(c("q1", "q2", "q3"), 4))
  )
  inherited <- JAGS_formula(
    formula = ~ 0 + f + diag(0 + f | id_treatment),
    parameter = "mu",
    data = data,
    prior_list = list(
      f = prior_factor(
        "normal",
        list(0, 1),
        contrast = "independent"
      )
    ),
    prior_random = prior_random(
      id_treatment = random_block(sd = prior("gamma", list(2, 2)))
    )
  )
  inherited_term <- inherited$formula_design$random_effects[[1L]]

  expect_equal(
    inherited$formula_design$contrast_matrices$f,
    contr.independent(levels(data$f))
  )
  expect_equal(
    inherited_term$contrast_matrices$f,
    contr.independent(levels(data$f))
  )
  expect_equal(inherited_term$n_columns, nlevels(data$f))

  result <- JAGS_formula(
    formula = ~ f +
      diag(0 + f | id_treatment) +
      diag(0 + f | id_meandif),
    parameter = "mu",
    data = data,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      f = prior_factor(
        "mnormal",
        list(0, 1),
        contrast = "orthonormal"
      )
    ),
    prior_random = prior_random(
      id_treatment = random_block(
        sd = prior("gamma", list(2, 2)),
        contrasts = c(f = "treatment")
      ),
      id_meandif = random_block(
        sd = prior("gamma", list(2, 2)),
        contrasts = c(f = "meandif")
      )
    )
  )

  fixed_matrix <- contr.orthonormal(levels(data$f))
  treatment_matrix <- stats::contr.treatment(levels(data$f))
  meandif_matrix <- contr.meandif(levels(data$f))
  treatment_term <- result$formula_design$random_effects[[1L]]
  meandif_term <- result$formula_design$random_effects[[2L]]

  expect_identical(result$formula_design$schema_version, 4L)
  expect_equal(result$formula_design$contrast_matrices$f, fixed_matrix)
  expect_equal(treatment_term$contrast_matrices$f, treatment_matrix)
  expect_equal(meandif_term$contrast_matrices$f, meandif_matrix)
  expect_identical(treatment_term$contrast_owner, "random_block")
  expect_identical(meandif_term$contrast_owner, "random_block")
  expect_false(isTRUE(all.equal(fixed_matrix, treatment_matrix)))
  expect_false(isTRUE(all.equal(fixed_matrix, meandif_matrix)))
  expect_identical(
    meandif_term$sd_parameter_names,
    paste0("mu__xREx__id_meandif_f[", 1:2, "]")
  )
  expect_identical(
    meandif_term$sd_leaves$leaf_terms_by_column,
    paste0("f[", 1:2, "]")
  )

  treatment_prediction <- .bt_random_effect_prediction_data(
    treatment_term,
    data
  )
  meandif_prediction <- .bt_random_effect_prediction_data(
    meandif_term,
    data
  )
  expect_equal(
    treatment_prediction$model_matrix,
    treatment_term$model_matrix
  )
  expect_equal(
    meandif_prediction$model_matrix,
    meandif_term$model_matrix
  )

  missing_basis <- treatment_term
  missing_basis$contrast_matrices <- NULL
  expect_error(
    .bt_random_effect_prediction_data(missing_basis, data),
    "missing its owner-scoped concrete factor basis",
    fixed = TRUE
  )
})

test_that("random-block contrast specifications validate their design scope", {

  expect_error(
    random_block(contrasts = "treatment"),
    "must be a named character vector",
    fixed = TRUE
  )
  expect_error(
    random_block(contrasts = c(f = "sum")),
    "Unknown random-block contrast value",
    fixed = TRUE
  )

  data <- data.frame(
    f = factor(rep(c("a", "b", "c"), 2)),
    id = factor(rep(c("s1", "s2"), each = 3))
  )
  expect_error(
    JAGS_formula(
      ~ 1 + diag(0 + f | id),
      parameter = "mu",
      data = data,
      prior_list = list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(
        id = random_block(
          sd = prior("gamma", list(2, 2)),
          contrasts = c(unknown = "treatment")
        )
      )
    ),
    "contrast overrides reference predictors outside this design",
    fixed = TRUE
  )
  expect_error(
    JAGS_formula(
      ~ 1 + cs(f | id),
      parameter = "mu",
      data = data,
      prior_list = list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(
        id = random_block(
          sd = prior("gamma", list(2, 2)),
          cor = prior("normal", list(0, 0.5)),
          contrasts = c(f = "orthonormal")
        )
      )
    ),
    "level basis is defined by the covariance structure",
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

  large_factor_df <- data.frame(
    f = factor(sprintf("level_%03d", seq_len(113L))),
    id = factor(rep(sprintf("group_%02d", seq_len(17L)), length.out = 113L))
  )
  large_result <- JAGS_formula(
      formula = ~ 1 + cs(f | id),
      parameter = "mu",
      data = large_factor_df,
      prior_list = list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(
        id = random_block(sd = sd_prior, cor = prior("normal", list(0, 0.5)))
      )
    )
  large_term <- large_result$formula_design$random_effects[[1]]
  expect_s3_class(
    large_term$latent_layout,
    "BayesTools_random_effect_structured_local_layout"
  )
  expect_equal(large_term$latent_layout$n_local, 113L)
  expect_equal(
    names(large_result$data),
    c("N_mu", "mu__xREx__id_xRE_MAPx", "mu__xREx__id_xRE_COLx")
  )
  expect_false(grepl("_xRE_CORx_L", large_result$formula_syntax, fixed = TRUE))
  expect_false(grepl("_xRE_CORx_R", large_result$formula_syntax, fixed = TRUE))

  coefficient_result <- JAGS_formula(
    formula = ~ 1 + cs(f | id),
    parameter = "mu",
    data = large_factor_df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(
        sd = sd_prior,
        cor = prior("normal", list(0, 0.5)),
        monitor = random_monitor(latent = FALSE, coefficients = TRUE)
      )
    )
  )
  coefficient_term <- coefficient_result$formula_design$random_effects[[1]]
  expect_null(coefficient_term$latent_layout)
  expect_true("mu__xREx__id_xRE_COEFx" %in% coefficient_result$add_parameters)
  expect_match(
    coefficient_result$formula_syntax,
    "mu__xREx__id_xRE_Zx[i,j] ~ dnorm(0, 1)",
    fixed = TRUE
  )

  old_multiplier <- getOption("BayesTools.random_effects_complexity_multiplier")
  on.exit(options(BayesTools.random_effects_complexity_multiplier = old_multiplier), add = TRUE)
  options(BayesTools.random_effects_complexity_multiplier = 0)
  expect_error(
    BayesTools:::.bt_random_effect_check_dense_complexity(
      random_term = list(block_name = "id"),
      structure = "cs",
      n_groups = 17L,
      n_columns = 113L,
      n_rows = 113L,
      monitor_policy = random_monitor()
    ),
    "must be a positive numeric scalar",
    fixed = TRUE
  )
  options(BayesTools.random_effects_complexity_multiplier = old_multiplier)

  options(BayesTools.random_effects_complexity_multiplier = Inf)
  expect_null(
    BayesTools:::.bt_random_effect_check_dense_complexity(
      random_term = list(block_name = "id"),
      structure = "cs",
      n_groups = 17L,
      n_columns = 113L,
      n_rows = 113L,
      monitor_policy = random_monitor()
    )
  )
  options(BayesTools.random_effects_complexity_multiplier = old_multiplier)

  default_cs <- JAGS_formula(
    formula = ~ 1 + cs(f | id),
    parameter = "mu",
    data = factor_df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(id = random_block(sd = sd_prior))
  )
  expect_equal(names(default_cs$prior_list), c(
    "mu_intercept", "mu__xREx__id_sd", "mu__xREx__id_rho"
  ))
  expect_equal(default_cs$prior_list$mu__xREx__id_rho$distribution, "uniform")
  expect_equal(
    default_cs$prior_list$mu__xREx__id_rho$parameters,
    list(a = -0.5, b = 1)
  )
  expect_equal(
    default_cs$formula_design$random_effects[[1]]$correlation$rho_scale,
    "rho"
  )

  fisher_cs <- JAGS_formula(
    formula = ~ 1 + cs(f | id),
    parameter = "mu",
    data = factor_df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(sd = sd_prior, cor = prior("normal", list(0, 0.5)))
    )
  )

  expect_equal(names(fisher_cs$prior_list), c("mu_intercept", "mu__xREx__id_sd", "mu__xREx__id_rho_z"))
  expect_equal(fisher_cs$prior_list$mu__xREx__id_rho_z$distribution, "normal")
  expect_equal(fisher_cs$prior_list$mu__xREx__id_rho_z$parameters, list(mean = 0, sd = 0.5))
  expect_true(fisher_cs$formula_design$random_effects[[1]]$homogeneous_sd)
  expect_equal(fisher_cs$formula_design$random_effects[[1]]$column_names,
               c("fa", "fb", "fc"))
  expect_null(fisher_cs$data$mu__xREx__id_xRE_DATAx)
  expect_equal(
    fisher_cs$formula_design$random_effects[[1]]$sd_parameter_names,
    rep("mu__xREx__id_sd", 3)
  )

  expect_error(
    JAGS_formula(
      formula = ~ 1 + cs(1 | id),
      parameter = "mu",
      data = continuous_df,
      prior_list = list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(
        id = random_block(sd = sd_prior, cor = prior("normal", list(0, 0.5)))
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
      id = random_block(sd = sd_prior, cor = prior("normal", list(0, 0.5)))
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

  wrapper_hcs_result <- JAGS_formula(
    formula = ~ 1 + random(f | id, covariance = "cs", hom = FALSE),
    parameter = "mu",
    data = factor_df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(sd = sd_prior, cor = prior("normal", list(0, 0.5)))
    )
  )
  expect_equal(wrapper_hcs_result$formula_design$random_effects[[1]]$structure, "hcs")
  expect_false(wrapper_hcs_result$formula_design$random_effects[[1]]$homogeneous_sd)

  hcs_composite <- JAGS_formula(
    formula = ~ 1 + hcs(f + g | id),
    parameter = "mu",
    data = factor_df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(sd = sd_prior, cor = prior("normal", list(0, 0.5)))
    )
  )
  expect_equal(hcs_composite$formula_design$random_effects[[1]]$structured_index$variables, c("f", "g"))
  expect_equal(hcs_composite$formula_design$random_effects[[1]]$structured_index$label, "f:g")
  expect_equal(hcs_composite$formula_design$random_effects[[1]]$n_columns, 6L)
  expect_null(hcs_composite$data$mu__xREx__id_xRE_DATAx)
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
  hcs_composite_raw <- .parameter_catalog_random_summary_samples(
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
  expect_true("(mu) id: sd(f:g[a.u])" %in% hcs_composite_raw_display)
  expect_true("(mu) id: cor" %in% hcs_composite_raw_display)
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
        target = "sd_component",
        sd = sd_prior,
        weights = prior("dirichlet", list(alpha = rep(1, 6)))
      ),
      id = random_block(cor = prior("normal", list(0, 0.5)))
    )
  )
  malformed_random_term <- hcs_composite_allocation$formula_design$random_effects[[1]]
  malformed_random_term$sd_binding$allocations[[1L]]$parent_factors <- NULL
  expect_error(
    BayesTools:::.bt_random_effect_sd_draws(
      random_term = malformed_random_term,
      n_columns = malformed_random_term$n_columns,
      posterior = matrix(
        rep(1, 7),
        nrow = 1,
        dimnames = list(
          NULL,
          c(
            "mu__xRE_ALLOCx_leaf_alloc__allocation_sd",
            paste0("prior_par_eta_mu__xRE_ALLOCx_leaf_alloc__weight[", 1:6, "]")
          )
        )
      ),
      prior_list = hcs_composite_allocation$prior_list
    ),
    "missing 'allocation\\$parent_factors'"
  )

  malformed_random_term <- hcs_composite_allocation$formula_design$random_effects[[1]]
  malformed_random_term$sd_binding$allocations[[1L]]$n_targets <- 6.5
  expect_error(
    BayesTools:::.bt_random_effect_sd_draws(
      random_term = malformed_random_term,
      n_columns = malformed_random_term$n_columns,
      posterior = matrix(
        rep(1, 7),
        nrow = 1,
        dimnames = list(
          NULL,
          c(
            "mu__xRE_ALLOCx_leaf_alloc__allocation_sd",
            paste0("prior_par_eta_mu__xRE_ALLOCx_leaf_alloc__weight[", 1:6, "]")
          )
        )
      ),
      prior_list = hcs_composite_allocation$prior_list
    ),
    "missing canonical 'allocation\\$n_targets'"
  )

  malformed_random_term <- hcs_composite_allocation$formula_design$random_effects[[1]]
  malformed_random_term$sd_binding$allocations[[1L]]$leaf_index_by_column <-
    malformed_random_term$sd_binding$allocations[[1L]]$leaf_index_by_column[-1L]
  expect_error(
    BayesTools:::.bt_random_effect_sd_draws(
      random_term = malformed_random_term,
      n_columns = malformed_random_term$n_columns,
      posterior = matrix(
        rep(1, 7),
        nrow = 1,
        dimnames = list(
          NULL,
          c(
            "mu__xRE_ALLOCx_leaf_alloc__allocation_sd",
            paste0("prior_par_eta_mu__xRE_ALLOCx_leaf_alloc__weight[", 1:6, "]")
          )
        )
      ),
      prior_list = hcs_composite_allocation$prior_list
    ),
    "leaf_index_by_column"
  )
  malformed_random_term <- hcs_composite_allocation$formula_design$random_effects[[1]]
  malformed_random_term$sd_binding$allocations[[1L]]$n_targets <- 1L
  malformed_random_term$sd_binding$allocations[[1L]]$leaf_index_by_column <-
    rep(1L, malformed_random_term$n_columns)
  malformed_random_term$sd_binding$allocations[[1L]]$leaf_names <- "only"
  malformed_random_term$sd_binding$allocations[[1L]]$leaf_terms <- "only"
  expect_error(
    BayesTools:::.bt_random_effect_sd_draws(
      random_term = malformed_random_term,
      n_columns = malformed_random_term$n_columns,
      posterior = matrix(
        rep(1, 7),
        nrow = 1,
        dimnames = list(
          NULL,
          c(
            "mu__xRE_ALLOCx_leaf_alloc__allocation_sd",
            paste0("prior_par_eta_mu__xRE_ALLOCx_leaf_alloc__weight[", 1:6, "]")
          )
        )
      ),
      prior_list = hcs_composite_allocation$prior_list
    ),
    "missing canonical 'allocation\\$n_targets'"
  )
  malformed_random_term <- hcs_composite_allocation$formula_design$random_effects[[1]]
  malformed_random_term$sd_binding$allocations[[1L]]$leaf_index_by_column <-
    rep(1L, malformed_random_term$n_columns)
  expect_error(
    BayesTools:::.bt_random_effect_sd_draws(
      random_term = malformed_random_term,
      n_columns = malformed_random_term$n_columns,
      posterior = matrix(
        rep(1, 7),
        nrow = 1,
        dimnames = list(
          NULL,
          c(
            "mu__xRE_ALLOCx_leaf_alloc__allocation_sd",
            paste0("prior_par_eta_mu__xRE_ALLOCx_leaf_alloc__weight[", 1:6, "]")
          )
        )
      ),
      prior_list = hcs_composite_allocation$prior_list
    ),
    "must cover every target",
    fixed = TRUE
  )
  malformed_random_term <- hcs_composite_allocation$formula_design$random_effects[[1]]
  malformed_random_term$sd_binding$allocations[[1L]]$leaf_names[[2L]] <-
    malformed_random_term$sd_binding$allocations[[1L]]$leaf_names[[1L]]
  expect_error(
    BayesTools:::.bt_random_effect_sd_draws(
      random_term = malformed_random_term,
      n_columns = malformed_random_term$n_columns,
      posterior = matrix(
        rep(1, 7),
        nrow = 1,
        dimnames = list(
          NULL,
          c(
            "mu__xRE_ALLOCx_leaf_alloc__allocation_sd",
            paste0("prior_par_eta_mu__xRE_ALLOCx_leaf_alloc__weight[", 1:6, "]")
          )
        )
      ),
      prior_list = hcs_composite_allocation$prior_list
    ),
    "missing canonical 'allocation\\$leaf_names'"
  )
  malformed_random_term <- hcs_composite_allocation$formula_design$random_effects[[1]]
  malformed_random_term$sd_binding$allocations[[1L]]$leaf_terms <- NULL
  expect_error(
    BayesTools:::.bt_random_effect_sd_draws(
      random_term = malformed_random_term,
      n_columns = malformed_random_term$n_columns,
      posterior = matrix(
        rep(1, 7),
        nrow = 1,
        dimnames = list(
          NULL,
          c(
            "mu__xRE_ALLOCx_leaf_alloc__allocation_sd",
            paste0("prior_par_eta_mu__xRE_ALLOCx_leaf_alloc__weight[", 1:6, "]")
          )
        )
      ),
      prior_list = hcs_composite_allocation$prior_list
    ),
    "missing canonical 'allocation\\$leaf_terms'"
  )

  ar1_result <- JAGS_formula(
    formula = ~ 1 + ar1(f | id),
    parameter = "mu",
    data = factor_df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(sd = sd_prior, cor = prior("normal", list(0, 0.5)))
    )
  )

  expect_equal(ar1_result$formula_design$random_effects[[1]]$column_names,
               c("fa", "fb", "fc"))
  expect_null(ar1_result$data$mu__xREx__id_xRE_DATAx)
  expect_equal(names(ar1_result$prior_list), c("mu_intercept", "mu__xREx__id_sd", "mu__xREx__id_rho_z"))
  expect_true(ar1_result$formula_design$random_effects[[1]]$homogeneous_sd)
  expect_equal(ar1_result$formula_design$random_effects[[1]]$structure, "ar1")
  expect_equal(
    ar1_result$formula_design$random_effects[[1]]$sd_parameter_names,
    rep("mu__xREx__id_sd", 3)
  )
  expect_match(
    ar1_result$formula_syntax,
    "mu__xREx__id_xRE_UNIT_COEFx[g,i - 1]",
    fixed = TRUE
  )

  ar1_implicit_levels <- JAGS_formula(
    formula = ~ 1 + ar1(f | id),
    parameter = "mu",
    data = factor_df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(sd = sd_prior, cor = prior("normal", list(0, 0.5)))
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

  wrapper_har_result <- JAGS_formula(
    formula = ~ 1 + re(f | id, covariance = "ar1", hom = FALSE),
    parameter = "mu",
    data = factor_df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(sd = sd_prior, cor = prior("normal", list(0, 0.5)))
    )
  )
  expect_equal(wrapper_har_result$formula_design$random_effects[[1]]$structure, "har")
  expect_false(wrapper_har_result$formula_design$random_effects[[1]]$homogeneous_sd)

  hcs_result <- JAGS_formula(
    formula = ~ 1 + hcs(f | id),
    parameter = "mu",
    data = factor_df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(sd = sd_prior, cor = prior("normal", list(0, 0.5)))
    )
  )

  expect_equal(hcs_result$formula_design$random_effects[[1]]$column_names,
               c("fa", "fb", "fc"))
  expect_null(hcs_result$data$mu__xREx__id_xRE_DATAx)
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
  hcs_implicit_levels <- JAGS_formula(
    formula = ~ 1 + hcs(f | id),
    parameter = "mu",
    data = factor_df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(sd = sd_prior, cor = prior("normal", list(0, 0.5)))
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
      id = random_block(sd = sd_prior, cor = prior("normal", list(0, 0.5)))
    )
  )
  expect_equal(har_implicit_levels$formula_design$random_effects[[1]]$column_names,
               c("fa", "fb", "fc"))
  expect_null(har_implicit_levels$data$mu__xREx__id_xRE_DATAx)
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
          cor = prior("normal", list(0, 0.5))
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
        cor = prior("normal", list(0, 0.5))
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
            cor = prior("normal", list(0, 0.5)),
            cor_scale = "cor"
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
  expect_match(
    raw_result$formula_syntax,
    "mu__xREx__id_xRE_AR_PHIX[i] * mu__xREx__id_xRE_UNIT_COEFx[g,i - 1]",
    fixed = TRUE
  )

  raw_inherited <- NULL
  expect_warning(
    raw_inherited <- JAGS_formula(
      formula = ~ 1 + har(f | id),
      parameter = "mu",
      data = factor_df,
      prior_list = list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(
        covariance = random_covariance(
          cor = prior("normal", list(0, 0.5)),
          cor_scale = "cor"
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
        covariance = random_covariance(cor = prior("normal", list(0, 0.5))),
        id = random_block(
          sd = sd_prior,
          covariance = random_covariance(structure = "har", cor_scale = "cor")
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
        covariance = random_covariance(cor_scale = "cor"),
        id = random_block(
          sd = sd_prior,
          cor = prior("normal", list(0, 0.5))
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
      id = random_block(sd = sd_prior, cor = prior("normal", list(0, 0.5)))
    )
  )
  car_term <- car_result$formula_design$random_effects[[1]]
  car_distance <- abs(outer(c(0, 0.5, 2), c(0, 0.5, 2), "-"))

  expect_equal(car_term$column_names, c("time_0", "time_0p5", "time_2"))
  expect_null(car_result$data$mu__xREx__id_xRE_DATAx)
  expect_equal(names(car_result$prior_list), c("mu_intercept", "mu__xREx__id_sd", "mu__xREx__id_rho_z"))
  expect_true(car_term$homogeneous_sd)
  expect_equal(car_term$structure, "car")
  expect_equal(car_term$car$time_variable, "time")
  expect_equal(car_term$car$time_values, c(0, 0.5, 2))
  expect_null(car_term$car$distance_matrix)
  expect_equal(car_term$correlation$bounds, c(lower = 0, upper = 1))
  expect_null(car_term$correlation$distance_matrix)
  expect_equal(car_term$correlation$time_values, c(0, 0.5, 2))
  expect_equal(car_result$prior_list$mu__xREx__id_rho_z$truncation$lower, 0)
  expect_equal(car_result$prior_list$mu__xREx__id_rho_z$truncation$upper, Inf)
  expect_match(
    car_result$formula_syntax,
    "mu__xREx__id_rho <- tanh(mu__xREx__id_rho_z)",
    fixed = TRUE
  )
  expect_match(
    car_result$formula_syntax,
    "tanh(mu__xREx__id_rho_z)",
    fixed = TRUE
  )
  expect_match(
    car_result$formula_syntax,
    paste0(
      "mu__xREx__id_xRE_CAR_LOG_PHIX[2] <- ",
      "0.5 * log(mu__xREx__id_rho)"
    ),
    fixed = TRUE
  )
  expect_match(
    car_result$formula_syntax,
    paste0(
      "mu__xREx__id_xRE_CAR_LOG_PHIX[3] <- ",
      "1.5 * log(mu__xREx__id_rho)"
    ),
    fixed = TRUE
  )
  expect_match(
    car_result$formula_syntax,
    paste0(
      "mu__xREx__id_xRE_CAR_PHIX[2] <- ",
      "exp(mu__xREx__id_xRE_CAR_LOG_PHIX[2])"
    ),
    fixed = TRUE
  )
  expect_match(
    car_result$formula_syntax,
    paste0(
      "mu__xREx__id_xRE_CAR_INNOV_VARx[2] <- ",
      "pexp(-2 * mu__xREx__id_xRE_CAR_LOG_PHIX[2], 1)"
    ),
    fixed = TRUE
  )
  expect_false(grepl(
    "pow(mu__xREx__id_rho",
    car_result$formula_syntax,
    fixed = TRUE
  ))
  expect_false(grepl(
    "1 - pow(mu__xREx__id_xRE_CAR_PHIX",
    car_result$formula_syntax,
    fixed = TRUE
  ))

  car_independent <- JAGS_formula(
    formula = ~ 1 + car(0 + time | id),
    parameter = "mu",
    data = car_df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(sd = sd_prior, cor = prior("point", list(location = 0)))
    )
  )
  car_independent_term <- car_independent$formula_design$random_effects[[1]]
  empty_posterior <- matrix(numeric(0), nrow = 2, ncol = 0)
  expect_equal(
    BayesTools:::.bt_random_effect_rho_draws(car_independent_term, empty_posterior),
    rep(
      BayesTools:::.bt_random_effect_representable_rho_bounds(
        car_independent_term$correlation$bounds,
        car_independent_term$structure
      )[["lower"]],
      2L
    )
  )
  expect_equal(
    bayestools_reference_random_effect_scalar_rho_support(
      empty_posterior,
      car_independent_term
    ),
    0
  )

  car_independent_raw <- JAGS_formula(
    formula = ~ 1 + car(0 + time | id),
    parameter = "mu",
    data = car_df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(
        sd = sd_prior,
        covariance = random_covariance(
          cor = prior("point", list(location = 0)),
          cor_scale = "cor"
        )
      )
    )
  )
  expect_equal(
    car_independent_raw$formula_design$random_effects[[1]]$correlation$sample_fixed,
    0
  )

  car_explicit_no_intercept <- JAGS_formula(
    formula = ~ 1 + car(0 + time | id),
    parameter = "mu",
    data = car_df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(sd = sd_prior, cor = prior("normal", list(0, 0.5)))
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
      id = random_block(sd = sd_prior, cor = prior("normal", list(0, 0.5)))
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
      id = random_block(sd = sd_prior, cor = prior("normal", list(0, 0.5)))
    )
  )
  expect_equal(fixed_car_scale$formula_design$random_effects[[1]]$car$time_values, c(0, 0.5, 2))
  expect_null(fixed_car_scale$formula_design$random_effects[[1]]$car$distance_matrix)
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
        id = random_block(sd = sd_prior, cor = prior("normal", list(0, 0.5)))
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
        id = random_block(sd = sd_prior, cor = prior("normal", list(0, 0.5)))
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
      id = random_block(sd = sd_prior, cor = prior("normal", list(0, 0.5)))
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
          cor = prior("normal", list(0, 0.5)),
          cor_scale = "logit"
        )
      )
    )
  )
  expect_equal(names(car_logit$prior_list), c("mu_intercept", "mu__xREx__id_sd", "mu__xREx__id_rho_logit"))
  expect_match(car_logit$formula_syntax, "mu__xREx__id_rho <- 0 +", fixed = TRUE)
  expect_match(car_logit$formula_syntax, "* ilogit(mu__xREx__id_rho_logit)", fixed = TRUE)

  ar_logit <- JAGS_formula(
    formula = ~ 1 + ar1(f | id),
    parameter = "mu",
    data = factor_df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(
        sd = sd_prior,
        covariance = random_covariance(
          cor = prior("normal", list(0, 0.5)),
          cor_scale = "logit"
        )
      )
    )
  )
  expect_equal(names(ar_logit$prior_list), c("mu_intercept", "mu__xREx__id_sd", "mu__xREx__id_rho_logit"))
  expect_match(ar_logit$formula_syntax, "mu__xREx__id_rho <- -1 + 2", fixed = TRUE)
  expect_match(ar_logit$formula_syntax, "* ilogit(mu__xREx__id_rho_logit)", fixed = TRUE)

  expect_error(
    JAGS_formula(
      formula = ~ 1 + car(0 + x + time | id),
      parameter = "mu",
      data = transform(car_df, x = seq_len(nrow(car_df))),
      prior_list = list(intercept = prior("normal", list(0, 1))),
      prior_random = prior_random(
        id = random_block(sd = sd_prior, cor = prior("normal", list(0, 0.5)))
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
        id = random_block(sd = sd_prior, cor = prior("normal", list(0, 0.5)))
      )
    ),
    "must be numeric or an ordered factor",
    fixed = TRUE
  )
})

test_that("structured correlation Cholesky syntax exposes intended covariance patterns", {

  sd_prior <- prior("normal", list(0, 1), truncation = list(lower = 0, upper = Inf))
  rho_prior <- prior("normal", list(0, 0.5))
  block_prior <- random_block(sd = sd_prior, cor = rho_prior)

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
  expect_error(
    BayesTools:::.bt_JAGS_structured_corr_cholesky(
      node_prefix = "mu__xREx__id",
      prior_prefix = "_xREx__id",
      K = 3,
      structure = "car",
      block_prior = block_prior,
      include_correlation = TRUE,
      distance_matrix = car_distance
    ),
    "Dense CAR Cholesky compilation is unsupported",
    fixed = TRUE
  )
  car_module <- BayesTools:::.bt_JAGS_structured_corr_direct(
    node_prefix = "mu__xREx__id",
    prior_prefix = "_xREx__id",
    K = 3,
    structure = "car",
    block_prior = block_prior,
    include_correlation = TRUE,
    distance_matrix = car_distance
  )
  expect_match(
    car_module$syntax,
    paste0(
      "mu__xREx__id_xRE_CORx_R[1,2] <- ",
      "exp(1.5 * log(mu__xREx__id_rho))"
    ),
    fixed = TRUE
  )
  expect_match(
    car_module$syntax,
    paste0(
      "mu__xREx__id_xRE_CORx_R[1,3] <- ",
      "exp(3 * log(mu__xREx__id_rho))"
    ),
    fixed = TRUE
  )
  expect_false(grepl(
    "pow(mu__xREx__id_rho",
    car_module$syntax,
    fixed = TRUE
  ))
  expect_false(grepl("_xRE_CORx_L", car_module$syntax, fixed = TRUE))
  expect_null(car_module$cholesky_name)
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
        cor = rho_prior,
        cor_scale = "logit"
      )
    ),
    include_correlation = TRUE
  )
  expect_match(logit_module$syntax, "mu__xREx__id_rho <- -0.5 + 1.5", fixed = TRUE)
  expect_match(logit_module$syntax, "* ilogit(mu__xREx__id_rho_logit)", fixed = TRUE)
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
        cor = prior("normal", list(0, 0.5))
      )
    )
  )
  ar1_lme4 <- lme4::lFormula(y ~ 1 + ar1(0 + f | id), data = df_ar1)

  expect_equal(ar1_result$formula_design$random_effects[[1]]$column_names, ar1_lme4$reTrms$cnms$id)
  expect_equal(ar1_result$data$mu__xREx__id_xRE_MAPx, unname(as.integer(ar1_lme4$reTrms$flist$id)))
  expect_equal(ar1_result$formula_design$random_effects[[1]]$column_names,
               c("fa", "fb", "fc"))
  expect_null(ar1_result$data$mu__xREx__id_xRE_DATAx)
  expect_equal(ar1_result$formula_design$random_effects[[1]]$structure, "ar1")
  expect_true(ar1_result$formula_design$random_effects[[1]]$homogeneous_sd)
})
