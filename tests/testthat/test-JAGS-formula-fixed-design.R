skip_if_not_test_profile("unit")

# ============================================================================ #
# TEST FILE: JAGS Formula - Fixed-effect formula design and validation
# ============================================================================ #
#
# PURPOSE:
#   JAGS_formula() design metadata against stats::model.matrix(), contrast resolution,
#   and the validation of fixed-effect formulas, priors and inputs.
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
  random_result <- NULL
  expect_no_warning(random_result <- JAGS_formula(
    formula   = ~ 1 + diag(0 + group | id),
    parameter = "mu",
    data      = prediction_data,
    prior_list = list(
      intercept = prior("normal", list(0, 1))
    ),
    prior_random = prior_random(
      id = random_block(
        sd = prior("gamma", list(2, 2)),
        contrasts = c(group = "orthonormal")
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
  attr(posterior, "formula_design") <- list(mu = formula_result$formula_design)
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
  attr(fixed_posterior, "formula_design") <- list(mu = valid_fixed$formula_design)
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
  attr(posterior, "formula_design") <- list(
    mu = JAGS_formula(~ x, "mu", data, prior_list)$formula_design
  )
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
  attr(fit, "formula_design") <- list(mu = valid_result$formula_design)
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
  replay_result <- JAGS_formula(
    log_formula,
    "mu",
    data,
    list(
      intercept = prior("point", list(1)),
      x = prior("normal", list(0, 1))
    )
  )
  attr(fitted, "formula_design") <- list(mu = replay_result$formula_design)
  invalid_replay_priors <- replay_result$prior_list
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

test_that("dot expansion is refused in the package's words and an empty formula_scale scales nothing", {

  data <- data.frame(x = c(-1, 0, 1, 2), g = factor(c("a", "b", "a", "b")))
  priors <- list(intercept = prior("normal", list(0, 1)),
                 x = prior("normal", list(0, 1)))
  for(formula in list(~ ., ~ x + .)){
    expect_error(
      JAGS_formula(formula, "mu", data, priors),
      "Unsupported fixed-formula term '.'. Dot expansion is not supported; list each data-frame column explicitly.",
      fixed = TRUE
    )
  }
  scaled_none <- JAGS_formula(~ x, "mu", data, priors, formula_scale = list())
  unscaled <- JAGS_formula(~ x, "mu", data, priors)
  expect_identical(scaled_none$data, unscaled$data)
  expect_length(scaled_none$formula_scale, 0L)
})
