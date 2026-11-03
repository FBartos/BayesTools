skip_if_not_test_profile("unit")

# ============================================================================ #
# TEST FILE: JAGS Formula - Formula expressions
# ============================================================================ #
#
# PURPOSE:
#   Expression terms in formulas: structural parsing, array and paired expressions,
#   offsets, and replay of sampled indexed parameters.
# ============================================================================ #

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

test_that("formula array expressions pair coordinates pointwise", {
  X <- matrix(1:6, nrow = 3L)
  column <- c(1, 2, 1)
  evaluate <- function(formula, model_data){
    draws <- JAGS_formula_draws(
      matrix(0, 2, 1, dimnames = list(NULL, "mu_intercept")),
      formula, "mu", data.frame(dummy = 1:3),
      list(intercept = prior("point", list(0))),
      model_data = model_data
    )
    unname(JAGS_evaluate_formula(draws, parameter = "mu"))
  }
  expect_equal(
    evaluate(~ 1 + expression(X[i, column[i]]), list(X = X, column = column)),
    matrix(c(1, 5, 3), 3, 2)
  )
  expect_equal(evaluate(~ 1 + expression(X[i, 2]), list(X = X)),
               matrix(c(4, 5, 6), 3, 2))
  expect_equal(
    evaluate(~ 1 + expression(X[2, column[i]]), list(X = X, column = column)),
    matrix(c(2, 5, 2), 3, 2)
  )
  expect_equal(evaluate(~ 1 + expression(X[2, 2]), list(X = X)),
               matrix(5, 3, 2))
  expect_equal(
    evaluate(~ 1 + expression(2 * X[i, column[i]] + column[i] - 1),
             list(X = X, column = column)),
    matrix(c(2, 11, 6), 3, 2)
  )
  expect_equal(
    evaluate(~ 1 + expression(X[i, column[i]]),
             list(X = X, column = rep(TRUE, 3))),
    matrix(c(1, 2, 3), 3, 2)
  )
  expect_equal(evaluate(~ 1 + expression(column[i]), list(column = column)),
               matrix(column, 3, 2))

  A <- array(seq_len(12), dim = c(3L, 2L, 2L))
  layer <- c(2L, 1L, 2L)
  expected <- vapply(seq_len(3L), function(row){
    A[row, column[row], layer[row]]
  }, numeric(1))
  expect_equal(
    evaluate(~ 1 + expression(A[i, column[i], layer[i]]),
             list(A = A, column = column, layer = layer)),
    matrix(expected, 3, 2)
  )
})

test_that("paired array expressions use explicit prediction and model row order", {
  formula <- ~ 1 + expression(X[i, column[i]])
  X <- matrix(1:6, nrow = 3L)
  column <- c(1L, 2L, 1L)
  samples <- matrix(0, 2, 1, dimnames = list(NULL, "mu_intercept"))
  priors <- list(intercept = prior("point", list(0)))
  draws <- JAGS_formula_draws(samples, formula, "mu", data.frame(dummy = 1:3),
    priors, model_data = list(X = X, column = column))
  rows <- c(3L, 1L, 3L, 2L)
  prediction_data <- data.frame(dummy = seq_along(rows),
    X = I(X[rows, , drop = FALSE]), column = column[rows])
  expected <- matrix(c(3, 1, 3, 5), 4, 2)
  expect_equal(
    unname(JAGS_evaluate_formula(draws, parameter = "mu", data = prediction_data)),
    expected
  )
  reordered <- JAGS_formula_draws(
    samples, formula, "mu", data.frame(dummy = seq_along(rows)), priors,
    model_data = list(X = X[rows, , drop = FALSE], column = column[rows])
  )
  expect_equal(unname(JAGS_evaluate_formula(reordered, parameter = "mu")), expected)
})

test_that("paired array expressions refuse invalid coordinates and geometry", {
  X <- matrix(1:6, nrow = 3L)
  invalid <- c("X[i, 0]", "X[i, -1]", "X[i, 1.5]", "X[i, 1 / 0]",
               "X[i, 0 / 0]", "X[i, 3]", "X[4, 1]")
  for(expression in invalid){
    expect_error(
      BayesTools:::.bt_formula_expression_row_values(
        expressions = list(expression), data = list(X = X), n_rows = 3L,
        context = "Test formula expression"),
      paste0("Test formula expression '", expression,
             "' could not be evaluated: Indexed expression coordinates must be finite positive integers within their array dimensions."),
      fixed = TRUE
    )
  }
  expect_error(
    BayesTools:::.bt_formula_expression_row_values(
      expressions = list("X[i, column[i - 1]]"),
      data = list(X = X, column = c(1L, 2L, 1L)), n_rows = 3L),
    "indices must have length 1 or a common nonzero length.", fixed = TRUE
  )
  expect_error(
    BayesTools:::.bt_formula_expression_row_values(
      expressions = list("X[i, 1]"), data = list(X = 1:3), n_rows = 3L),
    "must supply one index for every array dimension.", fixed = TRUE
  )
  expect_error(
    BayesTools:::.bt_formula_expression_row_values(
      expressions = list("X[i, 1]"),
      data = list(X = array(1:12, c(3L, 2L, 2L))), n_rows = 3L),
    "must supply one index for every array dimension.", fixed = TRUE
  )
  expect_error(
    BayesTools:::.bt_formula_expression_row_values(
      expressions = list("X[i, i[0]]"), data = list(X = X), n_rows = 3L),
    "indices must have length 1 or a common nonzero length.", fixed = TRUE
  )
  expect_error(
    BayesTools:::.bt_formula_expression_eval(
      list(label = "X[i, ]", parsed = quote(X[i, ])), data = list(X = X),
      n_rows = 3L, context = "Test formula expression"),
    "Test formula expression 'X[i, ]' could not be evaluated:", fixed = TRUE
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
  reconstructed <- bayestools_marglik_parameters_formula_design(
    samples = samples,
    design = formula_result$formula_design,
    formula_prior_list = formula_result$prior_list,
    prior_list_parameters = list(),
    log_intercept = FALSE
  )
  expect_equal(unname(reconstructed), c(1 + 2 * 1 + 10, 1 + 2 * 2 + 20))

  bridge_plan <- BayesTools:::.bt_dnode_linear_predictor_fixed_plan(
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

  reconstructed <- bayestools_marglik_parameters_formula_design(
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

  bridge_plan <- BayesTools:::.bt_dnode_linear_predictor_fixed_plan(
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
