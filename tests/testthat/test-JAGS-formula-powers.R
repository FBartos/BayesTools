skip_if_not_test_profile("unit")

test_that("R133 literal whole formula powers expand fixed interactions and replay", {
  data <- data.frame(x = c(-1, 0, 2, 3), z = c(2, 1, -1, 4), y = 0)
  priors <- setNames(rep(list(prior("normal", list(0, 1))), 4), c("intercept", "x", "z", "x:z"))
  explicit <- JAGS_formula(~ x + z + x:z, "mu", data, priors)
  reference <- cbind(1, data$x, data$z, data$x * data$z)
  for(formula in list(~ (x + z)^2, ~ (x + z)^3, y ~ (x + z)^2, ~ x*z)){
    result <- JAGS_formula(formula, "mu", data, priors)
    expect_equal(result$formula_design$model_matrix, explicit$formula_design$model_matrix)
    expect_equal(result$formula_design$model_matrix, reference, ignore_attr = TRUE)
    expect_identical(result$formula_design$model_terms, explicit$formula_design$model_terms)
    expect_identical(result$formula_design$name_map, explicit$formula_design$name_map)
    posterior <- coda::mcmc(matrix(c(1, 2, 3, 4), 1, dimnames = list(NULL, names(result$prior_list))))
    attr(posterior, "formula_design") <- list(mu = result$formula_design)
    replay <- JAGS_evaluate_formula(posterior, parameter = "mu", data = data, prior_list = result$prior_list)
    expect_equal(as.numeric(replay), as.numeric(reference %*% c(1, 2, 3, 4)))
  }
})

test_that("R133 random whole powers retain the explicit design and structural covariance", {
  data <- data.frame(x = c(-1, 0, 2, 3), z = c(2, 1, -1, 4), g = factor(c("a", "a", "b", "b")))
  compile <- function(formula) JAGS_formula(formula, "mu", data,
    list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(sd = prior("point", list(1))))
  explicit <- compile(~ diag(x + z + x:z | g))
  for(formula in list(~ diag((x + z)^2 | g), ~ diag((x + z)^3 | g))){
    result <- compile(formula)
    block <- result$formula_design$random_effects[[1L]]
    expect_equal(block$model_matrix, explicit$formula_design$random_effects[[1L]]$model_matrix)
    expect_equal(block$model_matrix, cbind(1, data$x, data$z, data$x * data$z), ignore_attr = TRUE)
    expect_identical(block$column_names, explicit$formula_design$random_effects[[1L]]$column_names)
    prediction <- .bt_random_effect_prediction_data(block, data)
    expect_equal(prediction$model_matrix, block$model_matrix)
    draws <- matrix(numeric(), nrow = 1, ncol = 0)
    covariance <- random_effects_marginal_vcov(result$formula_design, posterior_samples = draws,
      prior_list = result$prior_list)
    reference <- tcrossprod(block$model_matrix) * outer(data$g, data$g, `==`)
    expect_equal(as.numeric(covariance$samples[1, , ]), as.numeric(reference))
  }
})

test_that("R133 unsupported public powers and arithmetic calls retain refusal", {
  data <- data.frame(x = 1:4, z = 4:1)
  priors <- setNames(rep(list(prior("normal", list(0, 1))), 4), c("intercept", "x", "z", "x:z"))
  for(formula in list(~ (x + z)^0, ~ (x + z)^1, ~ (x + z)^(-1), ~ (x + z)^degree)){
    expect_error(JAGS_formula(formula, "mu", data, priors), "invalid power")
  }
  expect_warning(expect_error(JAGS_formula(~ (x + z)^Inf, "mu", data, priors), "invalid power"), "NAs introduced by coercion")
  expect_error(JAGS_formula(~ (x + z)^2.5, "mu", data, priors), "Unsupported fixed-formula expression")
  expect_error(JAGS_formula(~ I(x^2), "mu", data, priors), "Unsupported fixed-formula call")
  expect_true(.bt_validate_fixed_formula_grammar(~ (x + z)^0))
  expect_true(.bt_validate_fixed_formula_grammar(~ (x + z)^1))
})
