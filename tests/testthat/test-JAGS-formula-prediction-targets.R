skip_if_not_test_profile("unit")

.formula_prediction_data <- function(){
  data.frame(
    x = c(-1, 0, 1, 2),
    id = factor(c("a", "a", "b", "b"), levels = c("a", "b"))
  )
}

.formula_prediction_sd_prior <- function(){
  prior(
    "normal",
    list(mean = 0, sd = 1),
    truncation = list(lower = 0, upper = Inf)
  )
}

.formula_prediction_result <- function(random_effects_compile = NULL,
                                       new_levels = NULL){
  JAGS_formula(
    formula = ~ 1 + x + diag(1 + x | id),
    parameter = "mu",
    data = .formula_prediction_data(),
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    ),
    prior_random = prior_random(
      id = random_block(
        sd = .formula_prediction_sd_prior(),
        monitor = random_monitor(
          latent = FALSE,
          coefficients = TRUE,
          correlation = FALSE
        ),
        new_levels = new_levels
      )
    ),
    random_effects_compile = random_effects_compile
  )
}

.formula_prediction_fit <- function(result){
  posterior <- matrix(
    c(
      10, 1, 2, 3,  0.5, 0.1, -0.5, -0.2,
      20, 2, 2, 3,  0.6, 0.2, -0.6, -0.3
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(NULL, c(
      "mu_intercept",
      "mu_x",
      "mu__xREx__id_intercept",
      "mu__xREx__id_x",
      "mu__xREx__id_xRE_COEFx[1,1]",
      "mu__xREx__id_xRE_COEFx[1,2]",
      "mu__xREx__id_xRE_COEFx[2,1]",
      "mu__xREx__id_xRE_COEFx[2,2]"
    ))
  )
  fit <- coda::mcmc(posterior)
  attr(fit, "formula_design") <- list(mu = result$formula_design)
  fit
}

test_that("conditional target equals fixed target for fixed-only formulas", {

  df <- .formula_prediction_data()
  result <- JAGS_formula(
    formula = ~ 1 + x,
    parameter = "mu",
    data = df,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    )
  )
  posterior <- matrix(
    c(10, 1, 20, 2),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(NULL, c("mu_intercept", "mu_x"))
  )
  fit <- coda::mcmc(posterior)
  attr(fit, "formula_design") <- list(mu = result$formula_design)

  fixed <- JAGS_evaluate_formula(
    fit = fit,
    parameter = "mu",
    formula_target = "fixed"
  )
  conditional <- JAGS_evaluate_formula(
    fit = fit,
    parameter = "mu",
    formula_target = "conditional"
  )

  expect_equal(conditional, fixed)
})

test_that("formula_target fixed and conditional preserve explicit semantics", {

  result <- .formula_prediction_result()
  fit <- .formula_prediction_fit(result)
  df <- .formula_prediction_data()

  expect_error(
    JAGS_evaluate_formula(
      fit = fit,
      formula = ~ 1 + x,
      parameter = "mu",
      data = df,
      prior_list = result$prior_list
    ),
    "silently dropping group-level contributions",
    fixed = TRUE
  )

  fixed <- JAGS_evaluate_formula(
    fit = fit,
    parameter = "mu",
    formula_target = "fixed"
  )
  expect_equal(
    unname(fixed),
    cbind(10 + df$x, 20 + 2 * df$x)
  )

  conditional <- JAGS_evaluate_formula(
    fit = fit,
    parameter = "mu",
    formula_target = "conditional"
  )
  expect_equal(
    unname(conditional),
    cbind(
      c(9.4, 10.5, 10.3, 11.1),
      c(18.4, 20.6, 21.1, 22.8)
    )
  )

  default_prediction <- JAGS_predict_formula(
    fit = fit,
    parameter = "mu"
  )
  expect_equal(default_prediction$value, conditional)
  expect_equal(default_prediction$mean, fixed)
  expect_equal(default_prediction$random, conditional - fixed)

  expect_error(
    JAGS_evaluate_formula(
      fit = fit,
      parameter = "mu",
      formula_target = "fixed",
      new_levels = "zero"
    ),
    "'new_levels' can be used only",
    fixed = TRUE
  )
})

test_that("conditional target handles new levels by explicit policy", {

  result <- .formula_prediction_result()
  fit <- .formula_prediction_fit(result)
  new_data <- data.frame(
    x = c(0, 1, 2),
    id = factor(c("a", "c", "c"), levels = c("a", "b", "c"))
  )

  expect_error(
    JAGS_evaluate_formula(
      fit = fit,
      parameter = "mu",
      data = new_data,
      formula_target = "conditional"
    ),
    "New random-effect level",
    fixed = TRUE
  )

  zero <- JAGS_evaluate_formula(
    fit = fit,
    parameter = "mu",
    data = new_data,
    formula_target = "conditional",
    new_levels = "zero"
  )
  expect_equal(
    unname(zero),
    cbind(c(10.5, 11, 12), c(20.6, 22, 24))
  )

  repeated_new <- data.frame(
    x = c(0, 1, 2),
    id = factor(c("c", "c", "c"), levels = c("a", "b", "c"))
  )
  fixed <- JAGS_evaluate_formula(
    fit = fit,
    parameter = "mu",
    data = repeated_new,
    formula_target = "fixed"
  )
  set.seed(123)
  sampled <- JAGS_evaluate_formula(
    fit = fit,
    parameter = "mu",
    data = repeated_new,
    formula_target = "conditional",
    new_levels = "sample"
  )
  random <- sampled - fixed
  expect_equal(
    unname(random[3, ] - random[2, ]),
    unname(random[2, ] - random[1, ]),
    tolerance = 1e-12
  )
})

test_that("stored new-level policies are used when no override is supplied", {

  result <- .formula_prediction_result(
    new_levels = random_new_levels(method = "zero")
  )
  fit <- .formula_prediction_fit(result)
  new_data <- data.frame(
    x = c(0, 1, 2),
    id = factor(c("a", "c", "c"), levels = c("a", "b", "c"))
  )

  implicit_zero <- JAGS_evaluate_formula(
    fit = fit,
    parameter = "mu",
    data = new_data,
    formula_target = "conditional"
  )
  explicit_zero <- JAGS_evaluate_formula(
    fit = fit,
    parameter = "mu",
    data = new_data,
    formula_target = "conditional",
    new_levels = "zero"
  )
  expect_equal(implicit_zero, explicit_zero)

  marginal_zero <- JAGS_predict_formula(
    fit = fit,
    parameter = "mu",
    data = new_data,
    formula_target = "marginal",
    marginal_method = "covariance"
  )
  expect_true(all(marginal_zero$vcov$samples[, 2:3, ] == 0))
  expect_true(all(marginal_zero$vcov$samples[, , 2:3] == 0))
})

test_that("JAGS_predict_formula composes fixed means with marginal covariance", {

  result <- .formula_prediction_result()
  fit <- .formula_prediction_fit(result)
  new_data <- data.frame(
    x = c(0, 1, 2),
    id = factor(c("a", "c", "c"), levels = c("a", "b", "c"))
  )

  expect_error(
    JAGS_predict_formula(
      fit = fit,
      parameter = "mu",
      formula_target = "fixed",
      blocks = "id"
    ),
    "'blocks' can be used only",
    fixed = TRUE
  )

  expect_error(
    JAGS_predict_formula(
      fit = fit,
      parameter = "mu",
      formula_target = "fixed",
      new_levels = "zero"
    ),
    "'new_levels' can be used only",
    fixed = TRUE
  )

  expect_error(
    JAGS_predict_formula(
      fit = fit,
      parameter = "mu",
      formula_target = "conditional",
      marginal_method = "sample"
    ),
    "'marginal_method' can be used only",
    fixed = TRUE
  )

  expect_error(
    JAGS_predict_formula(
      fit = fit,
      parameter = "mu",
      data = new_data,
      formula_target = "marginal"
    ),
    "explicit new-level policy",
    fixed = TRUE
  )

  prediction <- JAGS_predict_formula(
    fit = fit,
    parameter = "mu",
    data = new_data,
    formula_target = "marginal",
    marginal_method = "covariance",
    new_levels = "sample"
  )
  direct_vcov <- random_effects_marginal_vcov(
    fit = fit,
    parameter = "mu",
    data = new_data,
    posterior_samples = as.matrix(fit),
    prior_list = result$prior_list,
    new_levels = "sample"
  )

  expect_s3_class(prediction, "BayesTools_formula_prediction")
  expect_equal(
    unname(prediction$value),
    unname(JAGS_evaluate_formula(
      fit = fit,
      parameter = "mu",
      data = new_data,
      formula_target = "fixed"
    ))
  )
  expect_equal(prediction$vcov$samples, direct_vcov$samples)
})

test_that("marginal target includes marginalized random-effect blocks", {

  result <- .formula_prediction_result(
    random_effects_compile = random_effects_compile(marginalized = "id")
  )
  random_term <- result$formula_design$random_effects[[1]]
  posterior <- matrix(
    c(10, 1, 2, 3),
    nrow = 1,
    dimnames = list(NULL, c(
      "mu_intercept",
      "mu_x",
      random_term$sd_parameter_names
    ))
  )
  fit <- coda::mcmc(posterior)
  attr(fit, "formula_design") <- list(mu = result$formula_design)

  expect_error(
    JAGS_evaluate_formula(
      fit = fit,
      parameter = "mu",
      formula_target = "conditional"
    ),
    "compiled as marginalized",
    fixed = TRUE
  )

  prediction <- JAGS_predict_formula(
    fit = fit,
    parameter = "mu",
    formula_target = "marginal",
    marginal_method = "covariance"
  )
  expect_equal(prediction$vcov$metadata$included_blocks, "id")
  expect_equal(prediction$vcov$metadata$blocks$id$compile_mode, "marginalized")
})
