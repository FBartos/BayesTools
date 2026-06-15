skip_if_not_test_profile("unit")

make_random_scale_table_fit <- function(formula_result, posterior){

  fit <- structure(
    list(
      mcmc = coda::mcmc.list(coda::mcmc(posterior)),
      sample = nrow(posterior),
      summary.pars = list(mutate = NULL),
      monitor = colnames(posterior)
    ),
    class = c("runjags", "BayesTools_fit", "list")
  )
  attr(fit, "prior_list") <- formula_result$prior_list
  attr(fit, "formula_design") <- list(mu = formula_result$formula_design)
  attr(fit, "formula_scale") <- list(mu = formula_result$formula_scale)

  fit
}

test_that("JAGS_estimates_table suppresses fixed warnings for random-only scaled slopes", {

  skip_if_not_installed("runjags")

  df <- data.frame(
    x = c(1, 2, 3, 4),
    id = factor(c("a", "a", "b", "b"))
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
        monitor = random_monitor(
          latent = FALSE,
          coefficients = FALSE,
          correlation = FALSE
        )
      )
    )
  )
  posterior <- matrix(
    rep(c(0, 2), 3),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(NULL, c("mu_intercept", "mu__xREx__id_x"))
  )
  fit <- make_random_scale_table_fit(formula_result, posterior)

  samples <- expect_silent(JAGS_estimates_table(
    fit,
    transform_scaled = TRUE,
    random_effects_summary = "standard",
    remove_diagnostics = TRUE,
    return_samples = TRUE
  ))

  expect_equal(
    unname(samples[, "(mu) sd(x | id)"]),
    rep(2 / formula_result$formula_scale$mu_x$sd, nrow(samples)),
    tolerance = 1e-12
  )
})

test_that("JAGS_estimates_table suppresses fixed warnings for homogeneous random-only slopes", {

  skip_if_not_installed("runjags")

  df <- data.frame(
    x = c(1, 2, 3, 4),
    id = factor(c("a", "a", "b", "b"))
  )
  formula_result <- JAGS_formula(
    formula = ~ 1 + id(0 + x | id),
    parameter = "mu",
    data = df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    formula_scale = list(x = TRUE),
    prior_random = prior_random(
      id = random_block(
        sd = prior("gamma", list(2, 2)),
        monitor = random_monitor(
          latent = FALSE,
          coefficients = FALSE,
          correlation = FALSE
        )
      )
    )
  )
  random_term <- formula_result$formula_design$random_effects[[1]]
  posterior <- matrix(
    rep(c(0, 2), 3),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(NULL, c("mu_intercept", random_term$sd_parameter_names))
  )
  fit <- make_random_scale_table_fit(formula_result, posterior)

  samples <- expect_silent(JAGS_estimates_table(
    fit,
    transform_scaled = TRUE,
    random_effects_summary = "standard",
    remove_diagnostics = TRUE,
    return_samples = TRUE
  ))

  sd_column <- grep("^\\(mu\\) sd\\(", colnames(samples), value = TRUE)
  expect_length(sd_column, 1L)
  expect_equal(
    unname(samples[, sd_column]),
    rep(2 / formula_result$formula_scale$mu_x$sd, nrow(samples)),
    tolerance = 1e-12
  )
})

test_that("JAGS_estimates_table treats sd as a valid random-only predictor name", {

  skip_if_not_installed("runjags")

  df <- data.frame(
    sd = c(1, 2, 3, 4),
    id = factor(c("a", "a", "b", "b"))
  )
  formula_result <- JAGS_formula(
    formula = ~ 1 + diag(0 + sd | id),
    parameter = "mu",
    data = df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    formula_scale = list(sd = TRUE),
    prior_random = prior_random(
      id = random_block(
        sd = prior("gamma", list(2, 2)),
        monitor = random_monitor(
          latent = FALSE,
          coefficients = FALSE,
          correlation = FALSE
        )
      )
    )
  )
  random_term <- formula_result$formula_design$random_effects[[1]]
  posterior <- matrix(
    rep(c(0, 2), 3),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(NULL, c("mu_intercept", random_term$sd_parameter_names))
  )
  fit <- make_random_scale_table_fit(formula_result, posterior)

  samples <- expect_silent(JAGS_estimates_table(
    fit,
    transform_scaled = TRUE,
    random_effects_summary = "standard",
    remove_diagnostics = TRUE,
    return_samples = TRUE
  ))

  sd_column <- grep("^\\(mu\\) sd\\(", colnames(samples), value = TRUE)
  expect_length(sd_column, 1L)
  expect_equal(
    unname(samples[, sd_column]),
    rep(2 / formula_result$formula_scale$mu_sd$sd, nrow(samples)),
    tolerance = 1e-12
  )
})

test_that("JAGS_estimates_table still warns about genuinely unused scale entries", {

  skip_if_not_installed("runjags")

  df <- data.frame(
    x = c(1, 2, 3, 4),
    id = factor(c("a", "a", "b", "b"))
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
        monitor = random_monitor(
          latent = FALSE,
          coefficients = FALSE,
          correlation = FALSE
        )
      )
    )
  )
  formula_result$formula_scale$mu_z <- list(mean = 0, sd = 1)
  random_term <- formula_result$formula_design$random_effects[[1]]
  posterior <- matrix(
    rep(c(0, 2), 3),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(NULL, c("mu_intercept", random_term$sd_parameter_names))
  )
  fit <- make_random_scale_table_fit(formula_result, posterior)

  warnings <- character()
  samples <- withCallingHandlers(
    JAGS_estimates_table(
      fit,
      transform_scaled = TRUE,
      random_effects_summary = "standard",
      remove_diagnostics = TRUE,
      return_samples = TRUE
    ),
    warning = function(w){
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  expect_length(warnings, 1L)
  expect_match(warnings, "mu_z", fixed = TRUE)
  expect_false(grepl("mu_x", warnings, fixed = TRUE))
  expect_equal(
    unname(samples[, "(mu) sd(x | id)"]),
    rep(2 / formula_result$formula_scale$mu_x$sd, nrow(samples)),
    tolerance = 1e-12
  )
})

test_that("JAGS_estimates_table unscales diagonal random intercept-slope blocks without correlations", {

  skip_if_not_installed("runjags")

  df <- data.frame(
    x = c(1, 2, 3, 4),
    id = factor(c("a", "a", "b", "b"))
  )
  formula_result <- JAGS_formula(
    formula = ~ 1 + diag(1 + x | id),
    parameter = "mu",
    data = df,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    formula_scale = list(x = TRUE),
    prior_random = prior_random(
      id = random_block(
        sd = prior("gamma", list(2, 2)),
        monitor = random_monitor(
          latent = FALSE,
          coefficients = FALSE,
          correlation = FALSE
        )
      )
    )
  )
  random_term <- formula_result$formula_design$random_effects[[1]]
  posterior <- matrix(
    rep(c(0, 1, 2), 3),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(
      NULL,
      c("mu_intercept", random_term$sd_parameter_names)
    )
  )
  fit <- make_random_scale_table_fit(formula_result, posterior)

  samples <- expect_silent(JAGS_estimates_table(
    fit,
    transform_scaled = TRUE,
    random_effects_summary = "standard",
    remove_diagnostics = TRUE,
    return_samples = TRUE
  ))

  scale_info <- formula_result$formula_scale$mu_x
  expected_intercept_sd <- sqrt(1^2 + (scale_info$mean / scale_info$sd)^2 * 2^2)
  expected_slope_sd <- 2 / scale_info$sd

  expect_equal(
    unname(samples[, "(mu) sd(intercept | id)"]),
    rep(expected_intercept_sd, nrow(samples)),
    tolerance = 1e-12
  )
  expect_equal(
    unname(samples[, "(mu) sd(x | id)"]),
    rep(expected_slope_sd, nrow(samples)),
    tolerance = 1e-12
  )
})

test_that("JAGS_estimates_table keeps fixed and random scaled slope transforms together", {

  skip_if_not_installed("runjags")

  df <- data.frame(
    x = c(1, 2, 3, 4),
    id = factor(c("a", "a", "b", "b"))
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
        monitor = random_monitor(
          latent = FALSE,
          coefficients = FALSE,
          correlation = FALSE
        )
      )
    )
  )
  random_term <- formula_result$formula_design$random_effects[[1]]
  posterior <- matrix(
    rep(c(10, 4, 2), 3),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(
      NULL,
      c("mu_intercept", "mu_x", random_term$sd_parameter_names)
    )
  )
  fit <- make_random_scale_table_fit(formula_result, posterior)

  samples <- expect_silent(JAGS_estimates_table(
    fit,
    transform_scaled = TRUE,
    random_effects_summary = "standard",
    remove_diagnostics = TRUE,
    return_samples = TRUE
  ))

  scale_info <- formula_result$formula_scale$mu_x
  expect_equal(
    unname(samples[, "(mu) intercept"]),
    rep(10 - scale_info$mean / scale_info$sd * 4, nrow(samples)),
    tolerance = 1e-12
  )
  expect_equal(
    unname(samples[, "(mu) x"]),
    rep(4 / scale_info$sd, nrow(samples)),
    tolerance = 1e-12
  )
  expect_equal(
    unname(samples[, "(mu) sd(x | id)"]),
    rep(2 / scale_info$sd, nrow(samples)),
    tolerance = 1e-12
  )
})
