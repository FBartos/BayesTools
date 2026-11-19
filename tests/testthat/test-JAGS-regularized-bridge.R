skip_if_not_test_profile("unit")

test_that("a centered shrinkage bridge counts each shared and local scale once", {
  model <- JAGS_formula(~ x + y, "mu", data.frame(x = 1:5, y = c(3, 7, 5, 1, 4)),
    prior_list = list(intercept = prior("point", list(0)),
      shrink = prior_regularized("horseshoe", scale = prior("gamma", list(2, 4)),
        slab_scale = 2, slab_df = 4)))
  priors <- model$prior_list
  x <- priors$mu_x$regularization
  y <- priors$mu_y$regularization
  draws <- cbind(mu_x = c(-0.8, 2), mu_y = c(0.6, -0.1),
    global = c(0.4, 0.9), slab = c(0.5, 2), local_x = c(1.2, 2), local_y = c(0.5, 2.3))
  colnames(draws)[3:6] <- c(x$global_name, x$slab_name, x$local_names, y$local_names)
  s <- function(lambda) 1 / sqrt(1 / (draws[, 3] * lambda)^2 + 1 / draws[, 4])
  expected <- stats::dnorm(draws[, 1], 0, s(draws[, 5]), log = TRUE) +
    stats::dnorm(draws[, 2], 0, s(draws[, 6]), log = TRUE) +
    stats::dgamma(draws[, 3], 2, 4, log = TRUE) +
    stats::dgamma(1 / draws[, 4], 2, 8, log = TRUE) - 2 * log(draws[, 4]) +
    2 * log(2) + stats::dcauchy(draws[, 5], log = TRUE) + stats::dcauchy(draws[, 6], log = TRUE)
  expect_equal(JAGS_marglik_priors_rows(draws, priors), expected)
  expect_equal(JAGS_marglik_priors(draws[1, ], priors), expected[1])
  expect_equal(JAGS_marglik_priors(as.list(draws[1, ]), priors), expected[1])
  bridge <- .bt_JAGS_bridge_compile_prior_list_evaluator(priors)
  expect_equal(bridge$log_prior(draws[2, ]), expected[2])
  expect_equal(bridge$parameters(draws[1, ])$mu_x, unname(draws[1, 1]))
  coordinates <- .JAGS_bridgesampling_posterior_info(priors)
  expect_setequal(as.character(coordinates), colnames(draws))
  expect_equal(unname(attr(coordinates, "lb")[colnames(draws)]), c(-Inf, -Inf, rep(0, 4)))
  expect_error(JAGS_marglik_priors_rows(draws[, -3], priors), class = "BayesTools_missing_monitored_columns")
  draws[1, 3] <- 0
  expect_equal(JAGS_marglik_priors_rows(draws, priors)[1], -Inf)
})

test_that("joint prior draws share the named scale across terms and mixture components", {
  model <- JAGS_formula(~ x + y, "mu", data.frame(x = 1:5, y = c(3, 7, 5, 1, 4)),
    prior_list = list(intercept = prior("point", list(0)),
      shrink = prior_regularized("lasso", scale = prior("lognormal", list(0, 1.5)))))
  priors <- model$prior_list
  global <- priors$mu_x$regularization$global_name
  draws <- .generate_prior_sample_matrix(priors, 3000L, seed = 832)
  # Conditional absolute Laplace draws / global are independent unit exponentials.
  expect_equal(unname(colMeans(abs(draws[, c("mu_x", "mu_y")]) / draws[, global])),
    c(1, 1), tolerance = 0.06)
  expect_gt(stats::cor(log(abs(draws[, "mu_x"])), log(abs(draws[, "mu_y"]))), 0.44)
  expect_equal(sum(colnames(draws) == global), 1L)

  null <- prior("point", list(0))
  priors$mu_x <- prior_mixture(list(null, priors$mu_x), is_null = c(TRUE, FALSE))
  draws <- .generate_prior_sample_matrix(priors, 3000L, seed = 832)
  expect_equal(draws[, "mu_x"], ifelse(draws[, "mu_x_indicator"] == 1, 0, draws[, "mu_x_component_2"]))
  expect_equal(mean(abs(draws[, "mu_x_component_2"]) / draws[, global]), 1, tolerance = 0.06)
  expect_equal(mean(abs(draws[, "mu_y"]) / draws[, global]), 1, tolerance = 0.06)
})

test_that("factor coefficients and their other group members use the same scale", {
  dat <- data.frame(x = 1:9, f = factor(rep(letters[1:3], each = 3)))
  model <- JAGS_formula(~ x + f, "mu", dat,
    prior_list = list(intercept = prior("point", list(0)),
      shrink = prior_regularized("ridge", scale = prior("lognormal", list(0, 1.5)))))
  priors <- model$prior_list
  global <- priors$mu_x$regularization$global_name
  draws <- .generate_prior_sample_matrix(priors, 4000L, seed = 37)
  columns <- c("mu_x", "mu_f[1]", "mu_f[2]")
  expect_equal(unname(colMeans((draws[, columns] / draws[, global])^2)), rep(1, 3), tolerance = 0.06)
  expect_setequal(.bt_parameter_coordinates_private_auxiliaries(priors), global)
  expect_setequal(colnames(.remove_auxiliary_parameters(draws, priors)$model_samples), c("mu_intercept", columns))
  standalone <- rng(priors$mu_f, 1000, transform_factor_samples = FALSE)
  expect_identical(dim(standalone), c(1000L, 2L))
  expect_gt(stats::cor(log(abs(standalone[, 1])), log(abs(standalone[, 2]))), 0.5)
})
