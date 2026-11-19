skip_if_not_test_profile("unit")

test_that("formula regularization scopes shared and local scales to their groups", {
  dat <- data.frame(x = 1:6, y = 6:1, f = factor(rep(1:3, 2)))
  global <- prior("cauchy", list(0, 0.4), list(0, Inf))
  intercept <- prior("normal", list(0, 2))
  pool <- prior_regularized("horseshoe", scale = global, slab_scale = 2, slab_df = 4)
  fit <- JAGS_formula(~ x + y + f, "mu", dat, list(intercept = intercept, h = pool))
  priors <- fit$prior_list
  syntax <- JAGS_add_priors("model {}", priors)
  expect_identical(sum(grepl("mu_h_global ~", strsplit(syntax, "\n")[[1]], fixed = TRUE)), 1L)
  expect_identical(sum(grepl("mu_h_slab2 ~", strsplit(syntax, "\n")[[1]], fixed = TRUE)), 1L)
  expect_match(syntax, "mu_h_slab2 ~ dbt_invgamma(2,8)", fixed = TRUE)
  expect_setequal(JAGS_to_monitor(priors), c("mu_intercept", "mu_x", "mu_y", "mu_f",
    "mu_h_global", "mu_h_slab2", "mu_x_local_1", "mu_y_local_1", "mu_f_local_1", "mu_f_local_2"))
  expect_identical(priors$mu_f$regularization$dimension, 2L)
  expect_equal(.bt_regularized_scale(priors$mu_x, 1e300, 1e300, 4), 2)
  expect_identical(priors$mu_intercept$distribution, "normal")

  separate <- list(intercept = intercept, f = prior_factor("normal", list(0, 1), contrast = "treatment"),
    a = prior_regularized("ridge", scale = global, terms = "x"),
    b = prior_regularized("lasso", scale = global, terms = "y"))
  split <- JAGS_formula(~ x + y + f, "mu", dat, separate)$prior_list
  expect_setequal(JAGS_to_monitor(split), c("mu_intercept", "mu_x", "mu_y", "mu_f", "mu_a_global", "mu_b_global"))
  expect_match(JAGS_add_priors("model {}", split), "mu_y ~ ddexp(", fixed = TRUE)

  # A named scalar is an ordinary coefficient prior; the same fixed-scale
  # constructor supplied as a named group fills its unclaimed slopes.
  scalar <- prior_regularized("ridge", scale = 0.7)
  fixed <- JAGS_formula(~ x + y, "mu", dat, list(intercept = intercept, x = scalar, shrink = scalar))
  expect_equal(fixed$prior_list$mu_y$parameters$sd, 0.7)
  expect_equal(JAGS_to_monitor(fixed$prior_list), c("mu_intercept", "mu_x", "mu_y"))

  independent <- prior_regularized("horseshoe", scale = global, contrast = "independent")
  cells <- JAGS_formula(~ 0 + f, "mu", dat, list(shrink = independent))
  expect_identical(cells$prior_list$mu_f$regularization$dimension, 3L)
  fixed_cells <- JAGS_formula(~ 0 + f, "mu", dat,
    list(shrink = prior_regularized("ridge", scale = .7, contrast = "independent")))
  expect_true(is.prior.independent(fixed_cells$prior_list$mu_f))
  direct_cells <- JAGS_formula(~ 0 + f, "mu", dat,
    list(f = prior_regularized("ridge", scale = .7, contrast = "independent")))
  expect_true(is.prior.independent(direct_cells$prior_list$mu_f))
  expect_error(prior_regularized("horseshoe", scale = global, contrast = "meandif"), "contrast")
})

test_that("regularization membership and node ownership survive product-space binding", {
  dat <- data.frame(x = 1:5, y = 5:1, h_global = 1:5)
  intercept <- prior("normal", list(0, 1))
  global <- prior("exp", list(1))
  pool <- prior_regularized("horseshoe", scale = global)
  fit <- JAGS_formula(~ x + y, "mu", dat, list(intercept = intercept, h = pool))
  gated <- lapply(fit$prior_list[c("mu_x", "mu_y")], function(p){
    prior_mixture(list(prior("point", list(0)), p), is_null = c(TRUE, FALSE))
  })
  refit <- JAGS_formula(~ x + y, "mu", dat, list(intercept = intercept, x = gated[[1]], y = gated[[2]]))
  syntax <- JAGS_add_priors("model {}", refit$prior_list)
  expect_match(syntax, "mu_x_component_2 ~ dnorm(", fixed = TRUE)
  expect_match(syntax, "mu_y_component_2 ~ dnorm(", fixed = TRUE)
  expect_setequal(JAGS_to_monitor(refit$prior_list), c("mu_intercept", "mu_x", "mu_y", "mu_x_indicator",
    "mu_y_indicator", "mu_h_global", "mu_x_local_1", "mu_y_local_1"))
  submodel <- JAGS_formula(~ y, "mu", dat, list(intercept = intercept, y = gated[[2]]))
  expect_match(JAGS_add_priors("model {}", submodel$prior_list), "mu_h_global ~ dexp(", fixed = TRUE)

  expect_error(JAGS_formula(~ x + y, "mu", dat, list(intercept = intercept, a = pool, b = pool)), "At most one")
  expect_error(JAGS_formula(~ x, "mu", dat, list(intercept = intercept, x = pool)), "distinct from")
  expect_error(JAGS_formula(~ x, "mu", dat, list(intercept = intercept, x = intercept,
    h = prior_regularized("ridge", terms = "x", scale = global))), "its own prior")
  expect_error(JAGS_formula(~ x + h_global, "mu", dat, list(intercept = intercept, h = pool)), "conflicts")
  expect_error(JAGS_add_priors("model {}", fit$prior_list["mu_y"]), "no owner")

  changed <- refit$prior_list
  changed$mu_y[[2L]]$regularization$scale <- prior("exp", list(10))
  expect_invisible(.bt_regularized_check_fitted_priors(refit$formula_design, refit$prior_list))
  expect_error(.bt_regularized_check_fitted_priors(refit$formula_design, changed),
               class = "BayesTools_refit_required")
})
