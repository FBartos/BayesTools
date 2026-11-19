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
  expect_null(priors$mu_f$regularization$unit_index)
  expect_null(priors$mu_f$regularization$variance_names)
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
  expect_error(prior_regularized("horseshoe", scale = global, contrast = "ordered"), "contrast")
})

test_that("centered factor kernels retain their coded scale and one scale per unit", {
  dat <- data.frame(f = factor(letters[1:4]))
  for(contrast in c("meandif", "orthonormal")){
    bind <- function(p) JAGS_formula(~ 0 + f, "mu", dat, list(sh = p))$prior_list$mu_f
    ridge <- bind(prior_regularized("ridge", .7, contrast = contrast))
    normal <- JAGS_formula(~ 0 + f, "mu", dat, list(f = prior_factor("mnormal", list(0, .7),
      contrast = contrast)))$prior_list$mu_f
    expect_identical(ridge, normal)
    lasso <- bind(prior_regularized("lasso", .7, contrast = contrast))
    expect_identical(lasso$regularization$unit_index, rep(1L, 3L))
    expect_identical(lasso$regularization$unit_dimension, 3L)
    expect_identical(lasso$regularization$unit_centered, TRUE)
    expect_identical(.bt_regularized_primitives(lasso)$mu_f_variance_1$parameters, list(shape = 2, rate = .5))
    syntax <- JAGS_add_priors("model {}", list(mu_f = lasso))
    expect_match(syntax, "mu_f_variance_1 ~ dgamma(2,0.5)", fixed = TRUE)
    for(j in 1:3){
      expect_match(syntax, paste0("mu_f[", j, "] ~ dnorm(0, exp(-2 * (log(", .prior_ordered_format_number(.7),
        ") + 0.5 * log(mu_f_variance_1)))"), fixed = TRUE)
    }
    hs <- bind(prior_regularized("horseshoe", .7, contrast = contrast))
    expect_identical(hs$regularization$local_names, "mu_f_local_1")
    expect_setequal(JAGS_to_monitor(list(mu_f = hs)), c("mu_f", "mu_f_local_1"))

    # A fixed shorthand placed directly on a factor has the same grouped law.
    direct <- JAGS_formula(~ 0 + f, "mu", dat,
      list(f = prior_regularized("lasso", .7, contrast = contrast)))$prior_list$mu_f
    expect_equal(.bt_regularized_jags(direct, "mu_f"), .bt_regularized_jags(lasso, "mu_f"))
    inherited <- prior_regularized("lasso", .7)
    attr(inherited, "regularized_contrast") <- contrast
    direct <- JAGS_formula(~ 0 + f, "mu", dat, list(f = inherited))$prior_list$mu_f
    expect_equal(.bt_regularized_jags(direct, "mu_f"), .bt_regularized_jags(lasso, "mu_f"))
  }
})

test_that("centered interactions share only within the fitted partner column", {
  dat <- expand.grid(f = factor(letters[1:3]), h = factor(LETTERS[1:4]))
  global <- prior("uniform", list(.1, 1))
  for(contrast in c("treatment", "independent", "meandif")){
    partners <- switch(contrast, treatment = 3L, independent = 4L, meandif = 1L)
    dimension <- if(contrast == "meandif") 6L else 2L
    h_prior <- if(contrast == "meandif") prior_factor("mnormal", list(0, 1), contrast = contrast) else
      prior_factor("normal", list(0, 1), contrast = contrast)
    for(formula in list(~ f * h, ~ h * f)){
      term <- tail(attr(terms(formula), "term.labels"), 1)
      term_name <- paste0("mu_", gsub(":", "__xXx__", term, fixed = TRUE))
      priors <- JAGS_formula(formula, "mu", dat, list(
        intercept = prior("point", list(0)), f = prior_factor("mnormal", list(0, 1), contrast = "orthonormal"),
        h = h_prior, pool = prior_regularized("lasso", global, terms = "f:h", contrast = "meandif")))$prior_list
      p <- priors[[term_name]]
      expected <- if(term == "f:h") rep(seq_len(partners), each = dimension) else rep(seq_len(partners), dimension)
      expect_identical(p$regularization$unit_index, expected)
      expect_identical(p$regularization$unit_dimension, rep(dimension, partners))
      variance <- .bt_regularized_primitives(p)[p$regularization$variance_names]
      expect_true(all(vapply(variance, function(v) identical(v$parameters, list(shape = (dimension + 1) / 2,
                                                                                 rate = .5)), logical(1))))
      # RoBMA rebinds already-assigned alternatives inside a term mixture.
      names(priors) <- sub("^mu_", "", names(priors))
      names(priors) <- gsub("__xXx__", ":", names(priors), fixed = TRUE)
      priors[[term]] <- prior_mixture(list(prior("point", list(0)), p), is_null = c(TRUE, FALSE))
      refit <- JAGS_formula(formula, "mu", dat, priors)$prior_list[[term_name]]
      expect_identical(refit[[2]]$regularization, p$regularization)
    }
  }
  centered <- prior_regularized("horseshoe", global, contrast = "orthonormal")
  subset <- JAGS_formula(~ f + f:h, "mu", dat,
    list(intercept = prior("point", list(0)), sh = centered))$prior_list$mu_f__xXx__h
  expect_identical(unname(attr(subset, "factor_contrasts")), c("contr.independent", "contr.orthonormal"))
  expect_identical(subset$regularization$unit_index, rep(1:3, 3))
  expect_length(subset$regularization$local_names, 3L)
  expect_error(JAGS_formula(~ f:h, "mu", dat, list(intercept = prior("point", list(0)), sh = centered)),
    "codes 'f' and 'h' by level indicators")
  ordinary <- JAGS_formula(~ f * h, "mu", dat, list(intercept = prior("point", list(0)),
    sh = prior_regularized("lasso", .7, contrast = "treatment")))$prior_list$mu_f__xXx__h
  expect_identical(ordinary$distribution, "laplace")
  expect_false(inherits(ordinary, "prior.regularized_term"))
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
