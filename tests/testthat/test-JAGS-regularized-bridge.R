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

test_that("centered factor bridge densities use a normalized spherical law for each shrinkage unit", {
  dat <- data.frame(f = factor(rep(letters[1:4], 2)))
  coefficients <- c(.5, -.6, .3)
  cn <- paste0("mu_f[", seq_len(3), "]")
  gamma <- prior("gamma", list(3, 2))
  for(contrast in c("meandif", "orthonormal")){
    for(family in c("ridge", "lasso", "horseshoe")){
      group <- prior_regularized(family, gamma, contrast = contrast,
        slab_scale = if(family == "horseshoe") 1.5)
      model <- JAGS_formula(~ f, "mu", dat, list(intercept = prior("point", list(0)), sh = group))
      priors <- model$prior_list
      spec <- priors$mu_f$regularization
      row <- c(stats::setNames(coefficients, cn), stats::setNames(.7, spec$global_name))
      log_hyper <- stats::dgamma(.7, 3, 2, log = TRUE)
      width <- .7
      if(family == "lasso"){
        expect_length(spec$variance_names, 1L)
        row[spec$variance_names] <- 1.8
        # The radial group-lasso law in dimension three uses shape two.
        log_hyper <- log_hyper + stats::dgamma(1.8, 2, .5, log = TRUE)
        width <- .7 * sqrt(1.8)
      }
      if(family == "horseshoe"){
        expect_length(spec$local_names, 1L)
        row[spec$local_names] <- 1.2
        row[spec$slab_name] <- 1.8
        log_hyper <- log_hyper + log(2) + stats::dcauchy(1.2, log = TRUE) +
          stats::dgamma(1 / 1.8, 2, 2 * 1.5^2, log = TRUE) - 2 * log(1.8)
        width <- 1 / sqrt(1 / (.7 * 1.2)^2 + 1 / 1.8)
      }
      expected <- log_hyper + sum(stats::dnorm(coefficients, 0, width, log = TRUE))
      expect_equal(JAGS_marglik_priors(row, priors), expected)
      expect_equal(JAGS_marglik_priors_rows(rbind(row, row), priors), rep(expected, 2))
      expect_equal(JAGS_marglik_parameters(row, priors)$mu_f, coefficients)
      info <- .JAGS_bridgesampling_posterior_info(priors)
      expect_setequal(as.character(info), names(row))
      private <- .bt_parameter_coordinates_private_auxiliaries(priors)
      expect_setequal(private, setdiff(names(row), cn))
    }
  }
})

test_that("group lasso prior draws have the radial Gamma law and preserve the fitted contrasts", {
  dat <- data.frame(f = factor(rep(letters[1:4], 2)))
  group <- prior_regularized("lasso", 0.7, terms = "f", contrast = "meandif")
  model <- JAGS_formula(~ f, "mu", dat, list(intercept = prior("point", list(0)), sh = group))
  p <- model$prior_list$mu_f
  n <- 5000L
  draws <- .generate_prior_sample_matrix(model$prior_list, n, seed = 183)
  columns <- .JAGS_prior_coefficient_names("mu_f", p)
  expect_length(columns, 3L)
  variance <- draws[, p$regularization$variance_names]
  z <- draws[, columns] / (0.7 * sqrt(variance))
  expect_lt(max(abs(crossprod(z) / n - diag(3))), .06)
  # The joint marginal density is proportional to exp(-||b||/scale).
  # Its standardized radius follows Gamma(d, 1).
  radius <- sqrt(rowSums(draws[, columns]^2)) / .7
  expect_equal(mean(radius < 3), stats::pgamma(3, 3, 1), tolerance = .025)
  set.seed(731)
  raw <- rng(p, 50, transform_factor_samples = FALSE)
  set.seed(731)
  levels <- rng(p, 50)
  expect_equal(unname(levels), unname(raw %*% t(attr(p, "factor_design"))))
  expect_lt(max(abs(rowSums(levels))), 1e-12)
})

test_that("interaction slices have separate shrinkage primitives and gate whole factor vectors", {
  dat <- expand.grid(f = factor(letters[1:3]), g = factor(LETTERS[1:3]))
  priors <- list(intercept = prior("point", list(0)),
    f = prior_factor("mnormal", list(0, 1), contrast = "orthonormal"),
    g = prior_factor("normal", list(0, 1), contrast = "treatment"))
  for(family in c("lasso", "horseshoe")){
    priors$shrink <- prior_regularized(family, .5, terms = "f:g")
    bound <- JAGS_formula(~ f*g, "mu", dat, priors)$prior_list
    name <- "mu_f__xXx__g"
    p <- bound[[name]]
    s <- p$regularization
    expect_identical(s$unit_index, c(1L, 1L, 2L, 2L))
    aux <- if(family == "lasso") s$variance_names else s$local_names
    expect_length(aux, 2)
    # Isolate this fixed-global term. The normalized conditionals of the
    # two contrast slices must read different positive coordinates.
    term <- bound[name]
    coef <- c(-.75, -.3, .25, .8)
    row <- c(stats::setNames(coef, .JAGS_prior_coefficient_names(name, p)),
             stats::setNames(c(.25, 1.5), aux))
    width <- if(family == "lasso") .5 * sqrt(c(.25, 1.5)) else .5 * c(.25, 1.5)
    hp <- if(family == "lasso") stats::dgamma(c(.25, 1.5), 1.5, .5, log = TRUE) else
      log(2) + stats::dcauchy(c(.25, 1.5), log = TRUE)
    expected <- sum(hp) + sum(stats::dnorm(coef, 0, rep(width, each = 2), log = TRUE))
    expect_equal(JAGS_marglik_priors(row, term), expected)
  }

  group <- prior_regularized("lasso", .7, terms = "f", contrast = "meandif")
  base <- JAGS_formula(~ f, "mu", dat, list(intercept = priors$intercept, shrink = group))$prior_list$mu_f
  spike <- prior_factor("point", list(0), contrast = "meandif")
  for(gate in list(prior_spike_and_slab(base), prior_mixture(list(spike, base), is_null = c(TRUE, FALSE)))){
    model <- JAGS_formula(~ f, "mu", dat, list(intercept = priors$intercept, f = gate))
    draws <- .generate_prior_sample_matrix(model$prior_list, 1000L, seed = 825)
    active <- draws[, "mu_f_indicator"] == if(is.prior.spike_and_slab(gate)) 1 else 2
    prefix <- if(is.prior.spike_and_slab(gate)) "mu_f_variable" else "mu_f_component_2"
    expected <- draws[, paste0(prefix, "[", 1:2, "]")] * active
    expect_equal(unname(draws[, paste0("mu_f[", 1:2, "]")]), unname(expected))
    # a mixture component is its own node and names its scales after it
    variance <- if(is.prior.spike_and_slab(gate)) "mu_f_variance_1" else "mu_f_component_2_variance_1"
    expect_equal(sum(colnames(draws) == variance), 1)
  }
})
