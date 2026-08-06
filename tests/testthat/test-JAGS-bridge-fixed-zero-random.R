skip_if_not_test_profile("unit")

.fixed_zero_bridge_fit <- function(formula_result, values){

  latent_names <- unlist(lapply(
    formula_result$formula_design$random_effects,
    function(random_term){
      as.vector(.bt_random_effect_latent_names(
        random_term = random_term,
        n_groups = random_term$n_groups,
        n_columns = random_term$n_columns
      ))
    }
  ))
  posterior <- matrix(
    rep(values[latent_names], each = 20L),
    nrow = 20L,
    dimnames = list(NULL, latent_names)
  )
  fit <- coda::mcmc(posterior)
  attr(fit, "prior_list") <- formula_result$prior_list
  attr(fit, "formula_design") <- list(mu = formula_result$formula_design)
  fit
}

test_that("fixed-zero random blocks leave the bridge target and dimension invariant", {

  data <- data.frame(
    id_zero = factor(c("a", "b", "c", "d")),
    id_active = factor(c("u", "u", "v", "v"))
  )
  fixed_prior <- list(intercept = prior("point", list(location = 0)))
  with_zero <- JAGS_formula(
    formula = ~ 1 + diag(1 | id_zero) + diag(1 | id_active),
    parameter = "mu",
    data = data,
    prior_list = fixed_prior,
    prior_random = prior_random(
      id_zero = random_block(sd = prior("point", list(location = 0))),
      id_active = random_block(sd = prior("point", list(location = 0.4)))
    )
  )
  active_only <- JAGS_formula(
    formula = ~ 1 + diag(1 | id_active),
    parameter = "mu",
    data = data,
    prior_list = fixed_prior,
    prior_random = prior_random(
      id_active = random_block(sd = prior("point", list(location = 0.4)))
    )
  )
  only_zero <- JAGS_formula(
    formula = ~ 1 + diag(1 | id_zero),
    parameter = "mu",
    data = data,
    prior_list = fixed_prior,
    prior_random = prior_random(
      id_zero = random_block(sd = prior("point", list(location = 0)))
    )
  )

  zero_term <- with_zero$formula_design$random_effects[[1L]]
  active_term <- with_zero$formula_design$random_effects[[2L]]
  zero_names <- as.vector(.bt_random_effect_latent_names(
    zero_term,
    zero_term$n_groups,
    zero_term$n_columns
  ))
  active_names <- as.vector(.bt_random_effect_latent_names(
    active_term,
    active_term$n_groups,
    active_term$n_columns
  ))
  values <- stats::setNames(
    c(rep(3, length(zero_names)), -0.25, 0.5),
    c(zero_names, active_names)
  )
  zero_fit <- .fixed_zero_bridge_fit(with_zero, values)
  active_fit <- .fixed_zero_bridge_fit(active_only, values)
  only_zero_fit <- .fixed_zero_bridge_fit(only_zero, values)

  bridge_dimensions <- integer()
  bridge_names <- list()
  testthat::local_mocked_bindings(
    bridge_sampler = function(...){
      arguments <- list(...)
      bridge_dimensions <<- c(bridge_dimensions, ncol(arguments$samples))
      bridge_names[[length(bridge_names) + 1L]] <<- colnames(arguments$samples)
      callback_names <- c(
        "data",
        "bridge_prior_evaluator",
        "bridge_formula_prior_evaluator",
        "bridge_formula_random_prior_evaluator",
        "bridge_formula_parameter_evaluator",
        "add_parameters",
        "fixed_random_latent",
        "bridge_context",
        "formula_design_list",
        "formula_data_list",
        "formula_prior_list"
      )
      logml <- do.call(
        arguments$log_posterior,
        c(
          list(samples.row = arguments$samples[1L, ]),
          arguments[callback_names]
        )
      )
      structure(
        list(logml = logml, niter = 1L, mcse_logml = 0, method = "normal"),
        class = "bridge"
      )
    },
    .package = "bridgesampling"
  )
  log_posterior <- function(parameters, data){
    sum(stats::dnorm(data$y, mean = parameters$mu, sd = 1, log = TRUE))
  }
  zero_result <- JAGS_bridgesampling(
    fit = zero_fit,
    log_posterior = log_posterior,
    data = list(y = c(-0.2, 0.1, 0.3, -0.1))
  )
  active_result <- JAGS_bridgesampling(
    fit = active_fit,
    log_posterior = log_posterior,
    data = list(y = c(-0.2, 0.1, 0.3, -0.1))
  )
  only_zero_result <- JAGS_bridgesampling(
    fit = only_zero_fit,
    log_posterior = log_posterior,
    data = list(y = c(-0.2, 0.1, 0.3, -0.1))
  )

  expect_equal(zero_result$logml, active_result$logml, tolerance = 1e-12)
  expect_identical(bridge_dimensions, c(2L, 2L))
  expect_identical(bridge_names[[1L]], active_names)
  expect_false(any(zero_names %in% bridge_names[[1L]]))
  expect_identical(
    only_zero_result$aggregation$rule,
    "exact_zero_dimensional"
  )
  expect_equal(
    only_zero_result$logml,
    sum(stats::dnorm(c(-0.2, 0.1, 0.3, -0.1), log = TRUE)),
    tolerance = 1e-12
  )
})

test_that("independent fixed-zero components prune only their latent columns", {

  data <- data.frame(
    x = c(0, 1, 0, 1),
    id = factor(c("a", "a", "b", "b"))
  )
  result <- JAGS_formula(
    formula = ~ 1 + diag(1 + x | id),
    parameter = "mu",
    data = data,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(
        sd = prior("gamma", list(2, 2)),
        terms = list(intercept = prior("point", list(location = 0)))
      )
    )
  )
  random_term <- result$formula_design$random_effects[[1L]]
  z_names <- .bt_random_effect_latent_names(
    random_term,
    random_term$n_groups,
    random_term$n_columns
  )
  bridge <- .bt_JAGS_formula_random_bridge_parameters(
    list(mu = result$formula_design)
  )

  expect_identical(names(bridge$fixed_latent), as.vector(z_names[, 1L]))
  expect_identical(bridge$parameters, as.vector(z_names[, 2L]))
})
