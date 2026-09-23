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
        "bridge_context_evaluator"
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

test_that("known group covariance blocks prune fixed-zero columns exactly", {

  data <- data.frame(
    id = factor(c("b", "a", "c", "b", "a", "c")),
    x = c(-1, 0, 1, 2, -0.5, 0.3)
  )
  K <- matrix(
    c(4, 1, 0.5, 1, 9, 2, 0.5, 2, 16),
    nrow = 3L,
    byrow = TRUE,
    dimnames = list(c("a", "b", "c"), c("a", "b", "c"))
  )
  result <- JAGS_formula(
    formula = random_effects_formula(
      ~ diag(1 + x | id),
      group_covariance = random_group_covariance(K, scale = "none")
    ),
    parameter = "mu",
    data = data,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      id = random_block(
        sd = prior("normal", list(0, 1), list(0, Inf)),
        terms = list(x = prior("point", list(location = 0)))
      )
    )
  )
  design <- list(mu = result$formula_design)
  random_term <- result$formula_design$random_effects[[1L]]
  z_names <- .bt_random_effect_latent_names(
    random_term,
    random_term$n_groups,
    random_term$n_columns
  )
  bridge <- .bt_JAGS_formula_random_bridge_parameters(design)
  expect_identical(names(bridge$fixed_latent), as.vector(z_names[, 2L]))

  # The remaining intercept column keeps its exact MVN(0, K) group density;
  # the omitted slope column integrates to one.
  evaluator <- .bt_JAGS_bridge_compile_formula_random_prior_evaluator(
    design,
    omitted_latent = bridge$omitted_latent
  )
  z_intercept <- c(0.4, -1.2, 0.7)
  samples <- stats::setNames(z_intercept, z_names[, 1L])
  precision <- solve(K)
  expected <- -0.5 * (
    3 * log(2 * pi) +
      as.numeric(determinant(K, logarithm = TRUE)$modulus) +
      as.numeric(crossprod(z_intercept, precision %*% z_intercept))
  )
  expect_equal(evaluator$log_prior(samples), expected, tolerance = 1e-12)

  # Partial omission inside a column is still not representable.
  expect_error(
    .bt_JAGS_bridge_compile_formula_random_prior_evaluator(
      design,
      omitted_latent = z_names[1L, 2L]
    ),
    "Bridge sampling can omit only a complete correlated random-effect",
    fixed = TRUE
  )

  # The public bridge reaches the likelihood with the pruned coordinates.
  set.seed(82)
  sd_name <- random_term$sd_parameter_names[[1L]]
  posterior <- cbind(
    mu_intercept = stats::rnorm(40L),
    stats::setNames(
      as.data.frame(matrix(abs(stats::rnorm(40L)) + 0.1, ncol = 1L)),
      sd_name
    ),
    as.data.frame(matrix(
      stats::rnorm(40L * length(z_names)),
      ncol = length(z_names),
      dimnames = list(NULL, as.vector(z_names))
    ))
  )
  fit <- coda::mcmc(as.matrix(posterior))
  attr(fit, "prior_list") <- result$prior_list
  attr(fit, "formula_design") <- design
  bridge_names <- NULL
  testthat::local_mocked_bindings(
    bridge_sampler = function(...){
      arguments <- list(...)
      bridge_names <<- colnames(arguments$samples)
      structure(
        list(logml = -1, niter = 1L, mcse_logml = 0, method = "normal"),
        class = "bridge"
      )
    },
    .package = "bridgesampling"
  )
  bridged <- JAGS_bridgesampling(
    fit = fit,
    log_posterior = function(parameters, data) 0,
    data = list()
  )
  expect_s3_class(bridged, "BayesTools_marglik")
  expect_false(any(as.vector(z_names[, 2L]) %in% bridge_names))
  expect_true(all(as.vector(z_names[, 1L]) %in% bridge_names))
})
