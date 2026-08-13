skip_if_not_test_profile("unit")

.bridge_marginal_random_fit <- function(formula_result, values, n_draws = 20L){

  posterior <- matrix(
    rep(values, each = n_draws),
    nrow = n_draws,
    dimnames = list(NULL, names(values))
  )
  fit <- coda::mcmc(posterior)
  attr(fit, "prior_list") <- formula_result$prior_list
  attr(fit, "formula_design") <- list(mu = formula_result$formula_design)
  fit
}

.bridge_marginal_random_mvn <- function(y, mean, covariance){

  chol_covariance <- chol(covariance)
  residual <- y - mean
  z <- backsolve(chol_covariance, residual, transpose = TRUE)
  -0.5 * (
    length(y) * log(2 * pi) +
      2 * sum(log(diag(chol_covariance))) +
      sum(z^2)
  )
}

.bridge_marginal_random_dense <- function(x){

  if(identical(x$representation, "factor_state")){
    x$factors <- Map(c, x$factor_plans, x$factor_states)
    x$representation <- "factor"
  }
  if(identical(x$representation, "dense")){
    return(x$covariance)
  }
  covariance <- matrix(0, nrow = x$dimension, ncol = x$dimension)
  for(index in x$row_blocks){
    block_covariance <- matrix(0, nrow = length(index), ncol = length(index))
    for(factor in x$factors){
      block_covariance <- block_covariance +
        .bt_JAGS_bridge_marginal_random_geometry_covariance(factor, index)
    }
    covariance[index, index] <- block_covariance
  }
  covariance
}

test_that("bridge-only random marginalization preserves the exact Gaussian target", {

  data <- data.frame(study = factor(c("a", "a", "b", "c")))
  formula_result <- JAGS_formula(
    formula = ~ 1 + diag(1 | study),
    parameter = "mu",
    data = data,
    prior_list = list(intercept = prior("point", list(location = 0))),
    prior_random = prior_random(
      study = random_block(sd = prior("point", list(location = 0.4)))
    )
  )
  random_term <- formula_result$formula_design$random_effects[[1L]]
  z_names <- as.vector(.bt_random_effect_latent_names(
    random_term,
    random_term$n_groups,
    random_term$n_columns
  ))
  fit <- .bridge_marginal_random_fit(
    formula_result,
    stats::setNames(c(-0.2, 0.3, 0.1), z_names)
  )
  y <- c(-0.1, 0.2, 0.4, -0.3)
  sampling_covariance <- diag(c(0.2, 0.3, 0.4, 0.5))
  same_study <- outer(data$study, data$study, "==")
  expected_random_covariance <- 0.4^2 * same_study
  seen <- new.env(parent = emptyenv())

  result <- JAGS_bridgesampling(
    fit = fit,
    data = list(y = y, sampling_covariance = sampling_covariance),
    bridge_context = "marginal",
    formula_random_effects_marginalize_list = list(
      mu = list(
        blocks = "study",
        row_blocks = list(1:2, 3L, 4L)
      )
    ),
    log_posterior = function(parameters, data, bridge_context){
      seen$parameters <- parameters
      seen$context <- bridge_context
      random_covariance <- .bridge_marginal_random_dense(
        bridge_context$marginalized_random$mu
      )
      .bridge_marginal_random_mvn(
        y = data$y,
        mean = parameters$mu,
        covariance = data$sampling_covariance +
          random_covariance
      )
    }
  )

  expect_identical(result$aggregation$rule, "exact_zero_dimensional")
  expect_equal(seen$parameters$mu, rep(0, nrow(data)), tolerance = 0)
  expect_s3_class(seen$context, "BayesTools_bridge_marginal_context")
  expect_identical(seen$context$marginalized_random$mu$blocks, "study")
  expect_identical(
    seen$context$marginalized_random$mu$representation,
    "factor"
  )
  expect_equal(
    unname(.bridge_marginal_random_dense(
      seen$context$marginalized_random$mu
    )),
    unname(expected_random_covariance),
    tolerance = 1e-12
  )
  expect_equal(
    result$logml,
    .bridge_marginal_random_mvn(
      y = y,
      mean = rep(0, length(y)),
      covariance = sampling_covariance + expected_random_covariance
    ),
    tolerance = 1e-12
  )
})

test_that("bridge-only random marginalization retains covariance parameters", {

  data <- data.frame(study = factor(c("a", "a", "b", "c")))
  formula_result <- JAGS_formula(
    formula = ~ 1 + diag(1 | study),
    parameter = "mu",
    data = data,
    prior_list = list(intercept = prior("point", list(location = 0))),
    prior_random = prior_random(
      study = random_block(sd = prior("gamma", list(2, 2)))
    )
  )
  random_term <- formula_result$formula_design$random_effects[[1L]]
  z_names <- as.vector(.bt_random_effect_latent_names(
    random_term,
    random_term$n_groups,
    random_term$n_columns
  ))
  sd_name <- random_term$sd_parameter_names[[1L]]
  values <- stats::setNames(c(0.6, -0.2, 0.3, 0.1), c(sd_name, z_names))
  fit <- .bridge_marginal_random_fit(formula_result, values)
  bridge_names <- NULL
  seen_covariance <- NULL

  testthat::local_mocked_bindings(
    bridge_sampler = function(...){
      arguments <- list(...)
      bridge_names <<- colnames(arguments$samples)
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

  JAGS_bridgesampling(
    fit = fit,
    data = list(),
    bridge_context = "marginal",
    formula_random_effects_marginalize_list = list(mu = "study"),
    log_posterior = function(parameters, data, bridge_context){
      seen_covariance <<- bridge_context$marginalized_random$mu$covariance
      0
    }
  )

  expect_identical(bridge_names, sd_name)
  expect_false(any(z_names %in% bridge_names))
  expect_equal(
    unname(seen_covariance),
    unname(values[[sd_name]]^2 * outer(data$study, data$study, "==")),
    tolerance = 1e-12
  )
})

test_that("marginalization removes exactly the selected latent coordinates", {

  data <- data.frame(
    study = factor(c("a", "a", "b", "c")),
    x = c(-1, 0, 1, 2)
  )
  formula_result <- JAGS_formula(
    formula = ~ 1 + us(1 + x | study),
    parameter = "mu",
    data = data,
    prior_list = list(intercept = prior("point", list(location = 0))),
    prior_random = prior_random(
      study = random_block(
        sd = prior("gamma", list(2, 2)),
        terms = list(x = prior("gamma", list(2, 2))),
        cor = prior_lkj(eta = 2, include_primitives = TRUE)
      )
    )
  )
  design <- list(mu = formula_result$formula_design)
  random_term <- formula_result$formula_design$random_effects[[1L]]
  latent_names <- as.vector(.bt_random_effect_latent_names(
    random_term,
    random_term$n_groups,
    random_term$n_columns
  ))
  joint <- .bt_JAGS_formula_random_bridge_parameters(design)
  marginalized <- .bt_JAGS_formula_random_bridge_parameters(
    formula_design_list = design,
    marginal_random_spec = list(
      mu = list(blocks = "study", row_blocks = NULL)
    )
  )

  expect_setequal(
    setdiff(joint$parameters, marginalized$parameters),
    latent_names
  )
  expect_length(setdiff(marginalized$parameters, joint$parameters), 0L)
  expect_true(all(random_term$correlation$primitive_names %in%
                    marginalized$parameters))
  expect_setequal(marginalized$omitted_latent, latent_names)
})

test_that("compiled bridge SD extraction safely reuses posterior positions", {

  data <- data.frame(
    id = factor(rep(c("a", "b"), each = 2L)),
    x = rep(c(0, 1), 2L)
  )
  formula_result <- JAGS_formula(
    formula = ~ 1 + diag(1 + x | id),
    parameter = "mu",
    data = data,
    prior_list = list(intercept = prior("point", list(location = 0))),
    prior_random = prior_random(
      id = random_block(
        sd = prior("gamma", list(2, 2)),
        terms = list(x = prior("gamma", list(3, 2)))
      )
    )
  )
  random_term <- formula_result$formula_design$random_effects[[1L]]
  evaluator <- .bt_JAGS_bridge_compile_random_sd_evaluator(
    random_term = random_term,
    prior_list = formula_result$prior_list
  )
  sd_names <- random_term$sd_parameter_names
  first <- matrix(
    c(9, 0.3, 0.5),
    nrow = 1L,
    dimnames = list(NULL, c("unrelated", sd_names))
  )
  reordered <- matrix(
    c(0.6, 8, 0.4),
    nrow = 1L,
    dimnames = list(NULL, c(sd_names[[2L]], "unrelated", sd_names[[1L]]))
  )
  override <- stats::setNames(list(0.7, 0.8), sd_names)

  expect_equal(
    evaluator$posterior_values(first),
    c(0.3, 0.5),
    tolerance = 0
  )
  expect_equal(
    evaluator$posterior_values(reordered),
    c(0.4, 0.6),
    tolerance = 0
  )
  expect_equal(
    evaluator$posterior_values(reordered, parameters = override),
    c(0.7, 0.8),
    tolerance = 0
  )
})

test_that("bridge-only random marginalization requires a covariance context", {

  data <- data.frame(study = factor(c("a", "b")))
  formula_result <- JAGS_formula(
    formula = ~ 1 + diag(1 | study),
    parameter = "mu",
    data = data,
    prior_list = list(intercept = prior("point", list(location = 0))),
    prior_random = prior_random(
      study = random_block(sd = prior("point", list(location = 0.4)))
    )
  )
  random_term <- formula_result$formula_design$random_effects[[1L]]
  z_names <- as.vector(.bt_random_effect_latent_names(
    random_term,
    random_term$n_groups,
    random_term$n_columns
  ))
  fit <- .bridge_marginal_random_fit(
    formula_result,
    stats::setNames(c(0.1, -0.1), z_names)
  )

  expect_error(
    JAGS_bridgesampling(
      fit = fit,
      data = list(),
      formula_random_effects_marginalize_list = list(mu = "study"),
      log_posterior = function(parameters, data) 0
    ),
    "requires bridge_context"
  )

  expect_error(
    JAGS_bridgesampling(
      fit = fit,
      data = list(),
      bridge_context = "marginal",
      formula_random_effects_marginalize_list = list(
        mu = list(
          blocks = "study",
          factor_state = TRUE
        )
      ),
      log_posterior = function(parameters, data, bridge_context) 0
    ),
    "factor_state.*requires exact 'row_blocks'"
  )
})

test_that("block covariance contract rejects separated random dependencies", {

  data <- data.frame(study = factor(c("a", "a", "b", "c")))
  formula_result <- JAGS_formula(
    formula = ~ 1 + diag(1 | study),
    parameter = "mu",
    data = data,
    prior_list = list(intercept = prior("point", list(location = 0))),
    prior_random = prior_random(
      study = random_block(sd = prior("point", list(location = 0.4)))
    )
  )
  random_term <- formula_result$formula_design$random_effects[[1L]]
  z_names <- as.vector(.bt_random_effect_latent_names(
    random_term,
    random_term$n_groups,
    random_term$n_columns
  ))
  fit <- .bridge_marginal_random_fit(
    formula_result,
    stats::setNames(c(0.1, -0.1, 0.2), z_names)
  )

  expect_error(
    JAGS_bridgesampling(
      fit = fit,
      data = list(),
      bridge_context = "marginal",
      formula_random_effects_marginalize_list = list(
        mu = list(
          blocks = "study",
          row_blocks = list(1L, 2L, 3:4)
        )
      ),
      log_posterior = function(parameters, data, bridge_context) 0
    ),
    "structurally nonzero covariance"
  )
})

test_that("bridge marginal evaluator supports every implemented covariance structure", {

  data <- data.frame(
    id = factor(rep(c("a", "b"), each = 3L)),
    x = rep(c(0, 1, 2), 2L),
    f = factor(rep(c("a", "b", "c"), 2L)),
    time = rep(c(0, 2, 5), 2L)
  )
  cases <- list(
    diag = list(
      formula = ~ 1 + diag(1 + x | id),
      prior_random = prior_random(
        id = random_block(
          sd = prior("point", list(location = 0.3)),
          terms = list(x = prior("point", list(location = 0.5)))
        )
      ),
      values = numeric()
    ),
    us = list(
      formula = ~ 1 + us(1 + x | id),
      prior_random = prior_random(
        id = random_block(
          sd = prior("point", list(location = 0.3)),
          terms = list(x = prior("point", list(location = 0.5))),
          cor = prior_lkj(eta = 2, include_primitives = TRUE)
        )
      ),
      values = NULL
    ),
    cs = list(
      formula = ~ 1 + cs(f | id),
      prior_random = prior_random(
        id = random_block(
          sd = prior("point", list(location = 0.4)),
          rho = prior("point", list(location = 0.2))
        )
      ),
      values = numeric()
    ),
    hcs = list(
      formula = ~ 1 + hcs(f | id),
      prior_random = prior_random(
        id = random_block(
          sd = prior("point", list(location = 0.3)),
          rho = prior("point", list(location = 0.2))
        )
      ),
      values = numeric()
    ),
    ar1 = list(
      formula = ~ 1 + ar1(f | id),
      prior_random = prior_random(
        id = random_block(
          sd = prior("point", list(location = 0.4)),
          rho = prior("point", list(location = 0.5))
        )
      ),
      values = numeric()
    ),
    har = list(
      formula = ~ 1 + har(f | id),
      prior_random = prior_random(
        id = random_block(
          sd = prior("point", list(location = 0.3)),
          rho = prior("point", list(location = 0.5))
        )
      ),
      values = numeric()
    ),
    car = list(
      formula = ~ 1 + car(time | id),
      prior_random = prior_random(
        id = random_block(
          sd = prior("point", list(location = 0.4)),
          rho = prior("point", list(location = 0.5))
        )
      ),
      values = numeric()
    )
  )

  for(structure in names(cases)){
    case <- cases[[structure]]
    formula_result <- JAGS_formula(
      formula = case$formula,
      parameter = "mu",
      data = data,
      prior_list = list(intercept = prior("point", list(location = 0))),
      prior_random = case$prior_random
    )
    random_term <- formula_result$formula_design$random_effects[[1L]]
    values <- case$values
    if(identical(structure, "us")){
      primitive_names <- random_term$correlation$primitive_names
      values <- stats::setNames(rep(0.5, length(primitive_names)), primitive_names)
    }
    posterior <- matrix(
      values,
      nrow = 1L,
      dimnames = list(NULL, names(values))
    )
    posterior <- .bt_random_effect_marginal_covariance_validate_posterior(
      posterior,
      allow_zero_columns = TRUE
    )
    reference <- .bt_random_effect_marginal_covariance_samples(
      design = formula_result$formula_design,
      posterior = posterior,
      prior_list = formula_result$prior_list,
      blocks = "id"
    )$samples[1L, , ]
    evaluator <- .bt_JAGS_bridge_compile_marginal_random_evaluator(
      formula_design_list = list(mu = formula_result$formula_design),
      marginal_random_spec = list(
        mu = list(blocks = "id", row_blocks = list(1:3, 4:6))
      ),
      formula_data_list = list(mu = data),
      formula_prior_list = list(mu = formula_result$prior_list),
      model_data = list()
    )
    actual_value <- evaluator$covariance(
      samples = values,
      prior_parameters = list(),
      formula_prior_parameters = list(mu = list()),
      formula_parameters = list(mu = rep(0, nrow(data)))
    )$mu
    actual <- .bridge_marginal_random_dense(actual_value)

    expect_identical(actual_value$representation, "factor", info = structure)
    expect_equal(
      unname(actual),
      unname(reference),
      tolerance = 1e-12,
      info = structure
    )
  }
})

test_that("bridge coefficient geometry reconstructs a structured factor once", {

  data <- data.frame(
    id = factor(rep(c("a", "b"), each = 3L)),
    f = factor(rep(c("a", "b", "c"), 2L))
  )
  formula_result <- JAGS_formula(
    formula = ~ 1 + cs(f | id),
    parameter = "mu",
    data = data,
    prior_list = list(intercept = prior("point", list(location = 0))),
    prior_random = prior_random(
      id = random_block(
        sd = prior("point", list(location = 0.4)),
        rho = prior("point", list(location = 0.2))
      )
    )
  )
  random_term <- formula_result$formula_design$random_effects[[1L]]
  posterior <- matrix(numeric(), nrow = 1L)
  rho <- .bt_random_effect_rho_draws(
    random_term = random_term,
    posterior = posterior,
    missing = "error",
    out_of_support = "error"
  )
  reconstruct <- .bt_random_effect_cholesky_draws
  reconstruction_count <- 0L
  testthat::local_mocked_bindings(
    .bt_random_effect_cholesky_draws = function(...){

      reconstruction_count <<- reconstruction_count + 1L
      reconstruct(...)
    },
    .package = "BayesTools"
  )

  actual <- .bt_JAGS_bridge_marginal_random_coefficient_geometry(
    random_term = random_term,
    posterior = posterior,
    column_scale = rep(0.4, 3L),
    covariance = TRUE,
    structure = "cs"
  )
  expected_correlation <- matrix(rho, nrow = 3L, ncol = 3L)
  diag(expected_correlation) <- 1

  expect_identical(reconstruction_count, 1L)
  expect_equal(
    tcrossprod(actual$factor),
    0.4^2 * expected_correlation,
    tolerance = 1e-15
  )
  expect_equal(actual$covariance, tcrossprod(actual$factor), tolerance = 0)
})

test_that("bridge marginal evaluator supports known group covariance", {

  data <- data.frame(id = factor(c("b", "a", "c", "b")))
  K <- matrix(
    c(4, 1, 0.5, 1, 9, 2, 0.5, 2, 16),
    nrow = 3L,
    byrow = TRUE,
    dimnames = list(c("a", "b", "c"), c("a", "b", "c"))
  )
  formula_result <- JAGS_formula(
    formula = random_effects_formula(
      ~ 1 | id,
      group_covariance = random_group_covariance(K, scale = "none")
    ),
    parameter = "mu",
    data = data,
    prior_list = list(intercept = prior("point", list(location = 0))),
    prior_random = prior_random(
      id = random_block(sd = prior("point", list(location = 0.5)))
    )
  )
  evaluator <- .bt_JAGS_bridge_compile_marginal_random_evaluator(
    formula_design_list = list(mu = formula_result$formula_design),
    marginal_random_spec = list(
      mu = list(
        blocks = "id",
        row_blocks = list(seq_len(nrow(data))),
        factor_state = TRUE
      )
    ),
    formula_data_list = list(mu = data),
    formula_prior_list = list(mu = formula_result$prior_list),
    model_data = list()
  )
  actual_value <- evaluator$covariance(
    samples = numeric(),
    prior_parameters = list(),
    formula_prior_parameters = list(mu = list()),
    formula_parameters = list(mu = rep(0, nrow(data)))
  )$mu
  actual <- .bridge_marginal_random_dense(actual_value)
  expected <- 0.5^2 * K[as.character(data$id), as.character(data$id)]

  expect_identical(actual_value$representation, "factor")
  expect_equal(unname(actual), unname(expected), tolerance = 1e-12)

  compact_value <- evaluator$covariance(
    samples = numeric(),
    prior_parameters = list(),
    formula_prior_parameters = list(mu = list()),
    formula_parameters = list(mu = rep(0, nrow(data))),
    factor_covariance = FALSE,
    factor_state = TRUE
  )$mu
  expect_identical(compact_value$representation, "factor_state")
  expect_true(is.environment(compact_value$contract_id))
  expect_identical(
    names(compact_value$factor_plans[[1L]]),
    c("type", "model_matrix", "group_map", "group_covariance")
  )
  expect_identical(
    names(compact_value$factor_states[[1L]]),
    "coefficient_factor"
  )
  expect_equal(
    unname(.bridge_marginal_random_dense(compact_value)),
    unname(expected),
    tolerance = 1e-12
  )
})

test_that("known group covariance carries the full coefficient covariance", {

  data <- data.frame(
    id = factor(c("b", "a", "c", "b")),
    x = c(-1, 0, 1, 2)
  )
  K <- matrix(
    c(4, 1, 0.5, 1, 9, 2, 0.5, 2, 16),
    nrow = 3L,
    byrow = TRUE,
    dimnames = list(c("a", "b", "c"), c("a", "b", "c"))
  )
  formula_result <- JAGS_formula(
    formula = random_effects_formula(
      ~ 1 + us(1 + x | id),
      group_covariance = random_group_covariance(K, scale = "none")
    ),
    parameter = "mu",
    data = data,
    prior_list = list(intercept = prior("point", list(location = 0))),
    prior_random = prior_random(
      id = random_block(
        sd = prior("point", list(location = 0.5)),
        terms = list(x = prior("point", list(location = 0.3))),
        cor = prior_lkj(eta = 2, include_primitives = TRUE)
      )
    )
  )
  random_term <- formula_result$formula_design$random_effects[[1L]]
  primitive_names <- random_term$correlation$primitive_names
  values <- stats::setNames(rep(0.7, length(primitive_names)), primitive_names)
  posterior <- matrix(values, nrow = 1L, dimnames = list(NULL, names(values)))
  reference <- .bt_random_effect_marginal_covariance_samples(
    design = formula_result$formula_design,
    posterior = posterior,
    prior_list = formula_result$prior_list,
    blocks = "id"
  )$samples[1L, , ]
  evaluator <- .bt_JAGS_bridge_compile_marginal_random_evaluator(
    formula_design_list = list(mu = formula_result$formula_design),
    marginal_random_spec = list(
      mu = list(blocks = "id", row_blocks = list(seq_len(nrow(data))))
    ),
    formula_data_list = list(mu = data),
    formula_prior_list = list(mu = formula_result$prior_list),
    model_data = list()
  )
  actual_value <- evaluator$covariance(
    samples = values,
    prior_parameters = list(),
    formula_prior_parameters = list(mu = list()),
    formula_parameters = list(mu = rep(0, nrow(data)))
  )$mu
  factor <- actual_value$factors[[1L]]

  expect_identical(factor$type, "known_group")
  expect_identical(dim(factor$coefficient_covariance), c(2L, 2L))
  expect_equal(
    factor$coefficient_covariance,
    tcrossprod(factor$coefficient_factor),
    tolerance = 0
  )
  expect_equal(
    unname(.bridge_marginal_random_dense(actual_value)),
    unname(reference),
    tolerance = 1e-12
  )

  compact_value <- evaluator$covariance(
    samples = values,
    prior_parameters = list(),
    formula_prior_parameters = list(mu = list()),
    formula_parameters = list(mu = rep(0, nrow(data))),
    factor_covariance = FALSE
  )$mu
  compact_factor <- compact_value$factors[[1L]]
  compact_coefficient_covariance <- tcrossprod(
    compact_factor$coefficient_factor
  )
  compact_covariance <- compact_factor$group_covariance[
    compact_factor$group_map,
    compact_factor$group_map,
    drop = FALSE
  ] * tcrossprod(
    compact_factor$model_matrix %*% compact_coefficient_covariance,
    compact_factor$model_matrix
  )

  expect_null(compact_factor$coefficient_covariance)
  expect_equal(
    unname(compact_covariance),
    unname(reference),
    tolerance = 1e-12
  )
})

test_that("row-indexed external SD sources remain covariance factors", {

  data <- data.frame(id = factor(c("a", "a", "b", "b")))
  formula_result <- JAGS_formula(
    formula = ~ 1 + random(1 | id, name = "id", covariance = "diag"),
    parameter = "mu",
    data = data,
    prior_list = list(intercept = prior("point", list(location = 0))),
    prior_random = prior_random(
      id = random_block(
        sd_source = random_sd_source("tau", shape = "row")
      )
    )
  )
  tau_names <- paste0("tau[", seq_len(nrow(data)), "]")
  values <- stats::setNames(c(0.2, 0.4, 0.3, 0.5), tau_names)
  posterior <- matrix(values, nrow = 1L, dimnames = list(NULL, names(values)))
  reference <- .bt_random_effect_marginal_covariance_samples(
    design = formula_result$formula_design,
    posterior = posterior,
    prior_list = formula_result$prior_list,
    blocks = "id"
  )$samples[1L, , ]
  evaluator <- .bt_JAGS_bridge_compile_marginal_random_evaluator(
    formula_design_list = list(mu = formula_result$formula_design),
    marginal_random_spec = list(
      mu = list(
        blocks = "id",
        row_blocks = list(1:2, 3:4),
        factor_state = TRUE
      )
    ),
    formula_data_list = list(mu = data),
    formula_prior_list = list(mu = formula_result$prior_list),
    model_data = list()
  )
  actual_value <- evaluator$covariance(
    samples = values,
    prior_parameters = list(),
    formula_prior_parameters = list(mu = list()),
    formula_parameters = list(mu = rep(0, nrow(data)))
  )$mu
  factor <- actual_value$factors[[1L]]

  expect_identical(factor$type, "row_group")
  expect_equal(factor$row_scale, unname(values), tolerance = 0)
  expect_equal(
    unname(.bridge_marginal_random_dense(actual_value)),
    unname(reference),
    tolerance = 1e-12
  )

  compact_value <- evaluator$covariance(
    samples = values,
    prior_parameters = list(),
    formula_prior_parameters = list(mu = list()),
    formula_parameters = list(mu = rep(0, nrow(data))),
    factor_covariance = FALSE,
    factor_state = TRUE
  )$mu
  compact_state <- compact_value$factor_states[[1L]]
  expect_identical(
    names(compact_state),
    c("coefficient_factor", "row_scale")
  )
  expect_equal(compact_state$row_scale, unname(values), tolerance = 0)
  expect_equal(
    unname(.bridge_marginal_random_dense(compact_value)),
    unname(reference),
    tolerance = 1e-12
  )
})
