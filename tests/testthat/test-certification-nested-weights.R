skip_if_not_test_profile("certification")

test_that("compiled nested weights and finite label densities use the declared orientation", {

  skip_if_not_installed("rjags")
  rjags::load.module("dic", quiet = TRUE)
  priors <- nested_weight_priors()
  syntax <- "model{\n y ~ dnorm(omega[2] + omega[4], 4)\n nuisance ~ dnorm(0, 1)\n}"
  fit_priors <- .bt_nested_weights_bind(priors, syntax, NULL, TRUE)$prior_list
  plan <- .bt_nested_weights_plan(fit_priors$bias)
  model <- JAGS_add_priors(syntax, fit_priors)
  a <- .35
  v <- .6
  u <- .4
  b <- a * v
  refinement <- a + (1 - a) * u
  expected_weights <- rbind(
    c(1, refinement, a, a, b),
    rep(1, 5),
    c(1, 1, 1, .35, .35),
    c(1, 1, a, a, b)
  )
  chart <- setNames(list(a, v, u), plan$nodes)
  log_density <- vapply(seq_len(nrow(expected_weights)), function(label){
    data <- c(list(y = .9, bias_indicator = label), chart)
    connection <- textConnection(model)
    on.exit(close(connection))
    compiled <- rjags::jags.model(connection, data = data, quiet = TRUE)
    draws <- rjags::jags.samples(compiled, c("omega", "deviance"), n.iter = 1L)
    expect_equal(as.numeric(draws$omega), expected_weights[label, ], tolerance = 1e-14)
    -as.numeric(draws$deviance) / 2
  }, numeric(1))
  spec <- selection_backend_spec(priors$bias, include_init = FALSE)
  expected <- stats::dnorm(.9, expected_weights[, 2] + expected_weights[, 4], .5, log = TRUE) +
    log(spec$prior_weights) + log(c(6 * a * (1 - a), 1, 1, 2 * a))
  expect_equal(log_density - log_density[2L], expected - expected[2L], tolerance = 1e-12)
  normalize <- function(log_values){
    values <- exp(log_values - max(log_values))
    values / sum(values)
  }
  expect_equal(normalize(log_density), normalize(expected), tolerance = 1e-12)

  expect_true(.BayesTools_load_module(quiet = TRUE))
  markers <- "BayesTools_finite_1 <- bt_finite_tag(bias_indicator)"
  connection <- textConnection(.bt_sampler_insert_markers(model, markers))
  on.exit(close(connection))
  compiled <- rjags::jags.model(connection, data = c(list(y = .9), chart), quiet = TRUE)
  readout <- "BayesTools_finite_readout_1"
  expect_identical(rjags::list.samplers(compiled)[["BayesTools::Finite"]], "bias_indicator")
  recorded <- as.matrix(rjags::coda.samples(compiled, readout, 4L, progress.bar = "none"))
  probs <- recorded[, sprintf("%s[%d]", readout, 1:4), drop = FALSE]
  expect_equal(unname(probs), matrix(rep(normalize(expected), each = 4L), nrow = 4L),
               tolerance = 1e-13)
  likelihood <- stats::dnorm(.9, expected_weights[, 2] + expected_weights[, 4], .5, log = TRUE)
  expected_pair_odds <- log(spec$prior_weights[1L] / spec$prior_weights[4L]) +
    likelihood[1L] - likelihood[4L] + log(3 * (1 - a))
  expect_equal(unname(log(probs[, 1L]) - log(probs[, 4L])),
               rep(expected_pair_odds, 4L), tolerance = 1e-13)
})

test_that("nested fits reconstruct from saved state and retain the original raw-monitor graph", {

  skip_if_not_installed("rjags")
  priors <- nested_weight_priors()
  priors$bias <- prior_mixture(c(as.list(priors$bias), list(prior_PET("normal", list(0, 1)))))
  priors$mu <- prior("normal", list(0, 1))
  syntax <- "model{ y ~ dnorm(mu + PET + omega[2] + omega[4], 4) }"
  run <- function(nested, add = NULL, probabilities = NULL){
    JAGS_fit(
      syntax, data = list(y = .9), prior_list = priors, add_parameters = add,
      sampler_gates = TRUE, sampler_nested_weights = nested,
      sampler_probabilities = probabilities,
      chains = 2, adapt = 100, burnin = 50, sample = 100, seed = 381
    )
  }
  fit <- run(TRUE)
  expect_s3_class(fit, "BayesTools_fit")
  layout <- attr(fit, "sampler_declarations")$nested_weights
  expect_true(layout$realized)
  plan <- .bt_nested_weights_plan(attr(fit, "prior_list")$bias)
  expect_setequal(
    fit$samplers$Node[fit$samplers$Sampler == "BayesTools::GatedSlice"],
    c(unname(plan$nodes), "PET_1")
  )
  expect_true(any(fit$samplers$Node == "bias_indicator" & fit$samplers$Sampler == "base::Finite"))
  draws <- as.matrix(fit$mcmc)
  expect_setequal(colnames(draws), c("bias_indicator", paste0("omega[", 1:5, "]"), "mu", "PET"))

  control <- list(max_Rhat = NULL, min_ESS = NULL, max_error = NULL, max_SD_error = NULL,
                  sample_extend = 20, max_extend = 1)
  reloaded <- unserialize(serialize(fit, NULL))
  extended <- JAGS_extend(reloaded, autofit_control = control)
  expect_identical(extended$model, fit$model)
  expect_identical(extended$data, fit$data)
  expect_identical(extended$samplers, fit$samplers)
  expect_identical(attr(extended, "sampler_declarations")$nested_weights, layout)
  for(i in seq_along(fit$mcmc)){
    expect_equal(nrow(extended$mcmc[[i]]), 120)
    expect_identical(extended$mcmc[[i]][1:100, ], fit$mcmc[[i]][1:100, ])
  }

  # An explicit raw request must return the same fitted raw nodes and draws.
  ordinary <- run(FALSE, "eta_component_1")
  fallback <- run(TRUE, "eta_component_1")
  declined <- attr(fallback, "sampler_declarations")$nested_weights
  expect_false(declined$realized)
  expect_match(declined$reason, "eta_component_1")
  expect_identical(fallback$model, ordinary$model)
  expect_identical(fallback$mcmc, ordinary$mcmc)
  expect_true(all(paste0("eta_component_1[", 1:4, "]") %in% colnames(as.matrix(fallback$mcmc))))

  # Readout requests also refer to the original graph, even if the requested
  # continuous node cannot be claimed by the finite sampler.
  request <- "eta_component_1 [ 2 ]"
  ordinary_readout <- run(FALSE, probabilities = request)
  fallback_readout <- run(TRUE, probabilities = request)
  declined <- attr(fallback_readout, "sampler_declarations")$nested_weights
  expect_false(declined$realized)
  expect_match(declined$reason, "eta_component_1")
  expect_length(fallback_readout$end_state[[1L]]$eta_component_1, 4L)
  expect_identical(fallback_readout$model, ordinary_readout$model)
  expect_identical(fallback_readout$mcmc, ordinary_readout$mcmc)
  expect_identical(fallback_readout$sampler_probabilities, list("eta_component_1[2]" = NULL))

  bias_readout <- run(TRUE, probabilities = "bias_indicator")
  expect_true(attr(bias_readout, "sampler_declarations")$nested_weights$realized)
  expect_s3_class(bias_readout$sampler_probabilities$bias_indicator, "mcmc.list")
})
