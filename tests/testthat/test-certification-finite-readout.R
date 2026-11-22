skip_if_not_test_profile("certification")

.finite_readout_model <- function(syntax, data, inits, threads = 1L){

  file <- textConnection(syntax)
  on.exit(close(file))
  model <- rjags::jags.model(file, data = data, inits = inits, n.chains = length(inits),
    n.threads = threads, n.adapt = 0L, quiet = TRUE)
  rjags::adapt(model, n.iter = 25L, end.adaptation = TRUE, progress.bar = "none")
  stats::update(model, n.iter = 11L, progress.bar = "none")
  model
}

.finite_readout_columns <- function(id, width){

  sprintf("BayesTools_finite_readout_%d[%d]", id, seq_len(width))
}

test_that("finite monitor records the actual stock-equivalent update in each chain", {

  skip_if_not_installed("rjags")
  expect_true(.BayesTools_load_module(quiet = TRUE))
  body <- c(
    "effect ~ dcat(c(1, 2))", "slab ~ dnorm(0, 1)",
    "mu <- (effect == 2) * slab", "y ~ dnorm(mu, 2)",
    "bias ~ dcat(c(2, 5, 2e-22))", "yb ~ dbern(rate[bias])",
    "het ~ dbern(0.33)", "yh ~ dbern(het_rate[het + 1])",
    "other ~ dbern(0.7)", "yo ~ dnorm(other, 2)"
  )
  declarations <- c(
    "BayesTools_finite_1 <- bt_finite_tag(bias)",
    "BayesTools_finite_2 <- bt_finite_witness_tag(effect, slab)",
    "BayesTools_finite_3 <- bt_finite_tag(het)"
  )
  data <- list(y = 0.86, yb = 1, yh = 1, yo = 0.3,
               rate = c(0.15, 0.77, 0.3), het_rate = c(0.2, 0.85))
  inits <- lapply(c(141235L, 442183L), function(seed){
    list(effect = 1L, slab = 0, bias = 1L, het = 0L, other = 0L,
         .RNG.name = "base::Super-Duper", .RNG.seed = seed)
  })
  stock <- .finite_readout_model(c("model {", body, "}"), data, inits)
  candidate <- .finite_readout_model(c("model {", body, declarations, "}"), data, inits, threads = 2L)
  other_data <- data
  other_data$yb <- 0
  other_model <- .finite_readout_model(c("model {", body, declarations, "}"), other_data, inits)

  roster <- rjags::list.samplers(stock)
  selected <- vapply(roster, function(nodes) any(nodes %in% c("bias", "effect", "het")), logical(1))
  expect_identical(unname(names(roster)[selected]), rep("base::Finite", 3L))
  names(roster)[selected] <- "BayesTools::Finite"
  expect_identical(rjags::list.samplers(candidate), roster)
  expect_identical(candidate$state(internal = TRUE), stock$state(internal = TRUE))

  ordinary <- c("bias", "effect", "het", "mu", "slab", "other")
  virtual <- paste0("BayesTools_finite_readout_", 1:3)
  expected_bias <- c(2, 5, 2e-22) * data$rate
  expected_bias <- expected_bias / sum(expected_bias)
  expected_het <- c(0.67, 0.33) * data$het_rate
  expected_het <- expected_het / sum(expected_het)
  for (n_iter in c(120L, 54L)) {
    reference <- rjags::coda.samples(stock, ordinary, n_iter,
      thin = 3L, progress.bar = "none")
    recorded <- rjags::coda.samples(candidate, c(ordinary, virtual), n_iter,
      thin = 3L, progress.bar = "none")
    same_names <- rjags::coda.samples(other_model, virtual[1], n_iter,
      thin = 3L, progress.bar = "none")
    expect_identical(candidate$state(internal = TRUE), stock$state(internal = TRUE))

    for (chain in seq_along(reference)) {
      draws <- as.matrix(recorded[[chain]])
      stock_draws <- as.matrix(reference[[chain]])
      expect_identical(draws[, colnames(stock_draws)], stock_draws)
      expect_identical(stats::time(recorded[[chain]]), stats::time(reference[[chain]]))
      expect_equal(nrow(draws), n_iter / 3)

      bias <- draws[, .finite_readout_columns(1L, 8L)]
      effect <- draws[, .finite_readout_columns(2L, 7L)]
      het <- draws[, .finite_readout_columns(3L, 6L)]
      expect_lt(max(abs(sweep(bias[, 1:3], 2L, expected_bias))), 1e-14)
      expect_lt(max(abs(sweep(het[, 1:2], 2L, expected_het))), 1e-14)
      expect_lt(max(abs(sweep(bias[, 4:6], 2L, log(expected_bias)))), 1e-13)
      expect_lt(max(abs(exp(effect[, 3:4]) - effect[, 1:2])), 1e-14)
      expect_identical(bias[, 7L], draws[, "bias"])
      expect_identical(effect[, 5L], draws[, "effect"])
      expect_identical(het[, 5L], draws[, "het"])
      expect_true(all(bias[, 8L] == 1 & effect[, 6L] == 1 & het[, 6L] == 0))
      expect_true(any(effect[, 7L] != draws[, "slab"]))
      log_odds <- log(2) + stats::dnorm(data$y, effect[, 7L], sqrt(0.5), log = TRUE) -
        stats::dnorm(data$y, 0, sqrt(0.5), log = TRUE)
      expect_lt(max(abs(effect[, 2L] - stats::plogis(log_odds))), 1e-14)
      expect_true(all(as.matrix(same_names[[chain]])[, 1] > bias[, 1]))
    }
  }

  saved_state <- unserialize(serialize(candidate$state(internal = TRUE), NULL))
  rm(candidate)
  invisible(gc())
  rebuilt <- .finite_readout_model(c("model {", body, declarations, "}"), data, saved_state)
  control <- .finite_readout_model(c("model {", body, "}"), data, saved_state)
  rebuilt_draws <- rjags::coda.samples(rebuilt, c(ordinary, virtual), 40L,
    thin = 4L, progress.bar = "none")
  control_draws <- rjags::coda.samples(control, ordinary, 40L, thin = 4L, progress.bar = "none")
  expect_identical(rebuilt$state(internal = TRUE), control$state(internal = TRUE))
  for (chain in seq_along(control_draws)) {
    recorded <- as.matrix(rebuilt_draws[[chain]])
    expected <- as.matrix(control_draws[[chain]])
    expect_identical(recorded[, colnames(expected)], expected)
    expect_identical(recorded[, .finite_readout_columns(1L, 8L)[7L]], expected[, "bias"])
  }
})

test_that("finite readout keeps underflowed log masses and refuses ineligible owners", {

  skip_if_not_installed("rjags")
  expect_true(.BayesTools_load_module(quiet = TRUE))
  syntax <- c(
    "model {", "label ~ dcat(c(1, 1, 0))", "y ~ dnorm(location[label], 2)",
    "continuous ~ dnorm(0, 1)", "z ~ dnorm(continuous, 1)",
    "BayesTools_finite_1 <- bt_finite_tag(label)",
    "BayesTools_finite_2 <- bt_finite_tag(continuous)", "}"
  )
  model <- .finite_readout_model(syntax, list(y = 0, z = 1, location = c(0, 40, 0)),
    list(list(label = 1L, continuous = 0, .RNG.name = "base::Super-Duper", .RNG.seed = 3001L)))
  draws <- as.matrix(rjags::coda.samples(model, c("label", "BayesTools_finite_readout_1"), 12L,
    thin = 2L, progress.bar = "none")[[1]])
  recorded <- draws[, .finite_readout_columns(1L, 8L)]
  expect_true(all(recorded[, 1:3] == matrix(c(1, 0, 0), 6, 3, byrow = TRUE)))
  expect_true(all(recorded[, 4L] == 0))
  expect_lt(max(abs(recorded[, 5L] + 1600)), 3e-12)
  expect_true(all(recorded[, 6L] == -Inf))
  expect_true(all(recorded[, 7L] == draws[, "label"]))
  expect_warning(
    rjags::coda.samples(model, c("label", "BayesTools_finite_readout_2"), 2L, progress.bar = "none"),
    "No available BayesTools::Finite sampler owns the marker", fixed = TRUE
  )
})

test_that("public finite readouts follow scientific draws and survive parallel saved-fit extension", {

  skip_if_not_installed("rjags")
  model <- paste(
    "model {", "gate[1] ~ dbern(0.4)", "label ~ dcat(c(2, 1, 1))",
    "y ~ dnorm(gate[1] + label, 1)", "z ~ dnorm(x, 1)", "}",
    sep = "\n"
  )
  args <- list(model, data = list(y = 2, z = 0.2), prior_list = list(x = prior("normal", list(0, 1))),
    chains = 3L, adapt = 50L, burnin = 50L, sample = 100L, thin = 2L, seed = 83251L,
    add_parameters = c("label", "gate[1]"))
  control <- do.call(JAGS_fit, args)
  fit <- do.call(JAGS_fit, c(args, list(sampler_probabilities = c("label", "x", "gate[1]"))))
  parallel <- do.call(JAGS_fit, c(args, list(
    sampler_probabilities = c("label", "x", "gate[1]"), parallel = TRUE, cores = 2L
  )))
  expect_s3_class(fit, "BayesTools_fit")
  expect_s3_class(parallel, "BayesTools_fit")
  expect_null(control$sampler_probabilities)
  expect_named(fit$sampler_probabilities, c("label", "x", "gate[1]"))
  expect_null(fit$sampler_probabilities[["x"]])
  expect_identical(fit$mcmc, control$mcmc)
  expect_identical(fit$end_state, control$end_state)
  expect_identical(fit$mcmc, parallel$mcmc)
  expect_identical(fit$sampler_probabilities, parallel$sampler_probabilities)
  expect_identical(fit$end_state, parallel$end_state)

  check_rows <- function(fit){
    for(node in c("label", "gate[1]")){
      readout <- fit$sampler_probabilities[[node]]
      categories <- if(node == "label") 1:3 else 0:1
      expect_s3_class(readout, "mcmc.list")
      expect_identical(colnames(readout[[1L]]),
        c(paste0("p[", categories, "]"), paste0("log_p[", categories, "]"), "selected", "lower"))
      for(i in seq_along(readout)){
        expect_identical(stats::time(readout[[i]]), stats::time(fit$mcmc[[i]]))
        expect_identical(as.numeric(readout[[i]][, "selected"]), as.numeric(fit$mcmc[[i]][, node]))
        expect_equal(rowSums(readout[[i]][, seq_along(categories), drop = FALSE]),
          rep(1, nrow(readout[[i]])), tolerance = 1e-14)
      }
    }
  }
  check_rows(fit)
  settings <- list(max_extend = 1L, sample_extend = 15L, max_Rhat = NULL, min_ESS = NULL,
    max_error = NULL, max_SD_error = NULL, max_time = NULL)
  extended <- JAGS_extend(unserialize(serialize(fit, NULL)), autofit_control = settings,
    parallel = TRUE, cores = 2L)
  expect_identical(extended$sample, 115)
  expect_identical(extended$end_state, JAGS_extend(control, autofit_control = settings)$end_state)
  expect_null(extended$sampler_probabilities[["x"]])
  expect_identical(names(extended$sampler_probabilities), names(fit$sampler_probabilities))
  check_rows(extended)
  expect_identical(lapply(extended$sampler_probabilities$label, function(draws){
    as.matrix(draws)[seq_len(100L), , drop = FALSE]
  }), lapply(fit$sampler_probabilities$label, as.matrix))
})

test_that("validated package declarations survive recompilation without replaying the callback", {

  skip_if_not_installed("rjags")
  called <- 0L
  fit <- JAGS_fit("model { y ~ dnorm(x + w, 1) }", data = list(y = 0.1),
    prior_list = list(x = prior("normal", list(0, 1)), w = prior("normal", list(0, 1))),
    chains = 2L, adapt = 50L, burnin = 50L, sample = 100L, seed = 43010L,
    sampler_joint = function(prior_list, formula_design, data){
      called <<- called + 1L
      list(markers = list(
        list(name = "package_marker_x", distribution = "dbt_block", arguments = list("x", 1)),
        list(name = "package_marker_w", distribution = "dbt_block", arguments = list("w", 1))
      ))
    }
  )
  expect_s3_class(fit, "BayesTools_fit")
  expect_identical(called, 1L)
  expect_setequal(fit$samplers$Node[fit$samplers$Sampler == "BayesTools::FactorSlice"], c("x", "w"))
  settings <- list(max_extend = 1L, sample_extend = 13L, max_Rhat = NULL, min_ESS = NULL,
    max_error = NULL, max_SD_error = NULL, max_time = NULL)
  extended <- JAGS_extend(unserialize(serialize(fit, NULL)), autofit_control = settings)
  expect_identical(extended$sample, 113)
  expect_identical(called, 1L)
  expect_identical(extended$samplers, fit$samplers)
  expect_identical(attr(extended, "sampler_declarations"), attr(fit, "sampler_declarations"))
})
