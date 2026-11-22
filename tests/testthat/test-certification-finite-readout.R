skip_if_not_test_profile("certification")

.finite_readout_model <- function(syntax, data, inits, threads = 1L){

  file <- textConnection(syntax)
  on.exit(close(file))
  .bt_rjags_prepare_finite_readout()
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
  unavailable <- rjags::coda.samples(model, c("label", "BayesTools_finite_readout_2"),
    2L, progress.bar = "none")
  expect_true(all(unavailable[[1]][, "BayesTools_finite_readout_2"] == -1))
})

test_that("recording preserves competing stock joint updates and the physical RNG", {

  skip_if_not_installed("rjags")
  model <- paste(c(
    "model {", "x ~ dbern(0.5)", "w ~ dbern(0.5)", "total ~ dsum(x, w)",
    "recorded ~ dbern(0.4)", "ordinary ~ dcat(c(0.3, 0.7))",
    "y ~ dnorm(recorded + ordinary + offset, 1)", "unused ~ dbern(0.5)", "}"
  ), collapse = "\n")
  args <- list(model, data = list(total = 1L, y = 1),
    prior_list = list(offset = prior("normal", list(0, 1))),
    add_parameters = c("x", "w", "recorded", "ordinary"),
    chains = 2L, adapt = 500L, burnin = 50L, sample = 100L, thin = 3L, seed = 71123L)
  stock <- do.call(JAGS_fit, args)
  fit <- do.call(JAGS_fit, c(args, list(sampler_probabilities = c("x", "recorded", "unused"))))
  expect_s3_class(fit, "BayesTools_fit")
  expect_identical(fit$mcmc, stock$mcmc)
  expect_identical(fit$end_state, stock$end_state)
  roster <- stock$samplers
  expect_setequal(roster$Node[roster$Sampler == "bugs::DiscreteDSum"], c("x", "w"))
  roster$Sampler[roster$Node == "recorded"] <- "BayesTools::Finite"
  expect_identical(fit$samplers, roster)
  expect_null(fit$sampler_probabilities$x)
  expect_null(fit$sampler_probabilities$unused)
  expect_s3_class(fit$sampler_probabilities$recorded, "mcmc.list")
  for(chain in fit$mcmc){
    expect_setequal(chain[, "x"], 0:1)
    expect_true(all(chain[, "x"] + chain[, "w"] == 1))
  }
})

test_that("a JAGS module loaded after registration keeps priority over the recorder", {

  skip_if_not_installed("rjags")
  expect_true(.BayesTools_load_module(quiet = TRUE))
  .bt_rjags_prepare_finite_readout()
  rjags::load.module("mix", quiet = TRUE)
  factories <- rjags::list.factories("sampler")$factory
  expect_lt(match("mix::LDA", factories), match("BayesTools::Finite", factories))
  expect_equal(match("BayesTools::Finite", factories) + 1L, match("base::Finite", factories))
  body <- c(
    "topic_prob[1:2] ~ ddirch(alpha[])",
    "for(j in 1:3) {", "topic[j] ~ dcat(topic_prob[])",
    "word[j] ~ dcat(word_prob[topic[j],])", "}",
    "for(k in 1:2) { word_prob[k,1:3] ~ ddirch(beta[]) }"
  )
  data <- list(alpha = c(1, 1), beta = c(1, 1, 1), word = c(1, 2, 3))
  inits <- list(list(.RNG.name = "base::Super-Duper", .RNG.seed = 19343L))
  stock <- .finite_readout_model(c("model {", body, "}"), data, inits)
  fit <- .finite_readout_model(c("model {", body, "BayesTools_finite_1 <- bt_finite_tag(topic[1])", "}"),
    data, inits)
  roster <- rjags::list.samplers(stock)
  expect_true("topic[1]" %in% roster[["mix::LDA"]])
  expect_identical(rjags::list.samplers(fit), roster)
  readout <- "BayesTools_finite_readout_1"
  expected <- rjags::coda.samples(stock, "topic", 40L, thin = 2L, progress.bar = "none")
  recorded <- rjags::coda.samples(fit, c("topic", readout), 40L, thin = 2L, progress.bar = "none")
  expect_identical(as.matrix(recorded[[1L]])[, colnames(expected[[1L]])], as.matrix(expected[[1L]]))
  expect_true(all(recorded[[1L]][, readout] == -1))
  expect_identical(fit$state(internal = TRUE), stock$state(internal = TRUE))
})

test_that("wide categorical indicators are recorded only with active stock Finite", {

  skip_if_not_installed("rjags")
  expect_true(.BayesTools_load_module(quiet = TRUE))
  body <- c("z ~ dcat(weights[])", "y ~ dnorm(z / 100, 2)")
  tagged <- c("model {", body, "BayesTools_finite_1 <- bt_finite_tag(z)", "}")
  weights <- seq_len(101L)
  data <- list(y = 0.23, weights = weights)
  inits <- list(list(z = 30L, .RNG.name = "base::Super-Duper", .RNG.seed = 12234L))
  stock <- .finite_readout_model(c("model {", body, "}"), data, inits)
  candidate <- .finite_readout_model(tagged, data, inits)
  expect_identical(names(rjags::list.samplers(stock)), "base::Finite")
  expect_identical(names(rjags::list.samplers(candidate)), "BayesTools::Finite")
  virtual <- "BayesTools_finite_readout_1"
  expected <- rjags::coda.samples(stock, "z", 6L, progress.bar = "none")
  retained <- rjags::coda.samples(candidate, c("z", virtual), 6L, progress.bar = "none")
  values <- retained[[1L]][, .finite_readout_columns(1L, 204L)]
  exact <- weights * stats::dnorm(data$y, seq_along(weights) / 100, sqrt(0.5))
  expect_lt(max(abs(sweep(values[, seq_along(weights)], 2L, exact / sum(exact)))), 1e-14)
  expect_identical(retained[[1L]][, "z"], expected[[1L]][, "z"])
  expect_identical(candidate$state(internal = TRUE), stock$state(internal = TRUE))

  rjags::set.factory("base::Finite", "sampler", FALSE)
  withr::defer(rjags::set.factory("base::Finite", "sampler", TRUE))
  stock <- .finite_readout_model(c("model {", body, "}"), data, inits)
  candidate <- .finite_readout_model(tagged, data, inits)
  expect_false("BayesTools::Finite" %in% names(rjags::list.samplers(candidate)))
  expect_identical(rjags::list.samplers(candidate), rjags::list.samplers(stock))
  expected <- rjags::coda.samples(stock, "z", 6L, progress.bar = "none")
  retained <- rjags::coda.samples(candidate, c("z", virtual), 6L, progress.bar = "none")
  expect_true(all(retained[[1L]][, virtual] == -1))
  expect_identical(retained[[1L]][, "z"], expected[[1L]][, "z"])
  expect_identical(candidate$state(internal = TRUE), stock$state(internal = TRUE))
})

test_that("compiled aliases share finite readouts and real scientific names remain usable", {

  skip_if_not_installed("rjags")
  model <- paste(c(
    "model {", "z ~ dbern(0.5)", "a <- z", "g[1] ~ dbern(0.5)",
    "y ~ dnorm(z + g[1] + offset, 1)", "}"
  ), collapse = "\n")
  args <- list(model, data = list(y = 1), prior_list = list(offset = prior("normal", list(0, 1))),
    add_parameters = c("z", "g"), chains = 1L, adapt = 50L, burnin = 50L, sample = 100L, seed = 3109L)
  stock <- do.call(JAGS_fit, args)
  fit <- do.call(JAGS_fit, c(args, list(sampler_probabilities = c("z", "a", "g", "g[1]"))))
  expect_s3_class(fit, "BayesTools_fit")
  expect_identical(fit$mcmc, stock$mcmc)
  expect_identical(fit$end_state, stock$end_state)
  expect_s3_class(fit$sampler_probabilities$z, "mcmc.list")
  expect_s3_class(fit$sampler_probabilities$g, "mcmc.list")
  expect_identical(fit$sampler_probabilities$z, fit$sampler_probabilities$a)
  expect_identical(fit$sampler_probabilities$g, fit$sampler_probabilities[["g[1]"]])

  scientific <- c("BayesTools_finite_readout_user", "BayesTools_finite_readout_1")
  syntax <- paste0("model { ", scientific[1], " ~ dbern(0.5) ", scientific[2], " ~ dbern(0.5)",
    " y ~ dnorm(", scientific[1], " + ", scientific[2], " + offset, 1) }")
  ordinary <- JAGS_fit(syntax, data = list(y = 1), prior_list = args$prior_list,
    add_parameters = scientific, chains = 1L, adapt = 50L, burnin = 50L, sample = 100L, seed = 3109L)
  expect_s3_class(ordinary, "BayesTools_fit")
  expect_true(all(scientific %in% colnames(ordinary$mcmc[[1L]])))
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
