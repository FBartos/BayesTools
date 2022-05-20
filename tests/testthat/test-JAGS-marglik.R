context("JAGS marginal likelihood functions")

test_that("JAGS model functions work (simple)", {

  skip_if_not_installed("rjags")
  all_priors  <- list(
    p1  = prior("normal", list(0, 1)),
    p2  = prior("normal", list(0, 1), list(1, Inf)),
    p3  = prior("lognormal", list(0, .5)),
    p4  = prior("t", list(0, .5, 5)),
    p5  = prior("Cauchy", list(1, 0.1), list(-10, 0)),
    p6  = prior("gamma", list(2, 1)),
    p7  = prior("invgamma", list(3, 2), list(1, 3)),
    p8  = prior("exp", list(1.5)),
    p9  = prior("beta", list(3, 2)),
    p10 = prior("uniform", list(1, 5)),
    PET = prior_PET("normal", list(0, 1)),
    PEESE = prior_PEESE("gamma", list(1, 1))
    #p12 = prior("bernoulli", list(0.75)) discrete priors are not supported with bridgesampling
  )
  log_posterior <- function(parameters, data){
    return(0)
  }


  for(i in seq_along(all_priors)){
    prior_list   <- all_priors[i]
    model_syntax <- JAGS_add_priors("model{}", prior_list)
    monitor      <- JAGS_to_monitor(prior_list)
    inits        <- JAGS_get_inits(prior_list, chains = 2, seed = 1)

    set.seed(1)
    model   <- rjags::jags.model(file = textConnection(model_syntax), inits = inits, n.chains = 2, quiet = TRUE)
    samples <- rjags::coda.samples(model = model, variable.names = monitor, n.iter = 5000, quiet = TRUE, progress.bar = "none")
    marglik <- JAGS_bridgesampling(samples, prior_list = prior_list, data = list(), log_posterior = log_posterior)
    expect_equal(marglik$logml, 0, tolerance = 1e-2)
  }

})

test_that("JAGS model functions work (vector)", {

  skip_if_not_installed("rjags")
  all_priors  <- list(
    p1  = prior("mnormal", list(mean = 0, sd = 1, K = 3),),
    p2  = prior("mcauchy", list(location = 0, scale = 1.5, K = 2)),
    p3  = prior("mt",      list(location = 2, scale = 0.5, df = 5, K = 2))
  )
  log_posterior <- function(parameters, data){
    return(0)
  }


  for(i in seq_along(all_priors)){
    prior_list   <- all_priors[i]
    model_syntax <- JAGS_add_priors("model{}", prior_list)
    monitor      <- JAGS_to_monitor(prior_list)
    inits        <- JAGS_get_inits(prior_list, chains = 2, seed = 1)

    set.seed(1)
    model   <- rjags::jags.model(file = textConnection(model_syntax), inits = inits, n.chains = 2, quiet = TRUE)
    samples <- rjags::coda.samples(model = model, variable.names = monitor, n.iter = 10000, quiet = TRUE, progress.bar = "none")
    marglik <- JAGS_bridgesampling(samples, prior_list = prior_list, data = list(), log_posterior = log_posterior)
    expect_equal(marglik$logml, 0, tolerance = 5*1e-2) # the mCauchy is a bit more variable
  }

})

test_that("JAGS model functions work (factor)", {

  skip_if_not_installed("rjags")
  all_priors   <- list(
    p1  = prior_factor("mnorm", list(mean = 0, sd = 1),    contrast = "orthonormal"),
    p2  = prior_factor("beta",  list(alpha = 1, beta = 1), contrast = "dummy"),
    p3  = prior_factor("beta",  list(alpha = 2, beta = 2), contrast = "dummy")
  )

  # add levels
  attr(all_priors[[1]], "levels") <- 3
  attr(all_priors[[2]], "levels") <- 2
  attr(all_priors[[3]], "levels") <- 3
  log_posterior <- function(parameters, data){
    return(0)
  }


  for(i in seq_along(all_priors)){
    prior_list   <- all_priors[i]
    model_syntax <- JAGS_add_priors("model{}", prior_list)
    monitor      <- JAGS_to_monitor(prior_list)
    inits        <- JAGS_get_inits(prior_list, chains = 2, seed = 1)

    set.seed(1)
    model   <- rjags::jags.model(file = textConnection(model_syntax), inits = inits, n.chains = 2, quiet = TRUE)
    samples <- rjags::coda.samples(model = model, variable.names = monitor, n.iter = 10000, quiet = TRUE, progress.bar = "none")
    marglik <- JAGS_bridgesampling(samples, prior_list = prior_list, data = list(), log_posterior = log_posterior)
    expect_equal(marglik$logml, 0, tolerance = 1e-2)
  }

})

test_that("JAGS model functions work (spike and slab)", {

  skip_if_not_installed("rjags")
  all_priors   <- list(
    p1  = prior_spike_and_slab("normal", list(mean = 0, sd = 1))
  )

  log_posterior <- function(parameters, data){
    return(0)
  }


  for(i in seq_along(all_priors)){
    prior_list   <- all_priors[i]
    model_syntax <- JAGS_add_priors("model{}", prior_list)
    monitor      <- JAGS_to_monitor(prior_list)
    inits        <- JAGS_get_inits(prior_list, chains = 2, seed = 1)

    set.seed(1)
    model   <- rjags::jags.model(file = textConnection(model_syntax), inits = inits, n.chains = 2, quiet = TRUE)
    samples <- rjags::coda.samples(model = model, variable.names = monitor, n.iter = 10000, quiet = TRUE, progress.bar = "none")
    marglik <- JAGS_bridgesampling(samples, prior_list = prior_list, data = list(), log_posterior = log_posterior)
    expect_equal(marglik$logml, 0, tolerance = 1e-2)
  }

})

test_that("JAGS model functions work (weightfunctions)", {

  skip_if_not_installed("rjags")
  all_priors  <- list(
    prior_weightfunction("one.sided", list(c(.05), c(1, 1))),
    prior_weightfunction("one.sided", list(c(.05, 0.10), c(1, 2, 3))),
    prior_weightfunction("one.sided", list(c(.05, 0.60), c(1, 1), c(1, 5))),
    prior_weightfunction("two.sided", list(c(.05), c(1, 1)))
  )
  log_posterior <- function(parameters, data){
    return(0)
  }


  for(i in seq_along(all_priors)){
    prior_list   <- all_priors[i]
    model_syntax <- JAGS_add_priors("model{}", prior_list)
    monitor      <- JAGS_to_monitor(prior_list)
    inits        <- JAGS_get_inits(prior_list, chains = 2, seed = 1)

    set.seed(1)
    model   <- rjags::jags.model(file = textConnection(model_syntax), inits = inits, n.chains = 2, quiet = TRUE)
    samples <- rjags::coda.samples(model = model, variable.names = monitor, n.iter = 5000, quiet = TRUE, progress.bar = "none")
    marglik <- JAGS_bridgesampling(samples, prior_list = prior_list, data = list(), log_posterior = log_posterior)
    expect_equal(marglik$logml, 0, tolerance = 1e-2)
  }

})

test_that("JAGS model functions work (complex scenario)", {

  skip_if_not_installed("rjags")
  # tests different model estimation techniques and passing additional arguments
  set.seed(1)
  data <- list(
    x = rnorm(50, 0, .5),
    N = 50
  )
  priors1  <- list(
    m  = prior("normal", list(0, 1)),
    s  = prior("normal", list(0, 1), list(0, Inf))
  )
  priors2  <- list(
    m  = prior("normal", list(0, 1)),
    s  = prior("spike",  list(1))
  )
  log_posterior <- function(parameters, data, return3){
    if(return3){
      return(3)
    }else{
      return(sum(stats::dnorm(data$x, mean = parameters[["m"]], sd = parameters[["s"]], log = TRUE)))
    }
  }
  model_syntax <-
  "model{
    for(i in 1:N){
      x[i] ~ dnorm(m, pow(s, -2))
    }
  }"



  model1   <- rjags::jags.model(
    file     = textConnection(JAGS_add_priors(model_syntax, priors1)),
    inits    = JAGS_get_inits(priors1, chains = 2, seed = 1),
    n.chains = 2,
    data     = data,
    quiet    = TRUE)
  samples1 <- rjags::jags.samples(
    model          = model1,
    variable.names = JAGS_to_monitor(priors1),
    data           = data,
    n.iter         = 5000,
    quiet          = TRUE,
    progress.bar   = "none")
  marglik1 <- JAGS_bridgesampling(
    samples1,
    prior_list    = priors1,
    data          = data,
    log_posterior = log_posterior,
    return3       = FALSE)

  runjags::runjags.options(silent.jags = TRUE, silent.runjags = TRUE)
  fit2 <- runjags::run.jags(
    model           = JAGS_add_priors(model_syntax, priors2),
    data            = data,
    inits           = JAGS_get_inits(priors2, chains = 2, seed = 1),
    monitor         = JAGS_to_monitor(priors2),
    n.chains        = 2,
    sample          = 5000,
    burnin          = 1000,
    adapt           = 500,
    summarise       = FALSE
  )
  marglik2 <- JAGS_bridgesampling(
    fit2,
    data          = data,
    prior_list    = priors2,
    log_posterior = log_posterior,
    return3       = FALSE)

  marglik3 <- JAGS_bridgesampling(
    fit2,
    data          = data,
    prior_list    = priors2,
    log_posterior = log_posterior,
    return3       = TRUE)


  expect_equal(marglik1$logml, -31.944, tolerance = 1e-2)
  expect_equal(marglik2$logml, -52.148, tolerance = 1e-2)
  expect_equal(marglik3$logml,   1.489, tolerance = 1e-2)
})

test_that("JAGS model functions work (formula)",{

  set.seed(1)

  data_formula <- data.frame(
    x_cont1 = rnorm(300),
    x_fac2t = factor(rep(c("A", "B"), 150), levels = c("A", "B")),
    x_fac3t = factor(rep(c("A", "B", "C"), 100), levels = c("A", "B", "C"))
  )
  data <- list(
    y = rnorm(300, .4 * data_formula$x_cont1 + ifelse(data_formula$x_fac3t == "A", 0.0, ifelse(data_formula$x_fac3t == "B", -0.2, 0.4)), ifelse(data_formula$x_fac2t == "A", 0.5, 1)),
    N = 300
  )

  # create an empty model ----
  formula_list0 <- list(
    mu    = ~ x_cont1 + x_fac3t
  )
  formula_data_list0 <- list(
    mu    = data_formula
  )
  formula_prior_list0 <- list(
    mu    = list(
      "intercept"       = prior("normal", list(0, 5)),
      "x_cont1"         = prior("normal", list(0, 1)),
      "x_fac3t"         = prior_factor("normal", contrast = "treatment", list(0, 1))
    )
  )
  prior_list0 <- list(
    sigma = prior("lognormal", list(0, 1))
  )
  model_syntax0 <- "model{}"

  fit0     <- JAGS_fit(
    model_syntax = model_syntax0, data = list(), prior_list = prior_list0,
    formula_list = formula_list0, formula_data_list = formula_data_list0, formula_prior_list = formula_prior_list0)

  log_posterior0 <- function(parameters, data){
    return(0)
  }

  marglik0 <- JAGS_bridgesampling(
    fit                = fit0,
    log_posterior      = log_posterior0,
    data               = list(),
    prior_list         = prior_list0,
    formula_list       = formula_list0,
    formula_data_list  = formula_data_list0,
    formula_prior_list = formula_prior_list0
  )

  expect_equal(marglik0$logml,   0, tolerance = 1e-3)


  # create model with mix of a formula and free parameters ---
  formula_list1 <- list(
    mu    = ~ x_cont1 + x_fac3t
  )
  formula_data_list1 <- list(
    mu    = data_formula
  )
  formula_prior_list1 <- list(
    mu    = list(
      "intercept"       = prior("normal", list(0, 5)),
      "x_cont1"         = prior("normal", list(0, 1)),
      "x_fac3t"         = prior_factor("normal", contrast = "treatment", list(0, 1))
    )
  )
  prior_list1 <- list(
    sigma = prior("lognormal", list(0, 1))
  )
  model_syntax1 <- paste0(
    "model{\n",
    "for(i in 1:N){\n",
    "  y[i] ~ dnorm(mu[i], 1/pow(sigma, 2))\n",
    "}\n",
    "}"
  )

  fit1     <- JAGS_fit(
    model_syntax = model_syntax1, data = data, prior_list = prior_list1,
    formula_list = formula_list1, formula_data_list = formula_data_list1, formula_prior_list = formula_prior_list1)

  log_posterior1 <- function(parameters, data){
    return(sum(stats::dnorm(data$y, mean = parameters[["mu"]], sd = parameters[["sigma"]], log = TRUE)))
  }

  marglik1 <- JAGS_bridgesampling(
    fit                = fit1,
    log_posterior      = log_posterior1,
    data               = data,
    prior_list         = prior_list1,
    formula_list       = formula_list1,
    formula_data_list  = formula_data_list1,
    formula_prior_list = formula_prior_list1)

  # more of a consistency test
  expect_equal(marglik1$logml, -370.87, tolerance = 1e-2)


  # create model with mix of a formula and free scaled parameters ---
  prior_list1s         <- prior_list1
  prior_list1s$scale3  <- prior("point", parameters = list(location = 1/3))
  formula_prior_list1s <- list(
    mu    = list(
      "intercept"       = prior("normal", list(0, 5)),
      "x_cont1"         = prior("normal", list(0, 1/2)),
      "x_fac3t"         = prior_factor("normal", contrast = "treatment", list(0, 1*3))
    )
  )
  attr(formula_prior_list1s$mu$x_cont1, "multiply_by") <- 2
  attr(formula_prior_list1s$mu$x_fac3t, "multiply_by") <- "scale3"

  fit1s     <- JAGS_fit(
    model_syntax = model_syntax1, data = data, prior_list = prior_list1s,
    formula_list = formula_list1, formula_data_list = formula_data_list1, formula_prior_list = formula_prior_list1s)

  log_posterior1s <- function(parameters, data){
    return(sum(stats::dnorm(data$y, mean = parameters[["mu"]], sd = parameters[["sigma"]], log = TRUE)))
  }

  marglik1s <- JAGS_bridgesampling(
    fit                = fit1s,
    log_posterior      = log_posterior1s,
    data               = data,
    prior_list         = prior_list1s,
    formula_list       = formula_list1,
    formula_data_list  = formula_data_list1,
    formula_prior_list = formula_prior_list1s)

  # more of a consistency test
  expect_equal(marglik1$logml, marglik1s$logml, tolerance = 1e-2)


  # create model with two formulas ---
  formula_list2 <- list(
    mu    = ~ x_cont1 + x_fac3t,
    sigma = ~ x_fac2t
  )

  formula_data_list2 <- list(
    mu    = data_formula,
    sigma = data_formula
  )

  formula_prior_list2 <- list(
    mu    = list(
      "intercept"       = prior("normal", list(0, 5)),
      "x_cont1"         = prior("normal", list(0, 1)),
      "x_fac3t"         = prior_factor("normal", contrast = "treatment", list(0, 1))
    ),
    sigma = list(
      "intercept"       = prior("normal", list(0, 1)),
      "x_fac2t"         = prior_factor("normal",  contrast = "treatment",   list(0, 1))
    )
  )
  model_syntax2 <- paste0(
    "model{\n",
    "for(i in 1:N){\n",
    "  y[i] ~ dnorm(mu[i], 1/pow(exp(sigma[i]), 2))\n",
    "}\n",
    "}"
  )

  fit2 <- JAGS_fit(
    model_syntax = model_syntax2, data = data, prior_list = NULL,
    formula_list = formula_list2, formula_data_list = formula_data_list2, formula_prior_list = formula_prior_list2)

  log_posterior2 <- function(parameters, data){
    return(sum(stats::dnorm(data$y, mean = parameters[["mu"]], sd = exp(parameters[["sigma"]]), log = TRUE)))
  }

  marglik2 <- JAGS_bridgesampling(
    fit                = fit2,
    log_posterior      = log_posterior2,
    data               = data,
    prior_list         = NULL,
    formula_list       = formula_list2,
    formula_data_list  = formula_data_list2,
    formula_prior_list = formula_prior_list2)

  # more of a consistency test
  expect_equal(marglik2$logml, -351.43, tolerance = 1e-2)
})
