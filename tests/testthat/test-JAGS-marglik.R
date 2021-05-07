context("JAGS marginal likelihood functions")

test_that("JAGS model functions work (simple)", {


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

test_that("JAGS model functions work (weightfunctions)", {


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
