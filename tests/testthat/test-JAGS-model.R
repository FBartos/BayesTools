context("JAGS model functions")

test_that("JAGS model functions work (simple)", {

  model_syntax <- "model{}"
  priors       <- list(
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
    p11 = prior("point", list(1)),
    PET = prior_PET("normal", list(0, 1)),
    PEESE = prior_PEESE("gamma", list(1, 1))
  )

  model_syntax <- JAGS_add_priors(model_syntax, priors)
  monitor      <- JAGS_to_monitor(priors)
  inits        <- JAGS_get_inits(priors, chains = 2, seed = 1)

  set.seed(1)
  model   <- rjags::jags.model(file = textConnection(model_syntax), inits = inits, n.chains = 2, quiet = TRUE)
  samples <- rjags::coda.samples(model = model, variable.names = monitor, n.iter = 5000, quiet = TRUE, progress.bar = "none")
  samples <- do.call(rbind, samples)


  for(i in seq_along(priors)){
    expect_doppelganger(paste0("JAGS-model-prior-",i), function(){
      hist(samples[,names(priors)[i]], breaks = 50, main = print(priors[[i]], plot = TRUE), freq = FALSE)
      lines(density(priors[[i]]))
    })
  }
})

test_that("JAGS model functions work (weightfunctions)", {

  priors       <- list(
    prior_weightfunction("one.sided", list(c(.05), c(1, 1))),
    prior_weightfunction("one.sided", list(c(.05, 0.10), c(1, 2, 3))),
    prior_weightfunction("one.sided", list(c(.05, 0.60), c(1, 1), c(1, 5))),
    prior_weightfunction("two.sided", list(c(.05), c(1, 1))),
    prior_weightfunction("one.sided.fixed", list(c(.05), c(1, .5))),
    prior_weightfunction("two.sided.fixed", list(c(.05, 0.10), c(1, .2, .5)))
  )

  for(i in 1:length(priors)){
    model_syntax <- "model{}"
    model_syntax <- JAGS_add_priors(model_syntax, priors[i])
    monitor      <- JAGS_to_monitor(priors[i])
    inits        <- JAGS_get_inits(priors[i], chains = 2, seed = 1)

    set.seed(1)
    model   <- rjags::jags.model(file = textConnection(model_syntax), inits = inits, n.chains = 2, quiet = TRUE)
    samples <- rjags::coda.samples(model = model, variable.names = monitor, n.iter = 5000, quiet = TRUE, progress.bar = "none")
    samples <- do.call(rbind, samples)


    expect_doppelganger(paste0("JAGS-model-weightfunction-",i), function(){
      densities <- density(priors[[i]], individual = TRUE)
      par(mfrow = c(1, length(densities)))
      for(j in seq_along(densities)){
        hist(samples[,paste0("omega[",j,"]")], breaks = 50, freq = FALSE)
        lines(densities[[j]])
      }
    })
  }
})

