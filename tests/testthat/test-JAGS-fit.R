context("JAGS fit functions")

test_that("JAGS model functions work (simple)", {

  skip_if_not_installed("rjags")
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
    PEESE = prior_PEESE("gamma", list(1, 1)),
    p12 = prior("bernoulli", list(0.75))
  )

  model_syntax <- JAGS_add_priors(model_syntax, priors)
  monitor      <- JAGS_to_monitor(priors)
  inits        <- JAGS_get_inits(priors, chains = 2, seed = 1)

  set.seed(1)
  model   <- rjags::jags.model(file = textConnection(model_syntax), inits = inits, n.chains = 2, quiet = TRUE)
  samples <- rjags::coda.samples(model = model, variable.names = monitor, n.iter = 5000, quiet = TRUE, progress.bar = "none")
  samples <- do.call(rbind, samples)


  for(i in seq_along(priors)){
    vdiffr::expect_doppelganger(paste0("JAGS-model-prior-",i), function(){
      if(is.prior.discrete(priors[[i]])){
        barplot(table(samples[,names(priors)[i]])/length(samples[,names(priors)[i]]), main = print(priors[[i]], plot = T), width = 1/(max(samples[,names(priors)[i]])+1), space = 0, xlim = c(-0.25, max(samples[,names(priors)[i]])+0.25))
      }else{
        hist(samples[,names(priors)[i]], breaks = 50, main = print(priors[[i]], plot = TRUE), freq = FALSE)
      }
      lines(priors[[i]], individual = TRUE)
    })
  }
})

# skip the rest as it takes too long
skip_on_cran()

test_that("JAGS model functions work (vector)", {

  skip_if_not_installed("rjags")
  model_syntax <- "model{}"
  priors       <- list(
    p1  = prior("mnormal", list(mean = 0, sd = 1, K = 3),),
    p2  = prior("mcauchy", list(location = 0, scale = 1.5, K = 2)),
    p3  = prior("mt",      list(location = 2, scale = 0.5, df = 5, K = 2))
  )


  model_syntax <- JAGS_add_priors(model_syntax, priors)
  monitor      <- JAGS_to_monitor(priors)
  inits        <- JAGS_get_inits(priors, chains = 2, seed = 1)

  set.seed(1)
  model   <- rjags::jags.model(file = textConnection(model_syntax), inits = inits, n.chains = 2, quiet = TRUE)
  samples <- rjags::coda.samples(model = model, variable.names = monitor, n.iter = 5000, quiet = TRUE, progress.bar = "none")
  samples <- do.call(rbind, samples)

  vdiffr::expect_doppelganger("JAGS-model-prior-vector-1", function(){

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfcol = oldpar[["mfcol"]]))
    par(mfcol = c(1, 2))

    hist(samples[,"p1[1]"], breaks = 50, main = print(priors[[1]], plot = TRUE), freq = FALSE)
    lines(prior("normal", list(0, 1)))

    plot(samples[,"p1[1]"], samples[,"p1[2]"], pch = 19, xlim = c(-3, 3), ylim = c(-3, 3), asp = 1,
         xlab = "X1", ylab = "X2", main = print(priors[[1]], plot = TRUE))
  })

  vdiffr::expect_doppelganger("JAGS-model-prior-vector-2", function(){

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfcol = oldpar[["mfcol"]]))
    par(mfcol = c(1, 2))

    hist(samples[,"p2[1]"][abs(samples[,"p2[1]"]) < 5], breaks = 20, main = print(priors[[2]], plot = TRUE), freq = FALSE)
    lines(prior("cauchy", list(0, 1.5)))

    plot(samples[,"p2[1]"], samples[,"p2[2]"], pch = 19, xlim = c(-5, 5), ylim = c(-5, 5), asp = 1,
         xlab = "X1", ylab = "X2", main = print(priors[[2]], plot = TRUE))
  })

  vdiffr::expect_doppelganger("JAGS-model-prior-vector-3", function(){

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfcol = oldpar[["mfcol"]]))
    par(mfcol = c(1, 2))

    hist(samples[,"p3[1]"], breaks = 50, main = print(priors[[3]], plot = TRUE), freq = FALSE)
    lines(prior("t", list(2, 0.5, 5)))

    plot(samples[,"p3[1]"], samples[,"p3[2]"], pch = 19, xlim = c(-3, 7), ylim = c(-3, 7), asp = 1,
         xlab = "X1", ylab = "X2", main = print(priors[[3]], plot = TRUE))
  })

})

test_that("JAGS model functions work (factor)", {

  skip_if_not_installed("rjags")
  skip_on_os(c("mac", "linux", "solaris")) # multivariate sampling does not exactly match across OSes

  model_syntax <- "model{}"
  priors       <- list(
    p1  = prior_factor("mnorm",   list(mean = 0, sd = 1),    contrast = "orthonormal"),
    p2  = prior_factor("beta",    list(alpha = 1, beta = 1), contrast = "treatment"),
    p3  = prior_factor("beta",    list(alpha = 2, beta = 2), contrast = "treatment"),
    p4  = prior_factor("gamma",   list(shape = 2, rate = 3), contrast = "independent"),
    p5  = prior_factor("uniform", list(a = -0.5, b = 1.5),   contrast = "independent"),
    p6  = prior_factor("mnorm",   list(mean = 0, sd = 0.5),  contrast = "meandif"),
    p7  = prior_factor("spike",   list(location = 1), contrast = "treatment"),
    p8  = prior_factor("spike",   list(location = 2), contrast = "independent"),
    p9  = prior_factor("spike",   list(location = 0), contrast = "orthonormal"),
    p10 = prior_factor("spike",   list(location = 0), contrast = "meandif")
  )

  # add levels
  attr(priors[[1]],  "levels") <- 3
  attr(priors[[2]],  "levels") <- 2
  attr(priors[[3]],  "levels") <- 3
  attr(priors[[4]],  "levels") <- 1
  attr(priors[[5]],  "levels") <- 3
  attr(priors[[6]],  "levels") <- 3
  attr(priors[[7]],  "levels") <- 2
  attr(priors[[8]],  "levels") <- 3
  attr(priors[[9]],  "levels") <- 3
  attr(priors[[10]], "levels") <- 3


  model_syntax <- JAGS_add_priors(model_syntax, priors)
  monitor      <- JAGS_to_monitor(priors)
  inits        <- JAGS_get_inits(priors, chains = 2, seed = 1)

  set.seed(1)
  model   <- rjags::jags.model(file = textConnection(model_syntax), inits = inits, n.chains = 2, quiet = TRUE)
  samples <- rjags::coda.samples(model = model, variable.names = monitor, n.iter = 5000, quiet = TRUE, progress.bar = "none")
  samples <- do.call(rbind, samples)

  expect_equal(colnames(samples), c("p1[1]", "p1[2]", "p10[1]", "p10[2]", "p2", "p3[1]", "p3[2]", "p4", "p5[1]", "p5[2]", "p5[3]", "p6[1]", "p6[2]", "p7", "p8[1]", "p8[2]", "p8[3]", "p9[1]", "p9[2]"))

  vdiffr::expect_doppelganger("JAGS-model-prior-factor-1", function(){

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfcol = oldpar[["mfcol"]]))
    par(mfcol = c(1, 2))

    hist(samples[,"p1[1]"], breaks = 50, main = print(priors[[1]], plot = TRUE), freq = FALSE)
    lines(prior("normal", list(0, 1)))

    plot(samples[,"p1[1]"], samples[,"p1[2]"], pch = 19, xlim = c(-3, 3), ylim = c(-3, 3), asp = 1,
         xlab = "X1", ylab = "X2", main = print(priors[[1]], plot = TRUE))
  })

  vdiffr::expect_doppelganger("JAGS-model-prior-factor-2", function(){

    hist(samples[,"p2"], breaks = 20, main = print(priors[[2]], plot = TRUE), freq = FALSE)
    lines(prior("beta", list(1, 1)))
  })

  vdiffr::expect_doppelganger("JAGS-model-prior-factor-3", function(){

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfcol = oldpar[["mfcol"]]))
    par(mfcol = c(1, 2))

    hist(samples[,"p3[1]"], breaks = 50, main = print(priors[[3]], plot = TRUE), freq = FALSE)
    lines(prior("beta", list(2, 2)))

    plot(samples[,"p3[1]"], samples[,"p3[2]"], pch = 19, xlim = c(0, 1), ylim = c(0, 1), asp = 1,
         xlab = "X1", ylab = "X2", main = print(priors[[3]], plot = TRUE), cex = .25)
  })

  vdiffr::expect_doppelganger("JAGS-model-prior-factor-4", function(){

    hist(samples[,"p4"], breaks = 20, main = print(priors[[4]], plot = TRUE), freq = FALSE)
    lines(prior("gamma", list(shape = 2, rate = 3)))
  })

  vdiffr::expect_doppelganger("JAGS-model-prior-factor-5", function(){

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfcol = oldpar[["mfcol"]]))
    par(mfcol = c(1, 3))

    hist(samples[,"p5[1]"], breaks = 50, main = print(priors[[5]], plot = TRUE), freq = FALSE)
    lines(prior("uniform", list(a = -0.5, b = 1.5)))

    hist(samples[,"p5[2]"], breaks = 50, main = print(priors[[5]], plot = TRUE), freq = FALSE)
    lines(prior("uniform", list(a = -0.5, b = 1.5)))

    hist(samples[,"p5[3]"], breaks = 50, main = print(priors[[5]], plot = TRUE), freq = FALSE)
    lines(prior("uniform", list(a = -0.5, b = 1.5)))
  })

  vdiffr::expect_doppelganger("JAGS-model-prior-factor-6", function(){

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfcol = oldpar[["mfcol"]]))
    par(mfcol = c(1, 2))

    hist(samples[,"p6[1]"], breaks = 50, main = print(priors[[6]], plot = TRUE), freq = FALSE)
    lines(prior("normal", list(mean = 0, sd = 0.5)))

    hist(samples[,"p6[2]"], breaks = 50, main = print(priors[[6]], plot = TRUE), freq = FALSE)
    lines(prior("normal", list(mean = 0, sd = 0.5)))
  })

  vdiffr::expect_doppelganger("JAGS-model-prior-factor-7", function(){

    hist(samples[,"p7"], breaks = 50, main = print(priors[[7]], plot = TRUE), freq = FALSE)
    lines(prior_factor("spike", list(location = 1), contrast = "treatment"))

  })

  vdiffr::expect_doppelganger("JAGS-model-prior-factor-8", function(){

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfcol = oldpar[["mfcol"]]))
    par(mfcol = c(1, 3))

    hist(samples[,"p8[1]"], breaks = 50, main = print(priors[[5]], plot = TRUE), freq = FALSE)
    lines(prior_factor("spike", list(location = 2), contrast = "independent"))

    hist(samples[,"p8[2]"], breaks = 50, main = print(priors[[5]], plot = TRUE), freq = FALSE)
    lines(prior_factor("spike", list(location = 2), contrast = "independent"))

    hist(samples[,"p8[3]"], breaks = 50, main = print(priors[[5]], plot = TRUE), freq = FALSE)
    lines(prior_factor("spike", list(location = 2), contrast = "independent"))
  })

  vdiffr::expect_doppelganger("JAGS-model-prior-factor-9", function(){

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfcol = oldpar[["mfcol"]]))
    par(mfcol = c(1, 2))

    hist(samples[,"p9[1]"], breaks = 50, main = print(priors[[5]], plot = TRUE), freq = FALSE)
    lines(prior_factor("spike", list(location = 0), contrast = "orthonormal"))

    hist(samples[,"p9[2]"], breaks = 50, main = print(priors[[5]], plot = TRUE), freq = FALSE)
    lines(prior_factor("spike", list(location = 0), contrast = "orthonormal"))
  })

  vdiffr::expect_doppelganger("JAGS-model-prior-factor-10", function(){

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfcol = oldpar[["mfcol"]]))
    par(mfcol = c(1, 2))

    hist(samples[,"p10[1]"], breaks = 50, main = print(priors[[5]], plot = TRUE), freq = FALSE)
    lines(prior_factor("spike", list(location = 0), contrast = "meandif"))

    hist(samples[,"p10[2]"], breaks = 50, main = print(priors[[5]], plot = TRUE), freq = FALSE)
    lines(prior_factor("spike", list(location = 0), contrast = "meandif"))
  })

})

test_that("JAGS model functions work (weightfunctions)", {

  skip_if_not_installed("rjags")
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


    vdiffr::expect_doppelganger(paste0("JAGS-model-weightfunction-",i), function(){
      densities <- density(priors[[i]], individual = TRUE)
      oldpar <- graphics::par(no.readonly = TRUE)
      on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))
      par(mfrow = c(1, length(densities)))
      for(j in seq_along(densities)){
        hist(samples[,paste0("omega[",j,"]")], breaks = 50, freq = FALSE)
        lines(densities[[j]])
      }
    })
  }
})

test_that("JAGS model functions work (spike and slab)", {

  skip_if_not_installed("rjags")
  priors       <- list(
    "mu"   = prior_spike_and_slab(prior("normal", list(0, 1)), prior_inclusion = prior("beta", list(1,1))),
    "beta" = prior_spike_and_slab(prior_factor("mnormal", list(0, 1), contrast = "orthonormal"),
                                  prior_inclusion = prior("beta", list(1,1)))
  )

  attr(priors$beta$variable, "levels") <- 3

  for(i in 1:length(priors)){
    model_syntax <- "model{}"
    model_syntax <- JAGS_add_priors(model_syntax, priors[i])
    monitor      <- JAGS_to_monitor(priors[i])
    inits        <- JAGS_get_inits(priors[i], chains = 2, seed = 1)

    set.seed(1)
    model   <- rjags::jags.model(file = textConnection(model_syntax), inits = inits, n.chains = 2, quiet = TRUE)
    samples <- rjags::coda.samples(model = model, variable.names = monitor, n.iter = 5000, quiet = TRUE, progress.bar = "none")
    samples <- do.call(rbind, samples)


    if(i == 1){
      vdiffr::expect_doppelganger(paste0("JAGS-model-prior_spike-and-slab-",i), function(){
        temp_samples <- samples[,names(priors)[i]]
        hs <- hist(temp_samples[temp_samples != 0], breaks = 50, plot = FALSE)
        hs$density <- hs$density * mean(temp_samples != 0)
        plot(hs, main = print(priors[[i]], plot = TRUE), freq = FALSE, ylim = range(c(0, max(hs$density), mean(temp_samples == 0))))
        lines(priors[[i]], individual = TRUE)
      })
    }else{
      vdiffr::expect_doppelganger(paste0("JAGS-model-prior_spike-and-slab-",i), function(){

        oldpar <- graphics::par(no.readonly = TRUE)
        on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))
        par(mfrow = c(1, 3))

        temp_samples <- samples[,paste0(names(priors)[i], "[", 1:2, "]")]  %*% t(contr.orthonormal(1:3))

        hs1 <- hist(temp_samples[temp_samples[,1] != 0, 1], breaks = 50, plot = FALSE)
        hs1$density <- hs1$density * mean(temp_samples[,1] != 0)
        plot(hs1, main = print(priors[[i]], plot = TRUE), freq = FALSE, ylim = range(c(0, max(hs1$density), mean(temp_samples == 0))))
        lines(priors[[i]], individual = TRUE)

        hs2 <- hist(temp_samples[temp_samples[,2] != 0, 2], breaks = 50, plot = FALSE)
        hs2$density <- hs2$density * mean(temp_samples[,1] != 0)
        plot(hs2, main = print(priors[[i]], plot = TRUE), freq = FALSE, ylim = range(c(0, max(hs2$density), mean(temp_samples == 0))))
        lines(priors[[i]], individual = TRUE)

        hs3 <- hist(temp_samples[temp_samples[,3] != 0, 3], breaks = 50, plot = FALSE)
        hs3$density <- hs3$density * mean(temp_samples[,1] != 0)
        plot(hs3, main = print(priors[[i]], plot = TRUE), freq = FALSE, ylim = range(c(0, max(hs3$density), mean(temp_samples == 0))))
        lines(priors[[i]], individual = TRUE)
      })
    }
  }
})

test_that("JAGS model functions work (mixture)", {

  skip_if_not_installed("rjags")
  priors       <- list(
    "mu" = prior_mixture(
      list(
        prior("normal", list(0,  1), prior_weights = 1),
        prior("normal", list(-3, 1), prior_weights = 5),
        prior("gamma",  list(5, 10), prior_weights = 1)
      ),
      is_null = c(T, F, T)
    ),
    "beta" = prior_mixture(
      list(
        prior("normal", list(0,  1), prior_weights = 1),
        prior("normal", list(-3, 1), prior_weights = 5)
      ),
      components = c("b", "a")
    ),
    "gamma" = prior_mixture(
      list(
        prior("spike", list(2)),
        prior("normal", list(-3, 1))
      )
    ),
    "bias" = prior_mixture(list(
      prior_weightfunction(distribution = "two.sided", parameters = list(alpha = c(1, 1), steps = c(0.05)), prior_weights = 1/12),
      prior_weightfunction(distribution = "two.sided", parameters = list(alpha = c(1, 1, 1), steps = c(0.05, 0.10)), prior_weights = 1/12),
      prior_weightfunction(distribution = "one.sided", parameters = list(alpha = c(1, 1), steps = c(0.05)), prior_weights = 1/12),
      prior_weightfunction(distribution = "one.sided", parameters = list(alpha = c(1, 1, 1), steps = c(0.025, 0.05)), prior_weights = 1/12),
      prior_weightfunction(distribution = "one.sided", parameters = list(alpha = c(1, 1, 1), steps = c(0.05, 0.5)), prior_weights = 1/12),
      prior_weightfunction(distribution = "one.sided", parameters = list(alpha = c(1, 1, 1, 1), steps = c(0.025, 0.05, 0.5)), prior_weights = 1/12),
      prior_PET(distribution = "Cauchy", parameters = list(0,1), truncation = list(0, Inf), prior_weights = 1/4),
      prior_PEESE(distribution = "Cauchy", parameters = list(0,5), truncation = list(0, Inf), prior_weights = 1/4)
    ))
  )



  for(i in 1:length(priors)){
    model_syntax <- "model{}"
    model_syntax <- JAGS_add_priors(model_syntax, priors[i])
    monitor      <- JAGS_to_monitor(priors[i])
    inits        <- JAGS_get_inits(priors[i], chains = 2, seed = 1)

    if(i == 4){
      if("RoBMA" %in% rownames(installed.packages())){
        require("RoBMA")
      }else{
        next
      }
    }

    set.seed(1)
    model   <- rjags::jags.model(file = textConnection(model_syntax), inits = inits, n.chains = 2, quiet = TRUE)
    samples <- rjags::coda.samples(model = model, variable.names = monitor, n.iter = 5000, quiet = TRUE, progress.bar = "none")
    samples <- do.call(rbind, samples)


    if(i != 4){
      vdiffr::expect_doppelganger(paste0("JAGS-model-prior_mixture-",i), function(){
        temp_samples <- samples[,names(priors)[i]]
        hist(temp_samples, breaks = 100, freq = FALSE, main = print(priors[[i]], plot = TRUE))
        lines(density(rng(priors[[i]], 1000000)))
      })
    }else{
      vdiffr::expect_doppelganger(paste0("JAGS-model-prior_mixture-",i), function(){

        oldpar <- graphics::par(no.readonly = TRUE)
        on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))
        par(mfrow = c(3, 3))

        samples_PET   <- samples[,"PET"]
        samples_PEESE <- samples[,"PEESE"]
        samples_omega <- samples[,paste0("omega[",1:6,"]")]
        samples_bias  <- samples[,"bias_indicator"]

        barplot(table(samples_bias)/length(samples_bias), main = "Bias componenets")

        hist(samples_PET[samples_PET != 0 & samples_PET < 10], breaks = 50, main = "PET", freq = FALSE)
        lines(priors$bias[[7]], individual = TRUE)

        hist(samples_PEESE[samples_PEESE != 0 & samples_PEESE < 20], breaks = 50, main = "PEESE", freq = FALSE)
        lines(priors$bias[[8]], individual = TRUE)

        hist(samples_omega[samples_bias == 2, 1], breaks = 50, main = "omega[2:1]", freq = FALSE)
        hist(samples_omega[samples_bias == 2, 2], breaks = 50, main = "omega[2:2]", freq = FALSE)
        hist(samples_omega[samples_bias == 2, 3], breaks = 50, main = "omega[2:3]", freq = FALSE)
        hist(samples_omega[samples_bias == 2, 4], breaks = 50, main = "omega[2:4]", freq = FALSE)
        hist(samples_omega[samples_bias == 2, 5], breaks = 50, main = "omega[2:5]", freq = FALSE)
        hist(samples_omega[samples_bias == 2, 6], breaks = 50, main = "omega[2:6]", freq = FALSE)

        })
    }
  }
})

test_that("JAGS fit function works" , {

  set.seed(1)
  data <- list(
    x = rnorm(50, 0, .5),
    N = 50
  )
  priors_list <- list(
    m  = prior("normal", list(0, 1)),
    s  = prior("normal", list(0, 1), list(0, Inf))
  )
  model_syntax <-
    "model
    {
      for(i in 1:N){
        x[i] ~ dnorm(m, pow(s, -2))
      }
    }"

  ### checking the default settings
  set.seed(1)
  fit <- JAGS_fit(model_syntax, data, priors_list)
  expect_equal(length(fit$mcmc), 4)
  expect_true(all(sapply(fit$mcmc, function(mcmc)dim(mcmc) == c(4000, 2))))
  vdiffr::expect_doppelganger("JAGS-fit-posterior", function(){
    samples <- do.call(rbind, fit$mcmc)
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))
    par(mfrow = c(1, 2))
    for(i in seq_along(priors_list)){
      hist(samples[,i], breaks = 50, freq = FALSE)
    }
  })

  ### checking control the main control arguments
  fit1 <- JAGS_fit(model_syntax, data, priors_list, chains = 1, adapt = 100, burnin = 150, sample = 175, thin = 3, seed = 2)
  expect_equal(length(fit1$mcmc), 1)
  expect_true(all(sapply(fit1$mcmc, function(mcmc)dim(mcmc) == c(175, 2))))
  expect_equal(fit1$burnin, 250) # adapt + burnin
  expect_equal(fit1$sample, 175)
  expect_equal(fit1$thin, 3)

  ### adding custom parameters
  model_syntax3 <-
    "model
    {
      g ~ dnorm(0, 1)
      for(i in 1:N){
        x[i] ~ dnorm(m, pow(s, -2))
      }
    }"
  fit3 <- JAGS_fit(model_syntax3, data, priors_list, add_parameters = "g",
                   chains = 1, adapt = 100, burnin = 100, sample = 100, seed = 3)
  expect_equal(colnames(fit3$mcmc[[1]]), c("m", "s", "g"))

  ### checking mcmc_error autofit
  priors_list4 <- list(
    m  = prior("normal", list(0, 1))
  )
  data4 <- list(
    x = c(-500),
    N = 1
  )
  model_syntax4 <-
    "model
    {
      l = 1
      for(i in 1:N){
        x[i] ~ dt(m, pow(.3, -2), 1)
      }
    }"
  runjags::runjags.options(silent.jags = TRUE, silent.runjags = TRUE)
  fit4 <- JAGS_fit(model_syntax4, data4, priors_list4, autofit = FALSE,
                   chains = 2, adapt = 100, burnin = 50, sample = 100, seed = 4)
  summary_4 <- suppressWarnings(summary(fit4))
  expect_true(summary_4[1,"MCerr"]   > 0.069)
  expect_true(summary_4[1,"MC%ofSD"] > 8)
  expect_true(summary_4[1,"SSeff"]   < 150)
  expect_true(summary_4[1,"psrf"]    > 1.007)

  convergence <- JAGS_check_convergence(fit4, prior_list = priors_list4, max_Rhat = 1.001, min_ESS = 500, max_error = 0.01, max_SD_error = 0.05)
  expect_true(!convergence)
  expect_equal(attr(convergence, "errors")[1], "R-hat 1.053 is larger than the set target (1.001).")
  expect_equal(attr(convergence, "errors")[2], "ESS 149 is lower than the set target (500).")
  expect_equal(attr(convergence, "errors")[3], "MCMC error 0.07422 is larger than the set target (0.01).")
  expect_equal(attr(convergence, "errors")[4], "MCMC SD error 0.087 is larger than the set target (0.05).")

  fit4b <- JAGS_fit(model_syntax4, data4, priors_list4, autofit = TRUE, autofit_control = list(max_error = 0.05, sample_extend = 100),
                    chains = 2, adapt = 100, burnin = 50, sample = 100, seed = 4)
  summary_4b <- summary(fit4b)
  expect_true(summary_4b[1,"MCerr"]   < 0.05)

  fit4c <- JAGS_fit(model_syntax4, data4, priors_list4, autofit = TRUE, autofit_control = list(max_Rhat = 1.001, sample_extend = 100),
                    chains = 2, adapt = 100, burnin = 50, sample = 100, seed = 4)
  summary_4c <- summary(fit4c)
  expect_true(summary_4c[1,"psrf"]    < 1.001)

  fit4d <- JAGS_fit(model_syntax4, data4, priors_list4, autofit = TRUE, autofit_control = list(min_ESS = 200, sample_extend = 100),
                    chains = 2, adapt = 100, burnin = 50, sample = 100, seed = 4)
  summary_4d <- summary(fit4d)
  expect_true(summary_4d[1,"SSeff"]   > 200)

  fit4e <- JAGS_fit(model_syntax4, data4, priors_list4, autofit = TRUE, autofit_control = list(max_SD_error = 0.06, sample_extend = 100),
                    chains = 2, adapt = 100, burnin = 50, sample = 100, seed = 4)
  summary_4e <- summary(fit4e)
  expect_true(summary_4e[1,"MC%ofSD"] < 6)

  fit4f <- JAGS_fit(model_syntax4, data4, priors_list4, autofit = TRUE, autofit_control = list(max_error = 0.0001, sample_extend = 100, max_time = list(time = 5, unit = "secs")),
                    chains = 2, adapt = 100, burnin = 50, sample = 100, seed = 4)
  summary_4f <- summary(fit4f)
  expect_true(summary_4f[1,"MCerr"] > 0.0001)
  expect_true(fit4f$timetaken < 5)

  # test extending the fit
  fite <- JAGS_extend(fit)
  expect_equal(length(fite$mcmc), 4)
  expect_true(all(sapply(fite$mcmc, function(mcmc)dim(mcmc) == c(5000, 2))))
})

test_that("JAGS fit function integration with formula works" , {

  skip_on_os(c("mac", "linux", "solaris")) # multivariate sampling does not exactly match across OSes
  skip_on_cran()

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

  fit1 <- JAGS_fit(
    model_syntax = model_syntax1, data = data, prior_list = prior_list1,
    formula_list = formula_list1, formula_data_list = formula_data_list1, formula_prior_list = formula_prior_list1)

  posterior1 <- suppressWarnings(coda::as.mcmc(fit1))

  lm_1 <- stats::lm(y ~ x_cont1 + x_fac3t, data = cbind(data_formula, y = data$y))

  # verify against the frequentist fit
  vdiffr::expect_doppelganger("JAGS-fit-formula-1", function(){

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfcol = oldpar[["mfcol"]]))
    par(mfrow = c(2, 2))

    hist(posterior1[,"mu_intercept"], freq = FALSE, main = "Intercept")
    curve(dnorm(x, mean = coef(lm_1)["(Intercept)"], sd = summary(lm_1)$coefficients["(Intercept)", "Std. Error"]), add = TRUE, lwd = 2)

    hist(posterior1[,"mu_x_cont1"], freq = FALSE, main = "mu_x_cont1")
    curve(dnorm(x, mean = coef(lm_1)["x_cont1"], sd = summary(lm_1)$coefficients["x_cont1", "Std. Error"]), add = TRUE, lwd = 2)

    hist(posterior1[,"mu_x_fac3t[1]"], freq = FALSE, main = "mu_x_fac3t")
    curve(dnorm(x, mean = coef(lm_1)["x_fac3tB"], sd = summary(lm_1)$coefficients["x_fac3tB", "Std. Error"]), add = TRUE, lwd = 2)

    hist(posterior1[,"mu_x_fac3t[2]"], freq = FALSE, main = "mu_x_fac3t")
    curve(dnorm(x, mean = coef(lm_1)["x_fac3tC"], sd = summary(lm_1)$coefficients["x_fac3tC", "Std. Error"]), add = TRUE, lwd = 2)
  })

  # create model with two formulas
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

  posterior2 <- suppressWarnings(coda::as.mcmc(fit2))

  # verify against the true values
  vdiffr::expect_doppelganger("JAGS-fit-formula-2", function(){

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfcol = oldpar[["mfcol"]]))
    par(mfrow = c(2, 3))

    hist(posterior2[,"mu_intercept"], freq = FALSE, main = "Intercept")
    abline(v = 0, lwd = 3, col = "blue")

    hist(posterior2[,"mu_x_cont1"], freq = FALSE, main = "mu_x_cont1")
    abline(v = .4, lwd = 3, col = "blue")

    hist(posterior2[,"mu_x_fac3t[1]"], freq = FALSE, main = "mu_x_fac3t")
    abline(v = -0.2, lwd = 3, col = "blue")

    hist(posterior2[,"mu_x_fac3t[2]"], freq = FALSE, main = "mu_x_fac3t")
    abline(v = 0.4, lwd = 3, col = "blue")

    hist(exp(posterior2[,"sigma_intercept"]), freq = FALSE, main = "sigma_intercept")
    abline(v = 0.5, lwd = 3, col = "blue")

    hist(exp(posterior2[,"sigma_intercept"] + posterior2[,"sigma_x_fac2t"]), freq = FALSE, main = "sigma_x_fac2t")
    abline(v = 1, lwd = 3, col = "blue")
  })

})

test_that("JAGS fit function integration with formula and spike and slab works" , {

  skip_on_os(c("mac", "linux", "solaris")) # multivariate sampling does not exactly match across OSes
  skip_on_cran()

  set.seed(1)

  data_formula <- data.frame(
    x_cont1 = rnorm(300),
    x_fac2t = factor(rep(c("A", "B"), 150), levels = c("A", "B")),
    x_fac3t = factor(rep(c("A", "B", "C"), 100), levels = c("A", "B", "C"))
  )
  data <- list(
    y = rnorm(300, .20 * data_formula$x_cont1 + ifelse(data_formula$x_fac3t == "A", 0.0, ifelse(data_formula$x_fac3t == "B", -0.2, 0.2)), ifelse(data_formula$x_fac2t == "A", 0.5, 1)),
    N = 300
  )


  # create model with mix of a formula and free parameters ---
  formula_list1 <- list(
    mu    = ~ x_cont1 + x_fac3t
  )
  formula_data_list1 <- list(
    mu    = data_formula
  )
  formula_prior_list1 <- list(
    mu    = list(
      "intercept"  = prior("normal", list(0, 5)),
      "x_cont1"    = prior_spike_and_slab(prior("normal", list(0, 1)), prior_inclusion = prior("beta", list(1,1))),
      "x_fac3t"    = prior_spike_and_slab(prior_factor("mnormal", list(0, 1), contrast = "orthonormal"),
                                          prior_inclusion = prior("spike", list(0.5)))
    )
  )
  attr(formula_prior_list1$mu$x_fac3t$variable, "multiply_by") <- "sigma"
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

  fit1 <- JAGS_fit(
    model_syntax = model_syntax1, data = data, prior_list = prior_list1,
    formula_list = formula_list1, formula_data_list = formula_data_list1, formula_prior_list = formula_prior_list1)

  posterior1 <- suppressWarnings(coda::as.mcmc(fit1))

  vdiffr::expect_doppelganger("JAGS-fit-formula-spike-and-slab-1", function(){

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfcol = oldpar[["mfcol"]]))
    par(mfrow = c(1, 3))

    hist(posterior1[,"mu_x_cont1"], freq = FALSE, main = "x_cont1")

    hist(posterior1[,"mu_x_cont1_variable"], freq = FALSE, main = "x_cont1_variable")

    hist(posterior1[,"mu_x_cont1_inclusion"], freq = FALSE, main = "x_cont1_inclusion")
  })

  vdiffr::expect_doppelganger("JAGS-fit-formula-spike-and-slab-2", function(){

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfcol = oldpar[["mfcol"]]))
    par(mfrow = c(2, 3))

    temp_samples          <- posterior1[,paste0("mu_x_fac3t[", 1:2, "]")]          %*% t(contr.orthonormal(1:3))
    temp_samples_variable <- posterior1[,paste0("mu_x_fac3t_variable[", 1:2, "]")]  %*% t(contr.orthonormal(1:3))

    hist(temp_samples[,1], freq = FALSE, main = "x_fac3t[A]")
    hist(temp_samples[,2], freq = FALSE, main = "x_fac3t[B]")
    hist(temp_samples[,3], freq = FALSE, main = "x_fac3t[C]")

    hist(temp_samples_variable[,1], freq = FALSE, main = "x_fac3t_variable[A]")
    hist(temp_samples_variable[,2], freq = FALSE, main = "x_fac3t_variable[B]")
    hist(temp_samples_variable[,3], freq = FALSE, main = "x_fac3t_variable[C]")
  })

})

test_that("JAGS fit function integration with formula, spike and slab works, and mixture works" , {

  skip_on_os(c("mac", "linux", "solaris")) # multivariate sampling does not exactly match across OSes
  skip_on_cran()

  set.seed(1)

  data_formula <- data.frame(
    x_cont1 = rnorm(300),
    x_fac2t = factor(rep(c("A", "B"), 150), levels = c("A", "B")),
    x_fac3t = factor(rep(c("A", "B", "C"), 100), levels = c("A", "B", "C"))
  )
  data <- list(
    y = rnorm(300, -0.15 + 0.20 * data_formula$x_cont1 + ifelse(data_formula$x_fac3t == "A", 0.0, ifelse(data_formula$x_fac3t == "B", -0.2, 0.2)), ifelse(data_formula$x_fac2t == "A", 0.5, 1)),
    N = 300
  )

  # create model with mix of a formula and free parameters ---
  formula_list1 <- list(
    mu    = ~ x_cont1 + x_fac3t
  )
  formula_data_list1 <- list(
    mu    = data_formula
  )
  formula_prior_list1 <- list(
    mu    = list(
      "intercept"  = prior_mixture(
        list(
          prior("spike",   list(0),       prior_weights = 2),
          prior("normal",  list(-1, 0.5), prior_weights = 1),
          prior("normal",  list( 1, 0.5), prior_weights = 1)
        ),
        is_null = c(T, F, F)
      ),
      "x_cont1"    = prior_mixture(
        list(
          prior("spike",   list(0),    prior_weights = 1),
          prior("normal",  list(0, 1), prior_weights = 1)
        ),
        is_null = c(T, F)
      ),
      "x_fac3t"    = prior_spike_and_slab(prior_factor("mnormal", list(0, 1), contrast = "orthonormal"),
                                          prior_inclusion = prior("spike", list(0.5)))
    )
  )
  attr(formula_prior_list1$mu$x_cont1, "multiply_by") <- "sigma"
  prior_list1 <- list(
    "sigma" = prior_mixture(
      list(
        prior("normal",    list(0, 1), truncation = list(0, Inf)),
        prior("lognormal", list(0, 1))
      ),
      is_null = c(T, F)
    )
  )
  model_syntax1 <- paste0(
    "model{\n",
    "for(i in 1:N){\n",
    "  y[i] ~ dnorm(mu[i], 1/pow(sigma, 2))\n",
    "}\n",
    "}"
  )

  fit1 <- JAGS_fit(
    model_syntax = model_syntax1, data = data, prior_list = prior_list1,
    formula_list = formula_list1, formula_data_list = formula_data_list1, formula_prior_list = formula_prior_list1)

  posterior1 <- suppressWarnings(coda::as.mcmc(fit1))

  vdiffr::expect_doppelganger("JAGS-fit-formula-mixture-1", function(){

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfcol = oldpar[["mfcol"]]))
    par(mfrow = c(3, 3))

    barplot(table(posterior1[,"mu_intercept_indicator"]) / nrow(posterior1), main = "Intercept indicator")
    barplot(table(posterior1[,"mu_x_cont1_indicator"]) / nrow(posterior1), main = "x_cont1 indicator")
    barplot(table(posterior1[,"sigma_indicator"]) / nrow(posterior1), main = "sigma indicator")

    hist(posterior1[,"mu_intercept"], freq = FALSE, main = "mu_intercept")
    hist(posterior1[,"mu_x_cont1"], freq = FALSE, main = "x_cont1")
    hist(posterior1[,"sigma"], freq = FALSE, main = "sigma")
  })

  vdiffr::expect_doppelganger("JAGS-fit-formula-mixture-2", function(){

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfcol = oldpar[["mfcol"]]))
    par(mfrow = c(2, 3))

    temp_samples          <- posterior1[,paste0("mu_x_fac3t[", 1:2, "]")]          %*% t(contr.orthonormal(1:3))
    temp_samples_variable <- posterior1[,paste0("mu_x_fac3t_variable[", 1:2, "]")]  %*% t(contr.orthonormal(1:3))

    hist(temp_samples[,1], freq = FALSE, main = "x_fac3t[A]")
    hist(temp_samples[,2], freq = FALSE, main = "x_fac3t[B]")
    hist(temp_samples[,3], freq = FALSE, main = "x_fac3t[C]")

    hist(temp_samples_variable[,1], freq = FALSE, main = "x_fac3t_variable[A]")
    hist(temp_samples_variable[,2], freq = FALSE, main = "x_fac3t_variable[B]")
    hist(temp_samples_variable[,3], freq = FALSE, main = "x_fac3t_variable[C]")
  })

})

test_that("JAGS fit with priors expressions", {

  skip_if_not_installed("rjags")

  # a simple prior
  model_syntax <- "model{}"
  priors       <- list(
    x        = prior("normal",   list(0, expression(x_sigma))),
    x_sigma  = prior("invgamma", list(1/2, 1/2))
  )

  model_syntax <- JAGS_add_priors(model_syntax, priors)
  monitor      <- JAGS_to_monitor(priors)
  inits        <- JAGS_get_inits(priors, chains = 2, seed = 1)

  set.seed(1)
  model   <- rjags::jags.model(file = textConnection(model_syntax), inits = inits, n.chains = 2, quiet = TRUE)
  samples <- rjags::coda.samples(model = model, variable.names = monitor, n.iter = 5000, quiet = TRUE, progress.bar = "none")
  samples <- do.call(rbind, samples)

  vdiffr::expect_doppelganger("JAGS-model-prior-e1", function(){
    x_samples <- samples[,"x"]
    hist(x_samples[abs(x_samples) < 10], breaks = 50, main = print(priors[[1]], plot = TRUE), freq = FALSE)
    lines(prior("Cauchy", list(0, 1), list(-10, 10)))
  })

  # a spike and slab prior
  model_syntax <- "model{}"
  priors       <- list(
    x        = prior_spike_and_slab(
      prior("normal",   list(0, expression(x_sigma)))
    ),
    x_sigma  = prior("invgamma", list(1/2, 1/2))
  )

  model_syntax <- JAGS_add_priors(model_syntax, priors)
  monitor      <- JAGS_to_monitor(priors)
  inits        <- JAGS_get_inits(priors, chains = 2, seed = 1)

  set.seed(1)
  model   <- rjags::jags.model(file = textConnection(model_syntax), inits = inits, n.chains = 2, quiet = TRUE)
  samples <- rjags::coda.samples(model = model, variable.names = monitor, n.iter = 5000, quiet = TRUE, progress.bar = "none")
  samples <- do.call(rbind, samples)

  vdiffr::expect_doppelganger("JAGS-model-prior-e2", function(){
    x_samples <- samples[,"x"]
    x_samples <- x_samples[x_samples != 0]
    hist(x_samples[abs(x_samples) < 10], breaks = 50, main = print(priors[[1]], plot = TRUE), freq = FALSE)
    lines(prior("Cauchy", list(0, 1), list(-10, 10)))
  })

  # a mixture prior
  model_syntax <- "model{}"
  priors       <- list(
    x        = prior_mixture(list(
      prior("normal",   list(0, expression(x_sigma))),
      prior("cauchy",   list(0, 1))
    ), is_null = c(T, F)),
    x_sigma  = prior("invgamma", list(1/2, 1/2))
  )

  model_syntax <- JAGS_add_priors(model_syntax, priors)
  monitor      <- JAGS_to_monitor(priors)
  inits        <- JAGS_get_inits(priors, chains = 2, seed = 1)

  set.seed(1)
  model   <- rjags::jags.model(file = textConnection(model_syntax), inits = inits, n.chains = 2, quiet = TRUE)
  samples <- rjags::coda.samples(model = model, variable.names = monitor, n.iter = 5000, quiet = TRUE, progress.bar = "none")
  samples <- do.call(rbind, samples)

  vdiffr::expect_doppelganger("JAGS-model-prior-e3", function(){
    x_samples <- samples[,"x"]
    hist(x_samples[abs(x_samples) < 10], breaks = 50, main = print(priors[[1]], plot = TRUE), freq = FALSE)
    lines(prior("Cauchy", list(0, 1), list(-10, 10)))
  })

})

test_that("JAGS fit function integration with formula and priors expressions", {

  skip_on_os(c("mac", "linux", "solaris")) # multivariate sampling does not exactly match across OSes
  skip_if_not_installed("rjags")
  skip_on_cran()

  set.seed(1)

  data_formula <- data.frame(
    x_cont1 = rnorm(100),
    x_fac2t = factor(rep(c("A", "B"), 50), levels = c("A", "B"))
  )
  data <- list(
    y = rnorm(100, .4 * data_formula$x_cont1 + ifelse(data_formula$x_fac2t == "A", 0.25, 0.50)),
    N = 100
  )

  # create model with mix of a formula and free parameters ---
  formula_list <- list(
    mu    = ~ x_cont1 + x_fac2t
  )
  formula_data_list <- list(
    mu    = data_formula
  )
  formula_prior_list1 <- list(
    mu    = list(
      "intercept"   = prior("normal", list(0, 5)),
      "x_cont1"     = prior("normal", list(0, 1)),
      "x_fac2t"     = prior_spike_and_slab(prior_factor("cauchy", contrast = "treatment", list(0, 1)))
    )
  )
  formula_prior_list2 <- list(
    mu    = list(
      "intercept"   = prior("normal", list(0, 5)),
      "x_cont1"     = prior("normal", list(0, 1)),
      "x_fac2t"     = prior_spike_and_slab(prior_factor("normal", contrast = "treatment", list(0, expression(tau))))
    )
  )
  prior_list1 <- list(
    sigma = prior("lognormal", list(0, 1))
  )
  prior_list2 <- list(
    sigma = prior("lognormal", list(0, 1)),
    tau   = prior("invgamma", list(1/2, 1/2))
  )

  model_syntax <- paste0(
    "model{\n",
    "for(i in 1:N){\n",
    "  y[i] ~ dnorm(mu[i], 1/pow(sigma, 2))\n",
    "}\n",
    "}"
  )

  # fit1 <- JAGS_fit(
  #   model_syntax = model_syntax, data = data, prior_list = prior_list1,
  #   formula_list = formula_list, formula_data_list = formula_data_list, formula_prior_list = formula_prior_list1, adapt = 1000)
  fit2 <- JAGS_fit(
    model_syntax = model_syntax, data = data, prior_list = prior_list2,
    formula_list = formula_list, formula_data_list = formula_data_list, formula_prior_list = formula_prior_list2, adapt = 1000)

  # runjags_estimates_table(fit1)
  # verified against the simpler model directly sampling cauchy
  expect_equal(capture_output_lines(print(runjags_estimates_table(fit2, conditional = FALSE, remove_parameters = "tau")), width = 150),  c(
    "                          Mean    SD    lCI Median   uCI error(MCMC) error(MCMC)/SD   ESS R-hat",
    "(mu) intercept           0.304 0.125  0.024  0.312 0.527     0.00197          0.016  4128 1.001",
    "(mu) x_cont1             0.385 0.111  0.169  0.383 0.604     0.00091          0.008 14631 1.000",
    "(mu) x_fac2t (inclusion) 0.280    NA     NA     NA    NA          NA             NA    NA    NA",
    "(mu) x_fac2t[B]          0.068 0.148 -0.002  0.000 0.494     0.00309          0.021  2299 1.002",
    "sigma                    0.980 0.071  0.854  0.975 1.133     0.00074          0.010  9500 1.001"
  ))

})

test_that("JAGS fit function integration with formula expressions", {

  skip_on_os(c("mac", "linux", "solaris")) # multivariate sampling does not exactly match across OSes
  skip_if_not_installed("rjags")
  skip_on_cran()

  set.seed(1)

  data_formula <- data.frame(
    x_cont1 = rnorm(200),
    x_fac2t = factor(rep(LETTERS[1:10], 20))
  )
  x_fac2t_values        <- rnorm(10, 0, 0.3)
  names(x_fac2t_values) <- LETTERS[1:10]
  data <- list(
    y = rnorm(200, .4 * data_formula$x_cont1 + x_fac2t_values[data_formula[["x_fac2t"]]]),
    N = 200
  )

  # add id mapping
  data[["mapping_id"]] <- as.integer(data_formula$x_fac2t)

  # create model with mix of a formula and free parameters ---
  formula_list <- list(
    mu    = ~ x_cont1 + expression(mu_id[mapping_id[i]])
  )
  formula_data_list <- list(
    mu    = data_formula
  )
  formula_prior_list <- list(
    mu    = list(
      "intercept"   = prior("normal", list(0, 5)),
      "x_cont1"     = prior("normal", list(0, 1))
    )
  )
  prior_list <- list(
    sigma = prior("lognormal", list(0, 1)),
    mu_id = prior_factor("normal", list(0, expression(tau)), contrast = "independent"),
    tau   = prior("normal", list(0, 10), list(0, Inf))
  )

  attr(prior_list$mu_id, "levels") <- 10

  model_syntax <- paste0(
    "model{\n",
    "for(i in 1:N){\n",
    "  y[i] ~ dnorm(mu[i], 1/pow(sigma, 2))\n",
    "}\n",
    "}"
  )

  fit <- JAGS_fit(
    model_syntax = model_syntax, data = data, prior_list = prior_list,
    formula_list = formula_list, formula_data_list = formula_data_list, formula_prior_list = formula_prior_list)

  fit_sumary <- runjags_estimates_table(fit, formula_prefix = FALSE)
  expect_true(cor(fit_sumary[grepl("id", rownames(fit_sumary)),"Mean"], x_fac2t_values) > 0.8)
  expect_equal(capture_output_lines(print(fit_sumary), width = 150),  c(
    "            Mean    SD    lCI Median   uCI error(MCMC) error(MCMC)/SD   ESS R-hat",
    "intercept  0.179 0.173 -0.169  0.179 0.525     0.00476          0.028  1319 1.003",
    "x_cont1    0.423 0.077  0.272  0.424 0.575     0.00064          0.008 14330 1.000",
    "sigma      0.996 0.051  0.903  0.994 1.103     0.00055          0.011  8587 1.000",
    "id[1]      0.107 0.246 -0.376  0.106 0.602     0.00470          0.019  2763 1.002",
    "id[2]      0.359 0.250 -0.118  0.351 0.867     0.00479          0.019  2724 1.003",
    "id[3]      0.213 0.250 -0.272  0.207 0.726     0.00466          0.019  2873 1.001",
    "id[4]     -0.340 0.249 -0.853 -0.331 0.127     0.00481          0.019  2699 1.001",
    "id[5]     -0.476 0.257 -1.001 -0.468 0.009     0.00492          0.019  2804 1.001",
    "id[6]      0.643 0.260  0.152  0.634 1.180     0.00507          0.019  2692 1.002",
    "id[7]     -0.213 0.247 -0.723 -0.208 0.260     0.00462          0.019  2882 1.001",
    "id[8]     -0.069 0.245 -0.551 -0.069 0.413     0.00465          0.019  2781 1.001",
    "id[9]     -0.362 0.249 -0.873 -0.354 0.107     0.00489          0.020  2656 1.001",
    "id[10]     0.125 0.245 -0.347  0.122 0.624     0.00446          0.018  3012 1.002",
    "tau        0.474 0.172  0.229  0.444 0.903     0.00319          0.019  2915 1.001"
  ))

})

test_that("JAGS parallel fit function works", {

  skip("requires parallel processing")
  skip_on_cran()
  skip_on_travis()
  skip_on_ci()

  priors_list <- list(
    m  = prior("normal", list(0, 1))
  )
  data <- list(
    x = c(-500),
    N = 1
  )
  model_syntax <-
    "model
    {
      l = 1
      for(i in 1:N){
        x[i] ~ dt(m, pow(.3, -2), 1)
      }
    }"


  fit <- JAGS_fit(model_syntax, data, priors_list, autofit = FALSE, parallel = TRUE,
                   chains = 4, adapt = 100, burnin = 50, sample = 100, seed = 4)
  expect_equal(length(fit$mcmc), 4)
  expect_true(all(sapply(fit$mcmc, function(mcmc)dim(mcmc) == c(100, 1))))


  ### checking mcmc_error autofit
  runjags::runjags.options(silent.jags = TRUE, silent.runjags = TRUE)
  fit1 <- JAGS_fit(model_syntax, data, priors_list, parallel = TRUE,
                   autofit = TRUE, autofit_control = list(max_error = 0.05, sample_extend = 100),
                   chains = 4, adapt = 100, burnin = 50, sample = 100, seed = 4)
  expect_equal(length(fit1$mcmc), 4)
  expect_true(all(sapply(fit1$mcmc, function(mcmc)dim(mcmc) == c(200, 1))))

})

test_that("JAGS fit function with JASP works" , {

  set.seed(1)
  data <- list(
    x = rnorm(50, 0, .5),
    N = 50
  )
  priors_list <- list(
    m  = prior("normal", list(0, 1)),
    s  = prior("normal", list(0, 1), list(0, Inf))
  )
  model_syntax <-
    "model
    {
      for(i in 1:N){
        x[i] ~ dnorm(m, pow(s, -2))
      }
    }"

  ### checking the default settings
  set.seed(1)

  ### checking control the main control arguments
  fit <- capture.output(JAGS_fit(model_syntax, data, priors_list, chains = 1, adapt = 100, burnin = 150, sample = 175, thin = 3, seed = 2, is_JASP = TRUE))
  expect_equal(fit, c(
    "Adapting and burnin the model(1)"                                                                                                   ,
    ".Sampling the model(5)"                                                                                                             ,
    "....."                                                                                                                              ,
    "JAGS model with 176 samples (thin = 3; adapt+burnin = 250)"                                                                         ,
    ""                                                                                                                                   ,
    "Full summary statistics have not been pre-calculated - use either the summary method or add.summary to calculate summary statistics",
    ""
    ))

})
