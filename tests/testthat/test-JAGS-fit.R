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
      lines(priors[[i]], individual = TRUE)
    })
  }
})

test_that("JAGS model functions work (vector)", {

  skip_if_not_installed("rjags")
  model_syntax <- "model{}"
  priors       <- list(
    p1  = prior("mnormal", list(mean = 0, sd = 2, K = 3),),
    p2  = prior("mcauchy", list(location = 0, scale = 1.5, K = 2)),
    p3  = prior("mt",      list(location = 0, scale = 0.5, df = 5, K = 5))
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
      lines(priors[[i]], individual = TRUE)
    })
  }
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


    expect_doppelganger(paste0("JAGS-model-weightfunction-",i), function(){
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
  expect_doppelganger("JAGS-fit-posterior", function(){
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
  summary_4 <- summary(fit4)
  expect_true(summary_4[1,"MCerr"]   > 0.069)
  expect_true(summary_4[1,"MC%ofSD"] > 8)
  expect_true(summary_4[1,"SSeff"]   < 150)
  expect_true(summary_4[1,"psrf"]    > 1.007)

  convergence <- JAGS_check_convergence(fit4, prior_list = priors_list4, max_Rhat = 1.001, min_ESS = 500, max_error = 0.01, max_SD_error = 0.05)
  expect_true(!convergence)
  expect_equal(attr(convergence, "errors")[1], "R-hat 1.007 is larger than the set target (1.001).")
  expect_equal(attr(convergence, "errors")[2], "ESS 149 is lower than the set target (500).")
  expect_equal(attr(convergence, "errors")[3], "MCMC error 0.0698 is larger than the set target (0.01).")
  expect_equal(attr(convergence, "errors")[4], "MCMC SD error 0.082 is larger than the set target (0.05).")

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
})

test_that("JAGS parallel fit function works" , {

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
