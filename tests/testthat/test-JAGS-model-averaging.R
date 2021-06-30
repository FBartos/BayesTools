context("JAGS model-averaging functions")

test_that("JAGS model-averaging functions work (simple)",{

  set.seed(1)
  data <- list(
    x = rnorm(20, 0, 1),
    N = 20
  )
  priors_list0 <- list(
    m  = prior("spike", list(0)),
    s  = prior("normal", list(0, 1), list(0, Inf))
  )
  priors_list1 <- list(
    m  = prior("normal", list(0, .3)),
    s  = prior("normal", list(0, 1), list(0, Inf))
  )
  model_syntax <-
    "model
    {
      for(i in 1:N){
        x[i] ~ dnorm(m, pow(s, -2))
      }
    }"
  log_posterior <- function(parameters, data){
    sum(stats::dnorm(data$x, parameters[["m"]], parameters[["s"]], log = TRUE))
  }
  # fit the models
  fit0 <- JAGS_fit(model_syntax, data, priors_list0, chains = 1, adapt = 100, burnin = 150, sample = 500, seed = 0)
  fit1 <- JAGS_fit(model_syntax, data, priors_list1, chains = 1, adapt = 100, burnin = 150, sample = 500, seed = 1)
  # get marginal likelihoods
  marglik0 <- JAGS_bridgesampling(fit0, data, priors_list0, log_posterior)
  marglik1 <- JAGS_bridgesampling(fit1, data, priors_list1, log_posterior)

  # make parameter inference
  inference_m             <- compute_inference(c(1, 1), c(marglik0$logml, marglik1$logml), c(T, F))
  inference_m_conditional <- compute_inference(c(1, 1), c(marglik0$logml, marglik1$logml), c(T, F), conditional = T)

  # manually mix posteriors
  mixed_posterior <- BayesTools:::.mix_posteriors.simple(list(fit0, fit1), list(priors_list0[["m"]], priors_list1[["m"]]), "m", inference_m$post_prob)
  mixed_posterior_conditional <- BayesTools:::.mix_posteriors.simple(list(fit0, fit1), list(priors_list0[["m"]], priors_list1[["m"]]), "m", inference_m_conditional$post_prob)

  expect_equal(mean(mixed_posterior == 0), inference_m$post_probs[1], tolerance = 1e-4)
  expect_equal(mean(attr(mixed_posterior, "models_ind") == 1), inference_m$post_probs[1], tolerance = 1e-4)

  expect_doppelganger("JAGS-model-averaging-1", function(){
    par(mfrow = c(1, 2))
    hist(mixed_posterior, main = "model-averaged")
    hist(mixed_posterior_conditional, main = "conditional")
  })

  # automatically mix posteriors
  models <- list(
    list(fit = fit0, marglik = marglik0, priors = priors_list0, prior_weights = 1),
    list(fit = fit1, marglik = marglik1, priors = priors_list1, prior_weights = 1)
  )
  inference <- ensemble_inference(model_list = models, parameters = c("m", "s"), is_null_list = list("m" = 1, "s" = 0), conditional = FALSE)
  inference_conditional <- ensemble_inference(model_list = models, parameters = c("m", "s"), is_null_list = list("m" = 1, "s" = 0), conditional = TRUE)

  mixed_posteriors <- mix_posteriors(model_list = models, parameters = c("m", "s"), is_null_list = list("m" = 1, "s" = 0), seed = 1)
  mixed_posteriors_conditional <- mix_posteriors(model_list = models, parameters = c("m", "s"), is_null_list = list("m" = 1, "s" = 0), conditional = TRUE)

  inference_s             <- compute_inference(c(1, 1), c(marglik0$logml, marglik1$logml), c(F, F))
  inference_s_conditional <- compute_inference(c(1, 1), c(marglik0$logml, marglik1$logml), c(F, F), conditional = T)

  expect_equal(inference$m, inference_m)
  expect_equal(inference_conditional$m, inference_m_conditional)
  expect_equal(inference$s, inference_s)
  expect_equal(inference_conditional$s, inference_s_conditional)
  expect_equal(mean(mixed_posteriors$m == 0), inference_m$post_probs[1], tolerance = 1e-4)
  expect_equal(mean(attr(mixed_posteriors$m, "models_ind") == 1), inference_m$post_probs[1], tolerance = 1e-4)
  expect_equal(mean(attr(mixed_posteriors$s, "models_ind") == 1), inference_m$post_probs[1], tolerance = 1e-4)
  expect_true(all(attr(mixed_posteriors_conditional$m, "models_ind") == 2))
  expect_equal(mean(attr(mixed_posteriors_conditional$s, "models_ind") == 1), inference_m$post_probs[1], tolerance = 1e-4)
  expect_doppelganger("JAGS-model-averaging-2", function(){
    par(mfrow = c(2, 2))
    hist(mixed_posteriors$m, main = "model-averaged (m)")
    hist(mixed_posteriors_conditional$m, main = "conditional (m)")
    hist(mixed_posteriors$s, main = "model-averaged (s)")
    hist(mixed_posteriors_conditional$s, main = "conditional = conditional (s)")
  })

  # dealing with missing unspecified null priors
  models2 <- list(
    list(fit = fit0, marglik = marglik0, priors = list(s  = prior("normal", list(0, 1), list(0, Inf))), prior_weights = 1),
    list(fit = fit1, marglik = marglik1, priors = priors_list1, prior_weights = 1)
  )
  mixed_posteriors2 <- mix_posteriors(model_list = models2, parameters = c("m", "s"), is_null_list = list("m" = 1, "s" = 0), seed = 1)
  expect_equal(mixed_posteriors, mixed_posteriors2)
})

test_that("JAGS model-averaging functions work (weightfunctions)",{

  set.seed(1)
  data <- list(
    x = rnorm(20, 0, 1),
    N = 20
  )
  priors_list0 <- list(
    m     =prior("normal", list(0, 1)),
    omega = prior_none()
  )
  priors_list1 <- list(
    m  = prior("normal", list(0, .5)),
    omega = prior_weightfunction("one.sided", list(c(0.05), c(1, 1)))
  )
  priors_list2 <- list(
    m  = prior("normal", list(0, .3)),
    omega = prior_weightfunction("one.sided", list(c(0.05, 0.50), c(1, 1, 1)))
  )
  model_syntax <-
    "model
    {
      for(i in 1:N){
        x[i] ~ dnorm(m, 1)
      }
    }"
  log_posterior <- function(parameters, data){
    return(0)
  }
  # fit the models
  fit0 <- JAGS_fit(model_syntax, data, priors_list0, chains = 1, adapt = 100, burnin = 150, sample = 500, seed = 0)
  fit1 <- JAGS_fit(model_syntax, data, priors_list1, chains = 1, adapt = 100, burnin = 150, sample = 500, seed = 1)
  fit2 <- JAGS_fit(model_syntax, data, priors_list2, chains = 1, adapt = 100, burnin = 150, sample = 500, seed = 1)
  # get marginal likelihoods
  marglik0 <- JAGS_bridgesampling(fit0, data, priors_list0, log_posterior)
  marglik1 <- JAGS_bridgesampling(fit1, data, priors_list1, log_posterior)
  marglik2 <- JAGS_bridgesampling(fit2, data, priors_list2, log_posterior)

  # check coefficient mapping
  expect_equal(weightfunctions_mapping(list(priors_list0$omega, priors_list1$omega, priors_list2$omega)), list(NULL, c(2, 1, 1), c(3, 2, 1)))
  expect_equal(weightfunctions_mapping(list(
    prior_weightfunction("two.sided", list(c(0.05), c(1, 1))),
    prior_weightfunction("one.sided", list(c(0.05, 0.50), c(1, 1, 1)))
  )), list(
    c(2, 1, 1, 1, 2),
    c(3, 3, 2, 1, 1))
  )
  expect_equal(weightfunctions_mapping(list(
    prior_weightfunction("two.sided", list(c(0.05), c(1, 1))),
    prior_weightfunction("one.sided", list(c(0.05, 0.50, .975), c(1, 1, 1), c(1, 1)))
  )), list(
    c(2, 1, 1, 1, 2),
    c(4, 4, 3, 2, 1))
  )

  models <- list(
    list(fit = fit0, marglik = marglik0, priors = priors_list0, prior_weights = 1),
    list(fit = fit1, marglik = marglik1, priors = priors_list1, prior_weights = 1),
    list(fit = fit2, marglik = marglik2, priors = priors_list2, prior_weights = 1)
  )

  # get models inference &  mix posteriors
  models    <- models_inference(models)
  inference <- ensemble_inference(model_list = models, parameters = c("m", "omega"), is_null_list = list("m" = 0, "omega" = 1), conditional = FALSE)
  inference_conditional <- ensemble_inference(model_list = models, parameters = c("m", "omega"), is_null_list = list("m" = 0, "omega" = 1), conditional = TRUE)
  mixed_posteriors <- mix_posteriors(model_list = models, parameters = c("m", "omega"), is_null_list = list("m" = 0, "omega" = 1), seed = 1)
  mixed_posteriors_conditional <-mix_posteriors(model_list = models, parameters = c("m", "omega"), is_null_list = list("m" = 0, "omega" = 1), seed = 1, conditional = TRUE)

  # checking posteriors and inferences
  expect_equal(names(models[[1]]$inference), c("m_number", "marglik", "prior_prob", "post_prob", "inclusion_BF"))
  expect_equal(unname(unlist(models[[1]]$inference)), c(1.0000000, -1.1023042, 0.3333333, 0.1998118, 0.4994120), tolerance = 1e-4)
  expect_equal(mean(mixed_posteriors$omega[,-1] == 1), inference$omega$post_probs[1], tolerance = 1e-4)
  expect_true(all(mixed_posteriors$omega[1,] == 1))
  expect_true(all(colnames(mixed_posteriors$omega[1,]) == c("omega[0,0.05]", "omega[0.05,0.5]", "omega[0.5,1]")))
  expect_equal(mean(attr(mixed_posteriors$omega, "models_ind") == 2), inference$omega$post_probs[2], tolerance = 1e-4)
  expect_equal(mean(attr(mixed_posteriors$omega, "models_ind") == 3), inference$omega$post_probs[3], tolerance = 1e-4)
  expect_equal(mean(attr(mixed_posteriors_conditional$omega, "models_ind") == 2), inference_conditional$omega$post_probs[2], tolerance = 1e-4)
  expect_equal(mean(attr(mixed_posteriors_conditional$omega, "models_ind") == 3), inference_conditional$omega$post_probs[3], tolerance = 1e-4)
  expect_doppelganger("JAGS-model-averaging-weightfunctions-1", function(){
    par(mfrow = c(2, 3))
    sapply(1:3, function(i)hist(mixed_posteriors$omega[,i], main = "model-averaged (omega)", xlab = colnames(mixed_posteriors$omega)[i]))
    sapply(1:3, function(i)hist(mixed_posteriors_conditional$omega[,i], main = "conditional (omega)", xlab = colnames(mixed_posteriors$omega)[i]))
  })


  ### checking fixed weightfunctions
  priors_list3 <- list(
    m  = prior("normal", list(0, .3)),
    omega = prior_weightfunction("two.sided.fixed", list(0.20, c(.3, 1)))
  )
  fit3     <- JAGS_fit(model_syntax, data, priors_list3, chains = 1, adapt = 100, burnin = 150, sample = 500, seed = 1)
  marglik3 <- JAGS_bridgesampling(fit3, data, priors_list3, log_posterior)
  models <- list(
    list(fit = fit0, marglik = marglik0, priors = priors_list0, prior_weights = 1),
    list(fit = fit1, marglik = marglik1, priors = priors_list1, prior_weights = 1),
    list(fit = fit2, marglik = marglik2, priors = priors_list2, prior_weights = 1),
    list(fit = fit3, marglik = marglik3, priors = priors_list3, prior_weights = 1)
  )

  inference <- ensemble_inference(model_list = models, parameters = c("m", "omega"), is_null_list = list("m" = 0, "omega" = 1), conditional = FALSE)
  mixed_posteriors <- mix_posteriors(model_list = models, parameters = c("m", "omega"), is_null_list = list("m" = 0, "omega" = 1), seed = 1)
  mixed_posteriors_conditional <- mix_posteriors(model_list = models, parameters = c("m", "omega"), is_null_list = list("m" = 0, "omega" = 1), seed = 1, conditional = TRUE)

  expect_equal(mean(mixed_posteriors$omega[,1] == .3), inference$omega$post_probs[4], tolerance = 1e-4)
  expect_equal(mean(mixed_posteriors$omega[,3] == 1), inference$omega$post_probs[4] + inference$omega$post_probs[1], tolerance = 1e-4)
  expect_doppelganger("JAGS-model-averaging-weightfunctions-2", function(){
    par(mfrow = c(2, 5))
    sapply(1:5, function(i)hist(mixed_posteriors$omega[,i], main = "model-averaged (omega)", xlab = colnames(mixed_posteriors$omega)[i]))
    sapply(1:5, function(i)hist(mixed_posteriors_conditional$omega[,i], main = "conditional (omega)", xlab = colnames(mixed_posteriors$omega)[i]))
  })
})
