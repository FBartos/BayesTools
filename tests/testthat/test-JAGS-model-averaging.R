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
  marglik0 <- JAGS_bridgesampling(fit0, log_posterior = log_posterior, data = data, prior_list = priors_list0)
  marglik1 <- JAGS_bridgesampling(fit1, log_posterior = log_posterior, data = data, prior_list = priors_list1)

  # make parameter inference
  inference_m             <- compute_inference(c(1, 1), c(marglik0$logml, marglik1$logml), c(T, F))
  inference_m_conditional <- compute_inference(c(1, 1), c(marglik0$logml, marglik1$logml), c(T, F), conditional = T)

  # manually mix posteriors
  mixed_posterior <- BayesTools:::.mix_posteriors.simple(list(fit0, fit1), list(priors_list0[["m"]], priors_list1[["m"]]), "m", inference_m$post_prob)
  mixed_posterior_conditional <- BayesTools:::.mix_posteriors.simple(list(fit0, fit1), list(priors_list0[["m"]], priors_list1[["m"]]), "m", inference_m_conditional$post_prob)

  expect_equal(mean(mixed_posterior == 0), inference_m$post_probs[1], tolerance = 1e-4)
  expect_equal(mean(attr(mixed_posterior, "models_ind") == 1), inference_m$post_probs[1], tolerance = 1e-4)

  vdiffr::expect_doppelganger("JAGS-model-averaging-1", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))
    par(mfrow = c(1, 2))
    hist(mixed_posterior, main = "model-averaged")
    hist(mixed_posterior_conditional, main = "conditional")
  })

  # automatically mix posteriors
  models <- list(
    list(fit = fit0, marglik = marglik0, prior_weights = 1),
    list(fit = fit1, marglik = marglik1, prior_weights = 1)
  )
  inference <- ensemble_inference(model_list = models, parameters = c("m", "s"), is_null_list = list("m" = 1, "s" = 0), conditional = FALSE)
  inference_conditional <- ensemble_inference(model_list = models, parameters = c("m", "s"), is_null_list = list("m" = 1, "s" = 0), conditional = TRUE)

  mixed_posteriors <- mix_posteriors(model_list = models, parameters = c("m", "s"), is_null_list = list("m" = 1, "s" = 0), seed = 1)
  mixed_posteriors_conditional <- mix_posteriors(model_list = models, parameters = c("m", "s"), is_null_list = list("m" = 1, "s" = 0), conditional = TRUE)

  inference_s             <- compute_inference(c(1, 1), c(marglik0$logml, marglik1$logml), c(F, F))
  inference_s_conditional <- compute_inference(c(1, 1), c(marglik0$logml, marglik1$logml), c(F, F), conditional = T)

  expect_equal(inference$m[c("prior_probs", "post_probs", "BF")], inference_m[c("prior_probs", "post_probs", "BF")])
  expect_equal(inference_conditional$m[c("prior_probs", "post_probs", "BF")], inference_m_conditional[c("prior_probs", "post_probs", "BF")])
  expect_equal(inference$s[c("prior_probs", "post_probs", "BF")], inference_s[c("prior_probs", "post_probs", "BF")])
  expect_equal(inference_conditional$s[c("prior_probs", "post_probs", "BF")], inference_s_conditional[c("prior_probs", "post_probs", "BF")])
  expect_equal(mean(mixed_posteriors$m == 0), inference_m$post_probs[1], tolerance = 1e-4)
  expect_equal(mean(attr(mixed_posteriors$m, "models_ind") == 1), inference_m$post_probs[1], tolerance = 1e-4)
  expect_equal(mean(attr(mixed_posteriors$s, "models_ind") == 1), inference_m$post_probs[1], tolerance = 1e-4)
  expect_true(all(attr(mixed_posteriors_conditional$m, "models_ind") == 2))
  expect_equal(mean(attr(mixed_posteriors_conditional$s, "models_ind") == 1), inference_m$post_probs[1], tolerance = 1e-4)
  vdiffr::expect_doppelganger("JAGS-model-averaging-2", function(){
    par(mfrow = c(2, 2))
    hist(mixed_posteriors$m, main = "model-averaged (m)")
    hist(mixed_posteriors_conditional$m, main = "conditional (m)")
    hist(mixed_posteriors$s, main = "model-averaged (s)")
    hist(mixed_posteriors_conditional$s, main = "conditional = conditional (s)")
  })

  # dealing with missing unspecified null priors
  models2 <- list(
    list(fit = fit0, marglik = marglik0, prior_weights = 1),
    list(fit = fit1, marglik = marglik1, prior_weights = 1)
  )
  mixed_posteriors2 <- mix_posteriors(model_list = models2, parameters = c("m", "s"), is_null_list = list("m" = 1, "s" = 0), seed = 1)
  expect_equal(mixed_posteriors, mixed_posteriors2)
})

# skip the rest as it takes too long
skip_on_cran()

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
  marglik0 <- JAGS_bridgesampling(fit0, log_posterior = log_posterior, data = data, prior_list = priors_list0)
  marglik1 <- JAGS_bridgesampling(fit1, log_posterior = log_posterior, data = data, prior_list = priors_list1)
  marglik2 <- JAGS_bridgesampling(fit2, log_posterior = log_posterior, data = data, prior_list = priors_list2)

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
    list(fit = fit0, marglik = marglik0, prior_weights = 1),
    list(fit = fit1, marglik = marglik1, prior_weights = 1),
    list(fit = fit2, marglik = marglik2, prior_weights = 1)
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
  vdiffr::expect_doppelganger("JAGS-model-averaging-weightfunctions-1", function(){
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
  marglik3 <- JAGS_bridgesampling(fit3, log_posterior = log_posterior, data = data, prior_list = priors_list3)
  models <- list(
    list(fit = fit0, marglik = marglik0, prior_weights = 1),
    list(fit = fit1, marglik = marglik1, prior_weights = 1),
    list(fit = fit2, marglik = marglik2, prior_weights = 1),
    list(fit = fit3, marglik = marglik3, prior_weights = 1)
  )

  inference <- ensemble_inference(model_list = models, parameters = c("m", "omega"), is_null_list = list("m" = 0, "omega" = 1), conditional = FALSE)
  mixed_posteriors <- mix_posteriors(model_list = models, parameters = c("m", "omega"), is_null_list = list("m" = 0, "omega" = 1), seed = 1)
  mixed_posteriors_conditional <- mix_posteriors(model_list = models, parameters = c("m", "omega"), is_null_list = list("m" = 0, "omega" = 1), seed = 1, conditional = TRUE)

  expect_equal(mean(mixed_posteriors$omega[,1] == .3), inference$omega$post_probs[4], tolerance = 1e-4)
  expect_equal(mean(mixed_posteriors$omega[,3] == 1), inference$omega$post_probs[4] + inference$omega$post_probs[1], tolerance = 1e-4)
  vdiffr::expect_doppelganger("JAGS-model-averaging-weightfunctions-2", function(){
    par(mfrow = c(2, 5))
    sapply(1:5, function(i)hist(mixed_posteriors$omega[,i], main = "model-averaged (omega)", xlab = colnames(mixed_posteriors$omega)[i]))
    sapply(1:5, function(i)hist(mixed_posteriors_conditional$omega[,i], main = "conditional (omega)", xlab = colnames(mixed_posteriors$omega)[i]))
  })
})

test_that("JAGS model-averaging functions work (formula + factors)",{

  skip_on_os(c("mac", "linux", "solaris")) # multivariate sampling does not exactly match across OSes
  skip_on_cran()

  set.seed(1)

  data_formula <- data.frame(
    x_cont1 = rnorm(60),
    x_fac2t = factor(rep(c("A", "B"), 30), levels = c("A", "B")),
    x_fac3o = factor(rep(c("A", "B", "C"), 20), levels = c("A", "B", "C")),
    x_fac3t = factor(rep(c("A", "B", "C"), 20), levels = c("A", "B", "C"))
  )
  data <- list(
    y = rnorm(60, .4 * data_formula$x_cont1 + ifelse(data_formula$x_fac3t == "A", 0.0, ifelse(data_formula$x_fac3t == "B", -0.2, 0.4)), 1),
    N = 60
  )

  # create model with mix of a formula and free parameters ---
  formula_list0 <- list(mu = ~ x_fac2t)
  formula_list1 <- list(mu = ~ x_cont1 + x_fac3t)
  formula_list2 <- list(mu = ~ x_fac3o)
  formula_list3 <- list(mu = ~ x_cont1 * x_fac3o)

  formula_prior_list0 <- list(
    mu    = list(
      "intercept"       = prior("normal", list(0, 5)),
      "x_fac2t"         = prior_factor("normal",  contrast = "treatment", list(0, 1))
    )
  )
  formula_prior_list1 <- list(
    mu    = list(
      "intercept"       = prior("normal", list(0, 5)),
      "x_cont1"         = prior("normal", list(0, 1)),
      "x_fac3t"         = prior_factor("normal", contrast = "treatment", list(0, 1))
    )
  )
  formula_prior_list2 <- list(
    mu    = list(
      "intercept"       = prior("normal", list(0, 5)),
      "x_fac3o"         = prior_factor("mnormal", contrast = "orthonormal", list(0, 1))
    )
  )
  formula_prior_list3 <- list(
    mu    = list(
      "intercept"       = prior("normal", list(0, 5)),
      "x_cont1"         = prior("normal", list(0, 1)),
      "x_fac3o"         = prior_factor("mnormal", contrast = "orthonormal", list(0, 1)),
      "x_cont1:x_fac3o" = prior_factor("mnormal", contrast = "orthonormal", list(0, 1))
    )
  )

  prior_list        <- list(sigma = prior("lognormal", list(0, 1)))
  formula_data_list <- list(mu = data_formula)

  model_syntax <- paste0(
    "model{\n",
    "for(i in 1:N){\n",
    "  y[i] ~ dnorm(mu[i], 1/pow(sigma, 2))\n",
    "}\n",
    "}"
  )

  log_posterior <- function(parameters, data){
    sum(stats::dnorm(data$y, parameters[["mu"]], parameters[["sigma"]], log = TRUE))
  }

  fit0 <- JAGS_fit(
    model_syntax = model_syntax, data = data, prior_list = prior_list,
    formula_list = formula_list0, formula_data_list = formula_data_list, formula_prior_list = formula_prior_list0, seed = 1)
  fit1 <- JAGS_fit(
    model_syntax = model_syntax, data = data, prior_list = prior_list,
    formula_list = formula_list1, formula_data_list = formula_data_list, formula_prior_list = formula_prior_list1, seed = 2)
  fit2 <- JAGS_fit(
    model_syntax = model_syntax, data = data, prior_list = prior_list,
    formula_list = formula_list2, formula_data_list = formula_data_list, formula_prior_list = formula_prior_list2, seed = 3)
  fit3 <- JAGS_fit(
    model_syntax = model_syntax, data = data, prior_list = prior_list,
    formula_list = formula_list3, formula_data_list = formula_data_list, formula_prior_list = formula_prior_list3, seed = 4)

  marglik0 <- JAGS_bridgesampling(
    fit0, log_posterior = log_posterior, data = data, prior_list = prior_list,
    formula_list = formula_list0, formula_data_list = formula_data_list, formula_prior_list = formula_prior_list0)
  marglik1 <- JAGS_bridgesampling(
    fit1, log_posterior = log_posterior, data = data, prior_list = prior_list,
    formula_list = formula_list1, formula_data_list = formula_data_list, formula_prior_list = formula_prior_list1)
  marglik2 <- JAGS_bridgesampling(
    fit2, log_posterior = log_posterior, data = data, prior_list = prior_list,
    formula_list = formula_list2, formula_data_list = formula_data_list, formula_prior_list = formula_prior_list2)
  marglik3 <- JAGS_bridgesampling(
    fit3, log_posterior = log_posterior, data = data, prior_list = prior_list,
    formula_list = formula_list3, formula_data_list = formula_data_list, formula_prior_list = formula_prior_list3)


  # mix posteriors
  models <- list(
    list(fit = fit0, marglik = marglik0, prior_weights = 1),
    list(fit = fit1, marglik = marglik1, prior_weights = 1),
    list(fit = fit2, marglik = marglik2, prior_weights = 1),
    list(fit = fit3, marglik = marglik3, prior_weights = 1)
  )


  inference <- ensemble_inference(
    model_list   = models,
    parameters   = c("mu_x_cont1", "mu_x_fac2t", "mu_x_fac3t", "mu_x_fac3o", "mu_x_cont1__xXx__x_fac3o"),
    is_null_list = list(
      "mu_x_cont1"         = c(TRUE,  FALSE, TRUE,  FALSE),
      "mu_x_fac2t"         = c(FALSE, TRUE,  TRUE,  TRUE),
      "mu_x_fac2t"         = c(TRUE,  FALSE, TRUE,  TRUE),
      "mu_x_fac2t"         = c(TRUE,  TRUE,  FALSE, FALSE),
      "mu_x_cont1:x_fac3o" = c(TRUE,  TRUE,  TRUE,  FALSE)
      ),
    conditional = FALSE)

  mixed_posteriors <- mix_posteriors(
    model_list   = models,
    parameters   = c("mu_x_cont1", "mu_x_fac2t", "mu_x_fac3t", "mu_x_fac3o", "mu_x_cont1__xXx__x_fac3o"),
    is_null_list = list(
      "mu_x_cont1"         = c(TRUE,  FALSE, TRUE,  FALSE),
      "mu_x_fac2t"         = c(FALSE, TRUE,  TRUE,  TRUE),
      "mu_x_fac2t"         = c(TRUE,  FALSE, TRUE,  TRUE),
      "mu_x_fac2t"         = c(TRUE,  TRUE,  FALSE, FALSE),
      "mu_x_cont1:x_fac3o" = c(TRUE,  TRUE,  TRUE,  FALSE)
    ),
    seed = 1, n_samples = 10000)

  mixed_posteriors_c <- mix_posteriors(
    model_list   = models,
    parameters   = c("mu_x_cont1", "mu_x_fac2t", "mu_x_fac3t", "mu_x_fac3o", "mu_x_cont1__xXx__x_fac3o"),
    is_null_list = list(
      "mu_x_cont1"         = c(TRUE,  FALSE, TRUE,  FALSE),
      "mu_x_fac2t"         = c(FALSE, TRUE,  TRUE,  TRUE),
      "mu_x_fac2t"         = c(TRUE,  FALSE, TRUE,  TRUE),
      "mu_x_fac2t"         = c(TRUE,  TRUE,  FALSE, FALSE),
      "mu_x_cont1:x_fac3o" = c(TRUE,  TRUE,  TRUE,  FALSE)
    ),
    seed = 1, n_samples = 10000, conditional = TRUE)


  expect_true(is.numeric(inference$mu_x_cont1$BF))
  expect_true(is.numeric(inference$mu_x_fac2t$BF))
  expect_true(is.numeric(inference$mu_x_fac3t$BF))
  expect_true(is.numeric(inference$mu_x_fac3o$BF))
  expect_true(is.numeric(inference$mu_x_cont1__xXx__x_fac3o$BF))

  expect_equal(length(mixed_posteriors$mu_x_cont1), 10000)
  expect_equal(length(mixed_posteriors$mu_x_fac2t), 10000)
  expect_equal(dim(mixed_posteriors$mu_x_fac3t),               c(10000, 2))
  expect_equal(dim(mixed_posteriors$mu_x_fac3o),               c(10000, 2))
  expect_equal(dim(mixed_posteriors$mu_x_cont1__xXx__x_fac3o), c(10000, 2))

  vdiffr::expect_doppelganger("JAGS-model-averaging-3", function(){
    par(mfrow = c(2, 3))
    hist(mixed_posteriors$mu_x_fac2t,       main = "averaged x_fac2t")
    hist(mixed_posteriors_c$mu_x_fac2t,     main = "conditiona x_fac2t")
    hist(mixed_posteriors_c$mu_x_fac3t[,1], main = "conditional mu_x_fac3t[1]")
    hist(mixed_posteriors_c$mu_x_fac3t[,2], main = "conditional mu_x_fac3t[2]")
    hist(mixed_posteriors_c$mu_x_cont1__xXx__x_fac3o[,1], main = "conditional mu_x_cont1__xXx__x_fac3o[1]")
    hist(mixed_posteriors_c$mu_x_cont1__xXx__x_fac3o[,2], main = "conditional mu_x_cont1__xXx__x_fac3o[2]")
  })

})

test_that("JAGS model-averaging functions work (formula + spike factors)",{

  skip_on_os(c("mac", "linux", "solaris")) # multivariate sampling does not exactly match across OSes
  skip_on_cran()

  set.seed(1)

  data_formula <- data.frame(
    x_fac3md = factor(rep(c("A", "B", "C"), 20), levels = c("A", "B", "C"))
  )
  data <- list(
    y = rnorm(60, ifelse(data_formula$x_fac3md == "A", 0.0, ifelse(data_formula$x_fac3md == "B", -0.2, 0.4)), 1),
    N = 60
  )

  # create model with mix of a formula and free parameters ---
  formula_list0a <- list(mu = ~ 1)
  formula_list0b <- list(mu = ~ x_fac3md)
  formula_list1  <- list(mu = ~ x_fac3md)


  formula_prior_list0a <- list(
    mu    = list(
      "intercept" = prior("normal", list(0, 5))
    )
  )
  formula_prior_list0b <- list(
    mu    = list(
      "intercept" = prior("normal", list(0, 5)),
      "x_fac3md"  = prior_factor("spike",  contrast = "meandif", list(0))
    )
  )
  formula_prior_list1 <- list(
    mu    = list(
      "intercept" = prior("normal", list(0, 5)),
      "x_fac3md"  = prior_factor("mnormal",  contrast = "meandif", list(0, 0.25))
    )
  )

  prior_list        <- list(sigma = prior("lognormal", list(0, 1)))
  formula_data_list <- list(mu = data_formula)

  model_syntax <- paste0(
    "model{\n",
    "for(i in 1:N){\n",
    "  y[i] ~ dnorm(mu[i], 1/pow(sigma, 2))\n",
    "}\n",
    "}"
  )

  log_posterior <- function(parameters, data){
    sum(stats::dnorm(data$y, parameters[["mu"]], parameters[["sigma"]], log = TRUE))
  }

  fit0a <- JAGS_fit(
    model_syntax = model_syntax, data = data, prior_list = prior_list,
    formula_list = formula_list0a, formula_data_list = formula_data_list, formula_prior_list = formula_prior_list0a, seed = 1)
  fit0b <- JAGS_fit(
    model_syntax = model_syntax, data = data, prior_list = prior_list,
    formula_list = formula_list0b, formula_data_list = formula_data_list, formula_prior_list = formula_prior_list0b, seed = 2)
  fit1 <- JAGS_fit(
    model_syntax = model_syntax, data = data, prior_list = prior_list,
    formula_list = formula_list1, formula_data_list = formula_data_list, formula_prior_list = formula_prior_list1, seed = 3)


  marglik0a <- JAGS_bridgesampling(
    fit0a, log_posterior = log_posterior, data = data, prior_list = prior_list,
    formula_list = formula_list0a, formula_data_list = formula_data_list, formula_prior_list = formula_prior_list0a)
  marglik0b <- JAGS_bridgesampling(
    fit0b, log_posterior = log_posterior, data = data, prior_list = prior_list,
    formula_list = formula_list0b, formula_data_list = formula_data_list, formula_prior_list = formula_prior_list0b)
  marglik1 <- JAGS_bridgesampling(
    fit1, log_posterior = log_posterior, data = data, prior_list = prior_list,
    formula_list = formula_list1, formula_data_list = formula_data_list, formula_prior_list = formula_prior_list1)


  # mix posteriors
  modelsA <- list(
    list(fit = fit0a, marglik = marglik0a, prior_weights = 1),
    list(fit = fit1, marglik = marglik1, prior_weights = 1)
  )
  modelsB <- list(
    list(fit = fit0b, marglik = marglik0b, prior_weights = 1),
    list(fit = fit1, marglik = marglik1, prior_weights = 1)
  )


  inferenceA <- ensemble_inference(
    model_list   = modelsA,
    parameters   = c("mu_x_fac3md"),
    is_null_list = list(
      "mu_x_fac3md" = c(TRUE, FALSE)
    ),
    conditional = FALSE)
  inferenceB <- ensemble_inference(
    model_list   = modelsB,
    parameters   = c("mu_x_fac3md"),
    is_null_list = list(
      "mu_x_fac3md" = c(TRUE, FALSE)
    ),
    conditional = FALSE)

  mixed_posteriorsA <- mix_posteriors(
    model_list   = modelsA,
    parameters   = c("mu_x_fac3md"),
    is_null_list = list(
      "mu_x_fac3md" = c(TRUE, FALSE)
    ),
    seed = 1, n_samples = 10000)
  mixed_posteriorsB <- mix_posteriors(
    model_list   = modelsB,
    parameters   = c("mu_x_fac3md"),
    is_null_list = list(
      "mu_x_fac3md" = c(TRUE, FALSE)
    ),
    seed = 1, n_samples = 10000)


  expect_equivalent(inferenceA, inferenceB, tolerance = 1e-2)

  common_attributes <- names(attributes(mixed_posteriorsB$mu_x_fac3md))
  common_attributes <- common_attributes[!common_attributes %in% c("sample_ind", "models_ind", "prior_list")]

  expect_equal(attributes(mixed_posteriorsA$mu_x_fac3md)[common_attributes], attributes(mixed_posteriorsB$mu_x_fac3md)[common_attributes])

})
