context("Model-averaging plot functions")
set.seed(1)

test_that("helper functions work", {

  # join duplicate
  prior_list       <- list(
    p1  = prior("normal", list(0, 1)),
    p2  = prior("lognormal", list(0, .5)),
    p3  = prior("point", list(1)),
    p4  = prior("normal", list(0, 1))
  )

  simplified_list <- .simplify_prior_list(prior_list)

  expect_equal(simplified_list, list(
    p1  = prior("normal", list(0, 1), prior_weights = 2),
    p2  = prior("lognormal", list(0, .5)),
    p3  = prior("point", list(1))
  ))


  # no duplicate
  prior_list       <- list(
    p1  = prior("normal", list(0, 1)),
    p2  = prior("lognormal", list(0, .5)),
    p3  = prior("point", list(1))
  )

  simplified_list <- .simplify_prior_list(prior_list)

  expect_equal(simplified_list, list(
    p1  = prior("normal", list(0, 1)),
    p2  = prior("lognormal", list(0, .5)),
    p3  = prior("point", list(1))
  ))


  # multiple duplicates
  prior_list       <- list(
    p1  = prior("normal", list(0, 1)),
    p2  = prior("lognormal", list(0, .5)),
    p3  = prior("point", list(1)),
    p1  = prior("normal", list(0, 1)),
    p4  = prior("normal", list(0, 1)),
    p2  = prior("lognormal", list(0, .5))
  )

  simplified_list <- .simplify_prior_list(prior_list)

  expect_equal(simplified_list, list(
    p1  = prior("normal", list(0, 1), prior_weights = 3),
    p2  = prior("lognormal", list(0, .5), prior_weights = 2),
    p3  = prior("point", list(1))
  ))
})


test_that("prior plot functions (simple) work", {

  ### simple cases
  # continuous
  prior_list       <- list(
    p1 = prior("normal", list(0, 1))
  )

  expect_doppelganger("model-averaging-plot-prior-simple-1", function(){
    plot_prior_list(prior_list, col = "red", lwd = 4)
    lines(prior_list$p1, col = "blue", lwd = 3, lty = 2)
  })
  expect_doppelganger("model-averaging-plot-prior-simple-2", {
    plot_prior_list(prior_list, plot_type = "ggplot", col = "red", lwd = 4) + geom_prior(prior_list$p1, col = "blue", lwd = 3, lty = 2)
  })

  # spike
  prior_list       <- list(
    p1 = prior("spike", list(.5))
  )

  expect_doppelganger("model-averaging-plot-prior-simple-3", function(){
    plot_prior_list(prior_list, col = "red", lwd = 4)
    lines(prior_list$p1, col = "blue", lwd = 3, lty = 2)
  })
  expect_doppelganger("model-averaging-plot-prior-simple-4", {
    plot_prior_list(prior_list, plot_type = "ggplot", col = "red", lwd = 4) + geom_prior(prior_list$p1, col = "blue", lwd = 3, lty = 2)
  })


  ### the prior joining should give the same prior (+ check truncation)
  prior_list       <- list(
    p1 = prior("normal", list(0, 1), truncation = list(0, Inf)),
    p2 = prior("normal", list(0, 1.001), truncation = list(0, Inf))
  )
  expect_doppelganger("model-averaging-plot-prior-simple-5", function(){
    plot_prior_list(prior_list, col = "red", lwd = 4)
    lines(prior_list$p1, col = "blue", lwd = 3, lty = 2)
  })
  expect_doppelganger("model-averaging-plot-prior-simple-6", function(){
    plot_prior_list(prior_list, plot_type = "ggplot", col = "red", lwd = 4) + geom_prior(prior_list$p1, col = "blue", lwd = 3, lty = 2)
  })


  ### mixtures
  prior_list       <- list(
    p1 = prior("normal", list(0, 1)),
    p2 = prior("normal", list(0, 1), list(1, Inf)),
    p3 = prior("spike", list(.5))
  )
  expect_doppelganger("model-averaging-plot-prior-simple-7", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))
    par(mar = c(4, 4, 1, 4))
    plot_prior_list(prior_list)
  })
  expect_doppelganger("model-averaging-plot-prior-simple-8", function(){
    plot_prior_list(prior_list, plot_type = "ggplot")
  })

  # with additional settings
  expect_doppelganger("model-averaging-plot-prior-simple-9", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))
    par(mar = c(4, 4, 1, 4))
    plot_prior_list(prior_list, xlab = "xlab", ylab = "ylab", ylab2 = "ylab2", main = "main")
  })
  expect_doppelganger("model-averaging-plot-prior-simple-10", function(){
    plot_prior_list(prior_list, plot_type = "ggplot", xlab = "xlab", ylab = "ylab", ylab2 = "ylab2", main = "main")
  })

  # and more spikes
  prior_list       <- list(
    p1 = prior("normal", list(0, 1)),
    p2 = prior("normal", list(0, 1), list(1, Inf)),
    p3 = prior("spike", list(.5)),
    p4 = prior("spike", list(-5))
  )
  expect_doppelganger("model-averaging-plot-prior-simple-11", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))
    par(mar = c(4, 4, 1, 4))
    plot_prior_list(prior_list)
  })
  expect_doppelganger("model-averaging-plot-prior-simple-12", function(){
    plot_prior_list(prior_list, plot_type = "ggplot")
  })

  # verify aggregation
  prior_list       <- list(
    p1 = prior("normal", list(0, 1)),
    p2 = prior("normal", list(0, 1)),
    p3 = prior("spike", list(.5)),
    p4 = prior("spike", list(.5))
  )
  expect_doppelganger("model-averaging-plot-prior-simple-13", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))
    par(mar = c(4, 4, 1, 4))
    plot_prior_list(prior_list, lwd = 2)
    lines_prior_list(prior_list, lty = 2, col = "red", lwd = 2)
  })
})


test_that("prior plot functions (PET-PEESE) work", {

  ### simple cases
  # continuous
  prior_list       <- list(
    p1 = prior_PET("cauchy",   list(0, 1))
  )
  prior_list_mu   <- list(
    m1 = prior("spike", list(0))
  )
  expect_doppelganger("model-averaging-plot-prior-PETPEESE-1", function(){
    plot_prior_list(prior_list, col = "red", lwd = 4, col.fill = ggplot2::alpha("red", .20), n_samples = 1000, n_points = 50, prior_list_mu = prior_list_mu)
    lines(prior_list$p1, col = "blue", lwd = 3, lty = 2, col.fill = ggplot2::alpha("blue", .20))
  })
  prior_list       <- list(
    p1 = prior_PEESE("cauchy",   list(0, 2))
  )
  expect_doppelganger("model-averaging-plot-prior-PETPEESE-2", {
    plot_prior_list(prior_list, plot_type = "ggplot", col = "red", lwd = 4, col.fill = scales::alpha("red", .20), n_samples = 1000, n_points = 50, prior_list_mu = prior_list_mu) + geom_prior(prior_list$p1, col = "blue", lwd = 3, lty = 2, col.fill = scales::alpha("blue", .20), n_samples = 1000, n_points = 50, prior_list_mu = prior_list_mu)
  })

  # spike
  prior_list       <- list(
    p1 = prior_PET("point", list(.1))
  )
  expect_doppelganger("model-averaging-plot-prior-PETPEESE-3", function(){
    plot_prior_list(prior_list, col = "red", lwd = 4, n_samples = 1000, n_points = 50, prior_list_mu = prior_list_mu)
    lines(prior_list$p1, col = "blue", lwd = 3, lty = 2)
  })
  prior_list       <- list(
    p1 = prior_PEESE("point", list(.05))
  )
  expect_doppelganger("model-averaging-plot-prior-PETPEESE-4", {
    plot_prior_list(prior_list, plot_type = "ggplot", col = "red", lwd = 4, n_samples = 1000, n_points = 50, prior_list_mu = prior_list_mu) + geom_prior(prior_list$p1, col = "blue", lwd = 3, lty = 2, n_samples = 1000, n_points = 50, prior_list_mu = prior_list_mu)
  })


  ### the prior joining should give the same prior
  prior_list       <- list(
    PET1  = prior_PET("cauchy",   list(0, 1)),
    PET2  = prior_PET("cauchy",   list(0, 1.001))
  )
  prior_list_mu   <- list(
    m1 = prior("spike", list(0)),
    m2 = prior("spike", list(0))
  )
  expect_doppelganger("model-averaging-plot-prior-PETPEESE-5", function(){
    plot_prior_list(prior_list, col = "red", lwd = 4, col.fill = scales::alpha("red", .20), n_samples = 1000, n_points = 50, prior_list_mu = prior_list_mu)
    lines(prior_list$PET1, col = "blue", lwd = 3, lty = 2, col.fill = scales::alpha("blue", .20))
  })
  ### the prior joining should give the same prior
  prior_list       <- list(
    PEESE1 = prior_PEESE("cauchy",   list(0, 1)),
    PEESE2 = prior_PEESE("cauchy",   list(0, 1.001))
  )
  expect_doppelganger("model-averaging-plot-prior-PETPEESE-6", function(){
    plot_prior_list(prior_list, plot_type = "ggplot", col = "red", col.fill = scales::alpha("red", .20), lwd = 4, n_samples = 1000, n_points = 50, prior_list_mu = prior_list_mu) + geom_prior_list(prior_list, col = "blue", col.fill = scales::alpha("blue", .20), lwd = 3, lty = 2, n_samples = 1000, n_points = 50, prior_list_mu = prior_list_mu)
  })


  ### mixtures
  prior_list       <- list(
    p1 = prior_PET("cauchy",     list(0, 1)),
    p2 = prior_PEESE("cauchy",   list(0, 5))
  )
  expect_doppelganger("model-averaging-plot-prior-PETPEESE-7", function(){
    plot_prior_list(prior_list, col = "red", lwd = 4, col.fill = scales::alpha("red", .20), n_samples = 1000, n_points = 50, prior_list_mu = prior_list_mu)
    lines_prior_list(prior_list, col = "blue", lwd = 3, lty = 2, col.fill = scales::alpha("blue", .20), n_samples = 1000, n_points = 50, prior_list_mu = prior_list_mu)
  })
  expect_doppelganger("model-averaging-plot-prior-PETPEESE-8", function(){
    plot_prior_list(prior_list, plot_type = "ggplot", col = "red", lwd = 4, col.fill = scales::alpha("red", .20), n_samples = 1000, n_points = 50, prior_list_mu = prior_list_mu) + geom_prior_list(prior_list, col = "blue", lwd = 3, lty = 2, col.fill = scales::alpha("blue", .20), n_samples = 1000, n_points = 50, prior_list_mu = prior_list_mu)
  })

  # with additional settings
  expect_doppelganger("model-averaging-plot-prior-PETPEESE-9", function(){
    plot_prior_list(prior_list, n_samples = 1000, n_points = 50, xlab = "xlab", ylab = "ylab", main = "main", prior_list_mu = prior_list_mu)
  })
  expect_doppelganger("model-averaging-plot-prior-PETPEESE-10", function(){
    plot_prior_list(prior_list, n_samples = 1000, n_points = 50, plot_type = "ggplot", xlab = "xlab", ylab = "ylab", main = "main", prior_list_mu = prior_list_mu)
  })


  ### dealing with other type of priors
  prior_list       <- list(
    p1 = prior_PET("cauchy",     list(0, 1)),
    p2 = prior_PEESE("cauchy",   list(0, 5)),
    p3 = prior_none(prior_weights = 4)
  )
  prior_list_mu    <- list(
    m1 = prior("spike",  list(0)),
    m2 = prior("spike",  list(0)),
    m3 = prior("normal", list(.3, .15))
  )
  expect_doppelganger("model-averaging-plot-prior-PETPEESE-11", function(){
    plot_prior_list(prior_list, col = "red", lwd = 4, col.fill = scales::alpha("red", .20), n_samples = 1000, n_points = 50, ylim = c(0, .5), prior_list_mu = prior_list_mu)
    lines_prior_list(prior_list, col = "blue", lwd = 3, lty = 2, col.fill = scales::alpha("blue", .20), n_samples = 1000, n_points = 50, prior_list_mu = prior_list_mu)
  })
  expect_doppelganger("model-averaging-plot-prior-PETPEESE-12", function(){
    plot_prior_list(prior_list, plot_type = "ggplot", col = "red", lwd = 4, col.fill = scales::alpha("red", .20), n_samples = 1000, n_points = 50, ylim = c(0, .5), prior_list_mu = prior_list_mu) + geom_prior_list(prior_list, col = "blue", lwd = 3, lty = 2, col.fill = scales::alpha("blue", .20), n_samples = 1000, n_points = 50, prior_list_mu = prior_list_mu)
  })
})


test_that("prior plot functions (weightfunctions) work", {

  ### simple cases
  # continuous
  prior_list       <- list(
    p1 = prior_weightfunction("one.sided", list(c(.05, 0.10), c(1, 1, 1)))
  )

  expect_doppelganger("model-averaging-plot-prior-wf-1", function(){
    plot_prior_list(prior_list, col = "red", lwd = 4, col.fill = scales::alpha("red", .20))
    lines(prior_list$p1, col = "blue", lwd = 3, lty = 2, col.fill = scales::alpha("blue", .20))
  })
  expect_doppelganger("model-averaging-plot-prior-wf-2", {
    plot_prior_list(prior_list, plot_type = "ggplot", col = "red", lwd = 4, col.fill = scales::alpha("red", .20)) + geom_prior(prior_list$p1, col = "blue", lwd = 3, lty = 2, col.fill = scales::alpha("blue", .20))
  })

  # spike
  prior_list       <- list(
    p1 = prior_weightfunction("one.sided.fixed", list(c(.05), c(1, .5)))
  )

  expect_doppelganger("model-averaging-plot-prior-wf-3", function(){
    plot_prior_list(prior_list, col = "red", lwd = 4)
    lines(prior_list$p1, col = "blue", lwd = 3, lty = 2)
  })
  expect_doppelganger("model-averaging-plot-prior-wf-4", {
    plot_prior_list(prior_list, plot_type = "ggplot", col = "red", lwd = 4) + geom_prior(prior_list$p1, col = "blue", lwd = 3, lty = 2)
  })


  ### the prior joining should give the same prior
  prior_list       <- list(
    p1 = prior_weightfunction("one.sided", list(c(.05, 0.10), c(1, 1, 1))),
    p2 = prior_weightfunction("one.sided", list(c(.05, 0.10), c(1, 1, 1.0001)))
  )
  expect_doppelganger("model-averaging-plot-prior-wf-5", function(){
    plot_prior_list(prior_list, col = "red", lwd = 4, col.fill = scales::alpha("red", .20))
    lines(prior_list$p1, col = "blue", lwd = 3, lty = 2, col.fill = scales::alpha("blue", .20))
  })
  expect_doppelganger("model-averaging-plot-prior-wf-6", function(){
    plot_prior_list(prior_list, plot_type = "ggplot", col = "red", lwd = 4, col.fill = scales::alpha("red", .20)) + geom_prior_list(prior_list, col = "blue", lwd = 3, lty = 2, col.fill = scales::alpha("blue", .20))
  })


  ### mixtures
  prior_list       <- list(
    p1 = prior_weightfunction("one.sided", list(c(.025), c(1, 1))),
    p2 = prior_weightfunction("two.sided", list(c(.05),  c(1, 1))),
    p3 = prior_weightfunction("one.sided.fixed", list(c(.05), c(1, .5)), prior_weights = 10)
  )
  expect_doppelganger("model-averaging-plot-prior-wf-7", function(){
    plot_prior_list(prior_list, col = "red", lwd = 4, col.fill = scales::alpha("red", .20), rescale_x = TRUE)
    lines_prior_list(prior_list, col = "blue", lwd = 3, lty = 2, col.fill = scales::alpha("blue", .20), rescale_x = TRUE)
  })
  expect_doppelganger("model-averaging-plot-prior-wf-8", function(){
    plot_prior_list(prior_list, plot_type = "ggplot", col = "red", lwd = 4, col.fill = scales::alpha("red", .20), rescale_x = TRUE) + geom_prior_list(prior_list, col = "blue", lwd = 3, lty = 2, col.fill = scales::alpha("blue", .20), rescale_x = TRUE)
  })

  # with additional settings
  expect_doppelganger("model-averaging-plot-prior-wf-9", function(){
    plot_prior_list(prior_list, xlab = "xlab", ylab = "ylab", main = "main")
  })
  expect_doppelganger("model-averaging-plot-prior-wf-10", function(){
    plot_prior_list(prior_list, plot_type = "ggplot", xlab = "xlab", ylab = "ylab", main = "main")
  })


  ### dealing with other type of priors
  prior_list       <- list(
    p1 = prior_weightfunction("one.sided", list(c(.5), c(1, 1))),
    p2 = prior_none(),
    p3 = prior_none()
  )
  expect_doppelganger("model-averaging-plot-prior-wf-11", function(){
    plot_prior_list(prior_list, col = "red", lwd = 4, col.fill = scales::alpha("red", .20), rescale_x = TRUE)
    lines_prior_list(prior_list, col = "blue", lwd = 3, lty = 2, col.fill = scales::alpha("blue", .20), rescale_x = TRUE)
  })
  expect_doppelganger("model-averaging-plot-prior-wf-12", function(){
    plot_prior_list(prior_list, plot_type = "ggplot", col = "red", lwd = 4, col.fill = scales::alpha("red", .20), rescale_x = TRUE) + geom_prior_list(prior_list, col = "blue", lwd = 3, lty = 2, col.fill = scales::alpha("blue", .20), rescale_x = TRUE)
  })

})


test_that("posterior plot functions (simple) work", {

  set.seed(1)
  data <- NULL
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
    }"
  log_posterior <- function(parameters, data){
    return(0)
  }
  # fit the models
  fit0 <- suppressWarnings(JAGS_fit(model_syntax, data, priors_list0, chains = 1, adapt = 100, burnin = 150, sample = 500, seed = 0))
  fit1 <- suppressWarnings(JAGS_fit(model_syntax, data, priors_list1, chains = 1, adapt = 100, burnin = 150, sample = 500, seed = 1))
  # get marginal likelihoods
  marglik0 <- JAGS_bridgesampling(fit0, log_posterior = log_posterior, data = data, prior_list = priors_list0)
  marglik1 <- JAGS_bridgesampling(fit1, log_posterior = log_posterior, data = data, prior_list = priors_list1)
  # automatically mix posteriors
  models <- list(
    list(fit = fit0, marglik = marglik0, priors = priors_list0, prior_weights = 1),
    list(fit = fit1, marglik = marglik1, priors = priors_list1, prior_weights = 1)
  )
  mixed_posteriors <- mix_posteriors(model_list = models, parameters = c("m", "s"), is_null_list = list("m" = 1, "s" = 0), seed = 1)


  expect_doppelganger("model-averaging-plot-posterior-simple-1", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))
    par(mar = c(4, 4, 1, 4))
    plot_posterior(mixed_posteriors, "m", lwd = 2, col = "red", par_name = expression(mu))
    lines_prior_list(attr(mixed_posteriors$m, "prior_list"), col = "blue")
  })
  expect_doppelganger("model-averaging-plot-posterior-simple-2", {
    plot_posterior(mixed_posteriors, "m", plot_type = "ggplot", lwd = 2, col = "red") + geom_prior_list(attr(mixed_posteriors$m, "prior_list"), col = "blue")
  })

  # checks truncation
  expect_doppelganger("model-averaging-plot-posterior-simple-3", function(){
    plot_posterior(mixed_posteriors, "s")
    lines_prior_list(attr(mixed_posteriors$s, "prior_list"), col = "blue")
  })
  expect_doppelganger("model-averaging-plot-posterior-simple-4", {
    plot_posterior(mixed_posteriors, "s", plot_type = "ggplot")
  })

  # check transformation
  expect_doppelganger("model-averaging-plot-posterior-simple-5", function(){
    plot_posterior(mixed_posteriors, "s", transformation = "exp")
    lines_prior_list(attr(mixed_posteriors$s, "prior_list"), col = "blue", transformation = "exp")
  })

  # prior and posterior
  expect_doppelganger("model-averaging-plot-posterior-simple-6", function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mar = oldpar[["mar"]]))
    par(mar = c(4, 4, 1, 4))
    plot_posterior(mixed_posteriors, "m", lwd = 2, col = "red", prior = TRUE, dots_prior = list(col = "blue", lty = 2))
  })

  expect_doppelganger("model-averaging-plot-posterior-simple-7", function(){
    plot_posterior(mixed_posteriors, "m", plot_type = "ggplot", lwd = 2, col = "red", prior = TRUE, dots_prior = list(col = "blue", lty = 2))
  })

  # check transformation
  expect_doppelganger("model-averaging-plot-posterior-simple-8", function(){
    plot_posterior(mixed_posteriors, "s", transformation = "exp", lwd = 2, col = "red", prior = TRUE, dots_prior = list(col = "blue", lty = 2))
  })
})


test_that("posterior plot functions (PET-PEESE) work", {

  set.seed(1)
  data <- NULL
  priors_list0 <- list(
    mu    = prior("spike", list(0)),
    PET   = prior_PET("normal", list(0, .2))
  )
  priors_list1 <- list(
    mu    = prior("spike", list(0)),
    PEESE = prior_PEESE("normal", list(0, .8))
  )
  model_syntax <-
    "model
    {
    }"
  log_posterior <- function(parameters, data){
    return(0)
  }
  # fit the models
  fit0 <- suppressWarnings(JAGS_fit(model_syntax, data, priors_list0, chains = 1, adapt = 100, burnin = 150, sample = 2000, seed = 0))
  fit1 <- suppressWarnings(JAGS_fit(model_syntax, data, priors_list1, chains = 1, adapt = 100, burnin = 150, sample = 2000, seed = 1))
  # get marginal likelihoods
  marglik0 <- JAGS_bridgesampling(fit0, log_posterior = log_posterior, data = data, prior_list = priors_list0)
  marglik1 <- JAGS_bridgesampling(fit1, log_posterior = log_posterior, data = data, prior_list = priors_list1)
  # automatically mix posteriors
  models <- list(
    list(fit = fit0, marglik = marglik0, priors = priors_list0, prior_weights = 1),
    list(fit = fit1, marglik = marglik1, priors = priors_list1, prior_weights = 1)
  )
  mixed_posteriors <- mix_posteriors(model_list = models, parameters = c("mu", "PET", "PEESE"), is_null_list = list("mu" = c(T, T), "PET" = c(F,T),  "PEESE" = c(T,F)), seed = 1)



  expect_doppelganger("model-averaging-plot-posterior-PETPEESE-1", function(){
    plot_posterior(mixed_posteriors, "PETPEESE", lwd = 2, col = "red", col.fill = scales::alpha("red", .20), par_name = "PET-PEESE", n_points = 50, ylim = c(0, 1))
    lines_prior_list(list(priors_list0$PET, priors_list1$PEESE), n_points = 50, n_samples = 1000, col = "blue", col.fill = scales::alpha("blue", .20), prior_list_mu = list(priors_list0$mu, priors_list1$mu))
  })
  expect_doppelganger("model-averaging-plot-posterior-PETPEESE-2", {
    plot_posterior(mixed_posteriors, "PETPEESE", plot_type = "ggplot", lwd = 2, col = "red", col.fill = scales::alpha("red", .20), ylim = c(0, .5)) + geom_prior_list(list(priors_list0$PET, priors_list1$PEESE), n_points = 50, n_samples = 1000, col = "blue", col.fill = scales::alpha("blue", .20), prior_list_mu = list(priors_list0$mu, priors_list1$mu))
  })

  # check transformation
  expect_doppelganger("model-averaging-plot-posterior-PETPEESE-5", function(){
    plot_posterior(mixed_posteriors, "PETPEESE", transformation = "lin", transformation_arguments = list(a = 0, b = 0.5), main = "PET-PEESE (1/2x)")
    lines_prior_list(list(priors_list0$PET, priors_list1$PEESE), n_points = 50, n_samples = 1000, col = "blue", col.fill = scales::alpha("blue", .20), transformation = "lin", transformation_arguments = list(a = 0, b = 0.5), prior_list_mu = list(priors_list0$mu, priors_list1$mu))
  })

  # prior and posterior
  expect_doppelganger("model-averaging-plot-posterior-PETPEESE-6", function(){
    plot_posterior(mixed_posteriors, "PETPEESE", lwd = 2, col = "red", col.fill = scales::alpha("red", .20), prior = TRUE, n_points = 50, n_samples = 1000, dots_prior = list(col = "blue", col.fill = scales::alpha("blue", .20), lty = 2))
  })

  expect_doppelganger("model-averaging-plot-posterior-PETPEESE-7", function(){
    plot_posterior(mixed_posteriors, "PETPEESE", plot_type = "ggplot", lwd = 2, col = "red", col.fill = scales::alpha("red", .20), n_points = 50, n_samples = 1000, prior = TRUE, dots_prior = list(col = "blue", col.fill = scales::alpha("blue", .20), lty = 2))
  })

  # check transformation
  expect_doppelganger("model-averaging-plot-posterior-PETPEESE-8", function(){
    plot_posterior(mixed_posteriors, "PETPEESE", transformation = "lin", transformation_arguments = list(a = 0, b = 0.5), lwd = 2, col = "red", n_points = 50, n_samples = 1000, col.fill = scales::alpha("red", .20), prior = TRUE, dots_prior = list(col = "blue", col.fill = scales::alpha("blue", .20), lty = 2))
  })

  # add an overhelming missing model
  priors_list2 <- list(
    mu = prior("normal", list(.2, .2), prior_weights = 4)
  )
  # fit the models
  fit2 <- suppressWarnings(JAGS_fit(model_syntax, data, priors_list2, chains = 1, adapt = 100, burnin = 150, sample = 2000, seed = 1))
  # get marginal likelihoods
  marglik2 <- JAGS_bridgesampling(fit2, log_posterior = log_posterior, data = data, prior_list = priors_list2)
  # automatically mix posteriors
  models <- list(
    list(fit = fit0, marglik = marglik0, priors = priors_list0, prior_weights = 1),
    list(fit = fit1, marglik = marglik1, priors = priors_list1, prior_weights = 1),
    list(fit = fit2, marglik = marglik2, priors = priors_list2, prior_weights = 4)
  )
  mixed_posteriors <- mix_posteriors(model_list = models, parameters = c("mu" ,"PET", "PEESE"), is_null_list = list("mu" = c(T, T, F),"PET" = c(F,T,F),  "PEESE" = c(T,F,F)), seed = 1)
  expect_doppelganger("model-averaging-plot-posterior-PETPEESE-9", function(){
    plot_posterior(mixed_posteriors, "PETPEESE", ylim = c(0, 3), lwd = 2, col = "red", col.fill = scales::alpha("red", .20), n_points = 50, n_samples = 1000, prior = TRUE, dots_prior = list(col = "blue", col.fill = scales::alpha("blue", .20), lty = 2))
  })

})


test_that("posterior plot functions (weightfunctions) work", {

  set.seed(1)
  data <- NULL
  priors_list0 <- list(
    omega = prior_weightfunction("one.sided", list(c(.025), c(1, 1)))
  )
  priors_list1 <- list(
    omega = prior_weightfunction("two.sided", list(c(.05),  c(1, 1)))
  )
  model_syntax <-
    "model
    {
    }"
  log_posterior <- function(parameters, data){
    return(0)
  }
  # fit the models
  fit0 <- suppressWarnings(JAGS_fit(model_syntax, data, priors_list0, chains = 1, adapt = 100, burnin = 150, sample = 2000, seed = 0))
  fit1 <- suppressWarnings(JAGS_fit(model_syntax, data, priors_list1, chains = 1, adapt = 100, burnin = 150, sample = 2000, seed = 1))
  # get marginal likelihoods
  marglik0 <- JAGS_bridgesampling(fit0, log_posterior = log_posterior, data = data, prior_list = priors_list0)
  marglik1 <- JAGS_bridgesampling(fit1, log_posterior = log_posterior, data = data, prior_list = priors_list1)
  # automatically mix posteriors
  models <- list(
    list(fit = fit0, marglik = marglik0, priors = priors_list0, prior_weights = 1),
    list(fit = fit1, marglik = marglik1, priors = priors_list1, prior_weights = 1)
  )
  mixed_posteriors <- mix_posteriors(model_list = models, parameters = "omega", is_null_list = list("omega" = c(F,F)), seed = 1)



  expect_doppelganger("model-averaging-plot-posterior-wf-1", function(){
    plot_posterior(mixed_posteriors, "omega", lwd = 2, col = "red", col.fill = scales::alpha("red", .20), par_name = "Selection Models", n_points = 50, ylim = c(0, 1))
    lines_prior_list(list(priors_list0$omega, priors_list1$omega), n_points = 50, n_samples = 1000, col = "blue", col.fill = scales::alpha("blue", .20))
  })
  expect_doppelganger("model-averaging-plot-posterior-wf-2", {
    plot_posterior(mixed_posteriors, "omega", plot_type = "ggplot", lwd = 2, col = "red", col.fill = scales::alpha("red", .20)) + geom_prior_list(list(priors_list0$omega, priors_list1$omega), n_points = 50, n_samples = 1000, col = "blue", col.fill = scales::alpha("blue", .20))
  })

  # rescale-x
  expect_doppelganger("model-averaging-plot-posterior-wf-3", function(){
    plot_posterior(mixed_posteriors, "omega", lwd = 2, rescale_x = TRUE, col = "red", col.fill = scales::alpha("red", .20), par_name = "Selection Models", n_points = 50, ylim = c(0, 1))
    lines_prior_list(list(priors_list0$omega, priors_list1$omega), rescale_x = TRUE, n_points = 50, n_samples = 1000, col = "blue", col.fill = scales::alpha("blue", .20))
  })
  expect_doppelganger("model-averaging-plot-posterior-wf-4", {
    plot_posterior(mixed_posteriors, "omega", rescale_x = TRUE, plot_type = "ggplot", lwd = 2, col = "red", col.fill = scales::alpha("red", .20)) + geom_prior_list(list(priors_list0$omega, priors_list1$omega), rescale_x = TRUE, n_points = 50, n_samples = 1000, col = "blue", col.fill = scales::alpha("blue", .20))
  })

  # prior and posterior
  expect_doppelganger("model-averaging-plot-posterior-wf-6", function(){
    plot_posterior(mixed_posteriors, "omega", lwd = 2, col = "red", col.fill = scales::alpha("red", .20), prior = TRUE, n_points = 50, n_samples = 1000, dots_prior = list(col = "blue", col.fill = scales::alpha("blue", .20), lty = 2))
  })

  expect_doppelganger("model-averaging-plot-posterior-wf-7", function(){
    plot_posterior(mixed_posteriors, "omega", plot_type = "ggplot", lwd = 2, col = "red", col.fill = scales::alpha("red", .20), n_points = 50, n_samples = 1000, prior = TRUE, dots_prior = list(col = "blue", col.fill = scales::alpha("blue", .20), lty = 2))
  })

  # rescale-x
  expect_doppelganger("model-averaging-plot-posterior-wf-8", function(){
    plot_posterior(mixed_posteriors, "omega", rescale_x = TRUE, lwd = 2, col = "red", col.fill = scales::alpha("red", .20), prior = TRUE, n_points = 50, n_samples = 1000, dots_prior = list(col = "blue", col.fill = scales::alpha("blue", .20), lty = 2))
  })

  expect_doppelganger("model-averaging-plot-posterior-wf-9", function(){
    plot_posterior(mixed_posteriors, "omega", rescale_x = TRUE, plot_type = "ggplot", lwd = 2, col = "red", col.fill = scales::alpha("red", .20), n_points = 50, n_samples = 1000, prior = TRUE, dots_prior = list(col = "blue", col.fill = scales::alpha("blue", .20), lty = 2))
  })

  # add an overhelming missing model
  priors_list2 <- list(
    mu = prior("normal", list(0, .8), prior_weights = 4)
  )
  # fit the models
  fit2 <- suppressWarnings(JAGS_fit(model_syntax, data, priors_list2, chains = 1, adapt = 100, burnin = 150, sample = 2000, seed = 1))
  # get marginal likelihoods
  marglik2 <- JAGS_bridgesampling(fit2, log_posterior = log_posterior, data = data, prior_list = priors_list2)
  # automatically mix posteriors
  models <- list(
    list(fit = fit0, marglik = marglik0, priors = priors_list0, prior_weights = 1),
    list(fit = fit1, marglik = marglik1, priors = priors_list1, prior_weights = 1),
    list(fit = fit2, marglik = marglik2, priors = priors_list2, prior_weights = 5)
  )
  mixed_posteriors <- mix_posteriors(model_list = models, parameters = "omega", is_null_list = list("omega" = c(F,F,F)), seed = 1)
  expect_doppelganger("model-averaging-plot-posterior-wf-10", function(){
    plot_posterior(mixed_posteriors, "omega", lwd = 2, col = "red", col.fill = scales::alpha("red", .20), n_points = 50, n_samples = 1000, prior = TRUE, dots_prior = list(col = "blue", col.fill = scales::alpha("blue", .20), lty = 2))
  })

})


test_that("models plot functions work", {

  ### prior distribution related functions
  p0 <- prior(distribution = "point",  parameters = list(location = 0))
  p1 <- prior(distribution = "normal", parameters = list(mean = 0, sd = 1))
  p2 <- prior(distribution = "normal", parameters = list(mean = 0, sd = 1), truncation = list(0, Inf))

  data <- list(
    x = rnorm(10),
    N = 10
  )

  ## create and fit models
  priors_list0 <- list(mu  = p0)
  priors_list1 <- list(mu  = p1)
  priors_list2 <- list(tau = p2)

  # define likelihood for the data
  model_syntax <-
    "model{
    for(i in 1:N){
      x[i] ~ dnorm(mu, 1)
    }
  }"
  model_syntax2 <-
    "model{
    for(i in 1:N){
      x[i] ~ dnorm(0, pow(tau, -2))
    }
  }"

  # define log posterior for bridge sampling
  log_posterior <- function(parameters, data){
    sum(dnorm(data$x, parameters$mu, 1, log = TRUE))
  }
  log_posterior2 <- function(parameters, data){
    sum(dnorm(data$x, 0, parameters$tau, log = TRUE))
  }
  # fit the models
  fit0 <- JAGS_fit(model_syntax,  data, priors_list0, chains = 1, adapt = 100, burnin = 200, sample = 1000, seed = 0)
  fit1 <- JAGS_fit(model_syntax,  data, priors_list1, chains = 1, adapt = 100, burnin = 200, sample = 1000, seed = 1)
  fit2 <- JAGS_fit(model_syntax2, data, priors_list2, chains = 1, adapt = 100, burnin = 200, sample = 1000, seed = 2)
  # get marginal likelihoods
  marglik0 <- list(
    logml = sum(dnorm(data$x, mean(p0), 1, log = TRUE))
  )
  class(marglik0) <- "bridge"
  marglik1 <- JAGS_bridgesampling(fit1, log_posterior = log_posterior, data = data, prior_list = priors_list1)
  marglik2 <- JAGS_bridgesampling(fit2, log_posterior = log_posterior2, data = data, prior_list = priors_list2)
  ## create an object with the models
  models <- list(
    list(fit = fit0, marglik = marglik0, priors = priors_list0, prior_weights = 1, fit_summary = runjags_estimates_table(fit0, priors_list0)),
    list(fit = fit1, marglik = marglik1, priors = priors_list1, prior_weights = 1, fit_summary = runjags_estimates_table(fit1, priors_list1)),
    list(fit = fit2, marglik = marglik2, priors = priors_list2, prior_weights = 1, fit_summary = runjags_estimates_table(fit2, priors_list2))
  )
  # compare and summarize the models
  models            <- models_inference(models)
  inference         <- ensemble_inference(model_list = models, parameters = c("mu", "tau"), is_null_list = list("mu" = c(1, 3), "tau" = c(1, 2)))
  mixed_posteriors  <- mix_posteriors(model_list = models, parameters = c("mu", "tau"), is_null_list = list("mu" = c(1, 3), "tau" = c(1, 2)), seed = 1)


  expect_doppelganger("model-averaging-plot-models-1", function(){
    plot_models(model_list = models, samples = mixed_posteriors, inference = inference, parameter = "mu")
  })
  expect_doppelganger("model-averaging-plot-models-2", function(){
    plot_models(model_list = models, samples = mixed_posteriors, inference = inference, parameter = "tau")
  })
  expect_doppelganger("model-averaging-plot-models-3", function(){
    plot_models(model_list = models, samples = mixed_posteriors, inference = inference, parameter = "mu", prior = TRUE)
  })
  expect_doppelganger("model-averaging-plot-models-4", function(){
    plot_models(model_list = models, samples = mixed_posteriors, inference = inference, parameter = "tau", prior = TRUE)
  })
  expect_doppelganger("model-averaging-plot-models-5", function(){
    plot_models(model_list = models, samples = mixed_posteriors, inference = inference, parameter = "mu", conditional = TRUE)
  })
  expect_doppelganger("model-averaging-plot-models-6", function(){
    plot_models(model_list = models, samples = mixed_posteriors, inference = inference, parameter = "tau", prior = TRUE, conditional = TRUE)
  })
  expect_doppelganger("model-averaging-plot-models-7", {
    plot_models(model_list = models, samples = mixed_posteriors, inference = inference, parameter = "mu", plot_type = "ggplot")
  })
  expect_doppelganger("model-averaging-plot-models-8", {
    plot_models(model_list = models, samples = mixed_posteriors, inference = inference, parameter = "tau", prior = TRUE, plot_type = "ggplot")
  })
  expect_doppelganger("model-averaging-plot-models-9", {
    plot_models(model_list = models, samples = mixed_posteriors, inference = inference, parameter = "mu", plot_type = "ggplot", y_axis2 = FALSE)
  })
  expect_doppelganger("model-averaging-plot-models-10", {
    plot_models(model_list = models, samples = mixed_posteriors, inference = inference, parameter = "mu", plot_type = "ggplot", show_estimates = FALSE)
  })

})
