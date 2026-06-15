skip_if_not_test_profile("fit")

expect_nonlocal_prior_only_samples <- function(prior, samples, tolerance = .08) {
  expect_true(all(is.finite(samples)))
  expect_true(all(samples >= prior$truncation[["lower"]]))
  expect_true(all(samples <= prior$truncation[["upper"]]))

  probs <- c(.1, .25, .5, .75, .9)
  quantiles <- quant(prior, probs)
  sample_cdf <- vapply(quantiles, function(x) mean(samples <= x), numeric(1))
  expect_equal(sample_cdf, probs, tolerance = tolerance)

  prior_mean <- mean(prior)
  prior_sd <- sd(prior)
  if (is.finite(prior_mean) && is.finite(prior_sd)) {
    expect_equal(mean(samples), prior_mean, tolerance = tolerance)
    expect_equal(stats::sd(samples), prior_sd, tolerance = tolerance)
  }
}

test_that("BayesTools JAGS module initializes truncated nonlocal priors", {
  skip_if_not_installed("rjags")
  skip_on_cran()

  skip_if_not(isTRUE(BayesTools_load_JAGS_module(quiet = TRUE, warn = FALSE)))

  priors <- list(
    moment = prior(
      "moment",
      list(mode = .5),
      truncation = list(lower = -.1, upper = .1)
    ),
    invmoment = prior(
      "invmoment",
      list(mode = .5, df = 6),
      truncation = list(lower = -.1, upper = .1)
    )
  )

  for(prior_i in priors){
    syntax <- JAGS_add_priors("model{}", list(theta = prior_i))
    expect_silent(local({
      con <- textConnection(syntax)
      on.exit(close(con), add = TRUE)
      rjags::jags.model(
        file     = con,
        data     = list(),
        n.chains = 1,
        n.adapt  = 0,
        quiet    = TRUE
      )
    }))
  }
})

test_that("BayesTools JAGS module samples nonlocal priors", {
  skip_if_not_installed("rjags")
  skip_if_not_installed("runjags")
  skip_if_not_installed("bridgesampling")
  skip_on_cran()

  skip_if_not(isTRUE(BayesTools_load_JAGS_module(quiet = TRUE, warn = FALSE)))

  priors <- list(
    moment = prior("moment", list(mode = .5)),
    invmoment = prior("invmoment", list(mode = .5, df = 6)),
    moment_truncated = prior("moment", list(mode = .5), truncation = list(lower = -Inf, upper = 0)),
    invmoment_truncated = prior("invmoment", list(mode = .5, df = 6), truncation = list(lower = -Inf, upper = 0))
  )

  for(prior_name in names(priors)){
    prior_list <- list(theta = priors[[prior_name]])
    fit <- suppressWarnings(JAGS_fit(
      model_syntax = "model{}",
      data = NULL,
      prior_list = prior_list,
      chains = 2,
      adapt = 250,
      burnin = 250,
      sample = 4000,
      silent = TRUE,
      seed = 1
    ))

    expect_s3_class(fit, "BayesTools_fit")
    expect_true("BayesTools" %in% attr(fit, "jags_modules"))
    samples <- as.matrix(fit$mcmc)
    expect_true("theta" %in% colnames(samples))
    expect_nonlocal_prior_only_samples(priors[[prior_name]], samples[, "theta"])

    marglik <- JAGS_bridgesampling(
      fit = fit,
      log_posterior = STANDARD_LOG_POSTERIOR,
      data = list(),
      prior_list = prior_list,
      maxiter = 2000
    )
    expect_s3_class(marglik, "bridge")
    expect_equal(marglik$logml, 0, tolerance = .08)
  }
})
