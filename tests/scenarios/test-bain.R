if (file.exists("helper-scenarios.R")) source("helper-scenarios.R") else source("tests/scenarios/helper-scenarios.R")
scenario_start("bain", width = 220)
# testthat::test_file("tests/scenarios/test-bain.R")

### Description
# Overall ratings in the 30 departments of the attitude survey. Compare
# the printed BayesTools and bain tables for the same Gaussian likelihood and
# restrictions. Match bain's default fractional normal prior on the tested
# coefficients. For these lm fits its covariance is vcov(lm) * n / rank,
# where rank counts independent rows of the full constraint matrix including
# its right-hand side. The interval contributes two distinct bounds.
# JAGS updates the matched prior under the Gaussian likelihood. bain uses
# a normal posterior approximation with the lm estimates and covariance.
# For a point, bain's Fit and Com are posterior and prior densities. For an
# inequality or interval they are posterior and prior probabilities. bain
# estimates interval probabilities by simulation, so even the matched prior
# can give slightly different printed values.

testthat::test_that("Hypothesis reporting for the mean department rating", {

  # One unit is ten percentage points. A mean within 0.1 of zero on this
  # scale corresponds to an overall rating between 64% and 66%.
  departments <- transform(attitude, rating = (rating - 65) / 10)
  mean_lm <- stats::lm(rating ~ 1, data = departments)
  fit_mean <- scenario_fit("mean_rating", {
    # The common boundary (and interval midpoint) is zero. The rank is two.
    mean_sd <- sqrt(stats::vcov(mean_lm)[1, 1] * stats::nobs(mean_lm) / 2)
    JAGS_fit("model { for (i in 1:N) { rating[i] ~ dnorm(mu, pow(sigma, -2)) } }",
             data = list(rating = departments$rating, N = nrow(departments)),
             prior_list = list(mu = prior("normal", list(0, mean_sd)),
                               sigma = prior("normal", list(0, 1), truncation = list(0, Inf))),
             chains = 4, adapt = 1000, burnin = 1000, sample = 5000, seed = 1)
  })

  # Point on the left of vs makes the numerator agree with bain's H1.
  # The two interval rows correspond to bain's BF.u and BF.c for H3.
  hypotheses <- c("mu = 0 vs unconstrained", "mu > 0",
                  "-0.1 < mu < 0.1 vs unconstrained", "-0.1 < mu < 0.1")
  scenario_text("mean-compact", {
    hypothesis_BF(fit_mean, hypotheses)
  })
  scenario_text("mean-reporting", {
    cat("BayesTools\n")
    print(hypothesis_BF(fit_mean, hypotheses, include_densities = TRUE,
                        include_probabilities = TRUE, include_method = TRUE))
    cat("\nbain\n")
    set.seed(1)
    bain::bain(mean_lm, "Intercept = 0; Intercept > 0; -0.1 < Intercept < 0.1",
               standardize = FALSE)
  })
})

testthat::test_that("Hypothesis reporting for adjusted department rating slopes", {

  # Slopes are changes in units of ten rating points per ten predictor points,
  # holding the other predictor fixed. Centering keeps the intercept at the
  # departments' average predictor values.
  departments <- transform(attitude, rating = (rating - 65) / 10,
                           complaints = (complaints - mean(complaints)) / 10,
                           raises = (raises - mean(raises)) / 10)
  regression_lm <- stats::lm(rating ~ complaints + raises, data = departments)
  regression <- scenario_fit("rating_regression", {
    # The hypotheses give rank three and center (complaints, raises) at (0.5, 0).
    # Factor their joint normal prior as p(raises) * p(complaints | raises).
    slopes <- c("complaints", "raises")
    V <- stats::vcov(regression_lm)[slopes, slopes] * stats::nobs(regression_lm) / 3
    loading <- V["complaints", "raises"] / V["raises", "raises"]
    conditional_sd <- sqrt(V["complaints", "complaints"] - loading * V["complaints", "raises"])
    fit <- JAGS_fit("model { for (i in 1:N) { rating[i] ~ dnorm(mu[i], pow(sigma, -2)) } }",
             data = list(rating = departments$rating, N = nrow(departments), loading = loading),
             prior_list = list(sigma = prior("normal", list(0, 1), truncation = list(0, Inf))),
             formula_list = list(mu = ~ complaints + raises),
             formula_data_list = list(mu = departments),
             formula_prior_list = list(mu = list(
               intercept = prior("normal", list(0, 1)),
               raises = prior("normal", list(0, sqrt(V["raises", "raises"]))),
               complaints = prior("normal", list(expression(0.5 + loading * mu_raises), conditional_sd)))),
             chains = 4, adapt = 1000, burnin = 1000, sample = 5000, seed = 1)
    list(fit = fit, complaints_prior = prior("normal", list(0.5, sqrt(V["complaints", "complaints"]))))
  })

  # H1 asks whether a ten-point improvement in complaint handling predicts
  # more than five rating points. H2 tests a zero slope for performance-based
  # raises. H3 allows up to one rating point per ten points on the raises item.
  # Both priors put half their mass above 0.5 for H1. Use the marginal prior
  # for this test because JAGS declares the complaints prior conditionally.
  # Both packages use continuous priors, so the interval endpoints have zero mass.
  scenario_text("regression-reporting", {
    cat("BayesTools\n")
    complaints <- parameter_catalog_resolve(parameter_catalog(regression$fit), "mu_complaints")
    draws <- as.numeric(as.matrix(parameter_draws(regression$fit, complaints)))
    print(hypothesis_BF(draws, regression$complaints_prior, "complaints > 0.5", parameter = "complaints",
                        include_probabilities = TRUE, include_method = TRUE))
    print(hypothesis_BF(regression$fit,
                        c("raises = 0 vs unconstrained",
                          "-0.1 <= raises <= 0.1", "-0.1 <= raises <= 0.1 vs unconstrained"),
                        include_densities = TRUE, include_probabilities = TRUE, include_method = TRUE))
    cat("\nbain\n")
    set.seed(1)
    bain::bain(regression_lm, "complaints > 0.5; raises = 0; -0.1 < raises < 0.1",
               standardize = FALSE)
  })
})
