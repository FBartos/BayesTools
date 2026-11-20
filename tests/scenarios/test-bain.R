if (file.exists("helper-scenarios.R")) source("helper-scenarios.R") else source("tests/scenarios/helper-scenarios.R")
scenario_start("bain", width = 220)
# testthat::test_file("tests/scenarios/test-bain.R")

### Description
# Overall ratings in the 30 departments of the attitude survey. Compare
# the printed BayesTools and bain tables for the same Gaussian likelihood and
# restrictions. BayesTools uses the priors specified below. bain uses its
# default fractional prior and the lm estimates, so BF magnitudes can differ.
# For a point, bain's Fit and Com are posterior and prior densities. For an
# inequality or interval they are posterior and prior probabilities.

testthat::test_that("Hypothesis reporting for the mean department rating", {

  # One unit is ten percentage points. A mean within 0.1 of zero on this
  # scale corresponds to an overall rating between 64% and 66%.
  departments <- transform(attitude, rating = (rating - 65) / 10)
  fit_mean <- scenario_fit("mean_rating", {
    JAGS_fit("model { for (i in 1:N) { rating[i] ~ dnorm(mu, pow(sigma, -2)) } }",
             data = list(rating = departments$rating, N = nrow(departments)),
             prior_list = list(mu = prior("normal", list(0, 1)),
                               sigma = prior("normal", list(0, 1), truncation = list(0, Inf))),
             chains = 4, adapt = 1000, burnin = 1000, sample = 5000, seed = 1)
  })
  mean_lm <- stats::lm(rating ~ 1, data = departments)

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
  fit_regression <- scenario_fit("rating_regression", {
    coefficient <- prior("normal", list(0, 1))
    JAGS_fit("model { for (i in 1:N) { rating[i] ~ dnorm(mu[i], pow(sigma, -2)) } }",
             data = list(rating = departments$rating, N = nrow(departments)),
             prior_list = list(sigma = prior("normal", list(0, 1), truncation = list(0, Inf))),
             formula_list = list(mu = ~ complaints + raises),
             formula_data_list = list(mu = departments),
             formula_prior_list = list(mu = list(intercept = coefficient, complaints = coefficient, raises = coefficient)),
             chains = 4, adapt = 1000, burnin = 1000, sample = 5000, seed = 1)
  })
  regression_lm <- stats::lm(rating ~ complaints + raises, data = departments)

  # H1 asks whether a ten-point improvement in complaint handling predicts
  # more than five rating points. H2 tests a zero slope for performance-based
  # raises. H3 allows up to one rating point per ten points on the raises item.
  # H1 illustrates the prior difference: bain puts half its prior above the
  # boundary, while the BayesTools N(0, 1) slope prior puts about 31% above 0.5.
  # Both packages use continuous priors, so the interval endpoints have zero mass.
  scenario_text("regression-reporting", {
    cat("BayesTools\n")
    print(hypothesis_BF(fit_regression,
                        c("complaints > 0.5", "raises = 0 vs unconstrained",
                          "-0.1 <= raises <= 0.1", "-0.1 <= raises <= 0.1 vs unconstrained"),
                        include_densities = TRUE, include_probabilities = TRUE, include_method = TRUE))
    cat("\nbain\n")
    set.seed(1)
    bain::bain(regression_lm, "complaints > 0.5; raises = 0; -0.1 < raises < 0.1",
               standardize = FALSE)
  })
})
