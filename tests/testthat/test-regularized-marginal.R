test_that("regularization declarations distinguish fixed widths and learned scales", {
  expect_true(is.prior.simple(prior_regularized("lasso", 0.4)))
  expect_equal(pdf(prior_regularized("ridge", 0.4), 1), dnorm(1, sd = 0.4))
  hs <- prior_regularized("horseshoe", 0.5, terms = "wt", slab_scale = 2)
  expect_true(is.prior.regularized(hs))
  expect_match(print(hs, silent = TRUE), "Regularized horseshoe.*slab df = 4.*wt")
  expect_error(prior_regularized("horseshoe", prior("normal", list(0, 1))), "positive support")
  expect_error(prior_regularized("ridge", 1, terms = "intercept"), "non-intercept")
  expect_error(prior_regularized("ridge", Inf), "finite")
  expect_error(prior_regularized("lasso", 1, slab_scale = 1), "only to the horseshoe")
  expect_error(prior_regularized("horseshoe", 1, slab_df = 7), "Specify 'slab_scale'")
  expect_error(rng(hs, 3), "named 'prior_list'")
})

test_that("learned scale marginal ordinates agree with independent integration", {
  bound <- function(family, scale, ...){
    p <- JAGS_formula(~ wt, "mu", mtcars, list(
      intercept = prior("point", list(0)), sh = prior_regularized(family, scale, ...)))
    p$prior_list[["mu_wt"]]
  }
  global <- prior("gamma", list(2, 1))
  ridge <- bound("ridge", global)
  lasso <- bound("lasso", global)
  for(value in c(0, 0.3, 2)){
    normal <- integrate(function(t) dnorm(value, sd = t) * dgamma(t, 2), 0, Inf)
    laplace <- integrate(function(t) exp(-value/t) / (2*t) * dgamma(t, 2), 0, Inf)
    expect_equal(pdf(ridge, value), normal$value, tolerance = 1e-6)
    expect_equal(pdf(lasso, value), laplace$value, tolerance = 1e-6)
  }
  expect_equal(quant(ridge, c(0, 0.5, 1)), c(-Inf, 0, Inf))
  expect_equal(cdf(lasso, quant(lasso, .82)), .82, tolerance = 1e-6)
  expect_equal(var(ridge), 6)
  expect_equal(var(lasso), 12)
  expect_equal(prior_density_ordinate(bound("ridge", prior("exp", list(1))), 0)$behavior,
               "infinite")
  hs <- bound("horseshoe", 0.35)
  expect_equal(pdf(hs, 0), Inf)
  expected <- integrate(function(t) dnorm(1, sd = 0.35 * t) * 2 * dcauchy(t), 0, Inf)
  expect_equal(pdf(hs, 1), expected$value, tolerance = 1e-6)
  finite <- bound("horseshoe", 0.35, slab_scale = 1)
  expect_equal(prior_density_ordinate(finite, 0)$behavior, "infinite")
  # Integrate the local scale and the inverse-chi-square slab independently.
  expected <- integrate(function(u) vapply(u, function(v){
    integrate(function(t) pnorm(-0.8 * sqrt(1 / (0.35*t)^2 + v/4)) * 2 * dcauchy(t),
              0, Inf)$value * dchisq(v, 4)
  }, numeric(1)), 0, Inf)
  expect_equal(ccdf(finite, 0.8), expected$value, tolerance = 1e-5)
})
