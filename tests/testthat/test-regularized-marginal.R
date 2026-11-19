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

test_that("random-global finite-slab horseshoe agrees with the direct scale-mixture integral", {
  global <- prior("cauchy", list(0, 0.2), list(0, Inf))
  hs <- JAGS_formula(~ wt, "mu", mtcars, list(
    intercept = prior("point", list(0)),
    shrink = prior_regularized("horseshoe", global, slab_scale = 2, slab_df = 4)
  ))$prior_list$mu_wt

  # For independent standard half-Cauchys the density of the log product
  # z is 2*z/(pi^2*sinh(z)). Mix the Normal over this coordinate and the
  # independent chi-square(4) coordinate of c^2 = 16 / chi-square(4).
  reference <- function(x, density){
    integrand <- function(z){
      mass <- ifelse(z == 0, 2 / pi^2, 2 * abs(z) / (pi^2 * sinh(abs(z))))
      vapply(seq_along(z), function(i){
        if(mass[i] == 0) return(0)
        stats::integrate(function(v){
          sd <- 1 / sqrt(25 * exp(-2 * z[i]) + v / 16)
          conditional <- if(density) stats::dnorm(x, sd = sd) else stats::pnorm(-x, sd = sd)
          conditional * stats::dchisq(v, 4)
        }, 0, Inf, rel.tol = 1e-8)$value * mass[i]
      }, numeric(1))
    }
    cuts <- sort(unique(c(-Inf, log(x / .2), 0, Inf)))
    sum(vapply(seq_len(length(cuts) - 1L), function(i){
      stats::integrate(integrand, cuts[i], cuts[i + 1L], rel.tol = 1e-7)$value
    }, numeric(1)))
  }
  x <- c(1e-10, 0.001, 1.4)
  expect_equal(pdf(hs, x), vapply(x, reference, numeric(1), density = TRUE), tolerance = 2e-6)
  expect_equal(ccdf(hs, x), vapply(x, reference, numeric(1), density = FALSE), tolerance = 2e-6)
  # Integrating the two half-Cauchys gives c² E[tau/(tau+c)].
  variance <- stats::integrate(function(v){
    c2 <- 16 / v
    a <- sqrt(c2) / 0.2
    c2 * (1 + 2 * a * log(a) / pi) / (1 + a^2) * stats::dchisq(v, 4)
  }, 0, Inf, rel.tol = 1e-10)$value
  expect_equal(var(hs), variance, tolerance = 1e-6)
  expect_equal(ccdf(hs, quant(hs, 0.99)), .01, tolerance = 1e-6)
})

test_that("ridge and lasso resolve the peak under a learned half-Cauchy width", {
  # In u=log(tau/0.2) coordinates the Gaussian or Laplace kernel's 1/tau
  # cancels the Jacobian. Integrate on both sides of the kernel transition.
  x <- 1e-10
  cuts <- c(-Inf, log(x / .2), 0, Inf)
  for(family in c("ridge", "lasso")){
    model <- JAGS_formula(~ wt, "mu", mtcars, list(intercept = prior("point", list(0)),
      sh = prior_regularized(family, prior("cauchy", list(0, .2), list(0, Inf)))))
    integrand <- function(u){
      z <- exp(log(x / .2) - u)
      kernel <- if(family == "ridge") stats::dnorm(z) else exp(-z) / 2
      kernel * (2 / (.2 * pi)) * stats::plogis(-2 * u)
    }
    exact <- sum(vapply(seq_len(3), function(i){
      stats::integrate(integrand, cuts[i], cuts[i + 1L], rel.tol = 1e-10)$value
    }, numeric(1)))
    expect_equal(pdf(model$prior_list$mu_wt, x), exact, tolerance = 1e-7)
  }
})
test_that("multi-column shrinkage terms distinguish joint and marginal methods", {
  group <- prior_regularized("ridge", prior("gamma", list(2, 1)))
  data <- data.frame(g = factor(rep(letters[1:3], 2)))
  term <- JAGS_formula(~ g, "mu", data, list(intercept = prior("point", list(0)),
                                            shrinkage = group))$prior_list$mu_g
  expect_error(pdf(term, c(0, 1)), class = "BayesTools_regularized_joint_unavailable")
  expect_equal(mpdf(term, 0), dnorm(0))
  expect_equal(dim(rng(term, 5, transform_factor_samples = FALSE)), c(5L, 2L))
  levels <- rng(term, 5)
  expect_equal(dim(levels), c(5L, 3L))
  expect_true(all(levels[, 1L] == 0))
})

