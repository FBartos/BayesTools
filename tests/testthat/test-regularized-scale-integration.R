skip_if_not_test_profile("unit")

test_that("half-t global shrinkage retains both scale tails in horseshoe marginal methods", {
  ht <- prior("t", list(0, .2, 3), list(0, Inf))
  x <- c(.01, .5, 2)
  # Independent Normal scale mixture. For W = tau*lambda, k = .2*sqrt(3)
  # and q = (w/k)^2, the density of W is 4/(k*pi²)*(q-1-log(q))/(q-1)²
  # by integration of the half-Cauchy density against the half-t(3).
  # The slab additionally mixes the Normal over 16/c² ~ chi-square(4).
  density <- list(c(5.330492525886, .1306424147123, .01245856909167),
                  c(5.350535395482, .1411936102545, .01278937283732))
  tail <- list(c(.4083532640487, .08718056005989, .02667760231420),
               c(.4081521671100, .07973305868269, .01290036659426))
  for(i in 1:2){
    group <- prior_regularized("horseshoe", ht, slab_scale = if(i == 2) 2)
    priors <- JAGS_formula(~ wt, "mu", mtcars, list(intercept = prior("point", list(0)), sh = group))$prior_list
    p <- priors$mu_wt
    expect_equal(pdf(p, x), density[[i]], tolerance = 1e-7)
    expect_equal(ccdf(p, x), tail[[i]], tolerance = 1e-7)
    expect_equal(ccdf(p, quant(p, .99)), .01, tolerance = 1e-7)
    curve <- ggplot2::ggplot_build(plot(p, plot_type = "ggplot", n_points = 64))$data[[1L]]
    expect_true(sum(is.finite(curve$y) & curve$y > 0) >= 60)
  }
})

test_that("shifted truncated global scales retain their probability and first moment", {
  scales <- list(prior("t", list(-.1, .2, 3), list(0, Inf)),
    prior("normal", list(.1, .2), list(0, Inf)),
    prior("laplace", list(-.1, .2), list(0, Inf)),
    prior("t", list(.1, .2, 3), list(.05, 1)))
  for(s in scales){
    rule <- .regularized_global_scale_rule(s, 4L, .regularized_probability_rule(4L), 6)
    expect_true(all(rule$tau > s$truncation$lower & rule$tau < s$truncation$upper))
    expect_equal(sum(exp(rule$log_weight)), 1, tolerance = 1e-8)
    expected <- stats::integrate(function(tau) tau * pdf(s, tau),
      s$truncation$lower, s$truncation$upper, rel.tol = 1e-9)$value
    expect_equal(sum(rule$tau * exp(rule$log_weight)), expected, tolerance = 1e-7)
  }
})
