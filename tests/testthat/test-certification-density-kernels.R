skip_if_not_test_profile("certification")

# Expensive reference sweeps of prior densities, ordinates, region grids, and
# quantile estimators, moved out of the unit lane. Each test here is the full
# sweep behind a cheaper representative check that stays in the unit test of
# the same behaviour (named in the comments); the sweeps compare with analytic
# or independently integrated references at the tolerances of those tests.

# The product sweeps of "scalar affine ordinates at a bound of the mapped
# support are the one-sided limit inside it" (test-prior-density-ordinate.R),
# whose unit test keeps one (w, offset) per source, two log combinations, the
# three named rounded combinations, the ulp loop, and the mixture.
test_that("scalar affine ordinates at a bound of the mapped support are the one-sided limit inside it (full sweep)", {

  # offset + w S with a scalar prior S and a point term o: the mapped bounds
  # are rounded and their inverse image (o + w - o) / w can round just outside
  # the support of S (above 1 for S ~ U(0, 1) with (w, o) = (.1, .2), (.3, .1),
  # (.3, .7)), where the ordinate was 0. The ordinate at a bound is the
  # one-sided limit inside the support, the density of S at its bound divided
  # by |w| (analytic references); outside the support it stays 0
  ulp <- function(x) 2^(floor(log2(abs(x))) - 52)
  affine <- function(source, w, o, source_transform = NULL){
    BayesTools:::.prior_linear_combination_density(
      prior_list        = list(s = source, p = prior("point", list(o))),
      weights           = c(s = w, p = 1),
      source_transforms = if(!is.null(source_transform)) c(s = source_transform, p = NA_character_),
      n_grid            = 64
    )
  }
  # the density of S at each finite source bound (the log source Z = log(S)
  # has the density f_S(s) s at log(s))
  normal_mass <- stats::pnorm(2) - stats::pnorm(-1)
  sources <- list(
    list(prior = prior("uniform", list(0, 1)), bounds = c(0, 1), limits = c(1, 1)),
    list(prior = prior("uniform", list(1, 3)), bounds = c(1, 3), limits = c(.5, .5)),
    list(prior = prior("beta", list(1, 3)), bounds = c(0, 1), limits = c(3, 0)),
    list(prior = prior("beta", list(3, 1)), bounds = c(0, 1), limits = c(0, 3)),
    list(prior = prior("normal", list(0, 1), list(-1, 2)), bounds = c(-1, 2),
         limits = stats::dnorm(c(-1, 2)) / normal_mass),
    list(prior = prior("exp", list(2)), bounds = 0, limits = 2)
  )
  expect_bound_limits <- function(source, w, o, log_source = FALSE){
    density <- affine(source$prior, w, o, if(log_source) "log")
    positive <- source$bounds > 0
    bounds <- if(log_source) log(source$bounds[positive]) else source$bounds
    limits <- source$limits[if(log_source) positive else TRUE]
    if(log_source) limits <- limits * source$bounds[positive]
    mapped <- o + w * bounds
    for(i in seq_along(bounds)){
      value <- mapped[i]
      label <- sprintf("%s, w = %g, o = %g, bound %g", source$prior$distribution, w, o, bounds[i])
      ordinate <- prior_density_ordinate(density, value)
      expect_identical(ordinate$method, "scalar_affine", info = label)
      expect_true(ordinate$exact, info = label)
      if(limits[i] > 0){
        expect_identical(ordinate$behavior, "regular", info = label)
        expect_equal(exp(ordinate$log_density), limits[i] / abs(w), tolerance = 1e-12, info = label)
      }else{
        expect_identical(ordinate$behavior, "zero", info = label)
        expect_identical(ordinate$log_density, -Inf, info = label)
      }
      # outside the support (4 ulps of the largest operand of the mapping, which
      # the rounding of the inverse image is in, and 1e-9 beyond the mapped
      # bound) the ordinate is 0; a one-sided bound has its outside on one side
      # only
      outward <- if(length(bounds) == 2L) {
        if(value == max(mapped)) 1 else -1
      }else{
        if(w > 0) -1 else 1
      }
      operands <- max(abs(c(value, o, w * bounds[i])))
      for(beyond in value + outward * c(1e-9 * max(1, abs(value)), 4 * ulp(operands))){
        outside <- prior_density_ordinate(density, beyond)
        expect_identical(outside$behavior, "zero", info = label)
        expect_identical(outside$log_density, -Inf, info = label)
        expect_true(outside$exact, info = label)
      }
    }
  }

  # S ~ U(0, 1): both bounds of 32 combinations (the three rounded ones among
  # them), and each source on a subset
  for(w in c(.1, .3, .7, 1.3, -.1, -.3, -.7, -1.3)){
    for(o in c(.1, .2, .7, -.35)){
      expect_bound_limits(sources[[1L]], w, o)
    }
  }
  for(source in sources[-1L]){
    for(w in c(.1, .3, 1.3, -.3, -.7)){
      for(o in c(.15, .2, .7)){
        expect_bound_limits(source, w, o)
      }
    }
  }

  # a log source term: the bounds of S ~ U(5, 9) and U(1, 3) are rounded by
  # exp(log(s)) as well (it is below 5 and above 9 and 3)
  for(source in list(
    list(prior = prior("uniform", list(5, 9)), bounds = c(5, 9), limits = c(.25, .25)),
    sources[[2L]]
  )){
    for(w in c(.1, .3, .7, 1.3, -.1, -.3, -.7, -1.3)){
      for(o in c(.1, .2, .7, -.35)){
        expect_bound_limits(source, w, o, log_source = TRUE)
      }
    }
  }
})

# The 201-threshold sweep of "a divided linear statement compiles to the
# correctly rounded threshold" (test-hypothesis-BF.R), whose unit test keeps the
# 13 thresholds whose quotient does not multiply back to the numerator.
.hypothesis_factor_marginal_for_test <- function(sd = 1, independent = FALSE){
  formula <- JAGS_formula(if(independent) ~0+fac else ~fac, "mu", data.frame(fac = factor(c("A", "B", "C"))),
    prior_list = if(independent) list(fac = prior_factor("normal", list(0, sd), contrast = "independent")) else
      list(intercept = prior("normal", list(0, 1)), fac = prior_factor("normal", list(0, sd), contrast = "treatment")))
  columns <- .JAGS_prior_factor_names("mu_fac", formula$prior_list$mu_fac)
  posterior <- matrix(unlist(lapply(seq_along(columns), function(i) seq(-i, i, length.out = 201L))), 201L, dimnames = list(NULL, columns))
  if(!independent) posterior <- cbind(mu_intercept = seq(-1, 1, length.out = 201L), posterior)
  fit <- structure(coda::mcmc(posterior), class = c("mcmc", "BayesTools_fit"))
  attr(fit, "prior_list") <- formula$prior_list
  fit <- attach_test_parameter_map(fit)
  marginal_posterior(as_mixed_posteriors(fit, "mu_fac"), "mu_fac", use_formula = FALSE, prior_samples = TRUE)
}

test_that("a divided linear statement compiles to the correctly rounded threshold (full sweep)", {

  levels <- .hypothesis_factor_marginal_for_test(independent = TRUE)
  shifted <- function(threshold) paste0("(mu_fac[A] + mu_fac[B]) / 3 + 0.25 > ", threshold)

  # every two-decimal threshold of the statement compiles (the round-trip check
  # refused 13 of these 201 ordinary statements)
  thresholds <- round(seq(-1, 1, by = 0.01), 2)
  refused <- vapply(thresholds, function(threshold){
    inherits(tryCatch(hypothesis_linear_target(levels, shifted(threshold), "mu_fac"), error = identity), "error")
  }, logical(1))
  expect_identical(thresholds[refused], numeric())
})

# The 10,000-draw variant of the pre-registered accuracy study of
# test-harrell-davis-quantile.R, whose unit test keeps the 1,000-draw variant.
# Pre-registered accuracy study (plan hd-bands): pointwise bands of lines
# y = b0 + b1 * x with b0 ~ N(0, 1) and b1 ~ N(0, .5^2), whose pointwise
# quantile is analytic, at p = .025 over a uniform grid of 41 points. The
# roughness of the error is the RMS of its second differences. Medians of the
# paired ratios over 200 seeded replicates: roughness of the empirical
# quantile over that of the HD estimate (at least 4 for 1,000 and 2 for
# 10,000 draws) and RMSE of the HD estimate over that of the empirical
# quantile (at most 1.02).
hd_accuracy_study <- function(m, minimum_roughness_ratio) {

  roughness <- function(error) sqrt(mean(diff(error, differences = 2)^2))
  rmse      <- function(error) sqrt(mean(error^2))
  grid      <- seq(-1, 1, length.out = 41)
  p         <- .025
  truth     <- stats::qnorm(p, 0, sqrt(1 + .25 * grid^2))

  set.seed(31, kind = "Mersenne-Twister", normal.kind = "Inversion")
  ratios <- vapply(seq_len(200), function(i) {
    b0  <- stats::rnorm(m)
    b1  <- stats::rnorm(m, 0, .5)
    y   <- outer(b0, rep(1, length(grid))) + outer(b1, grid)
    emp <- apply(y, 2, stats::quantile, probs = p, names = FALSE)
    hd  <- harrell_davis_quantile(y, p)[1, ]
    c(roughness = roughness(emp - truth) / roughness(hd - truth),
      rmse      = rmse(hd - truth) / rmse(emp - truth))
  }, numeric(2))

  expect_gte(stats::median(ratios["roughness", ]), minimum_roughness_ratio)
  expect_lte(stats::median(ratios["rmse", ]), 1.02)
}

test_that("harrell_davis_quantile attenuates the kinks of the empirical quantile without losing accuracy (10,000 draws)", {

  hd_accuracy_study(10000, 2)
})

# The default display grid of 1,000 values of the ordered direct displays of
# test-prior-ordered-direct-display.R.
test_that("ordered direct displays of a total away from zero follow the structural route (default grid)", {

  p <- ordered_plot_test_fixture(prior("uniform", list(1, 2)), levels = c("a", "b", "c", "d"))$prior
  # the default grid has no zero height inside the support and each level keeps its unit mass
  default <- density(p)
  for(level in 2:3){
    inside <- default[[level]]$x > 0 & default[[level]]$x < 2
    expect_true(all(default[[level]]$y[inside] > 0))
    grid <- default[[level]]$x
    height <- default[[level]]$y
    expect_equal(sum(diff(grid) * (head(height, -1L) + tail(height, -1L)) / 2), 1, tolerance = 1e-4)
  }
})

test_that("unresolved ordered display values are reported, not dropped silently (default grid)", {

  # a lognormal total with a Dirichlet(1, 1, .05) allocation: the quadrature of the third level does not
  # converge and its structural route resolves 1 of 200 display values; the other levels are complete
  p <- ordered_plot_test_fixture(prior("lognormal", list(0, 1)), prior("dirichlet", list(c(1, 1, .05))),
    levels = c("a", "b", "c", "d"))$prior
  # the default display of 1000 values resolves every value of the level and does not warn
  default <- NULL
  expect_no_warning(default <- density(p))
  expect_length(default[[3L]]$x, 1000L)
})

# The cross-checks of the vectorized PET-PEESE bands with the scalar quantile at
# the grid points 2, 37 and 101 (test-model-averaging-plots-edge-cases.R keeps
# points 1 and 11 and the independent oracles).
test_that("vectorized PET-PEESE prior bands reproduce the scalar quantiles (full sweep)", {

  probs <- c(.5, .025, .975)
  check_band <- function(context, x, check){
    band <- .petpeese_prior_quantile_band(context, x, probs)
    expect_identical(dim(band), c(3L, length(x)))
    expect_true(all(is.finite(band)))
    for(i in check){
      expect_equal(band[, i], .petpeese_prior_cdf_quantile(context, x[i], probs), tolerance = 1e-4)
      if(x[i] > 0){
        # The exact CDF at the band quantile meets the band's 1e-7 criterion.
        for(j in seq_along(probs)){
          expect_lte(abs(as.numeric(.petpeese_prior_cdf(context, band[j, i], x[i])) - probs[j]), 2e-7)
        }
      }
    }
    band
  }

  # RoBMA-like ensemble: no-bias, PET and PEESE models over a spike-and-normal
  # effect, with half-Cauchy bias priors.
  mu <- prior_mixture(list(prior("point", list(0), prior_weights = 1), prior("normal", list(0, .7))),
    is_null = c(TRUE, FALSE))
  bias <- list(prior_none(prior_weights = 2), prior_PET("cauchy", list(0, 1), prior_weights = 1),
    prior_PEESE("cauchy", list(0, 5), prior_weights = 1))
  context <- .petpeese_prior_cdf_context(bias, list(mu, mu, mu), c(2, 1, 1) / 4)
  x <- seq(0, 1, length.out = 101)
  check_band(context, x, c(2L, 37L, 101L))
})

# The ten two-term convolution builds that the unit test of "two-term
# convolution regions integrate over the term without a 1/x^2 tail"
# (test-hypothesis-prior-region-grid.R) leaves out: it builds the six
# gamma-partner terms, the Cauchy term with a beta partner (heavy first), and
# the Cauchy term with a uniform partner (heavy second).
test_that("two-term convolution regions integrate over the term without a 1/x^2 tail (remaining builds)", {

  # a 1/x^2-tail term (Cauchy, half-Cauchy, inverse gamma with shape 1) with a
  # gamma, beta or uniform partner: the region probability integrates the
  # partner's density times the heavy term's exact region probability in both
  # list orders (the heavy-tailed first term hit the QUADPACK roundoff flag)
  heavy <- list(
    cauchy = list(
      prior = prior("cauchy", list(0, 1)),
      lower = function(x) stats::pcauchy(x),
      upper = function(x) stats::pcauchy(x, lower.tail = FALSE)
    ),
    halfcauchy = list(
      prior = prior("cauchy", list(0, 1), list(0, Inf)),
      lower = function(x) ifelse(x <= 0, 0, 2 * stats::pcauchy(x) - 1),
      upper = function(x) ifelse(x <= 0, 1, 2 * stats::pcauchy(x, lower.tail = FALSE))
    ),
    invgamma = list(
      prior = prior("invgamma", list(1, .5)),
      lower = function(x) ifelse(x <= 0, 0, stats::pgamma(1 / pmax(x, 1e-300), 1, rate = .5, lower.tail = FALSE)),
      upper = function(x) ifelse(x <= 0, 1, stats::pgamma(1 / pmax(x, 1e-300), 1, rate = .5))
    )
  )
  partners <- list(
    gamma   = list(prior = prior("gamma", list(2, 1)), density = function(g) stats::dgamma(g, 2, 1), range = c(0, Inf)),
    beta    = list(prior = prior("beta", list(2, 2)), density = function(g) stats::dbeta(g, 2, 2), range = c(0, 1)),
    uniform = list(prior = prior("uniform", list(0, 1)), density = function(g) stats::dunif(g, 0, 1), range = c(0, 1))
  )
  regions <- list(
    "theta > 0" = function(term, g) term$upper(-g),
    "theta < 1" = function(term, g) term$lower(1 - g),
    "theta < -1 | theta > 1" = function(term, g) term$lower(-1 - g) + term$upper(1 - g)
  )
  reference <- function(term, partner, region){
    knots <- unique(c(partner$range[1L], if(partner$range[1L] < 1 && 1 < partner$range[2L]) 1, partner$range[2L]))
    sum(vapply(seq_len(length(knots) - 1L), function(i){
      stats::integrate(function(g) partner$density(g) * region(term, g),
                       knots[i], knots[i + 1L], rel.tol = 1e-13, subdivisions = 1000L)$value
    }, numeric(1)))
  }

  for(term_name in names(heavy)){
    for(partner_name in names(partners)){
      for(heavy_first in c(TRUE, FALSE)){
        in_unit <- partner_name == "gamma" ||
          (term_name == "cauchy" && partner_name == "beta" && heavy_first) ||
          (term_name == "cauchy" && partner_name == "uniform" && !heavy_first)
        if(in_unit) next
        terms <- list(a = heavy[[term_name]]$prior, b = partners[[partner_name]]$prior)
        density <- .prior_linear_combination_density(
          if(heavy_first) terms else terms[c("b", "a")], c(a = 1, b = 1)
        )
        route <- .prior_density_route_from_adaptive(attr(density, "adaptive_evaluation", exact = TRUE))
        expect_identical(route$type, "convolution")
        for(hypothesis in names(regions)){
          label <- paste(term_name, partner_name, if(heavy_first) "heavy first" else "heavy second", hypothesis)
          probability <- tryCatch(
            as.numeric(.prior_linear_density_region_probability(
              density, .hypothesis_prior_region(str2lang(hypothesis), "theta")
            )),
            error = function(e) conditionMessage(e)
          )
          expect_equal(probability, reference(heavy[[term_name]], partners[[partner_name]], regions[[hypothesis]]),
                       tolerance = 1e-8, label = label)
        }
      }
    }
  }
})
