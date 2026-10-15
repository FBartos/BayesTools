skip_if_not_test_profile("unit")

test_that("factor joint quantiles fail clearly instead of returning a function", {

  p <- prior_factor("mnormal", list(mean = 0, sd = 1), contrast = "orthonormal")
  expect_error(
    quant(p, .5),
    "Joint quantiles are unavailable for factor priors. Use 'mquant()' for marginal quantiles.",
    fixed = TRUE
  )
})

test_that("unsupported prior classes stop explicitly in every distribution method", {

  selection <- prior_weightfunction("one-sided", c(.025), wf_fixed(c(1, .5)))
  phacking  <- prior_phacking(form = "linear")
  priors <- list(
    ordered  = prior_ordered(prior("normal", list(0, 1))),
    none     = prior_none(),
    phacking = phacking,
    bias     = prior_bias(selection, phacking)
  )
  kinds <- c(ordered = "ordered priors", none = "'prior_none()' priors")
  methods <- list(
    "random generation"          = function(p) rng(p, 3),
    "cdf"                        = function(p) cdf(p, .5),
    "ccdf"                       = function(p) ccdf(p, .5),
    "lpdf"                       = function(p) lpdf(p, .5),
    "lpdf"                       = function(p) pdf(p, .5),
    "quantile function"          = function(p) quant(p, .5),
    "mcdf"                       = function(p) mcdf(p, .5),
    "mccdf"                      = function(p) mccdf(p, .5),
    "mlpdf"                      = function(p) mlpdf(p, .5),
    "mlpdf"                      = function(p) mpdf(p, .5),
    "marginal quantile function" = function(p) mquant(p, .5),
    "mean"                       = function(p) mean(p),
    "variance"                   = function(p) var(p),
    "variance"                   = function(p) sd(p),
    "range"                      = function(p) range(p),
    "density"                    = function(p) density(p)
  )
  supported <- list(
    ordered  = c("random generation", "range", "density"),
    phacking = "random generation",
    bias     = "random generation"
  )

  for(prior_name in names(priors)){
    for(i in seq_along(methods)){
      method <- names(methods)[i]
      if(method %in% supported[[prior_name]]){
        next
      }
      message <- if(prior_name %in% names(kinds)){
        paste0("No ", method, " is implemented for ", kinds[[prior_name]], ".")
      }else if(identical(prior_name, "phacking")){
        paste0("No ", if(identical(method, "quantile function")) "quantile functions" else method,
               " is implemented for p-hacking priors")
      }else{
        paste0("No ", if(identical(method, "quantile function")) "quantile functions" else method,
               " is implemented for composed bias priors")
      }
      expect_error(methods[[i]](priors[[prior_name]]), message, fixed = TRUE)
    }
  }

  mixtures <- list(
    spike_and_slab = prior_spike_and_slab(prior("normal", list(0, 1))),
    mixture        = prior_mixture(list(prior("point", list(0)), prior("normal", list(0, 1))))
  )
  expect_error(range(mixtures$spike_and_slab),
               "No range is implemented for spike and slab priors.", fixed = TRUE)
  expect_error(range(mixtures$mixture),
               "No range is implemented for prior mixtures.", fixed = TRUE)
  expect_error(density(mixtures$mixture),
               "No density is implemented for prior mixtures.", fixed = TRUE)
  expect_error(density(prior("mnormal", list(mean = 0, sd = 1, K = 2))),
               "No density is implemented for this prior distribution.", fixed = TRUE)
})

# ============================================================================ #
# TEST FILE: Prior Distribution Coverage Tests
# ============================================================================ #
#
# PURPOSE:
#   Targeted tests to ensure code coverage for edge cases and error paths
#   in R/priors.R that are not covered by main test-priors.R
#
# DEPENDENCIES:
#   - No external packages required beyond testthat
#
# SKIP CONDITIONS:
#   - None (fast, pure R tests)
#
# MODELS/FIXTURES:
#   - None required
#
# TAGS: @coverage, @edge-cases, @priors
# ============================================================================ #


# ============================================================================ #
# SECTION: prior() constructor errors
# ============================================================================ #

test_that("prior() rejects unknown distribution names", {
  expect_error(prior("unknown_dist", list(0, 1)),
               "The specified distribution name")
})


# ============================================================================ #
# SECTION: prior_factor() contrast validation
# ============================================================================ #

test_that("prior_factor() requires multivariate prior for orthonormal/meandif contrasts", {
  # orthonormal contrast requires multivariate distribution
  expect_error(prior_factor("normal", list(0, 1), contrast = "orthonormal"),
               "contrasts require multivariate prior")

  # meandif contrast requires multivariate distribution
  expect_error(prior_factor("normal", list(0, 1), contrast = "meandif"),
               "contrasts require multivariate prior")

  # bernoulli is not a valid multivariate distribution
  expect_error(prior_factor("bernoulli", list(0.5), contrast = "orthonormal"),
               "contrasts require multivariate prior")
})

test_that("mean-difference and orthonormal factor priors require exact zero centering", {
  for(contrast in c("meandif", "orthonormal")){
    expect_error(
      prior_factor("mnormal", list(mean = 1e-12, sd = 1), contrast = contrast),
      "centered exactly at zero",
      fixed = TRUE
    )
    expect_error(
      prior_factor("mt", list(location = -1, scale = 1, df = 3), contrast = contrast),
      "centered exactly at zero",
      fixed = TRUE
    )
    expect_error(
      prior_factor("point", list(location = 1), contrast = contrast),
      "centered exactly at zero",
      fixed = TRUE
    )
  }

  expect_no_error(
    prior_factor("mnormal", list(mean = 0, sd = 1), contrast = "orthonormal")
  )
  expect_no_error(
    prior_factor("normal", list(mean = 1, sd = 1), contrast = "treatment")
  )
  expect_no_error(
    prior_factor("point", list(location = 1), contrast = "independent")
  )
})


test_that("prior_factor() requires univariate prior for treatment contrast", {
  expect_error(prior_factor("mnormal", list(0, 1, 2), contrast = "treatment"),
               "contrasts require univariate prior")
})


# ============================================================================ #
# SECTION: spike_and_slab prior construction and helpers
# ============================================================================ #

test_that("prior_spike_and_slab() works with factor priors", {
  p_factor <- prior_factor("mnormal", list(mean = 0, sd = 1), contrast = "orthonormal")
  p_factor$parameters[["K"]] <- 2

  p_ss <- prior_spike_and_slab(
    prior_parameter = p_factor,
    prior_inclusion = prior("beta", list(1, 1))
  )

  expect_true(is.prior.spike_and_slab(p_ss))
  expect_s3_class(p_ss, "prior.spike_and_slab")
})


test_that("prior_spike_and_slab() preserves model prior weights", {
  p_ss <- prior_spike_and_slab(
    prior_parameter = prior("normal", list(0, 1)),
    prior_inclusion = prior("point", list(0.25)),
    prior_weights = 7
  )
  p_plain <- prior("normal", list(0, 1), prior_weights = 2)

  expect_equal(BayesTools:::.prior_model_weight(p_ss), 7)
  expect_equal(sapply(list(p_ss, p_plain), BayesTools:::.prior_model_weight), c(7, 2))

  p_ss <- BayesTools:::.set_prior_model_weight(p_ss, 3)
  expect_equal(BayesTools:::.prior_model_weight(p_ss), 3)
  expect_false("prior_weights" %in% names(p_ss))
})

test_that("prior_spike_and_slab() requires scalar inclusion probability priors", {
  expect_error(
    prior_spike_and_slab(
      prior_parameter = prior("normal", list(0, 1)),
      prior_inclusion = prior("mpoint", list(0.5, 2))
    ),
    "scalar probability prior"
  )
  expect_error(
    prior_spike_and_slab(
      prior_parameter = prior("normal", list(0, 1)),
      prior_inclusion = prior_weightfunction("one-sided", c(.05), wf_fixed(c(1, 1)))
    ),
    "scalar probability prior"
  )
})


test_that(".get_spike_and_slab_variable() requires spike_and_slab prior", {
  expect_error(BayesTools:::.get_spike_and_slab_variable(prior("normal", list(0, 1))),
               "only works with spike_and_slab priors")
})


test_that(".get_spike_and_slab_inclusion() requires spike_and_slab prior", {
  expect_error(BayesTools:::.get_spike_and_slab_inclusion(prior("normal", list(0, 1))),
               "only works with spike_and_slab priors")
})


# ============================================================================ #
# SECTION: prior_mixture() construction
# ============================================================================ #

test_that("prior_mixture() creates factor_mixture with spike component", {
  p1 <- prior("spike", list(0))
  p2 <- prior_factor("mnormal", list(0, 1), contrast = "orthonormal")
  p2$parameters[["K"]] <- 2

  p_mix <- prior_mixture(list(p1, p2), components = c("null", "alt"))

  expect_s3_class(p_mix, "prior.factor_mixture")
})


test_that("prior_mixture() handles prior_none in factor mixtures", {
  p1 <- prior_none()
  p2 <- prior_factor("mnormal", list(0, 1), contrast = "orthonormal")
  p2$parameters[["K"]] <- 2

  p_mix <- prior_mixture(list(p1, p2), components = c("null", "alt"))

  expect_s3_class(p_mix, "prior.factor_mixture")
})


test_that("prior_mixture() creates bias_mixture for PET/PEESE/weightfunction", {
  p_pet <- prior_PET("normal", list(0, 1))
  p_wf  <- prior_weightfunction("one-sided", c(0.05), wf_cumulative(c(1, 1)))

  p_mix <- prior_mixture(list(p_pet, p_wf), components = c("a", "b"))

  expect_s3_class(p_mix, "prior.bias_mixture")
})


# ============================================================================ #
# SECTION: Distribution parameter validation
# ============================================================================ #

test_that("uniform prior requires a < b", {
  expect_error(prior("uniform", list(a = 5, b = 1)),
               "lower than")
})


test_that("uniform prior honors truncation", {
  p <- prior("uniform", list(a = 0, b = 10), truncation = list(2, 5))

  expect_equal(p$truncation, list(lower = 2, upper = 5))
  expect_equal(pdf(p, c(1, 3, 6)), c(0, 1 / 3, 0))
  expect_equal(cdf(p, c(1, 2, 3.5, 5, 6)), c(0, 0, .5, 1, 1))
  samples <- rng(p, 100)
  expect_true(all(samples >= 2 & samples <= 5))

  # JAGS cannot bound dunif: the syntax is the uniform distribution on the
  # intersection of the range and the truncation
  syntax <- JAGS_add_priors("model{}", list(theta = p))
  expect_match(syntax, "theta ~ dunif(2,5)\n", fixed = TRUE)
  expect_false(grepl("T(", syntax, fixed = TRUE))

  posterior <- matrix(c(2.5, 3.5), ncol = 1)
  colnames(posterior) <- "theta"
  prepared <- JAGS_bridgesampling_posterior(posterior, list(theta = p))
  expect_equal(attr(prepared, "lb"), c(theta = 2))
  expect_equal(attr(prepared, "ub"), c(theta = 5))
})


test_that("a truncated uniform prior emits dunif on the intersection without T()", {

  emitted <- function(prior_list) JAGS_add_priors("model{}", prior_list)

  # the truncation narrows the range on one or both sides
  expect_match(emitted(list(x = prior("uniform", list(-1, 2), truncation = list(0, 1)))), "x ~ dunif(0,1)\n", fixed = TRUE)
  expect_match(emitted(list(x = prior("uniform", list(-1, 2), truncation = list(-1, 1)))), "x ~ dunif(-1,1)\n", fixed = TRUE)
  expect_match(emitted(list(x = prior("uniform", list(-1, 2), truncation = list(.5, 2)))), "x ~ dunif(0.5,2)\n", fixed = TRUE)
  expect_match(emitted(list(x = prior("uniform", list(-1, 2), truncation = list(-Inf, 1)))), "x ~ dunif(-1,1)\n", fixed = TRUE)
  expect_match(emitted(list(x = prior("uniform", list(-1, 2), truncation = list(0, Inf)))), "x ~ dunif(0,2)\n", fixed = TRUE)

  # an untruncated uniform prior is unchanged
  expect_identical(emitted(list(x = prior("uniform", list(-1, 2)))), "model{\n x ~ dunif(-1,2)\n\n}")
  expect_identical(emitted(list(x = prior("uniform", list(0, 1)))), "model{\n x ~ dunif(0,1)\n\n}")

  # other distributions keep T()
  expect_match(emitted(list(x = prior("normal", list(0, 1), truncation = list(0, 1)))), "x ~ dnorm(0,1)T(0,1)\n", fixed = TRUE)

  # the distribution of the emitted syntax is the one of the prior
  p <- prior("uniform", list(-1, 2), truncation = list(0, 1))
  x <- c(-.5, 0, .25, .75, 1, 1.5)
  expect_equal(stats::dunif(x, 0, 1), pdf(p, x))

  # also inside spike-and-slab, mixture and factor priors
  slab <- emitted(list(x = prior_spike_and_slab(p)))
  expect_match(slab, "x_variable ~ dunif(0,1)\n", fixed = TRUE)
  expect_false(grepl("T(", slab, fixed = TRUE))
  mixture <- emitted(list(x = prior_mixture(list(p, prior("normal", list(0, 1))))))
  expect_match(mixture, "x_component_1 ~ dunif(0,1)\n", fixed = TRUE)
  expect_false(grepl("T(", mixture, fixed = TRUE))
  factor_syntax <- emitted(list(x = structure(
    prior_factor("uniform", list(-1, 2), truncation = list(0, 1), contrast = "treatment"),
    levels = 2L
  )))
  expect_match(factor_syntax, "x[i] ~ dunif(0,1)\n", fixed = TRUE)
  expect_false(grepl("T(", factor_syntax, fixed = TRUE))
})


test_that("near-support truncation remains an exact truncation", {
  lower <- 1e-12
  p <- prior("exp", list(rate = 1), truncation = list(lower, Inf))

  expect_false(.is_prior_default_range(p))
  expect_identical(pdf(p, lower / 2), 0)
  expect_equal(
    pdf(p, 1),
    stats::dexp(1) / stats::pexp(lower, lower.tail = FALSE)
  )
})


test_that("vector priors reject unsupported truncation", {
  expect_error(
    prior("mnormal", list(mean = 0, sd = 1, K = 2), truncation = list(0, Inf)),
    "Vector priors do not support truncation"
  )
  expect_error(
    prior("mcauchy", list(location = 0, scale = 1, K = 2), truncation = list(-1, 1)),
    "Vector priors do not support truncation"
  )
  expect_error(
    prior("mt", list(location = 0, scale = 1, df = 3, K = 2), truncation = list(-1, Inf)),
    "Vector priors do not support truncation"
  )

  p <- prior("mnormal", list(mean = 0, sd = 1, K = 2))
  p$truncation <- list(lower = 0, upper = Inf)

  expect_error(lpdf(p, c(0, 0)), "Vector priors do not support truncation")
  expect_error(JAGS_add_priors("model{}", list(theta = p)), "Vector priors do not support truncation")

  posterior <- matrix(c(0, 0, 1, 1), ncol = 2, byrow = TRUE)
  colnames(posterior) <- c("theta[1]", "theta[2]")
  expect_error(
    JAGS_bridgesampling_posterior(posterior, list(theta = p)),
    "Vector priors do not support truncation"
  )
})


test_that("exchangeable vector priors expose exact numeric and JAGS APIs", {
  mvt_lpdf <- function(x, location, scale, df) {
    x <- matrix(x, ncol = length(location))
    K <- length(location)
    centered <- sweep(x, 2, location, "-")
    mahalanobis <- rowSums((centered / scale)^2)

    lgamma((df + K) / 2) - lgamma(df / 2) -
      K / 2 * log(df * pi) - K * log(scale) -
      (df + K) / 2 * log1p(mahalanobis / df)
  }

  p_mnormal <- prior("mnormal", list(mean = 1, sd = 2, K = 3))
  p_mcauchy <- prior("mcauchy", list(location = -1, scale = .5, K = 2))
  p_mt <- prior("mt", list(location = 2, scale = .5, df = 5, K = 2))

  expect_s3_class(p_mnormal, "prior.vector")
  expect_s3_class(p_mcauchy, "prior.vector")
  expect_s3_class(p_mt, "prior.vector")
  expect_equal(p_mcauchy$distribution, "mt")
  expect_equal(p_mcauchy$parameters$df, 1)

  x_mnormal <- rbind(c(1, 1, 1), c(3, 1, -1))
  expect_equal(
    lpdf(p_mnormal, x_mnormal),
    rowSums(stats::dnorm(x_mnormal, mean = 1, sd = 2, log = TRUE)),
    tolerance = 1e-12
  )
  expect_equal(pdf(p_mnormal, x_mnormal), exp(lpdf(p_mnormal, x_mnormal)), tolerance = 1e-12)

  x_mcauchy <- rbind(c(-1, -1), c(-.5, 0))
  expect_equal(
    lpdf(p_mcauchy, x_mcauchy),
    mvt_lpdf(x_mcauchy, location = c(-1, -1), scale = .5, df = 1),
    tolerance = 1e-12
  )

  x_mt <- rbind(c(2, 2), c(2.5, 3))
  expect_equal(
    lpdf(p_mt, x_mt),
    mvt_lpdf(x_mt, location = c(2, 2), scale = .5, df = 5),
    tolerance = 1e-12
  )

  q <- rbind(c(1, 2, 3), c(4, 5, 6))
  expect_equal(
    mlpdf(p_mnormal, q),
    matrix(stats::dnorm(as.vector(q), mean = 1, sd = 2, log = TRUE), nrow = nrow(q)),
    tolerance = 1e-12
  )
  expect_equal(
    mcdf(p_mnormal, q),
    matrix(stats::pnorm(as.vector(q), mean = 1, sd = 2), nrow = nrow(q)),
    tolerance = 1e-12
  )

  probs <- c(.25, .5, .75)
  expect_equal(
    mquant(p_mnormal, probs),
    matrix(rep(stats::qnorm(probs, mean = 1, sd = 2), times = 3), ncol = 3),
    tolerance = 1e-12
  )
  expect_equal(mean(p_mnormal), rep(1, 3))
  expect_equal(var(p_mnormal), rep(4, 3))
  expect_equal(sd(p_mnormal), rep(2, 3))
  expect_true(all(is.nan(mean(p_mcauchy))))

  expect_equal(JAGS_to_monitor(list(theta = p_mnormal, beta = p_mt)), c("theta", "beta"))
  expect_equal(
    JAGS_to_monitor(list(
      fixed_scalar = prior("point", list(0)),
      fixed_vector = prior("mpoint", list(1, 2)),
      fixed_factor = prior_factor("point", list(0), contrast = "treatment")
    )),
    c("fixed_scalar", "fixed_vector", "fixed_factor")
  )

  syntax_mnormal <- JAGS_add_priors("model{}", list(theta = p_mnormal))
  expect_match(syntax_mnormal, "prior_par1_theta = rep(1,3)", fixed = TRUE)
  expect_match(syntax_mnormal, "prior_par2_theta[i,i] <- 0.25", fixed = TRUE)
  expect_match(syntax_mnormal, "theta ~ dmnorm(prior_par1_theta,prior_par2_theta)", fixed = TRUE)

  syntax_mt <- JAGS_add_priors("model{}", list(beta = p_mt))
  expect_match(syntax_mt, "prior_par_s_beta ~ dgamma(2.5, 2.5)", fixed = TRUE)
  expect_match(syntax_mt, "prior_par2_beta[i,i] <- 4", fixed = TRUE)
  expect_match(syntax_mt, "beta[i] <- prior_par_z_beta[i]/sqrt(prior_par_s_beta) + 2", fixed = TRUE)

  posterior <- matrix(
    c(1, 1, 1, 3, 1, -1),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(NULL, paste0("theta[", 1:3, "]"))
  )
  prepared <- JAGS_bridgesampling_posterior(posterior, list(theta = p_mnormal))
  expect_equal(as.vector(prepared), as.vector(posterior))
  expect_equal(colnames(prepared), colnames(posterior))
  expect_equal(attr(prepared, "lb"), c("theta[1]" = -Inf, "theta[2]" = -Inf, "theta[3]" = -Inf))
  expect_equal(attr(prepared, "ub"), c("theta[1]" = Inf, "theta[2]" = Inf, "theta[3]" = Inf))

  samples <- posterior[1, ]
  expect_equal(
    JAGS_marglik_priors(samples, list(theta = p_mnormal)),
    sum(stats::dnorm(samples, mean = 1, sd = 2, log = TRUE)),
    tolerance = 1e-12
  )
  expect_equal(JAGS_marglik_parameters(samples, list(theta = p_mnormal))$theta, samples)
})


test_that("D6 numeric Dirichlet declarations enforce the inclusive concentration minimum", {

  message <- "The 'alpha' must be equal or higher than 0.01."
  expect_error(prior("dirichlet", list(alpha = c(.0099, 1))), message, fixed = TRUE)
  expect_error(prior("dirichlet", list(alpha = c(.01, .0099))), message, fixed = TRUE)
  expect_error(prior("simplex", list(concentration = c(.0099, 1))), message, fixed = TRUE)
  expect_identical(prior("dirichlet", list(alpha = c(.01, 1)))$parameters$alpha, c(.01, 1))

  # The established type, positivity, finiteness and length checks run first.
  expect_error(prior("dirichlet", list(alpha = "small")), "The 'alpha' must be a numeric vector.", fixed = TRUE)
  expect_error(prior("dirichlet", list(alpha = c(0, .0099))), "The 'alpha' must be positive.", fixed = TRUE)
  expect_error(prior("dirichlet", list(alpha = c(Inf, .0099))), "The 'alpha' must be finite.", fixed = TRUE)
  expect_error(prior("dirichlet", list(alpha = .0099)), "The Dirichlet 'alpha' concentration vector must contain at least two values.", fixed = TRUE)

  symbolic <- prior("dirichlet", list(alpha = expression(.001, alpha2)))
  expect_identical(symbolic$parameters$alpha, expression(.001, alpha2))
  expect_match(JAGS_add_priors("model{}", list(w = symbolic)),
               "prior_par_eta_w\\[1\\] ~ dgamma\\(0.001, 1\\)")
  expect_error(rng(symbolic, 10), "The 'prior' argument must not contain parameter expressions.", fixed = TRUE)
})

test_that("D6 supported Dirichlet draws retain the extraDistr values and RNG stream", {

  for(alpha in list(c(.01, .01), c(.01, 1), c(1, 1))){
    set.seed(1)
    expected <- extraDistr::rdirichlet(10, alpha = alpha)
    expected_seed <- .Random.seed
    set.seed(1)
    actual <- rng(prior("dirichlet", list(alpha = alpha)), 10)
    expect_identical(unname(actual), expected)
    expect_identical(.Random.seed, expected_seed)
    expect_identical(colnames(actual), c("V1", "V2"))
  }
})

test_that("D6 degenerate Dirichlet sampling has a stable availability condition", {

  set.seed(1)
  error <- tryCatch(rng(prior("dirichlet", list(alpha = c(1e308, 1e308))), 10), error = identity)
  expect_s3_class(error, "BayesTools_dirichlet_rng_unavailable")
  if(inherits(error, "error")){
    expect_identical(class(error), c("BayesTools_dirichlet_rng_unavailable", "BayesTools_prior_rng_unavailable", "error", "condition"))
    expect_identical(conditionMessage(error), "Dirichlet draws are unavailable because sampling produced non-finite or degenerate simplex draws. Use a less extreme 'alpha' specification.")
    expect_null(conditionCall(error))
  }
})

test_that("Dirichlet simplex priors expose joint, marginal, JAGS, and bridge APIs", {
  p <- prior("dirichlet", list(alpha = c(2, 3, 5)))
  p_alias <- prior("simplex", list(concentration = c(2, 3, 5)))

  expect_true(is.prior.simplex(p))
  expect_s3_class(p, "prior.vector")
  expect_equal(p_alias$parameters$alpha, p$parameters$alpha)
  expect_equal(mean(p), c(.2, .3, .5))
  expect_equal(
    var(p),
    c(2, 3, 5) * (10 - c(2, 3, 5)) / (10^2 * 11)
  )
  expect_equal(sd(p), sqrt(var(p)))

  expect_error(prior("dirichlet", list(alpha = 1)), "at least two")
  expect_error(prior("dirichlet", list(alpha = c(1, 0))), "positive")
  expect_error(prior("dirichlet", list(alpha = c(1, Inf))), "finite")
  expect_error(prior("dirichlet", list(alpha = expression(alpha1))), "at least two")

  p_expression <- prior("dirichlet", list(alpha = expression(alpha1, alpha2)))
  expect_true(BayesTools:::.is_prior_expression(p_expression))
  expect_equal(p_expression$parameters$K, 2)
  syntax_expression <- JAGS_add_priors("model{}", list(w = p_expression))
  expect_match(syntax_expression, "prior_par_eta_w\\[1\\] ~ dgamma\\(alpha1, 1\\)")
  expect_match(syntax_expression, "prior_par_eta_w\\[2\\] ~ dgamma\\(alpha2, 1\\)")

  d <- density(p, x_seq = c(.25, .5), truncate_end = FALSE)
  expect_s3_class(d, "density.prior.simplex")
  expect_equal(names(d), c("V1", "V2", "V3"))
  expect_s3_class(d[[1]], "density.prior.simple")
  expect_equal(
    d[[1]]$y,
    stats::dbeta(c(.25, .5), shape1 = 2, shape2 = 8),
    tolerance = 1e-12
  )
  expect_equal(
    d[[2]]$y,
    stats::dbeta(c(.25, .5), shape1 = 3, shape2 = 7),
    tolerance = 1e-12
  )
  expect_s3_class(
    plot(p, plot_type = "ggplot", n_points = 16, show_figures = 1),
    "ggplot"
  )

  x <- c(.2, .3, .5)
  expected_lpdf <- lgamma(10) - sum(lgamma(c(2, 3, 5))) +
    sum((c(2, 3, 5) - 1) * log(x))
  expect_equal(lpdf(p, x), expected_lpdf, tolerance = 1e-12)
  expect_equal(pdf(p, x), exp(expected_lpdf), tolerance = 1e-12)
  roundoff_x <- c(.2, .3, .5 + .Machine$double.eps)
  canonical_x <- roundoff_x / sum(roundoff_x)
  expect_equal(
    lpdf(p, roundoff_x),
    lpdf(p, canonical_x),
    tolerance = 1e-12
  )
  expect_equal(lpdf(p, c(.2, .3, .5 + 1e-10)), -Inf)
  expect_equal(lpdf(p, c(.2, .3, .4)), -Inf)
  expect_equal(lpdf(p, c(.2, -.3, 1.1)), -Inf)

  q <- rbind(c(.2, .3, .5), c(.6, .2, .2))
  expect_equal(
    mlpdf(p, q),
    matrix(
      c(
        stats::dbeta(q[, 1], 2, 8, log = TRUE),
        stats::dbeta(q[, 2], 3, 7, log = TRUE),
        stats::dbeta(q[, 3], 5, 5, log = TRUE)
      ),
      nrow = 2
    ),
    tolerance = 1e-12
  )
  expect_equal(
    mcdf(p, q),
    matrix(
      c(
        stats::pbeta(q[, 1], 2, 8),
        stats::pbeta(q[, 2], 3, 7),
        stats::pbeta(q[, 3], 5, 5)
      ),
      nrow = 2
    ),
    tolerance = 1e-12
  )
  expect_equal(
    mquant(p, c(.25, .75)),
    cbind(
      stats::qbeta(c(.25, .75), 2, 8),
      stats::qbeta(c(.25, .75), 3, 7),
      stats::qbeta(c(.25, .75), 5, 5)
    ),
    tolerance = 1e-12
  )

  set.seed(1)
  draws <- rng(p, 1000)
  expect_equal(dim(draws), c(1000L, 3L))
  expect_true(all(draws >= 0 & draws <= 1))
  expect_equal(rowSums(draws), rep(1, 1000), tolerance = 1e-12)

  syntax <- JAGS_add_priors("model{}", list(w = p))
  expect_match(syntax, "prior_par_eta_w\\[1\\] ~ dgamma\\(2, 1\\)")
  expect_match(syntax, "w\\[3\\] <- prior_par_eta_w\\[3\\] / sum\\(prior_par_eta_w\\[1:3\\]\\)")

  inits <- JAGS_get_inits(list(w = p), chains = 1, seed = 1)[[1]]
  expect_true("prior_par_eta_w" %in% names(inits))
  expect_equal(length(inits$prior_par_eta_w), 3L)
  expect_true(all(inits$prior_par_eta_w > 0))
  expect_equal(JAGS_to_monitor(list(w = p)), c("w", "prior_par_eta_w"))

  expect_error(prior("dirichlet", list(alpha = c(1e-300, 1e-300))),
               "The 'alpha' must be equal or higher than 0.01.", fixed = TRUE)

  posterior <- matrix(
    c(1, 2, 3, 4, 5, 6),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(NULL, paste0("prior_par_eta_w[", 1:3, "]"))
  )
  prepared <- JAGS_bridgesampling_posterior(posterior, list(w = p))
  expect_equal(colnames(prepared), colnames(posterior))
  expect_equal(attr(prepared, "lb"), stats::setNames(rep(0, 3), colnames(posterior)))
  expect_equal(attr(prepared, "ub"), stats::setNames(rep(Inf, 3), colnames(posterior)))

  samples <- posterior[1, ]
  expect_equal(
    JAGS_marglik_priors(samples, list(w = p)),
    sum(stats::dgamma(c(1, 2, 3), shape = c(2, 3, 5), rate = 1, log = TRUE)),
    tolerance = 1e-12
  )
  expect_equal(JAGS_marglik_parameters(samples, list(w = p))$w, c(1, 2, 3) / 6)

  public_posterior <- matrix(
    c(.2, .3, .5, 1, 2, 3),
    nrow = 1,
    dimnames = list(NULL, c("w[1]", "w[2]", "w[3]", paste0("prior_par_eta_w[", 1:3, "]")))
  )
  cleaned <- BayesTools:::.remove_auxiliary_parameters(public_posterior, list(w = p))
  expect_equal(colnames(cleaned$model_samples), paste0("w[", 1:3, "]"))
})


test_that("JAGS prior APIs reject duplicate parameter names", {
  p <- prior("normal", list(0, 1))
  duplicate_prior_list <- list(theta = p, theta = p)

  expect_error(
    JAGS_to_monitor(duplicate_prior_list),
    "must not contain duplicate names"
  )
  expect_error(
    JAGS_add_priors("model{}", duplicate_prior_list),
    "must not contain duplicate names"
  )
  expect_error(
    JAGS_get_inits(duplicate_prior_list, chains = 1, seed = 1),
    "must not contain duplicate names"
  )

  posterior <- matrix(c(0, 1), nrow = 1)
  colnames(posterior) <- c("theta", "theta")
  expect_error(
    JAGS_bridgesampling_posterior(posterior, duplicate_prior_list),
    "must not contain duplicate names"
  )
})

test_that("JAGS prior APIs require complete parameter names", {
  p <- prior("normal", list(0, 1))
  malformed_prior_lists <- list(
    list(p),
    structure(list(p), names = ""),
    structure(list(p), names = NA_character_),
    structure(list(theta = p, p), names = c("theta", ""))
  )

  for(prior_list in malformed_prior_lists){
    expect_error(
      JAGS_to_monitor(prior_list),
      "must be a fully named list",
      fixed = TRUE
    )
    expect_error(
      JAGS_add_priors("model{}", prior_list),
      "must be a fully named list",
      fixed = TRUE
    )
    expect_error(
      JAGS_get_inits(prior_list, chains = 1, seed = 1),
      "must be a fully named list",
      fixed = TRUE
    )
  }
})


test_that("JAGS prior APIs reject malformed factor and vector attributes", {
  p_factor <- prior_factor("normal", list(0, 1), contrast = "treatment")
  attr(p_factor, "levels") <- 0
  expect_error(
    JAGS_add_priors("model{}", list(theta = p_factor)),
    "levels.*equal or higher than 1"
  )

  p_vector <- prior("mnormal", list(0, 1, 2))
  p_vector$parameters[["K"]] <- 0
  expect_error(
    JAGS_add_priors("model{}", list(theta = p_vector)),
    "K.*equal or higher than 1"
  )
  expect_error(
    JAGS_bridgesampling_posterior(
      matrix(numeric(0), nrow = 1, dimnames = list(NULL, character())),
      list(theta = p_vector)
    ),
    "K.*equal or higher than 1"
  )

  p_vector$parameters[["K"]] <- 1.5
  expect_error(
    JAGS_add_priors("model{}", list(theta = p_vector)),
    "K.*integer vector"
  )
})


test_that("multivariate prior APIs reject malformed vector prior internals", {
  p_vector <- prior("mnormal", list(mean = 0, sd = 1, K = 2))
  expect_error(
    mlpdf(p_vector, matrix(c(0, 1, 2), nrow = 1)),
    "x.*2 columns"
  )
  expect_error(
    mcdf(p_vector, matrix(c(0, 1, 2), nrow = 1)),
    "q.*2 columns"
  )

  p_vector$parameters[["mean"]] <- c(0, 1)
  expect_error(
    mlpdf(p_vector, matrix(c(0, 1), nrow = 1)),
    "Only exchangeable vector prior distributions"
  )
  expect_error(
    mcdf(p_vector, matrix(c(0, 1), nrow = 1)),
    "Only exchangeable vector prior distributions"
  )
  expect_error(
    mquant(p_vector, 0.5),
    "Only exchangeable vector prior distributions"
  )
  expect_error(
    mean(p_vector),
    "Only exchangeable vector prior distributions"
  )
})


test_that("density.prior rejects bad transformation lists", {
  p <- prior("normal", list(0, 1))

  expect_error(
    density(p, x_seq = 0, transformation = list(fun = exp, inv = log)),
    "The 'jac' objects are missing in the 'transformation' argument.",
    fixed = TRUE
  )
  expect_error(
    density(p, x_seq = 0, transformation = list(fun = exp, inv = log, bad = exp)),
    "'jac' objects are missing"
  )
  expect_error(
    density(p, x_seq = 0, transformation = 1),
    "Uknown format"
  )
})


# ============================================================================ #
# SECTION: rng() function with sample_components
# ============================================================================ #

test_that("rng() spike_and_slab returns component indicators", {
  p_ss <- prior_spike_and_slab(
    prior_parameter = prior("normal", list(0, 1)),
    prior_inclusion = prior("beta", list(1, 1))
  )

  set.seed(1)
  components <- rng(p_ss, 10, sample_components = TRUE)

  expect_true(all(components %in% c(0, 1)))
  expect_length(components, 10)
})


test_that("rng() mixture returns component indicators", {
  p_mix <- prior_mixture(
    list(prior("normal", list(0, 1)), prior("normal", list(3, 1))),
    components = c("a", "b")
  )

  set.seed(1)
  components <- rng(p_mix, 10, sample_components = TRUE)

  expect_true(all(components %in% 1:2))
  expect_length(components, 10)
})


test_that("rng() factor_mixture with transform_factor_samples", {
  p_mix <- prior_mixture(
    list(
      prior("spike", list(0)),
      prior_factor("mnormal", list(0, 1), contrast = "orthonormal")
    ),
    components = c("null", "alt")
  )
  p_mix <- prior_factor_levels(p_mix, 3)

  set.seed(1)
  samples <- rng(p_mix, 10, transform_factor_samples = FALSE)
  expect_true(is.matrix(samples))

  samples2 <- rng(p_mix, 10, transform_factor_samples = TRUE)
  expect_true(is.matrix(samples2))
  expect_equal(ncol(samples2), 3)  # K+1 columns
})


test_that("rng() orthonormal prior with transform_factor_samples", {
  p <- prior_factor("mnormal", list(0, 1), contrast = "orthonormal")
  p$parameters[["K"]] <- 2

  set.seed(1)
  samples <- rng(p, 10, transform_factor_samples = TRUE)

  expect_true(is.matrix(samples))
  expect_equal(ncol(samples), 3)  # K+1 columns
})

test_that("rng() factor spike-and-slab honours transform_factor_samples", {

  slab <- prior_factor("mnormal", list(0, 1), contrast = "meandif")
  p <- prior_spike_and_slab(slab, prior_inclusion = prior("spike", list(.5)))
  p <- prior_factor_levels(p, 3)

  set.seed(4103)
  raw <- rng(p, 50, transform_factor_samples = FALSE)
  set.seed(4103)
  levels <- rng(p, 50, transform_factor_samples = TRUE)

  expect_equal(dim(raw), c(50L, 2L))
  expect_equal(dim(levels), c(50L, 3L))
  expect_identical(attr(raw, "inclusion"), attr(levels, "inclusion"))
  expect_equal(as.vector(levels),
               as.vector(unclass(raw) %*% t(contr.meandif(1:3))))
  expect_true(all(raw[attr(raw, "inclusion") == 0, ] == 0))
})


# ============================================================================ #
# SECTION: cdf() function edge cases
# ============================================================================ #

test_that("cdf() not implemented for spike_and_slab", {
  p_ss <- prior_spike_and_slab(
    prior_parameter = prior("normal", list(0, 1)),
    prior_inclusion = prior("beta", list(1, 1))
  )

  expect_error(cdf(p_ss, 0), "No cdfs are implemented for spike and slab")
})


test_that("cdf() handles truncated priors correctly", {
  p <- prior("normal", list(0, 1), truncation = list(-2, 2))

  expect_true(cdf(p, 0) > 0)
  expect_true(cdf(p, -3) == 0)  # Below truncation
  expect_true(cdf(p, 3) >= 1 - 1e-6)  # Above truncation
  expect_equal(cdf(p, c(NA_real_, -3, 0, 3)), c(NA_real_, 0, cdf(p, 0), 1))
})


# ============================================================================ #
# SECTION: ccdf() function edge cases
# ============================================================================ #

test_that("ccdf() not implemented for spike_and_slab or mixture", {
  p_ss <- prior_spike_and_slab(
    prior_parameter = prior("normal", list(0, 1)),
    prior_inclusion = prior("beta", list(1, 1))
  )
  expect_error(ccdf(p_ss, 0), "No ccdf are implemented for spike and slab")

  p_mix <- prior_mixture(
    list(prior("normal", list(0, 1)), prior("normal", list(3, 1))),
    components = c("a", "b")
  )
  expect_error(ccdf(p_mix, 0), "No ccdf are implemented for prior mixtures")
})


test_that("ccdf() handles truncated priors correctly", {
  p <- prior("normal", list(0, 1), truncation = list(-2, 2))
  expect_true(ccdf(p, 0) > 0)
  expect_true(ccdf(p, 3) == 0)  # Above truncation
  expect_equal(ccdf(p, c(NA_real_, -3, 0, 3)), c(NA_real_, 1, ccdf(p, 0), 0))
})


# ============================================================================ #
# SECTION: lpdf() function edge cases
# ============================================================================ #

test_that("lpdf() not implemented for spike_and_slab or mixture", {
  p_ss <- prior_spike_and_slab(
    prior_parameter = prior("normal", list(0, 1)),
    prior_inclusion = prior("beta", list(1, 1))
  )
  expect_error(lpdf(p_ss, 0), "No lpdf are implemented for spike and slab")

  p_mix <- prior_mixture(
    list(prior("normal", list(0, 1)), prior("normal", list(3, 1))),
    components = c("a", "b")
  )
  expect_error(lpdf(p_mix, 0), "No lpdf are implemented for prior mixtures")
})


# ============================================================================ #
# SECTION: quant() function edge cases
# ============================================================================ #

test_that("quant() not implemented for spike_and_slab or mixture", {
  p_ss <- prior_spike_and_slab(
    prior_parameter = prior("normal", list(0, 1)),
    prior_inclusion = prior("beta", list(1, 1))
  )
  expect_error(quant(p_ss, 0.5), "No quant(ile)? functions? are implemented for spike and slab")

  p_mix <- prior_mixture(
    list(prior("normal", list(0, 1)), prior("normal", list(3, 1))),
    components = c("a", "b")
  )
  expect_error(quant(p_mix, 0.5), "No quant(ile)? functions? are implemented for prior mixtures")
})


test_that("quant() handles truncated priors with optimization", {
  p <- prior("normal", list(0, 1), truncation = list(0.5, 2))

  q <- quant(p, 0.5)
  expect_true(q > 0.5 && q < 2)

  q_low <- quant(p, 0.01)
  q_high <- quant(p, 0.99)
  expect_true(q_low >= 0.5)
  expect_true(q_high <= 2)
})


test_that("extreme finite truncated normal methods remain tail-stable", {
  p_left  <- prior("normal", list(0, 1), truncation = list(-40, -39))
  p_right <- prior("normal", list(0, 1), truncation = list(39, 40))
  probs   <- c(0, .01, .25, .5, .75, .99, 1)

  q_left  <- quant(p_left, probs)
  q_right <- quant(p_right, probs)

  expect_true(all(is.finite(q_left)))
  expect_true(all(is.finite(q_right)))
  expect_true(all(diff(q_left) > 0))
  expect_true(all(diff(q_right) > 0))
  expect_equal(q_left[c(1, length(q_left))], c(-40, -39))
  expect_equal(q_right[c(1, length(q_right))], c(39, 40))
  expect_equal(q_left, -rev(q_right), tolerance = 1e-12)

  expect_equal(cdf(p_left, q_left), probs, tolerance = 1e-10)
  expect_equal(cdf(p_right, q_right), probs, tolerance = 1e-10)
  expect_equal(ccdf(p_left, q_left), 1 - probs, tolerance = 1e-10)
  expect_equal(ccdf(p_right, q_right), 1 - probs, tolerance = 1e-10)

  x_left    <- c(-40, -39.5, -39)
  x_right   <- c(39, 39.5, 40)
  pdf_left  <- pdf(p_left, x_left)
  pdf_right <- pdf(p_right, x_right)

  expect_true(all(is.finite(pdf_left)))
  expect_true(all(is.finite(pdf_right)))
  expect_true(all(pdf_left > 0))
  expect_true(all(pdf_right > 0))
  expect_equal(pdf_left, rev(pdf_right), tolerance = 1e-12)
  expect_equal(log(pdf_left), lpdf(p_left, x_left), tolerance = 1e-12)
  expect_equal(log(pdf_right), lpdf(p_right, x_right), tolerance = 1e-12)

  expect_equal(
    stats::integrate(function(x) pdf(p_left, x), -40, -39)$value,
    1,
    tolerance = 1e-10
  )
  expect_equal(
    stats::integrate(function(x) pdf(p_right, x), 39, 40)$value,
    1,
    tolerance = 1e-10
  )
})

test_that("far-tail truncated normal and t draws stay finite and inside the truncation", {

  # With n = 2000 draws the median of cdf(draws) ~ U(0, 1) has SE 0.011;
  # 0.05 is about 4.5 SE.
  set.seed(4101)
  for(truncation in list(c(-40, -39), c(-Inf, -39), c(39, 40))){
    p <- prior("normal", list(0, 1), truncation = as.list(truncation))
    x <- rng(p, 2000)
    expect_true(all(is.finite(x)))
    expect_true(all(x >= truncation[1] & x <= truncation[2]))
    expect_lt(abs(stats::median(cdf(p, x)) - .5), .05)
  }
  # Ordinary truncations keep the lower-tail inversion bit for bit.
  set.seed(4102)
  half_normal <- rng(prior("normal", list(0, 1), truncation = list(0, Inf)), 5)
  set.seed(4102)
  expect_identical(half_normal, stats::qnorm(stats::runif(5, .5, 1)))

  # Truncation masses below 1e-12 need the exact upper tail of stats::pt.
  far_t <- prior("t", list(0, 1, 3), truncation = list(1e5, 2e5))
  upper <- stats::pt(c(1e5, 2e5), 3, lower.tail = FALSE)
  probs <- c(.1, .5, .9)
  reference <- stats::qt(upper[1] - probs * (upper[1] - upper[2]), 3, lower.tail = FALSE)
  expect_equal(quant(far_t, probs), reference, tolerance = 1e-10)
  expect_equal(cdf(far_t, reference), probs, tolerance = 1e-8)
  expect_equal(stats::integrate(function(x) pdf(far_t, x), 1e5, 2e5, rel.tol = 1e-10)$value,
               1, tolerance = 1e-8)
  x <- rng(far_t, 2000)
  expect_true(all(x >= 1e5 & x <= 2e5))

  # Cauchy: S(x) = atan(1 / x) / pi, so the conditional median beyond x0 is
  # 1 / tan(pi * S(x0) / 2).
  far_cauchy <- prior("cauchy", list(0, 1), truncation = list(1e13, Inf))
  expect_equal(quant(far_cauchy, .5), 1 / tan(pi * atan(1 / 1e13) / pi / 2), tolerance = 1e-10)
  expect_equal(ccdf(far_cauchy, 2e13), atan(1 / 2e13) / atan(1 / 1e13), tolerance = 1e-10)
  expect_true(all(rng(far_cauchy, 2000) >= 1e13))
})


# ============================================================================ #
# SECTION: Multivariate distribution functions (mcdf, mccdf, mlpdf, mquant)
# ============================================================================ #

test_that("mcdf() works for orthonormal and meandif priors", {
  p_orth <- prior_factor("mnormal", list(0, 1), contrast = "orthonormal")
  p_orth$parameters[["K"]] <- 2

  cdf_val <- mcdf(p_orth, 0)
  expect_true(cdf_val >= 0 && cdf_val <= 1)

  p_md <- prior_factor("mnormal", list(0, 1), contrast = "meandif")
  p_md$parameters[["K"]] <- 2

  cdf_val2 <- mcdf(p_md, 0)
  expect_true(cdf_val2 >= 0 && cdf_val2 <= 1)
})


test_that("mccdf() works for orthonormal and meandif priors", {
  p_orth <- prior_factor("mnormal", list(0, 1), contrast = "orthonormal")
  p_orth$parameters[["K"]] <- 2

  ccdf_val <- mccdf(p_orth, 0)
  expect_true(ccdf_val >= 0 && ccdf_val <= 1)

  p_md <- prior_factor("mnormal", list(0, 1), contrast = "meandif")
  p_md$parameters[["K"]] <- 2

  ccdf_val2 <- mccdf(p_md, 0)
  expect_true(ccdf_val2 >= 0 && ccdf_val2 <= 1)
})


test_that("mlpdf() works for orthonormal and meandif priors", {
  p_orth <- prior_factor("mnormal", list(0, 1), contrast = "orthonormal")
  p_orth$parameters[["K"]] <- 2

  lpdf_val <- mlpdf(p_orth, 0)
  expect_true(is.finite(lpdf_val))

  p_md <- prior_factor("mnormal", list(0, 1), contrast = "meandif")
  p_md$parameters[["K"]] <- 2

  lpdf_val2 <- mlpdf(p_md, 0)
  expect_true(is.finite(lpdf_val2))
})


test_that("mquant() works for orthonormal and meandif priors", {
  p_orth <- prior_factor("mnormal", list(0, 1), contrast = "orthonormal")
  p_orth$parameters[["K"]] <- 2

  q <- mquant(p_orth, 0.5)
  expect_true(is.numeric(q))

  p_md <- prior_factor("mnormal", list(0, 1), contrast = "meandif")
  p_md$parameters[["K"]] <- 2

  q2 <- mquant(p_md, 0.5)
  expect_true(is.numeric(q2))
})


# ============================================================================ #
# SECTION: S3 dispatch and generic functions
# ============================================================================ #

test_that("pdf() S3 dispatch works for prior objects", {
  p <- prior("normal", list(0, 1))
  expect_true(is.numeric(pdf(p, 0)))
})


# ============================================================================ #
# SECTION: mean() function edge cases
# ============================================================================ #

test_that("mean() works for spike_and_slab priors", {
  p_ss <- prior_spike_and_slab(
    prior_parameter = prior("normal", list(1, 1)),
    prior_inclusion = prior("point", list(0.5))
  )

  m <- mean(p_ss)
  expect_true(is.numeric(m))
  expect_equal(m, 0.5)  # mean(normal(1,1)) * 0.5
})


test_that("mean() handles truncated distributions and undefined moments", {
  # Truncated normal
  p <- prior("normal", list(0, 1), truncation = list(0, Inf))
  m <- mean(p)
  expect_true(m > 0)

  # Bounded truncation makes otherwise undefined heavy-tailed means finite
  p_t <- prior("t", list(0, 1, 1), truncation = list(-1, 1))
  m_t <- mean(p_t)
  expect_false(is.nan(m_t))

  p_ig <- prior("invgamma", list(.5, 1), truncation = list(.1, 10))
  expect_false(is.nan(mean(p_ig)))

  # Unbounded problematic tails still have undefined means
  expect_true(is.nan(mean(prior("t", list(0, 1, 1), truncation = list(0, Inf)))))
  expect_true(is.nan(mean(prior("invgamma", list(.5, 1)))))
})

test_that("truncated normal moments are analytic for narrow, offset, and far-tail priors", {

  # References: closed forms for the half-normal and the effectively
  # untruncated / shifted cases; the others from the closed form evaluated in
  # 80-digit mpmath arithmetic (erfc normaliser), confirmed by split mpmath
  # quadrature. Integration previously returned 0 for the narrow and offset
  # priors.
  cases <- list(
    list(prior = prior("normal", list(0, 1), list(0, Inf)),
         mean = sqrt(2 / pi), var = 1 - 2 / pi, tol_var = 1e-12),
    list(prior = prior("normal", list(0, 1), list(-1, 2)),
         mean = 0.2296371790913289686, var = 0.5197625392115339359, tol_var = 1e-12),
    list(prior = prior("normal", list(0, .02), list(-10, 10)),
         mean = 0, var = 4e-4, tol_var = 1e-12),
    list(prior = prior("normal", list(1e4, 1), list(0, Inf)),
         mean = 1e4, var = 1, tol_var = 1e-12),
    list(prior = prior("normal", list(0, 1), list(8, 9)),
         mean = 8.121188992979797123, var = 0.01414854278274811104, tol_var = 1e-10),
    # 39 SD into the tail the variance cancels terms of order 39^2 whose
    # log-space ratios carry about 1e-13 relative error: the closed form is
    # still accurate to 1e-6 (actual error about 1e-7), but the conservative
    # rounding-error bound of var() exceeds 1e-6 there, so var() and sd()
    # stop; the closed form itself is checked below.
    list(prior = prior("normal", list(0, 1), list(39, 40)),
         mean = 39.02560741993010845, var = 6.548827702932774827e-4, tol_var = 1e-6,
         bounded = FALSE),
    list(prior = prior("normal", list(0, 1), list(-Inf, -39)),
         mean = -39.02560741993010846, var = 6.548827702932843032e-4, tol_var = 1e-6,
         bounded = FALSE)
  )
  for(case in cases){
    if(case$mean == 0){
      expect_lt(abs(mean(case$prior)), 1e-15)
    }else{
      expect_equal(mean(case$prior), case$mean, tolerance = 1e-12)
    }
    if(isFALSE(case$bounded)){
      expect_equal(
        BayesTools:::.prior_normal_truncated_variance(case$prior)$value,
        case$var,
        tolerance = case$tol_var
      )
      expect_error(var(case$prior), "is unavailable in double-precision arithmetic", fixed = TRUE)
      expect_error(sd(case$prior), "is unavailable in double-precision arithmetic", fixed = TRUE)
      next
    }
    expect_equal(var(case$prior), case$var, tolerance = case$tol_var)
    expect_equal(sd(case$prior), sqrt(case$var), tolerance = case$tol_var)
  }

  expect_error(
    var(prior("normal", list(0, 1), list(5000, Inf))),
    "The variance of the truncated normal prior is unavailable in double-precision arithmetic",
    fixed = TRUE
  )
})

test_that("truncated normal means stop when their rounding error may exceed 1e-6 of the interval width", {

  # References: the closed form evaluated in 80-digit mpmath arithmetic from the
  # exact double bounds (erfc normaliser). A narrow interval cancels two ratios
  # of about 1 / width, so the double-precision mean left its support: N(0, 100)
  # on [1, 1 + 1e-6] gave 1.0000184 and N(0, 1) on [-3, -3 + 1e-12] gave -3.0023.
  accepted <- list(
    list(mu = 0, sigma = 1,   lower = 1,  upper = 1.01, mean = 1.004991625042014735311553),
    list(mu = 0, sigma = 1,   lower = 1,  upper = 1.001, mean = 1.000499916625004115091039),
    list(mu = 0, sigma = 1,   lower = 5,  upper = 5.1,  mean = 5.045810831656684276787259),
    list(mu = 0, sigma = 1,   lower = 0,  upper = 0.1,  mean = 0.04995834723791308296934359),
    list(mu = 0, sigma = 1,   lower = -3, upper = -2.9, mean = -2.947546040865694359270231),
    list(mu = 2, sigma = 0.5, lower = 2.5, upper = 2.6, mean = 2.548170581517787731710686),
    list(mu = 0, sigma = 100, lower = 1,  upper = 1.5,  mean = 1.249997395835503641118044),
    list(mu = 0, sigma = 1,   lower = 39, upper = 40,   mean = 39.02560741993010844958296)
  )
  for(case in accepted){
    p <- prior("normal", list(case$mu, case$sigma), truncation = list(case$lower, case$upper))
    m <- mean(p)
    expect_gte(m, case$lower)
    expect_lte(m, case$upper)
    expect_lte(abs(m - case$mean), 1e-6 * (case$upper - case$lower))
    expect_equal(m, case$mean, tolerance = 1e-9)
  }

  # intervals with an infinite bound have a single ratio and are not gated
  expect_equal(mean(prior("normal", list(0, 1), list(8, Inf))), 8.12136811223611268065352, tolerance = 1e-12)
  expect_equal(mean(prior("normal", list(0, 1), list(-Inf, -39))), -39.02560741993010846, tolerance = 1e-12)

  # narrow intervals refuse like the variance of the same prior
  refused <- list(
    prior("normal", list(0, 100), truncation = list(1, 1 + 1e-6)),
    prior("normal", list(0, 1),   truncation = list(1, 1 + 1e-8)),
    prior("normal", list(0, 1),   truncation = list(-3, -3 + 1e-12))
  )
  for(p in refused){
    expect_error(
      mean(p),
      "The mean of the truncated normal prior is unavailable in double-precision arithmetic",
      fixed = TRUE)
    expect_error(
      var(p),
      "The variance of the truncated normal prior is unavailable in double-precision arithmetic",
      fixed = TRUE)
  }
  expect_error(
    mean(prior("normal", list(0, 1), truncation = list(1, 1 + 1e-8))),
    paste0(
      "The mean of the truncated normal prior is unavailable in double-precision ",
      "arithmetic: its estimated rounding error is [0-9.e+]+ times the width of the ",
      "truncation interval for an interval of width 1e-08 prior standard deviations ",
      "starting 1 prior standard deviations from the mean[.] Estimate the mean from ",
      "rng[(][)] draws of the prior, which sample the truncated distribution exactly, ",
      "or use a wider truncation interval[.]"
    )
  )
  expect_error(
    mean(prior_spike_and_slab(prior("normal", list(0, 1), truncation = list(1, 1 + 1e-8)))),
    "The mean of the truncated normal prior is unavailable in double-precision arithmetic",
    fixed = TRUE)

  # a returned mean never leaves the support
  for(lower in c(-3, 0, 1, 5, 20)){
    for(width in 10^(-12:1)){
      p <- prior("normal", list(0, 1), truncation = list(lower, lower + width))
      m <- tryCatch(mean(p), error = function(e) NA_real_)
      if(!is.na(m)){
        expect_gte(m, p$truncation[["lower"]])
        expect_lte(m, p$truncation[["upper"]])
      }
    }
  }
})

test_that("truncated normal variances stop when their rounding error may exceed 1e-6", {

  # References: the closed form evaluated in 50-digit mpmath arithmetic from
  # the exact double bounds (agreeing with an 80-digit evaluation to < 1e-21).
  # The returned values must lie within the estimated rounding-error bound,
  # which is at most 1e-6 whenever a value is returned. Validation over 4610
  # far-tail and narrow truncations from 1 to 5000 SD found no accepted value
  # with a larger error and no error above the bound.
  accepted <- list(
    list(prior = prior("normal", list(0, 1), list(10, Inf)),
         var = 0.00944537782565626116413681765037),
    list(prior = prior("normal", list(0, 1), list(-Inf, -20)),
         var = 0.00246326161505216359968528619983),
    list(prior = prior("normal", list(0, 1), list(20, 21)),
         var = 0.0024632604300108293913512820229),
    list(prior = prior("normal", list(0, 1), list(30, Inf)),
         var = 0.0011037715118900910011367413855),
    list(prior = prior("normal", list(0.3, 0.7), list(0.3 + 0.7 * 25, Inf)),
         var = 0.000776572320492515523112347483473),
    list(prior = prior("normal", list(0, 1), list(5, 5.1)),
         var = 0.000822546087919772895415146023238),
    list(prior = prior("normal", list(0, 1), list(1, 1.01)),
         var = 0.00000833326347178311381650572698637)
  )
  for(case in accepted){
    bound <- BayesTools:::.prior_normal_truncated_variance(case$prior)$relative_error
    expect_lte(bound, 1e-6)
    expect_lte(abs(var(case$prior) - case$var) / case$var, bound)
    expect_lte(abs(sd(case$prior) - sqrt(case$var)) / sqrt(case$var), bound)
  }

  # 1000 SD into a tail the double-precision value is 49 times too large
  # (4.90e-5 instead of 1.0e-6).
  expect_error(
    var(prior("normal", list(0, 1), list(1000, Inf))),
    paste0(
      "The variance of the truncated normal prior is unavailable in double-precision ",
      "arithmetic: its estimated rounding error exceeds the computed variance for a ",
      "truncation from 1000 to Inf prior standard deviations from the mean. Estimate ",
      "the variance from rng() draws of the prior, which sample the truncated ",
      "distribution exactly, or use a truncation closer to the prior mean or a wider ",
      "truncation interval."
    ),
    fixed = TRUE
  )
  expect_error(
    sd(prior("normal", list(2, 0.5), list(-Inf, 2 - 0.5 * 1000))),
    "exceeds the computed variance for a truncation from -Inf to -1000 prior standard deviations",
    fixed = TRUE
  )
  expect_error(
    var(prior("normal", list(0, 1), list(50, Inf))),
    paste0(
      "The variance of the truncated normal prior is unavailable in double-precision ",
      "arithmetic: its estimated relative rounding error is 1.2e-05 for a ",
      "truncation from 50 to Inf prior standard deviations from the mean. Estimate ",
      "the variance from rng() draws of the prior, which sample the truncated ",
      "distribution exactly, or use a truncation closer to the prior mean or a wider ",
      "truncation interval."
    ),
    fixed = TRUE
  )
  # The suggested alternative: the variance of rng() draws matches the mpmath
  # reference 9.99994000049999482e-07 at 1000 SD within Monte Carlo error
  # (far-tail draws are exponential-like, kurtosis 9: relative SE sqrt(8 / n),
  # 0.9% for n = 1e5; tolerance 5 SE).
  set.seed(1)
  far_tail_draws <- rng(prior("normal", list(0, 1), list(1000, Inf)), 1e5)
  expect_true(all(far_tail_draws >= 1000))
  expect_equal(stats::var(far_tail_draws), 9.99994000049999482e-07, tolerance = 5 * sqrt(8 / 1e5))
  # A 1e-4 SD wide interval: the value is 1e-3 relative off (mpmath reference
  # 8.333333330387e-10).
  expect_error(
    var(prior("normal", list(0, 1), list(0.2, 0.2001))),
    "its estimated relative rounding error is 0.075 for a truncation from 0.2 to 0.2001 prior",
    fixed = TRUE
  )
  # Spike-and-slab variances use the same slab variance.
  expect_error(
    var(prior_spike_and_slab(prior("normal", list(0, 1), list(1000, Inf)))),
    "is unavailable in double-precision arithmetic",
    fixed = TRUE
  )
})

test_that("native inverse-gamma helpers match gamma-transform identities", {
  p <- prior("invgamma", list(shape = 3, scale = 2))
  x <- c(.25, 1, 2)
  expected_lpdf <- 3 * log(2) - lgamma(3) - 4 * log(x) - 2 / x

  expect_equal(lpdf(p, x), expected_lpdf, tolerance = 1e-12)
  expect_equal(pdf(p, x), exp(expected_lpdf), tolerance = 1e-12)
  expect_equal(
    cdf(p, x),
    stats::pgamma(1 / x, shape = 3, rate = 2, lower.tail = FALSE),
    tolerance = 1e-12
  )
  expect_equal(
    ccdf(p, x),
    stats::pgamma(1 / x, shape = 3, rate = 2, lower.tail = TRUE),
    tolerance = 1e-12
  )

  probs <- c(.05, .5, .95)
  expected_q <- 1 / stats::qgamma(probs, shape = 3, rate = 2, lower.tail = FALSE)
  expect_equal(quant(p, probs), expected_q, tolerance = 1e-12)
  expect_equal(cdf(p, quant(p, probs)), probs, tolerance = 1e-12)
  expect_equal(
    BayesTools:::.pinvgamma_prior(x, 3, 2, log.p = TRUE),
    log(cdf(p, x)),
    tolerance = 1e-12
  )
  expect_equal(
    BayesTools:::.qinvgamma_prior(log(probs), 3, 2, log.p = TRUE),
    expected_q,
    tolerance = 1e-12
  )

  set.seed(11)
  draws <- rng(p, 100)
  expect_length(draws, 100)
  expect_true(all(is.finite(draws)))
  expect_true(all(draws > 0))
})


test_that("mean() returns NaN for multivariate t with df <= 1", {
  p_mt <- prior_factor("mt", list(0, 1, 1), contrast = "orthonormal")
  p_mt$parameters[["K"]] <- 2

  m <- mean(p_mt)
  expect_true(is.nan(m))
})


# ============================================================================ #
# SECTION: var() function edge cases
# ============================================================================ #

test_that("var() S3 dispatch works for numeric vectors", {
  x <- c(1, 2, 3, 4, 5)
  expect_equal(var(x), stats::var(x))
})


test_that("var() works for spike_and_slab priors", {
  # spike_and_slab with beta inclusion
  p_ss <- prior_spike_and_slab(
    prior_parameter = prior("normal", list(0, 1)),
    prior_inclusion = prior("beta", list(1, 1))
  )
  v <- var(p_ss)
  expect_equal(v, 0.5)

  # spike_and_slab with point inclusion
  p_ss2 <- prior_spike_and_slab(
    prior_parameter = prior("normal", list(0, 1)),
    prior_inclusion = prior("point", list(0.5))
  )
  v2 <- var(p_ss2)
  expect_equal(v2, 0.5)
})


test_that("spike-and-slab moments include nonzero slab mean", {
  p_ss <- prior_spike_and_slab(
    prior_parameter = prior("normal", list(mean = 1, sd = 2)),
    prior_inclusion = prior("point", list(location = 0.25))
  )

  expect_equal(mean(p_ss), 0.25)
  expect_equal(var(p_ss), 0.25 * (2^2 + 1^2) - 0.25^2)
})


test_that("var() returns NaN for distributions with undefined variance", {
  expect_equal(
    var(prior("invgamma", list(4, 2))),
    2^2 / ((4 - 1)^2 * (4 - 2))
  )

  # Bounded truncation makes otherwise undefined heavy-tailed variances finite
  expect_false(is.nan(var(prior("t", list(0, 1, 2), truncation = list(-1, 1)))))
  expect_false(is.nan(var(prior("invgamma", list(2, 1), truncation = list(0.1, 10)))))

  # Unbounded problematic tails still have undefined variances
  expect_true(is.nan(var(prior("t", list(0, 1, 2), truncation = list(0, Inf)))))
  expect_true(is.nan(var(prior("invgamma", list(2, 1)))))
})


test_that("Bernoulli truncation conditions on included atoms", {
  p <- prior("bernoulli", list(0.25), truncation = list(0.5, 1))

  expect_equal(pdf(p, 1), 1)
  expect_equal(pdf(p, 0), 0)
  expect_equal(cdf(p, 0.75), 0)
  expect_equal(cdf(p, 1), 1)
  expect_equal(mean(p), 1)
  expect_equal(var(p), 0)
  expect_equal(quant(p, c(0, .5, 1)), c(1, 1, 1))
  expect_equal(unique(rng(p, 20)), 1)
  expect_error(
    prior("bernoulli", list(0.25), truncation = list(0.25, 0.75)),
    "support point"
  )
  expect_error(
    prior("bernoulli", list(0), truncation = list(0.5, 1)),
    "positive probability mass"
  )
})


test_that("var() works for orthonormal and meandif priors", {
  # orthonormal with mpoint returns 0
  p_mp <- prior_factor("mpoint", list(0), contrast = "orthonormal")
  p_mp$parameters[["K"]] <- 2
  v <- var(p_mp)
  expect_equal(v, 0)

  # orthonormal with mt and df <= 2 returns NaN
  p_mt <- prior_factor("mt", list(0, 1, 2), contrast = "orthonormal")
  p_mt$parameters[["K"]] <- 2
  v_mt <- var(p_mt)
  expect_true(is.nan(v_mt))

  # meandif with mnormal returns positive variance
  p_md <- prior_factor("mnormal", list(0, 1), contrast = "meandif")
  p_md$parameters[["K"]] <- 2
  v_md <- var(p_md)
  expect_true(v_md > 0)
})


test_that("marginal methods and var() of bound factor priors use the declared factor design", {
  design_prior <- function(contrast, formula, data, prior_factor_args = NULL, term){
    p <- do.call(prior_factor, c(list(contrast = contrast), prior_factor_args))
    terms <- attr(stats::terms(formula), "term.labels")
    priors <- stats::setNames(rep(list(p), length(terms)), terms)
    compiled <- JAGS_formula(formula, "mu", data, c(list(intercept = prior("normal", list(0, 1))), priors))
    compiled$prior_list[[term]]
  }
  data_2x3 <- expand.grid(a = factor(c("a1", "a2")), b = factor(c("b1", "b2", "b3")))
  data_3x3 <- expand.grid(g = factor(c("a", "b", "c")), h = factor(c("u", "v", "w")))

  for(contrast in c("orthonormal", "meandif")){
    # the interaction of a 2 x 3 design: every level has the marginal scale of its design row
    bound  <- design_prior(contrast, ~ a * b, data_2x3,
                           list(distribution = "mnormal", parameters = list(mean = 0, sd = 1.5)), "mu_a__xXx__b")
    design <- attr(bound, "factor_design")
    expect_identical(dim(design), c(6L, 2L))
    expect_equal(max(abs(rowSums(design^2) - rowSums(design^2)[1L])), 0, tolerance = 1e-12)
    scale <- 1.5 * sqrt(sum(design[1L, ]^2))
    expect_equal(var(bound),  scale^2)
    expect_equal(sd(bound),   scale)
    expect_equal(mcdf(bound,  c(-.5, .5)), stats::pnorm(c(-.5, .5), 0, scale))
    expect_equal(mccdf(bound, c(-.5, .5)), stats::pnorm(c(-.5, .5), 0, scale, lower.tail = FALSE))
    expect_equal(mlpdf(bound, c(-.5, .5)), stats::dnorm(c(-.5, .5), 0, scale, log = TRUE))
    expect_equal(mpdf(bound,  c(-.5, .5)), stats::dnorm(c(-.5, .5), 0, scale))
    expect_equal(mquant(bound, c(.025, .975)), stats::qnorm(c(.025, .975), 0, scale))

    # multivariate t levels keep their degrees of freedom
    bound_t <- design_prior(contrast, ~ a * b, data_2x3,
                            list(distribution = "mt", parameters = list(location = 0, scale = 1.5, df = 5)), "mu_a__xXx__b")
    design_t <- attr(bound_t, "factor_design")
    scale_t  <- 1.5 * sqrt(sum(design_t[1L, ]^2))
    expect_equal(var(bound_t), scale_t^2 * 5 / 3)
    expect_equal(mcdf(bound_t, .5), extraDistr::plst(.5, df = 5, mu = 0, sigma = scale_t))
    expect_equal(mquant(bound_t, .975), extraDistr::qlst(.975, df = 5, mu = 0, sigma = scale_t))
    expect_equal(mlpdf(bound_t, .5), extraDistr::dlst(.5, df = 5, mu = 0, sigma = scale_t, log = TRUE))
  }

  # the orthonormal 2 x 3 interaction: 1/3 (the single-factor contrast of the coefficient count gave 2/3)
  orthonormal <- design_prior("orthonormal", ~ a * b, data_2x3,
                              list(distribution = "mnormal", parameters = list(mean = 0, sd = 1)), "mu_a__xXx__b")
  expect_equal(var(orthonormal), 1 / 3)
  expect_equal(mcdf(orthonormal, .5), stats::pnorm(.5, 0, sqrt(1 / 3)))

  # the 3 x 3 interaction: 4/9
  orthonormal <- design_prior("orthonormal", ~ g * h, data_3x3,
                              list(distribution = "mnormal", parameters = list(mean = 0, sd = 1)), "mu_g__xXx__h")
  expect_equal(var(orthonormal), 4 / 9)
  expect_equal(mquant(orthonormal, .975), stats::qnorm(.975, 0, 2 / 3))

  # single-factor priors, bound or sized by K, keep the single-factor contrast
  for(contrast in c("orthonormal", "meandif")){
    contrast_row <- if(contrast == "orthonormal") contr.orthonormal(1:4)[1L, ] else contr.meandif(1:4)[1L, ]
    scale <- 2 * sqrt(sum(contrast_row^2))
    sized <- prior_factor("mnormal", list(0, 2), contrast = contrast)
    sized$parameters[["K"]] <- 3
    bound <- prior_factor_levels(prior_factor("mnormal", list(0, 2), contrast = contrast), c("a", "b", "c", "d"))
    for(single in list(sized, bound)){
      expect_equal(var(single), scale^2)
      expect_equal(mcdf(single, .5), stats::pnorm(.5, 0, scale))
      expect_equal(mccdf(single, .5), stats::pnorm(.5, 0, scale, lower.tail = FALSE))
      expect_equal(mlpdf(single, .5), stats::dnorm(.5, 0, scale, log = TRUE))
      expect_equal(mquant(single, .975), stats::qnorm(.975, 0, scale))
    }
  }

  # levels with different marginal scales have no common marginal
  uneven <- prior_factor_levels(prior_factor("mnormal", list(0, 1), contrast = "orthonormal"), c("a", "b", "c"))
  design <- attr(uneven, "factor_design")
  design[1L, ] <- 2 * design[1L, ]
  attr(uneven, "factor_design") <- design
  message <- "The marginal distribution of the factor prior is unavailable: the levels of the declared factor design have different marginal scales. Use 'rng()' to obtain the level samples."
  expect_error(var(uneven),         message, fixed = TRUE)
  expect_error(sd(uneven),          message, fixed = TRUE)
  expect_error(mcdf(uneven, .5),    message, fixed = TRUE)
  expect_error(mccdf(uneven, .5),   message, fixed = TRUE)
  expect_error(mlpdf(uneven, .5),   message, fixed = TRUE)
  expect_error(mpdf(uneven, .5),    message, fixed = TRUE)
  expect_error(mquant(uneven, .5),  message, fixed = TRUE)
})


test_that("var() not implemented for mixture priors", {
  p_mix <- prior_mixture(
    list(prior("normal", list(0, 1)), prior("normal", list(3, 1))),
    components = c("a", "b")
  )
  expect_error(var(p_mix), "No var is implemented for prior mixtures")
})


# ============================================================================ #
# SECTION: weightfunction log-omega marginal semantics
# ============================================================================ #

test_that("weightfunction marginal functions accept the matrix that mquant() returns", {

  # mquant() returns a 1 x J matrix; the other marginal functions took only
  # vectors ("The 'q' argument must be a numeric vector")
  wfs <- list(
    cumulative  = prior_weightfunction("one-sided", c(.05), wf_cumulative(c(5, .5))),
    fixed       = prior_weightfunction("one-sided", c(.05), wf_fixed(c(1, .5))),
    independent = prior_weightfunction("one-sided", c(.05, .6), wf_independent(prior("beta", list(1, 1)))),
    log_omega   = prior_weightfunction("one-sided", c(.05), wf_independent(prior("normal", list(0, 1)), scale = "log_omega"))
  )
  for(wf in wfs){
    q <- mquant(wf, .5)
    expect_identical(dim(q), c(1L, length(wf$steps) + 1L))
    expect_identical(mcdf(wf, q), mcdf(wf, as.vector(q)))
    expect_identical(mccdf(wf, q), mccdf(wf, as.vector(q)))
    expect_identical(mlpdf(wf, q), mlpdf(wf, as.vector(q)))
    expect_identical(mpdf(wf, q), mpdf(wf, as.vector(q)))
    expect_identical(dim(mcdf(wf, q)), c(length(q), length(q)))
  }
})

test_that("log_omega weightfunction priors expose exact transformed marginals", {
  log_prior <- prior("normal", list(mean = 0.25, sd = 0.5))
  p_wf <- prior_weightfunction(
    side    = "one-sided",
    steps   = c(.025, .05),
    weights = wf_independent(log_prior, scale = "log_omega")
  )

  q <- c(0, .5, 1, 2)
  expected_cdf <- cbind(
    ppoint(q, location = 1),
    stats::pnorm(log(q), mean = 0.25, sd = 0.5),
    stats::pnorm(log(q), mean = 0.25, sd = 0.5)
  )
  expected_cdf[q <= 0, 2:3] <- 0

  expect_equal(mcdf(p_wf, q), expected_cdf, tolerance = 1e-12)
  expect_equal(mccdf(p_wf, q), 1 - expected_cdf, tolerance = 1e-12)
  expect_true(all(is.na(mcdf(p_wf, NA_real_))))
  expect_true(all(is.na(mccdf(p_wf, NA_real_))))

  x <- c(0, .5, 1, 2)
  expected_lpdf <- cbind(
    dpoint(x, location = 1, log = TRUE),
    stats::dnorm(log(x), mean = 0.25, sd = 0.5, log = TRUE) - log(x),
    stats::dnorm(log(x), mean = 0.25, sd = 0.5, log = TRUE) - log(x)
  )
  expected_lpdf[x <= 0, 2:3] <- -Inf

  expect_equal(mlpdf(p_wf, x), expected_lpdf, tolerance = 1e-12)
  expect_equal(mpdf(p_wf, x), exp(expected_lpdf), tolerance = 1e-12)

  probs <- c(.1, .5, .9)
  expect_equal(
    mquant(p_wf, probs),
    cbind(
      qpoint(probs, location = 1),
      exp(stats::qnorm(probs, mean = 0.25, sd = 0.5)),
      exp(stats::qnorm(probs, mean = 0.25, sd = 0.5))
    ),
    tolerance = 1e-12
  )
  expect_equal(mquant(p_wf, 0)[1, 1], 1)

  lognormal_mean <- exp(0.25 + 0.5^2 / 2)
  lognormal_var  <- (exp(0.5^2) - 1) * exp(2 * 0.25 + 0.5^2)

  expect_equal(mean(p_wf), c(1, lognormal_mean, lognormal_mean), tolerance = 1e-8)
  expect_equal(var(p_wf), c(0, lognormal_var, lognormal_var), tolerance = 1e-8)
})

test_that("log_omega weight moments follow the log-weight prior's exponential moments", {
  moments <- function(p){
    p_wf <- prior_weightfunction("one-sided", .05, wf_independent(p, scale = "log_omega"))
    c(mean = mean(p_wf)[2], var = var(p_wf)[2])
  }
  lognormal_moments <- function(mu, s) c(exp(mu + s^2 / 2), exp(2 * mu + s^2) * expm1(s^2))
  # exact lognormal identities, also for concentrated and wide priors
  for(case in list(c(2, .01), c(log(1.5), .001), c(0, 10), c(.25, .5))){
    expect_equal(unname(moments(prior("normal", list(case[1], case[2])))), lognormal_moments(case[1], case[2]), tolerance = 1e-8)
  }

  reference <- function(density, lower, upper, anchor = 0){
    m <- vapply(1:2, function(k){
      stats::integrate(function(x) exp(k * (x - anchor)) * density(x), lower, upper, rel.tol = 1e-13)$value
    }, numeric(1))
    c(mean = exp(anchor) * m[1], var = exp(2 * anchor) * (m[2] - m[1]^2))
  }
  # truncated normal, bounded priors and a truncated heavy tail: integration
  trunc <- function(x) stats::dnorm(x) / (stats::pnorm(2) - stats::pnorm(-1))
  expect_equal(unname(moments(prior("normal", list(0, 1), list(-1, 2)))), unname(reference(trunc, -1, 2)), tolerance = 1e-8)
  expect_equal(unname(moments(prior("beta", list(2, 3)))), unname(reference(function(x) stats::dbeta(x, 2, 3), 0, 1)), tolerance = 1e-8)
  expect_equal(unname(moments(prior("uniform", list(-1, 1)))), c(sinh(1), (exp(2) - exp(-2)) / 4 - sinh(1)^2), tolerance = 1e-8)
  expect_equal(unname(moments(prior("uniform", list(-30, -29)))), c(exp(-29) * (1 - exp(-1)), exp(-58) * ((1 - exp(-2)) / 2 - (1 - exp(-1))^2)), tolerance = 1e-8)
  t_trunc <- prior("t", list(0, .5, 30), list(-Inf, 2))
  expect_equal(unname(moments(t_trunc)), unname(reference(function(x) pdf(t_trunc, x), -Inf, 2, anchor = 1)), tolerance = 1e-8)
  # gamma and exponential with enough rate: finite exponential moments
  expect_equal(unname(moments(prior("gamma", list(2, 3)))), c(2.25, (3 / (3 - 2))^2 - 2.25^2), tolerance = 1e-8)
  expect_equal(unname(moments(prior("exp", list(3)))), c(3 / 2, 3 / 1 - (3 / 2)^2), tolerance = 1e-8)

  # right tails without (enough) exponential moments give infinite moments
  expect_identical(unname(moments(prior("t", list(0, .5, 30)))), c(Inf, Inf))
  expect_identical(unname(moments(prior("lognormal", list(0, .5)))), c(Inf, Inf))
  expect_identical(unname(moments(prior("invgamma", list(3, 2)))), c(Inf, Inf))
  gamma_moments <- moments(prior("gamma", list(2, 1.5)))
  expect_equal(unname(gamma_moments[1]), 9, tolerance = 1e-8)
  expect_identical(unname(gamma_moments[2]), Inf)
  expect_identical(unname(moments(prior("exp", list(1)))), c(Inf, Inf))
  # a truncated right tail has all moments
  expect_true(all(is.finite(moments(t_trunc))))

  # a moment that no integration resolves, or that no double holds, is refused
  inverse_moment <- prior("invmoment", list(location = 0, tau = 1, order = 1, df = 3))
  component <- list(type = "prior", scale = "log_omega", prior = inverse_moment)
  expect_error(.prior_weightfunction_component_mean(component), class = "BayesTools_numerical_unavailable")
  expect_error(.prior_weightfunction_component_var(component), class = "BayesTools_numerical_unavailable")
  wide <- list(type = "prior", scale = "log_omega", prior = prior("normal", list(0, 40)))
  expect_error(.prior_weightfunction_component_mean(wide), class = "BayesTools_numerical_range_limit")
  expect_error(.prior_weightfunction_component_var(wide), class = "BayesTools_numerical_range_limit")

  # the plotted mean is the same quantity
  p_wf <- prior_weightfunction("one-sided", .05, wf_independent(prior("normal", list(2, .01)), scale = "log_omega"))
  plotted <- .weightfunction_prior_marginal_components(.weightfunction_prior_list_context(list(p_wf)), 2)[[1]]
  expect_equal(.weightfunction_component_mean(plotted), mean(p_wf)[[2]])
  heavy <- list(type = "prior", scale = "log_omega", prior = prior("t", list(0, .5, 30)), weight = 1)
  expect_identical(.weightfunction_component_mean(heavy), Inf)
  expect_identical(.weightfunction_mixture_mean(list(heavy, list(type = "point", location = 1, weight = 0))), Inf)
  expect_identical(.weightfunction_mixture_mean(list(modifyList(heavy, list(weight = 0)), list(type = "point", location = 1, weight = 1))), 1)
})

test_that("cumulative weightfunction priors expose exact beta marginals", {
  alpha <- c(1, 2, 3)
  p_wf <- prior_weightfunction(
    side    = "one-sided",
    steps   = c(.025, .05),
    weights = wf_cumulative(alpha)
  )

  q <- c(0, .2, .8, 1)
  expected_cdf <- cbind(
    ppoint(q, location = 1),
    stats::pbeta(q, shape1 = 5, shape2 = 1),
    stats::pbeta(q, shape1 = 3, shape2 = 3)
  )
  expected_lpdf <- cbind(
    dpoint(q, location = 1, log = TRUE),
    stats::dbeta(q, shape1 = 5, shape2 = 1, log = TRUE),
    stats::dbeta(q, shape1 = 3, shape2 = 3, log = TRUE)
  )

  expect_equal(mcdf(p_wf, q), expected_cdf, tolerance = 1e-12)
  expect_equal(mccdf(p_wf, q), 1 - expected_cdf, tolerance = 1e-12)
  expect_equal(mlpdf(p_wf, q), expected_lpdf, tolerance = 1e-12)
  expect_equal(mpdf(p_wf, q), exp(expected_lpdf), tolerance = 1e-12)

  probs <- c(.25, .5, .75)
  expect_equal(
    mquant(p_wf, probs),
    cbind(
      qpoint(probs, location = 1),
      stats::qbeta(probs, shape1 = 5, shape2 = 1),
      stats::qbeta(probs, shape1 = 3, shape2 = 3)
    ),
    tolerance = 1e-12
  )
  expect_equal(mean(p_wf), c(1, 5 / 6, 1 / 2), tolerance = 1e-12)
  expect_equal(var(p_wf), c(0, 5 / 252, 9 / 252), tolerance = 1e-12)
})

test_that("weightfunction marginal CCDFs preserve extreme upper tails", {
  p_cumulative <- prior_weightfunction(
    side    = "one-sided",
    steps   = .05,
    weights = wf_cumulative(c(100, 1))
  )
  beta_tail <- mccdf(p_cumulative, .9)[1, 2]
  beta_expected <- stats::pbeta(
    .9, shape1 = 1, shape2 = 100, lower.tail = FALSE
  )

  expect_gt(beta_tail, 0)
  expect_equal(beta_tail / beta_expected, 1, tolerance = 1e-12)

  p_omega <- prior_weightfunction(
    side    = "one-sided",
    steps   = .05,
    weights = wf_independent(prior("beta", list(1, 100)))
  )
  omega_tail <- mccdf(p_omega, .9)[1, 2]

  expect_gt(omega_tail, 0)
  expect_equal(omega_tail / beta_expected, 1, tolerance = 1e-12)

  p_log_omega <- prior_weightfunction(
    side    = "one-sided",
    steps   = .05,
    weights = wf_independent(
      prior("normal", list(mean = 0, sd = 1)),
      scale = "log_omega"
    )
  )
  log_omega_tail <- mccdf(p_log_omega, exp(10))[1, 2]
  log_omega_expected <- stats::pnorm(10, lower.tail = FALSE)

  expect_gt(log_omega_tail, 0)
  expect_equal(
    log_omega_tail / log_omega_expected,
    1,
    tolerance = 1e-12
  )
})

test_that("the literal of a numeric constant is the shortest text that parses back", {

  # the rule: the as.character() text where it parses back, otherwise the
  # shortest of 16 and 17 significant digits that does
  expect_identical(.JAGS_prior_number(0.5), "0.5")
  expect_identical(.JAGS_prior_number(c(0, 1, -2.5, 100, 1e5, 1e-5, 0.1)), c("0", "1", "-2.5", "100", "1e+05", "1e-05", "0.1"))
  expect_identical(.JAGS_prior_number(1/3), "0.3333333333333333")
  expect_identical(.JAGS_prior_number(0.1 + 0.2), "0.30000000000000004")
  expect_identical(.JAGS_prior_number(1e8 + 1e-7), "100000000.0000001")
  expect_identical(.JAGS_prior_number("expression"), "expression")
  expect_identical(.JAGS_prior_number(numeric()), character())
  hard <- c(1/3, 2/3, 0.1 + 0.2, 1e8 + 1e-7, pi, -exp(1), sqrt(2) * 1e10, 1e-300 / 3, 1e22 + 1e6, 123456789.123456789, 5L)
  expect_identical(as.numeric(.JAGS_prior_number(hard)), as.numeric(hard))
})

test_that("numeric constants of the emitted JAGS syntax parse back to the declared doubles", {

  # `number` is the first capture of `pattern` in the syntax of `prior_list`
  number <- "(-?[0-9][0-9.]*(?:e[+-]?[0-9]+)?)"
  captured <- function(prior_list, pattern){
    syntax  <- JAGS_add_priors("model{}", prior_list)
    pattern <- gsub("<num>", number, pattern, fixed = TRUE)
    matches <- regmatches(syntax, regexec(pattern, syntax, perl = TRUE))[[1]]
    expect_gt(length(matches), 1L)
    as.numeric(matches[-1])
  }

  # a spike at 1/3 is sampled as the declared atom; 0.5 keeps its text
  for(location in c(1/3, 2/3, 0.1 + 0.2, 1e8 + 1e-7, pi, -1/7, 1e-300, 0.5, 2, 100)){
    expect_identical(captured(list(x = prior("spike", list(location))), "x = <num>\n"), location)
  }
  expect_identical(
    JAGS_add_priors("model{}", list(x = prior("spike", list(.5)))),
    "model{\n x = 0.5\n\n}"
  )

  # means, precisions and truncation bounds
  expect_identical(
    captured(list(x = prior("normal", list(1e8 + 1e-7, 1e-8))), "x ~ dnorm\\(<num>,<num>\\)"),
    c(1e8 + 1e-7, 1 / (1e-8)^2)
  )
  expect_identical(
    JAGS_add_priors("model{}", list(x = prior("normal", list(0, 1)))),
    "model{\n x ~ dnorm(0,1)\n\n}"
  )
  expect_identical(
    captured(list(x = prior("normal", list(0, 1), truncation = list(1/3, 2/3))), "x ~ dnorm\\(0,1\\)T\\(<num>,<num>\\)"),
    c(1/3, 2/3)
  )
  expect_identical(
    captured(list(x = prior("gamma", list(1/3, 2/7))), "x ~ dgamma\\(<num>,<num>\\)"),
    c(1/3, 2/7)
  )
  expect_identical(
    captured(list(x = prior("t", list(1/3, 1/7, 2/3))), "x ~ dt\\(<num>,<num>,<num>\\)"),
    c(1/3, 1 / (1/7)^2, 2/3)
  )

  # vector priors
  expect_identical(
    captured(list(x = prior("mnormal", list(mean = 1/3, sd = 1/7, K = 2))), "prior_par1_x = rep\\(<num>,2\\)"),
    1/3
  )
  expect_identical(
    captured(list(x = prior("mnormal", list(mean = 1/3, sd = 1/7, K = 2))), "prior_par2_x\\[i,i\\] <- <num>"),
    1 / (1/7)^2
  )
  mt <- prior("mt", list(location = 1/3, scale = 1/7, df = 3.3, K = 2))
  expect_identical(captured(list(x = mt), "prior_par_s_x ~ dgamma\\(<num>, <num>\\)"), rep(3.3 / 2, 2))
  expect_identical(captured(list(x = mt), "prior_par2_x\\[i,i\\] <- <num>"), 1 / (1/7)^2)
  expect_identical(captured(list(x = mt), "\\) \\+ <num> \n"), 1/3)
  expect_identical(
    captured(list(x = prior("dirichlet", list(alpha = c(1/3, 2/3, 1)))), "prior_par_eta_x\\[2\\] ~ dgamma\\(<num>, 1\\)"),
    2/3
  )
  expect_identical(
    captured(list(x = prior("mpoint", list(location = 1/3, K = 2))), "x\\[i\\] = <num> \n"),
    1/3
  )

  # the weights of a mixture
  weights <- c(1/3, 2/3)
  normals <- list(
    prior("normal", list(0, 1), prior_weights = weights[1]),
    prior("normal", list(1, 1), prior_weights = weights[2])
  )
  expect_identical(
    captured(list(x = prior_mixture(normals)), "x_indicator ~ dcat\\(c\\(<num>, <num>\\)\\)"),
    weights
  )
  pet_peese <- prior_mixture(list(
    prior_PET("normal", list(0, 1), prior_weights = weights[1]),
    prior_PEESE("normal", list(0, 1), prior_weights = weights[2])
  ))
  expect_identical(captured(list(bias = pet_peese), "bias_indicator ~ dcat\\(c\\(<num>, <num>\\)\\)"), weights)
  selection <- prior_mixture(list(
    prior_none(prior_weights = weights[1]),
    prior_weightfunction("one-sided", c(.05), wf_fixed(c(1, .5)), prior_weights = weights[2])
  ))
  expect_identical(captured(list(bias = selection), "bias_indicator ~ dcat\\(c\\(<num>, <num>\\)\\)"), weights)

  # fixed weights and the shapes of cumulative weights
  fixed <- prior_weightfunction("one-sided", c(.05, .1), wf_fixed(c(1, 1/3, 0.2)))
  expect_identical(captured(list(omega = fixed), "omega\\[2\\] <- <num>\n"), 1/3)
  expect_identical(captured(list(omega = fixed), "omega\\[3\\] <- <num>\n"), 0.2)
  cumulative <- prior_weightfunction("one-sided", c(.05, .1), wf_cumulative(c(1/3, 0.7, 1)))
  expect_identical(captured(list(omega = cumulative), "eta\\[1\\] ~ dgamma\\(<num>, 1\\)"), 1/3)
  binary <- prior_weightfunction("one-sided", c(.05), wf_cumulative(c(1/3, 2/3)))
  beta   <- .weightfunction_alpha_marginal(c(1/3, 2/3), 2L)
  expect_identical(
    captured(list(omega = binary), "omega_ratio ~ dbeta\\(<num>, <num>\\)"),
    c(beta$alpha, beta$beta)
  )
})

test_that("expression-valued priors are found in mixtures, factor priors, and ordered totals", {

  expression_normal <- prior("normal", list(0, expression(tau)))
  expression_point  <- prior("point", list(location = expression(2 * theta)))
  numeric_normal    <- prior("normal", list(0, 1))
  numeric_point     <- prior("point", list(location = 0.25))

  # simple priors read their own parameters
  expect_true(.bt_prior_has_expression(expression_normal))
  expect_true(.bt_prior_has_expression(expression_point))
  expect_false(.bt_prior_has_expression(numeric_normal))
  expect_false(.bt_prior_has_expression(numeric_point))
  expect_false(.bt_prior_has_expression(prior_none()))
  expect_false(.bt_prior_has_expression(NULL))
  expect_false(.bt_prior_has_expression(1))

  # mixture, spike-and-slab, and model lists look at every component
  expect_true(.bt_prior_has_expression(prior_mixture(
    list(numeric_normal, expression_point), is_null = c(FALSE, FALSE))))
  expect_false(.bt_prior_has_expression(prior_mixture(
    list(numeric_normal, numeric_point), is_null = c(FALSE, FALSE))))
  expect_true(.bt_prior_has_expression(prior_spike_and_slab(
    expression_normal, prior("spike", list(0.5)))))
  expect_false(.bt_prior_has_expression(prior_spike_and_slab(
    numeric_normal, prior("spike", list(0.5)))))
  expect_true(.bt_prior_has_expression(list(numeric_normal, expression_normal)))
  expect_false(.bt_prior_has_expression(list(numeric_normal, numeric_point)))

  # factor priors, factor mixtures, and factor spike-and-slab priors
  factor_expression <- prior_factor("mnormal", list(0, expression(tau)), contrast = "meandif")
  factor_numeric    <- prior_factor("mnormal", list(0, 1), contrast = "meandif")
  expect_true(.bt_prior_has_expression(factor_expression))
  expect_false(.bt_prior_has_expression(factor_numeric))
  expect_true(.bt_prior_has_expression(prior_mixture(
    list(factor_expression, factor_numeric), is_null = c(FALSE, FALSE))))
  expect_true(.bt_prior_has_expression(prior_spike_and_slab(
    factor_expression, prior("spike", list(0.5)))))
  expect_false(.bt_prior_has_expression(prior_mixture(
    list(factor_numeric, factor_numeric), is_null = c(FALSE, FALSE))))

  # an ordered prior has an expression exactly when its total has one
  expect_true(.bt_prior_has_expression(prior_ordered(expression_point)))
  expect_true(.bt_prior_has_expression(prior_ordered(expression_normal)))
  expect_true(.bt_prior_has_expression(prior_ordered(prior_spike_and_slab(
    expression_normal, prior("spike", list(0.5))))))
  expect_false(.bt_prior_has_expression(prior_ordered(numeric_normal)))
  expect_false(.bt_prior_has_expression(prior_ordered(numeric_point, allocation = c(.25, .25, .5))))
})

test_that("only a numeric zero point is a spike at zero", {

  expect_true(.is_prior_point_zero(prior("point", list(location = 0))))
  expect_true(.is_prior_point_zero(prior("spike", list(location = 0))))
  expect_false(.is_prior_point_zero(prior("point", list(location = 0.25))))
  expect_false(.is_prior_point_zero(prior("point", list(location = expression(2 * theta)))))
  expect_false(.is_prior_point_zero(prior("point", list(location = expression(0)))))
  expect_false(.is_prior_point_zero(prior("normal", list(0, 1))))
  expect_false(.is_prior_point_zero(prior_none()))
  expect_true(.is_prior_point_zero(prior("mpoint", list(0, 2))))
  expect_false(.is_prior_point_zero(prior("mpoint", list(1, 2))))
})
