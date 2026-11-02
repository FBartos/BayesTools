skip_if_not_test_profile("unit")
source(testthat::test_path("common-functions.R"))

# ============================================================================ #
# TEST FILE: Hypotheses on the levels of ordered terms in interactions
# ============================================================================ #
#
# The fits are synthetic: their posterior draws are draws from the priors of
# the model (totals, indicators, allocations and the coefficients they give,
# from .prior_ordered_draws(), which does not use any prior-density code), so
# the posterior is the prior and every Bayes factor of a statement is 1 up to
# Monte Carlo and kernel-density error.
#
# ============================================================================ #

N01 <- prior("normal", list(0, 1))
treatment <- prior_factor("normal", list(0, 1), contrast = "treatment")

.interaction_fit_data <- function(){

  data.frame(
    f = ordered(rep(c("lo", "mid", "hi"), 6), levels = c("lo", "mid", "hi")),
    o = ordered(rep(c("p", "q", "r"), each = 6), levels = c("p", "q", "r")),
    g = factor(rep(c("A", "B", "C"), each = 6)),
    h = factor(rep(c("U", "V"), 9)),
    x = sin(seq_len(18))
  )
}

# A fit whose draws are prior draws; 'formula_scale' marks scaled predictors.
.interaction_fit <- function(formula, priors, n = 20000L, seed = 1L, formula_scale = NULL){

  info <- JAGS_formula(formula, "mu", .interaction_fit_data(), priors, formula_scale = formula_scale)
  bound <- info$prior_list
  set.seed(seed)
  columns <- list()
  shares <- list()
  for(name in names(bound)){
    prior <- bound[[name]]
    if(!is.prior.ordered(prior)){
      names_ <- if(is.prior.factor(prior)) .JAGS_prior_factor_names(name, prior) else name
      columns[[name]] <- matrix(stats::rnorm(n * length(names_)), n, length(names_),
                                dimnames = list(NULL, names_))
      next
    }
    draws <- .prior_ordered_draws(prior, n)
    coefficients <- draws$coefficients
    colnames(coefficients) <- .JAGS_prior_factor_names(name, prior)
    columns[[name]] <- cbind(coefficients, .prior_ordered_total_samples(draws, name))
    for(record in .prior_ordered_dirichlet_records(prior)){
      eta_name <- .JAGS_prior_dirichlet_eta_name(record$node)
      if(eta_name %in% names(shares)) next
      shares[[eta_name]] <- draws$allocation_samples[[record$key]]
      columns[[eta_name]] <- matrix(shares[[eta_name]] * stats::rexp(n), n, record$dim,
        dimnames = list(NULL, paste0(eta_name, "[", seq_len(record$dim), "]")))
    }
  }
  values <- do.call(cbind, unname(columns))
  fit <- structure(
    list(mcmc = coda::mcmc.list(coda::mcmc(values)), sample = n,
         summary.pars = list(mutate = NULL), monitor = colnames(values)),
    class = c("runjags", "BayesTools_fit", "list")
  )
  attr(fit, "prior_list")     <- bound
  attr(fit, "formula_design") <- list(mu = info$formula_design)
  if(length(info$formula_scale) > 0L) attr(fit, "formula_scale") <- list(mu = info$formula_scale)
  attach_test_parameter_map(fit)
}

.interaction_bf <- function(fit, statement, ...){

  attr(suppressWarnings(hypothesis_BF(fit, statement, seed = 1, ...)), "raw_BF", exact = TRUE)
}

test_that("the levels of an ordered interaction are tested within and across slices", {

  slab <- prior_spike_and_slab(N01, prior("point", list(.5)))
  for(total in list(N01, slab)){
    conditional <- is.prior.spike_and_slab(total)
    fit <- .interaction_fit(~ f * g, list(intercept = N01, f = prior_ordered(total), g = treatment,
                                          "f:g" = prior_ordered(total)), n = 40000L, seed = 31L)
    for(statement in c(
      # cells and differences within a slice
      "f:g[f=hi, g=B] > 0", "f:g[f=mid, g=C] > 0.1", "f:g[f=hi, g=C] > f:g[f=mid, g=C]",
      "f:g[f=hi, g=B] - f:g[f=mid, g=B] > 0.2", "f:g[f=mid, g=B] = 0.2", "f:g[f=hi, g=C] = 0.3",
      # across slices: independent totals and allocations
      "f:g[f=mid, g=B] > f:g[f=mid, g=C]", "f:g[f=hi, g=B] - f:g[f=mid, g=C] = 0.1",
      "f:g[f=hi, g=B] = f:g[f=hi, g=C]")){
      expect_equal(.interaction_bf(fit, statement, conditional = conditional), 1, tolerance = 0.1,
                   info = paste(statement, conditional))
    }
  }
})

test_that("the levels of other interaction layouts are tested from their prior", {

  slab <- prior_spike_and_slab(N01, prior("point", list(.5)))
  # the ordinary factor first, and two ordinary factors (one slice per coordinate pair)
  reversed <- .interaction_fit(~ g * f, list(intercept = N01, g = treatment, f = prior_ordered(N01),
                                             "g:f" = prior_ordered(N01)), n = 30000L, seed = 32L)
  for(statement in c("g:f[g=B, f=hi] > 0", "g:f[g=B, f=hi] > g:f[g=C, f=hi]",
                     "g:f[g=B, f=hi] - g:f[g=B, f=mid] = 0.15")){
    expect_equal(.interaction_bf(reversed, statement), 1, tolerance = 0.1, info = statement)
  }
  threeway <- .interaction_fit(~ f * g * h, list(intercept = N01, f = prior_ordered(N01), g = treatment,
    h = treatment, "f:g" = prior_ordered(N01), "f:h" = prior_ordered(N01), "g:h" = treatment,
    "f:g:h" = prior_ordered(N01)), n = 30000L, seed = 33L)
  for(statement in c("f:g:h[f=hi, g=B, h=V] > f:g:h[f=hi, g=C, h=V]",
                     "f:g:h[f=hi, g=B, h=V] - f:g:h[f=mid, g=B, h=V] = 0.2")){
    expect_equal(.interaction_bf(threeway, statement), 1, tolerance = 0.1, info = statement)
  }
  # two ordered factors: the cell is a total times two shares
  crossed <- .interaction_fit(~ f * o, list(intercept = N01, f = prior_ordered(N01), o = prior_ordered(N01),
                                            "f:o" = prior_ordered(N01)), n = 30000L, seed = 34L)
  for(statement in c("f:o[f=mid, o=q] > 0", "f:o[f=mid, o=q] = 0.2", "f:o[f=hi, o=q] > f:o[f=mid, o=q]",
                     "f:o[f=hi, o=q] - f:o[f=mid, o=q] > 0.05", "f:o[f=hi, o=q] - f:o[f=hi, o=p] = 0.1",
                     "f:o[f=mid, o=q] - f:o[f=lo, o=r] > 0")){
    expect_equal(.interaction_bf(crossed, statement), 1, tolerance = 0.1, info = statement)
  }
  # a shared allocation: the slices are one set of coordinates
  shared <- .interaction_fit(~ f * g, list(intercept = N01, f = prior_ordered(N01), g = treatment,
    "f:g" = prior_ordered(N01, id = "shape")), n = 30000L, seed = 35L)
  for(statement in c("f:g[f=mid, g=B] > f:g[f=mid, g=C]", "f:g[f=hi, g=B] - f:g[f=mid, g=C] = 0.1")){
    expect_equal(.interaction_bf(shared, statement), 1, tolerance = 0.1, info = statement)
  }
  # a spike-and-slab total with the shared indicator, tested within its inclusion event
  slab_fit <- .interaction_fit(~ f * g, list(intercept = N01, f = prior_ordered(slab), g = treatment,
    "f:g" = prior_ordered(slab, id = "shape")), n = 30000L, seed = 36L)
  expect_equal(.interaction_bf(slab_fit, "f:g[f=hi, g=B] - f:g[f=mid, g=C] = 0.1", conditional = TRUE), 1,
               tolerance = 0.1)
})

test_that("combinations without a route are refused with the measure condition, and tables complete", {

  heavy <- prior("t", list(0, 1, 3))
  sliced <- .interaction_fit(~ f * g, list(intercept = N01, f = prior_ordered(N01), g = treatment,
    "f:g" = prior_ordered(heavy)), n = 4000L)
  # several totals that are not normal (two slices of a t total)
  expect_error(.interaction_bf(sliced, "f:g[f=mid, g=B] > f:g[f=mid, g=C]"),
               class = "BayesTools_formula_measure_unavailable")
  # one t total is a scale product
  expect_true(is.finite(.interaction_bf(sliced, "f:g[f=mid, g=B] > 0")))
  # the table of statements completes on layouts with and without routes
  crossed <- .interaction_fit(~ f * o, list(intercept = N01, f = prior_ordered(N01), o = prior_ordered(N01),
                                            "f:o" = prior_ordered(N01)), n = 4000L)
  shared <- .interaction_fit(~ f * g, list(intercept = N01, f = prior_ordered(N01), g = treatment,
    "f:g" = prior_ordered(N01, id = "shape")), n = 4000L)
  for(case in list(list(sliced, "f:g"), list(crossed, "f:o"), list(shared, "f:g"))){
    quantities <- hypothesis_quantities(case[[1L]])
    expect_s3_class(quantities, "data.frame")
    expect_true(case[[2L]] %in% quantities$term)
  }

  # transformed prior densities of a model with scaled predictors
  scaled <- .interaction_fit(~ f * o + x, list(intercept = N01, f = prior_ordered(N01),
    o = prior_ordered(N01), "f:o" = prior_ordered(N01), x = N01), n = 2000L,
    formula_scale = list(x = TRUE))
  mixed <- as_mixed_posteriors(scaled, parameters = "mu_f__xXx__o", transform_scaled = TRUE)
  densities <- posterior_metadata(mixed[["mu_f__xXx__o"]], "prior_densities")
  expect_true(all(vapply(densities, inherits, logical(1), what = "prior_linear_density")))
  heavy_scaled <- .interaction_fit(~ f * g + x, list(intercept = N01, f = prior_ordered(N01), g = treatment,
    "f:g" = prior_ordered(heavy), x = N01), n = 2000L, formula_scale = list(x = TRUE))
  expect_no_error(as_mixed_posteriors(heavy_scaled, parameters = "mu_f__xXx__g", transform_scaled = TRUE))
})
