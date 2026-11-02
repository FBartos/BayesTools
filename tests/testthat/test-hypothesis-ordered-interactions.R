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

test_that("the levels of an ordered interaction are tested within their slice", {

  slab <- prior_spike_and_slab(N01, prior("point", list(.5)))
  for(total in list(N01, slab)){
    conditional <- is.prior.spike_and_slab(total)
    fit <- .interaction_fit(~ f * g, list(intercept = N01, f = prior_ordered(total), g = treatment,
                                          "f:g" = prior_ordered(total)), n = 40000L, seed = 31L)
    for(statement in c("f:g[f=hi, g=B] > 0", "f:g[f=mid, g=C] > 0.1", "f:g[f=hi, g=C] > f:g[f=mid, g=C]",
                       "f:g[f=hi, g=B] - f:g[f=mid, g=B] > 0.2", "f:g[f=mid, g=B] = 0.2",
                       "f:g[f=hi, g=C] = 0.3")){
      expect_equal(.interaction_bf(fit, statement, conditional = conditional), 1, tolerance = 0.1,
                   info = paste(statement, conditional))
    }
  }
})

test_that("ordered interactions degrade to measure refusals instead of aborting", {

  sliced <- .interaction_fit(~ f * g, list(intercept = N01, f = prior_ordered(N01), g = treatment,
                                           "f:g" = prior_ordered(N01)), n = 4000L)
  crossed <- .interaction_fit(~ f * o, list(intercept = N01, f = prior_ordered(N01),
    o = prior_ordered(N01), "f:o" = prior_ordered(N01)), n = 4000L)
  # a statement across slices, and a level of two ordered factors, are refused with the typed condition
  expect_error(.interaction_bf(sliced, "f:g[f=mid, g=B] > f:g[f=mid, g=C]"),
               class = "BayesTools_formula_measure_unavailable")
  expect_error(.interaction_bf(crossed, "f:o[f=mid, o=q] > 0"),
               class = "BayesTools_formula_measure_unavailable")
  # the table of the statements completes
  for(case in list(list(sliced, "f:g"), list(crossed, "f:o"))){
    quantities <- hypothesis_quantities(case[[1L]])
    expect_s3_class(quantities, "data.frame")
    expect_true(case[[2L]] %in% quantities$term)
  }

  # transformed prior densities of a model with scaled predictors report the
  # unavailable levels instead of failing
  scaled <- .interaction_fit(~ f * o + x, list(intercept = N01, f = prior_ordered(N01),
    o = prior_ordered(N01), "f:o" = prior_ordered(N01), x = N01), n = 2000L,
    formula_scale = list(x = TRUE))
  expect_no_error(as_mixed_posteriors(scaled, parameters = "mu_f__xXx__o", transform_scaled = TRUE))
})
