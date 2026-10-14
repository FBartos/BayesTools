skip_if_not_test_profile("fixture")

# ============================================================================ #
# TEST FILE: Priors with expression parameters on cached fits
# ============================================================================ #
#
# PURPOSE:
#   Priors whose parameters are expression()s are evaluated only by JAGS: they
#   have no R-side prior law. Every consumer of the prior density either returns
#   what does not need it (draws, support, atoms) with the prior density
#   recorded as unavailable (cause "expression_prior_law_unavailable") or
#   refuses with the typed condition of that reason; none raises an unclassed
#   error.
#
# MODELS/FIXTURES:
#   - fit_dnode_point_expression (expression point coefficient, formula mu = ~ x)
#   - fit_expression_simple, fit_expression_mixture, fit_expression_spike_slab
#     (simple expression priors) from test-00-model-fits.R
#
# TAGS: @fixture, @expression-priors
# ============================================================================ #

source(testthat::test_path("common-functions.R"))

.expression_fit <- function(name){

  skip_if_not_installed("rjags")
  skip_if_missing_fits(name)
  readRDS(file.path(temp_fits_dir, paste0(name, ".RDS")))
}

.expression_measure_row <- function(x, measure){

  table <- posterior_metadata(x, "measure_unavailable")
  table[table$measure == measure, , drop = FALSE]
}

test_that("formula marginals of an expression point record its prior law as unavailable with its cause", {

  fit <- .expression_fit("fit_dnode_point_expression")
  samples <- as_mixed_posteriors(fit, c("mu_intercept", "mu_x", "theta"))

  # the slope enters the marginal at x = 1: only its prior law (not the draws) is unavailable
  for(prior_samples in c(FALSE, TRUE)){
    marginal <- marginal_posterior(samples, "mu_intercept", formula = ~ x, at = list(x = 1),
                                   prior_samples = prior_samples)[[1L]]
    row <- .expression_measure_row(marginal, "prior_density")
    expect_identical(row$cause, "expression_prior_law_unavailable")
    expect_identical(row$reason, "The prior of 'mu_x' has parameters given as expressions, which only JAGS evaluates")
    refusal <- tryCatch(posterior_measure_check(marginal, "prior_density"), error = identity)
    expect_s3_class(refusal, "BayesTools_formula_measure_unavailable")
    expect_identical(refusal$reason, "expression_prior_law_unavailable")
    expect_null(posterior_metadata(marginal, "prior_density"))
    expect_length(as.numeric(marginal), nrow(as.matrix(.fit_to_posterior(fit))))
    # the contribution refusals of the atoms and the support are not the prior law
    expect_identical(.expression_measure_row(marginal, "atoms")$reason, "unsupported_contribution_measure")
    expect_identical(.expression_measure_row(marginal, "support")$reason, "unsupported_contribution_measure")
  }

  # at x = 0 the slope has no weight: the law of the intercept is available
  marginal <- marginal_posterior(samples, "mu_intercept", formula = ~ x, at = list(x = 0),
                                 prior_samples = TRUE)[[1L]]
  expect_null(posterior_metadata(marginal, "measure_unavailable"))
  expect_s3_class(posterior_metadata(marginal, "prior_density"), "prior_linear_density")

  # the formula prior-density interface refuses with the same typed cause
  for(call in list(quote(JAGS_formula_prior_density(fit, "mu", "mu_x")),
                   quote(JAGS_formula_prior_density(fit, "mu", weights = c(mu_intercept = 1, mu_x = 1))))){
    refusal <- tryCatch(eval(call), error = identity)
    expect_s3_class(refusal, "BayesTools_formula_prior_density_unavailable")
    expect_s3_class(refusal, "BayesTools_formula_measure_unavailable")
    expect_identical(refusal$reason, "expression_prior_law_unavailable")
  }
})
