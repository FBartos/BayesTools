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

.expression_cases <- list(
  list(fit = "fit_dnode_point_expression", parameter = "mu_x"),
  list(fit = "fit_expression_simple", parameter = "x"),
  list(fit = "fit_expression_mixture", parameter = "x"),
  list(fit = "fit_expression_spike_slab", parameter = "x")
)

test_that("the catalog and the mixed posterior of an expression prior keep the draws and omit the prior law", {

  for(case in .expression_cases){
    fit <- .expression_fit(case$fit)
    draws <- unname(as.matrix(.fit_to_posterior(fit))[, case$parameter])
    selection <- parameter_catalog_resolve(parameter_catalog(fit), case$parameter)

    expect_null(parameter_prior_density(fit, selection), info = case$fit)
    mixed <- parameter_mixed_posterior(fit, selection)
    expect_equal(as.numeric(mixed), draws, tolerance = 0, info = case$fit)
    expect_null(posterior_metadata(mixed, "prior_density"), info = case$fit)
    expect_null(posterior_metadata(mixed, "measure_unavailable"), info = case$fit)
    # the declared atoms are the point masses of the draws
    if(identical(case$fit, "fit_expression_spike_slab")){
      atoms <- .posterior_atoms_get(mixed)
      expect_identical(as.numeric(atoms$locations), 0, info = case$fit)
      expect_equal(atoms$mass, mean(draws == 0), tolerance = 1e-12, info = case$fit)
    }else if(!identical(case$fit, "fit_dnode_point_expression")){
      expect_length(.posterior_atoms_get(mixed)$mass, 0L)
    }

    # the control parameter keeps its density
    control <- parameter_catalog_resolve(parameter_catalog(fit), if(identical(case$parameter, "mu_x")) "theta" else "x_sigma")
    expect_s3_class(parameter_prior_density(fit, control), "prior_linear_density")
  }
})

test_that("marginal posteriors and transformed draws of an expression prior record the unavailable law", {

  for(case in .expression_cases){
    fit <- .expression_fit(case$fit)
    parameter <- case$parameter
    samples <- as_mixed_posteriors(fit, parameter)
    draws <- as.numeric(samples[[parameter]])
    expect_unavailable <- function(x, info){
      row <- .expression_measure_row(x, "prior_density")
      expect_identical(row$column, parameter, info = info)
      expect_identical(row$cause, "expression_prior_law_unavailable", info = info)
      expect_null(posterior_metadata(x, "prior_density"), info = info)
      refusal <- tryCatch(posterior_measure_check(x, "prior_density"), error = identity)
      expect_s3_class(refusal, "BayesTools_formula_measure_unavailable")
      expect_identical(refusal$reason, "expression_prior_law_unavailable", info = info)
    }

    marginal <- marginal_posterior(samples, parameter, prior_samples = TRUE, use_formula = FALSE)
    expect_equal(as.numeric(marginal), draws, tolerance = 0, info = case$fit)
    expect_unavailable(marginal, case$fit)
    marginal <- marginal_posterior(samples, parameter, use_formula = FALSE, transformation = "lin",
      transformation_arguments = list(a = 1, b = 2))
    expect_equal(as.numeric(marginal), 1 + 2 * draws, tolerance = 1e-12, info = case$fit)
    expect_unavailable(marginal, case$fit)
    transformed <- posterior_transform(samples[[parameter]], "lin", list(a = 1, b = 2))
    expect_equal(as.numeric(transformed), 1 + 2 * draws, tolerance = 1e-12, info = case$fit)
    expect_unavailable(transformed, case$fit)

    # the samples of an unscaled predictor are their own original-scale samples
    scaled <- as_mixed_posteriors(fit, parameter, transform_scaled = TRUE)
    expect_equal(as.numeric(scaled[[parameter]]), draws, tolerance = 0, info = case$fit)
  }
})

test_that("plot_posterior() draws the posterior of an expression prior without its prior curve", {

  for(case in .expression_cases){
    fit <- .expression_fit(case$fit)
    parameter <- case$parameter
    samples <- as_mixed_posteriors(fit, parameter)
    plot <- NULL
    expect_warning(plot <- plot_posterior(samples, parameter, prior = TRUE, plot_type = "ggplot"),
      class = "BayesTools_prior_curve_unavailable")
    expect_no_warning(posterior_only <- plot_posterior(samples, parameter, prior = FALSE, plot_type = "ggplot"))
    expect_identical(ggplot2::ggplot_build(plot)$data, ggplot2::ggplot_build(posterior_only)$data, info = case$fit)

    # the base graphics path draws too
    grDevices::pdf(NULL)
    expect_warning(plot_posterior(samples, parameter, prior = TRUE), class = "BayesTools_prior_curve_unavailable")
    grDevices::dev.off()

    # draws whose only refusal is the unavailable prior law are drawn, too
    transformed <- samples
    transformed[[parameter]] <- posterior_transform(samples[[parameter]], "lin", list(a = 1, b = 2))
    expect_warning(plot_posterior(transformed, parameter, prior = TRUE, plot_type = "ggplot"),
      class = "BayesTools_prior_curve_unavailable")
  }
})

test_that("a hypothesis on the prior of an expression prior refuses with its typed cause", {

  for(case in .expression_cases){
    fit <- .expression_fit(case$fit)
    samples <- as_mixed_posteriors(fit, case$parameter)
    prior_object <- attr(fit, "prior_list")[[case$parameter]]
    for(hypothesis in c("x = 0", "x > 0")){
      refusal <- tryCatch(hypothesis_BF(as.numeric(samples[[case$parameter]]), prior = prior_object,
        hypothesis = hypothesis), error = identity)
      expect_s3_class(refusal, "BayesTools_formula_prior_density_unavailable")
      expect_s3_class(refusal, "BayesTools_formula_measure_unavailable")
      expect_identical(refusal$reason, "expression_prior_law_unavailable", info = case$fit)
    }
  }
})

test_that("hypotheses on a marginal posterior without a prior law refuse with the recorded cause", {

  fit <- .expression_fit("fit_dnode_point_expression")
  samples <- as_mixed_posteriors(fit, c("mu_intercept", "mu_x", "theta"))
  levels <- marginal_posterior(samples, "mu_intercept", formula = ~ x, at = list(x = 1), prior_samples = TRUE)
  expect_null(posterior_metadata(levels[[1L]], "prior_density"))
  for(hypothesis in c("mu_intercept = 0.1", "mu_intercept > 0", "mu_intercept > 0 vs mu_intercept < 0")){
    refusal <- tryCatch(hypothesis_BF(levels, hypothesis = hypothesis, parameter = "mu_intercept"), error = identity)
    expect_s3_class(refusal, "BayesTools_formula_measure_unavailable")
    expect_identical(refusal$reason, "expression_prior_law_unavailable", info = hypothesis)
    expect_false(grepl("is required", conditionMessage(refusal), fixed = TRUE))
  }
  # a marginal without any record of its prior keeps the plain message
  plain <- levels
  plain[[1L]] <- .bt_meta_set(plain[[1L]], "measure_unavailable", NULL)
  expect_error(hypothesis_BF(plain, hypothesis = "mu_intercept = 0.1", parameter = "mu_intercept"),
               "Prior density is required for point hypotheses.", fixed = TRUE)
  expect_error(hypothesis_BF(plain, hypothesis = "mu_intercept > 0", parameter = "mu_intercept"),
               "Prior information is required for region hypotheses.", fixed = TRUE)
})

test_that("summary tables list the expression point of a fitted model", {

  fit <- .expression_fit("fit_dnode_point_expression")
  model <- list(fit = fit, inference = list(m_number = 1L, marglik = 0, prior_prob = 1, post_prob = 1, inclusion_BF = 1))
  text <- paste(utils::capture.output(print(model_summary_table(model))), collapse = "\n")
  expect_match(text, "(mu) x ~ Spike(2 * theta + delta)", fixed = TRUE)
  expect_false(grepl("intercept", text, fixed = TRUE))
  table <- ensemble_summary_table(list(model), c(mu_x = "mu_x"))
  expect_identical(as.data.frame(table)[["(mu) x"]], "Spike(2 * theta + delta)")
})

test_that("the estimates of a scaled expression point are reported on both scales", {

  fit <- .expression_fit("fit_dnode_point_expression_scaled")
  expect_true(.bt_formula_scale_has_scaled(attr(fit, "formula_scale")))
  posterior <- as.matrix(.fit_to_posterior(fit))
  x <- attr(fit, "formula_design")$mu$source_data$x
  fitted <- JAGS_estimates_table(fit)
  original <- JAGS_estimates_table(fit, transform_scaled = TRUE)

  expect_equal(fitted["(mu) x", "Mean"], mean(posterior[, "mu_x"]), tolerance = 1e-8)
  expect_equal(fitted["(mu) intercept", "Mean"], mean(posterior[, "mu_intercept"]), tolerance = 1e-8)
  # the slope and the intercept on the original predictor scale
  slope <- posterior[, "mu_x"] / stats::sd(x)
  intercept <- posterior[, "mu_intercept"] - posterior[, "mu_x"] * mean(x) / stats::sd(x)
  expect_equal(original["(mu) x", "Mean"], mean(slope), tolerance = 1e-8)
  expect_equal(original["(mu) x", "SD"], stats::sd(slope), tolerance = 1e-8)
  expect_equal(original["(mu) intercept", "Mean"], mean(intercept), tolerance = 1e-8)
  expect_equal(original["(mu) intercept", "SD"], stats::sd(intercept), tolerance = 1e-8)
  expect_equal(original["(mu) x", "Mean"], mean(2 * posterior[, "theta"]) / stats::sd(x), tolerance = 1e-8)

  # the original-scale prior laws of the intercept and the slope are unavailable with their cause
  samples <- as_mixed_posteriors(fit, c("mu_intercept", "mu_x", "theta"), transform_scaled = TRUE)
  for(parameter in c("mu_intercept", "mu_x")){
    expect_identical(.expression_measure_row(samples[[parameter]], "prior_density")$cause, "expression_prior_law_unavailable")
  }
  expect_null(posterior_metadata(samples$theta, "measure_unavailable"))
  for(target in c("mu_intercept", "mu_x")){
    refusal <- tryCatch(JAGS_formula_prior_density(fit, "mu", target), error = identity)
    expect_s3_class(refusal, "BayesTools_formula_measure_unavailable")
    expect_identical(refusal$reason, "expression_prior_law_unavailable")
  }
  # the model table of the fit
  model <- list(fit = fit, inference = list(m_number = 1L, marglik = 0, prior_prob = 1, post_prob = 1, inclusion_BF = 1))
  expect_match(paste(utils::capture.output(print(model_summary_table(model))), collapse = "\n"),
    "(mu) x ~ Spike(2 * theta)", fixed = TRUE)
})
