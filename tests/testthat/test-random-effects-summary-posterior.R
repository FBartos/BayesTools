skip_if_not_test_profile("unit")

.random_effects_mean_variance_allocation_fit <- function(alpha = c(2, 3)){

  data <- data.frame(
    study = factor(c("s1", "s1", "s2", "s2")),
    x = c(-1, 0, 1, 2)
  )
  formula_result <- JAGS_formula(
    formula = ~ 1 + random(1 + x | study, name = "study", covariance = "diag"),
    parameter = "mu",
    data = data,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      allocation = random_variance_allocation(name = "allocation",
        terms = "study",
        target = "sd_component",
        scale = "mean_variance",
        sd = prior("gamma", list(2, 2)),
        weights = prior("dirichlet", list(alpha = alpha))
      )
    )
  )
  allocation <- formula_result$formula_design$random_effects[[1]]$sd_binding$allocations[[1]]
  samples <- matrix(
    c(
      2, 0.25, 0.75,
      2, 0.75, 0.25
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      NULL,
      c(
        allocation$source_node,
        paste0(allocation$weight_name, "[1]"),
        paste0(allocation$weight_name, "[2]")
      )
    )
  )
  fit <- structure(
    list(mcmc = coda::mcmc.list(coda::mcmc(samples)), sample = nrow(samples)),
    class = c("runjags", "BayesTools_fit", "list")
  )
  attr(fit, "prior_list") <- formula_result$prior_list
  attr(fit, "formula_design") <- list(mu = formula_result$formula_design)

  attach_test_parameter_map(fit)
}

.random_effects_total_variance_allocation_fit <- function(){

  data <- data.frame(
    study = factor(c("s1", "s1", "s2", "s2")),
    drug = factor(c("a", "b", "a", "b"))
  )
  formula_result <- JAGS_formula(
    formula = ~ 1 +
      random(1 | study, name = "study", covariance = "diag") +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = data,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      allocation = random_variance_allocation(name = "allocation",
        terms = c("study", "drug"),
        sd = prior("gamma", list(2, 2)),
        weights = prior("dirichlet", list(alpha = c(2, 3)))
      )
    )
  )
  allocation <- formula_result$formula_design$random_allocations[[1]]
  samples <- matrix(
    c(
      2, 0.25, 0.75,
      2, 0.75, 0.25
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      NULL,
      c(
        allocation$source_node,
        paste0(allocation$weight_name, "[1]"),
        paste0(allocation$weight_name, "[2]")
      )
    )
  )
  fit <- structure(
    list(mcmc = coda::mcmc.list(coda::mcmc(samples)), sample = nrow(samples)),
    class = c("runjags", "BayesTools_fit", "list")
  )
  attr(fit, "prior_list") <- formula_result$prior_list
  attr(fit, "formula_design") <- list(mu = formula_result$formula_design)

  attach_test_parameter_map(fit)
}

test_that("random-effect summary posterior extracts mean-variance ratios", {

  skip_if_not_installed("runjags")

  fit <- .random_effects_mean_variance_allocation_fit()
  testthat::local_mocked_bindings(
    .bt_random_effect_summary_samples = function(...) {
      stop("parallel random summary path was used", call. = FALSE)
    },
    .package = "BayesTools"
  )
  ratios <- random_effects_summary_posterior(fit, summary = "var_ratio")
  ratio_name <- "(mu) allocation: var_ratio(x)"

  expect_s3_class(ratios, "mixed_posteriors")
  expect_true(ratio_name %in% names(ratios))
  expect_s3_class(ratios[[ratio_name]], "marginal_posterior")
  expect_equal(unname(as.numeric(ratios[[ratio_name]])), c(1.5, 0.5), tolerance = 1e-12)
  ratio_selection <- parameter_catalog_resolve(
    parameter_catalog(fit),
    ratio_name,
    namespace = "mu"
  )
  expect_identical(
    parameter_transform(fit, ratio_selection),
    list(type = "affine", offset = 0, scale = 2)
  )
  supplied_samples <- as.matrix(fit$mcmc[[1L]])
  supplied_samples[, grep("weight\\[1\\]$", colnames(supplied_samples))] <-
    c(.1, .2)
  supplied_samples[, grep("weight\\[2\\]$", colnames(supplied_samples))] <-
    c(.9, .8)
  supplied_draws <- parameter_draws(
    fit,
    ratio_selection,
    model_samples = supplied_samples
  )
  expect_equal(
    unname(as.numeric(supplied_draws[[1L]][, 1L])),
    c(1.8, 1.6),
    tolerance = 1e-12
  )

  prior_density <- attr(ratios[[ratio_name]], "prior_density", exact = TRUE)
  expect_s3_class(prior_density, "prior_linear_density")
  expect_equal(attr(prior_density, "support", exact = TRUE), c(0, 2))
  expect_equal(
    BayesTools:::.prior_linear_density_height(prior_density, 1),
    stats::dbeta(0.5, 3, 2) / 2,
    tolerance = 1e-8
  )

  intercept_ratio <- random_effects_summary_posterior(
    fit,
    summary = "var_ratio",
    component = "intercept"
  )
  expect_equal(names(intercept_ratio), "(mu) allocation: var_ratio(intercept)")

  expect_error(
    random_effects_summary_posterior(fit, summary = "var_prop"),
    "Mean-variance SD-component allocations are returned by summary = \"var_ratio\"",
    fixed = TRUE
  )
  expect_s3_class(
    plot_posterior(ratios, ratio_name, prior = TRUE, plot_type = "ggplot"),
    "ggplot"
  )
})

test_that("random-effect summary posterior extracts SD ratios", {

  skip_if_not_installed("runjags")

  fit <- .random_effects_mean_variance_allocation_fit()
  multipliers <- random_effects_summary_posterior(fit, summary = "sd_ratio")
  multiplier_name <- "(mu) allocation: sd_ratio(x)"

  expect_true(multiplier_name %in% names(multipliers))
  expect_equal(
    unname(as.numeric(multipliers[[multiplier_name]])),
    sqrt(c(1.5, 0.5)),
    tolerance = 1e-12
  )

  prior_density <- attr(multipliers[[multiplier_name]], "prior_density", exact = TRUE)
  expect_s3_class(prior_density, "prior_linear_density")
  expect_equal(attr(prior_density, "support", exact = TRUE), c(0, sqrt(2)))
  expect_equal(
    BayesTools:::.prior_linear_density_height(prior_density, 1),
    stats::dbeta(0.5, 3, 2),
    tolerance = 1e-8
  )
})

test_that("full estimates summaries add SD ratios to standard quantities", {

  skip_if_not_installed("runjags")

  fit <- .random_effects_mean_variance_allocation_fit()
  standard <- JAGS_estimates_table(
    fit,
    random_effects_summary = "standard",
    return_samples = TRUE
  )
  full <- JAGS_estimates_table(
    fit,
    random_effects_summary = "full",
    return_samples = TRUE
  )

  expect_setequal(
    colnames(standard),
    c(
      "(mu) allocation: sd_common",
      "(mu) allocation: var_ratio(intercept)",
      "(mu) allocation: var_ratio(x)"
    )
  )
  expect_false(any(grepl(": sd_ratio\\(", colnames(standard))))
  expect_false(any(grepl("(^|: )sd\\(", colnames(standard))))
  expect_false(any(grepl("var_common", colnames(standard), fixed = TRUE)))
  expect_true("(mu) allocation: var_ratio(x)" %in% colnames(full))
  expect_true("(mu) allocation: sd_ratio(x)" %in% colnames(full))
  expect_true("(mu) allocation: var_common" %in% colnames(full))
  expect_true("(mu) sd(x)" %in% colnames(full))
})

test_that("known group covariance retains its fitted kernel scale", {

  skip_if_not_installed("runjags")

  data <- data.frame(study = factor(c("s1", "s1", "s2", "s2")))
  kernel <- matrix(
    c(1, .25, .25, 1),
    nrow = 2,
    dimnames = list(c("s1", "s2"), c("s1", "s2"))
  )
  random_formula <- random_effects_formula(
    ~ 1 | study,
    group_covariance = random_group_covariance(kernel, scale = "none")
  )
  formula_result <- JAGS_formula(
    formula = random_formula,
    parameter = "mu",
    data = data,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      study = random_block(sd = prior("gamma", list(2, 2)))
    )
  )
  random_term <- formula_result$formula_design$random_effects[[1L]]
  samples <- matrix(
    c(1, 2),
    ncol = 1L,
    dimnames = list(NULL, random_term$sd_parameter_names)
  )
  fit <- structure(
    list(mcmc = coda::mcmc.list(coda::mcmc(samples)), sample = nrow(samples)),
    class = c("runjags", "BayesTools_fit", "list")
  )
  attr(fit, "prior_list") <- formula_result$prior_list
  attr(fit, "formula_design") <- list(mu = formula_result$formula_design)
  fit <- attach_test_parameter_map(fit)

  standard <- JAGS_estimates_table(
    fit,
    random_effects_summary = "standard",
    return_samples = TRUE
  )
  full <- JAGS_estimates_table(
    fit,
    random_effects_summary = "full",
    return_samples = TRUE
  )

  expect_identical(colnames(standard), "(mu) sd(intercept)")
  expect_true("(mu) sd(intercept)" %in% colnames(full))
  expect_true("(mu) var(intercept)" %in% colnames(full))

  simplified <- JAGS_estimates_table(
    fit,
    random_effects_summary = "standard",
    simplify_names = TRUE,
    return_samples = TRUE
  )
  expect_identical(colnames(simplified), "(mu) sd")
  expect_equal(simplified[, 1L], samples[, 1L])
})

test_that("random-effect summary posterior extracts total-variance proportions", {

  skip_if_not_installed("runjags")

  fit <- .random_effects_total_variance_allocation_fit()
  proportions <- random_effects_summary_posterior(fit, summary = "var_prop")
  proportion_name <- "(mu) allocation: var_prop(drug)"

  expect_true(proportion_name %in% names(proportions))
  expect_equal(
    unname(as.numeric(proportions[[proportion_name]])),
    c(0.75, 0.25),
    tolerance = 1e-12
  )

  prior_density <- attr(proportions[[proportion_name]], "prior_density", exact = TRUE)
  expect_s3_class(prior_density, "prior_linear_density")
  expect_equal(attr(prior_density, "support", exact = TRUE), c(0, 1))
  expect_equal(
    BayesTools:::.prior_linear_density_height(prior_density, 0.5),
    stats::dbeta(0.5, 3, 2),
    tolerance = 1e-7
  )

  expect_error(
    random_effects_summary_posterior(fit, summary = "var_ratio"),
    "Variance-ratio summaries are created only",
    fixed = TRUE
  )
})

test_that("random-effect summary posterior handles singular Dirichlet boundaries", {

  skip_if_not_installed("runjags")

  fit <- .random_effects_mean_variance_allocation_fit(alpha = c(0.5, 2))
  ratios <- random_effects_summary_posterior(
    fit,
    summary = "var_ratio",
    component = "intercept",
    n_prior_points = 64
  )
  prior_density <- attr(ratios[[1]], "prior_density", exact = TRUE)

  expect_equal(attr(prior_density, "support", exact = TRUE), c(0, 2))
  expect_true(all(is.finite(prior_density$density$x)))
  expect_true(all(is.finite(prior_density$density$y)))
  expect_true(min(prior_density$density$x) > 0)
  expect_identical(
    BayesTools:::.prior_linear_density_height(prior_density, 0),
    Inf
  )
  expect_equal(
    BayesTools:::.prior_linear_density_height(prior_density, 1),
    stats::dbeta(.5, .5, 2) / 2,
    tolerance = 1e-12
  )
  expect_equal(
    attr(prior_density, "singular_boundaries", exact = TRUE),
    0
  )

  prior_plot_data <- BayesTools:::.prior_linear_density_to_plot_data(
    prior_density,
    n_points = 32
  )
  expect_true(all(is.finite(prior_plot_data$density$y)))
})

test_that("scaled-Beta analytic evaluators preserve square-root endpoint limits", {

  prior_density <- BayesTools:::.bt_random_effect_summary_posterior_scaled_beta_density(
    alpha = .5,
    beta = 1,
    scale = 4,
    transform = "sqrt",
    n_grid = 64
  )

  expect_equal(
    BayesTools:::.prior_linear_density_height(prior_density, 0),
    1 / sqrt(4),
    tolerance = 1e-12
  )
  expect_equal(
    BayesTools:::.prior_linear_density_height(prior_density, 2),
    2 * .5 / 2,
    tolerance = 1e-12
  )
  expect_equal(
    BayesTools:::.prior_linear_density_height(prior_density, -1),
    0
  )
})
