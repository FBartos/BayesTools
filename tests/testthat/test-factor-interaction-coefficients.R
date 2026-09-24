skip_if_not_test_profile("unit")

test_that("as_mixed_posteriors handles treatment factor-continuous interaction coefficients", {

  df <- data.frame(
    alloc = factor(
      rep(c("alternate", "random", "systematic"), each = 4),
      levels = c("alternate", "random", "systematic")
    ),
    year = seq(1960, 1971, length.out = 12)
  )
  formula_result <- JAGS_formula(
    formula       = ~ alloc * year,
    parameter     = "mu",
    data          = df,
    formula_scale = list(year = TRUE),
    prior_list    = list(
      intercept    = prior("normal", list(0, 1)),
      alloc        = prior_factor("normal", list(0, 1), contrast = "treatment"),
      year         = prior("normal", list(0, 1)),
      "alloc:year" = prior_factor("normal", list(0, 1), contrast = "treatment")
    )
  )
  interaction_prior <- formula_result$prior_list$mu_alloc__xXx__year

  posterior <- matrix(seq_len(20), nrow = 10, ncol = 2)
  colnames(posterior) <- paste0("mu_alloc__xXx__year[", 1:2, "]")
  fit <- coda::mcmc(posterior)
  class(fit) <- c("mcmc", "BayesTools_fit")
  attr(fit, "prior_list") <- formula_result$prior_list[
    "mu_alloc__xXx__year"
  ]
  attr(fit, "formula_scale") <- list(mu = formula_result$formula_scale)
  fit <- attach_test_parameter_map(fit)

  samples <- as_mixed_posteriors(
    model            = fit,
    parameters       = "mu_alloc__xXx__year",
    transform_scaled = TRUE,
    n_prior_samples  = 64
  )

  expect_equal(
    colnames(samples$mu_alloc__xXx__year),
    c(
      "mu_alloc[random]__xXx__year",
      "mu_alloc[systematic]__xXx__year"
    )
  )
  expect_equal(attr(samples$mu_alloc__xXx__year, "factor_design"), attr(interaction_prior, "factor_design"))

  plot_data <- BayesTools:::.plot_data_samples.factor(
    samples                  = samples,
    parameter                = "mu_alloc__xXx__year",
    n_points                 = 32,
    transformation           = NULL,
    transformation_arguments = NULL,
    transformation_settings  = FALSE
  )
  density_entries <- plot_data[vapply(plot_data, inherits, logical(1), what = "density.prior.factor")]

  expect_equal(length(density_entries), 2L)
  expect_equal(
    unname(vapply(density_entries, attr, character(1), which = "level_name")),
    colnames(samples$mu_alloc__xXx__year)
  )

  prior_plot_data <- BayesTools:::.plot_data_prior_factor_density_transformed(
    prior_density_context = attr(samples, "prior_density_context"),
    samples               = samples,
    parameter             = "mu_alloc__xXx__year",
    prior_list            = attr(samples$mu_alloc__xXx__year, "prior_list"),
    n_points              = 32
  )

  expect_equal(length(prior_plot_data), 2L)
  expect_equal(
    unname(vapply(prior_plot_data, attr, character(1), which = "level_name")),
    colnames(samples$mu_alloc__xXx__year)
  )
})


test_that("marginal estimates unscale fixed-factor interaction summaries", {

  df <- data.frame(
    alloc = factor(
      rep(c("alternate", "random", "systematic"), each = 4),
      levels = c("alternate", "random", "systematic")
    ),
    year = seq(10, 32, length.out = 12)
  )
  formula_result <- JAGS_formula(
    formula       = ~ alloc * year,
    parameter     = "mu",
    data          = df,
    formula_scale = list(year = TRUE),
    prior_list    = list(
      intercept    = prior("normal", list(0, 1)),
      alloc        = prior_factor("normal", list(0, 1), contrast = "treatment"),
      year         = prior("normal", list(0, 1)),
      "alloc:year" = prior_factor("normal", list(0, 1), contrast = "treatment")
    )
  )

  parameter <- "mu_alloc__xXx__year"
  interaction_prior <- formula_result$prior_list[[parameter]]
  posterior <- matrix(seq_len(12), nrow = 6, ncol = 2)
  colnames(posterior) <- paste0(parameter, "[", 1:2, "]")
  fit <- coda::mcmc(posterior)
  class(fit) <- c("mcmc", "BayesTools_fit")
  attr(fit, "prior_list") <- formula_result$prior_list
  attr(fit, "formula_scale") <- list(mu = formula_result$formula_scale)
  attr(fit, "formula_design") <- list(mu = formula_result$formula_design)

  scaled_samples <- as_mixed_posteriors(
    model = fit,
    parameters = parameter
  )
  marginal <- marginal_posterior(
    samples = scaled_samples,
    parameter = parameter,
    use_formula = FALSE
  )
  samples <- stats::setNames(list(marginal), parameter)
  inference <- stats::setNames(
    list(stats::setNames(as.list(rep(1, length(marginal))), names(marginal))),
    parameter
  )

  table <- marginal_estimates_table(
    samples = samples,
    inference = inference,
    parameters = parameter,
    probs = c(0.25, 0.75),
    transform_scaled = TRUE,
    formula_scale = attr(fit, "formula_scale")
  )

  expected_samples <- posterior %*% t(attr(interaction_prior, "factor_design"))
  expected_samples <- expected_samples / formula_result$formula_scale$mu_year$sd
  expected_summary <- cbind(
    Mean = colMeans(expected_samples),
    SD = apply(expected_samples, 2, stats::sd),
    "0.25" = apply(expected_samples, 2, stats::quantile, probs = 0.25),
    "0.75" = apply(expected_samples, 2, stats::quantile, probs = 0.75)
  )

  expect_equal(
    unname(as.matrix(table[, colnames(expected_summary)])),
    unname(expected_summary),
    tolerance = 1e-12
  )
})


test_that("marginal_posterior handles treatment factor-continuous interaction coefficients", {

  df <- data.frame(
    alloc = factor(
      rep(c("alternate", "random", "systematic"), each = 4),
      levels = c("alternate", "random", "systematic")
    ),
    year = seq(1960, 1971, length.out = 12)
  )
  formula_result <- JAGS_formula(
    formula    = ~ alloc * year,
    parameter  = "mu",
    data       = df,
    prior_list = list(
      intercept    = prior("normal", list(0, 1)),
      alloc        = prior_factor("normal", list(0, 1), contrast = "treatment"),
      year         = prior("normal", list(0, 1)),
      "alloc:year" = prior_factor("normal", list(0, 1), contrast = "treatment")
    )
  )
  interaction_prior <- formula_result$prior_list$mu_alloc__xXx__year

  posterior <- matrix(seq_len(20), nrow = 10, ncol = 2)
  colnames(posterior) <- paste0("mu_alloc__xXx__year[", 1:2, "]")
  fit <- coda::mcmc(posterior)
  class(fit) <- c("mcmc", "BayesTools_fit")
  attr(fit, "prior_list") <- formula_result$prior_list

  samples <- as_mixed_posteriors(
    model      = fit,
    parameters = "mu_alloc__xXx__year"
  )
  marginal <- marginal_posterior(
    samples       = samples,
    parameter     = "mu_alloc__xXx__year",
    prior_samples = TRUE,
    use_formula   = FALSE,
    n_samples     = 32
  )
  expected <- posterior %*% t(attr(interaction_prior, "factor_design"))

  expect_equal(names(marginal), c("alternate", "random", "systematic"))
  expect_equal(as.numeric(marginal$alternate), as.numeric(expected[, 1]))
  expect_equal(as.numeric(marginal$random), as.numeric(expected[, 2]))
  expect_equal(as.numeric(marginal$systematic), as.numeric(expected[, 3]))
  expect_equal(
    BayesTools:::.prior_linear_density_point_mass(attr(marginal$alternate, "prior_density"), 0),
    1
  )
  expect_equal(
    BayesTools:::.prior_linear_density_point_mass(attr(marginal$random, "prior_density"), 0),
    0
  )
})


# A fitted object for a formula whose factor interaction has no main effect of
# its other component, `~ g + g:x` or `~ g + g:h`: the interaction design then
# codes `g` by level indicators, so the term has one coordinate per level of
# `g` (including the first) whatever the contrast. The draws are deterministic
# normal quantiles; interaction coordinate j has mean `slope_means[j]`.
full_rank_interaction_fit <- function(levels, contrast, formula = ~ g + g:x,
                                      slope_means = c(0.3, -0.2, 0.5)) {

  data <- data.frame(
    g = factor(rep(levels, each = 4), levels = levels),
    h = factor(rep(c("u", "v"), 6), levels = c("u", "v")),
    x = seq(-1, 1, length.out = 12)
  )
  distribution <- if (contrast == "meandif") "mnormal" else "normal"
  term <- attr(stats::terms(formula), "term.labels")[2]
  prior_list <- list(
    intercept = prior("normal", list(0, 1)),
    g         = prior_factor(distribution, list(0, 1), contrast = contrast)
  )
  prior_list[[term]] <- prior_factor(distribution, list(0, 1), contrast = contrast)
  formula_result <- JAGS_formula(
    formula    = formula,
    parameter  = "mu",
    data       = data,
    prior_list = prior_list
  )

  parameter <- paste0("mu_", gsub(":", "__xXx__", term, fixed = TRUE))
  n_g       <- BayesTools:::.get_prior_factor_levels(formula_result$prior_list$mu_g)
  n_term    <- BayesTools:::.get_prior_factor_levels(formula_result$prior_list[[parameter]])
  draws     <- stats::qnorm(stats::ppoints(200))
  posterior <- cbind(
    draws,
    matrix(rev(draws), nrow = length(draws), ncol = n_g),
    vapply(seq_len(n_term), function(j) slope_means[j] + 0.3 * draws, numeric(length(draws)))
  )
  colnames(posterior) <- c(
    "mu_intercept",
    paste0("mu_g[", seq_len(n_g), "]"),
    paste0(parameter, "[", seq_len(n_term), "]")
  )

  fit <- structure(
    list(
      mcmc         = coda::mcmc.list(coda::mcmc(posterior)),
      sample       = nrow(posterior),
      summary.pars = list(mutate = NULL),
      monitor      = colnames(posterior)
    ),
    class = c("runjags", "BayesTools_fit", "list")
  )
  attr(fit, "prior_list")     <- formula_result$prior_list
  attr(fit, "formula_design") <- list(mu = formula_result$formula_design)
  fit <- attach_test_parameter_map(fit)

  list(fit = fit, posterior = posterior, parameter = parameter)
}


test_that("full-rank factor-by-continuous interaction slopes are one quantity per level", {

  for (contrast in c("treatment", "meandif")) {
    for (levels in list(c("a", "b", "c"), c("1", "2", "3"))) {

      info      <- paste0(contrast, ": ", paste0(levels, collapse = ", "))
      synthetic <- full_rank_interaction_fit(levels, contrast)
      fit       <- synthetic$fit
      parameter <- synthetic$parameter
      slopes    <- unname(synthetic$posterior[, paste0(parameter, "[", 1:3, "]")])

      # the fitted design gives every level of `g` its own slope coordinate
      expect_equal(
        attr(attr(fit, "prior_list")[[parameter]], "factor_design"),
        diag(3),
        info = info
      )

      samples <- as_mixed_posteriors(fit, c("mu_intercept", "mu_g", parameter))
      expect_identical(
        colnames(samples[[parameter]]),
        if (contrast == "treatment") {
          paste0("mu_g[", levels, "]__xXx__x")
        } else {
          paste0(parameter, "{", 1:3, "}")
        },
        info = info
      )

      marginal <- marginal_posterior(
        samples       = samples,
        parameter     = parameter,
        use_formula   = FALSE,
        prior_samples = TRUE
      )
      expect_identical(names(marginal), levels, info = info)

      catalog <- parameter_catalog(fit)
      for (i in seq_along(levels)) {

        level      <- levels[i]
        level_info <- paste0(info, "; level ", level)

        # independent reference: the level's slope is its JAGS column
        selection <- parameter_catalog_resolve(
          catalog, paste0("g:x[", level, "]"), namespace = "mu"
        )
        expect_identical(
          as.numeric(as.matrix(parameter_draws(fit, selection)[[1]])),
          slopes[, i],
          info = level_info
        )
        expect_identical(as.numeric(marginal[[level]]), slopes[, i], info = level_info)

        # The level's point hypothesis is the scalar Savage-Dickey test of its
        # column: each slope coordinate has a N(0, 1) prior (dnorm for a
        # treatment coordinate, an identity-precision dmnorm for the
        # mean-difference coefficients), and both routes estimate the
        # posterior ordinate from the same draws.
        level_test <- hypothesis_BF(
          marginal,
          hypothesis = paste0("`", parameter, "[", level, "]` = 0"),
          columns    = "all"
        )
        column_test <- hypothesis_BF(
          slopes[, i],
          prior      = prior("normal", list(0, 1)),
          hypothesis = "theta = 0",
          parameter  = "theta",
          columns    = "all"
        )
        expect_equal(level_test$prior, stats::dnorm(0), tolerance = 1e-12, info = level_info)
        expect_equal(level_test$posterior, column_test$posterior, tolerance = 1e-10, info = level_info)
        expect_equal(level_test$BF, column_test$BF, tolerance = 1e-10, info = level_info)
      }

      if (contrast == "treatment") {
        summary_table <- runjags_estimates_table(fit)
        expect_identical(
          rownames(summary_table),
          c(
            "(mu) intercept",
            paste0("(mu) g[", levels[-1], "]"),
            paste0("(mu) g[", levels, "]:x")
          ),
          info = info
        )
        expect_equal(
          unname(summary_table[paste0("(mu) g[", levels, "]:x"), "Mean"]),
          colMeans(slopes),
          tolerance = 1e-12,
          info = info
        )
      }

      model_list <- list(
        list(fit = fit, marglik = bridgesampling_object(0), prior_weights = 1),
        list(fit = fit, marglik = bridgesampling_object(0), prior_weights = 1)
      )
      mixed <- mix_posteriors(
        model_list   = model_list,
        parameters   = parameter,
        is_null_list = stats::setNames(list(c(FALSE, FALSE)), parameter),
        seed         = 1,
        n_samples    = 50
      )
      expect_identical(colnames(mixed[[parameter]]), colnames(samples[[parameter]]), info = info)
    }
  }
})


test_that("partially full-rank factor interactions are named by their level cells", {

  synthetic <- full_rank_interaction_fit(c("a", "b", "c"), "treatment", ~ g + g:h)
  fit       <- synthetic$fit
  parameter <- synthetic$parameter
  cells     <- paste0("mu_g[", c("a", "b", "c"), "]__xXx__h[v]")

  # `g` is coded by level indicators and `h` by treatment contrasts: the three
  # coordinates are the cells (a, v), (b, v), and (c, v)
  expect_equal(
    attr(attr(fit, "prior_list")[[parameter]], "factor_design"),
    rbind(matrix(0, 3, 3), diag(3))
  )

  samples <- as_mixed_posteriors(fit, parameter)
  expect_identical(colnames(samples[[parameter]]), cells)

  renamed <- BayesTools:::.rename_factor_levels(
    synthetic$posterior,
    attr(fit, "prior_list")
  )
  expect_identical(
    colnames(renamed),
    c("mu_intercept", "mu_g[b]", "mu_g[c]", cells)
  )
  expect_identical(
    rownames(runjags_estimates_table(fit)),
    c("(mu) intercept", "(mu) g[b]", "(mu) g[c]", paste0("(mu) g[", c("a", "b", "c"), "]:h[v]"))
  )

  selection <- parameter_catalog_resolve(
    parameter_catalog(fit), "mu_g__xXx__h[g=a, h=v]"
  )
  expect_identical(
    as.numeric(as.matrix(parameter_draws(fit, selection)[[1]])),
    unname(synthetic$posterior[, paste0(parameter, "[1]")])
  )
})
