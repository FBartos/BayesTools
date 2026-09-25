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

  # mixed columns are the canonical selectors of the level cells
  expect_equal(
    colnames(samples$mu_alloc__xXx__year),
    c(
      "mu_alloc__xXx__year[random]",
      "mu_alloc__xXx__year[systematic]"
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
    prior_density_context = .bt_meta_get(samples, "prior_context"),
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
  fit <- coda::mcmc(complete_test_posterior(posterior, formula_result$prior_list))
  class(fit) <- c("mcmc", "BayesTools_fit")
  attr(fit, "prior_list") <- formula_result$prior_list
  attr(fit, "formula_scale") <- list(mu = formula_result$formula_scale)
  attr(fit, "formula_design") <- list(mu = formula_result$formula_design)
  fit <- attach_test_parameter_map(fit)

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
  fit <- coda::mcmc(complete_test_posterior(posterior, formula_result$prior_list))
  class(fit) <- c("mcmc", "BayesTools_fit")
  attr(fit, "prior_list") <- formula_result$prior_list
  fit <- attach_test_parameter_map(fit)

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
    BayesTools:::.prior_linear_density_point_mass(.bt_meta_get(marginal$alternate, "prior_density"), 0),
    1
  )
  expect_equal(
    BayesTools:::.prior_linear_density_point_mass(.bt_meta_get(marginal$random, "prior_density"), 0),
    0
  )
})


# A fitted object for a formula with a factor interaction, such as `~ g + g:x`,
# `~ g + g:h`, or `~ x + g:x`. Without the main effect of its other component
# (`x` or `h`), the interaction codes `g` by level indicators, so the term has
# one coordinate per level of `g` (including the first); mean-difference and
# orthonormal priors are unavailable for such a term. The main-effect factor
# terms get a factor prior with `contrast` and the interaction (the last term)
# one with `interaction_contrast`. The draws are deterministic normal
# quantiles; interaction coordinate j has mean `slope_means[j]`.
full_rank_interaction_fit <- function(levels, contrast, formula = ~ g + g:x,
                                      slope_means = c(0.3, -0.2, 0.5),
                                      interaction_contrast = contrast) {

  data <- data.frame(
    g = factor(rep(levels, each = 4), levels = levels),
    h = factor(rep(c("u", "v"), 6), levels = c("u", "v")),
    x = seq(-1, 1, length.out = 12)
  )
  factor_prior <- function(term_contrast) {
    distribution <- if (term_contrast %in% c("meandif", "orthonormal")) "mnormal" else "normal"
    prior_factor(distribution, list(0, 1), contrast = term_contrast)
  }
  term_labels <- attr(stats::terms(formula), "term.labels")
  prior_list  <- c(
    list(intercept = prior("normal", list(0, 1))),
    lapply(seq_along(term_labels), function(i) {
      if (identical(term_labels[i], "x")) {
        prior("normal", list(0, 1))
      } else if (i == length(term_labels)) {
        factor_prior(interaction_contrast)
      } else {
        factor_prior(contrast)
      }
    })
  )
  names(prior_list) <- c("intercept", term_labels)
  formula_result <- JAGS_formula(
    formula    = formula,
    parameter  = "mu",
    data       = data,
    prior_list = prior_list
  )

  parameters <- paste0("mu_", gsub(":", "__xXx__", c("intercept", term_labels), fixed = TRUE))
  parameter  <- parameters[length(parameters)]
  draws      <- stats::qnorm(stats::ppoints(200))
  columns    <- list()
  for (term_parameter in parameters) {
    term_prior <- formula_result$prior_list[[term_parameter]]
    n_columns  <- if (is.prior.factor(term_prior)) {
      BayesTools:::.get_prior_factor_levels(term_prior)
    } else {
      1L
    }
    column_draws <- if (identical(term_parameter, "mu_intercept")) {
      matrix(draws, ncol = 1)
    } else if (identical(term_parameter, parameter)) {
      vapply(seq_len(n_columns), function(j) slope_means[j] + 0.3 * draws, numeric(length(draws)))
    } else {
      matrix(rev(draws), nrow = length(draws), ncol = n_columns)
    }
    colnames(column_draws) <- if (is.prior.factor(term_prior) && n_columns > 1L) {
      paste0(term_parameter, "[", seq_len(n_columns), "]")
    } else {
      term_parameter
    }
    columns[[term_parameter]] <- column_draws
  }
  posterior <- do.call(cbind, unname(columns))

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

  list(fit = fit, posterior = posterior, parameter = parameter, parameters = parameters)
}


# Contrasts of the main effect of `g` and of the indicator-coded `g:x` term.
# The term does not use the contrast of `g`, so an independent prior on it
# combines with any contrast of the main effect.
full_rank_interaction_designs <- list(
  c(main = "treatment",   interaction = "treatment"),
  c(main = "independent", interaction = "independent"),
  c(main = "meandif",     interaction = "independent"),
  c(main = "treatment",   interaction = "independent")
)


test_that("full-rank factor-by-continuous interaction slopes are one quantity per level", {

  for (design in full_rank_interaction_designs) {
    for (levels in list(c("a", "b", "c"), c("1", "2", "3"))) {

      info      <- paste0(design[["main"]], " + ", design[["interaction"]], ": ", paste0(levels, collapse = ", "))
      synthetic <- full_rank_interaction_fit(
        levels, design[["main"]], interaction_contrast = design[["interaction"]]
      )
      fit       <- synthetic$fit
      parameter <- synthetic$parameter
      slopes    <- unname(synthetic$posterior[, paste0(parameter, "[", 1:3, "]")])

      # the fitted design gives every level of `g` its own slope coordinate,
      # recorded as the independent (indicator) coding of `g` in the term,
      # while the main effect keeps the contrast of `g`
      expect_equal(
        attr(attr(fit, "prior_list")[[parameter]], "factor_design"),
        diag(3),
        info = info
      )
      expect_identical(
        attr(attr(fit, "prior_list")[[parameter]], "factor_contrasts"),
        c(g = "contr.independent"),
        info = info
      )
      expect_identical(
        attr(attr(fit, "prior_list")$mu_g, "factor_contrasts"),
        c(g = paste0("contr.", design[["main"]])),
        info = info
      )

      samples <- as_mixed_posteriors(fit, synthetic$parameters)
      expect_identical(
        colnames(samples[[parameter]]),
        paste0("mu_g__xXx__x[", levels, "]"),
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
        # column: each treatment or independent slope coordinate has a
        # N(0, 1) prior (one dnorm per coordinate), and both routes estimate
        # the posterior ordinate from the same draws.
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

      main_effect_rows <- switch(
        design[["main"]],
        treatment   = paste0("(mu) g[", levels[-1], "]"),
        independent = paste0("(mu) g[", levels, "]"),
        meandif     = paste0("(mu) g{", 1:2, "}")
      )
      summary_table <- runjags_estimates_table(fit)
      expect_identical(
        rownames(summary_table),
        c(
          "(mu) intercept",
          main_effect_rows,
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
  cells     <- paste0("mu_g__xXx__h[g=", c("a", "b", "c"), ", h=v]")

  # `g` is coded by level indicators and `h` by treatment contrasts: the three
  # coordinates are the cells (a, v), (b, v), and (c, v)
  expect_equal(
    attr(attr(fit, "prior_list")[[parameter]], "factor_design"),
    rbind(matrix(0, 3, 3), diag(3))
  )

  samples <- as_mixed_posteriors(fit, parameter)
  expect_identical(colnames(samples[[parameter]]), cells)

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


test_that("mean-difference and orthonormal priors are unavailable for indicator-coded factor terms", {

  data <- data.frame(
    g = factor(rep(c("a", "b", "c"), each = 4), levels = c("a", "b", "c")),
    h = factor(rep(c("u", "v"), 6), levels = c("u", "v")),
    x = seq(-1, 1, length.out = 12)
  )
  normal <- prior("normal", list(0, 1))
  message_for <- function(contrast, term, factor, missing_term) {
    paste0(
      "The '", contrast, "' prior of the factor term '", term,
      "' is unavailable: the formula has no term '", missing_term, "', so '",
      term, "' codes '", factor, "' by level indicators and has one ",
      "coefficient per level instead of '", contrast, "' contrast ",
      "coefficients. Add '", missing_term, "' to the formula ",
      "to keep the '", contrast, "' contrast, or use ",
      "prior_factor(contrast = \"independent\") for one independent ",
      "coefficient per level."
    )
  }

  for (contrast in c("meandif", "orthonormal")) {
    factor_prior <- prior_factor("mnormal", list(0, 1), contrast = contrast)

    expect_error(
      JAGS_formula(~ g + g:x, "mu", data, list(
        intercept = normal, g = factor_prior, "g:x" = factor_prior
      )),
      message_for(contrast, "g:x", "g", "x"),
      fixed = TRUE
    )
    # default factor priors (as RoBMA's model-averaging constructors supply)
    expect_error(
      JAGS_formula(~ g + g:x, "mu", data, list(
        intercept = normal, "__default_factor" = factor_prior
      )),
      message_for(contrast, "g:x", "g", "x"),
      fixed = TRUE
    )
    expect_error(
      JAGS_formula(~ g + g:x, "mu", data, list(
        intercept = normal, g = factor_prior,
        "g:x" = prior_spike_and_slab(factor_prior)
      )),
      message_for(contrast, "g:x", "g", "x"),
      fixed = TRUE
    )
    expect_error(
      JAGS_formula(~ g:x, "mu", data, list(
        intercept = normal, "g:x" = factor_prior
      )),
      message_for(contrast, "g:x", "g", "x"),
      fixed = TRUE
    )
    expect_error(
      JAGS_formula(~ g + g:h, "mu", data, list(
        intercept = normal, g = factor_prior, "g:h" = factor_prior
      )),
      message_for(contrast, "g:h", "g", "h"),
      fixed = TRUE
    )

    # with the lower-order term, the interaction keeps the contrast basis
    reduced <- JAGS_formula(~ g * x, "mu", data, list(
      intercept = normal, g = factor_prior, x = normal, "g:x" = factor_prior
    ))
    reduced_prior <- reduced$prior_list$mu_g__xXx__x
    contrast_matrix <- if (contrast == "meandif") {
      contr.meandif(c("a", "b", "c"))
    } else {
      contr.orthonormal(c("a", "b", "c"))
    }
    expect_equal(BayesTools:::.get_prior_factor_levels(reduced_prior), 2)
    expect_equal(attr(reduced_prior, "factor_design"), unname(contrast_matrix), tolerance = 1e-12)
    expect_equal(
      unname(reduced$data$mu_data_g__xXx__x),
      unname(contrast_matrix[as.integer(data$g), ] * data$x),
      tolerance = 1e-12
    )
    # so does a factor that enters only an interaction with a present margin
    expect_equal(
      BayesTools:::.get_prior_factor_levels(
        JAGS_formula(~ x + g:x, "mu", data, list(
          intercept = normal, x = normal, "x:g" = factor_prior
        ))$prior_list$mu_x__xXx__g
      ),
      2
    )
  }

  # treatment and independent priors keep one coefficient per level
  for (contrast in c("treatment", "independent")) {
    factor_prior <- prior_factor("normal", list(0, 1), contrast = contrast)
    full_rank <- JAGS_formula(~ g + g:x, "mu", data, list(
      intercept = normal, g = factor_prior, "g:x" = factor_prior
    ))
    expect_equal(BayesTools:::.get_prior_factor_levels(full_rank$prior_list$mu_g__xXx__x), 3, info = contrast)
    expect_equal(attr(full_rank$prior_list$mu_g__xXx__x, "factor_design"), diag(3), info = contrast)
  }

  # an indicator-coded term does not use the contrast of `g`, so its
  # independent prior combines with a mean-difference main effect; a
  # reduced-rank interaction uses that contrast, so conflicting priors on it
  # and on the main effect are still rejected
  meandif_prior     <- prior_factor("mnormal", list(0, 1), contrast = "meandif")
  independent_prior <- prior_factor("normal", list(0, 1), contrast = "independent")
  mixed <- JAGS_formula(~ g + g:x, "mu", data, list(
    intercept = normal, g = meandif_prior, "g:x" = independent_prior
  ))
  expect_identical(attr(mixed$prior_list$mu_g, "factor_contrasts"), c(g = "contr.meandif"))
  expect_identical(attr(mixed$prior_list$mu_g__xXx__x, "factor_contrasts"), c(g = "contr.independent"))
  expect_equal(
    unname(mixed$data$mu_data_g),
    unname(contr.meandif(c("a", "b", "c"))[as.integer(data$g), ]),
    tolerance = 1e-12
  )
  expect_equal(
    unname(mixed$data$mu_data_g__xXx__x),
    unname(diag(3)[as.integer(data$g), ] * data$x),
    tolerance = 1e-12
  )
  expect_error(
    JAGS_formula(~ g * x, "mu", data, list(
      intercept = normal, g = meandif_prior, x = normal, "g:x" = independent_prior
    )),
    "Factor predictor 'g' has conflicting contrast priors across formula terms.",
    fixed = TRUE
  )

  # a mean-difference point mass at zero (the null-hypothesis prior of
  # RoBMA's model-averaged terms) fixes every coefficient at zero in any basis
  null_term <- JAGS_formula(~ g + g:x, "mu", data, list(
    intercept = normal, g = meandif_prior,
    "g:x" = prior_factor("spike", list(0), contrast = "meandif")
  ))
  expect_equal(attr(null_term$prior_list$mu_g__xXx__x, "factor_design"), diag(3))
  expect_identical(attr(null_term$prior_list$mu_g__xXx__x, "factor_contrasts"), c(g = "contr.independent"))
})


test_that("formula marginal posteriors accept predictors that enter only interactions", {

  for (design in full_rank_interaction_designs) {

    contrast  <- paste(design, collapse = " + ")
    synthetic <- full_rank_interaction_fit(
      c("a", "b", "c"), design[["main"]], interaction_contrast = design[["interaction"]]
    )
    fit       <- synthetic$fit
    posterior <- synthetic$posterior
    parameter <- synthetic$parameter
    samples   <- as_mixed_posteriors(fit, synthetic$parameters)

    # independent reference: the linear predictor from the fitted term designs
    g_columns     <- grep("^mu_g\\[", colnames(posterior), value = TRUE)
    level_effects <- posterior[, g_columns, drop = FALSE] %*%
      t(attr(attr(fit, "prior_list")$mu_g, "factor_design"))
    slopes <- posterior[, paste0(parameter, "[", 1:3, "]")]

    # `x` has no main-effect term: it is continuous, varied over -1, 0, 1
    # for the interaction and held at 0 for the main effect of `g`
    marginal <- marginal_posterior(
      samples       = samples,
      parameter     = parameter,
      formula       = ~ g + g:x,
      prior_samples = TRUE,
      n_samples     = 100
    )
    expect_identical(
      names(marginal),
      paste0(rep(c("a", "b", "c"), 3), ", ", rep(c("-1", "0", "1"), each = 3), "SD"),
      info = contrast
    )
    for (level_i in 1:3) {
      for (x_i in 1:3) {
        x <- c(-1, 0, 1)[x_i]
        expect_equal(
          as.numeric(marginal[[(x_i - 1) * 3 + level_i]]),
          unname(posterior[, "mu_intercept"] + level_effects[, level_i] + slopes[, level_i] * x),
          tolerance = 1e-12,
          info = paste(contrast, level_i, x)
        )
      }
    }

    main_effect <- marginal_posterior(
      samples       = samples,
      parameter     = "mu_g",
      formula       = ~ g + g:x,
      prior_samples = TRUE,
      n_samples     = 100
    )
    expect_identical(names(main_effect), c("a", "b", "c"), info = contrast)
    for (level_i in 1:3) {
      expect_equal(
        as.numeric(main_effect[[level_i]]),
        unname(posterior[, "mu_intercept"] + level_effects[, level_i]),
        tolerance = 1e-12,
        info = paste(contrast, level_i)
      )
    }
  }

  # a mean-difference factor `g` without a main-effect term takes its fitted
  # levels and contrast from `x:g` (coded by the contrast, since `x` is in
  # the formula); it is omitted (coded as zero) for the main effect of `x`
  synthetic <- full_rank_interaction_fit(c("a", "b", "c"), "meandif", ~ x + g:x)
  posterior <- synthetic$posterior
  parameter <- synthetic$parameter
  samples   <- as_mixed_posteriors(synthetic$fit, synthetic$parameters)
  design    <- attr(attr(synthetic$fit, "prior_list")[[parameter]], "factor_design")
  expect_equal(design, unname(contr.meandif(c("a", "b", "c"))), tolerance = 1e-12)
  level_slopes <- posterior[, paste0(parameter, "[", 1:2, "]")] %*% t(design)
  slope_cells <- marginal_posterior(
    samples       = samples,
    parameter     = parameter,
    formula       = ~ x + g:x,
    prior_samples = TRUE,
    n_samples     = 100
  )
  expect_identical(
    names(slope_cells),
    paste0(rep(c("-1", "0", "1"), 3), "SD, ", rep(c("a", "b", "c"), each = 3))
  )
  expect_equal(
    as.numeric(slope_cells[["1SD, c"]]),
    unname(posterior[, "mu_intercept"] + posterior[, "mu_x"] + level_slopes[, 3]),
    tolerance = 1e-12
  )
  main_effect <- marginal_posterior(
    samples       = samples,
    parameter     = "mu_x",
    formula       = ~ x + g:x,
    prior_samples = TRUE,
    n_samples     = 100
  )
  expect_identical(names(main_effect), c("-1SD", "0SD", "1SD"))
  expect_equal(
    as.numeric(main_effect[["1SD"]]),
    unname(posterior[, "mu_intercept"] + posterior[, "mu_x"]),
    tolerance = 1e-12
  )

  # a treatment factor `h` without a main-effect term takes its fitted levels
  # and contrast from `g:h`: it is held at its reference level `u`
  synthetic <- full_rank_interaction_fit(c("a", "b", "c"), "treatment", ~ g + g:h)
  posterior <- synthetic$posterior
  samples   <- as_mixed_posteriors(synthetic$fit, synthetic$parameters)
  cells <- marginal_posterior(
    samples       = samples,
    parameter     = synthetic$parameter,
    formula       = ~ g + g:h,
    prior_samples = TRUE,
    n_samples     = 100
  )
  expect_identical(names(cells), c("a, u", "b, u", "c, u", "a, v", "b, v", "c, v"))
  expect_equal(
    as.numeric(cells[["c, v"]]),
    unname(posterior[, "mu_intercept"] + posterior[, "mu_g[2]"] + posterior[, "mu_g__xXx__h[3]"]),
    tolerance = 1e-12
  )
  main_effect <- marginal_posterior(
    samples       = samples,
    parameter     = "mu_g",
    formula       = ~ g + g:h,
    prior_samples = TRUE,
    n_samples     = 100
  )
  expect_identical(
    as.character(attr(main_effect, "data")$h),
    rep("u", 3)
  )
  expect_equal(
    as.numeric(main_effect[["b"]]),
    unname(posterior[, "mu_intercept"] + posterior[, "mu_g[1]"]),
    tolerance = 1e-12
  )
})


test_that("JAGS_evaluate_formula codes indicator-coded terms as fitted, also for draws without a fit", {

  # In `~ x + g:z + x:g`, `g:z` precedes `x:g` and codes `g` by level
  # indicators, so it records the independent coding of `g`; `x:g` codes `g`
  # by its fitted contrast. Evaluation through the fitted design (attached to
  # the draws, or built by JAGS_formula_draws() for draws without a fit)
  # reproduces the fitted JAGS data columns.
  data <- data.frame(
    x = c(-2, -1, 1, 2, 3, 4, 0.5, -0.5, 1.5),
    z = c(0.3, -1, 2, 1, -0.5, 0.7, 1.1, -0.2, 0.4),
    g = factor(rep(c("a", "b", "c"), 3), levels = c("a", "b", "c"))
  )
  normal <- prior("normal", list(0, 1))
  term_priors <- list(
    treatment = list(
      "g:z" = prior_factor("normal", list(0, 1), contrast = "treatment"),
      "x:g" = prior_factor("normal", list(0, 1), contrast = "treatment")
    ),
    meandif = list(
      "g:z" = prior_factor("normal", list(0, 1), contrast = "independent"),
      "x:g" = prior_factor("mnormal", list(0, 1), contrast = "meandif")
    )
  )

  for (contrast in names(term_priors)) {
    formula_result <- JAGS_formula(
      ~ x + g:z + x:g, "mu", data,
      c(list(intercept = normal, x = normal), term_priors[[contrast]])
    )
    prior_list <- formula_result$prior_list
    expect_identical(
      attr(prior_list$mu_g__xXx__z, "factor_contrasts"),
      c(g = "contr.independent"),
      info = contrast
    )
    expect_identical(
      attr(prior_list$mu_x__xXx__g, "factor_contrasts"),
      c(g = paste0("contr.", contrast)),
      info = contrast
    )

    # reference: the fitted JAGS data columns times the coefficients
    coefficients <- list(
      mu_intercept = 0.5,
      mu_x         = -0.25,
      mu_g__xXx__z = c(0.3, -0.6, 0.9),
      mu_x__xXx__g = c(0.4, -0.8)
    )
    expected <- coefficients$mu_intercept +
      coefficients$mu_x * formula_result$data$mu_data_x +
      drop(formula_result$data$mu_data_g__xXx__z %*% coefficients$mu_g__xXx__z) +
      drop(formula_result$data$mu_data_x__xXx__g %*% coefficients$mu_x__xXx__g)
    posterior <- coda::mcmc(matrix(
      unlist(coefficients, use.names = FALSE),
      nrow = 1,
      dimnames = list(NULL, c(
        "mu_intercept", "mu_x",
        paste0("mu_g__xXx__z[", 1:3, "]"),
        paste0("mu_x__xXx__g[", 1:2, "]")
      ))
    ))

    attr(posterior, "formula_design") <- list(mu = formula_result$formula_design)
    prediction <- JAGS_evaluate_formula(posterior, ~ x + g:z + x:g, "mu", data, prior_list)
    expect_equal(unname(drop(prediction)), unname(expected), tolerance = 1e-12, info = contrast)

    # draws without a fit evaluate through the design JAGS_formula_draws() builds
    formula_draws <- JAGS_formula_draws(
      as.matrix(posterior), ~ x + g:z + x:g, "mu", data,
      c(list(intercept = normal, x = normal), term_priors[[contrast]])
    )
    draws_prediction <- JAGS_evaluate_formula(formula_draws, ~ x + g:z + x:g, "mu", data, prior_list)
    expect_equal(unname(drop(draws_prediction)), unname(expected), tolerance = 1e-12, info = contrast)
  }
})
