skip_if_not_test_profile("unit")

test_that("R2 convergence and block declarations distinguish values, gates and allocation copies", {

  data <- data.frame(x = seq_len(8), f = factor(rep(letters[1:4], 2)))
  probability <- prior("spike", list(0.5))
  total <- prior_spike_and_slab(prior("beta", list(2, 3)), prior_inclusion = prior("beta", list(2, 2)))
  compile <- function(formula, budget) {
    JAGS_formula(formula, "mu", data, list(intercept = prior("normal", list(0, 1)), Signal = budget))
  }
  info <- compile(~ x + f, prior_R2(
    total = total, inclusion = prior("beta", list(2, 2)), model = R2_model(contrast = "treatment", sigma = 1)
  ))
  priors <- info$prior_list
  core <- priors$mu_Signal
  nodes <- .bt_R2_core_nodes(core)
  roles <- .bt_convergence_role_prior("mu_Signal", core)
  role <- function(coordinates) unname(roles$role[match(coordinates, roles$name)])
  expect_identical(role(nodes$total), "sampled")
  expect_true(all(role(c(nodes$gates, nodes$indicator)) == "indicator"))
  expect_true(all(role(c(nodes$open, nodes$allocation)) == "derived"))
  expect_true(all(role(c(nodes$rate, nodes$rate_u, nodes$inclusion, nodes$variable, nodes$eta)) == "auxiliary"))

  set.seed(785)
  chains <- lapply(seq_len(2L), function(i) coda::mcmc(.generate_prior_sample_matrix(priors, 100L)))
  fit <- .parameter_catalog_test_fit(coda::mcmc.list(chains), priors, formula_design = list(mu = info$formula_design))
  requested <- function(indicators) {
    diagnostic <- attr(JAGS_check_convergence(fit, check_indicators = indicators), "diagnostics")
    diagnostic$parameter[diagnostic$state != "not_requested"]
  }
  default <- requested(FALSE)
  expect_true(all(c("mu_x", nodes$total) %in% default))
  expect_false(any(c(nodes$gates, nodes$indicator, nodes$open, nodes$allocation, nodes$rate) %in% default))
  expect_setequal(setdiff(requested(TRUE), default), c(nodes$gates, nodes$indicator))

  share <- parameter_catalog_resolve(parameter_catalog(fit), "Signal[share: x]")
  provider <- hypothesis_provider(fit)
  entries <- provider$metadata()$entries
  root <- entries$parameter[match(share$quantity_id, entries$quantity_id)]
  samples <- provider$random$mixed_posterior(root, TRUE, FALSE, NULL)
  expect_named(samples, root)
  expect_named(attr(samples, "prior_list"), root)
  columns <- posterior_metadata(samples[[root]], "quantities")
  expect_identical(columns$column, root)
  expect_identical(columns$quantity_id, share$quantity_id)

  # The sampled gate and its realized open flag differ on draws in the spike.
  expect_true(any(chains[[1]][, nodes$gates[1]] != chains[[1]][, nodes$open[1]]))
  for (target in c("mu_x", "mu_Signal")) {
    event <- if (target == "mu_x") nodes$open[1] else nodes$indicator
    for (statistic in c("inclusion", "component", "switching")) {
      trace <- JAGS_diagnostics_trace(fit, target, plot_type = "ggplot", statistic = statistic, batch_size = 50)
      layers <- ggplot2::ggplot_build(trace)$data
      for (i in seq_along(chains)) {
        batches <- split(as.numeric(chains[[i]][, event]), rep(1:2, each = 50L))
        expected <- vapply(batches, function(z) {
          if (statistic == "switching") mean(diff(z) != 0) else mean(z)
        }, numeric(1))
        expect_equal(layers[[i]]$y, unname(expected))
      }
    }
  }

  # A truncated, one-term slab always opens its raw gate. The model's spike
  # can still turn off the realized inclusion.
  one <- compile(~ x, prior_R2(total = total, inclusion = probability, model = R2_model(sigma = 1)))
  roles <- .bt_convergence_role_prior("mu_Signal", one$prior_list$mu_Signal)
  single <- .bt_R2_core_nodes(one$prior_list$mu_Signal)
  expect_identical(role(single$gates), "structural")
  expect_identical(role(single$indicator), "indicator")
  expect_identical(role(single$open), "derived")

  ordinary <- compile(~ x + f, prior_R2(model = R2_model(sigma = 1)))
  allocation <- .bt_convergence_role_prior("mu_Signal_allocation", ordinary$prior_list$mu_Signal_allocation)
  expect_identical(allocation$role, "derived")
  targets <- c("mu_x", "mu_f", "mu_Signal")
  sources <- JAGS_sampler_sources(priors, targets)
  expect_identical(sources$prior, targets)
  expect_identical(sources$node, rep(NA_character_, 3L))
  expect_identical(sources$status, rep("declined", 3L))
  expect_identical(sources$reason, c(rep("its bound R-squared coefficient is deterministic", 2L),
                                    "its realized R-squared budget total is deterministic"))
  for (i in seq_along(targets)) {
    expect_error(.bt_sampler_block_declarations(list(targets[[i]]), priors),
                 sources$reason[[i]], fixed = TRUE,
                 class = "BayesTools_sampler_declaration_error")
  }
  ordinary_sources <- JAGS_sampler_sources(ordinary$prior_list, c("mu_x", "mu_Signal"))
  expect_identical(ordinary_sources$node, c(NA_character_, "mu_Signal"))
  expect_identical(ordinary_sources$status, c("declined", "eligible"))
})

test_that("an unavailable kernel factor projection only refuses hypotheses about that level", {

  set.seed(784)
  data <- data.frame(f = factor(rep(letters[1:3], each = 4)))
  info <- JAGS_formula(~ f, "mu", data, list(
    intercept = prior("normal", list(0, 1)),
    Signal = prior_R2(model = R2_model(law = "normal", sigma = 1))
  ))
  samples <- .generate_prior_sample_matrix(info$prior_list, 100L)
  fit <- .parameter_catalog_test_fit(coda::mcmc.list(coda::mcmc(samples)), info$prior_list,
                                     formula_design = list(mu = info$formula_design))
  # Limit the listing to the unsupported projection and an unaffected total.
  # The supported coefficient marginal has its own numerical R2 tests.
  provider <- hypothesis_provider(fit)
  metadata <- provider$metadata()
  expect_identical(provider$random$select("mu_f[a]", FALSE)$spec$quantity, "coefficient")
  expect_identical(provider$random$select("mu_Signal", FALSE)$spec$quantity, "total")
  expect_identical(metadata$entries$source[metadata$entries$parameter == "mu_Signal"], "BayesTools::r2_budget")
  metadata$entries <- metadata$entries[
    metadata$entries$parameter %in% c("mu_f[a]", "mu_Signal"), , drop = FALSE]
  provider$metadata <- function() metadata
  eligibility <- .hypothesis_model_quantities_model(provider, transform_scaled = FALSE)
  idx <- match(c("(mu) f[a]", "(mu) Signal{total}"), eligibility$alias)
  expect_false(anyNA(idx))
  expect_identical(eligibility$direction_test[idx], c(FALSE, TRUE))
  expect_false(eligibility$point_test[idx[1L]])
  expect_match(eligibility$reason[idx[1L]], "combination of several coded coefficients", fixed = TRUE)

  expect_error(hypothesis_BF_model(fit, "f[a] > 0.1"), class = "BayesTools_R2_density_unavailable")
  expect_s3_class(hypothesis_BF_model(fit, "Signal{total} > 0.2"), "BayesTools_hypothesis_BF")
})

test_that("a never-open term reports its local refusal without hiding the budget total", {

  set.seed(786)
  data <- data.frame(x = seq_len(8), z = rep(c(-1, 1), 4))
  info <- JAGS_formula(~ x + z, "mu", data, list(
    intercept = prior("normal", list(0, 1)),
    Signal = prior_R2(inclusion = list(x = prior("spike", list(0.5))), model = R2_model(sigma = 1))
  ))
  samples <- .generate_prior_sample_matrix(info$prior_list, 200L)
  closed <- samples[samples[, "mu_Signal_open"] == 0, , drop = FALSE]
  fit <- .parameter_catalog_test_fit(coda::mcmc.list(coda::mcmc(closed)), info$prior_list,
                                     formula_design = list(mu = info$formula_design))
  provider <- hypothesis_provider(fit)
  metadata <- provider$metadata()
  wanted <- c("(mu) Signal[share: x]", "(mu) Signal{total}")
  labels <- vapply(metadata$entries$aliases, `[[`, character(1), 1L)
  metadata$entries <- metadata$entries[labels %in% wanted, , drop = FALSE]
  provider$metadata <- function() metadata
  eligibility <- .hypothesis_model_quantities_model(provider, transform_scaled = FALSE)
  ids <- match(wanted, eligibility$alias)
  expect_false(anyNA(ids))
  expect_identical(eligibility$direction_test[ids], c(FALSE, TRUE))
  expect_match(eligibility$reason[ids[1L]], "no posterior draw defines the quantity", fixed = TRUE)
  expect_error(hypothesis_BF_model(fit, "Signal[share: x] > 0.2"), class = "BayesTools_no_defined_posterior_draws")
})

test_that("plot samples in original predictor units retain the fitted projection and prior", {

  set.seed(7933)
  data <- data.frame(x = rep(c(-4, -2, 0, 2, 4), 3), f = factor(rep(letters[1:3], each = 5)))
  fit_formula <- function(formula, priors) {
    info <- JAGS_formula(formula, "mu", data, priors, formula_scale = TRUE)
    samples <- .generate_prior_sample_matrix(info$prior_list, 150L)
    .parameter_catalog_test_fit(coda::mcmc.list(coda::mcmc(samples)), info$prior_list,
                                formula_scale = list(mu = info$formula_scale),
                                formula_design = list(mu = info$formula_design))
  }
  original_factor <- 1 / stats::sd(data$x)
  check_projection <- function(fit, alias) {
    selection <- parameter_catalog_resolve(parameter_catalog(fit), alias, namespace = "mu")
    key <- selection$quantities$extraction_key[[1L]]
    weights <- if(identical(key$type, "coordinate")) 1 else key$weights
    out <- parameter_plot_samples(fit, selection)
    name <- selection$quantities$canonical_name
    q <- posterior_metadata(out[[name]], "quantities")
    # sigma = 2.5 scales the linear predictor. The selected coefficients
    # remain in residual SDs per original predictor unit.
    expect_identical(q$quantity_id, selection$quantity_id)
    expect_identical(q$dependencies[[1L]], key$dependencies)
    expect_identical(q$weights[[1L]], as.numeric(original_factor * weights))
    coordinates <- as.matrix(JAGS_materialize_draws(fit, parameters = key$dependencies))
    expect_equal(as.numeric(out[[name]]), as.vector(coordinates %*% q$weights[[1L]]), tolerance = 1e-13)
    expect_identical(.bt_meta_get(out, "transform_scaled"), TRUE)
    densities <- posterior_metadata(out, "prior_densities")
    expect_named(densities, name)
    expect_identical(densities[[name]], posterior_metadata(out[[name]], "prior_density"))
    out
  }

  scalar <- fit_formula(~ x, list(intercept = prior("normal", list(0, 1)),
    Signal = prior_R2(total = prior("beta", list(1, 3)), model = R2_model(sigma = 2.5))))
  original <- check_projection(scalar, "x")
  sel <- parameter_catalog_resolve(parameter_catalog(scalar), "x")
  fitted <- parameter_plot_samples(scalar, sel, transform_scaled = FALSE)
  expect_identical(posterior_metadata(fitted[[1L]], "quantities")$weights[[1L]], 1)
  expect_null(.bt_meta_get(fitted, "transform_scaled"))
  expect_null(posterior_metadata(fitted, "prior_densities"))
  total <- parameter_catalog_resolve(parameter_catalog(scalar), "Signal{total}")
  expect_null(.bt_meta_get(parameter_plot_samples(scalar, total), "transform_scaled"))
  # The rendered prior scales by the same Jacobian as its analytic density.
  graph <- function(out, limits) {
    ggplot2::ggplot_build(plot_posterior(out, "mu_x", plot_type = "ggplot",
      prior = TRUE, n_points = 40L, xlim = limits))$data[[2L]]
  }
  limits <- c(-1, 1)
  before <- graph(fitted, limits)
  after <- graph(original, original_factor * limits)
  expect_equal(after$x, original_factor * before$x, tolerance = 1e-10)
  expect_equal(after$y, before$y / original_factor, tolerance = 1e-4)

  # This level is a projection of both indexed interaction coefficients.
  interaction <- fit_formula(~ x * f, list(
    intercept = prior("normal", list(0, 1)), x = prior("normal", list(0, 1)),
    f = prior_factor("mnormal", list(0, 1), contrast = "orthonormal"),
    Signal = prior_R2(terms = "x:f", total = prior("beta", list(1, 3)),
                      model = R2_model(sigma = 2.5, contrast = "orthonormal"))
  ))
  cell <- check_projection(interaction, "x:f[a]")
  q <- posterior_metadata(cell[[1L]], "quantities")
  expect_equal(length(q$dependencies[[1L]]), 2L)
  expect_true(all(q$weights[[1L]] != 0))
})
