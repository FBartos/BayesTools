skip_if_not_test_profile("fixture")

# ============================================================================ #
# TEST FILE: Random-effect summary posteriors of fitted models
# ============================================================================ #
#
# PURPOSE:
#   Atoms, gate states, point-free status, and prior densities of catalog
#   quantities of real fits: composite and allocated original-scale
#   correlations, monitored linear predictors, and the point components of
#   mixture and spike-and-slab formula and random-effect SD priors.
#
# MODELS/FIXTURES:
#   - fit_re_summary_* models from test-00-model-fits.R
#
# TAGS: @fixture, @JAGS, @random-effects
# ============================================================================ #

source(testthat::test_path("common-functions.R"))

.re_summary_cached_fit <- function(name){

  skip_if_not_installed("rjags")
  skip_if_missing_fits(name)
  readRDS(file.path(temp_fits_dir, paste0(name, ".RDS")))
}

test_that("parameter_mixed_posterior declares atoms from structure without a prior density", {

  # the original-scale correlation of a us() block with a scaled slope is a
  # composite of the LKJ primitive and both SDs: no prior density
  fit <- .re_summary_cached_fit("fit_re_summary_composite")
  catalog <- parameter_catalog(fit)
  correlation <- parameter_catalog_resolve(catalog, "(mu) cor(intercept,x)")
  expect_identical(correlation$quantities$source_type, "composite")
  expect_null(parameter_prior_density(fit, correlation))

  # none of its fitted coordinates takes a point mass: declared atom-free
  mixed <- parameter_mixed_posterior(fit, correlation)
  expect_true(posterior_atoms_free(mixed))
  expect_null(parameter_gate_states(fit, correlation))
  samples <- list(rho = mixed)
  class(samples) <- c("as_mixed_posteriors", "mixed_posteriors", "list")
  expect_s3_class(plot_posterior(samples, "rho", plot_type = "ggplot"), "ggplot")
  expect_no_error(plot_posterior(samples, "rho"))

  # region hypotheses are the posterior over prior odds of the draws, as
  # before
  set.seed(2)
  prior_draws <- stats::runif(2000, -1, 1)
  region <- hypothesis_BF(
    posterior  = data.frame(rho = as.numeric(mixed)),
    prior      = data.frame(rho = prior_draws),
    hypothesis = "rho > 0"
  )
  odds <- function(p) p / (1 - p)
  expect_equal(
    as.numeric(region$BF),
    odds(mean(as.numeric(mixed) > 0)) / odds(mean(prior_draws > 0)),
    tolerance = 1e-12
  )
})

test_that("correlations of allocated blocks are point-free up to common factors of their SDs", {

  # Four us() blocks with a scaled slope: g splits a continuous scale prior
  # into its SDs (no gate); s splits the SD of a gate-only allocation with an
  # inclusion gate (the gate multiplies both SDs of the block); t splits a
  # scale prior with a spike at 0 (the source multiplies both SDs); p has a
  # spike-and-slab prior on each SD (a gate per SD).
  fit <- .re_summary_cached_fit("fit_re_summary_allocated")
  draws <- as.matrix(BayesTools:::.fit_to_posterior(fit))
  coordinates <- parameter_coordinates(fit)
  nodes <- JAGS_deterministic_nodes(fit)
  source <- function(allocation) paste0("mu__xRE_ALLOCx_", allocation, "__allocation_sd")
  gate <- "mu__xRE_ALLOCx_s_gate__include_s_indicator"

  # The allocated SDs are derived coordinates of registered 'random_sd' nodes
  # computed from the scale prior, the Dirichlet weights, and the gates.
  sd_names <- paste0("mu__xREx__", c("g", "s", "t"), "_x")
  expect_identical(
    coordinates$convergence_role[match(sd_names, coordinates$coordinate_name)],
    rep("derived", 3L)
  )
  expect_identical(nodes$family[match(sd_names, nodes$node)], rep("random_sd", 3L))
  expect_identical(
    nodes$dependencies[[match("mu__xREx__g_x", nodes$node)]],
    c(source("g_split"),
      "mu__xRE_ALLOCx_g_split__weight[1]", "mu__xRE_ALLOCx_g_split__weight[2]")
  )
  expect_true(gate %in% nodes$dependencies[[match("mu__xREx__s_x", nodes$node)]])
  # a node is point-free when every dependency is: a continuous scale prior
  # and Dirichlet weights are; an inclusion gate or a scale prior with a point
  # component propagates to the SD
  context <- BayesTools:::.bt_parameter_point_free_context(fit)
  point_free <- vapply(sd_names, function(name){
    BayesTools:::.bt_parameter_coordinate_point_free(context, name)
  }, logical(1), USE.NAMES = FALSE)
  expect_identical(point_free, c(TRUE, FALSE, FALSE))

  # The original-scale correlation of each block mixes the LKJ primitive and
  # both SDs: no prior density and no gate plan.
  resolve <- function(block, object = fit){
    parameter_catalog_resolve(parameter_catalog(object),
                              paste0("(mu) ", block, ": cor(intercept,x)"))
  }
  for(block in c("g", "s", "t", "p")){
    selection <- resolve(block)
    expect_identical(selection$quantities$source_type, "composite", info = block)
    expect_null(parameter_prior_density(fit, selection), info = block)
    expect_null(parameter_gate_states(fit, selection), info = block)
  }
  # The scale source and the gate of s enter both SDs of their block only as
  # one common factor; the gates of p act on one SD each.
  common_factors <- function(block){
    BayesTools:::.bt_parameter_correlation_common_factors(
      BayesTools:::.bt_parameter_point_free_context(fit),
      resolve(block)$quantities
    )
  }
  expect_identical(common_factors("g"), source("g_split"))
  expect_identical(common_factors("s"), c(source("s_gate"), gate))
  expect_identical(common_factors("t"), source("t_split"))
  expect_identical(common_factors("p"), character())

  # A correlation is scale-invariant: a common factor of both SDs leaves it
  # unchanged where it is positive and makes it undefined where it is 0 (all
  # SDs 0; the draw is dropped). With continuous scales (g), a whole-block
  # gate (s), or a spike scale source (t), the correlation is atom-free on the
  # defined draws and plots; without a prior density, prior = TRUE draws the
  # posterior alone with the classed warning.
  defined <- list(
    g = rep(TRUE, nrow(draws)),
    s = draws[, gate] == 1,
    t = draws[, source("t_split")] > 0
  )
  expect_true(all(vapply(defined[c("s", "t")], function(x) any(!x) && any(x), logical(1))))
  for(block in c("g", "s", "t")){
    name <- paste0("(mu) ", block, ": cor(intercept,x)")
    mixed <- parameter_mixed_posterior(fit, resolve(block))
    expect_length(mixed, sum(defined[[block]]))
    expect_true(posterior_atoms_free(mixed), info = block)
    expect_identical(posterior_metadata(mixed, "atoms")$source, "parameter_structure",
                     info = block)
    samples <- stats::setNames(list(mixed), name)
    posterior_only <- plot_posterior(samples, name, prior = FALSE, plot_type = "ggplot")
    expect_s3_class(posterior_only, "ggplot")
    with_prior <- NULL
    expect_warning(
      with_prior <- plot_posterior(samples, name, prior = TRUE, plot_type = "ggplot"),
      class = "BayesTools_prior_curve_unavailable"
    )
    expect_identical(ggplot2::ggplot_build(with_prior)$data,
                     ggplot2::ggplot_build(posterior_only)$data, info = block)
  }

  # Gates of single SDs (p): a draw with one SD 0 has a singular
  # original-scale covariance and a correlation of -1 or 1 (up to rounding),
  # so the atom status stays undeclared and plots stop.
  name <- "(mu) p: cor(intercept,x)"
  mixed <- parameter_mixed_posterior(fit, resolve("p"))
  expect_true(any(abs(as.numeric(mixed)) > 1 - 1e-8))
  expect_null(posterior_metadata(mixed, "atoms"))
  expect_error(
    plot_posterior(stats::setNames(list(mixed), name), name, plot_type = "ggplot"),
    "Posterior atom status is unknown",
    fixed = TRUE
  )

  # On the fitted scale the correlations have the exact LKJ marginal prior
  # density and no atoms, and prior = TRUE draws the prior.
  fitted <- fit
  attr(fitted, "formula_scale") <- list()
  for(block in c("g", "s", "t", "p")){
    name <- paste0("(mu) ", block, ": cor(intercept,x)")
    mixed <- parameter_mixed_posterior(fitted, resolve(block, fitted))
    expect_true(posterior_atoms_free(mixed), info = block)
    expect_false(is.null(posterior_metadata(mixed, "prior_density")), info = block)
    plot <- plot_posterior(stats::setNames(list(mixed), name), name, prior = TRUE,
                           plot_type = "ggplot")
    expect_length(plot$layers, 2L)
  }
})

test_that("monitored linear predictors are not point-free through their registry dependencies", {

  # mu[i] = mu_x * x[i] + mu_z * z[i] is a registered node whose dependencies
  # are continuous, yet it is 0 in every draw of the row with x = z = 0: its
  # atom status stays undeclared rather than atom-free.
  fit <- .re_summary_cached_fit("fit_re_summary_linear_predictor")
  nodes <- JAGS_deterministic_nodes(fit)
  expect_identical(nodes$family[match("mu", nodes$node)], "linear_predictor")
  expect_identical(nodes$dependencies[[match("mu", nodes$node)]], c("mu_x", "mu_z"))
  context <- BayesTools:::.bt_parameter_point_free_context(fit)
  expect_true(BayesTools:::.bt_parameter_coordinate_point_free(context, "mu_x"))
  expect_false(BayesTools:::.bt_parameter_coordinate_point_free(context, "mu[1]"))

  catalog  <- parameter_catalog(fit)
  constant <- parameter_mixed_posterior(fit, parameter_catalog_resolve(catalog, "mu[1]", "mu"))
  expect_length(unique(as.numeric(constant)), 1L)
  expect_null(posterior_metadata(constant, "atoms"))
  expect_false(posterior_atoms_free(constant))
})

test_that("parameter_mixed_posterior declares the point components of mixture priors", {

  fit <- .re_summary_cached_fit("fit_re_summary_point_components")
  draws <- as.matrix(BayesTools:::.fit_to_posterior(fit))
  catalog <- parameter_catalog(fit)
  resolve <- function(name) parameter_catalog_resolve(catalog, name)
  atoms_of <- function(name){
    atoms <- posterior_metadata(parameter_mixed_posterior(fit, resolve(name)), "atoms")
    list(x = as.numeric(atoms$locations[, 1L]), mass = atoms$mass)
  }

  # the spike of each spike-and-slab prior (indicator 0) and the spike(0.2)
  # component of the mixture (component 2) are atoms with the posterior share
  # of the draws whose indicator selects them
  spike_x  <- mean(draws[, "mu_x_indicator"] == 0)
  spike_f  <- mean(draws[, "mu_f_indicator"] == 0)
  spike_sd <- mean(draws[, "mu__xREx__g_intercept_indicator"] == 0)
  spike_intercept <- mean(draws[, "mu_intercept_indicator"] == 2)
  expect_true(all(c(spike_x, spike_f, spike_sd, spike_intercept) > 0 &
                    c(spike_x, spike_f, spike_sd, spike_intercept) < 1))
  expect_identical(atoms_of("mu_x"), list(x = 0, mass = spike_x))
  expect_identical(atoms_of("mu_f[b]"), list(x = 0, mass = spike_f))
  expect_identical(atoms_of("mu_f[c]"), list(x = 0, mass = spike_f))
  expect_identical(atoms_of("mu_intercept"), list(x = 0.2, mass = spike_intercept))
  expect_identical(atoms_of("(mu) sd(intercept)"), list(x = 0, mass = spike_sd))
  expect_identical(atoms_of("(mu) var(intercept)"), list(x = 0, mass = spike_sd))
  # the atoms are where the draws are exactly on the point components
  expect_equal(mean(draws[, "mu_x"] == 0), spike_x)
  expect_equal(mean(draws[, "mu_intercept"] == 0.2), spike_intercept)

  # the per-draw states, from the posterior draws or from given draws
  states <- parameter_gate_states(fit, resolve("mu_x"))
  expect_true(states$known)
  expect_identical(states$atom, ifelse(draws[, "mu_x_indicator"] == 0, 0, NA_real_))
  expect_identical(states$continuous, draws[, "mu_x_indicator"] == 1)
  expect_null(states$event)
  row <- parameter_gate_states(fit, resolve("(mu) var(intercept)"), draws = draws[5L, ])
  expect_identical(row$atom, if(draws[5L, "mu__xREx__g_intercept_indicator"] == 0) 0 else NA_real_)
  expect_error(
    parameter_gate_states(fit, resolve("mu_x"), draws = draws[, "mu_x", drop = FALSE]),
    "The draws do not contain the inclusion indicators 'mu_x_indicator' of the quantity.",
    fixed = TRUE
  )

  # remaining undeclared: auxiliary nodes of mixture priors
  expect_null(posterior_metadata(parameter_mixed_posterior(fit, resolve("mu_x_variable")), "atoms"))

  # the inclusion indicator of the spike-and-slab SD prior is Bernoulli on
  # {0, 1} with the prior inclusion probability 0.5 (the location of its
  # spike(0.5) inclusion prior), not the SD's prior; its atoms are the shares
  # of the indicator draws
  inclusion <- resolve("(mu) inclusion(sd(intercept))")
  support <- inclusion$quantities$support[[1L]]
  expect_identical(support$bounds, c(0, 1))
  expect_identical(support$points, c(0, 1))
  expect_identical(support$type, "points")
  expect_true(support$exact)
  inclusion_prior <- parameter_prior_density(fit, inclusion)
  expect_equal(inclusion_prior$points$x, c(0, 1))
  expect_equal(inclusion_prior$points$p, c(.5, .5), tolerance = 1e-12)
  expect_null(inclusion_prior$density)
  indicator <- draws[, "mu__xREx__g_intercept_indicator"]
  expect_identical(as.numeric(table(indicator)), c(158, 42))
  expect_identical(atoms_of("(mu) inclusion(sd(intercept))"),
                   list(x = c(0, 1), mass = c(158, 42) / 200))
  inclusion_draws <- parameter_mixed_posterior(fit, inclusion)
  expect_identical(as.numeric(inclusion_draws), as.numeric(indicator))
  expect_identical(parameter_gate_states(fit, inclusion)$atom, as.numeric(indicator))

  # the inclusion indicators of the components of a mixture SD prior are
  # Bernoulli with the components' prior weights 1:2:1
  mixture_fit <- .re_summary_cached_fit("fit_re_summary_mixture_sd")
  mixture_catalog <- parameter_catalog(mixture_fit)
  mixture_draws <- as.matrix(BayesTools:::.fit_to_posterior(mixture_fit))
  for(component in c("null", "narrow", "wide")){
    name <- paste0("(mu) inclusion(sd(intercept)[", component, "])")
    selection <- parameter_catalog_resolve(mixture_catalog, name)
    expect_identical(selection$quantities$support[[1L]]$points, c(0, 1), info = name)
    density <- parameter_prior_density(mixture_fit, selection)
    probability <- c(null = .25, narrow = .5, wide = .25)[[component]]
    expect_equal(density$points$p, c(1 - probability, probability), tolerance = 1e-12,
                 info = name)
    selected <- mixture_draws[, "mu__xREx__g_intercept_indicator"] ==
      match(component, c("null", "narrow", "wide"))
    atoms <- posterior_metadata(parameter_mixed_posterior(mixture_fit, selection), "atoms")
    expect_equal(atoms$mass[as.numeric(atoms$locations[, 1L]) == 1], mean(selected),
                 info = name)
    expect_equal(sum(atoms$mass), 1, info = name)
  }
})
