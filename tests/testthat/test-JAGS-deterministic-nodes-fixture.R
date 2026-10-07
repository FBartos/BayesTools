skip_if_not_test_profile("fixture")

# ============================================================================ #
# TEST FILE: Deterministic nodes of fitted models
# ============================================================================ #
#
# PURPOSE:
#   The parity tests remove the monitored nodes from the posterior draws of
#   small real fits and rebuild them with JAGS_evaluate_deterministic(): the R
#   evaluators must reproduce the JAGS monitors (bit-identical where R and
#   JAGS perform the same floating-point operations, otherwise within 1e-14).
#
# MODELS/FIXTURES:
#   - fit_dnode_* models from test-00-model-fits.R
#
# TAGS: @fixture, @JAGS, @deterministic-nodes
# ============================================================================ #

source(testthat::test_path("common-functions.R"))

.dnode_cached_fit <- function(name){

  skip_if_not_installed("rjags")
  skip_if_missing_fits(name)
  readRDS(file.path(temp_fits_dir, paste0(name, ".RDS")))
}

test_that("expression point coefficients replay genuine JAGS monitors through their registered family", {
  fit <- .dnode_cached_fit("fit_dnode_point_expression")
  posterior <- as.matrix(.fit_to_posterior(fit))
  reduced <- posterior[, setdiff(colnames(posterior), c("mu_x", paste0("mu[", 1:24, "]"))), drop = FALSE]
  rebuilt <- JAGS_evaluate_deterministic(fit, reduced, nodes = c("mu_x", "mu"))
  expect_equal(as.numeric(rebuilt[, "mu_x"]), 2 * posterior[, "theta"] + .25, tolerance = 1e-14)
  expect_equal(as.vector(rebuilt), as.vector(posterior[, colnames(rebuilt), drop = FALSE]), tolerance = 1e-14)
  expect_equal(as.vector(JAGS_evaluate_formula(fit, parameter = "mu")),
    as.vector(t(posterior[, paste0("mu[", 1:24, "]"), drop = FALSE])), tolerance = 1e-14)
})

# Removes every monitored coordinate of the family's nodes from the posterior,
# rebuilds them, and returns the rebuilt and monitored columns.
.dnode_rebuild <- function(fit, family){

  posterior <- as.matrix(BayesTools:::.fit_to_posterior(fit))
  nodes <- JAGS_deterministic_nodes(fit)
  nodes <- nodes[nodes$family == family & nodes$monitored, , drop = FALSE]
  coordinates <- unique(unlist(nodes$coordinates))
  reduced <- posterior[, setdiff(colnames(posterior), coordinates), drop = FALSE]
  rebuilt <- JAGS_evaluate_deterministic(fit, draws = reduced, nodes = nodes$node)

  list(
    nodes = nodes,
    rebuilt = rebuilt[, coordinates, drop = FALSE],
    monitored = posterior[, coordinates, drop = FALSE]
  )
}

test_that("allocation SD nodes reproduce the JAGS monitors of gated, nested, and external-source allocations", {

  # A gated total-variance root allocation over the blocks g and d, with the
  # correlated block g split into SD components by a mean-variance child.
  fit <- .dnode_cached_fit("fit_dnode_allocation")
  syntax <- attr(fit, "model_syntax")

  result <- .dnode_rebuild(fit, "random_sd")
  expect_setequal(
    result$nodes$node,
    c("mu__xREx__g_intercept", "mu__xREx__g_x", "mu__xREx__d_intercept")
  )
  expect_identical(
    unname(unlist(result$nodes$dependencies[match("mu__xREx__d_intercept", result$nodes$node)])),
    c(
      "mu__xRE_ALLOCx_tot__allocation_sd",
      "mu__xRE_ALLOCx_tot__weight[1]", "mu__xRE_ALLOCx_tot__weight[2]",
      "mu__xRE_ALLOCx_tot__include_d_indicator"
    )
  )
  # The factors are multiplied in the order of the model syntax and the gates
  # are 0/1, so the rebuilt SDs are bit-identical to the monitors.
  expect_identical(unname(result$rebuilt), unname(result$monitored))

  # The consumed parent component is a node of its own (not monitored), and
  # every node is emitted into the model syntax from its specification.
  nodes <- JAGS_deterministic_nodes(fit)
  component <- nodes[nodes$node == "mu__xRE_ALLOCx_tot__component_g_sd", , drop = FALSE]
  expect_identical(nrow(component), 1L)
  expect_false(component$monitored)
  registry <- BayesTools:::.bt_deterministic_nodes_fit(fit)
  for(node in registry[nodes$node[nodes$family == "random_sd"]]){
    expect_true(grepl(BayesTools:::.bt_deterministic_node_emit(node), syntax, fixed = TRUE))
  }
  expect_true(grepl(
    "mu__xREx__g_x = mu__xRE_ALLOCx_tot__component_g_sd * sqrt(2 * mu__xRE_ALLOCx_gc__weight[2])",
    syntax,
    fixed = TRUE
  ))

  # Prior draws evaluate the same chain.
  prior_draws <- transform_prior_samples(fit, n_samples = 300, seed = 5, formula_scale = list())
  expect_identical(
    unname(prior_draws[, "mu__xREx__g_x"]),
    unname(
      prior_draws[, "mu__xRE_ALLOCx_tot__allocation_sd"] *
        sqrt(prior_draws[, "mu__xRE_ALLOCx_tot__weight[1]"]) *
        sqrt(2 * prior_draws[, "mu__xRE_ALLOCx_gc__weight[2]"])
    )
  )
  expect_identical(
    unname(prior_draws[, "mu__xREx__d_intercept"]),
    unname(
      prior_draws[, "mu__xRE_ALLOCx_tot__allocation_sd"] *
        sqrt(prior_draws[, "mu__xRE_ALLOCx_tot__weight[2]"]) *
        prior_draws[, "mu__xRE_ALLOCx_tot__include_d_indicator"]
    )
  )

  # An external scalar SD source is the root of the chain.
  external <- .dnode_cached_fit("fit_dnode_allocation_external")
  external_result <- .dnode_rebuild(external, "random_sd")
  expect_identical(
    unname(unlist(external_result$nodes$dependencies[1L])),
    c("tau", "mu__xRE_ALLOCx_tot__weight[1]", "mu__xRE_ALLOCx_tot__weight[2]")
  )
  expect_identical(unname(external_result$rebuilt), unname(external_result$monitored))
})

test_that("scalar correlation nodes reproduce the JAGS monitors of every structure", {

  fit <- .dnode_cached_fit("fit_dnode_correlation")
  syntax <- attr(fit, "model_syntax")

  result <- .dnode_rebuild(fit, "random_rho")
  expect_setequal(
    result$nodes$node,
    paste0("mu__xREx__", c("g", "s", "d", "p"), "_rho")
  )
  expect_identical(
    unname(unlist(result$nodes$dependencies[match("mu__xREx__g_rho", result$nodes$node)])),
    "mu__xREx__g_rho_logit"
  )
  expect_identical(
    unname(unlist(result$nodes$dependencies[match("mu__xREx__s_rho", result$nodes$node)])),
    "mu__xREx__s_rho_z"
  )
  # tanh() and the logistic map are the same floating-point operations in R and
  # JAGS, so the rebuilt correlations are bit-identical to the monitors.
  expect_identical(unname(result$rebuilt), unname(result$monitored))

  # The model syntax is the emitted node definition.
  registry <- BayesTools:::.bt_deterministic_nodes_fit(fit)
  for(node in registry[result$nodes$node]){
    expect_true(grepl(
      BayesTools:::.bt_deterministic_node_emit(node),
      syntax,
      fixed = TRUE
    ))
  }
  expect_true(grepl(
    "mu__xREx__g_rho <- -0.33333333333333331 + 1.3333333333333333 * ilogit(mu__xREx__g_rho_logit)",
    syntax,
    fixed = TRUE
  ))

  # A fixed Fisher-z correlation is evaluated from its point prior.
  expect_identical(
    unname(result$rebuilt[, "mu__xREx__p_rho"]),
    rep(tanh(0.4), nrow(result$rebuilt))
  )

  # Convergence roles read the declared parents of the generated nodes: the
  # fixed correlation depends only on its point prior (structural), the others
  # on sampled coordinates (derived).
  parents <- BayesTools:::.bt_deterministic_node_parent_bases(registry)
  expect_identical(parents[["mu__xREx__g_rho"]], "mu__xREx__g_rho_logit")
  roles <- parameter_coordinates(fit)
  role <- function(name) roles$convergence_role[roles$coordinate_name == name]
  expect_identical(role("mu__xREx__p_rho"), "structural")
  expect_identical(role("mu__xREx__g_rho"), "derived")
  expect_identical(role("mu__xREx__s_rho"), "derived")

  # Prior draws carry the same node definition.
  prior_draws <- transform_prior_samples(fit, n_samples = 200, seed = 3, formula_scale = list())
  expect_identical(
    unname(prior_draws[, "mu__xREx__s_rho"]),
    unname(tanh(prior_draws[, "mu__xREx__s_rho_z"]))
  )
  expect_identical(
    unname(prior_draws[, "mu__xREx__d_rho"]),
    unname(0 + 1 * stats::plogis(prior_draws[, "mu__xREx__d_rho_logit"]))
  )
})

test_that("LKJ Cholesky, correlation, and partial-correlation nodes reproduce the JAGS monitors", {

  fit <- .dnode_cached_fit("fit_dnode_lkj")
  syntax <- attr(fit, "model_syntax")

  result <- .dnode_rebuild(fit, "lkj")
  expect_setequal(result$nodes$node, c("mu__xREx__g_xRE_CORx", "mu__xREx__s_xRE_CORx"))
  g <- match("mu__xREx__g_xRE_CORx", result$nodes$node)
  expect_identical(
    unname(unlist(result$nodes$dependencies[g])),
    paste0("mu__xREx__g_xRE_CORx_lkj_u[", 1:3, "]")
  )
  # K = 3: 9 Cholesky cells, 9 correlation cells, and 3 partial correlations.
  expect_length(unlist(result$nodes$coordinates[g]), 21L)
  # The Cholesky factor comes from the module's own kernel, R = L L' sums the
  # products in the module's order, and cpc = 2 u - 1: all bit-identical.
  expect_identical(unname(result$rebuilt), unname(result$monitored))

  registry <- BayesTools:::.bt_deterministic_nodes_fit(fit)
  for(node in registry[result$nodes$node]){
    expect_true(grepl(
      paste(BayesTools:::.bt_deterministic_node_emit(node), collapse = "\n"),
      syntax,
      fixed = TRUE
    ))
  }

  # Prior draws evaluate the same node on the drawn primitives.
  prior_draws <- transform_prior_samples(fit, n_samples = 200, seed = 4, formula_scale = list())
  u_names <- paste0("mu__xREx__g_xRE_CORx_lkj_u[", 1:3, "]")
  L <- BayesTools:::.bt_lkj_cholesky_cpc_u_to_L(prior_draws[, u_names], K = 3L)
  expect_identical(unname(prior_draws[, "mu__xREx__g_xRE_CORx_L[3,2]"]), L[, 3, 2])
  expect_identical(
    unname(prior_draws[, "mu__xREx__g_xRE_CORx_lkj_cpc[2]"]),
    unname(2 * prior_draws[, u_names[2]] - 1)
  )
  expect_true(all(prior_draws[, "mu__xREx__s_xRE_CORx_R[2,2]"] == 1))
})

test_that("publication-weight nodes reproduce the JAGS monitors of every weight function", {

  fits <- lapply(c(
    cumulative_two_sided = "fit_dnode_omega_cumulative",
    binary               = "fit_dnode_omega_binary",
    log_independent      = "fit_dnode_omega_log_independent",
    fixed                = "fit_dnode_omega_fixed",
    bias_mixture         = "fit_dnode_omega_bias_mixture"
  ), .dnode_cached_fit)

  expect_identical(
    unname(unlist(JAGS_deterministic_nodes(fits$cumulative_two_sided)$dependencies)),
    paste0("eta[", 1:3, "]")
  )
  expect_identical(
    unname(unlist(JAGS_deterministic_nodes(fits$binary)$dependencies)),
    "omega[2]"
  )
  bias_nodes <- JAGS_deterministic_nodes(fits$bias_mixture)
  expect_identical(
    unname(unlist(bias_nodes$dependencies[bias_nodes$node == "omega"])),
    c("bias_indicator", paste0("eta_component_2[", 1:3, "]"), "omega_ratio_component_3",
      paste0("log_omega_component_4[", 2:3, "]"))
  )

  for(name in names(fits)){
    fit <- fits[[name]]
    posterior <- as.matrix(BayesTools:::.fit_to_posterior(fit))
    node <- JAGS_deterministic_nodes(fit)
    node <- node[node$family == "omega", , drop = FALSE]
    expect_identical(node$node, "omega")
    coordinates <- unlist(node$coordinates)
    # The sampled free weights stay; every other bin is rebuilt.
    removed <- setdiff(coordinates, unlist(node$dependencies))
    reduced <- posterior[, setdiff(colnames(posterior), removed), drop = FALSE]
    rebuilt <- JAGS_evaluate_deterministic(fit, draws = reduced, nodes = "omega")
    # Mapped, binary, and independent weights are copies of their free
    # coordinates (bit-identical). The cumulative weights normalize and sum
    # the gamma auxiliaries in a different order than JAGS, and fixed weights
    # are emitted as 15-digit literals (1/3): both differ by a few ulps.
    expect_equal(
      unname(rebuilt[, coordinates]),
      unname(posterior[, coordinates]),
      tolerance = 1e-14
    )
    if(name %in% c("binary", "log_independent")){
      expect_identical(unname(rebuilt[, coordinates]), unname(posterior[, coordinates]))
    }

    # The registered emitter is the syntax of the fitted JAGS model.
    model <- JAGS_add_priors(attr(fit, "model_syntax"), attr(fit, "prior_list"))
    registry <- BayesTools:::.bt_deterministic_nodes_fit(fit)
    for(piece in BayesTools:::.bt_deterministic_node_emit(registry$omega)){
      expect_true(grepl(piece, model, fixed = TRUE))
    }

    # The marginal-likelihood parameters of a draw are the node's values.
    if(name != "bias_mixture"){
      prior_list <- attr(fit, "prior_list")
      expect_identical(
        JAGS_marglik_parameters(posterior[7, ], prior_list["omega"])$omega,
        unname(JAGS_evaluate_deterministic(fit, draws = posterior[7, ], nodes = "omega")[1, ])
      )
    }
  }
})

test_that("spike-and-slab and mixture nodes reproduce the JAGS monitors and give the marginal-likelihood parameters", {

  fit <- .dnode_cached_fit("fit_dnode_mixture")
  posterior <- as.matrix(BayesTools:::.fit_to_posterior(fit))
  nodes <- JAGS_deterministic_nodes(fit)
  expect_setequal(nodes$node, c("a", "b", "omega", "PET", "PEESE"))
  expect_identical(unname(unlist(nodes$dependencies[nodes$node == "a"])), c("a_variable", "a_indicator"))
  expect_identical(
    unname(unlist(nodes$dependencies[nodes$node == "b"])),
    c("b_indicator", "b_component_1", "b_component_3")
  )

  mixture_nodes <- nodes$node[nodes$family == "prior_mixture"]
  reduced <- posterior[, setdiff(colnames(posterior), mixture_nodes), drop = FALSE]
  rebuilt <- JAGS_evaluate_deterministic(fit, draws = reduced, nodes = mixture_nodes)
  # Products with a 0/1 indicator and the active component: bit-identical.
  expect_identical(unname(rebuilt[, mixture_nodes]), unname(posterior[, mixture_nodes]))
  model <- JAGS_add_priors(attr(fit, "model_syntax"), attr(fit, "prior_list"))
  registry <- BayesTools:::.bt_deterministic_nodes_fit(fit)
  for(node in registry[mixture_nodes]){
    expect_true(grepl(BayesTools:::.bt_deterministic_node_emit(node), model, fixed = TRUE))
  }

  # JAGS_marglik_parameters() evaluates the same nodes for one draw, including
  # spike-and-slab, mixture, and publication-bias mixture priors.
  prior_list <- attr(fit, "prior_list")
  for(row in c(3L, 40L, 77L)){
    parameters <- JAGS_marglik_parameters(posterior[row, ], prior_list)
    expect_identical(parameters$a, unname(posterior[row, "a"]))
    expect_identical(parameters$b, unname(posterior[row, "b"]))
    expect_identical(parameters$PET, unname(posterior[row, "PET"]))
    expect_identical(parameters$PEESE, unname(posterior[row, "PEESE"]))
    expect_equal(parameters$omega, unname(posterior[row, paste0("omega[", 1:3, "]")]), tolerance = 1e-14)
  }
  row <- posterior[3L, setdiff(colnames(posterior), "a_variable")]
  expect_error(
    JAGS_marglik_parameters(row, prior_list["a"]),
    "'samples' does not contain all monitored spike-and-slab parameters of 'a'.",
    fixed = TRUE
  )
  expect_error(
    JAGS_marglik_parameters(posterior[3L, setdiff(colnames(posterior), "bias_indicator")], prior_list["bias"]),
    "'samples' does not contain all monitored bias-mixture parameters of 'bias'.",
    fixed = TRUE
  )
  # Bridge sampling still refuses the discrete mixture indicators.
  expect_error(
    JAGS_marglik_priors(posterior[3L, ], prior_list["a"]),
    "prior mixture priors is not implemented"
  )

  # Factor spike-and-slab and mixture formula priors, per coefficient; the
  # point component of the mixture is a constant.
  # Factor spike-and-slab and mixture formula priors, per coefficient; the
  # point component of the mixture is a constant.
  factor_fit <- .dnode_cached_fit("fit_dnode_mixture_factor")
  factor_nodes <- JAGS_deterministic_nodes(factor_fit)
  expect_identical(
    unname(unlist(factor_nodes$coordinates[factor_nodes$node == "mu_t"])),
    paste0("mu_t[", 1:3, "]")
  )
  factor_posterior <- as.matrix(BayesTools:::.fit_to_posterior(factor_fit))
  coordinates <- c("mu_x", paste0("mu_t[", 1:3, "]"))
  factor_rebuilt <- JAGS_evaluate_deterministic(
    factor_fit,
    draws = factor_posterior[, setdiff(colnames(factor_posterior), coordinates), drop = FALSE],
    nodes = c("mu_x", "mu_t")
  )
  expect_identical(unname(factor_rebuilt[, coordinates]), unname(factor_posterior[, coordinates]))
})

test_that("linear predictor nodes reproduce the monitored formula output", {

  fit <- .dnode_cached_fit("fit_dnode_linear_predictor")
  nodes <- JAGS_deterministic_nodes(fit)
  node <- nodes[nodes$family == "linear_predictor", , drop = FALSE]
  expect_identical(node$node, "mu")
  expect_true(node$monitored)
  dependencies <- unlist(node$dependencies)
  expect_true(all(c("mu_intercept", "mu_x", "mu_d[1]", "mu_d[2]", "b_scale",
                    "mu__xREx__g_xRE_Zx[1,1]", "mu__xREx__g_intercept",
                    "mu__xREx__g_xRE_CORx_L[2,1]", "mu__xREx__s_rho") %in% dependencies))

  # The formula output of the model syntax is the emitted node.
  registry <- BayesTools:::.bt_deterministic_nodes_fit(fit)
  expect_true(grepl(
    BayesTools:::.bt_deterministic_node_emit(registry$mu),
    attr(fit, "model_syntax"),
    fixed = TRUE
  ))

  # Rebuilt from the coefficients, the multiplier, the expression, and the
  # latent random effects. The data JAGS reads are serialized with 15
  # significant digits and the random-effect contributions are summed in
  # another order, so the rebuilt predictor agrees to a few ulps.
  posterior <- as.matrix(BayesTools:::.fit_to_posterior(fit))
  coordinates <- unlist(node$coordinates)
  reduced <- posterior[, setdiff(colnames(posterior), coordinates), drop = FALSE]
  rebuilt <- JAGS_evaluate_deterministic(fit, draws = reduced, nodes = "mu")
  expect_lte(max(abs(rebuilt[, coordinates] - posterior[, coordinates])), 1e-14)

  # Without the latent effects the node is unavailable.
  latent <- grepl("_xRE_Zx", colnames(reduced), fixed = TRUE)
  expect_error(
    JAGS_evaluate_deterministic(fit, draws = reduced[, !latent, drop = FALSE], nodes = "mu"),
    "Deterministic node 'mu' is unavailable from 'draws'",
    fixed = TRUE
  )

  # A formula design without the priors of its terms defines no linear
  # predictor node.
  design <- attr(fit, "formula_design")
  design$mu$prior_list <- NULL
  attr(fit, "formula_design") <- design
  expect_false("mu" %in% JAGS_deterministic_nodes(fit)$node)
})

test_that("linear predictor nodes of mean-centered random intercepts are evaluated from the latent effects", {

  fit <- .dnode_cached_fit("fit_dnode_mean_centered")
  posterior <- as.matrix(BayesTools:::.fit_to_posterior(fit))
  nodes <- JAGS_deterministic_nodes(fit)
  node <- nodes[nodes$node == "mu", , drop = FALSE]

  # The model syntax adds the group locations in place of the intercept; the
  # evaluator adds the latent deviations to the intercept, so the declared
  # dependencies are the monitored latent effects and SD, not the locations.
  registry <- BayesTools:::.bt_deterministic_nodes_fit(fit)
  expect_true(grepl(
    BayesTools:::.bt_deterministic_node_emit(registry$mu),
    attr(fit, "model_syntax"),
    fixed = TRUE
  ))
  dependencies <- unlist(node$dependencies)
  expect_true(all(dependencies %in% colnames(posterior)))
  sd_name <- unique(attr(fit, "formula_design")$mu$random_effects[[1L]]$sd_parameter_names)
  expect_true(all(c("mu_intercept", "mu_x", "mu__xREx__g_xRE_Zx[1,1]", sd_name) %in% dependencies))
  expect_false(any(grepl("_xRE_MEANx", dependencies, fixed = TRUE)))

  coordinates <- unlist(node$coordinates)
  reduced <- posterior[, setdiff(colnames(posterior), coordinates), drop = FALSE]
  rebuilt <- JAGS_evaluate_deterministic(fit, draws = reduced, nodes = "mu")
  expect_lte(max(abs(rebuilt[, coordinates] - posterior[, coordinates])), 1e-14)
})

test_that("JAGS_evaluate_deterministic() validates its node selection", {

  fit <- .dnode_cached_fit("fit_dnode_ar1")
  nodes <- JAGS_deterministic_nodes(fit)
  expect_identical(
    names(nodes),
    c("node", "family", "parameter", "block", "coordinates", "dependencies", "monitored")
  )
  # The scalar correlation of block s and the (unmonitored) linear predictor.
  expect_identical(nodes$node, c("mu__xREx__s_rho", "mu"))
  expect_identical(nodes$family, c("random_rho", "linear_predictor"))
  expect_identical(nodes$parameter, c("mu", "mu"))
  expect_identical(nodes$block, c("s", NA_character_))
  expect_identical(nodes$monitored, c(TRUE, FALSE))

  posterior <- as.matrix(BayesTools:::.fit_to_posterior(fit))
  expect_error(
    JAGS_evaluate_deterministic(fit, nodes = "mu__xREx__s_sd"),
    "'nodes' contains names that are not generated deterministic nodes of 'fit': 'mu__xREx__s_sd'. See JAGS_deterministic_nodes(fit) for the available nodes.",
    fixed = TRUE
  )
  reduced <- posterior[, setdiff(colnames(posterior), c("mu__xREx__s_rho", "mu__xREx__s_rho_z")), drop = FALSE]
  expect_error(
    JAGS_evaluate_deterministic(fit, draws = reduced, nodes = "mu__xREx__s_rho"),
    "Deterministic node 'mu__xREx__s_rho' is unavailable from 'draws': its dependencies 'mu__xREx__s_rho_z' are not all available as columns of 'draws' or as point priors.",
    fixed = TRUE
  )
  # Without a selection, unavailable nodes are skipped.
  expect_identical(dim(JAGS_evaluate_deterministic(fit, draws = reduced)), c(nrow(reduced), 0L))

  # A single draw is a named vector.
  single <- JAGS_evaluate_deterministic(fit, draws = posterior[5, ])
  expect_identical(unname(single[1, "mu__xREx__s_rho"]), unname(posterior[5, "mu__xREx__s_rho"]))
  expect_error(
    JAGS_evaluate_deterministic(fit, draws = unname(posterior[5, ])),
    "A single draw in 'draws' must be a named numeric vector.",
    fixed = TRUE
  )
})

test_that("linear predictor nodes of row-indexed external SD sources declare the source rows and allocation factors", {

  fit <- .dnode_cached_fit("fit_dnode_row_source")
  posterior <- as.matrix(BayesTools:::.fit_to_posterior(fit))
  nodes <- JAGS_deterministic_nodes(fit)
  node <- nodes[nodes$node == "mu", , drop = FALSE]
  coordinates <- unlist(node$coordinates)
  dependencies <- unlist(node$dependencies)
  factor_dependencies <- c(
    "mu__xRE_ALLOCx_tot__weight[1]", "mu__xRE_ALLOCx_tot__weight[2]",
    "mu__xRE_ALLOCx_tot__include_g_indicator",
    "mu__xRE_ALLOCx_sc__weight[1]", "mu__xRE_ALLOCx_sc__weight[2]"
  )
  expect_true(all(c(paste0("tau[", 1:24, "]"), factor_dependencies) %in% dependencies))
  expect_true(all(dependencies %in% colnames(posterior)))

  # Every declared dependency available: the default call rebuilds 'mu'.
  reduced <- posterior[, setdiff(colnames(posterior), coordinates), drop = FALSE]
  rebuilt <- JAGS_evaluate_deterministic(fit, draws = reduced)
  expect_lte(max(abs(rebuilt[, coordinates] - posterior[, coordinates])), 1e-14)

  # Without the source rows or an allocation weight, the node is unavailable:
  # skipped by default, and an error when requested.
  for(drop in c("^tau\\[", "^mu__xRE_ALLOCx_sc__weight\\[")){
    partial <- reduced[, !grepl(drop, colnames(reduced)), drop = FALSE]
    expect_false(any(coordinates %in% colnames(JAGS_evaluate_deterministic(fit, draws = partial))))
    expect_error(
      JAGS_evaluate_deterministic(fit, draws = partial, nodes = "mu"),
      "Deterministic node 'mu' is unavailable from 'draws'",
      fixed = TRUE
    )
  }

  # The convergence roles read 'tau' as a parent of the linear predictor.
  registry <- BayesTools:::.bt_deterministic_nodes_fit(fit)
  expect_true("tau" %in% BayesTools:::.bt_deterministic_node_parent_bases(registry)[["mu"]])

  # A source with a 'values' function is reconstructed by that function from
  # the inputs it declares: they replace the source rows as dependencies, and
  # the allocation factors still are dependencies.
  values_fit <- .dnode_cached_fit("fit_dnode_row_source_values")
  values_posterior <- as.matrix(BayesTools:::.fit_to_posterior(values_fit))
  values_nodes <- JAGS_deterministic_nodes(values_fit)
  values_dependencies <- unlist(values_nodes$dependencies[values_nodes$node == "mu"])
  expect_false(any(grepl("^tau\\[", values_dependencies)))
  expect_true(all(c("tau_scale", factor_dependencies) %in% values_dependencies))
  values_reduced <- values_posterior[, !grepl("^tau\\[", colnames(values_posterior)) &
                                       !colnames(values_posterior) %in% coordinates, drop = FALSE]
  values_rebuilt <- JAGS_evaluate_deterministic(values_fit, draws = values_reduced)
  expect_lte(max(abs(values_rebuilt[, coordinates] - values_posterior[, coordinates])), 1e-14)

  # Without a declared input, the node is skipped by default and reported
  # unavailable when requested, instead of failing inside the function.
  no_input <- values_reduced[, colnames(values_reduced) != "tau_scale", drop = FALSE]
  expect_false(any(coordinates %in% colnames(JAGS_evaluate_deterministic(values_fit, draws = no_input))))
  expect_error(
    JAGS_evaluate_deterministic(values_fit, draws = no_input, nodes = "mu"),
    "Deterministic node 'mu' is unavailable from 'draws': its dependencies",
    fixed = TRUE
  )
  expect_error(
    JAGS_evaluate_deterministic(values_fit, draws = no_input, nodes = "mu"),
    "'tau_scale'",
    fixed = TRUE
  )
  registry <- BayesTools:::.bt_deterministic_nodes_fit(values_fit)
  expect_true("tau_scale" %in% BayesTools:::.bt_deterministic_node_parent_bases(registry)[["mu"]])
})

test_that("marginal posteriors of formula parameters evaluate the linear predictor node's fixed part", {

  fit <- .dnode_cached_fit("fit_dnode_multiplied")
  samples <- as_mixed_posteriors(fit, parameters = c("mu_intercept", "mu_x", "mu_z", "b_scale"))
  marginal <- marginal_posterior(samples, parameter = "mu_x", formula = ~ x + z, at = list(z = 0.3))
  evaluated <- JAGS_evaluate_formula(
    fit, parameter = "mu", data = data.frame(x = c(-1, 0, 1), z = 0.3)
  )

  # The same terms, multipliers, and floating-point operations as the formula
  # evaluator: (multiplier * coefficient) * data, summed in model order.
  expect_identical(names(marginal), c("-1SD", "0SD", "1SD"))
  for(i in 1:3){
    expect_identical(as.vector(marginal[[i]]), unname(evaluated[i, ]))
  }
})
