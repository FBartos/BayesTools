skip_if_not_test_profile("unit")

# ============================================================================ #
# TEST FILE: Generated deterministic node registry
# ============================================================================ #
#
# PURPOSE:
#   Every deterministic node that BayesTools generates in a JAGS model is
#   defined once, by a registered family that emits its JAGS syntax and
#   evaluates it in R. The parity tests remove the monitored nodes from the
#   posterior draws of small real fits and rebuild them with
#   JAGS_evaluate_deterministic(): the R evaluators must reproduce the JAGS
#   monitors (bit-identical where R and JAGS perform the same floating-point
#   operations, otherwise within 1e-14).
#
# TAGS: @unit, @JAGS, @deterministic-nodes
# ============================================================================ #

.dnode_data <- function(){
  set.seed(11)
  n <- 24L
  data.frame(
    x = stats::rnorm(n, 5, 3),
    z = stats::rnorm(n, -2, 0.5),
    t = factor(rep(c("t1", "t2", "t3", "t4"), 6), levels = c("t1", "t2", "t3", "t4")),
    time = rep(c(0, 1, 2.5, 4), 6),
    g = factor(rep(sprintf("g%d", 1:6), each = 4)),
    s = factor(rep(sprintf("s%d", 1:4), each = 6)),
    d = factor(rep(c("a", "b", "c"), 8)),
    p = factor(rep(sprintf("p%d", 1:3), each = 8))
  )
}

.dnode_fit <- function(formula, prior_list, prior_random = NULL,
                       formula_scale = NULL, extra_prior = NULL,
                       add_parameters = NULL, seed = 1L){

  skip_if_not_installed("rjags")
  skip_if_not_installed("runjags")
  data <- .dnode_data()
  set.seed(seed)
  y <- stats::rnorm(nrow(data), 0.1 * data$x, 1)
  suppressWarnings(JAGS_fit(
    model_syntax = "model{\n  for(i in 1:N_mu){\n    y[i] ~ dnorm(mu[i], 1)\n  }\n}",
    data = list(y = y),
    prior_list = extra_prior,
    formula_list = list(mu = formula),
    formula_data_list = list(mu = data),
    formula_prior_list = list(mu = prior_list),
    formula_scale_list = if(!is.null(formula_scale)) list(mu = formula_scale),
    formula_random_prior_list = if(!is.null(prior_random)) list(mu = prior_random),
    add_parameters = add_parameters,
    chains = 1, adapt = 50, burnin = 50, sample = 100, seed = seed, silent = TRUE
  ))
}

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
  fit <- .dnode_fit(
    ~ 1 + x + us(1 + x | g) + random(1 | d, name = "d", covariance = "diag"),
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1))
    ),
    prior_random = prior_random(
      random_variance_allocation(
        name = "tot", terms = c(g = "g", d = "d"), scale = "total_variance",
        sd = prior("gamma", list(2, 2)),
        weights = prior("dirichlet", list(alpha = c(1.5, 2.5))),
        inclusion = list(d = prior("beta", list(3, 2)))
      ),
      random_variance_allocation(
        name = "gc", parent = allocation_ref("tot", "g"), terms = "g",
        target = "sd_component", scale = "mean_variance",
        weights = prior("dirichlet", list(alpha = c(1, 2)))
      )
    ),
    formula_scale = list(x = TRUE)
  )
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
  external <- .dnode_fit(
    ~ 1 + random(1 | g, name = "g", covariance = "diag") +
      random(1 | d, name = "d", covariance = "diag"),
    prior_list = list(intercept = prior("normal", list(0, 1))),
    extra_prior = list(tau = prior("normal", list(0, 1), list(0, Inf))),
    prior_random = prior_random(random_variance_allocation(
      name = "tot", terms = c(g = "g", d = "d"),
      sd_source = random_sd_source("tau"),
      weights = prior("dirichlet", list(alpha = c(2, 2)))
    )),
    seed = 2L
  )
  external_result <- .dnode_rebuild(external, "random_sd")
  expect_identical(
    unname(unlist(external_result$nodes$dependencies[1L])),
    c("tau", "mu__xRE_ALLOCx_tot__weight[1]", "mu__xRE_ALLOCx_tot__weight[2]")
  )
  expect_identical(unname(external_result$rebuilt), unname(external_result$monitored))
})

test_that("scalar correlation nodes reproduce the JAGS monitors of every structure", {

  sd_prior <- prior("normal", list(0, 1), list(0, Inf))
  fit <- .dnode_fit(
    ~ 1 + cs(t | g) + ar1(t | s) + car(time | d) + hcs(t | p),
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      g = random_block(sd = sd_prior, covariance = random_covariance(
        cor = prior("normal", list(0, 1.5)), cor_scale = "logit")),
      s = random_block(sd = sd_prior, cor = prior("normal", list(0, 0.5))),
      d = random_block(sd = sd_prior, covariance = random_covariance(
        cor = prior("normal", list(0, 1)), cor_scale = "logit")),
      p = random_block(sd = sd_prior, cor = prior("spike", list(location = 0.4)),
                       monitor = random_monitor(correlation = TRUE))
    )
  )
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

  sd_prior <- prior("normal", list(0, 1), list(0, Inf))
  fit <- .dnode_fit(
    ~ 1 + x + z + us(1 + x + z | g) + us(1 + x | s),
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior("normal", list(0, 1)),
      z = prior("normal", list(0, 1))
    ),
    prior_random = prior_random(
      g = random_block(sd = sd_prior, cor = prior_lkj(eta = 1.5, include_primitives = TRUE),
                       monitor = random_monitor(lkj_primitives = TRUE)),
      s = random_block(sd = sd_prior, cor = prior_lkj(eta = 2))
    ),
    formula_scale = list(x = TRUE)
  )
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

.dnode_prior_fit <- function(prior_list, add_parameters = NULL, seed = 1L){

  skip_if_not_installed("rjags")
  skip_if_not_installed("runjags")
  set.seed(seed)
  suppressWarnings(JAGS_fit(
    model_syntax = "model{\n  for(i in 1:N){\n    x[i] ~ dnorm(m, 1)\n  }\n}",
    data = list(x = stats::rnorm(10), N = 10L),
    prior_list = c(list(m = prior("normal", list(0, 1))), prior_list),
    add_parameters = add_parameters,
    chains = 1, adapt = 50, burnin = 50, sample = 100, seed = seed, silent = TRUE
  ))
}

test_that("publication-weight nodes reproduce the JAGS monitors of every weight function", {

  fits <- list(
    cumulative_two_sided = .dnode_prior_fit(list(omega = prior_weightfunction(
      "two-sided", c(0.05, 0.1), wf_cumulative(c(1, 2, 1))))),
    binary = .dnode_prior_fit(list(omega = prior_weightfunction(
      "two-sided", c(0.05), wf_cumulative(c(1, 1))))),
    log_independent = .dnode_prior_fit(list(omega = prior_weightfunction(
      "one-sided", c(0.025, 0.5), wf_independent(prior("normal", list(0, 1)), scale = "log_omega")))),
    fixed = .dnode_prior_fit(list(omega = prior_weightfunction(
      "one-sided", c(0.025, 0.5), wf_fixed(c(1, 1/3, 0.2))))),
    bias_mixture = .dnode_prior_fit(
      list(bias = prior_mixture(list(
        prior_PET("normal", list(0, 1)),
        prior_weightfunction("one-sided", c(0.025, 0.05), wf_cumulative(c(1, 1, 1))),
        prior_weightfunction("two-sided", c(0.05), wf_cumulative(c(1, 1))),
        prior_weightfunction("two-sided", c(0.05, 0.1), wf_independent(prior("normal", list(0, 1)), scale = "log_omega"))
      ), is_null = c(FALSE, FALSE, FALSE, FALSE))),
      add_parameters = c("eta_component_2", "omega_ratio_component_3", "log_omega_component_4")
    )
  )

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

  fit <- .dnode_prior_fit(
    list(
      a = prior_spike_and_slab(prior("normal", list(0, 1)), prior_inclusion = prior("beta", list(1, 1))),
      b = prior_mixture(list(prior("normal", list(0, 1), list(0, Inf)), prior("spike", list(0)),
                             prior("normal", list(2, 0.5))), is_null = c(FALSE, TRUE, FALSE)),
      bias = prior_mixture(list(
        prior_PET("normal", list(0, 1)),
        prior_PEESE("normal", list(0, 1)),
        prior_weightfunction("one-sided", c(0.025, 0.05), wf_cumulative(c(1, 1, 1)))
      ), is_null = c(FALSE, FALSE, FALSE))
    ),
    add_parameters = c("b_component_1", "b_component_3", "PET_1", "PEESE_1", "eta_component_3"),
    seed = 3L
  )
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
  factor_fit <- .dnode_fit(
    ~ 1 + x + t,
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = prior_spike_and_slab(prior("normal", list(0, 1)), prior_inclusion = prior("beta", list(1, 1))),
      t = prior_mixture(list(prior_factor("spike", list(0), contrast = "treatment"),
                             prior_factor("normal", list(0, 1), contrast = "treatment")),
                        is_null = c(TRUE, FALSE))
    ),
    add_parameters = "mu_t_component_2"
  )
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

test_that("JAGS_evaluate_deterministic() validates its node selection", {

  sd_prior <- prior("normal", list(0, 1), list(0, Inf))
  fit <- .dnode_fit(
    ~ 1 + ar1(t | s),
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      s = random_block(sd = sd_prior, cor = prior("normal", list(0, 0.5)))
    )
  )
  nodes <- JAGS_deterministic_nodes(fit)
  expect_identical(
    names(nodes),
    c("node", "family", "parameter", "block", "coordinates", "dependencies", "monitored")
  )
  expect_identical(nodes$parameter, "mu")
  expect_identical(nodes$block, "s")
  expect_true(nodes$monitored)

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
