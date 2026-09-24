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
