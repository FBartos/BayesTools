skip_if_not_test_profile("unit")

# ============================================================================ #
# TEST FILE: Generated deterministic node registry
# ============================================================================ #
#
# PURPOSE:
#   Every deterministic node that BayesTools generates in a JAGS model is
#   defined once, by a registered family that emits its JAGS syntax and
#   evaluates it in R. The parity tests against the JAGS monitors of small
#   real fits are in test-JAGS-deterministic-nodes-fixture.R.
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

# A fitted object with a gated total-variance allocation of the SDs of two
# random intercepts over synthetic draws (no JAGS run).
.dnode_gated_fit <- function(n = 200L){

  result <- JAGS_formula(
    formula = ~ 1 + random(1 | study, name = "study", covariance = "diag") +
      random(1 | esid, name = "esid", covariance = "diag"),
    parameter = "mu",
    data = data.frame(study = factor(c("a", "a", "b", "b")), esid = factor(1:4)),
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      sd = prior("gamma", list(2, 2)),
      allocation = list(random_variance_allocation(
        name = "split", terms = c(study = "study", esid = "esid"),
        sd = prior("normal", list(0, 1), list(0, Inf)),
        inclusion = list(study = prior("spike", list(location = .5)))
      ))
    )
  )
  share <- seq(0.02, 0.98, length.out = n)
  samples <- cbind(
    mu_intercept = rep(c(-0.1, 0.1), length.out = n),
    mu__xRE_ALLOCx_split__allocation_sd = 0.3 + 0.4 * share,
    "mu__xRE_ALLOCx_split__weight[1]" = share,
    "mu__xRE_ALLOCx_split__weight[2]" = 1 - share,
    mu__xRE_ALLOCx_split__include_study_indicator = rep(c(0, 1), length.out = n)
  )
  fit <- coda::mcmc.list(coda::mcmc(samples))
  class(fit) <- c("BayesTools_fit", class(fit))
  attr(fit, "prior_list") <- result[["prior_list"]]
  attr(fit, "formula_design") <- list(mu = result[["formula_design"]])
  fit <- BayesTools:::.bt_attach_parameter_map(fit)
  fit <- BayesTools:::.bt_attach_draw_geometry(fit)
  BayesTools:::.bt_attach_fit_contract(fit)
}

test_that("JAGS_deterministic_evaluator() resolves the nodes once and reproduces their evaluators", {

  fit <- .dnode_gated_fit()
  draws <- as.matrix(fit[[1L]])
  nodes <- c("mu__xREx__study_intercept", "mu__xREx__esid_intercept")
  registry <- BayesTools:::.bt_deterministic_nodes_fit(fit)
  prior_list <- attr(fit, "prior_list")
  # the node evaluator every other consumer uses
  reference <- function(draws){
    lookup <- BayesTools:::.bt_deterministic_lookup(draws, prior_list)
    do.call(cbind, lapply(registry[nodes], BayesTools:::.bt_deterministic_node_evaluate,
                          lookup = lookup))
  }

  # the registry is built once, when the evaluator is created
  builds <- 0L
  build_registry <- BayesTools:::.bt_deterministic_nodes_fit
  local_mocked_bindings(.bt_deterministic_nodes_fit = function(fit){
    builds <<- builds + 1L
    build_registry(fit)
  })
  evaluate <- JAGS_deterministic_evaluator(fit, nodes = nodes)
  all_rows <- evaluate(draws)
  by_row <- t(vapply(seq_len(nrow(draws)), function(i) evaluate(draws[i, ]), numeric(2)))
  evaluate(draws[1:3, rev(colnames(draws))])
  expect_identical(builds, 1L)

  expect_identical(colnames(all_rows), nodes)
  expect_equal(unname(all_rows), unname(reference(draws)), tolerance = 1e-14)
  expect_identical(all_rows, reference(draws))
  expect_identical(unname(by_row), unname(all_rows))
  expect_identical(all_rows, JAGS_evaluate_deterministic(fit, draws, nodes = nodes))
  # the SD of the gated study intercept is T * gate * sqrt(w[1])
  expect_identical(
    unname(all_rows[, 1L]),
    unname(draws[, "mu__xRE_ALLOCx_split__allocation_sd"] *
             sqrt(draws[, "mu__xRE_ALLOCx_split__weight[1]"]) *
             draws[, "mu__xRE_ALLOCx_split__include_study_indicator"])
  )
  # columns in another order are located again
  reordered <- draws[, rev(colnames(draws))]
  expect_identical(evaluate(reordered), all_rows)

  # the Dirichlet weights from their gamma auxiliaries
  eta_names <- paste0("prior_par_eta_mu__xRE_ALLOCx_split__weight[", 1:2, "]")
  expect_identical(
    BayesTools:::.JAGS_prior_dirichlet_eta_name("mu__xRE_ALLOCx_split__weight"),
    "prior_par_eta_mu__xRE_ALLOCx_split__weight"
  )
  eta <- draws[, setdiff(colnames(draws), paste0("mu__xRE_ALLOCx_split__weight[", 1:2, "]"))]
  eta <- cbind(eta, 3 * draws[, paste0("mu__xRE_ALLOCx_split__weight[", 1:2, "]")])
  colnames(eta)[ncol(eta) - 1:0] <- eta_names
  expect_identical(evaluate(eta), reference(eta))

  # unavailable dependencies and invalid gates stop as in the node evaluator
  expect_error(
    evaluate(draws[, colnames(draws) != "mu__xRE_ALLOCx_split__allocation_sd"]),
    "Deterministic node 'mu__xREx__study_intercept' is unavailable from 'draws'",
    fixed = TRUE
  )
  expect_identical(
    dim(JAGS_deterministic_evaluator(fit)(draws[, colnames(draws) != "mu__xRE_ALLOCx_split__allocation_sd"])),
    c(nrow(draws), 0L)
  )
  expect_error(
    evaluate(draws[, colnames(draws) != "mu__xRE_ALLOCx_split__include_study_indicator"]),
    "Random-effect allocation inclusion samples are missing Bernoulli indicator 'mu__xRE_ALLOCx_split__include_study_indicator'.",
    fixed = TRUE
  )
  invalid <- draws
  invalid[1L, "mu__xRE_ALLOCx_split__include_study_indicator"] <- .5
  expect_error(evaluate(invalid), class = "BayesTools_random_effect_allocation_out_of_support")
  expect_error(JAGS_deterministic_evaluator(fit, nodes = "sd"),
               "'nodes' contains names that are not generated deterministic nodes of 'fit': 'sd'.",
               fixed = TRUE)

  # draws without column names stop, also on the first call of an evaluator
  # (before any column names were checked) and without requested nodes
  unnamed_message <- paste0(
    "'draws' must be a numeric matrix, 'mcmc', or 'mcmc.list' object with ",
    "unique column names, or a named numeric vector."
  )
  expect_error(JAGS_deterministic_evaluator(fit)(unname(draws)), unnamed_message, fixed = TRUE)
  expect_error(JAGS_evaluate_deterministic(fit, unname(draws)), unnamed_message, fixed = TRUE)
  expect_error(evaluate(unname(draws)), unnamed_message, fixed = TRUE)
})

test_that("JAGS_formula() emits the formula syntax from the linear predictor node", {

  sd_prior <- prior("normal", list(0, 1), list(0, Inf))
  x_prior <- prior("normal", list(0, 1))
  attr(x_prior, "multiply_by") <- "b_scale"
  formula_args <- list(
    formula = ~ 1 + x + d + expression(0.25 * z[i]) + us(1 + x | g) + ar1(t | s),
    parameter = "mu",
    data = .dnode_data(),
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      x = x_prior,
      d = prior_factor("mnormal", list(0, 1), contrast = "meandif")
    ),
    prior_random = prior_random(
      g = random_block(sd = sd_prior, cor = prior_lkj(eta = 1)),
      s = random_block(sd = sd_prior, cor = prior("normal", list(0, 0.5)))
    )
  )
  output <- do.call(JAGS_formula, formula_args)
  node <- BayesTools:::.bt_dnode_linear_predictor(output$formula_design)
  emitted <- BayesTools:::.bt_deterministic_node_emit(node)
  expect_identical(
    emitted,
    paste0(
      "for(i in 1:N_mu){\n",
      "  mu[i] = mu_intercept + b_scale * mu_x * mu_data_x[i] + ",
      "inprod(mu_d, mu_data_d[i,]) + 0.25 * z[i] + mu__xREx__g[i] + mu__xREx__s[i]\n",
      "}\n"
    )
  )
  expect_identical(substr(output$formula_syntax, 1L, nchar(emitted)), emitted)

  # The node is the only definition of the formula syntax: a different node
  # emitter gives a different model, with the random-effect syntax unchanged.
  local_mocked_bindings(
    .bt_dnode_linear_predictor_emit = function(node) paste0("<", node$node, " node>\n"),
    .package = "BayesTools"
  )
  mocked <- do.call(JAGS_formula, formula_args)
  expect_identical(
    mocked$formula_syntax,
    sub(emitted, "<mu node>\n", output$formula_syntax, fixed = TRUE)
  )
})

test_that("selection_backend_spec() emits the weights from the omega node", {

  priors <- list(
    two_sided = prior_weightfunction("two-sided", c(0.05, 0.1), wf_cumulative(c(1, 2, 1))),
    phacking_only = prior_bias(phacking = prior_phacking()),
    mixture = prior_mixture(list(
      prior_none(),
      prior_bias(phacking = prior_phacking()),
      prior_weightfunction("one-sided", c(0.025, 0.05), wf_cumulative(c(1, 1, 1))),
      prior_weightfunction("two-sided", c(0.05), wf_independent(prior("normal", list(0, 1)), scale = "log_omega"))
    ), is_null = c(TRUE, FALSE, FALSE, FALSE))
  )
  for(name in c("two_sided", "mixture")){
    spec <- selection_backend_spec(priors[[name]], include_init = FALSE)
    node <- BayesTools:::.bt_dnode_omega("omega", priors[[name]])
    expect_identical(node$coordinates, spec$step$coefficient_ids)
    code <- paste(spec$prior_code, spec$transform_code, sep = "\n")
    for(piece in BayesTools:::.bt_deterministic_node_emit(node)){
      expect_true(grepl(piece, code, fixed = TRUE))
    }
  }
  # A p-hacking-only prior has unit weights and no registered node.
  expect_null(BayesTools:::.bt_dnode_omega("omega", priors$phacking_only))
  expect_true(grepl(
    "omega[1] <- 1",
    selection_backend_spec(priors$phacking_only, include_init = FALSE)$prior_code,
    fixed = TRUE
  ))

  # Renamed weights on a finer global grid are the node of that specification.
  one_sided <- prior_weightfunction("one-sided", c(0.025, 0.5), wf_cumulative(c(1, 1, 1)))
  breaks <- c(0, 0.01, 0.025, 0.5, 1)
  custom <- selection_backend_spec(
    one_sided, names = list(omega = "w"), global_breaks = breaks, include_init = FALSE
  )
  custom_node <- BayesTools:::.bt_dnode_omega(
    "w", one_sided,
    spec = BayesTools:::.bt_dnode_omega_spec(one_sided, global_breaks = breaks, name = "w")
  )
  expect_identical(custom_node$coordinates, custom$step$coefficient_ids)
  expect_identical(custom_node$coordinates, paste0("w[", 1:4, "]"))
  for(piece in BayesTools:::.bt_deterministic_node_emit(custom_node)){
    expect_true(grepl(piece, custom$prior_code, fixed = TRUE))
  }

  # The node specification is the only definition of the weight syntax.
  local_mocked_bindings(
    .bt_dnode_omega_emit_branch = function(spec, k) paste0("<", spec$name, " branch ", k, ">\n"),
    .bt_dnode_omega_emit_composition = function(spec){
      if(spec$uses_indicator) paste0("<", spec$name, " composition>") else character()
    },
    .package = "BayesTools"
  )
  mocked <- selection_backend_spec(priors$mixture, include_init = FALSE)
  for(k in 1:4){
    expect_true(grepl(paste0("<omega branch ", k, ">"), mocked$prior_code, fixed = TRUE))
  }
  expect_false(grepl("omega_component_", mocked$prior_code, fixed = TRUE))
  expect_true(grepl("<omega composition>", mocked$transform_code, fixed = TRUE))
  expect_true(grepl(
    "<w branch 1>",
    selection_backend_spec(one_sided, names = list(omega = "w"), include_init = FALSE)$prior_code,
    fixed = TRUE
  ))
})

test_that("parameter_source() values functions declare and receive their inputs", {

  values <- function(parameters, data, n_rows) rep(parameters[["a"]], n_rows)
  declared <- parameter_source("tau", shape = "row", values = values, inputs = c("a", "b[2]"))
  expect_identical(declared$inputs, c("a", "b[2]"))
  # a values function must declare its inputs
  undeclared <- tryCatch(
    parameter_source("tau", shape = "row", values = values),
    error = function(e) e
  )
  expect_s3_class(undeclared, c("BayesTools_missing_source_inputs", "BayesTools_parameter_source"))
  expect_identical(
    conditionMessage(undeclared),
    paste0(
      "'inputs' must declare the posterior coordinates that the 'values' ",
      "function reads from 'parameters' (character() for a function of ",
      "'data' alone)."
    )
  )
  expect_false("inputs" %in% names(parameter_source("tau", shape = "row")))
  expect_identical(parameter_source("tau", shape = "row", values = values, inputs = character())$inputs, character())
  expect_error(
    parameter_source("tau", shape = "row", inputs = "a"),
    "'inputs' is supported only for parameter sources with a 'values' function.",
    fixed = TRUE
  )
  expect_error(
    parameter_source("tau", shape = "row", values = values, inputs = c("a", NA)),
    "'inputs' must be a character vector of posterior coordinate names.",
    fixed = TRUE
  )
  broken <- declared
  broken$values <- NULL
  expect_error(
    BayesTools:::.bt_check_parameter_source(broken),
    "'source$inputs' is supported only for parameter sources with a 'values' function.",
    fixed = TRUE
  )
  expect_true(any(grepl("inputs: a, b[2]", capture.output(print(declared)), fixed = TRUE)))

  # The function receives only its declared inputs; reading another parameter
  # or a draw without a declared input stops.
  posterior <- matrix(c(2, 3, 5, 7), nrow = 1, dimnames = list(NULL, c("a", "b[1]", "b[2]", "c")))
  seen <- NULL
  recorder <- parameter_source("tau", shape = "row", inputs = c("a", "b[2]"),
    values = function(parameters, data, n_rows){
      seen <<- names(parameters)
      rep(parameters[["a"]] * parameters[["b[2]"]], n_rows)
    })
  expect_identical(
    unname(BayesTools:::.bt_parameter_source_value_draws(recorder, n_rows = 2, posterior = posterior)[1, ]),
    c(10, 10)
  )
  expect_identical(seen, c("a", "b[2]"))
  reader <- parameter_source("tau", shape = "row", inputs = "a",
    values = function(parameters, data, n_rows) rep(parameters$c, n_rows))
  expect_error(
    BayesTools:::.bt_parameter_source_value_draws(reader, n_rows = 2, posterior = posterior),
    "reads 'c', which is not among its declared 'inputs'.",
    fixed = TRUE
  )
  expect_error(
    BayesTools:::.bt_parameter_source_value_draws(recorder, n_rows = 2, posterior = posterior[, c("a", "c"), drop = FALSE]),
    "Parameter source 'tau[row]' is missing its declared input(s) 'b[2]'.",
    fixed = TRUE
  )
  # A function of the data alone receives no parameters, and a source object
  # without declared inputs is refused where it is used.
  seen <- NULL
  BayesTools:::.bt_parameter_source_value_draws(
    parameter_source("tau", shape = "row", inputs = character(), values = function(parameters, data, n_rows){
      seen <<- names(parameters)
      rep(1, n_rows)
    }),
    n_rows = 2, posterior = posterior
  )
  expect_identical(seen, character())
  undeclared_object <- declared
  undeclared_object$inputs <- NULL
  expect_error(
    BayesTools:::.bt_parameter_source_value_draws(undeclared_object, n_rows = 2, posterior = posterior),
    "'source$inputs' must declare the posterior coordinates",
    fixed = TRUE,
    class = "BayesTools_missing_source_inputs"
  )
})

test_that("JAGS_marglik_parameters() builds the weight node of a prior once across draws", {

  cumulative <- list(omega = prior_weightfunction("one-sided", c(0.025, 0.5), wf_cumulative(c(1, 2, 3))))
  fixed <- list(
    a = list(omega = prior_weightfunction("one-sided", c(0.025, 0.5), wf_fixed(c(1, 0.5, 0.2)))),
    b = list(omega = prior_weightfunction("one-sided", c(0.025, 0.5), wf_fixed(c(1, 0.4, 0.2))))
  )
  cache <- BayesTools:::.bt_dnode_omega_cache
  cache$entries <- NULL
  build <- BayesTools:::.bt_dnode_omega
  built <- 0L
  local_mocked_bindings(
    .bt_dnode_omega = function(...){
      built <<- built + 1L
      build(...)
    },
    .package = "BayesTools"
  )

  # 20 draws of one prior list build its node once.
  omega <- lapply(1:20, function(i){
    JAGS_marglik_parameters(c("eta[1]" = i, "eta[2]" = 2, "eta[3]" = 3), cumulative)$omega
  })
  expect_identical(built, 1L)
  expect_identical(omega[[5L]], c(1, 0.5, 0.3))

  # Priors that differ in their weights are different nodes.
  for(i in 1:3){
    expect_identical(JAGS_marglik_parameters(numeric(), fixed$a)$omega, c(1, 0.5, 0.2))
    expect_identical(JAGS_marglik_parameters(numeric(), fixed$b)$omega, c(1, 0.4, 0.2))
  }
  expect_identical(built, 3L)
})

test_that("prior draws carry the auxiliary nodes of mixture and spike-and-slab priors", {

  formula_result <- JAGS_formula(
    ~ f, "mu", data.frame(f = factor(c("a", "b", "c", "a"))),
    prior_list = list(
      intercept = prior("normal", list(0, 1)),
      f = prior_spike_and_slab(prior_factor("mnormal", list(0, 1), contrast = "meandif"))
    )
  )
  prior_list <- c(
    formula_result$prior_list,
    list(
      m = prior_spike_and_slab(prior("normal", list(0, 1)),
                               prior_inclusion = prior("beta", list(1, 1))),
      t = prior_mixture(list(prior("spike", list(0)), prior("normal", list(0, 1))),
                        is_null = c(TRUE, FALSE))
    )
  )
  columns <- c("mu_intercept", "mu_f[1]", "mu_f[2]", "mu_f_indicator", "mu_f_inclusion",
               "mu_f_variable[1]", "mu_f_variable[2]",
               "m", "m_indicator", "m_inclusion", "m_variable", "t", "t_indicator")
  posterior <- matrix(.5, nrow = 4, ncol = length(columns), dimnames = list(NULL, columns))
  fit <- coda::mcmc(posterior)
  class(fit) <- c("mcmc", "BayesTools_fit")
  attr(fit, "prior_list") <- prior_list
  attr(fit, "formula_design") <- list(mu = formula_result$formula_design)
  fit <- attach_test_parameter_map(fit)

  draws <- transform_prior_samples(fit, n_samples = 2000, seed = 11)
  expect_identical(colnames(draws), columns)
  # the values follow the fitted node definitions
  expect_identical(unname(draws[, "m"]), unname(draws[, "m_variable"] * draws[, "m_indicator"]))
  expect_identical(
    unname(draws[, c("mu_f[1]", "mu_f[2]")]),
    unname(draws[, c("mu_f_variable[1]", "mu_f_variable[2]")] * draws[, "mu_f_indicator"])
  )
  expect_true(all(draws[, "t"][draws[, "t_indicator"] == 1] == 0))
  expect_true(all(draws[, "t_indicator"] %in% 1:2))
  expect_true(all(draws[, "mu_f_inclusion"] == .5))

  # the value columns keep rng()'s random-number stream (as before the
  # auxiliary columns were added)
  generated <- BayesTools:::.generate_prior_sample_matrix(
    prior_list[c("m", "t")], n_samples = 500, seed = 3
  )
  set.seed(3)
  expected_m <- rng(prior_list$m, 500)
  expected_t <- rng(prior_list$t, 500)
  expect_identical(unname(generated[, "m"]), as.numeric(expected_m))
  expect_identical(unname(generated[, "m_indicator"]), as.numeric(attr(expected_m, "inclusion")))
  expect_identical(unname(generated[, "t"]), as.numeric(expected_t))
  expect_identical(unname(generated[, "t_indicator"]), as.numeric(attr(expected_t, "components")))

  # catalog quantities of the auxiliary nodes resolve on the prior draws
  catalog <- parameter_catalog(fit)
  indicator <- parameter_draws(
    fit,
    parameter_catalog_resolve(catalog, "m_indicator"),
    model_samples = draws
  )
  expect_identical(as.numeric(as.matrix(indicator)), unname(draws[, "m_indicator"]))
})

test_that("prior draws carry the fitted nodes of ordered-prior totals", {

  data <- data.frame(
    f = ordered(rep(c("low", "mid", "high"), 2), levels = c("low", "mid", "high")),
    g = factor(rep(c("a", "b", "c"), each = 2))
  )
  prior_columns <- function(prior_list){
    unlist(lapply(names(prior_list), function(name){
      prior <- prior_list[[name]]
      if(BayesTools:::.bt_prior_is_factor_family(prior)){
        BayesTools:::.JAGS_prior_factor_names(name, prior)
      }else{
        name
      }
    }), use.names = FALSE)
  }
  prior_fit <- function(formula, prior_list, total_columns){
    formula_result <- JAGS_formula(formula, "mu", data, prior_list)
    columns <- c(prior_columns(formula_result$prior_list), total_columns)
    fit <- coda::mcmc(matrix(.5, nrow = 4, ncol = length(columns), dimnames = list(NULL, columns)))
    class(fit) <- c("mcmc", "BayesTools_fit")
    attr(fit, "prior_list") <- formula_result$prior_list
    attr(fit, "formula_design") <- list(mu = formula_result$formula_design)
    attach_test_parameter_map(fit)
  }
  totals <- list(
    spike = prior_spike_and_slab(prior("normal", list(0, 1)),
                                 prior_inclusion = prior("beta", list(1, 1))),
    mixture = prior_mixture(list(prior("spike", list(0)), prior("normal", list(0, 1))),
                            is_null = c(TRUE, FALSE))
  )
  total_nodes <- list(
    spike   = paste0("mu_f_ordered_total", c("", "_indicator", "_inclusion", "_variable")),
    mixture = paste0("mu_f_ordered_total", c("", "_indicator"))
  )
  n <- 2000
  for(name in names(totals)){
    fit <- prior_fit(~ f, list(intercept = prior("normal", list(0, 1)),
                               f = prior_ordered(totals[[name]])), total_nodes[[name]])
    prior_list <- attr(fit, "prior_list")
    draws <- transform_prior_samples(fit, n_samples = n, seed = 11)
    expect_identical(colnames(draws),
                     c("mu_intercept", "mu_f[1]", "mu_f[2]", total_nodes[[name]]), info = name)

    # the existing columns keep rng()'s random-number stream, and the total
    # and its indicator are the draws of the total's own rng() stream
    set.seed(11)
    intercept <- rng(prior_list$mu_intercept, n)
    total <- rng(totals[[name]], n)
    set.seed(11)
    rng(prior_list$mu_intercept, n)
    coefficients <- rng(prior_list$mu_f, n, transform_factor_samples = FALSE)
    expect_identical(unname(draws[, "mu_intercept"]), as.numeric(intercept), info = name)
    expect_identical(unname(draws[, c("mu_f[1]", "mu_f[2]")]), unname(unclass(coefficients)[, 1:2]),
                     info = name)
    expect_identical(unname(draws[, "mu_f_ordered_total"]), as.numeric(total), info = name)
    indicator <- if(name == "spike") attr(total, "inclusion") else attr(total, "components")
    expect_identical(unname(draws[, "mu_f_ordered_total_indicator"]), as.numeric(indicator),
                     info = name)
    # the coefficients allocate the total
    expect_equal(unname(draws[, "mu_f[1]"] + draws[, "mu_f[2]"]),
                 unname(draws[, "mu_f_ordered_total"]), tolerance = 1e-12, info = name)
    if(name == "spike"){
      expect_identical(
        unname(draws[, "mu_f_ordered_total"]),
        unname(draws[, "mu_f_ordered_total_variable"] * draws[, "mu_f_ordered_total_indicator"])
      )
      expect_true(all(draws[, "mu_f_ordered_total_inclusion"] > 0 &
                        draws[, "mu_f_ordered_total_inclusion"] < 1))
    }else{
      expect_true(all(draws[, "mu_f_ordered_total"][draws[, "mu_f_ordered_total_indicator"] == 1] == 0))
    }

    # catalog quantities of the total's nodes resolve on the prior draws
    catalog <- parameter_catalog(fit)
    for(node in total_nodes[[name]]){
      values <- parameter_draws(fit, parameter_catalog_resolve(catalog, node), model_samples = draws)
      expect_identical(as.numeric(as.matrix(values)), unname(draws[, node]), info = node)
    }
  }

  # an interaction with two slices: the total has one node per slice, and a
  # spike-and-slab total one indicator and inclusion probability shared by
  # the slices, as in the fitted model
  simple <- prior_fit(~ g * f, list(
    intercept = prior("normal", list(0, 1)),
    g = prior_factor("normal", list(0, 1), contrast = "treatment"),
    f = prior_ordered(prior("normal", list(0, 1))),
    "g:f" = prior_ordered(prior("normal", list(0, 1)))
  ), c("mu_f_ordered_total", "mu_g__xXx__f_ordered_total[1]", "mu_g__xXx__f_ordered_total[2]"))
  draws <- transform_prior_samples(simple, n_samples = n, seed = 12)
  expect_true(all(c("mu_g__xXx__f_ordered_total[1]", "mu_g__xXx__f_ordered_total[2]") %in% colnames(draws)))
  interaction <- attr(simple, "prior_list")$mu_g__xXx__f
  coordinates <- BayesTools:::.JAGS_prior_factor_names("mu_g__xXx__f", interaction)
  slices <- unlist(attr(interaction, "ordered_metadata")$slice_index)
  # the coefficients of each slice allocate its total
  for(slice in 1:2){
    expect_equal(
      unname(rowSums(draws[, coordinates[slices == slice], drop = FALSE])),
      unname(draws[, paste0("mu_g__xXx__f_ordered_total[", slice, "]")]),
      tolerance = 1e-12
    )
  }

  spike <- prior_fit(~ g * f, list(
    intercept = prior("normal", list(0, 1)),
    g = prior_factor("normal", list(0, 1), contrast = "treatment"),
    f = prior_ordered(prior("normal", list(0, 1))),
    "g:f" = prior_ordered(totals$spike)
  ), c("mu_f_ordered_total", paste0("mu_g__xXx__f_ordered_total",
         c("_indicator", "_inclusion", "[1]", "[2]", "_variable[1]", "_variable[2]"))))
  draws <- transform_prior_samples(spike, n_samples = n, seed = 12)
  slice_nodes <- paste0("mu_g__xXx__f_ordered_total",
                        c("_indicator", "_inclusion", "[1]", "[2]", "_variable[1]", "_variable[2]"))
  expect_true(all(slice_nodes %in% colnames(draws)))
  for(slice in 1:2){
    expect_identical(
      unname(draws[, paste0("mu_g__xXx__f_ordered_total[", slice, "]")]),
      unname(draws[, paste0("mu_g__xXx__f_ordered_total_variable[", slice, "]")] *
               draws[, "mu_g__xXx__f_ordered_total_indicator"])
    )
  }
  interaction <- attr(spike, "prior_list")$mu_g__xXx__f
  coordinates <- BayesTools:::.JAGS_prior_factor_names("mu_g__xXx__f", interaction)
  excluded <- draws[, "mu_g__xXx__f_ordered_total_indicator"] == 0
  expect_true(any(excluded) && any(!excluded))
  expect_true(all(draws[excluded, coordinates] == 0))
  expect_true(all(draws[!excluded, coordinates] != 0))
  # all six nodes are catalog quantities that resolve on the prior draws
  catalog <- parameter_catalog(spike)
  expect_true(all(slice_nodes %in% catalog$quantities$canonical_name))
  for(node in slice_nodes){
    values <- parameter_draws(spike, parameter_catalog_resolve(catalog, node), model_samples = draws)
    expect_identical(as.numeric(as.matrix(values)), unname(draws[, node]), info = node)
  }
})

