skip_if_not_test_profile("unit")

# Sampler declarations of JAGS_fit() ('sampler_gates', 'sampler_blocks'): the
# expansion of priors into marker nodes, their eligibility and guard rules, and
# the classed errors. No model is compiled here; the samplers themselves are
# certified in the 'fit-assertions' case.

.declare <- function(prior_list, model_syntax = "model{\n  y ~ dnorm(0, 1)\n}", data = list(y = 0),
                     sampler_gates = FALSE, sampler_blocks = NULL, add_parameters = NULL){
  BayesTools:::.bt_sampler_declarations(
    model_text = JAGS_add_priors(model_syntax, prior_list),
    data = data, prior_list = prior_list, add_parameters = add_parameters,
    sampler_gates = sampler_gates, sampler_blocks = sampler_blocks
  )
}

.gates <- function(declarations){
  gates <- declarations$declarations$gates
  stats::setNames(gates$active, gates$node)
}

.wf_cumulative2 <- function() prior_weightfunction("one-sided", c(0.05), wf_cumulative(c(1, 1)))
.wf_cumulative3 <- function() prior_weightfunction("one-sided", c(0.05, 0.5), wf_cumulative(c(1, 1, 1)))

test_that("default arguments leave the model text and data unchanged", {

  prior_list <- list(
    mu = prior("normal", list(0, 1)),
    bias = prior_mixture(list(prior_none(), prior_PET("normal", list(0, 1))), is_null = c(TRUE, FALSE))
  )
  data <- list(y = 1)
  text <- JAGS_add_priors("model{\n  y ~ dnorm(mu + PET, 1)\n}", prior_list)
  out <- BayesTools:::.bt_sampler_declarations(text, data, prior_list, NULL)
  expect_identical(out$model, text)
  expect_identical(out$data, data)
  expect_null(out$declarations)
  expect_false(out$uses_module)
})

test_that("JAGS_fit() hands the driver the same model text and data with default arguments", {

  skip_if_not_installed("rjags")
  # the backend call is captured before any sampling
  capture <- function(...){
    call <- NULL
    modules <- NULL
    testthat::local_mocked_bindings(
      .JAGS_run_backend = function(fun, model_call, parallel){
        call <<- model_call
        stop(structure(class = c("captured_call", "condition"), list(message = "", call = NULL)))
      },
      .JAGS_load_modules = function(jags_modules, ...){
        modules <<- jags_modules
        invisible(TRUE)
      },
      .package = "BayesTools"
    )
    tryCatch(JAGS_fit(...), captured_call = function(e) NULL)
    list(call = call, modules = modules)
  }

  set.seed(1)
  formula_data <- data.frame(x = stats::rnorm(10), f = factor(rep(c("a", "b"), 5)))
  model_set <- list(
    list(syntax = "model{\n  y ~ dnorm(mu, 1)\n}", data = list(y = 0.3),
         prior_list = list(mu = prior("normal", list(0, 1)))),
    list(syntax = "model{\n  y ~ dnorm(mu + PET + PEESE, 1)\n}", data = list(y = 0.3),
         prior_list = list(
           mu = prior_spike_and_slab(prior("normal", list(0, 1))),
           bias = prior_mixture(list(prior_none(), prior_PET("normal", list(0, 1)), prior_PEESE("normal", list(0, 1))),
                                is_null = c(TRUE, FALSE, FALSE))
         )),
    list(syntax = "model{\n  y ~ dnorm(mu * omega[1], 1)\n}", data = list(y = 0.3),
         prior_list = list(
           mu = prior_mixture(list(prior("normal", list(0, 1)), prior("point", list(0)))),
           bias = prior_mixture(list(prior_none(), .wf_cumulative3(), prior_phacking(form = "linear")),
                                is_null = c(TRUE, FALSE, FALSE))
         )),
    list(syntax = "model{\n  for(i in 1:N_mu){ y[i] ~ dnorm(mu[i], 1) }\n}", data = list(y = stats::rnorm(10)),
         prior_list = NULL, formula_list = list(mu = ~ x + f), formula_data_list = list(mu = formula_data),
         formula_prior_list = list(mu = list(
           intercept = prior("normal", list(0, 1)),
           x = prior("normal", list(0, 1)),
           f = prior_factor("normal", list(0, 1), contrast = "treatment")
         )))
  )
  for(model in model_set){
    arguments <- list(
      model_syntax = model$syntax, data = model$data, prior_list = model$prior_list,
      formula_list = model$formula_list, formula_data_list = model$formula_data_list,
      formula_prior_list = model$formula_prior_list,
      chains = 1, adapt = 50, burnin = 50, sample = 100, seed = 1
    )
    default <- do.call(capture, arguments)
    explicit <- do.call(capture, c(arguments, list(sampler_gates = FALSE, sampler_blocks = NULL)))
    expect_identical(default$call$model, explicit$call$model)
    expect_identical(default$call$data, explicit$call$data)
    if(is.null(model$formula_list)){
      prior_list <- BayesTools:::.complete_factor_metadata_prior_list(model$prior_list)
      expect_identical(default$call$model, JAGS_add_priors(model$syntax, prior_list))
      expect_identical(default$call$data, model$data)
    }
    expect_false(any(grepl("BayesTools_(gate|block)_", default$call$model)))
    expect_false(any(startsWith(names(default$call$data), "BayesTools_")))
  }

  # with declarations the markers and their data reach the driver and the module is required
  declared <- do.call(capture, list(
    model_syntax = "model{\n  y ~ dnorm(mu + PET, 1)\n}", data = list(y = 0.3),
    prior_list = list(
      mu = prior("normal", list(0, 1)), tau = prior("normal", list(0, 1), list(0, Inf)),
      bias = prior_mixture(list(prior_none(), prior_PET("normal", list(0, 1))), is_null = c(TRUE, FALSE))
    ),
    chains = 1, adapt = 50, burnin = 50, sample = 100, seed = 1,
    sampler_gates = TRUE, sampler_blocks = list(c("mu", "tau"))
  ))
  expect_true(grepl("BayesTools_gate_1 ~ dbt_gate(PET_1, equals(bias_indicator, 2))", declared$call$model, fixed = TRUE))
  expect_true(grepl("BayesTools_block_1 ~ dbt_block(mu, 1)", declared$call$model, fixed = TRUE))
  expect_true(grepl("BayesTools_block_2 ~ dbt_block(tau, 1)", declared$call$model, fixed = TRUE))
  expect_identical(declared$call$data[c("BayesTools_gate_1", "BayesTools_block_1", "BayesTools_block_2")],
                   list(BayesTools_gate_1 = 0, BayesTools_block_1 = 0, BayesTools_block_2 = 0))
  expect_true("BayesTools" %in% declared$modules)
})

test_that("gates are declared for PET and PEESE in every branch position", {

  PET <- prior_PET("normal", list(0, 1))
  PEESE <- prior_PEESE("normal", list(0, 1))
  orders <- list(
    list(prior_none(), PET, PEESE),
    list(PEESE, prior_none(), PET),
    list(PET, PEESE, prior_none()),
    list(PEESE, PET)
  )
  for(branches in orders){
    is_null <- vapply(branches, is.prior.none, logical(1))
    prior_list <- list(bias = prior_mixture(branches, is_null = is_null))
    out <- .declare(prior_list, "model{\n  y ~ dnorm(PET + PEESE, 1)\n}", sampler_gates = TRUE)
    position <- function(test) which(vapply(branches, test, logical(1)))
    expect_identical(.gates(out), c(
      PET_1 = paste0("equals(bias_indicator, ", position(is.prior.PET), ")"),
      PEESE_1 = paste0("equals(bias_indicator, ", position(is.prior.PEESE), ")")
    )[order(c(position(is.prior.PET), position(is.prior.PEESE)))])
    expect_identical(nrow(out$declarations$declined), 0L)
    expect_true(out$uses_module)
  }
})

test_that("gates cover every weight-function, p-hacking, and combined branch leaf", {

  wf_independent <- prior_weightfunction("one-sided", c(0.025, 0.05), wf_independent(prior("gamma", list(9, 3))))
  wf_log <- prior_weightfunction("one-sided", c(0.025), wf_independent(prior("normal", list(0, 1)), "log_omega"))
  wf_fixed <- prior_weightfunction("two-sided", c(0.05), wf_fixed(c(1, 0.4)))
  phacking <- prior_phacking(form = "linear")
  prior_list <- list(bias = prior_mixture(
    list(prior_none(), .wf_cumulative2(), .wf_cumulative3(), wf_independent, wf_log, wf_fixed,
         phacking, prior_bias(.wf_cumulative2(), phacking), prior_PET("normal", list(0, 1))),
    is_null = c(TRUE, rep(FALSE, 8))
  ))
  out <- .declare(prior_list, "model{\n  y ~ dnorm(PET * omega[1] * alpha, 1)\n}", sampler_gates = TRUE)
  active <- function(k) paste0("equals(bias_indicator, ", k, ")")
  expect_identical(.gates(out), c(
    omega_ratio_component_2 = active(2),
    `eta_component_3[1]` = active(3), `eta_component_3[2]` = active(3), `eta_component_3[3]` = active(3),
    `omega_local_component_4[2]` = active(4), `omega_local_component_4[3]` = active(4),
    `log_omega_component_5[2]` = active(5),
    alpha_component_7 = active(7),
    omega_ratio_component_8 = active(8), alpha_component_8 = active(8),
    PET_1 = active(9)
  ))
  expect_identical(out$declarations$gates$branch, c(2L, 3L, 3L, 3L, 4L, 4L, 5L, 7L, 8L, 8L, 9L))
  expect_identical(out$declarations$gates$prior, rep("bias", 11L))
  expect_identical(out$declarations$gates$marker, paste0("BayesTools_gate_", 1:11))
  expect_identical(nrow(out$declarations$declined), 0L)

  # each leaf statement and its marker are in the model text, markers before the closing brace
  for(i in seq_len(nrow(out$declarations$gates))){
    gate <- out$declarations$gates[i, ]
    expect_true(grepl(paste0(gate$marker, " ~ dbt_gate(", gate$node, ", ", gate$active, ")"), out$model, fixed = TRUE))
  }
  expect_true(grepl("BayesTools_gate_11 ~ dbt_gate(PET_1, equals(bias_indicator, 9))\n}", out$model, fixed = TRUE))
  expect_identical(unlist(out$data[out$declarations$gates$marker], use.names = FALSE), rep(0, 11L))
  expect_identical(out$data$y, 0)
})

test_that("gates are not declared without a bias mixture of several branches or for point leaves", {

  cases <- list(
    list(bias = prior_mixture(list(prior_PET("normal", list(0, 1))), is_null = FALSE)),
    list(omega = .wf_cumulative3()),
    list(mu = prior_mixture(list(prior("normal", list(0, 1)), prior("normal", list(1, 1))))),
    list(mu = prior_spike_and_slab(prior("normal", list(0, 1)))),
    list(bias = prior_mixture(list(prior_none(), prior_PET("point", list(0.5))), is_null = c(TRUE, FALSE))),
    list(bias = prior_mixture(list(prior_none(), prior_phacking(form = "linear", alpha = prior("point", list(0.3)))),
                              is_null = c(TRUE, FALSE)))
  )
  for(prior_list in cases){
    out <- .declare(prior_list, sampler_gates = TRUE)
    expect_identical(nrow(out$declarations$gates), 0L)
    expect_false(out$uses_module)
    expect_false(grepl("BayesTools_gate_", out$model, fixed = TRUE))
  }
})

test_that("discrete and point leaves of a bias-mixture branch get no gate and no declined row", {

  cases <- list(
    list(bias = prior_mixture(list(prior_none(), prior_PET("bernoulli", list(0.5))), is_null = c(TRUE, FALSE))),
    list(bias = prior_mixture(list(prior_none(), prior_PEESE("point", list(0.5))), is_null = c(TRUE, FALSE))),
    list(bias = prior_mixture(list(prior_none(), prior_phacking(form = "linear", alpha = prior("point", list(0.3)))),
                              is_null = c(TRUE, FALSE)))
  )
  for(prior_list in cases){
    out <- .declare(prior_list, "model{\n  y ~ dnorm(0, 1)\n}", sampler_gates = TRUE)
    expect_identical(nrow(out$declarations$gates), 0L)
    expect_identical(nrow(out$declarations$declined), 0L)
    expect_false(out$uses_module)
  }
  # a continuous branch next to a discrete one keeps its gate
  prior_list <- list(bias = prior_mixture(
    list(prior_none(), prior_PET("bernoulli", list(0.5)), prior_PEESE("normal", list(0, 1))),
    is_null = c(TRUE, FALSE, FALSE)
  ))
  out <- .declare(prior_list, "model{\n  y ~ dnorm(PET + PEESE, 1)\n}", sampler_gates = TRUE)
  expect_identical(.gates(out), c(PEESE_1 = "equals(bias_indicator, 3)"))
  expect_identical(nrow(out$declarations$declined), 0L)
})

test_that("gates are declined when a branch node occurs outside its gated composition", {

  bias <- prior_mixture(list(prior_none(), prior_PET("normal", list(0, 1)), .wf_cumulative3()),
                        is_null = c(TRUE, FALSE, FALSE))
  wf_gates <- paste0("eta_component_3[", 1:3, "]")

  # composed nodes and monitoring are not occurrences
  out <- .declare(list(bias = bias), "model{\n  y ~ dnorm(PET + omega[1], 1)\n}",
                  sampler_gates = TRUE, add_parameters = c("PET_1", "eta_component_3"))
  expect_identical(names(.gates(out)), c("PET_1", wf_gates))

  # caller syntax naming a leaf or an intermediate of a branch
  out <- .declare(list(bias = bias), "model{\n  y ~ dnorm(PET_1, 1)\n}", sampler_gates = TRUE)
  expect_identical(names(.gates(out)), wf_gates)
  expect_identical(out$declarations$declined$kind, "gate")
  expect_identical(out$declarations$declined$branch, 2L)
  expect_identical(out$declarations$declined$block, NA_integer_)
  expect_identical(out$declarations$declined$nodes, "PET_1")
  expect_match(out$declarations$declined$reason, "'PET_1' occurs outside the statements of the mixture", fixed = TRUE)

  out <- .declare(list(bias = bias), "model{\n  y ~ dnorm(omega_component_3[2] + PET, 1)\n}", sampler_gates = TRUE)
  expect_identical(names(.gates(out)), "PET_1")
  expect_identical(out$declarations$declined$branch, 3L)

  out <- .declare(list(bias = bias), "model{\n  y ~ dnorm(std_eta_component_3[1], 1)\n}", sampler_gates = TRUE)
  expect_identical(names(.gates(out)), "PET_1")

  # a comment is not an occurrence
  out <- .declare(list(bias = bias), "model{\n  y ~ dnorm(PET, 1) # not PET_1\n}", sampler_gates = TRUE)
  expect_identical(names(.gates(out)), c("PET_1", wf_gates))

  # another prior's expression
  prior_list <- list(
    bias = bias,
    mu = prior("normal", list(mean = expression(PET_1), sd = 1))
  )
  out <- .declare(prior_list, "model{\n  y ~ dnorm(mu + PET, 1)\n}", sampler_gates = TRUE)
  expect_identical(names(.gates(out)), wf_gates)
  expect_identical(out$declarations$declined$branch, 2L)
})

test_that("reserved marker names and invalid arguments stop with the declaration error class", {

  prior_list <- list(mu = prior("normal", list(0, 1)), tau = prior("normal", list(0, 1)))
  expect_error(
    .declare(prior_list, data = list(y = 0, BayesTools_gate_1 = 0), sampler_gates = TRUE),
    "The model uses reserved BayesTools sampler node name 'BayesTools_gate_1'.",
    fixed = TRUE, class = "BayesTools_sampler_declaration_error"
  )
  expect_error(
    .declare(prior_list, "model{\n  BayesTools_block_7 <- 1\n  y ~ dnorm(0, 1)\n}", sampler_blocks = list(c("mu", "tau"))),
    class = "BayesTools_sampler_declaration_error"
  )
  expect_error(
    .declare(prior_list, sampler_gates = TRUE, add_parameters = "BayesTools_gate_x"),
    class = "BayesTools_sampler_declaration_error"
  )
  # without declarations the names are not reserved
  expect_identical(.declare(prior_list, data = list(y = 0, BayesTools_gate_1 = 0))$data$BayesTools_gate_1, 0)

  for(value in list(NA, "yes", c(TRUE, FALSE), 1)){
    expect_error(.declare(prior_list, sampler_gates = value),
                 "The 'sampler_gates' argument must be TRUE or FALSE.", fixed = TRUE,
                 class = "BayesTools_sampler_declaration_error")
  }
  for(value in list(list(), "mu", list(c("mu", NA)), list(character()), list(1:2), list(c("mu", "")))){
    expect_error(.declare(prior_list, sampler_blocks = value),
                 paste0("The 'sampler_blocks' argument must be NULL, a list of character vectors of prior names, ",
                        "or a function of the combined prior list that returns such a list or NULL."), fixed = TRUE,
                 class = "BayesTools_sampler_declaration_error")
  }
  condition <- tryCatch(.declare(prior_list, sampler_gates = NA), error = identity)
  expect_identical(class(condition), c("BayesTools_sampler_declaration_error", "error", "condition"))
  expect_null(conditionCall(condition))
})

test_that("blocks expand every supported prior to its scalar leaves", {

  set.seed(1)
  formula_data <- data.frame(x = stats::rnorm(24), g = factor(rep(1:4, 6)), d = factor(rep(1:3, 8)))
  allocation <- JAGS_formula(
    formula = ~ 1 + x + (1 | g) + random(1 | d, name = "d", covariance = "diag"),
    parameter = "mu", data = formula_data,
    prior_list = list(intercept = prior("normal", list(0, 1)), x = prior("normal", list(0, 1))),
    prior_random = prior_random(random_variance_allocation(
      name = "tot", terms = c(g = "g", d = "d"), scale = "total_variance",
      sd = prior("gamma", list(2, 2)),
      weights = prior("dirichlet", list(alpha = c(1.5, 2.5)))
    ))
  )$prior_list
  sd_name <- "mu__xRE_ALLOCx_tot__allocation_sd"
  weight_name <- "mu__xRE_ALLOCx_tot__weight"

  prior_list <- c(allocation, list(
    tau = prior("normal", list(0, 1), list(0, Inf)),
    rho = prior("beta", list(1, 1)),
    f_t = prior_factor("normal", list(0, 1), contrast = "treatment"),
    f_i = prior_factor("normal", list(0, 1), contrast = "independent"),
    omega = .wf_cumulative3(),
    PET = prior_PET("normal", list(0, 1)),
    PEESE = prior_PEESE("normal", list(0, 1))
  ))
  prior_list$f_t <- prior_factor_levels(prior_list$f_t, 3)
  prior_list$f_i <- prior_factor_levels(prior_list$f_i, 2)
  out <- .declare(prior_list, sampler_blocks = list(
    c("mu_intercept", "mu_x", sd_name, weight_name, "tau", "rho"),
    c("f_t", "f_i"),
    c("omega", "PET", "PEESE")
  ))
  blocks <- out$declarations$blocks
  expect_identical(blocks$node, c(
    "mu_intercept", "mu_x", sd_name,
    paste0("prior_par_eta_", weight_name, "[", 1:2, "]"), "tau", "rho",
    paste0("f_t[", 1:2, "]"), paste0("f_i[", 1:2, "]"),
    paste0("eta[", 1:3, "]"), "PET", "PEESE"
  ))
  expect_identical(blocks$block, rep(1:3, c(7L, 4L, 5L)))
  expect_identical(blocks$prior, c("mu_intercept", "mu_x", sd_name, weight_name, weight_name, "tau", "rho",
                                   "f_t", "f_t", "f_i", "f_i", "omega", "omega", "omega", "PET", "PEESE"))
  expect_identical(blocks$marker, paste0("BayesTools_block_", 1:16))
  expect_true(grepl("BayesTools_block_4 ~ dbt_block(prior_par_eta_mu__xRE_ALLOCx_tot__weight[1], 1)", out$model, fixed = TRUE))
  expect_true(grepl("BayesTools_block_16 ~ dbt_block(PEESE, 3)\n}", out$model, fixed = TRUE))
  expect_identical(nrow(out$declarations$gates), 0L)
  expect_true(out$uses_module)

  # single weight functions: binary cumulative ratio, independent weights on both scales
  wf_blocks <- function(omega){
    .declare(list(omega = omega, mu = prior("normal", list(0, 1))), sampler_blocks = list(c("omega", "mu")))$declarations$blocks$node
  }
  expect_identical(wf_blocks(.wf_cumulative2()), c("omega_ratio", "mu"))
  expect_identical(wf_blocks(prior_weightfunction("one-sided", c(0.025, 0.05), wf_independent(prior("gamma", list(9, 3))))),
                   c("omega[2]", "omega[3]", "mu"))
  expect_identical(wf_blocks(prior_weightfunction("two-sided", c(0.05), wf_independent(prior("normal", list(0, 1)), "log_omega"))),
                   c("log_omega[2]", "mu"))
})

test_that("sampler sources report emitted nodes and explain priors without eligible leaves", {

  priors <- list(
    mu = prior("normal", list(0, 1)),
    fixed = prior_weightfunction("two-sided", c(0.05), wf_fixed(c(1, 0.4))),
    cumulative = .wf_cumulative3(),
    point = prior("point", list(0))
  )
  sources <- JAGS_sampler_sources(priors, c("cumulative", "fixed", "point", "mu"))

  expect_identical(sources$prior, c(rep("cumulative", 3L), "fixed", "point", "mu"))
  expect_identical(sources$node, c(paste0("eta[", 1:3, "]"), NA_character_, NA_character_, "mu"))
  expect_identical(sources$status, c(rep("eligible", 3L), "empty", "declined", "eligible"))
  expect_identical(sources$reason, c(rep(NA_character_, 3L), "it has no supported scalar stochastic nodes",
                                    "it is a point prior", NA_character_))
  expect_identical(JAGS_sampler_sources(list()), sources[0L, ])
  expect_error(JAGS_sampler_sources(priors, "unknown"), class = "BayesTools_sampler_declaration_error")

  declared <- .declare(priors[c("cumulative", "mu")], sampler_blocks = list(c("cumulative", "mu")))
  expect_identical(declared$declarations$blocks$node, sources$node[sources$status == "eligible"])
})

test_that("ordinary mixture sources keep gate order, literal points and raw coefficient shape", {

  slab <- prior("normal", list(0, 1), prior_weights = 3)
  point <- prior("point", list(1.25), prior_weights = 7)
  priors <- list(mu = prior_mixture(list(slab = slab, point = point)))
  sources <- JAGS_mixture_sources(priors, "mu")
  expect_identical(sources$gate, "mu_indicator")
  expect_equal(sources$probabilities, c(0.3, 0.7))
  expect_identical(sources$branches[[1L]]$family, "scalar")
  expect_identical(sources$branches[[1L]]$source, "mu_component_1")
  expect_true(sources$branches[[1L]]$stochastic)
  expect_identical(sources$branches[[2L]]$family, "point")
  expect_identical(sources$branches[[2L]]$source, "mu_component_2")
  expect_identical(sources$branches[[2L]]$location, 1.25)
  expect_false(sources$branches[[2L]]$stochastic)

  rho <- list(rho = prior_mixture(list(prior("beta", list(2, 3)))))
  rho_sources <- JAGS_mixture_sources(rho, "rho")
  expect_identical(rho_sources$gate, "rho_indicator")
  expect_identical(rho_sources$probabilities, 1)
  expect_identical(rho_sources$branches[[1L]]$source, "rho_component_1")
  expect_true(rho_sources$branches[[1L]]$stochastic)
  expect_match(JAGS_add_priors("model{ y ~ dnorm(0, 1) }", rho),
               "rho_component_1 ~ dbeta(", fixed = TRUE)

  group <- data.frame(group = factor(rep(letters[1:4], 2)))
  factor_slab <- prior_factor("mnormal", list(0, 1), contrast = "meandif")
  compiled <- JAGS_formula(~ group, "mu", group, list(
    intercept = prior("normal", list(0, 1)),
    group = prior_mixture(list(prior("point", list(0)), factor_slab), is_null = c(TRUE, FALSE))
  ))
  factor_sources <- JAGS_mixture_sources(compiled$prior_list, "mu_group")
  block <- factor_sources$branches[[2L]]
  expect_identical(factor_sources$gate, "mu_group_indicator")
  expect_identical(block$family, "vector")
  expect_identical(block$coordinates, paste0("mu_group_component_2[", 1:3, "]"))
  expect_true(block$stochastic)
  syntax <- JAGS_add_priors("model{ y ~ dnorm(0, 1) }", compiled$prior_list)
  expect_match(syntax, paste0(block$source, " ~ dmnorm("), fixed = TRUE)
})

test_that("bias mixture sources identify raw slopes and physical weights on the common grid", {

  priors <- list(bias = prior_mixture(list(
    pet = prior_PET("normal", list(0, 1), prior_weights = 2),
    wide = prior_weightfunction("one-sided", c(0.025, 0.5), wf_cumulative(c(1, 1, 1)), prior_weights = 3),
    none = prior_none(prior_weights = 7),
    peese = prior_PEESE("normal", list(0, 1), prior_weights = 4),
    coarse = prior_weightfunction("one-sided", 0.025, wf_cumulative(c(1, 1)), prior_weights = 5)
  )))
  sources <- JAGS_mixture_sources(priors, "bias")
  expect_identical(sources$gate, "bias_indicator")
  expect_equal(sources$probabilities, c(2, 3, 7, 4, 5) / 21)
  expect_identical(vapply(sources$branches, `[[`, character(1), "family"),
                   c("PET", "selection", "none", "PEESE", "selection"))
  expect_identical(vapply(sources$branches[c(1L, 4L)], `[[`, character(1), "source"), c("PET_1", "PEESE_1"))
  expect_true(all(vapply(sources$branches[c(1L, 4L)], `[[`, logical(1), "stochastic")))
  expect_identical(sources$branches[[3L]]$omega, character())
  syntax <- JAGS_add_priors("model{ y ~ dnorm(0, 1) }", priors)
  for(i in c(2L, 5L)){
    expect_identical(sources$branches[[i]]$omega, paste0("omega_component_", i, "[", 1:3, "]"))
    for(node in sources$branches[[i]]$omega){
      expect_match(syntax, paste0(node, " <- "), fixed = TRUE)
    }
  }
  expect_error(JAGS_mixture_sources(list(bias = prior_mixture(list(.wf_cumulative2()))), "bias"),
               "at least two branches", fixed = TRUE, class = "BayesTools_sampler_declaration_error")
})

test_that("nested weight markers follow both active branches and their private dependencies", {

  model <- "model{\n  y ~ dnorm(omega[2], 1)\n}"
  priors <- list(bias = prior_mixture(list(
    prior_none(),
    prior_weightfunction("one-sided", c(0.05, 0.5), wf_cumulative(c(1, 1, 1))),
    prior_PET("normal", list(0, 1)),
    prior_weightfunction("one-sided", c(0.025, 0.05, 0.5), wf_cumulative(c(1, 1, 1, 1)))
  )))
  bound <- .bt_nested_weights_bind(priors, model, add_parameters = NULL,
                                   sampler_nested_weights = TRUE)$prior_list
  plan <- .bt_nested_weights_plan(bound$bias)
  expect_identical(c(plan$coarse, plan$fine), c(2L, 4L))
  out <- .declare(bound, model_syntax = model, sampler_gates = TRUE)
  markers <- out$declarations$gates
  expect_identical(markers$node, c("PET_1", unname(plan$nodes)))
  expect_identical(markers$active, c("equals(bias_indicator, 3)",
    rep("equals(bias_indicator, 2) + equals(bias_indicator, 4)", 2L), "equals(bias_indicator, 4)"))
  expect_identical(markers$branch, c(3L, NA_integer_, NA_integer_, 4L))
  expect_identical(nrow(out$declarations$declined), 0L)
  expect_true(all(vapply(unname(plan$nodes), function(node) {
    grepl(paste0("dbt_gate(", node, ","), out$model, fixed = TRUE)
  }, logical(1))))

  bypass <- paste0("model{\n  y ~ dnorm(omega[2] + ", plan$nodes[["a"]], ", 1)\n}")
  guarded <- .declare(bound, model_syntax = bypass, sampler_gates = TRUE)$declarations
  expect_identical(guarded$gates$node, "PET_1")
  expect_identical(nrow(guarded$declined), 2L)
  expect_true(all(grepl("occurs outside the statements of the mixture", guarded$declined$reason, fixed = TRUE)))
})

test_that("gate claims use only emitted copies to resolve reported JAGS aliases", {

  priors <- list(bias = prior_mixture(list(
    prior_none(),
    prior_weightfunction("one-sided", 0.025, wf_cumulative(c(1, 1))),
    prior_weightfunction("one-sided", c(0.025, 0.5), wf_cumulative(c(1, 1, 1))),
    prior_PET("normal", list(0, 1))
  )))
  gates <- .declare(priors, sampler_gates = TRUE)$declarations$gates
  samplers <- data.frame(Index.No = 1:4, Sampler = c("BayesTools::GatedSlice", rep("base::RealSlicer", 3L)),
                        Node = c("omega_component_2[2]", "omega_component_3[2]", "eta_component_3[2]", "PET_1"))
  claims <- .bt_sampler_gate_claims(gates, priors, samplers)
  expect_identical(claims$marker, gates$marker)
  expect_identical(claims$node, gates$node)
  expect_identical(claims$reported_node, c("omega_component_2[2]", NA_character_, "eta_component_3[2]",
                                         NA_character_, "PET_1"))
  expect_identical(claims$assignment, c(1L, NA_integer_, 3L, NA_integer_, 4L))
  expect_identical(claims$status, c("claimed", "unresolved", "other", "unresolved", "other"))

  ambiguous <- rbind(samplers, data.frame(Index.No = 5L, Sampler = "BayesTools::GatedSlice",
                                          Node = "omega_local_component_2[2]"))
  expect_identical(.bt_sampler_gate_claims(gates, priors, ambiguous)$status[[1L]], "unresolved")
})

test_that("blocks refuse guarded priors with the declaration error class", {

  set.seed(1)
  formula_data <- data.frame(x = stats::rnorm(24), g = factor(rep(1:4, 6)), d = factor(rep(1:3, 8)))
  inclusion <- JAGS_formula(
    formula = ~ 1 + (1 | g) + random(1 | d, name = "d", covariance = "diag"),
    parameter = "mu", data = formula_data,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(random_variance_allocation(
      name = "tot", terms = c(g = "g", d = "d"), scale = "total_variance",
      sd = prior_spike_and_slab(prior("gamma", list(2, 2))),
      weights = prior("dirichlet", list(alpha = c(1.5, 2.5))),
      inclusion = list(g = prior("beta", list(1, 1)))
    ))
  )$prior_list
  # random-effect priors that are not simple continuous scalars: a factor SD,
  # point, spike-and-slab and mixture SDs, a point correlation, a vector prior
  # with random-effect metadata, and a summary prior
  sd_prior <- prior("normal", list(0, 1), list(0, Inf))
  factor_sd <- JAGS_formula(
    formula = ~ 1 + diag(1 + d | g), parameter = "mu", data = formula_data,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(g = random_block(sd = sd_prior))
  )$prior_list
  random_prior <- function(prior, type = "random_sd"){
    attr(prior, type) <- TRUE
    attr(prior, "random_factor") <- "g"
    prior
  }
  summary_sd <- random_prior(sd_prior)
  attr(summary_sd, "random_summary") <- "sd"
  # simple scalar priors of a variance allocation that are neither its plain SD
  # nor its Dirichlet weights, with only the allocation metadata or also that of
  # a direct SD: only the direct SDs and correlations outside an allocation join
  allocation_only <- sd_prior
  attr(allocation_only, "random_allocation") <- "tot"
  allocation_sd_metadata <- random_prior(sd_prior)
  attr(allocation_sd_metadata, "random_allocation") <- "tot"
  random_guarded <- list(
    re_allocation = allocation_only,
    re_allocation_sd_metadata = allocation_sd_metadata,
    re_point = random_prior(prior("point", list(0.3))),
    re_s_s = random_prior(prior_spike_and_slab(sd_prior)),
    re_mix = random_prior(prior_mixture(list(sd_prior, prior("gamma", list(2, 2))))),
    re_cor_point = random_prior(prior("spike", list(0.3)), "random_correlation"),
    re_vector = random_prior(prior("mnormal", list(mean = 0, sd = 1, K = 2))),
    re_summary = summary_sd
  )

  prior_list <- c(inclusion, factor_sd[setdiff(names(factor_sd), "mu_intercept")], random_guarded, list(
    mu = prior("normal", list(0, 1)),
    s_s = prior_spike_and_slab(prior("normal", list(0, 1))),
    mix = prior_mixture(list(prior("normal", list(0, 1)), prior("normal", list(1, 1)))),
    bias = prior_mixture(list(prior_none(), prior_PET("normal", list(0, 1))), is_null = c(TRUE, FALSE)),
    ord = prior_ordered(prior("normal", list(0, 1), list(0, Inf))),
    point = prior("point", list(0)),
    discrete = prior("bernoulli", list(0.5)),
    orth = prior_factor("mnormal", list(0, 1), contrast = "orthonormal"),
    mnorm = prior("mnormal", list(mean = 0, sd = 1, K = 2)),
    dir = prior("dirichlet", list(alpha = c(1, 1))),
    phack = prior_phacking(form = "linear")
  ))
  prior_list$orth <- prior_factor_levels(prior_list$orth, 3)
  # the declarations of blocks alone: the guarded priors need not be emitted here
  .blocks <- function(sampler_blocks){
    BayesTools:::.bt_sampler_block_declarations(sampler_blocks, prior_list)
  }
  expect_block_error <- function(name, message){
    expect_error(
      .blocks(list(c("mu", name))),
      paste0("Prior '", name, "' cannot join a block of 'sampler_blocks': ", message, "."),
      fixed = TRUE, class = "BayesTools_sampler_declaration_error"
    )
  }
  expect_block_error("mu__xRE_ALLOCx_tot__allocation_sd", "it is a spike-and-slab prior")
  expect_block_error("mu__xRE_ALLOCx_tot__include_g_prob", "it is an inclusion prior")
  random_other <- paste0("it is a random-effect prior other than a simple continuous scalar prior, ",
                         "a plain allocation SD or Dirichlet allocation weights")
  expect_block_error("mu__xREx__g_d", random_other)
  expect_block_error("re_vector", random_other)
  expect_block_error("re_summary", random_other)
  expect_block_error("re_allocation", random_other)
  expect_block_error("re_allocation_sd_metadata", random_other)
  expect_block_error("re_point", "it is a point prior")
  expect_block_error("re_s_s", "it is a spike-and-slab prior")
  expect_block_error("re_mix", "it is a mixture prior")
  expect_block_error("re_cor_point", "it is a point prior")
  expect_block_error("s_s", "it is a spike-and-slab prior")
  expect_block_error("mix", "it is a mixture prior")
  expect_block_error("bias", "it is a mixture prior")
  expect_block_error("ord", "it is an ordered prior")
  expect_block_error("point", "it is a point prior")
  expect_block_error("discrete", "it is a discrete prior")
  expect_block_error("orth", "it is a vector prior")
  expect_block_error("mnorm", "it is a vector prior")
  expect_block_error("dir", "it is a vector prior")
  expect_block_error("phack", "its prior type is not supported in a block")

  # the allocation weights stay allowed next to a spike-and-slab allocation SD,
  # and the scalar SD of a factor block's intercept next to its factor SD
  expect_identical(.blocks(list(c("mu", "mu__xRE_ALLOCx_tot__weight")))$declared$node,
                   c("mu", paste0("prior_par_eta_mu__xRE_ALLOCx_tot__weight[", 1:2, "]")))
  expect_identical(.blocks(list(c("mu", "mu__xREx__g_intercept")))$declared$node,
                   c("mu", "mu__xREx__g_intercept"))

  expect_error(.blocks(list(c("mu", "nope"))),
               "Prior 'nope' in 'sampler_blocks' is not in the combined prior list of the model.", fixed = TRUE,
               class = "BayesTools_sampler_declaration_error")
  expect_error(.blocks(list(c("mu", "s_s"), "mu")),
               class = "BayesTools_sampler_declaration_error")
  expect_error(.blocks(list(c("mu", "x1"), c("mu", "x2"))),
               "Prior 'mu' appears more than once", fixed = TRUE,
               class = "BayesTools_sampler_declaration_error")
  expect_error(.blocks(list(c("mu", "mu"))),
               "Prior 'mu' appears more than once in 'sampler_blocks'; a prior may be in one block only.", fixed = TRUE,
               class = "BayesTools_sampler_declaration_error")
})

test_that("blocks take random-effect SDs and scalar correlations with simple continuous priors", {

  set.seed(1)
  formula_data <- data.frame(x = stats::rnorm(24), g = factor(rep(1:4, 6)), s = factor(rep(1:6, each = 4)),
                             t = rep(1:3, 8))
  sd_prior <- prior("normal", list(0, 1), list(0, Inf))
  prior_list <- JAGS_formula(
    formula = ~ 1 + x + (1 + x | g) + cs(t | s), parameter = "mu", data = formula_data,
    prior_list = list(intercept = prior("normal", list(0, 1)), x = prior("normal", list(0, 1))),
    prior_random = prior_random(
      g = random_block(sd = sd_prior, cor = prior_lkj(eta = 1)),
      s = random_block(sd = prior("gamma", list(2, 2)), covariance = random_covariance(
        cor = prior("uniform", list(0, 0.9)), cor_scale = "cor"))
    )
  )$prior_list
  random_names <- c("mu__xREx__g_intercept", "mu__xREx__g_x", "mu__xREx__s_sd", "mu__xREx__s_rho")
  expect_true(all(random_names %in% names(prior_list)))
  expect_true(all(vapply(prior_list[random_names], BayesTools:::.bt_is_random_effect_prior, logical(1))))
  expect_true(isTRUE(attr(prior_list[["mu__xREx__s_rho"]], "random_correlation")))

  # each prior is a leaf under its own name, in a block with fixed coefficients
  out <- .declare(prior_list, sampler_blocks = list(c("mu_intercept", random_names[1:2]), c(random_names[3:4], "mu_x")))
  blocks <- out$declarations$blocks
  expect_identical(blocks$node, c("mu_intercept", random_names, "mu_x"))
  expect_identical(blocks$prior, blocks$node)
  expect_identical(blocks$block, rep(1:2, each = 3L))
  expect_identical(nrow(out$declarations$declined), 0L)
  expect_true(grepl("BayesTools_block_5 ~ dbt_block(mu__xREx__s_rho, 2)", out$model, fixed = TRUE))
  expect_true(out$uses_module)

  # a block of a single random-effect SD is declined like any other single leaf
  single <- .declare(prior_list, sampler_blocks = list("mu__xREx__s_sd"))
  expect_identical(nrow(single$declarations$blocks), 0L)
  expect_identical(single$declarations$declined$nodes, "mu__xREx__s_sd")
})

test_that("a block of fewer than two sampled nodes is declined, not declared", {

  prior_list <- list(
    mu = prior("normal", list(0, 1)),
    tau = prior("normal", list(0, 1), list(0, Inf)),
    omega = prior_weightfunction("two-sided", c(0.05), wf_fixed(c(1, 0.4)))
  )
  model_syntax <- "model{\n  y ~ dnorm(mu, 1 / tau^2)\n}"
  plain <- .declare(prior_list, model_syntax)
  for(sampler_blocks in list(list("mu"), list("omega"))){
    out <- .declare(prior_list, model_syntax, sampler_blocks = sampler_blocks)
    expect_identical(out$model, plain$model)
    expect_identical(out$data, plain$data)
    expect_false(out$uses_module)
    expect_identical(nrow(out$declarations$blocks), 0L)
    expect_identical(out$declarations$declined, data.frame(
      kind = "block", prior = sampler_blocks[[1L]], branch = NA_integer_, block = 1L,
      nodes = if(identical(sampler_blocks[[1L]], "mu")) "mu" else "",
      reason = "the block expands to fewer than two sampled nodes", stringsAsFactors = FALSE
    ))
  }
  # a declined block next to a declared one: block ids keep their position
  out <- .declare(prior_list, model_syntax, sampler_blocks = list("omega", c("mu", "tau")))
  expect_identical(out$declarations$blocks$node, c("mu", "tau"))
  expect_identical(out$declarations$blocks$block, c(2L, 2L))
  expect_identical(out$declarations$blocks$marker, paste0("BayesTools_block_", 1:2))
  expect_identical(out$declarations$declined$block, 1L)
  expect_true(grepl("BayesTools_block_1 ~ dbt_block(mu, 2)", out$model, fixed = TRUE))
  expect_true(out$uses_module)
  # guarded priors stay errors in a block that would be declined
  expect_error(.declare(c(prior_list, list(point = prior("point", list(0)))), sampler_blocks = list("point")),
               "it is a point prior", class = "BayesTools_sampler_declaration_error")
})

test_that("a function of the combined prior list may select the blocks", {

  set.seed(1)
  formula_data <- data.frame(x = stats::rnorm(24), g = factor(rep(1:4, 6)), d = factor(rep(1:3, 8)))
  prior_list <- c(JAGS_formula(
    formula = ~ 1 + x + (1 | g) + random(1 | d, name = "d", covariance = "diag"),
    parameter = "mu", data = formula_data,
    prior_list = list(intercept = prior("normal", list(0, 1)), x = prior("normal", list(0, 1))),
    prior_random = prior_random(random_variance_allocation(
      name = "tot", terms = c(g = "g", d = "d"), scale = "total_variance",
      sd = prior("gamma", list(2, 2)),
      weights = prior("dirichlet", list(alpha = c(1.5, 2.5)))
    ))
  )$prior_list, list(
    tau = prior("normal", list(0, 1), list(0, Inf)),
    bias = prior_mixture(list(prior_none(), prior_PET("normal", list(0, 1))), is_null = c(TRUE, FALSE))
  ))
  # by attribute: the plain allocation SD and the Dirichlet allocation weights;
  # by class: the simple scalar priors
  calls <- 0L
  by_metadata <- function(priors){
    calls <<- calls + 1L
    allocation <- vapply(priors, function(prior){
      isTRUE(attr(prior, "random_allocation_sd", exact = TRUE)) ||
        (!is.null(attr(prior, "random_allocation", exact = TRUE)) && identical(prior[["distribution"]], "dirichlet"))
    }, logical(1))
    simple <- vapply(priors, function(prior){
      is.prior.simple(prior) && !is.prior.point(prior) && !is.prior.factor(prior) &&
        !is.prior.PET(prior) && !is.prior.PEESE(prior) && is.null(attributes(prior)[["random_allocation"]])
    }, logical(1))
    list(names(priors)[allocation], names(priors)[simple])
  }
  explicit <- list(
    c("mu__xRE_ALLOCx_tot__allocation_sd", "mu__xRE_ALLOCx_tot__weight"),
    c("mu_intercept", "mu_x", "tau")
  )
  from_function <- .declare(prior_list, sampler_gates = TRUE, sampler_blocks = by_metadata)
  from_list <- .declare(prior_list, sampler_gates = TRUE, sampler_blocks = explicit)
  expect_identical(calls, 1L)
  expect_identical(from_function, from_list)
  expect_identical(from_function$declarations$blocks$prior[c(1L, 4L)],
                   c("mu__xRE_ALLOCx_tot__allocation_sd", "mu_intercept"))

  # NULL declares nothing: text and data unchanged
  plain <- .declare(prior_list)
  none <- .declare(prior_list, sampler_blocks = function(priors) NULL)
  expect_identical(none$model, plain$model)
  expect_identical(none$data, plain$data)
  expect_false(none$uses_module)
  expect_identical(nrow(none$declarations$blocks), 0L)

  # invalid values and guarded priors
  for(value in list(list(), "mu_x", list(c("mu_x", NA)), list(1:2), 1, prior("normal", list(0, 1)))){
    expect_error(.declare(prior_list, sampler_blocks = function(priors) value),
                 "The 'sampler_blocks' function returned an invalid value; it must return NULL or a list of character vectors of prior names.",
                 fixed = TRUE, class = "BayesTools_sampler_declaration_error")
  }
  expect_error(.declare(prior_list, sampler_blocks = function(priors) list(c("mu_x", "bias"))),
               "Prior 'bias' cannot join a block of 'sampler_blocks': it is a mixture prior.", fixed = TRUE,
               class = "BayesTools_sampler_declaration_error")
})

test_that("blocks and gates of one fit are numbered by kind in emission order", {

  prior_list <- list(
    mu = prior("normal", list(0, 1)),
    tau = prior("normal", list(0, 1), list(0, Inf)),
    bias = prior_mixture(list(prior_none(), prior_PEESE("normal", list(0, 1)), prior_PET("normal", list(0, 1))),
                         is_null = c(TRUE, FALSE, FALSE))
  )
  out <- .declare(prior_list, "model{\n  y ~ dnorm(mu + PET + PEESE, 1/tau^2)\n}",
                  sampler_gates = TRUE, sampler_blocks = list(c("tau", "mu")))
  expect_identical(out$declarations$gates$marker, paste0("BayesTools_gate_", 1:2))
  expect_identical(out$declarations$gates$node, c("PEESE_1", "PET_1"))
  expect_identical(out$declarations$blocks$node, c("tau", "mu"))
  expect_identical(names(out$data), c("y", "BayesTools_gate_1", "BayesTools_gate_2", "BayesTools_block_1", "BayesTools_block_2"))
  markers <- regmatches(out$model, gregexpr("BayesTools_(gate|block)_[0-9]+ ~ [^\n]+", out$model))[[1L]]
  expect_identical(markers, c(
    "BayesTools_gate_1 ~ dbt_gate(PEESE_1, equals(bias_indicator, 2))",
    "BayesTools_gate_2 ~ dbt_gate(PET_1, equals(bias_indicator, 3))",
    "BayesTools_block_1 ~ dbt_block(tau, 1)",
    "BayesTools_block_2 ~ dbt_block(mu, 1)"
  ))
  # the text before the markers is the text without declarations
  plain <- .declare(prior_list, "model{\n  y ~ dnorm(mu + PET + PEESE, 1/tau^2)\n}")$model
  expect_identical(sub("  BayesTools_gate_1[^}]*", "", out$model), plain)
})
