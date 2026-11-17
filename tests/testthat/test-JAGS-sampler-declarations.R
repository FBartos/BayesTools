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

test_that("JAGS_fit() hands runjags the same model text and data with default arguments", {

  skip_if_not_installed("runjags")
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

  # with declarations the markers and their data reach runjags and the module is required
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
  expect_identical(out$declarations$declined$branch, 2L)
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
                 "The 'sampler_blocks' argument must be NULL or a list of character vectors of prior names.", fixed = TRUE,
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
  block_sd <- JAGS_formula(
    formula = ~ 1 + (1 + x | g), parameter = "mu", data = formula_data,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(g = random_block(sd = prior("normal", list(0, 1), list(0, Inf)), cor = prior_lkj(eta = 1)))
  )$prior_list

  prior_list <- c(inclusion, block_sd[setdiff(names(block_sd), "mu_intercept")], list(
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
  expect_block_error("mu__xREx__g_intercept",
                     "it is a random-effect prior other than a plain allocation SD or Dirichlet allocation weights")
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

  # the allocation weights stay allowed next to a spike-and-slab allocation SD
  expect_identical(.blocks(list(c("mu", "mu__xRE_ALLOCx_tot__weight")))$node,
                   c("mu", paste0("prior_par_eta_mu__xRE_ALLOCx_tot__weight[", 1:2, "]")))

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
  expect_error(.blocks(list("mu")),
               "Block 1 of 'sampler_blocks' expands to fewer than two sampled nodes.", fixed = TRUE,
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
