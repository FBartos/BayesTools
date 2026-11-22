skip_if_not_test_profile("unit")

test_that("JAGS_fit binds package sampler declarations before creating the backend", {

  skip_if_not_installed("rjags")
  called <- 0L
  supplied <- NULL
  backend <- NULL
  testthat::local_mocked_bindings(
    .JAGS_load_modules = function(...) invisible(TRUE),
    .bt_rjags_fit = function(...){
      backend <<- list(...)
      simpleError("Captured the backend input")
    },
    .package = "BayesTools"
  )
  declare <- function(data, formula_design, prior_list){
    called <<- called + 1L
    supplied <<- list(priors = prior_list, design = formula_design, data = data)
    list(markers = list(list(name = "Package_joint_marker", distribution = "dpackage_marker",
      arguments = list(c("mu_intercept", "mu_x"), "Package_joint_data", c(1, 2)))),
      data = list(Package_joint_data = matrix(1:6, 2)))
  }
  observed <- c(1, 2, 4)
  normal <- prior("normal", list(0, 1))
  fit <- JAGS_fit(
    "model { for(i in 1:N_mu) { y[i] ~ dnorm(mu[i], tau) } }",
    data = list(y = observed), prior_list = list(tau = prior("gamma", list(1, 1))),
    formula_list = list(mu = ~ x), formula_data_list = list(mu = data.frame(x = c(0, 1, 3))),
    formula_prior_list = list(mu = list(intercept = normal, x = normal)),
    sampler_joint = declare, jags_modules = "Package", seed = 1L
  )
  expect_identical(called, 1L)
  expect_named(supplied$priors, c("mu_intercept", "mu_x", "tau"))
  expect_named(supplied$design, "mu")
  expect_identical(supplied$data$y, observed)
  expect_true("N_mu" %in% names(supplied$data))
  expect_identical(backend$data$Package_joint_marker, 0)
  expect_identical(backend$data$Package_joint_data, matrix(1:6, 2))
  expect_match(backend$model,
    "Package_joint_marker ~ dpackage_marker(c(mu_intercept, mu_x), Package_joint_data, c(1, 2))",
    fixed = TRUE)
  expect_false(grepl("Package_joint_marker", attr(fit, "model_syntax"), fixed = TRUE))
  declarations <- attr(fit, "sampler_declarations", exact = TRUE)
  expect_identical(declarations$joint$data, "Package_joint_data")
  expect_identical(declarations$joint$markers[[1L]]$name, "Package_joint_marker")
})

test_that("package markers leave the scientific gate scan intact", {

  priors <- list(bias = prior_mixture(
    list(prior_none(), prior_PET("normal", list(0, 1))),
    is_null = c(TRUE, FALSE)
  ))
  data <- list(y = 0)
  syntax <- JAGS_add_priors("model { y ~ dnorm(PET, 1) }", priors)
  legacy <- .bt_sampler_declarations(syntax, data, priors, NULL, sampler_gates = TRUE)
  plan <- .bt_sampler_apply_plan(legacy, NULL, priors, NULL, data,
    sampler_joint = function(prior_list, formula_design, data){
      list(markers = list(list(name = "Package_joint_marker", distribution = "dbt_block",
        arguments = list("PET_1", 1))))
    },
    sampler_probabilities = "bias_indicator")
  expect_identical(plan$declarations$gates, legacy$declarations$gates)
  expect_identical(plan$declarations$gates$node, "PET_1")
  expect_match(plan$model, "BayesTools_finite_1 <- bt_finite_tag(bias_indicator)", fixed = TRUE)
  expect_identical(plan$declarations$probabilities$readout, "BayesTools_finite_readout_1")
})

test_that("invalid sampler declarations fail before JAGS compilation", {

  expect_error(.bt_sampler_check_plan_arguments(1, NULL), class = "BayesTools_sampler_declaration_error")
  for(nodes in list(c("g", "g "), "g[1:2]", "g[3 2]", "g\n", NA_character_)){
    expect_error(.bt_sampler_check_plan_arguments(NULL, nodes), class = "BayesTools_sampler_declaration_error")
  }
  expect_identical(.bt_sampler_check_plan_arguments(NULL, c("a", "b [ 2, 1 ]")), c("a", "b[2,1]"))

  marker <- list(name = "package_marker", distribution = "dpackage_marker", arguments = list("x", 1))
  check <- function(fields){
    marker[names(fields)] <- fields
    .bt_sampler_joint_markers(list(markers = list(marker)), c("x", "y"), "Package")
  }
  for(bad in list("y", "BayesTools_finite_readout_1", "mark\n", "x <- 1")){
    expect_error(check(list(name = bad)), class = "BayesTools_sampler_declaration_error")
  }
  for(bad in list("equals(x, 1)", "unknown", x ~ y, Inf, list(x = 2))){
    expect_error(check(list(arguments = list(bad))), class = "BayesTools_sampler_declaration_error")
  }
  expect_error(check(list(distribution = "dnorm")),
    "must be a BayesTools marker or use the prefix", fixed = TRUE, class = "BayesTools_sampler_declaration_error")
  expect_error(.bt_sampler_joint_markers(list(markers = list(marker, marker)), c("x", "y"), "Package"),
    "used more than once", fixed = TRUE, class = "BayesTools_sampler_declaration_error")
  expect_error(.bt_sampler_joint_markers(list(markers = list(marker), data = list(x = 0)), c("x", "y"), "Package"),
    "already used or reserved", fixed = TRUE, class = "BayesTools_sampler_declaration_error")
  expect_error(.bt_sampler_joint_markers(list(markers = list(marker)), c("x", "y"), "Package", retired = "x"),
    "use 'sampler_nested_weights' = FALSE", fixed = TRUE, class = "BayesTools_sampler_declaration_error")
  expect_error(.bt_sampler_joint_markers(list(markers = list(marker), data = list(spare = 1)), c("x", "y"), "Package"),
    "data 'spare' is not referenced by any marker", fixed = TRUE, class = "BayesTools_sampler_declaration_error")
})

test_that("sampler references come from model relations and inputs", {

  syntax <- paste(
    "data { n ~ dpois(4) }", "model {",
    "for(i in 1:2) { logit(prob[i]) <- theta }",
    "y ~ dbin(prob[1], n)", "theta ~ dnorm(0, 1)", "}", sep = "\n"
  )
  data <- list(y = 1, known_constant = 9)
  prior_list <- list(bias = prior("normal", list(0, 1)))
  legacy <- list(model = syntax, data = data)
  declare <- function(node){
    function(prior_list, formula_design, data){
      list(markers = list(list(name = "trace_marker", distribution = "dpackage_marker",
        arguments = list(node, "n", "known_constant"))))
    }
  }
  plan <- function(nodes = NULL, joint = NULL){
    .bt_sampler_apply_plan(legacy, "undeclared", prior_list, NULL, data,
      sampler_joint = joint, sampler_probabilities = nodes, jags_modules = "Package")
  }
  for(node in c("bias", "dnorm", "i", "logit", "undeclared")){
    expect_error(plan(nodes = node), "unavailable in the model", fixed = TRUE,
      class = "BayesTools_sampler_declaration_error")
    expect_error(plan(joint = declare(node)), paste0("unknown node '", node, "'"), fixed = TRUE,
      class = "BayesTools_sampler_declaration_error")
  }
  accepted <- plan(nodes = c("prob[1]", "n", "theta"), joint = declare("prob[1:2]"))
  expect_identical(accepted$declarations$probabilities$node, c("prob[1]", "n", "theta"))
  expect_match(accepted$model, "trace_marker ~ dpackage_marker(prob[1:2], n, known_constant)", fixed = TRUE)
})

test_that("callback and ordinary FactorSlice block markers cannot share one fit", {

  priors <- stats::setNames(rep(list(prior("normal", list(0, 1))), 3L), c("x", "w", "z"))
  data <- list(y = 1)
  syntax <- JAGS_add_priors("model { alias <- x\n y ~ dnorm(alias + w + z, 1) }", priors)
  legacy <- .bt_sampler_declarations(syntax, data, priors, NULL, sampler_blocks = list(c("x", "w")))
  declare <- function(target, id, distribution = "dbt_block"){
    function(prior_list, formula_design, data){
      list(markers = list(list(name = "package_marker", distribution = distribution,
        arguments = list(target, id))))
    }
  }
  plan <- function(joint) .bt_sampler_apply_plan(legacy, NULL, priors, NULL, data,
    sampler_joint = joint, sampler_probabilities = NULL, jags_modules = "Package")
  expect_error(plan(declare("alias", 2L)), "blocks from 'sampler_blocks'", fixed = TRUE,
    class = "BayesTools_sampler_declaration_error")
  expect_error(plan(declare("z", 1L)), "blocks from 'sampler_blocks'", fixed = TRUE,
    class = "BayesTools_sampler_declaration_error")
  expect_true(plan(declare("z", 1L, "dpackage_marker"))$uses_module)
})

test_that("typed callbacks cannot bypass the scientific gate check", {

  callback <- function(prior_list, formula_design, data){
    list(markers = list(list(name = "manual_gate", distribution = "dbt_gate",
      arguments = list("x", "open"))))
  }
  expect_error(JAGS_fit(
    "model { open ~ dbern(0.99)\n y ~ dnorm(x, 100) }", data = list(y = 10),
    prior_list = list(x = prior("normal", list(0, 1))),
    sampler_joint = callback, jags_modules = "bt",
    runtime_setup = function(context) stop("Unexpected backend setup")
  ), "distribution 'dbt_gate' is reserved for gates verified from the model. Use 'sampler_gates' to enable those gates.",
    fixed = TRUE, class = "BayesTools_sampler_declaration_error")
})

test_that("parallel merge checks probability schema by indicator and reports mixed gate ownership", {

  make_worker <- function(ids, start, gate_sampler){
    draws <- function(values) coda::mcmc.list(lapply(ids, function(i){
      coda::mcmc(values(i), start = start, thin = 2L)
    }))
    list(
      mcmc = draws(function(i) cbind(x = rep(i, 3L))),
      end_state = lapply(ids, function(id) list(x = id)),
      samplers = data.frame(Index.No = c(1L, 2L, 3L),
        Sampler = c(gate_sampler, "BayesTools::GatedSlice", "BayesTools::Finite"),
        Node = c("x", "w", "b")),
      sampler_probabilities = list(b = draws(function(i){
        p <- rep(i / 4, 3L)
        cbind("p[0]" = p, "p[1]" = 1 - p, "log_p[0]" = log(p),
          "log_p[1]" = log1p(-p), selected = 0, lower = 0)
      }))
    )
  }
  groups <- list(c(1L, 3L), 2L)
  workers <- list(make_worker(groups[[1L]], 101L, "BayesTools::GatedSlice"),
                  make_worker(groups[[2L]], 51L, "base::RealSlicer"))
  plan <- data.frame(node = "b")
  fit <- .bt_rjags_merge_workers(workers, groups, 2L, plan)
  expect_identical(fit$samplers, workers[[1L]]$samplers)
  expect_identical(lapply(fit$samplers_by_worker, "[[", "chains"), groups)
  expect_identical(lapply(fit$samplers_by_worker, "[[", "samplers"), lapply(workers, "[[", "samplers"))
  expect_identical(vapply(fit$mcmc, function(chain) chain[1L, "x"], integer(1)), 1:3)
  expect_identical(vapply(fit$sampler_probabilities$b, function(chain) 4 * chain[1L, "p[0]"], numeric(1)), c(1, 2, 3))
  expect_identical(lapply(fit$mcmc, stats::time), lapply(fit$sampler_probabilities$b, stats::time))

  gates <- data.frame(marker = c("first", "second"), node = c("x", "w"),
    branch = NA_integer_, prior = c("x", "w"))
  claims <- .bt_sampler_fit_gate_claims(gates, NULL, fit)
  expect_identical(claims$status, c("mixed", "claimed"))
  expect_identical(claims$sampler, c(NA_character_, "BayesTools::GatedSlice"))
  fit_without_readouts <- .bt_rjags_merge_workers(workers, groups, 2L, NULL)
  expect_identical(.bt_sampler_fit_gate_claims(gates, NULL, fit_without_readouts), claims)

  incompatible <- workers
  colnames(incompatible[[2L]]$sampler_probabilities$b[[1L]])[1L] <- "p[3]"
  expect_error(.bt_rjags_merge_workers(incompatible, groups, 2L, plan),
    "incompatible probability traces for 'b'", fixed = TRUE, class = "BayesTools_sampler_probabilities_error")
  incompatible[[2L]]$sampler_probabilities["b"] <- list(NULL)
  expect_error(.bt_rjags_merge_workers(incompatible, groups, 2L, plan),
    "incompatible probability traces for 'b'", fixed = TRUE, class = "BayesTools_sampler_probabilities_error")
})

test_that("readout decoding distinguishes an unowned node from an incomplete trace", {

  readout <- c(alias = "BayesTools_finite_readout_1", declined = "BayesTools_finite_readout_2")
  row <- c(0.25, 0.75, log(0.25), log(0.75), 1, 0, -1)
  values <- matrix(rep(row, 3L), nrow = 3L, byrow = TRUE,
    dimnames = list(NULL, c(paste0(readout[[1L]], "[", 1:6, "]"), readout[[2L]])))
  draws <- coda::mcmc.list(coda::mcmc(values, start = 20L, thin = 2L))
  decoded <- .bt_rjags_probability_draws(draws, readout)
  expect_identical(names(decoded), names(readout))
  expect_s3_class(decoded$alias, "mcmc.list")
  expect_null(decoded$declined)
  expect_identical(stats::time(decoded$alias[[1L]]), stats::time(draws[[1L]]))
  values[2L, 7L] <- 1
  invalid <- coda::mcmc.list(coda::mcmc(values))
  expect_error(.bt_rjags_probability_draws(invalid, readout),
    "incomplete probability trace for 'declined'", fixed = TRUE, class = "BayesTools_sampler_probabilities_error")
})
