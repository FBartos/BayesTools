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
})

test_that("finite trace ownership handles JAGS names of one-cell arrays", {

  plan <- data.frame(node = c("gate[1]", "other", "array[1,2]", "declined[1]"),
    readout = letters[1:4])
  assignments <- data.frame(Node = c("array[1,2]", "gate", "declined"), Sampler = c(
    "BayesTools::Finite", "BayesTools::Finite", "base::Finite"))
  expect_identical(.bt_rjags_finite_readouts(plan, assignments), c("gate[1]" = "a", "array[1,2]" = "c"))
})
