skip_if_not_test_profile("unit")

# A balanced three-level factor f and continuous x and z orthogonal to f and
# to each other, so that the term blocks of ~ x + f + z are orthogonal.
.R2_test_data <- function(){
  set.seed(11)
  f <- factor(rep(c("a", "b", "c"), each = 8))
  centered_f <- stats::model.matrix(~ f)[, -1] - rep(colMeans(stats::model.matrix(~ f)[, -1]), each = 24)
  orthogonal <- function(v, M) stats::residuals(stats::lm(v ~ M))
  x <- orthogonal(stats::rnorm(24), centered_f)
  z <- orthogonal(stats::rnorm(24), cbind(centered_f, x))
  data.frame(f = f, x = x, z = z, w = stats::rnorm(24))
}

.R2_beta <- function(a, b) prior("beta", list(alpha = a, beta = b))
.R2_dirichlet <- function(...) prior("dirichlet", list(alpha = c(...)))

.R2_test_budget <- function(...){
  arguments <- list(
    terms = c("x", "f", "z"), total = .R2_beta(1, 2), allocation = .R2_dirichlet(1, 2, 1),
    law = "exact", convention = "common", concentration = "inherited", sigma = "sigma",
    contrast = "orthonormal"
  )
  overrides <- list(...)
  arguments[names(overrides)] <- overrides
  do.call(.R2_prior, arguments)
}

test_that("prior_R2() has defaults, takes its fixed choices from R2_model() and validates its arguments", {

  # the defaults: every formula term, Beta(1, 1), a flat allocation, the
  # exact law with the common convention and inherited concentrations
  default <- prior_R2()
  expect_null(default$terms)
  expect_identical(default[c("a", "b", "law", "convention", "concentration", "sigma", "contrast")],
    list(a = 1, b = 1, law = "exact", convention = "common", concentration = "inherited", sigma = "sigma",
         contrast = "orthonormal"))
  expect_output(print(default), "every formula term without its own prior (resolved when the formula binds the budget)", fixed = TRUE)
  expect_output(print(default), "allocation: flat Dirichlet", fixed = TRUE)
  expect_error(prior_R2(model = list(law = "exact")), "'model' must be a specification from 'R2_model()'", fixed = TRUE)
  expect_error(prior_R2(law = "exact"), "unused argument")
  expect_error(R2_model(law = "student"), "law", fixed = TRUE)
  expect_error(R2_model(sigma = -1), "sigma")
  expect_error(R2_model(contrast = "treatment2"), "contrast")
  model <- R2_model(law = "normal", convention = "linked", sigma = 2)
  expect_identical(check_R2_model(model), model)
  expect_output(print(model), "convention:    linked", fixed = TRUE)
  expect_error(prior_R2(allocation = prior("dirichlet", list(alpha = c(1, 1)))), "'R2_dirichlet()'", fixed = TRUE)
  expect_error(R2_dirichlet(0), "alpha")
  expect_error(R2_dirichlet(c(x = 1, x = 2)), "names each term once", fixed = TRUE)
  expect_output(print(R2_dirichlet(0.5)), "symmetric Dirichlet (concentration 0.5 for every term)", fixed = TRUE)
  expect_error(.R2_test_budget(law = "student"), "law", fixed = TRUE)
  expect_error(.R2_test_budget(terms = c(1, 2)), "'terms' must be a character vector", fixed = TRUE)
  expect_error(.R2_test_budget(terms = c("x:f", "f:x", "z")), "duplicated: 'f:x'", fixed = TRUE)
  expect_error(.R2_test_budget(terms = c("x", "intercept", "z")), "The intercept cannot be part", fixed = TRUE)
  expect_error(.R2_test_budget(total = prior("normal", list(0, 1))), "'total' must be a Beta prior", fixed = TRUE)
  expect_error(.R2_test_budget(total = prior_mixture(list(.R2_beta(1, 1), prior("spike", list(0))), is_null = c(FALSE, TRUE))),
    "other mixtures are not supported", fixed = TRUE)
  expect_error(.R2_test_budget(total = prior("beta", list(1, 1), list(0, 0.5))), "cannot be truncated", fixed = TRUE)
  expect_error(.R2_test_budget(total = prior("beta", list(0.01, 1), list(1e-9, 1))), "cannot be truncated", fixed = TRUE)
  expect_equal(.R2_test_budget(total = prior("beta", list(alpha = c(shape = 2), beta = 3)))[c("a", "b")], list(a = 2, b = 3))
  expect_error(.R2_test_budget(terms = "x", allocation = .R2_dirichlet(1, 1)),
    "'allocation' has 2 concentrations, but 'terms' has 1 term(s)", fixed = TRUE)
  expect_error(.R2_test_budget(allocation = .R2_dirichlet(1, 1)), "has 2 concentrations, but 'terms' has 3", fixed = TRUE)
  # allocations: symmetric, named (matched to the terms, covering them), one term
  expect_equal(prior_R2(terms = c("x", "f"), allocation = R2_dirichlet(0.5))$alpha, c(x = 0.5, f = 0.5))
  expect_equal(prior_R2(terms = c("x", "f"), allocation = R2_dirichlet(c(f = 3, x = 1)))$alpha, c(x = 1, f = 3))
  expect_error(prior_R2(terms = c("x", "f"), allocation = R2_dirichlet(c(x = 1))), "each once", fixed = TRUE)
  expect_equal(prior_R2(terms = "x", allocation = R2_dirichlet(2))$alpha, c(x = 2))
  # a one-term linked budget has its own concentration as the full-model total
  linked_one <- prior_R2(terms = "x", total = .R2_beta(2, 3), allocation = R2_dirichlet(2), model = R2_model(convention = "linked"))
  expect_equal(linked_one$a, 2)
  expect_error(prior_R2(terms = "x", allocation = R2_dirichlet(2), model = R2_model(convention = "linked")),
    "must equal the total concentration of the full model, 2", fixed = TRUE)
  expect_error(.R2_test_budget(sigma = -1), "sigma")
  expect_error(.R2_test_budget(convention = "conditioning", total = .R2_beta(1, 1), terms = c("x", "f"),
    allocation = .R2_dirichlet(1, 2)), "a - sum(alpha) + min(alpha) > 0", fixed = TRUE)
  # a linked total has the full-model total concentration as its first shape
  expect_error(.R2_test_budget(convention = "linked"),
    "must equal the total concentration of the full model, 4 for this budget; use total = prior(\"beta\", list(alpha = 4, beta = 2))",
    fixed = TRUE, class = "BayesTools_R2_invalid")

  budget <- .R2_test_budget(convention = "linked", total = .R2_beta(4, 2), concentration = "fixed_total")
  expect_true(is.R2_budget(budget))
  expect_equal(budget$alpha, c(x = 1, f = 2, z = 1))
  expect_output(print(budget), "Beta(4, 2) in the full model; Beta(A_M, 2) in a model with retained concentration A_M", fixed = TRUE)
  expect_output(print(budget), "retained concentrations rescaled to sum to 4 in every model", fixed = TRUE)
  expect_output(print(.R2_test_budget(sigma = 1)), "coefficients in units of a fixed residual SD of 1", fixed = TRUE)

  # the flat allocation is the default
  flat <- .R2_test_budget(allocation = NULL)
  expect_equal(unname(flat$alpha), c(1, 1, 1))
  expect_output(print(flat), "allocation: Dirichlet(1, 1, 1), retained terms keep their concentrations", fixed = TRUE)
  expect_equal(unname(.R2_test_budget(terms = "x", allocation = NULL, contrast = NULL)$alpha), 1)
})

test_that("Beta shapes follow the cross-model conventions", {

  budget <- .R2_test_budget(total = .R2_beta(2, 3))
  expect_equal(.bt_R2_model_shapes(budget, c(1, 2), 4), list(a = 2, b = 3, alpha = c(1, 2), A = 3))

  linked <- .R2_test_budget(convention = "linked", total = .R2_beta(4, 3))
  expect_equal(.bt_R2_model_shapes(linked, c(1, 2), 4)[c("a", "b")], list(a = 3, b = 3))

  # a - A and a + b stay at their full-model values (A_ref = 4)
  conditioning <- .R2_test_budget(convention = "conditioning", total = .R2_beta(5, 3))
  shapes <- .bt_R2_model_shapes(conditioning, c(1, 2), 4)
  expect_equal(c(shapes$a - shapes$A, shapes$a + shapes$b), c(5 - 4, 5 + 3))

  # fixed_total keeps the sum of the declared concentrations (6)
  fixed <- .R2_test_budget(concentration = "fixed_total", allocation = .R2_dirichlet(2, 3, 1))
  expect_equal(.bt_R2_model_shapes(fixed, c(1, 2), NA)$alpha, c(2, 4))
  expect_equal(.bt_R2_model_shapes(fixed, 2, NA)$alpha, 6)

  # kernel laws count coded coefficients, so feasibility is checked per model
  infeasible <- .R2_test_budget(law = "normal", convention = "conditioning", total = .R2_beta(1, 1),
    terms = c("x", "f"), allocation = .R2_dirichlet(1, 1))
  expect_error(.bt_R2_model_shapes(infeasible, 1, 4), class = "BayesTools_R2_invalid")
})

test_that("JAGS_formula() binds a budget with deletion, redistribution and exact design metrics", {

  data <- .R2_test_data()
  info <- JAGS_formula(~ x + f, "mu", data,
    list(intercept = prior("normal", list(0, 1)), R2 = .R2_test_budget()))
  priors <- info$prior_list
  expect_identical(names(priors), c("mu_intercept", "mu_x", "mu_f", "mu_R2", "mu_R2_allocation"))
  expect_true(is.prior.R2(priors$mu_x) && is.prior.simple(priors$mu_x) && !is.prior.factor(priors$mu_x))
  expect_true(is.prior.R2(priors$mu_f) && is.prior.factor(priors$mu_f) && is.prior.orthonormal(priors$mu_f))
  # z is declared but absent: the budget is redistributed over x and f
  expect_identical(priors$mu_R2$distribution, "beta")
  expect_equal(unname(unlist(priors$mu_R2$parameters)), c(1, 2))
  expect_equal(priors$mu_R2_allocation$parameters$alpha, c(1, 2))
  expect_identical(c(priors$mu_x$units, priors$mu_f$units), 1:2)
  expect_identical(attr(priors$mu_x, "multiply_by"), "sigma")

  # the coefficient map whitens the centered empirical metric of each block
  for(prior in priors[c("mu_x", "mu_f")]){
    expect_equal(t(prior$map) %*% prior$metric %*% prior$map, diag(prior$d), tolerance = 1e-12)
  }
  X <- info$data$mu_data_f
  expect_equal(priors$mu_f$metric, crossprod(sweep(X, 2, colMeans(X))) / nrow(X), tolerance = 1e-14)
})

test_that("kernel laws allocate over coded coefficients and scale them by column SD", {

  data <- .R2_test_data()
  info <- JAGS_formula(~ x + f, "mu", data,
    list(intercept = prior("normal", list(0, 1)), R2 = .R2_test_budget(law = "laplace")))
  priors <- info$prior_list
  expect_equal(priors$mu_R2_allocation$parameters$alpha, c(1, 2, 2))
  expect_identical(priors$mu_f$units, 2:3)
  expect_equal(diag(priors$mu_f$map), unname(1 / sqrt(diag(priors$mu_f$metric))))
  expect_equal(priors$mu_f$map[upper.tri(priors$mu_f$map)], c(0))

  syntax <- .JAGS_add_priors.fun(priors)
  expect_match(syntax, "prior_par_z_mu_f[1] ~ ddexp(0, ", fixed = TRUE)
  expect_match(syntax, "mu_f[2] <- sqrt(mu_R2 / (1 - mu_R2) * mu_R2_allocation[3]) * ", fixed = TRUE)
})

test_that("excluded, conflicting and unsupported budget terms are handled explicitly", {

  data <- .R2_test_data()
  budget <- .R2_test_budget()

  # a spike at zero excludes a term from the model and from the budget
  info <- JAGS_formula(~ x + f, "mu", data,
    list(intercept = prior("normal", list(0, 1)), x = prior("spike", list(0)), R2 = budget))
  expect_true(is.prior.point(info$prior_list$mu_x))
  expect_null(info$prior_list$mu_R2_allocation)
  expect_identical(info$prior_list$mu_f$n_units, 1L)
  info_none <- JAGS_formula(~ x + f, "mu", data,
    list(intercept = prior("normal", list(0, 1)), x = prior_none(), R2 = budget))
  expect_true(is.prior.point(info_none$prior_list$mu_x))

  expect_error(JAGS_formula(~ x + f, "mu", data,
    list(intercept = prior("normal", list(0, 1)), x = prior("normal", list(0, 1)), R2 = budget)),
    class = "BayesTools_R2_invalid")
  # the contrast of the factor terms defaults to orthonormal
  orthonormal <- JAGS_formula(~ x + f, "mu", data,
    list(intercept = prior("normal", list(0, 1)), R2 = .R2_test_budget(contrast = NULL)))$prior_list
  expect_true(is.prior.orthonormal(orthonormal$mu_f))
  other <- .R2_test_budget(terms = "x", allocation = NULL)
  expect_error(JAGS_formula(~ x + f, "mu", data,
    list(intercept = prior("normal", list(0, 1)), R2 = budget, other = other)),
    "only one R-squared budget", fixed = TRUE)
  expect_error(JAGS_formula(~ x + f, "mu", data,
    list(intercept = prior("normal", list(0, 1)), R2 = budget, R2 = budget)),
    "distinct names", fixed = TRUE)
  expect_error(JAGS_formula(~ x + f, "mu", data,
    list(intercept = prior("normal", list(0, 1)), R2 = list(budget, other))),
    "Several R-squared budgets are separate entries", fixed = TRUE)

  # the entry name names the budget: a JAGS node name that is not a term
  expect_error(JAGS_formula(~ x + f, "mu", data, list(intercept = prior("normal", list(0, 1)), x = budget)),
    "is named like a formula term", fixed = TRUE, class = "BayesTools_R2_invalid")
  expect_error(JAGS_formula(~ x + f, "mu", data, list(intercept = prior("normal", list(0, 1)), z = budget)),
    "is named like a formula term", fixed = TRUE)
  expect_error(JAGS_formula(~ x + f, "mu", data, list(intercept = prior("normal", list(0, 1)), my_intercept = budget)),
    "cannot contain 'intercept'", fixed = TRUE)
  expect_error(JAGS_formula(~ x + f, "mu", data, list(intercept = prior("normal", list(0, 1)), "R 2" = budget)),
    "must be a valid JAGS node name", fixed = TRUE)
  expect_error(JAGS_formula(~ x + f, "mu", data, list(intercept = prior("normal", list(0, 1)),
    R2 = .R2_test_budget(terms = "f", allocation = NULL), R2_allocation = other)),
    "cannot end in", fixed = TRUE)

  # terms outside every budget keep their own priors; defaults fill only them
  info_mixed <- JAGS_formula(~ x + f + w, "mu", data,
    list(intercept = prior("normal", list(0, 1)), "__default_continuous" = prior("normal", list(0, 2)),
         R2 = budget))
  expect_false(is.prior.R2(info_mixed$prior_list$mu_w))
  expect_true(is.prior.R2(info_mixed$prior_list$mu_x))
})

test_that("several budgets are bound independently (compulsory covariates)", {

  data <- .R2_test_data()
  compulsory <- .R2_test_budget(terms = "z", allocation = NULL, contrast = NULL)
  optional <- .R2_test_budget(terms = c("x", "f"), allocation = NULL, total = .R2_beta(1.5, 1))
  info <- JAGS_formula(~ x + f + z, "mu", data,
    list(intercept = prior("normal", list(0, 1)), covariates = compulsory, effects = optional))
  priors <- info$prior_list
  expect_true(all(c("mu_covariates", "mu_effects", "mu_effects_allocation") %in% names(priors)))
  expect_null(priors$mu_covariates_allocation)
  expect_identical(priors$mu_z$total_name, "mu_covariates")
  expect_identical(priors$mu_f$total_name, "mu_effects")
  expect_equal(unname(unlist(priors$mu_effects$parameters)), c(1.5, 1))

  # the compulsory-only model has no signal from the optional budget
  info_c <- JAGS_formula(~ z, "mu", data,
    list(intercept = prior("normal", list(0, 1)), covariates = compulsory, effects = optional))
  expect_false(any(grepl("effects", names(info_c$prior_list))))
})

test_that("the conditioning convention of a kernel law uses the full declared dimensions", {

  data <- .R2_test_data()
  budget <- .R2_test_budget(law = "normal", convention = "conditioning", total = .R2_beta(6.5, 1))
  info <- JAGS_formula(~ x, "mu", data, list(intercept = prior("normal", list(0, 1)), R2 = budget))
  # full model units: x (1), f (2 coded coefficients with alpha 2), z (1) => A_ref = 6
  expect_equal(unname(unlist(info$prior_list$mu_R2$parameters)), c(6.5 - (6 - 1), 1 + (6 - 1)))
  expect_error(JAGS_formula(~ x, "mu", data, list(intercept = prior("normal", list(0, 1)),
    R2 = .R2_test_budget(law = "normal", convention = "conditioning", total = .R2_beta(4, 1)))),
    "gives the Beta shapes (-1, 6) in this model", fixed = TRUE)
  budget_missing <- .R2_test_budget(law = "normal", convention = "conditioning", total = .R2_beta(7, 1),
    terms = c("x", "f", "z", "v"), allocation = .R2_dirichlet(1, 2, 1, 1))
  expect_error(JAGS_formula(~ x, "mu", data, list(intercept = prior("normal", list(0, 1)), R2 = budget_missing)),
    "the data lack 'v'", fixed = TRUE)
})

test_that("the linked convention of a kernel law checks the full-model concentration when binding", {

  data <- .R2_test_data()
  # full model units: x (1), f (2 coded coefficients with alpha 2), z (1) => A = 6
  expect_error(JAGS_formula(~ x, "mu", data, list(intercept = prior("normal", list(0, 1)),
    R2 = .R2_test_budget(law = "normal", convention = "linked", total = .R2_beta(4, 1)))),
    "the full model, 6 for this budget 'R2'; use total = prior(\"beta\", list(alpha = 6, beta = 1))",
    fixed = TRUE, class = "BayesTools_R2_invalid")
  budget <- .R2_test_budget(law = "normal", convention = "linked", total = .R2_beta(6, 1))
  expect_equal(unname(unlist(JAGS_formula(~ x, "mu", data, list(intercept = prior("normal", list(0, 1)),
    R2 = budget))$prior_list$mu_R2$parameters)), c(1, 1))
  expect_equal(unname(unlist(JAGS_formula(~ x + f + z, "mu", data, list(intercept = prior("normal", list(0, 1)),
    R2 = budget))$prior_list$mu_R2$parameters)), c(6, 1))
  # with fixed_total, the total concentration is that of the allocation
  expect_error(.R2_test_budget(law = "normal", convention = "linked", concentration = "fixed_total",
    total = .R2_beta(6, 1)), "the full model, 4 for this budget", fixed = TRUE)
})

test_that("budgets refuse random effects and log-intercept formulas", {

  data <- .R2_test_data()
  data$g <- factor(rep(1:4, 6))
  expect_error(JAGS_formula(~ x + (1 | g), "mu", data,
    list(intercept = prior("normal", list(0, 1)), R2 = .R2_test_budget(terms = "x", allocation = NULL)),
    prior_random = prior_random(g = random_block(sd = prior("gamma", list(2, 2))))),
    class = "BayesTools_R2_unavailable")
  formula <- ~ x
  attr(formula, "log(intercept)") <- TRUE
  expect_error(JAGS_formula(formula, "mu", data,
    list(intercept = prior("normal", list(1, 1), list(0, Inf)), R2 = .R2_test_budget(terms = "x", allocation = NULL))),
    class = "BayesTools_R2_unavailable")
})

test_that("coefficient nodes reproduce the coefficient laws and the exact signal identity", {

  data <- .R2_test_data()
  for(law in c("exact", "normal", "laplace")){
    info <- JAGS_formula(~ x + f + z, "mu", data,
      list(intercept = prior("normal", list(0, 1)), R2 = .R2_test_budget(law = law)))
    priors <- info$prior_list
    draws <- .generate_prior_sample_matrix(priors, 200L, seed = 3)
    W <- draws[, "mu_R2"] / (1 - draws[, "mu_R2"])
    for(term in c("mu_x", "mu_f", "mu_z")){
      prior <- priors[[term]]
      beta <- draws[, .bt_R2_coefficient_names(prior, term), drop = FALSE]
      z <- draws[, .bt_R2_direction_coordinates(prior, term), drop = FALSE]
      phi <- draws[, paste0("mu_R2_allocation[", prior$units, "]"), drop = FALSE]
      expected <- if(law == "exact"){
        sqrt(W * phi[, 1]) * (z / sqrt(rowSums(z^2))) %*% t(prior$map)
      }else{
        sqrt(W * phi) * sweep(z, 2, diag(prior$map), "*")
      }
      expect_equal(unname(beta), unname(expected), tolerance = 1e-13)
      if(law == "exact"){
        signal <- rowSums((beta %*% prior$metric) * beta)
        expect_equal(signal, W * phi[, 1], tolerance = 1e-10)
      }
    }
    # orthogonal blocks: under the exact law the realized R-squared equals r
    if(law == "exact"){
      fitted <- draws[, "mu_x"] %o% info$data$mu_data_x +
        draws[, "mu_z"] %o% info$data$mu_data_z +
        draws[, c("mu_f[1]", "mu_f[2]")] %*% t(info$data$mu_data_f)
      signal <- apply(fitted, 1, function(row) mean((row - mean(row))^2))
      expect_equal(signal / (1 + signal), unname(draws[, "mu_R2"]), tolerance = 1e-10)
    }
    syntax <- .JAGS_add_priors.fun(priors)
    for(term in c("mu_x", "mu_f", "mu_z")){
      node <- .bt_dnode_R2_coefficient(term, priors[[term]])
      expect_true(all(vapply(.bt_deterministic_node_emit(node), grepl, logical(1), x = syntax, fixed = TRUE)))
    }
  }
})

test_that("factor terms are excluded by scalar or factor zero spikes", {

  data <- .R2_test_data()
  budget <- .R2_test_budget()
  for(exclusion in list(prior("spike", list(0)), prior_none(), prior_factor("spike", list(0), contrast = "orthonormal"))){
    info <- JAGS_formula(~ x + f, "mu", data,
      list(intercept = prior("normal", list(0, 1)), f = exclusion, R2 = budget))
    expect_true(is.prior.point(info$prior_list$mu_f) && is.prior.factor(info$prior_list$mu_f))
    expect_true(is.prior.orthonormal(info$prior_list$mu_f))
    expect_identical(info$prior_list$mu_x$n_units, 1L)
  }
  # a scalar zero spike excludes a factor term with the budget's contrast
  # (orthonormal by default)
  default_contrast <- JAGS_formula(~ x + f, "mu", data,
    list(intercept = prior("normal", list(0, 1)), f = prior("spike", list(0)), R2 = .R2_test_budget(contrast = NULL)))
  expect_true(is.prior.orthonormal(default_contrast$prior_list$mu_f))
})

test_that("the coefficient law does not depend on the units of the predictors", {

  data <- .R2_test_data()
  small <- data
  small$x <- small$x * 1e-9
  priors <- list(intercept = prior("normal", list(0, 1)), R2 = .R2_test_budget())
  original <- JAGS_formula(~ x + f, "mu", data, priors)$prior_list$mu_x
  rescaled <- JAGS_formula(~ x + f, "mu", small, priors)$prior_list$mu_x
  expect_equal(rescaled$map * 1e-9, original$map, tolerance = 1e-12)

  constant <- data
  constant$x <- 1
  expect_error(JAGS_formula(~ x + f, "mu", constant, priors), "is constant after centering", fixed = TRUE)
})

test_that("coefficient nodes replay from exactly their declared dependencies", {

  data <- .R2_test_data()
  info <- JAGS_formula(~ x + f, "mu", data,
    list(intercept = prior("normal", list(0, 1)), R2 = .R2_test_budget(law = "normal")))
  priors <- info$prior_list
  draws <- .generate_prior_sample_matrix(priors, 5L, seed = 7)
  for(term in c("mu_x", "mu_f")){
    node <- .bt_dnode_R2_coefficient(term, priors[[term]])
    expect_true(all(grepl("prior_par_eta_mu_R2_allocation", node$dependencies[grepl("allocation", node$dependencies)])))
    reduced <- draws[, node$dependencies, drop = FALSE]
    values <- .bt_deterministic_node_evaluate(node, .bt_deterministic_lookup(reduced, priors))
    expect_equal(unname(values), unname(draws[, node$coordinates, drop = FALSE]), tolerance = 1e-13)

    # changed gamma coordinates with stale normalized weights: the replay uses the coordinates
    changed <- draws
    eta <- grep("prior_par_eta", colnames(changed))
    changed[, eta] <- changed[, eta] * rep(c(9, 1, 1), each = nrow(changed))
    replay <- .bt_deterministic_node_evaluate(node, .bt_deterministic_lookup(changed, priors))
    phi <- changed[, eta] / rowSums(changed[, eta])
    W <- changed[, "mu_R2"] / (1 - changed[, "mu_R2"])
    z <- changed[, .bt_R2_direction_coordinates(priors[[term]], term), drop = FALSE]
    expected <- sqrt(W * phi[, priors[[term]]$units, drop = FALSE]) * sweep(z, 2, diag(priors[[term]]$map), "*")
    expect_equal(unname(replay), unname(expected), tolerance = 1e-13)
  }
})

test_that("a budget without terms takes the formula's terms, equivalently to explicit terms", {

  data <- .R2_test_data()
  bind <- function(prior_list, formula = ~ x + z + f){
    JAGS_formula(formula, "mu", data, c(list(intercept = prior("normal", list(0, 1))), prior_list))$prior_list
  }
  syntax <- function(priors) .JAGS_add_priors.fun(priors)

  # every term without its own prior; a zero spike keeps a term, excluded
  implicit <- bind(list(x = prior("normal", list(0, 1)), z = prior("spike", list(0)),
                        R2 = prior_R2(allocation = R2_dirichlet(c(f = 2, z = 1)))))
  explicit <- bind(list(x = prior("normal", list(0, 1)), z = prior("spike", list(0)),
                        R2 = prior_R2(terms = c("z", "f"), allocation = R2_dirichlet(c(1, 2)))))
  expect_identical(syntax(implicit), syntax(explicit))
  expect_identical(implicit$mu_f$budget$terms, c("z", "f"))
  expect_equal(implicit$mu_f$budget$alpha, explicit$mu_f$budget$alpha)
  expect_true(is.prior.point(implicit$mu_z))
  expect_match(print(implicit$mu_f, silent = TRUE), "(budget 'R2': z, f)", fixed = TRUE)
  expect_output(print(implicit$mu_f$budget), "z, f (resolved from the formula)", fixed = TRUE)

  # the allocation forms are equivalent; a budget without terms needs named
  # or single concentrations, and its named lists are checked at binding
  flat <- syntax(bind(list(R2 = prior_R2())))
  expect_identical(flat, syntax(bind(list(R2 = prior_R2(allocation = R2_dirichlet(1))))))
  expect_identical(flat, syntax(bind(list(R2 = prior_R2(terms = c("x", "z", "f"), allocation = R2_dirichlet(c(1, 1, 1)))))))
  expect_error(bind(list(R2 = prior_R2(allocation = R2_dirichlet(c(1, 2, 3))))), "need the explicit 'terms'", fixed = TRUE)
  expect_error(bind(list(R2 = prior_R2(allocation = R2_dirichlet(c(x = 1, z = 2))))), "each once", fixed = TRUE)
  expect_error(bind(list(R2 = prior_R2(inclusion = list(w = prior("spike", list(0.5)))))),
    "'inclusion' names terms that are not terms of the budget: 'w'", fixed = TRUE)
  gated <- bind(list(R2 = prior_R2(inclusion = list(z = prior("spike", list(0.5))))))
  expect_identical(unname(gated$mu_R2$term_gate), c(0L, 1L, 0L))

  # a budget without terms takes the terms that no other budget names; at most
  # one, and not an empty one
  two <- bind(list(cov = prior_R2(terms = "x"), R2 = prior_R2()))
  expect_identical(two$mu_z$budget$terms, c("z", "f"))
  expect_identical(two$mu_x$budget$name, "cov")
  expect_error(bind(list(R2 = prior_R2(), other = prior_R2())), "At most one R-squared budget", fixed = TRUE)
  expect_error(bind(list(x = prior("normal", list(0, 1)), R2 = prior_R2()), ~ x), "has no terms", fixed = TRUE)

  # linked: the default total names the required first shape at binding
  expect_error(bind(list(R2 = prior_R2(model = R2_model(convention = "linked")))),
    "must equal the total concentration of the full model, 3 for this budget 'R2'", fixed = TRUE)
})
