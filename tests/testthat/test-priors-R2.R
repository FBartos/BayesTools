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

.R2_test_budget <- function(...){
  arguments <- utils::modifyList(list(
    alpha = c(x = 1, f = 2, z = 1), a = 1, b = 2, law = "exact", convention = "common",
    concentration = "inherited", sigma = "sigma", contrast = "orthonormal"
  ), list(...), keep.null = TRUE)
  do.call(prior_R2, arguments)
}

test_that("prior_R2() requires every undecided choice and validates its arguments", {

  expect_error(prior_R2(alpha = c(x = 1), a = 1, b = 1, law = "exact", sigma = "sigma"),
    "The 'convention', 'concentration' arguments are required", fixed = TRUE)
  expect_error(.R2_test_budget(law = "student"), "law", fixed = TRUE)
  expect_error(.R2_test_budget(convention = "linked"),
    "'a' must be NULL with convention = \"linked\"", fixed = TRUE)
  expect_error(.R2_test_budget(a = NULL), "The 'a' argument is required", fixed = TRUE)
  expect_error(.R2_test_budget(concentration = "fixed_total"), "The 'total' argument is required", fixed = TRUE)
  expect_error(.R2_test_budget(total = 2), "'total' applies only", fixed = TRUE)
  expect_error(.R2_test_budget(alpha = c(1, 2)), "named numeric vector", fixed = TRUE)
  expect_error(.R2_test_budget(alpha = c("x:f" = 1, "f:x" = 1)), "duplicated: 'f:x'", fixed = TRUE)
  expect_error(.R2_test_budget(alpha = c(x = 0.001)), "alpha")
  expect_error(.R2_test_budget(sigma = -1), "sigma")
  expect_error(.R2_test_budget(name = "my_intercept"), "cannot contain 'intercept'", fixed = TRUE)
  expect_error(.R2_test_budget(convention = "conditioning", a = 1, alpha = c(x = 1, f = 2)),
    "a - sum(alpha) + min(alpha) > 0", fixed = TRUE)

  budget <- .R2_test_budget(convention = "linked", a = NULL, concentration = "fixed_total", total = 3)
  expect_true(is.R2_budget(budget))
  expect_null(budget$a)
  expect_output(print(budget), "R2 ~ Beta(A_M, 2), A_M the retained concentration", fixed = TRUE)
  expect_output(print(budget), "rescaled to total 3", fixed = TRUE)
})

test_that("Beta shapes follow the cross-model conventions", {

  budget <- .R2_test_budget(a = 2, b = 3)
  expect_equal(.bt_R2_model_shapes(budget, c(1, 2), 4), list(a = 2, b = 3, alpha = c(1, 2), A = 3))

  linked <- .R2_test_budget(convention = "linked", a = NULL, b = 3)
  expect_equal(.bt_R2_model_shapes(linked, c(1, 2), 4)[c("a", "b")], list(a = 3, b = 3))

  # a - A and a + b stay at their full-model values (A_ref = 4)
  conditioning <- .R2_test_budget(convention = "conditioning", a = 5, b = 3)
  shapes <- .bt_R2_model_shapes(conditioning, c(1, 2), 4)
  expect_equal(c(shapes$a - shapes$A, shapes$a + shapes$b), c(5 - 4, 5 + 3))

  fixed <- .R2_test_budget(concentration = "fixed_total", total = 6)
  expect_equal(.bt_R2_model_shapes(fixed, c(1, 2), NA)$alpha, c(2, 4))
  expect_equal(.bt_R2_model_shapes(fixed, 2, NA)$alpha, 6)

  # kernel laws count coded coefficients, so feasibility is checked per model
  infeasible <- .R2_test_budget(law = "normal", convention = "conditioning", a = 1, b = 1, alpha = c(x = 1, f = 1))
  expect_error(.bt_R2_model_shapes(infeasible, 1, 4), class = "BayesTools_R2_invalid")
})

test_that("JAGS_formula() binds a budget with deletion, redistribution and exact design metrics", {

  data <- .R2_test_data()
  info <- JAGS_formula(~ x + f, "mu", data,
    list(intercept = prior("normal", list(0, 1)), "__R2" = .R2_test_budget()))
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
    list(intercept = prior("normal", list(0, 1)), "__R2" = .R2_test_budget(law = "laplace")))
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
    list(intercept = prior("normal", list(0, 1)), x = prior("spike", list(0)), "__R2" = budget))
  expect_true(is.prior.point(info$prior_list$mu_x))
  expect_null(info$prior_list$mu_R2_allocation)
  expect_identical(info$prior_list$mu_f$n_units, 1L)
  info_none <- JAGS_formula(~ x + f, "mu", data,
    list(intercept = prior("normal", list(0, 1)), x = prior_none(), "__R2" = budget))
  expect_true(is.prior.point(info_none$prior_list$mu_x))

  expect_error(JAGS_formula(~ x + f, "mu", data,
    list(intercept = prior("normal", list(0, 1)), x = prior("normal", list(0, 1)), "__R2" = budget)),
    class = "BayesTools_R2_invalid")
  expect_error(JAGS_formula(~ x + f, "mu", data,
    list(intercept = prior("normal", list(0, 1)), "__R2" = .R2_test_budget(contrast = NULL))),
    "needs a 'contrast'", fixed = TRUE)
  expect_error(JAGS_formula(~ x + f, "mu", data,
    list(intercept = prior("normal", list(0, 1)), "__R2" = list(budget, .R2_test_budget(alpha = c(x = 1), name = "other")))),
    "only one R-squared budget", fixed = TRUE)
  expect_error(JAGS_formula(~ x + f, "mu", data,
    list(intercept = prior("normal", list(0, 1)), "__R2" = list(budget, budget))),
    "distinct names", fixed = TRUE)
  expect_error(JAGS_formula(~ x + f, "mu", data,
    list(intercept = prior("normal", list(0, 1)), "__R2" = "budget")),
    "must be a prior_R2() budget", fixed = TRUE)

  # terms outside every budget keep their own priors; defaults fill only them
  info_mixed <- JAGS_formula(~ x + f + w, "mu", data,
    list(intercept = prior("normal", list(0, 1)), "__default_continuous" = prior("normal", list(0, 2)),
         "__R2" = budget))
  expect_false(is.prior.R2(info_mixed$prior_list$mu_w))
  expect_true(is.prior.R2(info_mixed$prior_list$mu_x))
})

test_that("several budgets are bound independently (compulsory covariates)", {

  data <- .R2_test_data()
  compulsory <- .R2_test_budget(alpha = c(z = 1), name = "covariates", contrast = NULL)
  optional <- .R2_test_budget(alpha = c(x = 1, f = 1), a = 1.5, b = 1, name = "effects")
  info <- JAGS_formula(~ x + f + z, "mu", data,
    list(intercept = prior("normal", list(0, 1)), "__R2" = list(compulsory, optional)))
  priors <- info$prior_list
  expect_true(all(c("mu_covariates", "mu_effects", "mu_effects_allocation") %in% names(priors)))
  expect_null(priors$mu_covariates_allocation)
  expect_identical(priors$mu_z$total_name, "mu_covariates")
  expect_identical(priors$mu_f$total_name, "mu_effects")
  expect_equal(unname(unlist(priors$mu_effects$parameters)), c(1.5, 1))

  # the compulsory-only model has no signal from the optional budget
  info_c <- JAGS_formula(~ z, "mu", data,
    list(intercept = prior("normal", list(0, 1)), "__R2" = list(compulsory, optional)))
  expect_false(any(grepl("effects", names(info_c$prior_list))))
})

test_that("the conditioning convention of a kernel law uses the full declared dimensions", {

  data <- .R2_test_data()
  budget <- .R2_test_budget(law = "normal", convention = "conditioning", a = 6.5, b = 1)
  info <- JAGS_formula(~ x, "mu", data, list(intercept = prior("normal", list(0, 1)), "__R2" = budget))
  # full model units: x (1), f (2 coded coefficients with alpha 2), z (1) => A_ref = 6
  expect_equal(unname(unlist(info$prior_list$mu_R2$parameters)), c(6.5 - (6 - 1), 1 + (6 - 1)))
  expect_error(JAGS_formula(~ x, "mu", data, list(intercept = prior("normal", list(0, 1)),
    "__R2" = .R2_test_budget(law = "normal", convention = "conditioning", a = 4, b = 1))),
    "gives the Beta shapes (-1, 6) in this model", fixed = TRUE)
  budget_missing <- .R2_test_budget(law = "normal", convention = "conditioning", a = 7, b = 1,
    alpha = c(x = 1, f = 2, z = 1, v = 1))
  expect_error(JAGS_formula(~ x, "mu", data, list(intercept = prior("normal", list(0, 1)), "__R2" = budget_missing)),
    "the data lack 'v'", fixed = TRUE)
})

test_that("budgets refuse random effects and log-intercept formulas", {

  data <- .R2_test_data()
  data$g <- factor(rep(1:4, 6))
  expect_error(JAGS_formula(~ x + (1 | g), "mu", data,
    list(intercept = prior("normal", list(0, 1)), "__R2" = .R2_test_budget(alpha = c(x = 1))),
    prior_random = prior_random(g = random_block(sd = prior("gamma", list(2, 2))))),
    class = "BayesTools_R2_unavailable")
  formula <- ~ x
  attr(formula, "log(intercept)") <- TRUE
  expect_error(JAGS_formula(formula, "mu", data,
    list(intercept = prior("normal", list(1, 1), list(0, Inf)), "__R2" = .R2_test_budget(alpha = c(x = 1)))),
    class = "BayesTools_R2_unavailable")
})

test_that("coefficient nodes reproduce the coefficient laws and the exact signal identity", {

  data <- .R2_test_data()
  for(law in c("exact", "normal", "laplace")){
    info <- JAGS_formula(~ x + f + z, "mu", data,
      list(intercept = prior("normal", list(0, 1)), "__R2" = .R2_test_budget(law = law)))
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
      list(intercept = prior("normal", list(0, 1)), f = exclusion, "__R2" = budget))
    expect_true(is.prior.point(info$prior_list$mu_f) && is.prior.factor(info$prior_list$mu_f))
    expect_true(is.prior.orthonormal(info$prior_list$mu_f))
    expect_identical(info$prior_list$mu_x$n_units, 1L)
  }
  expect_error(JAGS_formula(~ x + f, "mu", data,
    list(intercept = prior("normal", list(0, 1)), f = prior("spike", list(0)),
         "__R2" = .R2_test_budget(contrast = NULL))),
    "is excluded by a scalar point prior", fixed = TRUE)
})

test_that("the coefficient law does not depend on the units of the predictors", {

  data <- .R2_test_data()
  small <- data
  small$x <- small$x * 1e-9
  priors <- list(intercept = prior("normal", list(0, 1)), "__R2" = .R2_test_budget())
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
    list(intercept = prior("normal", list(0, 1)), "__R2" = .R2_test_budget(law = "normal")))
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
