skip_if_not_test_profile("unit")

.point_replay_carrier <- function(location = expression(2 * theta),
                                  draws = cbind(theta = c(.25, .75), mu_x = c(.5, 1.5)),
                                  parents = list(theta = prior("normal", list(0, 1))),
                                  model_data = NULL){
  attr(draws, "prior_list") <- parents
  JAGS_formula_draws(draws, ~ x, "mu", data.frame(x = c(-2, -1, 1, 2)),
    list(intercept = prior("point", list(0)), x = prior("point", list(location = location))),
    model_data = model_data)
}

.point_replay_fit <- function(carrier){
  fit <- coda::mcmc.list(coda::mcmc(as.matrix(carrier)))
  class(fit) <- c("BayesTools_fit", class(fit))
  for(name in c("prior_list", "formula_design", "formula_scale")) attr(fit, name) <- attr(carrier, name)
  fit <- .bt_attach_parameter_map(fit)
  fit <- .bt_attach_draw_geometry(fit)
  .bt_attach_fit_contract(fit)
}

test_that("point coefficients prefer valid monitors and replay declared parents when absent", {
  monitored <- .point_replay_carrier()
  expect_equal(unname(JAGS_evaluate_formula(monitored, parameter = "mu")), outer(c(-2, -1, 1, 2), c(.5, 1.5)))
  parent_only <- .point_replay_carrier(draws = cbind(theta = c(.25, .75)))
  expect_equal(JAGS_evaluate_formula(parent_only, parameter = "mu"), JAGS_evaluate_formula(monitored, parameter = "mu"))
  monitor_only <- .point_replay_carrier(draws = cbind(mu_x = c(.5, 1.5)))
  expect_equal(JAGS_evaluate_formula(monitor_only, parameter = "mu"), JAGS_evaluate_formula(monitored, parameter = "mu"))
  fit <- .point_replay_fit(monitored)
  stale <- cbind(theta = c(1, 2), mu_x = c(.5, 1.5), "mu[1]" = 99, "mu[2]" = 99, "mu[3]" = 99, "mu[4]" = 99)
  expect_equal(as.vector(JAGS_evaluate_formula(fit, parameter = "mu")), as.vector(JAGS_evaluate_formula(monitored, parameter = "mu")))
  rebuilt <- JAGS_evaluate_deterministic(fit, stale, nodes = c("mu_x", "mu"))
  expect_equal(as.numeric(rebuilt[, "mu_x"]), c(2, 4))
  expect_equal(unname(rebuilt[, paste0("mu[", 1:4, "]")]), outer(c(2, 4), c(-2, -1, 1, 2)))
  expect_equal(unname(JAGS_evaluate_deterministic(fit, as.matrix(monitor_only), nodes = "mu")), outer(c(.5, 1.5), c(-2, -1, 1, 2)))
  expect_error(JAGS_evaluate_deterministic(fit, as.matrix(monitor_only), nodes = "mu_x"), class = "BayesTools_formula_point_unavailable")
  table <- JAGS_deterministic_nodes(fit)
  expect_identical(table$dependencies[[match("mu_x", table$node)]], "theta")
  expect_true("mu_x" %in% table$dependencies[[match("mu", table$node)]])
  expect_identical(parameter_coordinates(fit)$convergence_role[parameter_coordinates(fit)$coordinate_name == "mu_x"], "derived")
})

test_that("point replay owns scalar data and array dimensions independently of prediction rows", {
  parents <- list(theta = prior("point", list(2)), delta = prior("point", list(location = expression(theta + 1))))
  carrier <- .point_replay_carrier(expression(delta * model_x[2] + A[2, 1]),
    draws = cbind(unrelated = c(1, 2)), parents = parents,
    model_data = list(model_x = c(7, 11, 13), A = matrix(1:6, 3, 2)))
  expected <- matrix(c(-2, -1, 1, 2) * 35, 4, 2)
  expect_equal(unname(JAGS_evaluate_formula(carrier, parameter = "mu")), expected)
  expect_equal(unname(JAGS_predict_formula(carrier, "mu", data = data.frame(x = c(3, 3, -1)), formula_target = "fixed")$value), matrix(c(105, 105, -35), 3, 2))
  expect_identical(dim(attr(carrier, "formula_design")$mu$point_expression_owner$data$A), c(3L, 2L))
  serialized <- unserialize(serialize(carrier, NULL))
  expect_identical(JAGS_evaluate_formula(serialized, parameter = "mu"), JAGS_evaluate_formula(carrier, parameter = "mu"))
  constant <- .point_replay_carrier(expression(2), cbind(unrelated = 1:2), parents = list())
  expect_equal(unname(JAGS_evaluate_formula(constant, parameter = "mu")), matrix(c(-4, -2, 2, 4), 4, 2))
  fixed <- .point_replay_fit(.point_replay_carrier(expression(2), cbind(mu_x = c(2, 2)), parents = list()))
  expect_identical(parameter_coordinates(fixed)$convergence_role[parameter_coordinates(fixed)$coordinate_name == "mu_x"], "derived")
  indexed <- .point_replay_carrier(expression(theta[2] + delta[2]),
    cbind("theta[2]" = c(2, 4)), parents = list(theta = prior("mnormal", list(mean = 0, sd = 1, K = 3)),
      delta = prior_factor_levels(prior_factor("point", list(2), contrast = "treatment"), 3)))
  expect_equal(unname(JAGS_evaluate_formula(indexed, parameter = "mu")), outer(c(-2, -1, 1, 2), c(4, 6)))
})

test_that("point source, syntax, shape, cycle, and owner refusals remain distinct", {
  missing <- .point_replay_carrier(draws = cbind(unrelated = 1:2))
  condition <- tryCatch(JAGS_evaluate_formula(missing, parameter = "mu"), error = identity)
  expect_s3_class(condition, "BayesTools_formula_point_unavailable")
  expect_identical(condition$reason, "missing_point_parent")
  for(location in list(expression(sin(theta)), expression(theta[i]))){
    monitored <- .point_replay_carrier(location)
    expect_equal(unname(JAGS_evaluate_formula(monitored, parameter = "mu")), outer(c(-2, -1, 1, 2), c(.5, 1.5)))
    fit <- .point_replay_fit(monitored)
    expect_equal(unname(JAGS_evaluate_deterministic(fit, as.matrix(monitored), nodes = "mu")), outer(c(.5, 1.5), c(-2, -1, 1, 2)))
    expect_error(JAGS_evaluate_deterministic(fit, as.matrix(monitored), nodes = "mu_x"), class = "BayesTools_formula_point_unavailable")
    unmonitored <- .point_replay_carrier(location, cbind(theta = c(.25, .75)))
    expect_error(JAGS_evaluate_formula(unmonitored, parameter = "mu"), class = "BayesTools_formula_point_unavailable")
  }
  unresolved <- .point_replay_carrier(expression(syntax_parent), cbind(unrelated = 1:2))
  condition <- tryCatch(JAGS_evaluate_formula(unresolved, parameter = "mu"), error = identity)
  expect_identical(condition$reason, "unresolved_point_parent")
  cycle <- .point_replay_carrier(expression(delta), cbind(unrelated = 1:2),
    list(delta = prior("point", list(location = expression(mu_x)))))
  condition <- tryCatch(JAGS_evaluate_formula(cycle, parameter = "mu"), error = identity)
  expect_identical(condition$reason, "cyclic_point_parent")
  for(location in list(expression(model_x), expression(model_x[0]), expression(model_x[1.5]), expression(model_x[4]), expression(sqrt(-1)))){
    carrier <- .point_replay_carrier(location, cbind(unrelated = 1:2), parents = list(), model_data = list(model_x = 1:3))
    expect_error(suppressWarnings(JAGS_evaluate_formula(carrier, parameter = "mu")), class = "BayesTools_formula_point_unavailable")
  }
  broken <- .point_replay_carrier()
  design <- attr(broken, "formula_design")
  design$mu$point_expression_owner <- NULL
  attr(broken, "formula_design") <- design
  expect_error(JAGS_evaluate_formula(broken, parameter = "mu"), class = "BayesTools_refit_required")
  expect_error(.point_replay_carrier(expression(theta), model_data = list(theta = 3)), class = "BayesTools_formula_point_unavailable")
  expect_error(.point_replay_carrier(expression(1, 2)), class = "BayesTools_formula_point_unavailable")
  unsafe <- .point_replay_carrier()
  design <- attr(unsafe, "formula_design")
  design$mu$point_expression_owner$points$mu_x$spec$parsed <- quote(system("forbidden"))
  attr(unsafe, "formula_design") <- design
  expect_error(JAGS_evaluate_formula(unsafe, parameter = "mu"), class = "BayesTools_formula_point_unavailable")
  expect_error(.bt_JAGS_bridge_compile_parameter_values(prior("point", list(location = expression(theta))), "mu_x"), class = "BayesTools_formula_point_unavailable")
  expect_error(JAGS_bridgesampling_posterior(cbind(mu_x = 1:20), list(mu_x = prior("point", list(location = expression(theta))))), "expression")
})
