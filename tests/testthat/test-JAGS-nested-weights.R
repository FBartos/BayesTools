skip_if_not_test_profile("unit")

test_that("nested coordinates require an opt-in and a unique guarded pair", {

  priors <- nested_weight_priors()
  syntax <- "model{ y ~ dnorm(omega[3], 1) }"
  resolve <- function(prior_list = priors, model = syntax, add = NULL, on = TRUE, data = NULL){
    .bt_nested_weights_bind(prior_list, model, add, on, data)
  }
  bind <- function(...) resolve(...)$prior_list
  expect_identical(bind(on = FALSE), priors)
  plan <- .bt_nested_weights_plan(bind()$bias)
  expect_identical(plan$coarse, 4L)
  expect_identical(plan$fine, 1L)
  expect_identical(resolve()$status$branches, c(coarse = 4L, fine = 1L))
  expect_true(resolve()$status$realized)
  expect_null(resolve()$status$reason)
  expect_null(resolve(on = FALSE)$status)
  expect_null(.bt_nested_weights_plan(priors$bias))
  expect_identical(bind(nested_weight_priors(alpha = c(1, 1, 2, 1))), nested_weight_priors(alpha = c(1, 1, 2, 1)))
  changed_split <- nested_weight_priors(fine_steps = c(.05, .1, .50))
  expect_identical(bind(changed_split), changed_split)
  ambiguous <- list(bias = prior_mixture(c(as.list(priors$bias), list(priors$bias[[4L]]))))
  expect_identical(bind(ambiguous), ambiguous)
  expect_match(resolve(ambiguous)$status$reason, "multiple eligible nested pairs")
  expect_identical(bind(add = "eta_component_1[2]"), priors)
  expect_match(resolve(add = "eta_component_1[2]")$status$reason, "eta_component_1")
  expect_identical(bind(add = "std_eta_component_4"), priors)
  expect_identical(bind(model = "model{ y ~ dnorm(omega_local_component_4[2], 1) }"), priors)
  with_expression <- c(priors, list(aux = prior("point", list(expression(eta_component_1[2])))))
  expect_identical(bind(with_expression), with_expression)
  expect_match(resolve(with_expression)$status$reason, "eta_component_1")
  expression_slope <- prior_PET("normal", list(expression(eta_component_1[2]), 1))
  with_slope <- list(bias = prior_mixture(c(as.list(priors$bias), list(expression_slope))))
  expect_identical(bind(with_slope), with_slope)
  expect_match(resolve(with_slope)$status$reason, "eta_component_1")
  from_data <- c(priors, list(aux = prior("point", list(expression(location_data)))))
  expect_true(resolve(from_data, data = list(location_data = .5))$status$realized)
  expect_identical(bind(data = setNames(list(.5), plan$nodes[["a"]])), priors)
  expect_identical(bind(bind(), on = FALSE), priors)
  expect_identical(
    selection_backend_spec(bind()$bias, include_init = FALSE)$monitor,
    c("bias_indicator", "omega")
  )
})

test_that("nested initial values preserve active physical weights and the label", {

  priors <- nested_weight_priors()
  fit_priors <- .bt_nested_weights_bind(priors, "model{}", NULL, TRUE)$prior_list
  plan <- .bt_nested_weights_plan(fit_priors$bias)
  branches <- lapply(priors$bias, .selection_branch_info)
  init <- list(eta_component_1 = c(.2, 2, .3, .4),
               eta_component_4 = c(.7, .2, 2))
  expanded <- .bt_dnode_omega("omega", fit_priors$bias)
  original <- .bt_dnode_omega("omega", priors$bias)

  for(label in seq_along(branches)){
    old <- c(init, list(bias_indicator = label))
    new <- .bt_nested_weights_init(plan, old, branches)
    expect_identical(new$bias_indicator, label)
    expect_true(all(unlist(new[plan$nodes]) > 0 & unlist(new[plan$nodes]) < 1))
    expect_false(any(names(init) %in% names(new)))
    old <- unlist(old)
    new <- unlist(new)
    names(old) <- sub("(eta_component_[14])([[:digit:]]+)$", "\\1[\\2]", names(old))
    expect_equal(
      .bt_dnode_omega_evaluate(expanded, .bt_deterministic_row_lookup(new)),
      .bt_dnode_omega_evaluate(original, .bt_deterministic_row_lookup(old)),
      tolerance = 1e-14
    )
  }
})
