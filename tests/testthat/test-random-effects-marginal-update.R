skip_if_not_test_profile("unit")

.random_update_test_fit <- function(formula, data){

  result <- JAGS_formula(
    formula = formula,
    parameter = "mu",
    data = data,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(sd = prior("gamma", list(2, 2)))
  )
  term <- result$formula_design$random_effects[[1L]]
  correlation <- term$correlation
  correlation_columns <- if(is.null(correlation)){
    character()
  }else if(identical(correlation$type, "lkj")){
    correlation$primitive_names
  }else{
    correlation$rho_name
  }
  columns <- unique(c(
    "mu_intercept",
    term$sd_parameter_names,
    correlation_columns
  ))
  draws <- matrix(
    0.5,
    nrow = 2L,
    ncol = length(columns),
    dimnames = list(NULL, columns)
  )
  fit <- coda::mcmc.list(coda::mcmc(draws))
  class(fit) <- c("BayesTools_fit", class(fit))
  attr(fit, "prior_list") <- result$prior_list
  attr(fit, "formula_design") <- list(mu = result$formula_design)
  attr(fit, "parameter_map") <- .bt_build_parameter_map(
    columns = columns,
    prior_list = result$prior_list,
    formula_design = list(mu = result$formula_design)
  )
  fit <- .bt_attach_draw_geometry(fit)
  .bt_attach_fit_contract(fit)
}


.random_update_test_plan <- function(fit, role, component = NULL){

  catalog <- parameter_catalog(fit)
  quantities <- catalog$quantities[
    catalog$quantities$role == role,
    ,
    drop = FALSE
  ]
  if(!is.null(component)){
    quantities <- quantities[quantities$component == component, , drop = FALSE]
  }
  stopifnot(nrow(quantities) == 1L)
  selection <- parameter_catalog_resolve(
    catalog,
    quantities$canonical_name[[1L]],
    "mu"
  )
  random_effects_marginal_update_plan(fit, selection)
}


test_that("random covariance updates are classified from formula metadata", {

  data <- data.frame(
    id = factor(rep(c("a", "b"), each = 3L)),
    f = factor(rep(c("x", "y", "z"), 2L))
  )
  diagonal <- .random_update_test_fit(~ 1 + diag(1 + f | id), data)
  hcs_fit  <- .random_update_test_fit(~ 1 + hcs(f | id), data)
  ar1_fit  <- .random_update_test_fit(~ 1 + ar1(f | id), data)
  ar1_two_fit <- .random_update_test_fit(
    ~ 1 + ar1(f | id),
    droplevels(data[data$f != "z", , drop = FALSE])
  )

  diagonal_sd <- .random_update_test_plan(
    diagonal,
    "random_sd",
    "intercept"
  )
  hcs_sd <- .random_update_test_plan(hcs_fit, "random_sd", "f[x]")
  hcs_cor <- .random_update_test_plan(hcs_fit, "random_correlation")
  ar1_cor <- .random_update_test_plan(ar1_fit, "random_correlation")
  ar1_two_cor <- .random_update_test_plan(
    ar1_two_fit,
    "random_correlation"
  )

  expect_identical(diagonal_sd$family, "affine")
  expect_identical(diagonal_sd$update, "column_scale")
  expect_identical(diagonal_sd$coefficient_input, "source")
  expect_identical(hcs_sd$family, "factor")
  expect_identical(hcs_cor$family, "affine")
  expect_identical(hcs_cor$update, "correlation")
  expect_identical(hcs_cor$coefficient_input, "source")
  expect_identical(ar1_cor$family, "markov")
  expect_identical(ar1_two_cor$family, "affine")
})


test_that("single direct random scale exposes its exact row covariance basis", {

  data <- data.frame(id = factor(c("a", "a", "b", "c")))
  fit <- .random_update_test_fit(~ 1 + (1 | id), data)
  plan <- .random_update_test_plan(fit, "random_sd")

  expect_identical(plan$family, "affine")
  expect_identical(plan$coefficient_transform, list(type = "square"))
  expect_equal(
    plan$invariant_covariance$update_covariance,
    outer(as.integer(data$id), as.integer(data$id), "==") * 1
  )
  expect_equal(plan$invariant_covariance$base_covariance, matrix(0, 4L, 4L))
})


test_that("shared independent-coefficient scale exposes an exact dense basis", {

  data <- data.frame(
    id = factor(c("a", "a", "b", "b")),
    x = c(-1, 0, 0.5, 2)
  )
  fit <- .random_update_test_fit(~ 1 + id(1 + x | id), data)
  plan <- .random_update_test_plan(fit, "random_sd")
  design <- attr(fit, "formula_design")$mu$random_effects[[1L]]
  same_group <- outer(as.integer(data$id), as.integer(data$id), "==") * 1
  expected <- same_group * tcrossprod(design$model_matrix)

  expect_identical(plan$family, "affine")
  expect_equal(plan$invariant_covariance$update_covariance, unname(expected))
})


test_that("known group covariance remains an exact affine scale basis", {

  data <- data.frame(id = factor(c("b", "a", "c", "b")))
  kernel <- matrix(
    c(4, 1, 0.5, 1, 9, 2, 0.5, 2, 16),
    nrow = 3L,
    byrow = TRUE,
    dimnames = list(c("a", "b", "c"), c("a", "b", "c"))
  )
  formula <- random_effects_formula(
    ~ 1 | id,
    group_covariance = random_group_covariance(kernel, scale = "none")
  )
  fit <- .random_update_test_fit(formula, data)
  plan <- .random_update_test_plan(fit, "random_sd_ratio")
  expected <- kernel[as.character(data$id), as.character(data$id)]

  expect_identical(plan$family, "affine")
  expect_equal(
    plan$invariant_covariance$update_covariance,
    unname(expected)
  )
})


test_that("LKJ correlation updates distinguish scalar and composite paths", {

  data <- data.frame(
    id = factor(rep(c("a", "b"), each = 3L)),
    x = c(-1, 0, 1, -0.5, 0.5, 2),
    z = c(0, 1, 0, 1, 0, 1)
  )
  two <- .random_update_test_fit(~ 1 + us(1 + x | id), data)
  three <- .random_update_test_fit(~ 1 + us(1 + x + z | id), data)
  two_plan <- .random_update_test_plan(two, "random_correlation")

  catalog <- parameter_catalog(three)
  correlations <- catalog$quantities[
    catalog$quantities$role == "random_correlation",
    ,
    drop = FALSE
  ]
  selection <- parameter_catalog_resolve(
    catalog,
    correlations$canonical_name[[1L]],
    "mu"
  )
  three_plan <- random_effects_marginal_update_plan(three, selection)

  expect_identical(two_plan$family, "affine")
  expect_identical(three_plan$family, "unsupported")
  expect_identical(three_plan$reason, "non_scalar_covariance_path")
})
