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


.random_update_test_allocation_fit <- function(){

  data <- data.frame(
    id = factor(c("a", "a", "b", "b")),
    x = c(-1, 0, 1, 2)
  )
  result <- JAGS_formula(
    formula = ~ 1 + random(1 + x | id, name = "study",
                           covariance = "diag"),
    parameter = "mu",
    data = data,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      allocation = random_variance_allocation(
        name = "allocation",
        terms = "study",
        target = "sd_component",
        scale = "total_variance",
        sd_source = random_sd_source("tau"),
        weights = prior("dirichlet", list(alpha = c(1, 1)))
      )
    )
  )
  weight <- "mu__xRE_ALLOCx_allocation__weight"
  columns <- c(
    "mu_intercept", "tau",
    result$formula_design$random_effects[[1L]]$sd_parameter_names,
    paste0(weight, "[", 1:2, "]"),
    paste0("prior_par_eta_", weight, "[", 1:2, "]")
  )
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


test_that("variance allocations declare exact scalar covariance inputs", {

  fit <- .random_update_test_allocation_fit()
  total_sd <- .random_update_test_plan(fit, "random_sd_total")
  total_var <- .random_update_test_plan(fit, "random_var_total")
  proportion <- .random_update_test_plan(
    fit,
    "random_var_prop",
    "intercept"
  )
  component_sd <- .random_update_test_plan(
    fit,
    "random_sd",
    "intercept"
  )

  expect_identical(total_sd$family, "affine")
  expect_identical(total_sd$coefficient_input, "source")
  expect_identical(total_sd$coefficient_transform, list(type = "square"))
  expect_identical(total_var$coefficient_transform, list(type = "square"))
  expect_identical(proportion$update, "allocation")
  expect_identical(proportion$coefficient_input, "source")
  expect_identical(proportion$coefficient_transform, list(type = "identity"))
  expect_identical(component_sd$update, "scale")
  expect_identical(component_sd$coefficient_input, "quantity")
  expect_identical(component_sd$coefficient_transform, list(type = "square"))
})


test_that("nested aggregate allocations use their public covariance scale", {

  data <- data.frame(
    study = factor(c("s1", "s1", "s2", "s2")),
    paper = factor(c("p1", "p2", "p1", "p2")),
    drug = factor(c("a", "b", "a", "b"))
  )
  result <- JAGS_formula(
    formula = ~ 1 +
      random(1 | study, name = "study", covariance = "diag") +
      random(1 | paper, name = "paper", covariance = "diag") +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = data,
    prior_list = list(intercept = prior("normal", list(0, 1))),
    prior_random = prior_random(
      random_variance_allocation(
        name = "total_re",
        terms = c(nested = "nested", drug = "drug"),
        sd = prior("gamma", list(2, 2)),
        weights = prior("dirichlet", list(alpha = c(1, 1)))
      ),
      random_variance_allocation(
        name = "nested_split",
        parent = allocation_ref("total_re", "nested"),
        terms = c(study = "study", paper = "paper"),
        weights = prior("dirichlet", list(alpha = c(3, 2)))
      )
    )
  )
  allocation_columns <- function(name){
    weight <- paste0("mu__xRE_ALLOCx_", name, "__weight")
    c(
      paste0(weight, "[", 1:2, "]"),
      paste0("prior_par_eta_", weight, "[", 1:2, "]")
    )
  }
  columns <- c(
    "mu_intercept",
    "mu__xRE_ALLOCx_total_re__allocation_sd",
    allocation_columns("total_re"),
    allocation_columns("nested_split")
  )
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
  fit <- .bt_attach_fit_contract(fit)

  nested_sd <- .random_update_test_plan(
    fit,
    "random_sd_total",
    "nested_split"
  )
  nested_var <- .random_update_test_plan(
    fit,
    "random_var_total",
    "nested_split"
  )
  outer_sd <- .random_update_test_plan(
    fit,
    "random_sd_total",
    "total_re"
  )

  expect_identical(outer_sd$blocks, c("study", "paper", "drug"))
  expect_identical(nested_sd$blocks, c("study", "paper"))
  expect_identical(nested_sd$coefficient_input, "quantity")
  expect_identical(nested_sd$coefficient_transform, list(type = "square"))
  expect_identical(nested_var$coefficient_input, "quantity")
  expect_identical(nested_var$coefficient_transform, list(type = "identity"))
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
