skip_if_not_test_profile("unit")

.bridge_formula_validation_inputs <- function(){

  list(
    formula_list = list(mu = ~ 1, sigma = ~ 1),
    formula_data_list = list(mu = list(), sigma = list()),
    formula_prior_list = list(mu = list(), sigma = list()),
    formula_scale_list = NULL,
    formula_random_prior_list = NULL,
    formula_random_effects_compile_list = NULL
  )
}

.bridge_formula_validation_design <- function(){

  formula_data <- data.frame(.row = 1)
  out <- lapply(c("mu", "sigma"), function(parameter){
    JAGS_formula(
      formula = ~ 1,
      parameter = parameter,
      data = formula_data,
      prior_list = list(intercept = prior("point", list(0)))
    )$formula_design
  })
  stats::setNames(out, c("mu", "sigma"))
}

.bridge_formula_validation_context <- function(inputs, fitted_design = NULL){

  fit <- list()
  if(is.null(fitted_design)){
    fitted_design <- .bridge_formula_validation_design()
  }
  attr(fit, "formula_design") <- fitted_design

  do.call(
    BayesTools:::.bt_JAGS_bridge_formula_context,
    c(list(fit = fit), inputs)
  )
}

test_that("bridge formula inputs require fully named unique required lists", {

  invalid_names <- list(
    unnamed = NULL,
    missing = c(NA_character_, "sigma"),
    empty = c("", "sigma"),
    duplicate = c("mu", "mu")
  )

  for(argument in c(
    "formula_list",
    "formula_data_list",
    "formula_prior_list"
  )){
    for(case in names(invalid_names)){
      inputs <- .bridge_formula_validation_inputs()
      if(is.null(invalid_names[[case]])){
        inputs[[argument]] <- unname(inputs[[argument]])
      }else{
        names(inputs[[argument]]) <- invalid_names[[case]]
      }

      expected <- if(identical(case, "duplicate")){
        paste0(
          "The '", argument,
          "' argument must not contain duplicate names"
        )
      }else{
        paste0(
          "The '", argument,
          "' argument must be a fully named list."
        )
      }
      expect_error(
        .bridge_formula_validation_context(inputs),
        expected,
        fixed = TRUE,
        info = paste(argument, case)
      )
    }
  }
})

test_that("bridge formula inputs validate explicitly supplied optional names", {

  invalid_names <- list(
    unnamed = NULL,
    missing = c(NA_character_, "sigma"),
    empty = c("", "sigma"),
    duplicate = c("mu", "mu"),
    unrecognized = c("mu", "tau")
  )

  for(argument in c(
    "formula_scale_list",
    "formula_random_prior_list",
    "formula_random_effects_compile_list"
  )){
    for(case in names(invalid_names)){
      inputs <- .bridge_formula_validation_inputs()
      inputs[[argument]] <- list(NULL, NULL)
      if(!is.null(invalid_names[[case]])){
        names(inputs[[argument]]) <- invalid_names[[case]]
      }

      expected <- if(identical(case, "duplicate")){
        paste0(
          "The '", argument,
          "' argument must not contain duplicate names"
        )
      }else if(identical(case, "unrecognized")){
        paste0(
          "The 'tau' objects are not recognized by the '",
          argument,
          "' argument."
        )
      }else{
        paste0(
          "The '", argument,
          "' argument must be a fully named list."
        )
      }
      expect_error(
        .bridge_formula_validation_context(inputs),
        expected,
        fixed = TRUE,
        info = paste(argument, case)
      )
    }
  }
})

test_that("bridge formula inputs preserve empty optional lists", {

  inputs <- .bridge_formula_validation_inputs()
  inputs$formula_scale_list <- list()
  inputs$formula_random_prior_list <- list()
  inputs$formula_random_effects_compile_list <- list()

  rebuildability <- do.call(
    BayesTools:::.bt_JAGS_bridge_formula_inputs_rebuildability,
    inputs
  )

  expect_identical(
    rebuildability,
    list(rebuildable = TRUE, detail = NULL)
  )
})

test_that("bridge formula name errors retain fitted-design mismatch context", {

  inputs <- .bridge_formula_validation_inputs()
  names(inputs$formula_prior_list) <- c("mu", "mu")

  expect_error(
    .bridge_formula_validation_context(inputs),
    paste0(
      "JAGS_bridgesampling() supplied formula-related inputs do not fully ",
      "match the fitted formula design. First mismatch: The ",
      "'formula_prior_list' argument must not contain duplicate names"
    ),
    fixed = TRUE
  )
})

test_that("bridge formula inputs cannot replace missing fitted metadata", {

  expect_error(
    do.call(
      BayesTools:::.bt_JAGS_bridge_formula_context,
      c(list(fit = list()), .bridge_formula_validation_inputs())
    ),
    paste0(
      "the fitted formula-design metadata are missing. Refit the model with ",
      "this version of BayesTools; supplied formula inputs cannot replace ",
      "fitted replay metadata."
    ),
    fixed = TRUE
  )
})

.bridge_formula_prior_context <- function(multiply_by = NULL){

  formula_data <- data.frame(x = c(-1, 0, 1))
  formula_priors <- list(
    intercept = prior("normal", list(0, 1)),
    x = prior("normal", list(0, 1))
  )
  if(!is.null(multiply_by)){
    attr(formula_priors$x, "multiply_by") <- multiply_by
  }
  formula_output <- JAGS_formula(
    formula = ~ 1 + x,
    parameter = "mu",
    data = formula_data,
    prior_list = formula_priors
  )

  list(
    design = formula_output$formula_design,
    inputs = list(
      formula_list = list(mu = formula_output$formula),
      formula_data_list = list(mu = formula_data),
      formula_prior_list = list(mu = formula_priors),
      formula_scale_list = NULL,
      formula_random_prior_list = NULL,
      formula_random_effects_compile_list = NULL
    )
  )
}

test_that("bridge formula context compares semantic prior multiplier metadata", {

  multiplier_cases <- list(
    added = list(fitted = NULL, supplied = "sigma"),
    removed = list(fitted = "sigma", supplied = NULL),
    changed_name = list(fitted = "sigma", supplied = "tau"),
    changed_value = list(fitted = 2, supplied = 3)
  )

  for(case in names(multiplier_cases)){
    fitted <- .bridge_formula_prior_context(
      multiplier_cases[[case]]$fitted
    )
    supplied <- .bridge_formula_prior_context(
      multiplier_cases[[case]]$supplied
    )

    expect_error(
      .bridge_formula_validation_context(
        supplied$inputs,
        fitted_design = list(mu = fitted$design)
      ),
      paste0(
        "formula prior metadata differ for parameter 'mu', ",
        "prior 'mu_x'"
      ),
      fixed = TRUE,
      info = case
    )
  }
})

test_that("bridge formula context ignores nonsemantic prior attributes", {

  fitted <- .bridge_formula_prior_context("sigma")
  supplied <- .bridge_formula_prior_context("sigma")
  attr(fitted$design$prior_list$mu_x, "bridge_note") <- "fitted only"
  attr(supplied$inputs$formula_prior_list$mu$x, "bridge_note") <- "supplied only"

  expect_silent(
    context <- .bridge_formula_validation_context(
      supplied$inputs,
      fitted_design = list(mu = fitted$design)
    )
  )
  expect_s3_class(
    context$formula_design_list$mu,
    "BayesTools_formula_design"
  )
})

test_that("bridge formula context treats fitted source semantics as authoritative", {

  formula_data <- data.frame(
    id = factor(c("a", "b", "a"), levels = c("a", "b")),
    tau_factor = c(1, 2, 3)
  )
  fitted_source <- parameter_source(
    "tau",
    shape = "row",
    values = function(parameters, data, n_rows){
      data$tau_factor[seq_len(n_rows)]
    }
  )
  fitted_prior_random <- prior_random(
    id = random_block(sd_source = random_sd_source(fitted_source))
  )
  formula_priors <- list(intercept = prior("normal", list(0, 1)))
  formula <- ~ 1 + random(1 | id, name = "id", covariance = "diag")
  fitted <- JAGS_formula(
    formula = formula,
    parameter = "mu",
    data = formula_data,
    prior_list = formula_priors,
    prior_random = fitted_prior_random
  )
  fit <- list()
  attr(fit, "formula_design") <- list(mu = fitted$formula_design)
  context_args <- list(
    fit = fit,
    formula_list = list(mu = formula),
    formula_data_list = list(mu = formula_data),
    formula_prior_list = list(mu = formula_priors),
    formula_scale_list = NULL,
    formula_random_prior_list = list(mu = fitted_prior_random)
  )

  context <- do.call(
    BayesTools:::.bt_JAGS_bridge_formula_context,
    context_args
  )
  expect_identical(
    context$formula_design_list$mu,
    fitted$formula_design
  )
  expect_null(context$formula_data_list)

  changed_callback <- parameter_source(
    "tau",
    shape = "row",
    values = function(parameters, data, n_rows){
      2 * data$tau_factor[seq_len(n_rows)]
    }
  )
  callback_args <- context_args
  callback_args$formula_random_prior_list <- list(mu = prior_random(
    id = random_block(sd_source = random_sd_source(changed_callback))
  ))
  expect_error(
    do.call(
      BayesTools:::.bt_JAGS_bridge_formula_context,
      callback_args
    ),
    "scale/allocation metadata differ",
    fixed = TRUE
  )

  source_data_args <- context_args
  source_data_args$formula_data_list$mu$tau_factor <- c(3, 2, 1)
  expect_error(
    do.call(
      BayesTools:::.bt_JAGS_bridge_formula_context,
      source_data_args
    ),
    "original formula source data differ for parameter 'mu'",
    fixed = TRUE
  )
})
