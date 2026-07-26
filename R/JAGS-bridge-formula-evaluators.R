.bt_JAGS_bridge_compile_formula_parameter_evaluator <- function(formula_list,
                                                                formula_data_list,
                                                                formula_prior_list,
                                                                formula_design_list,
                                                                model_data){

  if(length(formula_prior_list) == 0L){
    return(list(parameters = function(samples, prior_list_parameters,
                                      formula_prior_parameters = list()) list()))
  }

  fixed_plans <- list()
  random_plans <- list()

  for(parameter in names(formula_prior_list)){
    formula_parameter <- if(!is.null(formula_list)) formula_list[[parameter]] else NULL
    formula_data <- if(!is.null(formula_data_list)) formula_data_list[[parameter]] else NULL
    log_intercept <- if(!is.null(formula_parameter)){
      isTRUE(attr(formula_parameter, "log(intercept)"))
    }else{
      FALSE
    }
    parameter_prior_list <- formula_prior_list[[parameter]]
    design <- if(!is.null(formula_design_list)) formula_design_list[[parameter]] else NULL
    log_intercept <- isTRUE(log_intercept) || isTRUE(design$log_intercept)
    if(log_intercept){
      .bt_validate_formula_log_intercept_prior(
        parameter_prior_list,
        parameter = parameter
      )
    }

    if(.bt_formula_design_has_any_random_effects(design)){
      fixed_prior_list <- .bt_JAGS_marglik_formula_fixed_priors(
        parameter_prior_list,
        parameter
      )
    }else{
      fixed_prior_list <- parameter_prior_list
    }
    if(.bt_formula_design_has_sampled_random_effects(design)){
      design <- .bt_JAGS_bridge_prepare_random_effect_allocation_design(design)
      source_data <- .bt_JAGS_marglik_parameter_source_data(
        model_data = model_data,
        formula_data = formula_data,
        design = design
      )
      random_plans[[parameter]] <- list(
        parameter = parameter,
        design = design,
        formula_prior_list = parameter_prior_list,
        source_data = source_data
      )
    }

    fixed_plans[[parameter]] <- .bt_JAGS_bridge_compile_formula_fixed_plan(
      parameter = parameter,
      formula_data = formula_data,
      formula_prior_list = fixed_prior_list,
      design = design,
      log_intercept = log_intercept
    )
  }

  list(
    parameters = function(samples, prior_list_parameters,
                          formula_prior_parameters = list()){
      parameters <- list()
      for(parameter in names(fixed_plans)){
        parameters[[parameter]] <- fixed_plans[[parameter]]$value(
          samples = samples,
          prior_list_parameters = prior_list_parameters
        )
      }

      if(length(random_plans) > 0L){
        source_base <- .bt_JAGS_bridge_formula_source_base_parameters(
          samples = samples,
          prior_list_parameters = prior_list_parameters,
          formula_prior_parameters = formula_prior_parameters
        )
        for(random_plan in random_plans){
          source_parameters <- .bt_JAGS_bridge_formula_source_parameters(
            source_base = source_base,
            formula_parameters = parameters
          )
          parameters[[random_plan$parameter]] <- parameters[[random_plan$parameter]] +
            .bt_JAGS_marglik_random_effects_value(
              samples = samples,
              design = random_plan$design,
              formula_prior_list = random_plan$formula_prior_list,
              data = random_plan$source_data,
              parameters = source_parameters
            )
        }
      }

      parameters
    }
  )
}

.bt_JAGS_bridge_prepare_random_effect_allocation_design <- function(design){

  if(!.bt_formula_design_has_sampled_random_effects(design)){
    return(design)
  }

  for(random_i in seq_along(design$random_effects)){
    if(!identical(.bt_random_effect_term_compile_mode(design$random_effects[[random_i]]),
                  "sampled")){
      next
    }
    design$random_effects[[random_i]] <-
      .bt_JAGS_bridge_prepare_random_effect_allocation_term(
        design$random_effects[[random_i]]
      )
  }

  design
}

.bt_JAGS_bridge_prepare_random_effect_allocation_term <- function(random_term){

  binding <- random_term$sd_binding
  if(is.null(binding)){
    return(random_term)
  }

  .bt_check_random_sd_binding(binding)
  sd_component_allocation <- isTRUE(binding$true_allocation) &&
    length(binding$allocations) > 0L &&
    is.list(binding$allocations[[1L]]) &&
    identical(binding$allocations[[1L]]$target, "sd_component")

  if(!isTRUE(sd_component_allocation)){
    binding$factors <- .bt_random_effect_allocation_factors_with_plan(
      binding$factors
    )
    if(is.list(binding$factors_by_column) && length(binding$factors_by_column) > 0L){
      for(column in seq_along(binding$factors_by_column)){
        binding$factors_by_column[[column]] <-
          .bt_random_effect_allocation_factors_with_plan(
            binding$factors_by_column[[column]]
          )
      }
    }
  }
  if(is.list(binding$allocations) && length(binding$allocations) > 0L){
    for(allocation_i in seq_along(binding$allocations)){
      allocation <- binding$allocations[[allocation_i]]
      if(is.list(allocation$factors) &&
         !identical(allocation$target, "sd_component")){
        allocation$factors <- .bt_random_effect_allocation_factors_with_plan(
          allocation$factors
        )
      }
      if(is.list(allocation$parent_factors) &&
         !identical(allocation$target, "sd_component")){
        allocation$parent_factors <- .bt_random_effect_allocation_factors_with_plan(
          allocation$parent_factors
        )
      }
      binding$allocations[[allocation_i]] <- allocation
    }
  }

  random_term$sd_binding <- binding
  random_term
}

.bt_JAGS_bridge_compile_formula_fixed_plan <- function(parameter,
                                                       formula_data,
                                                       formula_prior_list,
                                                       design,
                                                       log_intercept){

  force(parameter)
  force(formula_data)
  force(formula_prior_list)
  force(log_intercept)

  if(.bt_JAGS_formula_design_can_reconstruct(design)){
    return(.bt_JAGS_bridge_compile_formula_design_plan(
      design = design,
      formula_prior_list = formula_prior_list,
      log_intercept = log_intercept
    ))
  }

  list(
    value = function(samples, prior_list_parameters){
      .JAGS_marglik_parameters_formula_get(
        samples = samples,
        parameter = parameter,
        formula_data_list = formula_data,
        formula_prior_list = formula_prior_list,
        prior_list_parameters = prior_list_parameters,
        log_intercept = log_intercept
      )
    }
  )
}

.bt_JAGS_bridge_compile_formula_design_plan <- function(design,
                                                        formula_prior_list,
                                                        log_intercept){

  force(design)
  force(formula_prior_list)

  parameter <- design$parameter
  n_rows <- nrow(design$model_matrix)
  intercept_name <- paste0(parameter, "_intercept")
  intercept_plan <- NULL
  if(intercept_name %in% names(formula_prior_list)){
    intercept_plan <- .bt_JAGS_bridge_compile_formula_intercept_plan(
      prior_object = formula_prior_list[[intercept_name]],
      parameter_name = intercept_name,
      log_intercept = isTRUE(log_intercept) || isTRUE(design$log_intercept)
    )
  }

  term_names <- setdiff(names(formula_prior_list), intercept_name)
  term_plans <- lapply(term_names, function(term_name){
    term_prior <- formula_prior_list[[term_name]]
    model_term <- sub(paste0("^", parameter, "_"), "", term_name)
    columns <- .bt_JAGS_formula_design_term_columns(design, model_term)
    .bt_JAGS_bridge_compile_formula_term_plan(
      term_name = term_name,
      term_prior = term_prior,
      term_data = design$model_matrix[, columns, drop = FALSE]
    )
  })

  list(
    value = function(samples, prior_list_parameters){
      output <- rep(0, n_rows)
      if(!is.null(intercept_plan)){
        output <- output + intercept_plan$value(
          samples = samples,
          prior_list_parameters = prior_list_parameters
        )
      }
      for(term_plan in term_plans){
        output <- output + term_plan$value(
          samples = samples,
          prior_list_parameters = prior_list_parameters
        )
      }
      as.vector(output)
    }
  )
}

.bt_JAGS_bridge_compile_formula_intercept_plan <- function(prior_object,
                                                           parameter_name,
                                                           log_intercept){

  force(log_intercept)
  .bt_validate_formula_reconstruction_prior(
    prior_object,
    parameter_name
  )
  value_evaluator <- .bt_JAGS_bridge_compile_parameter_values(
    prior_object = prior_object,
    parameter_names = parameter_name
  )
  multiply_by_evaluator <- .bt_JAGS_bridge_compile_prior_multiply_by(prior_object)

  list(
    value = function(samples, prior_list_parameters){
      value <- value_evaluator(samples)
      if(isTRUE(log_intercept)){
        value <- log(value)
      }
      multiply_by_evaluator(prior_list_parameters) * value
    }
  )
}

.bt_JAGS_bridge_compile_formula_term_plan <- function(term_name,
                                                      term_prior,
                                                      term_data){

  force(term_name)
  force(term_prior)
  force(term_data)

  .bt_validate_formula_reconstruction_prior(term_prior, term_name)
  multiply_by_evaluator <- .bt_JAGS_bridge_compile_prior_multiply_by(term_prior)

  if(is.prior.point(term_prior) && !is.prior.factor(term_prior)){
    location <- term_prior[["parameters"]][["location"]]
    term_vector <- as.vector(term_data)
    return(list(
      value = function(samples, prior_list_parameters){
        multiply_by_evaluator(prior_list_parameters) * location * term_vector
      }
    ))
  }

  if(is.prior.point(term_prior) && is.prior.factor(term_prior)){
    location <- term_prior[["parameters"]][["location"]]
    levels <- .get_prior_factor_levels(term_prior)
    if(levels == 1L){
      term_vector <- as.vector(term_data)
      return(list(
        value = function(samples, prior_list_parameters){
          multiply_by_evaluator(prior_list_parameters) * location * term_vector
        }
      ))
    }

    term_values <- rep(location, levels)
    return(list(
      value = function(samples, prior_list_parameters){
        multiply_by_evaluator(prior_list_parameters) *
          as.vector(term_data %*% term_values)
      }
    ))
  }

  if(is.prior.factor(term_prior)){
    levels <- .get_prior_factor_levels(term_prior)
    if(levels == 1L){
      value_evaluator <- .bt_JAGS_bridge_compile_parameter_values(
        prior_object = term_prior,
        parameter_names = term_name
      )
      term_vector <- as.vector(term_data)
      return(list(
        value = function(samples, prior_list_parameters){
          multiply_by_evaluator(prior_list_parameters) *
            value_evaluator(samples) * term_vector
        }
      ))
    }

    parameter_names <- paste0(term_name, "[", seq_len(levels), "]")
    value_evaluator <- .bt_JAGS_bridge_compile_parameter_values(
      prior_object = term_prior,
      parameter_names = parameter_names
    )
    return(list(
      value = function(samples, prior_list_parameters){
        multiply_by_evaluator(prior_list_parameters) *
          as.vector(term_data %*% value_evaluator(samples))
      }
    ))
  }

  if(is.prior.simple(term_prior)){
    value_evaluator <- .bt_JAGS_bridge_compile_parameter_values(
      prior_object = term_prior,
      parameter_names = term_name
    )
    term_vector <- as.vector(term_data)
    return(list(
      value = function(samples, prior_list_parameters){
        multiply_by_evaluator(prior_list_parameters) *
          value_evaluator(samples) * term_vector
      }
    ))
  }

  stop(
    "Internal formula reconstruction prior dispatch failed for '",
    term_name, "'.",
    call. = FALSE
  )
}

.bt_JAGS_bridge_compile_prior_multiply_by <- function(prior_object){

  multiply_by <- attr(prior_object, "multiply_by")
  if(is.null(multiply_by)){
    return(function(prior_list_parameters) 1)
  }
  if(is.numeric(multiply_by)){
    force(multiply_by)
    return(function(prior_list_parameters) multiply_by)
  }

  force(multiply_by)
  function(prior_list_parameters){
    .bt_JAGS_marglik_resolve_named_multiply_by(
      multiply_by = multiply_by,
      prior_list_parameters = prior_list_parameters
    )
  }
}

.bt_JAGS_bridge_formula_source_base_parameters <- function(samples,
                                                           prior_list_parameters,
                                                           formula_prior_parameters = list()){

  out <- as.list(samples)
  if(length(prior_list_parameters) > 0L){
    out[names(prior_list_parameters)] <- prior_list_parameters
  }
  if(length(formula_prior_parameters) > 0L){
    out[names(formula_prior_parameters)] <- formula_prior_parameters
  }

  out
}

.bt_JAGS_bridge_formula_source_parameters <- function(source_base,
                                                      formula_parameters){

  out <- source_base
  if(length(formula_parameters) > 0L){
    out[names(formula_parameters)] <- formula_parameters
  }

  out
}
