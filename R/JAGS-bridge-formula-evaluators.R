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
  state_constants <- .bt_formula_bridge_state_constants(formula_design_list)

  for(parameter in names(formula_prior_list)){
    formula_parameter <- if(!is.null(formula_list)) formula_list[[parameter]] else NULL
    if(!is.null(formula_parameter)){
      .bt_validate_formula_replay_grammar(formula_parameter)
    }
    formula_data <- if(!is.null(formula_data_list)) formula_data_list[[parameter]] else NULL
    log_intercept <- if(!is.null(formula_parameter)){
      isTRUE(attr(formula_parameter, "log(intercept)"))
    }else{
      FALSE
    }
    parameter_prior_list <- formula_prior_list[[parameter]]
    design <- if(!is.null(formula_design_list)) formula_design_list[[parameter]] else NULL
    if(length(design$expression_specs) > 0L ||
       length(design$transformed_terms) > 0L){
      .bt_JAGS_marglik_parameter_source_data(
        model_data = model_data,
        formula_data = formula_data,
        design = design
      )
    }
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
        source_data = source_data,
        value_plans = lapply(.bt_formula_design_sampled_random_effects(design),
          .bt_JAGS_bridge_compile_random_value_plan,
          prior_list = parameter_prior_list)
      )
    }

    fixed_plans[[parameter]] <- .bt_JAGS_bridge_compile_formula_fixed_plan(
      parameter = parameter,
      formula_prior_list = fixed_prior_list,
      design = design,
      log_intercept = log_intercept
    )
  }

  list(
    parameters = function(samples, prior_list_parameters,
                          formula_prior_parameters = list()){
      parameters <- list()
      if(length(formula_prior_parameters) == 0L){
        formula_prior_parameters <- .bt_formula_decode_source_parameters(samples,
          formula_prior_list, formula_design_list)
      }
      source_base <- .bt_JAGS_bridge_formula_source_base_parameters(samples,
        prior_list_parameters, formula_prior_parameters, state_constants)
      for(parameter in names(fixed_plans)){
        parameters[[parameter]] <- fixed_plans[[parameter]]$value(
          samples = samples,
          prior_list_parameters = source_base
        )
      }

      if(length(random_plans) > 0L){
        for(random_plan in random_plans){
          source_parameters <- .bt_JAGS_bridge_formula_source_parameters(
            source_base = source_base,
            formula_parameters = parameters
          )
          source_parameters <- .bt_parameter_source_forbid_formula_parameters(
            source_parameters,
            names(random_plans)
          )
          parameters[[random_plan$parameter]] <- parameters[[random_plan$parameter]] +
            .bt_JAGS_marglik_random_effects_value(
              samples = samples,
              design = random_plan$design,
              formula_prior_list = random_plan$formula_prior_list,
              data = random_plan$source_data,
              parameters = source_parameters,
              value_plans = random_plan$value_plans
            )
        }
      }

      parameters
    }
  )
}

.bt_JAGS_bridge_compile_random_value_plan <- function(random_term, prior_list){

  if(.bt_random_effect_has_row_indexed_external_sd(random_term) ||
     inherits(random_term$latent_layout,
              "BayesTools_random_effect_structured_local_layout")){
    return(NULL)
  }
  structure <- .bt_random_effect_structure(
    random_term,
    context = "Bridge-sampling random-effect metadata"
  )
  if(random_term$n_columns > 1L &&
     structure %in% c("cs", "hcs", "ar1", "car", "har")){
    return(NULL)
  }
  sd_evaluator <- .bt_JAGS_bridge_compile_random_sd_evaluator(random_term, prior_list)
  if(is.null(sd_evaluator)){
    return(NULL)
  }
  list(
    structure = structure,
    sd_evaluator = sd_evaluator,
    latent_names = .bt_random_effect_latent_names(
      random_term = random_term,
      n_groups = random_term$n_groups,
      n_columns = random_term$n_columns
    )
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
                                                       formula_prior_list,
                                                       design,
                                                       log_intercept){

  if(!.bt_JAGS_formula_design_can_reconstruct(design)){
    .bt_stop_refit_required(
      "JAGS_bridgesampling() cannot reconstruct formula parameter '", parameter,
      "' because its fitted formula-design metadata are missing. Refit the ",
      "model with this version of BayesTools."
    )
  }

  .bt_JAGS_bridge_compile_formula_design_plan(
    design = design,
    formula_prior_list = formula_prior_list,
    log_intercept = log_intercept
  )
}

# The fixed linear predictor of a fitted design is the registered
# 'linear_predictor' node's fixed part.
.bt_JAGS_bridge_compile_formula_design_plan <- function(design,
                                                        formula_prior_list,
                                                        log_intercept){

  .bt_dnode_linear_predictor_fixed_plan(
    design = design,
    formula_prior_list = formula_prior_list,
    log_intercept = log_intercept
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
                                                           formula_prior_parameters = list(),
                                                           state_constants = list()){

  out <- as.list(samples)
  for(name in intersect(names(prior_list_parameters), names(formula_prior_parameters))){
    if(!identical(prior_list_parameters[[name]], formula_prior_parameters[[name]])){
      .bt_formula_transform_stop("Natural formula state has contradictory owners.", reason = "contradictory_state_constants", state = name)
    }
  }
  if(length(prior_list_parameters) > 0L){
    out[names(prior_list_parameters)] <- prior_list_parameters
  }
  if(length(formula_prior_parameters) > 0L){
    out[names(formula_prior_parameters)] <- formula_prior_parameters
  }
  for(name in names(state_constants)){
    declared <- c(prior_list_parameters, formula_prior_parameters)
    if(name %in% names(declared) && (!is.numeric(declared[[name]]) ||
       length(declared[[name]]) != 1L || !is.finite(declared[[name]]) ||
       declared[[name]] != state_constants[[name]])){
      .bt_formula_transform_stop("Natural formula state disagrees with its retained scalar declaration.", reason = "contradictory_state_constants", state = name)
    }
    out[[name]] <- state_constants[[name]]
  }

  out
}

# Direct formula helpers need only referenced natural sources and declared
# numeric points. Unused nuisance coordinates must not preempt the fixed-plan
# missing-input checks or the allocation reader's cached auxiliary route.
.bt_formula_decode_source_parameters <- function(samples, formula_prior_list,
                                                  designs){

  priors <- do.call(c, unname(formula_prior_list))
  needed <- names(priors)[vapply(priors, function(prior){
    !is.null(.bt_formula_numeric_point(prior))
  }, logical(1))]
  for(design in designs){
    needed <- union(needed, intersect(names(design$prior_list),
      paste0(design$parameter, "_", design$model_terms)))
    needed <- union(needed, .bt_formula_predictor_multiplier_dependencies(design))
    random_terms <- .bt_formula_design_random_effects(design)
    needed <- union(needed, unlist(lapply(random_terms, `[[`, "sd_parameter_names"), use.names = FALSE))
    for(term in random_terms){
      if(.bt_random_effect_has_row_indexed_external_sd(term)){
        source <- .bt_random_effect_row_indexed_source(term)
        needed <- union(needed, as.character(.bt_parameter_source_inputs(source)))
      }
    }
    owner <- design$point_expression_owner
    if(!is.null(owner)) needed <- union(needed,
      unlist(lapply(owner$points, `[[`, "parameter_dependencies"), use.names = FALSE))
  }
  needed <- unique(sub("\\[[0-9]+\\]$", "", needed[!is.na(needed)]))
  priors <- priors[intersect(needed, names(priors))]
  sample_names <- names(samples)
  priors <- Filter(function(prior){
    !.is_prior_expression(prior)
  }, priors)
  available <- vapply(names(priors), function(name){
    prior <- priors[[name]]
    if(!is.null(.bt_formula_numeric_point(prior))) return(TRUE)
    if(is.prior.simplex(prior) && identical(prior$distribution, "dirichlet")){
      eta_names <- paste0(.JAGS_prior_dirichlet_eta_name(name),
        "[", seq_len(prior$parameters[["K"]]), "]")
      return(all(eta_names %in% sample_names))
    }
    all(JAGS_to_monitor(stats::setNames(list(prior), name)) %in% sample_names)
  }, logical(1))
  JAGS_marglik_parameters(samples, priors[available])
}

.bt_JAGS_bridge_formula_source_parameters <- function(source_base,
                                                      formula_parameters){

  out <- source_base
  if(length(formula_parameters) > 0L){
    out[names(formula_parameters)] <- formula_parameters
  }

  out
}
