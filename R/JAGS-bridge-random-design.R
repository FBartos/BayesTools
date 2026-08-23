.bt_JAGS_bridge_validate_formula_random_designs <- function(fitted_formula_design,
                                                            rebuilt_formula_design){

  if(is.null(fitted_formula_design) || length(fitted_formula_design) == 0L){
    return(invisible(TRUE))
  }

  fitted_formula_design <- .bt_JAGS_bridge_formula_design_list(fitted_formula_design)
  rebuilt_formula_design <- .bt_JAGS_bridge_formula_design_list(rebuilt_formula_design)
  if(length(fitted_formula_design) == 0L){
    return(invisible(TRUE))
  }

  fitted_random_parameters <- names(fitted_formula_design)[
    vapply(fitted_formula_design, .bt_formula_design_has_any_random_effects, logical(1))
  ]
  rebuilt_random_parameters <- names(rebuilt_formula_design)[
    vapply(rebuilt_formula_design, .bt_formula_design_has_any_random_effects, logical(1))
  ]
  if(length(fitted_random_parameters) == 0L &&
     length(rebuilt_random_parameters) == 0L){
    return(invisible(TRUE))
  }
  if(!setequal(fitted_random_parameters, rebuilt_random_parameters)){
    stop(
      "JAGS_bridgesampling() rebuilt formula random-effect design does not match the fitted design. ",
      "Random-effect formula parameter(s) differ; fitted: ",
      paste(fitted_random_parameters, collapse = ", "),
      "; rebuilt: ",
      paste(rebuilt_random_parameters, collapse = ", "),
      ". Supply the same formula, data, scaling, and prior_random() metadata used to fit the model.",
      call. = FALSE
    )
  }

  for(parameter in fitted_random_parameters){
    .bt_JAGS_bridge_validate_formula_random_design(
      parameter = parameter,
      fitted = fitted_formula_design[[parameter]],
      rebuilt = rebuilt_formula_design[[parameter]]
    )
  }

  invisible(TRUE)
}

.bt_JAGS_bridge_formula_design_list <- function(formula_design){

  if(inherits(formula_design, "BayesTools_formula_design")){
    parameter <- formula_design$parameter
    if(is.null(parameter) || length(parameter) != 1L || !nzchar(parameter)){
      parameter <- ""
    }
    out <- list(formula_design)
    names(out) <- parameter
    return(out)
  }
  if(!is.list(formula_design)){
    return(list())
  }

  design_names <- names(formula_design)
  if(is.null(design_names)){
    design_names <- rep("", length(formula_design))
  }
  for(i in seq_along(formula_design)){
    if(!nzchar(design_names[i]) &&
       inherits(formula_design[[i]], "BayesTools_formula_design") &&
       !is.null(formula_design[[i]]$parameter) &&
       length(formula_design[[i]]$parameter) == 1L){
      design_names[i] <- formula_design[[i]]$parameter
    }
  }
  names(formula_design) <- design_names
  formula_design
}

.bt_JAGS_bridge_design_parameter_name <- function(design, fallback = ""){

  if(inherits(design, "BayesTools_formula_design") &&
     is.character(design$parameter) &&
     length(design$parameter) == 1L &&
     !is.na(design$parameter) &&
     nzchar(design$parameter)){
    return(design$parameter)
  }
  if(is.character(fallback) && length(fallback) == 1L &&
     !is.na(fallback) && nzchar(fallback)){
    return(fallback)
  }

  ""
}

.bt_JAGS_bridge_validate_formula_random_design <- function(parameter, fitted,
                                                          rebuilt){

  if(!identical(fitted$random_effects_interface, rebuilt$random_effects_interface)){
    .bt_JAGS_bridge_random_design_mismatch(
      parameter,
      "random-effect interface differs"
    )
  }

  fitted_all <- .bt_formula_design_random_effects(fitted)
  rebuilt_all <- .bt_formula_design_random_effects(rebuilt)
  fitted_blocks <- .bt_JAGS_bridge_random_block_names(fitted_all)
  rebuilt_blocks <- .bt_JAGS_bridge_random_block_names(rebuilt_all)
  if(!identical(fitted_blocks, rebuilt_blocks)){
    .bt_JAGS_bridge_random_design_mismatch(
      parameter,
      "random-effect block names or order differ"
    )
  }
  if(!identical(
    .bt_JAGS_bridge_random_compile_metadata(fitted),
    .bt_JAGS_bridge_random_compile_metadata(rebuilt)
  )){
    .bt_JAGS_bridge_random_design_mismatch(
      parameter,
      "random-effect compile metadata differ"
    )
  }

  for(block in fitted_blocks){
    fitted_term <- fitted_all[[match(block, fitted_blocks)]]
    rebuilt_term <- rebuilt_all[[match(block, rebuilt_blocks)]]
    .bt_JAGS_bridge_validate_random_term(
      parameter = parameter,
      block = block,
      fitted = fitted_term,
      rebuilt = rebuilt_term
    )
  }

  .bt_JAGS_bridge_validate_formula_random_compile(
    parameter = parameter,
    design = fitted,
    label = "fitted"
  )
  .bt_JAGS_bridge_validate_formula_random_compile(
    parameter = parameter,
    design = rebuilt,
    label = "rebuilt"
  )

  invisible(TRUE)
}

.bt_JAGS_bridge_validate_formula_random_compile <- function(parameter, design,
                                                            label){

  random_effects <- .bt_formula_design_random_effects(design)
  blocks <- .bt_JAGS_bridge_random_block_names(random_effects)
  modes <- .bt_random_effects_compile_modes_from_terms(random_effects)
  if(!identical(names(modes), blocks)){
    .bt_JAGS_bridge_random_design_mismatch(
      parameter,
      paste0(label, " random-effect compile metadata is missing block names")
    )
  }

  expected <- .bt_JAGS_bridge_random_compile_metadata_from_terms(random_effects)
  stored <- .bt_JAGS_bridge_random_compile_metadata_from_policy(
    design$random_effects_compile
  )
  if(!is.null(stored) &&
     !identical(expected, stored)){
    .bt_JAGS_bridge_random_design_mismatch(
      parameter,
      paste0(label, " random-effect compile metadata is inconsistent with random-effect terms")
    )
  }

  for(random_i in seq_along(random_effects)){
    .bt_JAGS_bridge_validate_random_term(
      parameter = parameter,
      block = blocks[[random_i]],
      fitted = random_effects[[random_i]],
      rebuilt = random_effects[[random_i]]
    )
  }

  invisible(TRUE)
}

.bt_JAGS_bridge_random_compile_metadata <- function(design){

  .bt_JAGS_bridge_random_compile_metadata_from_terms(
    .bt_formula_design_random_effects(design)
  )
}

.bt_JAGS_bridge_random_compile_metadata_from_terms <- function(random_effects){

  modes <- .bt_random_effects_compile_modes_from_terms(random_effects)
  blocks <- names(modes)

  list(
    sampled = blocks[modes == "sampled"],
    marginalized = blocks[modes == "marginalized"],
    mode = modes
  )
}

.bt_JAGS_bridge_random_compile_metadata_from_policy <- function(policy){

  if(is.null(policy)){
    return(NULL)
  }
  .bt_check_random_effects_compile(policy)
  if(is.null(policy$mode)){
    return(NULL)
  }

  list(
    sampled = if(is.null(policy$sampled)) character() else policy$sampled,
    marginalized = if(is.null(policy$marginalized)) character() else policy$marginalized,
    mode = policy$mode
  )
}

.bt_JAGS_bridge_random_block_names <- function(random_effects){

  if(length(random_effects) == 0L){
    return(character())
  }

  vapply(random_effects, function(random_term){
    block_name <- random_term$block_name
    if(is.null(block_name) || length(block_name) != 1L){
      return("")
    }
    as.character(block_name)
  }, character(1))
}

.bt_JAGS_bridge_validate_random_term <- function(parameter, block, fitted,
                                                rebuilt){

  if(!identical(.bt_JAGS_bridge_random_term_structure(fitted),
                .bt_JAGS_bridge_random_term_structure(rebuilt))){
    .bt_JAGS_bridge_random_design_mismatch(
      parameter,
      "random-effect covariance structure differs",
      block = block
    )
  }
  if(!identical(fitted$group_label, rebuilt$group_label)){
    .bt_JAGS_bridge_random_design_mismatch(
      parameter,
      "random-effect grouping label differs",
      block = block
    )
  }
  if(!identical(as.character(fitted$group_levels),
                as.character(rebuilt$group_levels))){
    .bt_JAGS_bridge_random_design_mismatch(
      parameter,
      "random-effect group levels or their order differ",
      block = block
    )
  }
  if(!identical(fitted$n_groups, rebuilt$n_groups) ||
     !identical(fitted$n_columns, rebuilt$n_columns)){
    .bt_JAGS_bridge_random_design_mismatch(
      parameter,
      "random-effect group or column counts differ",
      block = block
    )
  }
  if(!identical(fitted$column_names, rebuilt$column_names) ||
     !identical(fitted$raw_column_names, rebuilt$raw_column_names)){
    .bt_JAGS_bridge_random_design_mismatch(
      parameter,
      "random-effect design column names differ",
      block = block
    )
  }
  if(!identical(fitted$contrast_owner, rebuilt$contrast_owner) ||
     !identical(
       fitted$contrast_matrices,
       rebuilt$contrast_matrices
     )){
    .bt_JAGS_bridge_random_design_mismatch(
      parameter,
      "random-effect factor basis owner or concrete matrices differ",
      block = block
    )
  }
  if(!identical(dim(fitted$model_matrix), dim(rebuilt$model_matrix)) ||
     !identical(colnames(fitted$model_matrix), colnames(rebuilt$model_matrix))){
    .bt_JAGS_bridge_random_design_mismatch(
      parameter,
      "random-effect model matrix shape or columns differ",
      block = block
    )
  }
  if(!identical(
    unname(fitted$model_matrix),
    unname(rebuilt$model_matrix)
  )){
    .bt_JAGS_bridge_random_design_mismatch(
      parameter,
      "random-effect model matrix values differ",
      block = block
    )
  }
  if(!identical(as.integer(fitted$group_map),
                as.integer(rebuilt$group_map))){
    .bt_JAGS_bridge_random_design_mismatch(
      parameter,
      "random-effect group map differs",
      block = block
    )
  }
  if(!identical(
    .bt_JAGS_bridge_scale_metadata(fitted),
    .bt_JAGS_bridge_scale_metadata(rebuilt)
  )){
    .bt_JAGS_bridge_random_design_mismatch(
      parameter,
      "random-effect scale/allocation metadata differ",
      block = block
    )
  }
  if(!identical(
    .bt_JAGS_bridge_structured_index_metadata(fitted),
    .bt_JAGS_bridge_structured_index_metadata(rebuilt)
  )){
    .bt_JAGS_bridge_random_design_mismatch(
      parameter,
      "structured random-effect index metadata differ",
      block = block
    )
  }
  if(!identical(
    .bt_JAGS_bridge_car_metadata(fitted),
    .bt_JAGS_bridge_car_metadata(rebuilt)
  )){
    .bt_JAGS_bridge_random_design_mismatch(
      parameter,
      "CAR random-effect index metadata differ",
      block = block
    )
  }
  if(!identical(
    .bt_JAGS_bridge_correlation_metadata(fitted),
    .bt_JAGS_bridge_correlation_metadata(rebuilt)
  )){
    .bt_JAGS_bridge_random_design_mismatch(
      parameter,
      "random-effect correlation metadata differ",
      block = block
    )
  }
  if(!identical(fitted$parameterization_requested,
                rebuilt$parameterization_requested) ||
     !identical(fitted$parameterization_resolved,
                rebuilt$parameterization_resolved) ||
     !identical(fitted$parameterization_reason,
                rebuilt$parameterization_reason) ||
     !identical(fitted$parameterization_policy,
                rebuilt$parameterization_policy)){
    .bt_JAGS_bridge_random_design_mismatch(
      parameter,
      "random-effect parameterization differs",
      block = block
    )
  }
  if(!identical(
    .bt_JAGS_bridge_group_covariance_metadata(fitted),
    .bt_JAGS_bridge_group_covariance_metadata(rebuilt)
  )){
    .bt_JAGS_bridge_random_design_mismatch(
      parameter,
      "random-effect group covariance metadata differ",
      block = block
    )
  }
  if(!identical(
    .bt_JAGS_bridge_latent_layout_metadata(fitted),
    .bt_JAGS_bridge_latent_layout_metadata(rebuilt)
  )){
    .bt_JAGS_bridge_random_design_mismatch(
      parameter,
      "random-effect latent layout differs",
      block = block
    )
  }

  invisible(TRUE)
}

.bt_JAGS_bridge_latent_layout_metadata <- function(random_term){

  layout <- random_term$latent_layout
  if(is.null(layout)){
    return(NULL)
  }
  list(
    type = layout$type,
    structure = layout$structure,
    global_n_columns = layout$global_n_columns,
    n_groups = layout$n_groups,
    n_local = layout$n_local,
    group_columns = layout$group_columns,
    row_column = layout$row_column,
    row_local = layout$row_local,
    column_coordinates = layout$column_coordinates,
    node_names = layout$node_names
  )
}

.bt_JAGS_bridge_random_term_structure <- function(random_term){

  .bt_random_effect_structure(
    random_term,
    context = "Bridge sampling random-effect metadata"
  )
}

.bt_JAGS_bridge_scale_metadata <- function(random_term){

  list(
    sd_parameter_names = random_term$sd_parameter_names,
    sd_binding = .bt_JAGS_bridge_sd_binding_metadata(
      random_term$sd_binding,
      n_columns = random_term$n_columns
    ),
    homogeneous_sd = .bt_random_effect_homogeneous_sd_metadata(
      random_term,
      context = "Bridge sampling random-effect metadata"
    )
  )
}

.bt_JAGS_bridge_sd_binding_metadata <- function(binding, n_columns = NULL){

  if(is.null(binding)){
    return(NULL)
  }

  .bt_check_random_sd_binding(binding)
  if(isTRUE(binding$true_allocation) &&
     length(binding$allocations) > 0L &&
     identical(binding$allocations[[1L]]$target, "sd_component")){
    .bt_check_random_sd_component_binding(
      binding = binding,
      n_columns = n_columns,
      context = "Bridge sampling random-effect metadata"
    )
  }
  out <- binding[intersect(
    names(binding),
    c(
      "application", "factors", "factors_by_column",
      "true_allocation", "allocations", "sd_component_names",
      "sd_component_terms", "sd_component_index_by_column"
    )
  )]
  out$source <- .bt_JAGS_bridge_parameter_source_metadata(binding$source)
  out$sources_by_column <- lapply(
    binding$sources_by_column,
    .bt_JAGS_bridge_parameter_source_metadata
  )
  if(!is.null(out$factors)){
    out$factors <- lapply(out$factors, .bt_JAGS_bridge_allocation_factor_metadata)
  }
  if(!is.null(out$factors_by_column)){
    out$factors_by_column <- lapply(out$factors_by_column, function(factors){
      lapply(factors, .bt_JAGS_bridge_allocation_factor_metadata)
    })
  }
  if(!is.null(out$allocations)){
    out$allocations <- lapply(out$allocations, .bt_JAGS_bridge_allocation_metadata)
  }

  out
}

.bt_JAGS_bridge_allocation_metadata <- function(allocation){

  if(is.null(allocation)){
    return(NULL)
  }

  allocation <- allocation[intersect(
    names(allocation),
    c(
      "label", "terms", "index", "target", "scale",
      "parent", "source", "factors", "parent_factors", "n_targets", "scale_name",
      "weight_name", "leaf_names", "leaf_terms",
      "leaf_index_by_column", "source_node", "weight_suffix",
      "scale_suffix", "sd_component_names", "sd_component_terms",
      "sd_component_index_by_column", "inclusion"
    )
  )]
  if(!is.null(allocation$source)){
    allocation$source <- .bt_JAGS_bridge_parameter_source_metadata(allocation$source)
  }
  if(!is.null(allocation$factors)){
    allocation$factors <- lapply(allocation$factors, .bt_JAGS_bridge_allocation_factor_metadata)
  }
  if(!is.null(allocation$parent_factors)){
    allocation$parent_factors <- lapply(allocation$parent_factors, .bt_JAGS_bridge_allocation_factor_metadata)
  }
  if(!is.null(allocation$inclusion)){
    allocation$inclusion <- lapply(
      allocation$inclusion,
      .bt_JAGS_bridge_allocation_inclusion_metadata
    )
  }

  allocation
}

.bt_JAGS_bridge_allocation_factor_metadata <- function(factor){

  factor[intersect(
    names(factor),
    c("weight_name", "index", "scale", "n_targets", "inclusion_name")
  )]
}

.bt_JAGS_bridge_allocation_inclusion_metadata <- function(inclusion){

  if(is.null(inclusion)){
    return(NULL)
  }

  inclusion[intersect(
    names(inclusion),
    c("component", "index", "prob_suffix", "prob_name", "indicator_name")
  )]
}

.bt_JAGS_bridge_parameter_source_metadata <- function(source){

  if(is.null(source)){
    return(NULL)
  }
  if(!is.list(source)){
    return(source)
  }

  source <- source[intersect(
    names(source),
    c(
      "name", "shape", "kind", "owned",
      "scale_name", "scale_suffix", "source", "values"
    )
  )]
  if(!is.null(source$source)){
    source$source <- .bt_JAGS_bridge_parameter_source_metadata(source$source)
  }

  source
}

.bt_JAGS_bridge_structured_index_metadata <- function(random_term){

  index <- random_term$structured_index
  if(is.null(index)){
    return(NULL)
  }

  index[intersect(names(index), c("variables", "name", "label", "structure"))]
}

.bt_JAGS_bridge_car_metadata <- function(random_term){

  car <- random_term$car
  if(is.null(car)){
    return(NULL)
  }

  car[intersect(names(car), c("time_variable", "time_values", "distance_matrix"))]
}

.bt_JAGS_bridge_correlation_metadata <- function(random_term){

  structure <- .bt_JAGS_bridge_random_term_structure(random_term)
  correlation <- .bt_random_effect_correlation_metadata(
    random_term,
    structure = structure,
    context = "Bridge sampling random-effect metadata"
  )
  if(is.null(correlation)){
    return(NULL)
  }

  correlation[intersect(
    names(correlation),
    c(
      "type", "eta", "primitive_names", "primitive_bounds",
      "structure", "rho_name", "sample_name", "prior_name", "parameter_name",
      "sample_fixed", "rho_scale", "bounds", "cholesky_name", "correlation_name",
      "time_variable", "time_values", "distance_matrix"
    )
  )]
}

.bt_JAGS_bridge_group_covariance_metadata <- function(random_term){

  .bt_random_effect_group_covariance_metadata(random_term)
}

.bt_JAGS_bridge_random_design_mismatch <- function(parameter, detail,
                                                   block = NULL){

  stop(
    "JAGS_bridgesampling() rebuilt formula random-effect design does not match the fitted design for parameter '",
    parameter,
    "'",
    if(!is.null(block)) paste0(", block '", block, "'") else "",
    ": ",
    detail,
    ". Supply the same formula, data, scaling, and prior_random() metadata used to fit the model.",
    call. = FALSE
  )
}

.JAGS_formula_scale_list_from_fit <- function(fit, formula_parameters){

  formula_scale <- attr(fit, "formula_scale")
  if(is.null(formula_scale) || length(formula_scale) == 0L){
    return(NULL)
  }

  scale_list <- vector("list", length(formula_parameters))
  names(scale_list) <- formula_parameters

  for(parameter in intersect(formula_parameters, names(formula_scale))){
    parameter_scale <- formula_scale[[parameter]]
    if(is.null(parameter_scale) || length(parameter_scale) == 0L){
      next
    }

    scaled_terms <- names(parameter_scale)
    if(is.null(scaled_terms) || length(scaled_terms) == 0L){
      next
    }

    parameter_prefix <- paste0(parameter, "_")
    predictor_terms  <- ifelse(
      startsWith(scaled_terms, parameter_prefix),
      substring(scaled_terms, nchar(parameter_prefix) + 1L),
      scaled_terms
    )

    scale_list[[parameter]] <- as.list(stats::setNames(rep(TRUE, length(predictor_terms)), predictor_terms))
  }

  scale_list <- scale_list[!vapply(scale_list, is.null, logical(1))]
  if(length(scale_list) == 0L){
    return(NULL)
  }

  return(scale_list)
}
