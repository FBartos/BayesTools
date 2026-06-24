.bt_JAGS_bridge_check_row_indexed_external_sd_sources <- function(formula_design_list,
                                                                  bridgesampling_posterior){

  if(length(formula_design_list) == 0L){
    return(invisible(TRUE))
  }
  if(!is.matrix(bridgesampling_posterior)){
    stop("'bridgesampling_posterior' must be a matrix.", call. = FALSE)
  }

  bridge_parameters <- colnames(bridgesampling_posterior)
  bridge_lb <- attr(bridgesampling_posterior, "lb")
  bridge_ub <- attr(bridgesampling_posterior, "ub")

  for(parameter in names(formula_design_list)){
    design <- formula_design_list[[parameter]]
    if(!.bt_formula_design_has_any_random_effects(design)){
      next
    }
    for(random_term in .bt_formula_design_random_effects(design)){
      if(.bt_random_effect_has_row_indexed_external_sd(random_term)){
        source <- .bt_random_effect_row_indexed_source(random_term)
        source_names <- .bt_parameter_source_row_names(
          source$source,
          nrow(random_term$model_matrix)
        )
        if(.bt_parameter_source_has_values(source)){
          present <- intersect(source_names, bridge_parameters)
          if(length(present) > 0L){
            stop(
              "JAGS_bridgesampling() row-indexed external SD source '",
              .bt_random_effect_external_sd_source_label(random_term),
              "' for random-effect block '",
              random_term$block_name,
              "' provides a parameter_source() values function and also appears ",
              "as bridge parameter column(s). Use one row-source reconstruction ",
              "path only. Conflicting source column(s): ",
              paste0("'", present[seq_len(min(3L, length(present)))], "'",
                     collapse = ", "),
              if(length(present) > 3L) ", ..." else "",
              ".",
              call. = FALSE
            )
          }
          next
        }
        if(all(source_names %in% bridge_parameters)){
          .bt_JAGS_bridge_check_row_indexed_external_sd_bounds(
            random_term = random_term,
            source_names = source_names,
            bridge_lb = bridge_lb,
            bridge_ub = bridge_ub
          )
          next
        }
        missing <- setdiff(source_names, bridge_parameters)
        if(length(missing) > 0L){
          stop(
            "JAGS_bridgesampling() cannot reconstruct row-indexed external SD source '",
            .bt_random_effect_external_sd_source_label(random_term),
            "' for random-effect block '",
            random_term$block_name,
            "'. Provide a parameter_source(..., values = function(parameters, data, n_rows) ...) ",
            "or include the row-wise source column(s) as bridge parameters with appropriate priors or add_parameters/add_bounds. ",
            "Missing source column(s): ",
            paste0("'", missing[seq_len(min(3L, length(missing)))], "'", collapse = ", "),
            if(length(missing) > 3L) ", ..." else "",
            ".",
            call. = FALSE
          )
        }
      }
    }
  }

  invisible(TRUE)
}

.bt_JAGS_bridge_check_row_indexed_external_sd_bounds <- function(random_term,
                                                                 source_names,
                                                                 bridge_lb,
                                                                 bridge_ub){

  if(is.null(bridge_lb) || is.null(bridge_ub) ||
     !all(source_names %in% names(bridge_lb)) ||
     !all(source_names %in% names(bridge_ub))){
    stop(
      "JAGS_bridgesampling() row-indexed external SD source '",
      .bt_random_effect_external_sd_source_label(random_term),
      "' for random-effect block '",
      random_term$block_name,
      "' is missing lower or upper bridge bounds.",
      call. = FALSE
    )
  }

  lower <- bridge_lb[source_names]
  upper <- bridge_ub[source_names]
  bad_bounds <- is.na(lower) | is.na(upper) | lower < 0 | lower >= upper
  if(any(bad_bounds)){
    bad_names <- source_names[bad_bounds]
    stop(
      "JAGS_bridgesampling() row-indexed external SD source '",
      .bt_random_effect_external_sd_source_label(random_term),
      "' for random-effect block '",
      random_term$block_name,
      "' is included as bridge parameter(s), but its bridge bounds must ",
      "be valid standard-deviation bounds with non-negative lower bounds. ",
      "Invalid source column(s): ",
      paste0("'", bad_names[seq_len(min(3L, length(bad_names)))], "'",
             collapse = ", "),
      if(length(bad_names) > 3L) ", ..." else "",
      ".",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

.bt_JAGS_bridge_formula_allocation_inclusion_names <- function(formula_design_list){

  formula_design_list <- .bt_JAGS_bridge_formula_design_list(formula_design_list)
  .bt_random_variance_allocation_inclusion_indicator_names(formula_design_list)
}

.bt_JAGS_bridge_check_no_allocation_inclusion <- function(formula_design_list){

  indicator_names <- .bt_JAGS_bridge_formula_allocation_inclusion_names(
    formula_design_list
  )
  if(length(indicator_names) == 0L){
    return(invisible(TRUE))
  }

  stop(
    "Bridge sampling for variance allocation inclusion gates is not implemented ",
    "because the gates introduce discrete Bernoulli indicators. ",
    "Affected indicator node(s): ",
    paste0("'", indicator_names[seq_len(min(3L, length(indicator_names)))], "'",
           collapse = ", "),
    if(length(indicator_names) > 3L) ", ..." else "",
    ".",
    call. = FALSE
  )
}

.bt_JAGS_formula_random_bridge_parameters <- function(formula_design_list){

  if(length(formula_design_list) == 0L){
    return(list(
      parameters = character(),
      bounds = list(lb = numeric(), ub = numeric())
    ))
  }

  parameters <- character()
  lb <- numeric()
  ub <- numeric()

  for(parameter in names(formula_design_list)){
    design <- formula_design_list[[parameter]]
    if(!.bt_formula_design_has_any_random_effects(design)){
      next
    }
    .bt_JAGS_bridge_validate_formula_random_compile(
      parameter = parameter,
      design = design,
      label = "stored"
    )
    if(!identical(design$random_effects_interface, "prior_random")){
      stop(
        "JAGS_bridgesampling() supports formula random effects only through the 'prior_random()' interface for parameter '",
        parameter,
        "'.",
        call. = FALSE
      )
    }
    for(random_term in .bt_formula_design_sampled_random_effects(design)){
      n_groups <- random_term$n_groups
      n_columns <- random_term$n_columns
      z_names <- as.vector(.bt_random_effect_latent_names(
        random_term = random_term,
        n_groups = n_groups,
        n_columns = n_columns
      ))
      parameters <- c(parameters, z_names)
      lb <- c(lb, stats::setNames(rep(-Inf, length(z_names)), z_names))
      ub <- c(ub, stats::setNames(rep( Inf, length(z_names)), z_names))

    }
    for(random_term in .bt_formula_design_random_effects(design)){
      n_columns <- random_term$n_columns
      if(identical(.bt_JAGS_bridge_random_term_structure(random_term), "us") &&
         n_columns > 1L){
        u_names <- .bt_random_effect_lkj_primitive_names(
          random_term,
          n_columns,
          context = "Bridge sampling random-effect metadata"
        )
        parameters <- c(parameters, u_names)
        lb <- c(lb, stats::setNames(rep(0, length(u_names)), u_names))
        ub <- c(ub, stats::setNames(rep(1, length(u_names)), u_names))
      }
    }
  }

  keep <- !duplicated(parameters)
  parameters <- parameters[keep]
  lb <- lb[parameters]
  ub <- ub[parameters]

  list(
    parameters = parameters,
    bounds = list(lb = lb, ub = ub)
  )
}

.bt_JAGS_bridge_merge_add_parameters <- function(add_parameters, add_bounds,
                                                 bridge_parameters,
                                                 bridge_bounds){

  if((is.null(add_parameters) || length(add_parameters) == 0L) &&
     !is.null(add_bounds)){
    stop("'add_bounds' requires at least one 'add_parameters' entry.", call. = FALSE)
  }

  if(length(bridge_parameters) == 0L){
    return(list(add_parameters = add_parameters, add_bounds = add_bounds))
  }

  user_parameters <- if(is.null(add_parameters)) character() else add_parameters
  if(length(user_parameters) > 0L){
    user_bounds <- .bt_JAGS_bridge_validate_add_bounds(user_parameters, add_bounds)
    user_lb <- user_bounds$lb
    user_ub <- user_bounds$ub
  }else{
    user_lb <- numeric()
    user_ub <- numeric()
  }

  bridge_lb <- bridge_bounds$lb
  bridge_ub <- bridge_bounds$ub
  overlapping_bridge <- intersect(user_parameters, bridge_parameters)
  if(length(overlapping_bridge) > 0L){
    conflicting_bridge <- overlapping_bridge[
      vapply(overlapping_bridge, function(parameter){
        !isTRUE(all.equal(unname(user_lb[[parameter]]), unname(bridge_lb[[parameter]]))) ||
          !isTRUE(all.equal(unname(user_ub[[parameter]]), unname(bridge_ub[[parameter]])))
      }, logical(1))
    ]
    if(length(conflicting_bridge) > 0L){
      stop(
        "User-supplied bounds conflict with automatically inferred formula random-effect bridge bounds for parameter(s): ",
        paste(conflicting_bridge, collapse = ", "),
        ".",
        call. = FALSE
      )
    }
  }
  missing_bridge <- setdiff(bridge_parameters, names(user_lb))
  combined_parameters <- unique(c(user_parameters, bridge_parameters))
  combined_lb <- c(user_lb, bridge_lb[missing_bridge])
  combined_ub <- c(user_ub, bridge_ub[missing_bridge])
  combined_lb <- combined_lb[combined_parameters]
  combined_ub <- combined_ub[combined_parameters]

  list(
    add_parameters = combined_parameters,
    add_bounds = list(lb = combined_lb, ub = combined_ub)
  )
}

.bt_JAGS_bridge_validate_add_bounds <- function(add_parameters, add_bounds){

  if(!is.character(add_parameters)){
    stop("'add_parameters' must be a character vector.", call. = FALSE)
  }
  if(anyDuplicated(add_parameters)){
    stop("'add_parameters' must be unique.", call. = FALSE)
  }
  if(!is.list(add_bounds)){
    stop("'add_bounds' must be a list.", call. = FALSE)
  }
  if(length(add_bounds) != 2L || !all(c("lb", "ub") %in% names(add_bounds)) ||
     !all(names(add_bounds) %in% c("lb", "ub"))){
    stop("'add_bounds' must contain lower and upper bounds ('lb' and 'ub').", call. = FALSE)
  }

  lb <- add_bounds[["lb"]]
  ub <- add_bounds[["ub"]]
  if(length(lb) != length(add_parameters) || length(ub) != length(add_parameters)){
    stop("'lb' and 'ub' must have the same length as 'add_parameters'.", call. = FALSE)
  }
  if(!is.numeric(lb) || !is.numeric(ub)){
    stop("'lb' and 'ub' must be numeric vectors.", call. = FALSE)
  }

  lb <- .bt_JAGS_bridge_normalize_bound_names(lb, add_parameters, "lb")
  ub <- .bt_JAGS_bridge_normalize_bound_names(ub, add_parameters, "ub")
  if(any(is.na(lb)) || any(is.na(ub))){
    stop("'add_bounds' must not contain NA values.", call. = FALSE)
  }
  if(any(lb >= ub)){
    stop("Lower bounds in 'add_bounds' must be smaller than upper bounds.", call. = FALSE)
  }

  list(lb = lb, ub = ub)
}

.bt_JAGS_bridge_normalize_bound_names <- function(bounds, add_parameters,
                                                  bound_name){

  bound_names <- names(bounds)
  if(is.null(bound_names)){
    stop(
      "'add_bounds$", bound_name,
      "' names must be unique and match 'add_parameters'.",
      call. = FALSE
    )
  }

  if(any(!nzchar(bound_names)) || anyDuplicated(bound_names)){
    stop(
      "'add_bounds$", bound_name,
      "' names must be unique and match 'add_parameters'.",
      call. = FALSE
    )
  }

  missing_names <- setdiff(add_parameters, bound_names)
  unknown_names <- setdiff(bound_names, add_parameters)
  if(length(missing_names) > 0L || length(unknown_names) > 0L){
    stop(
      "'add_bounds$", bound_name,
      "' names must match 'add_parameters'.",
      call. = FALSE
    )
  }

  bounds[add_parameters]
}

.bt_JAGS_bridge_check_random_posterior <- function(posterior, bridge_parameters){

  if(length(bridge_parameters) == 0L){
    return(invisible(TRUE))
  }
  missing <- setdiff(bridge_parameters, colnames(posterior))
  if(length(missing) == 0L){
    return(invisible(TRUE))
  }

  stop(
    "Bridge sampling for formula random effects requires posterior samples of standardized latent random effects and LKJ primitive coordinates. ",
    "Refit with 'prior_random()' and 'random_monitor(latent = TRUE)' for the affected blocks. Missing parameter(s): ",
    paste(utils::head(missing, 8L), collapse = ", "),
    if(length(missing) > 8L) ", ..." else "",
    ".",
    call. = FALSE
  )
}

.bt_JAGS_bridge_apply_random_scalar_rho_bounds <- function(bridgesampling_posterior,
                                                           formula_design_list){

  rho_bridge <- .bt_JAGS_formula_random_scalar_rho_bridge_parameters(
    formula_design_list
  )
  if(length(rho_bridge$parameters) == 0L){
    return(bridgesampling_posterior)
  }

  missing <- setdiff(rho_bridge$parameters, colnames(bridgesampling_posterior))
  if(length(missing) > 0L){
    stop(
      "Bridge sampling for formula random effects requires posterior samples of scalar random-effect correlation coordinates. ",
      "Missing parameter(s): ",
      paste(utils::head(missing, 8L), collapse = ", "),
      if(length(missing) > 8L) ", ..." else "",
      ".",
      call. = FALSE
    )
  }

  lb <- attr(bridgesampling_posterior, "lb")
  ub <- attr(bridgesampling_posterior, "ub")
  for(parameter in rho_bridge$parameters){
    if(!parameter %in% names(lb) || !parameter %in% names(ub)){
      stop(
        "Bridge sampling scalar random-effect correlation coordinate '",
        parameter,
        "' is missing lower or upper bridge bounds.",
        call. = FALSE
      )
    }
    next_lb <- max(lb[[parameter]], rho_bridge$bounds$lb[[parameter]])
    next_ub <- min(ub[[parameter]], rho_bridge$bounds$ub[[parameter]])
    if(is.na(next_lb) || is.na(next_ub) || next_lb >= next_ub){
      stop(
        "Bridge sampling scalar random-effect correlation bounds conflict for parameter '",
        parameter,
        "'.",
        call. = FALSE
      )
    }
    lb[[parameter]] <- next_lb
    ub[[parameter]] <- next_ub
  }

  attr(bridgesampling_posterior, "lb") <- lb
  attr(bridgesampling_posterior, "ub") <- ub

  bridgesampling_posterior
}

.bt_JAGS_formula_random_scalar_rho_bridge_parameters <- function(formula_design_list){

  if(length(formula_design_list) == 0L){
    return(list(
      parameters = character(),
      bounds = list(lb = numeric(), ub = numeric())
    ))
  }

  parameters <- character()
  lb <- numeric()
  ub <- numeric()

  design_names <- names(formula_design_list)
  if(is.null(design_names)){
    design_names <- rep("", length(formula_design_list))
  }
  for(design_i in seq_along(formula_design_list)){
    design <- formula_design_list[[design_i]]
    if(!.bt_formula_design_has_any_random_effects(design)){
      next
    }
    .bt_JAGS_bridge_validate_formula_random_compile(
      parameter = .bt_JAGS_bridge_design_parameter_name(
        design,
        fallback = design_names[[design_i]]
      ),
      design = design,
      label = "stored"
    )
    for(random_term in .bt_formula_design_random_effects(design)){
      rho_parameter <- .bt_JAGS_bridge_scalar_rho_parameter(random_term)
      if(is.null(rho_parameter)){
        next
      }
      parameters <- c(parameters, rho_parameter$parameter)
      lb <- c(lb, stats::setNames(rho_parameter$lower, rho_parameter$parameter))
      ub <- c(ub, stats::setNames(rho_parameter$upper, rho_parameter$parameter))
    }
  }

  if(length(parameters) == 0L){
    return(list(
      parameters = character(),
      bounds = list(lb = numeric(), ub = numeric())
    ))
  }

  unique_parameters <- unique(parameters)
  out_lb <- vapply(unique_parameters, function(parameter){
    max(lb[parameters == parameter])
  }, numeric(1))
  out_ub <- vapply(unique_parameters, function(parameter){
    min(ub[parameters == parameter])
  }, numeric(1))
  names(out_lb) <- unique_parameters
  names(out_ub) <- unique_parameters

  if(any(out_lb >= out_ub)){
    stop(
      "Bridge sampling scalar random-effect correlation bounds are inconsistent.",
      call. = FALSE
    )
  }

  list(
    parameters = unique_parameters,
    bounds = list(lb = out_lb, ub = out_ub)
  )
}

.bt_JAGS_bridge_scalar_rho_parameter <- function(random_term){

  structure <- .bt_JAGS_bridge_random_term_structure(random_term)
  if(!structure %in% c("cs", "hcs", "ar1", "car", "har") ||
     random_term$n_columns <= 1L){
    return(NULL)
  }

  correlation <- .bt_random_effect_correlation_metadata(
    random_term,
    structure = structure,
    context = "Bridge sampling random-effect metadata"
  )
  if(is.null(correlation) || !identical(correlation$type, "rho")){
    stop(
      "Bridge sampling random-effect metadata",
      .bt_random_effect_metadata_block_detail(random_term),
      " are missing canonical scalar 'random_term$correlation'.",
      call. = FALSE
    )
  }
  if(!is.null(.bt_random_effect_rho_fixed_sample_metadata(
    correlation,
    random_term
  ))){
    return(NULL)
  }

  parameter <- .bt_JAGS_bridge_scalar_rho_sample_name(correlation, random_term)
  bounds <- .bt_JAGS_bridge_scalar_rho_sample_bounds(correlation, random_term)

  list(
    parameter = parameter,
    lower = bounds[["lower"]],
    upper = bounds[["upper"]]
  )
}

.bt_JAGS_bridge_scalar_rho_sample_name <- function(correlation, random_term){

  rho_scale <- .bt_random_effect_rho_scale_metadata(correlation, random_term)
  parameter <- if(identical(rho_scale, "rho")){
    correlation$rho_name
  }else{
    correlation$sample_name
  }

  if(!is.character(parameter) || length(parameter) != 1L ||
     is.na(parameter) || !nzchar(parameter)){
    stop(
      "Bridge sampling random-effect metadata",
      .bt_random_effect_metadata_block_detail(random_term),
      " are missing canonical scalar correlation sample name.",
      call. = FALSE
    )
  }

  parameter
}

.bt_JAGS_bridge_scalar_rho_sample_bounds <- function(correlation,
                                                     random_term){

  .bt_random_effect_rho_sample_bounds(
    correlation = correlation,
    random_term = random_term,
    context = "Bridge sampling random-effect metadata"
  )
}

