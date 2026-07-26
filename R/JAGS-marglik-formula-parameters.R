# .JAGS_marglik_parameters.spike_and_slab <- function(samples, prior, parameter_name){
#
#   .check_prior(prior)
#   if(!is.prior.spike_and_slab(prior))
#     stop("improper prior provided")
#   check_char(parameter_name, "parameter_name")
#
#   parameter <- list()
#   parameter[paste0(parameter_name, "_variable")]  <- .JAGS_marglik_parameters.simple(samples, prior[["variable"]],  paste0(parameter_name, "_variable"))
#   if(!is.prior.point(prior[[parameter_name]][["inclusion"]])){
#     parameter[paste0(parameter_name, "_inclusion")] <- .JAGS_marglik_parameters.simple(samples, prior[["inclusion"]], paste0(parameter_name, "_inclusion"))
#   }
#
#   return(parameter)
# }

#' @rdname JAGS_marglik_parameters
JAGS_marglik_parameters_formula      <- function(samples, formula_list, formula_data_list, formula_prior_list, prior_list_parameters,
                                                 formula_design_list = NULL,
                                                 model_data = NULL){

  # return empty list in case that no prior was specified
  if(length(formula_prior_list) == 0){
    return(list())
  }

  parameters <- list()
  random_parameters <- character()

  for(parameter in names(formula_prior_list)){
    # check for log(intercept) attribute on the formula
    formula_parameter <- if(!is.null(formula_list)) formula_list[[parameter]] else NULL
    if(!is.null(formula_parameter)){
      .bt_validate_formula_replay_grammar(formula_parameter)
    }
    log_intercept <- if(!is.null(formula_parameter)) isTRUE(attr(formula_parameter, "log(intercept)")) else FALSE
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
      parameter_prior_list <- .bt_JAGS_marglik_formula_fixed_priors(parameter_prior_list, parameter)
    }
    if(.bt_formula_design_has_sampled_random_effects(design)){
      random_parameters <- c(random_parameters, parameter)
    }
    if(.bt_JAGS_formula_design_can_reconstruct(design)){
      parameters[[parameter]] <- .bt_JAGS_marglik_parameters_formula_design(
        samples = samples,
        design = design,
        formula_prior_list = parameter_prior_list,
        prior_list_parameters = prior_list_parameters,
        log_intercept = log_intercept
      )
    }else{
      parameters[[parameter]] <- .JAGS_marglik_parameters_formula_get(samples, parameter, formula_data_list[[parameter]], parameter_prior_list, prior_list_parameters, log_intercept)
    }
  }

  for(parameter in random_parameters){
    design <- if(!is.null(formula_design_list)) formula_design_list[[parameter]] else NULL
    if(.bt_formula_design_has_sampled_random_effects(design)){
      source_parameters <- .bt_JAGS_marglik_parameter_source_parameters(
        samples = samples,
        prior_list_parameters = prior_list_parameters,
        formula_parameters = parameters
      )
      source_parameters <- .bt_parameter_source_forbid_formula_parameters(
        source_parameters,
        unique(random_parameters)
      )
      source_formula_data <- if(!is.null(formula_data_list)){
        formula_data_list[[parameter]]
      }else{
        NULL
      }
      source_data <- .bt_JAGS_marglik_parameter_source_data(
        model_data = model_data,
        formula_data = source_formula_data,
        design = design
      )
      parameters[[parameter]] <- parameters[[parameter]] +
        .bt_JAGS_marglik_random_effects_value(
          samples = samples,
          design = design,
          formula_prior_list = formula_prior_list[[parameter]],
          data = source_data,
          parameters = source_parameters
        )
    }
  }

  return(parameters)
}

.bt_JAGS_marglik_parameter_source_parameters <- function(samples,
                                                         prior_list_parameters,
                                                         formula_parameters){

  out <- as.list(samples)
  if(length(prior_list_parameters) > 0L){
    out[names(prior_list_parameters)] <- prior_list_parameters
  }
  if(length(formula_parameters) > 0L){
    out[names(formula_parameters)] <- formula_parameters
  }

  out
}

.bt_JAGS_marglik_parameter_source_data <- function(model_data,
                                                   formula_data,
                                                   design){

  out <- list()
  out <- .bt_JAGS_marglik_merge_source_data(out, formula_data)
  out <- .bt_JAGS_marglik_merge_source_data(out, model_data)
  if(inherits(design, "BayesTools_formula_design") &&
     !is.null(design$source_data)){
    out <- .bt_JAGS_marglik_fill_source_data(out, design$source_data)
  }

  out
}

.bt_JAGS_marglik_fill_source_data <- function(out, data){

  if(is.null(data)){
    return(out)
  }
  data_list <- if(is.data.frame(data)){
    as.list(data)
  }else if(is.list(data)){
    data
  }else{
    return(out)
  }
  data_names <- names(data_list)
  if(is.null(data_names)){
    return(out)
  }
  keep <- !is.na(data_names) & nzchar(data_names)
  if(!any(keep)){
    return(out)
  }
  data_list <- data_list[keep]
  missing <- setdiff(names(data_list), names(out))
  out[missing] <- data_list[missing]

  out
}

.bt_JAGS_marglik_merge_source_data <- function(out, data){

  if(is.null(data)){
    return(out)
  }
  data_list <- if(is.data.frame(data)){
    as.list(data)
  }else if(is.list(data)){
    data
  }else{
    return(out)
  }
  data_names <- names(data_list)
  if(is.null(data_names)){
    return(out)
  }
  keep <- !is.na(data_names) & nzchar(data_names)
  if(!any(keep)){
    return(out)
  }
  data_list <- data_list[keep]
  overlap <- intersect(names(out), names(data_list))
  for(name in overlap){
    if(!.bt_JAGS_marglik_source_data_equal(out[[name]], data_list[[name]])){
      stop(
        "JAGS_bridgesampling() row-indexed source reconstruction received ",
        "conflicting data for variable '", name, "'. Supply formula data, fitted ",
        "formula metadata, and model data with matching row-aligned values, or ",
        "remove the duplicate variable from the non-formula data.",
        call. = FALSE
      )
    }
  }
  out[names(data_list)] <- data_list

  out
}

.bt_JAGS_marglik_source_data_equal <- function(x, y){

  identical(x, y)
}

.bt_JAGS_formula_design_can_reconstruct <- function(design){

  inherits(design, "BayesTools_formula_design") &&
    !is.null(design$parameter) &&
    !is.null(design$model_matrix) &&
    !is.null(design$assign) &&
    !is.null(design$model_terms)
}

.bt_JAGS_marglik_parameters_formula_design <- function(samples, design,
                                                       formula_prior_list,
                                                       prior_list_parameters,
                                                       log_intercept = FALSE){

  parameter <- design$parameter
  output <- rep(0, nrow(design$model_matrix))
  if(length(formula_prior_list) == 0L){
    return(output)
  }

  intercept_name <- paste0(parameter, "_intercept")
  if(intercept_name %in% names(formula_prior_list)){
    intercept_prior <- formula_prior_list[[intercept_name]]
    .bt_validate_formula_reconstruction_prior(
      intercept_prior,
      intercept_name
    )
    intercept_value <- .JAGS_marglik_parameter_values(
      samples,
      intercept_prior,
      intercept_name
    )
    if(isTRUE(log_intercept) || isTRUE(design$log_intercept)){
      intercept_value <- log(intercept_value)
    }
    output <- output + .bt_JAGS_marglik_prior_multiply_by(
      intercept_prior,
      prior_list_parameters
    ) * intercept_value
  }

  remaining_terms <- setdiff(names(formula_prior_list), intercept_name)
  for(term in remaining_terms){
    term_prior <- formula_prior_list[[term]]
    .bt_validate_formula_reconstruction_prior(term_prior, term)
    model_term <- sub(paste0("^", parameter, "_"), "", term)
    columns <- .bt_JAGS_formula_design_term_columns(design, model_term)
    term_data <- design$model_matrix[, columns, drop = FALSE]
    multiply_by <- .bt_JAGS_marglik_prior_multiply_by(
      term_prior,
      prior_list_parameters
    )

    if(is.prior.point(term_prior) && !is.prior.factor(term_prior)){
      output <- output + multiply_by * term_prior[["parameters"]][["location"]] * as.vector(term_data)

    }else if(is.prior.point(term_prior) && is.prior.factor(term_prior)){
      if(.get_prior_factor_levels(term_prior) == 1){
        output <- output + multiply_by * term_prior[["parameters"]][["location"]] * as.vector(term_data)
      }else{
        output <- output + multiply_by * as.vector(term_data %*% rep(term_prior[["parameters"]][["location"]], .get_prior_factor_levels(term_prior)))
      }

    }else if(is.prior.factor(term_prior)){
      if(.get_prior_factor_levels(term_prior) == 1){
        term_value <- .JAGS_marglik_parameter_values(samples, term_prior, term)
        output <- output + multiply_by * term_value * as.vector(term_data)
      }else{
        term_names <- paste0(term, "[", 1:.get_prior_factor_levels(term_prior), "]")
        term_values <- .JAGS_marglik_parameter_values(samples, term_prior, term_names)
        output <- output + multiply_by * as.vector(term_data %*% term_values)
      }

    }else if(is.prior.simple(term_prior)){
      term_value <- .JAGS_marglik_parameter_values(samples, term_prior, term)
      output <- output + multiply_by * term_value * as.vector(term_data)
    }else{
      stop(
        "Internal formula reconstruction prior dispatch failed for '",
        term, "'.",
        call. = FALSE
      )
    }
  }

  as.vector(output)
}

.bt_JAGS_formula_design_term_columns <- function(design, model_term){

  term_index <- match(model_term, design$model_terms)
  if(is.na(term_index)){
    stop(
      "Stored formula design for parameter '",
      design$parameter,
      "' is missing model term '",
      model_term,
      "'.",
      call. = FALSE
    )
  }

  columns <- which(design$assign == (term_index - 1L))
  if(length(columns) == 0L){
    stop(
      "Stored formula design for parameter '",
      design$parameter,
      "' is missing model-matrix columns for term '",
      model_term,
      "'.",
      call. = FALSE
    )
  }

  columns
}

.bt_JAGS_marglik_prior_multiply_by <- function(prior, prior_list_parameters){

  multiply_by <- attr(prior, "multiply_by")
  if(is.null(multiply_by)){
    return(1)
  }
  if(is.numeric(multiply_by)){
    return(multiply_by)
  }

  .bt_JAGS_marglik_resolve_named_multiply_by(
    multiply_by = multiply_by,
    prior_list_parameters = prior_list_parameters
  )
}

.bt_JAGS_marglik_resolve_named_multiply_by <- function(multiply_by,
                                                        prior_list_parameters){

  value <- prior_list_parameters[[multiply_by]]
  if(is.null(value)){
    stop(
      "Formula prior 'multiply_by' parameter '",
      multiply_by,
      "' is missing from 'prior_list_parameters'.",
      call. = FALSE
    )
  }

  value
}

.bt_JAGS_marglik_formula_fixed_priors <- function(formula_prior_list, parameter){

  if(length(formula_prior_list) == 0L){
    return(formula_prior_list)
  }

  random_prefix <- paste0(parameter, "__xREx__")
  is_random <- startsWith(names(formula_prior_list), random_prefix) |
    vapply(formula_prior_list, function(prior){
      .bt_is_random_effect_prior(prior, include_summary = FALSE)
    }, logical(1))

  formula_prior_list[!is_random]
}

.bt_JAGS_marglik_random_effects_value <- function(samples, design,
                                                  formula_prior_list,
                                                  data = NULL,
                                                  parameters = NULL){

  output <- rep(0, nrow(design$model_matrix))
  for(random_term in .bt_formula_design_sampled_random_effects(design)){
    output <- output + .bt_JAGS_marglik_random_effect_value(
      samples = samples,
      random_term = random_term,
      prior_list = formula_prior_list,
      data = data,
      parameters = parameters
    )
  }

  output
}

.bt_JAGS_marglik_random_effect_value <- function(samples, random_term,
                                                 prior_list,
                                                 data = NULL,
                                                 parameters = NULL){

  .bt_JAGS_marglik_check_random_effect_dirichlet_samples(
    samples = samples,
    random_term = random_term,
    prior_list = prior_list
  )

  if(.bt_random_effect_has_row_indexed_external_sd(random_term)){
    return(.bt_JAGS_marglik_random_effect_row_indexed_value(
      samples = samples,
      random_term = random_term,
      prior_list = prior_list,
      data = data,
      parameters = parameters
    ))
  }

  if(inherits(
    random_term$latent_layout,
    "BayesTools_random_effect_structured_local_layout"
  )){
    return(.bt_JAGS_marglik_random_effect_structured_local_value(
      samples = samples,
      random_term = random_term,
      prior_list = prior_list
    ))
  }

  model_matrix <- random_term$model_matrix
  group_map <- random_term$group_map
  n_groups <- random_term$n_groups
  n_columns <- random_term$n_columns
  structure <- .bt_random_effect_structure(
    random_term,
    context = "Bridge-sampling random-effect metadata"
  )
  if(n_columns > 1L &&
     structure %in% c("cs", "hcs", "ar1", "car", "har")){
    posterior <- .bt_JAGS_marglik_random_effect_posterior_row(samples)
    sd_values <- .bt_JAGS_marglik_random_effect_sd_values(
      samples = samples,
      random_term = random_term,
      prior_list = prior_list
    )
    contribution <- .bt_random_effect_structured_contribution_from_latent(
      random_term = random_term,
      model_matrix = model_matrix,
      group_map = group_map,
      posterior = posterior,
      scale_draws = matrix(sd_values, nrow = 1L)
    )
    if(!is.null(contribution)){
      return(as.vector(contribution[, 1L]))
    }
  }

  z_names <- .bt_random_effect_latent_names(
    random_term = random_term,
    n_groups = n_groups,
    n_columns = n_columns
  )
  if(!all(as.vector(z_names) %in% names(samples))){
    stop(
      "Bridge samples are missing standardized latent random effects for block '",
      random_term$block_name,
      "'.",
      call. = FALSE
    )
  }
  z <- matrix(
    unname(samples[as.vector(z_names)]),
    nrow = n_groups,
    ncol = n_columns
  )
  z_draws <- lapply(seq_len(n_columns), function(column){
    matrix(z[, column], nrow = 1L)
  })

  sd_values <- .bt_JAGS_marglik_random_effect_sd_values(
    samples = samples,
    random_term = random_term,
    prior_list = prior_list
  )
  if(structure %in% c("diag", "id")){
    contribution <- .bt_random_effect_independent_contribution_from_latent(
      random_term = random_term,
      model_matrix = model_matrix,
      group_map = group_map,
      posterior = .bt_JAGS_marglik_random_effect_posterior_row(samples),
      column_scale_draws = matrix(sd_values, nrow = 1L)
    )
    if(is.null(contribution)){
      stop(
        "Bridge samples are missing standardized latent random effects for block '",
        random_term$block_name,
        "'.",
        call. = FALSE
      )
    }
    return(as.vector(contribution[, 1L]))
  }
  L <- .bt_JAGS_marglik_random_effect_cholesky(
    samples = samples,
    random_term = random_term
  )

  contribution <- .bt_random_effect_contribution_from_latent_draws(
    model_matrix = model_matrix,
    group_map = group_map,
    z_draws = z_draws,
    sd_draws = matrix(sd_values, nrow = 1L),
    cholesky = array(L, dim = c(1L, n_columns, n_columns))
  )
  as.vector(contribution[, 1L])
}

.bt_JAGS_marglik_random_effect_structured_local_value <- function(
    samples, random_term, prior_list){

  layout <- random_term$latent_layout
  if(!all(layout$node_names %in% names(samples))){
    stop(
      "Bridge samples are missing group-local standardized latent random effects for block '",
      random_term$block_name, "'.",
      call. = FALSE
    )
  }
  sd_values <- .bt_JAGS_marglik_random_effect_sd_values(
    samples = samples,
    random_term = random_term,
    prior_list = prior_list
  )
  rho <- .bt_JAGS_marglik_random_effect_rho(samples, random_term)
  coefficients <- vector("list", layout$n_groups)
  for(group in seq_len(layout$n_groups)){
    columns <- layout$group_columns[[group]]
    names <- .bt_random_effect_structured_local_node_names(
      parameter_stem = random_term$parameter_stem,
      group = rep(group, length(columns)),
      column = columns
    )
    unit <- .bt_random_effect_structured_subset_transform(
      structure = layout$structure,
      columns = columns,
      latent = unname(samples[names]),
      rho = rho,
      global_n_columns = layout$global_n_columns,
      column_coordinates = layout$column_coordinates
    )
    coefficients[[group]] <- unit * sd_values[columns]
  }

  vapply(seq_along(layout$row_column), function(row){
    coefficients[[random_term$group_map[row]]][layout$row_local[row]]
  }, numeric(1))
}

.bt_JAGS_marglik_check_random_effect_dirichlet_samples <- function(
    samples, random_term, prior_list){

  parameter_names <- .bt_JAGS_marglik_random_effect_allocation_parameters(
    random_term
  )
  for(parameter_name in parameter_names){
    if(!parameter_name %in% names(prior_list)){
      next
    }
    prior <- prior_list[[parameter_name]]
    if(!is.prior.simplex(prior) || !identical(prior$distribution, "dirichlet")){
      next
    }
    K <- prior$parameters[["K"]]
    weight_names <- paste0(parameter_name, "[", seq_len(K), "]")
    eta_names <- paste0(.JAGS_prior_dirichlet_eta_name(parameter_name), "[", seq_len(K), "]")
    if(all(weight_names %in% names(samples)) && all(eta_names %in% names(samples))){
      stop(
        "Bridge samples contain both normalized Dirichlet allocation coordinates ",
        "and auxiliary eta coordinates for parameter '",
        parameter_name,
        "'. Use the auxiliary eta bridge coordinates only.",
        call. = FALSE
      )
    }
  }

  invisible(TRUE)
}

.bt_JAGS_marglik_random_effect_allocation_parameters <- function(random_term){

  binding <- random_term$sd_binding
  if(is.null(binding)){
    return(character())
  }
  .bt_check_random_sd_binding(binding)

  out <- character()
  collect_factors <- function(factors){
    if(!is.list(factors) || length(factors) == 0L){
      return(character())
    }
    vapply(factors, function(factor){
      if(is.list(factor) &&
         is.character(factor$weight_name) &&
         length(factor$weight_name) == 1L &&
         !is.na(factor$weight_name)){
        return(factor$weight_name)
      }
      NA_character_
    }, character(1))
  }
  out <- c(out, collect_factors(binding$factors))
  if(is.list(binding$factors_by_column)){
    for(factors in binding$factors_by_column){
      out <- c(out, collect_factors(factors))
    }
  }
  if(is.list(binding$allocations)){
    for(allocation in binding$allocations){
      if(!is.list(allocation)){
        next
      }
      if(is.character(allocation$weight_name) &&
         length(allocation$weight_name) == 1L &&
         !is.na(allocation$weight_name)){
        out <- c(out, allocation$weight_name)
      }
      out <- c(out, collect_factors(allocation$factors))
      out <- c(out, collect_factors(allocation$parent_factors))
    }
  }

  unique(out[!is.na(out) & nzchar(out)])
}

.bt_JAGS_marglik_random_effect_row_indexed_value <- function(samples,
                                                             random_term,
                                                             prior_list,
                                                             data = NULL,
                                                             parameters = NULL){

  model_matrix <- random_term$model_matrix
  group_map <- random_term$group_map
  n_groups <- random_term$n_groups
  n_columns <- random_term$n_columns

  posterior <- .bt_JAGS_marglik_random_effect_posterior_row(samples)
  structure <- .bt_random_effect_structure(
    random_term,
    context = "Bridge-sampling random-effect metadata"
  )
  unit_contribution <- NULL
  unit_columns <- NULL
  independent <- structure %in% c("diag", "id")
  if(n_columns > 1L &&
     structure %in% c("cs", "hcs", "ar1", "car", "har")){
    unit_contribution <- .bt_random_effect_structured_contribution_from_latent(
      random_term = random_term,
      model_matrix = model_matrix,
      group_map = group_map,
      posterior = posterior,
      scale_draws = matrix(1, nrow = 1L, ncol = n_columns)
    )
  }
  if(is.null(unit_contribution) && !isTRUE(independent)){
    z_names <- .bt_random_effect_latent_names(
      random_term = random_term,
      n_groups = n_groups,
      n_columns = n_columns
    )
    if(!all(as.vector(z_names) %in% names(samples))){
      stop(
        "Bridge samples are missing standardized latent random effects for block '",
        random_term$block_name,
        "'.",
        call. = FALSE
      )
    }
    z <- matrix(
      unname(samples[as.vector(z_names)]),
      nrow = n_groups,
      ncol = n_columns
    )
    z_draws <- lapply(seq_len(n_columns), function(column){
      matrix(z[, column], nrow = 1L)
    })

    L <- .bt_JAGS_marglik_random_effect_cholesky(
      samples = samples,
      random_term = random_term
    )
    unit_columns <- .bt_random_effect_column_contributions_from_latent_draws(
      model_matrix = model_matrix,
      group_map = group_map,
      z_draws = z_draws,
      sd_draws = matrix(1, nrow = 1L, ncol = n_columns),
      cholesky = array(L, dim = c(1L, n_columns, n_columns))
    )
  }

  source_draws <- .bt_JAGS_marglik_row_indexed_external_sd_source_draws(
    random_term = random_term,
    n_rows = nrow(model_matrix),
    posterior = posterior,
    data = data,
    parameters = parameters,
    context = "Bridge-sampling reconstruction"
  )
  if(any(!is.finite(source_draws) | source_draws < 0)){
    .bt_JAGS_marglik_out_of_support(
      "Bridge samples contain out-of-support row-indexed external SD source values for block '",
      random_term$block_name,
      "'."
    )
  }
  column_allocation_draws <- .bt_JAGS_marglik_row_indexed_column_allocation_draws(
    random_term = random_term,
    posterior = posterior,
    prior_list = prior_list,
    n_columns = n_columns
  )
  if(!is.null(column_allocation_draws)){
    if(any(!is.finite(column_allocation_draws) | column_allocation_draws < 0)){
      .bt_JAGS_marglik_out_of_support(
        "Bridge samples contain out-of-support row-indexed external SD allocation values for block '",
        random_term$block_name,
        "'."
      )
    }
    if(isTRUE(independent)){
      contribution <- .bt_random_effect_independent_contribution_from_latent(
        random_term = random_term,
        model_matrix = model_matrix,
        group_map = group_map,
        posterior = posterior,
        column_scale_draws = column_allocation_draws,
        row_scale_draws = source_draws
      )
      if(is.null(contribution)){
        stop(
          "Bridge samples are missing standardized latent random effects for block '",
          random_term$block_name,
          "'.",
          call. = FALSE
        )
      }
      return(as.vector(contribution[, 1L]))
    }
    if(!is.null(unit_contribution)){
      row_column <- .bt_random_effect_structured_indicator_columns(
        model_matrix,
        context = "Bridge reconstruction for a row-indexed scalar-structured random effect"
      )
      return(as.vector(
        unit_contribution[, 1L] *
          source_draws[1L, ] *
          column_allocation_draws[1L, row_column]
      ))
    }
    return(as.vector(.bt_random_effect_apply_row_indexed_source_to_unit_columns(
      unit_columns = unit_columns,
      source_draws = source_draws,
      allocation_draws = column_allocation_draws
    )[, 1L]))
  }
  allocation_draws <- .bt_JAGS_marglik_row_indexed_allocation_draws(
    random_term = random_term,
    posterior = posterior,
    prior_list = prior_list
  )
  if(any(!is.finite(allocation_draws) | allocation_draws < 0)){
    .bt_JAGS_marglik_out_of_support(
      "Bridge samples contain out-of-support row-indexed external SD allocation values for block '",
      random_term$block_name,
      "'."
    )
  }

  if(isTRUE(independent)){
    contribution <- .bt_random_effect_independent_contribution_from_latent(
      random_term = random_term,
      model_matrix = model_matrix,
      group_map = group_map,
      posterior = posterior,
      row_scale_draws = source_draws,
      draw_scale = allocation_draws
    )
    if(is.null(contribution)){
      stop(
        "Bridge samples are missing standardized latent random effects for block '",
        random_term$block_name,
        "'.",
        call. = FALSE
      )
    }
    return(as.vector(contribution[, 1L]))
  }

  if(!is.null(unit_contribution)){
    return(as.vector(
      unit_contribution[, 1L] *
        source_draws[1L, ] *
        allocation_draws[1L]
    ))
  }
  as.vector(.bt_random_effect_apply_row_indexed_source_to_unit_columns(
    unit_columns = unit_columns,
    source_draws = source_draws,
    allocation_draws = allocation_draws
  )[, 1L])
}

.bt_JAGS_marglik_row_indexed_column_allocation_draws <- function(
    random_term, posterior, prior_list, n_columns){

  tryCatch(
    .bt_random_effect_row_indexed_column_allocation_draws(
      random_term = random_term,
      posterior = posterior,
      prior_list = prior_list,
      n_columns = n_columns
    ),
    error = function(e){
      if(.bt_JAGS_marglik_random_effect_support_error(e)){
        .bt_JAGS_marglik_out_of_support(conditionMessage(e))
      }
      stop(e)
    }
  )
}

.bt_JAGS_marglik_row_indexed_allocation_draws <- function(random_term,
                                                          posterior,
                                                          prior_list){

  tryCatch(
    .bt_random_effect_row_indexed_allocation_draws(
      random_term = random_term,
      posterior = posterior,
      prior_list = prior_list
    ),
    error = function(e){
      if(.bt_JAGS_marglik_random_effect_support_error(e)){
        .bt_JAGS_marglik_out_of_support(conditionMessage(e))
      }
      stop(e)
    }
  )
}

.bt_JAGS_marglik_random_effect_support_error <- function(error){

  inherits(error, "BayesTools_random_effect_allocation_out_of_support")
}

.bt_JAGS_marglik_row_indexed_external_sd_source_draws <- function(random_term,
                                                                  n_rows,
                                                                  posterior,
                                                                  data = NULL,
                                                                  parameters = NULL,
  context = "Bridge-sampling reconstruction"){

  source <- .bt_random_effect_row_indexed_source(random_term)
  source_parameter <- source$source
  source_label <- .bt_random_sd_source_label(source)
  source_names <- .bt_parameter_source_row_names(source_parameter, n_rows)
  values_function <- .bt_parameter_source_values_function(source_parameter)
  if(is.null(values_function)){
    return(.bt_random_effect_row_indexed_source_draws(
      random_term = random_term,
      n_rows = n_rows,
      posterior = posterior,
      data = data,
      parameters = parameters,
      context = context
    ))
  }
  present <- intersect(source_names, colnames(posterior))
  if(length(present) > 0L){
    stop(
      context, " for source '", source_label,
      "' is ambiguous: the source provides a values function and the posterior ",
      "also contains sampled row-source column(s): ",
      paste0("'", present[seq_len(min(3L, length(present)))], "'",
             collapse = ", "),
      if(length(present) > 3L) ", ..." else "",
      ". Use one row-source reconstruction path only.",
      call. = FALSE
    )
  }
  if(!is.matrix(posterior)){
    stop("'posterior' must be a matrix.", call. = FALSE)
  }
  check_int(n_rows, "n_rows", lower = 1, allow_NA = FALSE)

  out <- matrix(NA_real_, nrow = nrow(posterior), ncol = n_rows)
  for(draw in seq_len(nrow(posterior))){
    draw_parameters <- .bt_JAGS_marglik_row_indexed_source_parameters(
      posterior = posterior,
      draw = draw,
      parameters = parameters
    )
    draw_parameters <- .bt_parameter_source_guard_parameters(
      draw_parameters,
      source_parameter
    )
    values <- tryCatch(
      values_function(
        parameters = draw_parameters,
        data = data,
        n_rows = n_rows
      ),
      error = function(e)e
    )
    if(inherits(values, "error")){
      stop(
        context, " for source '", source_label,
        "' failed: ", conditionMessage(values),
        call. = FALSE
      )
    }
    if(!is.numeric(values)){
      stop(
        context, " for source '", source_label,
        "' must return a numeric vector.",
        call. = FALSE
      )
    }
    values <- as.numeric(values)
    if(length(values) != n_rows){
      stop(
        context, " for source '", source_label,
        "' must return a numeric vector of length ", n_rows, ".",
        call. = FALSE
      )
    }
    if(anyNA(values)){
      .bt_JAGS_marglik_out_of_support(
        context, " for source '", source_label,
        "' returned missing values."
      )
    }
    out[draw, ] <- values
  }

  colnames(out) <- source_names
  out
}

.bt_JAGS_marglik_row_indexed_source_parameters <- function(posterior,
                                                           draw,
                                                           parameters = NULL){

  if(nrow(posterior) == 1L &&
     !is.null(parameters) &&
     is.list(parameters) &&
     all(colnames(posterior) %in% names(parameters))){
    return(parameters)
  }

  .bt_parameter_source_draw_parameters(
    posterior = posterior,
    draw = draw,
    parameters = parameters
  )
}

.bt_JAGS_marglik_random_effect_sd_values <- function(samples, random_term,
                                                     prior_list){

  .bt_JAGS_marglik_check_random_effect_dirichlet_samples(
    samples = samples,
    random_term = random_term,
    prior_list = prior_list
  )

  if(.bt_random_effect_has_row_indexed_external_sd(random_term)){
    .bt_JAGS_marglik_row_indexed_external_sd_stop(random_term)
  }

  posterior <- .bt_JAGS_marglik_random_effect_posterior_row(samples)
  sd_draws <- tryCatch(
    .bt_random_effect_sd_draws(
      random_term = random_term,
      n_columns = random_term$n_columns,
      posterior = posterior,
      prior_list = prior_list
    ),
    error = function(e){
      if(.bt_JAGS_marglik_random_effect_support_error(e)){
        .bt_JAGS_marglik_out_of_support(conditionMessage(e))
      }
      stop(e)
    }
  )
  if(is.null(sd_draws)){
    .bt_JAGS_marglik_explain_random_effect_sd_missing(samples, random_term, prior_list)
    stop(
      "Random-effect SD metadata are incomplete for block '",
      random_term$block_name,
      "'.",
      call. = FALSE
    )
  }

  unname(sd_draws[1L, ])
}

.bt_JAGS_marglik_row_indexed_external_sd_stop <- function(random_term){

  stop(
    "Random-effect block '",
    random_term$block_name,
    "' uses row-indexed external SD source '",
    .bt_random_effect_external_sd_source_label(random_term),
    "'; scalar random-effect SD values are not defined for row-indexed sources.",
    call. = FALSE
  )
}

.bt_JAGS_marglik_explain_random_effect_sd_missing <- function(samples,
                                                              random_term,
                                                              prior_list){

  binding <- random_term$sd_binding
  if(!is.null(binding) && isTRUE(binding$true_allocation) &&
     length(binding$allocations) > 0L){
    allocation <- binding$allocations[[1L]]
    if(!.bt_JAGS_marglik_random_effect_parameter_available(
      samples = samples,
      parameter_name = .bt_random_sd_binding_source_name(allocation$source),
      prior_list = prior_list
    )){
      stop("'posterior' does not contain all monitored formula prior parameters.", call. = FALSE)
    }
    factors <- if(identical(allocation$target, "sd_component")){
      .bt_random_effect_allocation_parent_factors_metadata(allocation)
    }else{
      .bt_random_effect_allocation_factors_metadata(allocation)
    }
    missing_factor <- vapply(factors, function(factor){
      !.bt_JAGS_marglik_random_effect_dirichlet_available(
        samples = samples,
        parameter_name = factor$weight_name,
        prior_list = prior_list
      )
    }, logical(1))
    if(any(missing_factor)){
      stop(
        "Bridge samples are missing Dirichlet allocation coordinates for parameter '",
        factors[[which(missing_factor)[1L]]]$weight_name,
        "'.",
        call. = FALSE
      )
    }
    if(identical(allocation$target, "sd_component") &&
       !.bt_JAGS_marglik_random_effect_dirichlet_available(
         samples = samples,
         parameter_name = allocation$weight_name,
         prior_list = prior_list
       )){
      stop(
        "Bridge samples are missing Dirichlet allocation coordinates for parameter '",
        allocation$weight_name,
        "'.",
        call. = FALSE
      )
    }
  }else{
    sd_names <- random_term$sd_parameter_names
    if(is.null(sd_names) || length(sd_names) != random_term$n_columns || any(is.na(sd_names))){
      return(invisible(FALSE))
    }
    missing_sd <- !vapply(sd_names, function(parameter_name){
      .bt_JAGS_marglik_random_effect_parameter_available(
        samples = samples,
        parameter_name = parameter_name,
        prior_list = prior_list
      )
    }, logical(1))
    if(any(missing_sd)){
      stop("'posterior' does not contain all monitored formula prior parameters.", call. = FALSE)
    }
  }

  invisible(FALSE)
}

.bt_JAGS_marglik_random_effect_parameter_available <- function(samples,
                                                               parameter_name,
                                                               prior_list){

  if(parameter_name %in% names(samples)){
    return(TRUE)
  }

  prior_name <- sub("\\[[0-9]+\\]$", "", parameter_name)
  if(!prior_name %in% names(prior_list)){
    return(FALSE)
  }

  is.prior.point(prior_list[[prior_name]])
}

.bt_JAGS_marglik_random_effect_dirichlet_available <- function(samples,
                                                               parameter_name,
                                                               prior_list){

  if(!parameter_name %in% names(prior_list)){
    return(FALSE)
  }
  prior <- prior_list[[parameter_name]]
  if(!is.prior.simplex(prior) || !identical(prior$distribution, "dirichlet")){
    return(FALSE)
  }

  K <- prior$parameters[["K"]]
  weight_names <- paste0(parameter_name, "[", seq_len(K), "]")
  eta_names <- paste0(.JAGS_prior_dirichlet_eta_name(parameter_name), "[", seq_len(K), "]")
  all(weight_names %in% names(samples)) || all(eta_names %in% names(samples))
}

.bt_JAGS_marglik_random_effect_cholesky <- function(samples, random_term){

  n_columns <- random_term$n_columns
  L <- .bt_random_effect_cholesky_draws(
    random_term = random_term,
    n_columns = n_columns,
    posterior = .bt_JAGS_marglik_random_effect_posterior_row(samples)
  )
  if(is.null(L)){
    if(.bt_JAGS_marglik_random_effect_correlation_sample_available(samples, random_term)){
      .bt_JAGS_marglik_out_of_support(
        "Bridge samples contain out-of-support random-effect correlation coordinates for block '",
        random_term$block_name,
        "'."
      )
    }
    stop(
      "Bridge samples are missing random-effect correlation coordinates for block '",
      random_term$block_name,
      "'.",
      call. = FALSE
    )
  }

  L[1L, , ]
}

.bt_JAGS_marglik_random_effect_correlation_sample_available <- function(samples,
                                                                        random_term){

  structure <- .bt_random_effect_structure(
    random_term,
    context = "Bridge sampling random-effect metadata"
  )
  if(structure %in% c("diag", "id")){
    return(FALSE)
  }
  correlation <- .bt_random_effect_correlation_metadata(
    random_term,
    structure = structure,
    context = "Bridge sampling random-effect metadata"
  )
  if(is.null(correlation)){
    return(FALSE)
  }
  sample_names <- character()
  if(identical(correlation$type, "rho")){
    sample_names <- c(correlation$rho_name, correlation$sample_name)
  }else if(identical(correlation$type, "lkj")){
    sample_names <- .bt_random_effect_lkj_primitive_names(
      random_term,
      random_term$n_columns,
      context = "Bridge sampling random-effect metadata"
    )
  }
  sample_names <- sample_names[!is.na(sample_names) & nzchar(sample_names)]

  any(sample_names %in% names(samples))
}

.bt_JAGS_marglik_out_of_support <- function(...){

  stop(structure(
    list(message = paste0(...), call = NULL),
    class = c("BayesTools_marglik_out_of_support", "error", "condition")
  ))
}

.bt_JAGS_marglik_random_effect_rho <- function(samples, random_term){

  rho <- .bt_random_effect_rho_draws(
    random_term = random_term,
    posterior = .bt_JAGS_marglik_random_effect_posterior_row(samples),
    context = "Bridge sampling random-effect metadata"
  )
  if(is.null(rho)){
    stop(
      "Bridge samples are missing or invalid scalar correlation coordinates for block '",
      random_term$block_name,
      "'.",
      call. = FALSE
    )
  }

  unname(rho[1L])
}

.bt_JAGS_marglik_random_effect_posterior_row <- function(samples){

  cached <- attr(samples, "BayesTools_marglik_posterior_row", exact = TRUE)
  if(!is.null(cached)){
    return(cached)
  }

  matrix(
    unname(samples),
    nrow = 1L,
    dimnames = list(NULL, names(samples))
  )
}

.JAGS_marglik_parameters_formula_get <- function(samples, parameter, formula_data_list, formula_prior_list, prior_list_parameters, log_intercept = FALSE){

  formula_terms            <- names(formula_prior_list)
  names(formula_data_list) <- sub(paste0("^", parameter, "_data_"), paste0(parameter, "_"), names(formula_data_list))

  # start with intercept
  if(sum(formula_terms == paste0(parameter, "_intercept")) == 1){

    intercept_prior <- formula_prior_list[[paste0(parameter, "_intercept")]]
    .bt_validate_formula_reconstruction_prior(
      intercept_prior,
      paste0(parameter, "_intercept")
    )
    multiply_by <- .bt_JAGS_marglik_prior_multiply_by(
      intercept_prior,
      prior_list_parameters
    )
    intercept_value <- .JAGS_marglik_parameter_values(samples, intercept_prior, paste0(parameter, "_intercept"))
    # apply log transformation if log(intercept) attribute is set
    if(log_intercept){
      intercept_value <- log(intercept_value)
    }
    output <- multiply_by * rep(intercept_value, formula_data_list[[paste0("N_", parameter)]])

  }else{
    output <- rep(0, formula_data_list[[paste0("N_", parameter)]])
  }

  # add the remaining terms
  remaining_terms <- formula_terms[formula_terms != paste0(parameter, "_intercept")]
  if(length(remaining_terms) > 0){
    for(term in remaining_terms){

      .bt_validate_formula_reconstruction_prior(
        formula_prior_list[[term]],
        term
      )
      multiply_by <- .bt_JAGS_marglik_prior_multiply_by(
        formula_prior_list[[term]],
        prior_list_parameters
      )

      if(is.prior.point(formula_prior_list[[term]]) && !is.prior.factor(formula_prior_list[[term]])){

        output <- output + multiply_by * formula_prior_list[[term]][["parameters"]][["location"]] * formula_data_list[[term]]

      }else if(is.prior.point(formula_prior_list[[term]]) && is.prior.factor(formula_prior_list[[term]])){

        if(.get_prior_factor_levels(formula_prior_list[[term]]) == 1){
          output <- output + multiply_by * formula_prior_list[[term]][["parameters"]][["location"]] * formula_data_list[[term]]
        }else{
          output <- output + multiply_by * formula_data_list[[term]] %*% rep(formula_prior_list[[term]][["parameters"]][["location"]], .get_prior_factor_levels(formula_prior_list[[term]]))
        }

      }else if(is.prior.factor(formula_prior_list[[term]])){

        if(.get_prior_factor_levels(formula_prior_list[[term]]) == 1){
          term_value <- .JAGS_marglik_parameter_values(samples, formula_prior_list[[term]], term)
          output     <- output + multiply_by * term_value * formula_data_list[[term]]
        }else{
          term_names  <- paste0(term,"[", 1:.get_prior_factor_levels(formula_prior_list[[term]]), "]")
          term_values <- .JAGS_marglik_parameter_values(samples, formula_prior_list[[term]], term_names)
          output      <- output + multiply_by * formula_data_list[[term]] %*% term_values
        }


      }else if(is.prior.simple(formula_prior_list[[term]])){

        term_value <- .JAGS_marglik_parameter_values(samples, formula_prior_list[[term]], term)
        output     <- output + multiply_by * term_value * formula_data_list[[term]]

      }else{
        stop(
          "Internal formula reconstruction prior dispatch failed for '",
          term, "'.",
          call. = FALSE
        )
      }

    }
  }


  return(as.vector(output))
}
