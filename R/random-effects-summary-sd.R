.bt_random_effect_summary_validate_values <- function(values, name, label,
                                                      n_draws){

  if(is.null(values)){
    stop(
      "Random-effect summary '", label,
      "' could not be computed for derived column '", name, "'.",
      call. = FALSE
    )
  }
  if(length(values) != n_draws){
    stop(
      "Random-effect summary '", label,
      "' returned ", length(values), " draw(s), but expected ", n_draws, ".",
      call. = FALSE
    )
  }
  if(all(is.na(values))){
    stop(
      "Random-effect summary '", label,
      "' contains only missing values.",
      call. = FALSE
    )
  }

  as.numeric(values)
}

.bt_random_effect_summary_complete_scaled_samples <- function(random_term,
                                                              model_samples,
                                                              prior_list,
                                                              parameter = NULL,
                                                              formula_scale = NULL){

  sd_names <- unique(random_term$sd_parameter_names)
  sd_names <- sd_names[!is.na(sd_names)]
  if(length(sd_names) == 0L || is.null(parameter) ||
     is.null(formula_scale) || length(formula_scale) == 0L ||
     is.null(formula_scale[[parameter]]) || length(formula_scale[[parameter]]) == 0L){
    return(model_samples)
  }
  parameter_scale <- formula_scale[[parameter]]
  column_groups <- .random_sd_column_unscale_groups(
    random_sd_cols = sd_names,
    formula_scale = parameter_scale,
    prefix = parameter
  )
  term_map <- if(is.null(column_groups)){
    .random_sd_term_map(sd_names, parameter_scale, parameter)
  }else{
    character()
  }
  if((is.null(column_groups) || length(column_groups) == 0L) &&
     length(term_map) == 0L){
    return(model_samples)
  }

  sd_draws <- .bt_random_effect_sd_draws(
    random_term = random_term,
    n_columns = random_term$n_columns,
    posterior = model_samples,
    prior_list = prior_list
  )
  if(is.null(sd_draws)){
    .bt_random_effect_summary_missing_sd_stop(random_term)
  }

  fallback_columns <- .bt_random_effect_summary_sd_columns(random_term, sd_names)

  completed_sd <- sd_draws[, fallback_columns, drop = FALSE]
  colnames(completed_sd) <- sd_names
  completed <- model_samples
  existing <- sd_names %in% colnames(completed)
  if(any(existing)){
    completed[, sd_names[existing]] <- completed_sd[, existing, drop = FALSE]
  }
  if(any(!existing)){
    completed <- cbind(completed, completed_sd[, !existing, drop = FALSE])
  }
  completed <- .bt_random_effect_summary_complete_correlation_samples(
    random_term = random_term,
    model_samples = completed
  )
  correlation_required_groups <- if(.bt_random_effect_summary_requires_correlation(random_term)){
    random_term$block_name
  }else{
    NULL
  }

  .apply_random_sd_unscale(
    posterior = completed,
    random_sd_cols = sd_names,
    formula_scale = parameter_scale,
    prefix = parameter,
    correlation_required_groups = correlation_required_groups
  )
}

.bt_random_effect_summary_complete_correlation_samples <- function(random_term,
                                                                   model_samples){

  if(!.bt_random_effect_summary_requires_correlation(random_term)){
    return(model_samples)
  }

  L_names <- .bt_random_effect_cholesky_names(
    random_term = random_term,
    n_columns = random_term$n_columns
  )
  L_vector_names <- as.vector(L_names)
  if(all(L_vector_names %in% colnames(model_samples))){
    return(model_samples)
  }

  cholesky <- .bt_random_effect_cholesky_draws(
    random_term = random_term,
    n_columns = random_term$n_columns,
    posterior = model_samples
  )
  if(is.null(cholesky)){
    .bt_random_effect_summary_rho_samples(random_term, model_samples)
    .bt_random_effect_summary_missing_correlation_stop(random_term)
  }

  completed_L <- matrix(NA_real_, nrow = nrow(model_samples), ncol = length(L_vector_names))
  colnames(completed_L) <- L_vector_names
  for(row in seq_len(random_term$n_columns)){
    for(column in seq_len(random_term$n_columns)){
      completed_L[, L_names[row, column]] <- cholesky[, row, column]
    }
  }

  completed <- model_samples
  existing <- L_vector_names %in% colnames(completed)
  if(any(existing)){
    completed[, L_vector_names[existing]] <- completed_L[, existing, drop = FALSE]
  }
  if(any(!existing)){
    completed <- cbind(completed, completed_L[, !existing, drop = FALSE])
  }

  completed
}

.bt_random_effect_summary_requires_correlation <- function(random_term){

  structure <- .bt_random_effect_summary_term_structure(random_term)
  is.numeric(random_term$n_columns) &&
    length(random_term$n_columns) == 1L &&
    !is.na(random_term$n_columns) &&
    random_term$n_columns > 1L &&
    structure %in% c("us", "cs", "hcs", "ar1", "car", "har")
}

.bt_random_effect_summary_allocation_target <- function(allocation){

  target <- allocation$target
  if(is.character(target) && length(target) == 1L && !is.na(target) &&
     target %in% c("block", "sd_component")){
    return(target)
  }

  stop(
    "Random-effect allocation metadata are missing canonical 'allocation$target'.",
    call. = FALSE
  )
}

.bt_random_effect_summary_sd_samples <- function(random_term, model_samples,
                                                 prior_list,
                                                 parameter = NULL,
                                                 formula_scale = NULL){

  if(!is.null(random_term$sd_binding) &&
     .bt_random_sd_binding_has_external_source(random_term$sd_binding) &&
     !isTRUE(random_term$sd_binding$true_allocation)){
    return(list(names = character(), components = character(),
                values = matrix(nrow = nrow(model_samples), ncol = 0L)))
  }

  sd_names <- unique(random_term$sd_parameter_names)
  sd_names <- sd_names[!is.na(sd_names)]
  if(length(sd_names) == 0L){
    return(list(names = character(), components = character(),
                values = matrix(nrow = nrow(model_samples), ncol = 0L)))
  }
  model_samples <- .bt_random_effect_summary_complete_scaled_samples(
    random_term = random_term,
    model_samples = model_samples,
    prior_list = prior_list,
    parameter = parameter,
    formula_scale = formula_scale
  )

  components <- .bt_random_effect_summary_sd_components(random_term, sd_names)
  values <- matrix(NA_real_, nrow = nrow(model_samples), ncol = length(sd_names))
  colnames(values) <- sd_names

  fallback <- NULL
  fallback_used <- rep(FALSE, length(sd_names))
  for(i in seq_along(sd_names)){
    if(sd_names[i] %in% colnames(model_samples)){
      values[, i] <- model_samples[, sd_names[i]]
    }else{
      if(is.null(fallback)){
        fallback <- .bt_random_effect_sd_draws(
          random_term = random_term,
          n_columns = random_term$n_columns,
          posterior = model_samples,
          prior_list = prior_list
        )
      }
      if(is.null(fallback)){
        .bt_random_effect_summary_missing_sd_stop(random_term)
      }else{
        column <- match(sd_names[i], random_term$sd_parameter_names)
        if(is.na(column)){
          .bt_random_effect_summary_inconsistent_sd_stop(random_term)
        }
        values[, i] <- fallback[, column]
        fallback_used[i] <- TRUE
      }
    }
  }
  if(!is.null(fallback) && any(fallback_used)){
    fallback_columns <- .bt_random_effect_summary_sd_columns(random_term, sd_names)
    if(all(fallback_used)){
      values <- fallback[, fallback_columns, drop = FALSE]
      colnames(values) <- sd_names
      values <- .bt_random_effect_summary_unscale_sd_fallback(
        values = values,
        sd_names = sd_names,
        parameter = parameter,
        formula_scale = formula_scale
      )
    }else{
      fallback_values <- fallback[, fallback_columns[fallback_used], drop = FALSE]
      colnames(fallback_values) <- sd_names[fallback_used]
      values[, fallback_used] <- .bt_random_effect_summary_unscale_sd_fallback(
        values = fallback_values,
        sd_names = sd_names[fallback_used],
        parameter = parameter,
        formula_scale = formula_scale
      )
    }
  }

  missing <- colSums(!is.na(values)) == 0L
  if(any(missing)){
    .bt_random_effect_summary_missing_sd_values_stop(
      random_term = random_term,
      sd_names = sd_names[missing]
    )
  }

  list(
    names = sd_names,
    components = components,
    values = values
  )
}

.bt_random_effect_summary_inclusion_samples <- function(random_term,
                                                        model_samples,
                                                        prior_list,
                                                        parameter){

  prior_names <- .bt_random_effect_summary_inclusion_prior_names(
    random_term,
    prior_list
  )
  if(length(prior_names) == 0L){
    return(list(
      names = character(),
      labels = character(),
      component_labels = character(),
      types = character(),
      components = character(),
      values = matrix(nrow = nrow(model_samples), ncol = 0L)
    ))
  }

  names <- labels <- component_labels <- types <- components <- character()
  values <- list()
  owner <- .bt_random_effect_public_name(random_term)

  for(prior_name in prior_names){
    prior <- prior_list[[prior_name]]
    indicator_name <- paste0(prior_name, "_indicator")
    if(!indicator_name %in% colnames(model_samples)){
      next
    }
    indicator <- as.numeric(model_samples[, indicator_name])
    prior_components <- attr(prior, "components")
    if(is.prior.spike_and_slab(prior)){
      alternative_index <- which(prior_components == "alternative")
      if(length(alternative_index) != 1L){
        next
      }
      names <- c(names, .bt_random_effect_summary_name(
        parameter = parameter,
        type = "inclusion",
        parts = c(random_term$block_name, .bt_random_effect_summary_safe_label(prior_name))
      ))
      display_effect <- .bt_random_effect_summary_inclusion_effect_label(
        random_term,
        prior_name
      )
      effect <- if(identical(display_effect, owner)){
        "sd"
      }else{
        paste0("sd(", display_effect, ")")
      }
      labels <- c(labels, .bt_random_effect_semantic_name(
        parameter = "",
        owner = owner,
        quantity = "inclusion",
        arguments = effect,
        formula_prefix = FALSE
      ))
      component_labels <- c(
        component_labels,
        paste0("inclusion(", effect, ")")
      )
      types <- c(types, "inclusion")
      components <- c(components, "alternative")
      values[[length(values) + 1L]] <- as.numeric(indicator %in% alternative_index)
    }else if(is.prior.mixture(prior)){
      unique_components <- unique(prior_components)
      for(component in unique_components){
        component_index <- which(prior_components == component)
        names <- c(names, .bt_random_effect_summary_name(
          parameter = parameter,
          type = "inclusion",
          parts = c(
            random_term$block_name,
            .bt_random_effect_summary_safe_label(prior_name),
            component
          )
        ))
        display_effect <- .bt_random_effect_summary_inclusion_effect_label(
          random_term,
          prior_name
        )
        effect <- if(identical(display_effect, owner)){
          "sd"
        }else{
          paste0("sd(", display_effect, ")")
        }
        argument <- paste0(effect, "[", component, "]")
        labels <- c(labels, .bt_random_effect_semantic_name(
          parameter = "",
          owner = owner,
          quantity = "inclusion",
          arguments = argument,
          formula_prefix = FALSE
        ))
        component_labels <- c(
          component_labels,
          paste0("inclusion(", argument, ")")
        )
        types <- c(types, "inclusion")
        components <- c(components, component)
        values[[length(values) + 1L]] <- as.numeric(indicator %in% component_index)
      }
    }
  }

  if(length(values) == 0L){
    value_matrix <- matrix(nrow = nrow(model_samples), ncol = 0L)
  }else{
    value_matrix <- do.call(cbind, values)
    colnames(value_matrix) <- names
  }

  list(
    names = names,
    labels = labels,
    component_labels = component_labels,
    types = types,
    components = components,
    values = value_matrix
  )
}

.bt_random_effect_summary_inclusion_effect_label <- function(random_term,
                                                            prior_name){

  sd_names <- unique(random_term$sd_parameter_names)
  sd_names <- sd_names[!is.na(sd_names)]
  sd_prior_names <- sub("[[][^][]+[]]$", "", sd_names)
  component_terms <- NULL
  if(!is.null(random_term$sd_binding) &&
     !is.null(random_term$sd_binding$sd_component_terms)){
    component_terms <- random_term$sd_binding$sd_component_terms
  }

  if(!is.null(component_terms) && length(component_terms) > 0L){
    matched <- component_terms[sd_prior_names == prior_name]
    matched <- unname(matched[!is.na(matched) & nzchar(matched)])
    if(length(matched) > 0L){
      base_terms <- unique(sub("[[][^][]+[]]$", "", matched))
      if(length(base_terms) == 1L && nzchar(base_terms)){
        return(base_terms)
      }
      return(paste(matched, collapse = ", "))
    }
  }

  .bt_random_effect_public_name(random_term)
}

.bt_random_effect_summary_inclusion_prior_names <- function(random_term,
                                                            prior_list){

  sd_names <- unique(random_term$sd_parameter_names)
  sd_names <- sd_names[!is.na(sd_names)]
  if(length(sd_names) == 0L || length(prior_list) == 0L){
    return(character())
  }

  prior_names <- unique(sub("[[][^][]+[]]$", "", sd_names))
  prior_names <- prior_names[prior_names %in% names(prior_list)]
  prior_names[vapply(prior_names, function(prior_name){
    prior <- prior_list[[prior_name]]
    is.prior.spike_and_slab(prior) || is.prior.mixture(prior)
  }, logical(1))]
}

.bt_random_effect_summary_sd_columns <- function(random_term, sd_names){

  columns <- match(sd_names, random_term$sd_parameter_names)
  if(anyNA(columns)){
    .bt_random_effect_summary_inconsistent_sd_stop(random_term)
  }

  columns
}

.bt_random_effect_summary_missing_sd_stop <- function(random_term){

  stop(
    "Random-effect summary samples are missing canonical SD coordinates for block '",
    random_term$block_name,
    "'. Expected monitored SD, point-prior, or allocation coordinates.",
    call. = FALSE
  )
}

.bt_random_effect_summary_missing_sd_values_stop <- function(random_term,
                                                             sd_names){

  stop(
    "Random-effect summary SD samples contain only missing values for block '",
    random_term$block_name,
    "' coordinate(s): ",
    paste0("'", sd_names, "'", collapse = ", "),
    ".",
    call. = FALSE
  )
}

.bt_random_effect_summary_missing_allocation_sd_stop <- function(allocation){

  label <- allocation
  if(!is.character(label) || length(label) != 1L || is.na(label) || !nzchar(label)){
    label <- "<unknown>"
  }

  stop(
    "Random-effect allocation-scale summary samples are missing canonical SD coordinates for allocation '",
    label,
    "'. Expected monitored allocation SD or point-prior coordinates.",
    call. = FALSE
  )
}

.bt_random_effect_summary_inconsistent_sd_stop <- function(random_term){

  stop(
    "Random-effect summary metadata are inconsistent for block '",
    random_term$block_name,
    "': semantic SD names do not match canonical 'random_term$sd_parameter_names'.",
    call. = FALSE
  )
}
