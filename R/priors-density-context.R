.prior_density_context <- function(prior_list, column_names, formula_scale = NULL,
                                   n_grid = .prior_linear_density_default_grid(),
                                   tail_prob = .prior_linear_density_tail_prob()){

  check_list(prior_list, "prior_list")
  check_char(column_names, "column_names", check_length = FALSE)
  check_list(formula_scale, "formula_scale", allow_NULL = TRUE)
  check_int(n_grid, "n_grid", lower = 16)
  check_real(tail_prob, "tail_prob", lower = 0, upper = 0.5, allow_bound = FALSE)

  if(!is.null(formula_scale) && length(formula_scale) > 0){
    .check_formula_scale_info(formula_scale)
  }

  transforms <- list()
  if(!is.null(formula_scale) && length(formula_scale) > 0){
    for(param_name in names(formula_scale)){
      affected_cols <- grep(paste0("^", param_name, "_"), column_names, value = TRUE)
      if(length(affected_cols) == 0){
        next
      }
      transforms[[param_name]] <- list(
        columns       = affected_cols,
        matrix        = .build_unscale_matrix(affected_cols, formula_scale[[param_name]], param_name),
        log_intercept = isTRUE(attr(formula_scale[[param_name]], "log_intercept")),
        intercept     = paste0(param_name, "_intercept")
      )
    }
  }

  out <- list(
    prior_list    = prior_list,
    column_names  = column_names,
    formula_scale = formula_scale,
    transforms    = transforms,
    n_grid        = n_grid,
    tail_prob     = tail_prob
  )
  class(out) <- "prior_density_context"
  return(out)
}

.prior_density_context_standardized_weights <- function(context, weights){

  if(!inherits(context, "prior_density_context")){
    stop("'context' must be a prior density context.", call. = FALSE)
  }
  if(is.null(names(weights))){
    stop("'weights' must be a named numeric vector.", call. = FALSE)
  }

  weights <- weights[intersect(names(weights), context$column_names)]
  out <- rep(0, length(context$column_names))
  names(out) <- context$column_names

  used <- rep(FALSE, length(weights))
  names(used) <- names(weights)

  for(transform in context$transforms){
    cols <- intersect(transform$columns, names(weights))
    if(length(cols) == 0){
      next
    }

    if(transform$log_intercept &&
       transform$intercept %in% cols &&
       abs(weights[[transform$intercept]]) > .prior_linear_density_zero_tol()){
      stop(
        "Linear-combination prior densities with log-intercept scaling are only available ",
        "for the transformed intercept coefficient itself.",
        call. = FALSE
      )
    }

    temp_weights <- rep(0, length(transform$columns))
    names(temp_weights) <- transform$columns
    temp_weights[cols] <- weights[cols]
    out[transform$columns] <- out[transform$columns] +
      as.numeric(temp_weights %*% transform$matrix)
    used[cols] <- TRUE
  }

  remaining <- names(weights)[!used]
  if(length(remaining) > 0){
    out[remaining] <- out[remaining] + weights[remaining]
  }

  out[abs(out) > .prior_linear_density_zero_tol()]
}

.prior_density_context_density <- function(context, weights,
                                           source_transforms = NULL,
                                           output_transformation = NULL,
                                           output_transformation_arguments = NULL){

  standardized_weights <- .prior_density_context_standardized_weights(context, weights)

  if(!is.null(source_transforms)){
    source_transforms <- source_transforms[names(standardized_weights)]
  }

  .prior_linear_combination_density(
    prior_list                       = context$prior_list,
    weights                          = standardized_weights,
    n_grid                           = context$n_grid,
    tail_prob                        = context$tail_prob,
    source_transforms                = source_transforms,
    output_transformation            = output_transformation,
    output_transformation_arguments  = output_transformation_arguments
  )
}

.prior_density_model_mixture_context <- function(prior_list, column_names,
                                                  n_grid = .prior_linear_density_default_grid(),
                                                  tail_prob = .prior_linear_density_tail_prob()){

  check_list(prior_list, "prior_list")
  check_char(column_names, "column_names", check_length = FALSE)

  prior_weights <- do.call(cbind, lapply(prior_list, function(parameter_priors){
    if(is.prior(parameter_priors)){
      return(.prior_model_weight(parameter_priors))
    }
    sapply(parameter_priors, .prior_model_weight)
  }))

  if(!all(prior_weights[, 1] == prior_weights)){
    stop("The model prior distributions are not aligned across parameters.", call. = FALSE)
  }

  model_weights <- prior_weights[, 1]
  model_weights <- model_weights / sum(model_weights)

  out <- list(
    prior_list   = prior_list,
    column_names = column_names,
    model_weights = model_weights,
    n_grid       = n_grid,
    tail_prob    = tail_prob
  )
  class(out) <- "prior_density_model_mixture_context"
  return(out)
}

.prior_density_condition_component <- function(prior){

  if(is.prior.spike_and_slab(prior)){
    components <- attr(prior, "components")
    if(!all(components %in% c("null", "alternative"))){
      stop("conditional mixture posterior distributions are available only for 'null' and 'alternative' components", call. = FALSE)
    }

    inclusion <- mean(.get_spike_and_slab_inclusion(prior))
    inclusion <- min(max(inclusion, 0), 1)
    probabilities <- ifelse(components == "alternative", inclusion, 1 - inclusion)

    return(lapply(seq_along(prior), function(i){
      list(
        prior       = prior[[i]],
        probability = probabilities[i],
        alternative = components[i] == "alternative"
      )
    }))
  }

  if(is.prior.mixture(prior)){
    components <- attr(prior, "components")
    if(!all(components %in% c("null", "alternative"))){
      stop("conditional mixture posterior distributions are available only for 'null' and 'alternative' components", call. = FALSE)
    }

    prior_weights <- attr(prior, "prior_weights")
    prior_weights <- prior_weights / sum(prior_weights)

    return(lapply(seq_along(prior), function(i){
      list(
        prior       = prior[[i]],
        probability = prior_weights[i],
        alternative = components[i] == "alternative"
      )
    }))
  }

  list(list(
    prior       = prior,
    probability = 1,
    alternative = TRUE
  ))
}

.prior_density_copy_parent_attributes <- function(component, parent){

  parent_attributes <- attributes(parent)
  skip <- c("class", "names", "components", "prior_weights", "inclusion_prior")

  for(attribute in setdiff(names(parent_attributes), skip)){
    if(is.null(attr(component, attribute, exact = TRUE))){
      attr(component, attribute) <- parent_attributes[[attribute]]
    }
  }

  component
}

.prior_density_condition_models <- function(prior_list, conditional, conditional_rule,
                                            condition_event = NULL){

  if(is.null(condition_event)){
    condition_event <- .condition_event(
      prior_list        = prior_list,
      conditional       = conditional,
      conditional_rule  = conditional_rule
    )
  }

  .condition_event_model_options(prior_list, condition_event)
}

.prior_density_conditional_context <- function(prior_list, column_names, conditional,
                                               conditional_rule = "AND", formula_scale = NULL,
                                               n_grid = .prior_linear_density_default_grid(),
                                               tail_prob = .prior_linear_density_tail_prob(),
                                               condition_event = NULL){

  if(is.null(condition_event)){
    condition_event <- .condition_event(
      prior_list        = prior_list,
      conditional       = conditional,
      conditional_rule  = conditional_rule
    )
  }
  condition_models <- .prior_density_condition_models(
    prior_list        = prior_list,
    conditional       = conditional,
    conditional_rule  = conditional_rule,
    condition_event   = condition_event
  )
  if(is.null(condition_models)){
    return(.prior_density_context(prior_list, column_names, formula_scale, n_grid, tail_prob))
  }

  out <- list(
    prior_list      = prior_list,
    column_names    = column_names,
    formula_scale   = formula_scale,
    conditional     = condition_event[["conditional"]],
    conditional_rule = condition_event[["conditional_rule"]],
    condition_event = condition_event,
    condition_key   = condition_event[["condition_key"]],
    prior_lists     = condition_models$prior_lists,
    model_weights   = condition_models$weights,
    n_grid          = n_grid,
    tail_prob       = tail_prob
  )
  class(out) <- "prior_density_conditional_context"
  return(out)
}

.prior_density_model_mixture_density <- function(context, weights,
                                                  output_transformation = NULL,
                                                  output_transformation_arguments = NULL){

  if(!inherits(context, "prior_density_model_mixture_context")){
    stop("'context' must be a prior density model-mixture context.", call. = FALSE)
  }

  dists <- vector("list", length(context$model_weights))
  for(model_i in seq_along(context$model_weights)){
    model_prior_list <- lapply(context$prior_list, function(parameter_priors){
      if(is.prior(parameter_priors)){
        return(parameter_priors)
      }
      parameter_priors[[model_i]]
    })
    names(model_prior_list) <- names(context$prior_list)

    for(parameter in names(model_prior_list)){
      if(is.null(model_prior_list[[parameter]])){
        model_prior_list[[parameter]] <- prior("point", list(location = 0))
      }
    }

    dists[[model_i]] <- .prior_linear_combination_density(
      prior_list = model_prior_list,
      weights    = weights,
      n_grid     = context$n_grid,
      tail_prob  = context$tail_prob
    )
  }

  dx <- min(vapply(dists, function(dist){
    if(!is.null(dist$density) && length(dist$density$x) > 1){
      return(dist$density$x[2] - dist$density$x[1])
    }
    Inf
  }, numeric(1)))
  if(!is.finite(dx)){
    dx <- NA_real_
  }

  dist <- .prior_linear_density_mix(
    dists   = dists,
    weights = context$model_weights,
    dx      = dx,
    n_grid  = context$n_grid
  )

  .prior_linear_density_transform(dist, output_transformation,
                                  output_transformation_arguments,
                                  n_grid = context$n_grid)
}

.prior_density_conditional_context_density <- function(context, weights,
                                                       source_transforms = NULL,
                                                       output_transformation = NULL,
                                                       output_transformation_arguments = NULL){

  if(!inherits(context, "prior_density_conditional_context")){
    stop("'context' must be a conditional prior density context.", call. = FALSE)
  }

  if(length(context$prior_lists) == 0){
    stop("No prior models remain after applying the conditional event.", call. = FALSE)
  }

  dists <- lapply(context$prior_lists, function(prior_list){
    if(!is.null(context$formula_scale) && length(context$formula_scale) > 0){
      component_context <- .prior_density_context(
        prior_list    = prior_list,
        column_names  = context$column_names,
        formula_scale = context$formula_scale,
        n_grid        = context$n_grid,
        tail_prob     = context$tail_prob
      )
      return(.prior_density_context_density(
        context           = component_context,
        weights           = weights,
        source_transforms = source_transforms
      ))
    }

    .prior_linear_combination_density(
      prior_list        = prior_list,
      weights           = weights,
      n_grid            = context$n_grid,
      tail_prob         = context$tail_prob,
      source_transforms = source_transforms
    )
  })

  dx <- min(vapply(dists, function(dist){
    if(!is.null(dist$density) && length(dist$density$x) > 1){
      return(dist$density$x[2] - dist$density$x[1])
    }
    Inf
  }, numeric(1)))
  if(!is.finite(dx)){
    dx <- NA_real_
  }

  dist <- .prior_linear_density_mix(
    dists   = dists,
    weights = context$model_weights,
    dx      = dx,
    n_grid  = context$n_grid
  )

  .prior_linear_density_transform(dist, output_transformation,
                                  output_transformation_arguments,
                                  n_grid = context$n_grid)
}

.prior_density_build_context <- function(prior_list, column_names, formula_scale = NULL,
                                         n_grid = .prior_linear_density_default_grid(),
                                         tail_prob = .prior_linear_density_tail_prob(),
                                         conditional = NULL,
                                         conditional_rule = "AND",
                                         condition_event = NULL){

  if(is.null(condition_event)){
    condition_event <- .condition_event(
      prior_list        = prior_list,
      conditional       = conditional,
      conditional_rule  = conditional_rule
    )
  }

  if(length(condition_event[["conditional"]]) > 0){
    return(.prior_density_conditional_context(
      prior_list       = prior_list,
      column_names     = column_names,
      conditional      = condition_event[["conditional"]],
      conditional_rule = conditional_rule,
      formula_scale    = formula_scale,
      n_grid           = n_grid,
      tail_prob        = tail_prob,
      condition_event  = condition_event
    ))
  }

  if(all(vapply(prior_list, is.prior, logical(1)))){
    return(.prior_density_context(prior_list, column_names, formula_scale, n_grid, tail_prob))
  }

  if(!is.null(formula_scale) && length(formula_scale) > 0){
    stop("Formula-scale prior densities for model-list mixtures are not implemented.", call. = FALSE)
  }

  .prior_density_model_mixture_context(prior_list, column_names, n_grid, tail_prob)
}

.prior_density_from_context <- function(context, weights,
                                        source_transforms = NULL,
                                        output_transformation = NULL,
                                        output_transformation_arguments = NULL){

  if(inherits(context, "prior_density_context")){
    return(.prior_density_context_density(
      context                         = context,
      weights                         = weights,
      source_transforms               = source_transforms,
      output_transformation           = output_transformation,
      output_transformation_arguments = output_transformation_arguments
    ))
  }

  if(inherits(context, "prior_density_model_mixture_context")){
    if(!is.null(source_transforms)){
      stop("Source transformations are not supported for model-list prior mixtures.", call. = FALSE)
    }
    return(.prior_density_model_mixture_density(
      context                         = context,
      weights                         = weights,
      output_transformation           = output_transformation,
      output_transformation_arguments = output_transformation_arguments
    ))
  }

  if(inherits(context, "prior_density_conditional_context")){
    return(.prior_density_conditional_context_density(
      context                         = context,
      weights                         = weights,
      source_transforms               = source_transforms,
      output_transformation           = output_transformation,
      output_transformation_arguments = output_transformation_arguments
    ))
  }

  stop("Unknown prior density context.", call. = FALSE)
}

.prior_density_from_context_rows <- function(context, weights,
                                             source_transforms = NULL,
                                             output_transformation = NULL,
                                             output_transformation_arguments = NULL){

  if(is.null(dim(weights))){
    return(.prior_density_from_context(
      context                         = context,
      weights                         = weights,
      source_transforms               = source_transforms,
      output_transformation           = output_transformation,
      output_transformation_arguments = output_transformation_arguments
    ))
  }

  weights <- as.matrix(weights)
  if(nrow(weights) == 0){
    return(.prior_linear_density_point(0))
  }
  if(nrow(weights) == 1){
    return(.prior_density_from_context(
      context                         = context,
      weights                         = weights[1, ],
      source_transforms               = source_transforms,
      output_transformation           = output_transformation,
      output_transformation_arguments = output_transformation_arguments
    ))
  }

  row_keys <- apply(signif(weights, 15), 1, paste0, collapse = "\r")
  unique_keys <- unique(row_keys)
  row_counts <- tabulate(match(row_keys, unique_keys), nbins = length(unique_keys))
  row_indices <- match(unique_keys, row_keys)

  dists <- lapply(row_indices, function(row_i){
    .prior_density_from_context(
      context                         = context,
      weights                         = weights[row_i, ],
      source_transforms               = source_transforms,
      output_transformation           = output_transformation,
      output_transformation_arguments = output_transformation_arguments
    )
  })

  dx <- min(vapply(dists, function(dist){
    if(!is.null(dist$density) && length(dist$density$x) > 1){
      return(dist$density$x[2] - dist$density$x[1])
    }
    Inf
  }, numeric(1)))
  if(!is.finite(dx)){
    dx <- NA_real_
  }

  .prior_linear_density_mix(
    dists   = dists,
    weights = row_counts,
    dx      = dx,
    n_grid  = if(!is.null(context$n_grid)) context$n_grid else NULL
  )
}

.prior_density_coefficient_weights <- function(column_names, parameter){

  weights <- rep(0, length(column_names))
  names(weights) <- column_names
  if(parameter %in% names(weights)){
    weights[[parameter]] <- 1
  }
  weights
}
