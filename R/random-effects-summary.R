# Internal random-effect summary helpers.

.bt_random_effect_summary_samples <- function(model_samples, prior_list,
                                              formula_design = NULL,
                                              mode = c("standard", "full", "raw", "none"),
                                              formula_scale = NULL){

  mode <- match.arg(mode)
  if(is.null(formula_design) || length(formula_design) == 0L || identical(mode, "raw")){
    return(list(model_samples = model_samples, prior_list = prior_list))
  }

  random_design <- .bt_random_effect_summary_designs(formula_design)
  if(length(random_design) == 0L){
    return(list(model_samples = model_samples, prior_list = prior_list))
  }

  model_samples <- as.matrix(model_samples)

  if(!identical(mode, "none")){
    derived <- .bt_random_effect_summary_derived_samples(
      model_samples = model_samples,
      prior_list = prior_list,
      random_design = random_design,
      mode = mode,
      formula_scale = formula_scale
    )
    if(ncol(derived$model_samples) > 0L){
      model_samples <- cbind(model_samples, derived$model_samples)
      prior_list <- c(prior_list, derived$prior_list)
    }
  }

  .bt_random_effect_summary_remove_raw(
    model_samples = model_samples,
    prior_list = prior_list,
    random_design = random_design
  )
}

.bt_random_effect_summary_designs <- function(formula_design){

  if(inherits(formula_design, "BayesTools_formula_design")){
    formula_design <- list(formula_design)
  }
  if(!is.list(formula_design)){
    return(list())
  }

  formula_design[vapply(formula_design, .bt_formula_design_has_random_effects, logical(1))]
}

.bt_random_effect_summary_derived_samples <- function(model_samples, prior_list,
                                                      random_design, mode,
                                                      formula_scale = NULL){

  columns <- list()
  summary_priors <- list()
  used_names <- character()

  add_summary <- function(name, values, parameter, type, label,
                          block = NULL, grouping = NULL,
                          structure = NULL, effect_label = NULL,
                          allocation = NULL,
                          component = NULL){
    if(is.null(values) || length(values) != nrow(model_samples) || all(is.na(values))){
      return(invisible(NULL))
    }
    name <- .bt_random_effect_summary_unique_name(name, used_names)
    used_names <<- c(used_names, name)
    columns[[name]] <<- as.numeric(values)
    summary_priors[[name]] <<- .bt_random_effect_summary_prior(
      parameter = parameter,
      type = type,
      label = label,
      block = block,
      grouping = grouping,
      structure = structure,
      effect_label = effect_label,
      allocation = allocation,
      component = component
    )
    invisible(NULL)
  }

  for(prior_name in names(prior_list)){
    prior <- prior_list[[prior_name]]
    if(isTRUE(attr(prior, "random_sd_total"))){
      allocation <- attr(prior, "random_allocation")
      values <- .bt_random_effect_parameter_draws(prior_name, model_samples, prior_list)
      add_summary(
        name = .bt_random_effect_summary_name(
          parameter = attr(prior, "parameter"),
          type = "sd_total",
          parts = allocation
        ),
        values = values,
        parameter = attr(prior, "parameter"),
        type = "sd_total",
        label = paste0("sd_total(", allocation, ")"),
        allocation = allocation
      )
    }
  }

  seen_allocations <- character()
  add_allocation_summary <- function(allocation, parameter, block = NULL,
                                     grouping = NULL, random_term = NULL){
    if(is.null(allocation) || allocation$weight_name %in% seen_allocations){
      return(invisible(NULL))
    }
    allocation_summary <- .bt_random_effect_summary_allocation_samples(
      allocation = allocation,
      random_term = random_term,
      model_samples = model_samples,
      prior_list = prior_list,
      include_multipliers = identical(mode, "full")
    )
    allocation_target <- .bt_random_effect_summary_allocation_target(allocation)
    for(i in seq_along(allocation_summary$names)){
      add_summary(
        name = allocation_summary$names[i],
        values = allocation_summary$values[, i],
        parameter = parameter,
        type = allocation_summary$types[i],
        label = allocation_summary$labels[i],
        block = if(identical(allocation_target, "sd")) block else NULL,
        grouping = if(identical(allocation_target, "sd")) grouping else NULL,
        structure = if(identical(allocation_target, "sd") && !is.null(random_term)) {
          .bt_random_effect_summary_term_structure(random_term)
        }else{
          NULL
        },
        effect_label = if(identical(allocation_target, "sd") && !is.null(random_term)) .bt_random_effect_public_name(random_term) else NULL,
        allocation = allocation$label,
        component = allocation_summary$components[i]
      )
    }
    seen_allocations <<- c(seen_allocations, allocation$weight_name)
    invisible(NULL)
  }

  for(design in random_design){
    parameter <- design$parameter
    for(random_term in design$random_effects){
      display_group <- .bt_random_effect_summary_group_label(random_term)
      display_structure <- .bt_random_effect_summary_term_structure(random_term)
      sd_summary <- .bt_random_effect_summary_sd_samples(
        random_term = random_term,
        model_samples = model_samples,
        prior_list = prior_list,
        parameter = parameter,
        formula_scale = formula_scale
      )
      correlation_model_samples <- .bt_random_effect_summary_complete_scaled_samples(
        random_term = random_term,
        model_samples = model_samples,
        prior_list = prior_list,
        parameter = parameter,
        formula_scale = formula_scale
      )
      for(i in seq_along(sd_summary$names)){
        add_summary(
          name = .bt_random_effect_summary_name(
            parameter = parameter,
            type = "sd",
            parts = c(random_term$block_name, sd_summary$components[i])
          ),
          values = sd_summary$values[, i],
          parameter = parameter,
          type = "sd",
          label = paste0("sd(", sd_summary$components[i], " | ", display_group, ")"),
          block = random_term$block_name,
          grouping = random_term$group_label,
          structure = display_structure,
          effect_label = .bt_random_effect_public_name(random_term),
          component = sd_summary$components[i]
        )
      }

      rho <- .bt_random_effect_rho_draws(random_term, model_samples)
      if(!is.null(rho)){
        add_summary(
          name = .bt_random_effect_summary_name(
            parameter = parameter,
            type = "rho",
            parts = random_term$block_name
          ),
          values = rho,
          parameter = parameter,
          type = "rho",
          label = paste0("rho(", display_group, ")"),
          block = random_term$block_name,
          grouping = random_term$group_label,
          structure = display_structure,
          effect_label = .bt_random_effect_public_name(random_term)
        )
      }

      correlation_summary <- .bt_random_effect_summary_correlation_samples(
        random_term = random_term,
        model_samples = correlation_model_samples
      )
      for(i in seq_along(correlation_summary$labels)){
        add_summary(
          name = .bt_random_effect_summary_name(
            parameter = parameter,
            type = "cor",
            parts = c(random_term$block_name, correlation_summary$parts[[i]])
          ),
          values = correlation_summary$values[, i],
          parameter = parameter,
          type = "cor",
          label = paste0("cor(", correlation_summary$labels[i], " | ", display_group, ")"),
          block = random_term$block_name,
          grouping = random_term$group_label,
          structure = display_structure,
          effect_label = .bt_random_effect_public_name(random_term),
          component = correlation_summary$labels[i]
        )
      }

      allocation <- random_term$allocation
      add_allocation_summary(
        allocation = allocation,
        parameter = parameter,
        block = random_term$block_name,
        grouping = random_term$group_label,
        random_term = random_term
      )
    }
    for(allocation in design$random_allocations){
      add_allocation_summary(
        allocation = allocation,
        parameter = parameter
      )
    }
  }

  if(length(columns) == 0L){
    summary_matrix <- matrix(nrow = nrow(model_samples), ncol = 0L)
  }else{
    summary_matrix <- do.call(cbind, columns)
  }
  list(model_samples = summary_matrix, prior_list = summary_priors)
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

  .apply_random_sd_unscale(
    posterior = completed,
    random_sd_cols = sd_names,
    formula_scale = formula_scale[[parameter]],
    prefix = parameter,
    correlation_required_groups = random_term$block_name
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

  components <- allocation$components
  if(is.character(components) && length(components) == 1L &&
     !is.na(components) && components %in% c("block", "sd")){
    return(components)
  }
  if(is.list(components)){
    target <- allocation$target
    if(is.character(target) && length(target) == 1L && !is.na(target) &&
       target %in% c("block", "sd")){
      return(target)
    }
    stop(
      "Random-effect allocation metadata are missing canonical 'allocation$target'.",
      call. = FALSE
    )
  }

  stop(
    "Random-effect allocation metadata are missing canonical 'allocation$components'.",
    call. = FALSE
  )
}

.bt_random_effect_summary_sd_samples <- function(random_term, model_samples,
                                                 prior_list,
                                                 parameter = NULL,
                                                 formula_scale = NULL){

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

  keep <- colSums(!is.na(values)) > 0L
  list(
    names = sd_names[keep],
    components = components[keep],
    values = values[, keep, drop = FALSE]
  )
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

.bt_random_effect_summary_inconsistent_sd_stop <- function(random_term){

  stop(
    "Random-effect summary metadata are inconsistent for block '",
    random_term$block_name,
    "': semantic SD names do not match canonical 'random_term$sd_parameter_names'.",
    call. = FALSE
  )
}

.bt_random_effect_summary_missing_correlation_stop <- function(random_term){

  stop(
    "Random-effect summary samples are missing or invalid canonical correlation coordinates for block '",
    random_term$block_name,
    "'. Expected monitored Cholesky, LKJ primitive, or scalar correlation coordinates.",
    call. = FALSE
  )
}

.bt_random_effect_summary_sd_components <- function(random_term, sd_names){

  leaves <- random_term$sd_leaves
  if(!is.null(leaves) && !is.null(leaves$leaf_terms)){
    out <- unname(leaves$leaf_terms[sd_names])
    missing <- is.na(out)
    if(any(missing)){
      out[missing] <- .bt_random_effect_summary_component_from_sd_name(
        random_term,
        sd_names[missing]
      )
    }
    return(.bt_random_effect_summary_display_components(random_term, out))
  }

  allocation <- random_term$allocation
  if(!is.null(allocation) && identical(allocation$components, "sd") &&
     !is.null(allocation$leaf_terms)){
    out <- unname(allocation$leaf_terms[sd_names])
    missing <- is.na(out)
    if(any(missing)){
      out[missing] <- .bt_random_effect_summary_component_from_sd_name(
        random_term,
        sd_names[missing]
      )
    }
    return(.bt_random_effect_summary_display_components(random_term, out))
  }

  .bt_random_effect_summary_display_components(
    random_term,
    .bt_random_effect_summary_component_from_sd_name(random_term, sd_names)
  )
}

.bt_random_effect_summary_component_from_sd_name <- function(random_term,
                                                            sd_names){

  prefix <- paste0(random_term$parameter_stem, "_")
  out <- sub(paste0("^", prefix), "", sd_names)
  out <- gsub("__xXx__", ":", out, fixed = TRUE)
  .bt_random_effect_summary_normalize_components(out)
}

.bt_random_effect_summary_normalize_components <- function(components){

  components <- gsub("__xXx__", ":", components, fixed = TRUE)
  components[components == "sd"] <- "shared"
  components[components == "(Intercept)"] <- "intercept"
  components
}

.bt_random_effect_summary_display_components <- function(random_term, components){

  components <- .bt_random_effect_summary_normalize_components(components)

  index <- random_term$structured_index
  if(is.null(index) || is.null(index$name) || is.null(index$label) ||
     identical(index$name, index$label)){
    return(components)
  }

  index_name <- as.character(index$name)
  index_label <- as.character(index$label)
  replace <- components == index_name | startsWith(components, paste0(index_name, "["))
  components[replace] <- paste0(
    index_label,
    substr(components[replace], nchar(index_name) + 1L, nchar(components[replace]))
  )
  components
}

.bt_random_effect_summary_unscale_sd_fallback <- function(values, sd_names,
                                                          parameter,
                                                          formula_scale){

  if(is.null(parameter) || is.null(formula_scale) || length(formula_scale) == 0L ||
     is.null(formula_scale[[parameter]]) || length(formula_scale[[parameter]]) == 0L){
    return(values)
  }

  transformed <- .apply_random_sd_unscale(
    posterior = values,
    random_sd_cols = sd_names,
    formula_scale = formula_scale[[parameter]],
    prefix = parameter
  )
  transformed[, sd_names, drop = FALSE]
}

.bt_random_effect_summary_correlation_samples <- function(random_term,
                                                          model_samples){

  out <- list(
    labels = character(),
    parts = list(),
    values = matrix(nrow = nrow(model_samples), ncol = 0L)
  )
  if(random_term$n_columns < 2L){
    return(out)
  }

  structure <- .bt_random_effect_summary_term_structure(random_term)
  correlation <- .bt_random_effect_correlation_metadata(
    random_term,
    structure = structure,
    context = "Random-effect summary metadata"
  )
  if(is.null(correlation) || !identical(correlation$type, "lkj")){
    return(out)
  }

  cholesky <- .bt_random_effect_cholesky_draws(
    random_term = random_term,
    n_columns = random_term$n_columns,
    posterior = model_samples
  )
  if(is.null(cholesky)){
    .bt_random_effect_summary_missing_correlation_stop(random_term)
  }

  pairs <- utils::combn(seq_len(random_term$n_columns), 2L)
  values <- matrix(NA_real_, nrow = nrow(model_samples), ncol = ncol(pairs))
  labels <- character(ncol(pairs))
  parts <- vector("list", ncol(pairs))
  for(i in seq_len(ncol(pairs))){
    first <- pairs[1L, i]
    second <- pairs[2L, i]
    first_values <- cholesky[, first, , drop = FALSE]
    second_values <- cholesky[, second, , drop = FALSE]
    dim(first_values) <- c(dim(cholesky)[1L], dim(cholesky)[3L])
    dim(second_values) <- c(dim(cholesky)[1L], dim(cholesky)[3L])
    values[, i] <- rowSums(first_values * second_values)
    pair <- .bt_random_effect_summary_column_components(random_term)[c(first, second)]
    labels[i] <- paste0(pair[1L], ",", pair[2L])
    parts[[i]] <- pair
  }

  list(labels = labels, parts = parts, values = values)
}

.bt_random_effect_summary_column_components <- function(random_term){

  components <- random_term$column_names
  leaves <- random_term$sd_leaves
  if(!is.null(leaves) && !is.null(leaves$leaf_terms_by_column) &&
     length(leaves$leaf_terms_by_column) == random_term$n_columns &&
     !identical(unique(unname(leaves$leaf_terms_by_column)), "sd")){
    components <- leaves$leaf_terms_by_column
  }

  .bt_random_effect_summary_display_components(random_term, components)
}

.bt_random_effect_summary_allocation_samples <- function(allocation,
                                                         random_term = NULL,
                                                         model_samples,
                                                         prior_list,
                                                         include_multipliers = FALSE){

  weights <- .bt_random_effect_dirichlet_draws(
    parameter_name = allocation$weight_name,
    posterior = model_samples,
    prior_list = prior_list
  )
  if(is.null(weights)){
    .bt_random_effect_summary_missing_allocation_stop(allocation)
  }

  components <- .bt_random_effect_summary_allocation_components(
    allocation = allocation,
    K = ncol(weights),
    random_term = random_term
  )
  names <- labels <- types <- character()
  values <- list()

  for(i in seq_len(ncol(weights))){
    names <- c(names, .bt_random_effect_summary_name(
      parameter = sub("__xRE_ALLOCx_.*$", "", allocation$weight_name),
      type = "var_frac",
      parts = c(allocation$label, components[i])
    ))
    labels <- c(labels, paste0("var_frac(", allocation$label, ": ", components[i], ")"))
    types <- c(types, "var_frac")
    values[[length(values) + 1L]] <- weights[, i]
  }

  allocation_target <- .bt_random_effect_summary_allocation_target(allocation)
  if(include_multipliers && identical(allocation_target, "sd")){
    for(i in seq_len(ncol(weights))){
      names <- c(names, .bt_random_effect_summary_name(
        parameter = sub("__xRE_ALLOCx_.*$", "", allocation$weight_name),
        type = "sd_mult",
        parts = c(allocation$label, components[i])
      ))
      labels <- c(labels, paste0("sd_mult(", allocation$label, ": ", components[i], ")"))
      types <- c(types, "sd_multiplier")
      values[[length(values) + 1L]] <- .bt_random_effect_allocation_multiplier(
        weights = weights[, i],
        scale = allocation$scale,
        n_targets = ncol(weights)
      )
    }
  }

  values <- do.call(cbind, values)
  list(
    names = names,
    labels = labels,
    types = types,
    components = rep(components, if(include_multipliers && identical(allocation_target, "sd")) 2L else 1L),
    values = values
  )
}

.bt_random_effect_summary_missing_allocation_stop <- function(allocation){

  label <- allocation$label
  if(!is.character(label) || length(label) != 1L || is.na(label) || !nzchar(label)){
    label <- allocation$weight_name
  }
  if(!is.character(label) || length(label) != 1L || is.na(label) || !nzchar(label)){
    label <- "<unknown>"
  }

  stop(
    "Random-effect allocation summary samples are missing Dirichlet allocation coordinates for allocation '",
    label,
    "'. Expected monitored simplex weights or Dirichlet auxiliary coordinates.",
    call. = FALSE
  )
}

.bt_random_effect_summary_allocation_components <- function(allocation, K,
                                                           random_term = NULL){

  allocation_target <- .bt_random_effect_summary_allocation_target(allocation)
  if(identical(allocation_target, "sd") && !is.null(allocation$leaf_terms)){
    components <- unname(allocation$leaf_terms)
    if(!is.null(random_term)){
      components <- .bt_random_effect_summary_display_components(
        random_term,
        components
      )
    }
  }else{
    terms <- allocation$terms
    components <- names(terms)
    if(is.null(components) || !all(nzchar(components))){
      components <- unname(terms)
    }
  }

  if(length(components) != K){
    components <- paste0("component_", seq_len(K))
  }
  components
}

.bt_random_effect_summary_prior <- function(parameter, type, label,
                                            block = NULL, grouping = NULL,
                                            structure = NULL,
                                            effect_label = NULL,
                                            allocation = NULL,
                                            component = NULL){

  out <- prior_none()
  attr(out, "parameter") <- parameter
  attr(out, "random_summary") <- type
  attr(out, "random_summary_label") <- label
  if(!is.null(block)){
    attr(out, "random_factor") <- block
  }
  if(!is.null(effect_label)){
    attr(out, "random_name") <- effect_label
  }else if(!is.null(block)){
    attr(out, "random_name") <- block
  }
  if(!is.null(grouping)){
    attr(out, "random_grouping_factor") <- grouping
  }
  if(!is.null(structure)){
    attr(out, "random_structure") <- structure
  }
  if(!is.null(allocation)){
    attr(out, "random_allocation") <- allocation
  }
  if(!is.null(component)){
    attr(out, "random_component") <- component
  }

  out
}

.bt_random_effect_summary_group_label <- function(random_term){

  group_label <- random_term$group_label
  if(!is.null(group_label) && length(group_label) == 1L && nzchar(group_label)){
    return(group_label)
  }

  random_term$block_name
}

.bt_random_effect_summary_name <- function(parameter, type, parts){

  paste0(
    parameter,
    "__xRE_SUMMARY__",
    type,
    "__",
    paste(vapply(parts, .bt_random_effect_summary_safe_label, character(1)),
          collapse = "__")
  )
}

.bt_random_effect_summary_safe_label <- function(x){

  x <- gsub("[^A-Za-z0-9_]", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  if(!nzchar(x)){
    x <- "component"
  }
  x
}

.bt_random_effect_summary_unique_name <- function(name, used_names){

  if(!name %in% used_names){
    return(name)
  }
  i <- 2L
  candidate <- paste0(name, "_", i)
  while(candidate %in% used_names){
    i <- i + 1L
    candidate <- paste0(name, "_", i)
  }

  candidate
}

.bt_random_effect_summary_metadata_table <- function(parameter_names,
                                                     prior_list,
                                                     formula_design = NULL){

  n_parameters <- length(parameter_names)
  out <- data.frame(
    "Random name" = rep("", n_parameters),
    "Random grouping" = rep("", n_parameters),
    "Random structure" = rep("", n_parameters),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  if(n_parameters == 0L){
    return(out)
  }

  for(i in seq_along(parameter_names)){
    metadata <- .bt_random_effect_summary_metadata_for_parameter(
      parameter_name = parameter_names[i],
      prior_list = prior_list,
      formula_design = formula_design
    )
    out[i, ] <- metadata
  }

  out
}

.bt_random_effect_summary_add_metadata_columns <- function(table,
                                                          parameter_names,
                                                          prior_list,
                                                          formula_design = NULL){

  metadata <- .bt_random_effect_summary_metadata_table(
    parameter_names = parameter_names,
    prior_list = prior_list,
    formula_design = formula_design
  )
  for(i in rev(seq_len(ncol(metadata)))){
    table <- add_column(
      table = table,
      column_title = colnames(metadata)[i],
      column_values = metadata[[i]],
      column_position = 1,
      column_type = "string"
    )
  }
  attr(table, "random_effects_metadata") <- TRUE

  table
}

.bt_random_effect_summary_metadata_for_parameter <- function(parameter_name,
                                                            prior_list,
                                                            formula_design = NULL){

  prior <- .bt_random_effect_summary_prior_for_column(parameter_name, prior_list)
  if(!is.null(prior) && isTRUE(.bt_random_effect_metadata(prior)$any)){
    effect_label <- .bt_random_effect_prior_name(prior)
    effect <- .bt_random_effect_prior_effect(prior)
    if(!nzchar(effect_label)){
      effect_label <- effect
    }
    return(c(
      "Random name" = effect_label,
      "Random grouping" = .bt_random_effect_prior_grouping(prior),
      "Random structure" = .bt_random_effect_prior_structure(prior)
    ))
  }

  raw_metadata <- .bt_random_effect_summary_raw_metadata_for_parameter(
    parameter_name = parameter_name,
    formula_design = formula_design
  )
  c(
    "Random name" = raw_metadata$name,
    "Random grouping" = raw_metadata$grouping,
    "Random structure" = raw_metadata$structure
  )
}

.bt_random_effect_summary_prior_for_column <- function(parameter_name,
                                                       prior_list){

  if(length(prior_list) == 0L){
    return(NULL)
  }
  if(parameter_name %in% names(prior_list)){
    return(prior_list[[parameter_name]])
  }

  base_name <- sub("\\[[^\\]]+\\]$", "", parameter_name)
  if(base_name %in% names(prior_list)){
    return(prior_list[[base_name]])
  }

  NULL
}

.bt_random_effect_summary_raw_metadata_for_parameter <- function(parameter_name,
                                                                formula_design = NULL){

  out <- list(name = "", grouping = "", structure = "")
  random_design <- .bt_random_effect_summary_designs(formula_design)
  if(length(random_design) == 0L){
    return(out)
  }

  for(design in random_design){
    for(random_term in design$random_effects){
      if(.bt_random_effect_summary_raw_parameter_matches(parameter_name, random_term)){
        return(list(
          name = .bt_random_effect_public_name(random_term),
          grouping = .bt_random_effect_summary_group_label(random_term),
          structure = .bt_random_effect_summary_term_structure(random_term)
        ))
      }
    }
  }

  out
}

.bt_random_effect_summary_raw_parameter_matches <- function(parameter_name,
                                                           random_term){

  stem <- random_term$parameter_stem
  if(is.null(stem) || length(stem) != 1L || !nzchar(stem)){
    return(FALSE)
  }

  startsWith(parameter_name, paste0(stem, "_"))
}

.bt_random_effect_summary_filter_raw_columns <- function(model_samples,
                                                         formula_design = NULL,
                                                         remove_random_effects = NULL,
                                                         keep_random_effects = NULL,
                                                         remove_random_structures = NULL,
                                                         keep_random_structures = NULL){

  column_names <- colnames(model_samples)
  if(length(column_names) == 0L){
    return(model_samples)
  }

  random_design <- .bt_random_effect_summary_designs(formula_design)
  if(length(random_design) == 0L){
    return(model_samples)
  }

  remove_columns <- rep(FALSE, length(column_names))
  for(design in random_design){
    for(random_term in design$random_effects){
      term_columns <- vapply(
        column_names,
        .bt_random_effect_summary_raw_parameter_matches,
        logical(1),
        random_term = random_term
      )
      if(!any(term_columns)){
        next
      }
      if(!is.null(remove_random_effects)){
        term_matches_remove_effect <- .bt_random_effect_summary_term_filter_matches(
          random_term = random_term,
          random_effects = remove_random_effects,
          random_structures = NULL
        )
        remove_columns <- remove_columns | (term_columns & term_matches_remove_effect)
      }
      if(!is.null(remove_random_structures)){
        term_matches_remove_structure <- .bt_random_effect_summary_term_filter_matches(
          random_term = random_term,
          random_effects = NULL,
          random_structures = remove_random_structures
        )
        remove_columns <- remove_columns | (term_columns & term_matches_remove_structure)
      }
      if(!is.null(keep_random_effects) || !is.null(keep_random_structures)){
        term_matches <- .bt_random_effect_summary_term_filter_matches(
          random_term = random_term,
          random_effects = keep_random_effects,
          random_structures = keep_random_structures
        )
        remove_columns <- remove_columns | (term_columns & !term_matches)
      }
    }
  }

  model_samples[, !remove_columns, drop = FALSE]
}

.bt_random_effect_summary_term_filter_matches <- function(random_term,
                                                         random_effects = NULL,
                                                         random_structures = NULL){

  matches <- TRUE
  if(!is.null(random_effects)){
    effect_names <- c(
      random_term$block_name,
      .bt_random_effect_public_name(random_term),
      random_term$group_label
    )
    matches <- matches && any(effect_names %in% random_effects)
  }
  if(!is.null(random_structures)){
    matches <- matches &&
      .bt_random_effect_summary_term_structure(random_term) %in% random_structures
  }

  matches
}

.bt_random_effect_summary_term_structure <- function(random_term){

  .bt_random_effect_structure(
    random_term,
    context = "Random-effect summary metadata"
  )
}

.bt_random_effect_summary_remove_raw <- function(model_samples, prior_list,
                                                 random_design){

  raw_prior_names <- names(prior_list)[vapply(
    prior_list,
    .bt_random_effect_summary_is_raw_prior,
    logical(1)
  )]

  raw_prefixes <- character()
  for(design in random_design){
    for(random_term in design$random_effects){
      raw_prefixes <- c(
        raw_prefixes,
        paste0(random_term$parameter_stem, "_")
      )
    }
  }
  raw_prefixes <- unique(raw_prefixes)

  raw_cols <- .bt_random_effect_summary_parameter_columns(
    colnames(model_samples),
    raw_prior_names
  )
  if(length(raw_prefixes) > 0L){
    raw_cols <- raw_cols | Reduce(
      "|",
      lapply(raw_prefixes, function(prefix) startsWith(colnames(model_samples), prefix))
    )
  }

  if(any(raw_cols)){
    model_samples <- model_samples[, !raw_cols, drop = FALSE]
  }
  if(length(raw_prior_names) > 0L){
    prior_list <- prior_list[!names(prior_list) %in% raw_prior_names]
  }

  list(model_samples = model_samples, prior_list = prior_list)
}

.bt_random_effect_summary_is_raw_prior <- function(prior){

  .bt_is_random_effect_prior(prior, include_summary = FALSE)
}

.bt_random_effect_summary_parameter_columns <- function(column_names,
                                                        parameter_names){

  if(length(parameter_names) == 0L || length(column_names) == 0L){
    return(rep(FALSE, length(column_names)))
  }

  out <- rep(FALSE, length(column_names))
  for(parameter_name in parameter_names){
    eta_name <- .JAGS_prior_dirichlet_eta_name(parameter_name)
    out <- out |
      column_names == parameter_name |
      startsWith(column_names, paste0(parameter_name, "[")) |
      column_names == eta_name |
      startsWith(column_names, paste0(eta_name, "["))
  }

  out
}

.bt_random_effect_summary_display_names <- function(names, raw_names,
                                                    prior_list,
                                                    formula_prefix = TRUE,
                                                    formula_design = NULL){

  if(length(raw_names) == 0L){
    return(names)
  }

  if(length(prior_list) > 0L){
    for(i in seq_along(raw_names)){
      prior <- prior_list[[raw_names[i]]]
      if(is.null(prior)){
        next
      }
      label <- attr(prior, "random_summary_label")
      if(is.null(label)){
        next
      }
      parameter <- attr(prior, "parameter")
      prefix <- .bt_random_effect_summary_formula_prefix(parameter, formula_prefix)
      names[i] <- paste0(prefix, label)
    }
  }

  .bt_random_effect_summary_raw_display_names(
    names = names,
    raw_names = raw_names,
    prior_list = prior_list,
    formula_prefix = formula_prefix,
    formula_design = formula_design
  )
}

.bt_random_effect_summary_formula_prefix <- function(parameter, formula_prefix){

  if(isTRUE(formula_prefix) && !is.null(parameter) &&
     length(parameter) == 1L && nzchar(parameter)){
    return(paste0("(", parameter, ") "))
  }

  ""
}

.bt_random_effect_summary_raw_display_names <- function(names, raw_names,
                                                        prior_list,
                                                        formula_prefix,
                                                        formula_design = NULL){

  random_design <- .bt_random_effect_summary_designs(formula_design)
  if(length(random_design) == 0L){
    return(names)
  }

  for(design in random_design){
    parameter <- design$parameter
    prefix <- .bt_random_effect_summary_formula_prefix(parameter, formula_prefix)
    for(random_term in design$random_effects){
      names <- .bt_random_effect_summary_raw_sd_display_names(
        names = names,
        raw_names = raw_names,
        prior_list = prior_list,
        random_term = random_term,
        prefix = prefix
      )
      names <- .bt_random_effect_summary_raw_rho_display_names(
        names = names,
        raw_names = raw_names,
        random_term = random_term,
        prefix = prefix
      )
      names <- .bt_random_effect_summary_raw_matrix_display_names(
        names = names,
        raw_names = raw_names,
        random_term = random_term,
        prefix = prefix
      )
    }
  }

  names
}

.bt_random_effect_summary_renamed_parameter_names <- function(parameter_names,
                                                              prior_list){

  if(length(parameter_names) == 0L || length(prior_list) == 0L){
    return(parameter_names)
  }

  dummy <- matrix(
    nrow = 0L,
    ncol = length(parameter_names),
    dimnames = list(NULL, parameter_names)
  )
  colnames(.rename_factor_levels(dummy, prior_list))
}

.bt_random_effect_summary_raw_sd_display_names <- function(names, raw_names,
                                                           prior_list,
                                                           random_term,
                                                           prefix){

  sd_names <- unique(random_term$sd_parameter_names)
  sd_names <- sd_names[!is.na(sd_names)]
  if(length(sd_names) == 0L){
    return(names)
  }

  display_sd_names <- .bt_random_effect_summary_renamed_parameter_names(
    parameter_names = sd_names,
    prior_list = prior_list
  )
  components <- .bt_random_effect_summary_sd_components(random_term, sd_names)
  group <- .bt_random_effect_summary_group_label(random_term)
  labels <- paste0(prefix, "sd(", components, " | ", group, ")")

  for(i in seq_along(sd_names)){
    matches <- raw_names %in% c(sd_names[i], display_sd_names[i])
    names[matches] <- labels[i]
  }

  names
}

.bt_random_effect_summary_raw_rho_display_names <- function(names, raw_names,
                                                            random_term,
                                                            prefix){

  stem <- random_term$parameter_stem
  if(is.null(stem) || length(stem) != 1L || !nzchar(stem)){
    return(names)
  }

  group <- .bt_random_effect_summary_group_label(random_term)
  rho_names <- c(
    rho = paste0(stem, "_rho"),
    rho_z = paste0(stem, "_rho_z")
  )
  for(rho_label in names(rho_names)){
    names[raw_names == rho_names[[rho_label]]] <- paste0(
      prefix,
      rho_label,
      "(",
      group,
      ")"
    )
  }

  names
}

.bt_random_effect_summary_raw_matrix_display_names <- function(names,
                                                               raw_names,
                                                               random_term,
                                                               prefix){

  stem <- random_term$parameter_stem
  if(is.null(stem) || length(stem) != 1L || !nzchar(stem)){
    return(names)
  }

  components <- .bt_random_effect_summary_column_components(random_term)
  group <- .bt_random_effect_summary_group_label(random_term)
  group_levels <- random_term$group_levels

  names <- .bt_random_effect_summary_raw_correlation_matrix_names(
    names = names,
    raw_names = raw_names,
    stem = stem,
    matrix = "_xRE_CORx_R",
    label = "cor",
    components = components,
    group = group,
    prefix = prefix
  )
  names <- .bt_random_effect_summary_raw_correlation_matrix_names(
    names = names,
    raw_names = raw_names,
    stem = stem,
    matrix = "_xRE_CORx_L",
    label = "cor_chol",
    components = components,
    group = group,
    prefix = prefix
  )
  names <- .bt_random_effect_summary_raw_effect_matrix_names(
    names = names,
    raw_names = raw_names,
    stem = stem,
    matrix = "_xRE_Zx",
    label = "z",
    components = components,
    group = group,
    group_levels = group_levels,
    prefix = prefix
  )
  names <- .bt_random_effect_summary_raw_effect_matrix_names(
    names = names,
    raw_names = raw_names,
    stem = stem,
    matrix = "_xRE_COEFx",
    label = "coef",
    components = components,
    group = group,
    group_levels = group_levels,
    prefix = prefix
  )

  names
}

.bt_random_effect_summary_raw_correlation_matrix_names <- function(names,
                                                                   raw_names,
                                                                   stem,
                                                                   matrix,
                                                                   label,
                                                                   components,
                                                                   group,
                                                                   prefix){

  matrix_prefix <- paste0(stem, matrix)
  matches <- startsWith(raw_names, paste0(matrix_prefix, "["))
  for(i in which(matches)){
    index <- .bt_random_effect_summary_matrix_index(raw_names[i], matrix_prefix)
    if(anyNA(index) || any(index < 1L) || any(index > length(components))){
      next
    }
    names[i] <- paste0(
      prefix,
      label,
      "(",
      components[index[1L]],
      ",",
      components[index[2L]],
      " | ",
      group,
      ")"
    )
  }

  names
}

.bt_random_effect_summary_raw_effect_matrix_names <- function(names,
                                                              raw_names,
                                                              stem,
                                                              matrix,
                                                              label,
                                                              components,
                                                              group,
                                                              group_levels,
                                                              prefix){

  matrix_prefix <- paste0(stem, matrix)
  matches <- startsWith(raw_names, paste0(matrix_prefix, "["))
  for(i in which(matches)){
    index <- .bt_random_effect_summary_matrix_index(raw_names[i], matrix_prefix)
    if(anyNA(index) || index[2L] < 1L || index[2L] > length(components)){
      next
    }
    group_level <- as.character(index[1L])
    if(!is.null(group_levels) && index[1L] >= 1L &&
       index[1L] <= length(group_levels)){
      group_level <- as.character(group_levels[index[1L]])
    }
    names[i] <- paste0(
      prefix,
      label,
      "(",
      group,
      "[",
      group_level,
      "], ",
      components[index[2L]],
      ")"
    )
  }

  names
}

.bt_random_effect_summary_matrix_index <- function(x, prefix){

  rest <- substring(x, nchar(prefix) + 1L)
  if(!grepl("^\\[[0-9]+,[0-9]+\\]$", rest)){
    return(c(NA_integer_, NA_integer_))
  }

  as.integer(strsplit(substr(rest, 2L, nchar(rest) - 1L), ",", fixed = TRUE)[[1]])
}
