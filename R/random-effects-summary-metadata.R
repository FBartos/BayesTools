.bt_random_effect_summary_rho_samples <- function(random_term, model_samples){

  structure <- .bt_random_effect_summary_term_structure(random_term)
  correlation <- .bt_random_effect_correlation_metadata(
    random_term,
    structure = structure,
    context = "Random-effect summary metadata"
  )
  if(is.null(correlation) || !identical(correlation$type, "rho")){
    return(NULL)
  }

  .bt_random_effect_rho_draws(
    random_term = random_term,
    posterior = model_samples,
    missing = "error",
    out_of_support = "error",
    context = "Random-effect summary metadata"
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
  allocation_type <- .bt_random_effect_summary_allocation_type(allocation)
  names <- labels <- types <- character()
  component_values <- character()
  component_indices <- integer()
  values <- list()

  for(i in seq_len(ncol(weights))){
    names <- c(names, .bt_random_effect_summary_name(
      parameter = sub("__xRE_ALLOCx_.*$", "", allocation$weight_name),
      type = allocation_type$name,
      parts = c(allocation$label, components[i])
    ))
    labels <- c(labels, paste0(allocation_type$label, "(", allocation$label, ": ", components[i], ")"))
    types <- c(types, allocation_type$summary)
    component_values <- c(component_values, components[i])
    component_indices <- c(component_indices, i)
    values[[length(values) + 1L]] <- .bt_random_effect_summary_allocation_values(
      weights = weights[, i],
      allocation = allocation,
      allocation_type = allocation_type,
      K = ncol(weights)
    )
  }

  allocation_target <- .bt_random_effect_summary_allocation_target(allocation)
  if(include_multipliers && identical(allocation_target, "sd_component")){
    allocation_scale <- .bt_random_effect_allocation_scale_metadata(
      allocation,
      context = "Random-effect summary metadata"
    )
    for(i in seq_len(ncol(weights))){
      names <- c(names, .bt_random_effect_summary_name(
        parameter = sub("__xRE_ALLOCx_.*$", "", allocation$weight_name),
        type = "sd_mult",
        parts = c(allocation$label, components[i])
      ))
      labels <- c(labels, paste0("sd_mult(", allocation$label, ": ", components[i], ")"))
      types <- c(types, "sd_multiplier")
      component_values <- c(component_values, components[i])
      component_indices <- c(component_indices, i)
      values[[length(values) + 1L]] <- .bt_random_effect_allocation_multiplier(
        weights = weights[, i],
        scale = allocation_scale,
        n_targets = .bt_random_effect_summary_allocation_n_targets(
          allocation,
          K = ncol(weights)
        )
      )
    }
  }

  values <- do.call(cbind, values)
  list(
    names = names,
    labels = labels,
    types = types,
    components = component_values,
    indices = component_indices,
    values = values
  )
}

.bt_random_effect_summary_allocation_inclusion_samples <- function(allocation,
                                                                   model_samples){

  inclusion <- allocation$inclusion
  if(is.null(inclusion) || length(inclusion) == 0L){
    return(list(
      names = character(),
      labels = character(),
      types = character(),
      components = character(),
      indices = integer(),
      values = matrix(nrow = nrow(model_samples), ncol = 0L)
    ))
  }

  names <- labels <- types <- character()
  component_values <- character()
  component_indices <- integer()
  values <- list()

  for(component_label in names(inclusion)){
    record <- inclusion[[component_label]]
    indicator_name <- record$indicator_name
    if(!is.character(indicator_name) || length(indicator_name) != 1L ||
       is.na(indicator_name) || !nzchar(indicator_name)){
      stop(
        "Random-effect allocation inclusion metadata are missing canonical 'indicator_name'.",
        call. = FALSE
      )
    }
    indicator <- .bt_random_effect_allocation_gate_draws(
      parameter_name = indicator_name,
      posterior = model_samples
    )
    if(is.null(indicator)){
      stop(
        "Random-effect allocation inclusion summary samples are missing Bernoulli indicator '",
        indicator_name,
        "'.",
        call. = FALSE
      )
    }

    names <- c(names, .bt_random_effect_summary_name(
      parameter = sub("__xRE_ALLOCx_.*$", "", allocation$weight_name),
      type = "inclusion",
      parts = c(allocation$label, component_label)
    ))
    labels <- c(labels, paste0("inclusion(", allocation$label, ": ", component_label, ")"))
    types <- c(types, "inclusion")
    component_values <- c(component_values, component_label)
    component_indices <- c(component_indices, record$index)
    values[[length(values) + 1L]] <- indicator
  }

  values <- do.call(cbind, values)
  list(
    names = names,
    labels = labels,
    types = types,
    components = component_values,
    indices = component_indices,
    values = values
  )
}

.bt_random_effect_summary_allocation_type <- function(allocation){

  allocation_target <- .bt_random_effect_summary_allocation_target(allocation)
  allocation_scale <- .bt_random_effect_allocation_scale_metadata(
    allocation,
    context = "Random-effect summary metadata"
  )

  if(identical(allocation_target, "sd_component") &&
     identical(allocation_scale, "mean_variance")){
    return(list(
      name = "var_ratio",
      label = "var_ratio",
      summary = "var_ratio",
      scale = allocation_scale
    ))
  }

  list(
    name = "var_frac",
    label = "var_frac",
    summary = "var_frac",
    scale = allocation_scale
  )
}

.bt_random_effect_summary_allocation_values <- function(weights, allocation,
                                                        allocation_type, K){

  if(identical(allocation_type$summary, "var_ratio")){
    return(.bt_random_effect_allocation_variance_ratio(
      weights = weights,
      n_targets = .bt_random_effect_summary_allocation_n_targets(
        allocation,
        K = K
      )
    ))
  }

  weights
}

.bt_random_effect_summary_allocation_n_targets <- function(allocation, K){

  n_targets <- allocation$n_targets
  if(is.numeric(n_targets) && length(n_targets) == 1L &&
     !is.na(n_targets) && n_targets == as.integer(n_targets) &&
     n_targets >= 2L && n_targets == K){
    return(as.integer(n_targets))
  }

  stop(
    "Random-effect allocation metadata are missing canonical 'allocation$n_targets'.",
    call. = FALSE
  )
}

.bt_random_effect_allocation_variance_ratio <- function(weights, n_targets){

  if(!is.numeric(n_targets) || length(n_targets) != 1L ||
     is.na(n_targets) || n_targets < 1L){
    stop(
      "Random-effect allocation metadata are missing canonical 'allocation$n_targets'.",
      call. = FALSE
    )
  }

  n_targets * weights
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
  if(identical(allocation_target, "sd_component") && !is.null(allocation$leaf_terms)){
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
                                            allocation_metadata = NULL,
                                            allocation_index = NULL,
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
  if(!is.null(allocation_metadata)){
    attr(out, "random_allocation_metadata") <- allocation_metadata
  }
  if(!is.null(allocation_index)){
    attr(out, "random_allocation_index") <- allocation_index
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
