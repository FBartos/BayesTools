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
  if(!inherits(leaves, "BayesTools_random_effect_sd_leaves") ||
     is.null(leaves$leaf_terms)){
    stop(
      "Random-effect summary metadata for block '",
      random_term$block_name,
      "' are missing canonical 'random_term$sd_leaves'.",
      call. = FALSE
    )
  }
  out <- unname(leaves$leaf_terms[sd_names])
  if(anyNA(out)){
    stop(
      "Random-effect summary metadata for block '",
      random_term$block_name,
      "' do not map every SD coordinate to a semantic component.",
      call. = FALSE
    )
  }
  if(.bt_random_effect_summary_term_structure(random_term) %in%
     c("id", "cs", "ar1", "car") && length(out) == 1L){
    out <- "shared"
  }
  .bt_random_effect_summary_display_components(random_term, out)
}

.bt_random_effect_summary_normalize_components <- function(components){

  components <- gsub("__xXx__", ":", components, fixed = TRUE)
  components[components == "(Intercept)"] <- "intercept"
  components
}

.bt_random_effect_semantic_name <- function(parameter, owner, quantity,
                                            arguments = character(),
                                            formula_prefix = TRUE){

  quantity_label <- .bt_random_effect_semantic_quantity_name(
    quantity,
    arguments
  )
  prefix <- .bt_random_effect_summary_formula_prefix(
    parameter,
    formula_prefix
  )
  owner_prefix <- if(nzchar(owner)) paste0(owner, ": ") else ""
  paste0(prefix, owner_prefix, quantity_label)
}

.bt_random_effect_allocation_public_name <- function(allocation){

  display_name <- allocation$display_name
  if(is.character(display_name) && length(display_name) == 1L &&
     !is.na(display_name)){
    return(display_name)
  }

  allocation$label
}

.bt_random_effect_semantic_quantity_name <- function(quantity,
                                                     arguments = character()){

  arguments <- as.character(arguments)
  arguments <- arguments[!is.na(arguments) & nzchar(arguments)]
  if(length(arguments) == 0L){
    return(quantity)
  }

  paste0(quantity, "(", paste(arguments, collapse = ","), ")")
}

.bt_random_effect_semantic_sd_arguments <- function(component){

  if(identical(component, "shared")) character() else component
}

.bt_random_effect_semantic_sd_display_arguments <- function(random_term,
                                                            component){

  arguments <- .bt_random_effect_semantic_sd_arguments(component)
  sd_names <- unique(random_term$sd_parameter_names)
  sd_names <- sd_names[!is.na(sd_names)]
  if(length(sd_names) == 1L && identical(arguments, "intercept")){
    return(character())
  }

  arguments
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

  index <- random_term$structured_index
  if(!is.null(index) && !is.null(index$name) && !is.null(index$label) &&
     length(index$name) == 1L && length(index$label) == 1L &&
     !is.null(random_term$xlevels[[index$name]]) &&
     length(random_term$xlevels[[index$name]]) == random_term$n_columns){
    return(paste0(
      index$label,
      "[",
      as.character(random_term$xlevels[[index$name]]),
      "]"
    ))
  }

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

  if(isTRUE(allocation$gate_only)){
    stop(
      "Gate-only random-effect allocations do not define variance weights.",
      call. = FALSE
    )
  }
  weights <- .bt_random_effect_dirichlet_draws(
    parameter_name = allocation$weight_name,
    posterior = model_samples,
    prior_list = prior_list
  )
  if(is.null(weights)){
    .bt_random_effect_summary_missing_allocation_stop(allocation)
  }

  realized <- .bt_random_effect_summary_realized_allocation(
    allocation = allocation,
    weights = weights,
    model_samples = model_samples
  )

  components <- .bt_random_effect_summary_allocation_components(
    allocation = allocation,
    K = ncol(weights),
    random_term = random_term
  )
  allocation_type <- .bt_random_effect_summary_allocation_type(allocation)
  allocation_owner <- .bt_random_effect_allocation_public_name(allocation)
  names <- labels <- types <- character()
  component_values <- character()
  component_indices <- integer()
  values <- list()

  for(i in seq_len(ncol(weights))){
    names <- c(names, .bt_random_effect_summary_name(
      parameter = .bt_random_effect_allocation_formula_parameter(allocation),
      type = allocation_type$name,
      parts = c(allocation$label, components[i])
    ))
    labels <- c(labels, .bt_random_effect_semantic_name(
      parameter = "",
      owner = allocation_owner,
      quantity = allocation_type$label,
      arguments = components[i],
      formula_prefix = FALSE
    ))
    types <- c(types, allocation_type$summary)
    component_values <- c(component_values, components[i])
    component_indices <- c(component_indices, i)
    values[[length(values) + 1L]] <- .bt_random_effect_summary_allocation_values(
      weights = if(identical(allocation_type$summary, "var_prop")){
        realized$proportions[, i]
      }else{
        weights[, i]
      },
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
        parameter = .bt_random_effect_allocation_formula_parameter(allocation),
        type = "sd_mult",
        parts = c(allocation$label, components[i])
      ))
      labels <- c(labels, .bt_random_effect_semantic_name(
        parameter = "",
        owner = allocation_owner,
        quantity = "sd_mult",
        arguments = components[i],
        formula_prefix = FALSE
      ))
      types <- c(types, "sd_mult")
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

.bt_random_effect_allocation_scale_role <- function(allocation){

  scale <- .bt_random_effect_allocation_scale_metadata(
    allocation,
    context = "Random-effect allocation semantic metadata"
  )
  if(identical(scale, "total_variance")) "total" else "common"
}

.bt_random_effect_allocation_sd_quantity <- function(allocation){

  paste0("sd_", .bt_random_effect_allocation_scale_role(allocation))
}

.bt_random_effect_allocation_var_quantity <- function(allocation){

  paste0("var_", .bt_random_effect_allocation_scale_role(allocation))
}

.bt_random_effect_summary_allocation_scale_samples <- function(
    allocation, model_samples, prior_list){

  source <- allocation$source
  if(!is.list(source) || !identical(source$shape, "scalar")){
    return(NULL)
  }
  source_name <- .bt_random_sd_binding_source_name(source)
  values <- .bt_random_effect_parameter_draws(
    parameter_name = source_name,
    posterior = model_samples,
    prior_list = prior_list
  )
  if(is.null(values)){
    return(NULL)
  }
  factors <- allocation$parent_factors
  if(is.null(factors)){
    factors <- list()
  }
  values <- .bt_random_effect_apply_allocation_factors(
    base = values,
    factors = factors,
    posterior = model_samples,
    prior_list = prior_list
  )
  if(is.null(values) ||
     !identical(.bt_random_effect_allocation_scale_metadata(
       allocation,
       context = "Random-effect allocation summary metadata"
     ), "total_variance") ||
     length(.bt_random_effect_summary_allocation_gate_names(
       allocation,
       include_parents = FALSE
     )) == 0L){
    return(values)
  }

  if(isTRUE(allocation$gate_only)){
    gate_fraction <- rep(1, nrow(model_samples))
    for(record in allocation$inclusion){
      gate <- .bt_random_effect_allocation_gate_draws(
        parameter_name = record$indicator_name,
        posterior = model_samples
      )
      if(is.null(gate)){
        stop(
          "Random-effect allocation inclusion samples are missing Bernoulli indicator '",
          record$indicator_name,
          "'.",
          call. = FALSE
        )
      }
      gate_fraction <- gate_fraction * gate
    }
    return(values * sqrt(gate_fraction))
  }

  weights <- .bt_random_effect_dirichlet_draws(
    parameter_name = allocation$weight_name,
    posterior = model_samples,
    prior_list = prior_list
  )
  if(is.null(weights)){
    .bt_random_effect_summary_missing_allocation_stop(allocation)
  }
  realized <- .bt_random_effect_summary_realized_allocation(
    allocation = allocation,
    weights = weights,
    model_samples = model_samples
  )

  values * sqrt(realized$total_fraction)
}

.bt_random_effect_summary_realized_allocation <- function(
    allocation, weights, model_samples){

  if(!is.matrix(weights) || !is.numeric(weights) ||
     nrow(weights) != nrow(model_samples)){
    stop(
      "Random-effect allocation weights do not align with posterior draws.",
      call. = FALSE
    )
  }
  K <- .bt_random_effect_summary_allocation_n_targets(
    allocation,
    K = ncol(weights)
  )
  scale <- .bt_random_effect_allocation_scale_metadata(
    allocation,
    context = "Random-effect allocation summary metadata"
  )
  if(!identical(scale, "total_variance")){
    return(list(
      total_fraction = rep(1, nrow(weights)),
      proportions = weights,
      defined = rep(TRUE, nrow(weights))
    ))
  }
  if(length(.bt_random_effect_summary_allocation_gate_names(
    allocation
  )) == 0L){
    return(list(
      total_fraction = rep(1, nrow(weights)),
      proportions = weights,
      defined = rep(TRUE, nrow(weights))
    ))
  }

  component_gates <- matrix(1, nrow = nrow(weights), ncol = K)
  inclusion <- allocation$inclusion
  if(is.null(inclusion)){
    inclusion <- list()
  }
  for(record in inclusion){
    index <- record$index
    if(!is.numeric(index) || length(index) != 1L || is.na(index) ||
       index != as.integer(index) || index < 1L || index > K){
      stop(
        "Random-effect allocation inclusion metadata reference an invalid component index.",
        call. = FALSE
      )
    }
    gate <- .bt_random_effect_allocation_gate_draws(
      parameter_name = record$indicator_name,
      posterior = model_samples
    )
    if(is.null(gate)){
      stop(
        "Random-effect allocation inclusion samples are missing Bernoulli indicator '",
        record$indicator_name,
        "'.",
        call. = FALSE
      )
    }
    component_gates[, as.integer(index)] <- gate
  }

  parent_active <- rep(TRUE, nrow(weights))
  parent_factors <- allocation$parent_factors
  if(is.null(parent_factors)){
    parent_factors <- list()
  }
  for(factor in parent_factors){
    if(is.null(factor$inclusion_name)){
      next
    }
    gate <- .bt_random_effect_allocation_gate_draws(
      parameter_name = factor$inclusion_name,
      posterior = model_samples
    )
    if(is.null(gate)){
      stop(
        "Random-effect allocation inclusion samples are missing Bernoulli indicator '",
        factor$inclusion_name,
        "'.",
        call. = FALSE
      )
    }
    parent_active <- parent_active & gate == 1
  }

  gated_weights <- weights * component_gates
  total_fraction <- rowSums(gated_weights) * as.numeric(parent_active)
  defined <- parent_active & total_fraction > 0
  proportions <- matrix(
    NA_real_,
    nrow = nrow(weights),
    ncol = K,
    dimnames = dimnames(weights)
  )
  if(any(defined)){
    proportions[defined, ] <-
      gated_weights[defined, , drop = FALSE] / total_fraction[defined]
  }

  list(
    total_fraction = total_fraction,
    proportions = proportions,
    defined = defined
  )
}

.bt_random_effect_summary_allocation_gate_names <- function(
    allocation, include_components = TRUE, include_parents = TRUE){

  names <- character()
  if(include_components && is.list(allocation$inclusion)){
    names <- c(names, unlist(lapply(allocation$inclusion, function(record){
      record$indicator_name
    }), use.names = FALSE))
  }
  if(include_parents && is.list(allocation$parent_factors)){
    names <- c(names, unlist(lapply(allocation$parent_factors, function(factor){
      factor$inclusion_name
    }), use.names = FALSE))
  }

  unique(names[!is.na(names) & nzchar(names)])
}

.bt_random_effect_allocation_formula_parameter <- function(allocation){

  parameter <- allocation$parameter
  if(is.character(parameter) && length(parameter) == 1L &&
     !is.na(parameter) && nzchar(parameter)){
    return(parameter)
  }
  weight_name <- allocation$weight_name
  if(is.character(weight_name) && length(weight_name) == 1L &&
     !is.na(weight_name) && nzchar(weight_name)){
    return(sub("__xRE_ALLOCx_.*$", "", weight_name))
  }

  stop(
    "Random-effect allocation metadata are missing canonical 'parameter'.",
    call. = FALSE
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
      name = "var_mult",
      label = "var_mult",
      summary = "var_mult",
      scale = allocation_scale
    ))
  }

  list(
    name = "var_prop",
    label = "var_prop",
    summary = "var_prop",
    scale = allocation_scale
  )
}

.bt_random_effect_summary_allocation_values <- function(weights, allocation,
                                                        allocation_type, K){

  if(identical(allocation_type$summary, "var_mult")){
    return(.bt_random_effect_allocation_variance_multiplier(
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

.bt_random_effect_allocation_variance_multiplier <- function(weights, n_targets){

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
    public_components <- allocation$component_names
    if(is.character(public_components) && length(public_components) == K &&
       !anyNA(public_components) && all(nzchar(public_components))){
      components <- unname(public_components)
    }
  }

  if(length(components) != K){
    components <- paste0("component_", seq_len(K))
  }
  components
}

.bt_random_effect_allocation_component_name <- function(allocation,
                                                        component_label){

  terms <- allocation$terms
  internal <- names(terms)
  public <- allocation$component_names
  if(is.character(internal) && is.character(public) &&
     length(internal) == length(public)){
    index <- match(component_label, internal)
    if(!is.na(index)){
      return(unname(public[index]))
    }
  }

  component_label
}

.bt_random_effect_summary_prior <- function(parameter, type, label,
                                             component_label = NULL,
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
  if(!is.null(component_label)){
    attr(out, "random_summary_component_label") <- component_label
  }
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
