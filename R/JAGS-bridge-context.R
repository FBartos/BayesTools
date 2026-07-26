.bt_JAGS_bridge_cache_posterior_row <- function(samples, cache){

  if(isTRUE(cache)){
    posterior_row <- matrix(
      unname(samples),
      nrow = 1L,
      dimnames = list(NULL, names(samples))
    )
    attr(
      posterior_row,
      "BayesTools_random_effect_dirichlet_draw_cache"
    ) <- new.env(parent = emptyenv())
    attr(samples, "BayesTools_marglik_posterior_row") <- posterior_row
  }

  samples
}

.bt_JAGS_bridge_call_log_posterior <- function(log_posterior, parameters,
                                               data, context,
                                               bridge_context, ...){

  if(isTRUE(bridge_context)){
    return(log_posterior(
      parameters = parameters,
      data = data,
      bridge_context = context,
      ...
    ))
  }

  log_posterior(parameters = parameters, data = data, ...)
}

.bt_JAGS_bridge_context <- function(samples, prior_parameters,
                                    formula_prior_parameters,
                                    formula_parameters,
                                    add_parameters,
                                    formula_design_list,
                                    formula_data_list,
                                    formula_prior_list,
                                    model_data){

  state <- .bt_JAGS_bridge_context_state(samples)
  state_matrix <- .bt_JAGS_bridge_context_state_matrix(samples)
  add_parameter_values <- .bt_JAGS_bridge_context_add_parameters(
    samples = samples,
    add_parameters = add_parameters
  )
  random <- .bt_JAGS_bridge_context_random(
    samples = samples,
    prior_parameters = prior_parameters,
    formula_prior_parameters = formula_prior_parameters,
    formula_parameters = formula_parameters,
    formula_design_list = formula_design_list,
    formula_data_list = formula_data_list,
    formula_prior_list = formula_prior_list,
    model_data = model_data
  )
  random_nodes <- .bt_JAGS_bridge_context_random_nodes(random)

  nodes <- .bt_JAGS_bridge_context_nodes(
    state = state,
    prior_parameters = prior_parameters,
    formula_prior_parameters = formula_prior_parameters,
    formula_parameters = formula_parameters,
    add_parameter_values = add_parameter_values,
    random_nodes = random_nodes
  )
  node_info <- .bt_JAGS_bridge_context_node_info(
    state = state,
    prior_parameters = prior_parameters,
    formula_prior_parameters = formula_prior_parameters,
    formula_parameters = formula_parameters,
    add_parameter_values = add_parameter_values,
    random = random
  )

  out <- list(
    state = state,
    state_matrix = state_matrix,
    nodes = nodes,
    prior_parameters = prior_parameters,
    formula_prior_parameters = formula_prior_parameters,
    formula_parameters = formula_parameters,
    add_parameters = add_parameter_values,
    random = random,
    node_info = node_info,
    metadata = list(
      source = "JAGS_bridgesampling",
      state_names = names(state),
      node_names = names(nodes),
      node_info_names = node_info$name,
      formula_parameters = names(formula_parameters),
      random_parameters = names(random)
    )
  )
  class(out) <- c("BayesTools_bridge_context", "list")
  out
}

.bt_JAGS_bridge_context_state <- function(samples){

  state <- as.numeric(samples)
  names(state) <- names(samples)
  state
}

.bt_JAGS_bridge_context_state_matrix <- function(samples){

  .bt_JAGS_marglik_random_effect_posterior_row(samples)
}

.bt_JAGS_bridge_context_add_parameters <- function(samples, add_parameters){

  if(is.null(add_parameters) || length(add_parameters) == 0L){
    return(list())
  }

  as.list(samples[add_parameters])
}

.bt_JAGS_bridge_context_nodes <- function(state, prior_parameters,
                                          formula_prior_parameters,
                                          formula_parameters,
                                          add_parameter_values,
                                          random_nodes){

  .bt_JAGS_bridge_merge_nodes(
    state,
    .bt_JAGS_bridge_flatten_parameter_nodes(prior_parameters),
    .bt_JAGS_bridge_flatten_parameter_nodes(formula_prior_parameters),
    .bt_JAGS_bridge_flatten_parameter_nodes(formula_parameters),
    .bt_JAGS_bridge_flatten_parameter_nodes(add_parameter_values),
    random_nodes
  )
}

.bt_JAGS_bridge_context_node_info <- function(state, prior_parameters,
                                              formula_prior_parameters,
                                              formula_parameters,
                                              add_parameter_values,
                                              random){

  .bt_JAGS_bridge_merge_node_info(
    .bt_JAGS_bridge_node_info(
      nodes = state,
      source = "state",
      owner = "bridge",
      role = "coordinate"
    ),
    .bt_JAGS_bridge_node_info(
      nodes = .bt_JAGS_bridge_flatten_parameter_nodes(prior_parameters),
      source = "prior_parameters",
      owner = "model",
      role = "prior_parameter"
    ),
    .bt_JAGS_bridge_node_info(
      nodes = .bt_JAGS_bridge_flatten_parameter_nodes(formula_prior_parameters),
      source = "formula_prior_parameters",
      owner = "formula",
      role = "formula_prior_parameter"
    ),
    .bt_JAGS_bridge_node_info(
      nodes = .bt_JAGS_bridge_flatten_parameter_nodes(formula_parameters),
      source = "formula_parameters",
      owner = "formula",
      role = "formula_value"
    ),
    .bt_JAGS_bridge_node_info(
      nodes = .bt_JAGS_bridge_flatten_parameter_nodes(add_parameter_values),
      source = "add_parameters",
      owner = "model",
      role = "additional_parameter"
    ),
    .bt_JAGS_bridge_context_random_node_info(random)
  )
}

.bt_JAGS_bridge_node_info <- function(nodes, source, owner, role,
                                      parameter = NA_character_,
                                      block_name = NA_character_){

  if(length(nodes) == 0L || is.null(names(nodes))){
    return(.bt_JAGS_bridge_empty_node_info())
  }

  node_names <- names(nodes)
  node_names <- node_names[!is.na(node_names) & nzchar(node_names)]
  if(length(node_names) == 0L){
    return(.bt_JAGS_bridge_empty_node_info())
  }

  data.frame(
    name = node_names,
    source = rep(source, length(node_names)),
    owner = rep(owner, length(node_names)),
    role = rep(role, length(node_names)),
    parameter = rep(parameter, length(node_names)),
    block_name = rep(block_name, length(node_names)),
    stringsAsFactors = FALSE
  )
}

.bt_JAGS_bridge_empty_node_info <- function(){

  data.frame(
    name = character(),
    source = character(),
    owner = character(),
    role = character(),
    parameter = character(),
    block_name = character(),
    stringsAsFactors = FALSE
  )
}

.bt_JAGS_bridge_merge_node_info <- function(...){

  pieces <- list(...)
  pieces <- pieces[vapply(pieces, nrow, integer(1L)) > 0L]
  if(length(pieces) == 0L){
    return(.bt_JAGS_bridge_empty_node_info())
  }

  out <- do.call(rbind, pieces)
  out <- out[!duplicated(out$name, fromLast = TRUE), , drop = FALSE]
  rownames(out) <- NULL
  out
}

.bt_JAGS_bridge_context_random_node_info <- function(random){

  out <- .bt_JAGS_bridge_empty_node_info()
  random_node_names <- character()
  allocation_names <- character()
  for(parameter in names(random)){
    for(block_name in names(random[[parameter]])){
      block <- random[[parameter]][[block_name]]
      block_nodes <- block$nodes
      duplicated_block_nodes <- names(block_nodes)[names(block_nodes) %in% random_node_names]
      if(length(duplicated_block_nodes) > 0L){
        out$block_name[out$name %in% duplicated_block_nodes] <- NA_character_
      }
      block_nodes <- block_nodes[!names(block_nodes) %in% random_node_names]
      random_node_names <- c(random_node_names, names(block_nodes))
      out <- .bt_JAGS_bridge_merge_node_info(
        out,
        .bt_JAGS_bridge_node_info(
          nodes = block_nodes,
          source = "random",
          owner = "formula",
          role = "random_effect",
          parameter = parameter,
          block_name = block_name
        )
      )
      allocations <- block$allocation$weights
      if(length(allocations) > 0L){
        for(allocation in allocations){
          allocation <- allocation[!names(allocation) %in% allocation_names]
          allocation_names <- c(allocation_names, names(allocation))
          out <- .bt_JAGS_bridge_merge_node_info(
            out,
            .bt_JAGS_bridge_node_info(
              nodes = allocation,
              source = "random",
              owner = "formula",
              role = "random_allocation",
              parameter = parameter,
              block_name = NA_character_
            )
          )
        }
      }
    }
  }

  out
}

.bt_JAGS_bridge_flatten_parameter_nodes <- function(parameters){

  if(length(parameters) == 0L){
    return(numeric())
  }

  out <- numeric()
  parameter_names <- names(parameters)
  if(is.null(parameter_names)){
    return(out)
  }

  for(parameter in parameter_names){
    if(!nzchar(parameter)){
      next
    }
    value <- parameters[[parameter]]
    out <- .bt_JAGS_bridge_merge_nodes(
      out,
      .bt_JAGS_bridge_flatten_parameter_value(
        parameter = parameter,
        value = value
      )
    )
  }

  out
}

.bt_JAGS_bridge_flatten_parameter_value <- function(parameter, value){

  if(is.null(value) || length(value) == 0L){
    return(numeric())
  }
  if(!is.numeric(value) && !is.logical(value)){
    return(numeric())
  }

  values <- as.numeric(value)
  dims <- dim(value)
  if(is.null(dims)){
    node_names <- if(length(values) == 1L){
      parameter
    }else{
      paste0(parameter, "[", seq_along(values), "]")
    }
  }else{
    indices <- do.call(expand.grid, lapply(dims, seq_len))
    node_names <- paste0(
      parameter,
      "[",
      apply(indices, 1L, paste, collapse = ","),
      "]"
    )
  }

  names(values) <- node_names
  values
}

.bt_JAGS_bridge_merge_nodes <- function(...){

  pieces <- list(...)
  out <- numeric()
  for(piece in pieces){
    if(length(piece) == 0L){
      next
    }
    if(is.null(names(piece))){
      next
    }
    out[names(piece)] <- as.numeric(piece)
  }

  out
}

.bt_JAGS_bridge_context_random <- function(samples, prior_parameters,
                                           formula_prior_parameters,
                                           formula_parameters,
                                           formula_design_list,
                                           formula_data_list,
                                           formula_prior_list,
                                           model_data){

  if(length(formula_design_list) == 0L){
    return(list())
  }

  out <- list()
  parameter_names <- names(formula_design_list)
  if(is.null(parameter_names)){
    parameter_names <- rep("", length(formula_design_list))
  }
  source_parameters <- .bt_JAGS_bridge_context_source_parameters(
    samples = samples,
    prior_parameters = prior_parameters,
    formula_prior_parameters = formula_prior_parameters,
    formula_parameters = formula_parameters
  )
  sampled_random_parameters <- vapply(
    seq_along(formula_design_list),
    function(parameter_i){
      design <- formula_design_list[[parameter_i]]
      if(!.bt_formula_design_has_sampled_random_effects(design)){
        return(NA_character_)
      }
      .bt_JAGS_bridge_design_parameter_name(
        design,
        fallback = parameter_names[[parameter_i]]
      )
    },
    character(1)
  )
  sampled_random_parameters <- unique(stats::na.omit(
    sampled_random_parameters
  ))
  source_parameters <- .bt_parameter_source_forbid_formula_parameters(
    source_parameters,
    sampled_random_parameters
  )

  for(parameter_i in seq_along(formula_design_list)){
    design <- formula_design_list[[parameter_i]]
    if(!.bt_formula_design_has_any_random_effects(design)){
      next
    }
    parameter <- .bt_JAGS_bridge_design_parameter_name(
      design,
      fallback = parameter_names[[parameter_i]]
    )
    if(.bt_formula_design_has_sampled_random_effects(design)){
      design <- .bt_JAGS_bridge_prepare_random_effect_allocation_design(design)
    }
    formula_data <- if(!is.null(formula_data_list)){
      formula_data_list[[parameter]]
    }else{
      NULL
    }
    parameter_prior_list <- if(!is.null(formula_prior_list)){
      formula_prior_list[[parameter]]
    }else{
      list()
    }
    if(is.null(parameter_prior_list)){
      parameter_prior_list <- list()
    }
    source_data <- .bt_JAGS_marglik_parameter_source_data(
      model_data = model_data,
      formula_data = formula_data,
      design = design
    )

    parameter_random <- list()
    for(random_term in .bt_formula_design_random_effects(design)){
      parameter_random[[random_term$block_name]] <-
        .bt_JAGS_bridge_context_random_block(
          samples = samples,
          random_term = random_term,
          prior_list = parameter_prior_list,
          formula_prior_parameters = formula_prior_parameters,
          data = source_data,
          parameters = source_parameters
        )
    }
    out[[parameter]] <- parameter_random
  }

  out
}

.bt_JAGS_bridge_context_source_parameters <- function(samples,
                                                      prior_parameters,
                                                      formula_prior_parameters,
                                                      formula_parameters){

  out <- as.list(samples)
  out[names(prior_parameters)] <- prior_parameters
  out[names(formula_prior_parameters)] <- formula_prior_parameters
  out[names(formula_parameters)] <- formula_parameters
  out
}

.bt_JAGS_bridge_context_random_block <- function(samples, random_term,
                                                prior_list,
                                                formula_prior_parameters,
                                                data,
                                                parameters){

  posterior <- .bt_JAGS_marglik_random_effect_posterior_row(samples)
  n_columns <- random_term$n_columns
  structure <- .bt_JAGS_bridge_random_term_structure(random_term)
  rho <- if(n_columns > 1L &&
            structure %in% c("cs", "hcs", "ar1", "car", "har")){
    .bt_JAGS_bridge_context_random_rho(
      random_term = random_term,
      posterior = posterior
    )
  }else{
    NULL
  }
  layout <- random_term$latent_layout
  group_local <- inherits(
    layout,
    "BayesTools_random_effect_structured_local_layout"
  )
  scalar_structure <- n_columns > 1L &&
    structure %in% c("cs", "hcs", "ar1", "car", "har") &&
    !is.null(rho)
  independent <- n_columns > 1L && structure %in% c("diag", "id")
  compact_correlation <- isTRUE(group_local) || isTRUE(scalar_structure) ||
    isTRUE(independent)
  if(isTRUE(compact_correlation)){
    cholesky <- NULL
    correlation <- NULL
    correlation_blocks <- NULL
  }else{
    cholesky <- .bt_JAGS_marglik_random_effect_cholesky(
      samples = samples,
      random_term = random_term
    )
    correlation <- tcrossprod(cholesky)
    correlation_blocks <- NULL
  }
  allocation <- .bt_JAGS_bridge_context_random_allocation(
    random_term = random_term,
    posterior = posterior,
    prior_list = prior_list,
    formula_prior_parameters = formula_prior_parameters
  )
  latent <- .bt_JAGS_bridge_context_random_latent(
    samples = samples,
    random_term = random_term
  )
  group_covariance <- .bt_random_effect_group_covariance_metadata(random_term)
  nodes <- numeric()
  column_names <- .bt_JAGS_bridge_context_random_column_names(random_term)

  out <- list(
    block_name = random_term$block_name,
    parameter_stem = random_term$parameter_stem,
    compile_mode = .bt_random_effect_term_compile_mode(random_term),
    structure = structure,
    dimensions = list(
      n_groups = random_term$n_groups,
      n_columns = n_columns
    ),
    levels = list(
      group = random_term$group_levels,
      column = column_names
    ),
    scale = list(
      type = NULL,
      column_sd = NULL,
      row_sd_source = NULL
    ),
    allocation = allocation,
    correlation = list(
      rho = rho,
      cholesky = cholesky,
      matrix = correlation,
      blocks = correlation_blocks,
      block_columns = NULL,
      structure = if(isTRUE(group_local)) layout$structure else if(
        isTRUE(scalar_structure)
      ) {
        structure
      } else if(isTRUE(independent)) {
        structure
      } else {
        NULL
      },
      column_coordinates = if(isTRUE(group_local)) {
        layout$column_coordinates
      } else if(isTRUE(scalar_structure) && identical(structure, "car")) {
        random_term$car$time_values
      } else if(isTRUE(scalar_structure)) {
        seq_len(n_columns)
      } else if(isTRUE(independent)) {
        seq_len(n_columns)
      } else {
        NULL
      }
    ),
    group_covariance = group_covariance,
    covariance = NULL,
    latent = latent,
    nodes = nodes
  )

  if(.bt_random_effect_has_row_indexed_external_sd(random_term)){
    source_draws <- .bt_JAGS_marglik_row_indexed_external_sd_source_draws(
      random_term = random_term,
      n_rows = nrow(random_term$model_matrix),
      posterior = posterior,
      data = data,
      parameters = parameters,
      context = "Bridge context"
    )
    out$scale$type <- "row_indexed"
    out$scale$row_sd_source <- unname(source_draws[1L, ])
    names(out$scale$row_sd_source) <- colnames(source_draws)
    out$nodes <- .bt_JAGS_bridge_merge_nodes(
      out$nodes,
      out$scale$row_sd_source
    )

    column_allocation <- .bt_JAGS_marglik_row_indexed_column_allocation_draws(
      random_term = random_term,
      posterior = posterior,
      prior_list = prior_list,
      n_columns = n_columns
    )
    if(!is.null(column_allocation)){
      out$allocation$column_multiplier <- unname(column_allocation[1L, ])
    }else{
      allocation_draws <- .bt_JAGS_marglik_row_indexed_allocation_draws(
        random_term = random_term,
        posterior = posterior,
        prior_list = prior_list
      )
      out$allocation$block_multiplier <- unname(allocation_draws[1L])
    }
    return(out)
  }

  sd_values <- .bt_JAGS_marglik_random_effect_sd_values(
    samples = samples,
    random_term = random_term,
    prior_list = prior_list
  )
  names(sd_values) <- column_names
  out$scale$type <- "column"
  out$scale$column_sd <- sd_values
  if(isTRUE(compact_correlation)){
    out$covariance <- NULL
  }else{
    out$covariance <- correlation * tcrossprod(sd_values)
  }
  sd_nodes <- .bt_JAGS_bridge_context_random_sd_nodes(
    random_term = random_term,
    sd_values = sd_values
  )
  out$nodes <- .bt_JAGS_bridge_merge_nodes(out$nodes, sd_nodes)

  out
}

.bt_JAGS_bridge_context_random_rho <- function(random_term, posterior){

  rho <- .bt_random_effect_rho_draws(
    random_term = random_term,
    posterior = posterior,
    missing = "null",
    out_of_support = "error",
    context = "Bridge context random-effect metadata"
  )
  if(is.null(rho)){
    return(NULL)
  }

  unname(rho[1L])
}

.bt_JAGS_bridge_context_random_allocation <- function(random_term, posterior,
                                                      prior_list,
                                                      formula_prior_parameters){

  parameter_names <- .bt_JAGS_marglik_random_effect_allocation_parameters(
    random_term
  )
  weights <- list()
  for(parameter_name in parameter_names){
    if(parameter_name %in% names(formula_prior_parameters)){
      value <- formula_prior_parameters[[parameter_name]]
      weights[[parameter_name]] <- as.numeric(value)
      names(weights[[parameter_name]]) <- paste0(
        parameter_name,
        "[",
        seq_along(weights[[parameter_name]]),
        "]"
      )
      next
    }

    value <- .bt_random_effect_dirichlet_draws(
      parameter_name = parameter_name,
      posterior = posterior,
      prior_list = prior_list
    )
    if(!is.null(value)){
      weights[[parameter_name]] <- unname(value[1L, ])
      names(weights[[parameter_name]]) <- paste0(
        parameter_name,
        "[",
        seq_along(weights[[parameter_name]]),
        "]"
      )
    }
  }

  list(weights = weights)
}

.bt_JAGS_bridge_context_random_latent <- function(samples, random_term){

  if(!identical(.bt_random_effect_term_compile_mode(random_term), "sampled")){
    return(NULL)
  }

  z_names <- .bt_random_effect_latent_names(
    random_term = random_term,
    n_groups = random_term$n_groups,
    n_columns = random_term$n_columns
  )
  if(!all(as.vector(z_names) %in% names(samples))){
    return(NULL)
  }

  layout <- random_term$latent_layout
  if(inherits(layout, "BayesTools_random_effect_structured_local_layout")){
    z <- list(
      values = unname(samples[as.vector(z_names)]),
      node_names = as.vector(z_names),
      group = layout$local_group,
      column = layout$local_column,
      dimensions = c(
        n_groups = random_term$n_groups,
        n_columns = random_term$n_columns
      )
    )
    class(z) <- c("BayesTools_group_local_latent", "list")
    return(z)
  }else{
    z <- matrix(
      unname(samples[as.vector(z_names)]),
      nrow = random_term$n_groups,
      ncol = random_term$n_columns
    )
  }
  dimnames(z) <- list(
    group = random_term$group_levels,
    column = .bt_JAGS_bridge_context_random_column_names(random_term)
  )
  z
}

.bt_JAGS_bridge_context_random_column_names <- function(random_term){

  column_names <- colnames(random_term$model_matrix)
  if(is.null(column_names) || length(column_names) != random_term$n_columns){
    column_names <- paste0("column", seq_len(random_term$n_columns))
  }

  column_names
}

.bt_JAGS_bridge_context_random_sd_nodes <- function(random_term, sd_values){

  sd_names <- random_term$sd_parameter_names
  if(is.null(sd_names) || length(sd_names) != length(sd_values) ||
     any(is.na(sd_names))){
    return(numeric())
  }

  out <- as.numeric(sd_values)
  names(out) <- sd_names
  out
}

.bt_JAGS_bridge_context_random_nodes <- function(random){

  nodes <- numeric()
  for(parameter in names(random)){
    for(block in names(random[[parameter]])){
      nodes <- .bt_JAGS_bridge_merge_nodes(
        nodes,
        random[[parameter]][[block]]$nodes
      )
      allocations <- random[[parameter]][[block]]$allocation$weights
      if(length(allocations) > 0L){
        for(allocation in allocations){
          nodes <- .bt_JAGS_bridge_merge_nodes(nodes, allocation)
        }
      }
    }
  }

  nodes
}
