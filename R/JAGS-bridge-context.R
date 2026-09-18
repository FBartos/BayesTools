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

.bt_JAGS_bridge_context_mode <- function(bridge_context){

  if(is.logical(bridge_context) && length(bridge_context) == 1L &&
     !is.na(bridge_context)){
    return(if(isTRUE(bridge_context)) "full" else "none")
  }
  if(is.character(bridge_context) && length(bridge_context) == 1L &&
     !is.na(bridge_context) &&
     bridge_context %in% c("none", "full", "nodes", "marginal")){
    return(bridge_context)
  }

  stop(
    "'bridge_context' must be FALSE, TRUE, or one of 'none', 'full', 'nodes', and 'marginal'.",
    call. = FALSE
  )
}

.bt_JAGS_bridge_compile_context_evaluator <- function(mode, add_parameters,
                                                       formula_design_list,
                                                       formula_data_list,
                                                       formula_prior_list,
                                                       model_data,
                                                       marginal_random_evaluator = NULL,
                                                       node_names = NULL){

  mode <- .bt_JAGS_bridge_context_mode(mode)
  if(identical(mode, "none")){
    return(list(
      context = function(samples, prior_parameters,
                         formula_prior_parameters, formula_parameters) NULL
    ))
  }

  random_node_names <- unique(unlist(lapply(formula_design_list, function(design){
    random_effects <- .bt_formula_design_random_effects(design)
    unlist(lapply(random_effects, `[[`, "sd_parameter_names"), use.names = FALSE)
  }), use.names = FALSE))
  random_node_names <- random_node_names[
    !is.na(random_node_names) & nzchar(random_node_names)
  ]
  selected_random_only <- !is.null(node_names) &&
    all(node_names %in% random_node_names)
  random_evaluator <- .bt_JAGS_bridge_compile_random_context_evaluator(
    formula_design_list = formula_design_list,
    formula_data_list = formula_data_list,
    formula_prior_list = formula_prior_list,
    model_data = model_data,
    node_names = if(selected_random_only) node_names else NULL
  )
  if(is.null(marginal_random_evaluator)){
    marginal_random_evaluator <- list(
      active = FALSE,
      covariance = function(samples, prior_parameters,
                            formula_prior_parameters, formula_parameters,
                            factor_covariance = TRUE,
                            factor_state = FALSE) list()
    )
  }
  # The node layout and the requested-node selection belong to this compiled
  # evaluator, i.e. to one bridge. A worker receives its own copy with the
  # closure, and a new bridge compiles a new evaluator with an empty cache.
  node_cache <- new.env(parent = emptyenv())
  if(identical(mode, "full")){
    return(list(
      context = function(samples, prior_parameters,
                         formula_prior_parameters, formula_parameters){
        .bt_JAGS_bridge_context(
          samples = samples,
          prior_parameters = prior_parameters,
          formula_prior_parameters = formula_prior_parameters,
          formula_parameters = formula_parameters,
          add_parameters = add_parameters,
          formula_design_list = formula_design_list,
          formula_data_list = formula_data_list,
          formula_prior_list = formula_prior_list,
          model_data = model_data,
          random_context_evaluator = random_evaluator,
          marginal_random_evaluator = marginal_random_evaluator,
          node_cache = node_cache
        )
      }
    ))
  }

  if(identical(mode, "marginal")){
    return(list(
      context = function(samples, prior_parameters,
                         formula_prior_parameters, formula_parameters){
        .bt_JAGS_bridge_marginal_context(
          samples = samples,
          prior_parameters = prior_parameters,
          formula_prior_parameters = formula_prior_parameters,
          formula_parameters = formula_parameters,
          add_parameters = add_parameters,
          random_context_evaluator = random_evaluator,
          marginal_random_evaluator = marginal_random_evaluator,
          node_names = node_names,
          random_only = selected_random_only,
          node_cache = node_cache
        )
      }
    ))
  }

  list(
    context = function(samples, prior_parameters,
                       formula_prior_parameters, formula_parameters){
      .bt_JAGS_bridge_nodes_context(
        samples = samples,
        prior_parameters = prior_parameters,
        formula_prior_parameters = formula_prior_parameters,
        formula_parameters = formula_parameters,
        add_parameters = add_parameters,
        random_context_evaluator = random_evaluator,
        node_names = node_names,
        random_only = selected_random_only,
        node_cache = node_cache
      )
    }
  )
}

.bt_JAGS_bridge_call_log_posterior <- function(log_posterior, parameters,
                                               data, context,
                                               bridge_context, ...){

  if(!identical(bridge_context, FALSE) &&
     !identical(bridge_context, "none")){
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
                                     model_data,
                                     random_context_evaluator = NULL,
                                     marginal_random_evaluator = NULL,
                                     node_cache = NULL){

  state <- .bt_JAGS_bridge_context_state(samples)
  state_matrix <- .bt_JAGS_bridge_context_state_matrix(samples)
  add_parameter_values <- .bt_JAGS_bridge_context_add_parameters(
    samples = samples,
    add_parameters = add_parameters
  )
  if(is.null(random_context_evaluator)){
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
  }else{
    random <- random_context_evaluator$random(
      samples = samples,
      prior_parameters = prior_parameters,
      formula_prior_parameters = formula_prior_parameters,
      formula_parameters = formula_parameters
    )
  }
  random_nodes <- .bt_JAGS_bridge_context_random_nodes(random)
  marginalized_random <- if(is.null(marginal_random_evaluator)){
    list()
  }else{
    marginal_random_evaluator$covariance(
      samples = samples,
      prior_parameters = prior_parameters,
      formula_prior_parameters = formula_prior_parameters,
      formula_parameters = formula_parameters
    )
  }

  nodes <- .bt_JAGS_bridge_context_nodes(
    state = state,
    prior_parameters = prior_parameters,
    formula_prior_parameters = formula_prior_parameters,
    formula_parameters = formula_parameters,
    add_parameter_values = add_parameter_values,
    random_nodes = random_nodes,
    node_cache = node_cache
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
    marginalized_random = marginalized_random,
    node_info = node_info,
    metadata = list(
      source = "JAGS_bridgesampling",
      state_names = names(state),
      node_names = names(nodes),
      node_info_names = node_info$name,
      formula_parameters = names(formula_parameters),
      random_parameters = names(random),
      marginalized_random_parameters = names(marginalized_random)
    )
  )
  class(out) <- c("BayesTools_bridge_context", "list")
  out
}

.bt_JAGS_bridge_marginal_context <- function(
    samples,
    prior_parameters,
    formula_prior_parameters,
    formula_parameters,
    add_parameters,
    random_context_evaluator,
    marginal_random_evaluator,
    node_names = NULL,
    random_only = FALSE,
    node_cache = NULL){

  # Both halves of this state read the same random-effect SD vectors, in
  # different shapes. The node half evaluates them first and leaves them here
  # for the marginal half.
  sd_cache <- new.env(parent = emptyenv())
  nodes <- if(!is.null(node_names) && length(node_names) == 0L){
    numeric()
  }else{
    .bt_JAGS_bridge_nodes_context(
      samples = samples,
      prior_parameters = prior_parameters,
      formula_prior_parameters = formula_prior_parameters,
      formula_parameters = formula_parameters,
      add_parameters = add_parameters,
      random_context_evaluator = random_context_evaluator,
      node_names = node_names,
      random_only = random_only,
      node_cache = node_cache,
      sd_cache = sd_cache
    )$nodes
  }
  marginalized_random <- marginal_random_evaluator$covariance(
    samples = samples,
    prior_parameters = prior_parameters,
    formula_prior_parameters = formula_prior_parameters,
    formula_parameters = formula_parameters,
    factor_covariance = FALSE,
    factor_state = TRUE,
    sd_cache = sd_cache
  )

  out <- list(
    nodes = nodes,
    marginalized_random = marginalized_random
  )
  class(out) <- c(
    "BayesTools_bridge_marginal_context",
    "BayesTools_bridge_context",
    "list"
  )
  out
}

.bt_JAGS_bridge_nodes_context <- function(samples, prior_parameters,
                                          formula_prior_parameters,
                                          formula_parameters,
                                          add_parameters,
                                          random_context_evaluator,
                                          node_names = NULL,
                                          random_only = FALSE,
                                          node_cache = NULL,
                                          sd_cache = NULL){

  state <- .bt_JAGS_bridge_context_state(samples)
  add_parameter_values <- .bt_JAGS_bridge_context_add_parameters(
    samples = samples,
    add_parameters = add_parameters
  )
  random_nodes <- random_context_evaluator$nodes(
    samples = samples,
    prior_parameters = prior_parameters,
    formula_prior_parameters = formula_prior_parameters,
    formula_parameters = formula_parameters,
    sd_cache = sd_cache
  )
  if(isTRUE(random_only)){
    nodes <- .bt_JAGS_bridge_select_nodes(random_nodes, node_names, node_cache)
    out <- list(nodes = nodes)
    class(out) <- c(
      "BayesTools_bridge_nodes_context",
      "BayesTools_bridge_context",
      "list"
    )
    return(out)
  }
  nodes <- .bt_JAGS_bridge_context_nodes(
    state = state,
    prior_parameters = prior_parameters,
    formula_prior_parameters = formula_prior_parameters,
    formula_parameters = formula_parameters,
    add_parameter_values = add_parameter_values,
    random_nodes = random_nodes,
    node_cache = node_cache
  )
  nodes <- .bt_JAGS_bridge_select_nodes(nodes, node_names, node_cache)

  out <- list(nodes = nodes)
  class(out) <- c(
    "BayesTools_bridge_nodes_context",
    "BayesTools_bridge_context",
    "list"
  )
  out
}

.bt_JAGS_bridge_select_nodes <- function(nodes, node_names, node_cache = NULL){

  if(is.null(node_names)){
    return(nodes)
  }
  # The requested names and the available ones are both fixed by the model, so
  # the availability check and the name lookup are resolved once per bridge and
  # replayed as positions while the available names stay the same. Name lookup
  # and `match()` take the first occurrence of a duplicated name alike.
  if(is.environment(node_cache)){
    available <- names(nodes)
    if(!identical(available, node_cache$select_names) ||
       !identical(node_names, node_cache$select_requested)){
      .bt_JAGS_bridge_check_selected_nodes(available, node_names)
      node_cache$select_names <- available
      node_cache$select_requested <- node_names
      node_cache$select_index <- match(node_names, available)
    }
    return(nodes[node_cache$select_index])
  }
  .bt_JAGS_bridge_check_selected_nodes(names(nodes), node_names)
  nodes[node_names]
}

.bt_JAGS_bridge_check_selected_nodes <- function(available, node_names){

  missing <- setdiff(node_names, available)
  if(length(missing) > 0L){
    stop(
      "Requested bridge context node(s) are unavailable: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(TRUE)
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
                                          random_nodes,
                                          node_cache = NULL){

  if(is.environment(node_cache)){
    entries <- .bt_JAGS_bridge_node_entries(
      prior_parameters,
      formula_prior_parameters,
      formula_parameters,
      add_parameter_values
    )
    layout <- .bt_JAGS_bridge_node_layout(
      node_cache = node_cache,
      state = state,
      entries = entries,
      random_nodes = random_nodes
    )
    if(!is.null(layout)){
      return(.bt_JAGS_bridge_fill_node_layout(
        layout = layout,
        state = state,
        entries = entries,
        random_nodes = random_nodes
      ))
    }
  }

  .bt_JAGS_bridge_merge_nodes(
    state,
    .bt_JAGS_bridge_flatten_parameter_nodes(prior_parameters),
    .bt_JAGS_bridge_flatten_parameter_nodes(formula_prior_parameters),
    .bt_JAGS_bridge_flatten_parameter_nodes(formula_parameters),
    .bt_JAGS_bridge_flatten_parameter_nodes(add_parameter_values),
    random_nodes
  )
}


# The pieces the reconstructed parameter lists contribute to the merged node
# vector, in the order .bt_JAGS_bridge_flatten_parameter_nodes() merges them,
# with the parameter name, value length and value dimensions that decide their
# node names. Named lookup keeps a duplicated parameter name resolving to the
# same value the flattening loop takes.
.bt_JAGS_bridge_node_entries <- function(...){

  parameter_lists <- list(...)
  n_max <- sum(vapply(parameter_lists, length, integer(1L)))
  parameter <- character(n_max)
  value_length <- integer(n_max)
  dimension <- vector("list", n_max)
  value <- vector("list", n_max)
  kept <- 0L

  for(parameters in parameter_lists){
    if(length(parameters) == 0L){
      next
    }
    parameter_names <- names(parameters)
    if(is.null(parameter_names)){
      next
    }
    for(name in parameter_names){
      if(!nzchar(name)){
        next
      }
      parameter_value <- parameters[[name]]
      if(is.null(parameter_value) || length(parameter_value) == 0L ||
         (!is.numeric(parameter_value) && !is.logical(parameter_value))){
        next
      }
      kept <- kept + 1L
      parameter[[kept]] <- name
      value_length[[kept]] <- length(parameter_value)
      # A vector has no dimensions; `[[<-` would drop the entry for NULL.
      dimension[kept] <- list(dim(parameter_value))
      value[[kept]] <- parameter_value
    }
  }

  kept_index <- seq_len(kept)

  list(
    parameter = parameter[kept_index],
    length = value_length[kept_index],
    dimension = dimension[kept_index],
    value = value[kept_index]
  )
}


# The node layout a bridge merges is a property of the model: the names,
# lengths and dimensions of the reconstructed pieces are fixed by the formula
# and the monitored coordinates, and only their values change between draws.
# The layout is therefore built once per bridge and replayed by position, which
# gives the vector .bt_JAGS_bridge_merge_nodes() builds by name -- first
# appearance decides the position, the last assignment decides the value, and
# one positional assignment applies the pieces in that same order -- without
# rebuilding the names. A draw whose pieces have another shape rebuilds the
# layout; a piece with a missing or empty node name has no position rule here
# and returns to the by-name merge.
.bt_JAGS_bridge_node_layout <- function(node_cache, state, entries,
                                        random_nodes){

  state_names  <- names(state)
  random_names <- names(random_nodes)
  shape <- list(
    state_names,
    entries$parameter,
    entries$length,
    entries$dimension,
    random_names
  )
  if(identical(shape, node_cache$nodes_shape)){
    return(node_cache$nodes_layout)
  }

  entry_names <- vector("list", length(entries$parameter))
  for(entry in seq_along(entry_names)){
    entry_names[[entry]] <- .bt_JAGS_bridge_flatten_parameter_names(
      parameter = entries$parameter[[entry]],
      value = entries$value[[entry]]
    )
  }
  all_names <- c(
    state_names,
    unlist(entry_names, use.names = FALSE),
    random_names
  )

  layout <- if(length(all_names) == 0L || anyNA(all_names) ||
               !all(nzchar(all_names))){
    NULL
  }else{
    merged <- unique(all_names)
    list(
      names = merged,
      index = match(all_names, merged),
      state = !is.null(state_names),
      random = !is.null(random_names)
    )
  }
  node_cache$nodes_shape <- shape
  node_cache$nodes_layout <- layout

  layout
}


.bt_JAGS_bridge_fill_node_layout <- function(layout, state, entries,
                                             random_nodes){

  values <- c(
    if(isTRUE(layout$state)) as.numeric(state) else NULL,
    unlist(entries$value, use.names = FALSE),
    if(isTRUE(layout$random)) as.numeric(random_nodes) else NULL
  )
  index <- layout$index
  if(length(values) != length(index)){
    stop("Bridge context node layout and values are inconsistent.",
         call. = FALSE)
  }

  out <- numeric(length(layout$names))
  names(out) <- layout$names
  out[index] <- values

  out
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
  names(values) <- .bt_JAGS_bridge_flatten_parameter_names(
    parameter = parameter,
    value = value
  )
  values
}


.bt_JAGS_bridge_flatten_parameter_names <- function(parameter, value){

  dims <- dim(value)
  if(is.null(dims)){
    if(length(value) == 1L){
      return(parameter)
    }

    return(paste0(parameter, "[", seq_along(value), "]"))
  }

  indices <- do.call(expand.grid, lapply(dims, seq_len))

  paste0(
    parameter,
    "[",
    apply(indices, 1L, paste, collapse = ","),
    "]"
  )
}

.bt_JAGS_bridge_merge_nodes <- function(...){

  # First appearance decides the position, the last assignment decides the
  # value. Assignment by name is the fast path for the short named vectors a
  # bridge merges for every draw: `[<-` matches and grows them in C, while a
  # single pass through `unlist()`, `duplicated()` and `match()` costs three
  # times as much on random-effects contexts (0.26 s against 0.75 s of a
  # 6 s marginal likelihood for an `us()` structure).
  out <- numeric()
  for(piece in list(...)){
    if(length(piece) == 0L || is.null(names(piece))){
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

  evaluator <- .bt_JAGS_bridge_compile_random_context_evaluator(
    formula_design_list = formula_design_list,
    formula_data_list = formula_data_list,
    formula_prior_list = formula_prior_list,
    model_data = model_data
  )
  evaluator$random(
    samples = samples,
    prior_parameters = prior_parameters,
    formula_prior_parameters = formula_prior_parameters,
    formula_parameters = formula_parameters
  )
}

.bt_JAGS_bridge_compile_random_context_evaluator <- function(
    formula_design_list,
    formula_data_list,
    formula_prior_list,
    model_data,
    node_names = NULL){

  plans <- list()
  sampled_random_parameters <- character()
  parameter_names <- names(formula_design_list)
  if(is.null(parameter_names)){
    parameter_names <- rep("", length(formula_design_list))
  }

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
      sampled_random_parameters <- c(sampled_random_parameters, parameter)
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
    random_terms <- .bt_formula_design_random_effects(design)
    term_plans <- lapply(random_terms, function(random_term){
      row_indexed <- .bt_random_effect_has_row_indexed_external_sd(random_term)
      list(
        random_term = random_term,
        row_indexed = row_indexed,
        row_source_evaluator = if(row_indexed){
          .bt_JAGS_bridge_compile_row_source_evaluator(
            random_term = random_term,
            n_rows = nrow(random_term$model_matrix),
            context = "Bridge context"
          )
        }else{
          NULL
        },
        allocation_parameters =
          .bt_JAGS_marglik_random_effect_allocation_parameters(random_term),
        sd_evaluator = .bt_JAGS_bridge_compile_random_sd_evaluator(
          random_term = random_term,
          prior_list = parameter_prior_list
        )
      )
    })
    if(!is.null(node_names)){
      term_plans <- Filter(function(term_plan){
        sd_names <- term_plan$random_term$sd_parameter_names
        any(sd_names[!is.na(sd_names)] %in% node_names)
      }, term_plans)
    }
    plans[[length(plans) + 1L]] <- list(
      parameter = parameter,
      term_plans = term_plans,
      prior_list = parameter_prior_list,
      source_data = .bt_JAGS_marglik_parameter_source_data(
        model_data = model_data,
        formula_data = formula_data,
        design = design
      )
    )
  }
  sampled_random_parameters <- unique(sampled_random_parameters)

  source_parameters <- function(samples, prior_parameters,
                                formula_prior_parameters,
                                formula_parameters){
    out <- .bt_JAGS_bridge_context_source_parameters(
      samples = samples,
      prior_parameters = prior_parameters,
      formula_prior_parameters = formula_prior_parameters,
      formula_parameters = formula_parameters
    )
    .bt_parameter_source_forbid_formula_parameters(
      out,
      sampled_random_parameters
    )
  }

  list(
    random = function(samples, prior_parameters,
                      formula_prior_parameters, formula_parameters){
      if(length(plans) == 0L){
        return(list())
      }
      parameter_sources <- source_parameters(
        samples = samples,
        prior_parameters = prior_parameters,
        formula_prior_parameters = formula_prior_parameters,
        formula_parameters = formula_parameters
      )
      posterior <- .bt_JAGS_marglik_random_effect_posterior_row(samples)
      out <- list()
      for(plan in plans){
        parameter_random <- list()
        for(term_plan in plan$term_plans){
          random_term <- term_plan$random_term
          parameter_random[[random_term$block_name]] <-
            .bt_JAGS_bridge_context_random_block(
              samples = samples,
              random_term = random_term,
              prior_list = plan$prior_list,
              formula_prior_parameters = formula_prior_parameters,
              data = plan$source_data,
              parameters = parameter_sources,
              posterior = posterior,
              row_indexed = term_plan$row_indexed,
              row_source_evaluator = term_plan$row_source_evaluator,
              allocation_parameters = term_plan$allocation_parameters,
              sd_evaluator = term_plan$sd_evaluator
            )
        }
        out[[plan$parameter]] <- parameter_random
      }
      out
    },
    nodes = function(samples, prior_parameters,
                     formula_prior_parameters, formula_parameters,
                     sd_cache = NULL){
      if(length(plans) == 0L){
        return(numeric())
      }
      parameter_sources <- source_parameters(
        samples = samples,
        prior_parameters = prior_parameters,
        formula_prior_parameters = formula_prior_parameters,
        formula_parameters = formula_parameters
      )
      posterior <- .bt_JAGS_marglik_random_effect_posterior_row(samples)
      nodes <- numeric()
      for(plan in plans){
        for(term_plan in plan$term_plans){
          random_term <- term_plan$random_term
          block <- .bt_JAGS_bridge_context_random_block_nodes(
            samples = samples,
            random_term = random_term,
            prior_list = plan$prior_list,
            data = plan$source_data,
            parameters = parameter_sources,
            posterior = posterior,
            row_indexed = term_plan$row_indexed,
            row_source_evaluator = term_plan$row_source_evaluator,
            sd_evaluator = term_plan$sd_evaluator,
            sd_cache = sd_cache,
            sd_share_key = .bt_JAGS_bridge_random_sd_share_key(
              plan$parameter,
              random_term$block_name
            )
          )
          nodes <- .bt_JAGS_bridge_merge_nodes(nodes, block$nodes)
          for(allocation in block$allocations){
            nodes <- .bt_JAGS_bridge_merge_nodes(nodes, allocation)
          }
        }
      }
      nodes
    }
  )
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

.bt_JAGS_bridge_compile_random_sd_evaluator <- function(random_term,
                                                        prior_list,
                                                        posterior_names = NULL){

  binding <- random_term$sd_binding
  if(is.null(binding)){
    return(NULL)
  }
  .bt_check_random_sd_binding(binding)
  if(.bt_random_effect_has_row_indexed_external_sd(random_term)){
    return(NULL)
  }

  allocation_parameters <-
    .bt_JAGS_marglik_random_effect_allocation_parameters(random_term)
  ambiguity <- lapply(allocation_parameters, function(parameter_name){
    prior <- prior_list[[parameter_name]]
    if(is.null(prior) || !is.prior.simplex(prior) ||
       !identical(prior$distribution, "dirichlet")){
      return(NULL)
    }
    K <- prior$parameters[["K"]]
    list(
      weight = paste0(parameter_name, "[", seq_len(K), "]"),
      eta = paste0(
        .JAGS_prior_dirichlet_eta_name(parameter_name),
        "[", seq_len(K), "]"
      ),
      parameter = parameter_name
    )
  })
  ambiguity <- Filter(Negate(is.null), ambiguity)

  check_ambiguity_names <- function(sample_names){
    for(spec in ambiguity){
      if(all(spec$weight %in% sample_names) &&
         all(spec$eta %in% sample_names)){
        stop(
          "Bridge samples contain both normalized Dirichlet allocation coordinates ",
          "and auxiliary eta coordinates for parameter '",
          spec$parameter,
          "'. Use the auxiliary eta bridge coordinates only.",
          call. = FALSE
        )
      }
    }
    invisible(TRUE)
  }
  check_ambiguity <- function(samples){
    sample_names <- if(is.matrix(samples)) colnames(samples) else names(samples)
    check_ambiguity_names(sample_names)
  }

  evaluate <- if(isTRUE(binding$true_allocation)){
    allocation <- binding$allocations[[1L]]
    target <- .bt_random_effect_allocation_target_metadata(allocation)
    .bt_random_effect_allocation_scale_metadata(allocation)
    factors <- if(identical(target, "sd_component")){
      .bt_random_effect_allocation_parent_factors_metadata(allocation)
    }else{
      .bt_random_effect_allocation_factors_metadata(allocation)
    }
    factor_plan <- .bt_random_effect_compile_allocation_factor_plan(factors)
    source_name <- .bt_random_sd_binding_source_name(allocation$source)
    source_evaluator <- .bt_JAGS_bridge_compile_parameter_draw_evaluator(
      parameter_name = source_name,
      prior_list = prior_list,
      posterior_names = posterior_names
    )
    factor_evaluator <- .bt_JAGS_bridge_compile_allocation_factor_evaluator(
      factor_plan = factor_plan,
      prior_list = prior_list,
      posterior_names = posterior_names
    )
    chain_coordinates <- c(
      source_name,
      vapply(factor_plan, function(factor){
        if(is.null(factor$weight_name)) NA_character_ else factor$weight_name
      }, character(1))
    )

    if(identical(target, "sd_component")){
      component <- .bt_check_random_sd_component_allocation(
        allocation = allocation,
        n_columns = random_term$n_columns,
        context = "Bridge-sampling random-effect allocation metadata"
      )
      leaf_index <- component$leaf_index_by_column
      n_targets <- component$n_targets
      weight_name <- allocation$weight_name
      scale <- allocation$scale
      weight_evaluator <- .bt_JAGS_bridge_compile_dirichlet_draw_evaluator(
        parameter_name = weight_name,
        prior_list = prior_list,
        posterior_names = posterior_names
      )
      multiplier <- .bt_JAGS_bridge_compile_allocation_multiplier(
        scale = scale,
        n_targets = n_targets
      )
      chain_coordinates <- c(chain_coordinates, weight_name)
    }

    function(posterior, parameters = NULL, prefer_weights = FALSE){
      base <- source_evaluator(posterior, parameters)
      if(is.null(base)){
        return(NULL)
      }
      base <- factor_evaluator(
        base = base,
        posterior = posterior,
        parameters = parameters,
        prefer_weights = prefer_weights
      )
      if(is.null(base)){
        return(NULL)
      }
      if(!identical(target, "sd_component")){
        return(matrix(
          base,
          nrow = nrow(posterior),
          ncol = random_term$n_columns
        ))
      }

      weights <- weight_evaluator(
        posterior,
        parameters,
        prefer_weights = prefer_weights
      )
      if(is.null(weights)){
        return(NULL)
      }
      if(ncol(weights) != n_targets){
        stop(
          "Random-effect allocation metadata for '", weight_name,
          "' expected ", n_targets, " Dirichlet coordinate(s), but found ",
          ncol(weights), ".",
          call. = FALSE
        )
      }
      out <- matrix(
        NA_real_,
        nrow = nrow(posterior),
        ncol = random_term$n_columns
      )
      for(column in seq_len(random_term$n_columns)){
        out[, column] <- base * multiplier(weights[, leaf_index[column]])
      }
      out
    }
  }else{
    sd_names <- random_term$sd_parameter_names
    sd_evaluators <- if(
      is.null(sd_names) || length(sd_names) != random_term$n_columns ||
      any(is.na(sd_names))
    ){
      NULL
    }else{
      lapply(
        sd_names,
        .bt_JAGS_bridge_compile_parameter_draw_evaluator,
        prior_list = prior_list,
        posterior_names = posterior_names
      )
    }
    chain_coordinates <- if(is.null(sd_evaluators)) character() else sd_names
    function(posterior, parameters = NULL, prefer_weights = FALSE){
      if(is.null(sd_evaluators)){
        return(NULL)
      }
      out <- matrix(
        NA_real_,
        nrow = nrow(posterior),
        ncol = random_term$n_columns
      )
      for(column in seq_len(random_term$n_columns)){
        values <- sd_evaluators[[column]](posterior, parameters)
        if(is.null(values)){
          return(NULL)
        }
        out[, column] <- values
      }
      out
    }
  }

  values <- function(samples, parameters = NULL){
    check_ambiguity(samples)
    posterior <- .bt_JAGS_marglik_random_effect_posterior_row(samples)
    sd_draws <- tryCatch(
      evaluate(posterior, parameters),
      error = function(e){
        if(.bt_JAGS_marglik_random_effect_support_error(e)){
          .bt_JAGS_marglik_out_of_support(conditionMessage(e))
        }
        stop(e)
      }
    )
    if(is.null(sd_draws)){
      .bt_JAGS_marglik_explain_random_effect_sd_missing(
        samples,
        random_term,
        prior_list
      )
      stop(
        "Random-effect SD metadata are incomplete for block '",
        random_term$block_name,
        "'.",
        call. = FALSE
      )
    }
    unname(sd_draws[1L, ])
  }

  direct_names <- if(!isTRUE(binding$true_allocation)){
    random_term$sd_parameter_names
  }else{
    NULL
  }
  direct_indices <- NULL
  posterior_values <- function(posterior, parameters = NULL){
    draws <- evaluate(
      posterior,
      parameters = parameters,
      prefer_weights = TRUE
    )
    if(is.null(draws)){
      return(NULL)
    }
    unname(draws[1L, ])
  }

  posterior_draws <- function(posterior, parameters = NULL){
    posterior <- as.matrix(posterior)
    posterior_names <- colnames(posterior)
    direct <- !is.null(direct_names) &&
      length(direct_names) == random_term$n_columns &&
      !anyNA(direct_names) && length(ambiguity) == 0L
    if(direct){
      if(is.list(parameters) &&
         all(direct_names %in% names(parameters))){
        parameter_values <- lapply(
          parameters[direct_names],
          as.numeric
        )
        parameter_lengths <- lengths(parameter_values)
        if(all(parameter_lengths %in% c(1L, nrow(posterior)))){
          return(unname(do.call(cbind, lapply(parameter_values, function(x){
            if(length(x) == 1L) rep(x, nrow(posterior)) else x
          }))))
        }
      }
      direct <- !is.null(posterior_names)
    }
    if(direct){
      indices_valid <- !is.null(direct_indices) &&
        length(direct_indices) == length(direct_names) &&
        all(direct_indices >= 1L) &&
        all(direct_indices <= length(posterior_names)) &&
        identical(posterior_names[direct_indices], direct_names)
      if(!indices_valid){
        direct_indices <<- match(direct_names, posterior_names)
      }
      if(!anyNA(direct_indices)){
        return(unname(posterior[, direct_indices, drop = FALSE]))
      }
    }

    evaluate(
      posterior,
      parameters = parameters,
      prefer_weights = TRUE
    )
  }

  # What the compiled chain reads, so a caller can decide - from the metadata
  # alone - whether two evaluators of the same term bind the same posterior
  # columns. 'coordinates' are the names the chain looks up in a supplied
  # parameter list; 'ambiguity' names the Dirichlet allocations whose
  # normalized weights and auxiliary eta coordinates are distinguished by
  # 'prefer_weights'.
  bindings <- list(
    block_name = random_term$block_name,
    n_columns = random_term$n_columns,
    coordinates = unique(chain_coordinates[!is.na(chain_coordinates)]),
    ambiguity = ambiguity,
    index_cache = new.env(parent = emptyenv())
  )

  list(
    values = values,
    posterior_values = posterior_values,
    posterior_draws = posterior_draws,
    bindings = bindings
  )
}


.bt_JAGS_bridge_random_sd_share_key <- function(parameter, block_name){

  if(!is.character(parameter) || length(parameter) != 1L ||
     !is.character(block_name) || length(block_name) != 1L ||
     is.na(parameter) || is.na(block_name)){
    return(NULL)
  }
  paste0(parameter, "\r", block_name)
}


# The two evaluators of one term differ only in 'prefer_weights', and that
# preference selects a posterior column only when a Dirichlet allocation
# carries both its normalized weights and its auxiliary eta coordinates.
# Nothing is evaluated here: the verdict comes from the compiled bindings and
# the bridge's own coordinate names.
.bt_JAGS_bridge_random_sd_shared_binding <- function(bindings, posterior_names){

  if(is.null(bindings) || is.null(posterior_names)){
    return(FALSE)
  }
  for(spec in bindings$ambiguity){
    if(all(spec$weight %in% posterior_names) &&
       all(spec$eta %in% posterior_names)){
      return(FALSE)
    }
  }
  TRUE
}


# What the posterior row alone gives for one coordinate of an SD chain: its own
# column, or, for a Dirichlet allocation parameter, the normalized weights the
# draw evaluator reconstructs from the coordinates the row carries. NULL when
# the row determines no value.
.bt_JAGS_bridge_random_sd_posterior_value <- function(coordinate, bindings,
                                                      posterior){

  posterior_names <- colnames(posterior)
  # A bridge replays thousands of states through the same column layout, so
  # the name lookup belongs to the compiled bindings, not to the state.
  cache <- bindings$index_cache
  resolved <- if(is.environment(cache) &&
                 exists(coordinate, envir = cache, inherits = FALSE)){
    entry <- get(coordinate, envir = cache, inherits = FALSE)
    if(identical(entry$names, posterior_names)) entry else NULL
  }else{
    NULL
  }
  if(is.null(resolved)){
    resolved <- list(names = posterior_names, kind = "none", index = integer())
    index <- match(coordinate, posterior_names)
    if(!is.na(index)){
      resolved$kind <- "column"
      resolved$index <- index
    }else{
      for(spec in bindings$ambiguity){
        if(!identical(spec$parameter, coordinate)){
          next
        }
        if(all(spec$weight %in% posterior_names)){
          resolved$kind <- "column"
          resolved$index <- match(spec$weight, posterior_names)
        }else if(all(spec$eta %in% posterior_names)){
          resolved$kind <- "eta"
          resolved$index <- match(spec$eta, posterior_names)
        }
        break
      }
    }
    if(is.environment(cache)){
      assign(coordinate, resolved, envir = cache)
    }
  }
  if(identical(resolved$kind, "column")){
    return(as.numeric(posterior[, resolved$index]))
  }
  if(identical(resolved$kind, "eta")){
    eta <- posterior[, resolved$index, drop = FALSE]
    return(as.numeric(eta / rowSums(eta)))
  }
  NULL
}


# A shared SD vector must also have been resolved from the same values. The
# node evaluator consults the state's parameter sources first, the marginal
# evaluator reads the posterior row, so a source that carries another value
# than the row does keeps the two evaluations apart.
.bt_JAGS_bridge_random_sd_sources_agree <- function(bindings, parameters,
                                                    posterior){

  coordinates <- bindings$coordinates
  if(length(coordinates) == 0L || !is.list(parameters)){
    return(TRUE)
  }
  coordinates <- coordinates[coordinates %in% names(parameters)]
  if(length(coordinates) == 0L){
    return(TRUE)
  }
  for(coordinate in coordinates){
    supplied <- as.numeric(parameters[[coordinate]])
    reference <- .bt_JAGS_bridge_random_sd_posterior_value(
      coordinate = coordinate,
      bindings = bindings,
      posterior = posterior
    )
    if(is.null(reference)){
      return(FALSE)
    }
    if(length(supplied) == 1L && length(reference) > 1L){
      supplied <- rep(supplied, length(reference))
    }
    if(!identical(supplied, reference)){
      return(FALSE)
    }
  }
  TRUE
}

.bt_JAGS_bridge_compile_parameter_draw_evaluator <- function(
    parameter_name, prior_list, posterior_names = NULL){

  prior_name <- sub("\\[[0-9]+\\]$", "", parameter_name)
  prior <- prior_list[[prior_name]]
  fixed <- if(!is.null(prior) && is.prior.point(prior)){
    location <- prior$parameters[["location"]]
    if(length(location) == 1L && !is.na(location)) location else NULL
  }else{
    NULL
  }
  fixed_posterior_names <- !is.null(posterior_names)
  cached_posterior_names <- posterior_names
  posterior_index <- if(fixed_posterior_names){
    match(parameter_name, posterior_names)
  }else{
    NA_integer_
  }
  force(parameter_name)
  force(fixed)

  function(posterior, parameters = NULL){
    if(is.list(parameters) && parameter_name %in% names(parameters)){
      value <- as.numeric(parameters[[parameter_name]])
      if(length(value) == 1L){
        return(rep(value, nrow(posterior)))
      }
      if(length(value) == nrow(posterior)){
        return(value)
      }
    }
    if(!fixed_posterior_names){
      current_names <- colnames(posterior)
      if(!identical(current_names, cached_posterior_names)){
        posterior_index <<- match(parameter_name, current_names)
        cached_posterior_names <<- current_names
      }
    }
    if(!is.na(posterior_index)){
      return(posterior[, posterior_index])
    }
    if(!is.null(fixed)){
      return(rep(fixed, nrow(posterior)))
    }
    NULL
  }
}

.bt_JAGS_bridge_compile_dirichlet_draw_evaluator <- function(
    parameter_name, prior_list, posterior_names = NULL){

  prior <- prior_list[[parameter_name]]
  if(is.null(prior) || !is.prior.simplex(prior) ||
     !identical(prior$distribution, "dirichlet")){
    return(function(posterior, parameters = NULL, prefer_weights = FALSE) NULL)
  }
  K <- prior$parameters[["K"]]
  eta_names <- paste0(
    .JAGS_prior_dirichlet_eta_name(parameter_name),
    "[", seq_len(K), "]"
  )
  weight_names <- paste0(parameter_name, "[", seq_len(K), "]")
  eta_cache_key <- .bt_random_effect_dirichlet_cache_key(
    parameter_name,
    K,
    "eta"
  )
  weight_cache_key <- .bt_random_effect_dirichlet_cache_key(
    parameter_name,
    K,
    "weights"
  )
  fixed_posterior_names <- !is.null(posterior_names)
  cached_posterior_names <- posterior_names
  eta_indices <- if(fixed_posterior_names){
    match(eta_names, posterior_names)
  }else{
    rep(NA_integer_, K)
  }
  weight_indices <- if(fixed_posterior_names){
    match(weight_names, posterior_names)
  }else{
    rep(NA_integer_, K)
  }
  force(parameter_name)
  force(K)

  function(posterior, parameters = NULL, prefer_weights = FALSE){
    if(is.list(parameters) && parameter_name %in% names(parameters)){
      value <- as.numeric(parameters[[parameter_name]])
      return(matrix(
        value,
        nrow = nrow(posterior),
        ncol = length(value),
        byrow = TRUE
      ))
    }
    if(!fixed_posterior_names){
      current_names <- colnames(posterior)
      if(!identical(current_names, cached_posterior_names)){
        eta_indices     <<- match(eta_names, current_names)
        weight_indices  <<- match(weight_names, current_names)
        cached_posterior_names <<- current_names
      }
    }
    cache <- .bt_random_effect_dirichlet_draw_cache(posterior)
    if(isTRUE(prefer_weights) && !anyNA(weight_indices)){
      if(!is.null(cache) &&
         exists(weight_cache_key, envir = cache, inherits = FALSE)){
        return(get(weight_cache_key, envir = cache, inherits = FALSE))
      }
      weights <- .bt_random_effect_validate_dirichlet_weights(
        weights = posterior[, weight_indices, drop = FALSE],
        parameter_name = parameter_name
      )
      .bt_random_effect_dirichlet_cache_assign(
        cache,
        weight_cache_key,
        weights
      )
      return(weights)
    }
    if(!anyNA(eta_indices)){
      if(!is.null(cache) &&
         exists(eta_cache_key, envir = cache, inherits = FALSE)){
        return(get(eta_cache_key, envir = cache, inherits = FALSE))
      }
      eta <- posterior[, eta_indices, drop = FALSE]
      eta_sum <- rowSums(eta)
      invalid <- !is.finite(eta) | eta <= 0
      invalid_row <- !is.finite(eta_sum) | eta_sum <= 0
      if(any(invalid) || any(invalid_row)){
        if(any(invalid)){
          invalid_column <- col(eta)[which(invalid)[1L]]
          detail <- paste0(" for '", eta_names[invalid_column], "'")
        }else{
          detail <- paste0(" at row ", which(invalid_row)[1L])
        }
        .bt_random_effect_allocation_out_of_support(
          "Random-effect Dirichlet allocation auxiliary samples must be finite and positive",
          detail,
          "."
        )
      }
      weights <- eta / eta_sum
      return(.bt_random_effect_dirichlet_cache_return(
        weights = weights,
        cache = cache,
        cache_key = eta_cache_key
      ))
    }
    if(!anyNA(weight_indices)){
      if(!is.null(cache) &&
         exists(weight_cache_key, envir = cache, inherits = FALSE)){
        return(get(weight_cache_key, envir = cache, inherits = FALSE))
      }
      weights <- .bt_random_effect_validate_dirichlet_weights(
        weights = posterior[, weight_indices, drop = FALSE],
        parameter_name = parameter_name
      )
      .bt_random_effect_dirichlet_cache_assign(
        cache,
        weight_cache_key,
        weights
      )
      return(weights)
    }
    NULL
  }
}

.bt_JAGS_bridge_compile_allocation_factor_evaluator <- function(
    factor_plan, prior_list, posterior_names = NULL){

  if(length(factor_plan) == 0L){
    return(function(base, posterior, parameters = NULL,
                    prefer_weights = FALSE) base)
  }
  plans <- lapply(factor_plan, function(factor){
    list(
      weight_evaluator = if(is.null(factor$weight_name)) NULL else
        .bt_JAGS_bridge_compile_dirichlet_draw_evaluator(
          parameter_name = factor$weight_name,
          prior_list = prior_list,
          posterior_names = posterior_names
        ),
      multiplier = if(is.null(factor$weight_name)) NULL else
        .bt_JAGS_bridge_compile_allocation_multiplier(
          scale = factor$scale,
          n_targets = factor$n_targets
        ),
      index = factor$index,
      n_targets = factor$n_targets,
      weight_name = factor$weight_name,
      inclusion_name = factor$inclusion_name,
      gate_evaluator = .bt_JAGS_bridge_compile_allocation_gate_evaluator(
        factor$inclusion_name,
        posterior_names = posterior_names
      )
    )
  })
  force(plans)

  function(base, posterior, parameters = NULL, prefer_weights = FALSE){
    out <- base
    for(plan in plans){
      multiplier <- 1
      if(!is.null(plan$weight_evaluator)){
        weights <- plan$weight_evaluator(
          posterior,
          parameters,
          prefer_weights = prefer_weights
        )
        if(is.null(weights)){
          return(NULL)
        }
        if(ncol(weights) != plan$n_targets || plan$index > ncol(weights)){
          stop(
            "Random-effect allocation factor metadata for '",
            plan$weight_name,
            "' do not match the reconstructed Dirichlet coordinates.",
            call. = FALSE
          )
        }
        multiplier <- plan$multiplier(weights[, plan$index])
      }
      gate <- plan$gate_evaluator(posterior)
      if(is.null(gate)){
        stop(
          "Random-effect allocation inclusion samples are missing Bernoulli indicator '",
          plan$inclusion_name,
          "'.",
          call. = FALSE
        )
      }
      out <- out * multiplier * gate
    }
    out
  }
}

.bt_JAGS_bridge_compile_allocation_gate_evaluator <- function(
    parameter_name, posterior_names = NULL){

  if(is.null(parameter_name)){
    return(function(posterior) rep(1, nrow(posterior)))
  }
  fixed_posterior_names <- !is.null(posterior_names)
  cached_posterior_names <- posterior_names
  posterior_index <- if(fixed_posterior_names){
    match(parameter_name, posterior_names)
  }else{
    NA_integer_
  }
  force(parameter_name)

  function(posterior){
    if(!fixed_posterior_names){
      current_names <- colnames(posterior)
      if(!identical(current_names, cached_posterior_names)){
        posterior_index <<- match(parameter_name, current_names)
        cached_posterior_names <<- current_names
      }
    }
    if(is.na(posterior_index)){
      return(NULL)
    }
    values <- as.numeric(posterior[, posterior_index])
    invalid <- !is.finite(values) | !(values %in% c(0, 1))
    if(any(invalid)){
      .bt_random_effect_allocation_out_of_support(
        "Random-effect allocation inclusion samples for '",
        parameter_name,
        "' must be Bernoulli indicators."
      )
    }
    values
  }
}

.bt_JAGS_bridge_compile_allocation_multiplier <- function(scale, n_targets){

  if(identical(scale, "mean_variance")){
    if(!is.numeric(n_targets) || length(n_targets) != 1L ||
       is.na(n_targets) || n_targets < 1L){
      stop(
        "Random-effect allocation metadata are missing canonical 'allocation$n_targets'.",
        call. = FALSE
      )
    }
    force(n_targets)
    return(function(weights) sqrt(n_targets * weights))
  }
  if(identical(scale, "total_variance")){
    return(function(weights) sqrt(weights))
  }
  stop(
    "Random-effect allocation metadata are missing canonical 'allocation$scale'.",
    call. = FALSE
  )
}

.bt_JAGS_bridge_context_random_block <- function(samples, random_term,
                                                 prior_list,
                                                 formula_prior_parameters,
                                                 data,
                                                 parameters,
                                                 posterior = NULL,
                                                 row_indexed = NULL,
                                                 allocation_parameters = NULL,
                                                 sd_evaluator = NULL,
                                                 row_source_evaluator = NULL){

  if(is.null(posterior)){
    posterior <- .bt_JAGS_marglik_random_effect_posterior_row(samples)
  }
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
    formula_prior_parameters = formula_prior_parameters,
    parameter_names = allocation_parameters
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

  if(is.null(row_indexed)){
    row_indexed <- .bt_random_effect_has_row_indexed_external_sd(random_term)
  }
  if(isTRUE(row_indexed)){
    source_draws <- if(is.null(row_source_evaluator)) {
      .bt_JAGS_marglik_row_indexed_external_sd_source_draws(
        random_term = random_term,
        n_rows = nrow(random_term$model_matrix),
        posterior = posterior,
        data = data,
        parameters = parameters,
        context = "Bridge context"
      )
    } else {
      row_source_evaluator(posterior, data = data, parameters = parameters)
    }
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

  sd_values <- if(is.null(sd_evaluator)) {
    .bt_JAGS_marglik_random_effect_sd_values(
      samples = samples,
      random_term = random_term,
      prior_list = prior_list
    )
  } else {
    sd_evaluator$values(samples, parameters = parameters)
  }
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

.bt_JAGS_bridge_context_random_block_nodes <- function(
    samples,
    random_term,
    prior_list,
    data,
    parameters,
    posterior = NULL,
    row_indexed = NULL,
    sd_evaluator = NULL,
    row_source_evaluator = NULL,
    sd_cache = NULL,
    sd_share_key = NULL){

  if(is.null(posterior)){
    posterior <- .bt_JAGS_marglik_random_effect_posterior_row(samples)
  }

  if(is.null(row_indexed)){
    row_indexed <- .bt_random_effect_has_row_indexed_external_sd(random_term)
  }
  if(isTRUE(row_indexed)){
    source_draws <- if(is.null(row_source_evaluator)) {
      .bt_JAGS_marglik_row_indexed_external_sd_source_draws(
        random_term = random_term,
        n_rows = nrow(random_term$model_matrix),
        posterior = posterior,
        data = data,
        parameters = parameters,
        context = "Bridge context"
      )
    } else {
      row_source_evaluator(posterior, data = data, parameters = parameters)
    }
    nodes <- unname(source_draws[1L, ])
    names(nodes) <- colnames(source_draws)
  }else{
    sd_values <- if(is.null(sd_evaluator)) {
      .bt_JAGS_marglik_random_effect_sd_values(
        samples = samples,
        random_term = random_term,
        prior_list = prior_list
      )
    } else {
      sd_evaluator$values(samples, parameters = parameters)
    }
    # The marginal block plan needs this state's SD vector in another shape.
    # Leave it where that plan can take it instead of evaluating the same
    # allocation chain a second time, but only while both would have read the
    # same values.
    if(!is.null(sd_evaluator) && is.environment(sd_cache) &&
       is.character(sd_share_key) && length(sd_share_key) == 1L &&
       .bt_JAGS_bridge_random_sd_sources_agree(
         bindings = sd_evaluator$bindings,
         parameters = parameters,
         posterior = posterior
       )){
      assign(sd_share_key, sd_values, envir = sd_cache)
    }
    nodes <- .bt_JAGS_bridge_context_random_sd_nodes(
      random_term = random_term,
      sd_values = sd_values
    )
  }

  # Normalized allocation parameters are already flattened by
  # .bt_JAGS_bridge_context_nodes(). Replaying them here would reconstruct and
  # merge the same exact nodes a second time.
  list(nodes = nodes, allocations = list())
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
                                                      formula_prior_parameters,
                                                      parameter_names = NULL){

  if(is.null(parameter_names)){
    parameter_names <- .bt_JAGS_marglik_random_effect_allocation_parameters(
      random_term
    )
  }
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
      prior_list = prior_list,
      prefer_weights = FALSE
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
