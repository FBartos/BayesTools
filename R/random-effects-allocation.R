# Internal random-effect variance-allocation helpers.

.bt_random_variance_allocation_components <- function(allocation){

  components <- allocation$components
  if(is.null(components)){
    components <- "block"
  }

  components
}

.bt_random_variance_allocation_scale <- function(allocation){

  scale <- allocation$scale
  if(is.null(scale)){
    scale <- "total_variance"
  }

  scale
}

.bt_random_variance_allocation_terms <- function(allocation, block_names,
                                                 used_blocks, n_allocations){

  terms <- allocation$terms
  components <- .bt_random_variance_allocation_components(allocation)
  if(is.null(terms)){
    if(!identical(components, "block") || n_allocations != 1L ||
       !is.null(allocation$parent)){
      stop(
        "A variance allocation prior without explicit 'terms' is supported only for a single root block allocation.",
        call. = FALSE
      )
    }
    terms <- setdiff(block_names, used_blocks)
    if(length(terms) < 2L){
      stop(
        "A variance allocation prior without explicit 'terms' requires at least two unallocated random-effect blocks.",
        call. = FALSE
      )
    }
  }

  if(identical(components, "sd") && length(terms) != 1L){
    stop("'components = \"sd\"' requires exactly one random-effect block in 'terms'.", call. = FALSE)
  }

  terms
}

.bt_random_variance_allocation_component_labels <- function(terms){

  labels <- names(terms)
  if(!is.null(labels) && any(nzchar(labels))){
    if(any(!nzchar(labels))){
      stop("Variance allocation 'terms' must be either all named or all unnamed.", call. = FALSE)
    }
    if(anyDuplicated(labels)){
      stop("Variance allocation term labels must be unique.", call. = FALSE)
    }
    return(labels)
  }

  labels <- vapply(terms, .bt_random_variance_allocation_label, character(1))
  if(anyDuplicated(labels)){
    stop("Variance allocation component labels must be unique after sanitization.", call. = FALSE)
  }

  labels
}

.bt_random_variance_allocation_prior <- function(allocation, K){

  allocation_prior <- allocation$allocation
  if(is.null(allocation_prior)){
    allocation_prior <- prior("dirichlet", list(alpha = rep(1, K)))
  }
  .bt_check_random_allocation_prior(allocation_prior)
  if(K != allocation_prior$parameters[["K"]]){
    stop(
      "The Dirichlet allocation dimension must match the number of targeted random-effect terms.",
      call. = FALSE
    )
  }

  allocation_prior
}

.bt_validate_random_variance_allocation_block_overrides <- function(terms,
                                                                    prior_random){

  for(term in terms){
    override <- prior_random$blocks[[term]]
    if(!is.null(override)){
      override_cov_sd <- if(!is.null(override$covariance)) override$covariance$sd else NULL
      if(!is.null(override$sd) || !is.null(override_cov_sd) || !is.null(override$terms)){
        stop(
          "Random-effect block '", term,
          "' cannot supply block-specific SD or term SD overrides while it is controlled by a variance allocation prior.",
          call. = FALSE
        )
      }
    }
  }

  invisible(TRUE)
}

.bt_random_variance_allocation_resolve_label <- function(allocation,
                                                         allocation_i,
                                                         allocations,
                                                         terms){

  label <- allocation$name
  allocation_names <- names(allocations)
  if(is.null(label) && !is.null(allocation_names) &&
     nzchar(allocation_names[allocation_i])){
    label <- allocation_names[allocation_i]
  }
  if(is.null(label)){
    label <- if(length(allocations) == 1L){
      "allocation"
    }else{
      paste(terms, collapse = "_")
    }
  }

  .bt_random_variance_allocation_label(label)
}

.bt_random_variance_allocation_names <- function(parameter, label){

  total_suffix <- paste0("_xRE_ALLOCx_", label, "_total_sd")
  weight_suffix <- paste0("_xRE_ALLOCx_", label, "_weight")

  list(
    total_suffix = total_suffix,
    weight_suffix = weight_suffix,
    total_name = paste0(parameter, "_", total_suffix),
    weight_name = paste0(parameter, "_", weight_suffix)
  )
}

.bt_random_variance_allocation_component_name <- function(parameter, label,
                                                          component_label){

  paste0(parameter, "__xRE_ALLOCx_", label, "_", component_label, "_sd")
}

.bt_random_variance_allocation_factor <- function(weight_name, index, scale,
                                                  n_targets){

  list(
    weight_name = weight_name,
    index = index,
    scale = scale,
    n_targets = n_targets
  )
}

.bt_random_variance_allocation_expression <- function(source_name, weight_name,
                                                      index, scale, n_targets){

  multiplier <- if(identical(scale, "mean_variance")){
    paste0(n_targets, " * ", weight_name, "[", index, "]")
  }else{
    paste0(weight_name, "[", index, "]")
  }

  paste0(source_name, " * sqrt(", multiplier, ")")
}

.bt_random_variance_allocation_context <- function(random_effects, prior_random,
                                                   parameter){

  empty_context <- list(
    prior_list = list(),
    syntax = character(),
    by_block = list(),
    allocations = list()
  )
  if(length(random_effects) == 0L || is.null(prior_random) ||
     is.null(prior_random$allocation)){
    return(empty_context)
  }

  .bt_check_prior_random(prior_random)
  check_char(parameter, "parameter", allow_NA = FALSE)

  allocations <- .bt_random_allocation_list(prior_random$allocation)
  block_names <- vapply(random_effects, function(term) term$block_name, character(1))

  prior_list <- list()
  syntax <- character()
  by_block <- list()
  allocation_meta <- list()
  used_blocks <- character()
  allocation_labels <- character(length(allocations))
  allocation_terms <- vector("list", length(allocations))
  allocation_component_labels <- vector("list", length(allocations))

  for(allocation_i in seq_along(allocations)){
    allocation <- allocations[[allocation_i]]
    terms <- .bt_random_variance_allocation_terms(
      allocation = allocation,
      block_names = block_names,
      used_blocks = character(),
      n_allocations = length(allocations)
    )
    label <- .bt_random_variance_allocation_resolve_label(
      allocation = allocation,
      allocation_i = allocation_i,
      allocations = allocations,
      terms = terms
    )
    if(label %in% allocation_labels){
      stop("Variance allocation labels must be unique.", call. = FALSE)
    }
    allocation_labels[allocation_i] <- label
    allocation_terms[[allocation_i]] <- terms
    allocation_component_labels[[allocation_i]] <- .bt_random_variance_allocation_component_labels(terms)
  }

  consumed_components <- character()
  for(allocation in allocations){
    if(!is.null(allocation$parent)){
      consumed_components <- c(
        consumed_components,
        paste(allocation$parent$allocation, allocation$parent$component, sep = "::")
      )
    }
  }
  if(anyDuplicated(consumed_components)){
    stop("A variance allocation parent component can be consumed by only one child allocation.", call. = FALSE)
  }

  for(allocation_i in seq_along(allocations)){
    allocation <- allocations[[allocation_i]]
    terms <- allocation_terms[[allocation_i]]
    component_labels <- allocation_component_labels[[allocation_i]]
    label <- allocation_labels[[allocation_i]]
    components <- .bt_random_variance_allocation_components(allocation)
    scale <- .bt_random_variance_allocation_scale(allocation)

    allocation_names <- .bt_random_variance_allocation_names(parameter, label)
    if(is.null(allocation$parent)){
      source_name <- allocation_names$total_name
      source_base_name <- allocation_names$total_name
      source_factors <- list()

      total_prior <- .bt_random_effect_force_nonnegative_prior(
        prior = allocation$sd,
        name = paste0("variance allocation '", label, "' total SD")
      )
      .bt_random_effect_check_scalar_sd_prior(
        total_prior,
        paste0("variance allocation '", label, "' total SD")
      )
      total_prior <- .bt_random_effect_set_total_sd_metadata(
        total_prior,
        allocation = label,
        terms = terms
      )

      prior_list[[allocation_names$total_suffix]] <- total_prior
    }else{
      parent_label <- allocation$parent$allocation
      parent_component <- allocation$parent$component
      if(!parent_label %in% names(allocation_meta)){
        stop(
          "Parent variance allocation '", parent_label,
          "' must be defined before child allocation '", label, "'.",
          call. = FALSE
        )
      }
      parent_info <- allocation_meta[[parent_label]]$components[[parent_component]]
      if(is.null(parent_info)){
        stop(
          "Parent variance allocation '", parent_label,
          "' does not contain component '", parent_component, "'.",
          call. = FALSE
        )
      }
      source_name <- parent_info$node_name
      source_base_name <- parent_info$base_name
      source_factors <- parent_info$factors
    }

    if(identical(components, "block")){
      unknown_terms <- setdiff(terms, block_names)
      unknown_unconsumed <- character()
      if(length(unknown_terms) > 0L){
        unknown_unconsumed <- unknown_terms[
          !(paste(label, component_labels[match(unknown_terms, terms)], sep = "::") %in% consumed_components)
        ]
      }
      if(length(unknown_unconsumed) > 0L){
        stop(
          "Variance allocation targets unknown random-effect block(s): ",
          paste(unknown_unconsumed, collapse = ", "),
          ". Unknown targets are allowed only when consumed by a child allocation.",
          call. = FALSE
        )
      }

      allocation_prior <- .bt_random_variance_allocation_prior(allocation, length(terms))
      allocation_prior <- .bt_random_effect_set_allocation_metadata(
        allocation_prior,
        allocation = label,
        terms = terms,
        parent = allocation$parent
      )
      prior_list[[allocation_names$weight_suffix]] <- allocation_prior

      component_meta <- list()
      for(term_i in seq_along(terms)){
        component_key <- paste(label, component_labels[term_i], sep = "::")
        expression <- .bt_random_variance_allocation_expression(
          source_name = source_name,
          weight_name = allocation_names$weight_name,
          index = term_i,
          scale = scale,
          n_targets = length(terms)
        )
        node_name <- .bt_random_variance_allocation_component_name(
          parameter = parameter,
          label = label,
          component_label = component_labels[term_i]
        )
        if(component_key %in% consumed_components && terms[term_i] %in% block_names &&
           !.bt_random_variance_allocation_component_has_sd_child(
             allocations = allocations,
             parent_label = label,
             component_label = component_labels[term_i],
             block = terms[term_i]
           )){
          stop(
            "Variance allocation component '", component_labels[term_i],
            "' in allocation '", label,
            "' is consumed by a child allocation but also names a random-effect block. ",
            "Use a symbolic component label that is not a block name, or allocate the block directly.",
            call. = FALSE
          )
        }
        if(component_key %in% consumed_components){
          syntax <- c(syntax, paste0(node_name, " = ", expression))
        }
        factor <- .bt_random_variance_allocation_factor(
          weight_name = allocation_names$weight_name,
          index = term_i,
          scale = scale,
          n_targets = length(terms)
        )
        component_meta[[component_labels[term_i]]] <- list(
          label = component_labels[term_i],
          term = terms[term_i],
          node_name = if(component_key %in% consumed_components) node_name else expression,
          expression = expression,
          base_name = source_base_name,
          factors = c(source_factors, list(factor)),
          index = term_i
        )
        if(!(component_key %in% consumed_components)){
          if(terms[term_i] %in% used_blocks){
            stop(
              "Random-effect block(s) cannot appear in more than one variance allocation prior: ",
              terms[term_i],
              ".",
              call. = FALSE
            )
          }
          .bt_validate_random_variance_allocation_block_overrides(terms[term_i], prior_random)
          by_block[[terms[term_i]]] <- list(
            label = label,
            terms = terms,
            index = term_i,
            components = "block",
            scale = scale,
            source_expression = expression,
            source_name = source_base_name,
            factors = c(source_factors, list(factor)),
            n_targets = length(terms),
            total_name = source_base_name,
            weight_name = allocation_names$weight_name,
            total_suffix = allocation_names$total_suffix,
            weight_suffix = allocation_names$weight_suffix
          )
          used_blocks <- c(used_blocks, terms[term_i])
        }
      }

      allocation_meta[[label]] <- list(
        label = label,
        terms = terms,
        component_labels = component_labels,
        components = component_meta,
        target = "block",
        scale = scale,
        total_name = allocation_names$total_name,
        weight_name = allocation_names$weight_name,
        total_suffix = allocation_names$total_suffix,
        weight_suffix = allocation_names$weight_suffix
      )
    }else if(identical(components, "sd")){
      block <- terms[[1L]]
      if(!block %in% block_names){
        stop(
          "Variance allocation targets unknown random-effect block(s): ",
          block,
          ".",
          call. = FALSE
        )
      }
      if(block %in% used_blocks){
        stop(
          "Random-effect block(s) cannot appear in more than one variance allocation prior: ",
          block,
          ".",
          call. = FALSE
        )
      }
      .bt_validate_random_variance_allocation_block_overrides(block, prior_random)
      by_block[[block]] <- list(
        label = label,
        terms = terms,
        index = NA_integer_,
        components = "sd",
        scale = scale,
        source_node = source_name,
        source_name = source_base_name,
        parent_factors = source_factors,
        allocation = allocation$allocation,
        weight_name = allocation_names$weight_name,
        weight_suffix = allocation_names$weight_suffix
      )
      allocation_meta[[label]] <- list(
        label = label,
        terms = terms,
        component_labels = character(),
        components = list(),
        target = "sd",
        scale = scale,
        total_name = allocation_names$total_name,
        weight_name = allocation_names$weight_name,
        total_suffix = allocation_names$total_suffix,
        weight_suffix = allocation_names$weight_suffix
      )
      used_blocks <- c(used_blocks, block)
    }

  }

  list(
    prior_list = prior_list,
    syntax = syntax,
    by_block = by_block,
    allocations = allocation_meta
  )
}

.bt_random_variance_allocation_component_has_sd_child <- function(allocations,
                                                                 parent_label,
                                                                 component_label,
                                                                 block){

  for(allocation in allocations){
    if(is.null(allocation$parent)){
      next
    }
    if(!identical(allocation$parent$allocation, parent_label) ||
       !identical(allocation$parent$component, component_label)){
      next
    }
    if(!identical(.bt_random_variance_allocation_components(allocation), "sd")){
      next
    }
    if(length(allocation$terms) == 1L &&
       identical(unname(allocation$terms), unname(block))){
      return(TRUE)
    }
  }

  FALSE
}

.bt_random_variance_allocation_label <- function(x){

  check_char(x, "name", allow_NA = FALSE)
  x <- gsub("[^A-Za-z0-9_]", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  if(!nzchar(x)){
    stop("Variance allocation labels must contain at least one letter, digit, or underscore.", call. = FALSE)
  }
  if(!grepl("^[A-Za-z]", x)){
    x <- paste0("allocation_", x)
  }

  x
}

.bt_random_variance_allocation_for_block <- function(allocation_context,
                                                     block_name){

  if(is.null(allocation_context) || length(allocation_context$by_block) == 0L){
    return(NULL)
  }

  allocation_context$by_block[[block_name]]
}
