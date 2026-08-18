.bt_random_sd_binding_context <- function(random_effects, prior_random,
                                          parameter){

  empty_context <- list(
    prior_list = list(),
    syntax = character(),
    by_block = list(),
    allocations = list(),
    add_parameters = character()
  )
  if(is.null(prior_random) || is.null(prior_random$allocation)){
    return(empty_context)
  }

  .bt_check_prior_random(prior_random)
  check_char(parameter, "parameter", allow_NA = FALSE)

  if(length(random_effects) == 0L){
    stop(
      "Variance allocation priors require formula random-effect terms.",
      call. = FALSE
    )
  }

  allocations <- .bt_random_allocation_list(prior_random$allocation)
  block_names <- vapply(random_effects, function(term) term$block_name, character(1))

  prior_list <- list()
  syntax <- character()
  add_parameters <- character()
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
    label <- .bt_random_variance_allocation_resolve_label(allocation)
    if(label %in% allocation_labels){
      stop("Variance allocation labels must be unique.", call. = FALSE)
    }
    allocation_labels[allocation_i] <- label
    allocation_terms[[allocation_i]] <- terms
    allocation_component_labels[[allocation_i]] <- .bt_random_variance_allocation_component_labels(terms)
    .bt_check_random_allocation_inclusion(
      allocation$inclusion,
      component_labels = allocation_component_labels[[allocation_i]],
      allocation_label = label
    )
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
    display_name <- allocation$display_name
    component_names <- allocation$component_names
    target <- .bt_random_variance_allocation_target(allocation)
    scale <- .bt_random_variance_allocation_scale(allocation)
    if(identical(target, "sd_component") && !is.null(allocation$inclusion)){
      stop(
        "Variance allocation inclusion currently supports target = \"block\".",
        call. = FALSE
      )
    }

    allocation_names <- .bt_random_variance_allocation_names(parameter, label)
    if(is.null(allocation$parent)){
      source <- .bt_random_variance_allocation_root_source(
        allocation = allocation,
        allocation_names = allocation_names,
        label = label,
        terms = terms
      )
      source_name <- .bt_random_variance_allocation_source_jags_expression(
        source,
        row_index = "i"
      )
      source_base_name <- source$name
      source_factors <- list()
      if(isTRUE(source$owned)){
        prior_list[[allocation_names$scale_suffix]] <- source$prior
      }
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
      source <- parent_info$source
    }

    if(identical(target, "block")){
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

      inclusion_info <- .bt_random_variance_allocation_inclusion_info(
        allocation = allocation,
        parameter = parameter,
        label = label,
        component_labels = component_labels
      )
      if(length(inclusion_info) > 0L){
        for(component_label in names(inclusion_info)){
          inclusion <- inclusion_info[[component_label]]
          prior_list[[inclusion$prob_suffix]] <- inclusion$prior
          syntax <- c(
            syntax,
            paste0(inclusion$indicator_name, " ~ dbern(", inclusion$prob_name, ")")
          )
          add_parameters <- c(add_parameters, inclusion$indicator_name)
        }
      }

      component_meta <- list()
      for(term_i in seq_along(terms)){
        component_key <- paste(label, component_labels[term_i], sep = "::")
        inclusion_name <- NULL
        component_inclusion <- inclusion_info[[component_labels[term_i]]]
        if(!is.null(component_inclusion)){
          inclusion_name <- component_inclusion$indicator_name
        }
        expression <- .bt_random_variance_allocation_expression(
          source_name = source_name,
          weight_name = allocation_names$weight_name,
          index = term_i,
          scale = scale,
          n_targets = length(terms),
          inclusion_name = inclusion_name
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
        row_indexed_source <- .bt_random_variance_allocation_source_is_row(source)
        if(component_key %in% consumed_components && !row_indexed_source){
          syntax <- c(syntax, paste0(node_name, " = ", expression))
        }
        factor <- .bt_random_variance_allocation_factor(
          weight_name = allocation_names$weight_name,
          index = term_i,
          scale = scale,
          n_targets = length(terms),
          inclusion_name = inclusion_name
        )
        component_meta[[component_labels[term_i]]] <- list(
          label = component_labels[term_i],
          term = terms[term_i],
          node_name = if(component_key %in% consumed_components && !row_indexed_source) node_name else expression,
          expression = expression,
          base_name = source_base_name,
          factors = c(source_factors, list(factor)),
          source = source,
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
          allocation_record <- list(
            label = label,
            display_name = display_name,
            component_names = component_names,
            terms = terms,
            index = term_i,
            target = "block",
            scale = scale,
            parent = allocation$parent,
            source_node = source_name,
            source = source,
            parent_factors = source_factors,
            factors = c(source_factors, list(factor)),
            n_targets = length(terms),
            scale_name = source$scale_name,
            weight_name = allocation_names$weight_name,
            scale_suffix = source$scale_suffix,
            weight_suffix = allocation_names$weight_suffix,
            inclusion = inclusion_info
          )
          by_block[[terms[term_i]]] <- .bt_random_sd_binding(
            source = source,
            application = "block",
            factors = allocation_record$factors,
            true_allocation = TRUE,
            allocations = list(allocation_record)
          )
          used_blocks <- c(used_blocks, terms[term_i])
        }
      }

      allocation_meta[[label]] <- list(
        label = label,
        display_name = display_name,
        component_names = component_names,
        terms = terms,
        component_labels = component_labels,
        components = component_meta,
        target = "block",
        scale = scale,
        parent = allocation$parent,
        source_node = source_name,
        source = source,
        parent_factors = source_factors,
        scale_name = source$scale_name,
        weight_name = allocation_names$weight_name,
        scale_suffix = source$scale_suffix,
        weight_suffix = allocation_names$weight_suffix,
        inclusion = inclusion_info
      )
    }else if(identical(target, "sd_component")){
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
      allocation_record <- list(
        label = label,
        display_name = display_name,
        component_names = component_names,
        terms = terms,
        index = NA_integer_,
        target = "sd_component",
        scale = scale,
        parent = allocation$parent,
        source_node = source_name,
        source = source,
        parent_factors = source_factors,
        weights = allocation$weights,
        weight_name = allocation_names$weight_name,
        weight_suffix = allocation_names$weight_suffix,
        inclusion = list()
      )
      by_block[[block]] <- .bt_random_sd_binding(
        source = source,
        application = "column",
        factors = source_factors,
        true_allocation = TRUE,
        allocations = list(allocation_record)
      )
      allocation_meta[[label]] <- list(
        label = label,
        display_name = display_name,
        component_names = component_names,
        terms = terms,
        component_labels = character(),
        components = list(),
        target = "sd_component",
        scale = scale,
        parent = allocation$parent,
        source_node = source_name,
        source = source,
        parent_factors = source_factors,
        weights = allocation$weights,
        allocation_record = allocation_record,
        scale_name = source$scale_name,
        weight_name = allocation_names$weight_name,
        scale_suffix = source$scale_suffix,
        weight_suffix = allocation_names$weight_suffix,
        inclusion = list()
      )
      used_blocks <- c(used_blocks, block)
    }

  }

  list(
    prior_list = prior_list,
    syntax = syntax,
    by_block = by_block,
    allocations = allocation_meta,
    add_parameters = unique(add_parameters)
  )
}

.bt_random_variance_allocation_inclusion_info <- function(allocation,
                                                          parameter,
                                                          label,
                                                          component_labels){

  if(is.null(allocation$inclusion)){
    return(list())
  }

  inclusion_info <- list()
  for(component_label in names(allocation$inclusion)){
    inclusion_names <- .bt_random_variance_allocation_inclusion_names(
      parameter = parameter,
      label = label,
      component_label = component_label
    )
    inclusion_prior <- allocation$inclusion[[component_label]]
    attr(inclusion_prior, "random_allocation") <- label
    attr(inclusion_prior, "random_allocation_inclusion") <- component_label
    attr(inclusion_prior, "random_component") <- component_label
    inclusion_info[[component_label]] <- list(
      component = component_label,
      index = match(component_label, component_labels),
      prob_suffix = inclusion_names$prob_suffix,
      prob_name = inclusion_names$prob_name,
      indicator_name = inclusion_names$indicator_name,
      prior = inclusion_prior
    )
  }

  inclusion_info
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
    if(!identical(.bt_random_variance_allocation_target(allocation), "sd_component")){
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
  .bt_validate_random_effect_reserved_name(
    x,
    context = "variance allocation labels"
  )
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

.bt_random_sd_binding_for_block <- function(binding_context, block_name){

  if(is.null(binding_context) || length(binding_context$by_block) == 0L){
    return(NULL)
  }

  binding_context$by_block[[block_name]]
}

.bt_random_variance_allocation_source_is_row <- function(source){

  !is.null(source) &&
    identical(source$kind, "external") &&
    identical(source$shape, "row")
}

.bt_random_variance_allocation_source_jags_expression <- function(source,
                                                                  row_index = NULL){

  if(!is.null(source) && identical(source$kind, "external")){
    return(.bt_random_sd_source_expression(source, row_index = row_index))
  }

  source$name
}

.bt_random_effect_has_row_indexed_external_sd <- function(random_term){

  .bt_random_sd_binding_has_row_external_source(random_term$sd_binding)
}

.bt_random_effect_external_sd_source_label <- function(random_term){

  .bt_random_sd_binding_external_source_label(random_term$sd_binding)
}
