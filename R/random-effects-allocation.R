# Internal random-effect variance-allocation helpers.

.bt_random_variance_allocation_target <- function(allocation){

  target <- allocation$target
  check_char(target, "allocation$target",
             allow_values = c("block", "sd_component"),
             allow_NA = FALSE)

  target
}

.bt_random_variance_allocation_scale <- function(allocation){

  scale <- allocation$scale
  check_char(scale, "allocation$scale",
             allow_values = c("total_variance", "mean_variance"),
             allow_NA = FALSE)

  scale
}

.bt_random_variance_allocation_terms <- function(allocation, block_names,
                                                 used_blocks, n_allocations){

  terms <- allocation$terms
  target <- .bt_random_variance_allocation_target(allocation)
  if(is.null(terms)){
    if(!identical(target, "block") || n_allocations != 1L ||
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

  if(identical(target, "block") && length(terms) < 2L){
    stop(
      "Block variance allocation requires at least two resolved random-effect blocks. ",
      "Use random_block(sd_source = ...) for a direct one-block SD source.",
      call. = FALSE
    )
  }

  if(identical(target, "sd_component") && length(terms) != 1L){
    stop("'target = \"sd_component\"' requires exactly one random-effect block in 'terms'.", call. = FALSE)
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
    .bt_validate_random_effect_reserved_name(
      labels,
      context = "variance allocation component labels"
    )
    bad <- !grepl("^[A-Za-z][A-Za-z0-9_]*$", labels)
    if(any(bad)){
      stop("Variance allocation term labels must start with a letter and contain only letters, numbers, and underscores.", call. = FALSE)
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

  allocation_prior <- allocation$weights
  if(is.null(allocation_prior)){
    allocation_prior <- prior("dirichlet", list(alpha = rep(1, K)))
  }
  .bt_check_random_allocation_prior(allocation_prior)
  if(K < 2L){
    stop(
      "Variance allocation requires at least two resolved targets.",
      call. = FALSE
    )
  }
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
      if(!is.null(override$sd) || !is.null(override$sd_source) ||
         !is.null(override_cov_sd) || !is.null(override$terms)){
        stop(
          "Random-effect block '", term,
          "' cannot supply block-specific SD, SD source, or term SD overrides while it is controlled by a variance allocation prior.",
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

  total_suffix <- paste0("_xRE_ALLOCx_", label, "__total_sd")
  weight_suffix <- paste0("_xRE_ALLOCx_", label, "__weight")

  list(
    total_suffix = total_suffix,
    weight_suffix = weight_suffix,
    total_name = paste0(parameter, "_", total_suffix),
    weight_name = paste0(parameter, "_", weight_suffix)
  )
}

.bt_random_variance_allocation_component_name <- function(parameter, label,
                                                          component_label){

  paste0(parameter, "__xRE_ALLOCx_", label, "__component_", component_label, "_sd")
}

.bt_random_variance_allocation_inclusion_names <- function(parameter, label,
                                                           component_label){

  stem <- paste0("_xRE_ALLOCx_", label, "__include_", component_label)
  list(
    prob_suffix = paste0(stem, "_prob"),
    prob_name = paste0(parameter, "_", stem, "_prob"),
    indicator_name = paste0(parameter, "_", stem, "_indicator")
  )
}

.bt_random_variance_allocation_inclusion_indicator_names <- function(formula_design){

  if(inherits(formula_design, "BayesTools_formula_design")){
    formula_design <- list(formula_design)
  }
  if(!is.list(formula_design)){
    return(character())
  }

  indicators <- character()
  collect_allocation <- function(allocation){
    if(!is.list(allocation) || !is.list(allocation$inclusion) ||
       length(allocation$inclusion) == 0L){
      return(invisible(NULL))
    }
    for(inclusion in allocation$inclusion){
      indicator_name <- inclusion$indicator_name
      if(is.character(indicator_name) && length(indicator_name) == 1L &&
         !is.na(indicator_name) && nzchar(indicator_name)){
        indicators <<- c(indicators, indicator_name)
      }
    }
    invisible(NULL)
  }

  for(design in formula_design){
    if(!inherits(design, "BayesTools_formula_design")){
      next
    }
    if(is.list(design$random_allocations)){
      for(allocation in design$random_allocations){
        collect_allocation(allocation)
      }
    }
    for(random_term in .bt_formula_design_random_effects(design)){
      binding <- random_term$sd_binding
      if(is.null(binding) || !is.list(binding$allocations)){
        next
      }
      for(allocation in binding$allocations){
        collect_allocation(allocation)
      }
    }
  }

  unique(indicators)
}

.bt_random_variance_allocation_factor <- function(weight_name, index, scale,
                                                  n_targets,
                                                  inclusion_name = NULL){

  list(
    weight_name = weight_name,
    index = index,
    scale = scale,
    n_targets = n_targets,
    inclusion_name = inclusion_name
  )
}

.bt_check_random_variance_allocation_factor <- function(factor){

  if(!is.list(factor)){
    stop("Random-effect allocation factor metadata are missing canonical fields.", call. = FALSE)
  }
  if(!is.character(factor$weight_name) || length(factor$weight_name) != 1L ||
     is.na(factor$weight_name) || !nzchar(factor$weight_name)){
    stop("Random-effect allocation factor metadata are missing canonical 'weight_name'.", call. = FALSE)
  }
  if(!is.numeric(factor$index) || length(factor$index) != 1L ||
     is.na(factor$index) || factor$index != as.integer(factor$index) ||
     factor$index < 1L){
    stop("Random-effect allocation factor metadata are missing canonical 'index'.", call. = FALSE)
  }
  check_char(factor$scale, "factor$scale",
             allow_values = c("total_variance", "mean_variance"),
             allow_NA = FALSE)
  if(!is.numeric(factor$n_targets) || length(factor$n_targets) != 1L ||
     is.na(factor$n_targets) || factor$n_targets != as.integer(factor$n_targets) ||
     factor$n_targets < 2L){
    stop("Random-effect allocation factor metadata are missing canonical 'n_targets'.", call. = FALSE)
  }
  if(factor$index > factor$n_targets){
    stop("Random-effect allocation factor metadata reference a coordinate outside 'n_targets'.", call. = FALSE)
  }
  if(!is.null(factor$inclusion_name)){
    check_char(factor$inclusion_name, "factor$inclusion_name",
               allow_NA = FALSE)
  }

  invisible(TRUE)
}

.bt_check_random_variance_allocation_factor_chain <- function(factors,
                                                              label = "factors"){

  if(!is.list(factors)){
    stop("Random-effect SD binding metadata are missing '", label, "'.", call. = FALSE)
  }
  for(factor_i in seq_along(factors)){
    .bt_check_random_variance_allocation_factor(factors[[factor_i]])
  }

  invisible(TRUE)
}

.bt_random_variance_allocation_multiplier_expression <- function(weight_name,
                                                                 index, scale,
                                                                 n_targets){

  multiplier <- if(identical(scale, "mean_variance")){
    paste0(n_targets, " * ", weight_name, "[", index, "]")
  }else{
    paste0(weight_name, "[", index, "]")
  }

  paste0("sqrt(", multiplier, ")")
}

.bt_random_variance_allocation_factor_expression <- function(factor){

  multiplier <- .bt_random_variance_allocation_multiplier_expression(
    weight_name = factor$weight_name,
    index = factor$index,
    scale = factor$scale,
    n_targets = factor$n_targets
  )

  if(!is.null(factor$inclusion_name)){
    multiplier <- paste0(factor$inclusion_name, " * ", multiplier)
  }

  multiplier
}

.bt_random_variance_allocation_factors_expression <- function(factors){

  if(length(factors) == 0L){
    return("1")
  }

  paste(
    vapply(factors, .bt_random_variance_allocation_factor_expression, character(1)),
    collapse = " * "
  )
}

.bt_random_variance_allocation_expression <- function(source_name, weight_name,
                                                      index, scale, n_targets,
                                                      inclusion_name = NULL){

  multiplier <- .bt_random_variance_allocation_factor_expression(
    .bt_random_variance_allocation_factor(
      weight_name = weight_name,
      index = index,
      scale = scale,
      n_targets = n_targets,
      inclusion_name = inclusion_name
    )
  )

  paste0(source_name, " * ", multiplier)
}

.bt_random_variance_allocation_root_source <- function(allocation,
                                                       allocation_names,
                                                       label,
                                                       terms){

  if(!is.null(allocation$sd_source)){
    .bt_check_random_sd_source(allocation$sd_source)
    source <- allocation$sd_source
    source$total_name <- source$name
    source$total_suffix <- NULL
    return(source)
  }

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

  list(
    kind = "prior",
    name = allocation_names$total_name,
    shape = "scalar",
    owned = TRUE,
    total_name = allocation_names$total_name,
    total_suffix = allocation_names$total_suffix,
    prior = total_prior
  )
}

