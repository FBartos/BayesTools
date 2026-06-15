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

.bt_random_variance_allocation_factor <- function(weight_name, index, scale,
                                                  n_targets){

  list(
    weight_name = weight_name,
    index = index,
    scale = scale,
    n_targets = n_targets
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

.bt_random_variance_allocation_expression <- function(source_name, weight_name,
                                                      index, scale, n_targets){

  multiplier <- .bt_random_variance_allocation_multiplier_expression(
    weight_name = weight_name,
    index = index,
    scale = scale,
    n_targets = n_targets
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

.bt_random_sd_binding <- function(source = NULL, sources_by_column = list(),
                                  application = c("block", "column"),
                                  factors = list(), factors_by_column = list(),
                                  true_allocation = FALSE,
                                  allocations = list(),
                                  sd_component_names = NULL,
                                  sd_component_terms = NULL,
                                  sd_component_index_by_column = NULL){

  application <- match.arg(application)
  if(!is.null(source)){
    .bt_check_random_sd_binding_source(source)
  }
  if(!is.list(sources_by_column)){
    stop("'sources_by_column' must be a list.", call. = FALSE)
  }
  for(source_i in seq_along(sources_by_column)){
    .bt_check_random_sd_binding_source(sources_by_column[[source_i]])
  }
  if(!is.list(factors)){
    stop("'factors' must be a list.", call. = FALSE)
  }
  if(!is.list(factors_by_column)){
    stop("'factors_by_column' must be a list.", call. = FALSE)
  }
  .bt_check_random_variance_allocation_factor_chain(factors)
  for(column_i in seq_along(factors_by_column)){
    .bt_check_random_variance_allocation_factor_chain(
      factors_by_column[[column_i]],
      label = "factors_by_column"
    )
  }
  if(!is.list(allocations)){
    stop("'allocations' must be a list.", call. = FALSE)
  }
  check_bool(true_allocation, "true_allocation")
  if(identical(application, "block") && is.null(source)){
    stop("Block SD bindings require a shared 'source'.", call. = FALSE)
  }
  if(identical(application, "block") && length(sources_by_column) > 0L){
    stop("Block SD bindings must not use 'sources_by_column'.", call. = FALSE)
  }
  if(identical(application, "block") && length(factors_by_column) > 0L){
    stop("Block SD bindings must not use 'factors_by_column'.", call. = FALSE)
  }
  if(identical(application, "column") && is.null(source) &&
     length(sources_by_column) == 0L){
    stop("Column SD bindings require shared or per-column sources.", call. = FALSE)
  }

  out <- list(
    source = source,
    sources_by_column = sources_by_column,
    application = application,
    factors = factors,
    factors_by_column = factors_by_column,
    true_allocation = true_allocation,
    allocations = allocations,
    sd_component_names = sd_component_names,
    sd_component_terms = sd_component_terms,
    sd_component_index_by_column = sd_component_index_by_column
  )
  class(out) <- c("random_sd_binding", "list")
  .bt_check_random_sd_binding(out)

  out
}

.bt_prior_owned_sd_source <- function(name, prior_name = NULL){

  check_char(name, "name", allow_NA = FALSE)
  if(!is.null(prior_name)){
    check_char(prior_name, "prior_name", allow_NA = FALSE)
  }

  out <- list(
    name = name,
    prior_name = prior_name,
    shape = "scalar",
    kind = "prior",
    owned = TRUE
  )
  class(out) <- c("prior_owned_sd_source", "list")

  out
}

.bt_check_random_sd_binding_source <- function(source){

  if(inherits(source, "random_sd_source")){
    .bt_check_random_sd_source(source)
    return(invisible(TRUE))
  }
  if(!inherits(source, "prior_owned_sd_source") && !is.list(source)){
    stop("Random SD binding sources must be source descriptor objects.", call. = FALSE)
  }
  if(!is.character(source$name) || length(source$name) != 1L ||
     is.na(source$name) || !nzchar(source$name)){
    stop("Random SD binding source metadata are missing 'name'.", call. = FALSE)
  }
  if(!is.character(source$shape) || length(source$shape) != 1L ||
     is.na(source$shape) || !source$shape %in% c("scalar", "row")){
    stop("Random SD binding source metadata are missing 'shape'.", call. = FALSE)
  }
  if(!is.character(source$kind) || length(source$kind) != 1L ||
     is.na(source$kind) || !source$kind %in% c("prior", "external")){
    stop("Random SD binding source metadata are missing 'kind'.", call. = FALSE)
  }
  if(!is.logical(source$owned) || length(source$owned) != 1L ||
     is.na(source$owned)){
    stop("Random SD binding source metadata are missing 'owned'.", call. = FALSE)
  }
  if(inherits(source, "prior_owned_sd_source")){
    if(!identical(source$shape, "scalar") ||
       !identical(source$kind, "prior") ||
       !identical(source$owned, TRUE)){
      stop("Prior-owned random SD binding source metadata are inconsistent.", call. = FALSE)
    }
    if(!is.null(source$prior_name)){
      check_char(source$prior_name, "source$prior_name", allow_NA = FALSE)
    }
    return(invisible(TRUE))
  }
  if(identical(source$kind, "prior") &&
     (!identical(source$shape, "scalar") || !identical(source$owned, TRUE))){
    stop("Prior random SD binding source metadata are inconsistent.", call. = FALSE)
  }
  if(identical(source$kind, "external") && !identical(source$owned, FALSE)){
    stop("External random SD binding source metadata are inconsistent.", call. = FALSE)
  }

  invisible(TRUE)
}

.bt_check_random_sd_binding <- function(binding, allow_NULL = FALSE){

  if(is.null(binding) && isTRUE(allow_NULL)){
    return(invisible(TRUE))
  }
  if(!inherits(binding, "random_sd_binding")){
    stop("Random-effect SD binding metadata are missing.", call. = FALSE)
  }
  if(!is.null(binding$source)){
    .bt_check_random_sd_binding_source(binding$source)
  }
  if(!is.list(binding$sources_by_column)){
    stop("Random-effect SD binding metadata are missing 'sources_by_column'.", call. = FALSE)
  }
  for(source_i in seq_along(binding$sources_by_column)){
    .bt_check_random_sd_binding_source(binding$sources_by_column[[source_i]])
  }
  check_char(binding$application, "binding$application",
             allow_values = c("block", "column"), allow_NA = FALSE)
  if(identical(binding$application, "block") && is.null(binding$source)){
    stop("Block SD bindings require a shared 'source'.", call. = FALSE)
  }
  if(identical(binding$application, "block") &&
     length(binding$sources_by_column) > 0L){
    stop("Block SD bindings must not use 'sources_by_column'.", call. = FALSE)
  }
  if(identical(binding$application, "column") && is.null(binding$source) &&
     length(binding$sources_by_column) == 0L){
    stop("Column SD bindings require shared or per-column sources.", call. = FALSE)
  }
  if(length(binding$sources_by_column) > 0L &&
     any(vapply(
       binding$sources_by_column,
       .bt_random_sd_binding_source_is_external_row,
       logical(1)
     ))){
    stop(
      "Row-indexed per-column SD sources are not supported. Use one shared row-shaped SD source with allocation factors.",
      call. = FALSE
    )
  }
  .bt_check_random_variance_allocation_factor_chain(binding$factors)
  if(!is.list(binding$factors_by_column)){
    stop("Random-effect SD binding metadata are missing 'factors_by_column'.", call. = FALSE)
  }
  for(column_i in seq_along(binding$factors_by_column)){
    .bt_check_random_variance_allocation_factor_chain(
      binding$factors_by_column[[column_i]],
      label = "factors_by_column"
    )
  }
  if(identical(binding$application, "block") &&
     length(binding$factors_by_column) > 0L){
    stop("Block SD bindings must not use 'factors_by_column'.", call. = FALSE)
  }
  if(!is.logical(binding$true_allocation) || length(binding$true_allocation) != 1L ||
     is.na(binding$true_allocation)){
    stop("Random-effect SD binding metadata are missing 'true_allocation'.", call. = FALSE)
  }
  if(!is.list(binding$allocations)){
    stop("Random-effect SD binding metadata are missing 'allocations'.", call. = FALSE)
  }
  if(isTRUE(binding$true_allocation)){
    if(length(binding$allocations) != 1L){
      stop(
        "True allocation SD bindings require exactly one allocation record.",
        call. = FALSE
      )
    }
    .bt_check_random_sd_binding_allocation_record(
      binding$allocations[[1L]],
      binding = binding
    )
  }else if(length(binding$allocations) != 0L){
    stop(
      "Non-allocation SD bindings must not contain allocation records.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

.bt_check_random_sd_binding_allocation_record <- function(allocation,
                                                         binding = NULL){

  if(!is.list(allocation)){
    stop(
      "Random-effect SD binding metadata are missing canonical allocation record.",
      call. = FALSE
    )
  }
  if(!is.character(allocation$target) || length(allocation$target) != 1L ||
     is.na(allocation$target) ||
     !allocation$target %in% c("block", "sd_component")){
    stop(
      "Random-effect SD binding metadata are missing canonical 'allocation$target'.",
      call. = FALSE
    )
  }
  if(!is.character(allocation$weight_name) ||
     length(allocation$weight_name) != 1L ||
     is.na(allocation$weight_name) || !nzchar(allocation$weight_name)){
    stop(
      "Random-effect SD binding metadata are missing canonical 'allocation$weight_name'.",
      call. = FALSE
    )
  }
  if(!is.character(allocation$scale) || length(allocation$scale) != 1L ||
     is.na(allocation$scale) ||
     !allocation$scale %in% c("total_variance", "mean_variance")){
    stop(
      "Random-effect SD binding metadata are missing canonical 'allocation$scale'.",
      call. = FALSE
    )
  }
  .bt_check_random_sd_binding_source(allocation$source)

  if(identical(allocation$target, "block")){
    .bt_check_random_variance_allocation_factor_chain(
      allocation$factors,
      label = "allocation$factors"
    )
    if(!is.numeric(allocation$index) || length(allocation$index) != 1L ||
       is.na(allocation$index) || allocation$index != as.integer(allocation$index) ||
       allocation$index < 1L){
      stop(
        "Random-effect SD binding metadata are missing canonical 'allocation$index'.",
        call. = FALSE
      )
    }
    if(!is.numeric(allocation$n_targets) || length(allocation$n_targets) != 1L ||
       is.na(allocation$n_targets) ||
       allocation$n_targets != as.integer(allocation$n_targets) ||
       allocation$n_targets < 2L){
      stop(
        "Random-effect SD binding metadata are missing canonical 'allocation$n_targets'.",
        call. = FALSE
      )
    }
    if(allocation$index > allocation$n_targets){
      stop(
        "Random-effect SD binding metadata reference an allocation coordinate outside 'allocation$n_targets'.",
        call. = FALSE
      )
    }
    if(!is.null(binding) && !identical(binding$factors, allocation$factors)){
      stop(
        "Block allocation SD bindings require 'binding$factors' to match 'allocation$factors'.",
        call. = FALSE
      )
    }
  }else{
    .bt_check_random_variance_allocation_factor_chain(
      allocation$parent_factors,
      label = "allocation$parent_factors"
    )
    if(!is.null(binding) && !identical(binding$factors, allocation$parent_factors)){
      stop(
        "SD-component allocation SD bindings require 'binding$factors' to match 'allocation$parent_factors'.",
        call. = FALSE
      )
    }
  }

  invisible(TRUE)
}

.bt_check_random_sd_component_allocation <- function(allocation,
                                                     n_columns = NULL,
                                                     context = "Random-effect allocation metadata"){

  if(!is.list(allocation) || !identical(allocation$target, "sd_component")){
    stop(
      context,
      " are missing canonical 'allocation$target'.",
      call. = FALSE
    )
  }
  if(!is.character(allocation$weight_name) ||
     length(allocation$weight_name) != 1L ||
     is.na(allocation$weight_name) || !nzchar(allocation$weight_name)){
    stop(
      context,
      " are missing canonical 'allocation$weight_name'.",
      call. = FALSE
    )
  }
  if(!is.character(allocation$scale) || length(allocation$scale) != 1L ||
     is.na(allocation$scale) ||
     !allocation$scale %in% c("total_variance", "mean_variance")){
    stop(
      context,
      " are missing canonical 'allocation$scale'.",
      call. = FALSE
    )
  }
  parent_factors <- allocation$parent_factors
  if(!is.list(parent_factors)){
    stop(
      context,
      " are missing canonical 'allocation$parent_factors'.",
      call. = FALSE
    )
  }
  .bt_check_random_variance_allocation_factor_chain(
    parent_factors,
    label = "allocation$parent_factors"
  )

  K <- allocation$n_targets
  if(!is.numeric(K) || length(K) != 1L || is.na(K) ||
     K != as.integer(K) || K < 2L){
    stop(
      context,
      " are missing canonical 'allocation$n_targets'.",
      call. = FALSE
    )
  }
  K <- as.integer(K)

  leaf_index <- allocation$leaf_index_by_column
  if(!is.numeric(leaf_index) || length(leaf_index) == 0L ||
     any(is.na(leaf_index)) || any(leaf_index != as.integer(leaf_index)) ||
     any(leaf_index < 1L) || any(leaf_index > K)){
    stop(
      context,
      " are missing canonical 'allocation$leaf_index_by_column'.",
      call. = FALSE
    )
  }
  leaf_index <- as.integer(leaf_index)
  if(!setequal(unique(leaf_index), seq_len(K))){
    stop(
      context,
      " 'allocation$leaf_index_by_column' must cover every target.",
      call. = FALSE
    )
  }
  if(!is.null(n_columns)){
    check_int(n_columns, "n_columns", lower = 1, allow_NA = FALSE)
    if(length(leaf_index) != n_columns){
      stop(
        context,
        " 'allocation$leaf_index_by_column' does not match the number of random-effect columns.",
        call. = FALSE
      )
    }
  }
  leaf_names <- allocation$leaf_names
  if(!is.character(leaf_names) || length(leaf_names) != K ||
     any(is.na(leaf_names)) || any(!nzchar(leaf_names)) ||
     anyDuplicated(leaf_names)){
    stop(
      context,
      " are missing canonical 'allocation$leaf_names'.",
      call. = FALSE
    )
  }
  leaf_terms <- allocation$leaf_terms
  if(!is.character(leaf_terms) || length(leaf_terms) != K ||
     any(is.na(leaf_terms)) || any(!nzchar(leaf_terms))){
    stop(
      context,
      " are missing canonical 'allocation$leaf_terms'.",
      call. = FALSE
    )
  }

  list(
    parent_factors = parent_factors,
    leaf_index_by_column = leaf_index,
    n_targets = K
  )
}

.bt_check_random_sd_component_binding <- function(binding, n_columns = NULL,
                                                  context = "Random-effect SD binding metadata"){

  .bt_check_random_sd_binding(binding)
  if(!isTRUE(binding$true_allocation)){
    return(invisible(TRUE))
  }

  allocation <- binding$allocations[[1L]]
  if(!is.list(allocation) || !identical(allocation$target, "sd_component")){
    return(invisible(TRUE))
  }

  metadata <- .bt_check_random_sd_component_allocation(
    allocation = allocation,
    n_columns = n_columns,
    context = context
  )
  if(!identical(binding$factors, metadata$parent_factors)){
    stop(
      context,
      " require 'binding$factors' to match 'allocation$parent_factors' for SD-component allocations.",
      call. = FALSE
    )
  }

  if(!.bt_random_sd_binding_has_row_external_source(binding)){
    return(invisible(TRUE))
  }

  if(!identical(binding$application, "column")){
    stop(
      context,
      " with row-indexed SD-component allocation must use column application.",
      call. = FALSE
    )
  }
  if(length(binding$factors_by_column) == 0L){
    stop(
      context,
      " are missing canonical 'binding$factors_by_column'.",
      call. = FALSE
    )
  }
  if(length(binding$factors_by_column) != length(metadata$leaf_index_by_column)){
    stop(
      context,
      " 'binding$factors_by_column' does not match 'allocation$leaf_index_by_column'.",
      call. = FALSE
    )
  }

  for(column in seq_along(binding$factors_by_column)){
    chain <- binding$factors_by_column[[column]]
    .bt_check_random_variance_allocation_factor_chain(
      chain,
      label = "factors_by_column"
    )
    if(length(chain) != length(metadata$parent_factors) + 1L){
      stop(
        context,
        " column factor chains must contain parent factors followed by one SD-component leaf factor.",
        call. = FALSE
      )
    }
    parent_part <- chain[seq_along(metadata$parent_factors)]
    if(!identical(parent_part, metadata$parent_factors)){
      stop(
        context,
        " column factor chains must start with 'allocation$parent_factors'.",
        call. = FALSE
      )
    }
    leaf_factor <- chain[[length(chain)]]
    expected_index <- metadata$leaf_index_by_column[[column]]
    if(!identical(leaf_factor$weight_name, allocation$weight_name) ||
       leaf_factor$index != expected_index ||
       !identical(leaf_factor$scale, allocation$scale) ||
       leaf_factor$n_targets != metadata$n_targets){
      stop(
        context,
        " column factor chains do not match the resolved SD-component allocation.",
        call. = FALSE
      )
    }
  }

  invisible(TRUE)
}

.bt_random_sd_binding_source_name <- function(source){
  .bt_check_random_sd_binding_source(source)
  source$name
}

.bt_random_sd_binding_source_shape <- function(source){
  .bt_check_random_sd_binding_source(source)
  source$shape
}

.bt_random_sd_binding_source_kind <- function(source){
  .bt_check_random_sd_binding_source(source)
  source$kind
}

.bt_random_sd_binding_source_owned <- function(source){
  .bt_check_random_sd_binding_source(source)
  source$owned
}

.bt_random_sd_binding_source_jags_expression <- function(source,
                                                         row_index = NULL){

  if(inherits(source, "random_sd_source")){
    return(.bt_random_sd_source_expression(source, row_index = row_index))
  }
  .bt_check_random_sd_binding_source(source)
  if(identical(source$shape, "row")){
    if(is.null(row_index)){
      stop(
        "Row-shaped SD source '", source$name,
        "' requires a row index in generated JAGS syntax.",
        call. = FALSE
      )
    }
    return(paste0(source$name, "[", row_index, "]"))
  }

  source$name
}

.bt_random_sd_binding_source_label <- function(source){

  if(inherits(source, "random_sd_source")){
    return(.bt_random_sd_source_label(source))
  }
  .bt_check_random_sd_binding_source(source)
  if(identical(source$shape, "row")){
    return(paste0(source$name, "[row]"))
  }

  source$name
}

.bt_random_sd_binding_source_is_external_row <- function(source){

  !is.null(source) &&
    identical(.bt_random_sd_binding_source_kind(source), "external") &&
    identical(.bt_random_sd_binding_source_shape(source), "row")
}

.bt_random_sd_binding_source_is_external <- function(source){

  !is.null(source) &&
    identical(.bt_random_sd_binding_source_kind(source), "external")
}

.bt_random_sd_binding_has_row_external_source <- function(binding){

  if(is.null(binding)){
    return(FALSE)
  }
  .bt_check_random_sd_binding(binding)
  if(.bt_random_sd_binding_source_is_external_row(binding$source)){
    return(TRUE)
  }
  any(vapply(
    binding$sources_by_column,
    .bt_random_sd_binding_source_is_external_row,
    logical(1)
  ))
}

.bt_random_sd_binding_has_external_source <- function(binding){

  if(is.null(binding)){
    return(FALSE)
  }
  .bt_check_random_sd_binding(binding)
  if(.bt_random_sd_binding_source_is_external(binding$source)){
    return(TRUE)
  }
  any(vapply(
    binding$sources_by_column,
    .bt_random_sd_binding_source_is_external,
    logical(1)
  ))
}

.bt_random_sd_binding_factors_expression <- function(factors){

  .bt_random_variance_allocation_factors_expression(factors)
}

.bt_random_sd_binding_shared_source_expression <- function(binding,
                                                           row_index = NULL){

  .bt_check_random_sd_binding(binding)
  .bt_random_sd_binding_source_jags_expression(binding$source, row_index = row_index)
}

.bt_random_sd_binding_external_source_label <- function(binding){

  if(is.null(binding)){
    return("<unknown>")
  }
  .bt_check_random_sd_binding(binding)
  sources <- c(list(binding$source), binding$sources_by_column)
  for(source in sources){
    if(.bt_random_sd_binding_source_is_external(source)){
      return(.bt_random_sd_binding_source_label(source))
    }
  }

  "<unknown>"
}

.bt_random_sd_binding_context <- function(random_effects, prior_random,
                                          parameter){

  empty_context <- list(
    prior_list = list(),
    syntax = character(),
    by_block = list(),
    allocations = list()
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
    target <- .bt_random_variance_allocation_target(allocation)
    scale <- .bt_random_variance_allocation_scale(allocation)

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
        prior_list[[allocation_names$total_suffix]] <- source$prior
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
        row_indexed_source <- .bt_random_variance_allocation_source_is_row(source)
        if(component_key %in% consumed_components && !row_indexed_source){
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
            terms = terms,
            index = term_i,
            target = "block",
            scale = scale,
            parent = allocation$parent,
            source_node = source_name,
            source = source,
            factors = c(source_factors, list(factor)),
            n_targets = length(terms),
            total_name = source$total_name,
            weight_name = allocation_names$weight_name,
            total_suffix = source$total_suffix,
            weight_suffix = allocation_names$weight_suffix
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
        terms = terms,
        component_labels = component_labels,
        components = component_meta,
        target = "block",
        scale = scale,
        parent = allocation$parent,
        source_node = source_name,
        source = source,
        total_name = source$total_name,
        weight_name = allocation_names$weight_name,
        total_suffix = source$total_suffix,
        weight_suffix = allocation_names$weight_suffix
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
        weight_suffix = allocation_names$weight_suffix
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
        total_name = source$total_name,
        weight_name = allocation_names$weight_name,
        total_suffix = source$total_suffix,
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

.bt_random_variance_allocation_factor_expression <- function(factor){

  .bt_random_variance_allocation_multiplier_expression(
    weight_name = factor$weight_name,
    index = factor$index,
    scale = factor$scale,
    n_targets = factor$n_targets
  )
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

.bt_random_effect_has_row_indexed_external_sd <- function(random_term){

  .bt_random_sd_binding_has_row_external_source(random_term$sd_binding)
}

.bt_random_effect_external_sd_source_label <- function(random_term){

  .bt_random_sd_binding_external_source_label(random_term$sd_binding)
}
