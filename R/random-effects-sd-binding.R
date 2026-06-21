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

