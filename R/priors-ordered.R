.prior_ordered_contrast_values <- c("cumulative", "cumulative_levels")

.prior_ordered_contrast_name <- function(contrast){
  switch(
    contrast,
    cumulative        = "contr.ordered_cumulative",
    cumulative_levels = "contr.ordered_cumulative_levels"
  )
}

.prior_ordered_is_contrast_name <- function(contrast){
  contrast %in% c("contr.ordered_cumulative", "contr.ordered_cumulative_levels")
}

.prior_ordered_check_total <- function(total){
  if(!is.prior(total)){
    stop("'total' must be a prior object.", call. = FALSE)
  }

  total_is_scalar <- is.prior.simple(total) &&
    !is.prior.vector(total) &&
    !is.prior.factor(total) &&
    !is.prior.simplex(total) &&
    !is.prior.weightfunction(total) &&
    !is.prior.PET(total) &&
    !is.prior.PEESE(total) &&
    !is_prior_phacking(total) &&
    !is_prior_bias(total)

  total_is_mixture <- is.prior.mixture(total) || is.prior.spike_and_slab(total)
  if(total_is_mixture){
    components <- vapply(total, is.prior, logical(1))
    total_is_scalar <- all(vapply(total[components], function(component){
      is.prior.simple(component) &&
        !is.prior.vector(component) &&
        !is.prior.factor(component) &&
        !is.prior.simplex(component) &&
        !is.prior.weightfunction(component) &&
        !is.prior.PET(component) &&
        !is.prior.PEESE(component) &&
        !is_prior_phacking(component) &&
        !is_prior_bias(component)
    }, logical(1)))
  }

  if(!total_is_scalar){
    stop("'total' must be a scalar prior, scalar prior mixture, or scalar spike-and-slab prior.", call. = FALSE)
  }

  invisible(TRUE)
}

.prior_ordered_allocation_spec <- function(allocation, name = "allocation"){

  if(is.null(allocation)){
    return(list(type = "default_dirichlet"))
  }

  if(is.prior(allocation)){
    if(!is.prior.simplex(allocation) || !identical(allocation[["distribution"]], "dirichlet")){
      stop(paste0("The '", name, "' prior must be a Dirichlet prior."), call. = FALSE)
    }
    return(list(
      type  = "dirichlet",
      alpha = allocation$parameters[["alpha"]],
      prior = allocation
    ))
  }

  if(is.list(allocation)){
    if(is.null(names(allocation)) || any(!nzchar(names(allocation)))){
      stop(paste0("The '", name, "' allocation list must be named."), call. = FALSE)
    }
    allocations <- lapply(seq_along(allocation), function(i){
      .prior_ordered_allocation_spec(allocation[[i]], paste0(name, "$", names(allocation)[i]))
    })
    names(allocations) <- names(allocation)
    return(list(type = "by_factor", allocations = allocations))
  }

  check_real(allocation, name, check_length = 0, allow_NA = FALSE)
  if(any(!is.finite(allocation))){
    stop(paste0("The '", name, "' fixed allocation must be finite."), call. = FALSE)
  }
  if(any(allocation < 0)){
    stop(paste0("The '", name, "' fixed allocation must be non-negative."), call. = FALSE)
  }
  if(!isTRUE(all.equal(sum(allocation), 1, tolerance = sqrt(.Machine$double.eps)))){
    stop(paste0("The '", name, "' fixed allocation must sum to one."), call. = FALSE)
  }

  list(type = "fixed", weights = as.numeric(allocation))
}

.prior_ordered_allocation_for_factor <- function(allocation, factor_term, n_ordered_factors){

  if(identical(allocation$type, "by_factor")){
    if(!factor_term %in% names(allocation$allocations)){
      stop(
        "The ordered allocation list is missing an entry for factor '",
        factor_term, "'.",
        call. = FALSE
      )
    }
    return(allocation$allocations[[factor_term]])
  }

  if(n_ordered_factors > 1L && !identical(allocation$type, "default_dirichlet")){
    stop(
      "Terms with multiple ordered factors require a named 'allocation' list.",
      call. = FALSE
    )
  }

  allocation
}

.prior_ordered_bind_allocation <- function(allocation, D, factor_term){

  if(D < 1L){
    stop("Ordered factors must have at least two ordered levels.", call. = FALSE)
  }

  if(identical(allocation$type, "default_dirichlet")){
    if(D == 1L){
      return(list(type = "fixed", weights = 1))
    }
    return(list(type = "dirichlet", alpha = rep(1, D), prior = NULL))
  }

  if(identical(allocation$type, "fixed")){
    if(length(allocation$weights) != D){
      stop(
        "The fixed allocation for ordered factor '", factor_term,
        "' has length ", length(allocation$weights), ", but ", D,
        " value(s) are required.",
        call. = FALSE
      )
    }
    return(allocation)
  }

  if(identical(allocation$type, "dirichlet")){
    if(length(allocation$alpha) != D){
      stop(
        "The Dirichlet allocation for ordered factor '", factor_term,
        "' has length ", length(allocation$alpha), ", but ", D,
        " value(s) are required.",
        call. = FALSE
      )
    }
    if(any(!is.finite(allocation$alpha)) || any(allocation$alpha <= 0)){
      stop("Dirichlet allocation concentrations must be finite and positive.", call. = FALSE)
    }
    return(allocation)
  }

  stop("Unsupported ordered allocation specification.", call. = FALSE)
}

.prior_ordered_id_for_factor <- function(id, factor_term, ordered_terms){

  if(is.null(id)){
    return(NULL)
  }

  if(length(ordered_terms) == 1L && is.null(names(id))){
    return(id[[1]])
  }

  if(is.null(names(id)) || any(!nzchar(names(id)))){
    stop("For terms with multiple ordered factors, 'id' must be named by factor.", call. = FALSE)
  }
  if(!factor_term %in% names(id)){
    return(NULL)
  }

  id[[factor_term]]
}

.prior_ordered_jags_name <- function(x){
  x <- gsub("[^A-Za-z0-9_]", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  if(!nzchar(x)){
    x <- "x"
  }
  if(grepl("^[0-9]", x)){
    x <- paste0("x_", x)
  }
  x
}

.prior_ordered_allocation_key <- function(parameter_name, factor_term, slice, id){

  if(!is.null(id)){
    return(paste0("id:", id, "|factor:", factor_term))
  }

  paste0("term:", parameter_name, "|factor:", factor_term, "|slice:", slice)
}

.prior_ordered_allocation_node <- function(parameter_name, factor_term, slice, id){

  if(!is.null(id)){
    return(paste0(
      "ordered_alloc_",
      .prior_ordered_jags_name(id),
      "_",
      .prior_ordered_jags_name(factor_term)
    ))
  }

  paste0(
    parameter_name,
    "_ordered_alloc_",
    .prior_ordered_jags_name(factor_term),
    "_",
    slice
  )
}

.prior_ordered_allocation_signature <- function(record){

  spec <- record$spec
  values <- switch(
    spec$type,
    fixed     = spec$weights,
    dirichlet = spec$alpha,
    stop("Unsupported ordered allocation specification.", call. = FALSE)
  )

  paste(
    c(
      "contrast", record$contrast,
      "dim", record$dim,
      "type", spec$type,
      format(values, digits = 17, scientific = TRUE, trim = TRUE)
    ),
    collapse = "|"
  )
}

.bt_validate_ordered_shared_allocations <- function(prior_list){

  registry <- list()
  for(i in seq_along(prior_list)){
    prior <- prior_list[[i]]
    if(!is.prior.ordered(prior)){
      next
    }
    metadata <- .prior_ordered_metadata(prior)
    for(record in metadata$allocations){
      signature <- .prior_ordered_allocation_signature(record)
      if(!is.null(registry[[record$key]])){
        if(!identical(registry[[record$key]]$signature, signature)){
          allocation_label <- if(is.null(record$id)){
            paste0("allocation key '", record$key, "'")
          }else{
            paste0("allocation id '", record$id, "' for ordered factor '", record$factor, "'")
          }
          stop(
            "Shared ordered ", allocation_label,
            " is used with incompatible allocation specifications.",
            call. = FALSE
          )
        }
        next
      }
      registry[[record$key]] <- list(
        signature = signature,
        parameter = names(prior_list)[i]
      )
    }
  }

  invisible(TRUE)
}

.prior_ordered_total_name <- function(parameter_name){
  paste0(parameter_name, "_ordered_total")
}

.prior_ordered_format_number <- function(x){
  if(length(x) != 1L || is.na(x) || !is.finite(x)){
    stop("Fixed ordered allocation weights must be finite.", call. = FALSE)
  }
  format(x, digits = 16, scientific = FALSE, trim = TRUE)
}

.prior_ordered_component_contrast_dims <- function(prior){

  level_names <- .factor_level_list(prior)
  factor_terms <- attr(prior, "factor_terms", exact = TRUE)
  factor_contrasts <- attr(prior, "factor_contrasts", exact = TRUE)

  contrast_matrices <- lapply(factor_terms, function(factor_term){
    .factor_contrast_matrix(level_names[[factor_term]], factor_contrasts[[factor_term]])
  })
  names(contrast_matrices) <- factor_terms

  contrast_matrices
}

.bt_bind_ordered_prior_metadata <- function(prior, parameter_name){

  if(!is.prior.ordered(prior)){
    return(prior)
  }

  level_names <- .factor_level_list(prior)
  factor_terms <- attr(prior, "factor_terms", exact = TRUE)
  factor_contrasts <- attr(prior, "factor_contrasts", exact = TRUE)
  factor_design <- attr(prior, "factor_design", exact = TRUE)

  if(is.null(level_names) || is.null(factor_terms) ||
     is.null(factor_contrasts) || is.null(factor_design)){
    stop("Ordered prior metadata are incomplete.", call. = FALSE)
  }

  ordered_terms <- factor_terms[
    vapply(factor_terms, function(factor_term){
      .prior_ordered_is_contrast_name(factor_contrasts[[factor_term]])
    }, logical(1))
  ]
  if(length(ordered_terms) == 0L){
    stop("A prior_ordered() term must contain at least one ordered factor.", call. = FALSE)
  }

  contrast_matrices <- .prior_ordered_component_contrast_dims(prior)
  coef_grid <- expand.grid(
    lapply(contrast_matrices, function(contrast_matrix) seq_len(ncol(contrast_matrix))),
    KEEP.OUT.ATTRS = FALSE
  )
  names(coef_grid) <- factor_terms

  ordinary_terms <- setdiff(factor_terms, ordered_terms)
  if(length(ordinary_terms) == 0L){
    slice_index <- rep(1L, nrow(coef_grid))
  }else{
    slice_index <- as.integer(interaction(coef_grid[ordinary_terms], drop = TRUE, lex.order = FALSE))
  }
  theta_dim <- max(slice_index)

  if(ncol(factor_design) != nrow(coef_grid)){
    stop(
      "The ordered factor design for '", parameter_name, "' has ", ncol(factor_design),
      " coefficient columns, but the ordered expansion implies ",
      nrow(coef_grid), ". This usually means that the formula expanded a ",
      "non-hierarchical factor interaction into full level indicators. Include ",
      "the lower-order ordinary factor terms so the contrast basis is preserved.",
      call. = FALSE
    )
  }

  allocation_records <- list()
  for(factor_term in ordered_terms){
    D <- ncol(contrast_matrices[[factor_term]])
    allocation <- .prior_ordered_allocation_for_factor(
      allocation = prior$allocation,
      factor_term = factor_term,
      n_ordered_factors = length(ordered_terms)
    )
    allocation <- .prior_ordered_bind_allocation(allocation, D, factor_term)
    id <- .prior_ordered_id_for_factor(prior$id, factor_term, ordered_terms)
    slices <- if(is.null(id)) seq_len(theta_dim) else 1L

    for(slice in slices){
      key <- .prior_ordered_allocation_key(parameter_name, factor_term, slice, id)
      allocation_records[[key]] <- list(
        key    = key,
        factor = factor_term,
        slice  = if(is.null(id)) slice else NA_integer_,
        id     = id,
        node   = .prior_ordered_allocation_node(parameter_name, factor_term, slice, id),
        contrast = factor_contrasts[[factor_term]],
        dim    = D,
        spec   = allocation
      )
    }
  }

  metadata <- list(
    parameter_name = parameter_name,
    factor_terms = factor_terms,
    ordered_terms = ordered_terms,
    ordinary_terms = ordinary_terms,
    factor_contrasts = factor_contrasts,
    coefficient_grid = coef_grid,
    slice_index = slice_index,
    theta_dim = theta_dim,
    coefficient_dim = nrow(coef_grid),
    allocations = allocation_records
  )

  attr(prior, "ordered_metadata") <- metadata
  attr(prior, "coefficient_dim") <- metadata$coefficient_dim

  prior
}

.prior_ordered_metadata <- function(prior){

  metadata <- attr(prior, "ordered_metadata", exact = TRUE)
  if(is.null(metadata)){
    stop("prior_ordered() must be bound to formula metadata before this operation.", call. = FALSE)
  }

  metadata
}

.prior_ordered_allocation_for_coefficient <- function(metadata, factor_term, slice){

  for(record in metadata$allocations){
    if(!identical(record$factor, factor_term)){
      next
    }
    if(is.na(record$slice) || identical(record$slice, slice)){
      return(record)
    }
  }

  stop("Internal ordered allocation lookup failed.", call. = FALSE)
}

.prior_ordered_coefficient_expression <- function(prior, parameter_name, coefficient_i){

  metadata <- .prior_ordered_metadata(prior)
  total_name <- .prior_ordered_total_name(parameter_name)
  slice <- metadata$slice_index[[coefficient_i]]

  expression <- if(metadata$theta_dim == 1L){
    total_name
  }else{
    paste0(total_name, "[", slice, "]")
  }

  for(factor_term in metadata$ordered_terms){
    allocation <- .prior_ordered_allocation_for_coefficient(metadata, factor_term, slice)
    increment <- metadata$coefficient_grid[[factor_term]][[coefficient_i]]
    term_expression <- if(identical(allocation$spec$type, "fixed")){
      .prior_ordered_format_number(allocation$spec$weights[[increment]])
    }else{
      paste0(allocation$node, "[", increment, "]")
    }
    expression <- paste0(expression, " * ", term_expression)
  }

  expression
}

.prior_ordered_dirichlet_records <- function(prior){
  metadata <- .prior_ordered_metadata(prior)
  records <- metadata$allocations[
    vapply(metadata$allocations, function(record){
      identical(record$spec$type, "dirichlet")
    }, logical(1))
  ]
  records
}

.prior_ordered_bridge_check <- function(prior){

  if(.is_prior_expression(prior$total)){
    stop(
      "Bridge sampling for prior_ordered() does not support parameter expressions in 'total'.",
      call. = FALSE
    )
  }
  if(is.prior.mixture(prior$total) || is.prior.spike_and_slab(prior$total)){
    stop(
      "Bridge sampling for prior_ordered() is only available when 'total' is a simple scalar prior.",
      call. = FALSE
    )
  }
  if(!is.prior.simple(prior$total) || is.prior.vector(prior$total) ||
     is.prior.factor(prior$total) || is.prior.simplex(prior$total)){
    stop(
      "Bridge sampling for prior_ordered() requires a simple scalar 'total' prior.",
      call. = FALSE
    )
  }
  if(is.prior.discrete(prior$total)){
    stop(
      "Bridge sampling for prior_ordered() requires a continuous or point-valued 'total' prior.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

.prior_ordered_total_monitor_names <- function(prior, parameter_name){

  metadata <- .prior_ordered_metadata(prior)
  total_name <- .prior_ordered_total_name(parameter_name)

  if(metadata$theta_dim == 1L){
    return(total_name)
  }

  paste0(total_name, "[", seq_len(metadata$theta_dim), "]")
}

.prior_ordered_allocation_rng <- function(spec, n){

  if(identical(spec$type, "fixed")){
    return(matrix(rep(spec$weights, each = n), nrow = n))
  }

  if(identical(spec$type, "dirichlet")){
    prior <- prior("dirichlet", list(alpha = spec$alpha))
    return(rng(prior, n))
  }

  stop("Unsupported ordered allocation specification.", call. = FALSE)
}

.prior_ordered_default_bound <- function(prior, parameter_name = ".ordered"){

  if(!is.null(attr(prior, "ordered_metadata", exact = TRUE))){
    return(prior)
  }

  n_levels <- attr(prior, "levels", exact = TRUE)
  if(is.null(n_levels) || is.na(n_levels)){
    n_levels <- if(identical(prior$contrast, "cumulative")) 2L else 1L
    warning(
      "Number of ordered factor levels was not specified; assuming ",
      n_levels, " level(s).",
      call. = FALSE
    )
  }

  level_names <- attr(prior, "level_names", exact = TRUE)
  if(is.null(level_names)){
    level_names <- seq_len(n_levels)
  }
  attr(prior, "level_names") <- level_names
  attr(prior, "levels") <- length(level_names)
  attr(prior, "factor_terms") <- ".factor"
  attr(prior, "factor_contrasts") <- stats::setNames(
    .prior_ordered_contrast_name(prior$contrast),
    ".factor"
  )
  design_info <- .factor_term_design_from_metadata(prior)
  attr(prior, "factor_design") <- design_info$design
  attr(prior, "factor_cell_names") <- design_info$cell_names

  .bt_bind_ordered_prior_metadata(prior, parameter_name)
}

.prior_ordered_rng <- function(prior, n, transform_factor_samples = TRUE, quantity = "level"){

  check_char(quantity, "quantity", allow_values = c("level", "coefficient", "increment", "total", "allocation"))
  prior <- .prior_ordered_default_bound(prior)
  metadata <- .prior_ordered_metadata(prior)

  theta <- if(metadata$theta_dim == 1L){
    matrix(rng(prior$total, n), nrow = n, ncol = 1L)
  }else{
    do.call(cbind, replicate(metadata$theta_dim, rng(prior$total, n), simplify = FALSE))
  }

  if(quantity == "total"){
    colnames(theta) <- if(metadata$theta_dim == 1L) "total" else paste0("total[", seq_len(metadata$theta_dim), "]")
    return(theta)
  }

  allocation_samples <- list()
  for(record in metadata$allocations){
    if(record$key %in% names(allocation_samples)){
      next
    }
    allocation_samples[[record$key]] <- .prior_ordered_allocation_rng(record$spec, n)
  }

  if(quantity == "allocation"){
    out <- do.call(cbind, allocation_samples)
    colnames(out) <- unlist(lapply(names(allocation_samples), function(key){
      paste0(key, "[", seq_len(ncol(allocation_samples[[key]])), "]")
    }), use.names = FALSE)
    return(out)
  }

  coefficients <- matrix(NA_real_, nrow = n, ncol = metadata$coefficient_dim)
  for(coefficient_i in seq_len(metadata$coefficient_dim)){
    slice <- metadata$slice_index[[coefficient_i]]
    value <- theta[, slice]
    for(factor_term in metadata$ordered_terms){
      record <- .prior_ordered_allocation_for_coefficient(metadata, factor_term, slice)
      increment <- metadata$coefficient_grid[[factor_term]][[coefficient_i]]
      value <- value * allocation_samples[[record$key]][, increment]
    }
    coefficients[, coefficient_i] <- value
  }
  colnames(coefficients) <- .JAGS_prior_factor_names(metadata$parameter_name, prior)

  if(quantity %in% c("coefficient", "increment") || !transform_factor_samples){
    return(coefficients)
  }

  .transform_factor_contrast_samples(
    coefficient_samples = coefficients,
    metadata            = prior,
    parameter           = metadata$parameter_name,
    transformed_class   = "mixed_posteriors.ordered_transformed"
  )
}
