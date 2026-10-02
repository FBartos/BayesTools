# Ordered primitives and deterministic recipes --------------------------------

.bt_ordered_stop <- function(message, class = "BayesTools_ordered_coordinates_unavailable"){

  stop(structure(list(message = message, call = NULL),
    class = c(class, "BayesTools_ordered_unavailable", "error", "condition")))
}

.bt_ordered_spec <- function(parameter, prior){

  metadata <- .prior_ordered_metadata(prior)
  total_names <- .prior_ordered_total_monitor_names(prior, parameter)
  allocations <- lapply(metadata$allocations, function(record){
    record$coordinates <- paste0(record$node, "[", seq_len(record$dim), "]")
    record$gamma_coordinates <- if(identical(record$spec$type, "dirichlet")){
      paste0(.JAGS_prior_dirichlet_eta_name(record$node), "[", seq_len(record$dim), "]")
    }else character()
    record
  })
  list(parameter = parameter, prior = prior, metadata = metadata,
    coefficient_names = .JAGS_prior_factor_names(parameter, prior),
    coefficient_grid = metadata$coefficient_grid, slice_index = metadata$slice_index,
    total_names = total_names, total_prior = prior$total,
    total_node = .bt_dnode_ordered_total(parameter, prior),
    allocations = allocations,
    label_parts = .bt_label_parts_term(parameter, prior))
}

#' @title Fitted ordered parameter specifications
#' @description Reads the persisted total, allocation and coefficient recipes
#' of ordered terms. These backend recipes do not add public catalog quantities.
#' @param fit a model fitted with [JAGS_fit()].
#' @param parameters optional names of ordered terms; code{NULL} selects all.
#' @param weights optional named numeric fitted-coordinate projection weights.
#' @param draws optional numeric source draws for the projection.
#' @return With no code{weights}, a named list of term specifications containing
#' the bound code{prior}, coefficient names and grid, slice indices, total names
#' and prior, deterministic total recipe, allocation records (key, factor, slice,
#' id, normalized and gamma coordinates, alpha or fixed weights), and label parts.
#' An empty named list is returned when no ordered terms are present. With
#' code{weights}, an ordered projection specification, including resolved values
#' and declared point states when code{draws} are supplied, is returned.
#' @details Missing fitted numeric provenance raises
#' code{BayesTools_ordered_metadata_unavailable} and
#' code{BayesTools_refit_required}: refit with this version. Missing source
#' draws raise code{BayesTools_ordered_coordinates_unavailable}; invalid component
#' states raise code{BayesTools_ordered_invalid_state}. Unclassified expression
#' structure raises code{BayesTools_ordered_expression_unavailable}. These
#' coordinate/expression conditions inherit code{BayesTools_ordered_unavailable}.
#' @export
JAGS_ordered_parameter_spec <- function(fit, parameters = NULL, weights = NULL, draws = NULL){

  .bt_require_fit_contract(fit, "fit")
  check_char(parameters, "parameters", check_length = 0, allow_NULL = TRUE, allow_NA = FALSE)
  priors <- attr(fit, "prior_list", exact = TRUE)
  ordered <- priors[vapply(priors, is.prior.ordered, logical(1))]
  if(!is.null(parameters)){
    if(any(!parameters %in% names(ordered))){
      stop("'parameters' must name ordered terms of 'fit'.", call. = FALSE)
    }
    ordered <- ordered[unique(parameters)]
  }
  for(prior in ordered){
    metadata <- attr(prior, "ordered_metadata", exact = TRUE)
    if(is.null(metadata) || is.null(metadata$numeric_literals)){
      .bt_stop_refit_required(
        "Fitted ordered numeric provenance is unavailable. Refit the model with this version of BayesTools.",
        class = "BayesTools_ordered_metadata_unavailable")
    }
  }
  specs <- Map(.bt_ordered_spec, names(ordered), ordered)
  if(length(specs) == 0L) names(specs) <- character()
  if(is.null(weights)) return(specs)
  .bt_ordered_projection(specs, weights, draws, prior_list = priors)
}

.bt_dnode_ordered_total <- function(parameter, prior){

  total <- prior$total
  if(!is.prior.mixture(total)) return(NULL)
  metadata <- .prior_ordered_metadata(prior)
  name <- .prior_ordered_total_name(parameter)
  if(metadata$theta_dim == 1L) return(.bt_dnode_prior_mixture(name, total))
  if(!is.prior.spike_and_slab(total)){
    .bt_ordered_stop("Multi-slice ordered interactions require a simple scalar 'total' prior.")
  }
  variable <- .get_spike_and_slab_variable(total)
  coordinates <- paste0(name, "_variable[", seq_len(metadata$theta_dim), "]")
  location <- if(is.prior.point(variable) && !.is_prior_expression(variable)) variable$parameters$location
  .bt_deterministic_node("prior_mixture", name,
    .prior_ordered_total_monitor_names(prior, parameter),
    dependencies = c(if(is.null(location)) coordinates, paste0(name, "_indicator")),
    parameter = parameter,
    spec = list(kind = "ordered_spike_and_slab", indicator = paste0(name, "_indicator"),
      components = list(list(name = paste0(name, "_variable"), coordinates = coordinates,
        location = location, dependencies = if(is.null(location)) coordinates else character()))))
}

.bt_dnode_ordered_allocation <- function(parameter, record){

  .bt_deterministic_node("ordered_allocation", record$node, record$coordinates,
    dependencies = record$gamma_coordinates, parameter = parameter, spec = record)
}

.bt_dnode_ordered_coefficients <- function(spec){

  .bt_deterministic_node("ordered_coefficient", spec$parameter, spec$coefficient_names,
    dependencies = c(spec$total_names,
      unlist(lapply(spec$allocations, `[[`, "coordinates"), use.names = FALSE)),
    parameter = spec$parameter, spec = spec)
}

.bt_ordered_allocation_values <- function(record, draws){

  if(identical(record$spec$type, "fixed")){
    return(matrix(rep(record$spec$weights, each = nrow(draws)), nrow = nrow(draws)))
  }
  if(all(record$gamma_coordinates %in% colnames(draws))){
    eta <- draws[, record$gamma_coordinates, drop = FALSE]
    if(any(!is.finite(eta)) || any(eta < 0) || any(rowSums(eta) <= 0)){
      .bt_ordered_stop("Ordered gamma coordinates must be finite and nonnegative with positive row sums.",
        "BayesTools_ordered_invalid_state")
    }
    return(eta / rowSums(eta))
  }
  if(all(record$coordinates %in% colnames(draws))){
    values <- draws[, record$coordinates, drop = FALSE]
    if(any(!is.finite(values)) || any(values < 0) ||
       any(abs(rowSums(values) - 1) > .Machine$double.eps * max(8, record$dim))){
      .bt_ordered_stop("Ordered normalized allocation coordinates must be finite nonnegative simplex weights.",
        "BayesTools_ordered_invalid_state")
    }
    return(values)
  }
  NULL
}

.bt_ordered_total_values <- function(spec, lookup){

  if(!is.null(spec$total_node)){
    values <- .bt_deterministic_node_evaluate(spec$total_node, lookup)
    if(!is.null(values)) return(values)
  }
  total <- spec$total_prior
  if(is.prior.point(total) && !.is_prior_expression(total)){
    return(matrix(total$parameters$location, lookup$n, length(spec$total_names)))
  }
  .bt_deterministic_lookup_values(lookup, spec$total_names)
}

.bt_dnode_ordered_allocation_evaluate <- function(node, lookup){
  .bt_ordered_allocation_values(node$spec, lookup$draws)
}

.bt_dnode_ordered_coefficient_evaluate <- function(node, lookup){

  spec <- node$spec
  totals <- .bt_ordered_total_values(spec, lookup)
  allocations <- lapply(spec$allocations, .bt_ordered_allocation_values, draws = lookup$draws)
  if(is.null(totals) || any(vapply(allocations, is.null, logical(1)))) return(NULL)
  values <- matrix(NA_real_, lookup$n, length(spec$coefficient_names))
  for(i in seq_along(spec$coefficient_names)){
    slice <- spec$slice_index[[i]]
    value <- totals[, slice]
    for(factor in spec$metadata$ordered_terms){
      record <- .prior_ordered_allocation_for_coefficient(spec$metadata, factor, slice)
      value <- value * allocations[[record$key]][, spec$coefficient_grid[[factor]][[i]]]
    }
    values[, i] <- value
  }
  values
}

.bt_dnode_ordered_allocation_emit <- function(node){

  record <- node$spec
  paste0(record$coordinates, " <- ", record$gamma_coordinates, " / sum(",
    .JAGS_prior_dirichlet_eta_name(record$node), "[1:", record$dim, "])")
}

.bt_dnode_ordered_coefficient_emit <- function(node){

  spec <- node$spec
  paste0(spec$coefficient_names, " <- ", vapply(seq_along(spec$coefficient_names), function(i){
    .prior_ordered_coefficient_expression(spec$prior, spec$parameter, i)
  }, character(1)))
}
