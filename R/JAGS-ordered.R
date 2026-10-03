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
#' @param parameters optional names of ordered terms; \code{NULL} selects all.
#' @param weights optional named numeric fitted-coordinate projection weights.
#' @param draws optional numeric source draws for sampled parameters or a projection.
#' @return With no \code{weights}, a named list of term specifications containing
#' the bound \code{prior}, coefficient names and grid, slice indices, total names
#' and prior, deterministic total recipe, allocation records (key, factor, slice,
#' id, normalized and gamma coordinates, alpha or fixed weights), and label parts.
#' \code{source_coordinates} combines total monitors, total primitive dependencies,
#' and Gamma monitors; \code{total_auxiliary_coordinates} lists the total's other
#' monitored coordinates, including indicators, inclusion probabilities and slab
#' sources. Some total-component dependencies are optional unmonitored sources;
#' the monitored total is their fitted snapshot alternative. No supplied draws
#' means no values field. With \code{weights = NULL} and supplied \code{draws},
#' each term also contains \code{sampled_parameters}: \code{values} is the fitted
#' total/share matrix, \code{label_parts} names its rendered labels, and
#' \code{allocation_names} identifies its sampled share columns. Fixed totals
#' materialize every slice; fixed allocations produce no sampled share columns.
#'
#' An empty named list is returned when no ordered terms are present. With
#' \code{weights}, the result is a \code{BayesTools_ordered_projection} containing
#' bound \code{specs}, fitted \code{weights}, original \code{requested_weights},
#' static tensor \code{contractions},
#' and \code{ordinary_weights}. With resolved draws it also contains row-aligned
#' \code{values}, \code{atom} (point location, otherwise \code{NA}), \code{state}
#' (\code{"point"}, \code{"continuous"}, or \code{"unavailable"}), and
#' \code{exact} (a full-simplex identity). \code{total_contraction_weights} holds
#' the per-term matrices of total slopes, and \code{allocation_contractions}
#' holds complete coordinate contractions per allocation key, including totals
#' and every retained other factor. \code{reason} is a typed condition when the
#' scalar measure is structurally unavailable.
#' Total-coordinate weights expand to the corresponding complete coefficient
#' slice, so a total target uses its declared full-simplex identity.
#' Available ordered semantic values in mixed posteriors, factor levels and
#' marginal views use these primitive projections, including continuous
#' cumulative effects. Declared view transformations and undefined rows are kept.
#' @details Missing fitted numeric provenance raises
#' \code{BayesTools_ordered_metadata_unavailable} and
#' \code{BayesTools_refit_required}: refit with this version. Missing source
#' draws raise \code{BayesTools_ordered_coordinates_unavailable}; invalid component
#' states raise \code{BayesTools_ordered_invalid_state}. Continuous total families
#' remain structurally continuous with expression parameters. Bernoulli totals
#' retain their declared support even with an expression probability. Expression
#' point locations can supply fitted snapshot values but have an unavailable
#' scalar measure with a \code{BayesTools_ordered_expression_unavailable} reason.
#' Without a registered certified ancestor recipe, expression totals cannot
#' supply numerical charts or correct replay after changing an ancestor. These
#' coordinate/expression conditions inherit \code{BayesTools_ordered_unavailable}.
#' Missing or inconsistent fitted provenance reports: "Fitted ordered numeric
#' provenance is unavailable. Refit the model with this version of BayesTools."
#' Missing projection totals report: "Ordered total sources for '<term>' are
#' unavailable in 'draws'. Include its declared source coordinates." Invalid
#' total states report: "Ordered total indicator '<coordinate>' does not select
#' a declared component."
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
    if(!.bt_ordered_metadata_valid(prior)){
      .bt_stop_refit_required(
        "Fitted ordered numeric provenance is unavailable. Refit the model with this version of BayesTools.",
        class = "BayesTools_ordered_metadata_unavailable")
    }
  }
  specs <- Map(.bt_ordered_spec, names(ordered), ordered)
  if(length(specs)){
    catalog <- parameter_catalog(fit)
    coordinates <- parameter_coordinates(fit)
    specs <- lapply(specs,function(spec){
      spec$total_quantities <- catalog$quantities[match(spec$total_names,catalog$quantities$canonical_name),,drop=FALSE]
      gamma_coordinates <- unlist(lapply(spec$allocations,`[[`,"gamma_coordinates"),use.names=FALSE)
      total_monitors <- setdiff(.JAGS_monitor.ordered(spec$prior,spec$parameter),
        c(spec$parameter,.bt_parameter_coordinates_base(gamma_coordinates)))
      total_coordinates <- coordinates$coordinate_name[coordinates$monitor_name %in% total_monitors]
      spec$total_auxiliary_coordinates <- setdiff(total_coordinates,spec$total_names)
      spec$source_coordinates <- unique(c(spec$total_names,spec$total_node$dependencies,gamma_coordinates))
      if(!is.null(draws) && is.null(weights)){
        values <- .bt_deterministic_draws_matrix(draws)
        source <- .bt_ordered_source_new(spec$parameter,list(spec),
          list(.bt_ordered_source_rows(spec,values,seq_len(nrow(values)))),
          rep(1L,nrow(values)),seq_len(nrow(values)))
        display <- .bt_ordered_source_display(source)
        spec$sampled_parameters <- list(values=display$samples,label_parts=display$parts,
          allocation_names=setdiff(colnames(display$samples),spec$total_names))
      }
      spec
    })
  }
  if(length(specs) == 0L) names(specs) <- character()
  if(is.null(weights)) return(specs)
  .bt_ordered_projection(specs, weights, draws, prior_list = priors)
}

.bt_ordered_metadata_valid <- function(prior){

  metadata <- attr(prior,"ordered_metadata",exact=TRUE)
  required <- c("parameter_name","factor_terms","ordered_terms","ordinary_terms","factor_contrasts",
    "coefficient_grid","slice_index","theta_dim","coefficient_dim","allocations","numeric_literals")
  if(!is.list(metadata) || !all(required %in% names(metadata))) return(FALSE)
  count <- function(x) is.numeric(x) && length(x)==1L && is.finite(x) && x>=1 && x==as.integer(x)
  if(!count(metadata$theta_dim) || !count(metadata$coefficient_dim) ||
     !is.data.frame(metadata$coefficient_grid) || nrow(metadata$coefficient_grid)!=metadata$coefficient_dim ||
     !identical(names(metadata$coefficient_grid),metadata$factor_terms) ||
     !is.numeric(metadata$slice_index) || length(metadata$slice_index)!=metadata$coefficient_dim ||
     anyNA(metadata$slice_index) || any(!metadata$slice_index %in% seq_len(metadata$theta_dim)) ||
     !is.list(metadata$allocations) || !length(metadata$allocations) || is.null(names(metadata$allocations))) return(FALSE)
  valid <- vapply(metadata$allocations,function(record){
    if(!is.list(record) || !all(c("key","factor","dim","spec") %in% names(record)) ||
       !is.character(record$factor) || length(record$factor)!=1L || !is.list(record$spec) ||
       !is.character(record$spec$type) || length(record$spec$type)!=1L) return(FALSE)
    values <- if(identical(record$spec$type,"fixed")) record$spec$weights else record$spec$alpha
    is.character(record$key) && length(record$key)==1L && record$factor %in% metadata$ordered_terms &&
      count(record$dim) && record$spec$type %in% c("fixed","dirichlet") && is.numeric(values) &&
      length(values)==record$dim && all(is.finite(values)) &&
      if(identical(record$spec$type,"fixed")) all(values>=0) && abs(sum(values)-1)<=.Machine$double.eps*max(8,length(values)) else all(values>0)
  },logical(1))
  all(valid) && identical(metadata$numeric_literals$total,.bt_ordered_numeric_provenance(prior$total)) &&
    identical(metadata$numeric_literals$total_syntax,.JAGS_prior.ordered_total(prior$total,
      .prior_ordered_total_name(metadata$parameter_name),metadata$theta_dim,.bt_dnode_ordered_total(metadata$parameter_name,prior))) &&
    identical(metadata$numeric_literals$allocations,lapply(metadata$allocations,function(record){
      values <- if(identical(record$spec$type,"fixed")) record$spec$weights else record$spec$alpha
      vapply(values,.prior_ordered_format_number,character(1))
    }))
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
    dependencies = c(record$gamma_coordinates,record$coordinates), parameter = parameter, spec = record)
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
    eta_sum <- rowSums(eta)
    if(any(!is.finite(eta)) || any(eta < 0) || any(!is.finite(eta_sum)) || any(eta_sum <= 0)){
      .bt_ordered_stop("Ordered gamma coordinates must be finite and nonnegative with finite positive row sums.",
        "BayesTools_ordered_invalid_state")
    }
    return(eta / eta_sum)
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
