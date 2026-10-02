.bt_ordered_allocation_shares <- function(parameter, prior, model_samples){

  metadata    <- .prior_ordered_metadata(prior)
  level_names <- .factor_level_list(prior)
  term_parts  <- .bt_label_parts_term(parameter, prior)
  term_label  <- paste(term_parts$components, collapse = ":")
  several_factors <- length(metadata$ordered_terms) > 1L

  columns <- list()
  parts   <- list()
  for(record in metadata$allocations){
    if(!identical(record$spec$type, "dirichlet")){
      next
    }
    eta_columns <- paste0(
      .JAGS_prior_dirichlet_eta_name(record$node),
      "[", seq_len(record$dim), "]"
    )
    if(!all(eta_columns %in% colnames(model_samples))){
      stop(
        "The allocation of the ordered prior '", parameter,
        "' was not monitored; refit the model with this version of BayesTools.",
        call. = FALSE
      )
    }
    eta <- model_samples[, eta_columns, drop = FALSE]
    # an increment of a cumulative contrast reaches the level after it; the
    # first increment of a 'cumulative_levels' contrast reaches the first level
    factor_levels <- level_names[[record$factor]]
    tokens <- factor_levels[length(factor_levels) - record$dim + seq_len(record$dim)]
    suffix <- paste0(
      "_ordered_allocation",
      if(several_factors) paste0("_", record$factor),
      if(metadata$theta_dim > 1L && is.null(record$id)) paste0("[", record$slice, "]")
    )
    share_names <- paste0(parameter, suffix, "[", tokens, "]")
    share <- eta / rowSums(eta)
    colnames(share) <- share_names
    columns[[length(columns) + 1L]] <- share
    for(i in seq_along(share_names)){
      parts[[share_names[i]]] <- .bt_label_parts(
        components        = paste0(term_label, suffix, "[", tokens[i], "]"),
        formula_parameter = term_parts$formula_parameter,
        selector          = share_names[i]
      )
    }
  }

  list(
    samples = if(length(columns) > 0L) do.call(cbind, columns) else model_samples[, 0L, drop = FALSE],
    parts   = parts
  )
}

# Retained source rows stay on the fitted scale, independently of transformations
# of effect columns. Models include the complete declared ensemble, even when no
# row was selected from a model.
.bt_ordered_source_model <- function(prior, parameter){

  if(is.prior.ordered(prior)) return(.bt_ordered_spec(parameter, prior))
  list(parameter = parameter, prior = prior,
    parameterization = if(is.null(prior) || .posterior_atoms_is_zero_point(prior)) "absent" else "unavailable")
}

.bt_ordered_source_rows <- function(spec, model_samples, rows){

  if(!is.null(spec$parameterization)) return(matrix(numeric(), length(rows), 0L))
  draws <- model_samples[rows, , drop = FALSE]
  totals <- .bt_ordered_total_values(spec, .bt_deterministic_lookup(draws))
  if(is.null(totals)) .bt_ordered_stop(paste0("Ordered total sources for '", spec$parameter,
    "' are unavailable. Refit the model with this version of BayesTools."))
  colnames(totals) <- spec$total_names
  allocations <- lapply(spec$allocations, function(record){
    values <- .bt_ordered_allocation_values(record, draws)
    if(is.null(values)) .bt_ordered_stop(paste0("Ordered allocation sources for '", spec$parameter,
      "' are unavailable. Refit the model with this version of BayesTools."))
    colnames(values) <- record$coordinates
    values
  })
  sources <- c(list(totals), allocations)
  if(!is.null(spec$total_node)){
    indicator <- spec$total_node$spec$indicator
    if(!indicator %in% colnames(draws)) .bt_ordered_stop(paste0("Ordered total indicator '", indicator,
      "' is unavailable. Refit the model with this version of BayesTools."))
    sources <- c(sources, list(draws[, indicator, drop = FALSE]))
  }
  do.call(cbind, sources)
}

.bt_ordered_source_new <- function(parameter, models, sources, model, draw_index){

  columns <- unique(unlist(lapply(sources, colnames), use.names = FALSE))
  values <- matrix(NA_real_, length(model), length(columns), dimnames = list(NULL, columns))
  for(i in seq_along(models)){
    rows <- which(model == i)
    if(!length(rows)) next
    source <- sources[[i]]
    values[rows, colnames(source)] <- source
  }
  structure(list(version = 1L, parameter = parameter, models = models,
    model = as.integer(model), draw_index = as.integer(draw_index), primitives = values),
    class = c("BayesTools_ordered_source", "list"))
}

.bt_ordered_source_validate <- function(value){

  if(!inherits(value, "BayesTools_ordered_source") || !is.list(value) ||
     !identical(value$version, 1L) || !is.character(value$parameter) || length(value$parameter) != 1L ||
     !is.list(value$models) || !.bt_meta_is_index(value$model) || !.bt_meta_is_index(value$draw_index) ||
     any(value$model > length(value$models)) || !is.matrix(value$primitives) || !is.numeric(value$primitives) ||
     is.null(colnames(value$primitives)) || anyDuplicated(colnames(value$primitives)) ||
     length(value$model) != nrow(value$primitives) || length(value$draw_index) != nrow(value$primitives)){
    return("it must contain validated row-aligned ordered primitive sources and model provenance")
  }
  for(i in seq_along(value$models)){
    spec <- value$models[[i]]
    rows <- which(value$model == i)
    if(!is.null(spec$parameterization)){
      if(!spec$parameterization %in% c("absent", "unavailable")) return("its model parameterization is invalid")
      next
    }
    if(!is.prior.ordered(spec$prior) || !identical(spec$parameter, value$parameter)){
      return("its models must contain authoritative bound ordered specifications")
    }
    if(!length(rows)) next
    if(!all(spec$total_names %in% colnames(value$primitives))) return("its declared total coordinates are missing")
    for(record in spec$allocations){
      if(!all(record$coordinates %in% colnames(value$primitives))) return("its declared allocation coordinates are missing")
      if(length(rows)) .bt_ordered_allocation_values(record, value$primitives[rows, , drop = FALSE])
    }
  }
  NULL
}

.bt_ordered_source_subset <- function(source, rows){

  source$model <- source$model[rows]
  source$draw_index <- source$draw_index[rows]
  source$primitives <- source$primitives[rows, , drop = FALSE]
  source
}

.bt_ordered_source_display <- function(source){

  models <- source$models
  ordered <- models[vapply(models, function(spec) is.null(spec$parameterization), logical(1))]
  if(!length(ordered) || any(vapply(models, function(spec) identical(spec$parameterization, "unavailable"), logical(1)))){
    return(NULL)
  }
  template <- ordered[[1L]]
  totals <- matrix(NA_real_, length(source$model), length(template$total_names),
    dimnames = list(NULL, template$total_names))
  for(i in seq_along(models)){
    rows <- which(source$model == i)
    if(!length(rows)) next
    spec <- models[[i]]
    if(identical(spec$parameterization, "absent")) totals[rows, ] <- 0 else{
      totals[rows, ] <- source$primitives[rows, spec$total_names, drop = FALSE]
    }
  }
  parts <- lapply(template$total_names, function(name){
    suffix <- substring(name, nchar(source$parameter) + 1L)
    .bt_label_parts(components = paste0(paste(template$label_parts$components, collapse = ":"), suffix),
      formula_parameter = template$label_parts$formula_parameter, selector = name)
  })
  names(parts) <- template$total_names
  columns <- list(totals)
  undefined <- character()
  levels <- .factor_level_list(template$prior)
  for(factor in template$metadata$ordered_terms){
    records <- unlist(lapply(ordered, function(spec){
      spec$allocations[vapply(spec$allocations, function(record) identical(record$factor, factor), logical(1))]
    }), recursive = FALSE)
    if(!any(vapply(records, function(record) identical(record$spec$type, "dirichlet"), logical(1)))) next
    sliced <- template$metadata$theta_dim > 1L && any(vapply(records, function(record) !is.na(record$slice), logical(1)))
    slices <- if(sliced) seq_len(template$metadata$theta_dim) else 1L
    for(slice in slices){
      record <- records[[1L]]
      tokens <- levels[[factor]][length(levels[[factor]]) - record$dim + seq_len(record$dim)]
      suffix <- paste0("_ordered_allocation", if(length(template$metadata$ordered_terms) > 1L) paste0("_", factor),
        if(sliced) paste0("[", slice, "]"))
      names <- paste0(source$parameter, suffix, "[", tokens, "]")
      values <- matrix(NA_real_, length(source$model), record$dim, dimnames = list(NULL, names))
      for(i in seq_along(models)){
        rows <- which(source$model == i)
        spec <- models[[i]]
        if(!length(rows) || !is.null(spec$parameterization)) next
        own <- .prior_ordered_allocation_for_coefficient(spec$metadata, factor, slice)
        coordinates <- spec$allocations[[own$key]]$coordinates
        values[rows, ] <- source$primitives[rows, coordinates, drop = FALSE]
      }
      columns[[length(columns) + 1L]] <- values
      undefined <- c(undefined, stats::setNames(rep("ordered_parameterization", length(names)), names))
      for(j in seq_along(names)){
        parts[[names[[j]]]] <- .bt_label_parts(
          components = paste0(paste(template$label_parts$components, collapse=":"), suffix, "[", tokens[[j]], "]"),
          formula_parameter = template$label_parts$formula_parameter, selector = names[[j]])
      }
    }
  }
  list(samples = do.call(cbind, columns), parts = parts, undefined_draws = undefined)
}

.bt_ordered_raw_table_samples <- function(x){

  source <- .bt_meta_get(x, "ordered_source")
  if(is.null(source)){
    .bt_ordered_stop("Raw ordered totals and allocations are unavailable without retained source draws. Recreate mixed posteriors from the source fits with this version of BayesTools.",
      "BayesTools_ordered_metadata_unavailable")
  }
  display <- .bt_ordered_source_display(source)
  if(is.null(display)) return(NULL)
  samples <- display$samples
  quantities <- .bt_draws_quantity_table(colnames(samples), rep("", ncol(samples)),
    rep(list(character()), ncol(samples)), rep(list(numeric()), ncol(samples)), unname(display$parts))
  samples <- .bt_meta_set(samples, "quantities", quantities)
  samples <- .bt_meta_set(samples, "undefined_draws", display$undefined_draws)
  samples
}

