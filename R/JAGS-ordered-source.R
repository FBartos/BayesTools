# Labels of the quantities of an ordered term follow the label grammar: `term{...}`
# marks a coefficient of the term that is not a level (`f{total}`), and
# `term[tag: level]` a per-level quantity (`f[share: mid]`, as `[dif: level]`
# marks a transformed level). The label is built from the term components in
# term order. The ordered components are the ordered factors; the slice
# components are the ordinary factors, whose contrast coordinates index the
# slices of a term with 'theta_dim' > 1: a total is marked on the last
# component that is not a slice component (`f{total}`, `f:x{total}`,
# `g{1}:f{total}`), a share on the component of its ordered factor
# (`f[share: mid]:o`), and the slice index '{s}' follows the last slice
# component (`f{total}:g:h{1}`). Only rendered labels follow this grammar; the
# selectors and JAGS nodes keep the backend names.
.bt_ordered_label_text <- function(components, metadata, quantity, factor = NULL,
                                   level = NULL, slice = NULL, inclusion = FALSE){

  slice_components <- components %in% metadata$ordinary_terms
  text <- components
  if(identical(quantity, "total")){
    marked <- max(c(1L, which(!slice_components)))
    text[marked] <- paste0(text[marked], "{total}")
  }else{
    marked <- match(factor, components)
    if(is.na(marked)) marked <- length(components)
    text[marked] <- paste0(text[marked], "[share: ", level, "]")
  }
  if(!is.null(slice) && any(slice_components)){
    last <- max(which(slice_components))
    text[last] <- paste0(text[last], "{", slice, "}")
  }
  paste0(paste(text, collapse = ":"), if(inclusion) " (inclusion)")
}

# The slice index of the totals of an ordered term: NULL without slices.
.bt_ordered_total_slices <- function(metadata){

  if(metadata$theta_dim > 1L) seq_len(metadata$theta_dim)
}

# Label parts of the shares of one allocation record, named by their selectors.
# With 'legacy', the labels the shares had before their labels followed the
# label grammar (`f_ordered_allocation[mid]`), which the catalog keeps as aliases.
.bt_ordered_allocation_parts <- function(spec, record, sliced = NULL, slice = record$slice,
                                         legacy = FALSE){

  if(is.null(sliced)) sliced <- spec$metadata$theta_dim > 1L && is.null(record$id)
  levels <- .factor_level_list(spec$prior)[[record$factor]]
  tokens <- levels[length(levels)-record$dim+seq_len(record$dim)]
  suffix <- paste0("_ordered_allocation",
    if(length(spec$metadata$ordered_terms)>1L) paste0("_",record$factor),
    if(sliced) paste0("[",slice,"]"))
  names <- paste0(spec$parameter,suffix,"[",tokens,"]")
  parts <- lapply(seq_along(names),function(i){
    text <- if(legacy){
      paste0(paste(spec$label_parts$components,collapse=":"),suffix,"[",tokens[[i]],"]")
    }else{
      .bt_ordered_label_text(spec$label_parts$components,spec$metadata,"share",
        factor=record$factor,level=tokens[[i]],slice=if(sliced) slice)
    }
    .bt_label_parts(components=text,
      formula_parameter=spec$label_parts$formula_parameter,selector=names[[i]])
  })
  structure(stats::setNames(parts,names), levels = tokens)
}

# Label parts of the total coordinates of the ordered prior 'prior' of
# 'parameter', named by the coordinates, and of the inclusion indicator of a
# mixture total (the label of the total without its slice index and with the
# inclusion suffix). With 'legacy', the labels the quantities had before their
# labels followed the label grammar (`f_ordered_total`, `f:g_ordered_total[1]`,
# `f_ordered_total_indicator`), which the catalog keeps as aliases.
.bt_ordered_total_label_parts <- function(parameter, prior,
                                          term = .bt_label_parts_term(parameter, prior),
                                          legacy = FALSE){

  metadata <- .prior_ordered_metadata(prior)
  names <- .prior_ordered_total_monitor_names(prior, parameter)
  slices <- .bt_ordered_total_slices(metadata)
  stem <- paste0(paste(term$components, collapse = ":"), "_ordered_total")
  text <- if(legacy){
    paste0(stem, if(!is.null(slices)) paste0("[", slices, "]"))
  }else{
    vapply(seq_along(names), function(i){
      .bt_ordered_label_text(term$components, metadata, "total", slice = slices[i])
    }, character(1))
  }
  if(is.prior.mixture(prior$total)){
    names <- c(names, paste0(.prior_ordered_total_name(parameter), "_indicator"))
    text <- c(text, if(legacy){
      paste0(stem, "_indicator")
    }else{
      .bt_ordered_label_text(term$components, metadata, "total", inclusion = TRUE)
    })
  }
  stats::setNames(lapply(seq_along(names), function(i){
    .bt_label_parts(components = text[[i]], formula_parameter = term$formula_parameter,
                    selector = names[[i]])
  }), names)
}

# The names of the bound ordered terms of 'prior_list' whose totals and shares
# are public quantities. The ordered SD components of a random block are
# backend coordinates and keep their backend labels.
.bt_ordered_fixed_terms <- function(prior_list){

  keep <- vapply(prior_list, function(prior){
    is.prior.ordered(prior) && !is.null(attr(prior, "ordered_metadata", exact = TRUE)) &&
      !.bt_is_random_effect_prior(prior)
  }, logical(1))

  names(prior_list)[keep]
}

# The label parts of the totals (and total indicators) of the ordered terms in
# 'prior_list', named by their coordinates.
.bt_ordered_total_label_parts_list <- function(prior_list){

  parts <- list()
  for(parameter in .bt_ordered_fixed_terms(prior_list)){
    parts <- c(parts, .bt_ordered_total_label_parts(parameter, prior_list[[parameter]]))
  }

  parts
}

# The labels the totals, total indicators, and shares of the ordered terms in
# 'prior_list' had before their labels followed the label grammar, as label
# parts named by their coordinates (the catalog keeps them as aliases).
.bt_ordered_legacy_label_parts <- function(prior_list){

  parts <- list()
  for(parameter in .bt_ordered_fixed_terms(prior_list)){
    prior <- prior_list[[parameter]]
    parts <- c(parts, .bt_ordered_total_label_parts(parameter, prior, legacy = TRUE))
    spec <- .bt_ordered_spec(parameter, prior)
    for(record in spec$allocations){
      if(identical(record$spec$type, "dirichlet")){
        parts <- c(parts, .bt_ordered_allocation_parts(spec, record, legacy = TRUE))
      }
    }
  }

  parts
}

.bt_ordered_allocation_shares <- function(parameter, prior, model_samples){

  spec <- .bt_ordered_spec(parameter,prior)

  columns <- list()
  parts   <- list()
  for(record in spec$allocations){
    if(!identical(record$spec$type, "dirichlet")){
      next
    }
    share <- .bt_ordered_allocation_values(record,model_samples)
    if(is.null(share)){
      .bt_stop_refit_required(
        "The allocation of the ordered prior '", parameter,
        "' was not monitored; refit the model with this version of BayesTools.",
        class = "BayesTools_ordered_metadata_unavailable"
      )
    }
    share_parts <- .bt_ordered_allocation_parts(spec,record)
    share_names <- names(share_parts)
    colnames(share) <- share_names
    columns[[length(columns) + 1L]] <- share
    parts <- c(parts,share_parts)
  }

  list(
    samples = if(length(columns) > 0L) do.call(cbind, columns) else model_samples[, 0L, drop = FALSE],
    parts   = parts
  )
}

# Retained source rows stay on the fitted scale, independently of transformations
# of effect columns. Models include the complete declared ensemble, even when no
# row was selected from a model.
.bt_ordered_source_model <- function(prior, parameter, fit = NULL){

  if(is.prior.ordered(prior)){
    return(if(is.null(fit)) .bt_ordered_spec(parameter, prior) else JAGS_ordered_parameter_spec(fit,parameter)[[parameter]])
  }
  list(parameter = parameter, prior = prior,
    parameterization = if(is.null(prior) || .posterior_atoms_is_zero_point(prior)) "absent" else "unavailable")
}

.bt_ordered_source_rows <- function(spec, model_samples, rows){

  if(!is.null(spec$parameterization)) return(matrix(numeric(), length(rows), 0L))
  draws <- model_samples[rows, , drop = FALSE]
  totals <- .bt_ordered_total_values(spec, .bt_deterministic_lookup(draws))
  if(is.null(totals)) .bt_ordered_stop(paste0("Ordered total sources for '", spec$parameter,
    "' are unavailable. Include its declared source coordinates."))
  colnames(totals) <- spec$total_names
  allocations <- lapply(spec$allocations, function(record){
    values <- .bt_ordered_allocation_values(record, draws)
    if(is.null(values)) .bt_ordered_stop(paste0("Ordered allocation sources for '", spec$parameter,
      "' are unavailable. Include its declared source coordinates."))
    colnames(values) <- record$coordinates
    values
  })
  sources <- c(list(totals), allocations)
  if(!is.null(spec$total_node)){
    indicator <- spec$total_node$spec$indicator
    if(!indicator %in% colnames(draws)) .bt_ordered_stop(paste0("Ordered total indicator '", indicator,
      "' is unavailable. Include its declared source coordinate."))
    sources <- c(sources, list(draws[, indicator, drop = FALSE]))
  }
  do.call(cbind, sources)
}

.bt_ordered_source_new <- function(parameter, models, sources, model, draw_index, posterior_pair = NULL){

  columns <- unique(unlist(lapply(sources, colnames), use.names = FALSE))
  values <- matrix(NA_real_, length(model), length(columns), dimnames = list(NULL, columns))
  for(i in seq_along(models)){
    rows <- which(model == i)
    if(!length(rows)) next
    source <- sources[[i]]
    values[rows, colnames(source)] <- source
  }
  if(is.null(posterior_pair)){
    if(length(models) != 1L) stop("Ordered model probabilities require an explicit producer owner.", call. = FALSE)
    posterior_pair <- .model_probability_pair(1, 0, "posterior", "ordinary")
  }
  .model_probability_validate(posterior_pair$probabilities, posterior_pair$logs, posterior_pair$declaration)
  structure(list(version = 2L, parameter = parameter, models = models,
    model = as.integer(model), draw_index = as.integer(draw_index), primitives = values,
    model_probabilities = posterior_pair$probabilities, model_log_probabilities = posterior_pair$logs,
    model_probability_declaration = posterior_pair$declaration),
    class = c("BayesTools_ordered_source", "list"))
}

.bt_ordered_source_validate <- function(value){

  if(!inherits(value, "BayesTools_ordered_source") || !is.list(value) ||
     !identical(value$version, 2L) || !is.character(value$parameter) || length(value$parameter) != 1L ||
     !is.list(value$models) || !.bt_meta_is_index(value$model) || !.bt_meta_is_index(value$draw_index) ||
     any(value$model > length(value$models)) || !is.matrix(value$primitives) || !is.numeric(value$primitives) ||
     (ncol(value$primitives) > 0L && (is.null(colnames(value$primitives)) || anyDuplicated(colnames(value$primitives)))) ||
     length(value$model) != nrow(value$primitives) || length(value$draw_index) != nrow(value$primitives)){
    parameter <- if(is.list(value) && is.character(value$parameter) &&
                    length(value$parameter) == 1L && !is.na(value$parameter)){
      paste0(" for parameter '", value$parameter, "'")
    }else{
      ""
    }
    return(paste0(
      "it must contain validated row-aligned ordered primitive sources and model provenance",
      parameter, "; recreate the draws from the fitted models"
    ))
  }
  .model_probability_validate(value$model_probabilities, value$model_log_probabilities, value$model_probability_declaration, normalized = TRUE)
  if(!is.null(value$conditioning)){
    declaration <- value$conditioning$model_probability_declaration
    .model_probability_validate(value$conditioning$prior_probs,value$conditioning$log_prior_probs,declaration$prior)
    .model_probability_validate(value$conditioning$post_probs,value$conditioning$log_post_probs,declaration$posterior)
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
    if(ncol(value$primitives) == 0L) return("its allocated ordered primitive coordinates are missing")
    if(!all(spec$total_names %in% colnames(value$primitives))) return("its declared total coordinates are missing")
    if(!is.null(spec$total_node) && !spec$total_node$spec$indicator %in% colnames(value$primitives)){
      return("its declared total indicator coordinate is missing")
    }
    for(record in spec$allocations){
      if(!all(record$coordinates %in% colnames(value$primitives))) return("its declared allocation coordinates are missing")
      if(length(rows)) .bt_ordered_allocation_values(record, value$primitives[rows, , drop = FALSE])
    }
  }
  if(!is.null(value$projection_design) && (!is.matrix(value$projection_design) || !is.numeric(value$projection_design) ||
     any(!is.finite(value$projection_design)) || is.null(colnames(value$projection_design)))){
    return("its projection design must contain finite named fitted-coordinate weights")
  }
  if(!is.null(value$projection_context) && (!is.list(value$projection_context) ||
     !is.matrix(value$projection_context$primitives) || nrow(value$projection_context$primitives)!=nrow(value$primitives))){
    return("its projection context must have one primitive row per retained draw")
  }
  NULL
}

.bt_ordered_source_subset <- function(source, rows){

  source$model <- source$model[rows]
  source$draw_index <- source$draw_index[rows]
  source$primitives <- source$primitives[rows, , drop = FALSE]
  if(!is.null(source$projection_context)) source$projection_context$primitives <- source$projection_context$primitives[rows,,drop=FALSE]
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
  parts <- .bt_ordered_total_label_parts(template$parameter, template$prior)[template$total_names]
  columns <- list(totals)
  undefined <- character()
  for(factor in template$metadata$ordered_terms){
    records <- unlist(lapply(ordered, function(spec){
      spec$allocations[vapply(spec$allocations, function(record) identical(record$factor, factor), logical(1))]
    }), recursive = FALSE)
    if(!any(vapply(records, function(record) identical(record$spec$type, "dirichlet"), logical(1)))) next
    sliced <- template$metadata$theta_dim > 1L && any(vapply(records, function(record) !is.na(record$slice), logical(1)))
    slices <- if(sliced) seq_len(template$metadata$theta_dim) else 1L
    for(slice in slices){
      record <- records[[1L]]
      share_parts <- .bt_ordered_allocation_parts(template,record,sliced=sliced,slice=slice)
      names <- names(share_parts)
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
      parts <- c(parts,share_parts)
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
  quantity_ids <- rep("",ncol(samples))
  ordered <- source$models[vapply(source$models,function(spec) is.null(spec$parameterization),logical(1))]
  # the totals and the shares of the models are catalog quantities of their fits
  for(spec in rev(ordered)){
    for(quantities in list(spec$total_quantities,spec$share_quantities)){
      if(is.null(quantities)) next
      matches <- match(quantities$canonical_name,colnames(samples))
      valid <- !is.na(matches) & !is.na(quantities$quantity_id)
      quantity_ids[matches[valid]] <- quantities$quantity_id[valid]
    }
  }
  quantities <- .bt_draws_quantity_table(colnames(samples), quantity_ids,
    rep(list(character()), ncol(samples)), rep(list(numeric()), ncol(samples)), unname(display$parts))
  samples <- .bt_meta_set(samples, "quantities", quantities)
  samples <- .bt_meta_set(samples, "undefined_draws", display$undefined_draws)
  samples
}

.bt_ordered_fitted_scale_footnote <- function(){

  paste0("Ordered-factor totals and allocations are summarized on the fitted ",
    "(standardized) scale. Use 'transform_factors = TRUE' for the original-scale level effects.")
}

