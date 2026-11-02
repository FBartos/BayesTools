# Ordered quantities of the parameter catalog -------------------------------
#
# The total of an ordered term is the catalog quantity of its total
# coordinate (one per slice): it has the support and the prior density of the
# `total` prior. Each share of a random (Dirichlet) allocation is a derived
# catalog quantity of the term, `<parameter>_ordered_allocation...`: the
# normalized gamma coordinates of its allocation, with the support [0, 1] and
# the Beta marginal of the Dirichlet prior. Shares of fixed allocations are
# constants of the recipe and have no quantity. The ordered SD components of a
# random block keep their backend coordinates.

# The term of 'prior_list' whose total coordinates include 'coordinate' (NULL
# when none does).
.bt_ordered_total_owner <- function(prior_list, coordinate){

  for(parameter in .bt_ordered_fixed_terms(prior_list)){
    if(coordinate %in% .prior_ordered_total_monitor_names(prior_list[[parameter]], parameter)){
      return(parameter)
    }
  }

  NULL
}

# The extraction key of the share 'index' of the allocation record 'record' of
# the ordered term 'parameter'. The record is recomputed from the prior list.
.bt_ordered_share_key <- function(parameter, record, index, level){

  list(
    type           = "ordered_allocation",
    dependencies   = record$gamma_coordinates,
    parameter      = parameter,
    allocation_key = record$key,
    index          = index,
    level          = level
  )
}

.bt_ordered_share_record <- function(prior_list, key){

  prior <- prior_list[[key$parameter]]
  record <- if(is.prior.ordered(prior) && !is.null(attr(prior, "ordered_metadata", exact = TRUE))){
    .bt_ordered_spec(key$parameter, prior)$allocations[[key$allocation_key]]
  }
  if(is.null(record) || !identical(record$spec$type, "dirichlet") ||
     !identical(record$gamma_coordinates, key$dependencies)){
    .bt_stop_refit_required(
      "The allocation of the ordered share quantity of '", key$parameter,
      "' is not part of the fitted prior list. Refit the model with this version of BayesTools.",
      class = "BayesTools_ordered_metadata_unavailable"
    )
  }

  record
}

# The catalog quantities of the shares of the random allocations of the ordered
# terms in 'prior_list': one per share, named by the selector of the display
# row (`mu_f_ordered_allocation[mid]`), derived from the gamma coordinates of
# its allocation (unavailable when a coordinate is). A shared ('id') allocation
# has the shares of every term that uses it. An allocation whose gamma
# coordinates are not coordinates of the fit (they were not monitored) has no
# share quantities.
.bt_ordered_share_quantities <- function(coordinates, prior_list){

  rows <- list()
  for(parameter in .bt_ordered_fixed_terms(prior_list)){
    spec <- .bt_ordered_spec(parameter, prior_list[[parameter]])
    for(record in spec$allocations){
      if(!identical(record$spec$type, "dirichlet")){
        next
      }
      gamma <- match(record$gamma_coordinates, coordinates$coordinate_name)
      if(anyNA(gamma)){
        next
      }
      parts <- .bt_ordered_allocation_parts(spec, record)
      status <- if(any(coordinates$monitor_status[gamma] == "unavailable")) "unavailable" else "derived"
      formula_parameter <- spec$label_parts$formula_parameter
      for(index in seq_along(parts)){
        rows[[length(rows) + 1L]] <- .bt_parameter_catalog_quantity(
          canonical_name    = names(parts)[[index]],
          namespace         = if(nzchar(formula_parameter)) formula_parameter else "model",
          role              = "parameter",
          formula_parameter = formula_parameter,
          label_parts       = parts[[index]],
          fitted_scale      = "unitless",
          display_scale     = "unitless",
          status            = status,
          source_type       = "composite",
          extraction_key    = .bt_ordered_share_key(
            parameter, record, index, attr(parts, "levels", exact = TRUE)[[index]]
          )
        )
      }
    }
  }
  if(length(rows) == 0L){
    return(.bt_parameter_catalog_empty_quantities())
  }

  .bt_parameter_catalog_bind_quantities(rows)
}

# The draws of a share: the normalized gamma coordinates of its allocation.
.bt_parameter_draw_ordered_allocation <- function(fit, key, model_samples){

  record <- .bt_ordered_share_record(attr(fit, "prior_list", exact = TRUE), key)
  values <- .bt_ordered_allocation_values(record, model_samples)
  if(is.null(values)){
    .bt_ordered_stop(paste0(
      "Ordered allocation sources for '", key$parameter,
      "' are unavailable in the draws. Include its declared source coordinates."
    ))
  }

  values[, key$index]
}

# The gate plan of a share: it has no atoms, and its inclusion event is that
# of the total of its term (the draws of the alternative components of a
# null/alternative mixture total); NULL for terms without such an event.
.bt_ordered_share_plan <- function(fit, key){

  prior_list <- attr(fit, "prior_list", exact = TRUE)
  total <- prior_list[[key$parameter]]$total
  if(!is.prior.mixture(total) ||
     !all(attr(total, "components", exact = TRUE) %in% c("null", "alternative"))){
    return(NULL)
  }

  list(
    kind        = "ordered_share",
    prior_list  = prior_list,
    chain_gates = character(),
    event_gates = key$parameter,
    event_rule  = "AND"
  )
}

# The per-draw states of the plan of a share: never on an atom, always
# defined. The inclusion event (the total of its term is in its alternative)
# is evaluated on request only, as it needs the indicator of the total: the
# draws and the prior of an unconditional share do not.
.bt_ordered_share_states <- function(plan, n, model_samples, with_event = FALSE){

  list(
    atom    = rep(NA_real_, n),
    defined = rep(TRUE, n),
    event   = if(isTRUE(with_event)){
      .condition_event_posterior_mask(
        .condition_event(plan$prior_list, plan$event_gates, "AND"),
        plan$prior_list, model_samples
      )
    },
    known   = TRUE
  )
}

# The declared support of an ordered quantity that is not an ordinary
# coordinate: the support of the total prior of a total coordinate, and [0, 1]
# for a share. NULL for other quantities.
.bt_ordered_quantity_support <- function(object, key){

  if(identical(key$type, "ordered_allocation")){
    return(.posterior_support_new(c(0, 1), source = "catalog", type = "interval"))
  }
  prior_list <- attr(object, "prior_list", exact = TRUE)
  owner <- if(identical(key$type, "coordinate")){
    .bt_ordered_total_owner(prior_list, key$dependencies)
  }
  if(is.null(owner)){
    return(NULL)
  }

  .posterior_support_from_prior(prior_list[[owner]]$total, source = "catalog")
}

# The declared support of an ordered quantity within its inclusion event: the
# catalog support of the quantity under the priors restricted to the event
# (.condition_event_model_options(): the continuous components of the total
# priors), united over the restricted models. The unconditional support keeps
# the bounds and the point of the null component, which the draws of the event
# never take.
.bt_ordered_conditioned_support <- function(fit, plan, quantity){

  prior_list <- plan$prior_list
  options <- .condition_event_model_options(
    prior_list, .condition_event(prior_list, plan$event_gates, "AND")
  )
  .posterior_support_union(lapply(options$prior_lists, function(restricted){
    restricted_fit <- fit
    attr(restricted_fit, "prior_list") <- restricted
    .bt_parameter_catalog_quantity_support(restricted_fit, quantity)
  }), source = "catalog")
}

# The prior density of the total coordinate of an ordered term (the density of
# its `total` prior, also for every slice) and of a share (the Beta marginal of
# its Dirichlet allocation). With 'conditional', the density of a total is
# that of the total prior restricted to its inclusion event: the mixture of the
# continuous ('alternative') components (parameter_mixed_posterior(conditional
# = TRUE)); a share is independent of the total.
.bt_parameter_prior_density_ordered <- function(object, quantity, key, n_grid,
                                                tail_prob, conditional = FALSE){

  prior_list <- attr(object, "prior_list", exact = TRUE)
  if(identical(key$type, "ordered_allocation")){
    record <- .bt_ordered_share_record(prior_list, key)
    source <- .bt_parameter_prior_density_simplex_marginal(
      prior("dirichlet", list(alpha = record$spec$alpha)), key$index
    )
  }else{
    source <- prior_list[[.bt_ordered_total_owner(prior_list, key$dependencies)]]$total
    if(isTRUE(conditional)){
      source <- .bt_ordered_total_conditioned(
        prior_list, .bt_ordered_total_owner(prior_list, key$dependencies),
        quantity$canonical_name
      )
    }
  }
  out <- .bt_parameter_prior_density_scalar(source, n_grid = n_grid, tail_prob = tail_prob)
  if(!is.null(out)){
    attr(out, "parameter_prior_density") <- c(
      list(quantity_id = quantity$quantity_id, source = "fitted_parameter_map"),
      if(isTRUE(conditional)) list(conditional = "inclusion")
    )
  }

  out
}

# The total prior of the ordered term 'parameter' within its inclusion event
# (.condition_event_model_options()): the continuous component, or the mixture
# of the continuous components with their prior probabilities renormalized.
.bt_ordered_total_conditioned <- function(prior_list, parameter, name){

  options <- .condition_event_model_options(
    prior_list, .condition_event(prior_list, parameter, "AND")
  )
  restricted <- lapply(options$prior_lists, function(restricted_list){
    restricted_list[[parameter]]$total
  })
  if(length(restricted) == 0L){
    stop(
      "The inclusion event of '", name, "' is unavailable: it has prior ",
      "probability zero. Use 'conditional = FALSE'.",
      call. = FALSE
    )
  }
  if(length(restricted) == 1L){
    return(restricted[[1L]])
  }

  prior_mixture(restricted, is_null = rep(FALSE, length(restricted)))
}

# The rows of a catalog quantities table that are totals or shares of ordered
# terms.
.bt_ordered_quantity_rows <- function(quantities, prior_list){

  vapply(seq_len(nrow(quantities)), function(i){
    key <- quantities$extraction_key[[i]]
    identical(quantities$provider[[i]], "BayesTools") &&
      (identical(key$type, "ordered_allocation") ||
         (identical(key$type, "coordinate") &&
            !is.null(.bt_ordered_total_owner(prior_list, key$dependencies))))
  }, logical(1))
}

# The ordered term that owns an ordered quantity (a total or a share).
.bt_ordered_quantity_owner <- function(quantity, prior_list){

  key <- quantity$extraction_key[[1L]]
  if(identical(key$type, "ordered_allocation")){
    return(key$parameter)
  }

  .bt_ordered_total_owner(prior_list, key$dependencies)
}

# The root with which hypotheses name an ordered quantity inside the engine:
# the canonical name with its trailing bracket (the level of a share, the
# slice of a total) in braces, so that the name is no level reference.
.bt_ordered_hypothesis_root <- function(quantity){

  name <- quantity$canonical_name
  key <- quantity$extraction_key[[1L]]
  root <- if(identical(key$type, "ordered_allocation")){
    paste0(substr(name, 1L, nchar(name) - nchar(key$level) - 2L), "{", key$level, "}")
  }else if(grepl("]$", name)){
    sub("\\[([^][]*)\\]$", "{\\1}", name)
  }else{
    name
  }
  if(!.bt_hypothesis_root_round_trips(root)){
    root <- paste0(name, "{}")
  }

  root
}

