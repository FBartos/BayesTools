# R-squared budget quantities of the parameter catalog -----------------------
#
# The total of a budget is the catalog quantity of its total coordinate
# (`<parameter>_<name>`), displayed as `<name>{total}`: it has the support and
# the prior density of its Beta prior in the model. Each retained term of a
# budget that retains two or more terms has a derived share quantity,
# `<parameter>_<name>_allocation[<term>]` displayed as `<name>[share: <term>]`:
# the normalized gamma coordinates of the allocation summed over the term's
# allocation units (its coded coefficients under the kernel laws), with the
# support [0, 1] and the Beta marginal of the Dirichlet allocation. The
# coordinates of the allocation are backend coordinates. A budget that retains
# one term allocates it the whole budget and has no share quantities.

# The budgets of a bound prior list: for every budget total, its name, the
# formula parameter, the total and allocation nodes, the gamma coordinates of
# the allocation, and the retained terms with their allocation units.
.bt_R2_budget_records <- function(prior_list){

  records <- list()
  if(length(prior_list) == 0L){
    return(records)
  }
  totals <- names(prior_list)[vapply(prior_list, function(prior){
    identical(attr(prior, "R2_budget", exact = TRUE)$role, "total")
  }, logical(1))]
  for(total_name in totals){
    terms <- names(prior_list)[vapply(prior_list, function(prior){
      is.prior.R2(prior) && isTRUE(prior$bound) && identical(prior$total_name, total_name)
    }, logical(1))]
    if(length(terms) == 0L){
      next
    }
    first <- prior_list[[terms[[1L]]]]
    allocation_name <- first$allocation_name
    records[[total_name]] <- list(
      name              = first$budget$name,
      parameter         = first$parameter,
      total_name        = total_name,
      allocation_name   = if(is.na(allocation_name)) NULL else allocation_name,
      n_units           = first$n_units,
      gamma_coordinates = if(!is.na(allocation_name)){
        paste0(.JAGS_prior_dirichlet_eta_name(allocation_name), "[", seq_len(first$n_units), "]")
      },
      terms             = lapply(stats::setNames(terms, terms), function(term){
        list(
          parameter = term,
          label     = gsub("__xXx__", ":", prior_list[[term]]$term, fixed = TRUE),
          units     = prior_list[[term]]$units
        )
      })
    )
  }

  records
}

# The allocation nodes of the budgets of a bound prior list (backend
# coordinates).
.bt_R2_allocation_names <- function(prior_list){

  names(prior_list)[vapply(prior_list, function(prior){
    identical(attr(prior, "R2_budget", exact = TRUE)$role, "allocation")
  }, logical(1))]
}

# Label parts of the totals of the budgets of a bound prior list, named by
# their coordinates.
.bt_R2_total_label_parts_list <- function(prior_list){

  records <- .bt_R2_budget_records(prior_list)
  stats::setNames(lapply(records, function(record){
    .bt_label_parts(
      components        = paste0(record$name, "{total}"),
      formula_parameter = record$parameter,
      selector          = record$total_name
    )
  }), names(records))
}

# Label parts of the prior-list entry 'parameter' when it is the total of a
# budget (NULL otherwise), as the catalog labels its quantity.
.bt_R2_total_prior_label_parts <- function(parameter, prior){

  budget <- attr(prior, "R2_budget", exact = TRUE)
  if(!identical(budget$role, "total") || is.null(budget$name)){
    return(NULL)
  }

  .bt_label_parts(
    components        = paste0(budget$name, "{total}"),
    formula_parameter = .bt_label_formula_parameter(prior),
    selector          = parameter
  )
}

# Canonical name and label parts of the share of 'term' in 'record'.
.bt_R2_share_name <- function(record, term){
  paste0(record$total_name, "_allocation[", term$label, "]")
}
.bt_R2_share_label_parts <- function(record, term){
  .bt_label_parts(
    components        = paste0(record$name, "[share: ", term$label, "]"),
    formula_parameter = record$parameter,
    selector          = .bt_R2_share_name(record, term)
  )
}

# The share draws of the terms of 'record': the normalized gamma coordinates
# summed over each term's units (NULL when the coordinates are not in
# 'model_samples').
.bt_R2_share_values <- function(record, model_samples){

  if(is.null(record$gamma_coordinates) ||
     !all(record$gamma_coordinates %in% colnames(model_samples))){
    return(NULL)
  }
  eta <- model_samples[, record$gamma_coordinates, drop = FALSE]
  phi <- eta / rowSums(eta)
  values <- vapply(record$terms, function(term){
    rowSums(phi[, term$units, drop = FALSE])
  }, numeric(nrow(phi)))
  matrix(values, nrow = nrow(phi), dimnames = list(NULL, vapply(record$terms, function(term){
    .bt_R2_share_name(record, term)
  }, character(1))))
}

# The catalog quantities of the shares of the budgets in 'prior_list'.
.bt_R2_share_quantities <- function(coordinates, prior_list){

  rows <- list()
  for(record in .bt_R2_budget_records(prior_list)){
    if(length(record$terms) < 2L){
      next
    }
    gamma <- match(record$gamma_coordinates, coordinates$coordinate_name)
    if(anyNA(gamma)){
      next
    }
    status <- if(any(coordinates$monitor_status[gamma] == "unavailable")) "unavailable" else "derived"
    for(term in record$terms){
      rows[[length(rows) + 1L]] <- .bt_parameter_catalog_quantity(
        canonical_name    = .bt_R2_share_name(record, term),
        namespace         = record$parameter,
        role              = "parameter",
        formula_parameter = record$parameter,
        label_parts       = .bt_R2_share_label_parts(record, term),
        fitted_scale      = "unitless",
        display_scale     = "unitless",
        status            = status,
        source_type       = "composite",
        extraction_key    = list(
          type         = "R2_share",
          dependencies = record$gamma_coordinates,
          total        = record$total_name,
          term         = term$parameter,
          units        = as.integer(term$units)
        )
      )
    }
  }
  if(length(rows) == 0L){
    return(.bt_parameter_catalog_empty_quantities())
  }

  .bt_parameter_catalog_bind_quantities(rows)
}

# The budget record of a share key, recomputed from the prior list.
.bt_R2_share_record <- function(prior_list, key){

  record <- .bt_R2_budget_records(prior_list)[[key$total]]
  term <- if(!is.null(record)) record$terms[[key$term]]
  if(is.null(term) || !identical(record$gamma_coordinates, key$dependencies) ||
     !identical(as.integer(term$units), as.integer(key$units))){
    .bt_stop_refit_required(
      "The R-squared budget of the share quantity of '", key$term,
      "' is not part of the fitted prior list. Refit the model with this version of BayesTools.",
      class = "BayesTools_R2_unavailable"
    )
  }

  list(record = record, term = term)
}

.bt_parameter_catalog_valid_R2_share_key <- function(key){

  is.character(key$dependencies) && length(key$dependencies) >= 2L &&
    is.character(key$total) && length(key$total) == 1L && !is.na(key$total) && nzchar(key$total) &&
    is.character(key$term) && length(key$term) == 1L && !is.na(key$term) && nzchar(key$term) &&
    is.numeric(key$units) && length(key$units) >= 1L && !anyNA(key$units) &&
    all(key$units >= 1L & key$units <= length(key$dependencies))
}

# The draws of a share.
.bt_parameter_draw_R2_share <- function(fit, key, model_samples){

  found <- .bt_R2_share_record(attr(fit, "prior_list", exact = TRUE), key)
  if(!all(key$dependencies %in% colnames(model_samples))){
    .bt_R2_unavailable(paste0(
      "The allocation coordinates of the R-squared budget share of '", key$term,
      "' are unavailable in the draws. Include its declared source coordinates."
    ))
  }
  eta <- model_samples[, key$dependencies, drop = FALSE]

  rowSums((eta / rowSums(eta))[, key$units, drop = FALSE])
}

# The prior density of a share: the Beta marginal of its units under the
# Dirichlet allocation of the model.
.bt_parameter_prior_density_R2_share <- function(object, quantity, key, n_grid, tail_prob){

  prior_list <- attr(object, "prior_list", exact = TRUE)
  found <- .bt_R2_share_record(prior_list, key)
  alpha <- prior_list[[found$record$allocation_name]]$parameters[["alpha"]]
  share <- sum(alpha[key$units])
  out <- .bt_parameter_prior_density_scalar(
    prior("beta", list(alpha = share, beta = sum(alpha) - share)),
    n_grid = n_grid, tail_prob = tail_prob
  )
  if(!is.null(out)){
    attr(out, "parameter_prior_density") <- list(
      quantity_id = quantity$quantity_id, source = "fitted_parameter_map"
    )
  }

  out
}

# The rows of a catalog quantities table that are totals or shares of
# R-squared budgets.
.bt_R2_quantity_rows <- function(quantities, prior_list){

  totals <- names(.bt_R2_budget_records(prior_list))
  vapply(seq_len(nrow(quantities)), function(i){
    key <- quantities$extraction_key[[i]]
    identical(quantities$provider[[i]], "BayesTools") &&
      (identical(key$type, "R2_share") ||
         (identical(key$type, "coordinate") && length(key$dependencies) == 1L &&
            key$dependencies %in% totals))
  }, logical(1))
}

# The root with which hypotheses name a budget quantity inside the engine:
# the canonical name of a total, and for a share the canonical name with its
# term in braces, so that the name is no level reference.
.bt_R2_hypothesis_root <- function(quantity){

  name <- quantity$canonical_name
  key <- quantity$extraction_key[[1L]]
  if(!identical(key$type, "R2_share")){
    return(name)
  }
  root <- sub("\\[([^][]*)\\]$", "{\\1}", name)
  if(!.bt_hypothesis_root_round_trips(root)){
    root <- paste0(sub("\\[([^][]*)\\]$", "", name), "{", gsub("[^A-Za-z0-9_.]", "_", key$term), "}")
  }

  root
}

# Estimates tables show the share of each retained term of a budget after its
# total, in place of the allocation coordinates. Returns the samples and the
# label parts of the created share columns.
.bt_JAGS_estimates_R2_shares <- function(model_samples, raw_model_samples, prior_list){

  created <- list()
  for(record in .bt_R2_budget_records(prior_list)){
    if(is.null(record$allocation_name)){
      next
    }
    allocation_columns <- paste0(record$allocation_name, "[", seq_len(record$n_units), "]")
    present <- colnames(model_samples) %in% allocation_columns
    if(!any(present)){
      next
    }
    position <- min(which(present))
    model_samples <- model_samples[, !present, drop = FALSE]
    if(length(record$terms) < 2L || !all(allocation_columns %in% colnames(raw_model_samples))){
      next
    }
    phi <- raw_model_samples[, allocation_columns, drop = FALSE]
    shares <- vapply(record$terms, function(term){
      rowSums(phi[, term$units, drop = FALSE])
    }, numeric(nrow(phi)))
    shares <- matrix(shares, nrow = nrow(phi))
    colnames(shares) <- vapply(record$terms, function(term) .bt_R2_share_name(record, term), character(1))
    if(record$total_name %in% colnames(model_samples)){
      position <- match(record$total_name, colnames(model_samples)) + 1L
    }
    model_samples <- cbind(
      if(position > 1L) model_samples[, seq_len(position - 1L), drop = FALSE],
      shares,
      if(position <= ncol(model_samples)) model_samples[, position:ncol(model_samples), drop = FALSE]
    )
    for(term in record$terms){
      created[[.bt_R2_share_name(record, term)]] <- .bt_R2_share_label_parts(record, term)
    }
  }

  list(model_samples = model_samples, created = created)
}
