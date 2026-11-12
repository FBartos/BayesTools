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
      core              = if(is.prior.R2_budget(prior_list[[total_name]])) prior_list[[total_name]],
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
          key       = prior_list[[term]]$term,
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

  c(names(prior_list)[vapply(prior_list, function(prior){
    identical(attr(prior, "R2_budget", exact = TRUE)$role, "allocation")
  }, logical(1))],
    unlist(lapply(prior_list[vapply(prior_list, is.prior.R2_budget, logical(1))], function(core){
      c(if(core$n_units > 1L) core$allocation_name, if(core$gates$K > 0L) paste0(core$total_name, c("_gate", "_open")))
    }), use.names = FALSE))
}

# Label parts of the totals of the budgets of a bound prior list, named by
# their coordinates.
.bt_R2_total_label_parts_list <- function(prior_list){

  records <- .bt_R2_budget_records(prior_list)
  parts <- stats::setNames(lapply(records, function(record){
    .bt_label_parts(
      components        = paste0(record$name, "{total}"),
      formula_parameter = record$parameter,
      selector          = record$total_name
    )
  }), names(records))
  # the budget indicator of a gated or spiked budget is the inclusion of its
  # total; a shared Beta rate of its gates is its inclusion rate
  for(record in records){
    nodes <- if(!is.null(record$core)) .bt_R2_core_nodes(record$core)
    if(!is.null(nodes$indicator)){
      parts[[nodes$indicator]] <- .bt_label_parts(
        components        = paste0(record$name, "{total} (inclusion)"),
        formula_parameter = record$parameter,
        selector          = nodes$indicator
      )
    }
    if(!is.null(nodes$rate)){
      parts[[nodes$rate]] <- .bt_label_parts(
        components        = paste0(record$name, "{inclusion rate}"),
        formula_parameter = record$parameter,
        selector          = nodes$rate
      )
    }
  }

  parts
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
    gating <- .bt_R2_record_gating(record)
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
          dependencies = c(record$gamma_coordinates, gating),
          total        = record$total_name,
          term         = term$parameter,
          units        = as.integer(term$units)
        )
      )
    }
  }
  # the inclusion of every gated term: its gate (times the budget indicator
  # of a spike-and-slab total)
  for(record in .bt_R2_budget_records(prior_list)){
    core <- record$core
    if(is.null(core) || core$gates$K == 0L){
      next
    }
    nodes <- .bt_R2_core_nodes(core)
    for(term in record$terms){
      gate <- core$term_gate[[term$key]]
      if(gate == 0L){
        next
      }
      dependencies <- nodes$open[[gate]]
      position <- match(dependencies, coordinates$coordinate_name)
      if(anyNA(position)){
        next
      }
      rows[[length(rows) + 1L]] <- .bt_parameter_catalog_quantity(
        canonical_name    = paste0(record$total_name, "_gate[", term$label, "]"),
        namespace         = record$parameter,
        role              = "parameter",
        formula_parameter = record$parameter,
        label_parts       = .bt_label_parts(
          components        = paste0(term$label, " (inclusion)"),
          formula_parameter = record$parameter,
          selector          = paste0(record$total_name, "_gate[", term$label, "]")
        ),
        fitted_scale      = "unitless",
        display_scale     = "unitless",
        status            = if(any(coordinates$monitor_status[position] == "unavailable")) "unavailable" else "derived",
        source_type       = "composite",
        extraction_key    = list(
          type         = "R2_open",
          dependencies = dependencies,
          total        = record$total_name,
          term         = term$parameter
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
  if(is.null(term) || !identical(c(record$gamma_coordinates, .bt_R2_record_gating(record)), key$dependencies) ||
     !identical(as.integer(term$units), as.integer(key$units))){
    .bt_stop_refit_required(
      "The R-squared budget of the share quantity of '", key$term,
      "' is not part of the fitted prior list. Refit the model with this version of BayesTools.",
      class = "BayesTools_R2_unavailable"
    )
  }

  list(record = record, term = term,
       alpha = if(is.null(record$core)) prior_list[[record$allocation_name]]$parameters[["alpha"]])
}

# The prior law of a share: the Beta marginal of its units under the
# Dirichlet allocation of the model, or for a gated budget the mixture over
# the open configurations that contain the term (given that it is open).
.bt_R2_share_prior_law <- function(found){

  record <- found$record
  if(!is.null(record$core)){
    return(.bt_R2_core_share_law(record$core, found$term$key))
  }
  alpha <- found$alpha
  # the complement is summed directly (no cancellation in sum(alpha) - share)
  prior("beta", list(alpha = sum(alpha[found$term$units]), beta = sum(alpha[-found$term$units])))
}

.bt_parameter_catalog_valid_R2_open_key <- function(key){

  is.character(key$dependencies) && length(key$dependencies) %in% 1:2 &&
    is.character(key$total) && length(key$total) == 1L && !is.na(key$total) &&
    is.character(key$term) && length(key$term) == 1L && !is.na(key$term)
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
  record <- found$record
  eta <- model_samples[, record$gamma_coordinates, drop = FALSE]
  open <- .bt_R2_record_open(record, model_samples)
  if(!is.null(record$core)){
    unit_open <- open[, vapply(record$core$unit_term, function(key){
      names(record$terms)[vapply(record$terms, function(term) identical(term$key, key), logical(1))]
    }, character(1)), drop = FALSE]
    eta <- eta * unit_open
  }
  total <- rowSums(eta)
  share <- rowSums(eta[, key$units, drop = FALSE]) / ifelse(total > 0, total, 1)
  # a share is undefined where its term is closed
  share[open[, key$term] == 0] <- NA_real_

  share
}

# The draws of the inclusion of a gated term.
.bt_parameter_draw_R2_open <- function(fit, key, model_samples){

  if(!all(key$dependencies %in% colnames(model_samples))){
    .bt_R2_unavailable(paste0("The gate of the R-squared budget term '", key$term,
      "' is unavailable in the draws. Include its declared source coordinates."))
  }

  apply(model_samples[, key$dependencies, drop = FALSE], 1, prod)
}

# The gate and indicator coordinates of a gated record (empty otherwise).
.bt_R2_record_gating <- function(record){

  core <- record$core
  if(is.null(core)){
    return(character())
  }
  nodes <- .bt_R2_core_nodes(core)
  c(nodes$gates, if(!is.null(core$gates$spike)) nodes$indicator)
}

# The prior density of a share: the Beta marginal of its units under the
# Dirichlet allocation of the model.
.bt_parameter_prior_density_R2_share <- function(object, quantity, key, n_grid, tail_prob){

  prior_list <- attr(object, "prior_list", exact = TRUE)
  found <- .bt_R2_share_record(prior_list, key)
  source <- .bt_R2_share_prior_law(found)
  out <- .bt_parameter_prior_density_scalar(source, n_grid = n_grid, tail_prob = tail_prob)
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

  records <- .bt_R2_budget_records(prior_list)
  # the totals and the shared inclusion rates of gated budgets
  totals <- c(names(records), unlist(lapply(records, function(record){
    if(!is.null(record$core)) .bt_R2_core_nodes(record$core)$rate
  }), use.names = FALSE))
  # the scalar coefficients of gated terms (their atom at zero needs the
  # catalog law) and the levels of R2 factor terms (the law of each level)
  gated <- unlist(lapply(names(prior_list), function(parameter){
    prior <- prior_list[[parameter]]
    if(is.prior.R2(prior) && .bt_R2_term_is_gated(prior) && prior$d == 1L && !is.prior.factor(prior)) parameter
  }), use.names = FALSE)
  factors <- names(prior_list)[vapply(prior_list, function(prior){
    is.prior.R2(prior) && is.prior.factor(prior)
  }, logical(1))]
  levels <- unlist(lapply(factors, function(parameter){
    .bt_R2_coefficient_names(prior_list[[parameter]], parameter)
  }), use.names = FALSE)
  vapply(seq_len(nrow(quantities)), function(i){
    key <- quantities$extraction_key[[i]]
    identical(quantities$provider[[i]], "BayesTools") &&
      (identical(key$type, "R2_share") ||
         (identical(key$type, "coordinate") && length(key$dependencies) == 1L &&
            key$dependencies %in% c(totals, gated)) ||
         (identical(key$type, "factor_level") && length(key$dependencies) > 0L &&
            all(key$dependencies %in% levels)))
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

# The open flag of every retained term of 'record' per draw (1 for a term of
# a budget without gates and spike): its gate times the budget indicator of a
# spike-and-slab total. NULL when a gate or indicator column is missing.
.bt_R2_record_open <- function(record, samples){

  n <- nrow(samples)
  core <- record$core
  out <- matrix(1, n, length(record$terms), dimnames = list(NULL, names(record$terms)))
  if(is.null(core)){
    return(out)
  }
  nodes <- .bt_R2_core_nodes(core)
  needed <- c(nodes$gates, if(!is.null(core$gates$spike)) nodes$indicator)
  if(!all(needed %in% colnames(samples))){
    return(NULL)
  }
  indicator <- if(!is.null(core$gates$spike)) samples[, nodes$indicator] else rep(1, n)
  for(term in names(record$terms)){
    gate <- core$term_gate[[record$terms[[term]]$key]]
    out[, term] <- indicator * (if(gate > 0L) samples[, nodes$gates[[gate]]] else 1)
  }

  out
}

# Whether the budget of 'record' can carry no signal (its null is a model):
# a spike-and-slab total or gates on every term.
.bt_R2_record_has_null <- function(record){

  core <- record$core
  !is.null(core) && (!is.null(core$gates$spike) || (core$gates$K > 0L && core$gates$all_gated))
}

# Estimates tables of R-squared budgets. The budget total is followed by the
# share of each retained term in place of the allocation coordinates (the
# budgets come from the complete fitted prior list 'prior_list' and the
# values from the unfiltered draws, so that the display filters act on the
# quantities: the shares follow their total, and backend coordinates are never
# shown). A share is undefined (NA) on the draws where its term is closed.
# Gated and spiked budgets show inclusion rows in the style of spike-and-slab
# priors: '<term> (inclusion)' before every gated term and
# '<name>{total} (inclusion)' before a total that can be zero; with
# 'conditional', the coefficients of a gated term (and a total) are summarized
# on the draws where the term (the budget) is open. Returns the samples, the
# label parts of the created columns, the inclusion columns and the warnings.
.bt_JAGS_estimates_R2_budgets <- function(model_samples, raw_model_samples, prior_list,
                                          conditional = FALSE, remove_inclusion = FALSE){

  created <- list()
  inclusion_columns <- character()
  warnings <- character()
  insert <- function(samples, columns, position){
    cbind(
      if(position > 1L) samples[, seq_len(position - 1L), drop = FALSE],
      columns,
      if(position <= ncol(samples)) samples[, position:ncol(samples), drop = FALSE]
    )
  }
  for(record in .bt_R2_budget_records(prior_list)){
    core <- record$core
    nodes <- if(!is.null(core)) .bt_R2_core_nodes(core)
    allocation_columns <- if(!is.null(record$allocation_name)){
      paste0(record$allocation_name, "[", seq_len(record$n_units), "]")
    }
    backend <- c(allocation_columns, nodes$gates, nodes$indicator, nodes$inclusion)
    model_samples <- model_samples[, !colnames(model_samples) %in% backend, drop = FALSE]
    open <- .bt_R2_record_open(record, raw_model_samples)
    if(is.null(open)){
      next
    }

    # conditional summaries: a gated term on the draws where it is open, a
    # total that can be zero on the draws where the budget is on
    if(conditional && !is.null(core)){
      for(term in names(record$terms)){
        if(core$term_gate[[record$terms[[term]]$key]] == 0L && is.null(core$gates$spike)){
          next
        }
        columns <- intersect(.bt_R2_coefficient_names(prior_list[[term]], term), colnames(model_samples))
        if(length(columns)){
          model_samples[open[, term] == 0, columns] <- NA_real_
          warnings <- c(warnings, .runjags_conditional_warning(columns, sum(open[, term] == 1)))
        }
      }
      if(.bt_R2_record_has_null(record) && record$total_name %in% colnames(model_samples)){
        on <- rowSums(open) > 0
        model_samples[!on, record$total_name] <- NA_real_
        warnings <- c(warnings, .runjags_conditional_warning(record$total_name, sum(on)))
      }
    }

    # the shares after the total
    if(length(record$terms) >= 2L && record$total_name %in% colnames(model_samples) &&
       all(allocation_columns %in% colnames(raw_model_samples))){
      phi <- raw_model_samples[, allocation_columns, drop = FALSE]
      shares <- vapply(names(record$terms), function(term){
        value <- rowSums(phi[, record$terms[[term]]$units, drop = FALSE])
        value[open[, term] == 0] <- NA_real_
        value
      }, numeric(nrow(phi)))
      shares <- matrix(shares, nrow = nrow(phi))
      colnames(shares) <- vapply(record$terms, function(term) .bt_R2_share_name(record, term), character(1))
      model_samples <- insert(model_samples, shares, match(record$total_name, colnames(model_samples)) + 1L)
      for(term in record$terms){
        created[[.bt_R2_share_name(record, term)]] <- .bt_R2_share_label_parts(record, term)
      }
    }

    # inclusion rows
    if(!remove_inclusion && !is.null(core)){
      if(.bt_R2_record_has_null(record) && record$total_name %in% colnames(model_samples)){
        column <- paste0(record$total_name, " (inclusion)")
        model_samples <- insert(model_samples, matrix(as.numeric(rowSums(open) > 0), ncol = 1L,
          dimnames = list(NULL, column)), match(record$total_name, colnames(model_samples)))
        created[[column]] <- .bt_label_parts_update(.bt_R2_total_label_parts_list(prior_list)[[record$total_name]],
          inclusion = "")[[1L]]
        inclusion_columns <- c(inclusion_columns, column)
      }
      for(term in names(record$terms)){
        if(core$term_gate[[record$terms[[term]]$key]] == 0L){
          next
        }
        columns <- .bt_R2_coefficient_names(prior_list[[term]], term)
        present <- which(colnames(model_samples) %in% columns)
        if(length(present) == 0L){
          next
        }
        column <- paste0(term, " (inclusion)")
        model_samples <- insert(model_samples, matrix(open[, term], ncol = 1L, dimnames = list(NULL, column)), min(present))
        created[[column]] <- .bt_label_parts_update(.bt_label_parts_inclusion_term(term, prior_list[[term]]),
          inclusion = "")[[1L]]
        inclusion_columns <- c(inclusion_columns, column)
      }
    }
  }

  list(model_samples = model_samples, created = created, inclusion_columns = inclusion_columns,
       warnings = warnings)
}

# The rows of runjags_inference_table() for the budgets of 'prior_list': the
# budget (when it can be zero) and every gated term, with the marginal prior
# inclusion probability, the posterior inclusion probability and the
# inclusion Bayes factor. Each row is list(parameter, label parts, prior,
# indicators per chain).
.bt_R2_inference_rows <- function(prior_list, model_samples_list){

  rows <- list()
  for(record in .bt_R2_budget_records(prior_list)){
    core <- record$core
    if(is.null(core)){
      next
    }
    marginals <- .bt_R2_gate_marginals(core$gates, core$terms)
    open_list <- lapply(model_samples_list, function(samples) .bt_R2_record_open(record, as.matrix(samples)))
    if(any(vapply(open_list, is.null, logical(1)))){
      next
    }
    if(.bt_R2_record_has_null(record)){
      rows[[length(rows) + 1L]] <- list(
        parameter  = record$total_name,
        parts      = .bt_R2_total_label_parts_list(prior_list)[[record$total_name]],
        prior_prob = marginals$budget,
        indicators = lapply(open_list, function(open) as.numeric(rowSums(open) > 0))
      )
    }
    for(term in names(record$terms)){
      key <- record$terms[[term]]$key
      if(core$term_gate[[key]] == 0L){
        next
      }
      rows[[length(rows) + 1L]] <- list(
        parameter  = term,
        parts      = .bt_label_parts_inclusion_term(term, prior_list[[term]]),
        prior_prob = marginals$open[[key]],
        indicators = lapply(open_list, function(open) open[, term])
      )
    }
  }

  rows
}

# The prior law of a catalog quantity of a gated or spiked budget (NULL for
# other quantities): the realized total (atom at zero), the budget indicator
# and the inclusion of a gated term (Bernoulli of the marginal probability).
.bt_R2_quantity_law <- function(object, key, conditional = FALSE){

  prior_list <- attr(object, "prior_list", exact = TRUE)
  if(identical(key$type, "R2_open")){
    record <- .bt_R2_budget_records(prior_list)[[key$total]]
    core <- record$core
    marginals <- .bt_R2_gate_marginals(core$gates, core$terms)
    return(.bt_R2_bernoulli_law(marginals$open[[record$terms[[key$term]]$key]]))
  }
  if(identical(key$type, "factor_level")){
    return(.bt_R2_level_law(prior_list, key, conditional))
  }
  if(!identical(key$type, "coordinate")){
    return(NULL)
  }
  if(length(key$dependencies) != 1L){
    return(.bt_R2_level_law(prior_list, key, conditional))
  }
  core <- .bt_R2_gated_total_core(prior_list, key$dependencies)
  if(!is.null(core)){
    return(.bt_R2_core_total_law(core, conditional = conditional))
  }
  core <- .bt_R2_indicator_core(prior_list, key$dependencies)
  if(!is.null(core)){
    return(.bt_R2_bernoulli_law(.bt_R2_gate_marginals(core$gates, core$terms)$budget))
  }
  core <- .bt_R2_rate_core(prior_list, key$dependencies)
  if(!is.null(core)){
    return(.bt_R2_core_rate_law(core))
  }

  .bt_R2_gated_coefficient_law(prior_list, key, conditional)
}

# The prior law of a factor level (or a coded coefficient of a factor) of an
# R2 term, gated or not (NULL for other quantities).
.bt_R2_level_law <- function(prior_list, key, conditional = FALSE){

  owner <- .bt_R2_term_owner(prior_list, key$dependencies)
  if(is.null(owner)){
    return(NULL)
  }
  if(.bt_R2_term_is_gated(prior_list[[owner]])){
    return(.bt_R2_gated_coefficient_law(prior_list, key, conditional))
  }

  .bt_R2_functional_law(prior_list[[owner]], owner, key)
}

# The point plan of the realized total of a gated or spiked budget: an atom at
# zero (the budget null) selected by the budget indicator, whose inclusion
# event is the budget being on. NULL for other quantities and budgets
# without a null.
.bt_R2_total_point_plan <- function(fit, quantity){

  key <- quantity$extraction_key[[1L]]
  if(!identical(key$type, "coordinate") || length(key$dependencies) != 1L){
    return(NULL)
  }
  prior_list <- attr(fit, "prior_list", exact = TRUE)
  core <- .bt_R2_gated_total_core(prior_list, key$dependencies)
  indicator <- if(!is.null(core)) .bt_R2_core_nodes(core)$indicator
  if(is.null(indicator)){
    return(NULL)
  }
  on <- .bt_R2_gate_marginals(core$gates, core$terms)$budget
  # the spike-and-slab form of the indicator: 0 the null, 1 the slab
  plan_prior <- prior_spike_and_slab(prior("beta", list(1, 1)), prior_inclusion = prior("spike", list(on)))

  list(
    kind        = "point",
    origin      = "prior",
    owner       = key$dependencies,
    prior       = plan_prior,
    locations   = vapply(attr(plan_prior, "components"), function(component) if(component == "null") 0 else NA_real_, numeric(1)),
    prior_probabilities = .bt_prior_component_probabilities(plan_prior),
    indicator   = indicator,
    unknown     = FALSE,
    chain_gates = character(),
    event_components = attr(plan_prior, "components") == "alternative",
    event_reason     = NULL,
    event_gates = indicator,
    event_rule  = "AND"
  )
}
