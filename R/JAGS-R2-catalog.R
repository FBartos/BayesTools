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

# Fully indexed nodes monitored by a bound budget, keyed by their fitted names.
# A declared source absent from the samples can still be registered as unavailable.
.bt_R2_declared_nodes <- function(prior_list){

  rows <- list()
  add <- function(coordinates, role, term, parameter){
    if(length(coordinates)){
      rows[[length(rows) + 1L]] <<- data.frame(
        coordinate = unname(coordinates), role = role, term = term,
        formula_parameter = parameter, stringsAsFactors = FALSE
      )
    }
  }
  for(record in .bt_R2_budget_records(prior_list)){
    nodes <- if(is.null(record$core)){
      list(total = record$total_name,
           allocation = if(!is.null(record$allocation_name)){
             paste0(record$allocation_name, "[", seq_len(record$n_units), "]")
           },
           eta = record$gamma_coordinates)
    }else{
      declared <- .bt_R2_core_nodes(record$core)
      monitors <- .bt_R2_core_monitor(record$core)
      lapply(declared, function(values){
        values[.bt_parameter_coordinates_base(values) %in% monitors]
      })
    }
    for(role in names(nodes)){
      add(nodes[[role]], role, record$name, record$parameter)
    }
    for(term in record$terms){
      prior <- prior_list[[term$parameter]]
      add(.bt_R2_coefficient_names(prior, term$parameter), "coefficient", term$label, record$parameter)
      add(.bt_R2_direction_coordinates(prior, term$parameter), "direction", term$label, record$parameter)
    }
  }
  if(!length(rows)){
    return(data.frame(coordinate = character(), role = character(), term = character(),
                      formula_parameter = character(), stringsAsFactors = FALSE))
  }
  do.call(rbind, rows)
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

# Normalize the shared Gamma draws over the open allocation units. Tables and
# deferred extraction use the same measure, including NA for a closed term.
.bt_R2_share_values <- function(record, samples, terms = names(record$terms),
                                open = .bt_R2_record_open(record, samples)){

  .bt_R2_require_columns(samples, record$gamma_coordinates)
  eta <- samples[, record$gamma_coordinates, drop = FALSE]
  if(any(!is.finite(eta)) || any(eta < 0)){
    .bt_R2_unavailable(paste0("The Gamma weights of R-squared budget '",
      record$name, "' must be finite and nonnegative on every supplied draw."))
  }
  if(!is.null(record$core)){
    term_keys <- vapply(record$terms, `[[`, character(1), "key")
    eta <- eta * open[, match(record$core$unit_term, term_keys), drop = FALSE]
  }
  total <- rowSums(eta)
  retained <- rowSums(open) > 0
  if(any(!is.finite(total[retained]) | total[retained] <= 0)){
    .bt_R2_unavailable(paste0("The open allocation units of R-squared budget '",
      record$name, "' require a positive finite sum of Gamma weights."))
  }
  shares <- matrix(NA_real_, nrow(samples), length(terms))
  for(i in seq_along(terms)){
    term <- record$terms[[terms[[i]]]]
    included <- open[, term$parameter] == 1
    numerator <- rowSums(eta[included, term$units, drop = FALSE])
    shares[included, i] <- numerator / total[included]
  }
  colnames(shares) <- vapply(record$terms[terms], function(term){
    .bt_R2_share_name(record, term)
  }, character(1))
  shares
}

# The catalog quantities of the shares of the budgets in 'prior_list'.
.bt_R2_share_quantities <- function(coordinates, prior_list){

  rows <- list()
  records <- .bt_R2_budget_records(prior_list)
  for(record in records){
    if(length(record$terms) < 2L){
      next
    }
    dependencies <- c(record$gamma_coordinates, .bt_R2_record_gating(record))
    sources <- match(dependencies, coordinates$coordinate_name)
    status <- if(anyNA(sources) || any(coordinates$monitor_status[sources] == "unavailable")){
      "unavailable"
    }else{
      "derived"
    }
    for(term in record$terms){
      rows[[length(rows) + 1L]] <- .bt_parameter_catalog_quantity(
        canonical_name    = .bt_R2_share_name(record, term),
        namespace         = record$parameter,
        role              = "parameter",
        formula_parameter = record$parameter,
        owner_type        = "r2_budget",
        owner_name        = record$name,
        quantity          = "share",
        term              = term$label,
        component         = term$label,
        arguments         = term$label,
        parent_quantity_id = .bt_parameter_catalog_quantity_id(record$total_name, record$parameter, "parameter"),
        label_parts       = .bt_R2_share_label_parts(record, term),
        fitted_scale      = "unitless",
        display_scale     = "unitless",
        status            = status,
        source_type       = "composite",
        extraction_key    = list(
          type         = "R2_share",
          dependencies = dependencies,
          total        = record$total_name,
          term         = term$parameter,
          units        = as.integer(term$units)
        )
      )
    }
  }
  # the inclusion of every gated term: its gate (times the budget indicator
  # of a spike-and-slab total)
  for(record in records){
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
      open_name <- paste0(record$total_name, "_open[", term$label, "]")
      rows[[length(rows) + 1L]] <- .bt_parameter_catalog_quantity(
        canonical_name    = open_name,
        namespace         = record$parameter,
        role              = "parameter",
        formula_parameter = record$parameter,
        owner_type        = "r2_budget",
        owner_name        = record$name,
        quantity          = "term_inclusion",
        term              = term$label,
        component         = term$label,
        arguments         = term$label,
        parent_quantity_id = .bt_parameter_catalog_quantity_id(record$total_name, record$parameter, "parameter"),
        label_parts       = .bt_label_parts(
          components        = paste0(term$label, " (inclusion)"),
          formula_parameter = record$parameter,
          selector          = open_name
        ),
        fitted_scale      = "unitless",
        display_scale     = "unitless",
        status            = if(anyNA(position) || any(coordinates$monitor_status[position] == "unavailable")) "unavailable" else "derived",
        source_type       = "identity",
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

# Add budget membership to the existing coordinate and factor-cell rows. The
# coefficient's public identity and factor coding belong to the formula.
.bt_R2_register_catalog_owners <- function(quantities, prior_list){

  for(record in .bt_R2_budget_records(prior_list)){
    formula_rows <- quantities$formula_parameter == record$parameter
    labels <- vapply(record$terms, `[[`, character(1), "label")
    coefficient_rows <- which(formula_rows & quantities$role == "fixed_coefficient" &
                              quantities$term %in% labels)
    nodes <- if(!is.null(record$core)) .bt_R2_core_nodes(record$core)
    public_nodes <- c(total = record$total_name, budget_inclusion = nodes$indicator,
                      inclusion_rate = nodes$rate)
    rows <- match(unname(public_nodes), quantities$canonical_name)
    if(anyNA(rows) || !all(formula_rows[rows])){
      .bt_stop_refit_required("The fitted parameter map is missing declared R-squared budget quantities. Refit the model.")
    }
    quantities$quantity[rows] <- names(public_nodes)
    quantities$quantity[coefficient_rows] <- "coefficient"
    members <- c(rows, coefficient_rows)
    quantities$owner_type[members] <- "r2_budget"
    quantities$owner_name[members] <- record$name
    total_row <- rows[[1L]]
    quantities$parent_quantity_id[setdiff(members, total_row)] <- quantities$quantity_id[total_row]
    quantities$arguments[coefficient_rows] <- as.list(quantities$term[coefficient_rows])
  }
  quantities
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
  as.numeric(.bt_R2_share_values(found$record, model_samples, terms = key$term))
}

# The draws of the inclusion of a gated term.
.bt_parameter_draw_R2_open <- function(fit, key, model_samples){

  open <- rep(1, nrow(model_samples))
  for(column in key$dependencies){
    open <- open * .bt_R2_event_draws(model_samples, column)
  }
  open
}

# A realized open flag includes the budget spike. An ungated member of a
# spiked budget reads the budget indicator; an always-open member has no input.
.bt_R2_record_open_columns <- function(record){

  core <- record$core
  nodes <- if(!is.null(core)) .bt_R2_core_nodes(core)
  lapply(record$terms, function(term) .bt_R2_core_open_column(core, term$key, nodes))
}

.bt_R2_record_gating <- function(record){

  unique(as.character(unlist(.bt_R2_record_open_columns(record), use.names = FALSE)))
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

# Registered quantities whose hypotheses use the budget law: total, rate,
# share, gated scalar coefficient, or factor functional with fitted sources.
.bt_R2_quantity_rows <- function(quantities, prior_list){

  owned <- quantities$provider == "BayesTools" & quantities$owner_type == "r2_budget"
  selected <- owned & quantities$quantity %in% c("total", "share", "inclusion_rate")
  for(i in which(owned & quantities$quantity == "coefficient")){
    key <- quantities$extraction_key[[i]]
    if(identical(key$type, "factor_level")){
      selected[i] <- length(key$dependencies) > 0L
    }else{
      prior <- prior_list[[quantities$canonical_name[i]]]
      selected[i] <- is.prior.R2(prior) && .bt_R2_term_is_gated(prior)
    }
  }
  selected
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

.bt_R2_require_columns <- function(samples, columns){

  missing <- setdiff(columns, colnames(samples))
  if(length(missing)){
    .bt_R2_unavailable(paste0(
      "The supplied draws lack R-squared budget source coordinates ",
      paste0("'", missing, "'", collapse = ", "), ". Include the declared coordinates."
    ), class = "BayesTools_R2_coordinates_unavailable")
  }
}

.bt_R2_event_draws <- function(samples, column){

  .bt_R2_require_columns(samples, column)
  values <- as.numeric(samples[, column])
  if(any(!values %in% c(0, 1))){
    .bt_R2_unavailable(paste0("The R-squared budget indicator '", column,
      "' must be zero or one on every supplied draw."))
  }
  values
}

# A retained term's effective state is read from its registered open node.
# Request only the terms whose states enter the calculation.
.bt_R2_record_open <- function(record, samples, terms = names(record$terms),
                               sources = .bt_R2_record_open_columns(record)){

  open <- matrix(1, nrow(samples), length(terms), dimnames = list(NULL, terms))
  for(term in terms){
    source <- sources[[term]]
    if(!is.null(source)){
      open[, term] <- .bt_R2_event_draws(samples, source)
    }
  }
  open
}

# Whether the budget of 'record' can carry no signal (its null is a model):
# a spike-and-slab total or gates on every term.
.bt_R2_record_has_null <- function(record){

  core <- record$core
  !is.null(core) && (!is.null(core$gates$spike) || (core$gates$K > 0L && core$gates$all_gated))
}

# The fitted catalog identifies the exact event displayed for a budget or one
# of its terms, independently of its rendering and backend coordinate name.
.bt_R2_event_quantity <- function(quantities, record, term = NULL){

  owner <- quantities$provider == "BayesTools" & quantities$owner_type == "r2_budget" &
    quantities$formula_parameter == record$parameter & quantities$owner_name == record$name
  event <- if(is.null(term)){
    quantities$quantity == "budget_inclusion"
  }else{
    quantities$quantity == "term_inclusion" & quantities$term == term$label
  }
  rows <- which(owner & event)
  if(length(rows) != 1L){
    .bt_stop_refit_required("The R-squared budget lacks a unique inclusion quantity in the fitted catalog. Refit the model.")
  }
  quantities[rows, , drop = FALSE]
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
.bt_JAGS_estimates_R2_budgets <- function(fit, model_samples, raw_model_samples, prior_list,
                                          conditional = FALSE, remove_inclusion = FALSE,
                                          raw = FALSE){

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
    if(!raw){
      backend <- c(allocation_columns, nodes$gates, nodes$open, nodes$indicator, nodes$inclusion)
      model_samples <- model_samples[, !colnames(model_samples) %in% backend, drop = FALSE]
    }
    total_present <- record$total_name %in% colnames(model_samples)
    show_shares <- !raw && total_present && length(record$terms) >= 2L
    term_columns <- lapply(names(record$terms), function(term){
      intersect(.bt_R2_coefficient_names(prior_list[[term]], term), colnames(model_samples))
    })
    names(term_columns) <- names(record$terms)
    needed_terms <- if(show_shares){
      names(record$terms)
    }else if(conditional || !remove_inclusion){
      names(term_columns)[lengths(term_columns) > 0L]
    }else{
      character()
    }
    open_columns <- .bt_R2_record_open_columns(record)
    need_total_event <- total_present && .bt_R2_record_has_null(record) &&
      (conditional || !remove_inclusion)
    sources <- unique(c(
      if(show_shares) record$gamma_coordinates,
      if(need_total_event) nodes$indicator,
      unlist(open_columns[needed_terms], use.names = FALSE)
    ))
    # Materialization honors the current draw view and gives its typed
    # missing-source diagnostic. It also supplies any structural coordinates.
    missing <- setdiff(sources, colnames(raw_model_samples))
    if(length(missing)){
      coordinates <- parameter_coordinates(fit)
      unavailable <- setdiff(missing, coordinates$coordinate_name[
        coordinates$monitor_status %in% c("sampled", "structural")
      ])
      if(length(unavailable)){
        .bt_stop_refit_required("The R-squared budget requires unmonitored coordinates ",
          paste0("'", unavailable, "'", collapse = ", "), ". Refit with these coordinates monitored.",
          class = "BayesTools_refit_monitoring")
      }
      raw_model_samples <- cbind(raw_model_samples, as.matrix(.bt_parameter_draw_dependencies(fit, missing)))
    }
    open <- .bt_R2_record_open(record, raw_model_samples, needed_terms, open_columns)
    if(need_total_event){
      budget_open <- .bt_R2_event_draws(raw_model_samples, nodes$indicator)
    }

    # Conditional summaries: a gated term on the draws where it is open, a
    # total that can be zero on the draws where the budget is on
    if(conditional && !is.null(core)){
      for(term in names(record$terms)){
        if(core$term_gate[[record$terms[[term]]$key]] == 0L && is.null(core$gates$spike)){
          next
        }
        columns <- term_columns[[term]]
        if(length(columns)){
          model_samples[open[, term] == 0, columns] <- NA_real_
          warnings <- c(warnings, .runjags_conditional_warning(columns, sum(open[, term] == 1)))
        }
      }
      if(need_total_event){
        model_samples[budget_open == 0, record$total_name] <- NA_real_
        warnings <- c(warnings, .runjags_conditional_warning(record$total_name, sum(budget_open)))
      }
    }

    # the shares after the total
    if(show_shares){
      shares <- .bt_R2_share_values(record, raw_model_samples, open = open)
      model_samples <- insert(model_samples, shares, match(record$total_name, colnames(model_samples)) + 1L)
      for(term in record$terms){
        created[[.bt_R2_share_name(record, term)]] <- .bt_R2_share_label_parts(record, term)
      }
    }

    # The canonical event IDs and labels are shared with direct selections.
    if(!raw && !remove_inclusion && !is.null(core)){
      quantities <- parameter_catalog(fit)$quantities
      if(need_total_event){
        event <- .bt_R2_event_quantity(quantities, record)
        column <- event$canonical_name
        model_samples <- insert(model_samples, matrix(budget_open, ncol = 1L,
          dimnames = list(NULL, column)), match(record$total_name, colnames(model_samples)))
        created[[column]] <- event$label_parts[[1L]]
        inclusion_columns <- c(inclusion_columns, column)
      }
      for(term in names(record$terms)){
        if(core$term_gate[[record$terms[[term]]$key]] == 0L){
          next
        }
        present <- which(colnames(model_samples) %in% term_columns[[term]])
        if(length(present) == 0L){
          next
        }
        event <- .bt_R2_event_quantity(quantities, record, record$terms[[term]])
        column <- event$canonical_name
        model_samples <- insert(model_samples, matrix(open[, term], ncol = 1L, dimnames = list(NULL, column)), min(present))
        created[[column]] <- event$label_parts[[1L]]
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
.bt_R2_inference_rows <- function(prior_list, model_samples_list, quantities){

  rows <- list()
  for(record in .bt_R2_budget_records(prior_list)){
    core <- record$core
    if(is.null(core)){
      next
    }
    marginals <- .bt_R2_gate_marginals(core$gates, core$terms)
    open_list <- lapply(model_samples_list, function(samples) .bt_R2_record_open(record, as.matrix(samples)))
    if(.bt_R2_record_has_null(record)){
      event <- .bt_R2_event_quantity(quantities, record)
      indicator <- event$extraction_key[[1L]]$dependencies
      rows[[length(rows) + 1L]] <- list(
        parameter  = record$total_name,
        event_id   = event$quantity_id,
        role       = "r2_budget_inclusion",
        parts      = .bt_R2_total_label_parts_list(prior_list)[[record$total_name]],
        prior_prob = marginals$budget,
        indicators = lapply(model_samples_list, .bt_R2_event_draws, column = indicator)
      )
    }
    for(term in names(record$terms)){
      key <- record$terms[[term]]$key
      if(core$term_gate[[key]] == 0L){
        next
      }
      event <- .bt_R2_event_quantity(quantities, record, record$terms[[term]])
      rows[[length(rows) + 1L]] <- list(
        parameter  = term,
        event_id   = event$quantity_id,
        role       = "r2_term_inclusion",
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

# The factor by which the scale transformation of the predictors multiplies
# the coefficients 'dependencies' of an R2 term of the formula parameter
# 'formula_parameter': 1 when it leaves them unchanged, the weight of a pure
# rescaling of a single coefficient, and NULL when an original-scale
# coefficient combines several fitted coefficients.
.bt_R2_original_scale_factor <- function(fit, formula_parameter, dependencies){

  if(is.null(attr(fit, "formula_scale", exact = TRUE)[[formula_parameter]])){
    return(1)
  }
  routes <- lapply(dependencies, function(dependency){
    JAGS_formula_coefficient_route(fit, formula_parameter, dependency)
  })
  types <- vapply(routes, function(route) route[["type"]], character(1))
  if(all(types == "identity")){
    return(1)
  }
  # every coefficient multiplied by its own constant: a common constant
  # multiplies every combination of them
  factors <- vapply(seq_along(dependencies), function(i){
    if(identical(types[[i]], "identity")) return(1)
    weights <- routes[[i]][["weights"]]
    if(identical(types[[i]], "affine") && length(weights) == 1L && identical(names(weights), dependencies[[i]]) &&
       is.finite(weights) && weights != 0){
      return(unname(weights))
    }
    NA_real_
  }, numeric(1))
  if(anyNA(factors) || max(abs(factors - factors[[1L]])) > 1e-12 * abs(factors[[1L]])){
    return(NULL)
  }
  factors[[1L]]
}

# The mixed posterior of a coefficient of an R2 term on the original scale of
# its predictor (see .bt_R2_original_scale_factor()): the draws multiplied by
# the factor and the prior density of the rescaled catalog law (the atom of a
# gate stays at zero). Other quantities are returned unchanged.
.bt_R2_original_scale_mixed_posterior <- function(fit, selection, values, conditional){

  quantity <- selection$quantities[1L, , drop = FALSE]
  key <- quantity$extraction_key[[1L]]
  prior_list <- attr(fit, "prior_list", exact = TRUE)
  if(!key$type %in% c("coordinate", "factor_level")){
    return(values)
  }
  owner <- .bt_R2_term_owner(prior_list, key$dependencies)
  if(is.null(owner)){
    return(values)
  }
  factor <- .bt_R2_original_scale_factor(fit, prior_list[[owner]]$parameter, key$dependencies)
  if(is.null(factor)){
    .bt_R2_unavailable(paste0("The original-scale coefficient '", quantity$canonical_name,
      "' combines several fitted coefficients of an R-squared budget; its prior law is unavailable."),
      class = "BayesTools_R2_density_unavailable")
  }
  if(factor == 1){
    return(values)
  }
  scaled <- selection
  weights <- if(is.null(key$weights)) rep(1, length(key$dependencies)) else key$weights
  scaled$quantities$extraction_key[[1L]]$weights <- weights * factor
  density <- .bt_parameter_prior_density_quantity_law(fit, scaled,
    n_grid = .prior_linear_density_default_grid(), tail_prob = .prior_linear_density_tail_prob(),
    conditional = conditional)
  # the draws are multiplied by the factor: the metadata that describe values
  # (support, atoms, prior density) are rescaled with them, the others kept
  meta <- .bt_meta_container(values)
  values[] <- as.numeric(values) * factor
  values <- .bt_meta_refresh(values)
  support <- meta$support
  if(!is.null(support)){
    support$bounds <- sort(support$bounds * factor)
    if(length(support$points)) support$points <- sort(support$points * factor)
  }
  atoms <- meta$atoms
  if(!is.null(atoms) && length(atoms$locations)){
    atoms$locations[] <- atoms$locations * factor
  }

  .bt_meta_assign(values, list(support = support, atoms = atoms, prior_density = density))
}

# The original-scale factor of the coefficient 'parameter' of the R2 term
# 'prior' from the declared predictor scales 'formula_scale' (1 when the
# scale transformation leaves it unchanged, the weight of a pure rescaling,
# NULL for a combination of several fitted coefficients).
.bt_R2_formula_scale_factor <- function(formula_scale, prior, parameter){

  prefix <- prior$parameter
  scale <- formula_scale[[prefix]]
  spec <- if(!is.null(scale)) attr(scale, "unscale_design", exact = TRUE)
  if(is.null(spec)){
    return(1)
  }
  transform <- .bt_formula_coefficient_transform(names(spec$multipliers), scale, prefix,
    log_intercept = isTRUE(attr(scale, "log_intercept", exact = TRUE)))
  recipe <- transform$prior_recipes[[parameter]]
  if(is.null(recipe) || !identical(recipe$type, "raw_affine")){
    return(NULL)
  }
  weights <- recipe$weights[recipe$weights != 0]
  if(length(weights) == 1L && identical(names(weights), parameter) && is.finite(weights)){
    return(unname(weights))
  }
  NULL
}

# The law of a coefficient multiplied by 'factor': its R2 components restricted
# to the multiplied combination, its atom at zero unchanged.
.bt_R2_scale_law <- function(law, factor){

  if(factor == 1){
    return(law)
  }
  scale_part <- function(part){
    if(is.prior.R2(part)){
      part$functional_L <- factor * (if(is.null(part$functional_L)) 1 else part$functional_L)
      class(part) <- c("prior", "prior.simple", "prior.R2")
    }
    part
  }
  if(!is.prior.mixture(law)){
    return(scale_part(law))
  }
  parts <- lapply(unclass(law), scale_part)
  weights <- vapply(parts, function(part) part$prior_weights, numeric(1))
  .bt_R2_mixture(parts, weights, vapply(parts, is.prior.point, logical(1)))
}

# The transformed prior densities of as_mixed_posteriors() with the scalar
# coefficients of R2 terms replaced by their own laws on the original scale
# (given the condition event), when the scale transformation rescales them;
# the other entries (and a coefficient combined with others) are kept.
.bt_R2_scaled_prior_densities <- function(prior_densities, prior_list, formula_scale, condition_event){

  for(parameter in names(prior_densities)){
    prior <- prior_list[[parameter]]
    if(!is.prior.R2(prior) || !isTRUE(prior$bound) || prior$d != 1L || is.prior.factor(prior)){
      next
    }
    factor <- .bt_R2_formula_scale_factor(formula_scale, prior, parameter)
    if(is.null(factor)){
      next
    }
    density <- tryCatch({
      law <- .bt_R2_conditioned_law(prior_list, condition_event, parameter)
      if(is.null(law)){
        key <- list(type = "coordinate", dependencies = parameter)
        law <- if(.bt_R2_term_is_gated(prior)){
          .bt_R2_gated_coefficient_law(prior_list, key)
        }else{
          .bt_R2_functional_law(prior, parameter, key)
        }
      }
      law <- .bt_R2_scale_law(law, factor)
      .bt_parameter_prior_density_scalar(.bt_R2_grid_law(law), n_grid = .prior_linear_density_default_grid(),
        tail_prob = .prior_linear_density_tail_prob())
    }, BayesTools_R2_density_unavailable = function(condition) NULL)
    if(!is.null(density)){
      prior_densities[[parameter]] <- density
    }
  }

  prior_densities
}
