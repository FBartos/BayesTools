#' Mixed posterior of one catalog quantity
#'
#' @description `parameter_mixed_posterior()` returns the posterior draws of one
#' semantic quantity of the [parameter_catalog()] as a mixed posterior (the
#' element form returned by [as_mixed_posteriors()]) together with the draw
#' metadata that inference and plots need ([posterior_metadata()]):
#' \describe{
#'   \item{`support`}{the exact support declared by the catalog.}
#'   \item{`prior_density`}{the canonical prior density,
#'   [parameter_prior_density()] (conditional on the inclusion event with
#'   `conditional = TRUE`).}
#'   \item{`atoms`}{the declared posterior point masses, from the structure
#'   of the quantity, never inferred from the draw values: the atoms on which
#'   the per-draw states of its inclusion gates and point components put it
#'   ([parameter_gate_states()]: the gates of a variance allocation, the
#'   mixture indicator of a scale prior with a point component, or the
#'   component indicator of a mixture or spike-and-slab prior of its fitted
#'   coordinates), with masses the shares of those draws and locations that
#'   match the point masses of the prior density when it is available. A
#'   quantity whose prior density has no point mass declares no atoms, and so
#'   does a quantity without a prior density none of whose fitted coordinates
#'   can take a point mass (each owned by a continuous prior or a continuous
#'   generated primitive, e.g. the correlations of an LKJ block). Atoms remain
#'   undeclared (and posterior plots and Savage-Dickey ratios stop) when the
#'   point states are not in the draws (an unmonitored component indicator)
#'   or the quantity combines several coordinates of which some can take a
#'   point mass (e.g. original-scale random-effect SDs and correlations of a
#'   block with spike-and-slab SD priors).}
#'   \item{`undefined_draws`}{for quantities that are undefined on some
#'   fitted draws (the catalog `definedness`, e.g. variance proportions when
#'   no allocation component is active), the reason; those draws are omitted,
#'   and the prior density and atom masses are conditional on definedness.}
#'   \item{`condition`}{the conditioning of the draws: unconditional
#'   (`averaged = TRUE`), or the inclusion event of the quantity.}
#'   \item{`quantities`}{the catalog quantity of the draws: its id and label
#'   parts, from which [parameter_labels()] and the summaries render its
#'   labels.}
#' }
#'
#' @param fit a model fitted with [JAGS_fit()].
#' @param selection a selection of one catalog quantity, created with
#'   [parameter_catalog_resolve()].
#' @param conditional whether to condition on the inclusion event of the
#'   quantity: the gates of its allocation chain on for an allocation-derived
#'   component SD or variance, its own inclusion gate on for a variance
#'   proportion, and at least one component active for an allocation total
#'   (`sd_total`, `var_total`) whose components all have inclusion gates.
#'   Only the draws in the event are kept, and the prior density is
#'   restricted to the event and renormalized (point masses inside the event,
#'   e.g. a variance proportion of one, keep their renormalized masses). The
#'   declared `condition` names the gates of the event. Quantities without an
#'   inclusion gate cannot be conditioned, nor can a total with an ungated
#'   component, for which some component is always active.
#'
#' @return A numeric vector of class `mixed_posteriors`,
#'   `mixed_posteriors.simple`, `marginal_posterior.simple`, and
#'   `marginal_posterior` (attribute `parameter`: the canonical name) with the
#'   draw metadata described above.
#'
#' @seealso [parameter_catalog()], [parameter_prior_density()],
#'   [posterior_metadata()], [posterior_atoms_free()]
#' @export
parameter_mixed_posterior <- function(fit, selection, conditional = FALSE){

  if(!inherits(fit, "BayesTools_fit")){
    stop("'fit' must be a 'BayesTools_fit' object.", call. = FALSE)
  }
  check_bool(conditional, "conditional", allow_NA = FALSE)
  catalog <- parameter_catalog(fit)
  .bt_validate_parameter_selection(selection, catalog = catalog)
  if(nrow(selection$quantities) != 1L){
    stop("'selection' must contain exactly one parameter quantity.",
         call. = FALSE)
  }

  .bt_parameter_mixed_posterior(
    fit         = fit,
    selection   = selection,
    conditional = conditional
  )
}

#' Per-draw gate and point states of a catalog quantity
#'
#' @description `parameter_gate_states()` returns, for each draw, the state
#' of the structural gates and point components that put one catalog quantity
#' on its atoms: the inclusion gates of random-effect variance allocations
#' (allocated SDs and variances, allocation totals, variance proportions, and
#' inclusion indicators), and the component indicator of a mixture or
#' spike-and-slab prior with point components (coefficients, factor levels,
#' and random-effect SDs and variances of such priors). These are the states
#' from which [parameter_mixed_posterior()] declares the posterior atoms.
#'
#' @param fit a model fitted with [JAGS_fit()].
#' @param selection a selection of one catalog quantity, created with
#'   [parameter_catalog_resolve()].
#' @param draws optional draws of the fitted coordinates to read the states
#'   from: a numeric matrix with one column per coordinate and one row per
#'   draw, an `mcmc` or `mcmc.list` object, or a named numeric vector holding
#'   one draw (as for [JAGS_evaluate_deterministic()]). It must contain the
#'   gate and indicator coordinates of the quantity. Defaults to the
#'   posterior draws of `fit`.
#'
#' @return `NULL` for a quantity without inclusion gates or point components.
#'   Otherwise a list with one element per draw in each of
#'   \describe{
#'     \item{`atom`}{the atom of the quantity in the draw, `NA` where the
#'     quantity is on its continuous part or undefined.}
#'     \item{`continuous`}{whether the quantity is on its continuous part
#'     (defined and not on an atom; `NA` where that is not determined).}
#'     \item{`defined`}{whether the quantity is defined in the draw (a
#'     variance proportion is undefined without an active component).}
#'     \item{`event`}{whether the draw is in the inclusion event of the
#'     quantity (see [parameter_mixed_posterior()]), or `NULL` when the
#'     quantity has no such event.}
#'   }
#'   and the scalar `known`: `FALSE` when the atom states are not determined
#'   by the draws (the component indicator of a prior with point components
#'   is not monitored); `atom` then holds only the gate atoms.
#'
#' @seealso [parameter_mixed_posterior()], [parameter_catalog()]
#' @export
parameter_gate_states <- function(fit, selection, draws = NULL){

  if(!inherits(fit, "BayesTools_fit")){
    stop("'fit' must be a 'BayesTools_fit' object.", call. = FALSE)
  }
  catalog <- parameter_catalog(fit)
  .bt_validate_parameter_selection(selection, catalog = catalog)
  if(nrow(selection$quantities) != 1L){
    stop("'selection' must contain exactly one parameter quantity.",
         call. = FALSE)
  }
  quantity <- selection$quantities[1L, , drop = FALSE]
  plan <- .bt_parameter_gate_plan(fit, quantity)
  if(is.null(plan)){
    return(NULL)
  }
  if(is.null(draws)){
    n <- nrow(as.matrix(.fit_to_posterior(fit)))
    states <- .bt_parameter_gate_states(fit, plan, n)
  }else{
    draws <- .bt_deterministic_draws_matrix(draws)
    states <- .bt_parameter_gate_states(fit, plan, nrow(draws), model_samples = draws)
  }

  continuous <- states$defined & is.na(states$atom)
  if(!isTRUE(states$known)){
    # a draw off the gate atoms may lie on an unmonitored point component
    continuous[continuous] <- NA
  }

  list(
    atom       = states$atom,
    continuous = continuous,
    defined    = states$defined,
    event      = states$event,
    known      = states$known
  )
}

.bt_parameter_mixed_posterior <- function(fit, selection, conditional = FALSE,
                                          n_grid = .prior_linear_density_default_grid(),
                                          tail_prob = .prior_linear_density_tail_prob(),
                                          simplify_label = FALSE){

  quantity <- selection$quantities[1L, , drop = FALSE]
  name <- quantity$canonical_name
  draws <- .bt_parameter_draws_from_quantities(fit, quantity)
  values <- unname(as.numeric(as.matrix(draws)[, 1L]))
  plan <- .bt_parameter_gate_plan(fit, quantity)
  states <- .bt_parameter_gate_states(fit, plan, length(values))

  keep <- !is.na(values)
  if(isTRUE(conditional)){
    if(is.null(plan) || length(plan$event_gates) == 0L){
      stop(
        "The inclusion event of '", name, "' is unavailable: the quantity has ",
        "no inclusion gate. Use 'conditional = FALSE'.",
        call. = FALSE
      )
    }
    if(is.null(states$event)){
      stop(
        "The inclusion event of '", name, "' is unavailable: it combines ",
        "parent-allocation and component inclusion gates. Use ",
        "'conditional = FALSE'.",
        call. = FALSE
      )
    }
    keep <- keep & states$event
  }
  if(!any(keep)){
    stop(
      "The posterior of '", name, "' is unavailable: no posterior draw ",
      if(isTRUE(conditional)) "lies in its inclusion event" else
        "defines the quantity", ".",
      call. = FALSE
    )
  }

  prior_density <- .bt_parameter_prior_density_quantity(
    object      = fit,
    selection   = selection,
    n_grid      = n_grid,
    tail_prob   = tail_prob,
    conditional = conditional
  )
  atoms <- .bt_parameter_mixed_posterior_atoms(
    fit           = fit,
    quantity      = quantity,
    prior_density = prior_density,
    plan          = plan,
    states        = states,
    keep          = keep
  )

  out <- values[keep]
  attr(out, "parameter")  <- name
  attr(out, "prior_list") <- prior_none()
  out <- .posterior_support_set(out, quantity$support[[1L]])
  if(!is.null(prior_density)){
    out <- .bt_meta_set(out, "prior_density", prior_density)
  }
  if(!is.null(atoms)){
    out <- .posterior_atoms_set(out, atoms)
  }
  if(!identical(quantity$definedness, "always")){
    out <- .bt_meta_set(out, "undefined_draws",
                        stats::setNames(quantity$definedness, name))
  }
  if(is.character(quantity$formula_parameter) &&
     !is.na(quantity$formula_parameter) &&
     nzchar(quantity$formula_parameter)){
    out <- .bt_meta_set(out, "formula_parameter", quantity$formula_parameter)
  }
  out <- .bt_meta_set(out, "condition", .bt_parameter_mixed_posterior_condition(
    plan, conditional
  ))
  # the draws are the catalog quantity: its id and label parts
  quantities <- .bt_catalog_quantity_table(
    quantity,
    column   = name,
    simplify = simplify_label
  )
  if(!is.null(quantities)){
    out <- .bt_meta_set(out, "quantities", quantities)
  }
  class(out) <- c(
    "mixed_posteriors",
    "mixed_posteriors.simple",
    "marginal_posterior.simple",
    "marginal_posterior"
  )
  out
}

# The 'condition' metadata of a catalog mixed posterior: unconditional, or
# the inclusion event of the quantity over its gate indicators (all gates of a
# component chain or a proportion's own gate: "AND"; any component of a
# total: "OR").
.bt_parameter_mixed_posterior_condition <- function(plan, conditional){

  if(!isTRUE(conditional)){
    return(list(
      conditional      = character(),
      conditional_rule = "AND",
      condition_key    = .condition_event_key(character()),
      averaged         = TRUE
    ))
  }

  list(
    conditional      = plan$event_gates,
    conditional_rule = plan$event_rule,
    condition_key    = .condition_event_key(plan$event_gates, plan$event_rule),
    averaged         = FALSE
  )
}

# Declared posterior atoms of a catalog quantity, from its structure:
# - a structural quantity is its fixed value;
# - a quantity with inclusion gates or point components (a gate plan) has
#   the atoms its per-draw gate and indicator states put it on, with masses
#   the shares of those draws (checked against the point masses of the prior
#   density when it is available);
# - a quantity whose prior density has no point mass, or, without a prior
#   density, none of whose fitted coordinates can take a point mass, has none;
# - otherwise (the point states are not in the draws, or a composite of
#   coordinates with point masses) the atom status is undeclared (NULL).
.bt_parameter_mixed_posterior_atoms <- function(fit, quantity, prior_density,
                                                plan, states, keep){

  name <- quantity$canonical_name
  if(identical(quantity$status, "structural") &&
     is.numeric(quantity$fixed_value) && length(quantity$fixed_value) == 1L &&
     is.finite(quantity$fixed_value)){
    return(.posterior_atoms_new(
      locations = matrix(quantity$fixed_value, ncol = 1L),
      mass      = 1,
      source    = "parameter_structure"
    ))
  }
  prior_atoms <- NULL
  if(!is.null(prior_density)){
    points <- prior_density$points
    prior_atoms <- if(is.data.frame(points) && nrow(points) > 0L){
      points$x[points$p > 0]
    }else{
      numeric()
    }
    if(length(prior_atoms) == 0L){
      return(.posterior_atoms_new(source = "parameter_structure"))
    }
  }
  if(is.null(plan)){
    if(is.null(prior_density) && .bt_parameter_point_free(fit, quantity)){
      return(.posterior_atoms_new(source = "parameter_structure"))
    }
    return(NULL)
  }
  if(is.null(states) || !isTRUE(states$known)){
    return(NULL)
  }

  locations <- states$atom[keep]
  on_atom <- !is.na(locations)
  if(!is.null(prior_atoms)){
    matched <- vapply(locations[on_atom], function(location){
      distance <- abs(prior_atoms - location)
      index <- which(distance <= 1e-12 * pmax(1, abs(prior_atoms)))
      if(length(index) == 0L) NA_real_ else prior_atoms[[index[[1L]]]]
    }, numeric(1))
    if(anyNA(matched)){
      stop(
        "The inclusion-gate atoms of '", name, "' do not match the point ",
        "masses of its prior density.",
        call. = FALSE
      )
    }
    locations[on_atom] <- matched
  }
  atom_locations <- sort(unique(locations[on_atom]))
  if(length(atom_locations) == 0L){
    return(.posterior_atoms_new(source = "parameter_gate_states"))
  }

  .posterior_atoms_new(
    locations = matrix(atom_locations, ncol = 1L),
    mass      = vapply(atom_locations, function(location){
      mean(on_atom & locations == location)
    }, numeric(1)),
    source    = "parameter_gate_states"
  )
}

# Whether none of the fitted coordinates a catalog quantity is computed from
# (the dependencies of its extraction key) can take a point mass: each is
# owned by a prior-list entry without a point component, or is a generated
# continuous primitive without a prior-list entry (the Beta primitives of an
# LKJ correlation and standardized random effects). Structural and indicator
# coordinates, and coordinates without either structure, can.
.bt_parameter_point_free <- function(fit, quantity){

  dependencies <- quantity$extraction_key[[1L]]$dependencies
  if(length(dependencies) == 0L){
    return(FALSE)
  }
  coordinates <- parameter_coordinates(fit)
  rows <- match(dependencies, coordinates$coordinate_name)
  if(anyNA(rows) ||
     any(!coordinates$convergence_role[rows] %in% c("sampled", "derived"))){
    return(FALSE)
  }
  prior_list <- attr(fit, "prior_list", exact = TRUE)
  for(i in seq_along(dependencies)){
    owner <- .bt_parameter_coordinate_owner(prior_list, dependencies[[i]])
    if(!is.null(owner)){
      if(!isFALSE(.bt_prior_has_point_component(prior_list[[owner]]))){
        return(FALSE)
      }
    }else if(!coordinates$role[rows[i]] %in%
             c("random_correlation_coordinate", "random_latent") ||
             !identical(coordinates$convergence_role[rows[i]], "sampled")){
      return(FALSE)
    }
  }

  TRUE
}

# The prior-list entry whose fitted coordinates include 'coordinate' (NULL
# when none does).
.bt_parameter_coordinate_owner <- function(prior_list, coordinate){

  for(parameter in names(prior_list)){
    prior <- prior_list[[parameter]]
    if(is.prior(prior) &&
       coordinate %in% .prior_linear_prior_columns(parameter, prior)){
      return(parameter)
    }
  }

  NULL
}

# Whether a prior puts point mass somewhere: TRUE for point priors, the spike
# of spike-and-slab priors, mixtures with such a component, and ordered
# priors with a zero total; FALSE for continuous priors; NA when the prior's
# structure is not classified here.
.bt_prior_has_point_component <- function(prior){

  if(!is.prior(prior)){
    return(NA)
  }
  if(is.prior.point(prior) || is.prior.spike_and_slab(prior)){
    return(TRUE)
  }
  if(is.prior.mixture(prior)){
    components <- vapply(seq_along(prior), function(i){
      .bt_prior_has_point_component(prior[[i]])
    }, logical(1))
    if(any(components %in% TRUE)){
      return(TRUE)
    }
    return(if(anyNA(components)) NA else FALSE)
  }
  if(is.prior.ordered(prior)){
    return(.posterior_atoms_is_ordered_zero_total(prior) ||
             .posterior_atoms_ordered_total_has_spike(prior$total))
  }
  if(is.prior.simple(prior) || is.prior.vector(prior) || is.prior.factor(prior)){
    return(FALSE)
  }

  NA
}

# The gate structure of a catalog quantity (NULL for quantities without
# inclusion gates or point components): what puts the quantity on an atom in
# a draw, and its inclusion event.
#   inclusion: an allocation inclusion indicator, on its own atom (0 or 1).
#   component: an allocation-derived component SD (or variance), zero when a
#     gate of its chain is off; 'degenerate' when every factor of the chain
#     is a gate without Dirichlet weights (a gate-only allocation).
#   total: an allocation total (sd_total, var_total, or the common SD and
#     variance), zero without an active component or with a parent gate off.
#   var_prop: a gated total-variance proportion, 0 when its component is
#     inactive while another is active and 1 when it is the only active one.
#   point: a quantity of the fitted coordinates of one mixture or
#     spike-and-slab prior with point components (a coefficient, a factor
#     level, or a random-effect SD or its variance), on the image of a point
#     component in the draws whose component indicator selects it.
.bt_parameter_gate_plan <- function(fit, quantity){

  key <- quantity$extraction_key[[1L]]
  if(!identical(key$type, "random_summary")){
    return(.bt_parameter_point_plan(fit, quantity))
  }
  gate_names <- function(records, field){
    names <- unlist(lapply(records, function(record) record[[field]]),
                    use.names = FALSE)
    names[!is.na(names) & nzchar(names)]
  }

  if(identical(key$evaluator, "allocation_inclusion")){
    # the indicator is its own gate state
    return(list(
      kind        = "inclusion",
      chain_gates = key$source_parameter,
      event_gates = character(),
      event_rule  = "AND"
    ))
  }
  if(key$evaluator %in% c("sd", "sd_variance") &&
     (isTRUE(key$allocation_derived) ||
        .bt_parameter_prior_density_gate_only_sd(fit, key))){
    chain <- .bt_parameter_prior_density_component_chain(fit, key)
    if(is.null(chain)){
      return(NULL)
    }
    gates <- gate_names(chain$factors, "inclusion_name")
    return(list(
      kind        = "component",
      square      = identical(key$evaluator, "sd_variance"),
      source      = chain$allocation$source,
      degenerate  = all(vapply(chain$factors, function(factor){
        is.null(factor$weight_name)
      }, logical(1))),
      chain_gates = gates,
      event_gates = gates,
      event_rule  = "AND"
    ))
  }

  if(!key$evaluator %in% c("allocation_sd", "allocation_var") &&
     !(identical(key$evaluator, "allocation") &&
       identical(quantity$quantity, "var_prop"))){
    return(.bt_parameter_point_plan(fit, quantity))
  }
  random_term <- if(nzchar(key$random_block)){
    .bt_parameter_catalog_find_random_term(fit, key)
  }else{
    NULL
  }
  allocation <- .bt_parameter_catalog_find_allocation(fit, key, random_term)
  total_variance <- identical(.bt_random_effect_allocation_scale_metadata(
    allocation,
    context = "Random-effect allocation metadata"
  ), "total_variance")
  component_gates <- if(total_variance){
    gate_names(allocation$inclusion, "indicator_name")
  }else{
    character()
  }
  parent_gates <- gate_names(allocation$parent_factors, "inclusion_name")

  if(identical(key$evaluator, "allocation")){
    if(length(component_gates) == 0L && length(parent_gates) == 0L){
      return(NULL)
    }
    own <- Filter(function(record){
      isTRUE(record$index == key$index)
    }, allocation$inclusion)
    own_gate <- gate_names(own, "indicator_name")
    return(list(
      kind            = "var_prop",
      allocation      = allocation,
      index           = as.integer(key$index),
      component_gates = component_gates,
      parent_gates    = parent_gates,
      event_gates     = own_gate,
      event_rule      = "AND"
    ))
  }

  # the event "some component active" is uncertain only when every
  # component has an inclusion gate; with an ungated component it is certain
  # and only the parent gates condition the total
  gated_indices <- unlist(lapply(allocation$inclusion, `[[`, "index"), use.names = FALSE)
  all_gated <- length(component_gates) > 0L && (
    isTRUE(allocation$gate_only) ||
      setequal(gated_indices, seq_len(allocation$n_targets))
  )
  list(
    kind              = "total",
    square            = identical(key$evaluator, "allocation_var"),
    source            = allocation$source,
    allocation        = allocation,
    total_variance    = total_variance,
    component_gates   = component_gates,
    parent_gates      = parent_gates,
    all_gated         = all_gated,
    parent_degenerate = all(vapply(allocation$parent_factors, function(factor){
      is.null(factor$weight_name)
    }, logical(1))),
    event_gates       = c(if(all_gated) component_gates, parent_gates),
    event_rule        = if(all_gated && length(parent_gates) == 0L) "OR" else "AND"
  )
}

# Per-draw states of a gate plan: 'atom' is the atom location of the quantity
# in each draw (NA on its continuous part and where it is undefined),
# 'defined' whether the quantity is defined in the draw (a variance
# proportion needs an active component), 'event' its inclusion event (NULL
# when unavailable), and 'known' whether the atom states are determined (the
# mixture indicator of a prior with a point component is monitored). The
# states are read from the posterior draws of 'fit', or from 'model_samples'
# (a matrix of fitted coordinates with one row per draw).
.bt_parameter_gate_states <- function(fit, plan, n, model_samples = NULL){

  if(is.null(plan)){
    return(NULL)
  }
  gate_columns <- unique(c(
    plan$chain_gates, plan$component_gates, plan$parent_gates
  ))
  source <- if(plan$kind %in% c("component", "total")){
    .bt_parameter_source_point_plan(fit, plan$source)
  }else{
    NULL
  }
  columns <- unique(c(gate_columns, source$indicator, plan$indicator))
  if(is.null(model_samples)){
    model_samples <- if(length(columns) > 0L){
      as.matrix(.bt_parameter_draw_dependencies(fit, columns))
    }else{
      matrix(numeric(), nrow = n, ncol = 0L)
    }
  }else{
    missing <- setdiff(columns, colnames(model_samples))
    if(length(missing) > 0L){
      stop(
        "The draws do not contain the inclusion indicators ",
        paste0("'", missing, "'", collapse = ", "), " of the quantity.",
        call. = FALSE
      )
    }
  }
  if(nrow(model_samples) != n){
    stop("Random-effect inclusion draws do not align with the quantity draws.",
         call. = FALSE)
  }
  defined <- rep(TRUE, n)
  gate_on <- function(name){
    gate <- .bt_random_effect_allocation_gate_draws(
      parameter_name = name,
      posterior = model_samples
    )
    if(is.null(gate)){
      stop(
        "Random-effect allocation inclusion samples are missing Bernoulli indicator '",
        name, "'.",
        call. = FALSE
      )
    }
    gate == 1
  }
  all_on <- function(names){
    out <- rep(TRUE, n)
    for(name in names){
      out <- out & gate_on(name)
    }
    out
  }

  if(identical(plan$kind, "inclusion")){
    return(list(
      atom    = as.numeric(gate_on(plan$chain_gates)),
      defined = defined,
      event   = NULL,
      known   = TRUE
    ))
  }
  if(identical(plan$kind, "point")){
    if(isTRUE(plan$unknown)){
      return(list(atom = rep(NA_real_, n), defined = defined, event = NULL,
                  known = FALSE))
    }
    component <- .bt_component_from_indicator(plan$prior, model_samples[, plan$indicator])
    return(list(
      atom    = plan$locations[component],
      defined = defined,
      event   = NULL,
      known   = TRUE
    ))
  }
  if(identical(plan$kind, "var_prop")){
    gates <- .bt_random_effect_summary_allocation_component_gates(
      allocation    = plan$allocation,
      K             = plan$allocation$n_targets,
      model_samples = model_samples
    )
    active <- gates$component_gates == 1 & gates$parent_active
    own <- active[, plan$index]
    others <- rowSums(active[, -plan$index, drop = FALSE]) > 0L
    return(list(
      atom    = ifelse(!own & others, 0, ifelse(own & !others, 1, NA_real_)),
      defined = own | others,
      event   = if(length(plan$event_gates) > 0L) all_on(plan$event_gates),
      known   = TRUE
    ))
  }

  point <- .bt_parameter_source_point_draws(source, model_samples, n)
  if(identical(plan$kind, "component")){
    on <- all_on(plan$chain_gates)
    continuous_multiplier <- !isTRUE(plan$degenerate)
    event <- if(length(plan$chain_gates) > 0L) on
  }else{
    parent_active <- all_on(plan$parent_gates)
    any_on <- full <- rep(TRUE, n)
    if(length(plan$component_gates) > 0L){
      if(isTRUE(plan$allocation$gate_only)){
        # a gate-only total is the source times its gates
        any_on <- full <- all_on(plan$component_gates)
      }else{
        # components without a gate are always active
        states <- .bt_random_effect_summary_allocation_component_gates(
          allocation    = plan$allocation,
          K             = plan$allocation$n_targets,
          model_samples = model_samples
        )$component_gates == 1
        any_on <- rowSums(states) > 0L
        full <- rowSums(!states) == 0L
      }
    }
    on <- parent_active & any_on
    continuous_multiplier <- !(isTRUE(plan$parent_degenerate) & full)
    # the declared event: the component gates (OR) when every component is
    # gated, the parent gates (AND) otherwise; both kinds together are not
    # one conjunction or disjunction of gates
    event <- if(!isTRUE(plan$all_gated)){
      parent_active
    }else if(length(plan$parent_gates) == 0L){
      any_on
    }
  }
  atom <- if(is.null(point)){
    ifelse(on, NA_real_, 0)
  }else{
    ifelse(!on | (!is.na(point) & point == 0), 0,
           ifelse(!is.na(point) & !continuous_multiplier, point, NA_real_))
  }
  if(isTRUE(plan$square)){
    atom <- atom^2
  }

  list(atom = atom, defined = defined, event = event, known = !is.null(point))
}

# The point plan of a quantity of the fitted coordinates of one mixture or
# spike-and-slab prior with point components (a coefficient or a factor level,
# the weighted sum of its coordinates, or a random-effect SD and its one-to-one
# transformations): the image of each component (NA for continuous
# components) and the monitored component indicator of the prior. NULL for
# other quantities.
.bt_parameter_point_plan <- function(fit, quantity){

  key <- quantity$extraction_key[[1L]]
  transform <- NULL
  if(key$type %in% c("coordinate", "factor_level")){
    dependencies <- key$dependencies
    weights <- if(identical(key$type, "coordinate")) rep(1, length(dependencies)) else key$weights
  }else if(identical(key$type, "random_summary") &&
           key$source_type %in% c("identity", "one_to_one_transform") &&
           is.character(key$source_parameter) && length(key$source_parameter) == 1L &&
           nzchar(key$source_parameter)){
    dependencies <- key$source_parameter
    weights <- 1
    if(identical(key$source_type, "one_to_one_transform")){
      transform <- .bt_parameter_transform_from_quantity(fit, quantity)
      if(is.null(transform)){
        return(NULL)
      }
    }
  }else{
    return(NULL)
  }
  if(length(dependencies) == 0L){
    return(NULL)
  }

  prior_list <- attr(fit, "prior_list", exact = TRUE)
  owners <- unique(vapply(dependencies, function(dependency){
    owner <- .bt_parameter_coordinate_owner(prior_list, dependency)
    if(is.null(owner)) NA_character_ else owner
  }, character(1)))
  if(length(owners) != 1L || is.na(owners)){
    return(NULL)
  }
  prior <- prior_list[[owners]]
  if(!.posterior_components_is_mixture(prior)){
    return(NULL)
  }
  columns <- .prior_linear_prior_columns(owners, prior)
  positions <- match(dependencies, columns)
  locations <- vapply(seq_along(prior), function(i){
    location <- .posterior_atoms_point_location(prior[[i]], length(columns))
    if(is.null(location)) NA_real_ else sum(weights * location[positions])
  }, numeric(1))
  if(all(is.na(locations))){
    return(NULL)
  }
  if(!is.null(transform)){
    locations[!is.na(locations)] <- parameter_transform_forward(
      locations[!is.na(locations)],
      transform
    )
  }
  indicator <- paste0(owners, "_indicator")
  coordinates <- parameter_coordinates(fit)
  status <- coordinates$monitor_status[match(indicator, coordinates$coordinate_name)]
  known <- !is.na(status) && status %in% c("sampled", "structural")

  list(
    kind        = "point",
    prior       = prior,
    locations   = locations,
    indicator   = if(known) indicator,
    unknown     = !known,
    chain_gates = character(),
    event_gates = character(),
    event_rule  = "AND"
  )
}

# How the point components of a scale prior enter the draws: the location of
# each component that is a point prior (NA for continuous components) and the
# name of the monitored mixture indicator; NULL for sources without a prior.
.bt_parameter_source_point_plan <- function(fit, source){

  prior <- source$prior
  if(!is.prior(prior)){
    return(NULL)
  }
  if(is.prior.point(prior)){
    return(list(locations = prior$parameters[["location"]], indicator = NULL,
                prior = prior))
  }
  if(is.prior.simple(prior)){
    return(list(locations = NA_real_, indicator = NULL, prior = prior))
  }
  if(!is.prior.mixture(prior) && !is.prior.spike_and_slab(prior)){
    return(NULL)
  }
  locations <- vapply(seq_along(prior), function(i){
    component <- prior[[i]]
    if(is.prior.point(component)) component$parameters[["location"]] else NA_real_
  }, numeric(1))
  indicator <- paste0(.bt_random_sd_binding_source_name(source), "_indicator")
  coordinates <- parameter_coordinates(fit)
  status <- coordinates$monitor_status[match(indicator, coordinates$coordinate_name)]
  if(all(is.na(locations))){
    return(list(locations = NA_real_, indicator = NULL, prior = prior))
  }
  if(is.na(status) || !status %in% c("sampled", "structural")){
    # the atom states of the source are not in the draws
    return(list(locations = locations, indicator = NULL, prior = prior,
                unknown = TRUE))
  }

  list(locations = locations, indicator = indicator, prior = prior)
}

# Per-draw point location of the scale source (NA on its continuous
# components), or NULL when the atom states of the source are unknown.
.bt_parameter_source_point_draws <- function(source_plan, model_samples, n){

  if(is.null(source_plan) || isTRUE(source_plan$unknown)){
    return(NULL)
  }
  if(is.null(source_plan$indicator)){
    return(rep(source_plan$locations[[1L]], n))
  }
  component <- .bt_component_from_indicator(
    source_plan$prior,
    model_samples[, source_plan$indicator]
  )

  source_plan$locations[component]
}
