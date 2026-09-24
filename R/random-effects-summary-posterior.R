#' Extract random-effect summary posterior distributions
#'
#' @description Extracts posterior draws for derived random-effect summary
#' quantities and attaches analytic marginal prior densities where available.
#' For allocation summaries, the helper keeps raw Dirichlet allocation weights
#' internal and exposes interpretable scalar summaries such as mean-variance
#' SD-component variance multipliers. Gated total-variance proportions are
#' normalized over active components and omit draws on which every component
#' is excluded because a variance share is undefined there.
#'
#' The posterior atoms of each summary are declared from its structure: none
#' when the quantity's prior density ([parameter_prior_density()]) has no
#' point mass; the inclusion-gate atoms of gated total-variance proportions (at
#' 0 and 1) and allocation totals (at 0), with masses from the posterior
#' inclusion-indicator draws. Other point masses (e.g. a scale prior with its
#' own spike) leave the atom status undeclared, and posterior plots of such
#' summaries stop instead of inferring point masses from the draws.
#'
#' @param fit model fit created by [JAGS_fit].
#' @param summary semantic quantity to extract: `"var_mult"`, `"var_prop"`,
#'   `"sd_mult"`, `"sd_total"`, `"var_total"`, `"sd_common"`, or
#'   `"var_common"`.
#' @param allocation optional allocation label filter.
#' @param component optional allocation component label filter.
#' @param formula_parameter optional formula parameter filter.
#' @param simplify_names whether returned quantities should use centrally
#'   generated simplified display labels. Defaults to `FALSE`.
#' @param n_prior_points number of grid points used for the attached analytic
#'   prior density.
#'
#' @return A named list of posterior vectors with class `mixed_posteriors` and
#'   `marginal_posterior`. The list can be passed to [plot_posterior], and the
#'   individual vectors carry a `prior_density` attribute.
#'
#' @export
random_effects_summary_posterior <- function(
    fit,
    summary = c(
      "var_mult", "var_prop", "sd_mult",
      "sd_total", "var_total", "sd_common", "var_common"
    ),
    allocation = NULL,
    component = NULL,
    formula_parameter = NULL,
    simplify_names = FALSE,
    n_prior_points = 4096){

  if(!inherits(fit, "runjags") || !inherits(fit, "BayesTools_fit")){
    stop("'fit' must be a BayesTools JAGS fit.", call. = FALSE)
  }
  summary <- .bt_random_effect_summary_posterior_type(summary)
  check_char(allocation, "allocation", check_length = 0, allow_NULL = TRUE, allow_NA = FALSE)
  check_char(component, "component", check_length = 0, allow_NULL = TRUE, allow_NA = FALSE)
  check_char(formula_parameter, "formula_parameter", check_length = 0, allow_NULL = TRUE, allow_NA = FALSE)
  check_bool(simplify_names, "simplify_names", allow_NA = FALSE)
  check_int(n_prior_points, "n_prior_points", lower = 16)

  prior_list <- attr(fit, "prior_list", exact = TRUE)
  check_list(prior_list, "prior_list")
  if(!all(vapply(prior_list, is.prior, logical(1)))){
    stop("'prior_list' must be a list of priors.", call. = FALSE)
  }

  catalog <- parameter_catalog(fit)
  quantities <- catalog$quantities
  keep <- startsWith(quantities$role, "random_") &
    !quantities$internal & quantities$status != "unavailable" &
    quantities$quantity == summary$summary
  keys <- quantities$extraction_key
  if(!is.null(allocation)){
    keep <- keep & vapply(keys, function(key){
      !is.null(key$allocation_label) && key$allocation_label %in% allocation
    }, logical(1))
  }
  if(!is.null(component)){
    keep <- keep & quantities$component %in% component
  }
  if(!is.null(formula_parameter)){
    keep <- keep & quantities$formula_parameter %in% formula_parameter
  }
  quantities <- quantities[keep, , drop = FALSE]
  if(nrow(quantities) == 0L && !any(startsWith(
    catalog$quantities$role,
    "random_"
  ))){
    stop("No random-effect summaries are available in 'fit'.", call. = FALSE)
  }
  if(nrow(quantities) == 0L){
    .bt_random_effect_summary_posterior_no_match_stop(
      summary = summary,
      allocation = allocation,
      component = component,
      formula_parameter = formula_parameter
    )
  }

  display_names <- if(simplify_names){
    quantities$display_label
  }else{
    quantities$canonical_name
  }

  out <- vector("list", nrow(quantities))
  names(out) <- display_names
  out_prior_list <- vector("list", nrow(quantities))
  names(out_prior_list) <- display_names

  for(i in seq_len(nrow(quantities))){
    quantity <- quantities[i, , drop = FALSE]
    key      <- quantity$extraction_key[[1L]]
    draws    <- .bt_parameter_draws_from_quantities(fit, quantity)
    values   <- unname(as.numeric(as.matrix(draws)[, 1L]))
    defined  <- rep(TRUE, length(values))
    if(identical(quantity$quantity, "var_prop")){
      defined <- !is.na(values)
      values  <- values[defined]
      if(length(values) == 0L){
        stop(
          "The selected variance proportion is undefined because no posterior draw has positive realized allocation variance.",
          call. = FALSE
        )
      }
    }
    atoms <- .bt_random_effect_summary_posterior_atoms(
      fit      = fit,
      catalog  = catalog,
      quantity = quantity,
      defined  = defined
    )
    values <- .bt_meta_set(values, "sample_ind", FALSE)
    values <- .bt_meta_set(values, "models_ind", rep(1, length(values)))
    attr(values, "parameter") <- display_names[i]
    attr(values, "summary_name") <- key$summary_name
    attr(values, "random_summary") <- summary$summary
    attr(values, "random_summary_label") <- display_names[i]
    attr(values, "random_allocation") <- key$allocation_label
    attr(values, "random_component") <- quantity$component
    values <- .bt_meta_set(values, "formula_parameter", quantity$formula_parameter)
    attr(values, "prior_list") <- prior_none()

    prior_density <- .bt_random_effect_summary_posterior_prior_density(
      fit      = fit,
      quantity = quantity,
      n_grid   = n_prior_points
    )
    if(!is.null(prior_density)){
      values <- .bt_meta_set(values, "prior_density", prior_density)
      values <- .posterior_support_set(
        values,
        .posterior_support_new(
          bounds = attr(prior_density, "support", exact = TRUE),
          exact = TRUE,
          source = "prior",
          type = "interval"
        )
      )
    }
    if(!is.null(atoms)){
      values <- .posterior_atoms_set(values, atoms)
    }

    class(values) <- c(
      "mixed_posteriors",
      "mixed_posteriors.simple",
      "marginal_posterior.simple",
      "marginal_posterior"
    )
    out[[i]] <- values
    out_prior_list[[i]] <- attr(values, "prior_list", exact = TRUE)
  }

  attr(out, "prior_list") <- out_prior_list
  attr(out, "summary_names") <- vapply(
    quantities$extraction_key,
    `[[`,
    character(1),
    "summary_name"
  )
  class(out) <- c("as_mixed_posteriors", "mixed_posteriors", "list")
  out
}

.bt_random_effect_summary_posterior_type <- function(summary){

  summary <- match.arg(
    summary,
    c(
      "var_mult", "var_prop", "sd_mult",
      "sd_total", "var_total", "sd_common", "var_common"
    )
  )
  list(input = summary, summary = summary)
}

.bt_random_effect_summary_posterior_no_match_stop <- function(summary,
                                                              allocation,
                                                              component,
                                                              formula_parameter){

  has_filters <- !is.null(allocation) || !is.null(component) ||
    !is.null(formula_parameter)
  filter_text <- if(has_filters) " matched the requested filters" else " are available"

  detail <- switch(
    summary$summary,
    "var_mult" = paste0(
      "Variance-multiplier summaries are created only for ",
      "random_variance_allocation(..., target = \"sd_component\", ",
      "scale = \"mean_variance\"). Total-variance allocations are returned ",
      "by summary = \"var_prop\"."
    ),
    "var_prop" = paste0(
      "Variance-proportion summaries are created for true total-variance ",
      "allocations. Mean-variance SD-component allocations are returned by ",
      "summary = \"var_mult\"."
    ),
    "sd_mult" = paste0(
      "SD-multiplier summaries are available only for SD-component ",
      "variance allocations."
    ),
    "Requested random-effect summaries are not available."
  )

  stop(
    "No random-effect ",
    summary$input,
    " summaries",
    filter_text,
    ". ",
    detail,
    call. = FALSE
  )
}

# Posterior atoms of a random-effect summary, declared from the structure of
# the quantity: no atoms when its canonical prior density
# (parameter_prior_density()) has no point mass; the inclusion-gate atoms of
# gated total-variance proportions and allocation totals, with masses from the
# posterior inclusion-indicator draws; undeclared (NULL) otherwise, so plots
# stop with the atom-status message. 'defined' marks the posterior draws kept
# in the summary (variance proportions omit draws without active components).
.bt_random_effect_summary_posterior_atoms <- function(fit, catalog, quantity,
                                                      defined){

  selection <- list(
    schema_version        = .bt_parameter_selection_version,
    parameter_map_version = catalog$schema_version,
    quantity_id           = quantity$quantity_id,
    quantities            = quantity
  )
  class(selection) <- c("BayesTools_parameter_selection", "list")
  prior_density <- parameter_prior_density(fit, selection)
  if(is.null(prior_density)){
    return(NULL)
  }
  prior_points <- prior_density$points
  if(is.null(prior_points) || !any(prior_points$p > 0)){
    return(.posterior_atoms_new(source = "random_summary_structure"))
  }

  locations <- .bt_random_effect_summary_gate_atom_locations(fit, quantity)
  if(is.null(locations)){
    return(NULL)
  }
  locations <- locations[defined]
  atom_locations <- sort(unique(locations[!is.na(locations)]))
  if(any(!atom_locations %in% prior_points$x[prior_points$p > 0])){
    stop(
      "Random-effect summary gate atoms do not match the point masses of the ",
      "prior density of '", quantity$canonical_name, "'.",
      call. = FALSE
    )
  }
  .posterior_atoms_new(
    locations = matrix(atom_locations, ncol = 1L),
    mass      = vapply(atom_locations, function(location){
      mean(!is.na(locations) & locations == location)
    }, numeric(1)),
    source    = "random_summary_inclusion_indicators"
  )
}

# Per-draw location of the inclusion-gate atom of a random-effect summary (NA
# on draws on its continuous part), from the fitted inclusion indicators:
# a total-variance proportion is 0 when its component is inactive while
# another component is active and 1 when it is the only active component; an
# allocation total without parent allocations is 0 when no component is
# active. NULL when the quantity's point masses are not all gate atoms (e.g.
# a scale prior with its own point mass, or nested allocations).
.bt_random_effect_summary_gate_atom_locations <- function(fit, quantity){

  key <- quantity$extraction_key[[1L]]
  gated_proportion <- identical(key$evaluator, "allocation") &&
    identical(quantity$quantity, "var_prop")
  gated_total <- key$evaluator %in% c("allocation_sd", "allocation_var") &&
    identical(key$source_type, "composite")
  if(!gated_proportion && !gated_total){
    return(NULL)
  }

  random_term <- if(nzchar(key$random_block)){
    .bt_parameter_catalog_find_random_term(fit, key)
  }else{
    NULL
  }
  allocation <- .bt_parameter_catalog_find_allocation(fit, key, random_term)
  if(!identical(.bt_random_effect_allocation_scale_metadata(
    allocation,
    context = "Random-effect summary atom metadata"
  ), "total_variance")){
    return(NULL)
  }
  if(gated_total){
    source_prior <- allocation$source$prior
    continuous_source <- is.prior(source_prior) && is.prior.simple(source_prior) &&
      !is.prior.point(source_prior) && !is.prior.discrete(source_prior) &&
      !is.prior.mixture(source_prior) && !is.prior.spike_and_slab(source_prior)
    if(!continuous_source || length(allocation$parent_factors) > 0L){
      return(NULL)
    }
  }

  model_samples <- as.matrix(.bt_parameter_draw_dependencies(
    fit,
    .bt_random_effect_summary_allocation_gate_names(allocation)
  ))
  gates <- .bt_random_effect_summary_allocation_component_gates(
    allocation    = allocation,
    K             = allocation$n_targets,
    model_samples = model_samples
  )
  active <- gates$component_gates == 1 & gates$parent_active

  if(gated_total){
    return(ifelse(rowSums(active) == 0L, 0, NA_real_))
  }
  index <- as.integer(key$index)
  others_active <- rowSums(active[, -index, drop = FALSE]) > 0L
  ifelse(
    !active[, index] & others_active, 0,
    ifelse(active[, index] & !others_active, 1, NA_real_)
  )
}

.bt_random_effect_summary_posterior_prior_density <- function(
    fit,
    quantity,
    n_grid){

  key <- quantity$extraction_key[[1L]]
  if(is.null(key$allocation_label) || is.null(key$index)){
    return(NULL)
  }
  random_term <- if(nzchar(key$random_block)){
    .bt_parameter_catalog_find_random_term(fit, key)
  }else{
    NULL
  }
  allocation <- .bt_parameter_catalog_find_allocation(fit, key, random_term)
  index <- key$index
  if(identical(quantity$quantity, "sd_mult") &&
     index > allocation$n_targets){
    index <- index - allocation$n_targets
  }
  if(is.null(allocation) || is.null(index) || length(index) != 1L || is.na(index)){
    return(NULL)
  }

  prior_list <- attr(fit, "prior_list", exact = TRUE)
  weight_prior <- prior_list[[allocation$weight_name]]
  if(is.null(weight_prior) ||
     !is.prior.simplex(weight_prior) ||
     !identical(weight_prior$distribution, "dirichlet")){
    return(NULL)
  }

  alpha <- weight_prior$parameters[["alpha"]]
  if(index < 1L || index > length(alpha)){
    return(NULL)
  }

  alpha_i <- alpha[index]
  beta_i <- sum(alpha) - alpha_i
  if(!is.finite(alpha_i) || !is.finite(beta_i) || alpha_i <= 0 || beta_i <= 0){
    return(NULL)
  }

  semantic_transform <- .bt_parameter_transform_from_quantity(fit, quantity)
  density_transform <- if(identical(
    semantic_transform,
    list(type = "identity")
  )){
    list(scale = 1, type = "linear")
  }else if(is.list(semantic_transform) &&
           identical(semantic_transform$type, "affine") &&
           identical(semantic_transform$offset, 0) &&
           semantic_transform$scale > 0){
    list(scale = semantic_transform$scale, type = "linear")
  }else if(is.list(semantic_transform) &&
           identical(semantic_transform$type, "sqrt_scale")){
    list(scale = semantic_transform$scale, type = "sqrt")
  }else{
    return(NULL)
  }

  .bt_random_effect_summary_posterior_scaled_beta_density(
    alpha = alpha_i,
    beta = beta_i,
    scale = density_transform$scale,
    transform = density_transform$type,
    n_grid = n_grid
  )
}

.bt_random_effect_summary_posterior_scaled_beta_density <- function(alpha, beta,
                                                                    scale,
                                                                    transform,
                                                                    n_grid){

  # scale * w (linear) or sqrt(scale * w) = exp(log(scale) / 2) w^(1/2)
  # (exp_lin) of w ~ Beta(alpha, beta): the recorded provenance evaluates
  # heights, probabilities and the plotting grid on the structural route.
  root <- identical(transform, "sqrt")
  tail_prob <- .prior_linear_density_tail_prob()
  provenance <- list(
    kind = "linear_combination",
    arguments = list(
      prior_list = list(source = prior("beta", list(alpha = alpha, beta = beta))),
      weights = c(source = if(root) 1 else scale),
      n_grid = n_grid,
      tail_prob = tail_prob,
      source_transforms = NULL,
      output_transformation = if(root) "exp_lin" else NULL,
      output_transformation_arguments = if(root) list(a = log(scale) / 2, b = 1 / 2) else NULL,
      grid_spacing = NULL
    )
  )
  route <- .prior_density_route_from_adaptive(provenance)

  if(root){
    support <- c(0, sqrt(scale))
    lower <- if(alpha < 0.5) sqrt(scale * stats::qbeta(tail_prob / 2, alpha, beta)) else support[1]
    upper <- if(beta < 1) sqrt(scale * stats::qbeta(1 - tail_prob / 2, alpha, beta)) else support[2]
  }else{
    support <- c(0, scale)
    lower <- if(alpha < 1) scale * stats::qbeta(tail_prob / 2, alpha, beta) else support[1]
    upper <- if(beta < 1) scale * stats::qbeta(1 - tail_prob / 2, alpha, beta) else support[2]
  }

  if(!is.finite(lower) || !is.finite(upper) || lower >= upper){
    stop(
      "Could not construct a finite plotting grid for the scaled-Beta prior.",
      call. = FALSE
    )
  }

  x <- seq(lower, upper, length.out = n_grid)
  # Keep the unit reference value on the grid only inside the plotted range;
  # outside it, x = 1 can be a singular support boundary (scale = 1, beta < 1).
  if(lower < 1 && 1 < upper){
    x <- sort(unique(c(x, 1)))
  }

  y <- .prior_density_route_density(route, x)
  if(anyNA(y) || any(!is.finite(y))){
    stop(
      "The scaled-Beta plotting grid includes a singular support boundary.",
      call. = FALSE
    )
  }

  out <- list(
    density = list(
      x = x,
      y = y,
      mass = 1
    ),
    points = .prior_linear_density_empty_points(),
    n_grid = length(x)
  )
  class(out) <- c("prior_linear_density", "prior_density")
  attr(out, "support") <- support
  attr(out, "adaptive_evaluation") <- provenance
  attr(out, "singular_boundaries") <- support[vapply(support, function(bound){
    identical(prior_density_ordinate(out, bound)$behavior, "infinite")
  }, logical(1))]
  out
}
