# Within a centered unit, the coefficient vector, conditional on its
# scales, is spherical Normal. A linear direction therefore uses one
# standard Normal and the Euclidean norm of its coded weights. For group
# lasso retain the dimension of the ORIGINAL unit in the Gamma mixing law.
.bt_regularized_direction <- function(prior, unit, weights){

  spec <- prior$regularization
  largest <- max(abs(weights))
  norm <- largest * sqrt(sum((weights / largest)^2))
  if(!is.finite(norm)) .regularized_numerical_unavailable("direction norm")
  spec$multiplier <- .prior_linear_fold_scale(norm, spec$multiplier)
  spec$dimension <- 1L
  spec$unit_index <- 1L
  spec$unit_dimension <- spec$unit_dimension[unit]
  spec$unit_centered <- spec$unit_centered[unit]
  if(length(spec$local_names)) spec$local_names <- spec$local_names[unit]
  if(length(spec$variance_names)) spec$variance_names <- spec$variance_names[unit]
  if(spec$family == "lasso" && spec$unit_centered && spec$unit_dimension > 1){
    spec$group_lasso_dimension <- spec$unit_dimension
  }
  prior$regularization <- spec
  multiplier <- attr(prior, "multiply_by", exact = TRUE)
  attributes(prior) <- list(names = names(prior), class = c("prior", "prior.regularized_term", "prior.simple"))
  attr(prior, "multiply_by") <- multiplier
  prior
}

.bt_has_regularized_centered <- function(prior){

  if(is.prior.mixture(prior)) return(any(vapply(prior, .bt_has_regularized_centered, logical(1))))
  inherits(prior, "prior.regularized_term") && is.prior.factor(prior) && any(prior$regularization$unit_centered)
}

# A factor mixture has one gate for the entire term. Project each component
# inside that gate. Splitting the term into independently gated directions
# would change the atom probability.
.bt_regularized_project_component <- function(prior, indices, weights){

  if(is.prior.spike_and_slab(prior)){
    slab <- .bt_regularized_project_component(.get_spike_and_slab_variable(prior), indices, weights)
    if(is.null(slab)) return(NULL)
    return(prior_spike_and_slab(slab, .get_spike_and_slab_inclusion(prior), .prior_model_weight(prior)))
  }
  if(is.prior.mixture(prior)){
    components <- lapply(prior, .bt_regularized_project_component, indices = indices, weights = weights)
    if(any(vapply(components, is.null, logical(1)))) return(NULL)
    for(i in seq_along(components)) components[[i]]$prior_weights <- .prior_model_weight(prior[[i]])
    return(prior_mixture(components, components = attr(prior, "components", exact = TRUE)))
  }
  if(is.prior.none(prior)) return(prior("point", list(0)))
  if(is.prior.point(prior)){
    location <- prior$parameters$location
    if(!is.numeric(location)) return(NULL)
    if(length(location) > 1L) location <- location[indices]
    return(prior("point", list(sum(weights * location))))
  }
  if(!inherits(prior, "prior.regularized_term")){
    if(is.prior.vector(prior) && prior$distribution %in% c("mnormal", "mt"))
      return(.prior_linear_vector_scalar_prior(prior, weights))
    return(NULL)
  }
  spec <- prior$regularization
  units <- unique(spec$unit_index[indices])
  if(length(units) != 1L && spec$family != "ridge") return(NULL)
  .bt_regularized_direction(prior, units[1L], weights)
}

.bt_regularized_centered_plan <- function(prior_list, entries, weights, source_transforms){

  active <- names(entries)[vapply(entries, function(entry) .bt_has_regularized_centered(entry$prior), logical(1))]
  if(!length(active)) return(NULL)
  for(name in active){
    entry <- entries[[name]]
    indices <- match(entry$present, entry$columns)
    w <- weights[entry$present]
    multiplier <- attr(entry$prior, "multiply_by", exact = TRUE)
    fixed_zero <- isTRUE(is.numeric(multiplier) && multiplier == 0) ||
      (is.character(multiplier) && is.prior.point(prior_list[[multiplier]]) &&
       isTRUE(prior_list[[multiplier]]$parameters$location == 0))
    projected <- if(fixed_zero) prior("point", list(0)) else if(!any(!is.na(source_transforms[entry$present]))){
      .bt_regularized_project_component(entry$prior, indices, w)
    }
    if(is.null(projected)){
      atom_free <- inherits(entry$prior, "prior.regularized_term") &&
        !any(!is.na(source_transforms)) && all(vapply(entries, function(e){
          !is.character(attr(e$prior, "multiply_by", exact = TRUE)) && !.bt_prior_has_expression(e$prior)
        }, logical(1)))
      return(list(reason = paste0("The prior density of the weighted direction of '", name,
        "' is unavailable for this combination of shrinkage units or transformations. Use joint prior draws."),
        sources = name, atom_free = atom_free))
    }
    if(!fixed_zero) attr(projected, "multiply_by") <- multiplier
    prior_list[[name]] <- projected
    weights <- c(weights[!names(weights) %in% entry$present], stats::setNames(1, name))
  }
  reduced <- .bt_regularized_linear_plan(prior_list, weights, source_transforms)
  if(is.null(reduced)) return(list(prior_list = prior_list, weights = weights))
  if(is.null(reduced$weights)) reduced$weights <- weights
  reduced
}

# Marginal methods of a centered factor describe its level effects. The
# stored design also covers interactions: different treatment slices may
# have independent locals but the same one-level marginal law.
.bt_regularized_level_prior <- function(prior){

  if(!.bt_has_regularized_centered(prior) || is.prior.mixture(prior)) return(prior)
  design <- attr(prior, "factor_design", exact = TRUE)
  if(!is.matrix(design)) .regularized_numerical_unavailable("unbound factor-level law")
  projections <- lapply(seq_len(nrow(design)), function(i){
    indices <- which(design[i, ] != 0)
    if(!length(indices)) return(NULL)
    law <- .bt_regularized_project_component(prior, indices, design[i, indices])
    if(is.null(law)) .regularized_numerical_unavailable("factor-level law across shrinkage units")
    law
  })
  projections <- Filter(Negate(is.null), projections)
  if(!length(projections)) return(prior("point", list(0)))
  first <- projections[[1L]]
  same <- vapply(projections, function(p){
    identical(p$regularization$group_lasso_dimension, first$regularization$group_lasso_dimension) &&
      abs(p$regularization$multiplier - first$regularization$multiplier) <=
        1e-12 * max(p$regularization$multiplier, first$regularization$multiplier)
  }, logical(1))
  if(!all(same)) .regularized_numerical_unavailable("common marginal across differently scaled factor levels")
  first
}
