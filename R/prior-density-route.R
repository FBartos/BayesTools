# Structural routes of prior densities.
#
# A route is the value-independent structure of a deterministic prior density
# (a linear combination of prior terms, optionally standardized by a density
# context, mixed over models, conditional components or design rows, and
# mapped by a named monotone output transformation). The prior-density
# ordinate (R/prior-density-ordinate.R) and the prior region probabilities
# (R/prior-density-region.R) evaluate the same route, so a combination is
# classified in one place. Routes are built from the prior definitions only
# (no numerical evaluation) and are transient: they are never stored on a
# density or in a provenance record.
#
# Leaves:
# * "atom": point masses (locations, probabilities);
# * "scalar": one scalar prior term (possibly a finite mixture) times a weight,
#   plus a deterministic offset, optionally through a log source;
# * "normal": a sum of normal terms (normal, mnormal, log-lognormal) and
#   points;
# * "conditional_normal": a Gaussian convolution or a conditional-normal
#   scale mixture (R/priors-linear-density-combinations.R);
# * "product": a general product term without a structural route; its value
#   at the deterministic offset is classified by the product singularity;
# * "unknown": a combination without a structural route.
# Internal nodes:
# * "mixture": a finite mixture of routes (mixture and spike-and-slab terms,
#   model and conditional mixtures, design rows);
# * "transform": a named monotone output transformation of a route.

# Mixture and spike-and-slab terms are expanded into their component
# combinations only while the number of resulting leaves is at most the number
# of initial quadrature rules the evaluation budget admits (the smallest rule
# has 15 evaluations), which bounds the combinatorial expansion before any
# numerical evaluation.
.prior_density_route_leaf_cap <- function(n_grid){

  floor(n_grid / 15)
}

# The prior list of model 'model_i' of a model-mixture prior list; a model
# without a term has the point prior at zero.
.prior_density_model_prior_list <- function(prior_list, model_i){

  model_prior_list <- lapply(prior_list, function(parameter_priors){
    if(is.prior(parameter_priors)) parameter_priors else parameter_priors[[model_i]]
  })
  names(model_prior_list) <- names(prior_list)
  for(parameter in names(model_prior_list)){
    if(is.null(model_prior_list[[parameter]])){
      model_prior_list[[parameter]] <- prior("point", list(location = 0))
    }
  }
  model_prior_list
}

# Bitwise-distinct rows of a weight matrix: the first row index and the count
# of each distinct row.
.prior_density_distinct_rows <- function(weights){

  row_keys <- apply(weights, 1L, function(row){
    paste(sprintf("%a", row), collapse = "\r")
  })
  unique_keys <- unique(row_keys)
  list(
    indices = match(unique_keys, row_keys),
    counts  = tabulate(match(row_keys, unique_keys), nbins = length(unique_keys)),
    n       = length(unique_keys)
  )
}

.prior_density_route_atom <- function(locations, probability, provenance){
  list(type = "atom", locations = locations, probability = probability,
       provenance = provenance)
}

.prior_density_route_unknown <- function(reason, provenance){
  list(type = "unknown", reason = reason, provenance = provenance)
}

.prior_density_route_mixture <- function(components, weights,
                                         provenance_extra = list()){
  list(type = "mixture", components = components, weights = weights,
       provenance_extra = provenance_extra)
}

# Number of positive-probability leaves of a (possibly nested) mixture or
# spike-and-slab prior; NULL when its weights are not exact numeric values.
.prior_density_route_prior_leaves <- function(prior){

  if(!is.prior.mixture(prior) && !is.prior.spike_and_slab(prior)){
    return(1)
  }
  probabilities <- .prior_density_ordinate_mixture_weights(prior)
  if(is.null(probabilities)){
    return(NULL)
  }
  count <- 0
  for(i in which(probabilities > 0)){
    leaves <- .prior_density_route_prior_leaves(prior[[i]])
    if(is.null(leaves)){
      return(NULL)
    }
    count <- count + leaves
  }
  count
}

# Expansion of the first mixture or spike-and-slab prior among 'parameters'
# into its positive-probability components (each a prior list with that
# component in place of the mixture); the other mixtures are expanded when the
# components are routed in turn. NULL when there is no expandable mixture or
# the expansion exceeds the leaf cap.
.prior_density_route_mixture_expansion <- function(prior_list, parameters, n_grid){

  parameters <- unique(parameters[parameters %in% names(prior_list)])
  mixtures <- parameters[vapply(prior_list[parameters], function(prior){
    is.prior.mixture(prior) || is.prior.spike_and_slab(prior)
  }, logical(1))]
  if(length(mixtures) == 0L){
    return(NULL)
  }
  leaves <- lapply(prior_list[mixtures], .prior_density_route_prior_leaves)
  if(any(vapply(leaves, is.null, logical(1))) ||
     prod(unlist(leaves)) > .prior_density_route_leaf_cap(n_grid)){
    return(NULL)
  }

  parameter <- mixtures[[1L]]
  parent <- prior_list[[parameter]]
  probabilities <- .prior_density_ordinate_mixture_weights(parent)
  indices <- which(probabilities > 0)
  list(
    parameter     = parameter,
    probabilities = probabilities[indices],
    prior_lists   = lapply(indices, function(i){
      component_priors <- prior_list
      component_priors[[parameter]] <- .prior_density_copy_parent_attributes(parent[[i]], parent)
      component_priors
    })
  )
}

.prior_density_route_expand <- function(expansion, weights, source_transforms, n_grid){

  .prior_density_route_mixture(
    components = lapply(expansion$prior_lists, function(component_priors){
      .prior_density_route_linear(component_priors, weights, source_transforms, n_grid)
    }),
    weights = expansion$probabilities
  )
}

# One scalar term plus point terms: the deterministic offset and the scalar
# group; NULL for other combinations and NA for an undefined point term.
.prior_density_route_scalar_parts <- function(prior_list, weights, source_transforms){

  groups <- tryCatch(
    .prior_linear_weight_groups(prior_list, weights),
    error = function(e) NULL
  )
  if(is.null(groups)){
    return(NULL)
  }

  offset <- 0
  random_group <- NULL
  for(group in groups){
    point_location <- .prior_density_ordinate_point_group_location(
      group,
      source_transforms
    )
    if(length(point_location) == 1L){
      if(is.na(point_location)){
        return(NULL)
      }
      offset <- offset + point_location
      next
    }
    if(length(group$weights) != 1L || !is.null(random_group) ||
       is.prior.vector(group$prior) || is.prior.ordered(group$prior)){
      return(NULL)
    }
    random_group <- group
  }
  list(offset = offset, group = random_group)
}

.prior_density_route_scalar <- function(parts, weights, source_transforms,
                                        record_weights){

  if(is.null(parts$group)){
    return(.prior_density_route_atom(
      locations   = parts$offset,
      probability = 1,
      provenance  = list(
        kind    = "scalar_affine",
        offset  = parts$offset,
        scale   = 0,
        weights = .prior_density_ordinate_compact(weights)
      )
    ))
  }
  parameter <- names(parts$group$weights)[1L]
  list(
    type             = "scalar",
    prior            = parts$group$prior,
    offset           = parts$offset,
    scale            = unname(parts$group$weights[[1L]]),
    source_transform = .prior_linear_source_transform(source_transforms[parameter]),
    weights          = if(isTRUE(record_weights)) .prior_density_ordinate_compact(weights)
  )
}

.prior_density_route_normal <- function(prior_list, weights, source_transforms){

  normal <- .prior_density_ordinate_linear_normal(
    prior_list, weights, source_transforms, 0
  )
  if(is.null(normal)){
    return(NULL)
  }
  list(type = "normal", prior_list = prior_list, weights = weights,
       source_transforms = source_transforms)
}

# Scalar or normal route of a product-free combination (NULL otherwise).
.prior_density_route_additive_factor <- function(prior_list, weights,
                                                 source_transforms,
                                                 record_weights = FALSE){

  weights <- weights[weights != 0]
  parts <- .prior_density_route_scalar_parts(prior_list, weights, source_transforms)
  if(!is.null(parts)){
    return(.prior_density_route_scalar(parts, weights, source_transforms, record_weights))
  }
  .prior_density_route_normal(prior_list, weights, source_transforms)
}

# Product-free combinations that neither the scalar nor the normal route
# covers: mixture and spike-and-slab terms are expanded into their component
# combinations, each routed on its own, so a density jump of one component
# (e.g. a truncated prior at its bound, using the one-sided limit inside its
# support) never reaches a numerical grid; a Gaussian term plus one other
# continuous scalar term is a positive-variance Gaussian convolution,
# evaluated by quadrature over that term's declared support. NULL otherwise.
.prior_density_route_additive_components <- function(prior_list, weights,
                                                     source_transforms, n_grid){

  weights <- weights[weights != 0]
  if(length(weights) == 0L){
    return(NULL)
  }
  groups <- tryCatch(
    .prior_linear_weight_groups(prior_list, weights),
    error = function(e) NULL
  )
  if(!is.null(groups)){
    expansion <- .prior_density_route_mixture_expansion(prior_list, names(groups), n_grid)
    if(!is.null(expansion)){
      return(.prior_density_route_expand(expansion, weights, source_transforms, n_grid))
    }
  }
  spec <- .prior_density_ordinate_gaussian_convolution_spec(
    prior_list, weights, source_transforms
  )
  if(is.null(spec)){
    return(NULL)
  }
  list(type = "conditional_normal", spec = spec, n_grid = n_grid)
}

# Route of the linear combination sum_j weights[j] * term_j of 'prior_list'
# (terms entering through a 'multiply_by' scale are products).
.prior_density_route_linear <- function(prior_list, weights, source_transforms,
                                        n_grid = .prior_linear_density_default_grid()){

  weights <- weights[weights != 0]
  if(length(weights) == 0L){
    return(.prior_density_route_atom(
      locations   = 0,
      probability = 1,
      provenance  = list(kind = "scalar_affine", offset = 0, scale = 0)
    ))
  }
  if(is.null(source_transforms)){
    source_transforms <- rep(NA_character_, length(weights))
    names(source_transforms) <- names(weights)
  }else{
    source_transforms <- source_transforms[names(weights)]
  }
  unsupported_transform <- !is.na(source_transforms) & source_transforms != "log"
  if(any(unsupported_transform)){
    return(.prior_density_route_unknown(
      reason     = "The source transformation is not supported structurally.",
      provenance = list(
        kind              = "unsupported_provenance",
        weights           = .prior_density_ordinate_compact(weights),
        source_transforms = .prior_density_ordinate_compact(source_transforms)
      )
    ))
  }

  split <- tryCatch(
    .prior_linear_split_multiply_groups(prior_list, weights),
    error = function(e) NULL
  )
  if(is.null(split)){
    return(.prior_density_route_unknown(
      reason     = NULL,
      provenance = list(kind = "unsupported_provenance")
    ))
  }

  if(length(split$product_groups) > 0L){
    if(length(split$product_groups) == 1L){
      product <- split$product_groups[[1L]]
      active <- unique(c(
        .prior_linear_active_parameters(prior_list, split$additive_weights),
        names(product$prior_list)
      ))
      mixtures <- active[vapply(prior_list[active], function(prior){
        is.prior.mixture(prior) || is.prior.spike_and_slab(prior)
      }, logical(1))]
      if(length(mixtures) > 0L){
        expansion <- .prior_conditional_normal_expansion(prior_list, split, source_transforms, n_grid)
        if(!is.null(expansion)){
          parent <- prior_list[[expansion$parameter]]
          probabilities <- .prior_density_ordinate_mixture_weights(parent)
          return(.prior_density_route_mixture(
            components = lapply(expansion$indices, function(i){
              component_priors <- prior_list
              component_priors[[expansion$parameter]] <-
                .prior_density_copy_parent_attributes(parent[[i]], parent)
              .prior_density_route_linear(component_priors, weights, source_transforms, n_grid)
            }),
            weights = probabilities[expansion$indices]
          ))
        }
      }
      product_constant <- .prior_density_ordinate_deterministic_offset(
        product$prior_list, product$weights, source_transforms
      )
      if(identical(product_constant, 0)){
        additive <- .prior_density_route_additive_factor(
          prior_list, split$additive_weights, source_transforms
        )
        if(is.null(additive)){
          additive <- .prior_density_route_additive_components(
            prior_list, split$additive_weights, source_transforms, n_grid
          )
        }
        if(!is.null(additive)){
          return(additive)
        }
      }
    }
    spec <- .prior_conditional_normal_spec(prior_list, split, source_transforms)
    if(!is.null(spec)){
      return(list(type = "conditional_normal", spec = spec, n_grid = n_grid))
    }
    return(list(
      type              = "product",
      prior_list        = prior_list,
      split             = split,
      source_transforms = source_transforms,
      weights           = weights
    ))
  }

  weights <- split$additive_weights
  weights <- weights[weights != 0]

  additive <- .prior_density_route_additive_factor(
    prior_list, weights, source_transforms, record_weights = TRUE
  )
  if(!is.null(additive)){
    return(additive)
  }
  components <- .prior_density_route_additive_components(
    prior_list, weights, source_transforms, n_grid
  )
  if(!is.null(components)){
    return(components)
  }

  .prior_density_route_unknown(
    reason     = "General numerical convolutions are not structurally classified.",
    provenance = list(
      kind              = "general_convolution",
      weights           = .prior_density_ordinate_compact(weights),
      source_families   = vapply(prior_list, function(prior){
        if(is.prior(prior) && !is.null(prior$distribution)) prior$distribution else "unknown"
      }, character(1)),
      source_transforms = .prior_density_ordinate_compact(source_transforms)
    )
  )
}

# A named monotone output transformation of 'source'. 'hull' returns the exact
# support hull of the source (NULL when unknown); region probabilities of
# 'exp_lin' transformations need it.
.prior_density_route_transform <- function(source, transformation, arguments, hull){

  if(is.null(transformation)){
    return(source)
  }
  list(type = "transform", source = source, transformation = transformation,
       arguments = arguments, hull = hull)
}

# Route of the linear combination recorded in the arguments of a
# 'linear_combination' density.
.prior_density_route_linear_arguments <- function(arguments){

  prior_list <- arguments$prior_list
  weights <- arguments$weights
  source_transforms <- arguments$source_transforms
  if(!is.list(prior_list) || !is.numeric(weights) || is.null(names(weights)) ||
     anyNA(weights) || any(!is.finite(weights))){
    return(.prior_density_route_unknown(
      reason     = NULL,
      provenance = list(kind = "unsupported_provenance")
    ))
  }
  n_grid <- if(is.null(arguments$n_grid)) .prior_linear_density_default_grid() else arguments$n_grid
  .prior_density_route_transform(
    source         = .prior_density_route_linear(prior_list, weights, source_transforms, n_grid),
    transformation = arguments$output_transformation,
    arguments      = arguments$output_transformation_arguments,
    hull           = function(){
      .prior_linear_combination_support_hull(prior_list, weights, source_transforms)
    }
  )
}

# Route of the combination 'weights' in a density context: standardized
# coefficients, model mixtures, or conditional mixtures; each model or
# conditional component is an independent route.
.prior_density_route_context <- function(context, weights, source_transforms,
                                         transformation, transformation_arguments){

  if(inherits(context, "prior_density_context")){
    standardized <- tryCatch(
      .prior_density_context_standardized_weights(context, weights),
      error = function(e) NULL
    )
    if(is.null(standardized)){
      return(.prior_density_route_unknown(
        reason     = NULL,
        provenance = list(kind = "density_context")
      ))
    }
    if(!is.null(source_transforms)){
      source_transforms <- source_transforms[names(standardized)]
    }
    route <- .prior_density_route_linear_arguments(list(
      prior_list                      = context$prior_list,
      n_grid                          = context$n_grid,
      weights                         = standardized,
      source_transforms               = source_transforms,
      output_transformation           = transformation,
      output_transformation_arguments = transformation_arguments
    ))
    route$context <- list(
      kind                 = "prior_density_context",
      requested_weights    = .prior_density_ordinate_compact(weights),
      standardized_weights = .prior_density_ordinate_compact(standardized)
    )
    return(route)
  }

  hull <- function(){
    .prior_linear_context_support_hull(context, weights, source_transforms)
  }
  if(inherits(context, "prior_density_model_mixture_context")){
    indices <- which(context$model_weights > 0)
    route <- .prior_density_route_mixture(
      components = lapply(indices, function(model_i){
        .prior_density_route_linear(
          .prior_density_model_prior_list(context$prior_list, model_i),
          weights, source_transforms, n_grid = context$n_grid
        )
      }),
      weights          = context$model_weights[indices],
      provenance_extra = list(context = "model_mixture")
    )
    return(.prior_density_route_transform(route, transformation, transformation_arguments, hull))
  }

  if(inherits(context, "prior_density_conditional_context")){
    indices <- which(context$model_weights > 0)
    route <- .prior_density_route_mixture(
      components = lapply(context$prior_lists[indices], function(prior_list){
        if(!is.null(context$formula_scale) && length(context$formula_scale) > 0L){
          component_context <- .prior_density_context(
            prior_list,
            context$column_names,
            context$formula_scale,
            context$n_grid,
            context$tail_prob
          )
          return(.prior_density_route_context(
            component_context, weights, source_transforms, NULL, NULL
          ))
        }
        .prior_density_route_linear(
          prior_list, weights, source_transforms, n_grid = context$n_grid
        )
      }),
      weights          = context$model_weights[indices],
      provenance_extra = list(context = "conditional_mixture")
    )
    return(.prior_density_route_transform(route, transformation, transformation_arguments, hull))
  }

  .prior_density_route_unknown(
    reason     = NULL,
    provenance = list(kind = "unknown_density_context")
  )
}

# Route of a recorded deterministic prior density ('adaptive_evaluation'
# attribute); NULL without such a record.
.prior_density_route_from_adaptive <- function(adaptive){

  if(!is.list(adaptive) || !is.character(adaptive$kind) ||
     length(adaptive$kind) != 1L || !is.list(adaptive$arguments)){
    return(NULL)
  }
  arguments <- adaptive$arguments
  if(identical(adaptive$kind, "linear_combination")){
    return(.prior_density_route_linear_arguments(arguments))
  }
  if(identical(adaptive$kind, "density_context")){
    return(.prior_density_route_context(
      arguments$context,
      arguments$weights,
      arguments$source_transforms,
      arguments$output_transformation,
      arguments$output_transformation_arguments
    ))
  }
  if(identical(adaptive$kind, "density_context_rows")){
    weights <- arguments$weights
    if(!is.numeric(weights) || anyNA(weights) || any(!is.finite(weights))){
      return(.prior_density_route_unknown(
        reason     = "The row-varying density weights are malformed.",
        provenance = list(kind = "density_context_rows")
      ))
    }
    if(is.null(dim(weights))){
      weights <- matrix(weights, nrow = 1L, dimnames = list(NULL, names(weights)))
    }
    weights <- as.matrix(weights)
    if(nrow(weights) == 0L){
      return(.prior_density_route_atom(
        locations   = 0,
        probability = 1,
        provenance  = list(kind = "density_context_rows", rows = 0L)
      ))
    }
    rows <- .prior_density_distinct_rows(weights)
    # Each distinct row is an independent route with the full evaluation
    # budget and its own convergence diagnostics.
    route <- .prior_density_route_mixture(
      components = lapply(rows$indices, function(row_i){
        .prior_density_route_context(
          arguments$context,
          weights[row_i, ],
          arguments$source_transforms,
          arguments$output_transformation,
          arguments$output_transformation_arguments
        )
      }),
      weights = rows$counts
    )
    route$rows <- list(
      weight_dimensions = as.integer(dim(weights)),
      unique_rows       = rows$n,
      total_rows        = nrow(weights),
      weights_hash      = .prior_density_ordinate_numeric_hash(weights)
    )
    return(route)
  }
  NULL
}

# ---- ordinates --------------------------------------------------------------

.prior_density_route_ordinate <- function(route, value){

  result <- switch(
    route$type,
    "atom" = .prior_density_ordinate_atom_result(
      value       = value,
      locations   = route$locations,
      probability = route$probability,
      method      = if(identical(route$provenance$kind, "density_context_rows")) "finite_mixture" else "scalar_affine",
      provenance  = route$provenance
    ),
    "scalar" = {
      scalar <- .prior_density_ordinate_prior_affine(
        prior            = route$prior,
        value            = value,
        offset           = route$offset,
        scale            = route$scale,
        source_transform = route$source_transform
      )
      if(!is.null(route$weights)){
        scalar$provenance$weights <- route$weights
      }
      scalar
    },
    "normal" = .prior_density_ordinate_linear_normal(
      route$prior_list, route$weights, route$source_transforms, value
    ),
    "conditional_normal" = .prior_conditional_normal_ordinate(route$spec, value, route$n_grid),
    "product" = {
      singularity <- .prior_density_ordinate_product_singularity(
        route$prior_list, route$split, route$source_transforms, value
      )
      if(is.null(singularity)){
        singularity <- .prior_density_ordinate_result(
          value       = value,
          behavior    = "unknown",
          log_density = NA_real_,
          exact       = FALSE,
          method      = "unsupported_provenance",
          reason      = "General products are not structurally classified.",
          provenance  = list(
            kind        = "general_product",
            weights     = .prior_density_ordinate_compact(route$weights),
            multipliers = names(route$split$product_groups)
          )
        )
      }
      singularity
    },
    "unknown" = .prior_density_ordinate_result(
      value       = value,
      behavior    = "unknown",
      log_density = NA_real_,
      exact       = FALSE,
      method      = "unsupported_provenance",
      reason      = route$reason,
      provenance  = route$provenance
    ),
    "mixture" = .prior_density_route_mixture_ordinate(route, value),
    "transform" = .prior_density_ordinate_named_transform(
      classifier        = function(source_value){
        .prior_density_route_ordinate(route$source, source_value)
      },
      source_provenance = .prior_density_route_provenance(route$source),
      transformation    = route$transformation,
      arguments         = route$arguments,
      value             = value
    )
  )
  if(!is.null(route$context)){
    result$provenance$context <- route$context
  }
  result
}

.prior_density_route_mixture_ordinate <- function(route, value){

  results <- lapply(route$components, .prior_density_route_ordinate, value = value)
  combined <- .prior_density_ordinate_combine(
    results, route$weights, value,
    provenance_extra = route$provenance_extra
  )
  if(is.null(route$rows)){
    return(combined)
  }

  continuous_behavior <- .prior_density_ordinate_continuous_behavior(combined)
  integration <- .prior_density_ordinate_integration(combined$provenance)
  combined$provenance <- c(
    list(kind = "finite_mixture", context = "density_context_rows"),
    route$rows[c("weight_dimensions", "unique_rows", "total_rows", "weights_hash")],
    list(row_classifications = Map(function(result, count){
      row <- list(
        count               = unname(count),
        behavior            = result$behavior,
        continuous_behavior =
          .prior_density_ordinate_continuous_behavior(result),
        point_mass          = result$point_mass,
        method              = result$method,
        source_kind         = result$provenance$kind
      )
      row_integration <- .prior_density_ordinate_integration(result$provenance)
      if(!is.null(row_integration)) row$integration <- row_integration
      row
    }, results, route$weights))
  )
  if(!is.null(integration)) combined$provenance$integration <- integration
  if(identical(combined$behavior, "point_mass")){
    combined$provenance$continuous_behavior <- continuous_behavior
  }
  combined
}

# Structural provenance of a route (the classification provenance without
# value-dependent entries), used by named transformations to classify support
# bounds, atoms and boundary limits of the source.
.prior_density_route_provenance <- function(route){

  switch(
    route$type,
    "atom" = route$provenance,
    "scalar" = {
      provenance <- list(
        kind             = "scalar_affine",
        offset           = unname(route$offset),
        scale            = unname(route$scale),
        source_transform = route$source_transform,
        source           = if(is.prior(route$prior)){
          .prior_density_ordinate_prior_definition_provenance(route$prior)
        }else{
          list(kind = "unsupported_provenance")
        }
      )
      if(!is.null(route$weights)){
        provenance$weights <- route$weights
      }
      provenance
    },
    "normal" = .prior_density_ordinate_linear_normal(
      route$prior_list, route$weights, route$source_transforms, 0
    )$provenance,
    "conditional_normal" = list(kind = "conditional_normal_mixture"),
    "product" = list(kind = "general_product"),
    "unknown" = route$provenance,
    "mixture" = {
      positive <- route$weights > 0
      weights <- route$weights[positive] / sum(route$weights[positive])
      c(
        list(
          kind       = "finite_mixture",
          weights    = unname(weights),
          components = Map(function(component, weight){
            list(weight = unname(weight),
                 provenance = .prior_density_route_provenance(component))
          }, route$components[positive], weights)
        ),
        route$provenance_extra
      )
    },
    "transform" = list(kind = "named_transform",
                       transformation = route$transformation,
                       source = .prior_density_route_provenance(route$source))
  )
}

# ---- region probabilities ----------------------------------------------------

.prior_density_route_region <- function(route, region){

  switch(
    route$type,
    "atom" = .prior_region_atoms(route$locations, route$probability, region),
    "scalar" = .prior_region_prior_affine(
      prior            = route$prior,
      region           = region,
      offset           = route$offset,
      scale            = route$scale,
      source_transform = route$source_transform
    ),
    "normal" = .prior_region_linear_normal(
      route$prior_list, route$weights, route$source_transforms, region
    ),
    "conditional_normal" = .prior_region_conditional_normal(route$spec, region, route$n_grid),
    "product" = .prior_region_unavailable(),
    "unknown" = .prior_region_unavailable(),
    "mixture" = .prior_region_combine(
      lapply(route$components, .prior_density_route_region, region = region),
      route$weights
    ),
    "transform" = .prior_region_transformed(
      transformation = route$transformation,
      arguments      = route$arguments,
      region         = region,
      evaluate       = function(source_region){
        .prior_density_route_region(route$source, source_region)
      },
      hull           = route$hull
    )
  )
}
