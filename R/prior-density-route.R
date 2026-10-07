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
#   scale mixture, including the pure scale mixture of a product without an
#   additive normal term (R/priors-linear-density-combinations.R);
# * "truncated_normal_convolution": a Gaussian convolution whose other term
#   is a truncated normal, in closed form (R/prior-density-truncated-normal.R);
# * "scale_product": the product of a non-normal scalar term and a scalar
#   multiplier (ordered levels with a non-normal total, 'multiply_by'
#   products of non-normal terms);
# * "convolution": the sum of two non-normal simple continuous scalar terms
#   (Cauchy terms first sum to one Cauchy term);
# * "log_scale_product": log(X) + G, one log-source term X and a Gaussian part
#   G, the log image of the scale product X exp(G)
#   (.prior_density_route_log_source_product());
# * "unknown": a combination without a structural route.
# Internal nodes:
# * "mixture": a finite mixture of routes (mixture and spike-and-slab terms,
#   model and conditional mixtures, design rows);
# * "transform": a named monotone output transformation of a route.
# Terms are rewritten before routing: ordered-prior levels as their total and
# share (.prior_density_route_ordered_terms()), and a linear combination of a
# multivariate t vector prior as its univariate t term
# (.prior_density_route_vector_t_terms(), recorded as 'multivariate_t').

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

# 'recipe' (prior list, weights, source transformations, grid size) builds
# the combination's numerical grid, the only grid use (plotted densities, and
# grid heights and probabilities of combinations of simple terms).
.prior_density_route_unknown <- function(reason, provenance, recipe = NULL){
  list(type = "unknown", reason = reason, provenance = provenance, recipe = recipe)
}

.prior_density_route_recipe <- function(prior_list, weights, source_transforms, n_grid){
  list(prior_list = prior_list, weights = weights,
       source_transforms = source_transforms, n_grid = n_grid)
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
       is.prior.vector(group$prior) || is.prior.ordered(group$prior) ||
       .prior_density_route_vector_mixture(group$prior)){
      return(NULL)
    }
    random_group <- group
  }
  list(offset = offset, group = random_group)
}

# Whether a mixture or spike-and-slab prior has a vector component (e.g. the
# mean-difference factor priors of a model-averaged factor term). A coordinate
# of such a term is not a scalar term: the mixture is expanded into its
# components, each routed as a vector prior (normal, the univariate t of a
# multivariate t, or an atom for a point component), and the ordinate is the
# weighted sum of theirs.
.prior_density_route_vector_mixture <- function(prior){

  if(!is.prior.mixture(prior) && !is.prior.spike_and_slab(prior)){
    return(FALSE)
  }
  any(vapply(prior, function(component){
    is.prior.vector(component) || .prior_density_route_vector_mixture(component)
  }, logical(1)))
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
# evaluated by quadrature over that term's declared support (in closed form
# when that term is a truncated normal). NULL otherwise.
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
  if(.prior_truncated_normal_convolution_eligible(spec)){
    return(list(type = "truncated_normal_convolution", spec = spec, n_grid = n_grid))
  }
  list(type = "conditional_normal", spec = spec, n_grid = n_grid)
}

# Ordered-prior levels (and allocation subsets) are the ordered total times an
# allocation share (.prior_ordered_linear_share()): a fixed share scales the
# total, and a Beta(alpha_1, alpha_2) share multiplies it. The ordered terms of
# a combination are rewritten as these terms: the total with the share's
# scale as weight and, for a Beta share, the share as its 'multiply_by' scale.
# NULL without ordered terms; a 'reason' when the ordered term has no such
# representation.
.prior_density_route_ordered_terms <- function(prior_list, weights, source_transforms){

  groups <- tryCatch(
    .prior_linear_weight_groups(prior_list, weights),
    error = function(e) NULL
  )
  if(is.null(groups)){
    return(NULL)
  }
  ordered <- names(groups)[vapply(groups, function(group){
    is.prior.ordered(group$prior)
  }, logical(1))]
  if(length(ordered) == 0L){
    return(NULL)
  }

  for(parameter in ordered){
    group <- groups[[parameter]]
    if(!is.null(attr(group$prior, "multiply_by", exact = TRUE))){
      return(list(reason = "Scaled ordered-prior terms are not structurally classified."))
    }
    share <- tryCatch(
      .prior_ordered_linear_share(group$prior, group$weights, group$indices),
      error = function(e) e
    )
    if(inherits(share, "error")){
      return(list(reason = conditionMessage(share)))
    }
    weights <- weights[setdiff(names(weights), names(group$weights))]
    if(identical(share$type, "point") && share$scale == 0){
      next
    }
    terms <- .prior_ordered_share_terms(
      total      = group$prior$total,
      share      = share,
      total_name = paste0(".ordered_total[", parameter, "]"),
      share_name = paste0(".ordered_share[", parameter, "]")
    )
    prior_list[names(terms$prior_list)] <- terms$prior_list
    weights[[terms$total_name]] <- terms$weight
    source_transforms[[terms$total_name]] <- NA_character_
  }
  list(prior_list = prior_list, weights = weights,
       source_transforms = source_transforms[names(weights)])
}

# Multivariate t terms (the 'mt' vector priors, 'mcauchy' with df = 1, e.g.
# mean-difference and orthonormal factor priors): X = mu 1 + z / sqrt(w) with
# z ~ N(0, S), S = s^2 I the scale matrix, and w ~ Gamma(nu / 2, nu / 2), as
# the JAGS emitter draws them (.JAGS_prior.vector()) and mvtnorm::rmvt() and
# dmvt() with sigma = S evaluate them. A linear combination a'X is therefore
# the univariate t with location mu sum(a), scale sqrt(a' S a) = s ||a|| and
# nu degrees of freedom (zero weights are dropped before, so a' S a > 0; a
# zero combination is the point mu sum(a) = 0). Each such group is rewritten
# as that scalar t term under the prior's own name (keeping its
# 'multiply_by' scale), and the combination is routed as any combination of
# scalar terms; 'provenance' records the rewritten groups. NULL without
# multivariate t groups, and when a group has no representable scalar t
# (non-numeric parameters, a truncation, which vector priors do not support,
# a source transformation, or a scale that underflows or overflows), whose
# combination keeps the general route.
.prior_density_route_vector_t_terms <- function(prior_list, weights, source_transforms){

  groups <- tryCatch(
    .prior_linear_weight_groups(prior_list, weights),
    error = function(e) NULL
  )
  if(is.null(groups)){
    return(NULL)
  }
  vector_t <- names(groups)[vapply(groups, function(group){
    prior <- group$prior
    is.prior.vector(prior) && identical(prior$distribution, "mt") &&
      !is.prior.mixture(prior) && !is.prior.spike_and_slab(prior)
  }, logical(1))]
  if(length(vector_t) == 0L){
    return(NULL)
  }

  provenance <- list()
  for(parameter in vector_t){
    group <- groups[[parameter]]
    parameters <- group$prior$parameters[c("location", "scale", "df")]
    numeric_parameters <- all(vapply(parameters, function(value){
      is.numeric(value) && length(value) == 1L && is.finite(value)
    }, logical(1)))
    full_support <- is.list(group$prior$truncation) &&
      identical(group$prior$truncation$lower, -Inf) &&
      identical(group$prior$truncation$upper, Inf)
    if(!numeric_parameters || !full_support ||
       any(!is.na(source_transforms[names(group$weights)]))){
      return(NULL)
    }
    location <- sum(group$weights) * parameters$location
    scale <- .prior_density_ordinate_stable_norm(group$weights) * parameters$scale
    if(!is.finite(location) || !is.finite(scale) || scale <= 0){
      return(NULL)
    }
    scalar <- prior("t", list(location = location, scale = scale, df = parameters$df))
    attr(scalar, "multiply_by") <- attr(group$prior, "multiply_by", exact = TRUE)
    prior_list[[parameter]] <- scalar
    weights <- weights[setdiff(names(weights), names(group$weights))]
    weights[[parameter]] <- 1
    source_transforms[[parameter]] <- NA_character_
    provenance[[length(provenance) + 1L]] <- list(
      parameter  = parameter,
      family     = "mt",
      weights    = .prior_density_ordinate_compact(group$weights),
      parameters = c(location = parameters$location, scale = parameters$scale,
                     df = parameters$df),
      t          = c(location = location, scale = scale, df = parameters$df)
    )
  }
  list(prior_list = prior_list, weights = weights,
       source_transforms = source_transforms[names(weights)],
       provenance = provenance)
}

# Route of a combination with 'multiply_by' products. One product term
# X = a + b * s (additive part a, multiplied part b, multiplier s):
# * mixture and spike-and-slab terms (also of the multiplier) are expanded
#   into component combinations;
# * a point multiplier k, or a deterministic multiplied part c, is folded
#   into an additive term (k b, or c s; nothing when zero);
# * a normal multiplied part with a normal or deterministic additive part is
#   the conditional-normal route (a pure scale mixture when a is
#   deterministic);
# * a single non-normal scalar multiplied term with a deterministic additive
#   part is a scale product (.prior_scale_product_ordinate()).
# Other products (several products, a non-normal additive term, several
# multiplied non-normal terms, a multiplier that also has its own weight) have
# no structural route.
.prior_density_route_product <- function(prior_list, weights, split,
                                         source_transforms, n_grid){

  general <- function(reason = "General products are not structurally classified."){
    .prior_density_route_unknown(
      reason     = reason,
      provenance = list(
        kind        = "general_product",
        weights     = .prior_density_ordinate_compact(weights),
        multipliers = names(split$product_groups)
      ),
      recipe     = .prior_density_route_recipe(prior_list, weights, source_transforms, n_grid)
    )
  }
  if(length(split$product_groups) != 1L){
    return(general())
  }
  product <- split$product_groups[[1L]]
  multiplier <- product$multiplier
  multiplier_prior <- prior_list[[multiplier]]
  additive_weights <- split$additive_weights[split$additive_weights != 0]
  if(!is.prior(multiplier_prior) || .prior_linear_prior_dimension(multiplier_prior) != 1L ||
     multiplier %in% names(additive_weights) ||
     !is.null(attr(multiplier_prior, "multiply_by", exact = TRUE))){
    return(general())
  }

  active <- unique(c(
    .prior_linear_active_parameters(prior_list, additive_weights),
    names(product$prior_list),
    multiplier
  ))
  mixtures <- active[vapply(prior_list[active], function(prior){
    is.prior.mixture(prior) || is.prior.spike_and_slab(prior)
  }, logical(1))]
  if(length(mixtures) > 0L){
    expansion <- .prior_density_route_mixture_expansion(prior_list, active, n_grid)
    if(is.null(expansion)){
      return(general("The mixture expansion of this product exceeds the leaf cap."))
    }
    return(.prior_density_route_expand(expansion, weights, source_transforms, n_grid))
  }

  # a point multiplier scales the multiplied terms
  if(is.prior.none(multiplier_prior) || is.prior.point(multiplier_prior)){
    k <- if(is.prior.none(multiplier_prior)) 0 else multiplier_prior$parameters$location
    folded_priors <- prior_list
    for(parameter in names(product$prior_list)){
      attr(folded_priors[[parameter]], "multiply_by") <- NULL
    }
    folded_weights <- weights
    folded <- tryCatch(.prior_linear_fold_scale(product$weights, k),
      BayesTools_numerical_condition = function(e) e)
    if(inherits(folded, "BayesTools_numerical_condition")) return(.prior_density_route_unknown(
      reason = folded$reason, provenance = list(kind = "numerical_scale_unavailable",
        numerical_condition = unclass(folded))))
    folded_weights[names(product$weights)] <- folded
    return(.prior_density_route_linear(folded_priors, folded_weights, source_transforms, n_grid))
  }
  if(!.prior_density_simple_continuous(multiplier_prior)){
    return(general())
  }

  # a deterministic multiplied part c adds c * s
  constant <- .prior_density_ordinate_deterministic_offset(
    product$prior_list, product$weights, source_transforms
  )
  if(!is.null(constant)){
    if(!is.finite(constant)){
      return(general())
    }
    folded_weights <- additive_weights
    if(constant != 0){
      folded_weights[[multiplier]] <- constant
    }
    folded_transforms <- source_transforms[names(folded_weights)]
    names(folded_transforms) <- names(folded_weights)
    return(.prior_density_route_linear(prior_list, folded_weights, folded_transforms, n_grid))
  }

  spec <- .prior_conditional_normal_spec(prior_list, split, source_transforms)
  if(!is.null(spec)){
    return(list(type = "conditional_normal", spec = spec, n_grid = n_grid))
  }

  # a single non-normal scalar multiplied term with a deterministic additive part
  offset <- .prior_density_ordinate_deterministic_offset(
    prior_list, additive_weights, source_transforms
  )
  multiplied <- tryCatch(
    .prior_linear_weight_groups(product$prior_list, product$weights),
    error = function(e) NULL
  )
  if(is.null(offset) || !is.finite(offset) || length(multiplied) != 1L ||
     length(multiplied[[1L]]$weights) != 1L ||
     !is.na(source_transforms[names(multiplied[[1L]]$weights)])){
    return(general())
  }
  factor <- multiplied[[1L]]$prior
  attr(factor, "multiply_by") <- NULL
  if(!.prior_density_simple_continuous(factor)){
    return(general())
  }
  list(
    type   = "scale_product",
    spec   = .prior_scale_product_spec(
      offset     = offset,
      scale      = unname(multiplied[[1L]]$weights[[1L]]),
      factor     = factor,
      multiplier = multiplier_prior,
      sources    = list(
        additive   = names(additive_weights),
        multiplied = names(multiplied),
        multiplier = multiplier
      )
    ),
    n_grid = n_grid
  )
}

# Cauchy terms (Cauchy or t with one degree of freedom, untruncated) of a
# product-free combination sum to one Cauchy term: sum_j w_j C(l_j, s_j) =
# C(sum_j w_j l_j, sum_j |w_j| s_j). NULL with fewer than two Cauchy terms.
.prior_density_route_cauchy_terms <- function(prior_list, weights, source_transforms){

  groups <- tryCatch(
    .prior_linear_weight_groups(prior_list, weights),
    error = function(e) NULL
  )
  if(is.null(groups)){
    return(NULL)
  }
  is_cauchy <- function(group){
    prior <- group$prior
    .prior_density_simple_continuous(prior) &&
      (identical(prior$distribution, "cauchy") ||
         (identical(prior$distribution, "t") && identical(prior$parameters$df, 1))) &&
      identical(prior$truncation$lower, -Inf) && identical(prior$truncation$upper, Inf) &&
      all(is.na(source_transforms[names(group$weights)]))
  }
  cauchy <- names(groups)[vapply(groups, is_cauchy, logical(1))]
  columns <- unlist(lapply(groups[cauchy], function(group) names(group$weights)), use.names = FALSE)
  if(length(columns) < 2L){
    return(NULL)
  }
  location <- 0
  scale <- 0
  for(parameter in cauchy){
    group <- groups[[parameter]]
    location <- location + sum(group$weights) * group$prior$parameters$location
    scale <- scale + sum(abs(group$weights)) * group$prior$parameters$scale
  }
  prior_list[[".cauchy_sum"]] <- prior("cauchy", list(location = location, scale = scale))
  weights <- weights[setdiff(names(weights), columns)]
  weights[[".cauchy_sum"]] <- 1
  source_transforms <- source_transforms[names(weights)]
  names(source_transforms) <- names(weights)
  list(prior_list = prior_list, weights = weights, source_transforms = source_transforms)
}

# Two simple continuous scalar terms (not both normal) plus points: the
# two-term convolution quadrature over the first term, the one with an
# infinite density at a finite bound when only one has such a bound. NULL for
# other combinations.
.prior_density_route_convolution <- function(prior_list, weights, source_transforms, n_grid){

  groups <- tryCatch(
    .prior_linear_weight_groups(prior_list, weights),
    error = function(e) NULL
  )
  if(is.null(groups)){
    return(NULL)
  }
  offset <- 0
  terms <- list()
  for(group in groups){
    point <- .prior_density_ordinate_point_group_location(group, source_transforms)
    if(length(point) == 1L){
      if(is.na(point)){
        return(NULL)
      }
      offset <- offset + point
      next
    }
    if(!.prior_density_simple_continuous(group$prior) ||
       any(!is.na(source_transforms[names(group$weights)]))){
      return(NULL)
    }
    for(weight in group$weights){
      terms[[length(terms) + 1L]] <- list(prior = group$prior, weight = unname(weight))
    }
  }
  if(length(terms) != 2L){
    return(NULL)
  }
  singular <- vapply(terms, function(term){
    bounds <- unlist(term$prior$truncation[c("lower", "upper")], use.names = FALSE)
    any(vapply(bounds[is.finite(bounds)], function(bound){
      identical(.prior_density_ordinate_continuous_behavior(
        .prior_density_ordinate_primitive(term$prior, bound)
      ), "infinite")
    }, logical(1)))
  }, logical(1))
  if(identical(singular, c(FALSE, TRUE))){
    terms <- terms[2:1]
  }
  list(
    type   = "convolution",
    spec   = .prior_convolution_spec(offset, terms, sources = names(groups)),
    n_grid = n_grid
  )
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

  ordered <- .prior_density_route_ordered_terms(prior_list, weights, source_transforms)
  if(!is.null(ordered)){
    if(!is.null(ordered$reason)){
      return(.prior_density_route_unknown(
        reason     = ordered$reason,
        provenance = list(kind = "unsupported_provenance", ordered = TRUE),
        recipe     = .prior_density_route_recipe(prior_list, weights, source_transforms, n_grid)
      ))
    }
    return(.prior_density_route_linear(
      ordered$prior_list, ordered$weights, ordered$source_transforms, n_grid
    ))
  }

  vector_t <- .prior_density_route_vector_t_terms(prior_list, weights, source_transforms)
  if(!is.null(vector_t)){
    route <- .prior_density_route_linear(
      vector_t$prior_list, vector_t$weights, vector_t$source_transforms, n_grid
    )
    route$multivariate_t <- vector_t$provenance
    return(route)
  }

  split <- tryCatch(
    .prior_linear_split_multiply_groups(prior_list, weights),
    BayesTools_numerical_condition = function(e) e,
    error = function(e) NULL
  )
  if(inherits(split, "BayesTools_numerical_condition")) return(.prior_density_route_unknown(
    reason = split$reason, provenance = list(kind = "numerical_scale_unavailable",
      numerical_condition = unclass(split))))
  if(is.null(split)){
    return(.prior_density_route_unknown(
      reason     = NULL,
      provenance = list(kind = "unsupported_provenance")
    ))
  }

  if(length(split$product_groups) > 0L){
    return(.prior_density_route_product(
      prior_list, weights, split, source_transforms, n_grid
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
  cauchy <- .prior_density_route_cauchy_terms(prior_list, weights, source_transforms)
  if(!is.null(cauchy)){
    return(.prior_density_route_linear(
      cauchy$prior_list, cauchy$weights, cauchy$source_transforms, n_grid
    ))
  }
  convolution <- .prior_density_route_convolution(prior_list, weights, source_transforms, n_grid)
  if(!is.null(convolution)){
    return(convolution)
  }
  product <- .prior_density_route_log_source_product(
    .prior_density_route_recipe(prior_list, weights, source_transforms, n_grid)
  )
  if(!is.null(product)){
    return(list(type = "log_scale_product", product = product))
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
    ),
    recipe     = .prior_density_route_recipe(prior_list, weights, source_transforms, n_grid)
  )
}

# A named monotone output transformation of 'source'. 'hull' returns the exact
# support hull of the source (NULL when unknown); region probabilities of
# 'exp_lin' transformations need it. The exp of a log-source term plus a
# Gaussian part is a scale product (.prior_density_route_exp_scale_product()).
.prior_density_route_transform <- function(source, transformation, arguments, hull){

  if(is.null(transformation)){
    return(source)
  }
  exp_product <- .prior_density_route_exp_scale_product(source, transformation,
                                                         arguments, hull)
  if(!is.null(exp_product)){
    return(exp_product)
  }
  list(type = "transform", source = source, transformation = transformation,
       arguments = arguments, hull = hull)
}

# The exp of X + G, with X one term through a log source (weight 1; a positive
# simple continuous prior, e.g. the unscaled intercept of a log-intercept
# formula scaling) and G ~ N(m, s), s > 0, the Gaussian part (normal terms and
# points), is the scale product Y = X W of X and the independent lognormal
# multiplier W = exp(G) ~ lognormal(m, s): the 'scale_product' leaf, with its
# exact ordinates, region probabilities, plotted values and the
# classification of its offset at 0. The log-scale sum X + G is routed as the
# log image of that product ('log_scale_product'), whose exp is the product
# itself. Mixture and spike-and-slab priors of X are expanded before (the exp
# of each component is rewritten on its own; a component that cannot be keeps
# the transformation node). NULL for other transformations and sources:
# log-source terms with another weight (X^w W, which the leaf's multiplier map
# does not represent), several log-source terms, a non-Gaussian other term, a
# Gaussian part without variance (the scalar route), and products.
.prior_density_route_exp_scale_product <- function(source, transformation,
                                                   arguments, hull){

  if(!identical(transformation, "exp") || length(arguments) > 0L){
    return(NULL)
  }
  if(identical(source$type, "log_scale_product")){
    return(source$product)
  }
  if(!identical(source$type, "mixture") || !is.null(source$rows)){
    return(NULL)
  }
  products <- lapply(source$components, .prior_density_route_exp_scale_product,
                     transformation = transformation, arguments = arguments,
                     hull = hull)
  rewritten <- !vapply(products, is.null, logical(1))
  if(!any(rewritten)){
    return(NULL)
  }
  # the exp of a mixture is the mixture of the components' exps
  source$components <- Map(function(component, product){
    if(!is.null(product)){
      return(product)
    }
    list(type = "transform", source = component, transformation = transformation,
         arguments = arguments, hull = hull)
  }, source$components, products)
  source
}

# The scale-product leaf of the exp of the combination in 'recipe' (prior
# list, weights, source transformations and grid size; see
# .prior_density_route_exp_scale_product()), or NULL.
.prior_density_route_log_source_product <- function(recipe){

  weights <- recipe$weights[recipe$weights != 0]
  source_transforms <- recipe$source_transforms
  if(length(weights) < 2L || is.null(source_transforms)){
    return(NULL)
  }
  source_transforms <- source_transforms[names(weights)]
  log_source <- !is.na(source_transforms) & source_transforms == "log"
  if(sum(log_source) != 1L || any(!is.na(source_transforms) & !log_source) ||
     unname(weights[log_source]) != 1){
    return(NULL)
  }
  groups <- tryCatch(
    .prior_linear_weight_groups(recipe$prior_list, weights),
    error = function(e) NULL
  )
  if(is.null(groups) || any(vapply(groups, function(group){
    !is.null(attr(group$prior, "multiply_by", exact = TRUE))
  }, logical(1)))){
    return(NULL)
  }
  log_column <- names(weights)[log_source]
  log_group <- groups[vapply(groups, function(group){
    log_column %in% names(group$weights)
  }, logical(1))]
  if(length(log_group) != 1L || length(log_group[[1L]]$weights) != 1L){
    return(NULL)
  }
  factor <- log_group[[1L]]$prior
  if(!.prior_density_simple_continuous(factor) || factor$truncation$lower < 0){
    return(NULL)
  }
  gaussian_weights <- weights[!log_source]
  gaussian <- .prior_density_ordinate_linear_normal(
    recipe$prior_list, gaussian_weights, source_transforms[names(gaussian_weights)], 0
  )
  if(is.null(gaussian) || !identical(gaussian$method, "linear_normal") ||
     !isTRUE(gaussian$exact) || !is.finite(gaussian$provenance$mean) ||
     !is.finite(gaussian$provenance$sd) || gaussian$provenance$sd <= 0){
    return(NULL)
  }
  multiplier <- prior("lognormal", list(meanlog = gaussian$provenance$mean,
                                         sdlog   = gaussian$provenance$sd))
  n_grid <- if(is.null(recipe$n_grid)) .prior_linear_density_default_grid() else recipe$n_grid
  list(
    type   = "scale_product",
    spec   = .prior_scale_product_spec(
      offset     = 0,
      scale      = 1,
      factor     = factor,
      multiplier = multiplier,
      sources    = list(
        additive   = character(),
        multiplied = names(log_group),
        multiplier = setdiff(names(groups), names(log_group))
      )
    ),
    n_grid = n_grid
  )
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
    standardized <- .prior_density_context_standardized_weights(context, weights, source_transforms)
    if(!is.null(source_transforms)){
      source_transforms <- source_transforms[names(standardized)]
    }
    canonical <- .bt_formula_context_canonical(context, standardized)
    offset <- attr(standardized, "formula_recipe_offset", exact = TRUE)
    attr(standardized, "formula_recipe_offset") <- NULL
    attr(standardized, "formula_contribution_transforms") <- NULL
    .bt_formula_require_multiplier_laws(canonical, standardized)
    route <- .prior_density_route_linear_arguments(list(
      prior_list                      = canonical$prior_list,
      n_grid                          = context$n_grid,
      weights                         = standardized,
      source_transforms               = source_transforms,
      output_transformation           = if(is.null(offset)) transformation else NULL,
      output_transformation_arguments = if(is.null(offset)) transformation_arguments else NULL
    ))
    if(!is.null(offset)){
      support <- .bt_formula_route_support(route)
      route <- .prior_density_route_transform(route, "lin", list(a = offset, b = 1),
        hull = function() if(is.null(support)) NULL else support$bounds)
      shifted_support <- .bt_formula_route_support(route)
      route <- .prior_density_route_transform(route, transformation, transformation_arguments,
        hull = function() if(is.null(shifted_support)) NULL else shifted_support$bounds)
    }
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
  if(identical(adaptive$kind, "density_mixture")){
    if(!is.list(arguments$dists) || !is.numeric(arguments$weights) ||
       length(arguments$dists) != length(arguments$weights) || any(!is.finite(arguments$weights)) ||
       any(arguments$weights < 0) || sum(arguments$weights) <= 0){
      stop("Stored leaf prior mixture metadata are malformed.", call. = FALSE)
    }
    components <- lapply(arguments$dists, function(dist){
      .prior_density_route_from_adaptive(attr(dist, "adaptive_evaluation", exact = TRUE))
    })
    if(any(vapply(components, is.null, logical(1)))){
      return(.prior_density_route_unknown("A leaf prior law has no supported recipe.", list(kind = "density_mixture")))
    }
    return(.prior_density_route_mixture(components, arguments$weights / sum(arguments$weights),
      provenance_extra = list(context = "model_specific_formula_laws")))
  }
  if(identical(adaptive$kind, "linear_combination")){
    return(.prior_density_route_linear_arguments(arguments))
  }
  if(identical(adaptive$kind, "allocation_product")){
    # a scale prior times a variance-allocation multiplier
    # (R/prior-density-allocation.R)
    return(.prior_density_route_allocation_product(arguments))
  }
  if(identical(adaptive$kind, "output_transformation")){
    # a monotone transformation of a recorded density whose builder cannot
    # take it (.prior_density_output_transform())
    source <- .prior_density_route_from_adaptive(arguments$source)
    if(is.null(source)){
      return(NULL)
    }
    return(.prior_density_route_transform(
      source         = source,
      transformation = arguments$output_transformation,
      arguments      = arguments$output_transformation_arguments,
      hull           = function(){
        .prior_linear_density_support_hull(arguments$source)
      }
    ))
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
    "conditional_normal" = .prior_numerical_ordinate(
      .prior_conditional_normal_ordinate(route$spec, value, route$n_grid), value,
      "conditional_normal_mixture", .prior_density_route_provenance(route)),
    "truncated_normal_convolution" = .prior_numerical_ordinate(
      .prior_truncated_normal_convolution_ordinate(route$spec, value), value,
      "truncated_normal_convolution", .prior_density_route_provenance(route)),
    "scale_product" = .prior_numerical_ordinate(
      .prior_scale_product_ordinate(route$spec, value, route$n_grid), value,
      "scale_mixture", .prior_density_route_provenance(route)),
    "log_scale_product" = .prior_numerical_ordinate(
      .prior_density_route_log_scale_product_ordinate(route, value), value,
      "scale_mixture", .prior_density_route_provenance(route)),
    "convolution" = .prior_numerical_ordinate(
      .prior_convolution_ordinate(route$spec, value, route$n_grid), value,
      "convolution", .prior_density_route_provenance(route)),
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
  if(!is.null(route$multivariate_t)){
    result$provenance$multivariate_t <- route$multivariate_t
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

# Whether the values 'x' are representable at full double precision: finite
# and not subnormal (|x| >= .Machine$double.xmin; zero is excluded as well).
# Below .Machine$double.xmin a result is rounded to a multiple of the smallest
# subnormal, so a density evaluated at (or through) it is shifted wherever it
# varies near zero, and a result that underflows to zero or overflows is not
# the argument at all. No ordinate is exact when its value is computed from
# such an intermediate: every route that evaluates a density at an argument
# derived from the requested value applies this rule to that argument (see
# .prior_density_affine_full_precision() and the audit in
# .agents/instructions/priors.md).
.prior_density_full_precision <- function(x){
  is.finite(x) & abs(x) >= .Machine$double.xmin
}

# The full-precision rule for an argument derived from the values 'x' by its
# distance from an exact structural point 'anchor' (the offset of an affine
# map, a product or a convolution; 0 for a primitive density) and that
# distance divided by each nonzero element of 'scales': TRUE where x equals
# the anchor (the route classifies that point structurally) or where the
# distance and all its standardizations are representable at full precision.
.prior_density_affine_full_precision <- function(x, anchor = 0, scales = 1){

  difference <- x - anchor
  out <- .prior_density_full_precision(difference)
  for(scale in scales[is.finite(scales) & scales != 0]){
    out <- out & .prior_density_full_precision(difference / scale)
  }
  (!is.na(difference) & difference == 0) | out
}

# Ordinate of the log image Z = log(Y) of a scale product Y at 'value':
# f_Z(z) = f_Y(e^z) e^z, with the product's classification at e^z > 0 (never
# its offset at 0) and its quadrature error scaled by e^z. A value whose
# exponential is not representable at full precision (zero, infinite, or
# subnormal: below .Machine$double.xmin, e^z is rounded, which moves f_Y(e^z)
# when f_Y varies near 0) has no ordinate value.
.prior_density_route_log_scale_product_ordinate <- function(route, value){

  y <- exp(value)
  provenance <- list(kind = "log_scale_product")
  if(!.prior_density_full_precision(y)){
    return(.prior_density_ordinate_imprecise(
      value, "The exponential of the requested value", "scale_mixture", provenance
    ))
  }
  source <- .prior_scale_product_ordinate(route$product$spec, y, route$product$n_grid)
  provenance$source <- source$provenance
  .prior_density_ordinate_wrap(
    source       = source,
    value        = value,
    log_jacobian = -value,
    method       = "scale_mixture",
    provenance   = provenance
  )
}

# Structural provenance of a route (the classification provenance without
# value-dependent entries), used by named transformations to classify support
# bounds, atoms and boundary limits of the source.
.prior_density_route_provenance <- function(route){

  provenance <- switch(
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
    "truncated_normal_convolution" = .prior_truncated_normal_convolution_provenance(route$spec),
    "scale_product" = .prior_scale_product_route_provenance(route$spec),
    "log_scale_product" = list(
      kind   = "log_scale_product",
      source = .prior_scale_product_route_provenance(route$product$spec)
    ),
    "convolution" = .prior_convolution_provenance(route$spec),
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
    "transform" = list(
      kind           = "named_transform",
      transformation = route$transformation,
      arguments      = .prior_density_ordinate_compact(
        .prior_density_ordinate_transform_arguments(route$transformation, route$arguments)
      ),
      source         = .prior_density_route_provenance(route$source)
    )
  )
  if(!is.null(route$multivariate_t)){
    provenance$multivariate_t <- route$multivariate_t
  }
  provenance
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
    "truncated_normal_convolution" = .prior_region_conditional_normal(route$spec, region, route$n_grid),
    "scale_product" = .prior_region_scale_product(route$spec, region, route$n_grid),
    "log_scale_product" = .prior_region_scale_product(
      route$product$spec,
      list(intervals = exp(region$intervals),
           indicator = function(values) region$indicator(log(values))),
      route$product$n_grid
    ),
    "convolution" = .prior_region_convolution(route$spec, region, route$n_grid),
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

# ---- plotted densities -------------------------------------------------------

# Continuous density of a route at the values 'x' (point masses excluded):
# closed forms vectorized over 'x', quadrature leaves by one batched
# quadrature over all values (.prior_density_route_quadrature_density(), with
# 'batch_singular' also for leaves with a singular integrand), and routes
# without a structural representation by linear interpolation of their own
# numerical grid, the only grid use. Non-finite and unavailable values are
# NA; a value where the density is infinite is Inf.
.prior_density_route_density <- function(route, x, batch_singular = FALSE){

  if(length(x) == 0L){
    return(numeric())
  }
  switch(
    route$type,
    "atom" = rep(0, length(x)),
    "scalar" = .prior_density_scalar_density(
      route$prior, x, route$offset, route$scale, route$source_transform
    ),
    "normal" = {
      normal <- .prior_density_ordinate_linear_normal(
        route$prior_list, route$weights, route$source_transforms, 0
      )
      if(!identical(normal$method, "linear_normal")){
        rep(NA_real_, length(x))
      }else{
        stats::dnorm(x, normal$provenance$mean, normal$provenance$sd)
      }
    },
    "truncated_normal_convolution" = .prior_truncated_normal_convolution_density(route$spec, x),
    "log_scale_product" = {
      y <- exp(x)
      representable <- .prior_density_full_precision(y)
      out <- rep(NA_real_, length(x))
      if(any(representable)){
        out[representable] <- .prior_density_route_quadrature_density(
          route$product, y[representable]
        ) * y[representable]
      }
      out
    },
    "mixture" = {
      positive <- route$weights > 0
      weights <- route$weights[positive] / sum(route$weights[positive])
      out <- numeric(length(x))
      for(i in seq_along(weights)){
        out <- out + weights[[i]] *
          .prior_density_route_density(route$components[positive][[i]], x, batch_singular)
      }
      out
    },
    "transform" = .prior_density_route_transform_density(route, x, batch_singular),
    "unknown" = .prior_density_route_grid_density(route, x),
    .prior_density_route_quadrature_density(route, x, batch_singular)
  )
}

# Density of a quadrature leaf (conditional-normal mixture, scale product or
# two-term convolution) at the values 'x'. Values outside a bounded support hull
# are zero; the leaf's special values (offsets, support bounds and meeting
# points, which its ordinate classifies structurally), all values of a leaf
# whose integrand has an integrable singularity, and values whose batched
# integral does not meet the acceptance criterion take the ordinate's own
# value; all other values share one batched quadrature of the ordinate's
# integrand over the ordinate's breakpoints
# (.prior_density_quadrature_batch()), within a relative error of 1e-8. With
# 'batch_singular', the values of a conditional-normal or scale-product leaf
# with a singular integrand, which the batched bisection cannot resolve to
# 1e-8, are batched as well, with the acceptance criterion of their ordinates
# (a relative error of 1e-4, .prior_linear_density_refinement_tolerance()),
# unless a term of the leaf has a strong singularity (an exponent below 0.1,
# .prior_density_strong_singularity(), where that error estimate is not
# reliable); the display grids of route products use it
# (.prior_linear_density_route_product()).
.prior_density_route_quadrature_density <- function(route, x, batch_singular = FALSE){

  plan <- switch(
    route$type,
    "conditional_normal" = .prior_conditional_normal_density_plan(route$spec, x, batch_singular),
    "scale_product"      = .prior_scale_product_density_plan(route$spec, x, batch_singular),
    "convolution"        = .prior_convolution_density_plan(route$spec, x)
  )
  zero <- if(is.null(plan$zero)) rep(FALSE, length(x)) else plan$zero
  out <- rep(NA_real_, length(x))
  out[zero] <- 0
  if(!is.null(plan$integrand)){
    tolerance <- .prior_density_quadrature_tolerance()
    if(!isTRUE(plan$batch)){
      tolerance$relative <- .prior_linear_density_refinement_tolerance()$relative
    }
    out[!zero & !plan$special] <- .prior_density_quadrature_batch(
      plan$integrand, plan$breakpoints, tolerance
    )
  }
  ordinate <- which(!zero & (plan$special | is.na(out)))
  out[ordinate] <- vapply(x[ordinate], function(value){
    result <- .prior_density_route_ordinate(route, value)
    height <- .prior_density_ordinate_height_value(result)
    if(is.na(height) && is.null(result$provenance$inverse_moment)){
      # display only: a quadrature that converged but misses the relative
      # acceptance criterion of ordinates (a far-tail density) is drawn; the
      # quadrature of an offset ordinate is its inverse moment E[1 / |s|],
      # not the density, so such an offset is not drawn
      estimate <- result$provenance$integration$estimate
      if(is.numeric(estimate) && length(estimate) == 1L && is.finite(estimate)){
        height <- estimate
      }
    }
    height
  }, numeric(1))
  out
}

# Whether a route has a leaf of one of the 'types'.
.prior_density_route_has_leaf <- function(route, types){

  if(is.null(route)){
    return(FALSE)
  }
  switch(
    route$type,
    "mixture"   = any(vapply(route$components, .prior_density_route_has_leaf,
                             logical(1), types = types)),
    "transform" = .prior_density_route_has_leaf(route$source, types),
    route$type %in% types
  )
}

.prior_density_route_has_quadrature <- function(route){

  .prior_density_route_has_leaf(route, c("conditional_normal", "scale_product",
                                         "log_scale_product", "convolution"))
}

# Values a plotted density of the route must include: 'points' (atoms,
# offsets of scale mixtures and products, where the density may peak or be
# infinite, normal means, and meeting points of convolution bounds) and the
# finite support bounds of scalar terms and of scale products and
# convolutions, where the density may jump: 'lower' bounds (the support of
# the term lies above them) and 'upper' bounds (below them). The density at
# such a bound is the one-sided limit inside the term's support.
.prior_density_route_display_points <- function(route){

  empty <- list(points = numeric(), lower = numeric(), upper = numeric())
  combine <- function(parts){
    list(points = unlist(lapply(parts, `[[`, "points"), use.names = FALSE),
         lower  = unlist(lapply(parts, `[[`, "lower"), use.names = FALSE),
         upper  = unlist(lapply(parts, `[[`, "upper"), use.names = FALSE))
  }
  finite <- function(values) values[is.finite(values)]
  hull <- function(points, hull){
    list(points = finite(points), lower = finite(hull[1L]), upper = finite(hull[2L]))
  }
  switch(
    route$type,
    "atom" = list(points = finite(route$locations), lower = numeric(), upper = numeric()),
    "scalar" = {
      bounds <- .prior_density_route_prior_bounds(route$prior)
      if(identical(route$source_transform, "log")){
        bounds <- lapply(bounds, function(values) log(values[values > 0]))
      }
      bounds <- lapply(bounds, function(values) finite(route$offset + route$scale * values))
      # a negative scale maps lower bounds of the term to upper bounds
      if(route$scale < 0){
        bounds <- list(lower = bounds$upper, upper = bounds$lower)
      }
      list(points = numeric(), lower = bounds$lower, upper = bounds$upper)
    },
    "normal" = {
      normal <- .prior_density_ordinate_linear_normal(
        route$prior_list, route$weights, route$source_transforms, 0
      )
      list(points = finite(normal$provenance$mean), lower = numeric(), upper = numeric())
    },
    "conditional_normal" = list(points = finite(route$spec$additive_mean),
                                lower = numeric(), upper = numeric()),
    "scale_product" = hull(route$spec$offset, .prior_scale_product_hull(route$spec)),
    "convolution" = hull(.prior_convolution_meeting_points(route$spec),
                         .prior_convolution_hull(route$spec)),
    "mixture" = combine(lapply(route$components[route$weights > 0],
                               .prior_density_route_display_points)),
    "transform" = {
      source <- .prior_density_route_display_points(route$source)
      arguments <- .prior_density_ordinate_transform_arguments(
        route$transformation, route$arguments
      )
      if(!is.character(route$transformation) || is.null(arguments)){
        empty
      }else{
        map <- function(values){
          finite(suppressWarnings(.density.prior_transformation_x(
            values, route$transformation, arguments
          )))
        }
        # 'lin' and 'exp_lin' with a negative 'b' are decreasing and map
        # lower bounds to upper bounds ('tanh' and 'exp' are increasing)
        decreasing <- route$transformation %in% c("lin", "exp_lin") && arguments$b < 0
        list(points = map(source$points),
             lower  = map(if(decreasing) source$upper else source$lower),
             upper  = map(if(decreasing) source$lower else source$upper))
      }
    },
    empty
  )
}

# Finite lower and upper truncation bounds of a simple prior or of the
# components of a mixture or spike-and-slab prior (point components are
# atoms, not bounds).
.prior_density_route_prior_bounds <- function(prior){

  if(is.prior.spike_and_slab(prior) || is.prior.mixture(prior)){
    bounds <- lapply(prior, .prior_density_route_prior_bounds)
    return(list(
      lower = unlist(lapply(bounds, `[[`, "lower"), use.names = FALSE),
      upper = unlist(lapply(bounds, `[[`, "upper"), use.names = FALSE)
    ))
  }
  empty <- list(lower = numeric(), upper = numeric())
  if(!is.prior.simple(prior) || is.prior.point(prior) || is.prior.discrete(prior)){
    return(empty)
  }
  bounds <- prior$truncation[c("lower", "upper")]
  if(!all(vapply(bounds, function(bound) is.numeric(bound) && length(bound) == 1L,
                 logical(1)))){
    return(empty)
  }
  lapply(bounds, function(bound) bound[is.finite(bound)])
}

# The route with the numerical grids of its leaves without a structural
# representation built once (plotted densities evaluate a route repeatedly).
.prior_density_route_with_grids <- function(route){

  if(is.null(route)){
    return(route)
  }
  switch(
    route$type,
    "unknown" = {
      if(!is.null(route$recipe) && is.null(route$grid)){
        route$grid <- .prior_density_route_recipe_grid(route$recipe)
      }
      route
    },
    "mixture" = {
      route$components <- lapply(route$components, .prior_density_route_with_grids)
      route
    },
    "transform" = {
      route$source <- .prior_density_route_with_grids(route$source)
      route
    },
    route
  )
}

# The continuous density of an ordinate result: its height when structural
# (Inf for an infinite density, 0 for a structural zero), NA otherwise.
.prior_density_ordinate_height_value <- function(ordinate){

  behavior <- .prior_density_ordinate_continuous_behavior(ordinate)
  if(identical(behavior, "infinite")){
    return(Inf)
  }
  if(identical(behavior, "zero")){
    return(0)
  }
  if(identical(behavior, "regular") && !is.na(ordinate$log_density)){
    return(exp(ordinate$log_density))
  }
  NA_real_
}

# Continuous density of offset + scale * S (S = log(T) for a log source) at
# 'x' for a scalar prior, a finite mixture of them, or a point prior (no
# continuous part).
.prior_density_scalar_density <- function(prior, x, offset, scale, source_transform = NULL){

  if(scale == 0 || is.prior.none(prior) || is.prior.point(prior) ||
     is.prior.discrete(prior)){
    return(rep(0, length(x)))
  }
  if(is.prior.mixture(prior) || is.prior.spike_and_slab(prior)){
    weights <- .prior_density_ordinate_mixture_weights(prior)
    if(is.null(weights)){
      return(rep(NA_real_, length(x)))
    }
    out <- numeric(length(x))
    for(i in which(weights > 0)){
      out <- out + weights[[i]] *
        .prior_density_scalar_density(prior[[i]], x, offset, scale, source_transform)
    }
    return(out)
  }
  if(!is.prior.simple(prior)){
    return(rep(NA_real_, length(x)))
  }
  source <- (x - offset) / scale
  if(identical(source_transform, "log")){
    original <- exp(source)
    out <- mpdf(prior, original) * original / abs(scale)
    out[original == 0] <- 0
    return(out)
  }
  mpdf(prior, source) / abs(scale)
}

# Density of a named monotone output transformation y = g(s) at 'x':
# f_s(g^-1(y)) / |g'(g^-1(y))|, zero outside the transformation's image.
.prior_density_route_transform_density <- function(route, x, batch_singular = FALSE){

  arguments <- .prior_density_ordinate_transform_arguments(
    route$transformation, route$arguments
  )
  if(!is.character(route$transformation) || is.null(arguments)){
    return(rep(NA_real_, length(x)))
  }
  if(route$transformation %in% c("lin", "exp_lin") && arguments$b == 0){
    return(rep(0, length(x)))
  }
  source <- suppressWarnings(.density.prior_transformation_inv_grid(
    x, route$transformation, arguments
  ))
  out <- rep(0, length(x))
  inside <- is.finite(source)
  if(any(inside)){
    source_density <- .prior_density_route_density(route$source, source[inside], batch_singular)
    out[inside] <- .density.prior_transformation_y(
      x[inside], source_density, route$transformation, arguments
    )
  }
  # the image of a zero source under exp_lin is a boundary whose density is
  # the ordinate's structural limit
  if(identical(route$transformation, "exp_lin") && any(x == 0)){
    out[x == 0] <- .prior_density_ordinate_height_value(
      .prior_density_route_ordinate(route, 0)
    )
  }
  out
}

# A route without a structural representation: linear interpolation of its
# own numerical grid (display only; no refinement), built from its recipe
# unless .prior_density_route_with_grids() stored it.
.prior_density_route_grid_density <- function(route, x){

  if(is.null(route$recipe)){
    return(rep(NA_real_, length(x)))
  }
  grid <- if(is.null(route$grid)) .prior_density_route_recipe_grid(route$recipe) else route$grid
  if(identical(grid, "unavailable")){
    return(rep(NA_real_, length(x)))
  }
  if(is.null(grid$density)){
    return(rep(0, length(x)))
  }
  stats::approx(grid$density$x, grid$density$y, xout = x, yleft = 0, yright = 0)$y *
    grid$density$mass
}

# Whether a leaf without a structural representation relies on a numerical
# grid with an unresolved product component (the grid's
# 'product_grid_resolution' record; .prior_linear_density_route_product()).
# Grids that .prior_density_route_with_grids() did not store are built.
.prior_density_route_unresolved_products <- function(route){

  if(is.null(route)){
    return(FALSE)
  }
  switch(
    route$type,
    "unknown" = {
      if(is.null(route$recipe)){
        return(FALSE)
      }
      grid <- if(is.null(route$grid)) .prior_density_route_recipe_grid(route$recipe) else route$grid
      is.list(grid) && isFALSE(attr(grid, "product_grid_resolution", exact = TRUE)$resolved)
    },
    "mixture"   = any(vapply(route$components, .prior_density_route_unresolved_products, logical(1))),
    "transform" = .prior_density_route_unresolved_products(route$source),
    FALSE
  )
}

# The numerical grid of a route recipe ("unavailable" when it cannot be built).
.prior_density_route_recipe_grid <- function(recipe){

  tryCatch(
    .prior_linear_combination_density(
      prior_list        = recipe$prior_list,
      weights           = recipe$weights,
      n_grid            = recipe$n_grid,
      source_transforms = recipe$source_transforms
    ),
    error = function(e) "unavailable"
  )
}

# ---- numerical grids ---------------------------------------------------------

# Whether a route contains a product without a structural route.
.prior_density_route_has_general_product <- function(route){

  if(is.null(route)){
    return(FALSE)
  }
  switch(
    route$type,
    "unknown"   = identical(route$provenance$kind, "general_product"),
    "mixture"   = any(vapply(route$components, .prior_density_route_has_general_product, logical(1))),
    "transform" = .prior_density_route_has_general_product(route$source),
    FALSE
  )
}

# Numerical grids stand in only for combinations of simple terms without a
# structural route: a product grid is capped in size and cannot be refined
# reliably, so a product without a structural route is unavailable for
# inference.
.prior_linear_density_check_grid <- function(route, quantity = NULL){

  if(.prior_density_route_has_general_product(route)){
    message <- paste0(
      "The prior density of this linear combination is unavailable: its ",
      "'multiply_by' product has no structural density route (a non-normal ",
      "additive term, several products, or several non-normal multiplied ",
      "terms), and numerical product grids are not used for inference. ",
      "Evaluate the terms separately."
    )
    if(identical(quantity, "probability")){
      stop(errorCondition(message, call = NULL,
        class = c("BayesTools_prior_region_route_unavailable",
                  "BayesTools_hypothesis_region")))
    }
    stop(message, call. = FALSE)
  }
  invisible(TRUE)
}
