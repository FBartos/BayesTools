# Prior region probabilities from deterministic prior-density provenance.
#
# A region is the continuous set of disjoint open intervals 'intervals' (a
# two-column matrix, possibly with infinite bounds) together with an exact
# 'indicator' of the region at given values, which decides point masses
# (including strict and inclusive boundaries). Its prior probability is
# evaluated with the structure that classifies prior-density ordinates
# (R/prior-density-ordinate.R), mirroring its routes:
# * point masses contribute their exact mass when the indicator includes them;
# * scalar (affine, optionally log-source) priors use their exact distribution
#   function, normal sums their normal distribution function;
# * Gaussian convolutions and conditional-normal scale mixtures use the
#   conditional-normal quadrature of the ordinate
#   (.prior_conditional_normal_region());
# * mixture, spike-and-slab, model and conditional mixture components, and
#   distinct design rows are expanded as for the ordinate and summed with
#   their probabilities;
# * named monotone output transformations map the region to the source scale.
# A combination for which the ordinate has no such representation (general
# convolutions and products, unsupported families, transformations or
# contexts) is unavailable here, and the caller keeps its grid evaluation.

.prior_region_result <- function(probability, integration = NULL){

  list(
    available      = TRUE,
    probability    = probability,
    absolute_error = if(is.null(integration)) 0 else integration$absolute_error,
    converged      = if(is.null(integration)) TRUE else isTRUE(integration$converged),
    messages       = if(is.null(integration) || isTRUE(integration$converged)){
      character()
    }else if(identical(integration$message, "OK")){
      "an absolute error above its error bound"
    }else{
      integration$message
    },
    evaluations    = if(is.null(integration)) 0 else integration$evaluations,
    quadratures    = if(is.null(integration)) 0L else 1L
  )
}

.prior_region_unavailable <- function(){
  list(available = FALSE)
}

.prior_region_combine <- function(results, weights){

  positive <- weights > 0
  results <- results[positive]
  weights <- weights[positive] / sum(weights[positive])
  if(length(results) == 0L ||
     !all(vapply(results, function(result) isTRUE(result$available), logical(1)))){
    return(.prior_region_unavailable())
  }
  list(
    available      = TRUE,
    probability    = sum(weights * vapply(results, `[[`, numeric(1), "probability")),
    absolute_error = sum(weights * vapply(results, `[[`, numeric(1), "absolute_error")),
    converged      = all(vapply(results, `[[`, logical(1), "converged")),
    messages       = unique(unlist(lapply(results, `[[`, "messages"), use.names = FALSE)),
    evaluations    = sum(vapply(results, `[[`, numeric(1), "evaluations")),
    quadratures    = sum(vapply(results, `[[`, numeric(1), "quadratures"))
  )
}

# Interval algebra on the continuous part of a region (closedness of the
# bounds is irrelevant there; point masses use the exact indicator).
.prior_region_intervals <- function(lower = numeric(), upper = numeric()){

  intervals <- cbind(as.numeric(lower), as.numeric(upper))
  intervals <- intervals[intervals[, 1L] < intervals[, 2L], , drop = FALSE]
  if(nrow(intervals) <= 1L){
    return(intervals)
  }
  intervals <- intervals[order(intervals[, 1L]), , drop = FALSE]
  out <- intervals[1L, , drop = FALSE]
  for(i in seq.int(2L, nrow(intervals))){
    last <- nrow(out)
    if(intervals[i, 1L] <= out[last, 2L]){
      out[last, 2L] <- max(out[last, 2L], intervals[i, 2L])
    }else{
      out <- rbind(out, intervals[i, , drop = FALSE])
    }
  }
  out
}

.prior_region_intervals_complement <- function(intervals){

  bounds <- c(-Inf, as.vector(t(intervals)), Inf)
  odd <- seq(1L, length(bounds), by = 2L)
  .prior_region_intervals(bounds[odd], bounds[odd + 1L])
}

.prior_region_intervals_union <- function(a, b){

  combined <- rbind(a, b)
  .prior_region_intervals(combined[, 1L], combined[, 2L])
}

.prior_region_intervals_intersect <- function(a, b){

  .prior_region_intervals_complement(.prior_region_intervals_union(
    .prior_region_intervals_complement(a),
    .prior_region_intervals_complement(b)
  ))
}

.prior_region_whole <- function(region){

  nrow(region$intervals) == 1L &&
    identical(unname(region$intervals[1L, ]), c(-Inf, Inf))
}

# Exact probability of point masses at 'locations'. A location for which the
# indicator is undefined (NA) makes the region unavailable.
.prior_region_atoms <- function(locations, probability, region){

  keep <- is.finite(locations) & is.finite(probability) & probability > 0
  locations <- locations[keep]
  probability <- probability[keep]
  if(length(locations) == 0L){
    return(.prior_region_result(0))
  }
  inside <- region$indicator(locations)
  if(anyNA(inside)){
    return(.prior_region_unavailable())
  }
  .prior_region_result(sum(probability[inside]))
}

# Probability of the source intervals under a simple continuous prior, from
# its exact (truncation-aware) distribution function; upper-tail
# probabilities are used above the median.
.prior_region_prior_mass <- function(prior, intervals){

  probability <- 0
  for(i in seq_len(nrow(intervals))){
    lower <- intervals[i, 1L]
    upper <- intervals[i, 2L]
    lower_cdf <- if(is.finite(lower)) cdf(prior, lower) else if(lower > 0) 1 else 0
    if(lower_cdf > .5){
      upper_ccdf <- if(is.finite(upper)) ccdf(prior, upper) else if(upper > 0) 0 else 1
      mass <- ccdf(prior, lower) - upper_ccdf
    }else{
      upper_cdf <- if(is.finite(upper)) cdf(prior, upper) else if(upper > 0) 1 else 0
      mass <- upper_cdf - lower_cdf
    }
    probability <- probability + max(0, as.numeric(mass))
  }
  probability
}

# Mirrors .prior_density_ordinate_prior_affine(): offset + scale * S with S
# from 'prior' (S = log(T) for a log source transformation).
.prior_region_prior_affine <- function(prior, region, offset, scale,
                                       source_transform = NULL){

  if(!is.finite(offset) || !is.finite(scale)){
    return(.prior_region_unavailable())
  }
  if(scale == 0){
    return(.prior_region_atoms(offset, 1, region))
  }
  if(is.prior.spike_and_slab(prior) || is.prior.mixture(prior)){
    weights <- .prior_density_ordinate_mixture_weights(prior)
    if(is.null(weights)){
      return(.prior_region_unavailable())
    }
    components <- lapply(prior, function(component){
      .prior_region_prior_affine(component, region, offset, scale, source_transform)
    })
    return(.prior_region_combine(components, weights))
  }
  if(is.prior.none(prior)){
    return(.prior_region_atoms(offset, 1, region))
  }
  if(is.prior.simple(prior) && prior$distribution %in% c("point", "bernoulli")){
    if(identical(prior$distribution, "point")){
      locations <- prior$parameters$location
      probability <- 1
    }else{
      discrete <- tryCatch(.prior_simple_truncated_discrete(prior),
                           error = function(e) NULL)
      if(is.null(discrete)){
        return(.prior_region_unavailable())
      }
      locations <- discrete$support
      probability <- discrete$prob
    }
    if(identical(source_transform, "log")){
      if(any(locations <= 0)){
        return(.prior_region_unavailable())
      }
      locations <- log(locations)
    }else if(!is.null(source_transform)){
      return(.prior_region_unavailable())
    }
    return(.prior_region_atoms(offset + scale * locations, probability, region))
  }
  if(!.prior_region_primitive_supported(prior)){
    return(.prior_region_unavailable())
  }
  if(!is.null(source_transform) &&
     (!identical(source_transform, "log") || prior$truncation$lower < 0)){
    return(.prior_region_unavailable())
  }
  if(.prior_region_whole(region)){
    return(.prior_region_result(1))
  }

  intervals <- (region$intervals - offset) / scale
  if(scale < 0){
    intervals <- intervals[, 2:1, drop = FALSE]
  }
  if(identical(source_transform, "log")){
    intervals <- exp(intervals)
  }
  intervals <- .prior_region_intervals(intervals[, 1L], intervals[, 2L])
  .prior_region_result(.prior_region_prior_mass(prior, intervals))
}

# Families and parameters whose ordinate .prior_density_ordinate_primitive()
# classifies, with a valid numeric truncation.
.prior_region_primitive_supported <- function(prior){

  supported <- c(
    "normal", "lognormal", "t", "gamma", "invgamma", "beta", "exp",
    "uniform", "moment", "invmoment"
  )
  truncation <- prior$truncation
  is.prior.simple(prior) && !is.prior.discrete(prior) && !is.prior.point(prior) &&
    is.character(prior$distribution) && length(prior$distribution) == 1L &&
    prior$distribution %in% supported &&
    .prior_density_ordinate_parameters_numeric(prior) &&
    is.list(truncation) && all(c("lower", "upper") %in% names(truncation)) &&
    all(vapply(truncation[c("lower", "upper")], function(bound){
      is.numeric(bound) && length(bound) == 1L && !is.na(bound)
    }, logical(1)))
}

# Mirrors .prior_density_ordinate_linear_scalar(): NULL when the combination
# is not one scalar term plus point terms.
.prior_region_linear_scalar <- function(prior_list, weights, source_transforms,
                                        region){

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
  if(is.null(random_group)){
    return(.prior_region_atoms(offset, 1, region))
  }
  parameter <- names(random_group$weights)[1L]
  .prior_region_prior_affine(
    prior            = random_group$prior,
    region           = region,
    offset           = offset,
    scale            = unname(random_group$weights[[1L]]),
    source_transform = .prior_linear_source_transform(source_transforms[parameter])
  )
}

# Mirrors .prior_density_ordinate_linear_normal(): NULL when a term is not a
# normal (or log-transformed lognormal) or point term.
.prior_region_linear_normal <- function(prior_list, weights, source_transforms,
                                        region){

  normal <- .prior_density_ordinate_linear_normal(
    prior_list, weights, source_transforms, 0
  )
  if(is.null(normal)){
    return(NULL)
  }
  if(!identical(normal$method, "linear_normal") || !isTRUE(normal$exact)){
    return(.prior_region_unavailable())
  }
  probability <- sum(vapply(seq_len(nrow(region$intervals)), function(i){
    .prior_normal_interval_probability(
      region$intervals[i, 1L], region$intervals[i, 2L],
      normal$provenance$mean, normal$provenance$sd
    )
  }, numeric(1)))
  .prior_region_result(probability)
}

.prior_region_conditional_normal <- function(spec, region, n_grid){

  if(nrow(region$intervals) == 0L){
    return(.prior_region_result(0))
  }
  if(.prior_region_whole(region)){
    return(.prior_region_result(1))
  }
  integral <- .prior_conditional_normal_region(spec, region$intervals, n_grid)
  .prior_region_result(integral$value, integral$integration)
}

# Mirrors .prior_density_ordinate_additive_factor().
.prior_region_additive_factor <- function(prior_list, weights, source_transforms,
                                          region){

  weights <- weights[weights != 0]
  scalar <- .prior_region_linear_scalar(prior_list, weights, source_transforms, region)
  if(!is.null(scalar)){
    return(scalar)
  }
  .prior_region_linear_normal(prior_list, weights, source_transforms, region)
}

# Mirrors .prior_density_ordinate_additive_components(): mixture expansion,
# then a Gaussian convolution; NULL when neither applies.
.prior_region_additive_components <- function(prior_list, weights,
                                              source_transforms, region,
                                              n_grid){

  weights <- weights[weights != 0]
  if(length(weights) == 0L){
    return(NULL)
  }
  groups <- tryCatch(
    .prior_linear_weight_groups(prior_list, weights),
    error = function(e) NULL
  )
  if(!is.null(groups)){
    plan <- .prior_density_ordinate_mixture_plan(prior_list, names(groups), n_grid)
    if(!is.null(plan)){
      components <- lapply(plan$prior_lists, function(component_priors){
        .prior_region_linear_base(
          component_priors, weights, source_transforms, region, n_grid
        )
      })
      combined <- .prior_region_combine(components, plan$probabilities)
      if(isTRUE(combined$available)){
        return(combined)
      }
    }
  }
  spec <- .prior_density_ordinate_gaussian_convolution_spec(
    prior_list, weights, source_transforms
  )
  if(is.null(spec)){
    return(NULL)
  }
  .prior_region_conditional_normal(spec, region, n_grid)
}

# Mirrors .prior_density_ordinate_linear_base() on the linear-predictor scale.
.prior_region_linear_base <- function(prior_list, weights, source_transforms,
                                      region,
                                      n_grid = .prior_linear_density_default_grid()){

  weights <- weights[weights != 0]
  if(length(weights) == 0L){
    return(.prior_region_atoms(0, 1, region))
  }
  if(is.null(source_transforms)){
    source_transforms <- rep(NA_character_, length(weights))
    names(source_transforms) <- names(weights)
  }else{
    source_transforms <- source_transforms[names(weights)]
  }
  if(any(!is.na(source_transforms) & source_transforms != "log")){
    return(.prior_region_unavailable())
  }
  split <- tryCatch(
    .prior_linear_split_multiply_groups(prior_list, weights),
    error = function(e) NULL
  )
  if(is.null(split)){
    return(.prior_region_unavailable())
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
          parameter <- expansion$parameter
          parent <- prior_list[[parameter]]
          probabilities <- .prior_density_ordinate_mixture_weights(parent)
          components <- lapply(expansion$indices, function(i){
            component_priors <- prior_list
            component_priors[[parameter]] <- .prior_density_copy_parent_attributes(parent[[i]], parent)
            .prior_region_linear_base(
              component_priors, weights, source_transforms, region, n_grid
            )
          })
          combined <- .prior_region_combine(components, probabilities[expansion$indices])
          if(isTRUE(combined$available)){
            return(combined)
          }
        }
      }
      product_constant <- .prior_density_ordinate_deterministic_offset(
        product$prior_list, product$weights, source_transforms
      )
      if(identical(product_constant, 0)){
        additive <- .prior_region_additive_factor(
          prior_list, split$additive_weights, source_transforms, region
        )
        if(is.null(additive)){
          additive <- .prior_region_additive_components(
            prior_list, split$additive_weights, source_transforms, region, n_grid
          )
        }
        if(!is.null(additive)){
          return(additive)
        }
      }
    }
    spec <- .prior_conditional_normal_spec(prior_list, split, source_transforms)
    if(!is.null(spec)){
      return(.prior_region_conditional_normal(spec, region, n_grid))
    }
    # general products (including the product singularity of the ordinate)
    return(.prior_region_unavailable())
  }

  weights <- split$additive_weights
  weights <- weights[weights != 0]
  scalar <- .prior_region_linear_scalar(prior_list, weights, source_transforms, region)
  if(!is.null(scalar)){
    return(scalar)
  }
  normal <- .prior_region_linear_normal(prior_list, weights, source_transforms, region)
  if(!is.null(normal)){
    return(normal)
  }
  components <- .prior_region_additive_components(
    prior_list, weights, source_transforms, region, n_grid
  )
  if(!is.null(components)){
    return(components)
  }
  .prior_region_unavailable()
}

# Maps an output-scale region through the inverse of a named monotone
# transformation to the source scale; 'evaluate' computes the source-scale
# probability. A constant transformation is a point mass. 'hull' returns the
# source support hull, required by 'exp_lin' (defined on a nonnegative
# source). Unsupported transformations are unavailable.
.prior_region_transformed <- function(transformation, arguments, region,
                                      evaluate, hull){

  if(is.null(transformation)){
    return(evaluate(region))
  }
  if(!is.character(transformation) || length(transformation) != 1L ||
     !transformation %in% c("lin", "exp", "exp_lin", "tanh")){
    return(.prior_region_unavailable())
  }
  arguments <- .prior_density_ordinate_transform_arguments(transformation, arguments)
  if(is.null(arguments)){
    return(.prior_region_unavailable())
  }
  if(transformation %in% c("lin", "exp_lin") && arguments$b == 0){
    location <- if(identical(transformation, "lin")) arguments$a else exp(arguments$a)
    if(!is.finite(location) || location == 0 && identical(transformation, "exp_lin")){
      return(.prior_region_unavailable())
    }
    return(.prior_region_atoms(location, 1, region))
  }
  if(identical(transformation, "exp_lin")){
    support <- hull()
    if(is.null(support) || anyNA(support) || support[1L] < 0){
      return(.prior_region_unavailable())
    }
  }

  lower <- region$intervals[, 1L]
  upper <- region$intervals[, 2L]
  source <- switch(
    transformation,
    "lin" = {
      mapped <- cbind((lower - arguments$a) / arguments$b,
                      (upper - arguments$a) / arguments$b)
      if(arguments$b < 0) mapped[, 2:1, drop = FALSE] else mapped
    },
    "exp" = cbind(
      ifelse(lower <= 0, -Inf, log(pmax(lower, 0))),
      ifelse(upper <= 0, -Inf, log(pmax(upper, 0)))
    ),
    "tanh" = cbind(
      atanh(pmin(pmax(lower, -1), 1)),
      atanh(pmin(pmax(upper, -1), 1))
    ),
    "exp_lin" = {
      inverse <- function(y){
        out <- exp((log(pmax(y, 0)) - arguments$a) / arguments$b)
        out[y <= 0] <- if(arguments$b > 0) 0 else Inf
        out
      }
      mapped <- cbind(inverse(lower), inverse(upper))
      if(arguments$b < 0) mapped[, 2:1, drop = FALSE] else mapped
    }
  )
  if(anyNA(source)){
    return(.prior_region_unavailable())
  }
  source_region <- list(
    intervals = .prior_region_intervals(source[, 1L], source[, 2L]),
    indicator = local({
      indicator <- region$indicator
      function(values){
        output <- suppressWarnings(.density.prior_transformation_x(
          values, transformation, arguments
        ))
        if(identical(transformation, "exp_lin")){
          output[values <= 0] <- NA_real_
        }
        out <- rep(NA, length(values))
        defined <- is.finite(output)
        out[defined] <- indicator(output[defined])
        out
      }
    })
  )
  evaluate(source_region)
}

# Mirrors .prior_density_ordinate_linear_arguments().
.prior_region_linear_arguments <- function(arguments, region){

  prior_list <- arguments$prior_list
  weights <- arguments$weights
  source_transforms <- arguments$source_transforms
  if(!is.list(prior_list) || !is.numeric(weights) || is.null(names(weights)) ||
     anyNA(weights) || any(!is.finite(weights))){
    return(.prior_region_unavailable())
  }
  n_grid <- if(is.null(arguments$n_grid)) .prior_linear_density_default_grid() else arguments$n_grid
  .prior_region_transformed(
    arguments$output_transformation,
    arguments$output_transformation_arguments,
    region,
    evaluate = function(source_region){
      .prior_region_linear_base(prior_list, weights, source_transforms,
                                source_region, n_grid)
    },
    hull = function(){
      .prior_linear_combination_support_hull(prior_list, weights, source_transforms)
    }
  )
}

# Mirrors .prior_density_ordinate_context_classifier(): each model or
# conditional component is evaluated with the full budget and weighted by its
# probability.
.prior_region_context <- function(context, weights, source_transforms,
                                  transformation, transformation_arguments,
                                  region){

  if(inherits(context, "prior_density_context")){
    standardized <- tryCatch(
      .prior_density_context_standardized_weights(context, weights),
      error = function(e) NULL
    )
    if(is.null(standardized)){
      return(.prior_region_unavailable())
    }
    if(!is.null(source_transforms)){
      source_transforms <- source_transforms[names(standardized)]
    }
    return(.prior_region_linear_arguments(list(
      prior_list                      = context$prior_list,
      n_grid                          = context$n_grid,
      weights                         = standardized,
      source_transforms               = source_transforms,
      output_transformation           = transformation,
      output_transformation_arguments = transformation_arguments
    ), region))
  }

  hull <- function(){
    .prior_linear_context_support_hull(context, weights, source_transforms)
  }
  if(inherits(context, "prior_density_model_mixture_context")){
    evaluate <- function(source_region){
      component_indices <- which(context$model_weights > 0)
      results <- lapply(component_indices, function(model_i){
        model_prior_list <- lapply(context$prior_list, function(parameter_priors){
          if(is.prior(parameter_priors)) parameter_priors else parameter_priors[[model_i]]
        })
        for(parameter in names(model_prior_list)){
          if(is.null(model_prior_list[[parameter]])){
            model_prior_list[[parameter]] <- prior("point", list(location = 0))
          }
        }
        .prior_region_linear_base(
          model_prior_list, weights, source_transforms, source_region,
          n_grid = context$n_grid
        )
      })
      .prior_region_combine(results, context$model_weights[component_indices])
    }
    return(.prior_region_transformed(transformation, transformation_arguments,
                                     region, evaluate, hull))
  }

  if(inherits(context, "prior_density_conditional_context")){
    evaluate <- function(source_region){
      component_indices <- which(context$model_weights > 0)
      results <- lapply(context$prior_lists[component_indices], function(prior_list){
        if(!is.null(context$formula_scale) && length(context$formula_scale) > 0L){
          component_context <- .prior_density_context(
            prior_list,
            context$column_names,
            context$formula_scale,
            context$n_grid,
            context$tail_prob
          )
          return(.prior_region_context(
            component_context, weights, source_transforms, NULL, NULL,
            source_region
          ))
        }
        .prior_region_linear_base(
          prior_list, weights, source_transforms, source_region,
          n_grid = context$n_grid
        )
      })
      .prior_region_combine(results, context$model_weights[component_indices])
    }
    return(.prior_region_transformed(transformation, transformation_arguments,
                                     region, evaluate, hull))
  }

  .prior_region_unavailable()
}

# Mirrors .prior_density_ordinate_from_adaptive(); NULL without recorded
# provenance.
.prior_region_from_adaptive <- function(adaptive, region){

  if(!is.list(adaptive) || !is.character(adaptive$kind) ||
     length(adaptive$kind) != 1L || !is.list(adaptive$arguments)){
    return(NULL)
  }
  arguments <- adaptive$arguments
  if(identical(adaptive$kind, "linear_combination")){
    return(.prior_region_linear_arguments(arguments, region))
  }
  if(identical(adaptive$kind, "density_context")){
    return(.prior_region_context(
      arguments$context,
      arguments$weights,
      arguments$source_transforms,
      arguments$output_transformation,
      arguments$output_transformation_arguments,
      region
    ))
  }
  if(identical(adaptive$kind, "density_context_rows")){
    weights <- arguments$weights
    if(!is.numeric(weights) || anyNA(weights) || any(!is.finite(weights))){
      return(.prior_region_unavailable())
    }
    if(is.null(dim(weights))){
      weights <- matrix(weights, nrow = 1L, dimnames = list(NULL, names(weights)))
    }
    weights <- as.matrix(weights)
    if(nrow(weights) == 0L){
      return(.prior_region_atoms(0, 1, region))
    }
    row_keys <- apply(weights, 1L, function(row){
      paste(sprintf("%a", row), collapse = "\r")
    })
    unique_keys <- unique(row_keys)
    row_counts <- tabulate(match(row_keys, unique_keys), nbins = length(unique_keys))
    row_indices <- match(unique_keys, row_keys)
    # each distinct row is an independent probability with the full budget
    results <- lapply(row_indices, function(row_i){
      .prior_region_context(
        arguments$context,
        weights[row_i, ],
        arguments$source_transforms,
        arguments$output_transformation,
        arguments$output_transformation_arguments,
        region
      )
    })
    return(.prior_region_combine(results, row_counts))
  }
  NULL
}

# Prior probability of 'region' from the deterministic provenance of a
# prior_linear_density; NULL when the density has no such representation
# (the caller then uses its grid). A conditional-normal quadrature that fails
# its acceptance criterion stops.
.prior_linear_density_region_probability <- function(x, region){

  result <- .prior_region_from_adaptive(
    attr(x, "adaptive_evaluation", exact = TRUE),
    region
  )
  if(is.null(result) || !isTRUE(result$available)){
    return(NULL)
  }
  if(!isTRUE(result$converged) || !is.finite(result$probability)){
    messages <- result$messages
    if(length(messages) == 0L){
      messages <- "a non-finite probability"
    }
    stop(
      "Conditional-normal prior probability was rejected by diagnostics: ",
      "integration reported '", paste(messages, collapse = "; "),
      "' with absolute error ", format(result$absolute_error),
      ". Inspect the prior specification and the region bounds.",
      call. = FALSE
    )
  }
  probability <- result$probability
  attr(probability, "numerical_diagnostics") <- list(
    method         = if(result$quadratures > 0) "conditional_normal_quadrature" else "exact",
    quadratures    = result$quadratures,
    absolute_error = result$absolute_error,
    evaluations    = result$evaluations
  )
  probability
}
