.plot_data_prior_list.PETPEESE       <- function(prior_list, x_seq, x_range, x_range_quant, n_points, n_samples,
                                                 transformation, transformation_arguments, transformation_settings, prior_list_mu,
                                                 effect_direction = "positive", force_samples = FALSE){

  pair <- .model_probability_petpeese_plot_check(prior_list, prior_list_mu)
  prior_list <- .model_probability_plot_priors(prior_list, validated_pair = pair)
  # The x-axis is the standard error; the effect-size transformation applies
  # only to the regression line (y-axis), so 'transformation_settings' does
  # not rescale x.
  if(is.null(x_seq)){
    x_seq <- seq(x_range[1], x_range[2], length.out = n_points)
  }

  check_bool(force_samples, "force_samples", allow_NA = FALSE)
  model_weights <- vapply(prior_list, .prior_model_weight, numeric(1))
  if(any(!is.finite(model_weights) | model_weights < 0) ||
     !is.finite(sum(model_weights)) || sum(model_weights) <= 0){
    stop("PET-PEESE prior model weights must be finite, nonnegative, and have positive total mass.", call. = FALSE)
  }
  supported <- function(prior){
    if(is.null(prior) || is.prior.none(prior)) return(TRUE)
    if(!is.prior(prior)) stop("PET-PEESE summaries require declared priors.", call. = FALSE)
    if(is.prior.spike_and_slab(prior)) return(supported(.get_spike_and_slab_variable(prior)))
    if(is.prior.mixture(prior)) return(all(vapply(prior, supported, logical(1))))
    is.prior.simple(prior) && (!is.prior.discrete(prior) || prior$distribution == "bernoulli")
  }
  eligible <- all(vapply(c(prior_list, prior_list_mu), supported, logical(1)))
  if(!force_samples && !is.list(transformation) && eligible){
    return(.plot_data_prior_list.PETPEESE_deterministic(
      prior_list               = prior_list,
      x_seq                    = x_seq,
      n_points                 = n_points,
      transformation           = transformation,
      transformation_arguments = transformation_arguments,
      prior_list_mu            = prior_list_mu,
      effect_direction         = effect_direction
    ))
  }

  .plot_data_prior_list.PETPEESE_sampled(
    prior_list               = prior_list,
    x_seq                    = x_seq,
    n_points                 = n_points,
    n_samples                = n_samples,
    transformation           = transformation,
    transformation_arguments = transformation_arguments,
    prior_list_mu            = prior_list_mu,
    effect_direction         = effect_direction
  )
}
.plot_data_prior_list.PETPEESE_deterministic <- function(prior_list, x_seq, n_points, transformation, transformation_arguments,
                                                         prior_list_mu, effect_direction = "positive"){

  pair <- .model_probability_petpeese_plot_check(prior_list, prior_list_mu)
  prior_list <- .model_probability_plot_priors(prior_list, validated_pair = pair)
  if(is.list(transformation)){
    stop("Custom transformations are handled by sampled PET-PEESE prior summaries.", call. = FALSE)
  }

  prior_weights <- sapply(prior_list, .prior_model_weight)
  keep <- is.finite(prior_weights) & prior_weights > 0
  prior_list    <- prior_list[keep]
  prior_list_mu <- prior_list_mu[keep]
  prior_weights <- prior_weights[keep]

  if(length(prior_list) == 0){
    stop("At least one PET-PEESE prior must have positive prior weight.", call. = FALSE)
  }

  model_weights <- prior_weights / sum(prior_weights)
  context <- .petpeese_prior_cdf_context(
    prior_list       = prior_list,
    prior_list_mu    = prior_list_mu,
    model_weights    = model_weights,
    effect_direction = effect_direction
  )

  quantiles <- vapply(x_seq, function(se){
    .petpeese_prior_cdf_quantile(context, se, c(.500, .025, .975))
  }, numeric(3))

  if(!is.null(transformation)){
    quantiles <- .petpeese_transform_quantiles(
      quantiles,
      transformation,
      transformation_arguments
    )
  }

  out <- list(
    call    = call("density", "PET-PEESE list"),
    bw      = NULL,
    n       = n_points,
    x       = x_seq,
    y       = quantiles[1,],
    y_lCI   = quantiles[2,],
    y_uCI   = quantiles[3,],
    samples = NULL
  )

  class(out) <- c("density", "density.prior", "density.prior.PETPEESE")
  attr(out, "x_range") <- range(x_seq)
  attr(out, "y_range") <- range(out$y, out$y_lCI, out$y_uCI)

  return(out)
}
.petpeese_prior_cdf_context <- function(prior_list, prior_list_mu, model_weights,
                                        effect_direction = "positive"){

  direction_sign <- if(effect_direction == "negative") -1 else 1

  models <- vector("list", length(prior_list))
  for(i in seq_along(prior_list)){
    bias_prior <- prior_list[[i]]
    if(is.prior.PET(bias_prior)){
      bias_type <- "PET"
      bias_components <- .petpeese_prior_components(bias_prior)
    }else if(is.prior.PEESE(bias_prior)){
      bias_type <- "PEESE"
      bias_components <- .petpeese_prior_components(bias_prior)
    }else{
      bias_type <- "none"
      bias_components <- .petpeese_prior_components(prior("point", list(location = 0)))
    }

    models[[i]] <- list(
      weight = model_weights[i],
      mu     = .petpeese_prior_components(prior_list_mu[[i]]),
      bias   = bias_components,
      type   = bias_type
    )
  }

  list(
    models         = models,
    direction_sign = direction_sign
  )
}
.petpeese_prior_components <- function(prior){

  if(is.null(prior) || is.prior.none(prior)){
    return(.petpeese_prior_components_normalize(list(
      list(weight = 1, type = "atom", x = 0)
    )))
  }

  if(is.prior.spike_and_slab(prior)){
    inclusion <- mean(.get_spike_and_slab_inclusion(prior))
    if(!is.finite(inclusion) || inclusion < 0 || inclusion > 1){
      stop("Spike-and-slab inclusion prior must have a finite mean in [0, 1].", call. = FALSE)
    }
    variable_components <- .petpeese_prior_components(.get_spike_and_slab_variable(prior))
    variable_components <- lapply(variable_components, function(component){
      component$weight <- component$weight * inclusion
      component
    })
    spike_component <- list(list(weight = 1 - inclusion, type = "atom", x = 0))
    return(.petpeese_prior_components_normalize(c(variable_components, spike_component)))
  }

  if(is.prior.mixture(prior)){
    weights <- attr(prior, "prior_weights")
    if(is.null(weights)){
      weights <- sapply(prior, .prior_model_weight)
    }
    weights <- weights / sum(weights)

    components <- list()
    for(i in seq_along(prior)){
      component_components <- .petpeese_prior_components(prior[[i]])
      component_components <- lapply(component_components, function(component){
        component$weight <- component$weight * weights[i]
        component
      })
      components <- c(components, component_components)
    }
    return(.petpeese_prior_components_normalize(components))
  }

  if(is.prior.point(prior)){
    return(.petpeese_prior_components_normalize(list(
      list(weight = 1, type = "atom", x = prior$parameters[["location"]])
    )))
  }

  if(is.prior.discrete(prior)){
    support <- switch(
      prior[["distribution"]],
      "bernoulli" = c(0, 1),
      stop("Unsupported discrete PET-PEESE prior distribution.", call. = FALSE)
    )
    probabilities <- mpdf(prior, support)
    keep <- is.finite(probabilities) & probabilities > 0
    if(!any(keep)){
      stop("Discrete PET-PEESE prior has zero probability mass.", call. = FALSE)
    }
    support <- support[keep]
    probabilities <- probabilities[keep] / sum(probabilities[keep])
    return(.petpeese_prior_components_normalize(lapply(seq_along(support), function(i){
      list(weight = probabilities[i], type = "atom", x = support[i])
    })))
  }

  if(is.prior.simple(prior)){
    prior_functions <- .petpeese_prior_simple_functions(prior)
    return(.petpeese_prior_components_normalize(list(
      list(
        weight = 1,
        type   = "continuous",
        prior  = prior,
        cdf    = prior_functions$cdf,
        ccdf   = prior_functions$ccdf,
        pdf    = prior_functions$pdf,
        quant  = prior_functions$quant,
        partition_cache = new.env(parent = emptyenv())
      )
    )))
  }

  stop("Unsupported PET-PEESE prior type for deterministic CDF plotting.", call. = FALSE)
}
.petpeese_prior_simple_functions <- function(prior){

  force(prior)
  list(
    cdf = function(q) .prior_simple_cdf(prior, q),
    ccdf = function(q) .prior_simple_ccdf(prior, q),
    pdf = function(x) .prior_simple_pdf(prior, x),
    quant = function(p) .prior_simple_quant(prior, p)
  )
}
.petpeese_prior_components_normalize <- function(components){

  weights <- vapply(components, `[[`, numeric(1), "weight")
  if(any(!is.finite(weights) | weights < 0) || !is.finite(sum(weights))){
    stop("PET-PEESE prior component weights must be finite and nonnegative.", call. = FALSE)
  }
  components <- components[vapply(components, function(component){
    is.finite(component$weight) && component$weight > 0
  }, logical(1))]

  if(length(components) == 0){
    stop("PET-PEESE prior components have zero total weight.", call. = FALSE)
  }

  total_weight <- sum(vapply(components, function(component) component$weight, numeric(1)))
  lapply(components, function(component){
    component$weight <- component$weight / total_weight
    component
  })
}
.petpeese_prior_cdf_quantile <- function(context, se, probs){

  vapply(probs, function(p){
    .petpeese_prior_cdf_one_quantile(context, se, p)
  }, numeric(1))
}
.petpeese_prior_cdf_one_quantile <- function(context, se, p){

  if(p <= 0){
    return(.petpeese_prior_cdf_range(context, se, tail_prob = .Machine$double.eps)[1])
  }
  if(p >= 1){
    return(.petpeese_prior_cdf_range(context, se, tail_prob = .Machine$double.eps)[2])
  }

  fast_quantile <- .petpeese_prior_fast_quantile(context, se, p)
  if(is.finite(fast_quantile)){
    return(fast_quantile)
  }

  bounds <- .petpeese_prior_cdf_range(context, se)
  lower <- bounds[1]
  upper <- bounds[2]

  if(!is.finite(lower) || !is.finite(upper) || lower >= upper){
    locations <- unlist(lapply(context$models, function(model){
      scale <- .petpeese_prior_scale(model$type, se, context$direction_sign)
      if(!all(vapply(model$mu, function(x) x$type == "atom", logical(1))) ||
         (scale != 0 && !all(vapply(model$bias, function(x) x$type == "atom", logical(1))))) return(NA_real_)
      unlist(lapply(model$mu, function(mu){
        if(scale == 0) return(mu$x)
        vapply(model$bias, function(bias) .petpeese_prior_sum_checked(mu$x, scale, bias$x), numeric(1))
      }), use.names = FALSE)
    }), use.names = FALSE)
    if(length(locations) > 0L && all(is.finite(locations)) && length(unique(locations)) == 1L) return(locations[[1L]])
    .petpeese_prior_numerical_stop("required_quantile_range", "Required interior quantiles have no finite distinct bracket")
  }

  cdf_lower <- .petpeese_prior_cdf(context, lower, se)
  cdf_upper <- .petpeese_prior_cdf(context, upper, se)
  width <- max(1, upper - lower, abs(lower), abs(upper))

  iter <- 0L
  while(is.finite(cdf_lower) && cdf_lower >= p && iter < 80L){
    if(!is.finite(width) || !is.finite(lower - width)){
      .petpeese_prior_numerical_stop("quantile bracket", "The required lower bracket lost representable range")
    }
    upper <- lower
    lower <- lower - width
    width <- width * 2
    cdf_lower <- .petpeese_prior_cdf(context, lower, se)
    iter <- iter + 1L
  }

  iter <- 0L
  while(is.finite(cdf_upper) && cdf_upper < p && iter < 80L){
    if(!is.finite(width) || !is.finite(upper + width)){
      .petpeese_prior_numerical_stop("quantile bracket", "The required upper bracket lost representable range")
    }
    lower <- upper
    upper <- upper + width
    width <- width * 2
    cdf_upper <- .petpeese_prior_cdf(context, upper, se)
    iter <- iter + 1L
  }

  if(!is.finite(cdf_lower) || !is.finite(cdf_upper) || cdf_lower >= p || cdf_upper < p){
    .petpeese_prior_numerical_stop("quantile bracket", "The unchanged bracket expansion budget was exhausted")
  }

  if(!.petpeese_prior_cdf_has_atoms(context, se)){
    root <- tryCatch(
      stats::uniroot(
        f        = function(q) .petpeese_prior_cdf(context, q, se) - p,
        interval = c(lower, upper),
        tol      = 1e-8 * max(1, abs(upper - lower))
      )$root,
      error = function(e){
        if(inherits(e, "BayesTools_numerical_condition")) stop(e)
        if(grepl("^convergence problem in zero finding", conditionMessage(e))) return(NA_real_)
        stop(e)
      }
    )
    if(is.finite(root)){
      return(root)
    }
  }

  converged <- FALSE
  for(iter in seq_len(100L)){
    mid <- lower / 2 + upper / 2
    cdf_mid <- .petpeese_prior_cdf(context, mid, se)
    if(!is.finite(cdf_mid)){
      stop("PET-PEESE prior CDF returned a non-finite value.", call. = FALSE)
    }

    if(cdf_mid >= p){
      upper <- mid
    }else{
      lower <- mid
    }

    if(is.finite(upper - lower) && abs(upper - lower) <= 1e-8 * max(1, abs(lower), abs(upper))){
      converged <- TRUE
      break
    }
  }

  if(!converged) .petpeese_prior_numerical_stop("quantile bisection", "The unchanged bisection budget was exhausted")
  upper
}
.petpeese_prior_numerical_stop <- function(operation, reason){

  .prior_numerical_signal(operation, "PET-PEESE", "full precision", integer(), reason, error = TRUE)
}

.petpeese_prior_scale <- function(type, se, direction_sign){

  if(!is.finite(se)) .petpeese_prior_numerical_stop("bias scale", "The standard error must be finite")
  out <- switch(type, "PET" = direction_sign * se,
    "PEESE" = direction_sign * se^2, "none" = 0)
  if(type != "none" && se != 0 && !.prior_density_full_precision(out)){
    .petpeese_prior_numerical_stop("bias scale", "A nonzero standard-error scale lost representable precision")
  }
  out
}

.petpeese_prior_sum_checked <- function(location, scale, bias){

  product <- scale * bias
  out <- location + product
  bad <- !is.finite(out) | (scale != 0 & bias != 0 &
    (!.prior_density_full_precision(product) | (location != 0 & out == location)))
  if(any(bad)) .petpeese_prior_numerical_stop("sum", "Required location and bias arithmetic lost representable precision")
  out
}

.petpeese_prior_cdf_has_atoms <- function(context, se){

  any(vapply(context$models, function(model){
    if(model$weight <= 0){
      return(FALSE)
    }
    scale <- .petpeese_prior_scale(model$type, se, context$direction_sign)

    any(vapply(model$mu, function(mu_component){
      any(vapply(model$bias, function(bias_component){
        .petpeese_prior_sum_has_atom(mu_component, bias_component, scale)
      }, logical(1)))
    }, logical(1)))
  }, logical(1)))
}
.petpeese_prior_sum_has_atom <- function(mu_component, bias_component, scale){

  if(scale == 0){
    return(mu_component$type == "atom")
  }

  mu_component$type == "atom" && bias_component$type == "atom"
}
.petpeese_prior_fast_quantile <- function(context, se, p){

  if(length(context$models) != 1L){
    return(NA_real_)
  }

  model <- context$models[[1]]
  if(length(model$mu) != 1L || length(model$bias) != 1L){
    return(NA_real_)
  }

  scale <- .petpeese_prior_scale(model$type, se, context$direction_sign)

  .petpeese_prior_sum_quantile(model$mu[[1]], model$bias[[1]], scale, p)
}
.petpeese_prior_sum_quantile <- function(mu_component, bias_component, scale, p){

  if(scale == 0){
    return(.petpeese_prior_component_quantile(mu_component, p))
  }

  if(mu_component$type == "atom" && bias_component$type == "atom"){
    return(.petpeese_prior_sum_checked(mu_component$x, scale, bias_component$x))
  }

  if(mu_component$type == "atom"){
    bias_p <- if(scale > 0) p else 1 - p
    return(.petpeese_prior_sum_checked(mu_component$x, scale, .petpeese_prior_component_quantile(bias_component, bias_p)))
  }

  if(bias_component$type == "atom"){
    return(.petpeese_prior_sum_checked(.petpeese_prior_component_quantile(mu_component, p), scale, bias_component$x))
  }

  NA_real_
}
.petpeese_prior_cdf <- function(context, q, se){

  values <- lapply(context$models, function(model){
    .petpeese_prior_model_cdf(model, q, se, context$direction_sign)
  })
  weights <- vapply(context$models, `[[`, numeric(1), "weight")
  cdf <- sum(weights * vapply(values, as.numeric, numeric(1)))
  error <- sum(weights * vapply(values, function(value){
    error <- attr(value, "absolute_error", exact = TRUE)
    if(is.null(error)) 0 else error
  }, numeric(1)))

  probability_bound <- max(error, 16 * .Machine$double.eps * max(1, abs(cdf)))
  if(!is.finite(cdf) || cdf < -probability_bound || cdf > 1 + probability_bound){
    .petpeese_prior_numerical_stop("mixture CDF", "The mixture probability is outside its numerical error criterion")
  }
  attr(cdf, "absolute_error") <- error
  attr(cdf, "numerical_diagnostics") <- lapply(values, attr, which = "numerical_diagnostics", exact = TRUE)
  cdf
}
.petpeese_prior_model_cdf <- function(model, q, se, direction_sign){

  scale <- .petpeese_prior_scale(model$type, se, direction_sign)

  cdf <- 0
  error <- 0
  diagnostics <- list()
  for(mu_component in model$mu){
    for(bias_component in model$bias){
      value <- .petpeese_prior_sum_cdf(mu_component, bias_component, scale, q)
      weight <- mu_component$weight * bias_component$weight
      cdf <- cdf + weight * as.numeric(value)
      pair_error <- attr(value, "absolute_error", exact = TRUE)
      if(!is.null(pair_error)) error <- error + weight * pair_error
      diagnostics[[length(diagnostics) + 1L]] <- list(weight = weight,
        integration = attr(value, "numerical_diagnostics", exact = TRUE))
    }
  }

  attr(cdf, "absolute_error") <- error
  attr(cdf, "numerical_diagnostics") <- diagnostics
  cdf
}
.petpeese_prior_sum_cdf <- function(mu_component, bias_component, scale, q){

  if(scale == 0){
    return(.petpeese_prior_component_cdf(mu_component, q))
  }

  if(mu_component$type == "atom" && bias_component$type == "atom"){
    return(as.numeric(q >= .petpeese_prior_sum_checked(mu_component$x, scale, bias_component$x)))
  }

  if(mu_component$type == "atom"){
    threshold <- .prior_region_inverse_affine(q, mu_component$x, scale)
    if(scale > 0){
      return(.petpeese_prior_component_cdf(bias_component, threshold))
    }else{
      return(.petpeese_prior_component_ccdf(bias_component, threshold))
    }
  }

  if(bias_component$type == "atom"){
    return(.petpeese_prior_component_cdf(mu_component, .petpeese_prior_sum_checked(q, -scale, bias_component$x)))
  }

  .petpeese_prior_probability_integral(mu_component, bias_component, scale, q)
}

.petpeese_prior_probability_integral <- function(mu, bias, scale, q){

  if(is.infinite(q)) return(as.numeric(q > 0))
  if(!is.finite(q)) .petpeese_prior_numerical_stop("pair CDF", "The event threshold is unavailable")
  omitted <- list()
  optional <- function(expression, operation){
    condition <- NULL
    value <- tryCatch(withCallingHandlers(expression,
      BayesTools_numerical_condition = function(e){
        if(inherits(e, "warning")){
          condition <<- e
          invokeRestart("muffleWarning")
        }
      }), BayesTools_numerical_condition = function(e){ condition <<- e; NULL })
    if(!is.null(condition)){
      omitted[[length(omitted) + 1L]] <<- list(operation = operation, condition = condition)
      return(NULL)
    }
    value
  }
  probabilities <- c(1e-6, 1e-3, .02, .25, .5, .75, .98, 1 - 1e-3, 1 - 1e-6)
  cache <- bias$partition_cache
  if(is.null(cache)) cache <- new.env(parent = emptyenv())
  if(!exists("quantiles", cache, inherits = FALSE)){
    quantiles <- lapply(probabilities, function(p){
      optional(.petpeese_prior_component_quantile(bias, p), paste0("bias quantile ", p))
    })
    cache$quantiles <- quantiles
    cache$omitted <- omitted
    omitted <- list()
  }
  omitted <- c(omitted, cache$omitted)
  support <- unlist(bias$prior$truncation[c("lower", "upper")], use.names = FALSE)
  finite_support <- support[is.finite(support)]
  transitions <- lapply(finite_support, function(bound){
    .petpeese_prior_sum_checked(q, -scale, bound)
  })
  quantile_transitions <- lapply(cache$quantiles, function(value){
    if(is.null(value)) return(NULL)
    optional(.petpeese_prior_sum_checked(q, -scale, value), "bias-quantile anchor")
  })
  knots <- c(0, .25, .5, .75, 1)
  add_anchor <- function(value){
    if(is.null(value)) return(invisible(NULL))
    probability <- optional(.petpeese_prior_component_cdf(mu, value), "location CDF anchor")
    if(!is.null(probability)) knots <<- c(knots, probability)
    invisible(NULL)
  }
  # Finite declared support transitions are required partition boundaries.
  # Only valid endpoint images (0/1) are redundant as additional knots.
  for(value in transitions){
    probability <- .petpeese_prior_component_cdf(mu, value)
    if(probability > 0 && probability < 1) knots <- c(knots, probability)
  }
  for(value in quantile_transitions) add_anchor(value)
  median <- optional(.petpeese_prior_component_quantile(mu, .5), "location median anchor")
  distance_transitions <- transitions
  if(length(finite_support) == 0L){
    bias_median <- cache$quantiles[[which(probabilities == .5)]]
    distance_transitions <- list(if(!is.null(bias_median)){
      optional(.petpeese_prior_sum_checked(q, -scale, bias_median), "bias-median transition")
    })
  }
  if(!is.null(median)){
    for(transition in distance_transitions){
      if(is.null(transition)) next
      distance <- optional({
        value <- abs(transition - median)
        if(!is.finite(value) || (value != 0 && !.prior_density_full_precision(value))){
          .petpeese_prior_numerical_stop("median distance", "An optional anchor distance lost representable precision")
        }
        value
      }, "median-distance anchor")
      if(is.null(distance)) next
      for(multiplier in c(.1, .5, 1, 2, 10)){
        for(direction in c(-1, 1)) add_anchor(optional(
          .petpeese_prior_sum_checked(median, direction * multiplier, distance), "median-distance image"))
      }
    }
  }
  knots <- sort(unique(knots))
  panels <- length(knots) - 1L
  if(panels < 1L || panels > 200L){
    .petpeese_prior_numerical_stop("pair partition", "The exact partition exceeds the subdivision budget")
  }
  remaining <- 200L
  results <- vector("list", panels)
  absolute_tolerance <- 1e-7 / (2 * panels)
  relative_tolerance <- 1e-7 / 2
  integrand <- function(u, upper_u = NULL){
    location <- .petpeese_prior_component_quantile(mu, u, upper_p = upper_u)
    threshold <- .prior_region_inverse_affine(q, location, scale)
    if(scale > 0) .petpeese_prior_component_cdf(bias, threshold) else .petpeese_prior_component_ccdf(bias, threshold)
  }
  for(i in seq_len(panels)){
    cap <- remaining - (panels - i)
    lower <- knots[[i]]
    upper <- knots[[i + 1L]]
    width <- upper - lower
    if(!is.finite(width) || width <= 0){
      .petpeese_prior_numerical_stop("pair integration", "A probability-panel width is nonfinite or collapsed")
    }
    panel_integrand <- function(t){
      u <- lower + width * t
      upper_u <- (1 - upper) + width * (1 - t)
      if(any(!is.finite(u) | u < lower | u > upper | u < 0 | u > 1 |
             !is.finite(upper_u) | upper_u < 0 | upper_u > 1)){
        .petpeese_prior_numerical_stop("pair integration", "A mapped probability-panel node lost its interior coordinate")
      }
      probability <- integrand(u, upper_u)
      value <- width * probability
      if(any(!is.finite(value) | (probability != 0 & value == 0))){
        .petpeese_prior_numerical_stop("pair integration", "A nonzero weighted probability-panel value was lost")
      }
      value
    }
    result <- stats::integrate(panel_integrand, 0, 1,
      subdivisions = cap, abs.tol = absolute_tolerance, rel.tol = relative_tolerance,
      stop.on.error = FALSE)
    result$physical_interval <- c(lower, upper)
    result$coordinate <- "affine probability panel"
    result$jacobian <- width
    consumed <- result$subdivisions
    local_criterion <- max(absolute_tolerance, relative_tolerance * abs(result$value))
    valid <- is.finite(result$value) && is.finite(result$abs.error) && result$abs.error >= 0 &&
      length(consumed) == 1L && is.finite(consumed) && consumed >= 1L && consumed <= cap
    false_maximum <- valid && identical(result$message, "maximum number of subdivisions reached") &&
      consumed == cap && result$abs.error <= local_criterion
    result$status <- if(false_maximum) "criterion_converged_at_cap" else result$message
    result$cap <- cap
    results[[i]] <- result
    if(!valid || (!identical(result$message, "OK") && !false_maximum)){
      .petpeese_prior_numerical_stop("pair integration", paste0("Integration reported '", result$message, "'"))
    }
    remaining <- remaining - consumed
  }
  value <- sum(vapply(results, `[[`, numeric(1), "value"))
  error <- sum(vapply(results, `[[`, numeric(1), "abs.error"))
  probability_bound <- max(error, 16 * .Machine$double.eps * max(1, abs(value)))
  if(!is.finite(value) || !is.finite(error) || error > max(1e-7, 1e-7 * abs(value)) ||
     value < -probability_bound || value > 1 + probability_bound){
    .petpeese_prior_numerical_stop("pair integration", "The total probability failed its unchanged numerical error criterion")
  }
  attr(value, "absolute_error") <- error
  attr(value, "numerical_diagnostics") <- list(method = "probability-coordinate quadrature",
    knots = knots, panels = results, subdivisions = 200L - remaining,
    optional_anchor_diagnostics = omitted, absolute_error = error, raw_value = value)
  value
}
.petpeese_prior_component_cdf <- function(component, q){

  if(component$type == "atom"){
    return(as.numeric(q >= component$x))
  }

  out <- component$cdf(q)
  if(any(!is.finite(out) | out < 0 | out > 1)) .petpeese_prior_numerical_stop("component CDF", "A primitive CDF did not return valid probabilities")
  out
}
.petpeese_prior_component_ccdf <- function(component, q){

  if(component$type == "atom"){
    return(as.numeric(q <= component$x))
  }

  out <- component$ccdf(q)
  if(any(!is.finite(out) | out < 0 | out > 1)) .petpeese_prior_numerical_stop("component CCDF", "A primitive CCDF did not return valid probabilities")
  out
}
.petpeese_prior_component_pdf <- function(component, x){

  if(component$type == "atom"){
    return(ifelse(x == component$x, Inf, 0))
  }

  component$pdf(x)
}
.petpeese_prior_component_quantile <- function(component, p, upper_p = NULL){

  if(component$type == "atom"){
    return(component$x)
  }

  if(is.null(upper_p)){
    out <- component$quant(p)
    required <- p > 0 & p < 1
  }else{
    if(length(p) != length(upper_p) || any(!is.finite(p) | !is.finite(upper_p) |
      p < 0 | p > 1 | upper_p < 0 | upper_p > 1)){
      .petpeese_prior_numerical_stop("component quantile", "A directly represented probability tail is invalid")
    }
    lower_tail <- p <= upper_p
    probability <- ifelse(lower_tail, p, upper_p)
    if(any(probability <= 0 | probability >= 1)){
      .petpeese_prior_numerical_stop("component quantile", "A required interior probability tail was lost")
    }
    if(.is_prior_default_range(component$prior)){
      out <- numeric(length(p))
      if(any(lower_tail)) out[lower_tail] <- .prior_simple_base_q(component$prior, p[lower_tail])
      if(any(!lower_tail)) out[!lower_tail] <- .prior_simple_base_q(component$prior,
        upper_p[!lower_tail], lower.tail = FALSE)
    }else{
      if(any(p <= 0 | p >= 1)){
        .petpeese_prior_numerical_stop("component quantile", "A required truncated lower probability lost its interior coordinate")
      }
      out <- component$quant(p)
    }
    required <- rep(TRUE, length(p))
  }
  if(any(required & (!is.finite(out) | (out != 0 & !.prior_density_full_precision(out))))){
    .petpeese_prior_numerical_stop("component quantile", "A required interior quantile lost representable precision")
  }
  out
}
.petpeese_prior_cdf_range <- function(context, se, tail_prob = .prior_linear_density_tail_prob()){

  ranges <- do.call(rbind, lapply(context$models, function(model){
    .petpeese_prior_model_range(model, se, context$direction_sign, tail_prob)
  }))

  if(any(!is.finite(ranges))) .petpeese_prior_numerical_stop("quantile range", "A required quantile range is nonfinite")
  range(ranges)
}
.petpeese_prior_model_range <- function(model, se, direction_sign, tail_prob){

  scale <- .petpeese_prior_scale(model$type, se, direction_sign)

  ranges <- list()
  for(mu_component in model$mu){
    mu_range <- .petpeese_prior_component_range(mu_component, tail_prob)
    for(bias_component in model$bias){
      if(scale == 0){
        ranges[[length(ranges) + 1L]] <- mu_range
      }else{
        bias_range <- .petpeese_prior_component_range(bias_component, tail_prob)
        scaled_bias_range <- sort(.petpeese_prior_sum_checked(0, scale, bias_range))
        ranges[[length(ranges) + 1L]] <- c(
          .petpeese_prior_sum_checked(mu_range[1], 1, scaled_bias_range[1]),
          .petpeese_prior_sum_checked(mu_range[2], 1, scaled_bias_range[2])
        )
      }
    }
  }

  ranges <- do.call(rbind, ranges)
  range(ranges)
}
.petpeese_prior_component_range <- function(component, tail_prob){

  if(component$type == "atom"){
    return(rep(component$x, 2))
  }

  .petpeese_prior_component_quantile(component, c(tail_prob, 1 - tail_prob))
}
.petpeese_transform_quantiles <- function(quantiles, transformation, transformation_arguments){

  transformed <- .density.prior_transformation_checked_x(
    as.vector(quantiles),
    transformation,
    transformation_arguments
  )
  transformed <- matrix(transformed, nrow = nrow(quantiles), ncol = ncol(quantiles))

  if(any(!is.finite(transformed))){
    stop("PET-PEESE transformed prior quantiles are non-finite.", call. = FALSE)
  }

  rbind(
    transformed[1,],
    pmin(transformed[2,], transformed[3,]),
    pmax(transformed[2,], transformed[3,])
  )
}
.plot_data_prior_list.PETPEESE_sampled <- function(prior_list, x_seq, n_points, n_samples,
                                                   transformation, transformation_arguments, prior_list_mu,
                                                   effect_direction = "positive"){

  pair <- .model_probability_petpeese_plot_check(prior_list, prior_list_mu)
  prior_list <- .model_probability_plot_priors(prior_list, validated_pair = pair)
  prior_weights  <- sapply(prior_list, .prior_model_weight)
  mixing_prop    <- prior_weights / sum(prior_weights)

  prior_list     <- prior_list[round(n_samples * mixing_prop) > 0]
  prior_list_mu  <- prior_list_mu[round(n_samples * mixing_prop) > 0]
  mixing_prop    <- mixing_prop[round(n_samples * mixing_prop) > 0]

  # get the samples
  samples_list <- list()
  for(i in seq_along(prior_list)){
    if(is.prior.PET(prior_list[[i]])){
      samples_list[[i]] <- cbind(rng(prior_list_mu[[i]], round(n_samples * mixing_prop[i])), rng(prior_list[[i]], round(n_samples * mixing_prop[i])), rep(0, length = round(n_samples * mixing_prop[i])))
    }else if(is.prior.PEESE(prior_list[[i]])){
      samples_list[[i]] <- cbind(rng(prior_list_mu[[i]], round(n_samples * mixing_prop[i])), rep(0, length = round(n_samples * mixing_prop[i])), rng(prior_list[[i]], round(n_samples * mixing_prop[i])))
    }else{
      samples_list[[i]] <- cbind(rng(prior_list_mu[[i]], round(n_samples * mixing_prop[i])), matrix(0, nrow = round(n_samples * mixing_prop[i]), ncol = 2))
    }
  }
  samples <- do.call(rbind, samples_list)

  summary <- .petpeese_line_summary_from_samples(
    samples                  = samples,
    x_seq                    = x_seq,
    transformation           = transformation,
    transformation_arguments = transformation_arguments,
    effect_direction         = effect_direction,
    quantile_method          = "empirical"
  )


  out <- list(
    call    = call("density", "PET-PEESE list"),
    bw      = NULL,
    n       = n_points,
    x       = x_seq,
    y       = summary$median,
    y_lCI   = summary$lCI,
    y_uCI   = summary$uCI,
    samples = summary$samples
  )


  class(out) <- c("density", "density.prior", "density.prior.PETPEESE")
  attr(out, "x_range") <- range(x_seq)
  attr(out, "y_range") <- range(out$y, out$y_lCI, out$y_uCI)

  return(out)
}
