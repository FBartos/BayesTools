.plot_data_prior_list.PETPEESE       <- function(prior_list, x_seq, x_range, x_range_quant, n_points, n_samples,
                                                 transformation, transformation_arguments, transformation_settings, prior_list_mu,
                                                 effect_direction = "positive"){

  if(is.null(x_seq)){
    x_seq <- seq(x_range[1], x_range[2], length.out = n_points)
  }

  # specify it on the transformed range if requested
  if(transformation_settings & !is.null(transformation)){
    x_seq   <- .density.prior_transformation_x(x_seq,   transformation, transformation_arguments)
    x_range <- .density.prior_transformation_x(x_range, transformation, transformation_arguments)
  }

  deterministic <- tryCatch(
    .plot_data_prior_list.PETPEESE_deterministic(
      prior_list               = prior_list,
      x_seq                    = x_seq,
      n_points                 = n_points,
      transformation           = transformation,
      transformation_arguments = transformation_arguments,
      prior_list_mu            = prior_list_mu,
      effect_direction         = effect_direction
    ),
    error = function(e) NULL
  )

  if(!is.null(deterministic)){
    return(deterministic)
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
        quant  = prior_functions$quant
      )
    )))
  }

  stop("Unsupported PET-PEESE prior type for deterministic CDF plotting.", call. = FALSE)
}
.petpeese_prior_simple_functions <- function(prior){

  default_range <- .is_prior_default_range(prior)
  if(default_range){
    return(list(
      cdf = function(q) .prior_simple_base_p(prior, q, lower.tail = TRUE),
      ccdf = function(q) .prior_simple_base_p(prior, q, lower.tail = FALSE),
      pdf = function(x) .prior_simple_base_d(prior, x, log = FALSE),
      quant = function(p) .prior_simple_base_q(prior, p)
    ))
  }

  C1 <- .prior_C1(prior)
  C2 <- .prior_C2(prior)
  C  <- .prior_C(prior)
  lower <- prior$truncation[["lower"]]
  upper <- prior$truncation[["upper"]]
  use_survival <- .prior_simple_use_survival_truncation(prior)
  if(use_survival){
    S1 <- .prior_simple_base_p(prior, lower, lower.tail = FALSE)
    S2 <- .prior_simple_base_p(prior, upper, lower.tail = FALSE)
  }

  list(
    cdf = function(q){
      p <- numeric(length(q))
      q_lower  <- q < lower
      q_higher <- q > upper
      q_inside <- !q_lower & !q_higher

      p[q_lower]  <- 0
      p[q_higher] <- 1
      if(any(q_inside)){
        if(use_survival){
          S_q         <- .prior_simple_base_p(prior, q[q_inside], lower.tail = FALSE)
          p[q_inside] <- (S1 - S_q) / C
        }else{
          p[q_inside] <- (.prior_simple_base_p(prior, q[q_inside], lower.tail = TRUE) - C1) / C
        }
      }
      p
    },
    ccdf = function(q){
      p <- numeric(length(q))
      q_lower  <- q < lower
      q_higher <- q > upper
      q_inside <- !q_lower & !q_higher

      p[q_lower]  <- 1
      p[q_higher] <- 0
      if(any(q_inside)){
        if(use_survival){
          p[q_inside] <- (.prior_simple_base_p(prior, q[q_inside], lower.tail = FALSE) - S2) / C
        }else{
          p[q_inside] <- (.prior_simple_base_p(prior, q[q_inside], lower.tail = FALSE) - (1 - C2)) / C
        }
      }
      p
    },
    pdf = function(x){
      y <- .prior_simple_base_d(prior, x, log = FALSE)
      y[x < lower | x > upper] <- 0
      y / C
    },
    quant = function(p){
      if(use_survival){
        return(.prior_simple_base_q(prior, S1 - p * (S1 - S2), lower.tail = FALSE))
      }
      .prior_simple_base_q(prior, C1 + p * C)
    }
  )
}
.petpeese_prior_components_normalize <- function(components){

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

  if(!is.finite(lower) || !is.finite(upper)){
    lower <- -1
    upper <-  1
  }

  if(isTRUE(all.equal(lower, upper))){
    return(lower)
  }

  cdf_lower <- .petpeese_prior_cdf(context, lower, se)
  cdf_upper <- .petpeese_prior_cdf(context, upper, se)
  width <- max(1, upper - lower, abs(lower), abs(upper))

  iter <- 0L
  while(is.finite(cdf_lower) && cdf_lower >= p && iter < 80L){
    upper <- lower
    lower <- lower - width
    width <- width * 2
    cdf_lower <- .petpeese_prior_cdf(context, lower, se)
    iter <- iter + 1L
  }

  iter <- 0L
  while(is.finite(cdf_upper) && cdf_upper < p && iter < 80L){
    lower <- upper
    upper <- upper + width
    width <- width * 2
    cdf_upper <- .petpeese_prior_cdf(context, upper, se)
    iter <- iter + 1L
  }

  if(!is.finite(cdf_lower) || !is.finite(cdf_upper) || cdf_lower >= p || cdf_upper < p){
    stop("Could not bracket PET-PEESE prior quantile.", call. = FALSE)
  }

  if(!.petpeese_prior_cdf_has_atoms(context, se)){
    root <- tryCatch(
      stats::uniroot(
        f        = function(q) .petpeese_prior_cdf(context, q, se) - p,
        interval = c(lower, upper),
        tol      = 1e-8 * max(1, abs(upper - lower))
      )$root,
      error = function(e) NA_real_
    )
    if(is.finite(root)){
      return(root)
    }
  }

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

    if(abs(upper - lower) <= 1e-8 * max(1, abs(lower), abs(upper))){
      break
    }
  }

  upper
}
.petpeese_prior_cdf_has_atoms <- function(context, se){

  any(vapply(context$models, function(model){
    if(model$weight <= 0){
      return(FALSE)
    }
    scale <- switch(
      model$type,
      "PET"   = context$direction_sign * se,
      "PEESE" = context$direction_sign * se^2,
      "none"  = 0
    )

    any(vapply(model$mu, function(mu_component){
      any(vapply(model$bias, function(bias_component){
        .petpeese_prior_sum_has_atom(mu_component, bias_component, scale)
      }, logical(1)))
    }, logical(1)))
  }, logical(1)))
}
.petpeese_prior_sum_has_atom <- function(mu_component, bias_component, scale){

  if(abs(scale) <= .prior_linear_density_zero_tol()){
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

  scale <- switch(
    model$type,
    "PET"   = context$direction_sign * se,
    "PEESE" = context$direction_sign * se^2,
    "none"  = 0
  )

  .petpeese_prior_sum_quantile(model$mu[[1]], model$bias[[1]], scale, p)
}
.petpeese_prior_sum_quantile <- function(mu_component, bias_component, scale, p){

  if(abs(scale) <= .prior_linear_density_zero_tol()){
    return(.petpeese_prior_component_quantile(mu_component, p))
  }

  if(mu_component$type == "atom" && bias_component$type == "atom"){
    return(mu_component$x + scale * bias_component$x)
  }

  if(mu_component$type == "atom"){
    bias_p <- if(scale > 0) p else 1 - p
    return(mu_component$x + scale * .petpeese_prior_component_quantile(bias_component, bias_p))
  }

  if(bias_component$type == "atom"){
    return(.petpeese_prior_component_quantile(mu_component, p) + scale * bias_component$x)
  }

  NA_real_
}
.petpeese_prior_cdf <- function(context, q, se){

  cdf <- sum(vapply(context$models, function(model){
    model$weight * .petpeese_prior_model_cdf(model, q, se, context$direction_sign)
  }, numeric(1)))

  pmin(pmax(cdf, 0), 1)
}
.petpeese_prior_model_cdf <- function(model, q, se, direction_sign){

  scale <- switch(
    model$type,
    "PET"   = direction_sign * se,
    "PEESE" = direction_sign * se^2,
    "none"  = 0
  )

  cdf <- 0
  for(mu_component in model$mu){
    for(bias_component in model$bias){
      cdf <- cdf + mu_component$weight * bias_component$weight *
        .petpeese_prior_sum_cdf(mu_component, bias_component, scale, q)
    }
  }

  cdf
}
.petpeese_prior_sum_cdf <- function(mu_component, bias_component, scale, q){

  if(abs(scale) <= .prior_linear_density_zero_tol()){
    return(.petpeese_prior_component_cdf(mu_component, q))
  }

  if(mu_component$type == "atom" && bias_component$type == "atom"){
    return(as.numeric(q >= mu_component$x + scale * bias_component$x))
  }

  if(mu_component$type == "atom"){
    threshold <- (q - mu_component$x) / scale
    if(scale > 0){
      return(.petpeese_prior_component_cdf(bias_component, threshold))
    }else{
      return(.petpeese_prior_component_ccdf(bias_component, threshold))
    }
  }

  if(bias_component$type == "atom"){
    return(.petpeese_prior_component_cdf(mu_component, q - scale * bias_component$x))
  }

  integration <- stats::integrate(
    f = function(b){
      .petpeese_prior_component_cdf(mu_component, q - scale * b) *
        .petpeese_prior_component_pdf(bias_component, b)
    },
    lower        = bias_component$prior$truncation[["lower"]],
    upper        = bias_component$prior$truncation[["upper"]],
    subdivisions = 200L,
    rel.tol      = 1e-7,
    stop.on.error = FALSE
  )

  if(!isTRUE(integration$message == "OK") && !is.finite(integration$value)){
    stop("PET-PEESE prior CDF integration failed.", call. = FALSE)
  }

  pmin(pmax(integration$value, 0), 1)
}
.petpeese_prior_component_cdf <- function(component, q){

  if(component$type == "atom"){
    return(as.numeric(q >= component$x))
  }

  component$cdf(q)
}
.petpeese_prior_component_ccdf <- function(component, q){

  if(component$type == "atom"){
    return(as.numeric(q <= component$x))
  }

  component$ccdf(q)
}
.petpeese_prior_component_pdf <- function(component, x){

  if(component$type == "atom"){
    return(ifelse(x == component$x, Inf, 0))
  }

  component$pdf(x)
}
.petpeese_prior_component_quantile <- function(component, p){

  if(component$type == "atom"){
    return(component$x)
  }

  component$quant(p)
}
.petpeese_prior_cdf_range <- function(context, se, tail_prob = .prior_linear_density_tail_prob()){

  ranges <- do.call(rbind, lapply(context$models, function(model){
    .petpeese_prior_model_range(model, se, context$direction_sign, tail_prob)
  }))

  range(ranges[,1], ranges[,2], finite = TRUE)
}
.petpeese_prior_model_range <- function(model, se, direction_sign, tail_prob){

  scale <- switch(
    model$type,
    "PET"   = direction_sign * se,
    "PEESE" = direction_sign * se^2,
    "none"  = 0
  )

  ranges <- list()
  for(mu_component in model$mu){
    mu_range <- .petpeese_prior_component_range(mu_component, tail_prob)
    for(bias_component in model$bias){
      if(abs(scale) <= .prior_linear_density_zero_tol()){
        ranges[[length(ranges) + 1L]] <- mu_range
      }else{
        bias_range <- .petpeese_prior_component_range(bias_component, tail_prob)
        scaled_bias_range <- sort(scale * bias_range)
        ranges[[length(ranges) + 1L]] <- c(
          mu_range[1] + scaled_bias_range[1],
          mu_range[2] + scaled_bias_range[2]
        )
      }
    }
  }

  ranges <- do.call(rbind, ranges)
  range(ranges[,1], ranges[,2], finite = TRUE)
}
.petpeese_prior_component_range <- function(component, tail_prob){

  if(component$type == "atom"){
    return(rep(component$x, 2))
  }

  component$quant(c(tail_prob, 1 - tail_prob))
}
.petpeese_transform_quantiles <- function(quantiles, transformation, transformation_arguments){

  transformed <- .density.prior_transformation_x(
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
    effect_direction         = effect_direction
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
