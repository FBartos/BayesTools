.plot_data_prior_list.weightfunction <- function(prior_list, x_seq, x_range, x_range_quant, n_points, n_samples){

  context <- .weightfunction_prior_list_context(prior_list)

  omega_cuts <- context$omega_cuts
  x_mean     <- numeric(length(omega_cuts) - 1)
  x_lCI      <- numeric(length(omega_cuts) - 1)
  x_uCI      <- numeric(length(omega_cuts) - 1)

  for(i in seq_len(length(omega_cuts) - 1)){
    components <- .weightfunction_prior_marginal_components(context, i)

    x_mean[i] <- .weightfunction_mixture_mean(components)
    x_lCI[i]  <- .weightfunction_mixture_quantile(components, .025)
    x_uCI[i]  <- .weightfunction_mixture_quantile(components, .975)
  }

  x_seq     <- omega_cuts
  if(length(x_seq) > 2){
    x_seq_rep <- c(1, sort(rep(2:(length(x_seq)-1), 2)) ,length(x_seq))
    x_val_rep <- sort(rep(1:(length(x_seq)-1), 2))
  }else{
    x_seq_rep <- c(1, 2)
    x_val_rep <- c(1, 1)
  }


  out <- list(
    call    = call("density", "weightfunction list"),
    bw      = NULL,
    n       = n_points,
    x       = x_seq[x_seq_rep],
    y       = x_mean[x_val_rep],
    y_lCI   = x_lCI[x_val_rep],
    y_uCI   = x_uCI[x_val_rep],
    samples = NULL
  )


  class(out) <- c("density", "density.prior", "density.prior.weightfunction")
  attr(out, "x_range") <- c(0, 1)
  attr(out, "y_range") <- c(0, max(1, x_mean, x_lCI, x_uCI, na.rm = TRUE))

  return(out)
}
.plot_data_prior_list.weightparameter<- function(prior_list, parameter, n_points, n_samples){

  context       <- .weightfunction_prior_list_context(prior_list)
  parameter_ind <- match(parameter, context$omega_names)

  if(is.na(parameter_ind)){
    stop(paste0("Parameter '", parameter, "' not found in the weightfunction prior."), call. = FALSE)
  }

  components <- .weightfunction_prior_marginal_components(context, parameter_ind)
  out        <- .plot_data_prior_weightparameter_components(components, parameter, n_points)

  return(out)
}
.weightfunction_prior_list_context <- function(prior_list, one_sided = NULL){

  omega_context <- attr(prior_list, "omega_context")
  if(is.null(one_sided)){
    one_sided <- if(!is.null(omega_context) && !is.null(omega_context$one_sided)){
      isTRUE(omega_context$one_sided)
    }else{
      .weightfunction_prior_context_uses_selection_mapping(prior_list)
    }
  }else{
    check_bool(one_sided, "one_sided")
  }

  prior_list <- .weightfunction_expand_bias_mixture_priors(prior_list)
  prior_list <- .simplify_prior_list(prior_list)
  prior_list <- .weightfunction_expand_bias_mixture_priors(prior_list)

  prior_weights <- sapply(prior_list, .prior_model_weight)
  keep          <- is.finite(prior_weights) & prior_weights > 0
  prior_list    <- prior_list[keep]
  prior_weights <- prior_weights[keep]

  if(length(prior_list) == 0){
    stop("At least one weightfunction prior must have positive prior weight.", call. = FALSE)
  }

  # Non-weightfunction bias components imply no selection adjustment and
  # therefore correspond to publication weights fixed at one.
  for(i in seq_along(prior_list)){
    if(.weightfunction_prior_has_selection(prior_list[[i]])){
      selection_priors <- .selection_prior_selection_priors(prior_list[[i]])
      selection_prior  <- selection_priors[[1L]]
      selection_prior$prior_weights <- prior_weights[i]
      prior_list[[i]] <- selection_prior
    }else if(!(is.prior.weightfunction(prior_list[[i]]) | is.prior.none(prior_list[[i]]))){
      prior_list[[i]] <- prior_none(prior_weights = prior_weights[i])
    }
  }

  prior_weights <- sapply(prior_list, .prior_model_weight)
  model_weights <- prior_weights / sum(prior_weights)
  omega_info    <- .weightfunction_mapping_info(prior_list, one_sided = one_sided)
  omega_mapping <- omega_info$mapping
  omega_cuts    <- omega_info$cuts
  omega_names   <- omega_info$names

  list(
    prior_list    = prior_list,
    model_weights = model_weights,
    omega_mapping = omega_mapping,
    omega_cuts    = omega_cuts,
    omega_names   = omega_names
  )
}
.weightfunction_prior_context_uses_selection_mapping <- function(prior_list){

  if(inherits(prior_list, "prior.bias_mixture") || is_prior_bias(prior_list)){
    return(TRUE)
  }
  if(is.prior(prior_list)){
    return(FALSE)
  }
  if(is.list(prior_list)){
    return(any(vapply(prior_list, function(prior){
      inherits(prior, "prior.bias_mixture") || is_prior_bias(prior)
    }, logical(1))))
  }

  FALSE
}
.weightfunction_expand_bias_mixture_priors <- function(prior_list){

  if(is.prior.mixture(prior_list)){
    class(prior_list) <- NULL
    return(prior_list)
  }
  if(is.prior(prior_list)){
    return(list(prior_list))
  }

  expanded <- list()
  for(i in seq_along(prior_list)){
    if(inherits(prior_list[[i]], "prior.bias_mixture")){
      prior_mixture_components <- prior_list[[i]]
      class(prior_mixture_components) <- NULL
      expanded <- c(expanded, prior_mixture_components)
    }else{
      expanded[[length(expanded) + 1L]] <- prior_list[[i]]
    }
  }

  expanded
}
.weightfunction_prior_has_selection <- function(prior){

  if(is.prior.weightfunction(prior) || is_prior_bias(prior) || inherits(prior, "prior.bias_mixture")){
    return(.selection_prior_has_selection(prior))
  }

  FALSE
}
.weightfunction_prior_marginal_components <- function(context, parameter_ind){

  components <- list()

  for(i in seq_along(context$prior_list)){
    prior <- context$prior_list[[i]]

    if(is.prior.weightfunction(prior)){
      component <- .weightfunction_prior_component(
        prior  = prior,
        index  = context$omega_mapping[[i]][parameter_ind],
        weight = context$model_weights[i]
      )
    }else{
      component <- list(
        type     = "point",
        weight   = context$model_weights[i],
        location = 1
      )
    }

    components[[length(components) + 1L]] <- component
  }

  components <- components[vapply(components, function(component){
    is.finite(component$weight) && component$weight > 0
  }, logical(1))]
  total_weight <- sum(vapply(components, function(component) component$weight, numeric(1)))

  lapply(components, function(component){
    component$weight <- component$weight / total_weight
    component
  })
}
.weightfunction_prior_component <- function(prior, index, weight){

  if(prior$weights$type == "fixed"){
    return(list(
      type     = "point",
      weight   = weight,
      location = prior$weights[["omega"]][index]
    ))
  }

  if(prior$weights$type == "cumulative"){
    return(.weightfunction_prior_component_cumdirichlet(prior$weights[["alpha"]], index, weight))
  }

  if(prior$weights$type == "independent"){
    if(index == 1L){
      return(list(
        type     = "point",
        weight   = weight,
        location = 1
      ))
    }
    return(list(
      type   = "prior",
      weight = weight,
      prior  = prior$weights$prior,
      scale  = prior$weights$scale
    ))
  }

  stop("Unsupported weightfunction prior specification.", call. = FALSE)
}
.weightfunction_prior_component_cumdirichlet <- function(alpha, index, weight){

  if(index <= 1L){
    return(list(
      type     = "point",
      weight   = weight,
      location = 1
    ))
  }

  list(
    type   = "beta",
    weight = weight,
    alpha  = sum(alpha[index:length(alpha)]),
    beta   = sum(alpha[seq_len(index - 1L)])
  )
}
.weightfunction_component_mean <- function(component){

  switch(
    component$type,
    "point" = component$location,
    "beta"  = component$alpha / (component$alpha + component$beta),
    "prior" = if(component$scale == "omega"){
      mean(component$prior)
    }else{
      stats::integrate(
        f     = function(x, prior) {
          y <- exp(x) * pdf(prior, x)
          y[!is.finite(y)] <- 0
          y
        },
        lower = component$prior$truncation[["lower"]],
        upper = component$prior$truncation[["upper"]],
        prior = component$prior
      )$value
    },
    "one_minus_product_beta" = {
      1 - (component$u_alpha / (component$u_alpha + component$u_beta)) *
        (component$v_alpha / (component$v_alpha + component$v_beta))
    }
  )
}
.weightfunction_component_cdf <- function(component, q){

  switch(
    component$type,
    "point" = as.numeric(q >= component$location),
    "beta"  = stats::pbeta(q, shape1 = component$alpha, shape2 = component$beta),
    "prior" = if(component$scale == "omega"){
      mcdf(component$prior, q)
    }else{
      p <- numeric(length(q))
      p[q <= 0] <- 0
      inside <- q > 0
      p[inside] <- mcdf(component$prior, log(q[inside]))
      p
    },
    "one_minus_product_beta" = .weightfunction_one_minus_product_beta_cdf(
      q       = q,
      u_alpha = component$u_alpha,
      u_beta  = component$u_beta,
      v_alpha = component$v_alpha,
      v_beta  = component$v_beta
    )
  )
}
.weightfunction_component_pdf <- function(component, x){

  switch(
    component$type,
    "point" = rep(0, length(x)),
    "beta"  = stats::dbeta(x, shape1 = component$alpha, shape2 = component$beta),
    "prior" = if(component$scale == "omega"){
      mpdf(component$prior, x)
    }else{
      y <- numeric(length(x))
      inside <- x > 0
      y[inside] <- mpdf(component$prior, log(x[inside])) / x[inside]
      y
    },
    "one_minus_product_beta" = .weightfunction_one_minus_product_beta_pdf(
      x       = x,
      u_alpha = component$u_alpha,
      u_beta  = component$u_beta,
      v_alpha = component$v_alpha,
      v_beta  = component$v_beta
    )
  )
}
.weightfunction_component_range <- function(component, quantiles = .005){

  switch(
    component$type,
    "point" = c(component$location, component$location),
    "beta"  = c(0, 1),
    "prior" = {
      if(component$scale == "omega"){
        lower <- component$prior$truncation[["lower"]]
        upper <- component$prior$truncation[["upper"]]

        lower <- if(is.infinite(lower)) mquant(component$prior, quantiles) else lower
        upper <- if(is.infinite(upper)) mquant(component$prior, 1 - quantiles) else upper
        c(lower, upper)
      }else{
        lower <- component$prior$truncation[["lower"]]
        upper <- component$prior$truncation[["upper"]]

        lower <- if(is.infinite(lower)) 0 else exp(lower)
        upper <- if(is.infinite(upper)) exp(mquant(component$prior, 1 - quantiles)) else exp(upper)
        c(lower, upper)
      }
    },
    "one_minus_product_beta" = c(0, 1)
  )
}
.weightfunction_components_range <- function(components, samples = NULL, quantiles = .005){

  ranges <- do.call(rbind, lapply(components, .weightfunction_component_range, quantiles = quantiles))
  values <- as.vector(ranges)
  if(!is.null(samples)){
    values <- c(values, samples)
  }

  x_range <- range(values, finite = TRUE)
  if(!all(is.finite(x_range))){
    x_range <- c(0, 1)
  }
  if(x_range[1] == x_range[2]){
    x_range <- range(c(0, 1, values), finite = TRUE)
  }
  x_range[1] <- max(0, x_range[1])

  x_range
}
.weightfunction_component_quantile <- function(component, p){

  switch(
    component$type,
    "point" = component$location,
    "beta"  = stats::qbeta(p, shape1 = component$alpha, shape2 = component$beta),
    "prior" = if(component$scale == "omega"){
      mquant(component$prior, p)
    }else{
      exp(mquant(component$prior, p))
    },
    "one_minus_product_beta" = {
      lower <- 0
      upper <- 1
      for(iter in seq_len(100L)){
        mid <- lower / 2 + upper / 2
        if(.weightfunction_component_cdf(component, mid) >= p){
          upper <- mid
        }else{
          lower <- mid
        }
        if(abs(upper - lower) <= 1e-8){
          break
        }
      }
      upper
    }
  )
}
.weightfunction_mixture_mean <- function(components){

  sum(vapply(components, function(component){
    component$weight * .weightfunction_component_mean(component)
  }, numeric(1)))
}
.weightfunction_mixture_cdf <- function(components, q){

  p <- sum(vapply(components, function(component){
    component$weight * .weightfunction_component_cdf(component, q)
  }, numeric(1)))

  pmin(pmax(p, 0), 1)
}
.weightfunction_mixture_quantile <- function(components, p){

  if(p <= 0){
    return(0)
  }
  if(p >= 1){
    return(max(.weightfunction_components_range(components)))
  }

  lower <- 0
  upper <- max(vapply(components, .weightfunction_component_quantile, numeric(1), p = p))
  if(!is.finite(upper) || upper <= lower){
    upper <- max(.weightfunction_components_range(components))
  }

  for(iter in seq_len(100L)){
    mid <- lower / 2 + upper / 2
    if(.weightfunction_mixture_cdf(components, mid) >= p){
      upper <- mid
    }else{
      lower <- mid
    }

    if(abs(upper - lower) <= 1e-8){
      break
    }
  }

  upper
}
.plot_data_prior_weightparameter_components <- function(components, parameter, n_points){

  x_range <- .weightfunction_components_range(components)
  x_den <- seq(x_range[1], x_range[2], length.out = n_points)
  y_den <- rep(0, length(x_den))

  for(component in components){
    if(component$type != "point"){
      y_component <- .weightfunction_component_pdf(component, x_den)
      y_component[!is.finite(y_component)] <- 0
      y_den <- y_den + component$weight * y_component
    }
  }

  point_components <- components[vapply(components, function(component){
    component$type == "point"
  }, logical(1))]

  x_points <- NULL
  y_points <- NULL
  if(length(point_components) > 0){
    point_locations <- vapply(point_components, function(component) component$location, numeric(1))
    point_keys      <- as.character(signif(point_locations, 15))
    point_groups    <- split(seq_along(point_components), point_keys)

    x_points <- unname(vapply(point_groups, function(ind) point_locations[ind[1]], numeric(1)))
    y_points <- unname(vapply(point_groups, function(ind){
      sum(vapply(point_components[ind], function(component) component$weight, numeric(1)))
    }, numeric(1)))
  }

  out <- list()

  if(any(y_den > 0)){
    out_den <- list(
      call    = call("density", "weightfunction prior"),
      bw      = NULL,
      n       = n_points,
      x       = x_den,
      y       = y_den,
      samples = NULL
    )

    class(out_den) <- c("density", "density.prior", "density.prior.simple")
    attr(out_den, "x_range") <- x_range
    attr(out_den, "y_range") <- c(0, max(y_den))
    attr(out_den, "parameter") <- parameter

    out[["density"]] <- out_den
  }

  if(!is.null(y_points)){
    for(i in seq_along(y_points)){
      temp_points <- list(
        call    = call("density", paste0("point", i)),
        bw      = NULL,
        n       = n_points,
        x       = x_points[i],
        y       = y_points[i],
        samples = NULL
      )

      class(temp_points) <- c("density", "density.prior", "density.prior.point")
      attr(temp_points, "x_range") <- x_range
      attr(temp_points, "y_range") <- c(0, max(y_points[i]))
      attr(temp_points, "parameter") <- parameter

      out[[paste0("points",i)]] <- temp_points
    }
  }

  out
}
.weightfunction_one_minus_product_beta_cdf <- function(q, u_alpha, u_beta, v_alpha, v_beta){

  vapply(q, function(q_i){
    if(q_i <= 0){
      return(0)
    }
    if(q_i >= 1){
      return(1)
    }

    product_lower <- 1 - q_i
    integration <- stats::integrate(
      f = function(u){
        stats::dbeta(u, shape1 = u_alpha, shape2 = u_beta) *
          stats::pbeta(product_lower / u, shape1 = v_alpha, shape2 = v_beta, lower.tail = FALSE)
      },
      lower         = product_lower,
      upper         = 1,
      subdivisions  = 200L,
      rel.tol       = 1e-7,
      stop.on.error = FALSE
    )

    if(!is.finite(integration$value)){
      stop("Weightfunction prior CDF integration failed.", call. = FALSE)
    }

    pmin(pmax(integration$value, 0), 1)
  }, numeric(1))
}
.weightfunction_one_minus_product_beta_pdf <- function(x, u_alpha, u_beta, v_alpha, v_beta){

  vapply(x, function(x_i){
    if(x_i <= 0 || x_i >= 1){
      return(0)
    }

    product_value <- 1 - x_i
    integration <- stats::integrate(
      f = function(u){
        stats::dbeta(u, shape1 = u_alpha, shape2 = u_beta) *
          stats::dbeta(product_value / u, shape1 = v_alpha, shape2 = v_beta) / u
      },
      lower         = product_value,
      upper         = 1,
      subdivisions  = 200L,
      rel.tol       = 1e-7,
      stop.on.error = FALSE
    )

    if(!is.finite(integration$value)){
      return(NA_real_)
    }

    pmax(integration$value, 0)
  }, numeric(1))
}
