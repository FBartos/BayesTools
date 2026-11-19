# A centered shrinkage term contributes the normalized conditional law of
# its coefficients. Its local scales, and its pool's shared scales at the
# owning term, are ordinary scalar bridge coordinates.
.bt_regularized_coefficient_names <- function(prior, parameter_name){

  if(is.prior.factor(prior)) .JAGS_prior_factor_names(parameter_name, prior) else parameter_name
}

.bt_regularized_bridge_info <- function(prior, parameter_name){

  coefficients <- .bt_regularized_coefficient_names(prior, parameter_name)
  auxiliary <- .JAGS_bridgesampling_posterior_info(.bt_regularized_primitives(prior))
  parameters <- c(coefficients, as.character(auxiliary))
  attr(parameters, "lb") <- stats::setNames(c(rep(-Inf, length(coefficients)), attr(auxiliary, "lb")), parameters)
  attr(parameters, "ub") <- stats::setNames(c(rep(Inf, length(coefficients)), attr(auxiliary, "ub")), parameters)
  parameters
}

.bt_regularized_bridge_row <- function(samples){

  if(is.list(samples)) samples <- unlist(samples, use.names = TRUE)
  .bt_JAGS_marglik_prior_rows(samples)
}

# Within a centered contrast block the coordinates are spherical conditional
# on their unit's scale. Treatment and level coordinates with historical
# metadata each own a separate local scale.
.bt_regularized_coefficient_scales <- function(prior, dimension, global, slab, lookup){

  spec <- prior$regularization
  units <- spec$unit_index
  if(is.null(units)) units <- seq_len(dimension)
  scales <- lapply(seq_len(max(units)), function(unit){
    local <- if(length(spec$local_names)){
      lookup(spec$local_names[[unit]])
    }else if(length(spec$variance_names)){
      sqrt(lookup(spec$variance_names[[unit]]))
    }else 1
    .bt_regularized_scale(prior, global, local, slab)
  })
  scales[units]
}

.bt_regularized_bridge_density <- function(prior, parameter_name){

  spec <- prior$regularization
  coefficients <- .bt_regularized_coefficient_names(prior, parameter_name)
  dependencies <- .bt_regularized_auxiliary_names(prior)
  primitive_density <- .bt_JAGS_marglik_compile_prior_rows_evaluator(.bt_regularized_primitives(prior))
  function(samples){
    if(!all(c(coefficients, dependencies) %in% colnames(samples))){
      .bt_JAGS_marglik_missing_columns("A regularized coefficient or its shrinkage scales are missing from 'samples'.")
    }
    out <- unname(primitive_density(samples))
    valid <- rep(TRUE, nrow(samples))
    for(name in dependencies){
      valid <- valid & is.finite(samples[, name]) & samples[, name] > 0
    }
    out[!valid] <- -Inf
    if(!any(valid)) return(out)

    column <- function(name) as.numeric(samples[valid, name, drop = TRUE])
    global <- if(is.null(spec$global_name)) spec$scale else column(spec$global_name)
    slab <- if(!is.null(spec$slab_name)) column(spec$slab_name)
    scales <- .bt_regularized_coefficient_scales(prior, length(coefficients), global, slab, column)
    for(j in seq_along(coefficients)){
      scale <- scales[[j]]
      invalid_scale <- which(!is.finite(scale) | scale <= 0)
      if(length(invalid_scale)){
        .prior_numerical_signal("bridge log density", spec$family, "positive finite",
          which(valid)[invalid_scale], "The conditional coefficient scale is unrepresentable", error = TRUE)
      }
      value <- column(coefficients[[j]])
      density <- if(identical(spec$family, "lasso") && !length(spec$variance_names)){
        -log(2) - log(scale) - abs(value) / scale
      }else{
        stats::dnorm(value, 0, scale, log = TRUE)
      }
      out[valid] <- out[valid] + density
    }
    out
  }
}
