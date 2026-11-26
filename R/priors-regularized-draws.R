.bt_has_regularized_term <- function(prior){

  is.prior.regularized_term(prior) ||
    (is.prior.mixture(prior) && any(vapply(prior, .bt_has_regularized_term, logical(1))))
}

.bt_regularized_auxiliary_names <- function(prior){

  if(is.prior.regularized_term(prior)){
    spec <- prior$regularization
    return(c(spec$global_name, spec$slab_name, spec$local_names, spec$variance_names))
  }
  if(is.prior.mixture(prior)){
    return(unique(unlist(lapply(prior, .bt_regularized_auxiliary_names), use.names = FALSE)))
  }
  character()
}

# The registry is local to one joint draw request. Store a full column per
# stochastic scale so that selecting or reordering terms preserves the law.
.bt_regularized_draw_columns <- function(prior, name, n, registry){

  if(is.prior.mixture(prior)){
    if(is.prior.spike_and_slab(prior)){
      inclusion <- rng(.get_spike_and_slab_inclusion(prior), n)
      gate <- stats::rbinom(n, 1, inclusion)
      draws <- cbind(gate, inclusion)
      colnames(draws) <- paste0(name, c("_indicator", "_inclusion"))
      parts <- .bt_regularized_draw_columns(.get_spike_and_slab_variable(prior),
        paste0(name, "_variable"), n, registry)
      draws <- cbind(draws, parts)
    }else{
      .check_mixture_prior_weights(attr(prior, "prior_weights"))
      gate <- sample.int(length(prior), n, replace = TRUE, prob = attr(prior, "prior_weights"))
      draws <- matrix(gate, ncol = 1L, dimnames = list(NULL, paste0(name, "_indicator")))
      for(i in seq_along(prior)){
        parts <- .bt_regularized_draw_columns(prior[[i]], paste0(name, "_component_", i), n, registry)
        draws <- cbind(draws, parts)
      }
    }
    node <- .bt_dnode_prior_mixture(name, prior)
    values <- .bt_deterministic_node_evaluate(node, .bt_deterministic_lookup(draws))
    return(cbind(values, draws))
  }
  if(!is.prior.regularized_term(prior)){
    return(.generate_prior_sample_matrix(stats::setNames(list(prior), name), n))
  }

  spec <- prior$regularization
  # A term sampled by itself still needs the pool's shared hyperpriors.
  primitives <- .bt_regularized_primitives(prior, include_shared = TRUE)
  for(key in names(primitives)){
    if(!exists(key, envir = registry, inherits = FALSE)){
      draws <- rng(primitives[[key]], n)
      failed <- which(!is.finite(draws) | draws <= 0)
      if(length(failed)){
        .prior_numerical_signal("scale sampling", spec$family, "positive finite", failed,
          paste0("The sampled scale '", key, "' is unrepresentable"), error = TRUE, rng = TRUE)
      }
      registry[[key]] <- draws
    }
  }
  coordinates <- .JAGS_prior_coefficient_names(name, prior)
  scales <- .bt_regularized_coefficient_scales(prior, length(coordinates),
    function(key) registry[[key]])
  out <- matrix(NA_real_, nrow = n, ncol = length(coordinates), dimnames = list(NULL, coordinates))
  for(j in seq_along(coordinates)){
    scale <- scales[[j]]
    failed <- which(!is.finite(scale) | scale <= 0)
    if(length(failed)){
      .prior_numerical_signal("sampling", spec$family, "positive finite", failed,
        "The conditional coefficient scale is unrepresentable", error = TRUE, rng = TRUE)
    }
    out[, j] <- if(identical(spec$family, "lasso") && !length(spec$variance_names)){
      z <- extraDistr::rlaplace(n, mu = 0, sigma = scale)
      if(!is.null(spec$contrast_ratio)) z <- z + spec$contrast_ratio * extraDistr::rlaplace(n, mu = 0, sigma = scale)
      z
    }else{
      stats::rnorm(n, 0, scale)
    }
  }
  .prior_numerical_finite(out, "joint regularized prior")
  out
}
