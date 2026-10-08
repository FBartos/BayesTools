# Prior-mixture deterministic node family.


# Prior mixtures ("prior_mixture") --------------------------------------------
#
# Spike-and-slab and mixture priors define their parameter as a deterministic
# node of the mixture components and the component indicator:
#   spike and slab:  p = p_variable * p_indicator
#   mixture:         p = p_component_1 * (p_indicator == 1) + ...
#                        + p_component_M * (p_indicator == M)
# and a publication-bias mixture defines its PET and PEESE terms as
#   PET <- PET_1 * equals(bias_indicator, k)
# for the branch k that holds them. Vector (factor) parameters apply the same
# definition to every coefficient. The R evaluator takes, per draw, the value
# of the active component (the other terms of the JAGS sum are exact zeros),
# so only the components that are active in some draw need to be available.
# Point components are constants and need no draws.

.bt_dnode_prior_mixture <- function(parameter, prior){

  if(is.prior.spike_and_slab(prior)){
    variable <- .get_spike_and_slab_variable(prior)
    components <- list(
      .bt_dnode_prior_mixture_component(paste0(parameter, "_variable"), variable)
    )
    coordinates <- .bt_dnode_prior_mixture_names(parameter, variable)
    return(.bt_deterministic_node(
      family = "prior_mixture",
      node = parameter,
      coordinates = coordinates,
      dependencies = c(components[[1L]]$dependencies, paste0(parameter, "_indicator")),
      parameter = parameter,
      spec = list(
        kind = "spike_and_slab",
        indicator = paste0(parameter, "_indicator"),
        components = components
      )
    ))
  }

  components <- lapply(seq_along(prior), function(k){
    .bt_dnode_prior_mixture_component(
      paste0(parameter, "_component_", k),
      prior[[k]]
    )
  })
  coordinates <- .bt_dnode_prior_mixture_names(parameter, prior[[1L]])
  .bt_deterministic_node(
    family = "prior_mixture",
    node = parameter,
    coordinates = coordinates,
    dependencies = c(
      paste0(parameter, "_indicator"),
      unlist(lapply(components, `[[`, "dependencies"), use.names = FALSE)
    ),
    parameter = parameter,
    spec = list(
      kind = "mixture",
      indicator = paste0(parameter, "_indicator"),
      components = components
    )
  )
}

# The PET or PEESE term of a publication-bias mixture.
.bt_dnode_prior_mixture_bias_term <- function(parameter, prior, term){

  is_term <- vapply(prior, if(identical(term, "PET")) is.prior.PET else is.prior.PEESE, logical(1))
  if(!any(is_term)){
    return(NULL)
  }
  source <- paste0(term, "_1")

  .bt_deterministic_node(
    family = "prior_mixture",
    node = term,
    coordinates = term,
    dependencies = c(source, "bias_indicator"),
    parameter = parameter,
    spec = list(
      kind = "bias_term",
      indicator = "bias_indicator",
      branch = which(is_term),
      source = source
    )
  )
}

.bt_dnode_prior_mixture_component <- function(name, prior){

  coordinates <- .bt_dnode_prior_mixture_names(name, prior)
  location <- if(is.prior.point(prior) && !.is_prior_expression(prior)){
    prior$parameters[["location"]]
  }else{
    NULL
  }

  list(
    name = name,
    coordinates = coordinates,
    location = location,
    dependencies = if(is.null(location)) coordinates else character()
  )
}

# Coordinate names of a (component) prior: a scalar, or one per coefficient of
# a factor prior.
.bt_dnode_prior_mixture_names <- function(name, prior){

  if(is.prior.factor(prior)){
    return(.JAGS_prior_factor_names(name, prior))
  }

  name
}

.bt_dnode_prior_mixture_emit <- function(node){

  spec <- node$spec
  if(identical(spec$kind, "ordered_spike_and_slab")){
    if(length(node$coordinates) == 1L){
      return(paste0(node$node, " = ", node$node, "_variable * ", spec$indicator))
    }
    component <- spec$components[[1L]]
    return(paste0(node$coordinates, " <- ", component$coordinates, " * ", spec$indicator))
  }
  if(identical(spec$kind, "spike_and_slab")){
    return(paste0(node$node, " = ", node$node, "_variable * ", spec$indicator))
  }
  if(identical(spec$kind, "bias_term")){
    return(paste0(node$node, " <- ", spec$source, " * equals(", spec$indicator, ", ", spec$branch, ")"))
  }

  paste0(
    node$node, " = ",
    paste0(
      vapply(spec$components, `[[`, character(1), "name"),
      " * (", spec$indicator, " == ", seq_along(spec$components), ")",
      collapse = " + "
    )
  )
}

.bt_dnode_prior_mixture_evaluate <- function(node, lookup){

  spec <- node$spec
  n <- lookup$n
  indicator <- .bt_deterministic_lookup_value(lookup, spec$indicator)
  if(is.null(indicator)){
    return(NULL)
  }

  if(identical(spec$kind, "bias_term")){
    values <- rep(0, n)
    active <- indicator == spec$branch
    if(any(active)){
      source <- .bt_deterministic_lookup_value(lookup, spec$source)
      if(is.null(source)){
        return(NULL)
      }
      values[active] <- source[active] * 1
    }
    return(values)
  }

  component_values <- function(component){
    if(!is.null(component$location)){
      return(matrix(component$location, nrow = n, ncol = length(component$coordinates)))
    }
    .bt_deterministic_lookup_values(lookup, component$coordinates)
  }

  if(spec$kind %in% c("spike_and_slab", "ordered_spike_and_slab")){
    if(any(!is.finite(indicator)) || any(!indicator %in% c(0, 1))){
      if(identical(spec$kind, "ordered_spike_and_slab")){
        .bt_ordered_stop("Ordered inclusion indicator draws must be zero or one.", "BayesTools_ordered_invalid_state")
      }
      stop("Inclusion indicator draws of '", node$node, "' must be zero or one.", call. = FALSE)
    }
    out <- matrix(0,n,length(node$coordinates))
    active <- indicator==1
    if(any(active)){
      variable <- component_values(spec$components[[1L]])
      if(is.null(variable)) return(NULL)
      out[active,] <- variable[active,,drop=FALSE]
    }
    return(out)
  }

  out <- matrix(NA_real_, nrow = n, ncol = length(node$coordinates))
  if(identical(spec$kind, "ordered_mixture") &&
     any(!is.finite(indicator) | !indicator %in% seq_along(spec$components))){
    .bt_ordered_stop(paste0("Ordered total indicator '", spec$indicator,
      "' does not select a declared component."), "BayesTools_ordered_invalid_state")
  }
  for(k in unique(indicator)){
    if(!k %in% seq_along(spec$components)){
      stop(
        "Mixture indicator draws of '", node$node,
        "' must index a mixture component.",
        call. = FALSE
      )
    }
    values <- component_values(spec$components[[k]])
    if(is.null(values)){
      return(NULL)
    }
    rows <- indicator == k
    out[rows, ] <- values[rows, , drop = FALSE]
  }

  out
}

# The value of a spike-and-slab or mixture parameter in one draw of the
# marginal-likelihood parameters.
.bt_dnode_prior_mixture_parameter_values <- function(samples, prior,
                                                     parameter_name){

  node <- .bt_dnode_prior_mixture(parameter_name, prior)
  values <- .bt_deterministic_node_evaluate(
    node,
    .bt_deterministic_row_lookup(samples)
  )
  if(is.null(values)){
    .bt_JAGS_marglik_missing_columns(paste0(
      "'samples' does not contain all monitored ",
      if(is.prior.spike_and_slab(prior)) "spike-and-slab" else "prior mixture",
      " parameters of '", parameter_name, "'."
    ))
  }

  as.vector(values)
}
