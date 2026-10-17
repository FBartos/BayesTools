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
      branch_count = length(prior),
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

# Bias indicator draws of a bias mixture are finite integers that index a
# declared branch, for the registered bias-term nodes and the bias-mixture
# parameters of the marginal likelihood alike.
.bt_dnode_bias_indicator_check <- function(indicator, branch_count, name){

  if(any(!is.finite(indicator))){
    stop("Bias indicator draws of '", name, "' must be finite.", call. = FALSE)
  }
  if(any(!indicator %in% seq_len(branch_count))){
    stop("Bias indicator draws of '", name, "' must index a declared bias branch.", call. = FALSE)
  }

  invisible(TRUE)
}

.bt_dnode_prior_mixture_evaluate <- function(node, lookup){

  spec <- node$spec
  n <- lookup$n
  indicator <- .bt_deterministic_lookup_value(lookup, spec$indicator)
  if(is.null(indicator)){
    return(NULL)
  }

  if(identical(spec$kind, "bias_term")){
    if(!is.numeric(spec$branch_count) || length(spec$branch_count) != 1L ||
       !is.finite(spec$branch_count) || spec$branch_count < 1 ||
       spec$branch_count != floor(spec$branch_count) ||
       !all(spec$branch %in% seq_len(spec$branch_count))){
      .bt_stop_refit_required("The fitted bias-term branch declaration is missing or inconsistent. Refit the model with this version of BayesTools.")
    }
    .bt_dnode_bias_indicator_check(indicator, spec$branch_count, node$node)
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
  if(any(!is.finite(indicator)) || any(!indicator %in% seq_along(spec$components))){
    stop("Mixture indicator draws of '", node$node,
      "' must index a mixture component.", call. = FALSE)
  }
  for(k in unique(indicator)){
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


# Ordinary vector priors ("prior_dirichlet", "prior_mt") ------------------------
#
# The Dirichlet and multivariate-t priors define their coordinates
# deterministically from sampled primitives (see .JAGS_prior.vector()):
#   Dirichlet:        theta[i] <- eta[i] / sum(eta[1:K]),  eta[i] ~ dgamma(alpha[i], 1)
#   multivariate t:   theta[i] <- z[i] / sqrt(s) + location,
#                     z ~ dmnorm(0, precision), s ~ dgamma(df / 2, df / 2)
# Orthonormal and meandif multivariate-t factors, and the components of
# mixtures and the slabs of spike-and-slabs that are such priors, share the
# syntax and therefore the nodes. The dependencies are the primitive
# coordinates 'eta', 'z' and 's'; the location is a constant of the prior. The
# Dirichlet weights of a variance allocation are coordinates of the allocation
# family (the SDs read their weights or gamma auxiliaries) and no node of this
# family. A multivariate-t location that is an expression is part of the
# syntax; its values are unavailable to the R evaluator.

# JAGS names a vector of length one by its scalar name.
.bt_dnode_prior_vector_names <- function(name, K){

  if(K == 1L) name else paste0(name, "[", seq_len(K), "]")
}

.bt_dnode_prior_dirichlet <- function(name, prior, parameter = name){

  K   <- prior$parameters[["K"]]
  eta <- .JAGS_prior_dirichlet_eta_name(name)

  .bt_deterministic_node(
    family = "prior_dirichlet",
    node = name,
    coordinates = .bt_dnode_prior_vector_names(name, K),
    dependencies = .bt_dnode_prior_vector_names(eta, K),
    parameter = parameter,
    spec = list(K = K, eta = eta)
  )
}

.bt_dnode_prior_mt <- function(name, prior, parameter = name){

  if(.is_prior_expression(prior)){
    prior <- .prior_expression_to_character(prior)
  }
  K <- if(is.prior.factor(prior)) .get_prior_factor_levels(prior) else prior$parameters[["K"]]
  z <- paste0("prior_par_z_", name)
  s <- paste0("prior_par_s_", name)

  .bt_deterministic_node(
    family = "prior_mt",
    node = name,
    coordinates = .bt_dnode_prior_vector_names(name, K),
    dependencies = c(.bt_dnode_prior_vector_names(z, K), s),
    parameter = parameter,
    spec = list(K = K, z = z, s = s, location = prior$parameters[["location"]])
  )
}

# The nodes of the ordinary vector priors of a (component of a) prior: the
# prior itself, the components of a mixture and the slab of a spike-and-slab
# (named as in .JAGS_prior.mixture() and .JAGS_prior.spike_and_slab()).
.bt_dnode_prior_vector_nodes <- function(name, prior, parameter = name){

  if(is.prior.spike_and_slab(prior)){
    return(.bt_dnode_prior_vector_nodes(
      paste0(name, "_variable"), .get_spike_and_slab_variable(prior), parameter
    ))
  }
  if(inherits(prior, "prior.bias_mixture")){
    return(list())
  }
  if(is.prior.mixture(prior)){
    return(unlist(lapply(seq_along(prior), function(k){
      .bt_dnode_prior_vector_nodes(paste0(name, "_component_", k), prior[[k]], parameter)
    }), recursive = FALSE))
  }
  if(is.prior.ordered(prior) || !is.prior.vector(prior)){
    return(list())
  }

  if(identical(prior[["distribution"]], "dirichlet") &&
     is.null(attr(prior, "random_allocation", exact = TRUE))){
    return(list(.bt_dnode_prior_dirichlet(name, prior, parameter)))
  }
  if(identical(prior[["distribution"]], "mt")){
    return(list(.bt_dnode_prior_mt(name, prior, parameter)))
  }

  list()
}

.bt_dnode_prior_dirichlet_emit <- function(node){

  paste0(node$node, "[", seq_len(node$spec$K), "] <- ", node$spec$eta, "[", seq_len(node$spec$K), "] / sum(",
         node$spec$eta, "[1:", node$spec$K, "])")
}

.bt_dnode_prior_mt_emit <- function(node){

  spec <- node$spec
  c(
    paste0("for(i in 1:", spec$K, "){"),
    paste0("  ", node$node, "[i] <- ", spec$z, "[i]/sqrt(", spec$s, ") + ",
           .JAGS_prior_number(spec$location), " "),
    "}"
  )
}

# The Dirichlet weights are the gamma auxiliaries over their sum, accumulated
# in the order of the JAGS 'sum()'.
.bt_dnode_prior_dirichlet_evaluate <- function(node, lookup){

  eta <- .bt_deterministic_lookup_values(lookup, node$dependencies)
  if(is.null(eta)){
    return(NULL)
  }

  total <- eta[, 1L]
  for(j in seq_len(node$spec$K)[-1L]){
    total <- total + eta[, j]
  }
  if(any(!is.finite(eta)) || any(eta < 0) || any(!is.finite(total)) || any(total <= 0)){
    stop(
      "Dirichlet gamma draws of '", node$node,
      "' must be finite and nonnegative with finite positive row sums.",
      call. = FALSE
    )
  }

  eta / total
}

# The coordinates are the normal primitives over the square root of the scale
# primitive plus the location.
.bt_dnode_prior_mt_evaluate <- function(node, lookup){

  spec <- node$spec
  if(!is.numeric(spec$location)){
    return(NULL)
  }
  z <- .bt_deterministic_lookup_values(lookup, .bt_dnode_prior_vector_names(spec$z, spec$K))
  s <- .bt_deterministic_lookup_value(lookup, spec$s)
  if(is.null(z) || is.null(s)){
    return(NULL)
  }
  if(any(!is.finite(z)) || any(!is.finite(s)) || any(s <= 0)){
    stop(
      "Multivariate-t primitive draws of '", node$node,
      "' must be finite, with positive scale draws '", spec$s, "'.",
      call. = FALSE
    )
  }

  z / sqrt(s) + spec$location
}
