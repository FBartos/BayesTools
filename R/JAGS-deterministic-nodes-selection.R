# Publication-weight deterministic node family.


# Publication weights ("omega") ------------------------------------------------
#
# A weight-function prior (alone, as the selection of a composed bias prior,
# or as branches of a publication-bias mixture) defines the publication
# weights 'omega' on the one-sided global p-value grid. One weight-function
# component with J local bins defines its local weights
#   cumulative, J > 2:  std_eta[j] <- eta[j] / sum(eta)
#                       omega_local[1] <- 1; omega_local[j] <- sum(std_eta[j:J])
#   cumulative, J = 2:  omega_local[1] <- 1; omega_local[2] <- omega_ratio
#   independent:        omega_local[1] <- 1, and omega_local[j] sampled on the
#                       omega scale or omega_local[j] <- exp(log_omega[j])
#   fixed:              omega_local[j] <- w[j]
# and maps them onto the global bins, 'omega[i] <- omega_local[index[i]]',
# which mirrors the bins of a two-sided weight function. A mixture defines
#   omega[i] <- sum_k omega_component_k[i] * equals(bias_indicator, k)
# with 'omega_component_k[i] <- 1' for branches without a selection.
# selection_backend_spec() writes this syntax into the model from the node
# specification (.bt_dnode_omega_spec(), .bt_dnode_omega_emit_branch(),
# .bt_dnode_omega_emit_composition()), placing each branch's weights next to
# its p-hacking syntax and the mixture composition among the transforms.
#
# The free coordinates of a single weight function are the ones bridge
# sampling uses: 'eta[j]', the monitored 'omega[2]' of a binary cumulative
# weight function (a copy of 'omega_ratio'), 'omega[j]' or 'log_omega[j]'
# (j > 1) of independent weights. Mixture components are read from their own
# nodes ('eta_component_k[j]', 'omega_ratio_component_k',
# 'omega_local_component_k[j]', 'log_omega_component_k[j]'), which models
# monitor only on request.

# JAGS syntax of one weight-function component, stochastic auxiliaries
# included. The model syntax and the initial values must agree on which node
# is stochastic (see .weightfunction_component_node_names()).
.bt_dnode_omega_component_syntax <- function(prior, component_id = NULL,
                                             global_cuts = NULL,
                                             force_one_sided = FALSE){

  node_names <- .weightfunction_component_node_names(
    prior           = prior,
    component_id    = component_id,
    global_cuts     = global_cuts,
    force_one_sided = force_one_sided
  )
  J <- node_names$n_bins
  syntax <- character()

  expansion     <- node_names$expansion
  all_cuts      <- node_names$all_cuts
  needs_mapping <- node_names$needs_mapping

  omega_local  <- node_names$omega_local
  omega_target <- node_names$omega_target

  if(prior$weights$type == "cumulative"){
    if(J == 2L){
      omega_ratio_name <- node_names$omega_ratio
      beta_parameters <- .weightfunction_alpha_marginal(prior$weights$alpha, 2L)
      syntax <- paste0(syntax,
        omega_ratio_name, " ~ dbeta(", beta_parameters$alpha, ", ", beta_parameters$beta, ")\n",
        omega_local, "[1] <- 1\n",
        omega_local, "[2] <- ", omega_ratio_name, "\n"
      )
    }else{
      eta_name <- node_names$eta
      std_eta_name <- node_names$std_eta

      for(i in seq_len(J)){
        syntax <- paste0(syntax, eta_name, "[", i, "] ~ dgamma(", prior$weights$alpha[i], ", 1)\n")
      }
      syntax <- paste0(syntax,
        "for(j in 1:", J, "){\n",
        "  ", std_eta_name, "[j] <- ", eta_name, "[j] / sum(", eta_name, ")\n",
        "}\n",
        omega_local, "[1] <- 1\n",
        "for(j in 2:", J, "){\n",
        "  ", omega_local, "[j] <- sum(", std_eta_name, "[j:", J, "])\n",
        "}\n"
      )
    }

  }else if(prior$weights$type == "fixed"){
    for(i in seq_len(J)){
      syntax <- paste0(syntax, omega_local, "[", i, "] <- ", prior$weights$omega[i], "\n")
    }

  }else if(prior$weights$type == "independent"){
    syntax <- paste0(syntax, omega_local, "[1] <- 1\n")
    if(J > 1L){
      for(i in 2:J){
        if(prior$weights$scale == "omega"){
          syntax <- paste0(syntax, .JAGS_prior.simple(prior$weights$prior, paste0(omega_local, "[", i, "]")))
        }else if(prior$weights$scale == "log_omega"){
          log_omega_name <- node_names$log_omega
          syntax <- paste0(
            syntax,
            .JAGS_prior.simple(prior$weights$prior, paste0(log_omega_name, "[", i, "]")),
            omega_local, "[", i, "] <- exp(", log_omega_name, "[", i, "])\n"
          )
        }
      }
    }
  }

  if(!is.null(component_id) || needs_mapping){
    global_bin_indices <- .weightfunction_global_bin_indices(all_cuts, expansion)
    for(i in seq_len(length(all_cuts) - 1L)){
      ind <- global_bin_indices[i]
      syntax <- paste0(syntax, omega_target, "[", i, "] <- ", omega_local, "[", expansion$index[ind], "]\n")
    }
  }

  syntax
}

# JAGS syntax of a mixture branch without a selection.
.bt_dnode_omega_none_component_syntax <- function(component_id, n_bins){

  syntax <- character()
  omega_target <- if(is.null(component_id)) "omega" else paste0("omega_component_", component_id)
  for(i in seq_len(n_bins)){
    syntax <- paste0(syntax, omega_target, "[", i, "] <- 1\n")
  }
  syntax
}

# JAGS syntax of the mixture composition of the global weights.
.bt_dnode_omega_mixture_syntax <- function(target, n_bins, indicator_terms){

  vapply(seq_len(n_bins), function(j){
    paste0(
      target, "[", j, "] <- ",
      paste0(
        "omega_component_", seq_along(indicator_terms), "[", j, "] * ",
        indicator_terms,
        collapse = " + "
      )
    )
  }, character(1))
}

# Free coordinates of one weight-function component: the single-weight-function
# coordinates of bridge sampling, or the component's own nodes in a mixture.
.bt_dnode_omega_free_names <- function(prior, component_id = NULL){

  J <- .weightfunction_n_bins(prior)
  suffix <- if(is.null(component_id)) "" else paste0("_component_", component_id)
  type <- prior$weights$type
  if(identical(type, "fixed")){
    return(character())
  }
  if(identical(type, "cumulative")){
    if(J == 2L){
      return(if(is.null(component_id)) "omega[2]" else paste0("omega_ratio", suffix))
    }
    return(paste0("eta", suffix, "[", seq_len(J), "]"))
  }
  if(J < 2L){
    return(character())
  }
  if(identical(prior$weights$scale, "omega")){
    return(paste0(if(is.null(component_id)) "omega" else paste0("omega_local", suffix), "[", 2:J, "]"))
  }

  paste0("log_omega", suffix, "[", 2:J, "]")
}

# Local weights (draws x J) of one weight-function component from its free
# coordinates, with the arithmetic of the marginal-likelihood parameters; NULL
# when a free coordinate is unavailable. Out-of-support auxiliaries signal the
# classed marginal-likelihood support condition. A node passes the component's
# free coordinates and number of bins it stores.
.bt_dnode_omega_local_values <- function(prior, lookup, component_id = NULL,
                                         free_names = .bt_dnode_omega_free_names(prior, component_id),
                                         J = .weightfunction_n_bins(prior)){

  n <- lookup$n
  type <- prior$weights$type
  if(identical(type, "fixed")){
    return(matrix(unname(prior$weights$omega), nrow = n, ncol = J, byrow = TRUE))
  }

  omega <- matrix(1, nrow = n, ncol = J)
  if(length(free_names) == 0L){
    return(omega)
  }
  free <- .bt_deterministic_lookup_values(lookup, free_names)
  if(is.null(free)){
    return(NULL)
  }

  if(identical(type, "cumulative") && J == 2L){
    invalid <- !is.finite(free[, 1L]) | free[, 1L] < 0 | free[, 1L] > 1
    if(any(invalid)){
      .bt_JAGS_marglik_out_of_support(
        "Bridge samples contain out-of-support binary cumulative weightfunction coordinate '",
        free_names,
        "'."
      )
    }
    omega[, 2L] <- free[, 1L]
    return(omega)
  }
  if(identical(type, "cumulative")){
    invalid <- !is.finite(free) | free <= 0
    if(any(invalid)){
      .bt_JAGS_marglik_out_of_support(
        "Bridge samples contain out-of-support positive auxiliary coordinate '",
        free_names[col(free)[which(invalid)[1L]]],
        "'."
      )
    }
    for(draw in seq_len(n)){
      eta <- free[draw, ]
      total <- sum(eta)
      if(!is.finite(total) || total <= 0){
        .prior_numerical_signal("cumulative weight normalization", "gamma", "natural", draw,
          paste0("The finite positive Gamma auxiliary total for '",
            paste(free_names, collapse = ", "), "' is not finite and positive"), error = TRUE)
      }
      std_eta <- eta / total
      omega[draw, ] <- rev(cumsum(rev(std_eta)))
    }
    return(omega)
  }
  if(identical(prior$weights$scale, "omega")){
    omega[, 2:J] <- free
  }else{
    omega[, 2:J] <- exp(free)
  }

  omega
}

# The definition of the weights of the selection priors 'priors' (a prior, a
# prior mixture, or a list of priors; see .selection_normalize_priors()): the
# branches, the global one-sided p-value grid (the union of the step-selection
# cuts, or the validated 'global_breaks' that contain them; c(0, 1) without a
# step selection), each branch's local bin of every global bin, and the name
# of the weight node. selection_backend_spec() and the registered 'omega' node
# are both built from it.
.bt_dnode_omega_spec <- function(priors, global_breaks = NULL, name = "omega"){

  branches <- lapply(.selection_normalize_priors(priors), .selection_branch_info)
  has_selection <- vapply(branches, function(branch) !is.null(branch$selection), logical(1))
  selections <- lapply(branches[has_selection], function(branch) branch$selection)
  if(is.null(global_breaks)){
    global_cuts <- if(length(selections) > 0L){
      weightfunctions_mapping(selections, cuts_only = TRUE, one_sided = TRUE)
    }else{
      c(0, 1)
    }
  }else{
    global_cuts <- .selection_validate_global_breaks(global_breaks)
    if(length(selections) > 0L){
      required_cuts <- weightfunctions_mapping(selections, cuts_only = TRUE, one_sided = TRUE)
      if(!all(vapply(required_cuts, function(x) any(abs(x - global_cuts) < sqrt(.Machine$double.eps)), logical(1)))){
        stop("'global_breaks' must contain all step-selection p-value breaks.", call. = FALSE)
      }
    }
  }

  global_index <- lapply(branches, function(branch){
    if(is.null(branch$selection)){
      return(NULL)
    }
    expansion <- .weightfunction_mapping_expansion(branch$selection, force_one_sided = TRUE)
    expansion$index[.weightfunction_global_bin_indices(global_cuts, expansion)]
  })

  list(
    name = name,
    branches = branches,
    has_selection = has_selection,
    uses_indicator = length(branches) > 1L,
    global_cuts = global_cuts,
    global_index = global_index,
    n_bins = length(global_cuts) - 1L
  )
}

# The registered 'omega' node of a prior with a step selection in some branch;
# NULL otherwise (the unit weights of p-hacking-only priors are not a node).
.bt_dnode_omega <- function(parameter, prior, spec = .bt_dnode_omega_spec(prior)){

  if(!any(spec$has_selection)){
    return(NULL)
  }
  branches <- spec$branches
  # The free coordinates and local bins of each branch, stored for the
  # evaluator, which one caller may run on thousands of single draws.
  spec$free_names <- lapply(seq_along(branches), function(k){
    if(is.null(branches[[k]]$selection)){
      return(character())
    }
    .bt_dnode_omega_free_names(
      branches[[k]]$selection,
      component_id = if(spec$uses_indicator) k else NULL
    )
  })
  spec$local_bins <- vapply(branches, function(branch){
    if(is.null(branch$selection)) 0L else .weightfunction_n_bins(branch$selection)
  }, integer(1))

  .bt_deterministic_node(
    family = "omega",
    node = spec$name,
    coordinates = if(spec$n_bins == 1L){
      spec$name
    }else{
      paste0(spec$name, "[", seq_len(spec$n_bins), "]")
    },
    dependencies = c(
      if(spec$uses_indicator) "bias_indicator",
      unlist(spec$free_names, use.names = FALSE)
    ),
    parameter = parameter,
    spec = spec
  )
}

# JAGS syntax of branch k of the weights: the local weights of its step
# selection mapped onto the global bins, or unit weights for a mixture branch
# (or a p-hacking-only prior) without a selection. A single branch writes the
# weight node directly.
.bt_dnode_omega_emit_branch <- function(spec, k){

  selection <- spec$branches[[k]]$selection
  if(spec$uses_indicator){
    if(is.null(selection)){
      return(.bt_dnode_omega_none_component_syntax(component_id = k, n_bins = spec$n_bins))
    }
    return(.bt_dnode_omega_component_syntax(
      prior           = selection,
      component_id    = k,
      global_cuts     = spec$global_cuts,
      force_one_sided = TRUE
    ))
  }

  syntax <- if(!is.null(selection)){
    .bt_dnode_omega_component_syntax(
      prior           = selection,
      component_id    = NULL,
      global_cuts     = spec$global_cuts,
      force_one_sided = TRUE
    )
  }else if(!is.null(spec$branches[[k]]$phacking)){
    .bt_dnode_omega_none_component_syntax(component_id = NULL, n_bins = spec$n_bins)
  }else{
    character()
  }

  .selection_rename_jags_node(syntax, "omega", spec$name)
}

# JAGS syntax of the mixture composition of the weights (none for a single
# branch).
.bt_dnode_omega_emit_composition <- function(spec){

  if(!spec$uses_indicator){
    return(character())
  }

  .bt_dnode_omega_mixture_syntax(
    target = spec$name,
    n_bins = spec$n_bins,
    indicator_terms = paste0("equals(bias_indicator, ", seq_along(spec$branches), ")")
  )
}

.bt_dnode_omega_emit <- function(node){

  spec <- node$spec
  pieces <- unlist(lapply(seq_along(spec$branches), function(k){
    .bt_dnode_omega_emit_branch(spec, k)
  }), use.names = FALSE)
  pieces <- sub("[\r\n]+$", "", pieces[nzchar(pieces)])

  c(pieces, .bt_dnode_omega_emit_composition(spec))
}

.bt_dnode_omega_evaluate <- function(node, lookup){

  spec <- node$spec
  branch_values <- function(k, branch_lookup){
    selection <- spec$branches[[k]]$selection
    if(is.null(selection)){
      return(matrix(1, nrow = branch_lookup$n, ncol = spec$n_bins))
    }
    local <- .bt_dnode_omega_local_values(
      prior = selection,
      lookup = branch_lookup,
      component_id = if(spec$uses_indicator) k else NULL,
      free_names = spec$free_names[[k]],
      J = spec$local_bins[[k]]
    )
    if(is.null(local)){
      return(NULL)
    }
    local[, spec$global_index[[k]], drop = FALSE]
  }

  if(!spec$uses_indicator){
    return(branch_values(1L, lookup))
  }

  # A mixture's weights are those of the active branch (the other terms of
  # the JAGS sum are exact zeros).
  indicator <- .bt_deterministic_lookup_value(lookup, "bias_indicator")
  if(is.null(indicator)){
    return(NULL)
  }
  out <- matrix(NA_real_, nrow = lookup$n, ncol = spec$n_bins)
  for(k in unique(indicator)){
    rows <- which(indicator == k)
    if(!k %in% seq_along(spec$branches)){
      stop("Bias-mixture indicator draws must index a mixture branch.", call. = FALSE)
    }
    branch_lookup <- .bt_deterministic_lookup(
      lookup$draws[rows, , drop = FALSE],
      lookup$prior_list
    )
    values <- branch_values(k, branch_lookup)
    if(is.null(values)){
      return(NULL)
    }
    out[rows, ] <- values
  }

  out
}

# The 'omega' nodes of the priors most recently evaluated one draw at a time:
# JAGS_marglik_parameters() is called once per draw with the same prior list
# (bridge sampling rows, RoBMA's IWMDE rows), so each node is built once. A
# node is a function of its parameter name and prior alone; entries are found
# with identical(), which returns immediately for the same prior object.
.bt_dnode_omega_cache <- new.env(parent = emptyenv())
.bt_dnode_omega_cache_size <- 16L

.bt_dnode_omega_cached <- function(parameter, prior){

  entries <- .bt_dnode_omega_cache$entries
  for(entry in entries){
    if(identical(entry$prior, prior) && identical(entry$parameter, parameter)){
      return(entry$node)
    }
  }
  node <- .bt_dnode_omega(parameter, prior)
  entries <- c(list(list(parameter = parameter, prior = prior, node = node)), entries)
  .bt_dnode_omega_cache$entries <- entries[seq_len(min(length(entries), .bt_dnode_omega_cache_size))]

  node
}

# The weights of a single weight-function prior on its own global bins (the
# 'omega' node of 'node', compiled once by the caller or cached), for one draw
# of the marginal-likelihood parameters and bridge sampling; unavailable free
# coordinates stop with the bridge-sampling message.
.bt_dnode_omega_prior_values <- function(prior, samples,
                                         node = .bt_dnode_omega_cached("omega", prior)){

  values <- .bt_deterministic_node_evaluate(
    node,
    .bt_deterministic_row_lookup(samples)
  )
  if(is.null(values)){
    .bt_dnode_omega_missing_stop(prior)
  }

  as.vector(values)
}

.bt_dnode_omega_missing_stop <- function(prior){

  J <- .weightfunction_n_bins(prior)
  if(identical(prior$weights$type, "cumulative") && J == 2L){
    .bt_JAGS_marglik_missing_columns(
      "'samples' does not contain the monitored binary cumulative weightfunction parameter."
    )
  }
  if(identical(prior$weights$type, "cumulative")){
    .bt_JAGS_marglik_missing_columns(
      "'samples' does not contain all monitored cumulative weightfunction parameters."
    )
  }

  .bt_JAGS_marglik_missing_columns(
    "'samples' does not contain all monitored independent weightfunction parameters."
  )
}
