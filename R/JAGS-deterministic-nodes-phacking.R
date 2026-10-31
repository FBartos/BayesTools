# P-hacking deterministic node family.


# P-hacking nodes ("phacking") ---------------------------------------------------
#
# A p-hacking branch (a prior_phacking(), the p-hacking of a prior_bias(), or a
# branch of a publication-bias mixture) with the severity alpha defines, from
# the constants of its form (phack_backend_constants()),
#   phack_kind <- k                          the form code
#   pi_null    <- alpha * a                  the null mass of the source bin
#   beta_null  <- alpha * b
#   phack_z_source[1:2], phack_z_dest[1:2]   the z geometry of the form
# A mixture writes these per branch ('<node>_component_k', zeros for the
# branches without p-hacking; the severities are the sampled
# 'alpha_component_k') and composes
#   <node> <- sum_k <node>_component_k * equals(bias_indicator, k)
# and the selection kernel reads the vector rule of its branches,
#   sel_vector_rule <- sum_k rule_k * equals(bias_indicator, k).
# Every one of these nodes is a registered node of the family with one R
# evaluator, .bt_dnode_phacking_quantity(): selection_backend_spec() writes the
# syntax from the node specification (.bt_dnode_phacking_spec()), and the
# bridge sampling and marginal-likelihood parameters of a draw evaluate the same
# quantities. The R evaluator takes, per draw, the value of the active branch
# (the other terms of the JAGS sums are exact zeros). 'pi_null' and 'beta_null'
# read the severity from the coordinate 'alpha' (sampled in a single branch,
# the composition in a mixture), so a fit that monitors 'alpha' rebuilds them
# without the component severities.

.bt_dnode_phacking_quantities <- c(
  "alpha", "phack_kind", "pi_null", "beta_null",
  "phack_z_source", "phack_z_dest", "sel_vector_rule"
)

.bt_dnode_phacking_vector_quantities <- c("phack_z_source", "phack_z_dest")

# The vector rule of a selection branch (0: none or the product rule; 1 and 2:
# the one- and two-sided vector rules).
.bt_dnode_phacking_vector_rule <- function(branch,
                                           model = if(!is.null(branch$selection)) selection_model_spec(branch$selection)){

  if(is.null(model) || identical(model[["weight_rule"]], "product")){
    return(0L)
  }

  if(identical(branch$selection[["side"]], "two-sided")) 2L else 1L
}

# The definition of the p-hacking nodes of the selection branches 'branches'
# (the branch information of .bt_dnode_omega_spec()): the public node names, the
# constants of each p-hacking branch, and the vector rules.
.bt_dnode_phacking_spec <- function(branches,
                                    names = .selection_backend_names(list(omega = "omega", alpha = "alpha")),
                                    vector_rule = NULL){

  uses_indicator <- length(branches) > 1L
  if(is.null(vector_rule)){
    vector_rule <- vapply(branches, .bt_dnode_phacking_vector_rule, integer(1))
  }

  list(
    names = names,
    branches = branches,
    has_phacking = vapply(branches, function(branch) !is.null(branch$phacking), logical(1)),
    uses_indicator = uses_indicator,
    indicator_terms = if(uses_indicator){
      paste0("equals(bias_indicator, ", seq_along(branches), ")")
    }else{
      rep("1", length(branches))
    },
    constants = lapply(branches, function(branch){
      phacking <- branch$phacking
      if(is.null(phacking)){
        return(NULL)
      }
      phack_backend_constants(
        phacking$form, phacking$source, phacking$destination, target = phacking$target
      )
    }),
    vector_rule = vector_rule
  )
}

# The JAGS node of a quantity inside branch k of a mixture, or the public node.
.bt_dnode_phacking_node_name <- function(spec, quantity, k = NULL){

  if(identical(quantity, "sel_vector_rule")){
    return(quantity)
  }

  .selection_phacking_node_name(quantity, if(spec$uses_indicator) k else NULL, spec$names)
}

.bt_dnode_phacking_coordinates <- function(spec, quantity){

  name <- .bt_dnode_phacking_node_name(spec, quantity)
  if(quantity %in% .bt_dnode_phacking_vector_quantities){
    return(paste0(name, "[", 1:2, "]"))
  }

  name
}

# The registered nodes of the p-hacking branches of a selection prior (a prior,
# a prior mixture, or a list of priors; see .selection_normalize_priors()): the
# quantities of .bt_dnode_phacking_quantities, without the composed severity
# 'alpha' of a single branch, which is the sampled severity itself. NULL
# without a p-hacking branch (the vector rule of a weight function alone is not
# a node).
.bt_dnode_phacking_nodes <- function(parameter, prior){

  branches <- lapply(.selection_normalize_priors(prior), .selection_branch_info)
  if(!any(vapply(branches, function(branch) !is.null(branch$phacking), logical(1)))){
    return(NULL)
  }
  spec <- .bt_dnode_phacking_spec(branches)

  quantities <- setdiff(
    .bt_dnode_phacking_quantities,
    if(!spec$uses_indicator) "alpha"
  )
  lapply(quantities, function(quantity){
    node_spec <- spec
    node_spec$quantity <- quantity
    .bt_deterministic_node(
      family = "phacking",
      node = .bt_dnode_phacking_node_name(spec, quantity),
      coordinates = .bt_dnode_phacking_coordinates(spec, quantity),
      dependencies = .bt_dnode_phacking_dependencies(spec, quantity),
      parameter = parameter,
      spec = node_spec
    )
  })
}

# The coordinates the evaluator of a quantity reads: the branch indicator of a
# mixture, and the severities of the p-hacking branches that are not point
# priors (the component severities for the composed 'alpha', the coordinate
# 'alpha' for the null masses).
.bt_dnode_phacking_dependencies <- function(spec, quantity){

  dependencies <- if(spec$uses_indicator) "bias_indicator" else character()
  if(quantity %in% c("alpha", "pi_null", "beta_null")){
    sampled <- which(spec$has_phacking)
    sampled <- sampled[!vapply(sampled, function(k){
      is.prior.point(spec$branches[[k]]$phacking$alpha)
    }, logical(1))]
    dependencies <- c(
      dependencies,
      if(length(sampled) > 0L){
        if(identical(quantity, "alpha")){
          vapply(sampled, function(k) .bt_dnode_phacking_node_name(spec, "alpha", k), character(1))
        }else{
          spec$names$alpha
        }
      }
    )
  }

  unique(dependencies)
}


# JAGS syntax --------------------------------------------------------------------

# The statements of branch k as a list by quantity (no trailing newline, except
# the prior statement of the severity): zeros for a branch without p-hacking in
# a mixture.
.bt_dnode_phacking_branch_lines <- function(spec, k){

  component_id <- if(spec$uses_indicator) k else NULL
  name <- function(quantity) .selection_phacking_node_name(quantity, component_id, spec$names)
  alpha_name <- name("alpha")
  prior <- spec$branches[[k]]$phacking

  if(is.null(prior)){
    return(list(
      alpha = paste0(alpha_name, " <- 0\n"),
      phack_kind = paste0(name("phack_kind"), " <- 0"),
      pi_null = paste0(name("pi_null"), " <- 0"),
      beta_null = paste0(name("beta_null"), " <- 0"),
      phack_z_source = paste0(name("phack_z_source"), "[", 1:2, "] <- 0"),
      phack_z_dest = paste0(name("phack_z_dest"), "[", 1:2, "] <- 0")
    ))
  }

  constants <- spec$constants[[k]]
  list(
    alpha = .JAGS_prior.simple(prior$alpha, alpha_name),
    phack_kind = paste0(name("phack_kind"), " <- ", constants$phack_kind),
    pi_null = paste0(
      name("pi_null"), " <- ", alpha_name, " * ",
      .selection_format_number(constants$pi_null_per_alpha)
    ),
    beta_null = paste0(
      name("beta_null"), " <- ", alpha_name, " * ",
      .selection_format_number(constants$beta_null_per_alpha)
    ),
    phack_z_source = paste0(
      name("phack_z_source"), "[", 1:2, "] <- ", .selection_format_number(constants$z_source)
    ),
    phack_z_dest = paste0(
      name("phack_z_dest"), "[", 1:2, "] <- ", .selection_format_number(constants$z_destination)
    )
  )
}

# JAGS syntax of branch k: its severity prior and its constants; nothing for a
# single branch without p-hacking.
.bt_dnode_phacking_emit_branch <- function(spec, k){

  if(!spec$uses_indicator && is.null(spec$branches[[k]]$phacking)){
    return(character())
  }

  lines <- .bt_dnode_phacking_branch_lines(spec, k)
  paste0(
    lines$alpha,
    paste0(unlist(lines[setdiff(names(lines), "alpha")], use.names = FALSE), "\n", collapse = "")
  )
}

# JAGS syntax of the vector rule of the selection kernel.
.bt_dnode_phacking_emit_vector_rule <- function(spec){

  paste0(
    "sel_vector_rule <- ",
    paste0(spec$vector_rule, " * ", spec$indicator_terms, collapse = " + ")
  )
}

# JAGS syntax of the mixture composition of the severity and the constants of
# the branches (none for a single branch), restricted to 'quantities'.
.bt_dnode_phacking_emit_composition <- function(spec, quantities = .bt_dnode_phacking_quantities){

  if(!spec$uses_indicator){
    return(character())
  }

  ids <- seq_along(spec$branches)
  scalar <- function(quantity){
    .selection_jags_active_scalar(
      paste0(quantity, "_component_"), spec$names[[quantity]], spec$indicator_terms, ids
    )
  }
  vector <- function(quantity, index){
    .selection_jags_active_vector(
      paste0(quantity, "_component_"), spec$names[[quantity]], index, spec$indicator_terms, ids
    )
  }

  syntax <- c(
    if("alpha" %in% quantities) scalar("alpha"),
    if("phack_kind" %in% quantities) scalar("phack_kind"),
    if("pi_null" %in% quantities) scalar("pi_null"),
    if("beta_null" %in% quantities) scalar("beta_null")
  )
  if(any(spec$has_phacking)){
    for(index in 1:2){
      syntax <- c(
        syntax,
        if("phack_z_source" %in% quantities) vector("phack_z_source", index),
        if("phack_z_dest" %in% quantities) vector("phack_z_dest", index)
      )
    }
  }

  syntax
}

# The statements that define one node: its branch statements and, in a mixture,
# its composition.
.bt_dnode_phacking_emit <- function(node){

  spec <- node$spec
  quantity <- spec$quantity
  if(identical(quantity, "sel_vector_rule")){
    return(.bt_dnode_phacking_emit_vector_rule(spec))
  }

  branch_lines <- if(identical(quantity, "alpha") && spec$uses_indicator){
    # the severities of the branches are prior statements
    character()
  }else{
    unlist(lapply(seq_along(spec$branches), function(k){
      .bt_dnode_phacking_branch_lines(spec, k)[[quantity]]
    }), use.names = FALSE)
  }

  c(branch_lines, .bt_dnode_phacking_emit_composition(spec, quantity))
}


# R evaluator -------------------------------------------------------------------

# The severity alpha of branch k on the rows of a lookup: the location of a point
# prior, or the draws of the sampled coordinate 'name' (checked against the
# support of an inverse-gamma prior); NULL when the coordinate is unavailable.
.bt_dnode_phacking_severity <- function(spec, k, lookup, name){

  alpha_prior <- spec$branches[[k]]$phacking$alpha
  if(is.prior.point(alpha_prior)){
    location <- .bt_formula_numeric_point(alpha_prior)
    if(is.null(location)){
      .bt_formula_point_stop("alpha", "missing_point_owner",
        "this generic prior reader has no expression replay owner")
    }
    return(matrix(location, nrow = lookup$n, ncol = 1L))
  }

  values <- .bt_deterministic_lookup_values(lookup, name)
  if(is.null(values)){
    return(NULL)
  }
  if(identical(alpha_prior[["distribution"]], "invgamma") &&
     any(!is.finite(values) | values <= 0)){
    .bt_JAGS_marglik_out_of_support(
      "Bridge samples contain out-of-support inverse-gamma coordinate 'alpha'."
    )
  }

  values
}

# The values of a quantity in branch k on the rows of a lookup (a matrix with
# one row per draw), or NULL when a coordinate is unavailable.
.bt_dnode_phacking_branch_values <- function(spec, quantity, k, lookup){

  n <- lookup$n
  vector_quantity <- quantity %in% .bt_dnode_phacking_vector_quantities
  if(identical(quantity, "sel_vector_rule")){
    return(matrix(as.numeric(spec$vector_rule[[k]]), nrow = n, ncol = 1L))
  }
  if(is.null(spec$branches[[k]]$phacking)){
    return(matrix(0, nrow = n, ncol = if(vector_quantity) 2L else 1L))
  }

  constants <- spec$constants[[k]]
  if(identical(quantity, "phack_kind")){
    return(matrix(as.numeric(constants$phack_kind), nrow = n, ncol = 1L))
  }
  if(identical(quantity, "phack_z_source")){
    return(matrix(constants$z_source, nrow = n, ncol = 2L, byrow = TRUE))
  }
  if(identical(quantity, "phack_z_dest")){
    return(matrix(constants$z_destination, nrow = n, ncol = 2L, byrow = TRUE))
  }

  # the composed severity reads the severity of its own branch, the null masses
  # the coordinate 'alpha'
  severity <- .bt_dnode_phacking_severity(
    spec, k, lookup,
    name = if(identical(quantity, "alpha")) .bt_dnode_phacking_node_name(spec, "alpha", k) else spec$names$alpha
  )
  if(is.null(severity) || identical(quantity, "alpha")){
    return(severity)
  }
  per_alpha <- if(identical(quantity, "pi_null")) constants$pi_null_per_alpha else constants$beta_null_per_alpha

  severity * per_alpha
}

# The values of a quantity on the draws of a lookup: the values of the active
# branch of every draw (one branch in a single prior); NULL when a coordinate
# is unavailable. This is the one evaluator of the family.
.bt_dnode_phacking_quantity <- function(spec, quantity, lookup, label = quantity){

  if(!spec$uses_indicator){
    return(.bt_dnode_phacking_branch_values(spec, quantity, 1L, lookup))
  }

  indicator <- .bt_deterministic_lookup_value(lookup, "bias_indicator")
  if(is.null(indicator)){
    return(NULL)
  }
  .bt_dnode_bias_indicator_check(indicator, length(spec$branches), label)

  out <- matrix(
    NA_real_, nrow = lookup$n,
    ncol = if(quantity %in% .bt_dnode_phacking_vector_quantities) 2L else 1L
  )
  for(k in unique(indicator)){
    rows <- which(indicator == k)
    branch_lookup <- .bt_deterministic_lookup(
      lookup$draws[rows, , drop = FALSE],
      lookup$prior_list
    )
    values <- .bt_dnode_phacking_branch_values(spec, quantity, k, branch_lookup)
    if(is.null(values)){
      return(NULL)
    }
    out[rows, ] <- values
  }

  out
}

.bt_dnode_phacking_evaluate <- function(node, lookup){

  .bt_dnode_phacking_quantity(node$spec, node$spec$quantity, lookup)
}

# The severity and the null masses (alpha, pi_null, beta_null) of one draw of
# bridge sampling or the marginal-likelihood parameters, with the arithmetic of
# the nodes; the composed severity of a mixture is evaluated from the component
# severities first, as the registered evaluator does. 'label' names the prior in
# the messages of invalid indicators, and 'missing' is the message of
# unavailable coordinates (default: that of the severity prior).
.bt_dnode_phacking_parameters <- function(spec, samples, label = "alpha", missing = NULL){

  lookup <- .bt_deterministic_row_lookup(samples)
  values <- list(alpha = .bt_dnode_phacking_quantity(spec, "alpha", lookup, label))
  if(!is.null(values$alpha)){
    if(spec$uses_indicator){
      draws <- lookup$draws[, setdiff(colnames(lookup$draws), spec$names$alpha), drop = FALSE]
      draws <- cbind(draws, matrix(values$alpha, ncol = 1L, dimnames = list(NULL, spec$names$alpha)))
      lookup <- .bt_deterministic_lookup(draws, lookup$prior_list)
    }
    for(quantity in c("pi_null", "beta_null")){
      values[quantity] <- list(.bt_dnode_phacking_quantity(spec, quantity, lookup, label))
    }
  }
  if(length(values) < 3L || any(vapply(values, is.null, logical(1)))){
    if(is.null(missing)){
      distribution <- spec$branches[[1L]]$phacking$alpha[["distribution"]]
      missing <- if(identical(distribution, "invgamma")){
        "'samples' does not contain all monitored inverse-gamma prior parameters."
      }else{
        "'samples' does not contain all monitored scalar prior parameters."
      }
    }
    .bt_JAGS_marglik_missing_columns(missing)
  }

  lapply(values, as.vector)
}

# The p-hacking specification of the most recent priors evaluated one draw at a
# time (see .bt_dnode_omega_cached()): found with identical().
.bt_dnode_phacking_cache <- new.env(parent = emptyenv())
.bt_dnode_phacking_cache_size <- 16L

.bt_dnode_phacking_cached <- function(prior){

  entries <- .bt_dnode_phacking_cache$entries
  for(entry in entries){
    if(identical(entry$prior, prior)){
      return(entry$spec)
    }
  }
  spec <- .bt_dnode_phacking_spec(lapply(.selection_normalize_priors(prior), .selection_branch_info))
  entries <- c(list(list(prior = prior, spec = spec)), entries)
  .bt_dnode_phacking_cache$entries <- entries[seq_len(min(length(entries), .bt_dnode_phacking_cache_size))]

  spec
}
