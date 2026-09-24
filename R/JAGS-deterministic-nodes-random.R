# Random-effect deterministic node families.


# Variance-allocation SDs ("random_sd") ----------------------------------------
#
# A variance allocation distributes a source SD (an allocation's own SD prior
# or an external scalar source) over random-effect blocks or over the SD
# components of one block. Each allocated SD is the source multiplied by a
# chain of allocation factors, from the root allocation to the target:
#   sd = source * [gate_1 *] sqrt(m_1) * [gate_2 *] sqrt(m_2) * ...
# with m = w[i] (total variance) or m = K * w[i] (mean variance) of a Dirichlet
# weight vector w of dimension K, and optional Bernoulli inclusion gates; a
# gate-only factor contributes its gate. The model defines
#   - block SDs as 'sd = source * <factor chain>',
#   - SD components of a block as 'leaf = <source node> * sqrt(m)', where the
#     source node is the root source or the consumed component node of the
#     parent allocation,
#   - consumed parent components as 'component = <source node> * <factor>'.
# The R evaluator applies the same chain to draws of the source, the weights,
# and the gates. Each factor multiplies by its multiplier and then by its gate;
# the gates are 0/1 indicators, so the product equals the JAGS product exactly.

.bt_dnode_allocation_multiplier_expression <- function(weight_name, index,
                                                       scale, n_targets){

  multiplier <- if(identical(scale, "mean_variance")){
    paste0(n_targets, " * ", weight_name, "[", index, "]")
  }else{
    paste0(weight_name, "[", index, "]")
  }

  paste0("sqrt(", multiplier, ")")
}

.bt_dnode_allocation_factor_expression <- function(factor){

  if(is.null(factor$weight_name)){
    return(factor$inclusion_name)
  }
  multiplier <- .bt_dnode_allocation_multiplier_expression(
    weight_name = factor$weight_name,
    index = factor$index,
    scale = factor$scale,
    n_targets = factor$n_targets
  )

  if(!is.null(factor$inclusion_name)){
    multiplier <- paste0(factor$inclusion_name, " * ", multiplier)
  }

  multiplier
}

.bt_dnode_allocation_factors_expression <- function(factors){

  if(length(factors) == 0L){
    return("1")
  }

  paste(
    vapply(factors, .bt_dnode_allocation_factor_expression, character(1)),
    collapse = " * "
  )
}

# 'source * factor' of one allocation component, with a JAGS expression of the
# source (a node name or a row-indexed source expression).
.bt_dnode_allocation_expression <- function(source_expression, factor){

  paste0(source_expression, " * ", .bt_dnode_allocation_factor_expression(factor))
}

.bt_dnode_allocation_multiplier <- function(weights, scale, n_targets){

  if(identical(scale, "mean_variance")){
    if(!is.numeric(n_targets) || length(n_targets) != 1L ||
       is.na(n_targets) || n_targets < 1L){
      stop(
        "Random-effect allocation metadata are missing canonical 'allocation$n_targets'.",
        call. = FALSE
      )
    }
    return(sqrt(n_targets * weights))
  }
  if(identical(scale, "total_variance")){
    return(sqrt(weights))
  }

  stop(
    "Random-effect allocation metadata are missing canonical 'allocation$scale'.",
    call. = FALSE
  )
}

# The factor chain applied to 'base' (draws of the source, or 1). The
# consumer supplies how a factor's Dirichlet weights (a draws x K matrix, or
# NULL when unavailable) and its gate (a 0/1 vector; ones without a gate) are
# read, so that draws, bridge rows, and supplied parameter values share the
# arithmetic. Returns NULL when weights are unavailable.
.bt_dnode_allocation_chain <- function(base, factors, weights_of, gate_of){

  out <- base
  for(factor in factors){
    multiplier <- 1
    if(!is.null(factor$weight_name)){
      weights <- weights_of(factor)
      if(is.null(weights)){
        return(NULL)
      }
      multiplier <- .bt_dnode_allocation_multiplier(
        weights = weights[, factor$index],
        scale = factor$scale,
        n_targets = factor$n_targets
      )
    }
    gate <- gate_of(factor)
    if(is.null(gate)){
      stop(
        "Random-effect allocation inclusion samples are missing Bernoulli indicator '",
        factor$inclusion_name,
        "'.",
        call. = FALSE
      )
    }
    out <- out * multiplier * gate
  }

  out
}

# The coordinates a chain of allocation factors reads: the Dirichlet weights
# and the inclusion gates of every factor.
.bt_dnode_allocation_factor_dependencies <- function(factors){

  unlist(lapply(factors, function(factor){
    c(
      if(!is.null(factor$weight_name)){
        paste0(factor$weight_name, "[", seq_len(factor$n_targets), "]")
      },
      factor$inclusion_name
    )
  }), use.names = FALSE)
}

# 'node' = 'emit_source * <emit_factors>' in the model syntax; its value is the
# source times the whole 'factors' chain from the root allocation.
.bt_dnode_random_sd <- function(name, source_name, factors, emit_source,
                                emit_factors, parameter = NA_character_,
                                block = NA_character_){

  .bt_deterministic_node(
    family = "random_sd",
    node = name,
    coordinates = name,
    dependencies = c(source_name, .bt_dnode_allocation_factor_dependencies(factors)),
    parameter = parameter,
    block = block,
    spec = list(
      source_name = source_name,
      factors = factors,
      emit_source = emit_source,
      emit_factors = emit_factors
    )
  )
}

.bt_dnode_random_sd_emit <- function(node){

  spec <- node$spec
  factor_expression <- .bt_dnode_allocation_factors_expression(spec$emit_factors)
  paste0(
    node$node,
    " = ",
    paste(
      c(spec$emit_source, if(!identical(factor_expression, "1")) factor_expression),
      collapse = " * "
    )
  )
}

.bt_dnode_random_sd_evaluate <- function(node, lookup){

  spec <- node$spec
  base <- .bt_deterministic_lookup_value(lookup, spec$source_name)
  if(is.null(base)){
    return(NULL)
  }

  .bt_dnode_allocation_chain(
    base = base,
    factors = spec$factors,
    weights_of = function(factor){
      .bt_deterministic_lookup_simplex(lookup, factor)
    },
    gate_of = function(factor){
      .bt_random_effect_allocation_gate_draws(
        parameter_name = factor$inclusion_name,
        posterior = lookup$draws
      )
    }
  )
}

# The allocated SD nodes of a random-effect block: one block SD, or one node
# per SD component of an 'sd_component' allocation. Blocks with a row-indexed
# external source define no SD nodes.
.bt_dnode_random_sd_from_random_term <- function(random_term,
                                                 parameter = NA_character_){

  binding <- random_term$sd_binding
  if(is.null(binding) || !isTRUE(binding$true_allocation) ||
     length(binding$allocations) == 0L ||
     .bt_random_sd_binding_has_row_external_source(binding)){
    return(list())
  }
  allocation <- binding$allocations[[1L]]
  target <- .bt_random_effect_allocation_target_metadata(allocation)
  source_name <- .bt_random_sd_binding_source_name(allocation$source)

  if(identical(target, "block")){
    name <- unique(random_term$sd_parameter_names)
    if(length(name) != 1L || is.na(name)){
      stop(
        "Random-effect allocation metadata",
        .bt_random_effect_metadata_block_detail(random_term),
        " must define one block SD parameter.",
        call. = FALSE
      )
    }
    factors <- .bt_random_effect_allocation_factor_plan(
      .bt_random_effect_allocation_factors_metadata(allocation)
    )
    return(list(.bt_dnode_random_sd(
      name = name,
      source_name = source_name,
      factors = factors,
      emit_source = .bt_random_sd_binding_source_jags_expression(allocation$source),
      emit_factors = factors,
      parameter = parameter,
      block = random_term$block_name
    )))
  }

  component <- .bt_check_random_sd_component_allocation(
    allocation = allocation,
    n_columns = random_term$n_columns,
    context = "Random-effect allocation metadata"
  )
  parent_factors <- .bt_random_effect_allocation_factor_plan(component$parent_factors)
  lapply(seq_len(component$n_targets), function(leaf){
    leaf_factor <- .bt_random_variance_allocation_factor(
      weight_name = allocation$weight_name,
      index = leaf,
      scale = allocation$scale,
      n_targets = component$n_targets
    )
    .bt_dnode_random_sd(
      name = allocation$leaf_names[[leaf]],
      source_name = source_name,
      factors = c(parent_factors, list(leaf_factor)),
      emit_source = allocation$source_node,
      emit_factors = list(leaf_factor),
      parameter = parameter,
      block = random_term$block_name
    )
  })
}

# The consumed components of a formula's variance allocations: parent
# allocation components that a child allocation splits further are defined as
# nodes of their own.
.bt_dnode_random_sd_components <- function(design){

  nodes <- list()
  for(allocation in design$random_allocations){
    if(!identical(allocation$target, "block")){
      next
    }
    for(component in allocation$components){
      if(identical(component$node_name, component$expression)){
        next
      }
      factors <- .bt_random_effect_allocation_factor_plan(component$factors)
      nodes[[length(nodes) + 1L]] <- .bt_dnode_random_sd(
        name = component$node_name,
        source_name = .bt_random_sd_binding_source_name(component$source),
        factors = factors,
        emit_source = allocation$source_node,
        emit_factors = factors[length(factors)],
        parameter = design$parameter
      )
    }
  }

  nodes
}


# Scalar correlations ("random_rho") ------------------------------------------
#
# A structured random-effect block (cs, hcs, ar1, har, car) with a correlation
# prior on the Fisher-z or logit scale samples 'sample_name' and defines
#   rho <- tanh(sample_name)                                   (Fisher-z)
#   rho <- lower + (upper - lower) * ilogit(sample_name)       (logit)
# with the structure's correlation bounds. A raw ('rho' scale) correlation is
# the sampled node itself and has no generated node.

.bt_dnode_rho <- function(rho_name, sample_name, rho_scale, bounds,
                          sample_fixed = NULL, parameter = NA_character_,
                          block = NA_character_){

  if(!rho_scale %in% c("fisher_z", "logit")){
    stop("Scalar correlation nodes are generated only on the Fisher-z or logit scale.",
         call. = FALSE)
  }

  .bt_deterministic_node(
    family = "random_rho",
    node = rho_name,
    coordinates = rho_name,
    dependencies = sample_name,
    parameter = parameter,
    block = block,
    spec = list(
      sample_name = sample_name,
      rho_scale = rho_scale,
      bounds = c(lower = bounds[["lower"]], upper = bounds[["upper"]]),
      sample_fixed = sample_fixed
    )
  )
}

.bt_dnode_rho_from_random_term <- function(random_term, parameter = NA_character_){

  correlation <- random_term$correlation
  if(!is.list(correlation) || !identical(correlation$type, "rho") ||
     !correlation$rho_scale %in% c("fisher_z", "logit")){
    return(NULL)
  }
  plan <- .bt_random_effect_compile_rho_draw_plan(
    random_term = random_term,
    context = "Random-effect deterministic node metadata"
  )

  .bt_dnode_rho(
    rho_name = correlation$rho_name,
    sample_name = correlation$sample_name,
    rho_scale = plan$rho_scale,
    bounds = plan$bounds,
    sample_fixed = plan$sample_fixed,
    parameter = parameter,
    block = random_term$block_name
  )
}

.bt_dnode_rho_emit <- function(node){

  spec <- node$spec
  if(identical(spec$rho_scale, "fisher_z")){
    return(paste0(node$node, " <- tanh(", spec$sample_name, ")"))
  }

  paste0(
    node$node, " <- ",
    .bt_JAGS_numeric_literal(spec$bounds[["lower"]]), " + ",
    .bt_JAGS_numeric_literal(spec$bounds[["upper"]] - spec$bounds[["lower"]]),
    " * ilogit(", spec$sample_name, ")"
  )
}

# The map from the sampled coordinate to the correlation.
.bt_dnode_rho_transform <- function(value, rho_scale, bounds){

  if(identical(rho_scale, "fisher_z")){
    return(tanh(value))
  }
  if(identical(rho_scale, "logit")){
    return(
      bounds[["lower"]] +
        (bounds[["upper"]] - bounds[["lower"]]) * stats::plogis(value)
    )
  }

  value
}

.bt_dnode_rho_evaluate <- function(node, lookup){

  spec <- node$spec
  value <- .bt_deterministic_lookup_value(lookup, spec$sample_name)
  if(is.null(value) && !is.null(spec$sample_fixed)){
    value <- rep(spec$sample_fixed, lookup$n)
  }
  if(is.null(value)){
    return(NULL)
  }

  .bt_dnode_rho_transform(value, spec$rho_scale, spec$bounds)
}


# LKJ Cholesky factors ("lkj") ----------------------------------------------------
#
# The LKJ-Cholesky module samples the canonical-partial-correlation primitives
# '<name>_lkj_u' and defines the Cholesky factor '<name>_L', the correlation
# matrix '<name>_R' = L L' (with an exact unit diagonal) through the
# BayesTools JAGS module, and the monitored partial correlations
# '<name>_lkj_cpc[p] <- 2 * u[p] - 1'. The R evaluator computes L with the
# module's own native kernel.

.bt_dnode_lkj <- function(name, K, include_correlation = TRUE,
                          include_primitives = FALSE,
                          parameter = NA_character_, block = NA_character_){

  n_pairs <- .bt_lkj_cholesky_n_pairs(K)
  L_name <- paste0(name, "_L")
  R_name <- paste0(name, "_R")
  u_name <- paste0(name, "_lkj_u")
  cpc_name <- paste0(name, "_lkj_cpc")
  primitive_names <- if(n_pairs > 0L){
    paste0(u_name, "[", seq_len(n_pairs), "]")
  }else{
    character()
  }
  cpc_names <- if(isTRUE(include_primitives) && n_pairs > 0L){
    paste0(cpc_name, "[", seq_len(n_pairs), "]")
  }else{
    character()
  }
  cells <- function(matrix_name){
    paste0(
      matrix_name,
      "[", rep(seq_len(K), each = K), ",", rep(seq_len(K), times = K), "]"
    )
  }

  .bt_deterministic_node(
    family = "lkj",
    node = name,
    coordinates = c(
      cells(L_name),
      if(isTRUE(include_correlation)) cells(R_name),
      cpc_names
    ),
    dependencies = primitive_names,
    parameter = parameter,
    block = block,
    spec = list(
      K = K,
      L_name = L_name,
      R_name = if(isTRUE(include_correlation)) R_name else NULL,
      u_name = u_name,
      primitive_names = primitive_names,
      cpc_name = cpc_name,
      cpc_names = cpc_names
    )
  )
}

.bt_dnode_lkj_from_random_term <- function(random_term, parameter = NA_character_){

  correlation <- random_term$correlation
  if(!is.list(correlation) || !identical(correlation$type, "lkj")){
    return(NULL)
  }
  K <- random_term$n_columns
  primitive_names <- .bt_random_effect_lkj_primitive_names(
    random_term,
    K,
    context = "Random-effect deterministic node metadata"
  )
  node <- .bt_dnode_lkj(
    name = sub("_L$", "", correlation$cholesky_name),
    K = K,
    include_correlation = !is.null(correlation$correlation_name),
    include_primitives = length(correlation$cpc_names) > 0L,
    parameter = parameter,
    block = random_term$block_name
  )
  if(!identical(node$spec$primitive_names, primitive_names) ||
     !identical(node$spec$cpc_names, as.character(correlation$cpc_names))){
    stop(
      "Stored LKJ primitive metadata",
      .bt_random_effect_metadata_block_detail(random_term),
      " do not match the random-effect dimension.",
      call. = FALSE
    )
  }

  node
}

.bt_dnode_lkj_emit <- function(node){

  spec <- node$spec
  K <- spec$K
  L_flat_name <- paste0(node$node, "_L_flat")
  R_flat_name <- paste0(node$node, "_R_flat")

  syntax <- character()
  if(length(spec$primitive_names) > 0L){
    syntax <- c(
      syntax,
      paste0(L_flat_name, "[1:", K * K, "] <- bt_lkj_cholesky(", spec$u_name, ", ", K, ")")
    )
    if(!is.null(spec$R_name)){
      syntax <- c(
        syntax,
        paste0(R_flat_name, "[1:", K * K, "] <- bt_lkj_corr(", spec$u_name, ", ", K, ")")
      )
    }
  }

  cell_syntax <- function(matrix_name, flat_name){
    out <- character()
    for(row in seq_len(K)){
      for(column in seq_len(K)){
        target <- paste0(matrix_name, "[", row, ",", column, "]")
        if(K == 1L){
          out <- c(out, paste0(target, " <- 1"))
        }else{
          flat_index <- .bt_lkj_cholesky_flat_index(row, column, K)
          out <- c(out, paste0(target, " <- ", flat_name, "[", flat_index, "]"))
        }
      }
    }
    out
  }
  syntax <- c(syntax, cell_syntax(spec$L_name, L_flat_name))
  if(!is.null(spec$R_name)){
    syntax <- c(syntax, cell_syntax(spec$R_name, R_flat_name))
  }
  for(p in seq_along(spec$cpc_names)){
    syntax <- c(
      syntax,
      paste0(spec$cpc_name, "[", p, "] <- 2 * ", spec$u_name, "[", p, "] - 1")
    )
  }

  syntax
}

.bt_dnode_lkj_evaluate <- function(node, lookup){

  spec <- node$spec
  K <- spec$K
  n <- lookup$n
  if(length(spec$primitive_names) == 0L){
    return(matrix(1, nrow = n, ncol = length(node$coordinates)))
  }
  u <- .bt_deterministic_lookup_values(lookup, spec$primitive_names)
  if(is.null(u)){
    return(NULL)
  }

  L <- .bt_lkj_cholesky_cpc_u_to_L(u, K = K)
  cells <- function(values){
    out <- matrix(NA_real_, nrow = n, ncol = K * K)
    for(row in seq_len(K)){
      for(column in seq_len(K)){
        out[, (row - 1L) * K + column] <- values[, row, column]
      }
    }
    out
  }
  values <- cells(L)
  if(!is.null(spec$R_name)){
    values <- cbind(values, cells(.bt_lkj_cholesky_L_to_R(L)))
  }
  if(length(spec$cpc_names) > 0L){
    values <- cbind(values, 2 * u - 1)
  }

  values
}

