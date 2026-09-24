# Random-effect deterministic node families.


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

