#' Generate a JAGS LKJ-Cholesky correlation module
#'
#' @description
#' Creates JAGS syntax for a Cholesky factor of a correlation matrix with the
#' same target distribution as Stan's `lkj_corr_cholesky(eta)`, using the
#' BayesTools compiled JAGS module.
#'
#' @param name character scalar. Prefix used for generated JAGS nodes.
#' @param K integer scalar. Dimension of the correlation matrix.
#' @param eta positive numeric scalar. LKJ concentration parameter.
#' @param include_correlation logical scalar. Whether to generate the
#'   deterministic correlation matrix `name_R`.
#' @param include_primitives logical scalar. Whether primitive beta/CPC nodes
#'   should be included in the returned monitor vector.
#'
#' @return A list with JAGS syntax, monitor names, primitive bridge coordinate
#'   names and bounds, and metadata for the generated LKJ-Cholesky block.
#'
#' @export
JAGS_lkj_corr_cholesky <- function(name, K, eta = 1,
                                   include_correlation = TRUE,
                                   include_primitives = FALSE){

  .bt_check_lkj_cholesky_module_inputs(name, K, eta, include_correlation, include_primitives)

  .bt_JAGS_lkj_corr_cholesky_module(
    name = name,
    K = K,
    eta = eta,
    include_correlation = include_correlation,
    include_primitives = include_primitives
  )
}

# `unit_columns`: the columns of an exchangeable unit of the block, whose
# within-unit primitives are fixed at 1/2 (see .bt_lkj_layout()). Only the free
# primitives are sampled and monitored; they keep the indices of their pairs.
.bt_JAGS_lkj_corr_cholesky_module <- function(name, K, eta,
                                              include_correlation,
                                              include_primitives,
                                              unit_columns = integer()){

  pairs <- .bt_lkj_cholesky_cpc_pairs(K = K, eta = eta)
  n_pairs <- nrow(pairs)
  layout <- .bt_lkj_layout(K, unit_columns)
  free_pairs <- layout$free_pairs
  n_free <- length(free_pairs)
  if(n_pairs > 0L && n_free == 0L){
    stop("An LKJ-Cholesky block needs at least one free primitive.", call. = FALSE)
  }

  L_name <- paste0(name, "_L")
  R_name <- paste0(name, "_R")
  u_name <- paste0(name, "_lkj_u")
  cpc_name <- paste0(name, "_lkj_cpc")
  alpha_name <- paste0(name, "_lkj_alpha")

  syntax <- c(paste0("# LKJ-Cholesky correlation module: ", name))

  if(n_free > 0L){
    for(p in free_pairs){
      syntax <- c(
        syntax,
        paste0(alpha_name, "[", p, "] <- ", .bt_JAGS_numeric_literal(pairs$alpha[p]))
      )
    }
    for(p in layout$fixed_pairs){
      syntax <- c(syntax, paste0(u_name, "[", p, "] <- 0.5"))
    }
    # the fixed primitives are the leading pairs, the free ones a contiguous tail
    syntax <- c(syntax, if(n_free == n_pairs){
      paste0(u_name, "[1:", n_pairs, "] ~ dbt_lkj_cpc(", alpha_name, ")")
    }else{
      free_range <- paste0("[", free_pairs[[1L]], ":", n_pairs, "]")
      paste0(u_name, free_range, " ~ dbt_lkj_cpc(", alpha_name, free_range, ")")
    })
  }

  # The Cholesky factor, correlation matrix, and partial correlations are the
  # registered 'lkj' deterministic nodes of the primitives.
  node <- .bt_dnode_lkj(
    name = name,
    K = K,
    include_correlation = include_correlation,
    include_primitives = include_primitives,
    unit_columns = unit_columns
  )
  syntax <- c(syntax, .bt_deterministic_node_emit(node))

  primitive_names <- node$spec$primitive_names
  cpc_names <- node$spec$cpc_names
  primitive_lb <- numeric(0)
  primitive_ub <- numeric(0)
  if(n_free > 0L){
    primitive_lb <- stats::setNames(rep(0, n_free), primitive_names)
    primitive_ub <- stats::setNames(rep(1, n_free), primitive_names)
  }

  monitor <- L_name
  if(include_correlation){
    monitor <- c(monitor, R_name)
  }
  if(include_primitives && n_free > 0L){
    monitor <- c(monitor, primitive_names, cpc_names)
  }

  out <- list(
    syntax = paste0(paste(syntax, collapse = "\n"), "\n"),
    monitor = monitor,
    name = name,
    K = K,
    eta = eta,
    cholesky_name = L_name,
    correlation_name = if(include_correlation) R_name else NULL,
    primitive_name = u_name,
    cpc_name = cpc_name,
    primitive_names = primitive_names,
    cpc_names = if(include_primitives) cpc_names else character(0),
    primitive_bounds = list(lb = primitive_lb, ub = primitive_ub),
    pairs = pairs,
    backend = "module",
    jags_module = "BayesTools",
    required_packages = "BayesTools",
    data = list()
  )

  class(out) <- c("BayesTools_JAGS_lkj_corr_cholesky", "list")
  out
}

.bt_check_lkj_cholesky_module_inputs <- function(name, K, eta,
                                                 include_correlation,
                                                 include_primitives){

  check_char(name, "name", allow_NA = FALSE)
  if(!grepl("^[A-Za-z][A-Za-z0-9_]*$", name)){
    stop("'name' must be a valid JAGS node prefix using letters, digits, and underscores.", call. = FALSE)
  }
  check_int(K, "K", lower = 1, allow_NA = FALSE)
  check_real(eta, "eta", lower = 0, allow_bound = FALSE, allow_NA = FALSE)
  if(!is.finite(eta)){
    stop("'eta' must be finite.", call. = FALSE)
  }
  check_bool(include_correlation, "include_correlation", allow_NA = FALSE)
  check_bool(include_primitives, "include_primitives", allow_NA = FALSE)

  invisible(TRUE)
}

.bt_lkj_cholesky_check_K <- function(K){

  if(!is.numeric(K) || length(K) != 1L || is.na(K) || !is.finite(K) ||
     K < 1L || K > .Machine$integer.max || K != floor(K)){
    stop("'K' must be a positive integer scalar.", call. = FALSE)
  }
  as.integer(K)
}

.bt_lkj_cholesky_check_eta <- function(eta){

  if(!is.numeric(eta) || length(eta) != 1L || is.na(eta) || !is.finite(eta) || eta <= 0){
    stop("'eta' must be a positive finite scalar.", call. = FALSE)
  }
  as.numeric(eta)
}

.bt_lkj_cholesky_n_pairs <- function(K){
  K * (K - 1L) / 2L
}

.bt_lkj_cholesky_cpc_pairs <- function(K, eta = 1){

  K <- .bt_lkj_cholesky_check_K(K)
  eta <- .bt_lkj_cholesky_check_eta(eta)

  if(K == 1L){
    return(data.frame(
      index = integer(0),
      i = integer(0),
      j = integer(0),
      alpha = numeric(0)
    ))
  }

  out <- vector("list", K * (K - 1L) / 2L)
  alpha <- .bt_lkj_cholesky_alpha(K = K, eta = eta)
  p <- 0L
  for(j in 2:K){
    for(i in 1:(j - 1L)){
      p <- p + 1L
      out[[p]] <- data.frame(
        index = p,
        i = i,
        j = j,
        alpha = alpha[p]
      )
    }
  }

  do.call(rbind, out)
}

.bt_lkj_cholesky_cpc_u_to_L <- function(u, K){

  K <- .bt_lkj_cholesky_check_K(K)
  n_pairs <- .bt_lkj_cholesky_n_pairs(K)
  .bt_lkj_cholesky_check_u_shape(u, n_pairs)
  .bt_lkj_cholesky_check_u_support(u)

  .BayesTools_require_native_lkj()
  .Call("BayesTools_lkj_cholesky_from_u", u, as.integer(K), PACKAGE = "BayesTools")
}

# The lower Cholesky factor in design order of the correlation matrix whose
# vine-order factor .bt_lkj_cholesky_cpc_u_to_L() builds: design column d is
# vine variable position[d]. It is the module's bt_lkj_cholesky_ordered()
# kernel, so JAGS and R agree to the last bit; the identity order returns the
# vine factor itself.
.bt_lkj_cholesky_cpc_u_to_L_ordered <- function(u, K, position){

  K <- .bt_lkj_cholesky_check_K(K)
  n_pairs <- .bt_lkj_cholesky_n_pairs(K)
  .bt_lkj_cholesky_check_u_shape(u, n_pairs)
  .bt_lkj_cholesky_check_u_support(u)
  if(!is.numeric(position) || length(position) != K || anyNA(position) ||
     !identical(sort(as.integer(position)), seq_len(K)) ||
     any(position != as.integer(position))){
    stop("'position' must be a permutation of 1, ..., K.", call. = FALSE)
  }

  .BayesTools_require_native_lkj()
  .Call("BayesTools_lkj_cholesky_from_u_ordered", u, as.integer(K), as.numeric(position),
        PACKAGE = "BayesTools")
}

# Correlation matrices R = L L' of lower Cholesky factors of correlation
# matrices, such as .bt_lkj_cholesky_cpc_u_to_L() output. The rows of such a
# factor have unit norm by construction, so the diagonal is exactly 1, as in
# the module's bt_lkj_corr(), instead of the rounded row sums of squares
# (1 +/- a few eps), which would make a constant diagonal vary between draws.
# The off-diagonal cells sum the products in the module's order, m = 1, ...,
# column, so they reproduce its monitors bit for bit; tcrossprod() would sum in
# the order of the BLAS that R uses.
# `L` is a K x K matrix or an n x K x K array; the result has the same shape.
# Draws with a missing factor (NA entries) keep missing diagonals.
.bt_lkj_cholesky_L_to_R <- function(L){

  single <- is.matrix(L)
  if(single){
    L <- array(L, dim = c(1L, dim(L)))
  }
  if(!is.array(L) || length(dim(L)) != 3L || dim(L)[2L] != dim(L)[3L]){
    stop("'L' must be a K x K matrix or an n x K x K array.", call. = FALSE)
  }

  K <- dim(L)[2L]
  R <- array(NA_real_, dim = dim(L))
  for(row in seq_len(K)){
    for(column in seq_len(row - 1L)){
      cell <- 0
      for(m in seq_len(column)){
        cell <- cell + L[, row, m] * L[, column, m]
      }
      R[, row, column] <- cell
      R[, column, row] <- cell
    }
    diagonal <- rowSums(L[, row, , drop = FALSE]^2)
    diagonal[is.finite(diagonal)] <- 1
    R[, row, row] <- diagonal
  }

  if(single){
    return(matrix(R[1L, , ], nrow = K, ncol = K))
  }
  R
}

# Exact inverse of .bt_lkj_cholesky_cpc_u_to_L(): recover the LKJ primitive
# coordinates from lower Cholesky factors of correlation matrices. For row i
# and column j < i, the canonical partial correlation is
#   z_ij = L_ij / sqrt(1 - sum_{k < j} L_ik^2) = L_ij / sqrt(sum_{k >= j} L_ik^2)
# and u_ij = (z_ij + 1) / 2, ordered like the native map (column-major over the
# upper-triangular pairs). `L` is a K x K matrix or an n x K x K array; the
# result is a vector or an n x K(K - 1)/2 matrix. Factors that are not valid
# lower Cholesky factors of a correlation matrix with an interior primitive
# (unit row norms, positive diagonal, |z| < 1) give NA rows.
.bt_lkj_cholesky_L_to_cpc_u <- function(L, K){

  K <- .bt_lkj_cholesky_check_K(K)
  n_pairs <- .bt_lkj_cholesky_n_pairs(K)
  single <- is.matrix(L)
  if(single){
    L <- array(L, dim = c(1L, dim(L)))
  }
  if(!is.array(L) || length(dim(L)) != 3L ||
     dim(L)[2L] != K || dim(L)[3L] != K){
    stop("'L' must be a K x K matrix or an n x K x K array.", call. = FALSE)
  }

  n_draws <- dim(L)[1L]
  u <- matrix(NA_real_, nrow = n_draws, ncol = n_pairs)
  if(n_pairs == 0L){
    return(if(single) numeric(0) else u)
  }

  valid <- rep(TRUE, n_draws)
  tolerance <- 1e-10
  for(row in seq_len(K)){
    row_values <- matrix(L[, row, , drop = FALSE], nrow = n_draws)
    valid <- valid & apply(is.finite(row_values), 1L, all) &
      abs(rowSums(row_values^2) - 1) <= tolerance &
      row_values[, row] > 0
    if(row < K){
      valid <- valid & rowSums(abs(row_values[, (row + 1L):K, drop = FALSE])) == 0
    }
  }

  for(row in seq_len(K)[-1L]){
    # For unit rows, 1 - sum_{k < j} L_ik^2 equals the tail sum
    # sum_{k >= j} L_ik^2; the tail sum avoids catastrophic cancellation when
    # earlier partial correlations are close to +/-1.
    remaining <- matrix(0, nrow = n_draws, ncol = row)
    tail_sum <- rep(0, n_draws)
    for(column in rev(seq_len(row))){
      tail_sum <- tail_sum + L[, row, column]^2
      remaining[, column] <- tail_sum
    }
    for(column in seq_len(row - 1L)){
      pair <- (row - 1L) * (row - 2L) / 2L + column
      cpc <- L[, row, column] / sqrt(remaining[, column])
      valid <- valid & is.finite(cpc) & abs(cpc) < 1
      u[, pair] <- (cpc + 1) / 2
    }
  }
  u[!valid, ] <- NA_real_

  if(single){
    return(as.numeric(u[1L, ]))
  }
  u
}

.bt_lkj_cholesky_check_u_shape <- function(u, n_pairs){

  if(is.matrix(u)){
    if(ncol(u) != n_pairs){
      stop("'u' must have K * (K - 1) / 2 columns.", call. = FALSE)
    }
  }else{
    if(length(u) != n_pairs){
      stop("'u' must have length K * (K - 1) / 2.", call. = FALSE)
    }
  }

  invisible(TRUE)
}

.bt_lkj_cholesky_u_in_support <- function(u){

  if(is.matrix(u)){
    if(ncol(u) == 0L){
      return(rep(TRUE, nrow(u)))
    }
    return(rowSums(!is.finite(u) | u <= 0 | u >= 1) == 0L)
  }

  if(length(u) == 0L){
    return(TRUE)
  }
  all(is.finite(u) & u > 0 & u < 1)
}

.bt_lkj_cholesky_check_u_support <- function(u){

  if(!all(.bt_lkj_cholesky_u_in_support(u))){
    stop("'u' values must be finite and strictly between 0 and 1.", call. = FALSE)
  }

  invisible(TRUE)
}

.bt_lkj_cholesky_cpc_u_log_prior <- function(u, K, eta = 1){

  K <- .bt_lkj_cholesky_check_K(K)
  eta <- .bt_lkj_cholesky_check_eta(eta)
  .bt_lkj_cholesky_u_log_density(u, alpha = .bt_lkj_cholesky_alpha(K = K, eta = eta))
}

# The log density of primitives u under independent Beta(alpha, alpha) laws
# (-Inf outside (0, 1)); `u` is a vector or a draws x length(alpha) matrix.
.bt_lkj_cholesky_u_log_density <- function(u, alpha){

  n_pairs <- length(alpha)

  if(is.matrix(u)){
    .bt_lkj_cholesky_check_u_shape(u, n_pairs)
    if(nrow(u) == 0L){
      return(numeric(0))
    }
    support <- .bt_lkj_cholesky_u_in_support(u)
    out <- rep(-Inf, nrow(u))
    if(any(support)){
      .BayesTools_require_native_lkj()
      out[support] <- .Call(
        "BayesTools_lkj_log_prior_u",
        u[support, , drop = FALSE],
        alpha,
        PACKAGE = "BayesTools"
      )
    }
    return(out)
  }

  .bt_lkj_cholesky_check_u_shape(u, n_pairs)
  if(!.bt_lkj_cholesky_u_in_support(u)){
    return(-Inf)
  }
  .BayesTools_require_native_lkj()
  as.numeric(.Call("BayesTools_lkj_log_prior_u", u, alpha, PACKAGE = "BayesTools"))
}

.bt_lkj_cholesky_alpha <- function(K, eta){

  K <- .bt_lkj_cholesky_check_K(K)
  eta <- .bt_lkj_cholesky_check_eta(eta)
  n_pairs <- .bt_lkj_cholesky_n_pairs(K)

  if(n_pairs == 0L){
    return(numeric(0))
  }

  alpha <- numeric(n_pairs)
  p <- 0L
  for(j in 2:K){
    for(i in 1:(j - 1L)){
      p <- p + 1L
      alpha[p] <- eta + (K - i - 1L) / 2
    }
  }

  alpha
}

.bt_lkj_cholesky_flat_index <- function(row, column, K){
  (row - 1L) * K + column
}

# LKJ vine layouts --------------------------------------------------------------
#
# A correlated block builds its C-vine in a vine order: design column d is vine
# variable position[d]. The columns of an exchangeable unit (`unit_columns`, in
# design order) come first, and their within-unit primitives, the leading pairs
# (i, j) with j <= k of the column-major pair order, are fixed at u = 1/2, so
# the unit block of the correlation matrix is the identity. The other (free)
# primitives keep their LKJ Beta laws, which gives the LKJ(eta) density
# det(R)^(eta - 1) restricted to that block (Lewandowski, Kurowicka & Joe, 2009:
# the vine partial correlations are independent, and those of the unit block
# determine it). Without a unit the vine order is the design order and every
# primitive is free (the identity layout).
.bt_lkj_layout <- function(K, unit_columns = integer()){

  K <- .bt_lkj_cholesky_check_K(K)
  if(!is.numeric(unit_columns) || anyNA(unit_columns) ||
     any(unit_columns != floor(unit_columns)) ||
     any(unit_columns < 1 | unit_columns > K) ||
     anyDuplicated(unit_columns) || length(unit_columns) == 1L){
    stop("'unit_columns' must be empty or hold at least two distinct columns of the block.", call. = FALSE)
  }
  unit_columns <- as.integer(unit_columns)
  vine_order <- c(unit_columns, setdiff(seq_len(K), unit_columns))
  n_pairs <- .bt_lkj_cholesky_n_pairs(K)
  n_fixed <- .bt_lkj_cholesky_n_pairs(length(unit_columns))

  list(
    K            = K,
    unit_columns = unit_columns,
    vine_order   = vine_order,
    position     = match(seq_len(K), vine_order),
    fixed_pairs  = seq_len(n_fixed),
    free_pairs   = seq_len(n_pairs)[seq_len(n_pairs) > n_fixed]
  )
}

.bt_lkj_layout_is_identity <- function(layout){
  length(layout$unit_columns) == 0L
}

.bt_lkj_layout_identity_order <- function(layout){
  identical(layout$position, seq_len(layout$K))
}

# The full primitives of the free primitives u (a vector, or a draws x free
# matrix): the fixed primitives are 1/2.
.bt_lkj_layout_full_u <- function(u, layout){

  n_free <- length(layout$free_pairs)
  if(if(is.matrix(u)) ncol(u) != n_free else length(u) != n_free){
    stop("'u' must hold one value per free LKJ primitive.", call. = FALSE)
  }
  if(length(layout$fixed_pairs) == 0L){
    return(u)
  }

  n_pairs <- .bt_lkj_cholesky_n_pairs(layout$K)
  if(is.matrix(u)){
    out <- matrix(0.5, nrow = nrow(u), ncol = n_pairs)
    out[, layout$free_pairs] <- u
  }else{
    out <- rep(0.5, n_pairs)
    out[layout$free_pairs] <- u
  }
  out
}

# The design-order lower Cholesky factors of the free primitives u (the
# module's bt_lkj_cholesky() in the vine order, bt_lkj_cholesky_ordered()
# otherwise): a K x K matrix, or a draws x K x K array.
.bt_lkj_layout_u_to_L <- function(u, layout){

  u_full <- .bt_lkj_layout_full_u(u, layout)
  if(.bt_lkj_layout_identity_order(layout)){
    return(.bt_lkj_cholesky_cpc_u_to_L(u_full, K = layout$K))
  }
  .bt_lkj_cholesky_cpc_u_to_L_ordered(u_full, K = layout$K, position = layout$position)
}

# The design-order correlation matrices of the free primitives u: the cells of
# the vine-order correlation matrix (the module's bt_lkj_corr()) permuted, so
# that the unit block is exactly the identity. R = L L' of the design-order
# factor only up to rounding.
.bt_lkj_layout_u_to_R <- function(u, layout){

  u_full <- .bt_lkj_layout_full_u(u, layout)
  R <- .bt_lkj_cholesky_L_to_R(.bt_lkj_cholesky_cpc_u_to_L(u_full, K = layout$K))
  if(.bt_lkj_layout_identity_order(layout)){
    return(R)
  }
  position <- layout$position
  if(is.matrix(R)){
    return(R[position, position, drop = FALSE])
  }
  R[, position, position, drop = FALSE]
}

.bt_lkj_layout_alpha <- function(layout, eta){

  .bt_lkj_cholesky_alpha(K = layout$K, eta = eta)[layout$free_pairs]
}

# The log prior density of the free primitives u (the fixed primitives have no
# density: they are constants of the vine, whose Jacobian factor is 1).
.bt_lkj_layout_log_prior <- function(u, layout, eta){

  eta <- .bt_lkj_cholesky_check_eta(eta)
  .bt_lkj_cholesky_u_log_density(u, alpha = .bt_lkj_layout_alpha(layout, eta))
}

# The free primitives of design-order lower Cholesky factors L (a K x K matrix
# or a draws x K x K array; the inverse of .bt_lkj_layout_u_to_L()): the
# primitives of the vine-order factor of the reordered correlation matrix. Rows
# are NA where L is not a valid factor (see .bt_lkj_cholesky_L_to_cpc_u()) or
# where its unit block is not the identity (the fixed primitives differ from
# 1/2 by more than 1e-10), which the layout cannot represent.
.bt_lkj_layout_L_to_u <- function(L, layout){

  K <- layout$K
  if(.bt_lkj_layout_is_identity(layout)){
    return(.bt_lkj_cholesky_L_to_cpc_u(L, K = K))
  }

  single <- is.matrix(L)
  if(single){
    L <- array(L, dim = c(1L, dim(L)))
  }
  if(!is.array(L) || length(dim(L)) != 3L || dim(L)[2L] != K || dim(L)[3L] != K){
    stop("'L' must be a K x K matrix or an n x K x K array.", call. = FALSE)
  }

  n_draws <- dim(L)[1L]
  # validity of the design-order factors (unit rows, positive diagonal)
  valid <- !is.na(.bt_lkj_cholesky_L_to_cpc_u(L, K = K)[, 1L])
  vine_L <- array(NA_real_, dim = dim(L))
  if(any(valid)){
    R <- .bt_lkj_cholesky_L_to_R(L[valid, , , drop = FALSE])
    vine_order <- layout$vine_order
    for(i in seq_len(sum(valid))){
      factor <- try(t(chol(R[i, vine_order, vine_order])), silent = TRUE)
      if(!inherits(factor, "try-error")){
        vine_L[which(valid)[i], , ] <- factor
      }
    }
  }
  u_full <- .bt_lkj_cholesky_L_to_cpc_u(vine_L, K = K)
  fixed <- u_full[, layout$fixed_pairs, drop = FALSE]
  u_full[rowSums(is.na(fixed) | abs(fixed - 0.5) > 1e-10) > 0L, ] <- NA_real_
  u <- u_full[, layout$free_pairs, drop = FALSE]

  if(single){
    return(as.numeric(u[1L, ]))
  }
  u
}
