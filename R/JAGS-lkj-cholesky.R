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

.bt_JAGS_lkj_corr_cholesky_module <- function(name, K, eta,
                                              include_correlation,
                                              include_primitives){

  pairs <- .bt_lkj_cholesky_cpc_pairs(K = K, eta = eta)
  n_pairs <- nrow(pairs)

  L_name <- paste0(name, "_L")
  R_name <- paste0(name, "_R")
  u_name <- paste0(name, "_lkj_u")
  cpc_name <- paste0(name, "_lkj_cpc")
  alpha_name <- paste0(name, "_lkj_alpha")

  syntax <- c(paste0("# LKJ-Cholesky correlation module: ", name))

  if(n_pairs > 0L){
    for(p in seq_len(n_pairs)){
      syntax <- c(
        syntax,
        paste0(alpha_name, "[", p, "] <- ", .bt_jags_number(pairs$alpha[p]))
      )
    }
    syntax <- c(
      syntax,
      paste0(u_name, "[1:", n_pairs, "] ~ dbt_lkj_cpc(", alpha_name, ")")
    )
  }

  # The Cholesky factor, correlation matrix, and partial correlations are the
  # registered 'lkj' deterministic nodes of the primitives.
  node <- .bt_dnode_lkj(
    name = name,
    K = K,
    include_correlation = include_correlation,
    include_primitives = include_primitives
  )
  syntax <- c(syntax, .bt_deterministic_node_emit(node))

  primitive_names <- node$spec$primitive_names
  cpc_names <- node$spec$cpc_names
  primitive_lb <- numeric(0)
  primitive_ub <- numeric(0)
  if(n_pairs > 0L){
    primitive_lb <- stats::setNames(rep(0, n_pairs), primitive_names)
    primitive_ub <- stats::setNames(rep(1, n_pairs), primitive_names)
  }

  monitor <- L_name
  if(include_correlation){
    monitor <- c(monitor, R_name)
  }
  if(include_primitives && n_pairs > 0L){
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

# Correlation matrices R = L L' of lower Cholesky factors of correlation
# matrices, such as .bt_lkj_cholesky_cpc_u_to_L() output. The rows of such a
# factor have unit norm by construction, so the diagonal is exactly 1, as in
# the module's bt_lkj_corr(), instead of the rounded row sums of squares
# (1 +/- a few eps), which would make a constant diagonal vary between draws.
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
  for(draw_i in seq_len(dim(L)[1L])){
    R[draw_i, , ] <- tcrossprod(matrix(L[draw_i, , ], nrow = K, ncol = K))
  }
  for(k in seq_len(K)){
    defined <- is.finite(R[, k, k])
    R[defined, k, k] <- 1
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
  alpha <- .bt_lkj_cholesky_alpha(K = K, eta = eta)
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

.bt_jags_number <- function(x){
  format(x, scientific = FALSE, digits = 17, trim = TRUE)
}

.bt_lkj_cholesky_flat_index <- function(row, column, K){
  (row - 1L) * K + column
}
