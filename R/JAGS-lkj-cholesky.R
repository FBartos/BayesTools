#' Generate a JAGS LKJ-Cholesky correlation module
#'
#' @description
#' Creates JAGS syntax for a Cholesky factor of a correlation matrix with the
#' same target distribution as Stan's `lkj_corr_cholesky(eta)`. The default
#' backend uses the BayesTools compiled JAGS module. The syntax backend keeps a
#' pure-JAGS fallback that samples canonical partial correlations through
#' independent beta nodes and deterministically builds the lower Cholesky factor.
#'
#' @param name character scalar. Prefix used for generated JAGS nodes.
#' @param K integer scalar. Dimension of the correlation matrix.
#' @param eta positive numeric scalar. LKJ concentration parameter.
#' @param include_correlation logical scalar. Whether to generate the
#'   deterministic correlation matrix `name_R`.
#' @param include_primitives logical scalar. Whether primitive beta/CPC nodes
#'   should be included in the returned monitor vector.
#' @param backend character scalar. Backend used for generated JAGS syntax.
#'   `"module"` uses the package-shipped compiled JAGS module, and `"syntax"`
#'   emits a pure-JAGS fallback.
#'
#' @return A list with JAGS syntax, monitor names, primitive bridge coordinate
#'   names and bounds, and metadata for the generated LKJ-Cholesky block.
#'
#' @export
JAGS_lkj_corr_cholesky <- function(name, K, eta = 1,
                                   include_correlation = TRUE,
                                   include_primitives = FALSE,
                                   backend = c("module", "syntax")){

  backend <- match.arg(backend)
  .bt_check_lkj_cholesky_module_inputs(name, K, eta, include_correlation, include_primitives, backend)

  if(backend == "module"){
    return(.bt_JAGS_lkj_corr_cholesky_module(
      name = name,
      K = K,
      eta = eta,
      include_correlation = include_correlation,
      include_primitives = include_primitives
    ))
  }

  .bt_JAGS_lkj_corr_cholesky_syntax(
    name = name,
    K = K,
    eta = eta,
    include_correlation = include_correlation,
    include_primitives = include_primitives
  )
}

.bt_JAGS_lkj_corr_cholesky_syntax <- function(name, K, eta,
                                              include_correlation,
                                              include_primitives){

  pairs <- .bt_lkj_cholesky_cpc_pairs(K = K, eta = eta)
  n_pairs <- nrow(pairs)
  pair_index <- matrix(NA_integer_, K, K)
  if(n_pairs > 0L){
    for(p in seq_len(n_pairs)){
      pair_index[pairs$i[p], pairs$j[p]] <- pairs$index[p]
    }
  }

  L_name <- paste0(name, "_L")
  R_name <- paste0(name, "_R")
  u_name <- paste0(name, "_lkj_u")
  cpc_name <- paste0(name, "_lkj_cpc")

  syntax <- c(paste0("# LKJ-Cholesky correlation module: ", name))

  if(n_pairs > 0L){
    for(p in seq_len(n_pairs)){
      alpha <- .bt_jags_number(pairs$alpha[p])
      syntax <- c(
        syntax,
        paste0(u_name, "[", p, "] ~ dbeta(", alpha, ", ", alpha, ")"),
        paste0(cpc_name, "[", p, "] <- 2 * ", u_name, "[", p, "] - 1")
      )
    }
  }

  for(row in seq_len(K)){
    for(column in seq_len(K)){
      target <- paste0(L_name, "[", row, ",", column, "]")
      if(column > row){
        syntax <- c(syntax, paste0(target, " <- 0"))
      }else if(row == 1L && column == 1L){
        syntax <- c(syntax, paste0(target, " <- 1"))
      }else if(column < row){
        syntax <- c(
          syntax,
          paste0(target, " <- ", .bt_lkj_cholesky_entry_expr(cpc_name, pair_index, row, column))
        )
      }else{
        syntax <- c(
          syntax,
          paste0(target, " <- ", .bt_lkj_cholesky_diag_expr(cpc_name, pair_index, row))
        )
      }
    }
  }

  if(include_correlation){
    for(row in seq_len(K)){
      for(column in seq_len(row)){
        target <- paste0(R_name, "[", row, ",", column, "]")
        syntax <- c(
          syntax,
          paste0(target, " <- inprod(", L_name, "[", row, ",1:", K, "], ", L_name, "[", column, ",1:", K, "])")
        )
        if(column != row){
          syntax <- c(
            syntax,
            paste0(R_name, "[", column, ",", row, "] <- ", target)
          )
        }
      }
    }
  }

  primitive_names <- character(0)
  cpc_names <- character(0)
  primitive_lb <- numeric(0)
  primitive_ub <- numeric(0)
  if(n_pairs > 0L){
    primitive_names <- paste0(u_name, "[", seq_len(n_pairs), "]")
    cpc_names <- paste0(cpc_name, "[", seq_len(n_pairs), "]")
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
    cpc_names = cpc_names,
    primitive_bounds = list(lb = primitive_lb, ub = primitive_ub),
    pairs = pairs,
    backend = "syntax",
    jags_module = NULL,
    required_packages = NULL,
    data = list()
  )

  class(out) <- c("BayesTools_JAGS_lkj_corr_cholesky", "list")
  out
}

.bt_JAGS_lkj_corr_cholesky_module <- function(name, K, eta,
                                              include_correlation,
                                              include_primitives){

  pairs <- .bt_lkj_cholesky_cpc_pairs(K = K, eta = eta)
  n_pairs <- nrow(pairs)

  L_name <- paste0(name, "_L")
  R_name <- paste0(name, "_R")
  L_flat_name <- paste0(name, "_L_flat")
  R_flat_name <- paste0(name, "_R_flat")
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
      paste0(u_name, "[1:", n_pairs, "] ~ dbt_lkj_cpc(", alpha_name, ")"),
      paste0(L_flat_name, "[1:", K * K, "] <- bt_lkj_cholesky(", u_name, ", ", K, ")")
    )
    if(include_correlation){
      syntax <- c(
        syntax,
        paste0(R_flat_name, "[1:", K * K, "] <- bt_lkj_corr(", u_name, ", ", K, ")")
      )
    }
  }

  for(row in seq_len(K)){
    for(column in seq_len(K)){
      target <- paste0(L_name, "[", row, ",", column, "]")
      if(K == 1L){
        syntax <- c(syntax, paste0(target, " <- 1"))
      }else{
        flat_index <- .bt_lkj_cholesky_flat_index(row, column, K)
        syntax <- c(syntax, paste0(target, " <- ", L_flat_name, "[", flat_index, "]"))
      }
    }
  }

  if(include_correlation){
    for(row in seq_len(K)){
      for(column in seq_len(K)){
        target <- paste0(R_name, "[", row, ",", column, "]")
        if(K == 1L){
          syntax <- c(syntax, paste0(target, " <- 1"))
        }else{
          flat_index <- .bt_lkj_cholesky_flat_index(row, column, K)
          syntax <- c(syntax, paste0(target, " <- ", R_flat_name, "[", flat_index, "]"))
        }
      }
    }
  }

  primitive_names <- character(0)
  cpc_names <- character(0)
  primitive_lb <- numeric(0)
  primitive_ub <- numeric(0)
  if(n_pairs > 0L){
    primitive_names <- paste0(u_name, "[", seq_len(n_pairs), "]")
    primitive_lb <- stats::setNames(rep(0, n_pairs), primitive_names)
    primitive_ub <- stats::setNames(rep(1, n_pairs), primitive_names)
    if(include_primitives){
      cpc_names <- paste0(cpc_name, "[", seq_len(n_pairs), "]")
      for(p in seq_len(n_pairs)){
        syntax <- c(syntax, paste0(cpc_name, "[", p, "] <- 2 * ", u_name, "[", p, "] - 1"))
      }
    }
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
                                                 include_primitives,
                                                 backend = c("module", "syntax")){

  backend <- match.arg(backend)
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
  check_char(backend, "backend", allow_values = c("module", "syntax"), allow_NA = FALSE)

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

.bt_lkj_cholesky_cpc_to_L <- function(cpc, K){

  K <- .bt_lkj_cholesky_check_K(K)
  n_pairs <- .bt_lkj_cholesky_n_pairs(K)
  if(length(cpc) != n_pairs){
    stop("'cpc' must have length K * (K - 1) / 2.", call. = FALSE)
  }
  if(n_pairs > 0L && (any(!is.finite(cpc)) || any(abs(cpc) >= 1))){
    stop("'cpc' values must be finite and strictly between -1 and 1.", call. = FALSE)
  }

  .bt_lkj_cholesky_cpc_u_to_L((cpc + 1) / 2, K = K)
}

.bt_lkj_cholesky_cpc_u_to_L <- function(u, K){

  K <- .bt_lkj_cholesky_check_K(K)
  n_pairs <- .bt_lkj_cholesky_n_pairs(K)
  .bt_lkj_cholesky_check_u_shape(u, n_pairs)
  .bt_lkj_cholesky_check_u_support(u)

  .BayesTools_require_native_lkj()
  .Call("BayesTools_lkj_cholesky_from_u", u, as.integer(K), PACKAGE = "BayesTools")
}

.bt_lkj_cholesky_rng <- function(n, K, eta = 1){

  if(!is.numeric(n) || length(n) != 1L || is.na(n) || n != as.integer(n) || n < 1L){
    stop("'n' must be a positive integer scalar.", call. = FALSE)
  }
  K <- .bt_lkj_cholesky_check_K(K)
  eta <- .bt_lkj_cholesky_check_eta(eta)
  alpha <- .bt_lkj_cholesky_alpha(K = K, eta = eta)
  n_pairs <- length(alpha)

  if(n_pairs == 0L){
    draws <- array(1, dim = c(n, 1L, 1L))
    return(draws)
  }

  u <- matrix(NA_real_, nrow = n, ncol = n_pairs)
  for(p in seq_len(n_pairs)){
    u[, p] <- stats::rbeta(n, alpha[p], alpha[p])
  }

  .bt_lkj_cholesky_cpc_u_to_L(u, K = K)
}

.bt_lkj_cholesky_cpc_u_to_R <- function(u, K){

  K <- .bt_lkj_cholesky_check_K(K)
  n_pairs <- .bt_lkj_cholesky_n_pairs(K)
  .bt_lkj_cholesky_check_u_shape(u, n_pairs)
  .bt_lkj_cholesky_check_u_support(u)

  .BayesTools_require_native_lkj()
  .Call("BayesTools_lkj_corr_from_u", u, as.integer(K), PACKAGE = "BayesTools")
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

.bt_lkj_cholesky_corr <- function(L){

  if(length(dim(L)) == 2L){
    return(L %*% t(L))
  }
  if(length(dim(L)) == 3L){
    out <- array(NA_real_, dim = dim(L))
    for(i in seq_len(dim(L)[1L])){
      out[i, , ] <- L[i, , ] %*% t(L[i, , ])
    }
    return(out)
  }
  stop("'L' must be a matrix or a draw-by-row-by-column array.", call. = FALSE)
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

.bt_lkj_cholesky_entry_expr <- function(cpc_name, pair_index, row, column){

  factors <- character(0)
  if(column > 1L){
    for(m in 1:(column - 1L)){
      factors <- c(factors, .bt_lkj_cholesky_sqrt_expr(cpc_name, pair_index[m, row]))
    }
  }
  factors <- c(paste0(cpc_name, "[", pair_index[column, row], "]"), factors)
  paste(factors, collapse = " * ")
}

.bt_lkj_cholesky_diag_expr <- function(cpc_name, pair_index, row){

  if(row == 1L){
    return("1")
  }
  factors <- character(0)
  for(m in 1:(row - 1L)){
    factors <- c(factors, .bt_lkj_cholesky_sqrt_expr(cpc_name, pair_index[m, row]))
  }
  paste(factors, collapse = " * ")
}

.bt_lkj_cholesky_sqrt_expr <- function(cpc_name, index){
  paste0("sqrt(1 - pow(", cpc_name, "[", index, "], 2))")
}

.bt_jags_number <- function(x){
  format(x, scientific = FALSE, digits = 17, trim = TRUE)
}

.bt_lkj_cholesky_flat_index <- function(row, column, K){
  (row - 1L) * K + column
}
