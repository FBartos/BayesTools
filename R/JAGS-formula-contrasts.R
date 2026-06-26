#' @title BayesTools Contrast Matrices
#'
#' @description BayesTools provides several contrast matrix functions for Bayesian factor analysis.
#' These functions create different types of contrast matrices that can be used with factor
#' variables in Bayesian models.
#'
#' @details
#' The package includes the following contrast functions:
#' \describe{
#'   \item{\code{contr.orthonormal}}{Return a matrix of orthonormal contrasts.
#'     Code is based on \code{stanova::contr.bayes} and corresponding to description
#'     by \insertCite{rouder2012default;textual}{BayesTools}. Returns a matrix with n rows and
#'     k columns, with k = n - 1 if \code{contrasts = TRUE} and k = n if \code{contrasts = FALSE}.}
#'   \item{\code{contr.meandif}}{Return a matrix of mean difference contrasts.
#'     This is an adjustment to the \code{contr.orthonormal} that ascertains that the prior
#'     distributions on difference between the gran mean and factor level are identical independent
#'     of the number of factor levels (which does not hold for the orthonormal contrast). Furthermore,
#'     the contrast is re-scaled so the specified prior distribution exactly corresponds to the prior
#'     distribution on difference between each factor level and the grand mean -- this is approximately
#'     twice the scale of \code{contr.orthonormal}. Returns a matrix with n rows and k columns,
#'     with k = n - 1 if \code{contrasts = TRUE} and k = n if \code{contrasts = FALSE}.}
#'   \item{\code{contr.independent}}{Return a matrix of independent contrasts -- a level for each term.
#'     Returns a matrix with n rows and k columns, with k = n if \code{contrasts = TRUE} and k = n
#'     if \code{contrasts = FALSE}.}
#'   \item{\code{contr.ordered_cumulative}}{Return a cumulative increment basis
#'     for ordered factors, with the first level fixed to zero and the last level
#'     equal to the full total effect.}
#'   \item{\code{contr.ordered_cumulative_levels}}{Return a cumulative level
#'     basis for ordered factors, where the first level receives the first
#'     allocation share and the last level equals the full total effect.}
#' }
#'
#' @param n a vector of levels for a factor, or the number of levels
#' @param contrasts logical indicating whether contrasts should be computed
#'
#' @examples
#' # Orthonormal contrasts
#' contr.orthonormal(c(1, 2))
#' contr.orthonormal(c(1, 2, 3))
#'
#' # Mean difference contrasts
#' contr.meandif(c(1, 2))
#' contr.meandif(c(1, 2, 3))
#'
#' # Independent contrasts
#' contr.independent(c(1, 2))
#' contr.independent(c(1, 2, 3))
#'
#' # Ordered cumulative contrasts
#' contr.ordered_cumulative(c("low", "mid", "high"))
#' contr.ordered_cumulative_levels(c("low", "mid", "high"))
#'
#' @references
#' \insertAllCited{}
#'
#' @aliases contr.orthonormal contr.meandif contr.independent
#'   contr.ordered_cumulative contr.ordered_cumulative_levels
#' @name contr.BayesTools
NULL

#' @rdname contr.BayesTools
#' @export
contr.orthonormal <- function(n, contrasts = TRUE){
  # based on: stanova::contr.bayes
  if(length(n) <= 1L){
    if(is.numeric(n) && length(n) == 1L && n > 1L){
      return(TRUE)
    }else{
      stop("Not enough degrees of freedom to define contrasts.")
    }
  }else{
    n <- length(n)
  }

  cont <- diag(n)
  if(contrasts){
    a       <- n
    I_a     <- diag(a)
    J_a     <- matrix(1, nrow = a, ncol = a)
    Sigma_a <- I_a - J_a/a
    cont    <- eigen(Sigma_a)$vectors[, seq_len(a - 1), drop = FALSE]
  }

  return(cont)
}

#' @rdname contr.BayesTools
#' @export
contr.meandif <- function(n, contrasts = TRUE){

  if(length(n) <= 1L){
    if(is.numeric(n) && length(n) == 1L && n > 1L){
      return(TRUE)
    }else{
      stop("Not enough degrees of freedom to define contrasts.")
    }
  }else{
    n <- length(n)
  }

  cont <- diag(n)
  if(contrasts){
    a       <- n
    I_a     <- diag(a)
    J_a     <- matrix(1, nrow = a, ncol = a)
    Sigma_a <- I_a - J_a/a
    cont    <- eigen(Sigma_a)$vectors[, seq_len(a - 1), drop = FALSE]
    cont    <- cont / (sqrt(1 - 1/n))
  }

  return(cont)
}


#' @rdname contr.BayesTools
#' @export
contr.independent <- function(n, contrasts = TRUE){

  if(length(n) <= 1L){
    if(is.numeric(n) && length(n) == 1L && n >= 1L){
      return(TRUE)
    }else{
      stop("Not enough degrees of freedom to define contrasts.")
    }
  }else{
    n <- length(n)
  }

  cont <- diag(x = 1, nrow = n, ncol = n)

  return(cont)
}

#' @rdname contr.BayesTools
#' @export
contr.ordered_cumulative <- function(n, contrasts = TRUE){

  if(length(n) <= 1L){
    if(is.numeric(n) && length(n) == 1L && n > 1L){
      return(TRUE)
    }else{
      stop("Not enough degrees of freedom to define contrasts.")
    }
  }else{
    n <- length(n)
  }

  cont <- matrix(0, nrow = n, ncol = n - 1L)
  for(i in seq_len(n)){
    if(i > 1L){
      cont[i, seq_len(i - 1L)] <- 1
    }
  }

  cont
}

#' @rdname contr.BayesTools
#' @export
contr.ordered_cumulative_levels <- function(n, contrasts = TRUE){

  if(length(n) <= 1L){
    if(is.numeric(n) && length(n) == 1L && n >= 1L){
      return(TRUE)
    }else{
      stop("Not enough degrees of freedom to define contrasts.")
    }
  }else{
    n <- length(n)
  }

  cont <- matrix(0, nrow = n, ncol = n)
  for(i in seq_len(n)){
    cont[i, seq_len(i)] <- 1
  }

  cont
}

