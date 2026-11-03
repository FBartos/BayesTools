#' @title BayesTools Contrast Matrices
#'
#' @description BayesTools provides several contrast matrix functions for Bayesian factor analysis.
#' These functions create different types of contrast matrices that can be used with factor
#' variables in Bayesian models.
#'
#' @details
#' The package includes the following contrast functions:
#' \describe{
#'   \item{\code{contr.orthonormal}}{Return a matrix of orthonormal contrasts: a normalised
#'     Helmert basis, whose columns have unit length and sum to zero and which is identical on
#'     every platform. Any orthonormal basis of the centred space implies the same exchangeable
#'     prior on the factor level effects under a spherical prior on the coefficients, which
#'     corresponds to the description by \insertCite{rouder2012default;textual}{BayesTools}.
#'     Returns a matrix with n rows and k columns, with k = n - 1 if \code{contrasts = TRUE}
#'     and k = n if \code{contrasts = FALSE}.}
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

.bt_contrast_level_count <- function(n, contrasts, minimum = 2L){

  check_bool(contrasts, "contrasts", allow_NA = FALSE)

  if(length(n) > 1L){
    return(length(n))
  }
  if(length(n) == 1L && !is.numeric(n) && minimum <= 1L){
    return(1L)
  }
  if(!is.numeric(n) || length(n) != 1L){
    stop("Not enough degrees of freedom to define contrasts.", call. = FALSE)
  }

  check_int(n, "n", allow_NA = FALSE)
  if(n < minimum){
    stop("Not enough degrees of freedom to define contrasts.", call. = FALSE)
  }

  as.integer(n)
}

# The normalised Helmert basis: column j holds -1/d in rows 1..j, j/d in row
# j + 1 and 0 below, with d = sqrt(j * (j + 1)). Its columns are orthonormal and
# sum to zero, so B B' = I - J/n. Every entry is one correctly rounded IEEE
# operation, so the matrix is bit-identical on every platform (a LAPACK
# eigenbasis of the same space is not).
.bt_helmert_contrast_basis <- function(n){

  basis <- matrix(0, nrow = n, ncol = n - 1L)
  for(j in seq_len(n - 1L)){
    d <- sqrt(j * (j + 1))
    basis[seq_len(j), j] <- -1 / d
    basis[j + 1L, j]     <- j / d
  }

  basis
}

#' @rdname contr.BayesTools
#' @export
contr.orthonormal <- function(n, contrasts = TRUE){

  n <- .bt_contrast_level_count(n, contrasts)

  cont <- diag(n)
  if(contrasts){
    cont <- .bt_helmert_contrast_basis(n)
  }

  return(cont)
}

#' @rdname contr.BayesTools
#' @export
contr.meandif <- function(n, contrasts = TRUE){

  n <- .bt_contrast_level_count(n, contrasts)

  cont <- diag(n)
  if(contrasts){
    cont <- .bt_helmert_contrast_basis(n) / (sqrt(1 - 1/n))
  }

  return(cont)
}


#' @rdname contr.BayesTools
#' @export
contr.independent <- function(n, contrasts = TRUE){

  n <- .bt_contrast_level_count(n, contrasts, minimum = 1L)

  cont <- diag(x = 1, nrow = n, ncol = n)

  return(cont)
}

#' @rdname contr.BayesTools
#' @export
contr.ordered_cumulative <- function(n, contrasts = TRUE){

  n <- .bt_contrast_level_count(n, contrasts)
  if(!contrasts){
    return(diag(n))
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

  n <- .bt_contrast_level_count(n, contrasts, minimum = 1L)
  if(!contrasts){
    return(diag(n))
  }

  cont <- matrix(0, nrow = n, ncol = n)
  for(i in seq_len(n)){
    cont[i, seq_len(i)] <- 1
  }

  cont
}
