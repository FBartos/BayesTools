#' @title Harrell-Davis quantiles of posterior draws
#'
#' @description Estimates quantiles of equally weighted draws with the
#' Harrell-Davis estimator \insertCite{harrell1982new}{BayesTools}, a weighted
#' average of all order statistics. Applied pointwise to the columns of a
#' draws-by-grid matrix, it gives smooth credible bands: the kinks that the
#' empirical quantile (`stats::quantile()`) shows when the draws that define the
#' quantile change from one grid point to the next are attenuated. The
#' estimand is the same as that of the empirical quantile; in simulations of
#' smooth curves the accuracy was similar (see the details).
#'
#' @param x a numeric vector of draws, or a numeric matrix with the draws in
#' rows and one quantity (e.g., a grid point) per column. Missing values are
#' not allowed.
#' @param probs probabilities in `[0, 1]`, in any order; duplicates are
#' allowed.
#' @param names whether the result is named by the probabilities (as
#' `stats::quantile()` names them). Defaults to `FALSE`.
#'
#' @details For `m` draws with order statistics `y(1), ..., y(m)` and a
#' probability `p` in `(0, 1)`, the estimate is the sum of `w_i * y(i)` with
#' weights `w_i = I_(i/m)(a, b) - I_((i-1)/m)(a, b)`, where `I` is the
#' regularized incomplete beta function, `a = p * (m + 1)`, and
#' `b = (1 - p) * (m + 1)`. The weights depend only on the number of draws and
#' `p` and are computed once for all columns that need them. The estimate is
#' exact for constant columns and does not overflow for finite draws. The
#' probabilities `0` and `1` return the minimum and the maximum, and a single
#' draw is returned as it is.
#'
#' Caveats of the estimator:
#' * The quantiles of the columns are pointwise: the result is not a
#'   simultaneous band of the underlying curve.
#' * The estimate is smooth, not more accurate. The slowly varying Monte Carlo
#'   error of the empirical quantile across neighboring columns remains, and
#'   only more draws reduce it. Extreme probabilities (e.g., `0.0005`) are not
#'   estimated more accurately than by the empirical quantile. In simulations
#'   of smooth curves (lines with normal coefficients, probability `0.025`,
#'   1,000 and 10,000 draws), the root mean squared error was about as large as
#'   that of the empirical quantile (median ratios of 0.92 and 0.97) and the
#'   roughness of the error along the grid was several times smaller; this is
#'   no guarantee for other draws or settings.
#' * The input draws have equal weights (there are no importance weights), but
#'   the estimator gives every draw a positive and unequal weight that depends
#'   on its rank and is largest for the draws closest to the target rank. A
#'   small number of draws (below roughly 200) and heavy tails can therefore
#'   move the estimate far from the empirical quantile.
#' * Draws with a point mass (an atom, e.g., a parameter fixed at zero in part
#'   of a model-averaged posterior) are not treated separately: within about
#'   two Monte Carlo standard errors (in probability) of the edge of a point
#'   mass, the estimate blends the point mass and the continuous draws.
#' * The weighted sum is computed in double precision. Its rounding error is
#'   of the order of the machine precision times the largest absolute draw.
#'   Where the sorted draws jump close to the target rank (e.g., at the edge of
#'   a point mass or in a sparse tail), the rounding of the cell boundaries
#'   `i/m` can make it larger by a factor that grows roughly as
#'   `sqrt(m * p / (1 - p))` (e.g., about 40 with 10,000 draws and
#'   `p = 0.975`, and about 1,400 with 100,000 draws and `p = 0.9995`, for
#'   draws of 0 and 1). Differences between draws that are much smaller than
#'   this error (e.g., the middle draw of `c(-1e16, 1, 1e16)`) are not
#'   resolved.
#'
#' An infinite draw has a positive weight at every probability strictly
#' between 0 and 1, which makes the weighted average infinite (or undefined for
#' infinite draws of both signs). A column that contains an infinite draw is
#' therefore summarized by the empirical quantile of that column
#' (`stats::quantile()`, type 7). If that quantile is `NaN` (an interpolation
#' between draws of `-Inf` and `Inf`), or if the weights of a probability
#' cannot be computed numerically (a probability below about `1e-310`), the
#' function stops with an error of class `BayesTools_harrell_davis_undefined`
#' (also of the family class `BayesTools_harrell_davis`). The weights are
#' computed only when a column with finite draws needs them, so columns of
#' infinite draws never raise the second error.
#'
#' @return a numeric vector of the same length as `probs` for a vector `x`,
#' and a matrix with `length(probs)` rows and one column per column of `x`
#' (keeping its column names) for a matrix `x`.
#'
#' @examples
#' set.seed(1)
#' draws <- rnorm(500)
#' harrell_davis_quantile(draws, c(.025, .5, .975))
#' quantile(draws, c(.025, .5, .975), names = FALSE)
#'
#' # pointwise credible band of a line from draws of its intercept and slope
#' intercept <- rnorm(500)
#' slope     <- rnorm(500, sd = .5)
#' x_grid    <- seq(-1, 1, length.out = 5)
#' line      <- outer(intercept, rep(1, 5)) + outer(slope, x_grid)
#' harrell_davis_quantile(line, c(.025, .5, .975), names = TRUE)
#'
#' @references
#' \insertAllCited{}
#'
#' @seealso [stats::quantile()]
#' @export
harrell_davis_quantile <- function(x, probs, names = FALSE){

  is_factor <- is.factor(x)
  x <- unclass(x)
  if(is_factor || !is.numeric(x) || !(is.null(dim(x)) || length(dim(x)) == 2L))
    stop("The 'x' argument must be a numeric vector or a numeric matrix with the draws in rows.", call. = FALSE)
  check_real(probs, "probs", lower = 0, upper = 1, check_length = 0, allow_NULL = TRUE, allow_NA = FALSE)
  check_bool(names, "names")
  probs <- as.numeric(probs)

  is_matrix <- is.matrix(x)
  if(!is_matrix){
    x <- matrix(as.numeric(x), ncol = 1L)
  }else if(!is.double(x)){
    storage.mode(x) <- "double"
  }
  n_draws <- nrow(x)
  n_cols  <- ncol(x)
  if(n_draws == 0L)
    stop("The 'x' argument must contain at least one draw.", call. = FALSE)
  if(anyNA(x))
    stop("The 'x' argument cannot contain NA/NaN values.", call. = FALSE)

  n_probs <- length(probs)
  result  <- matrix(NA_real_, nrow = n_probs, ncol = n_cols)

  if(n_probs > 0L && n_cols > 0L){
    if(n_draws == 1L){
      result[] <- rep(x[1L, ], each = n_probs)
    }else{
      probs_unique <- unique(probs)
      match_probs  <- match(probs, probs_unique)
      is_lower     <- probs_unique <= 0
      is_upper     <- probs_unique >= 1
      interior     <- which(!is_lower & !is_upper)
      anchors      <- pmin(pmax(round(probs_unique[interior] * (n_draws + 1)), 1), n_draws)
      weights      <- NULL

      for(j in seq_len(n_cols)){
        draws <- sort.int(x[, j], method = "radix")
        if(is.infinite(draws[1L]) || is.infinite(draws[n_draws])){
          # the weighted average is not defined with infinite draws
          result[, j] <- stats::quantile(draws, probs = probs, names = FALSE, type = 7)
          next
        }
        if(is.null(weights)){
          # only columns with finite draws need the weights
          weights <- lapply(probs_unique[interior], .harrell_davis_weights, m = n_draws)
        }
        values           <- numeric(length(probs_unique))
        values[is_lower] <- draws[1L]
        values[is_upper] <- draws[n_draws]
        values[interior] <- .harrell_davis_sorted(draws, weights, anchors)
        result[, j]      <- values[match_probs]
      }
    }
  }

  if(anyNA(result)){
    bad_column <- which(colSums(is.na(result)) > 0)[1L]
    .harrell_davis_stop_undefined(paste0(
      "the empirical quantile of ",
      if(is_matrix) paste0("column ", bad_column, " of 'x'") else "'x'",
      " is NaN because its infinite draws have opposite signs. Remove or recode the infinite draws."
    ))
  }

  if(names && n_probs > 0L){
    rownames(result) <- base::names(stats::quantile(0, probs = probs, names = TRUE))
  }
  if(!is_matrix){
    return(result[, 1L])
  }
  colnames(result) <- colnames(x)

  return(result)
}

# Weighted averages of the sorted draws in the shift form: the weighted
# differences from an anchor draw (the draw at the rounded rank of the
# probability) are added to the anchor, so that constant columns are exact.
# Draws whose range overflows are scaled by their largest absolute value.
.harrell_davis_sorted <- function(draws, weights, anchors){

  n_draws <- length(draws)
  scale   <- 1
  if(!is.finite(draws[n_draws] - draws[1L])){
    scale <- max(abs(draws[1L]), abs(draws[n_draws]))
    draws <- draws / scale
  }

  values <- numeric(length(anchors))
  for(i in seq_along(anchors)){
    anchor    <- draws[anchors[i]]
    values[i] <- anchor + sum(weights[[i]] * (draws - anchor))
  }

  return(values * scale)
}

# Harrell-Davis weights of m sorted draws at the probability p in (0, 1):
# differences of the Beta(p * (m + 1), (1 - p) * (m + 1)) distribution function
# over the cells ((i - 1) / m, i / m]. The cells up to the mode use the lower
# tail and the remaining cells the upper tail, so that the small weights in
# either tail keep their relative accuracy.
.harrell_davis_weights <- function(p, m){

  a     <- p * (m + 1)
  b     <- (1 - p) * (m + 1)
  grid  <- (0:m) / m
  mode  <- min(max((a - 1) / (m - 1), 0), 1)
  lower <- min(sum(grid <= mode), m)

  weights <- numeric(m)
  suppressWarnings({
    weights[seq_len(lower)] <- diff(stats::pbeta(grid[seq_len(lower + 1L)], a, b))
    if(lower < m){
      weights[(lower + 1L):m] <- -diff(stats::pbeta(grid[(lower + 1L):(m + 1L)], a, b, lower.tail = FALSE))
    }
  })
  if(!all(is.finite(weights))){
    .harrell_davis_stop_undefined(paste0(
      "the weights of the probability ", format(p), " cannot be computed numerically. ",
      "Use values of 'probs' further from 0 and 1."
    ))
  }

  return(weights)
}

# The error of a quantile that is not available (class
# BayesTools_harrell_davis_undefined; callers match the class, not the message).
.harrell_davis_stop_undefined <- function(reason){

  stop(structure(
    class = c("BayesTools_harrell_davis_undefined", "BayesTools_harrell_davis", "error", "condition"),
    list(message = paste0("The Harrell-Davis quantile is unavailable: ", reason), call = NULL)
  ))
}
