#' Fixed randomized quasi-Monte Carlo design
#'
#' @description
#' Constructs a deterministic shifted-Halton design in an explicitly supplied
#' integration dimension. Random shifts are generated locally and do not alter
#' R's random-number state. Returned points are strictly inside the open unit
#' hypercube; a design that lands on a boundary is rejected rather than
#' modified.
#'
#' @param dimensions positive integration dimension.
#' @param points number of points per scramble.
#' @param scrambles number of independently shifted designs.
#' @param seed non-negative integer seed used only to construct the shifts.
#'
#' @return Numeric array with dimensions `scrambles`, `points`, and
#'   `dimensions`.
#'
#' @export
selection_qmc_design <- function(dimensions, points, scrambles, seed = 1L){

  check_int(dimensions, "dimensions", lower = 1L, check_length = 1L,
            allow_NA = FALSE)
  check_int(points, "points", lower = 1L, check_length = 1L,
            allow_NA = FALSE)
  check_int(scrambles, "scrambles", lower = 2L, check_length = 1L,
            allow_NA = FALSE)
  check_int(seed, "seed", lower = 0L, check_length = 1L, allow_NA = FALSE)

  .bt_selection_shifted_halton_design(
    dimensions = as.integer(dimensions),
    points     = as.integer(points),
    scrambles  = as.integer(scrambles),
    seed       = as.double(seed)
  )
}


.bt_selection_shifted_halton_design <- function(dimensions, points, scrambles,
                                                 seed){

  bases <- .bt_selection_first_primes(dimensions)
  base_design <- vapply(bases, function(base){
    .bt_selection_radical_inverse(seq_len(points), base)
  }, numeric(points))
  if(dimensions == 1L){
    base_design <- matrix(base_design, ncol = 1L)
  }

  shifts <- .bt_selection_local_uniforms(
    n = scrambles * dimensions,
    seed = seed
  )
  shifts <- matrix(shifts, nrow = scrambles, ncol = dimensions, byrow = TRUE)
  design <- array(NA_real_, dim = c(scrambles, points, dimensions))
  for(scramble in seq_len(scrambles)){
    shifted <- sweep(base_design, 2L, shifts[scramble, ], "+") %% 1
    if(any(!is.finite(shifted)) || any(shifted <= 0) || any(shifted >= 1)){
      stop(
        "The shifted-Halton design reached the boundary of the unit ",
        "hypercube; use a different 'seed'.",
        call. = FALSE
      )
    }
    design[scramble, , ] <- shifted
  }
  design
}


.bt_selection_radical_inverse <- function(index, base){

  value <- numeric(length(index))
  factor <- 1 / base
  remaining <- as.integer(index)
  while(any(remaining > 0L)){
    digit <- remaining %% base
    value <- value + digit * factor
    remaining <- remaining %/% base
    factor <- factor / base
  }
  value
}


.bt_selection_first_primes <- function(n){

  primes <- integer()
  candidate <- 2L
  while(length(primes) < n){
    upper <- floor(sqrt(candidate))
    divisors <- if(upper >= 2L) 2L:upper else integer()
    is_prime <- length(divisors) == 0L || all(candidate %% divisors != 0L)
    if(is_prime){
      primes <- c(primes, candidate)
    }
    candidate <- candidate + 1L
  }
  primes
}


.bt_selection_local_uniforms <- function(n, seed){

  modulus <- 4294967296
  state <- (as.double(seed) + 1) %% modulus
  out <- numeric(n)
  for(i in seq_len(n)){
    state <- (1664525 * state + 1013904223) %% modulus
    out[[i]] <- (state + 0.5) / modulus
  }
  out
}


#' Check Gaussian selection-event support
#'
#' @description
#' Checks whether a Gaussian candidate law can reach a positive-weight event
#' for almost every retained latent context. The complete candidate-plus-retained
#' covariance must be positive definite. Retained covariance itself may be
#' singular. This is a support check, not a numerical integration diagnostic.
#'
#' @param prior a weightfunction prior specifying one publication group's rule.
#' @param candidate_basis finite numeric matrix whose column space is exactly
#'   the candidate variation in that group. Zero columns are allowed, including
#'   a matrix with no columns. Use compiler-declared source geometry, not
#'   posterior covariance estimates or numerical rank truncation.
#'
#' @details Positive weights require no candidate variation. With hard-zero
#'   bins, the best-result rule needs at least one varying result. Two-sided
#'   product selection, and one-sided product selection with a positive
#'   least-significant weight, need every result to vary. One-sided product
#'   selection with a zero least-significant weight instead needs a strictly
#'   positive vector in the candidate column space.
#'
#'   The positive-direction check uses sign-definite columns, then base-R QR
#'   and optimization as candidate finders. A computed direction is accepted
#'   only when its positive entries exceed a standard dot-product roundoff
#'   bound. Failed searches return unknown rather than proving infeasibility.
#'
#'   The caller must check every candidate span having positive prior
#'   probability, including zero-SD and inclusion branches. When only a
#'   guaranteed subspace is supplied, only a positive result is conclusive;
#'   other results cannot classify the complete candidate law. This helper does
#'   not derive scale guarantees or change weights, priors, or source geometry.
#'
#' @return A list with `feasible` (`TRUE`, `FALSE`, or `NA`) and a stable
#'   character `reason`. `FALSE` proves that some positive-probability retained
#'   contexts have zero acceptance probability. `NA` means that a positive
#'   direction could not be certified. Reasons are `positive_weights`,
#'   `varying_best_result`, `varying_two_tails`, `positive_direction`,
#'   `deterministic_candidate`, `deterministic_candidate_rows`,
#'   `opposed_candidate_direction`, or `positive_direction_unavailable`.
#'
#' @export
selection_event_support <- function(prior, candidate_basis){

  if(!is.prior.weightfunction(prior)){
    stop("'prior' must be a weightfunction prior.", call. = FALSE)
  }
  model <- selection_model_spec(prior)
  if(!is.numeric(candidate_basis) || !is.matrix(candidate_basis) ||
     nrow(candidate_basis) < 1L || any(!is.finite(candidate_basis))){
    stop("'candidate_basis' must be a finite numeric matrix with at least one row.",
         call. = FALSE)
  }
  result <- function(feasible, reason){

    list(feasible = feasible, reason = reason)
  }
  weights <- prior$weights
  if(!identical(weights$type, "fixed") || all(weights$omega > 0)){
    return(result(TRUE, "positive_weights"))
  }
  varying <- rowSums(candidate_basis != 0) > 0L
  if(identical(model$weight_rule, "best")){
    return(if(any(varying)) result(TRUE, "varying_best_result") else
      result(FALSE, "deterministic_candidate"))
  }
  if(!all(varying)){
    return(result(FALSE, "deterministic_candidate_rows"))
  }
  if(identical(prior$side, "two-sided") || tail(weights$omega, 1L) > 0){
    return(result(TRUE, "varying_two_tails"))
  }
  basis <- candidate_basis[, colSums(candidate_basis != 0) > 0L, drop = FALSE]
  sign_definite <- vapply(seq_len(ncol(basis)), function(column){
    all(basis[, column] >= 0) || all(basis[, column] <= 0)
  }, logical(1))
  if(all(rowSums(basis[, sign_definite, drop = FALSE] != 0) > 0L)){
    return(result(TRUE, "positive_direction"))
  }
  if(ncol(basis) == 1L){
    return(result(FALSE, "opposed_candidate_direction"))
  }
  for(row in seq_len(nrow(basis) - 1L)){
    if(any(vapply(seq.int(row + 1L, nrow(basis)), function(other){
      all(basis[row, ] == -basis[other, ])
    }, logical(1)))){
      return(result(FALSE, "opposed_candidate_direction"))
    }
  }
  positive_direction <- function(coefficients){

    if(length(coefficients) != ncol(basis) || any(!is.finite(coefficients))){
      return(FALSE)
    }
    direction <- as.vector(basis %*% coefficients)
    magnitude <- as.vector(abs(basis) %*% abs(coefficients))
    product_error <- ncol(basis) * .Machine$double.eps
    gamma <- product_error / (1 - product_error)
    bound <- gamma * magnitude / (1 - gamma)
    all(is.finite(direction)) && all(is.finite(bound)) && all(direction > bound)
  }
  coefficients <- tryCatch(qr.solve(basis, rep(1, nrow(basis))),
                           error = function(e) numeric())
  if(positive_direction(coefficients)){
    return(result(TRUE, "positive_direction"))
  }
  objective <- function(coefficients){

    deficit <- pmin(as.vector(basis %*% coefficients) - 1, 0)
    sum(deficit^2) / 2
  }
  gradient <- function(coefficients){

    deficit <- pmin(as.vector(basis %*% coefficients) - 1, 0)
    as.vector(crossprod(basis, deficit))
  }
  optimum <- tryCatch(stats::optim(rep(0, ncol(basis)), objective, gradient,
                                   method = "BFGS"), error = function(e) NULL)
  if(!is.null(optimum) && positive_direction(optimum$par)){
    return(result(TRUE, "positive_direction"))
  }
  result(NA, "positive_direction_unavailable")
}
