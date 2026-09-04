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
