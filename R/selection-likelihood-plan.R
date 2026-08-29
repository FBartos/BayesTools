#' Numerical design for a finite-vector selection likelihood
#'
#' @description
#' Creates fixed randomized quasi-Monte Carlo designs for deterministic
#' evaluation of latent-Gaussian selection normalizers. Random shifts are
#' generated locally and persisted in the returned plan, so repeated likelihood
#' evaluations use exactly the same integration points and do not alter R's
#' random-number state.
#'
#' @param block_sizes positive integer selection-block sizes.
#' @param points_per_scramble number of shifted Halton points per scramble.
#' @param scrambles number of independent random shifts.
#' @param seed non-negative integer seed used only to construct fixed shifts.
#' @param relative_tolerance requested relative normalizer error.
#'
#' @return A versioned `BayesTools_selection_likelihood_plan` object containing
#'   fixed integration designs and their lower-triangle covariance ordering.
#'
#' @export
selection_likelihood_plan <- function(
    block_sizes, points_per_scramble = 256L, scrambles = 8L, seed = 1L,
    relative_tolerance = 5e-4){

  check_int(
    block_sizes,
    "block_sizes",
    lower = 1L,
    check_length = 0,
    allow_NA = FALSE
  )
  check_int(
    points_per_scramble,
    "points_per_scramble",
    lower = 8L,
    check_length = 1L,
    allow_NA = FALSE
  )
  check_int(
    scrambles,
    "scrambles",
    lower = 2L,
    check_length = 1L,
    allow_NA = FALSE
  )
  check_int(seed, "seed", lower = 0L, check_length = 1L, allow_NA = FALSE)
  check_real(
    relative_tolerance,
    "relative_tolerance",
    lower = 0,
    upper = 1,
    check_length = 1L,
    allow_NA = FALSE
  )
  if(relative_tolerance == 0){
    stop("'relative_tolerance' must be positive.", call. = FALSE)
  }

  unique_sizes <- sort(unique(as.integer(block_sizes)))
  lower_pairs <- stats::setNames(lapply(unique_sizes, function(block_size){
    .bt_lower_triangle_pairs(seq_len(block_size))
  }), as.character(unique_sizes))
  designs <- stats::setNames(lapply(unique_sizes, function(block_size){
    block_seed <- (
      as.double(seed) + 104729 * as.double(block_size)
    ) %% 4294967296
    .bt_selection_shifted_halton_design(
      dimensions = 2L * block_size,
      points = as.integer(points_per_scramble),
      scrambles = as.integer(scrambles),
      seed = block_seed
    )
  }), as.character(unique_sizes))

  structure(list(
    schema_version = 2L,
    block_sizes = as.integer(block_sizes),
    points_per_scramble = as.integer(points_per_scramble),
    scrambles = as.integer(scrambles),
    relative_tolerance = as.numeric(relative_tolerance),
    seed = as.integer(seed),
    lower_pairs = lower_pairs,
    designs = designs,
    exactness = "E2",
    statistical_target = "finite_vector_product_selection"
  ), class = c("BayesTools_selection_likelihood_plan", "list"))
}


.bt_lower_triangle_pairs <- function(rows){

  row_1 <- integer()
  row_2 <- integer()
  for(column in seq_along(rows)){
    local_rows <- column:length(rows)
    row_1 <- c(row_1, rows[local_rows])
    row_2 <- c(row_2, rep.int(rows[[column]], length(local_rows)))
  }
  data.frame(row_1 = row_1, row_2 = row_2)
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
    shifted[shifted <= 0] <- .Machine$double.eps
    shifted[shifted >= 1] <- 1 - .Machine$double.eps
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
