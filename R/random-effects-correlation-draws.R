#' Scalar-structured random-effect correlation draws
#'
#' @description
#' Reconstructs dense correlation matrices for one scalar-structured formula
#' random-effect block from its compiled metadata and posterior draws. The
#' matrices are derived only when this function is called; fitted objects can
#' otherwise retain the compact scalar correlation coordinate.
#'
#' @param random_term compiled random-effect term metadata from a
#'   `BayesTools_formula_design` object.
#' @param posterior_samples posterior sample matrix, data frame, `mcmc`, or
#'   `mcmc.list`. Samples must contain the canonical scalar correlation
#'   coordinate unless it is fixed in `random_term` metadata.
#'
#' @return A dense numeric array with dimensions
#'   `draw x coefficient x coefficient`. Coefficient dimension names use the
#'   compiled random-effect column names when available.
#'
#' @details
#' Supported structures are compound symmetry (`"CS"` and `"HCS"`), discrete
#' autoregressive (`"AR1"` and `"HAR"`), and continuous-time autoregressive
#' (`"CAR"`). Correlation coordinates are transformed and checked against the
#' structure-specific support recorded in `random_term`. `"CAR"` matrices use
#' the compact ordered time coordinates stored by the formula compiler.
#'
#' The returned array is dense and therefore requires memory proportional to
#' the number of posterior draws times the squared coefficient count.
#'
#' @examples
#' dat <- data.frame(
#'   study = factor(rep(c("s1", "s2"), each = 3)),
#'   time = factor(rep(1:3, 2))
#' )
#' generated <- JAGS_formula(
#'   ~ 1 + ar1(time | study),
#'   parameter = "mu",
#'   data = dat,
#'   prior_list = list(intercept = prior("normal", list(0, 1))),
#'   prior_random = prior_random(
#'     sd = prior("point", list(location = 1)),
#'     rho = prior("normal", list(0, 0.5))
#'   )
#' )
#' random_term <- generated$formula_design$random_effects[[1]]
#' posterior <- matrix(c(-0.2, 0.4), ncol = 1)
#' colnames(posterior) <- random_term$correlation$rho_name
#' random_effects_correlation_draws(random_term, posterior)
#'
#' @seealso [prior_random()] [random_effects_marginal_vcov()]
#' @export
random_effects_correlation_draws <- function(random_term, posterior_samples){

  if(!is.list(random_term)){
    stop("'random_term' must be compiled random-effect term metadata.",
         call. = FALSE)
  }
  posterior <- .bt_random_effect_correlation_draws_posterior(posterior_samples)
  context <- "Random-effect correlation reconstruction metadata"
  structure <- .bt_random_effect_structure(random_term, context = context)
  supported <- c("cs", "hcs", "ar1", "har", "car")
  if(!structure %in% supported){
    stop(
      "Random-effect correlation reconstruction supports only scalar ",
      "CS, HCS, AR1, HAR, and CAR structures; block '",
      if(is.null(random_term$block_name)) "unknown" else random_term$block_name,
      "' uses '", toupper(structure), "'.",
      call. = FALSE
    )
  }

  n_columns <- .bt_random_effect_correlation_draws_n_columns(
    random_term = random_term,
    context = context
  )
  coefficient_names <- .bt_random_effect_correlation_draws_column_names(
    random_term = random_term,
    n_columns = n_columns
  )
  output_dimnames <- list(
    draw = rownames(posterior),
    row = coefficient_names,
    column = coefficient_names
  )
  if(n_columns == 1L){
    return(array(
      1,
      dim = c(nrow(posterior), 1L, 1L),
      dimnames = output_dimnames
    ))
  }

  correlation <- .bt_random_effect_correlation_metadata(
    random_term = random_term,
    structure = structure,
    context = context
  )
  if(!is.list(correlation) || !identical(correlation$type, "rho")){
    stop(
      context,
      .bt_random_effect_metadata_block_detail(random_term),
      " does not define canonical scalar correlation metadata.",
      call. = FALSE
    )
  }
  rho <- .bt_random_effect_rho_draws(
    random_term = random_term,
    posterior = posterior,
    missing = "error",
    out_of_support = "error",
    context = context
  )
  distance <- .bt_random_effect_correlation_draws_distance(
    random_term = random_term,
    correlation = correlation,
    structure = structure,
    n_columns = n_columns,
    context = context
  )
  values <- vapply(
    as.vector(distance),
    function(exponent) rho^exponent,
    numeric(length(rho))
  )

  array(
    values,
    dim = c(length(rho), n_columns, n_columns),
    dimnames = output_dimnames
  )
}


# Coerce and validate posterior samples used for correlation reconstruction.
.bt_random_effect_correlation_draws_posterior <- function(posterior_samples){

  supported <- is.matrix(posterior_samples) || is.data.frame(posterior_samples) ||
    inherits(posterior_samples, "mcmc") || inherits(posterior_samples, "mcmc.list")
  if(!supported){
    stop(
      "'posterior_samples' must be a matrix, data frame, mcmc, or mcmc.list object.",
      call. = FALSE
    )
  }

  posterior <- as.matrix(posterior_samples)
  .bt_random_effect_marginal_covariance_validate_posterior(posterior)
}


# Validate the compiled coefficient count used by a random-effect term.
.bt_random_effect_correlation_draws_n_columns <- function(random_term, context){

  n_columns <- random_term$n_columns
  if(!is.numeric(n_columns) || length(n_columns) != 1L || is.na(n_columns) ||
     !is.finite(n_columns) || n_columns < 1 || n_columns != floor(n_columns) ||
     n_columns > .Machine$integer.max){
    stop(
      context,
      .bt_random_effect_metadata_block_detail(random_term),
      " is missing canonical positive integer 'random_term$n_columns'.",
      call. = FALSE
    )
  }

  as.integer(n_columns)
}


# Resolve readable coefficient names without making them a metadata requirement.
.bt_random_effect_correlation_draws_column_names <- function(random_term,
                                                             n_columns){

  column_names <- random_term$column_names
  valid <- is.character(column_names) && length(column_names) == n_columns &&
    !anyNA(column_names) && all(nzchar(column_names)) && !anyDuplicated(column_names)
  if(!isTRUE(valid)){
    column_names <- as.character(seq_len(n_columns))
  }

  column_names
}


# Build the structure-specific exponent matrix after canonical validation.
.bt_random_effect_correlation_draws_distance <- function(
    random_term, correlation, structure, n_columns, context){

  if(structure %in% c("cs", "hcs")){
    return(1 - diag(n_columns))
  }
  if(structure %in% c("ar1", "har")){
    return(abs(outer(seq_len(n_columns), seq_len(n_columns), "-")))
  }

  time_values <- correlation$time_values
  if(is.null(time_values) && is.list(random_term$car)){
    time_values <- random_term$car$time_values
  }
  valid_time <- is.numeric(time_values) && length(time_values) == n_columns &&
    !anyNA(time_values) && all(is.finite(time_values)) &&
    !anyDuplicated(time_values) && all(diff(time_values) > 0)
  if(!isTRUE(valid_time)){
    stop(
      context,
      .bt_random_effect_metadata_block_detail(random_term),
      " is missing canonical ordered CAR time coordinates.",
      call. = FALSE
    )
  }

  .bt_random_effect_validate_car_distance_matrix(
    abs(outer(time_values, time_values, "-")),
    n_columns
  )
}
