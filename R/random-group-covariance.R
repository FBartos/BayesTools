# Known group-level covariance kernels for formula random effects.

#' Known group covariance for random effects
#'
#' @description
#' Creates a fixed group-level covariance kernel for use with
#' [random_effects_formula()]. The kernel is a random-effect design property:
#' it describes covariance across grouping levels and is separate from
#' column-level random-effect covariance structures such as
#' [random_covariance()].
#'
#' @param covariance numeric square matrix with row and column names identifying
#'   grouping levels.
#' @param scale how to scale the kernel before fitting. `"cor"` converts the
#'   reordered matrix to a correlation matrix, `"none"` uses the matrix as
#'   supplied, `"cor0"` applies `(stats::cov2cor(K) - min(K)) / (1 - min(K))`,
#'   and `"cov0"` applies `K - min(K)`.
#'
#' @return A `random_group_covariance` object.
#'
#' @seealso [random_effects_formula()]
#'
#' @export
random_group_covariance <- function(covariance,
                                    scale = c("cor", "none", "cor0", "cov0")){

  scale <- match.arg(scale)
  covariance <- .bt_random_group_covariance_matrix(
    covariance,
    name = "covariance"
  )
  triangle_source <- attr(covariance, "triangle_source", exact = TRUE)

  out <- list(
    covariance = covariance,
    scale = scale,
    triangle_source = triangle_source
  )
  class(out) <- c("random_group_covariance", "list")
  out
}

#' @export
print.random_group_covariance <- function(x, ...){

  cat("random_group_covariance()\n")
  cat("  scale: ", x$scale, "\n", sep = "")
  cat("  levels: ", nrow(x$covariance), "\n", sep = "")
  invisible(x)
}

.bt_is_random_group_covariance <- function(x){

  inherits(x, "random_group_covariance")
}

.bt_check_random_group_covariance <- function(x, allow_NULL = FALSE){

  if(is.null(x) && isTRUE(allow_NULL)){
    return(invisible(TRUE))
  }
  if(!.bt_is_random_group_covariance(x)){
    stop(
      "'group_covariance' must be created by random_group_covariance().",
      call. = FALSE
    )
  }
  .bt_random_group_covariance_matrix(x$covariance, name = "group_covariance$covariance")
  check_char(
    x$scale,
    "group_covariance$scale",
    allow_values = c("cor", "none", "cor0", "cov0")
  )

  invisible(TRUE)
}

.bt_as_random_group_covariance <- function(x){

  if(.bt_is_random_group_covariance(x)){
    .bt_check_random_group_covariance(x)
    return(x)
  }
  if(is.matrix(x) || is.data.frame(x)){
    return(random_group_covariance(x, scale = "cor"))
  }

  stop(
    "'group_covariance' entries must be random_group_covariance() objects or numeric matrices.",
    call. = FALSE
  )
}

.bt_random_group_covariance_matrix <- function(covariance, name = "covariance"){

  if(is.data.frame(covariance)){
    covariance <- as.matrix(covariance)
  }
  if(!is.matrix(covariance)){
    covariance <- tryCatch(
      as.matrix(covariance),
      error = function(e) NULL
    )
  }
  if(!is.matrix(covariance) || !is.numeric(covariance)){
    stop("'", name, "' must be coercible to a numeric matrix.", call. = FALSE)
  }
  if(nrow(covariance) != ncol(covariance)){
    stop("'", name, "' must be a square matrix.", call. = FALSE)
  }
  if(nrow(covariance) < 1L){
    stop("'", name, "' must contain at least one level.", call. = FALSE)
  }
  check_real(
    as.vector(covariance),
    name,
    check_length = 0,
    allow_NA = FALSE
  )
  if(any(!is.finite(covariance))){
    stop("'", name, "' must contain only finite values.", call. = FALSE)
  }

  dimension_names <- dimnames(covariance)
  if(is.null(dimension_names) ||
     is.null(dimension_names[[1L]]) ||
     is.null(dimension_names[[2L]])){
    stop("'", name, "' must have row and column names.", call. = FALSE)
  }
  row_names <- as.character(dimension_names[[1L]])
  column_names <- as.character(dimension_names[[2L]])
  if(anyNA(row_names) || anyNA(column_names) ||
     any(!nzchar(row_names)) || any(!nzchar(column_names))){
    stop("'", name, "' row and column names must be non-missing and non-empty.", call. = FALSE)
  }
  if(anyDuplicated(row_names)){
    stop("'", name, "' row names must be unique.", call. = FALSE)
  }
  if(anyDuplicated(column_names)){
    stop("'", name, "' column names must be unique.", call. = FALSE)
  }
  if(!setequal(row_names, column_names)){
    stop("'", name, "' row and column names must identify the same levels.", call. = FALSE)
  }

  storage.mode(covariance) <- "double"
  dimnames(covariance) <- list(row_names, column_names)
  lower_values <- covariance[lower.tri(covariance)]
  upper_values <- covariance[upper.tri(covariance)]
  lower_supplied <- any(lower_values != 0)
  upper_supplied <- any(upper_values != 0)
  triangle_source <- "complete"
  if(lower_supplied && !upper_supplied){
    covariance[upper.tri(covariance)] <- t(covariance)[upper.tri(covariance)]
    triangle_source <- "lower"
  }else if(upper_supplied && !lower_supplied){
    covariance[lower.tri(covariance)] <- t(covariance)[lower.tri(covariance)]
    triangle_source <- "upper"
  }else if(!isTRUE(all(covariance == t(covariance)))){
    stop(
      "'", name,
      "' must be exactly symmetric, or supply values in only one triangle.",
      call. = FALSE
    )
  }
  attr(covariance, "triangle_source") <- triangle_source
  covariance
}

.bt_prepare_group_covariance_kernel <- function(x, group_levels, block_name){

  x <- .bt_as_random_group_covariance(x)
  raw <- .bt_random_group_covariance_matrix(
    x$covariance,
    name = "group_covariance$covariance"
  )
  group_levels <- as.character(group_levels)
  if(length(group_levels) < 1L || anyNA(group_levels) ||
     any(!nzchar(group_levels))){
    stop(
      "Known group covariance for random-effect block '",
      block_name,
      "' requires non-empty fitted group levels.",
      call. = FALSE
    )
  }
  if(anyDuplicated(group_levels)){
    stop(
      "Known group covariance for random-effect block '",
      block_name,
      "' requires unique fitted group levels.",
      call. = FALSE
    )
  }

  raw_levels <- rownames(raw)
  missing_levels <- group_levels[!group_levels %in% raw_levels]
  if(length(missing_levels) > 0L){
    stop(
      "Known group covariance for random-effect block '",
      block_name,
      "' is missing fitted level(s): ",
      paste(missing_levels, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  dropped_levels <- raw_levels[!raw_levels %in% group_levels]
  kernel <- raw[group_levels, group_levels, drop = FALSE]
  kernel <- .bt_random_group_covariance_scale(
    kernel = kernel,
    scale = x$scale,
    block_name = block_name
  )
  kernel <- .bt_random_group_covariance_validate_kernel(
    kernel = kernel,
    block_name = block_name
  )
  chol_kernel <- chol(kernel)
  precision <- chol2inv(chol_kernel)
  log_det <- 2 * sum(log(diag(chol_kernel)))
  if(any(!is.finite(precision)) || !is.finite(log_det)){
    stop(
      "Known group covariance for random-effect block '",
      block_name,
      "' must produce finite precision and log-determinant metadata after scaling.",
      call. = FALSE
    )
  }

  out <- list(
    type = "known",
    scale = x$scale,
    levels = group_levels,
    kernel = kernel,
    precision = precision,
    log_det = log_det,
    dropped_levels = dropped_levels
  )
  class(out) <- c("random_group_covariance_kernel", "list")
  out
}

.bt_random_group_covariance_scale <- function(kernel, scale, block_name){

  if(scale %in% c("cor", "cor0")){
    if(any(diag(kernel) <= 0)){
      stop(
        "Known group covariance for random-effect block '",
        block_name,
        "' must have positive diagonal entries when scale = '",
        scale,
        "'.",
        call. = FALSE
      )
    }
    kernel <- stats::cov2cor(kernel)
  }
  if(identical(scale, "cor0")){
    denominator <- 1 - min(kernel)
    if(!is.finite(denominator) || denominator <= 0){
      stop(
        "Known group covariance for random-effect block '",
        block_name,
        "' cannot be rescaled with scale = 'cor0'.",
        call. = FALSE
      )
    }
    kernel <- (kernel - min(kernel)) / denominator
  }else if(identical(scale, "cov0")){
    kernel <- kernel - min(kernel)
  }

  if(any(!is.finite(kernel))){
    stop(
      "Known group covariance for random-effect block '",
      block_name,
      "' produced non-finite scaled values.",
      call. = FALSE
    )
  }

  kernel
}

.bt_random_group_covariance_validate_kernel <- function(kernel, block_name){

  if(!isTRUE(all(kernel == t(kernel)))){
    stop(
      "Known group covariance for random-effect block '",
      block_name,
      "' must be symmetric after scaling.",
      call. = FALSE
    )
  }
  chol_result <- tryCatch(
    chol(kernel),
    error = function(e) NULL
  )
  if(is.null(chol_result)){
    stop(
      "Known group covariance for random-effect block '",
      block_name,
      "' must be positive definite after scaling.",
      call. = FALSE
    )
  }

  kernel
}

.bt_random_effect_group_covariance_input <- function(random_term){

  random_term$group_covariance
}

.bt_random_effect_has_known_group_covariance <- function(random_term){

  group_covariance <- .bt_random_effect_group_covariance_input(random_term)
  is.list(group_covariance) && identical(group_covariance$type, "known")
}

.bt_random_effect_group_covariance_metadata <- function(random_term){

  if(!.bt_random_effect_has_known_group_covariance(random_term)){
    return(NULL)
  }
  group_covariance <- random_term$group_covariance
  group_covariance[intersect(
    names(group_covariance),
    c("type", "scale", "levels", "kernel", "precision", "log_det", "dropped_levels")
  )]
}

.bt_random_effect_known_group_covariance <- function(random_term, context){

  if(.bt_random_effect_has_known_group_covariance(random_term)){
    return(random_term$group_covariance)
  }

  stop(
    context,
    " metadata for random-effect block '",
    random_term$block_name,
    "' are missing known group covariance metadata.",
    call. = FALSE
  )
}

.bt_random_effect_allows_new_levels <- function(random_term){

  !.bt_random_effect_has_known_group_covariance(random_term)
}

.bt_random_effect_prepare_known_group_covariance <- function(random_term,
                                                            group_levels,
                                                            n_columns,
                                                            model_matrix,
                                                            random_structure,
                                                            compile_mode,
                                                            row_indexed_external_sd,
                                                            parameterization){

  group_covariance <- .bt_random_effect_group_covariance_input(random_term)
  if(is.null(group_covariance)){
    return(NULL)
  }

  .bt_random_effect_validate_group_covariance_supported(
    random_term = random_term,
    n_columns = n_columns,
    model_matrix = model_matrix,
    random_structure = random_structure,
    compile_mode = compile_mode,
    row_indexed_external_sd = row_indexed_external_sd,
    parameterization = parameterization
  )

  .bt_prepare_group_covariance_kernel(
    x = group_covariance,
    group_levels = group_levels,
    block_name = random_term$block_name
  )
}

.bt_random_effect_validate_group_covariance_supported <- function(
    random_term,
    n_columns,
    model_matrix,
    random_structure,
    compile_mode,
    row_indexed_external_sd,
    parameterization){

  if(!compile_mode %in% c("sampled", "marginalized")){
    stop(
      "Known group covariance for random-effect block '",
      random_term$block_name,
      "' received an unknown random-effect compile mode.",
      call. = FALSE
    )
  }
  if(!random_structure %in% c("id", "diag", "us")){
    stop(
      "Known group covariance for random-effect block '",
      random_term$block_name,
      "' only supports unstructured, diagonal, or identity random-intercept blocks.",
      call. = FALSE
    )
  }
  if(n_columns > 1L && identical(parameterization, "centered")){
    stop(
      "Known group covariance for random-effect block '",
      random_term$block_name,
      "' with multiple random-effect columns requires the exact noncentered parameterization.",
      call. = FALSE
    )
  }
  if(!is.matrix(model_matrix) || ncol(model_matrix) != n_columns ||
     any(!is.finite(model_matrix))){
    stop(
      "Known group covariance for random-effect block '",
      random_term$block_name,
      "' requires a finite random-effect model matrix.",
      call. = FALSE
    )
  }
  if(isTRUE(row_indexed_external_sd)){
    stop(
      "Known group covariance for random-effect block '",
      random_term$block_name,
      "' does not support row-indexed external SD sources.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

.bt_random_effect_latent_log_density <- function(random_term, samples){

  n_groups <- random_term$n_groups
  n_columns <- random_term$n_columns
  z_names <- as.vector(.bt_random_effect_latent_names(
    random_term = random_term,
    n_groups = n_groups,
    n_columns = n_columns
  ))
  if(!all(z_names %in% names(samples))){
    stop(
      "Bridge samples are missing standardized latent random effects for block '",
      random_term$block_name,
      "'.",
      call. = FALSE
    )
  }

  z_values <- samples[z_names]
  if(any(is.na(z_values))){
    return(-Inf)
  }

  if(.bt_random_effect_has_known_group_covariance(random_term)){
    group_covariance <- .bt_random_effect_known_group_covariance(
      random_term,
      context = "Bridge sampling"
    )
    z_values <- matrix(
      as.numeric(z_values),
      nrow = n_groups,
      ncol = n_columns
    )
    out <- 0
    for(column in seq_len(n_columns)){
      out <- out + .bt_mvn_zero_log_density(
        z = z_values[, column],
        precision = group_covariance$precision,
        log_det = group_covariance$log_det
      )
    }
    return(out)
  }

  marglik <- sum(stats::dnorm(z_values, mean = 0, sd = 1, log = TRUE))
  if(is.na(marglik)){
    return(-Inf)
  }
  marglik
}

.bt_mvn_zero_log_density <- function(z, precision, log_det){

  if(any(!is.finite(z))){
    return(-Inf)
  }
  n <- length(z)
  if(!is.matrix(precision) || nrow(precision) != n || ncol(precision) != n){
    stop("Known group covariance precision matrix has incompatible dimensions.", call. = FALSE)
  }
  quadratic <- as.numeric(crossprod(z, precision %*% z))
  if(!is.finite(quadratic)){
    return(-Inf)
  }

  -0.5 * (n * log(2 * pi) + log_det + quadratic)
}

.bt_random_effect_sd_is_multiplier <- function(random_term){

  if(!.bt_random_effect_has_known_group_covariance(random_term)){
    return(FALSE)
  }
  group_covariance <- random_term$group_covariance
  identical(group_covariance$scale, "none") ||
    any(abs(diag(group_covariance$kernel) - 1) > sqrt(.Machine$double.eps))
}

.bt_random_effect_sd_summary_label <- function(component, group, random_term){

  prefix <- if(.bt_random_effect_sd_is_multiplier(random_term)){
    "sd_multiplier"
  }else{
    "sd"
  }
  paste0(prefix, "(", component, " | ", group, ")")
}
