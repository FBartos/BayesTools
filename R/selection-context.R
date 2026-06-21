#' @title Selection context helpers
#'
#' @description Validate and subset row-wise selection-kernel contexts and
#' prepare small argument lists for native selected-normal backends.
#'
#' @details Context validation checks row-wise posterior fields, observation
#' fields, and compiled selection-bin metadata together. Custom row-wise
#' fields may be listed in a character vector named \code{row_fields}; row
#' subsetting rejects undeclared fields that otherwise look row-wise. Native
#' argument helpers accept either an augmented context or a bare
#' \code{selection_backend_spec()} object, using compiled \code{step},
#' \code{phacking}, and \code{data} fallbacks where available.
#'
#' @param context selection context list.
#' @param selection_spec selection backend specification or context list.
#' @param n_samples optional expected number of posterior/sample rows.
#' @param required character vector of fields that must be present.
#' @param rows row indices or logical row mask.
#' @param idx observation indices.
#' @param S number of rows for row-wise native arguments.
#' @param alpha optional p-hacking severity values.
#' @param phack_kind optional p-hacking form indicators.
#' @param kernel_mode optional selection-kernel mode indicators.
#' @param x scalar or row-wise vector.
#' @param n required row count.
#' @param name argument name used in errors.
#'
#' @return Validated or subset context lists, or lists of prepared native
#' arguments.
#'
#' @export
selection_context_validate <- function(context, n_samples = NULL,
                                       required = character()){

  check_list(context, "context")
  if(length(required) > 0L){
    check_char(required, "required", check_length = 0, allow_NA = FALSE)
  }
  required_unknown <- setdiff(required, .selection_context_known_fields())
  if(length(required_unknown) > 0L){
    stop(
      "Unknown selection context required field",
      if(length(required_unknown) > 1L) "s" else "",
      ": '", paste(required_unknown, collapse = "', '"), "'.",
      call. = FALSE
    )
  }

  if(is.null(n_samples)){
    n_samples <- .selection_context_n_samples(context)
  }
  check_int(n_samples, "n_samples", lower = 1, allow_NA = FALSE)

  out <- context
  n_bins <- .selection_context_n_bins(out)
  out <- .selection_context_validate_row_fields(out)

  if("omega" %in% names(out) || "omega" %in% required){
    if(is.null(out[["omega"]]) || !is.matrix(out[["omega"]]) ||
       nrow(out[["omega"]]) != n_samples ||
       (!is.numeric(out[["omega"]]) && !is.integer(out[["omega"]])) ||
       any(!is.finite(out[["omega"]])) ||
       any(out[["omega"]] < 0)){
      stop("Invalid selection context 'omega'.", call. = FALSE)
    }
    if(!is.null(n_bins) && ncol(out[["omega"]]) != n_bins){
      stop(
        "Invalid selection context 'omega': the number of columns must match the selection bin count.",
        call. = FALSE
      )
    }
  }

  if("alpha" %in% names(out) || "alpha" %in% required){
    if(is.null(out[["alpha"]])){
      stop("Missing selection context 'alpha'.", call. = FALSE)
    }
    out[["alpha"]] <- selection_row_arg(out[["alpha"]], n_samples, "alpha")
    if(!is.numeric(out[["alpha"]]) && !is.integer(out[["alpha"]])){
      stop("Invalid selection context 'alpha'.", call. = FALSE)
    }
    if(any(!is.finite(out[["alpha"]]))){
      stop("Invalid selection context 'alpha'.", call. = FALSE)
    }
    out[["alpha"]] <- as.numeric(out[["alpha"]])
    if(any(out[["alpha"]] < 0 | out[["alpha"]] >= 1)){
      stop("Invalid selection context 'alpha'.", call. = FALSE)
    }
  }

  for(field in c("phack_kind", "kernel_mode", "bias_indicator")){
    if(!field %in% names(out) && !field %in% required){
      next
    }
    if(is.null(out[[field]])){
      stop("Missing selection context '", field, "'.", call. = FALSE)
    }
    out[[field]] <- selection_row_arg(out[[field]], n_samples, field)
    if(!is.numeric(out[[field]]) && !is.integer(out[[field]])){
      stop("Invalid selection context '", field, "'.", call. = FALSE)
    }
    if(any(!is.finite(out[[field]])) ||
       any(abs(out[[field]] - round(out[[field]])) > sqrt(.Machine$double.eps))){
      stop("Invalid selection context '", field, "'.", call. = FALSE)
    }
    out[[field]] <- as.integer(round(out[[field]]))
  }

  if("use_normal" %in% names(out) || "use_normal" %in% required){
    if(is.null(out[["use_normal"]]) || !is.logical(out[["use_normal"]])){
      stop("Invalid selection context 'use_normal'.", call. = FALSE)
    }
    out[["use_normal"]] <- selection_row_arg(
      out[["use_normal"]],
      n_samples,
      "use_normal"
    )
    if(any(is.na(out[["use_normal"]]))){
      stop("Invalid selection context 'use_normal'.", call. = FALSE)
    }
  }

  if("kernel_mode" %in% names(out) &&
     any(!out[["kernel_mode"]] %in% 0:3)){
    stop("Invalid selection context 'kernel_mode'.", call. = FALSE)
  }
  if("phack_kind" %in% names(out) &&
     any(!out[["phack_kind"]] %in% 0:2)){
    stop("Invalid selection context 'phack_kind'.", call. = FALSE)
  }
  if("bias_indicator" %in% names(out) &&
     any(out[["bias_indicator"]] < 1L)){
    stop("Invalid selection context 'bias_indicator'.", call. = FALSE)
  }
  n_branches <- .selection_context_n_branches(out)
  if("bias_indicator" %in% names(out) &&
     !is.null(n_branches) &&
     any(out[["bias_indicator"]] > n_branches)){
    stop("Invalid selection context 'bias_indicator'.", call. = FALSE)
  }
  out <- .selection_context_validate_observations(out, required, n_bins)

  return(out)
}


#' @rdname selection_context_validate
#' @export
selection_context_subset_rows <- function(context, rows){

  check_list(context, "context")
  S <- .selection_context_n_samples(context)
  if(is.logical(rows)){
    check_bool(rows, "rows", check_length = 0, allow_NA = FALSE)
    if(length(rows) != S){
      stop("'rows' logical mask must have one value per selection context row.",
           call. = FALSE)
    }
    rows <- which(rows)
  }else{
    if(length(rows) == 0L){
      stop("'rows' must select at least one selection context row.",
           call. = FALSE)
    }
    check_int(rows, "rows", check_length = 0, lower = 1, allow_NA = FALSE)
  }
  if(length(rows) == 0L){
    stop("'rows' must select at least one selection context row.",
         call. = FALSE)
  }

  out <- context
  if(any(rows > S)){
    stop("'rows' contains indices outside the selection context.",
         call. = FALSE)
  }
  .selection_context_reject_undeclared_rowwise_fields(out, S)

  for(field in .selection_context_row_fields(out)){
    value <- out[[field]]
    if(is.null(value)){
      next
    }
    if(is.matrix(value) && nrow(value) == S){
      out[[field]] <- value[rows, , drop = FALSE]
    }else if(!is.matrix(value) && length(value) == S){
      out[[field]] <- value[rows]
    }
  }

  out <- selection_context_validate(out, n_samples = length(rows))

  return(.selection_context_reset_native_cache(out))
}


#' @rdname selection_context_validate
#' @export
selection_context_subset_observations <- function(context, idx){

  check_list(context, "context")
  if(length(idx) == 0L){
    stop("'idx' must select at least one selection context observation.",
         call. = FALSE)
  }
  check_int(idx, "idx", check_length = 0, lower = 1, allow_NA = FALSE)

  out <- context
  for(field in c("obs_bin", "yi", "sei")){
    value <- out[[field]]
    if(is.null(value)){
      next
    }
    if(max(idx) > length(value)){
      stop("'idx' contains indices outside selection context '", field, "'.",
           call. = FALSE)
    }
    out[[field]] <- value[idx]
  }
  out <- .selection_context_validate_observations(out, character(), .selection_context_n_bins(out))

  return(.selection_context_reset_native_cache(out))
}


#' @rdname selection_context_validate
#' @export
selection_native_static_args <- function(selection_spec){

  check_list(selection_spec, "selection_spec")

  cache <- selection_spec[["native_cache"]]
  if(is.environment(cache) &&
     exists("static", envir = cache, inherits = FALSE)){
    return(get("static", envir = cache, inherits = FALSE))
  }

  segments <- selection_spec[["segments"]]
  if(is.null(segments)){
    segments <- .selection_native_segments(selection_spec)
  }

  out <- list(
    z_lower       = as.numeric(.selection_spec_z_lower(selection_spec)),
    z_upper       = as.numeric(.selection_spec_z_upper(selection_spec)),
    sign          = as.integer(.selection_spec_sign(selection_spec)),
    phack_q       = as.integer(.selection_spec_phack_q(selection_spec)),
    phack_z_source = as.numeric(.selection_spec_phack_z_source(selection_spec)),
    phack_z_dest  = as.numeric(.selection_spec_phack_z_dest(selection_spec)),
    segment_bounds = as.numeric(.selection_null_default(segments[["bounds"]], numeric())),
    segment_step_bin = as.integer(.selection_null_default(segments[["step_bin"]], integer())),
    segment_phack_region = as.integer(.selection_null_default(segments[["phack_region"]], integer())),
    telescope_probabilities = isTRUE(selection_spec[["telescope_probabilities"]])
  )
  .selection_validate_native_static_args(
    out,
    kernel_mode = .selection_spec_kernel_mode(selection_spec)
  )

  if(is.environment(cache)){
    assign("static", out, envir = cache)
  }

  return(out)
}


#' @rdname selection_context_validate
#' @export
selection_native_kernel_args <- function(selection_spec, S, alpha = NULL,
                                         phack_kind = NULL,
                                         kernel_mode = NULL){

  check_list(selection_spec, "selection_spec")
  check_int(S, "S", lower = 1, allow_NA = FALSE)

  if(is.null(alpha)){
    alpha <- .selection_null_default(selection_spec[["alpha"]], rep(0, S))
  }
  if(is.null(phack_kind)){
    if(.selection_spec_mixed_phack_q(selection_spec)){
      stop(
        "'phack_kind' is required for mixed linear/quadratic p-hacking forms.",
        call. = FALSE
      )
    }
    phack_kind <- rep(
      if(.selection_spec_has_phack(selection_spec)) .selection_spec_phack_q(selection_spec) else 0L,
      S
    )
  }
  if(is.null(kernel_mode)){
    kernel_mode <- rep(
      .selection_spec_kernel_mode(selection_spec),
      S
    )
  }

  alpha <- selection_row_arg(alpha, S, "alpha")
  if(!is.numeric(alpha) && !is.integer(alpha)){
    stop("Invalid selection native argument 'alpha'.", call. = FALSE)
  }
  if(any(!is.finite(alpha)) || any(alpha < 0 | alpha >= 1)){
    stop("Invalid selection native argument 'alpha'.", call. = FALSE)
  }
  alpha <- as.numeric(alpha)

  phack_kind <- selection_row_arg(phack_kind, S, "phack_kind")
  if(!is.numeric(phack_kind) && !is.integer(phack_kind)){
    stop("Invalid selection native argument 'phack_kind'.", call. = FALSE)
  }
  if(any(!is.finite(phack_kind)) ||
     any(abs(phack_kind - round(phack_kind)) > sqrt(.Machine$double.eps))){
    stop("Invalid selection native argument 'phack_kind'.", call. = FALSE)
  }
  phack_kind <- as.integer(round(phack_kind))
  if(any(!phack_kind %in% 0:2)){
    stop("Invalid selection native argument 'phack_kind'.", call. = FALSE)
  }

  kernel_mode <- selection_row_arg(kernel_mode, S, "kernel_mode")
  if(!is.numeric(kernel_mode) && !is.integer(kernel_mode)){
    stop("Invalid selection native argument 'kernel_mode'.", call. = FALSE)
  }
  if(any(!is.finite(kernel_mode)) ||
     any(abs(kernel_mode - round(kernel_mode)) > sqrt(.Machine$double.eps))){
    stop("Invalid selection native argument 'kernel_mode'.", call. = FALSE)
  }
  kernel_mode <- as.integer(round(kernel_mode))
  if(any(!kernel_mode %in% 0:3)){
    stop("Invalid selection native argument 'kernel_mode'.", call. = FALSE)
  }

  return(list(
    alpha       = alpha,
    phack_kind  = phack_kind,
    kernel_mode = kernel_mode,
    static      = selection_native_static_args(selection_spec)
  ))
}


#' @rdname selection_context_validate
#' @export
selection_row_arg <- function(x, n, name){

  check_int(n, "n", lower = 1, allow_NA = FALSE)
  check_char(name, "name", check_length = 1, allow_NA = FALSE)

  if(length(x) == 1L){
    return(rep(x, n))
  }
  if(length(x) != n){
    stop("Selection argument '", name, "' must have length 1 or ", n, ".",
         call. = FALSE)
  }

  return(x)
}


.selection_context_known_fields <- function(){

  c(.selection_context_builtin_row_fields(),
    .selection_context_observation_fields())
}

.selection_context_builtin_row_fields <- function(){

  c("omega", "alpha", "phack_kind", "kernel_mode",
    "bias_indicator", "use_normal")
}

.selection_context_observation_fields <- function(){

  c("obs_bin", "yi", "sei")
}

.selection_context_row_fields <- function(context){

  unique(c(.selection_context_builtin_row_fields(),
           .selection_context_declared_row_fields(context)))
}

.selection_context_declared_row_fields <- function(context){

  row_fields <- context[["row_fields"]]
  if(is.null(row_fields)){
    return(character())
  }
  if(!is.character(row_fields) ||
     anyNA(row_fields) ||
     any(!nzchar(row_fields))){
    stop("Invalid selection context 'row_fields'.", call. = FALSE)
  }
  if(any(!row_fields %in% names(context))){
    stop(
      "Selection context 'row_fields' contains fields that are not present in the context.",
      call. = FALSE
    )
  }

  return(row_fields)
}

.selection_context_validate_row_fields <- function(context){

  .selection_context_declared_row_fields(context)
  return(context)
}

.selection_context_structural_fields <- function(){

  c(
    "mode", "family", "branch_type", "prior_weights",
    "jags_omega", "jags_alpha", "jags_pi_null", "jags_beta_null",
    "jags_phack_kind", "jags_phack_z_source", "jags_phack_z_dest",
    "jags_kernel_mode", "jags_kernel_mode_expr", "jags_code",
    "step", "phacking", "prior_code", "transform_code", "monitor", "init",
    "data", "backend_data", "jags_data", "native_cache", "row_fields",
    "z_lower", "z_upper", "sign", "n_bins", "p_rule", "p_cuts",
    "telescope_probabilities", "has_step", "has_phack", "phack_q",
    "phack_q_values", "mixed_phack_q", "phack_z_source", "phack_z_dest",
    "segments", "branch_kernel_mode", "fixed_omega",
    "jags_use_step_switch"
  )
}

.selection_context_reject_undeclared_rowwise_fields <- function(context, S){

  allowed <- unique(c(.selection_context_row_fields(context),
                      .selection_context_observation_fields(),
                      .selection_context_structural_fields()))
  candidates <- setdiff(names(context), allowed)
  candidates <- candidates[nzchar(candidates)]
  if(length(candidates) == 0L || S <= 1L){
    return(invisible(TRUE))
  }

  rowwise <- vapply(candidates, function(field){
    value <- context[[field]]
    if(is.matrix(value)){
      return(nrow(value) == S)
    }
    is.atomic(value) && length(value) == S
  }, logical(1))

  if(any(rowwise)){
    stop(
      "Selection context field",
      if(sum(rowwise) > 1L) "s" else "",
      " '", paste(candidates[rowwise], collapse = "', '"),
      "' appear row-wise but are not declared in 'row_fields'.",
      call. = FALSE
    )
  }

  return(invisible(TRUE))
}

.selection_context_n_bins <- function(context){

  bins <- integer()
  add_bin <- function(value){
    if(is.null(value) || length(value) != 1L || is.na(value)){
      return()
    }
    if(!is.numeric(value) && !is.integer(value)){
      return()
    }
    if(!is.finite(value) || abs(value - round(value)) > sqrt(.Machine$double.eps) || value < 1){
      stop("Invalid selection context bin count.", call. = FALSE)
    }
    bins <<- c(bins, as.integer(round(value)))
  }

  add_bin(context[["n_bins"]])
  if(!is.null(context[["step"]])){
    add_bin(context[["step"]][["n_bins"]])
  }
  if(!is.null(context[["data"]])){
    add_bin(context[["data"]][["sel_n_bins"]])
  }
  if(!is.null(context[["z_lower"]])){
    bins <- c(bins, length(context[["z_lower"]]))
  }
  if(!is.null(context[["z_upper"]])){
    bins <- c(bins, length(context[["z_upper"]]))
  }
  if(!is.null(context[["data"]][["sel_z_lower"]])){
    bins <- c(bins, length(context[["data"]][["sel_z_lower"]]))
  }
  if(!is.null(context[["data"]][["sel_z_upper"]])){
    bins <- c(bins, length(context[["data"]][["sel_z_upper"]]))
  }

  bins <- unique(bins)
  if(length(bins) == 0L){
    return(NULL)
  }
  if(length(bins) > 1L){
    stop("Selection context bin count fields are inconsistent.", call. = FALSE)
  }

  return(bins)
}

.selection_context_n_branches <- function(context){

  counts <- integer()
  if(!is.null(context[["branch_type"]])){
    counts <- c(counts, length(context[["branch_type"]]))
  }
  if(!is.null(context[["prior_weights"]])){
    counts <- c(counts, length(context[["prior_weights"]]))
  }
  counts <- unique(counts[counts > 0L])
  if(length(counts) == 0L){
    return(NULL)
  }
  if(length(counts) > 1L){
    stop("Selection context branch metadata are inconsistent.", call. = FALSE)
  }

  return(counts)
}

.selection_context_validate_observations <- function(context, required, n_bins){

  lengths <- integer()
  for(field in .selection_context_observation_fields()){
    if(!field %in% names(context) && !field %in% required){
      next
    }
    if(is.null(context[[field]])){
      stop("Missing selection context '", field, "'.", call. = FALSE)
    }
    value <- context[[field]]
    if(length(value) == 0L){
      stop("Invalid selection context '", field, "'.", call. = FALSE)
    }
    if(field == "obs_bin"){
      if((!is.numeric(value) && !is.integer(value)) ||
         any(!is.finite(value)) ||
         any(abs(value - round(value)) > sqrt(.Machine$double.eps)) ||
         any(value < 1)){
        stop("Invalid selection context 'obs_bin'.", call. = FALSE)
      }
      if(!is.null(n_bins) && any(value > n_bins)){
        stop("Invalid selection context 'obs_bin'.", call. = FALSE)
      }
      context[["obs_bin"]] <- as.integer(round(value))
    }else if(field == "yi"){
      if((!is.numeric(value) && !is.integer(value)) ||
         any(!is.finite(value))){
        stop("Invalid selection context 'yi'.", call. = FALSE)
      }
      context[["yi"]] <- as.numeric(value)
    }else if(field == "sei"){
      if((!is.numeric(value) && !is.integer(value)) ||
         any(!is.finite(value)) ||
         any(value <= 0)){
        stop("Invalid selection context 'sei'.", call. = FALSE)
      }
      context[["sei"]] <- as.numeric(value)
    }
    lengths <- c(lengths, length(context[[field]]))
  }

  lengths <- unique(lengths)
  if(length(lengths) > 1L){
    stop("Selection context observation fields must have the same length.",
         call. = FALSE)
  }

  return(context)
}

.selection_spec_data <- function(selection_spec){

  data <- selection_spec[["data"]]
  if(is.null(data)){
    data <- selection_spec[["backend_data"]]
  }
  if(is.null(data)){
    data <- list()
  }

  return(data)
}

.selection_spec_mode_code <- function(selection_spec){

  data <- .selection_spec_data(selection_spec)
  if(!is.null(selection_spec[["kernel_mode"]])){
    mode <- selection_spec[["kernel_mode"]]
    if(length(mode) == 1L){
      return(as.integer(mode))
    }
  }
  if(!is.null(data[["kernel_mode"]])){
    return(as.integer(data[["kernel_mode"]]))
  }
  if(!is.null(selection_spec[["mode"]])){
    return(.selection_mode_code(selection_spec[["mode"]]))
  }

  return(0L)
}

.selection_spec_kernel_mode <- function(selection_spec){

  mode <- .selection_spec_mode_code(selection_spec)
  if(length(mode) != 1L ||
     is.na(mode) ||
     !is.finite(mode) ||
     abs(mode - round(mode)) > sqrt(.Machine$double.eps) ||
     !as.integer(round(mode)) %in% 0:3){
    stop("Invalid selection specification 'kernel_mode'.", call. = FALSE)
  }

  return(as.integer(round(mode)))
}

.selection_spec_has_phack <- function(selection_spec){

  if(!is.null(selection_spec[["has_phack"]])){
    return(isTRUE(selection_spec[["has_phack"]]))
  }
  .selection_spec_kernel_mode(selection_spec) %in% c(2L, 3L)
}

.selection_spec_z_lower <- function(selection_spec){

  data <- .selection_spec_data(selection_spec)
  out <- selection_spec[["z_lower"]]
  if(is.null(out)){
    out <- selection_spec[["step"]][["z_lower"]]
  }
  if(is.null(out)){
    out <- data[["sel_z_lower"]]
  }
  if(is.null(out)){
    out <- numeric()
  }

  return(out)
}

.selection_spec_z_upper <- function(selection_spec){

  data <- .selection_spec_data(selection_spec)
  out <- selection_spec[["z_upper"]]
  if(is.null(out)){
    out <- selection_spec[["step"]][["z_upper"]]
  }
  if(is.null(out)){
    out <- data[["sel_z_upper"]]
  }
  if(is.null(out)){
    out <- numeric()
  }

  return(out)
}

.selection_spec_sign <- function(selection_spec){

  data <- .selection_spec_data(selection_spec)
  out <- selection_spec[["sign"]]
  if(is.null(out)){
    out <- data[["sel_sign"]]
  }
  if(is.null(out)){
    out <- 1L
  }

  return(out)
}

.selection_spec_phack_q <- function(selection_spec){

  if(!.selection_spec_has_phack(selection_spec)){
    return(1L)
  }
  q <- selection_spec[["phack_q"]]
  if(is.null(q)){
    q <- selection_spec[["phacking"]][["q"]]
  }
  if(is.null(q) || length(q) == 0L){
    return(1L)
  }
  q <- unique(q)
  if(length(q) > 1L){
    return(1L)
  }

  return(q)
}

.selection_spec_mixed_phack_q <- function(selection_spec){

  if(!is.null(selection_spec[["mixed_phack_q"]])){
    return(isTRUE(selection_spec[["mixed_phack_q"]]))
  }
  q <- selection_spec[["phacking"]][["q"]]
  !is.null(q) && length(unique(q)) > 1L
}

.selection_spec_phack_z_source <- function(selection_spec){

  out <- selection_spec[["phack_z_source"]]
  if(is.null(out)){
    out <- selection_spec[["phacking"]][["z_source"]]
  }
  if(is.null(out)){
    out <- c(0, 0)
  }
  if(is.matrix(out)){
    if(nrow(out) != 1L &&
       any(out != matrix(out[1L,], nrow = nrow(out), ncol = ncol(out), byrow = TRUE))){
      stop(
        "Selection specification requires explicit 'segments' for mixed p-hacking geometry.",
        call. = FALSE
      )
    }
    out <- as.numeric(out[1L,])
  }

  return(out)
}

.selection_spec_phack_z_dest <- function(selection_spec){

  out <- selection_spec[["phack_z_dest"]]
  if(is.null(out)){
    out <- selection_spec[["phacking"]][["z_destination"]]
  }
  if(is.null(out)){
    out <- c(0, 0)
  }
  if(is.matrix(out)){
    if(nrow(out) != 1L &&
       any(out != matrix(out[1L,], nrow = nrow(out), ncol = ncol(out), byrow = TRUE))){
      stop(
        "Selection specification requires explicit 'segments' for mixed p-hacking geometry.",
        call. = FALSE
      )
    }
    out <- as.numeric(out[1L,])
  }

  return(out)
}

.selection_spec_p_cuts <- function(selection_spec){

  data <- .selection_spec_data(selection_spec)
  out <- selection_spec[["p_cuts"]]
  if(is.null(out)){
    out <- selection_spec[["step"]][["breaks"]]
  }
  if(is.null(out)){
    out <- data[["sel_p_cuts"]]
  }
  if(is.null(out)){
    out <- c(0, 1)
  }

  return(.selection_validate_global_breaks(out))
}

.selection_native_segment_midpoint <- function(lower, upper){

  if(is.infinite(lower) && lower < 0){
    return(upper - 1)
  }
  if(is.infinite(upper) && upper > 0){
    return(lower + 1)
  }

  return((lower + upper) / 2)
}

.selection_native_step_bin_from_z <- function(z, p_cuts){

  p_value <- stats::pnorm(z, lower.tail = FALSE)
  close <- vapply(p_cuts, function(cut) abs(p_value - cut) <= 1e-12, logical(1))
  if(any(close)){
    p_value <- p_cuts[which(close)[1L]]
  }
  bin <- findInterval(p_value, p_cuts, rightmost.closed = TRUE, left.open = TRUE)
  bin <- pmin(pmax(bin, 1L), length(p_cuts) - 1L)

  return(as.integer(bin))
}

.selection_native_segments <- function(selection_spec){

  p_cuts <- .selection_spec_p_cuts(selection_spec)
  z_lower <- stats::qnorm(1 - p_cuts[-1])
  z_upper <- stats::qnorm(1 - p_cuts[-length(p_cuts)])
  bounds <- c(-Inf, Inf, z_lower[is.finite(z_lower)], z_upper[is.finite(z_upper)])

  has_phacking <- .selection_spec_has_phack(selection_spec)
  phack_z_source <- .selection_spec_phack_z_source(selection_spec)
  phack_z_dest <- .selection_spec_phack_z_dest(selection_spec)
  if(has_phacking){
    bounds <- c(bounds, phack_z_source, phack_z_dest)
  }

  bounds <- sort(unique(bounds))
  n_segments <- length(bounds) - 1L
  step_bin <- integer(n_segments)
  phack_region <- integer(n_segments)

  for(i in seq_len(n_segments)){
    mid <- .selection_native_segment_midpoint(bounds[i], bounds[i + 1L])
    step_bin[i] <- .selection_native_step_bin_from_z(mid, p_cuts)
    if(has_phacking){
      if(mid >= phack_z_source[1L] && mid <= phack_z_source[2L]){
        phack_region[i] <- 1L
      }else if(mid > phack_z_dest[1L] && mid <= phack_z_dest[2L]){
        phack_region[i] <- 2L
      }
    }
  }

  return(list(
    bounds       = bounds,
    step_bin     = step_bin,
    phack_region = phack_region
  ))
}

.selection_validate_native_static_args <- function(args, kernel_mode){

  if((!is.numeric(args[["z_lower"]]) && !is.integer(args[["z_lower"]])) ||
     (!is.numeric(args[["z_upper"]]) && !is.integer(args[["z_upper"]])) ||
     anyNA(args[["z_lower"]]) ||
     anyNA(args[["z_upper"]]) ||
     length(args[["z_lower"]]) != length(args[["z_upper"]])){
    stop("Invalid selection native static z bounds.", call. = FALSE)
  }
  if(kernel_mode != 0L && length(args[["z_lower"]]) == 0L){
    stop("Invalid selection native static z bounds.", call. = FALSE)
  }
  if(any(args[["z_lower"]] >= args[["z_upper"]])){
    stop("Invalid selection native static z bounds.", call. = FALSE)
  }

  if(length(args[["sign"]]) != 1L ||
     is.na(args[["sign"]]) ||
     !args[["sign"]] %in% c(-1L, 1L)){
    stop("Invalid selection native static sign.", call. = FALSE)
  }
  if(length(args[["phack_q"]]) != 1L ||
     is.na(args[["phack_q"]]) ||
     !args[["phack_q"]] %in% c(1L, 2L)){
    stop("Invalid selection native static phack_q.", call. = FALSE)
  }
  for(field in c("phack_z_source", "phack_z_dest")){
    if((!is.numeric(args[[field]]) && !is.integer(args[[field]])) ||
       length(args[[field]]) != 2L ||
       any(!is.finite(args[[field]])) ||
       args[[field]][1L] > args[[field]][2L]){
      stop("Invalid selection native static ", field, ".", call. = FALSE)
    }
  }

  bounds <- args[["segment_bounds"]]
  step_bin <- args[["segment_step_bin"]]
  phack_region <- args[["segment_phack_region"]]
  if((!is.numeric(bounds) && !is.integer(bounds)) ||
     anyNA(bounds) ||
     length(bounds) == 1L ||
     any(diff(bounds) <= 0)){
    stop("Invalid selection native static segment bounds.", call. = FALSE)
  }
  if(length(step_bin) != max(length(bounds) - 1L, 0L) ||
     length(phack_region) != max(length(bounds) - 1L, 0L)){
    stop("Invalid selection native static segments.", call. = FALSE)
  }
  n_bins <- length(args[["z_lower"]])
  if(length(step_bin) > 0L &&
     (any(is.na(step_bin)) || any(step_bin < 1L) || any(step_bin > n_bins))){
    stop("Invalid selection native static segment step bins.", call. = FALSE)
  }
  if(length(phack_region) > 0L &&
     (any(is.na(phack_region)) || any(!phack_region %in% 0:2))){
    stop("Invalid selection native static segment p-hacking regions.",
         call. = FALSE)
  }

  return(invisible(TRUE))
}

.selection_context_n_samples <- function(context){

  counts <- integer()
  if(!is.null(context[["omega"]]) && !is.null(nrow(context[["omega"]]))){
    counts <- c(counts, nrow(context[["omega"]]))
  }
  for(field in .selection_context_row_fields(context)){
    if(field == "omega"){
      next
    }
    if(!is.null(context[[field]])){
      counts <- c(counts, length(context[[field]]))
    }
  }
  counts <- unique(counts[counts > 0L])

  non_scalar <- counts[counts > 1L]
  if(length(non_scalar) > 1L){
    stop("Selection context row fields have inconsistent lengths.",
         call. = FALSE)
  }
  if(length(non_scalar) == 1L){
    return(non_scalar)
  }
  if(1L %in% counts){
    return(1L)
  }

  stop("Cannot infer the number of rows in selection context.",
       call. = FALSE)
}


.selection_context_reset_native_cache <- function(context){

  if(!is.null(context)){
    context[["native_cache"]] <- new.env(parent = emptyenv())
  }

  return(context)
}


