#' @title Selection context helpers
#'
#' @description Validate and subset row-wise selection-kernel contexts and
#' prepare small argument lists for native selected-normal backends.
#'
#' @details Context validation checks row-wise posterior fields, observation
#' fields, and compiled selection-bin metadata together. Custom row-wise
#' fields may be listed in a character vector named \code{row_fields}; row
#' subsetting rejects undeclared fields that otherwise look row-wise.
#' The built-in \code{vector_rule} field accepts the integer codes \code{0}
#' (product), \code{1} (best one-sided p-value), and \code{2} (best two-sided
#' p-value), with one value or one value per posterior row.
#' Native argument helpers accept either an augmented context or a bare
#' \code{selection_backend_spec()} object, using compiled \code{step},
#' \code{phacking}, and \code{data} fallbacks where available.
#' Spec-level \code{kernel_mode} is the union of branch kernels (a capability
#' flag). Row routing uses \code{branch_kernel_mode} or an explicit per-row
#' vector. Branches with kernel mode \code{0} carry no selection kernel and
#' never compete for a route, so a specification whose active branches all
#' share one kernel routes on that kernel. When two or more \emph{different}
#' active kernels coexist, their bitwise union names a third kernel rather
#' than a route, and passing it is an error unless \code{bias_indicator} can
#' map each row to a branch.
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
  out <- .selection_context_validate_row_fields(out, n_samples)

  if("omega" %in% names(out) || "omega" %in% required){
    # One pass over the weight matrix rejects what the separate finiteness and
    # sign scans rejected: a missing, infinite or NaN entry leaves the range
    # non-finite, and a negative entry lowers its minimum.
    omega_range <- if((is.numeric(out[["omega"]]) || is.integer(out[["omega"]])) &&
                      length(out[["omega"]]) > 0L){
      range(out[["omega"]])
    }else{
      # Nothing to scan; the type and shape clauses below decide the outcome.
      c(0, 0)
    }
    if(is.null(out[["omega"]]) || !is.matrix(out[["omega"]]) ||
       nrow(out[["omega"]]) != n_samples ||
       (!is.numeric(out[["omega"]]) && !is.integer(out[["omega"]])) ||
       !is.finite(omega_range[1L]) || !is.finite(omega_range[2L]) ||
       omega_range[1L] < 0){
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
    # One pass rejects the same values: a missing, infinite or NaN severity
    # leaves the range non-finite, and the bounds are the range's bounds.
    alpha_range <- range(out[["alpha"]])
    if(!is.finite(alpha_range[1L]) || !is.finite(alpha_range[2L]) ||
       alpha_range[1L] < 0 || alpha_range[2L] >= 1){
      stop("Invalid selection context 'alpha'.", call. = FALSE)
    }
    if(!is.double(out[["alpha"]])){
      out[["alpha"]] <- as.numeric(out[["alpha"]])
    }
  }

  for(field in c("phack_kind", "kernel_mode", "bias_indicator")){
    if(!field %in% names(out) && !field %in% required){
      next
    }
    if(is.null(out[[field]])){
      stop("Missing selection context '", field, "'.", call. = FALSE)
    }
    out[[field]] <- .selection_context_integer_field(
      selection_row_arg(out[[field]], n_samples, field),
      field
    )
  }

  if("vector_rule" %in% names(out) || "vector_rule" %in% required){
    if(is.null(out[["vector_rule"]])){
      stop("Missing selection context 'vector_rule'.", call. = FALSE)
    }
    vector_rule <- out[["vector_rule"]]
    if(!is.numeric(vector_rule) || !is.null(dim(vector_rule)) ||
       any(!vector_rule %in% 0:2)){
      stop("Invalid selection context 'vector_rule'.", call. = FALSE)
    }
    out[["vector_rule"]] <- as.integer(selection_row_arg(vector_rule, n_samples, "vector_rule"))
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
    if(anyNA(out[["use_normal"]])){
      stop("Invalid selection context 'use_normal'.", call. = FALSE)
    }
  }

  # The routing fields are validated integer vectors without missing values by
  # now, so their range decides membership exactly as the element-wise tests
  # did, in one pass instead of one logical vector per test.
  if("kernel_mode" %in% names(out)){
    kernel_mode_range <- range(out[["kernel_mode"]])
    if(kernel_mode_range[1L] < 0L || kernel_mode_range[2L] > 3L){
      stop("Invalid selection context 'kernel_mode'.", call. = FALSE)
    }
  }
  if("phack_kind" %in% names(out)){
    phack_kind_range <- range(out[["phack_kind"]])
    if(phack_kind_range[1L] < 0L || phack_kind_range[2L] > 2L){
      stop("Invalid selection context 'phack_kind'.", call. = FALSE)
    }
  }
  bias_indicator_range <- if("bias_indicator" %in% names(out)){
    range(out[["bias_indicator"]])
  }else{
    NULL
  }
  if(!is.null(bias_indicator_range) && bias_indicator_range[1L] < 1L){
    stop("Invalid selection context 'bias_indicator'.", call. = FALSE)
  }
  n_branches <- .selection_context_n_branches(out)
  if(!is.null(bias_indicator_range) &&
     !is.null(n_branches) &&
     bias_indicator_range[2L] > n_branches){
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

  # Bridge sampling calls this once per evaluated state with S = 1, so the
  # row-wise arguments of a call are the constants the previous call already
  # validated and expanded. The validated set is kept with the spec's native
  # cache -- the same per-object cache that already holds the static arguments
  # -- and is reused only while every input this function reads is unchanged,
  # 'S' included, so a reused set is one whose inputs were accepted here; any
  # other input revalidates and rejects exactly as before.
  cache <- selection_spec[["native_cache"]]
  cache_key <- if(is.environment(cache)){
    list(
      S = S,
      alpha = alpha,
      phack_kind = phack_kind,
      kernel_mode = kernel_mode,
      spec = selection_spec[.selection_native_kernel_args_fields()]
    )
  }else{
    NULL
  }
  if(!is.null(cache_key) &&
     exists("kernel_args_key", envir = cache, inherits = FALSE) &&
     identical(
       cache_key,
       get("kernel_args_key", envir = cache, inherits = FALSE)
     )){
    return(get("kernel_args", envir = cache, inherits = FALSE))
  }
  check_int(S, "S", lower = 1, allow_NA = FALSE)

  # A row-wise argument is usually one value shared by every posterior row.
  # Validating the supplied value and expanding it afterwards keeps every
  # rejection while the checks stop repeating over rows known to be identical.
  if(is.null(alpha)){
    alpha <- .selection_null_default(selection_spec[["alpha"]], 0)
  }
  if(is.null(phack_kind)){
    if(.selection_spec_mixed_phack_q(selection_spec)){
      stop(
        "'phack_kind' is required for mixed linear/quadratic p-hacking forms.",
        call. = FALSE
      )
    }
    phack_kind <-
      if(.selection_spec_has_phack(selection_spec)) .selection_spec_phack_q(selection_spec) else 0L
  }

  .selection_native_row_length(alpha, S, "alpha")
  if(!is.numeric(alpha) && !is.integer(alpha)){
    stop("Invalid selection native argument 'alpha'.", call. = FALSE)
  }
  alpha_range <- range(alpha)
  if(!is.finite(alpha_range[1L]) || !is.finite(alpha_range[2L]) ||
     alpha_range[1L] < 0 || alpha_range[2L] >= 1){
    stop("Invalid selection native argument 'alpha'.", call. = FALSE)
  }
  alpha <- selection_row_arg(as.numeric(alpha), S, "alpha")

  phack_kind <- .selection_native_integer_arg(
    phack_kind, S, "phack_kind", 2L
  )
  kernel_mode <- .selection_native_integer_arg(
    .selection_row_kernel_mode(selection_spec, kernel_mode, S),
    S, "kernel_mode", 3L
  )

  out <- list(
    alpha       = alpha,
    phack_kind  = phack_kind,
    kernel_mode = kernel_mode,
    static      = selection_native_static_args(selection_spec)
  )

  if(!is.null(cache_key)){
    assign("kernel_args_key", cache_key, envir = cache)
    assign("kernel_args", out, envir = cache)
  }

  return(out)
}


# Every specification field selection_native_kernel_args() reads, directly or
# through the defaults and the row routing it resolves. The cached arguments
# are reused only while all of them are unchanged.
.selection_native_kernel_args_fields <- function(){

  c("alpha", "phack_kind", "phack_q", "has_phack", "mixed_phack_q",
    "phacking", "kernel_mode", "branch_kernel_mode", "bias_indicator")
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


# The length rejection of selection_row_arg(), without the expansion: a native
# row argument is validated before it is repeated over the posterior rows.
.selection_native_row_length <- function(x, n, name){

  if(length(x) != 1L && length(x) != n){
    stop("Selection argument '", name, "' must have length 1 or ", n, ".",
         call. = FALSE)
  }

  return(invisible(TRUE))
}


# Validate one row-wise integer native argument before expanding it. An integer
# vector needs no rounding repair, and a real vector's finiteness scan and its
# membership test become range comparisons; the accepted values and every
# rejection are the ones the element-wise checks produced.
.selection_native_integer_arg <- function(x, n, name, upper){

  .selection_native_row_length(x, n, name)
  if(!is.numeric(x) && !is.integer(x)){
    stop("Invalid selection native argument '", name, "'.", call. = FALSE)
  }
  if(is.integer(x)){
    if(anyNA(x)){
      stop("Invalid selection native argument '", name, "'.", call. = FALSE)
    }
  }else{
    limits <- range(x)
    if(!is.finite(limits[1L]) || !is.finite(limits[2L]) ||
       any(abs(x - round(x)) > sqrt(.Machine$double.eps))){
      stop("Invalid selection native argument '", name, "'.", call. = FALSE)
    }
    x <- as.integer(round(x))
  }
  limits <- range(x)
  if(limits[1L] < 0L || limits[2L] > upper){
    stop("Invalid selection native argument '", name, "'.", call. = FALSE)
  }

  return(selection_row_arg(x, n, name))
}


# Row-wise integer fields are already integer vectors on every path that builds
# a context from posterior draws. Recognizing that skips the rounding repair
# such a vector cannot need, and the finiteness scan of a real vector becomes
# one range pass; the values and the rejections are the ones the element-wise
# checks produce.
.selection_context_integer_field <- function(x, field){

  if(!is.numeric(x) && !is.integer(x)){
    stop("Invalid selection context '", field, "'.", call. = FALSE)
  }
  if(is.integer(x)){
    if(anyNA(x)){
      stop("Invalid selection context '", field, "'.", call. = FALSE)
    }

    return(x)
  }

  limits <- range(x)
  if(!is.finite(limits[1L]) || !is.finite(limits[2L]) ||
     any(abs(x - round(x)) > sqrt(.Machine$double.eps))){
    stop("Invalid selection context '", field, "'.", call. = FALSE)
  }

  return(as.integer(round(x)))
}

.selection_context_known_fields <- function(){

  c(.selection_context_builtin_row_fields(),
    .selection_context_observation_fields())
}

.selection_context_builtin_row_fields <- function(){

  c("omega", "alpha", "phack_kind", "kernel_mode", "vector_rule",
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

.selection_context_validate_row_fields <- function(context, n_samples){

  row_fields <- .selection_context_declared_row_fields(context)
  for(field in row_fields){
    value <- context[[field]]
    field_rows <- if(is.matrix(value)) nrow(value) else length(value)
    if(!field_rows %in% c(1L, n_samples)){
      stop(
        "Selection context row field '", field,
        "' must have either one row/value or 'n_samples' rows/values.",
        call. = FALSE
      )
    }
  }

  return(context)
}

.selection_context_structural_fields <- function(){

  c(
    "mode", "family", "branch_type", "prior_weights",
    "jags_omega", "jags_alpha", "jags_pi_null", "jags_beta_null",
    "jags_phack_kind", "jags_phack_z_source", "jags_phack_z_dest",
    "jags_kernel_mode", "jags_kernel_mode_expr", "jags_vector_rule", "jags_code",
    "step", "phacking", "prior_code", "transform_code", "monitor", "init",
    "data", "backend_data", "jags_data", "native_cache", "row_fields",
    "z_lower", "z_upper", "sign", "n_bins", "p_rule", "p_cuts",
    "telescope_probabilities", "has_step", "has_phack", "phack_q",
    "phack_q_values", "mixed_phack_q", "phack_z_source", "phack_z_dest",
    "segments", "branch_kernel_mode", "branch_vector_rule", "branch_model", "fixed_omega",
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

.selection_spec_kernel_mode <- function(selection_spec){

  # An augmented context carries one kernel mode per posterior row, and the
  # native batches ask for it twice per call, so this is a scan over the whole
  # batch. A validated integer vector needs no rounding repair and no copy, and
  # membership in 0:3 is the range of an integer vector. The rejections are the
  # ones the element-wise tests made.
  mode <- selection_spec[["kernel_mode"]]
  if(length(mode) == 0L){
    stop("Invalid selection specification 'kernel_mode'.", call. = FALSE)
  }
  if(is.integer(mode)){
    if(anyNA(mode)){
      stop("Invalid selection specification 'kernel_mode'.", call. = FALSE)
    }
  }else{
    limits <- range(mode)
    if(!is.finite(limits[1L]) || !is.finite(limits[2L]) ||
       any(abs(mode - round(mode)) > sqrt(.Machine$double.eps))){
      stop("Invalid selection specification 'kernel_mode'.", call. = FALSE)
    }
    mode <- as.integer(round(mode))
  }
  limits <- range(mode)
  if(limits[1L] < 0L || limits[2L] > 3L){
    stop("Invalid selection specification 'kernel_mode'.", call. = FALSE)
  }

  return(mode)
}

# Spec-level kernel_mode is the union of branch kernels (capability flag).
# Row routing must use branch_kernel_mode (or an explicit per-row vector).
.selection_row_kernel_mode <- function(selection_spec, kernel_mode, S){

  branch_modes <- selection_spec[["branch_kernel_mode"]]
  if(!is.null(branch_modes)){
    branch_modes <- as.integer(branch_modes)
  }
  unique_branches <- unique(branch_modes)
  # Mode 0 is the absence of a selection kernel, not a competing route: a
  # branch without selection never enters the kernel, so it cannot disagree
  # with one that does. Only two or more distinct *active* kernels make a row
  # route ambiguous - and there the bitwise union names a third kernel
  # (step | phack_power == step_phack_power), which is why it must not be used
  # as a route.
  active_branches <- unique(branch_modes[branch_modes != 0L])
  ambiguous <- length(active_branches) > 1L
  route_mode <- if(length(unique_branches) == 0L){
    integer()
  }else if(length(active_branches) == 0L){
    0L
  }else if(!ambiguous){
    active_branches
  }else{
    integer()
  }
  union_mode <- if(length(unique_branches) > 0L){
    Reduce(function(a, b) bitwOr(as.integer(a), as.integer(b)), unique_branches)
  }else{
    NA_integer_
  }
  spec_union <- unique(.selection_spec_kernel_mode(selection_spec))
  if(length(spec_union) == 1L){
    union_mode <- spec_union
  }

  if(is.null(kernel_mode)){
    if(length(route_mode) == 1L){
      return(route_mode)
    }
    if(ambiguous){
      kernel_mode <- .selection_map_union_kernel_mode(
        selection_spec,
        union_mode,
        branch_modes,
        S
      )
      if(!is.null(kernel_mode)){
        return(kernel_mode)
      }
      stop(
        "Row kernel_mode is required when selection branches use different kernels. ",
        "The spec-level kernel_mode is a union capability flag, not a row route.",
        call. = FALSE
      )
    }
    return(.selection_spec_kernel_mode(selection_spec))
  }

  kernel_mode <- as.integer(round(kernel_mode))
  if(length(kernel_mode) == 1L &&
     ambiguous &&
     identical(kernel_mode, as.integer(union_mode))){
    mapped <- .selection_map_union_kernel_mode(
      selection_spec,
      union_mode,
      branch_modes,
      S
    )
    if(!is.null(mapped)){
      return(mapped)
    }
    stop(
      "Cannot route rows on the union kernel_mode. Pass per-row modes or ",
      "branch_kernel_mode indexed by bias_indicator.",
      call. = FALSE
    )
  }

  kernel_mode
}

.selection_map_union_kernel_mode <- function(selection_spec, union_mode,
                                             branch_modes, S){

  if(is.null(branch_modes) ||
     length(unique(branch_modes[branch_modes != 0L])) <= 1L){
    return(NULL)
  }
  indicator <- selection_spec[["bias_indicator"]]
  if(is.null(indicator)){
    return(NULL)
  }
  indicator <- selection_row_arg(indicator, S, "bias_indicator")
  if(any(indicator < 1L | indicator > length(branch_modes))){
    stop("Invalid selection context 'bias_indicator'.", call. = FALSE)
  }

  branch_modes[indicator]
}

.selection_spec_has_phack <- function(selection_spec){

  if(!is.null(selection_spec[["has_phack"]])){
    return(isTRUE(selection_spec[["has_phack"]]))
  }
  any(.selection_spec_kernel_mode(selection_spec) %in% c(2L, 3L))
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
    # Experimental p-hacking: missing q is not a geometry. Static native
    # arguments currently fall back to linear (bin 1).
    return(1L)
  }
  q <- unique(q)
  if(length(q) > 1L){
    # Experimental p-hacking: mixed linear/quadratic forms cannot share one
    # static phack_q. Callers that need mixed q must pass phack_kind per row;
    # this fallback is bin 1 and is not a geometry choice.
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
      # Experimental p-hacking: mixed source cuts cannot be expressed through
      # `segments`. Reject rather than inventing a third geometry path.
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
      # Experimental p-hacking: mixed destination cuts cannot be expressed
      # through `segments`. Reject rather than inventing a third geometry path.
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
  bin <- findInterval(p_value, p_cuts, rightmost.closed = TRUE, left.open = TRUE)
  bin <- pmin(pmax(bin, 1L), length(p_cuts) - 1L)

  return(as.integer(bin))
}

.selection_native_segments <- function(selection_spec){

  p_cuts <- .selection_spec_p_cuts(selection_spec)
  z_lower <- stats::qnorm(p_cuts[-1], lower.tail = FALSE)
  z_upper <- stats::qnorm(
    p_cuts[-length(p_cuts)],
    lower.tail = FALSE
  )
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
  if(any(kernel_mode != 0L) && length(args[["z_lower"]]) == 0L){
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
      count <- if(is.matrix(context[[field]])){
        nrow(context[[field]])
      }else{
        length(context[[field]])
      }
      counts <- c(counts, count)
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
