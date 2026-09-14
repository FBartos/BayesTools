# Marginal-likelihood result contract.

.bt_marglik_result_version <- 1L

.bt_marglik_empty_repetitions <- function(){

  data.frame(
    repetition = integer(),
    logml = numeric(),
    niter = numeric(),
    mcse = numeric(),
    finite = logical(),
    success = logical(),
    within_maxiter = logical(),
    warning = character(),
    error = character(),
    method = character(),
    n_chains = integer(),
    n_draws = integer(),
    stringsAsFactors = FALSE
  )
}

.bt_marglik_manual_result <- function(logml){

  if(!is.numeric(logml) || length(logml) != 1L || is.na(logml) ||
     is.nan(logml) || is.infinite(logml) && logml > 0){
    stop(
      "'logml' must be one numeric natural-log marginal likelihood and may only be infinite when it is -Inf.",
      call. = FALSE
    )
  }

  out <- list(
    schema_version = .bt_marglik_result_version,
    logml = as.numeric(logml),
    scale = "natural_log",
    repetitions = .bt_marglik_empty_repetitions(),
    aggregation = list(
      rule = "supplied_scalar",
      nonfinite_policy = "not_applicable",
      n_repetitions = 0L,
      n_included = 1L,
      n_failed = 0L
    ),
    diagnostics = list(
      upstream = NULL,
      upstream_warnings = character(),
      chains = list(
        count = NA_integer_,
        draws = NA_integer_,
        draws_per_chain = integer()
      )
    )
  )
  class(out) <- c("BayesTools_marglik", "list")
  out
}

.bt_marglik_exact_result <- function(logml, chain_metadata){

  out <- .bt_marglik_manual_result(logml)
  out[["aggregation"]] <- list(
    rule = "exact_zero_dimensional",
    nonfinite_policy = "not_applicable",
    n_repetitions = 0L,
    n_included = 1L,
    n_failed = 0L
  )
  out[["diagnostics"]][["chains"]] <- chain_metadata
  out
}

.bt_marglik_repetition_field <- function(x, n_repetitions, name,
                                         missing = NA_real_){

  if(is.null(x)){
    return(rep(missing, n_repetitions))
  }
  if(length(x) == 1L && n_repetitions > 1L){
    return(rep(x, n_repetitions))
  }
  if(length(x) != n_repetitions){
    stop(
      "The upstream bridge result has ", length(x), " '", name,
      "' value(s) for ", n_repetitions, " repetition(s).",
      call. = FALSE
    )
  }
  x
}

.bt_marglik_append_detail <- function(current, rows, detail){

  if(!any(rows)){
    return(current)
  }
  current[rows] <- vapply(current[rows], function(existing){
    if(is.na(existing) || !nzchar(existing)){
      detail
    }else{
      paste(existing, detail, sep = " ")
    }
  }, character(1))
  current
}

.bt_marglik_from_upstream <- function(upstream, maxiter, nonfinite,
                                      chain_metadata,
                                      upstream_warnings = character()){

  if(!is.list(upstream) || is.null(upstream[["logml"]]) ||
     !is.numeric(upstream[["logml"]]) ||
     length(upstream[["logml"]]) == 0L){
    stop(
      "The upstream bridge sampler did not return numeric repetition-level 'logml' values.",
      call. = FALSE
    )
  }

  logml <- as.numeric(upstream[["logml"]])
  n_repetitions <- length(logml)
  niter <- as.numeric(.bt_marglik_repetition_field(
    upstream[["niter"]],
    n_repetitions,
    "niter"
  ))
  mcse <- as.numeric(.bt_marglik_repetition_field(
    upstream[["mcse_logml"]],
    n_repetitions,
    "mcse_logml"
  ))
  method <- as.character(.bt_marglik_repetition_field(
    upstream[["method"]],
    n_repetitions,
    "method",
    missing = "unknown"
  ))

  finite <- is.finite(logml)
  within_maxiter <- is.na(niter) | niter <= maxiter
  iteration_limit_exceeded <- !is.na(niter) & niter > maxiter
  success <- finite & within_maxiter
  warning_details <- rep(NA_character_, n_repetitions)
  error_details <- rep(NA_character_, n_repetitions)

  maxiter_message <- paste(
    "Marginal likelihood could not be estimated within the maximum number",
    "of iterations and might be more variable than usual."
  )
  warning_details <- .bt_marglik_append_detail(
    warning_details,
    iteration_limit_exceeded,
    maxiter_message
  )
  error_details <- .bt_marglik_append_detail(
    error_details,
    !finite,
    "The bridge sampler returned a non-finite log marginal likelihood."
  )

  repetitions <- data.frame(
    repetition = seq_len(n_repetitions),
    logml = logml,
    niter = niter,
    mcse = mcse,
    finite = finite,
    success = success,
    within_maxiter = within_maxiter,
    warning = warning_details,
    error = error_details,
    method = method,
    n_chains = rep(chain_metadata$count, n_repetitions),
    n_draws = rep(chain_metadata$draws, n_repetitions),
    stringsAsFactors = FALSE
  )

  if(any(!finite) && identical(nonfinite, "error")){
    failed <- repetitions[!finite, , drop = FALSE]
    condition <- errorCondition(
      message = paste0(
        "Bridge sampling returned a non-finite natural-log marginal likelihood ",
        "for repetition(s) ",
        paste(failed$repetition, collapse = ", "),
        ". Set 'nonfinite = \"drop\"' only when excluding failed repetitions is scientifically intended."
      ),
      call = NULL,
      class = "BayesTools_marglik_repetition_failure",
      repetitions = failed
    )
    stop(condition)
  }
  if(!any(finite)){
    stop(
      "Bridge sampling did not return any finite natural-log marginal likelihoods.",
      call. = FALSE
    )
  }

  if(any(!finite)){
    warning(
      "Dropped non-finite bridge-sampling repetition(s) ",
      paste(which(!finite), collapse = ", "),
      " by explicit request.",
      call. = FALSE
    )
  }

  out <- list(
    schema_version = .bt_marglik_result_version,
    logml = stats::median(logml[finite]),
    scale = "natural_log",
    repetitions = repetitions,
    aggregation = list(
      rule = "median_finite_logml",
      nonfinite_policy = nonfinite,
      n_repetitions = n_repetitions,
      n_included = sum(finite),
      n_failed = sum(!finite)
    ),
    diagnostics = list(
      upstream = upstream,
      upstream_warnings = unique(upstream_warnings),
      chains = chain_metadata
    )
  )
  class(out) <- c("BayesTools_marglik", "list")
  out
}

.bt_validate_marglik_result <- function(x, name = "marglik"){

  if(!inherits(x, "BayesTools_marglik")){
    stop(
      "'", name, "' must be a 'BayesTools_marglik' object created by ",
      "JAGS_bridgesampling() or bridgesampling_object().",
      call. = FALSE
    )
  }
  if(!identical(x[["schema_version"]], .bt_marglik_result_version)){
    stop(
      "'", name, "' uses an unsupported marginal-likelihood schema. ",
      "Recompute the marginal likelihood with the current BayesTools version.",
      call. = FALSE
    )
  }
  if(!identical(x[["scale"]], "natural_log")){
    stop(
      "'", name, "' must declare scale = \"natural_log\".",
      call. = FALSE
    )
  }
  logml <- x[["logml"]]
  if(!is.numeric(logml) || length(logml) != 1L || is.na(logml) ||
     is.nan(logml) || is.infinite(logml) && logml > 0){
    stop(
      "'", name, "$logml' must be one natural-log marginal likelihood ",
      "and may only be infinite when it is -Inf.",
      call. = FALSE
    )
  }
  invisible(x)
}

.bt_marglik_value <- function(x, name = "marglik"){

  .bt_validate_marglik_result(x, name)
  as.numeric(x[["logml"]])
}

.bt_JAGS_bridge_chain_metadata <- function(fit, posterior){

  chain_draws <- integer()

  if(inherits(fit, "mcmc.list")){
    chain_draws <- vapply(fit, nrow, integer(1))
  }else if(inherits(fit, "mcmc")){
    chain_draws <- nrow(fit)
  }else if(inherits(fit, "runjags")){
    chains <- tryCatch(
      coda::as.mcmc.list(fit),
      error = function(e) NULL
    )
    if(!is.null(chains)){
      chain_draws <- vapply(chains, nrow, integer(1))
    }
  }else if(is.list(fit) && length(fit) > 0L &&
           all(vapply(fit, inherits, logical(1), what = "mcarray"))){
    dimensions <- dim(fit[[1L]])
    if(length(dimensions) >= 2L){
      chain_draws <- rep.int(
        as.integer(dimensions[length(dimensions) - 1L]),
        as.integer(dimensions[length(dimensions)])
      )
    }
  }

  total_draws <- if(length(chain_draws) > 0L){
    sum(chain_draws)
  }else{
    nrow(posterior)
  }
  count <- if(length(chain_draws) > 0L){
    length(chain_draws)
  }else{
    NA_integer_
  }

  list(
    count = as.integer(count),
    draws = as.integer(total_draws),
    draws_per_chain = as.integer(chain_draws)
  )
}

#' @export
print.BayesTools_marglik <- function(x, ...){

  .bt_validate_marglik_result(x)
  cat("BayesTools marginal-likelihood result\n")
  cat("  logml (natural log): ", format(x[["logml"]]), "\n", sep = "")
  cat("  aggregation: ", x[["aggregation"]][["rule"]], "\n", sep = "")
  if(x[["aggregation"]][["n_repetitions"]] > 0L){
    cat(
      "  repetitions: ",
      x[["aggregation"]][["n_included"]],
      " included, ",
      x[["aggregation"]][["n_failed"]],
      " failed\n",
      sep = ""
    )
  }
  invisible(x)
}
