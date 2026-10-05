# Draw geometry and coordinate-based deterministic materialization.

.bt_draw_geometry_version <- 1L
.bt_backend_anchor_name <- "BayesTools_backend_anchor"

#' JAGS draw geometry and deterministic materialization
#'
#' @description
#' `JAGS_draw_geometry()` returns the retained iteration and chain geometry
#' stored from the actual chains returned by JAGS. Geometry is independent of
#' which coordinates are public.
#'
#' `JAGS_materialize_draws()` reconstructs requested coordinates from the
#' fitted parameter map. Sampled coordinates are selected from the fitted
#' chains and structural coordinates are filled with their exact `fixed_value`.
#' Internal coordinates, including a private backend anchor, are excluded by
#' default.
#'
#' `JAGS_with_draws()` returns a derived-draw view with supplied draws and
#' matching geometry, plus the complete untouched fit in `original_fit`.
#' Views have class `BayesTools_draws_view` and `BayesTools_fit`, with only
#' `original_fit` and `mcmc` slots. Nonruntime analysis attributes are copied;
#' view-local changes do not modify the original. Replacing a view retains its
#' original fit and current local analysis attributes without nesting views.
#' Fit-level draw metadata and stored posterior densities/ordinates are cleared;
#' attach fresh estimates with [posterior_metadata()] after replacement.
#' A view supports descriptive estimates, diagnostic plots, mixed-posterior
#' extraction and marginal inference from its supplied draws. It does not
#' support original-model convergence, bridge sampling or sampling extension:
#' use `fit$original_fit` and regenerate the view after extending that fit.
#' These refusals have common class `BayesTools_draws_view_unavailable`, with
#' leaves `BayesTools_draws_view_sampling_unavailable` (also
#' `BayesTools_sampling_unavailable`) and
#' `BayesTools_draws_view_inference_unavailable`.
#' Missing requested sampled coordinates signal
#' `BayesTools_draws_view_coordinate_unavailable`, with a `missing` field.
#' Views need not contain every original sampled coordinate; prior-draw
#' transformations retain their existing intersection-of-columns policy.
#' `coda::as.mcmc.list()` preserves chain timing. `coda::as.mcmc()` pools in
#' chain-major order with indices starting at one and thinning one; pooled
#' indices do not describe a continuous MCMC trajectory.
#'
#' @param fit fitted object created by [JAGS_fit()].
#' @param parameters optional exact vector of `coordinate_name` values.
#'   `NULL` selects every available public coordinate in map order.
#' @param include_internal whether internal coordinates may be returned.
#' @param draws replacement draws coercible to a `coda::mcmc.list`.
#'
#' @return `JAGS_draw_geometry()` returns a `BayesTools_draw_geometry` list.
#' `JAGS_draw_geometry_schema()` returns field descriptions.
#' `JAGS_materialize_draws()` returns a `coda::mcmc.list`, including a valid
#' zero-column list when the fit has no public coordinate.
#' `JAGS_with_draws()` returns a `BayesTools_draws_view`.
#'
#' @export JAGS_draw_geometry
#' @export JAGS_draw_geometry_schema
#' @export JAGS_materialize_draws
#' @export JAGS_with_draws
#' @name JAGS_draw_geometry
NULL

#' @rdname JAGS_draw_geometry
JAGS_draw_geometry <- function(fit){

  if(!inherits(fit, "BayesTools_fit")){
    stop("'fit' must be a 'BayesTools_fit' object.", call. = FALSE)
  }
  JAGS_validate_fit_contract(fit, requires = "draw_geometry")
  geometry <- attr(fit, "draw_geometry", exact = TRUE)
  .bt_validate_draw_geometry(geometry)
  geometry
}

#' @rdname JAGS_draw_geometry
JAGS_draw_geometry_schema <- function(){

  data.frame(
    field = c(
      "schema_version", "chains", "chain", "iterations", "start", "end",
      "thin", "draw_start", "draw_end", "total_draws", "chain_order"
    ),
    type = c(
      "integer", "data.frame", rep("integer", 7L), "integer", "integer"
    ),
    description = c(
      "Draw-geometry schema version.",
      "One row per retained chain.",
      "One-based chain number.",
      "Number of retained iterations in the chain.",
      "First retained MCMC iteration.",
      "Last retained MCMC iteration.",
      "MCMC thinning interval.",
      "First row occupied by the chain after chain-major concatenation.",
      "Last row occupied by the chain after chain-major concatenation.",
      "Total retained draws across chains.",
      "Explicit chain-major ordering."
    ),
    stringsAsFactors = FALSE
  )
}

#' @rdname JAGS_draw_geometry
JAGS_materialize_draws <- function(fit, parameters = NULL,
                                   include_internal = FALSE){

  check_char(parameters, "parameters", check_length = FALSE, allow_NULL = TRUE,
             allow_NA = FALSE)
  check_bool(include_internal, "include_internal", allow_NA = FALSE)
  geometry <- JAGS_draw_geometry(fit)
  coordinates <- parameter_coordinates(fit)

  available <- coordinates$monitor_status %in% c("sampled", "structural")
  if(!include_internal){
    available <- available & !coordinates$internal
  }
  if(is.null(parameters)){
    selected <- which(available)
  }else{
    if(anyDuplicated(parameters)){
      stop("'parameters' must not contain duplicates.", call. = FALSE)
    }
    matches <- match(parameters, coordinates$coordinate_name)
    if(anyNA(matches)){
      stop(
        "Unknown parameter coordinate",
        if(sum(is.na(matches)) > 1L) "s: " else ": ",
        paste0("'", parameters[is.na(matches)], "'", collapse = ", "),
        ".",
        call. = FALSE
      )
    }
    if(any(!available[matches])){
      stop(
        "Requested parameter is unavailable or internal. Set 'include_internal = TRUE' only for explicit backend diagnostics.",
        call. = FALSE
      )
    }
    selected <- matches
  }
  selected_coordinates <- coordinates[selected, , drop = FALSE]
  sampled <- which(selected_coordinates$monitor_status == "sampled")
  structural <- which(selected_coordinates$monitor_status != "sampled")
  sampled_names <- selected_coordinates$coordinate_name[sampled]
  fixed_values <- selected_coordinates$fixed_value[structural]

  chains <- .extract_posterior_samples(fit, as_list = TRUE)
  if(length(chains) != nrow(geometry$chains)){
    .bt_stop_refit_required(
      "The fitted chains disagree with the stored draw geometry. Refit the model with this version of BayesTools."
    )
  }
  out <- vector("list", length(chains))
  for(chain_i in seq_along(chains)){
    chain <- as.matrix(chains[[chain_i]])
    chain_geometry <- geometry$chains[chain_i, , drop = FALSE]
    if(nrow(chain) != chain_geometry$iterations){
      .bt_stop_refit_required(
        "The fitted chains disagree with the stored draw geometry. Refit the model with this version of BayesTools."
      )
    }
    values <- matrix(
      numeric(nrow(chain) * nrow(selected_coordinates)),
      nrow = nrow(chain),
      ncol = nrow(selected_coordinates),
      dimnames = list(NULL, selected_coordinates$coordinate_name)
    )
    columns <- match(sampled_names, colnames(chain))
    if(anyNA(columns)){
      if(inherits(fit, "BayesTools_draws_view")){
        missing <- sampled_names[is.na(columns)]
        stop(.bt_draws_view_condition(
          "BayesTools_draws_view_coordinate_unavailable",
          paste0("The derived-draw view does not contain requested sampled coordinates: ",
                 paste0("'", missing, "'", collapse = ", "),
                 ". Regenerate the view with those coordinates or use 'fit$original_fit'."),
          missing = missing
        ))
      }
      .bt_stop_refit_required(
        "A sampled parameter coordinate is missing from the fitted chains. Refit the model with this version of BayesTools."
      )
    }
    if(length(sampled) > 0L){
      values[, sampled] <- chain[, columns, drop = FALSE]
    }
    if(length(structural) > 0L){
      values[, structural] <- rep(fixed_values, each = nrow(chain))
    }
    out[[chain_i]] <- coda::mcmc(
      values,
      start = chain_geometry$start,
      end = chain_geometry$end,
      thin = chain_geometry$thin
    )
  }
  coda::mcmc.list(out)
}

#' @rdname JAGS_draw_geometry
JAGS_with_draws <- function(fit, draws){

  if(!inherits(fit, "BayesTools_fit")){
    stop("'fit' must be a 'BayesTools_fit' object.", call. = FALSE)
  }
  draws <- coda::as.mcmc.list(draws)
  if(is.null(fit[["mcmc"]])){
    stop("'fit' has no replaceable 'mcmc' component.", call. = FALSE)
  }
  original <- if(inherits(fit, "BayesTools_draws_view")) fit[["original_fit"]] else fit
  if(!inherits(original, "BayesTools_fit") ||
     inherits(original, "BayesTools_draws_view") || is.null(original[["mcmc"]])){
    stop("'fit$original_fit' must be a 'BayesTools_fit' with an 'mcmc' component and must not be a derived-draw view.", call. = FALSE)
  }
  out <- structure(list(original_fit = original, mcmc = draws),
                    class = c("BayesTools_draws_view", "BayesTools_fit"))
  excluded <- c("names", "class", "runtime_setup", "runtime_cache", "runtime_state",
                 "bayestools_meta", "posterior_density", "posterior_densities",
                 "posterior_ordinate", "posterior_ordinates")
  analysis_attributes <- attributes(fit)
  for(name in setdiff(names(analysis_attributes), excluded)){
    attr(out, name) <- analysis_attributes[[name]]
  }
  attr(out, "draw_geometry") <- .bt_draw_geometry_from_chains(draws)
  out
}

.bt_is_jags_analysis_fit <- function(fit){

  inherits(fit, c("runjags", "BayesTools_draws_view"))
}

.bt_draws_view_condition <- function(class, message, ...){

  errorCondition(message, class = c(class, "BayesTools_draws_view_unavailable"),
                 call = NULL, ...)
}

#' @rdname JAGS_draw_geometry
#' @param x a derived-draw view.
#' @param ... additional arguments.
#' @exportS3Method coda::as.mcmc.list
as.mcmc.list.BayesTools_draws_view <- function(x, ...){

  coda::as.mcmc.list(x[["mcmc"]])
}

#' @rdname JAGS_draw_geometry
#' @exportS3Method coda::as.mcmc
as.mcmc.BayesTools_draws_view <- function(x, ...){

  chains <- coda::as.mcmc.list(x)
  coda::mcmc(do.call(rbind, lapply(chains, as.matrix)), start = 1, thin = 1)
}

#' @rdname JAGS_draw_geometry
#' @exportS3Method
print.BayesTools_draws_view <- function(x, ...){

  geometry <- JAGS_draw_geometry(x)
  cat("Derived-draw view:", nrow(geometry$chains), "chains,",
      geometry$total_draws, "draws.\n",
      "The complete original sampling fit is available as 'fit$original_fit'.\n")
  invisible(x)
}

.bt_draw_geometry_from_chains <- function(chains){

  chains <- coda::as.mcmc.list(chains)
  if(length(chains) == 0L){
    stop("Cannot record draw geometry without at least one retained chain.",
         call. = FALSE)
  }
  chain_rows <- lapply(seq_along(chains), function(chain_i){
    chain <- chains[[chain_i]]
    mcpar <- attr(chain, "mcpar", exact = TRUE)
    if(is.null(mcpar) || length(mcpar) != 3L){
      stop("Retained chains have missing or malformed MCMC timing metadata.",
           call. = FALSE)
    }
    data.frame(
      chain = as.integer(chain_i),
      iterations = as.integer(nrow(chain)),
      start = as.integer(mcpar[1L]),
      end = as.integer(mcpar[2L]),
      thin = as.integer(mcpar[3L]),
      stringsAsFactors = FALSE
    )
  })
  chain_table <- do.call(rbind, chain_rows)
  rownames(chain_table) <- NULL
  draw_end <- cumsum(chain_table$iterations)
  chain_table$draw_start <- as.integer(draw_end - chain_table$iterations + 1L)
  chain_table$draw_end <- as.integer(draw_end)
  geometry <- list(
    schema_version = .bt_draw_geometry_version,
    chains = chain_table,
    total_draws = as.integer(sum(chain_table$iterations)),
    chain_order = as.integer(chain_table$chain)
  )
  class(geometry) <- c("BayesTools_draw_geometry", "list")
  .bt_validate_draw_geometry(geometry)
  geometry
}

.bt_validate_draw_geometry <- function(geometry){

  fields <- c("schema_version", "chains", "total_draws", "chain_order")
  chain_fields <- c(
    "chain", "iterations", "start", "end", "thin", "draw_start", "draw_end"
  )
  valid <- inherits(geometry, "BayesTools_draw_geometry") &&
    is.list(geometry) && identical(names(geometry), fields) &&
    identical(geometry$schema_version, .bt_draw_geometry_version) &&
    is.data.frame(geometry$chains) &&
    identical(names(geometry$chains), chain_fields) &&
    nrow(geometry$chains) > 0L
  if(!valid){
    .bt_stop_refit_required(
      "Draw geometry is missing, malformed, or unsupported. Refit the model with this version of BayesTools."
    )
  }
  integer_fields <- vapply(geometry$chains, is.integer, logical(1))
  if(!all(integer_fields) || anyNA(geometry$chains) ||
     any(geometry$chains$iterations < 1L) || any(geometry$chains$start < 1L) ||
     any(geometry$chains$thin < 1L) ||
     !identical(geometry$chains$chain, seq_len(nrow(geometry$chains))) ||
     !identical(geometry$chain_order, geometry$chains$chain) ||
     any(geometry$chains$end != geometry$chains$start +
           (geometry$chains$iterations - 1L) * geometry$chains$thin) ||
     !identical(geometry$chains$draw_start,
                c(1L, head(geometry$chains$draw_end, -1L) + 1L)) ||
     !identical(geometry$chains$draw_end,
                as.integer(cumsum(geometry$chains$iterations))) ||
     !is.integer(geometry$total_draws) || length(geometry$total_draws) != 1L ||
     !identical(geometry$total_draws, sum(geometry$chains$iterations))){
    .bt_stop_refit_required(
      "Draw geometry contains inconsistent chain timing or ordering. Refit the model with this version of BayesTools."
    )
  }
  invisible(TRUE)
}

.bt_attach_draw_geometry <- function(fit){

  if(inherits(fit, "error")){
    return(fit)
  }
  chains <- .extract_posterior_samples(fit, as_list = TRUE)
  attr(fit, "draw_geometry") <- .bt_draw_geometry_from_chains(chains)
  fit
}

.bt_add_backend_anchor <- function(model_syntax, data, prior_list,
                                   add_parameters, monitor){

  monitor <- unique(monitor[nzchar(monitor)])
  if(length(monitor) > 0L){
    return(list(
      model_syntax = model_syntax,
      monitor = monitor,
      backend_anchor = NULL
    ))
  }
  anchor <- .bt_backend_anchor_name
  occupied <- unique(c(names(data), names(prior_list), add_parameters))
  token_pattern <- paste0(
    "(^|[^A-Za-z0-9_.])", JAGS_regex_escape(anchor),
    "([^A-Za-z0-9_.]|$)"
  )
  if(anchor %in% occupied || grepl(token_pattern, model_syntax, perl = TRUE)){
    stop(
      "The model uses reserved BayesTools backend node '", anchor, "'.",
      call. = FALSE
    )
  }
  opening_bracket <- regexpr("{", model_syntax, fixed = TRUE)[1L]
  if(opening_bracket < 1L){
    stop("The JAGS model syntax has no opening model brace.", call. = FALSE)
  }
  syntax_start <- substr(model_syntax, 1L, opening_bracket)
  syntax_end <- substr(model_syntax, opening_bracket + 1L, nchar(model_syntax))
  model_syntax <- paste0(
    syntax_start, "\n  ", anchor, " <- 0\n", syntax_end
  )
  list(
    model_syntax = model_syntax,
    monitor = anchor,
    backend_anchor = anchor
  )
}
