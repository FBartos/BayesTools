#' @title Assess convergence of a runjags model
#'
#' @description Checks whether the supplied \link[runjags]{runjags-package} model
#' satisfied convergence criteria.
#' @param fit a 'BayesTools_fit' object created by [JAGS_fit()]
#' @param prior_list named list of prior distribution
#' (names correspond to the parameter names). Retained for compatibility and
#' validated when supplied; the classification of the fitted parameters comes
#' from the parameter map stored with \code{fit} (see Details).
#' @param max_Rhat maximum R-hat error for the autofit function.
#'   Defaults to \code{1.05}. With one chain, this criterion is skipped with
#'   a warning; the remaining enabled criteria are still assessed.
#' @param min_ESS minimum effective sample size. Defaults to \code{500}.
#' @param max_error maximum MCMC error. Defaults to \code{0.01}.
#' @param max_SD_error maximum MCMC error as the proportion of standard
#'   deviation of the parameters. Defaults to \code{0.05}.
#' @param add_parameters vector of parameter names that are excluded
#' from the default selection. Parameters named in \code{monitor} are
#' checked even when they are listed here. Automatic fitting in
#' \code{JAGS_fit()} and \code{JAGS_extend()} excludes nothing through this
#' argument; it uses the default selection described in Details.
#' @param fail_fast whether the function should stop after the first failed convergence check.
#' @param check_indicators whether model indicator variables should be included
#' in convergence checks. Binary indicators are checked as Bernoulli
#' occupancies and categorical indicators are checked separately for every
#' observed state. When \code{monitor} is supplied, eligible indicators are
#' added to that selection. Defaults to \code{FALSE}.
#' @param monitor optional character vector selecting parameters for convergence
#' checks. A base name selects all of its indexed elements. Requests are
#' resolved against all monitored parameters, including derived and auxiliary
#' ones that the default selection leaves out.
#' \code{NULL} selects the default parameters; \code{character()} requests no
#' parameters.
#' @param allow_not_assessable whether requested sampled parameters with
#' undefined diagnostics may be ignored. Defaults to \code{FALSE}. A sampled
#' column that never changes, including a constant model indicator, is not
#' evidence of convergence and remains not assessable, whatever its draws.
#'
#' @details \code{JAGS_fit()} assigns every fitted coordinate one convergence
#' role, stored in the \code{convergence_role} column of
#' [parameter_coordinates()] and derived from the declared model, never from
#' the draws:
#' \describe{
#'   \item{\code{"sampled"}}{stochastic parameters, including the user's
#'   \code{add_parameters} that depend on a stochastic node; checked by
#'   default.}
#'   \item{\code{"indicator"}}{declared model indicators of mixture,
#'   spike-and-slab, and variance-allocation priors; checked with
#'   \code{check_indicators = TRUE}.}
#'   \item{\code{"structural"}}{declared constants: point priors, one-state
#'   indicators, reference and fixed publication-weight bins, a p-hacking kind
#'   shared by every mixture branch, the point total of an ordered prior and
#'   the coefficients it fixes, unit correlation diagonals and Cholesky
#'   constants, and monitored deterministic nodes (\code{<-} or \code{=})
#'   whose ancestors in the model syntax are all data or constants. They are
#'   reported as \code{"structural_constant"}.}
#'   \item{\code{"derived"}}{deterministic functions of sampled nodes that
#'   BayesTools generates for formulas (such as random-effect correlation
#'   matrices, their Cholesky factors, and derived latent effects or SDs) and
#'   point priors whose location is an expression; checked only when named in
#'   \code{monitor}.}
#'   \item{\code{"auxiliary"}}{inclusion probabilities, the backend anchor of
#'   models without monitored parameters, and private implementation nodes;
#'   checked only when named in \code{monitor}.}
#' }
#' The default selection (\code{monitor = NULL}) consists of the sampled and
#' structural parameters, and of the indicators when \code{check_indicators}
#' is set. Automatic fitting in \code{JAGS_fit()} and \code{JAGS_extend()}
#' uses the same selection.
#'
#' @examples \dontrun{
#' # simulate data
#' set.seed(1)
#' data <- list(
#'   x = rnorm(10),
#'   N = 10
#' )
#' data$x
#'
#' # define priors
#' priors_list <- list(mu = prior("normal", list(0, 1)))
#'
#' # define likelihood for the data
#' model_syntax <-
#'   "model{
#'     for(i in 1:N){
#'       x[i] ~ dnorm(mu, 1)
#'     }
#'   }"
#'
#' # fit the models
#' fit <- JAGS_fit(model_syntax, data, priors_list)
#' JAGS_check_convergence(fit)
#' }
#' @return \code{JAGS_check_convergence} returns a boolean indicating whether
#' all requested, assessable parameters satisfy the enabled criteria. An
#' explicitly empty \code{monitor} returns \code{logical(0)} rather than
#' claiming convergence. When no parameters remain available (empty sample
#' columns and no structural parameters) or the selection is empty, the
#' function returns \code{TRUE}: there is nothing assessable, which is
#' treated as vacuously satisfied. The \code{diagnostics} attribute contains
#' one row per available parameter and classifies it as
#' \code{"assessable"}, \code{"structural_constant"}, \code{"not_assessable"},
#' \code{"not_requested"}, or \code{"not_checked"} when \code{fail_fast = TRUE}
#' stops before reaching that parameter. The \code{errors} attribute carries
#' failed checks.
#'
#' @seealso [JAGS_fit()] [parameter_coordinates()]
#' @export
JAGS_check_convergence <- function(
    fit,
    prior_list = NULL,
    max_Rhat = 1.05,
    min_ESS = 500,
    max_error = 0.01,
    max_SD_error = 0.05,
    add_parameters = NULL,
    fail_fast = FALSE,
    check_indicators = FALSE,
    monitor = NULL,
    allow_not_assessable = FALSE){

  # check input
  if(!inherits(fit, "runjags"))
    stop("'fit' must be a runjags fit", call. = FALSE)
  check_list(prior_list, "prior_list", allow_NULL = TRUE)
  if(!is.null(prior_list) && any(!vapply(prior_list, is.prior, logical(1))))
    stop("'prior_list' must be a list of priors.", call. = FALSE)
  check_real(max_Rhat,     "max_Rhat",     lower = 1, allow_NULL = TRUE)
  check_real(min_ESS,      "min_ESS",      lower = 0, allow_NULL = TRUE)
  check_real(max_error,    "max_error",    lower = 0, allow_NULL = TRUE)
  check_real(max_SD_error, "max_SD_error", lower = 0, upper = 1, allow_NULL = TRUE)
  check_char(add_parameters, "add_parameters", check_length = 0, allow_NULL = TRUE)
  check_bool(fail_fast, "fail_fast", allow_NA = FALSE)
  check_bool(check_indicators, "check_indicators", allow_NA = FALSE)
  check_char(monitor, "monitor", check_length = 0, allow_NULL = TRUE,
             allow_NA = FALSE)
  check_bool(allow_not_assessable, "allow_not_assessable", allow_NA = FALSE)
  if(!inherits(fit, "BayesTools_fit")){
    stop(
      "'fit' must be a 'BayesTools_fit' created by JAGS_fit(). Refit the model with this version of BayesTools.",
      call. = FALSE
    )
  }

  .bt_check_convergence(
    fit = fit,
    coordinates = parameter_coordinates(fit),
    max_Rhat = max_Rhat,
    min_ESS = min_ESS,
    max_error = max_error,
    max_SD_error = max_SD_error,
    add_parameters = add_parameters,
    fail_fast = fail_fast,
    check_indicators = check_indicators,
    monitor = monitor,
    allow_not_assessable = allow_not_assessable
  )
}

# Convergence check of a fit against its coordinate table. Automatic fitting
# calls it with the coordinates it builds after the first sampling run, so
# that fitting, extension, and the public check share one selection.
.bt_check_convergence <- function(
    fit,
    coordinates,
    max_Rhat,
    min_ESS,
    max_error,
    max_SD_error,
    add_parameters = NULL,
    fail_fast = FALSE,
    check_indicators = FALSE,
    monitor = NULL,
    allow_not_assessable = FALSE){

  prepared <- .bt_convergence_prepare(
    fit = fit,
    coordinates = coordinates,
    add_parameters = add_parameters,
    monitor = monitor
  )
  mcmc_samples_list <- prepared$mcmc_samples_list
  targets           <- prepared$targets
  metadata          <- targets$metadata

  explicitly_empty <- !is.null(monitor) && length(monitor) == 0L
  if(explicitly_empty){
    diagnostics <- .bt_convergence_diagnostics(character())
    diagnostics[["assessable"]] <- NULL
    return(.bt_convergence_result(logical(), diagnostics, NULL))
  }

  if(is.null(monitor)){
    selected_parameters <- metadata$parameter[
      metadata$role %in% c("sampled", "structural") |
        (check_indicators & metadata$role == "indicator")
    ]
  }else{
    selected_parameters <- .bt_convergence_resolve_monitor(
      monitor,
      metadata$parameter,
      metadata$source
    )
    if(check_indicators){
      selected_parameters <- unique(c(
        selected_parameters,
        metadata$parameter[metadata$role == "indicator"]
      ))
    }
  }

  if(nrow(metadata) == 0L){
    diagnostics <- .bt_convergence_diagnostics(character())
    diagnostics[["assessable"]] <- NULL
    return(.bt_convergence_result(TRUE, diagnostics, NULL))
  }

  diagnostics <- .bt_convergence_diagnostics(metadata$parameter)
  selected_rows <- match(selected_parameters, diagnostics[["parameter"]])
  diagnostics[["state"]][selected_rows] <- "assessable"
  structural_rows <- metadata$role == "structural" &
    diagnostics[["parameter"]] %in% selected_parameters
  diagnostics[["state"]][structural_rows] <- "structural_constant"

  sample_rows <- selected_rows[
    diagnostics[["state"]][selected_rows] == "assessable"
  ]
  if(length(sample_rows) > 0L){
    diagnostics[["state"]][sample_rows] <- "not_checked"
    if(length(mcmc_samples_list) == 1L && !is.null(max_Rhat)){
      warning(
        "Only one chain was run. R-hat cannot be computed; checking the remaining enabled convergence criteria.",
        call. = FALSE
      )
      max_Rhat <- NULL
    }
    chain_lengths <- vapply(mcmc_samples_list, nrow, integer(1))
    chain_id <- rep.int(seq_along(chain_lengths), chain_lengths)
    for(row in sample_rows){
      diagnostics[["state"]][[row]] <- "assessable"
      parameter_diagnostics <- .bt_convergence_parameter_diagnostics(
        targets$samples[[row]],
        chain_id = chain_id,
        n_chains = length(mcmc_samples_list),
        assess_Rhat = !is.null(max_Rhat),
        assess_ESS = !is.null(min_ESS),
        assess_error = !is.null(max_error),
        assess_SD_error = !is.null(max_SD_error)
      )
      diagnostics[row, names(parameter_diagnostics)] <-
        parameter_diagnostics
      if(!isTRUE(parameter_diagnostics[["assessable"]])){
        diagnostics[["state"]][[row]] <- "not_assessable"
      }
      if(fail_fast && length(.bt_convergence_failures(
        diagnostics = diagnostics[row, , drop = FALSE],
        max_Rhat = max_Rhat, min_ESS = min_ESS,
        max_error = max_error, max_SD_error = max_SD_error,
        allow_not_assessable = allow_not_assessable
      )) > 0L){
        break
      }
    }
  }

  diagnostics[["assessable"]] <- NULL
  fails <- .bt_convergence_failures(
    diagnostics = diagnostics,
    max_Rhat = max_Rhat,
    min_ESS = min_ESS,
    max_error = max_error,
    max_SD_error = max_SD_error,
    allow_not_assessable = allow_not_assessable
  )
  if(fail_fast && length(fails) > 1L){
    fails <- fails[[1L]]
  }

  .bt_convergence_result(length(fails) == 0L, diagnostics, fails)
}

# One convergence target per displayed parameter: the visible posterior
# columns (labelled as in the summary tables) with the convergence role of
# their coordinate, followed by the structural coordinates without draws.
# Indicators are checked as state occupancies.
.bt_convergence_prepare <- function(fit, coordinates, add_parameters, monitor){

  mcmc_samples_list <- .extract_posterior_samples(fit, as_list = TRUE)
  mcmc_samples      <- do.call(rbind, mcmc_samples_list)
  columns <- colnames(mcmc_samples)
  if(is.null(columns)){
    columns <- character()
  }

  visible <- .bt_convergence_visible_columns(
    columns = columns,
    prior_list = attr(fit, "prior_list", exact = TRUE),
    remove_parameters = .bt_convergence_excluded_add_parameters(
      add_parameters,
      monitor
    )
  )
  roles <- coordinates$convergence_role[
    match(visible$column, coordinates$coordinate_name)
  ]
  roles[is.na(roles)] <- "sampled"
  unmonitored <- coordinates$coordinate_name[
    coordinates$convergence_role == "structural" &
      !coordinates$coordinate_name %in% columns
  ]

  samples  <- list()
  metadata <- list()
  add_target <- function(parameter, source, values, role){
    samples[[length(samples) + 1L]] <<- values
    metadata[[length(metadata) + 1L]] <<- data.frame(
      parameter = parameter,
      source = source,
      role = role,
      stringsAsFactors = FALSE
    )
  }

  for(i in seq_len(nrow(visible))){
    label  <- visible$label[[i]]
    values <- mcmc_samples[, visible$column[[i]]]
    if(!identical(roles[[i]], "indicator")){
      add_target(label, label, values, roles[[i]])
      next
    }

    if(any(!is.finite(values)) || any(values != round(values))){
      stop(
        "Model indicator '", label,
        "' must contain finite integer states.",
        call. = FALSE
      )
    }
    observed <- sort(unique(values))
    if(length(observed) <= 1L){
      add_target(label, label, values, "indicator")
    }else if(length(observed) == 2L){
      add_target(label, label, as.numeric(values == observed[[2L]]), "indicator")
    }else{
      for(state in observed){
        state_label <- format(
          state,
          digits = 17,
          scientific = FALSE,
          trim = TRUE
        )
        add_target(
          paste0(label, " (state ", state_label, ")"),
          label,
          as.numeric(values == state),
          "indicator"
        )
      }
    }
  }
  for(parameter in setdiff(unmonitored, visible$label)){
    add_target(parameter, parameter, NULL, "structural")
  }

  if(length(metadata) == 0L){
    metadata <- data.frame(
      parameter = character(),
      source = character(),
      role = character(),
      stringsAsFactors = FALSE
    )
  }else{
    metadata <- do.call(rbind, metadata)
    rownames(metadata) <- NULL
  }

  list(
    mcmc_samples_list = mcmc_samples_list,
    targets = list(samples = samples, metadata = metadata)
  )
}

.bt_convergence_monitor_base <- function(parameters){

  sub("\\[.*$", "", sub(" \\(state [^)]*\\)$", "", parameters))
}

.bt_convergence_excluded_add_parameters <- function(add_parameters, monitor){

  if(length(add_parameters) == 0L || length(monitor) == 0L){
    return(add_parameters)
  }

  requested <- .bt_convergence_monitor_base(monitor)
  add_parameters[!.bt_convergence_monitor_base(add_parameters) %in% requested]
}

# Resolve an explicit convergence monitor against the fitted columns before
# further sampling, so that an unknown request fails without discarding work.
.bt_convergence_validate_monitor <- function(fit, coordinates, monitor){

  if(length(monitor) == 0L){
    return(invisible(TRUE))
  }

  metadata <- .bt_convergence_prepare(
    fit = fit,
    coordinates = coordinates,
    add_parameters = NULL,
    monitor = monitor
  )$targets$metadata
  .bt_convergence_resolve_monitor(
    monitor,
    metadata$parameter,
    metadata$source
  )

  invisible(TRUE)
}

# Before the first sampling run only the monitored node names are known.
# Validate the requested base names against them; indexed elements are
# resolved against the fitted columns by JAGS_check_convergence().
.bt_convergence_validate_monitor_names <- function(monitor, monitored_names){

  if(length(monitor) == 0L){
    return(invisible(TRUE))
  }

  monitored_names <- monitored_names[!is.na(monitored_names) & nzchar(monitored_names)]
  available <- unique(.bt_convergence_monitor_base(monitored_names))
  for(parameter in unique(monitor)){
    if(!.bt_convergence_monitor_base(parameter) %in% available){
      stop(
        "The requested convergence monitor '", parameter,
        "' is not monitored by the model.",
        call. = FALSE
      )
    }
  }

  invisible(TRUE)
}

.bt_convergence_resolve_monitor <- function(monitor, available_parameters,
                                            available_sources = available_parameters){

  selected <- character()
  for(parameter in unique(monitor)){
    if(grepl("\\[", parameter, fixed = FALSE)){
      matches <- available_sources == parameter |
        available_parameters == parameter
    }else{
      matches <- sub("\\[.*$", "", available_sources) == parameter |
        available_parameters == parameter
    }
    if(!any(matches)){
      stop(
        "The requested convergence monitor '", parameter,
        "' is not available in the fitted model.",
        call. = FALSE
      )
    }
    selected <- c(selected, available_parameters[matches])
  }
  unique(selected)
}

.bt_convergence_diagnostics <- function(parameters){

  diagnostics <- data.frame(
    parameter = parameters,
    state = rep.int("not_requested", length(parameters)),
    Rhat = rep.int(NA_real_, length(parameters)),
    ESS = rep.int(NA_real_, length(parameters)),
    MCMC_error = rep.int(NA_real_, length(parameters)),
    MCMC_SD_error = rep.int(NA_real_, length(parameters)),
    assessable = rep.int(NA, length(parameters)),
    stringsAsFactors = FALSE
  )
  class(diagnostics) <- c(
    "BayesTools_convergence_diagnostics",
    "data.frame"
  )
  diagnostics
}

.bt_convergence_result <- function(converged, diagnostics, errors){

  if(length(errors) == 0L){
    errors <- NULL
  }
  attr(converged, "diagnostics") <- diagnostics
  attr(converged, "errors") <- errors
  converged
}

.bt_convergence_parameter_diagnostics <- function(
    samples,
    chain_id,
    n_chains,
    assess_Rhat,
    assess_ESS,
    assess_error,
    assess_SD_error){

  output <- list(
    Rhat = NA_real_,
    ESS = NA_real_,
    MCMC_error = NA_real_,
    MCMC_SD_error = NA_real_,
    assessable = TRUE
  )
  enabled <- c(
    Rhat = assess_Rhat && n_chains >= 2L,
    ESS = assess_ESS,
    MCMC_error = assess_error,
    MCMC_SD_error = assess_SD_error
  )
  if(!any(enabled)){
    return(output)
  }
  if(any(!is.finite(samples)) || length(unique(samples)) < 2L){
    output[["assessable"]] <- FALSE
    return(output)
  }

  chains <- lapply(seq_len(n_chains), function(chain){
    values <- samples[chain_id == chain]
    coda::as.mcmc(matrix(
      values,
      ncol = 1L,
      dimnames = list(NULL, "parameter")
    ))
  })
  mcmc_samples <- coda::as.mcmc.list(chains)

  if(assess_Rhat && n_chains >= 2L){
    output[["Rhat"]] <- tryCatch(
      {
        psrf <- suppressWarnings(coda::gelman.diag(
          mcmc_samples,
          multivariate = FALSE,
          autoburnin = FALSE
        )[["psrf"]])
        max(psrf[1L, ], na.rm = FALSE)
      },
      error = function(e) NA_real_
    )
  }
  if(assess_ESS){
    output[["ESS"]] <- tryCatch(
      suppressWarnings(as.numeric(coda::effectiveSize(mcmc_samples)[[1L]])),
      error = function(e) NA_real_
    )
    if(identical(output[["ESS"]], 0)){
      output[["ESS"]] <- NA_real_
    }
  }
  if(assess_error || assess_SD_error){
    statistics <- tryCatch(
      suppressWarnings(summary(mcmc_samples, quantiles = NULL)[["statistics"]]),
      error = function(e) NULL
    )
    if(!is.null(statistics)){
      if(is.null(dim(statistics))){
        statistics <- t(statistics)
      }
      time_series_se <- statistics[1L, "Time-series SE"]
      standard_deviation <- statistics[1L, "SD"]
      if(assess_error){
        output[["MCMC_error"]] <- time_series_se
      }
      if(assess_SD_error){
        output[["MCMC_SD_error"]] <- time_series_se / standard_deviation
      }
    }
  }

  enabled_names <- names(enabled)[enabled]
  output[["assessable"]] <- all(is.finite(unlist(output[enabled_names])))
  if(!output[["assessable"]]){
    for(name in enabled_names){
      if(!is.finite(output[[name]])){
        output[[name]] <- NA_real_
      }
    }
  }
  output
}

.bt_convergence_failures <- function(
    diagnostics,
    max_Rhat,
    min_ESS,
    max_error,
    max_SD_error,
    allow_not_assessable){

  fails <- character()
  not_assessable <- diagnostics[["state"]] == "not_assessable"
  if(any(not_assessable) && !allow_not_assessable){
    enabled <- c(
      Rhat = !is.null(max_Rhat),
      ESS = !is.null(min_ESS),
      MCMC_error = !is.null(max_error),
      MCMC_SD_error = !is.null(max_SD_error)
    )
    labels <- c(
      Rhat = "R-hat",
      ESS = "ESS",
      MCMC_error = "MCMC error",
      MCMC_SD_error = "MCMC SD error"
    )
    for(row in which(not_assessable)){
      unavailable <- names(enabled)[
        enabled & !is.finite(unlist(diagnostics[row, names(enabled)]))
      ]
      fails <- c(
        fails,
        paste0(
          paste(labels[unavailable], collapse = ", "),
          if(length(unavailable) == 1L) " diagnostic is" else " diagnostics are",
          " not assessable for '",
          diagnostics[["parameter"]][[row]],
          "'."
        )
      )
    }
  }
  assessable <- diagnostics[["state"]] == "assessable"
  checks <- list(
    list(
      column = "Rhat",
      threshold = max_Rhat,
      fails = function(x, target) x > target,
      message = function(x, target, parameter){
        paste0(
          "R-hat ", round(x, 3), " for '", parameter,
          "' is larger than the set target (", target, ")."
        )
      }
    ),
    list(
      column = "ESS",
      threshold = min_ESS,
      fails = function(x, target) x < target,
      message = function(x, target, parameter){
        paste0(
          "ESS ", round(x), " for '", parameter,
          "' is lower than the set target (", target, ")."
        )
      }
    ),
    list(
      column = "MCMC_error",
      threshold = max_error,
      fails = function(x, target) x > target,
      message = function(x, target, parameter){
        paste0(
          "MCMC error ", round(x, 5), " for '", parameter,
          "' is larger than the set target (", target, ")."
        )
      }
    ),
    list(
      column = "MCMC_SD_error",
      threshold = max_SD_error,
      fails = function(x, target) x > target,
      message = function(x, target, parameter){
        paste0(
          "MCMC SD error ", round(x, 3), " for '", parameter,
          "' is larger than the set target (", target, ")."
        )
      }
    )
  )
  for(check in checks){
    if(is.null(check[["threshold"]])){
      next
    }
    failed <- assessable & check[["fails"]](
      diagnostics[[check[["column"]]]],
      check[["threshold"]]
    )
    if(any(failed)){
      fails <- c(
        fails,
        mapply(
          FUN = check[["message"]],
          x = diagnostics[[check[["column"]]]][failed],
          parameter = diagnostics[["parameter"]][failed],
          MoreArgs = list(target = check[["threshold"]]),
          USE.NAMES = FALSE
        )
      )
    }
  }
  fails
}
