#' @title Assess convergence of a runjags model
#'
#' @description Checks whether the supplied \link[runjags]{runjags-package} model
#' satisfied convergence criteria.
#' @param fit a runjags model
#' @param prior_list named list of prior distribution
#' (names correspond to the parameter names)
#' @param max_Rhat maximum R-hat error for the autofit function.
#'   Defaults to \code{1.05}.
#' @param min_ESS minimum effective sample size. Defaults to \code{500}.
#' @param max_error maximum MCMC error. Defaults to \code{0.01}.
#' @param max_SD_error maximum MCMC error as the proportion of standard
#'   deviation of the parameters. Defaults to \code{0.05}.
#' @param add_parameters vector of additional parameter names that should be used
#' (only allows removing last, fixed, omega element if omega is tracked manually).
#' @param fail_fast whether the function should stop after the first failed convergence check.
#' @param check_indicators whether model indicator variables should be included
#' in convergence checks. Binary indicators are checked as Bernoulli
#' occupancies and categorical indicators are checked separately for every
#' observed state. Defaults to \code{FALSE}.
#' @param monitor optional character vector selecting parameters for convergence
#' checks. A base name selects all of its indexed elements. \code{NULL} selects
#' every eligible parameter; \code{character()} requests no parameters.
#' @param allow_not_assessable whether requested sampled parameters with
#' undefined diagnostics may be ignored. Defaults to \code{FALSE}.
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
#' JAGS_check_convergence(fit, priors_list)
#' }
#' @return \code{JAGS_check_convergence} returns a boolean indicating whether
#' all requested, assessable parameters satisfy the enabled criteria. An
#' explicitly empty \code{monitor} returns \code{logical(0)} rather than
#' claiming convergence. When no parameters remain available after cleaning
#' (empty sample columns and no structural priors), the function returns
#' \code{TRUE}: there is nothing assessable, which is treated as vacuously
#' satisfied rather than as an empty selection. The \code{diagnostics}
#' attribute contains one row per available parameter and classifies it as
#' \code{"assessable"}, \code{"structural_constant"}, \code{"not_assessable"},
#' or \code{"not_requested"}. The \code{errors} attribute carries failed checks.
#'
#' @seealso [JAGS_fit()]
#' @export
JAGS_check_convergence <- function(
    fit,
    prior_list,
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

  # extract samples and parameter information
  mcmc_samples_list <- .extract_posterior_samples(fit, as_list = TRUE)
  mcmc_samples      <- do.call(rbind, mcmc_samples_list)

  # Remove parameters that are intentionally excluded from automatic checks.
  # Structural point parameters are added back below from prior metadata.
  remove_params <- c(
    names(prior_list)[vapply(
      prior_list,
      .bt_convergence_is_structural_prior,
      logical(1)
    )],
    add_parameters
  )

  cleaned <- .remove_auxiliary_parameters(mcmc_samples, prior_list, remove_params)
  mcmc_samples <- cleaned$model_samples

  targets <- .bt_convergence_sample_targets(mcmc_samples, prior_list)
  sample_parameters <- targets$metadata$parameter
  structural_parameters <- .bt_convergence_structural_parameters(prior_list)
  available_parameters <- unique(c(sample_parameters, structural_parameters))
  available_sources <- c(
    targets$metadata$source,
    structural_parameters[!structural_parameters %in% sample_parameters]
  )

  explicitly_empty <- !is.null(monitor) && length(monitor) == 0L
  if(explicitly_empty){
    diagnostics <- .bt_convergence_diagnostics(character())
    diagnostics[["assessable"]] <- NULL
    return(.bt_convergence_result(logical(), diagnostics, NULL))
  }

  if(is.null(monitor)){
    selected_parameters <- targets$metadata$parameter[
      !(targets$metadata$is_inclusion |
          (!check_indicators & targets$metadata$is_indicator))
    ]
    selected_parameters <- unique(c(
      selected_parameters,
      structural_parameters
    ))
  }else{
    selected_parameters <- .bt_convergence_resolve_monitor(
      monitor,
      available_parameters,
      available_sources
    )
  }

  if(length(available_parameters) == 0L){
    diagnostics <- .bt_convergence_diagnostics(character())
    diagnostics[["assessable"]] <- NULL
    return(.bt_convergence_result(TRUE, diagnostics, NULL))
  }

  diagnostics <- .bt_convergence_diagnostics(available_parameters)
  selected_rows <- match(selected_parameters, diagnostics[["parameter"]])
  diagnostics[["state"]][selected_rows] <- "assessable"
  structural_rows <- diagnostics[["parameter"]] %in% structural_parameters &
    diagnostics[["parameter"]] %in% selected_parameters
  structural_rows <- structural_rows |
    diagnostics[["parameter"]] %in%
      targets$metadata$parameter[targets$metadata$structural] &
    diagnostics[["parameter"]] %in% selected_parameters
  diagnostics[["state"]][structural_rows] <- "structural_constant"

  sample_rows <- selected_rows[
    diagnostics[["state"]][selected_rows] == "assessable"
  ]
  if(length(sample_rows) > 0L){
    chain_lengths <- vapply(mcmc_samples_list, nrow, integer(1))
    chain_id <- rep.int(seq_along(chain_lengths), chain_lengths)
    for(row in sample_rows){
      parameter <- diagnostics[["parameter"]][[row]]
      target_row <- match(parameter, targets$metadata$parameter)
      parameter_samples <- targets$samples[[target_row]]
      parameter_diagnostics <- .bt_convergence_parameter_diagnostics(
        parameter_samples,
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

.bt_convergence_sample_targets <- function(mcmc_samples, prior_list){

  columns <- colnames(mcmc_samples)
  if(is.null(columns)){
    columns <- character()
  }
  supports <- .bt_convergence_indicator_supports(prior_list)
  samples  <- list()
  metadata <- list()

  add_target <- function(parameter, source, values, is_indicator,
                         is_inclusion, structural = FALSE){
    samples[[length(samples) + 1L]] <<- values
    metadata[[length(metadata) + 1L]] <<- data.frame(
      parameter = parameter,
      source = source,
      is_indicator = is_indicator,
      is_inclusion = is_inclusion,
      structural = structural,
      stringsAsFactors = FALSE
    )
  }

  for(column in columns){
    values       <- mcmc_samples[, column]
    is_indicator <- grepl("_indicator(\\[[^]]+\\])?$", column)
    is_inclusion <- grepl("_inclusion(\\[[^]]+\\])?$", column)
    if(!is_indicator){
      add_target(column, column, values, FALSE, is_inclusion)
      next
    }

    if(any(!is.finite(values)) || any(values != round(values))){
      stop(
        "Model indicator '", column,
        "' must contain finite integer states.",
        call. = FALSE
      )
    }
    support <- supports[[column]]
    if(!is.null(support) && any(!values %in% support)){
      stop(
        "Model indicator '", column,
        "' contains states outside its prior support.",
        call. = FALSE
      )
    }
    if(!is.null(support) && length(support) == 1L){
      add_target(column, column, values, TRUE, FALSE, structural = TRUE)
      next
    }

    observed <- sort(unique(values))
    if(length(observed) <= 1L){
      add_target(column, column, values, TRUE, FALSE)
    }else if(length(observed) == 2L){
      occupancy <- as.numeric(values == observed[[2L]])
      add_target(column, column, occupancy, TRUE, FALSE)
    }else{
      for(state in observed){
        state_label <- format(
          state,
          digits = 17,
          scientific = FALSE,
          trim = TRUE
        )
        parameter <- paste0(column, " (state ", state_label, ")")
        occupancy <- as.numeric(values == state)
        add_target(parameter, column, occupancy, TRUE, FALSE)
      }
    }
  }

  if(length(metadata) == 0L){
    metadata <- data.frame(
      parameter = character(),
      source = character(),
      is_indicator = logical(),
      is_inclusion = logical(),
      structural = logical(),
      stringsAsFactors = FALSE
    )
  }else{
    metadata <- do.call(rbind, metadata)
    rownames(metadata) <- NULL
  }

  list(samples = samples, metadata = metadata)
}

.bt_convergence_indicator_supports <- function(prior_list){

  supports <- list()
  if(length(prior_list) == 0L){
    return(supports)
  }
  prior_names <- names(prior_list)
  if(is.null(prior_names)){
    prior_names <- rep.int("", length(prior_list))
  }

  for(i in seq_along(prior_list)){
    prior     <- prior_list[[i]]
    parameter <- prior_names[[i]]
    if(is.na(parameter) || !nzchar(parameter)){
      next
    }
    if(is.prior.spike_and_slab(prior)){
      inclusion <- .get_spike_and_slab_inclusion(prior)
      supports[[paste0(parameter, "_indicator")]] <-
        .bt_convergence_binary_indicator_support(inclusion)
    }else if(is.prior.mixture(prior)){
      indicator <- if(inherits(prior, "prior.bias_mixture")){
        "bias_indicator"
      }else{
        paste0(parameter, "_indicator")
      }
      supports[[indicator]] <- seq_along(prior)
    }else if(is.prior.factor(prior) && is.prior.ordered(prior)){
      metadata <- .prior_ordered_metadata(prior)
      if(is.prior.spike_and_slab(prior$total) && metadata$theta_dim > 1L){
        inclusion <- .get_spike_and_slab_inclusion(prior$total)
        total_name <- .prior_ordered_total_name(parameter)
        supports[[paste0(total_name, "_indicator")]] <-
          .bt_convergence_binary_indicator_support(inclusion)
      }
    }
  }

  supports
}

.bt_convergence_binary_indicator_support <- function(inclusion){

  if(is.prior.point(inclusion)){
    probability <- as.numeric(inclusion$parameters[["location"]])
    if(probability %in% c(0, 1)){
      return(probability)
    }
  }
  c(0, 1)
}

.bt_convergence_is_structural_prior <- function(prior){

  if(is.prior.point(prior)){
    return(TRUE)
  }
  if(is.prior.mixture(prior) && length(prior) == 1L){
    return(.bt_convergence_is_structural_prior(prior[[1L]]))
  }
  FALSE
}

.bt_convergence_structural_parameters <- function(prior_list){

  if(length(prior_list) == 0L){
    return(character())
  }
  structural <- names(prior_list)[vapply(
    prior_list,
    .bt_convergence_is_structural_prior,
    logical(1)
  )]
  unlist(lapply(structural, function(parameter){
    prior <- prior_list[[parameter]]
    names <- tryCatch(
      {
        if(is.prior.vector(prior) || is.prior.factor(prior)){
          .JAGS_prior_factor_names(parameter, prior)
        }else{
          parameter
        }
      },
      error = function(e) parameter
    )
    as.character(names)
  }), use.names = FALSE)
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
    Rhat = assess_Rhat,
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
