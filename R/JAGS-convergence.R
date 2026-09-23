#' @title Assess convergence of a runjags model
#'
#' @description Checks whether the supplied \link[runjags]{runjags-package} model
#' satisfied convergence criteria.
#' @param fit a runjags model
#' @param prior_list named list of prior distribution
#' (names correspond to the parameter names)
#' @param max_Rhat maximum R-hat error for the autofit function.
#'   Defaults to \code{1.05}. With one chain, this criterion is skipped with
#'   a warning; the remaining enabled criteria are still assessed.
#' @param min_ESS minimum effective sample size. Defaults to \code{500}.
#' @param max_error maximum MCMC error. Defaults to \code{0.01}.
#' @param max_SD_error maximum MCMC error as the proportion of standard
#'   deviation of the parameters. Defaults to \code{0.05}.
#' @param add_parameters vector of additional parameter names that are excluded
#' from the default selection (only allows removing last, fixed, omega element
#' if omega is tracked manually). Parameters named in \code{monitor} are
#' checked even when they are listed here.
#' @param fail_fast whether the function should stop after the first failed convergence check.
#' @param check_indicators whether model indicator variables should be included
#' in convergence checks. Binary indicators are checked as Bernoulli
#' occupancies and categorical indicators are checked separately for every
#' observed state. When \code{monitor} is supplied, eligible indicators are
#' added to that selection. Auxiliary inclusion-probability coordinates remain
#' excluded unless named explicitly in \code{monitor}. Defaults to \code{FALSE}.
#' @param monitor optional character vector selecting parameters for convergence
#' checks. A base name selects all of its indexed elements. Requests are
#' resolved against all monitored columns, including \code{add_parameters}.
#' \code{NULL} selects every eligible parameter; \code{character()} requests no
#' parameters.
#' @param allow_not_assessable whether requested sampled parameters with
#' undefined diagnostics may be ignored. Defaults to \code{FALSE}. A sampled
#' column that never changes, including a constant model indicator, is not
#' evidence of convergence and remains not assessable. Only constants declared
#' by the prior are structural, such as point priors, reference and fixed
#' publication-weight bins, a p-hacking kind shared by every mixture branch,
#' and the point total of an ordered prior.
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
#' \code{"not_requested"}, or \code{"not_checked"} when \code{fail_fast = TRUE}
#' stops before reaching that parameter. The \code{errors} attribute carries
#' failed checks.
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

  prepared <- .bt_convergence_prepare(
    fit = fit,
    prior_list = prior_list,
    add_parameters = add_parameters,
    monitor = monitor
  )
  mcmc_samples_list     <- prepared$mcmc_samples_list
  targets               <- prepared$targets
  structural_parameters <- prepared$structural_parameters
  available_parameters  <- prepared$available_parameters
  available_sources     <- prepared$available_sources

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
    if(check_indicators){
      selected_parameters <- unique(c(
        selected_parameters,
        targets$metadata$parameter[
          targets$metadata$is_indicator & !targets$metadata$is_inclusion
        ]
      ))
    }
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

.bt_convergence_prepare <- function(fit, prior_list, add_parameters, monitor){

  # extract samples and parameter information
  mcmc_samples_list <- .extract_posterior_samples(fit, as_list = TRUE)
  mcmc_samples      <- do.call(rbind, mcmc_samples_list)

  # Remove parameters that are intentionally excluded from automatic checks.
  # Structural point parameters are added back below from prior metadata.
  # Additional monitors stay excluded unless 'monitor' requests them.
  remove_params <- c(
    names(prior_list)[vapply(
      prior_list,
      .bt_convergence_is_structural_prior,
      logical(1)
    )],
    .bt_convergence_excluded_add_parameters(add_parameters, monitor)
  )

  cleaned <- .remove_auxiliary_parameters(mcmc_samples, prior_list, remove_params)
  mcmc_samples <- cleaned$model_samples

  targets <- .bt_convergence_sample_targets(mcmc_samples, prior_list)
  sample_parameters <- targets$metadata$parameter
  structural_parameters <- .bt_convergence_structural_parameters(prior_list)

  list(
    mcmc_samples_list = mcmc_samples_list,
    targets = targets,
    structural_parameters = structural_parameters,
    available_parameters = unique(c(sample_parameters, structural_parameters)),
    available_sources = c(
      targets$metadata$source,
      structural_parameters[!structural_parameters %in% sample_parameters]
    )
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
.bt_convergence_validate_monitor <- function(fit, prior_list, add_parameters,
                                             monitor){

  if(length(monitor) == 0L){
    return(invisible(TRUE))
  }

  prepared <- .bt_convergence_prepare(
    fit = fit,
    prior_list = prior_list,
    add_parameters = add_parameters,
    monitor = monitor
  )
  .bt_convergence_resolve_monitor(
    monitor,
    prepared$available_parameters,
    prepared$available_sources
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

.bt_convergence_sample_targets <- function(mcmc_samples, prior_list){

  columns <- colnames(mcmc_samples)
  if(is.null(columns)){
    columns <- character()
  }
  supports <- .bt_convergence_indicator_supports(prior_list)
  structural_columns <- c(
    .bt_convergence_structural_omega_bins(prior_list),
    .bt_convergence_structural_phacking_columns(prior_list),
    .bt_convergence_structural_ordered_columns(prior_list)
  )
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
      add_target(
        column,
        column,
        values,
        FALSE,
        is_inclusion,
        structural = column %in% structural_columns
      )
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

.bt_convergence_structural_omega_bins <- function(prior_list){

  if(length(prior_list) == 0L){
    return(character())
  }

  structural <- character()
  for(prior in prior_list){
    if(is.prior.weightfunction(prior)){
      cuts <- weightfunctions_mapping(list(prior), cuts_only = TRUE)
      names <- if(length(cuts) >= 2L){
        paste0("omega[", cuts[-length(cuts)], ",", cuts[-1], "]")
      }else{
        character()
      }
      structural <- c(structural, "omega[1]")
      if(identical(prior$weights$type, "fixed")){
        structural <- c(
          structural,
          names,
          paste0("omega[", seq_len(max(length(cuts) - 1L, 0L)), "]")
        )
      }else if(length(names) > 0L){
        structural <- c(structural, names[[1L]])
      }
    }else if(is_prior_bias(prior) || is_prior_phacking(prior) ||
             inherits(prior, "prior.bias_mixture")){
      structural <- c(
        structural,
        .bt_convergence_structural_selection_bins(prior)
      )
    }
  }

  unique(structural)
}

# Composed bias priors, p-hacking priors, and publication-bias mixtures share
# one omega vector on the global one-sided cut grid of the selection backend.
# A global bin is constant by construction when every mixture branch fixes it
# to the same value: branches without a step selection contribute 1, the
# reference bin (local bin 1) of a selection is 1, and every bin of fixed
# weights is its declared weight. Mirrored two-sided bins map to their local
# bin through the component expansion. Single composed priors with a selection
# are summarized under their renamed global bins; other priors keep the raw
# monitored omega coordinates.
.bt_convergence_structural_selection_bins <- function(prior){

  branches <- .selection_normalize_priors(prior)
  branch_info <- lapply(branches, .selection_branch_info)
  has_selection <- vapply(branch_info, function(x) !is.null(x$selection), logical(1))
  has_phacking  <- vapply(branch_info, function(x) !is.null(x$phacking), logical(1))
  if(!any(has_selection) && !any(has_phacking)){
    return(character())
  }

  cuts <- if(any(has_selection)){
    weightfunctions_mapping(
      lapply(branch_info[has_selection], function(x) x$selection),
      cuts_only = TRUE,
      one_sided = TRUE
    )
  }else{
    c(0, 1)
  }
  n_bins <- length(cuts) - 1L

  branch_values <- vapply(branch_info, function(x){
    .bt_convergence_selection_bin_values(x$selection, cuts)
  }, numeric(n_bins))
  branch_values <- matrix(branch_values, nrow = n_bins)
  constant <- apply(branch_values, 1L, function(values){
    all(!is.na(values)) && all(values == values[[1L]])
  })

  bin_names <- if(n_bins == 1L){
    "omega"
  }else if(is_prior_bias(prior) && any(has_selection)){
    .weightfunction_omega_names(cuts)
  }else{
    paste0("omega[", seq_len(n_bins), "]")
  }

  bin_names[constant]
}

.bt_convergence_selection_bin_values <- function(selection, cuts){

  n_bins <- length(cuts) - 1L
  if(is.null(selection)){
    return(rep(1, n_bins))
  }

  expansion <- .weightfunction_mapping_expansion(selection, force_one_sided = TRUE)
  local_bins <- expansion$index[.weightfunction_global_bin_indices(cuts, expansion)]
  if(identical(selection$weights$type, "fixed")){
    return(as.numeric(selection$weights$omega[local_bins]))
  }

  ifelse(local_bins == 1L, 1, NA_real_)
}

# The monitored p-hacking kind is a declared constant of each mixture branch:
# the form code of a p-hacking branch and 0 for branches without p-hacking.
# It is structural when every branch declares the same code, and remains
# assessable when the kind varies with the mixture indicator.
.bt_convergence_structural_phacking_columns <- function(prior_list){

  if(length(prior_list) == 0L){
    return(character())
  }

  structural <- character()
  for(prior in prior_list){
    if(!(is_prior_bias(prior) || is_prior_phacking(prior) ||
         inherits(prior, "prior.bias_mixture"))){
      next
    }
    branch_info <- lapply(.selection_normalize_priors(prior), .selection_branch_info)
    has_phacking <- vapply(branch_info, function(x) !is.null(x$phacking), logical(1))
    if(!any(has_phacking)){
      next
    }
    kinds <- vapply(branch_info, function(x){
      if(is.null(x$phacking)) 0 else as.numeric(.phack_kind(x$phacking$form))
    }, numeric(1))
    if(all(kinds == kinds[[1L]])){
      structural <- c(structural, "phack_kind")
    }
  }

  unique(structural)
}

# Ordered priors with a point total monitor a constant total. Their level
# coefficients are constant as well when the total is zero or when every
# allocation is fixed.
.bt_convergence_structural_ordered_columns <- function(prior_list){

  prior_names <- names(prior_list)
  if(length(prior_list) == 0L || is.null(prior_names)){
    return(character())
  }

  structural <- character()
  for(i in seq_along(prior_list)){
    prior     <- prior_list[[i]]
    parameter <- prior_names[[i]]
    if(is.na(parameter) || !nzchar(parameter) || !is.prior.ordered(prior)){
      next
    }
    total_value <- .bt_convergence_point_value(prior$total)
    metadata <- attr(prior, "ordered_metadata", exact = TRUE)
    if(is.null(total_value) || is.null(metadata)){
      next
    }

    total_name <- .prior_ordered_total_name(parameter)
    structural <- c(
      structural,
      if(metadata$theta_dim == 1L){
        total_name
      }else{
        paste0(total_name, "[", seq_len(metadata$theta_dim), "]")
      }
    )

    fixed_allocation <- all(vapply(metadata$allocations, function(record){
      identical(record$spec$type, "fixed")
    }, logical(1)))
    if(total_value == 0 || fixed_allocation){
      structural <- c(
        structural,
        if(metadata$coefficient_dim == 1L){
          parameter
        }else{
          paste0(parameter, "[", seq_len(metadata$coefficient_dim), "]")
        }
      )
    }
  }

  structural
}

.bt_convergence_point_value <- function(prior){

  if(is.prior.mixture(prior) && length(prior) == 1L){
    return(.bt_convergence_point_value(prior[[1L]]))
  }
  if(!is.prior.point(prior) || .is_prior_expression(prior)){
    return(NULL)
  }

  location <- prior$parameters[["location"]]
  if(!is.numeric(location) || length(location) != 1L || !is.finite(location)){
    return(NULL)
  }

  location
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
