#' @title Check and list 'JAGS' fitting settings
#'
#' @description Checks and lists settings for the
#' [JAGS_fit] function.
#'
#' @param check_mins named list of minimal values for which
#' should some input be checked. Defaults to:
#' \describe{
#'   \item{chains}{\code{1}}
#'   \item{adapt}{\code{50}}
#'   \item{burnin}{\code{50}}
#'   \item{sample}{\code{100}}
#'   \item{thin}{\code{1}}
#' }
#' @param skip_sample_extend whether \code{sample_extend}
#' is allowed to be NULL and skipped in the check
#'
#' @inheritParams JAGS_fit
#' @inheritParams check_input
#'
#' @return \code{JAGS_check_and_list_fit_settings} invisibly returns a
#' list of checked fit settings. \code{JAGS_check_and_list_autofit_settings}
#' invisibly returns a list of checked autofit settings.
#' parameter names.
#'
#' @export JAGS_check_and_list_fit_settings
#' @export JAGS_check_and_list_autofit_settings
#' @name JAGS_check_and_list
NULL

#' @rdname JAGS_check_and_list
JAGS_check_and_list_fit_settings     <- function(chains, adapt, burnin, sample, thin, autofit, parallel, cores, silent, seed, check_mins = list(chains = 1, adapt = 50, burnin = 50, sample = 100, thin = 1), call = ""){

  check_int(chains, "chains", lower = check_mins[["chains"]], allow_NA = FALSE, call = call)
  check_int(adapt,  "adapt",  lower = check_mins[["adapt"]],  allow_NA = FALSE, call = call)
  check_int(burnin, "burnin", lower = check_mins[["burnin"]], allow_NA = FALSE, call = call)
  check_int(sample, "sample", lower = check_mins[["sample"]], allow_NA = FALSE, call = call)
  check_int(thin,   "thin",   lower = check_mins[["thin"]],   allow_NA = FALSE, call = call)
  check_bool(parallel, "parallel",                allow_NA = FALSE, call = call)
  check_int(cores,     "cores", lower = 1,        allow_NA = FALSE, call = call)
  check_bool(autofit,  "autofit",                 allow_NA = FALSE, call = call)
  check_bool(silent,   "silent",                  allow_NA = FALSE, call = call)
  check_int(seed,      "seed", allow_NULL = TRUE, allow_NA = FALSE, call = call)

  return(invisible(list(
    chains   = chains,
    adapt    = adapt,
    burnin   = burnin,
    sample   = sample,
    thin     = thin,
    autofit  = autofit,
    parallel = parallel,
    cores    = cores,
    silent   = silent,
    seed     = seed
  )))
}

#' @rdname JAGS_check_and_list
JAGS_check_and_list_autofit_settings <- function(autofit_control, skip_sample_extend = FALSE, call = ""){

  check_list(
    autofit_control,
    "autofit_control",
    check_names = c(
      "max_Rhat", "min_ESS", "max_error", "max_SD_error", "max_time",
      "sample_extend", "restarts", "max_extend", "check_indicators",
      "monitor", "allow_not_assessable"
    ),
    call = call
  )
  if(is.null(autofit_control[["check_indicators"]])){
    autofit_control[["check_indicators"]] <- FALSE
  }
  if(is.null(autofit_control[["allow_not_assessable"]])){
    autofit_control[["allow_not_assessable"]] <- FALSE
  }
  check_real(autofit_control[["max_Rhat"]],     "max_Rhat",     lower = 1, allow_NULL = TRUE, allow_NA = FALSE, call = call)
  check_real(autofit_control[["min_ESS"]],      "min_ESS",      lower = 0, allow_NULL = TRUE, allow_NA = FALSE, call = call)
  check_real(autofit_control[["max_error"]],    "max_error",    lower = 0, allow_NULL = TRUE, allow_NA = FALSE, call = call)
  check_real(autofit_control[["max_SD_error"]], "max_SD_error", lower = 0, upper = 1, allow_NULL = TRUE, allow_NA = FALSE, call = call)
  for(name in c("max_Rhat", "min_ESS", "max_error", "max_SD_error")){
    value <- autofit_control[[name]]
    if(!is.null(value) && any(!is.finite(value))){
      stop(
        paste0(call, "The '", name, "' argument must contain only finite values."),
        call. = FALSE
      )
    }
  }
  check_bool(autofit_control[["check_indicators"]], "check_indicators", allow_NA = FALSE, call = call)
  check_char(
    autofit_control[["monitor"]],
    "monitor",
    check_length = 0,
    allow_NULL = TRUE,
    allow_NA = FALSE,
    call = call
  )
  if(!is.null(autofit_control[["monitor"]]) &&
     length(autofit_control[["monitor"]]) == 0L){
    stop(
      paste0(
        call,
        "The 'monitor' argument must select at least one parameter for ",
        "automatic fitting."
      ),
      call. = FALSE
    )
  }
  check_bool(
    autofit_control[["allow_not_assessable"]],
    "allow_not_assessable",
    allow_NA = FALSE,
    call = call
  )
  check_list(autofit_control[["max_time"]],     "max_time", check_names = c("time", "unit"), check_length = 2, allow_NULL = TRUE, call = call)
  if(!is.null(autofit_control[["max_time"]])){
    if(is.null(names(autofit_control[["max_time"]]))){
      names(autofit_control[["max_time"]]) <- c("time", "unit")
    }
    check_real(autofit_control[["max_time"]][["time"]], "max_time:time", lower = 0, allow_NA = FALSE, call = call)
    if(any(!is.finite(autofit_control[["max_time"]][["time"]]))){
      stop(
        paste0(call, "The 'max_time:time' argument must contain only finite values."),
        call. = FALSE
      )
    }
    check_char(autofit_control[["max_time"]][["unit"]], "max_time:unit", allow_values = c("secs", "mins", "hours", "days", "weeks"), allow_NA = FALSE, call = call)
  }
  check_int(autofit_control[["sample_extend"]], "sample_extend", lower = 1, allow_NULL = skip_sample_extend, allow_NA = FALSE, call = call)
  check_int(autofit_control[["restarts"]], "restarts", lower = 1, allow_NULL = TRUE, allow_NA = FALSE, call = call)
  check_int(autofit_control[["max_extend"]], "max_extend", lower = 1, allow_NULL = TRUE, allow_NA = FALSE, call = call)

  return(invisible(autofit_control))
}
