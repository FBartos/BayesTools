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

  check_int(chains, "chains", lower = check_mins[["chains"]], call = call)
  check_int(adapt,  "adapt",  lower = check_mins[["adapt"]],  call = call)
  check_int(burnin, "burnin", lower = check_mins[["burnin"]], call = call)
  check_int(sample, "sample", lower = check_mins[["sample"]], call = call)
  check_int(thin,   "thin",   lower = check_mins[["thin"]],   call = call)
  check_bool(parallel, "parallel",                call = call)
  check_int(cores,     "cores", lower = 1,        call = call)
  check_bool(autofit,  "autofit",                 call = call)
  check_bool(silent,   "silent",                  call = call)
  check_int(seed,      "seed", allow_NULL = TRUE, call = call)

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

  check_list(autofit_control, "autofit_control", check_names = c("max_Rhat", "min_ESS", "max_error", "max_SD_error",  "max_time", "sample_extend", "restarts", "max_extend", "check_indicators"), call = call)
  if(is.null(autofit_control[["check_indicators"]])){
    autofit_control[["check_indicators"]] <- FALSE
  }
  check_real(autofit_control[["max_Rhat"]],     "max_Rhat",     lower = 1, allow_NULL = TRUE, call = call)
  check_real(autofit_control[["min_ESS"]],      "min_ESS",      lower = 0, allow_NULL = TRUE, call = call)
  check_real(autofit_control[["max_error"]],    "max_error",    lower = 0, allow_NULL = TRUE, call = call)
  check_real(autofit_control[["max_SD_error"]], "max_SD_error", lower = 0, upper = 1, allow_NULL = TRUE, call = call)
  check_bool(autofit_control[["check_indicators"]], "check_indicators", call = call)
  check_list(autofit_control[["max_time"]],     "max_time", check_names = c("time", "unit"), check_length = 2, allow_NULL = TRUE, call = call)
  if(!is.null(autofit_control[["max_time"]])){
    if(is.null(names(autofit_control[["max_time"]]))){
      names(autofit_control[["max_time"]]) <- c("time", "unit")
    }
    check_real(autofit_control[["max_time"]][["time"]], "max_time:time", lower = 0, call = call)
    check_char(autofit_control[["max_time"]][["unit"]], "max_time:unit", allow_values = c("secs", "mins", "hours", "days", "weeks"), call = call)
  }
  check_int(autofit_control[["sample_extend"]], "sample_extend", lower = 1, allow_NULL = skip_sample_extend, call = call)
  check_int(autofit_control[["restarts"]], "restarts", lower = 1, allow_NULL = TRUE, call = call)
  check_int(autofit_control[["max_extend"]], "max_extend", lower = 1, allow_NULL = TRUE, call = call)

  return(invisible(autofit_control))
}
