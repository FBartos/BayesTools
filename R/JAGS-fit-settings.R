#' @title Check and list 'JAGS' fitting settings
#'
#' @description Checks and lists settings for the
#' [JAGS_fit] function. \code{adapt}, \code{burnin} and \code{sample} count
#' draws at the thinning \code{thin}; settings whose last retained draw,
#' iteration \code{(adapt + burnin + sample - 1) * thin + 1}, would be
#' numbered beyond \code{.Machine$integer.max} are refused. The settings are
#' described as the arguments of [JAGS_fit()]; the defaults named there are
#' those of [JAGS_fit()], as these functions require every setting except
#' \code{check_mins}, \code{skip_sample_extend} and \code{call}.
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
#' @param adapt length of the adaptation phase in draws per chain at the
#' thinning \code{thin} (\code{adapt * thin} iterations where a sampler
#' adapts), defaults to \code{500}
#' @param skip_sample_extend whether \code{sample_extend}
#' is allowed to be NULL and skipped in the check
#'
#' @inheritParams JAGS_fit
#' @inheritParams check_input
#'
#' @return \code{JAGS_check_and_list_fit_settings} invisibly returns a
#' list of checked fit settings. \code{JAGS_check_and_list_autofit_settings}
#' invisibly returns a list of checked autofit settings, with numeric
#' \code{burnin_extend} and \code{sample_extend} as the whole numbers that
#' are run.
#'
#' @export JAGS_check_and_list_fit_settings
#' @export JAGS_check_and_list_autofit_settings
#' @name JAGS_check_and_list
NULL

#' @rdname JAGS_check_and_list
JAGS_check_and_list_fit_settings     <- function(chains, adapt, burnin, sample, thin, autofit, parallel, cores, silent, seed, check_mins = list(chains = 1, adapt = 50, burnin = 50, sample = 100, thin = 1), call = ""){

  check_int(chains, "chains", lower = check_mins[["chains"]], allow_NA = FALSE, call = call)
  .JAGS_check_adapt(adapt, lower = check_mins[["adapt"]], call = call)
  check_int(burnin, "burnin", lower = check_mins[["burnin"]], allow_NA = FALSE, call = call)
  check_int(sample, "sample", lower = check_mins[["sample"]], allow_NA = FALSE, call = call)
  check_int(thin,   "thin",   lower = check_mins[["thin"]],   allow_NA = FALSE, call = call)
  # 'adapt', 'burnin' and 'sample' count draws at the fit's thinning; the
  # iteration numbers of the draws are integers
  last_iteration <- .JAGS_last_iteration(round(adapt), round(burnin), round(sample), round(thin))
  if(last_iteration > .Machine$integer.max){
    stop(paste0(call, "The settings 'adapt', 'burnin', 'sample' and 'thin' would number the last retained draw ",
      format(last_iteration, scientific = FALSE), ", beyond the largest iteration number (", .Machine$integer.max,
      "): a fit runs ('adapt' + 'burnin' + 'sample') times 'thin' iterations. Reduce 'adapt', 'burnin', ",
      "'sample' or 'thin'."), call. = FALSE)
  }
  check_bool(parallel, "parallel",                allow_NA = FALSE, call = call)
  check_int(cores,     "cores", lower = 1,        allow_NA = FALSE, call = call)
  check_bool(autofit,  "autofit",                 allow_NA = FALSE, call = call)
  check_bool(silent,   "silent",                  allow_NA = FALSE, call = call)
  check_int(seed,      "seed", lower = -.Machine$integer.max, upper = .Machine$integer.max, allow_NULL = TRUE, allow_NA = FALSE, call = call)

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

# 'adapt' of JAGS_fit() and of an extension, in draws at the fit's thinning.
# rjags counts the 'adapt * thin' adaptation iterations as an integer, so
# larger values of 'adapt' alone cannot be run; the product is checked with
# the other settings of a fit (.JAGS_last_iteration()) and with the fit's
# 'thin' in JAGS_extend().
.JAGS_check_adapt <- function(adapt, lower = 50, call = ""){

  check_int(adapt, "adapt", lower = lower, upper = .Machine$integer.max, allow_NA = FALSE, call = call)
}

# The iteration number of the last retained draw of a fit: 'adapt', 'burnin'
# and 'sample' count draws at the thinning 'thin', so the fit runs
# 'adapt * thin' adaptation, 'burnin * thin' burn-in and 'sample * thin'
# sampling iterations, and its first draw is iteration
# '(adapt + burnin) * thin + 1'. In double arithmetic, which these whole
# numbers keep exact far beyond the integer range.
.JAGS_last_iteration <- function(adapt, burnin, sample, thin){

  (as.double(adapt) + burnin) * thin + 1 + (as.double(sample) - 1) * thin
}

# The whole number that a count accepted by check_int() stands for. The
# validator accepts values within a tolerance of a whole number, which rjags
# would truncate: a 'thin' of 2 - 1e-9 would thin by 1 and keep twice the
# draws. Integers and whole doubles are returned as given.
.JAGS_whole_count <- function(x){

  if(is.integer(x)) x else round(x)
}

#' @rdname JAGS_check_and_list
JAGS_check_and_list_autofit_settings <- function(autofit_control, skip_sample_extend = FALSE, call = ""){

  check_list(
    autofit_control,
    "autofit_control",
    check_names = c(
      "max_Rhat", "min_ESS", "max_error", "max_SD_error", "max_time",
      "sample_extend", "restarts", "max_extend", "check_indicators",
      "monitor",
      "burnin_max_Rhat", "burnin_extend", "burnin_max_extend"
    ),
    call = call
  )
  if(is.null(autofit_control[["check_indicators"]])){
    autofit_control[["check_indicators"]] <- FALSE
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
  # a positive whole number of draws, or "auto" (.bt_autofit_extension_size())
  if(is.character(autofit_control[["sample_extend"]])){
    check_char(autofit_control[["sample_extend"]], "sample_extend", allow_values = "auto", allow_NA = FALSE, call = call)
    # the plain string, without names, which the fitting loops compare with
    # identical()
    autofit_control[["sample_extend"]] <- "auto"
  }else{
    check_int(autofit_control[["sample_extend"]], "sample_extend", lower = 1, allow_NULL = skip_sample_extend, allow_NA = FALSE, call = call)
    # the whole number of draws that the validator accepted (.JAGS_whole_count())
    if(!is.null(autofit_control[["sample_extend"]])){
      autofit_control[["sample_extend"]] <- .JAGS_whole_count(autofit_control[["sample_extend"]])
    }
  }
  check_int(autofit_control[["restarts"]], "restarts", lower = 1, allow_NULL = TRUE, allow_NA = FALSE, call = call)
  check_int(autofit_control[["max_extend"]], "max_extend", lower = 1, allow_NULL = TRUE, allow_NA = FALSE, call = call)

  # The readiness check between burn-in and sampling: an omitted element takes
  # its default (.bt_burnin_check_settings()) and an explicit
  # 'burnin_max_Rhat = NULL' turns the check off, so these elements are
  # validated where present and never filled in.
  supplied <- function(name) name %in% names(autofit_control)
  # a finite number of at least 1, or "auto" (.bt_burnin_threshold()); NULL
  # turns the check off, an empty value is refused
  if(!is.null(autofit_control[["burnin_max_Rhat"]]) && length(autofit_control[["burnin_max_Rhat"]]) == 0L){
    stop(paste0(call, "The 'burnin_max_Rhat' argument must have length '1'."), call. = FALSE)
  }
  if(is.character(autofit_control[["burnin_max_Rhat"]])){
    check_char(autofit_control[["burnin_max_Rhat"]], "burnin_max_Rhat", allow_values = "auto", allow_NA = FALSE, call = call)
    autofit_control[["burnin_max_Rhat"]] <- "auto"
  }else{
    check_real(autofit_control[["burnin_max_Rhat"]], "burnin_max_Rhat", lower = 1, allow_NULL = TRUE, allow_NA = FALSE, call = call)
    if(!is.null(autofit_control[["burnin_max_Rhat"]]) && !is.finite(autofit_control[["burnin_max_Rhat"]])){
      stop(
        paste0(call, "The 'burnin_max_Rhat' argument must contain only finite values."),
        call. = FALSE
      )
    }
  }
  check_int(autofit_control[["burnin_extend"]], "burnin_extend", lower = 100, upper = .Machine$integer.max,
    allow_NULL = !supplied("burnin_extend"), allow_NA = FALSE, call = call)
  if(!is.null(autofit_control[["burnin_extend"]])){
    autofit_control[["burnin_extend"]] <- .JAGS_whole_count(autofit_control[["burnin_extend"]])
  }
  check_int(autofit_control[["burnin_max_extend"]], "burnin_max_extend", lower = 0,
    allow_NULL = !supplied("burnin_max_extend"), allow_NA = FALSE, call = call)

  return(invisible(autofit_control))
}

# The settings of the readiness check between burn-in and sampling. Unlike the
# convergence criteria, an omitted element takes its default; an explicit
# 'burnin_max_Rhat = NULL' turns the check off.
.bt_burnin_check_settings <- function(autofit_control){

  settings <- list(burnin_max_Rhat = "auto", burnin_extend = 500, burnin_max_extend = 10)
  for(name in names(settings)){
    if(name %in% names(autofit_control)){
      settings[name] <- list(autofit_control[[name]])
    }
  }
  settings
}

# The threshold of the readiness check for a fit with 'chains' chains: "auto"
# falls with the number of chains; a number, or NULL (check off), stays.
# "auto" is 1.7 + 0.9 / chains, computed as one division so that the result is
# the double nearest the exact value (1.85 for six chains; the sum gives the
# double below it).
.bt_burnin_threshold <- function(burnin_max_Rhat, chains){

  if(identical(burnin_max_Rhat, "auto")) (17 * chains + 9) / (10 * chains) else burnin_max_Rhat
}
