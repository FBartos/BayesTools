# Run the chains assigned to one R process. The caller owns process-level
# setup and reporting. Keep JAGS to one thread per process for both serial
# fits and PSOCK workers, including workers that own multiple chains.
.bt_rjags_group <- function(inits, model, data, monitor, adapt, burnin, sample,
                            thin, silent, keep_live, live_model = NULL){

  progress <- if(silent) "none" else "text"
  if(is.null(live_model)){
    uses_lecuyer <- vapply(inits, function(chain){
      identical(chain[[".RNG.name"]], "lecuyer::RngStream")
    }, logical(1))
    if(any(uses_lecuyer)){
      rjags::load.module("lecuyer", quiet = TRUE)
    }
    connection <- textConnection(model)
    on.exit(close(connection))
    live_model <- rjags::jags.model(connection, data = data, inits = inits,
      n.chains = length(inits), n.threads = 1L, n.adapt = 0L, quiet = silent)
    adapted <- rjags::adapt(live_model, n.iter = adapt, end.adaptation = TRUE,
      progress.bar = progress)
    if(!adapted){
      warning(warningCondition(
        paste0("The JAGS adaptation phase was incomplete after ", adapt, " iterations."),
        class = "BayesTools_adaptation_incomplete"
      ))
    }
  }
  if(burnin > 0){
    stats::update(live_model, n.iter = burnin, progress.bar = progress)
  }
  draws <- .bt_rjags_samples(live_model, monitor, sample, thin, progress)
  assignments <- rjags::list.samplers(live_model)
  samplers <- data.frame(
    Index.No = rep(seq_along(assignments), lengths(assignments)),
    Sampler = rep(as.character(names(assignments)), lengths(assignments)),
    Node = as.character(unlist(assignments, use.names = FALSE))
  )

  list(mcmc = draws, end_state = live_model$state(internal = TRUE), samplers = samplers,
    live_model = if(keep_live) live_model)
}

# rjags 5-5 treats names containing an internal "_word_" like its reserved
# monitors and passes their subscript to JAGS as part of the variable name.
# Read those arrays by name, retaining just the coordinates the user requested.
# Other indexed requests can still read a small slice of a large array.
.bt_rjags_samples <- function(model, monitor, sample, thin, progress){

  monitor <- unname(trimws(monitor))
  base <- trimws(.bt_jags_parameter_base(monitor))
  indexed <- base != monitor
  whole <- unique(c(base[!indexed], base[indexed & grepl("_[A-Za-z][A-Za-z0-9]*_", base)]))
  native <- ifelse(base %in% whole, base, monitor)
  draws <- rjags::coda.samples(model, variable.names = unique(unname(native)),
    n.iter = as.double(sample) * thin, thin = thin, na.rm = FALSE, progress.bar = progress)

  # JAGS returns variables in name order, and each array in column-major order.
  # Retain the original request order across variables and overlapping slices.
  columns <- colnames(draws[[1L]])
  by_name <- split(seq_along(columns), .bt_jags_parameter_base(columns))
  selection <- lapply(seq_along(monitor), function(i){
    available <- by_name[[base[i]]]
    if(!length(available)){
      stop("JAGS returned no draws for requested monitor '", monitor[i], "'.", call. = FALSE)
    }
    if(!indexed[i]) return(stats::setNames(available, columns[available]))

    spec <- .bt_fit_mcarray_parse_varname(monitor[i])
    if(is.null(spec$lower) || any(!is.finite(c(spec$lower, spec$upper))) ||
       any(spec$lower < 1 | spec$lower != floor(spec$lower) | spec$upper != floor(spec$upper))){
      stop("Invalid JAGS monitor subset '", monitor[i], "'.", call. = FALSE)
    }
    size <- spec$upper - spec$lower + 1
    if(prod(size) == 1){
      requested <- paste0(base[i], "[", paste(spec$lower, collapse = ","), "]")
    }else{
      if(prod(size) > length(available)){
        stop("JAGS did not return the full requested subset '", monitor[i], "'.", call. = FALSE)
      }
      requested <- rjags::expand.name(monitor[i], size)
    }
    position <- match(requested, columns)
    # A full one-element array is returned under its bare name. Keep its
    # explicit index when the caller selected the element by subscript.
    if(length(available) == 1L && identical(columns[available], base[i]) &&
       all(spec$lower == 1 & spec$upper == 1)) position <- available
    if(anyNA(position)){
      stop("JAGS did not return the full requested subset '", monitor[i], "'.", call. = FALSE)
    }
    stats::setNames(position, requested)
  })
  # Validate ranges against the declared array above, before applying coda's
  # omission of cells undefined in any chain (e.g. x[1] and x[3] are defined
  # while x[2] is not). A hole does not make x[1:3] an out-of-bounds request.
  if(any(vapply(draws, anyNA, logical(1)))){
    defined <- Reduce(`&`, lapply(draws, function(chain) colSums(is.na(chain)) == 0))
    selection <- lapply(selection, function(index) index[defined[index]])
  }
  selection <- unlist(selection, use.names = TRUE)
  selection <- selection[!duplicated(selection)]
  if(identical(names(selection), columns)) return(draws)
  coda::mcmc.list(lapply(draws, function(chain){
    chain <- chain[, unname(selection), drop = FALSE]
    colnames(chain) <- names(selection)
    chain
  }))
}

# Warnings emitted on a worker must reach the calling fit, including the
# incomplete-adaptation warning. Do not send compiled JAGS objects over sockets.
.bt_rjags_worker <- function(inits, ...){

  warnings <- list()
  value <- tryCatch(withCallingHandlers(
    .bt_rjags_group(inits, ..., keep_live = FALSE),
    warning = function(w){
      warnings[[length(warnings) + 1L]] <<- w
      invokeRestart("muffleWarning")
    }), error = identity)
  list(value = value, warnings = warnings)
}

.bt_rjags_run <- function(model, data, inits, monitor, adapt, burnin, sample,
                          thin, cl, silent, live_model = NULL){

  if(is.null(cl)){
    return(.bt_rjags_group(inits, model, data, monitor, adapt, burnin, sample,
      thin, silent, keep_live = TRUE, live_model = live_model))
  }

  # The runtime callbacks use the same round-robin chain assignment. Restore
  # the original chain order before returning draws and their matching states.
  groups <- split(seq_along(inits), rep_len(seq_along(cl), length(inits)))
  results <- parallel::clusterApply(cl, lapply(groups, function(ids) inits[ids]),
    .bt_rjags_worker, model = model, data = data, monitor = monitor, adapt = adapt,
    burnin = burnin, sample = sample, thin = thin, silent = silent)
  for(result in results){
    for(w in result$warnings) warning(w)
    if(inherits(result$value, "error")) stop(result$value)
  }
  out <- lapply(results, "[[", "value")
  chain_order <- order(unlist(groups, use.names = FALSE))
  mcmc <- unlist(lapply(out, function(group) unclass(group$mcmc)),
    recursive = FALSE, use.names = FALSE)
  end_state <- unlist(lapply(out, "[[", "end_state"), recursive = FALSE, use.names = FALSE)
  list(mcmc = coda::mcmc.list(mcmc[chain_order]), end_state = end_state[chain_order],
    samplers = out[[1L]]$samplers, live_model = NULL)
}

.bt_rjags_fit <- function(model, data, inits, monitor, adapt, burnin, sample, thin,
                          cl = NULL, silent = TRUE){

  out <- .bt_rjags_run(model, data, inits, monitor, adapt, burnin, sample, thin, cl, silent)
  # Public chain indices count the requested adaptation and burn-in, even
  # when JAGS selects only samplers that never require adaptation.
  out$mcmc <- coda::mcmc.list(lapply(out$mcmc, function(chain){
    coda::mcmc(as.matrix(chain), start = as.double(adapt) + burnin + 1, thin = thin)
  }))
  out$model <- model
  out["data"] <- list(data)
  out$monitor <- monitor
  out$sample <- sample
  out$thin <- thin
  out$burnin <- burnin
  out$backend <- list(jags = as.character(rjags::jags.version()),
    rjags = as.character(utils::packageVersion("rjags")))
  class(out) <- "BayesTools_fit"
  out
}

# Realize a BUGS data block once so the fit records exactly which elements
# are observed, including partly generated arrays. Re-running a stochastic
# data block on different workers or extensions would change the posterior.
.bt_rjags_realize_data <- function(model, data, seed, monitor){

  if(!grepl("data", model, fixed = TRUE)) return(list(model = model, data = data))
  masked <- .JAGS_mask_comments(model)
  data_at <- regexpr("(?<![A-Za-z0-9_.])data[[:space:]]*\\{", masked, perl = TRUE)[1L]
  if(data_at < 1L) return(list(model = model, data = data))
  model_at <- regexpr("(?<![A-Za-z0-9_.])model[[:space:]]*\\{", masked, perl = TRUE)[1L]
  if(model_at < data_at){
    stop("The BUGS data block must precede the model block.", call. = FALSE)
  }
  prefix <- substr(model, 1L, model_at - 1L)
  inputs <- .bt_sampler_tokens(prefix)
  symbols <- .bt_sampler_tokens(model)
  anchor <- utils::tail(make.unique(c(symbols, "BayesTools_data_anchor")), 1L)
  connection <- textConnection(paste0(prefix, "\nmodel { ", anchor, " <- 0 }"))
  on.exit(close(connection))
  if(grepl("~", substr(masked, data_at, model_at - 1L), fixed = TRUE)){
    rjags::jags.seed(seed)
  }
  generator <- rjags::jags.model(connection, data = data[names(data) %in% inputs],
    n.chains = 1L, n.threads = 1L, n.adapt = 0L, quiet = TRUE)
  generated <- generator$data()
  data[names(generated)] <- generated
  declarations <- substr(model, 1L, data_at - 1L)
  sampling_block <- substring(model, model_at)
  sampling_model <- paste0(declarations, sampling_block)
  # An observed array may be requested as a monitor even if it was referenced
  # solely in the original data graph. Keep its name in the rebuilt program:
  # rjags discards input data whose name does not occur in the model text.
  monitored_data <- intersect(.bt_jags_parameter_base(monitor), names(generated))
  extra <- setdiff(monitored_data, .bt_sampler_tokens(sampling_model))
  if(length(extra)){
    references <- paste0(anchor, "[", seq_along(extra), "] <- length(", extra, ")")
    sampling_model <- paste(declarations, "data {", paste(references, collapse = "\n"),
      "}", sampling_block, sep = "\n")
  }
  list(model = sampling_model, data = data[names(data) %in% .bt_sampler_tokens(sampling_model)])
}

.bt_rjags_extend <- function(fit, sample, adapt, cl = NULL, silent = TRUE){

  out <- .bt_rjags_run(fit$model, fit$data, fit$end_state, fit$monitor, adapt,
    0L, sample, fit$thin, cl, silent, live_model = fit$live_model)
  fit$mcmc <- coda::mcmc.list(Map(function(previous, additional){
    if(!identical(colnames(previous), colnames(additional))){
      stop("The rebuilt model returned different monitored coordinates.", call. = FALSE)
    }
    coda::mcmc(rbind(as.matrix(previous), as.matrix(additional)),
      start = stats::start(previous), thin = fit$thin)
  }, fit$mcmc, out$mcmc))
  fit$sample <- as.double(fit$sample) + sample
  fit$end_state <- out$end_state
  fit$samplers <- out$samplers
  fit$live_model <- out$live_model
  fit
}
