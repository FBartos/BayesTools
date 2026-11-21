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
  draws <- rjags::coda.samples(live_model, variable.names = monitor,
    n.iter = as.double(sample) * thin, thin = thin, progress.bar = progress)
  assignments <- rjags::list.samplers(live_model)
  samplers <- data.frame(
    Index.No = rep(seq_along(assignments), lengths(assignments)),
    Sampler = rep(as.character(names(assignments)), lengths(assignments)),
    Node = as.character(unlist(assignments, use.names = FALSE))
  )

  list(mcmc = draws, end_state = live_model$state(internal = TRUE), samplers = samplers,
    live_model = if(keep_live) live_model)
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

  prepared <- .bt_rjags_realize_data(model, data, inits[[1L]][[".RNG.seed"]])
  model <- prepared$model
  data <- prepared$data
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

# A stochastic BUGS data block makes one data set for the entire fit. Run
# that block once, otherwise workers and later extensions would each simulate
# their own observations and hence target different posteriors.
.bt_rjags_realize_data <- function(model, data, seed){

  if(!grepl("data", model, fixed = TRUE)) return(list(model = model, data = data))
  masked <- .JAGS_mask_comments(model)
  data_at <- regexpr("(?<![A-Za-z0-9_.])data[[:space:]]*\\{", masked, perl = TRUE)[1L]
  if(data_at < 1L) return(list(model = model, data = data))
  model_at <- regexpr("(?<![A-Za-z0-9_.])model[[:space:]]*\\{", masked, perl = TRUE)[1L]
  if(model_at < data_at){
    stop("The BUGS data block must precede the model block.", call. = FALSE)
  }
  if(!grepl("~", substr(masked, data_at, model_at - 1L), fixed = TRUE)){
    return(list(model = model, data = data))
  }
  prefix <- substr(model, 1L, model_at - 1L)
  inputs <- .bt_sampler_tokens(prefix)
  anchor <- utils::tail(make.unique(c(inputs, "BayesTools_data_anchor")), 1L)
  connection <- textConnection(paste0(prefix, "\nmodel { ", anchor, " <- 0 }"))
  on.exit(close(connection))
  rjags::jags.seed(seed)
  generator <- rjags::jags.model(connection, data = data[names(data) %in% inputs],
    n.chains = 1L, n.threads = 1L, n.adapt = 0L, quiet = TRUE)
  generated <- generator$data()
  data[names(generated)] <- generated
  sampling_model <- paste0(substr(model, 1L, data_at - 1L), substring(model, model_at))
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
