# Run the chains assigned to one R process. The caller owns process-level
# setup and reporting. Keep JAGS to one thread per process for both serial
# fits and PSOCK workers, including workers that own multiple chains.
.bt_rjags_group <- function(inits, model, data, monitor, adapt, burnin, sample,
                            thin, silent, keep_live, live_model = NULL, probability_plan = NULL){

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
  assignments <- rjags::list.samplers(live_model)
  samplers <- data.frame(
    Index.No = rep(seq_along(assignments), lengths(assignments)),
    Sampler = rep(as.character(names(assignments)), lengths(assignments)),
    Node = as.character(unlist(assignments, use.names = FALSE))
  )

  readouts <- if(!is.null(probability_plan)){
    .bt_rjags_finite_readouts(probability_plan, samplers)
  }
  draws <- .bt_rjags_samples(live_model, monitor, sample, thin, progress, readouts)
  out <- list(mcmc = draws, end_state = live_model$state(internal = TRUE), samplers = samplers,
    live_model = if(keep_live) live_model)
  if(!is.null(probability_plan)){
    probabilities <- stats::setNames(vector("list", nrow(probability_plan)), probability_plan$node)
    if(length(readouts)){
      probabilities[names(readouts)] <- attr(draws, "sampler_probabilities", exact = TRUE)
    }
    attr(out$mcmc, "sampler_probabilities") <- NULL
    out$sampler_probabilities <- probabilities
  }
  out
}

# JAGS reports a length-one array by its bare name, even when the supplied
# reference selects its only cell explicitly. The monitor factory checks
# the actual marker and sampler identity when it opens the trace.
.bt_rjags_finite_readouts <- function(plan, samplers){

  assigned <- samplers$Node[samplers$Sampler == "BayesTools::Finite"]
  scalar <- sub("\\[1(?:,1)*\\]$", "", plan$node, perl = TRUE)
  owned <- plan$node %in% assigned | scalar %in% assigned
  stats::setNames(plan$readout[owned], plan$node[owned])
}

.bt_rjags_probability_stop <- function(...){

  stop(errorCondition(paste0(...), class = "BayesTools_sampler_probabilities_error"))
}

.bt_rjags_probability_draws <- function(draws, readouts){

  columns <- colnames(draws[[1L]])
  by_name <- split(seq_along(columns), .bt_jags_parameter_base(columns))
  stats::setNames(lapply(seq_along(readouts), function(i){
    indices <- by_name[[readouts[[i]]]]
    width <- length(indices)
    if(width < 4L || width %% 2L != 0L ||
       !identical(columns[indices], paste0(readouts[[i]], "[", seq_len(width), "]"))){
      .bt_rjags_probability_stop("The compiled finite sampler returned an incomplete probability trace for '",
        names(readouts)[[i]], "'.")
    }
    lower <- draws[[1L]][1L, indices[[width]]]
    if(!is.finite(lower) || lower != floor(lower) ||
       any(vapply(draws, function(chain){
         support <- chain[, indices[[width]]]
         anyNA(support) || any(support != lower)
       }, logical(1)))){
      .bt_rjags_probability_stop("The compiled finite sampler returned inconsistent support labels for '",
        names(readouts)[[i]], "'.")
    }
    labels <- lower + seq_len((width - 2L) / 2L) - 1L
    coda::mcmc.list(lapply(draws, function(chain){
      values <- chain[, indices, drop = FALSE]
      if(anyNA(values) || !all(values[, width - 1L] %in% labels)){
        .bt_rjags_probability_stop("The compiled finite sampler returned invalid probability trace rows for '",
          names(readouts)[[i]], "'.")
      }
      colnames(values) <- c(paste0("p[", labels, "]"), paste0("log_p[", labels, "]"), "selected", "lower")
      values
    }))
  }), names(readouts))
}

# rjags 5-5 treats names containing an internal "_word_" like its reserved
# monitors and passes their subscript to JAGS as part of the variable name.
# Read those arrays by name, retaining just the coordinates the user requested.
# Other indexed requests can still read a small slice of a large array.
.bt_rjags_samples <- function(model, monitor, sample, thin, progress, readouts = NULL){

  monitor <- unname(trimws(monitor))
  base <- trimws(.bt_jags_parameter_base(monitor))
  indexed <- base != monitor
  whole <- unique(c(base[!indexed], base[indexed & grepl("_[A-Za-z][A-Za-z0-9]*_", base)]))
  native <- ifelse(base %in% whole, base, monitor)
  draws <- rjags::coda.samples(model, variable.names = unique(c(unname(native), unname(readouts))),
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
  if(!length(readouts) && identical(names(selection), columns)) return(draws)
  scientific <- coda::mcmc.list(lapply(draws, function(chain){
    chain <- chain[, unname(selection), drop = FALSE]
    colnames(chain) <- names(selection)
    chain
  }))
  if(length(readouts)){
    attr(scientific, "sampler_probabilities") <- .bt_rjags_probability_draws(draws, readouts)
  }
  scientific
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
                          thin, cl, silent, live_model = NULL, probability_plan = NULL){

  if(is.null(cl)){
    return(.bt_rjags_group(inits, model, data, monitor, adapt, burnin, sample,
      thin, silent, keep_live = TRUE, live_model = live_model, probability_plan = probability_plan))
  }

  # The runtime callbacks use the same round-robin chain assignment. Restore
  # the original chain order before returning draws and their matching states.
  groups <- split(seq_along(inits), rep_len(seq_along(cl), length(inits)))
  results <- parallel::clusterApply(cl, lapply(groups, function(ids) inits[ids]),
    .bt_rjags_worker, model = model, data = data, monitor = monitor, adapt = adapt,
    burnin = burnin, sample = sample, thin = thin, silent = silent, probability_plan = probability_plan)
  for(result in results){
    for(w in result$warnings) warning(w)
    if(inherits(result$value, "error")) stop(result$value)
  }
  out <- lapply(results, "[[", "value")
  chain_order <- order(unlist(groups, use.names = FALSE))
  mcmc <- unlist(lapply(out, function(group) unclass(group$mcmc)),
    recursive = FALSE, use.names = FALSE)
  end_state <- unlist(lapply(out, "[[", "end_state"), recursive = FALSE, use.names = FALSE)
  merged <- list(mcmc = coda::mcmc.list(mcmc[chain_order]), end_state = end_state[chain_order],
    samplers = out[[1L]]$samplers, live_model = NULL)
  if(!is.null(probability_plan)){
    if(!all(vapply(out, function(worker) identical(worker$samplers, merged$samplers), logical(1)))){
      .bt_rjags_probability_stop("Parallel JAGS workers selected different sampler assignments for the probability trace.")
    }
    merged$sampler_probabilities <- stats::setNames(lapply(probability_plan$node, function(node){
      values <- lapply(out, function(worker) worker$sampler_probabilities[[node]])
      if(all(vapply(values, is.null, logical(1)))) return(NULL)
      if(any(vapply(values, is.null, logical(1)))){
        .bt_rjags_probability_stop("Parallel JAGS workers returned incompatible probability traces for '", node, "'.")
      }
      chains <- unlist(lapply(values, unclass), recursive = FALSE, use.names = FALSE)
      coda::mcmc.list(chains[chain_order])
    }), probability_plan$node)
  }
  merged
}

.bt_rjags_fit <- function(model, data, inits, monitor, adapt, burnin, sample, thin,
                          cl = NULL, silent = TRUE, probability_plan = NULL){

  out <- .bt_rjags_run(model, data, inits, monitor, adapt, burnin, sample, thin, cl, silent,
    probability_plan = probability_plan)
  # Public chain indices count the requested adaptation and burn-in, even
  # when JAGS selects only samplers that never require adaptation.
  public_iteration <- function(draws) coda::mcmc.list(lapply(draws, function(chain){
    coda::mcmc(as.matrix(chain), start = as.double(adapt) + burnin + 1, thin = thin)
  }))
  out$mcmc <- public_iteration(out$mcmc)
  if(!is.null(probability_plan)){
    attr(out, "sampler_declarations") <- list(probabilities = probability_plan)
    out$sampler_probabilities <- lapply(out$sampler_probabilities, function(draws){
      if(!is.null(draws)) public_iteration(draws)
    })
  }
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

  probability_plan <- attr(fit, "sampler_declarations", exact = TRUE)[["probabilities"]]
  out <- .bt_rjags_run(fit$model, fit$data, fit$end_state, fit$monitor, adapt,
    0L, sample, fit$thin, cl, silent, live_model = fit$live_model, probability_plan = probability_plan)
  append_draws <- function(previous, additional, mismatch){
    coda::mcmc.list(Map(function(kept, fresh){
      if(!identical(colnames(kept), colnames(fresh))) mismatch()
      coda::mcmc(rbind(as.matrix(kept), as.matrix(fresh)), start = stats::start(kept), thin = fit$thin)
    }, previous, additional))
  }
  fit$mcmc <- append_draws(fit$mcmc, out$mcmc, function(){
    stop("The rebuilt model returned different monitored coordinates.", call. = FALSE)
  })
  if(!is.null(probability_plan)){
    fit$sampler_probabilities <- stats::setNames(lapply(probability_plan$node, function(node){
      kept <- fit$sampler_probabilities[[node]]
      fresh <- out$sampler_probabilities[[node]]
      mismatch <- function(){
        .bt_rjags_probability_stop("The rebuilt JAGS model returned an incompatible probability trace for '", node, "'.")
      }
      if(is.null(kept) != is.null(fresh)) mismatch()
      if(is.null(kept)) return(NULL)
      append_draws(kept, fresh, mismatch)
    }), probability_plan$node)
  }
  fit$sample <- as.double(fit$sample) + sample
  fit$end_state <- out$end_state
  fit$samplers <- out$samplers
  fit$live_model <- out$live_model
  fit
}
