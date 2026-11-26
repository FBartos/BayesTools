# The JAGS libraries whose linear-algebra threads BayesTools controls: the
# Windows builds of JAGS 5 and of later major versions.
.bt_jags_libraries <- paste0("libjags-", 5:9, ".dll")

# The OpenMP thread count of the loaded JAGS library. Sets it to 'threads'
# and returns the previous count; with 'threads' NA only returns the count;
# NA when no library of 'libraries' is loaded or it has no thread controls,
# and always off Windows (src/r-jags-threads.c gives the reasons).
.bt_jags_threads <- function(threads = NA_integer_, libraries = .bt_jags_libraries){

  .Call("BayesTools_jags_threads", libraries, threads, PACKAGE = "BayesTools")
}

# Evaluates 'expr' with JAGS's linear algebra on one thread and restores the
# previous count on exit. The Windows build of JAGS 5 otherwise starts a team
# of threads for every larger BLAS or LAPACK call, which makes models with
# multivariate nodes up to a thousand times slower; the draws were identical
# in every comparison.
# Used where JAGS compiles, adapts or samples in a process: the steps of the
# chain groups (serial and on PSOCK workers) and the data-block generator.
.bt_with_single_threaded_jags <- function(expr){

  # the count is read and its restoration registered before it is changed, so
  # that an interrupt in between cannot leave it changed
  previous <- .bt_jags_threads()
  if(!is.na(previous)){
    on.exit(.bt_jags_threads(previous), add = TRUE)
    .bt_jags_threads(1L)
  }
  expr
}

# Run the chains assigned to one R process. The caller owns process-level
# setup and reporting. Keep JAGS to one thread per process for both serial
# fits and PSOCK workers, including workers that own multiple chains.
# With 'prepare', the compiled model stays in 'session' for later calls;
# a positive 'block' then runs the last 'block * thin' iterations of burn-in
# with the monitors on and returns their 'block' draws (every 'thin'-th
# iteration) for the readiness check, or NULL for a model without samplers.
.bt_rjags_group <- function(inits, model, data, monitor, adapt, burnin, sample,
                            thin, silent, keep_live, live_model = NULL, probability_plan = NULL,
                            prepare = FALSE, block = 0L, session = .bt_rjags_worker_session){

  progress <- if(silent) "none" else "text"
  if(isTRUE(live_model)){
    # Only the coordinator's private worker cluster uses this marker.
    live_model <- session$model
    if(is.null(live_model)) stop("The worker lost its compiled JAGS model.", call. = FALSE)
  }
  if(is.null(live_model)){
    uses_lecuyer <- vapply(inits, function(chain){
      identical(chain[[".RNG.name"]], "lecuyer::RngStream")
    }, logical(1))
    if(any(uses_lecuyer)){
      rjags::load.module("lecuyer", quiet = TRUE)
    }
    if(!is.null(probability_plan)) .bt_rjags_prepare_finite_readout()
    connection <- textConnection(model)
    on.exit(close(connection))
    live_model <- rjags::jags.model(connection, data = data, inits = inits,
      n.chains = length(inits), n.threads = 1L, n.adapt = 0L, quiet = silent)
    adapted <- rjags::adapt(live_model, n.iter = adapt, end.adaptation = TRUE,
      progress.bar = progress)
    if(!adapted){
      warning(warningCondition(
        paste0("The JAGS adaptation phase was incomplete after ", adapt, " iterations."),
        class = c("BayesTools_adaptation_incomplete", "BayesTools_fit_diagnostic")
      ))
    }
  }
  # A model without samplers draws every monitored node afresh in each sweep:
  # its chains cannot disagree, and monitoring a block would change their
  # random-number stream. Its whole burn-in runs unmonitored, as with the
  # readiness check off.
  if(block > 0L && length(rjags::list.samplers(live_model)) == 0L){
    block <- 0L
  }
  # rjags refuses an update of zero iterations; the iterations of the block
  # are counted in double arithmetic
  unmonitored <- burnin - as.double(block) * thin
  if(unmonitored > 0){
    stats::update(live_model, n.iter = unmonitored, progress.bar = progress)
  }
  if(prepare){
    session$model <- live_model
    session$chunks <- list()
    if(block > 0L){
      return(.bt_rjags_burnin_block(block, thin, monitor, progress, session))
    }
    return(NULL)
  }
  .bt_rjags_batch(sample, monitor, thin, progress, keep_live, probability_plan,
    live_model = live_model, session = NULL)
}

# A block of 'draws' burn-in draws taken every 'thin' iterations ('draws *
# thin' iterations) with the fit's monitors on and without probability
# read-outs, for the readiness check. Its draws never enter the retained
# chunks. Monitoring does not change a sampler; JAGS forward-samples a
# monitored node that is otherwise not sampled (one without observed
# descendants) from the chain's own generator.
.bt_rjags_burnin_block <- function(draws, thin, monitor, progress, session = .bt_rjags_worker_session){

  .bt_rjags_samples(session$model, monitor, draws, thin, progress)
}

# The burn-in blocks of the chain groups, as one mcmc.list in the original
# chain order.
.bt_rjags_merge_blocks <- function(blocks, groups){

  chains <- unlist(lapply(blocks, unclass), recursive = FALSE, use.names = FALSE)
  chains <- chains[order(unlist(groups, use.names = FALSE))]
  coda::mcmc.list(lapply(chains, function(chain) coda::mcmc(as.matrix(chain))))
}

# Drop a worker's compiled model and retained chunks.
.bt_rjags_release_worker <- function(session = .bt_rjags_worker_session){

  session$model  <- NULL
  session$chunks <- NULL
  NULL
}

# Each JAGS_fit() creates its own PSOCK workers and stops them on exit.
# Serial progress uses an environment local to the fitting call instead.
.bt_rjags_worker_session <- new.env(parent = emptyenv())

.bt_rjags_batch <- function(sample, monitor, thin, progress, keep_live, probability_plan,
                            last = TRUE, session = .bt_rjags_worker_session,
                            live_model = session$model, release_model = FALSE){

  readouts <- if(!is.null(probability_plan)){
    stats::setNames(probability_plan$readout, probability_plan$node)
  }
  draws <- .bt_rjags_samples(live_model, monitor, sample, thin, progress, readouts, session, last)
  if(!last) return(NULL)
  assignments <- rjags::list.samplers(live_model)
  samplers <- data.frame(
    Index.No = rep(seq_along(assignments), lengths(assignments)),
    Sampler = rep(as.character(names(assignments)), lengths(assignments)),
    Node = as.character(unlist(assignments, use.names = FALSE))
  )

  out <- list(mcmc = draws, end_state = live_model$state(internal = TRUE), samplers = samplers,
    live_model = if(keep_live) live_model)
  if(release_model) session$model <- NULL
  if(!is.null(probability_plan)){
    out$sampler_probabilities <- attr(draws, "sampler_probabilities", exact = TRUE)
    attr(out$mcmc, "sampler_probabilities") <- NULL
  }
  out
}

# Put the recorder at the stock finite factory's priority before compilation.
.bt_rjags_prepare_finite_readout <- function(){

  invisible(.Call("BayesTools_prepare_finite_readout", PACKAGE = "BayesTools"))
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
    # Ownership is determined from the compiled marker by the monitor factory.
    # An unowned marker returns one sentinel column on the same schedule.
    if(width == 1L && identical(columns[indices], readouts[[i]]) &&
       all(vapply(draws, function(chain) isTRUE(all(chain[, indices] == -1)), logical(1)))){
      return(NULL)
    }
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
.bt_rjags_samples <- function(model, monitor, sample, thin, progress, readouts = NULL,
                              session = NULL, last = TRUE){

  monitor <- unname(trimws(monitor))
  base <- trimws(.bt_jags_parameter_base(monitor))
  indexed <- base != monitor
  whole <- unique(c(base[!indexed], base[indexed & grepl("_[A-Za-z][A-Za-z0-9]*_", base)]))
  native <- ifelse(base %in% whole, base, monitor)
  draws <- rjags::coda.samples(model, variable.names = unique(c(unname(native), unname(readouts))),
    n.iter = as.double(sample) * thin, thin = thin, na.rm = FALSE, progress.bar = progress)
  if(!is.null(session)){
    session$chunks[[length(session$chunks) + 1L]] <- draws
    if(!last) return(NULL)
    # Combine once, before checking undefined coordinates and probability traces.
    if(length(session$chunks) > 1L){
      draws <- coda::mcmc.list(lapply(seq_along(draws), function(i){
        chains <- lapply(session$chunks, `[[`, i)
        coda::mcmc(do.call(rbind, chains), start = stats::start(chains[[1L]]), thin = thin)
      }))
    }
    session$chunks <- NULL
  }

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
.bt_rjags_worker <- function(input, step, ...){

  warnings <- list()
  value <- tryCatch(withCallingHandlers(
    .bt_with_single_threaded_jags(step(input, ...)),
    warning = function(w){
      warnings[[length(warnings) + 1L]] <<- w
      invokeRestart("muffleWarning")
    }), error = identity)
  list(value = value, warnings = warnings)
}

# 'readiness' (initial fits only; NULL otherwise) asks for the readiness check
# between burn-in and sampling: the draws 'block' of the checked block of
# burn-in (taken every 'thin' iterations), and the settings and coordinator
# callbacks of .bt_burnin_readiness(). Its record is returned as the
# 'burnin_check' element.
.bt_rjags_run <- function(model, data, inits, monitor, adapt, burnin, sample,
                          thin, cl, silent, live_model = NULL, probability_plan = NULL,
                          jasp_prefix = NULL, retain_workers = FALSE, readiness = NULL){

  # The runtime callbacks use the same round-robin chain assignment. Restore
  # the original chain order before returning draws and their matching states.
  local <- is.null(cl)
  groups <- if(local) list(seq_along(inits)) else split(seq_along(inits), rep_len(seq_along(cl), length(inits)))
  session <- new.env(parent = emptyenv())
  # Whether the workers can be contacted: not while a worker call is in
  # flight (an interrupt leaves its reply pending) or after it failed. A
  # failure of the call itself, as opposed to an error a worker returns, is a
  # transport failure: new initial values cannot repair the connection.
  reachable <- TRUE
  run_groups <- function(step, inputs, ...){
    if(local) return(list(.bt_with_single_threaded_jags(step(inputs[[1L]], ..., session = session))))
    reachable <<- FALSE
    results <- tryCatch(parallel::clusterApply(cl, inputs, .bt_rjags_worker, step, ...),
      error = function(e) stop(.JAGS_worker_connection_error(e)))
    reachable <<- TRUE
    for(result in results){
      for(w in result$warnings) warning(w)
      if(inherits(result$value, "error")) stop(result$value)
    }
    lapply(results, "[[", "value")
  }
  jasp <- !is.null(jasp_prefix)
  progress <- if(silent) "none" else "text"
  block <- if(is.null(readiness)) 0L else readiness$block
  # The coordinator sees the end of burn-in of every chain before sampling
  # starts: the processes keep their compiled models between the calls of a
  # staged run, as the JASP progress batches do.
  staged <- jasp || block > 0L
  warmed_up <- FALSE
  # The caller catches sampler and worker errors. Let UI exceptions escape
  # that handler with their original condition, including JASP IPC failures.
  notify <- function(callback, ...){
    withCallingHandlers(callback(...), error = function(e){
      tryInvokeRestart("BayesTools_jasp_progress_error", e, warmed_up)
    })
  }
  run_staged <- function(){
    if(jasp) notify(.JASP_progress_bar_start, 1L, paste0(jasp_prefix, "Adapting and burnin the model"))
    out <- run_groups(.bt_rjags_group, lapply(groups, function(ids) inits[ids]),
      model = model, data = data, monitor = monitor, adapt = adapt, burnin = burnin, sample = sample,
      thin = thin, silent = silent, keep_live = local, live_model = live_model,
      probability_plan = probability_plan, prepare = staged, block = block)
    record <- NULL
    # A process without samplers returns no block; every process compiles the
    # same graph, so the processes agree.
    checked <- !vapply(out, is.null, logical(1))
    if(block > 0L && any(checked) && !all(checked)){
      stop("Internal error: the JAGS processes disagree on whether the model has samplers.", call. = FALSE)
    }
    if(block > 0L && all(checked)){
      record <- .bt_burnin_readiness(
        first = .bt_rjags_merge_blocks(out, groups),
        run_block = function(draws){
          .bt_rjags_merge_blocks(run_groups(.bt_rjags_burnin_block, rep(list(draws), length(groups)),
            thin = thin, monitor = monitor, progress = progress), groups)
        },
        assess          = readiness$assess,
        burnin_max_Rhat = readiness$burnin_max_Rhat,
        block           = block,
        burnin_extend   = readiness$burnin_extend,
        thin            = thin,
        max_extensions  = readiness$max_extensions,
        autofit         = readiness$autofit,
        time_exceeded   = readiness$time_exceeded,
        extension_room  = readiness$extension_room
      )
    }
    if(jasp){
      warmed_up <<- TRUE
      notify(.JASP_progress_bar_tick)
      batches <- min(5L, sample)
      sizes <- sample %/% batches + (seq_len(batches) <= sample %% batches)
      notify(.JASP_progress_bar_start, batches, paste0(jasp_prefix, "Sampling the model"))
      for(i in seq_along(sizes)){
        out <- run_groups(.bt_rjags_batch, rep(list(sizes[[i]]), length(groups)), monitor = monitor,
          thin = thin, progress = progress, keep_live = local,
          probability_plan = probability_plan, last = i == batches,
          release_model = !local && !retain_workers && i == batches)
        notify(.JASP_progress_bar_tick)
      }
    }else if(staged){
      out <- run_groups(.bt_rjags_batch, rep(list(sample), length(groups)), monitor = monitor,
        thin = thin, progress = progress, keep_live = local,
        probability_plan = probability_plan, last = TRUE,
        release_model = !local && !retain_workers)
    }
    list(out = out, record = record)
  }
  # A worker keeps its compiled model between the calls of a staged run. On
  # every exit before the run is complete (an error, an interrupt, a JASP
  # cancellation), release it while the workers can be contacted, so that a
  # model whose run failed is never reused or kept by a restart.
  completed <- FALSE
  if(!local && staged){
    on.exit(if(!completed && reachable){
      try(parallel::clusterCall(cl, .bt_rjags_release_worker), silent = TRUE)
    }, add = TRUE)
  }
  run <- run_staged()
  out <- run$out
  if(local){
    out <- out[[1L]]
  }else{
    out <- .bt_rjags_merge_workers(out, groups, thin, probability_plan)
    if((jasp && retain_workers) || isTRUE(live_model)) out$live_model <- TRUE
  }
  out$burnin_check <- run$record
  completed <- TRUE
  out
}

.bt_rjags_merge_workers <- function(out, groups, thin, probability_plan){

  chain_order <- order(unlist(groups, use.names = FALSE))
  merge_draws <- function(draws){
    chains <- unlist(lapply(draws, unclass), recursive = FALSE, use.names = FALSE)
    chains <- chains[chain_order]
    # Workers can spend different numbers of iterations in adaptation.
    # Merge by retained draw order, as the public fit indexes these rows.
    if(length(unique(vapply(chains, stats::start, numeric(1)))) > 1L){
      chains <- lapply(chains, function(chain) coda::mcmc(as.matrix(chain), start = 1, thin = thin))
    }
    coda::mcmc.list(chains)
  }
  end_state <- unlist(lapply(out, "[[", "end_state"), recursive = FALSE, use.names = FALSE)
  merged <- list(mcmc = merge_draws(lapply(out, "[[", "mcmc")), end_state = end_state[chain_order],
    samplers = out[[1L]]$samplers, live_model = NULL)
  if(!all(vapply(out, function(worker) identical(worker$samplers, merged$samplers), logical(1)))){
    merged$samplers_by_worker <- Map(function(chains, worker){
      list(chains = chains, samplers = worker$samplers)
    }, groups, out)
  }
  if(!is.null(probability_plan)){
    merged$sampler_probabilities <- stats::setNames(lapply(probability_plan$node, function(node){
      values <- lapply(out, function(worker) worker$sampler_probabilities[[node]])
      if(all(vapply(values, is.null, logical(1)))) return(NULL)
      columns <- colnames(values[[1L]][[1L]])
      if(!all(vapply(values, function(draws){
        !is.null(draws) && identical(colnames(draws[[1L]]), columns)
      }, logical(1)))){
        .bt_rjags_probability_stop("Parallel JAGS workers returned incompatible probability traces for '", node, "'.")
      }
      merge_draws(values)
    }), probability_plan$node)
  }
  merged
}

.bt_rjags_fit <- function(model, data, inits, monitor, adapt, burnin, sample, thin,
                          cl = NULL, silent = TRUE, probability_plan = NULL, jasp_prefix = NULL,
                          retain_workers = FALSE, readiness = NULL){

  out <- .bt_rjags_run(model, data, inits, monitor, adapt, burnin, sample, thin, cl, silent,
    probability_plan = probability_plan, jasp_prefix = jasp_prefix, retain_workers = retain_workers,
    readiness = readiness)
  record <- out$burnin_check
  out$burnin_check <- NULL
  # The burn-in actually run: the requested one and the extensions of the
  # readiness check, in the type the burn-in was given (the range check of
  # the extensions keeps an integer sum within the integer range).
  if(!is.null(record) && record$extension > 0){
    burnin <- burnin + if(is.integer(burnin)) as.integer(record$extension) else record$extension
  }
  # Public chain indices count the requested adaptation and the burn-in run,
  # even when JAGS selects only samplers that never require adaptation.
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
  attr(out, "burnin_check") <- record
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
  generator <- .bt_with_single_threaded_jags(rjags::jags.model(connection,
    data = data[names(data) %in% inputs], n.chains = 1L, n.threads = 1L, n.adapt = 0L, quiet = TRUE))
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
  fit$samplers_by_worker <- out[["samplers_by_worker"]]
  fit$live_model <- out$live_model
  fit
}
