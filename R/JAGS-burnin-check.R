# The readiness check between burn-in and sampling of JAGS_fit(): the chains
# start sampling once they agree over the last block of burn-in. The block runs
# with the fit's monitors on; its draws are checked and dropped.

# The split-chain point estimate of the potential scale reduction factor of one
# quantity over a block: each chain's block is split into its first and last
# half (the middle draw of an odd block is left out), and the Gelman-Rubin point
# estimate is computed over the 2 x chains half-chains. Degenerate blocks are
# decided before coda, which returns NaN for some of them: NA (not assessed)
# when all values are equal or a value is not finite, Inf when every
# half-chain is constant but they differ, and 1 when the half-chains vary and
# have equal means.
.bt_burnin_split_Rhat <- function(values, chain_id, n_chains){

  if(any(!is.finite(values)) || length(unique(values)) < 2L){
    return(NA_real_)
  }
  halves <- unlist(lapply(seq_len(n_chains), function(chain){
    draws <- values[chain_id == chain]
    half  <- length(draws) %/% 2L
    list(draws[seq_len(half)], draws[length(draws) - half + seq_len(half)])
  }), recursive = FALSE)
  if(all(vapply(halves, function(half) all(half == half[[1L]]), logical(1)))){
    return(Inf)
  }
  if(length(unique(vapply(halves, mean, numeric(1)))) == 1L){
    return(1)
  }
  estimate <- tryCatch(
    suppressWarnings(coda::gelman.diag(
      coda::mcmc.list(lapply(halves, coda::mcmc)),
      autoburnin = FALSE,
      multivariate = FALSE
    )[["psrf"]][1L, 1L]),
    error = function(e) NA_real_
  )
  if(is.nan(estimate)) NA_real_ else estimate
}

# The split-chain R-hat of every quantity that the fit's convergence check
# would assess (the same selection: 'monitor', 'check_indicators' and the
# convergence roles of the coordinates; indicators as state occupancies) in a
# block of draws (an mcmc.list of the chains). Quantities that cannot be
# assessed are left out.
.bt_burnin_block_Rhat <- function(block, coordinates, prior_list, monitor, check_indicators){

  block_fit <- structure(list(mcmc = block), class = "BayesTools_fit")
  attr(block_fit, "prior_list") <- prior_list
  prepared <- .bt_convergence_prepare(
    fit = block_fit,
    coordinates = coordinates,
    add_parameters = NULL,
    monitor = monitor
  )
  metadata <- prepared$targets$metadata
  selected <- .bt_convergence_select(metadata, monitor, check_indicators)
  rows <- which(metadata$parameter %in% selected & metadata$role != "structural")
  chain_lengths <- vapply(prepared$mcmc_samples_list, NROW, integer(1))
  chain_id <- rep.int(seq_along(chain_lengths), chain_lengths)
  values <- vapply(rows, function(row){
    .bt_burnin_split_Rhat(prepared$targets$samples[[row]], chain_id, length(chain_lengths))
  }, numeric(1))
  names(values) <- metadata$parameter[rows]
  values[!is.na(values)]
}

# The check and its extensions. 'first' holds the draws of the last 'block'
# sweeps of the requested burn-in; 'run_block(sweeps)' runs a further block
# with the monitors on and returns its draws; 'assess(draws)' returns the
# split-chain R-hat per assessable quantity. Without 'autofit' the chains are
# not extended. With it, blocks of 'burnin_extend' sweeps are added (each
# checked alone) until the chains are ready, 'max_extensions' blocks were
# added, or 'time_exceeded()' (read before a block is started). Returns the
# record stored in attr(fit, "burnin_check").
.bt_burnin_readiness <- function(first, run_block, assess, burnin_max_Rhat, block, burnin_extend,
                                 max_extensions, autofit, time_exceeded){

  blocks     <- list()
  draws      <- first
  sweeps     <- block
  extensions <- 0L
  repeat{
    values <- sort(assess(draws), decreasing = TRUE)
    blocks[[length(blocks) + 1L]] <- data.frame(
      block     = length(blocks),
      sweeps    = sweeps,
      Rhat      = if(length(values)) unname(values[[1L]]) else NA_real_,
      parameter = if(length(values)) names(values)[[1L]] else NA_character_,
      stringsAsFactors = FALSE
    )
    stop_reason <- if(length(values) == 0L){
      "not_assessable"
    }else if(values[[1L]] <= burnin_max_Rhat){
      "ready"
    }else if(!autofit){
      "no_autofit"
    }else if(extensions >= max_extensions){
      "max_extend"
    }else if(time_exceeded()){
      "max_time"
    }
    if(!is.null(stop_reason)){
      break
    }
    draws      <- run_block(burnin_extend)
    sweeps     <- as.integer(burnin_extend)
    extensions <- extensions + 1L
  }

  failing <- values[values > burnin_max_Rhat]
  list(
    burnin_max_Rhat = burnin_max_Rhat,
    block           = block,
    burnin_extend   = burnin_extend,
    blocks          = do.call(rbind, blocks),
    extension       = extensions * burnin_extend,
    status          = switch(stop_reason, ready = "ready", not_assessable = "not_assessable", "not_ready"),
    stop_reason     = stop_reason,
    failing         = data.frame(parameter = names(failing), Rhat = unname(failing), stringsAsFactors = FALSE)
  )
}

# The warning for chains that were not ready when sampling started, or NULL.
.bt_burnin_warning <- function(record){

  if(is.null(record) || !identical(record$status, "not_ready")){
    return(NULL)
  }
  failing <- record$failing
  listed  <- utils::head(failing, 5L)
  values  <- paste0(sprintf("%.2f", listed$Rhat), " for '", listed$parameter, "'", collapse = ", ")
  if(nrow(failing) > 5L){
    values <- paste0(values, " and ", nrow(failing) - 5L, " more")
  }
  extensions <- record$extension / record$burnin_extend
  last_block <- utils::tail(record$blocks$sweeps, 1L)
  paste0(
    "The chains had not reached agreement at the end of burn-in",
    if(extensions > 0) paste0(" after ", extensions, ngettext(extensions, " extension", " extensions"),
                              " of ", record$burnin_extend, " iterations"),
    ": the split-chain R-hat of the last ", last_block, " burn-in iterations was ", values,
    " (above 'burnin_max_Rhat' of ", record$burnin_max_Rhat, "). ",
    switch(record$stop_reason,
      no_autofit = "Increase 'burnin', or fit with 'autofit = TRUE'.",
      max_extend = "Increase 'burnin' or 'burnin_max_extend'.",
      max_time   = "Increase 'burnin' or 'max_time'."
    )
  )
}

# Signal the warning of chains that were not ready (whatever 'silent' is) and
# return its message for the fit's warnings.
.bt_burnin_warn <- function(record){

  message <- .bt_burnin_warning(record)
  if(!is.null(message)){
    warning(warningCondition(message,
      class = c("BayesTools_burnin_not_converged", "BayesTools_fit_diagnostic")))
  }
  message
}
