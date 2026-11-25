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
  targets  <- .bt_burnin_requested_states(prepared, monitor, prior_list)
  metadata <- targets$metadata
  selected <- .bt_convergence_select(metadata, monitor, check_indicators)
  rows <- which(metadata$parameter %in% selected & metadata$role != "structural")
  chain_lengths <- vapply(prepared$mcmc_samples_list, NROW, integer(1))
  chain_id <- rep.int(seq_along(chain_lengths), chain_lengths)
  values <- vapply(rows, function(row){
    .bt_burnin_split_Rhat(targets$samples[[row]], chain_id, length(chain_lengths))
  }, numeric(1))
  names(values) <- metadata$parameter[rows]
  values[!is.na(values)]
}

# The convergence targets of a block name an indicator's states by the states
# that the block visits (with two visited states, the occupancy carries the
# indicator's own name). A state-specific 'monitor' element such as
# "mu_indicator (state 3)" names a categorical indicator and one of its
# declared states whatever the block visits: a requested state without its
# own target is added as its occupancy in the block, which is constant (not
# assessable) when the block never visits it. A selector that names no
# categorical indicator or a state it cannot take is left to fail in the
# selection, as in the convergence check.
.bt_burnin_requested_states <- function(prepared, monitor, prior_list){

  targets   <- prepared$targets
  requested <- setdiff(unique(monitor), targets$metadata$parameter)
  if(length(requested) == 0L){
    return(targets)
  }
  draws <- do.call(rbind, prepared$mcmc_samples_list)
  declared <- .bt_burnin_indicator_states(prior_list)
  for(indicator in names(declared)){
    states <- declared[[indicator]]
    labels <- paste0(indicator, " (state ", states, ")")
    if(!any(labels %in% requested) || !indicator %in% colnames(draws) ||
       !any(targets$metadata$source == indicator & targets$metadata$role == "indicator")){
      next
    }
    for(i in which(labels %in% requested)){
      targets$samples[[length(targets$samples) + 1L]] <- as.numeric(draws[, indicator] == states[[i]])
      targets$metadata <- rbind(targets$metadata, data.frame(
        parameter = labels[[i]], source = indicator, role = "indicator", stringsAsFactors = FALSE))
    }
  }
  targets
}

# The declared states of the categorical model indicators, by indicator node:
# the component indicator of a mixture of three or more components (dcat over
# the components; the indicator of a mixture total of an ordered prior
# included). Binary indicators have no state-specific targets.
.bt_burnin_indicator_states <- function(prior_list){

  states <- list()
  add <- function(parameter, prior){
    if(is.prior.factor(prior) && is.prior.ordered(prior)){
      return(add(.prior_ordered_total_name(parameter), prior$total))
    }
    if(!is.prior.mixture(prior) || is.prior.spike_and_slab(prior) || length(prior) < 3L){
      return(NULL)
    }
    indicator <- if(inherits(prior, "prior.bias_mixture")) "bias_indicator" else paste0(parameter, "_indicator")
    states[[indicator]] <<- seq_along(prior)
    NULL
  }
  parameters <- names(prior_list)
  for(i in seq_along(parameters)){
    if(!is.na(parameters[[i]]) && nzchar(parameters[[i]])){
      add(parameters[[i]], prior_list[[i]])
    }
  }
  states
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
